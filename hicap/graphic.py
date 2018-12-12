import copy
import logging
import re
import tempfile
import xml.etree.ElementTree as ET


import Bio.Graphics.GenomeDiagram
import Bio.SeqFeature
import reportlab.lib.colors


COLOURS_COMPLETE = {
        'one':      '#7bcebe',
        'two':      '#f9958b',
        'three':    '#ffe589',
        'none':     '#d3d3d3'
                    }
COLOURS_BLAST = {
        'one':      '#b0e6db',
        'two':      '#ffcbbf',
        'three':    '#ffeea7'
                    }
COLOURS_BROKEN = {
        'one':      '#8fa8a3',
        'two':      '#d3a7a2',
        'three':    '#d8ceab'
        }

LABEL_SIZE = 12
TRACK_LABEL_SIZE = 18

# Groups show in brackets
# '(90) 56.775, (1709.957) 56.775, 1709.957 93.225, 90 93.225'
HPOINTS_RE = re.compile(r'^([0-9.]+)[^,]+, ([0-9.]+).+$')

# '90 56.775, 1709.957 56.775, 1709.957 93.225, 90 (93.225)'
VPOINTS_RE = re.compile(r'^.+ ([0-9.]+)$')

# 'M (1709.957101),56.775000 L (1709.957101),93.225000 Z'
PATH_RE = re.compile(r'^M ([0-9.]+).+?L ([0-9.]+).+Z$')

# Matches only inner brackets - outer are literal
# ' matrix((1.000000),0.000000,-0.000000,1.000000,(132.856235),(93.225000))'
LABEL_RE = re.compile(r'^ matrix\(([-10.]+).+? ?([0-9.]+), ?([0-9.]+)\)')

# Matches all instances of:
# '(172.493076 69.532500),' or '(172.493076 62.242500)$' (endline is not literal)
# in a points array for an arrow
SA_PATH_RE = re.compile(r'([0-9.]+ [0-9.]+)(?:,|$)')


TRANSFORM_TEMPLATE = ' matrix(1.000000, 0.000000, -0.000000, 1.000000, %s, %s)'


def prepare_genbank(records):
    # Check if we need to rotate any records
    for record in records:
        last_position = 0
        for i, feature in enumerate(record.features):
            if feature.location.start - last_position > 5000:
                # Must rotate
                break
            last_position = feature.location.end
        else:
            # If we don't break, then go to next record
            continue
        # Apply rotation
        rotate_locus(record, i)

    # Order records from longest to shortest
    return sorted(records, key=lambda r: len(r.seq), reverse=True)


def rotate_locus(record, index):
    features_lower = record.features[:index]
    features_upper = record.features[index:]
    upper_block_size = len(record.seq) - features_upper[0].location.start

    # Rotate features
    upper_block_offset = features_upper[0].location.start
    for feature in features_upper:
        feature.location._start = Bio.SeqFeature.ExactPosition(feature.location.start - upper_block_offset)
        feature.location._end = Bio.SeqFeature.ExactPosition(feature.location.end - upper_block_offset)
    for feature in features_lower:
        feature.location._start = Bio.SeqFeature.ExactPosition(feature.location.start + upper_block_size)
        feature.location._end = Bio.SeqFeature.ExactPosition(feature.location.end + upper_block_size)
    record.features = sorted(record.features, key=lambda f: f.location.start)

    # Rotate sequence and trim
    sequence = record.seq[upper_block_offset:] + record.seq[:upper_block_offset]
    record.seq = sequence[:record.features[-1].location.end]
    logging.warning('The contig "%s" has been rotated for the graphical output', record.name)


def create_graphic(records, prefix):
    graphic = Bio.Graphics.GenomeDiagram.Diagram(prefix)
    for record in records:
        feature_track = graphic.new_track(1, name=record.name, greytrack=True, start=0, end=len(record))
        track_features = feature_track.new_set()
        for feature in record.features:
            if feature.type != 'CDS':
                continue
            # TODO: can we clean this up a little?
            # Accept quals as list or single item for interop
            notes = process_notes(get_qualifier(feature.qualifiers['note']))
            if notes['region'] != 'none':
                gene_name = get_qualifier(feature.qualifiers['gene'])
                gene_border = reportlab.lib.colors.HexColor(0x000000) # black
                sigil = 'BIGARROW'
            else:
                gene_name = ''
                gene_border = reportlab.lib.colors.HexColor(0x808080) # gray
                sigil = 'ARROW'

            if notes['fragment']:
                # Truncated hit with ORF
                gene_colour = COLOURS_BROKEN[notes['region']]
            elif notes['no_orf']:
                # Blast hit without ORF
                gene_colour = COLOURS_BLAST[notes['region']]
                gene_border = reportlab.lib.colors.HexColor(0x00000000, hasAlpha=True) # remove gene border
            elif notes['is']:
                # IS hit
                gene_colour = reportlab.lib.colors.HexColor(0x7aa7cc) # dark(ish) pastel blue
                gene_border = reportlab.lib.colors.HexColor(0x000000) # black
            else:
                # Complete hit with ORF
                gene_colour = COLOURS_COMPLETE[notes['region']]
            label_position = 'start' if feature.strand == 1 else 'end'

            track_features.add_feature(feature, sigil=sigil, label=True, name=gene_name,
                                       border=gene_border, label_position=label_position,
                                       label_angle=0, label_size=LABEL_SIZE, color=gene_colour)

        # TODO: Scale width with size
        # TODO: add contig boundaries symbols
        height = len(records) * 150
        end = max(len(record) for record in records)
        graphic.draw(format='linear', pagesize=(1800, height), fragments=1, track_size=0.30, start=0, end=end)
    return graphic


def process_notes(note_str):
    notes = {'region': 'none', 'fragment': False, 'no_orf': False, 'is': False}
    for note_token in note_str.split(';'):
        if note_token.startswith('region_'):
            notes['region'] = note_token.replace('region_', '')
        elif note_token == 'fragment':
            notes['fragment'] = True
        elif note_token == 'no_orf':
            notes['no_orf'] = True
        elif note_token == 'insertion_sequence':
            notes['is'] = True
    return notes


def patch_graphic(graphic_data):
    # TODO: dodge labels for very short contigs
    svg_data = get_svg_data(graphic_data)
    svg_tree = ET.fromstring(svg_data)
    visual_parent = svg_tree.find('.//{http://www.w3.org/2000/svg}g[@transform=""]')

    # Get track bound sizes - must round to 3 significant digits for backwards compatibility
    track_style = 'stroke: rgb(96%,96%,96%); stroke-linecap: butt; stroke-width: 1; fill: rgb(96%,96%,96%);'
    track_backgrounds = svg_tree.findall('.//*[@style="%s"]' % track_style)
    track_hbounds = set()
    for track_background in track_backgrounds:
        bounds = HPOINTS_RE.match(track_background.get('points')).groups()
        bounds = tuple(round(float(bound), 3) for bound in bounds)
        track_hbounds.add(bounds)

    # Send the mid-lines backwards
    patch_track_midlines(visual_parent, track_hbounds, svg_tree)

    # Move small arrows to center of track
    patch_small_arrow_annotations(svg_tree)

    # Fix reversed, mirrored labels and pad other labels
    # Genes on the non-coding strand have their labels upside down which is difficult to read
    patch_element_labels(svg_tree)

    # Remove original track labels and add new ones
    patch_track_name_labels(visual_parent, track_hbounds, track_backgrounds, svg_tree, graphic_data)
    return ET.tostring(svg_tree, encoding='unicode')


def patch_track_midlines(visual_parent, track_hbounds, svg_tree):
    '''Move the track mid lines to the back (but in front of the gray background'''
    # Place them behind the gene symbols but in front track background shading
    paths = svg_tree.findall('.//{http://www.w3.org/2000/svg}g[@transform=""]/{http://www.w3.org/2000/svg}path')
    line_elements = list()
    for path in paths:
        # Round to 3 significant digits for backwards compatibility
        path_hbounds = tuple(round(float(bound), 3) for bound in PATH_RE.match(path.get('d')).groups())
        if path_hbounds in track_hbounds:
            # Remove mid lines and store for later insert at appropriate position
            line_elements.append(path)
            visual_parent.remove(path)
        elif path_hbounds[0] not in {p for hbounds in track_hbounds for p in hbounds}:
            # Remove x-axis ticks
            visual_parent.remove(path)

    # Insert track midpoint lines immediately after track background elements
    for line_element in line_elements:
        visual_parent.insert(len(track_hbounds), line_element)


def patch_small_arrow_annotations(svg_tree):
    '''Vertically center small arrow annotations on the track'''
    # Collect all elements to center
    element_styles = [
            'stroke: rgb(0%,0%,0%); stroke-linecap: round; stroke-width: 1; fill: rgb(47%,65%,80%);',
            'stroke: rgb(50%,50%,50%); stroke-linecap: round; stroke-width: 1; fill: rgb(82%,82%,82%);',
            ]
    elements = list()
    for style in element_styles:
        elements.extend(svg_tree.findall('.//*[@style="%s"]' % style))

    for element in elements:
        # Get arrow point coordinates
        coords_str = SA_PATH_RE.findall(element.get('points'))
        coord_gen = (coord.split(' ') for coord in coords_str)
        coords = [(float(x), float(y)) for x, y in coord_gen]
        # Calculate move distance - half the height of the arrow
        # Determine move direction by finding direction of arrow
        move_dist = (coords[2][1] - coords[4][1]) / 2
        if coords[0][0] > coords[1][0]:
            move_dist *= -1
        # Apply move distance and build points string
        points_strings = list()
        for x, y in coords:
            points_strings.append('%s %s' % (x, y-move_dist))
        points_string = ', '.join(points_strings)
        # Assign new points string
        element.attrib['points'] = points_string


def patch_element_labels(svg_tree):
    '''Mirror and move non-coding strand labels'''
    # Set the padding to be equal
    texts = svg_tree.findall('.//{http://www.w3.org/2000/svg}g[@transform=""]/{http://www.w3.org/2000/svg}g')
    for text in texts:
        a, x, y = LABEL_RE.match(text.get('transform')).groups()
        if float(a) == -1:
            y_adjusted = float(y) - LABEL_SIZE
            transform = TRANSFORM_TEMPLATE % (x, y_adjusted)
            text.set('transform', transform)
        else:
            y_adjusted = float(y) + 3
            transform = TRANSFORM_TEMPLATE % (x, y_adjusted)
            text.set('transform', transform)


def patch_track_name_labels(visual_parent, track_hbounds, track_backgrounds, svg_tree, graphic_data):
    # Remove original track name label - on short tracks, these don't even appear
    name_style = 'font-family: Helvetica; font-size: 8px; fill: rgb(60%,60%,60%);'
    names = svg_tree.findall('.//*[@style="%s"]..' % name_style)
    for name in names:
        visual_parent.remove(name)

    # Add contig names as track labels
    contig_names = [track.name for track in graphic_data.tracks.values()]
    track_vbounds = [VPOINTS_RE.match(tb.get('points')).group(1) for tb in track_backgrounds]

    # Get the correct position to insert these labels
    for insert_index, element in enumerate(visual_parent):
        if element.tag == '{http://www.w3.org/2000/svg}g':
            break

    # Create new elements and place them into the svg document
    text_format = 'font-family: Helvetica; font-size: %spx; fill: rgb(0%%,0%%,0%%);' % TRACK_LABEL_SIZE
    text_transform = 'translate(0,0) scale(1,-1)'
    text_attribs = {'style': text_format, 'transform': text_transform, 'x': '0', 'y': '0'}
    text_element = ET.Element('{http://www.w3.org/2000/svg}text', attrib=text_attribs)
    for contig_name, vbound, hbounds in zip(contig_names, track_vbounds, track_hbounds):
        x = hbounds[0]
        y = float(vbound) + TRACK_LABEL_SIZE
        name_group = ET.Element('{http://www.w3.org/2000/svg}g', attrib={'transform': TRANSFORM_TEMPLATE % (x, y)})
        name_text = copy.deepcopy(text_element)
        name_text.text = contig_name
        name_group.append(name_text)
        visual_parent.insert(insert_index, name_group)


def get_svg_data(graphic_data):
    # Writing to string fails under Python3, must write to disk
    with tempfile.TemporaryFile('w+') as fh:
        graphic_data.write(fh, 'SVG')
        fh.seek(0)
        return fh.read()


def get_qualifier(qualifier):
    return qualifier[0] if isinstance(qualifier, list) else qualifier
