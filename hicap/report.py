import logging
import pathlib


import Bio.Alphabet
import Bio.Seq
import Bio.SeqRecord
import Bio.SeqFeature
import Bio.SeqIO
import Bio.Graphics.GenomeDiagram
import reportlab.lib.colors


from . import database
from . import locus
from . import utility


SEQ_PADDING = 1000
COLOURS_COMPLETE = {'one': '#7bcebe',
                    'two': '#f9958b',
                    'three': '#ffe589',
                   }
COLOURS_BROKEN= {'one': '#8fa8a3',
                  'two': '#d3a7a2',
                  'three': '#d8ceab',
                 }


class SummaryData:

    def __init__(self):
        self.completeness = dict.fromkeys(database.SCHEME, tuple())
        self.duplciated = None
        self.multiple_contigs = None
        self.serotypes = None
        self.hits_by_contig = None


def create_report(region_groups, fasta_fp, prefix, output_dir):
    # Report
    output_report_fp = pathlib.Path(output_dir, '%s.tsv' % prefix)
    summary_data = create_summary(region_groups)
    write_summary(summary_data, output_report_fp)

    # Genbank
    output_gbk_fp = pathlib.Path(output_dir, '%s.gbk' % prefix)
    genbank_data = create_genbank_record(region_groups, fasta_fp)
    with output_gbk_fp.open('w') as fh:
        Bio.SeqIO.write(genbank_data, fh, 'genbank')

    # Graphic
    # TODO: clean up graphic. force addition of contig names to short tracks
    # TODO: legend
    output_pdf_fp = pathlib.Path(output_dir, '%s.pdf' % prefix)
    graphic_data = create_graphic(genbank_data, prefix)
    graphic_data.write(str(output_pdf_fp), 'PDF')


def create_summary(region_groups):
    summary_data = SummaryData()
    for region in ('one', 'two', 'three'):
        group = region_groups[region]
        # Completeness and good genes
        genes_found = {hit.sseqid for hit in group.hits if not hit.broken}
        if region in ('one', 'three'):
            genes_expected = database.SCHEME[region]
        else:
            genes_expected = set()
            for serotype in group.serotypes:
                genes_expected.update(database.SEROTYPES[serotype])
        genes_missing = genes_expected - genes_found
        summary_data.completeness[region] = (genes_missing, len(genes_found), len(genes_expected))
        # Duplication. Verbose for clarity
        if region in ('one', 'three') and is_duplicated(group.hits, 8000):
            summary_data.duplicated = True

    # Serotype
    summary_data.serotypes = region_groups['two'].serotypes

    # On multiple contigs
    hits_all_gen = (hit for group in region_groups.values() for hit in group.hits)
    contig_hits = locus.sort_hits_by_contig(hits_all_gen)
    if len(contig_hits) > 1:
        summary_data.mutliple_contigs = True

    # Store hits sorted by contig
    summary_data.hits_by_contig = contig_hits
    return summary_data


def write_summary(data, fp):
    with fp.open('w') as fh:
        print('#', ','.join(data.serotypes), sep='', file=fh)
        if any(region_complete[0] for region_complete in data.completeness.values()):
            print('#incomplete', file=fh)
        else:
            print('#complete', file=fh)
        for contig, contig_hits in data.hits_by_contig.items():
            hits_sorted = sorted(contig_hits, key=lambda k: k.orf.start)
            gene_names = get_gene_names(hits_sorted)
            print(contig, ','.join(gene_names), sep='\t', file=fh)


def get_gene_names(hits):
    gene_names = list()
    for hit in hits:
        if hit.broken:
            gene_names.append(hit.sseqid + '*')
        else:
            gene_names.append(hit.sseqid)
    return gene_names


def is_duplicated(hits, distance):
    hits_sorted = sorted(hits, key=lambda k: k.orf.start)
    hits_start = min(hits_sorted[0].orf.start, hits_sorted[0].orf.end)
    hits_end = max(hits_sorted[-1].orf.start, hits_sorted[-1].orf.end)
    return (hits_end - hits_start) >= distance


def create_genbank_record(region_groups, fasta_fp):
    logging.info('Creating genbank records')
    genbank_records = list()
    fasta = utility.read_fasta(fasta_fp)
    hit_gen = (hit for group in region_groups.values() for hit in group.hits)
    contig_hits = locus.sort_hits_by_contig(hit_gen)
    for i, (contig, contig_hits) in enumerate(contig_hits.items(), 1):
        position_delta, block_sequence = get_block_sequence(contig_hits, fasta[contig], SEQ_PADDING)
        block_sequence = Bio.Seq.Seq(block_sequence, Bio.Alphabet.IUPAC.unambiguous_dna)
        block_genbank = Bio.SeqRecord.SeqRecord(seq=block_sequence, name='locus_part_%s' % i, id=fasta_fp.stem[:15])
        for hit in sorted(contig_hits, key=lambda k: k.orf.start):
            # Get appropriate representation of gene name
            region = hit.region if hit.region else locus.get_gene_region(hit.sseqid)
            qualifiers = {'gene': hit.sseqid, 'region': region}
            if hit.broken:
                qualifiers['note'] = 'fragment'
            feature_start = hit.orf.start - position_delta
            feature_end = hit.orf.end - position_delta
            feature_loc = Bio.SeqFeature.FeatureLocation(start=feature_start, end=feature_end, strand=hit.orf.strand)
            feature = Bio.SeqFeature.SeqFeature(location=feature_loc, type='CDS', qualifiers=qualifiers)
            block_genbank.features.append(feature)
        genbank_records.append(block_genbank)
    return genbank_records


def get_block_sequence(hits, block_contig, margin):
    '''Extract sequence for loci with margin on either side'''
    hits_sorted = sorted(hits, key=lambda k: k.orf.start)
    orf_start = min(hits_sorted[0].orf.start, hits_sorted[0].orf.end)
    orf_end = min(hits_sorted[-1].orf.start, hits_sorted[-1].orf.end)

    # Calculate start end position for sequence with margin
    block_start = orf_start - margin if orf_start >= margin else 0
    block_end = orf_end + margin if (orf_end + margin) <= len(block_contig) else len(block_contig)

    # Return delta in position and sequence with margin
    return block_start, block_contig[block_start:block_end]


def create_graphic(records, prefix):
    graphic = Bio.Graphics.GenomeDiagram.Diagram(prefix)
    for record in records:
        feature_track = graphic.new_track(1, name=record.name, greytrack=True, start=0, end=len(record))
        track_features = feature_track.new_set()
        for feature in record.features:
            if feature.type != 'CDS':
                continue
            # Accept quals as list or single item for interop
            region = get_qualifier(feature.qualifiers['region'])
            gene_name = get_qualifier(feature.qualifiers['gene'])
            if 'note' in feature.qualifiers and get_qualifier(feature.qualifiers['note']) == 'fragment':
                gene_colour = COLOURS_BROKEN[region]
            else:
                gene_colour = COLOURS_COMPLETE[region]
            track_features.add_feature(feature, sigil='BIGARROW', label=True, name=gene_name,
                                       border=reportlab.lib.colors.black, label_position='start',
                                       label_angle=0, color=gene_colour)

        # TODO: check that height scaling is appropriate
        height = len(records) * 100
        end = max(len(record) for record in records)
        graphic.draw(format='linear', pagesize=(1800, height), fragments=1, start=0, end=end)
        return graphic


def get_qualifier(qualifier):
    return qualifier[0] if isinstance(qualifier, list) else qualifier
