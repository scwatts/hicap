import Bio.Alphabet
import Bio.Graphics.GenomeDiagram
import Bio.Seq
import Bio.SeqFeature
import Bio.SeqRecord


from . import locus


SEQ_PADDING = 1000


def create_genbank_record(locus_data, contig_sequences):
    # Create base records
    position_deltas, gb_records = create_base_records(contig_sequences)

    # Add ORF hits
    orf_hits = locus.get_all_orf_hits(locus_data.regions)
    for contig, contig_hits in locus.sort_hits_by_contig(orf_hits).items():
        add_region_hit_features(contig_hits, position_deltas, contig, gb_records)

    # Add blast hits
    blast_hits = locus.get_all_blast_hits(locus_data)
    for contig, contig_hits in locus.sort_hits_by_contig(blast_hits).items():
        add_region_hit_features(contig_hits, position_deltas, contig, gb_records)

    # Add IS1016 hits
    for contig, contig_hits in locus.sort_hits_by_contig(locus_data.is_hits).items():
        add_is_hit_features(contig_hits, position_deltas, contig, gb_records)

    # Add nearby ORFs
    orf_counter = 0
    for contig, orfs in locus.sort_orfs_by_contig(locus_data.nearby_orfs).items():
        orf_counter = add_misc_orf_features(orfs, position_deltas, orf_counter, contig, gb_records)

    # Sort features by location
    for contig in gb_records.keys():
        gb_records[contig].features = sorted(gb_records[contig].features, key=lambda f: f.location.start)
    return [record for record in gb_records.values()]


def create_base_records(contig_sequences):
    # For some programs (like EMBOSS seqret) to parse the output correctly, there must be a valid
    # entry between '^ORGANISM +\.$'. The 'COMMENT' field works for this purpose
    anno = {'comment': 'created by hicap'}
    position_deltas = dict()
    gb_records = dict()
    for contig, (position_delta, sequence) in contig_sequences.items():
        sequence_record = Bio.Seq.Seq(sequence, Bio.Alphabet.IUPAC.unambiguous_dna)
        gb_records[contig] = Bio.SeqRecord.SeqRecord(seq=sequence_record, name=contig, annotations=anno)
        position_deltas[contig] = position_delta
    return position_deltas, gb_records


def add_region_hit_features(contig_hits, position_deltas, contig, gb_records):
    position_delta = position_deltas[contig]
    for hit in sorted(contig_hits, key=lambda h: locus.get_hit_start(h)):
        # Get appropriate representation of gene name
        region = hit.region if hit.region else locus.get_gene_region(hit.sseqid)
        qualifiers = {'gene': hit.sseqid, 'note': 'region_%s' % region}
        if hit.broken:
            qualifiers['note'] += ';fragment'
        # Get object with info relating to input query sequence
        if getattr(hit, 'orf'):
            element = hit.orf
        elif getattr(hit, 'seq_section'):
            element = hit.seq_section
            qualifiers['note'] += ';no_orf'
        # Create feature record
        start, end = element.start, element.end
        feature = create_cds_feature(start, end, position_delta, element.strand, qualifiers)
        gb_records[contig].features.append(feature)


def add_is_hit_features(contig_hits, position_deltas, contig, gb_records):
    position_delta = position_deltas[contig]
    for hit in sorted(contig_hits, key=lambda h: locus.get_hit_start(h)):
        # Get appropriate representation of gene name
        # TODO: maybe change the use of the gene field here
        qualifiers = {'gene': 'IS1016', 'note': 'insertion_sequence'}
        # Create feature record
        start, end = hit.seq_section.start, hit.seq_section.end
        feature = create_cds_feature(start, end, position_delta, hit.seq_section.strand, qualifiers)
        gb_records[contig].features.append(feature)


def add_misc_orf_features(orfs, position_deltas, orf_counter, contig, gb_records):
    position_delta = position_deltas[contig]
    for orf in sorted(orfs, key=lambda o: o.start):
        orf_counter += 1
        qualifiers = {'gene': 'orf_%s' % orf_counter, 'note': 'misc_orf'}
        # Create feature record
        feature = create_cds_feature(orf.start, orf.end, position_delta, orf.strand, qualifiers)
        gb_records[contig].features.append(feature)
    return orf_counter


def add_locus_feature(gb_records):
    locus_features = 0
    for record in gb_records:
        # Skip records with no cap hits
        if not any('region' in feature.qualifiers['note'] for feature in record.features):
            continue
        # Create qualifiers
        locus_features += 1
        quals = {'note': 'locus_%s' % locus_features}
        # Location - only include features that are part of the locus
        flocs = (f.location for f in record.features if 'region' in f.qualifiers['note'])
        flocs = sorted(flocs, key=lambda l: l.start)
        start, end = flocs[0].start, flocs[-1].end
        floc = Bio.SeqFeature.FeatureLocation(start=start, end=end)
        # Feature
        feature = Bio.SeqFeature.SeqFeature(location=floc, type='misc_feature', qualifiers=quals)
        record.features.insert(0, feature)


def collect_contig_sequences(fasta, locus_data):
    # Sort all Orfs and SeqSections by contig
    hits = locus.get_all_hits(locus_data)
    contig_elements = dict()
    for hit in hits:
        # Get element
        if getattr(hit, 'orf'):
            element = hit.orf
        elif getattr(hit, 'seq_section'):
            element = hit.seq_section
        # Add
        try:
            contig_elements[element.contig].add(element)
        except:
            contig_elements[element.contig] = {element}

    # Add nearby ORF hits
    for orf in locus_data.nearby_orfs:
        try:
            contig_elements[orf.contig].add(orf)
        except:
            contig_elements[orf.contig] = {orf}

    contig_sequences = dict()
    for contig, elements in contig_elements.items():
        # Get the most left and most right element associated with a hit
        elements_sorted = sorted(elements, key=lambda e: e.start)
        element_start = None
        element_end = None
        for element in elements_sorted:
            if element in locus_data.nearby_orfs:
                continue
            element_end = element
            if not element_start:
                element_start = element

        # Apply sequencing padding - extend if we do not extend beyond nearby orfs
        start = element_start.start - SEQ_PADDING
        if start > elements_sorted[0].start:
            start = elements_sorted[0].start - 16
        start = max(start, 0)
        end = element_end.end + SEQ_PADDING
        if end < elements_sorted[-1].end:
            end = elements_sorted[-1].end + 16
        end = min(end, len(fasta[contig]))
        contig_sequences[contig] = (start, fasta[contig][start:end])
    return contig_sequences


def create_cds_feature(start, end, delta, strand, qualifiers):
    feature_start = start - delta if (start - delta) > 1 else 1
    feature_end = end - delta
    feature_loc = Bio.SeqFeature.FeatureLocation(start=feature_start-1, end=feature_end, strand=strand)
    return Bio.SeqFeature.SeqFeature(location=feature_loc, type='CDS', qualifiers=qualifiers)
