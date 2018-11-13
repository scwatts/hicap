import logging


import Bio.Alphabet
import Bio.Graphics.GenomeDiagram
import Bio.Seq
import Bio.SeqFeature
import Bio.SeqRecord


from . import locus


SEQ_PADDING = 1000


def create_genbank_record(hits_all, nearby_orfs, contig_sequences):
    # Create base records
    logging.info('Creating genbank records')
    position_deltas, gb_records = create_base_records(contig_sequences)

    # Add hits and misc_feature for locus - SeqRecord instances modified inplace within function
    for contig, contig_hits in locus.sort_hits_by_contig(hits_all).items():
        add_hit_features(contig_hits, position_deltas, contig, gb_records)

    # Add ORFs - SeqRecord instances modified inplace within function
    orf_counter = 0
    for contig, orfs in locus.sort_orfs_by_contig(nearby_orfs).items():
        orf_counter = add_orf_features(orfs, position_deltas, orf_counter, contig, gb_records)

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


def add_hit_features(contig_hits, position_deltas, contig, gb_records):
    position_delta = position_deltas[contig]
    for hit in sorted(contig_hits, key=lambda k: k.orf.start):
        # Get appropriate representation of gene name
        region = hit.region if hit.region else locus.get_gene_region(hit.sseqid)
        qualifiers = {'gene': hit.sseqid, 'note': 'region_%s' % region}
        if hit.broken:
            qualifiers['note'] += ';fragment'
        # Create feature record
        start, end = hit.orf.start, hit.orf.end
        feature = create_cds_feature(start, end, position_delta, hit.orf.strand, qualifiers)
        gb_records[contig].features.append(feature)


def add_orf_features(orfs, position_deltas, orf_counter, contig, gb_records):
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
        if all(feature.qualifiers['gene'].startswith('orf_') for feature in record.features):
            continue
        # Create qualifiers
        locus_features += 1
        quals = {'note': 'locus_%s' % locus_features}
        # Location
        flocs = (f.location for f in record.features if not f.qualifiers['gene'].startswith('orf_'))
        flocs = sorted(flocs, key=lambda l: l.start)
        start, end = flocs[0].start, flocs[-1].end
        floc = Bio.SeqFeature.FeatureLocation(start=start, end=end)
        # Feature
        feature = Bio.SeqFeature.SeqFeature(location=floc, type='misc_feature', qualifiers=quals)
        record.features.insert(0, feature)


def collect_contig_sequences(fasta, hits, nearby_orfs):
    # Sort all ORFs by contig
    hit_orfs = {hit.orf for hit in hits}
    contig_orfs = locus.sort_orfs_by_contig(hit_orfs | nearby_orfs)
    contig_sequences = dict()
    for contig, orfs in contig_orfs.items():
        # Get the most left and most right ORF associated with a hit
        orfs_sorted = sorted(orfs, key=lambda o: o.start)
        orf_start = None
        orf_end = None
        for orf in orfs_sorted:
            if orf in nearby_orfs:
                continue
            orf_end = orf
            if not orf_start:
                orf_start = orf

        # Apply sequencing padding - extend if we do not extend beyond nearby orfs
        start = orf_start.start - SEQ_PADDING
        if start > orfs_sorted[0].start:
            start = orfs_sorted[0].start - 16
        start = max(start, 0)
        end = orf_end.end + SEQ_PADDING
        if end < orfs_sorted[-1].end:
            end = orfs_sorted[-1].end + 16
        end = min(end, len(fasta[contig]))
        contig_sequences[contig] = (start, fasta[contig][start:end])
    return contig_sequences


def create_cds_feature(start, end, delta, strand, qualifiers):
    feature_start = start - delta if (start - delta) > 1 else 1
    feature_end = end - delta
    feature_loc = Bio.SeqFeature.FeatureLocation(start=feature_start-1, end=feature_end, strand=strand)
    return Bio.SeqFeature.SeqFeature(location=feature_loc, type='CDS', qualifiers=qualifiers)
