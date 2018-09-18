import logging
import pathlib


import Bio.Alphabet
import Bio.Seq
import Bio.SeqRecord
import Bio.SeqFeature


from . import locus
from . import utility


SEQ_PADDING = 1000


def create_report(region_groups, fasta_fp, prefix, output_dir):
    # TODO: create summary table
    # TODO: name broken genes with suffix of 'fragment_n'
    # Genbank
    output_gbk_fp = pathlib.Path(output_dir, '%s.gbk' % prefix)
    genbank_data = create_genbank_record(region_groups, fasta_fp)
    with output_gbk_fp.open('w') as fh:
        Bio.SeqIO.write(genbank_data, fh, 'genbank')


def create_genbank_record(region_groups, fasta_fp):
    logging.info('Creating genbank records')
    genbank_records = list()
    fasta = utility.read_fasta(fasta_fp)
    hit_gen = (hit for group in region_groups.values() for hit in group.hits)
    contig_hits = locus.sort_hits_by_contig(hit_gen)
    for i, (contig, contig_hits) in enumerate(contig_hits.items(), 1):
        position_delta, block_sequence = get_block_sequence(contig_hits, fasta[contig], SEQ_PADDING)
        block_genbank = Bio.SeqRecord.SeqRecord(
                seq=Bio.Seq.Seq(block_sequence, Bio.Alphabet.IUPAC.unambiguous_dna),
                name='locus_part_%s' % i,
                id=fasta_fp.stem[:15])
        for hit in sorted(contig_hits, key=lambda k: k.orf.start):
            # Get appropriate representation of gene name
            qualifiers = {'gene': hit.sseqid}
            feature_start = hit.orf.start - position_delta
            feature_end = hit.orf.end - position_delta
            feature_location = Bio.SeqFeature.FeatureLocation(start=feature_start, end=feature_end)
            feature = Bio.SeqFeature.SeqFeature(
                    location=feature_location,
                    type='CDS',
                    qualifiers=qualifiers)
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
