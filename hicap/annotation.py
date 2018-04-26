import logging
import pathlib
import re
import tempfile


import Bio.Alphabet
import Bio.Seq
import Bio.SeqRecord
import Bio.SeqFeature


from . import utility


PRODIGAL_RE = re.compile(r'^>[0-9]+_([0-9]+)_([0-9]+)_([-\+])$')


class Orf():

    def __init__(self, start, end, strand):
        self.start = start
        self.end = end
        self.strand = strand

        self.hit = None


def discover_orfs(loci_data):
    logging.info('Discovering ORFs and matching with BLAST alignments')
    for i, locus_data in enumerate(loci_data.values(), 1):
        with tempfile.TemporaryDirectory() as dh:
            fasta_fp = pathlib.Path(dh, '%s.fasta' % id(locus_data))
            utility.write_fasta(id(locus_data), locus_data.sequence, fasta_fp)
            prodigal_stdout = annotate(fasta_fp)
            locus_data.orfs = process_prodigal_stdout(prodigal_stdout)
            logging.info('Found %s ORFs for locus %s', len(locus_data.orfs), i)

        # Most inefficiently match orfs with BLAST alignments
        for orf in locus_data.orfs:
            for hit in locus_data.hits:
                hit_range = range(*sorted((hit.qstart, hit.qend)))
                orf_start = orf.start + locus_data.sequence_offset
                orf_end = orf.end + locus_data.sequence_offset
                orf_range = range(orf_start, orf_end)
                # Hmmm...
                small_range, large_range = sorted((hit_range, orf_range), key=lambda k: len(k))
                if max(small_range) in large_range or min(small_range) in large_range:
                    orf.hit = hit
        logging.info('Matched %s ORFs for locus %s', len([o for o in locus_data.orfs if o.hit]), i)


def generate_genbank(loci_data, query_name):
    # Create genbank records
    logging.info('Creating genbank records')
    for locus_data in loci_data.values():
        locus_data.genbank = Bio.SeqRecord.SeqRecord(
                seq=Bio.Seq.Seq(locus_data.sequence, Bio.Alphabet.IUPAC.unambiguous_dna),
                id=str(id(locus_data)),
                name=query_name[:15])
        for orf in locus_data.orfs:
            if orf.hit:
                if 'type' not in orf.hit.sseqid:
                    qualifiers = {'gene': orf.hit.sseqid}
                else:
                    qualifiers = {'locus_tag': orf.hit.sseqid}
            else:
                qualifiers = {}
            feature_location = Bio.SeqFeature.FeatureLocation(start=orf.start, end=orf.end)
            feature = Bio.SeqFeature.SeqFeature(
                    location=feature_location,
                    type='CDS',
                    qualifiers=qualifiers)
            locus_data.genbank.features.append(feature)


def annotate(query_fp):
    logging.debug('Annotating %s using Prodigal', query_fp)
    command = 'prodigal -c -f sco -i %s -m -p meta'
    result = utility.execute_command(command % query_fp)
    return result.stdout


def process_prodigal_stdout(prodigal_results):
    logging.debug('Parsing %s Prodgial results', len(prodigal_results))
    results_line_gen = (line for line in prodigal_results.split('\n') if not line.startswith('#'))
    results_group_gen =(PRODIGAL_RE.match(line).groups() for line in results_line_gen if line)
    return [Orf(int(s), int(e), d) for s, e, d in results_group_gen]
