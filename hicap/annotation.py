import logging
import re


from . import utility


PRODIGAL_RESULT_RE = re.compile(r'^>[0-9]+_([0-9]+)_([0-9]+)_([-\+])$')
PRODIGAL_CONTIG_RE = re.compile(r'^# Sequence.+?seqhdr="(.+?)"(?:;|$)')


class Orf():

    def __init__(self, contig, start, end, strand):
        self.contig = contig
        self.start = int(start)
        self.end = int(end)
        self.strand = strand

        self.sequence = str()
        self.hits = dict()
        self.broken = False


def collect_orfs(fasta_fp):
    logging.info('Collecting ORFs from FASTA file')
    prodigal_stdout = annotate(fasta_fp)
    orfs = process_prodigal_stdout(prodigal_stdout)

    logging.info('Extracting nucleotide sequence of ORFs')
    fasta = utility.read_fasta(fasta_fp)
    for orf in orfs:
        orf.sequence = fasta[orf.contig][orf.start-1:orf.end]

    logging.info('Found %s ORFs', len(orfs))
    return orfs


def annotate(query_fp):
    logging.debug('Annotating %s using Prodigal', query_fp)
    command = 'prodigal -c -f sco -i %s -m -p meta'
    result = utility.execute_command(command % query_fp)
    return result.stdout


def process_prodigal_stdout(prodigal_results):
    logging.debug('Parsing %s Prodgial results', len(prodigal_results))
    orfs = list()
    contig = str()
    for line in prodigal_results.rstrip().split('\n'):
        if line.startswith('# Sequence Data'):
            contig = PRODIGAL_CONTIG_RE.match(line).group(1)
        elif line.startswith('# Model Data'):
            continue
        else:
            result = PRODIGAL_RESULT_RE.match(line).groups()
            orfs.append(Orf(contig, *result))
    return orfs
