import logging
import re


from . import utility


PRODIGAL_RESULT_RE = re.compile(r'^>[0-9]+_([0-9]+)_([0-9]+)_([-+])$')
PRODIGAL_CONTIG_RE = re.compile(r'^# Sequence.+?seqhdr="(.+?)"(?:;|$)')


class Orf:
    '''Used to represent an ORF and is associated with a BlastResult instance'''

    def __init__(self, contig, start, end, strand):
        self.contig = contig
        self.start = int(start)
        self.end = int(end)
        self.strand = +1 if strand == '+' else -1
        self.sequence = str()


class SeqSection:
    '''Used to represent an section of query sequence for hits without ORFs'''

    def __init__(self, contig, start, end, strand):
        self.contig = contig
        self.start = int(start)
        self.end = int(end)
        self.strand = strand


def collect_orfs(fasta_fp, model_fp):
    logging.info('Collecting ORFs from input assembly')
    prodigal_stdout = annotate(fasta_fp, model_fp)
    orfs = process_prodigal_stdout(prodigal_stdout)

    fasta = utility.read_fasta(fasta_fp)
    for orf in orfs:
        orf.sequence = fasta[orf.contig][orf.start-1:orf.end]
    return orfs


def annotate(query_fp, model_fp):
    logging.debug('Annotating %s using Prodigal', query_fp)
    command = 'gzip -c -d %s | prodigal -f sco -i /dev/stdin -m -t %s' if query_fp.name.endswith('gz') else 'prodigal -f sco -i %s -m -t %s'
    result = utility.execute_command(command % (query_fp, model_fp))
    # Prodigal includes \r from FASTAs, causing problems with the output. Remove \r here
    return result.stdout.replace('\r', '')


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
