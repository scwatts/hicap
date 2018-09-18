import logging
import pathlib


from . import utility


BlastFormat = {'qseqid': str,
               'sseqid': str,
               'qlen': int,
               'slen': int,
               'qstart': int,
               'qend': int,
               'sstart': int,
               'send': int,
               'length': int,
               'evalue': float,
               'bitscore': float,
               'pident': float,
               'nident': int,
               'mismatch': int,
               'gaps': int}


class BlastResults:

    def __init__(self, *values):
        for (attr, attr_type), value in zip(BlastFormat.items(), values):
            setattr(self, attr, attr_type(value))
        self.region = None
        self.orf = None
        self.broken = False


def create_blast_database(input_fp, output_dir):
    logging.debug('Create BLAST database for %s', input_fp)
    output_fp = pathlib.Path(output_dir, input_fp.name)
    command = 'makeblastdb -in %s -out %s -dbtype nucl' % (input_fp, output_fp)
    utility.execute_command(command)
    return output_fp


def align(query_fp, blast_db_fp):
    logging.debug('Aligning %s to sequences in %s using BLAST', query_fp, blast_db_fp)
    command = 'blastn -task blastn -db %s -query %s -outfmt "6 %s"'
    result = utility.execute_command(command % (blast_db_fp, query_fp, ' '.join(BlastFormat)))
    return result.stdout


def parse_blast_stdout(blast_results):
    logging.debug('Parsing %s BLAST results', len(blast_results))
    line_token_gen = (line.split() for line in blast_results.split('\n'))
    return [BlastResults(*lts) for lts in line_token_gen if lts]
