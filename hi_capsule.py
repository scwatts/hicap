#!/usr/bin/env python3
'''
Copyright 2018 Stephen Watts
https://github.com/scwatts/hi_capsule

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''


import argparse
import collections
import pathlib
import logging
import subprocess
import sys
import tempfile


BLAST_FORMAT = ['qseqid', 'sseqid', 'qlen', 'slen', 'qstart', 'qend', 'sstart', 'send',
                'length', 'evalue', 'bitscore', 'pident', 'nident', 'mismatch', 'gaps']
REGION_ONE = ('bexA', 'bexB', 'bexC', 'bexD')
REGION_TWO = ('type_a', 'type_b', 'type_c', 'type_d', 'type_e', 'type_f')
REGION_THREE = ('hcsA', 'hcsB')


BlastResults = collections.namedtuple('BlastResults', BLAST_FORMAT)

class Gene():

    def __init__(self, name, region):
        self.name = name
        self.region = region

        self.database_fp = None
        self.blast_database_fp = None
        self.blast_results = None
        self.complete_hits = None


def get_arguments():
    parser = argparse.ArgumentParser()
    # Inputs
    # TODO: need to work on generalising structure and naming of database
    parser.add_argument('--database_dir', required=True, type=pathlib.Path,
            help='Directory containing locus database')
    parser.add_argument('--query_fp', required=True, type=pathlib.Path,
            help='Input FASTA query')

    # Options
    parser.add_argument('--debug', action='store_const', dest='log_level', const=logging.DEBUG,
            default=logging.WARNING, help='Print debug messages')
    parser.add_argument('--log_fp', type=pathlib.Path,
            help='Record logging messages to file')

    args = parser.parse_args()

    # Glob for database files
    args.database_fps = list(args.database_dir.glob('*fasta'))
    return args


def main():
    # Initialise
    args = get_arguments()
    initialise_logging(args.log_level, args.log_fp)
    check_arguments(args)

    # Get instances of Gene for all flanking genes (region I and III) and associate databases
    # TODO: must standardise naming and aggregate genes by region in a more robust way. this
    # could be done by compiling into a multi FASTA ~ region
    genes = dict()
    for region_genes, region in zip((REGION_ONE, REGION_THREE), ('one', 'two')):
        for gene in region_genes:
            genes[gene] = Gene(gene, region)

    for database_fp in args.database_fps:
        database_gene = database_fp.stem
        if database_gene in REGION_TWO:
            continue
        try:
            genes[database_gene].database_fp = database_fp
        except KeyError:
            logging.error('Can\'t match database %s to flanking gene', database_fp)
            sys.exit(1)

    # Run
    # Verbose for clarity
    with tempfile.TemporaryDirectory() as dh:
        # TODO: q parallelise, use asyncio
        for gene in genes.values():
            # Create blast databases
            gene.blast_database_fp = create_blast_database(gene.database_fp, dh)
            # Align flanking region genes to query
            blast_stdout = blast_query(args.query_fp, gene.blast_database_fp)
            gene.blast_results = parse_blast_stdout(blast_stdout)

        # TODO: get complete gene hits
        # TODO: separate into regions
        # TODO: catch loci which cannot be complete due to discontiguous sequences
        # TODO: note missing flanking genes and attempt to locate incomplete matches

        # TODO: QN: separate region II into single genes rather than aligning entire sequence

def initialise_logging(log_level, log_file):
    log_handles = list()
    log_handles.append(logging.StreamHandler())
    # Ensure that parent dir of log file exists, otherwise raise error during check_arguments
    if log_file and log_file.exists():
        log_handles.append(logging.FileHandler(log_file, mode='w'))

    log_message_format = '%(asctime)s %(levelname)s: %(message)s'
    log_formatter= logging.Formatter(fmt=log_message_format, datefmt='%d/%m/%Y %H:%M:%S')
    logger = logging.getLogger()
    for log_handle in log_handles:
        log_handle.setFormatter(log_formatter)
        logger.addHandler(log_handle)

    logger.setLevel(log_level)


def check_arguments(args):
    # Files
    if args.log_fp:
        check_filepath_exists(args.log_fp.parent, 'Directory %s for log filepath does not exist')
    check_filepath_exists(args.database_dir, 'Database directory %s does not exist')
    check_filepath_exists(args.query_fp, 'Input query %s does not exist')
    if not args.database_fps:
        msg = 'Could not find any database files (.fasta extension) in %s.'
        logging.error(msg, args.database_dir)
        sys.exit(1)


def check_filepath_exists(filepath, message_format):
    if not filepath.exists():
        logging.error(message_format, filepath)
        sys.exit(1)


def execute_command(command):
    logging.debug(command)
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True,
                            encoding='utf-8')
    if result.returncode != 0:
        logging.critical('Failed to run command: %s', result.args)
        logging.critical('stdout: %s', result.stdout)
        logging.critical('stderr: %s', result.stderr)
        sys.exit(1)
    return result


def create_blast_database(input_fp, output_dir):
    output_fp = pathlib.Path(output_dir, input_fp.name)
    command = 'makeblastdb -in %s -out %s -dbtype nucl' % (input_fp, output_fp)
    execute_command(command)
    return output_fp


def blast_query(query_fp, blast_db_fp):
    command = 'blastn -task blastn -db %s -query %s -outfmt "6 %s"'
    result = execute_command(command % (blast_db_fp, query_fp, ' '.join(BLAST_FORMAT)))
    return result.stdout


def parse_blast_stdout(blast_results):
    line_token_gen = (line.split() for line in blast_results.split('\n'))
    return [BlastResults(*lts) for lts in line_token_gen if lts]


if __name__ == '__main__':
    main()
