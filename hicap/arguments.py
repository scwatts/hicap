import argparse
import logging
import pathlib
import sys


from . import utility


def get_args():
    parser = argparse.ArgumentParser()
    # Inputs
    # TODO: need to work on generalising structure and naming of database
    parser.add_argument('--database_dir', required=True, type=pathlib.Path,
            help='Directory containing locus database')
    parser.add_argument('--query_fp', required=True, type=pathlib.Path,
            help='Input FASTA query')
    parser.add_argument('--output_dir', required=True, type=pathlib.Path,
            help='Output directory')

    # Options
    parser.add_argument('--gene_coverage', default=0.80, type=float,
            help='Minimum percentage coverage to consider a single gene complete. [default: 0.80]')
    parser.add_argument('--gene_identity', default=0.75, type=float,
            help='Minimum percentage identity to consider a single gene complete. [default: 0.75]')
    parser.add_argument('--type_coverage', default=0.90, type=float,
            help='Minimum percentage coverage to consider a locus type. [default: 0.90]')
    parser.add_argument('--type_identity', default=0.90, type=float,
            help='Minimum percentage identity to consider a locus type. [default: 0.90]')
    parser.add_argument('--broken_gene_length', default=60, type=int,
            help='Minimum length to consider a broken gene. [default: 60]')
    parser.add_argument('--broken_gene_identity', default=0.90, type=float,
            help='Minimum percentage identity to consider a broken gene. [default: 0.90]')
    parser.add_argument('--broken_type_length', default=400, type=int,
            help='Minimum length to consider a broken type. [default: 400]')
    parser.add_argument('--broken_type_identity', default=0.90, type=float,
            help='Minimum percentage identity to consider a broken type. [default: 0.90]')
    parser.add_argument('--debug', action='store_const', dest='log_level', const=logging.DEBUG,
            default=logging.INFO, help='Print debug messages')
    parser.add_argument('--log_fp', type=pathlib.Path,
            help='Record logging messages to file')

    args = parser.parse_args()

    # Glob for database files
    args.database_fps = list(args.database_dir.glob('*fasta'))
    return args


def check_args(args):
    # Files
    if args.log_fp:
        utility.check_filepath_exists(args.log_fp.parent, 'Directory %s for log filepath does not exist')
    utility.check_filepath_exists(args.database_dir, 'Database directory %s does not exist')
    utility.check_filepath_exists(args.query_fp, 'Input query %s does not exist')
    if not args.database_fps:
        msg = 'Could not find any database files (.fasta extension) in %s.'
        logging.error(msg, args.database_dir)
        sys.exit(1)

    # Directory
    if not args.output_dir.exists():
        logging.error('Output directory %s does not exist', args.output_dir)
        sys.exit(1)

    # Parameters
    if args.gene_coverage <= 0:
        logging.error('--gene_coverage must be greater than zero')
        sys.exit(1)
    if args.gene_coverage > 1:
        logging.error('--gene_coverage cannot be greater than 1.0')
        sys.exit(1)
