import argparse
import logging
import pathlib
import sys


from . import utility


def get_args():
    # Parsers
    parser_parent = argparse.ArgumentParser()
    parser_files = parser_parent.add_argument_group('File input and output')
    parser_params = parser_parent.add_argument_group('Search parameters')
    parser_other = parser_parent.add_argument_group('Other')

    # Inputs
    # TODO: need to work on generalising structure and naming of database
    database_dir = pathlib.Path(__file__).parent / 'database'
    parser_files.add_argument('--query_fp', required=True, type=pathlib.Path,
            help='Input FASTA query')
    parser_files.add_argument('--output_dir', required=True, type=pathlib.Path,
            help='Output directory')
    parser_files.add_argument('--database_dir', required=False, type=pathlib.Path, default=database_dir,
            help='Directory containing locus database. [default: %s]' % database_dir)

    # Parameters
    parser_params.add_argument('--gene_coverage', default=0.80, type=float,
            help='Minimum percentage coverage to consider a single gene complete. [default: 0.80]')
    parser_params.add_argument('--gene_identity', default=0.75, type=float,
            help='Minimum percentage identity to consider a single gene complete. [default: 0.75]')
    parser_params.add_argument('--broken_gene_length', default=60, type=int,
            help='Minimum length to consider a broken gene. [default: 60]')
    parser_params.add_argument('--broken_gene_identity', default=0.90, type=float,
            help='Minimum percentage identity to consider a broken gene. [default: 0.90]')

    # Other
    parser_other.add_argument('--debug', action='store_const', dest='log_level', const=logging.DEBUG,
            default=logging.INFO, help='Print debug messages')
    parser_other.add_argument('--log_fp', type=pathlib.Path,
            help='Record logging messages to file')
    parser_other.add_argument('--help_all', help='Display extended help')

    # Change behaviour of help display
    quick_help_args = ('--query_fp', '--output_dir', '--help_all')
    if '--help_all' in sys.argv[1:]:
        parser_parent.print_help()
        sys.exit(0)
    else:
        for arg in parser_parent._actions:
            if not any(qarg in arg.option_strings for qarg in quick_help_args):
                arg.help = argparse.SUPPRESS

    # Glob for database files
    args = parser_parent.parse_args()
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
