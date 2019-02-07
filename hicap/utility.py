from distutils.version import LooseVersion
import logging
import re
import shutil
import subprocess
import sys


import Bio.SeqIO.FastaIO


from . import database
from . import locus


def initialise_logging(log_level, log_file):
    log_handles = list()
    log_handles.append(logging.StreamHandler())
    # Ensure that parent dir of log file exists, otherwise raise error during check_arguments
    if log_file and log_file.parent.exists():
        log_handles.append(logging.FileHandler(log_file, mode='w'))

    log_message_format = '%(asctime)s %(levelname)s: %(message)s'
    log_formatter = logging.Formatter(fmt=log_message_format, datefmt='%d/%m/%Y %H:%M:%S')
    logger = logging.getLogger()
    for log_handle in log_handles:
        log_handle.setFormatter(log_formatter)
        logger.addHandler(log_handle)

    logger.setLevel(log_level)


def check_filepath_exists(filepath, message_format):
    if not filepath.exists():
        logging.error(message_format, filepath)
        sys.exit(1)


def execute_command(command, check=True):
    logging.debug('command: %s', command)
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True,
                            encoding='utf-8')
    if check and result.returncode != 0:
        logging.critical('Failed to run command: %s', result.args)
        logging.critical('stdout: %s', result.stdout)
        logging.critical('stderr: %s', result.stderr)
        sys.exit(1)
    return result


def check_dependencies():
    logging.info('Checking dependencies')
    dependencies = {
            'blastn': {
                'vcommand': 'blastn -version',
                'vregex': re.compile(r'^blastn: (.+)\n'),
                'vrequired': '2.2.28'},
            'makeblastdb': {
                'vcommand': 'makeblastdb -version',
                'vregex': re.compile(r'^makeblastdb: (.+)\n'),
                'vrequired': '2.2.28'},
            'prodigal': {
                'vcommand': 'prodigal -v 2>&1',
                'vregex': re.compile(r'.*Prodigal V(.+?):'),
                'vrequired': '2.6.1'}
            }
    for dependency, version_data in dependencies.items():
        if not shutil.which(dependency):
            logging.critical('Could not find dependency %s' % dependency)
            sys.exit(1)
        result = execute_command(version_data['vcommand'], check=False)
        try:
            version = version_data['vregex'].search(result.stdout).group(1)
        except AttributeError:
            # TODO: should we have an option to skip dependency check?
            logging.critical('Unable to determine version for %s' % dependency)
            sys.exit(1)
        if LooseVersion(version) < LooseVersion(version_data['vrequired']):
            msg = '%s version %s or high is required'
            logging.critical(msg % (dependency, version_data['vrequired']))
            sys.exit(1)
        else:
            logging.debug('Found %s version %s' % (dependency, version))


def read_fasta(filepath):
    with filepath.open('r') as fh:
        # TODO: cleaner way to do this?
        fasta = {desc: seq for desc, seq in Bio.SeqIO.FastaIO.SimpleFastaParser(fh)}
    if not fasta:
        logging.error('Could not parse any valid FASTA records from %s', filepath)
        sys.exit(1)
    return fasta
