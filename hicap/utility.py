from distutils.version import LooseVersion
import logging
import re
import shutil
import subprocess
import sys


import Bio.Alphabet
import Bio.Seq
import Bio.SeqRecord
import Bio.SeqFeature
import Bio.SeqIO.FastaIO


from . import database


def initialise_logging(log_level, log_file):
    log_handles = list()
    log_handles.append(logging.StreamHandler())
    # Ensure that parent dir of log file exists, otherwise raise error during check_arguments
    if log_file and log_file.exists():
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
    logging.info('Collecting nucleotide sequence')
    with filepath.open('r') as fh:
        # TODO: cleaner way to do this?
        fasta = {desc: seq for desc, seq in Bio.SeqIO.FastaIO.SimpleFastaParser(fh)}
    if not fasta:
        logging.error('Could not parse any valid FASTA records from %s', filepath)
        sys.exit(1)
    return fasta


def create_genbank_record(loci_blocks, fasta_fp):
    '''Construct a Genbank record from loci_blocks'''
    logging.info('Creating genbank records')
    sequence_margin = 1000
    genbank_records = list()
    fasta = read_fasta(fasta_fp)
    for i, loci_block in enumerate(loci_blocks, 1):
        position_delta, loci_block_sequence = get_block_sequence(loci_block, fasta, sequence_margin)
        loci_genbank = Bio.SeqRecord.SeqRecord(
                seq=Bio.Seq.Seq(loci_block_sequence, Bio.Alphabet.IUPAC.unambiguous_dna),
                name='locus_part_%s' % i,
                id=fasta_fp.stem[:15])
        for orf in loci_block.orfs:
            # Get appropriate representation of gene name
            name, orf_hit = get_orf_hit_and_name(orf)
            qualifiers = {'gene': name}
            feature_start = orf.start - position_delta
            feature_end = orf.end - position_delta
            feature_location = Bio.SeqFeature.FeatureLocation(start=feature_start, end=feature_end)
            feature = Bio.SeqFeature.SeqFeature(
                    location=feature_location,
                    type='CDS',
                    qualifiers=qualifiers)
            loci_genbank.features.append(feature)
        genbank_records.append(loci_genbank)
    return genbank_records


def write_summary_data(loci_blocks, output_fp):
    serotypes = {s for lb in loci_blocks for s in lb.serotypes}
    with output_fp.open('w') as fh:
        print('#', ','.join(serotypes), sep='', file=fh)
        for loci_block in loci_blocks:
            genes = list()
            for orf in loci_block.orfs:
                # Get appropriate representation of gene name
                name, orf_hit = get_orf_hit_and_name(orf)
                genes.append(name)
            start = loci_block.orfs[0].start
            end = loci_block.orfs[-1].end
            print(loci_block.contig, start, end, ','.join(genes), sep='\t', file=fh)


def get_orf_hit_and_name(orf):
    '''Collect the first ORF hit and appropriate name'''
    # There should only ever be one hit per ORF here
    database_name, [orf_hit] = list(orf.hits.items())[0]
    # TODO: is there an appreciable difference if we construct a reverse hash map
    for region, databases in database.SCHEME.items():
        if database_name in databases:
            break
    name = orf_hit.sseqid if region == 'two' else database_name
    return name, orf_hit


def get_block_sequence(loci_block, contigs, margin):
    '''Extract sequence for loci with margin on either side'''
    loci_contig = contigs[loci_block.contig]
    orf_start = loci_block.orfs[0].start
    orf_end = loci_block.orfs[-1].end

    # Calculate start end position for sequence with margin
    loci_start = orf_start - margin if orf_start >= margin else 0
    loci_end = orf_end + margin if (orf_end + margin) <= len(loci_contig) else len(loci_contig)

    # Return delta in position and sequence with margin
    return loci_start, loci_contig[loci_start:loci_end]
