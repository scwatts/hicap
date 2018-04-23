#!/usr/bin/env python3
'''
Copyright 2018 Stephen Watts
https://github.com/scwatts/hicap

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


# TODO: examine corner cases (GCA_002987225; bexD on a different contig)
# TODO: organise functions
# TODO: refactor main function into discrete units
# TODO: many more logging messages
# TODO: many more unit tests
# TODO: software version checking
# TODO: place all region II into a single flat database
# TODO: QN: separate region II into single genes rather than aligning entire sequence


import argparse
import collections
from distutils.version import LooseVersion
import pathlib
import logging
import re
import shutil
import statistics
import subprocess
import sys
import tempfile


from Bio.SeqIO.FastaIO import SimpleFastaParser


BLAST_FORMAT = ['qseqid', 'sseqid', 'qlen', 'slen', 'qstart', 'qend', 'sstart', 'send',
                'length', 'evalue', 'bitscore', 'pident', 'nident', 'mismatch', 'gaps']
REGION_ONE = ('bexA', 'bexB', 'bexC', 'bexD')
REGION_TWO = ('type_a', 'type_b', 'type_c', 'type_d', 'type_e', 'type_f')
REGION_THREE = ('hcsA', 'hcsB')
FLANK_ORDER = REGION_ONE + REGION_THREE


BlastResults = collections.namedtuple('BlastResults', BLAST_FORMAT)


class GeneData():

    def __init__(self, name, region):
        self.name = name
        self.region = region

        self.database_fp = pathlib.Path()
        self.blast_database_fp = pathlib.Path()
        self.blast_results = list()
        self.complete_hits = list()
        self.incomplete_hits = list()


class Locus():

    def __init__(self, identifier, contig=None, genes=None):
        self.identifier = identifier
        if contig:
            self.contigs = [contig]
        if genes:
            self.genes = {gene.sseqid: gene for gene in genes}

        self.type_hits = list()
        self.broken_genes = list()
        self._is_complete = bool()

    @property
    def is_complete(self):
        if self._is_complete is bool():
            self._is_complete = ( len(self.genes) + len(self.broken_genes) ) == len(FLANK_ORDER)
        return self._is_complete

    @property
    def discontiguous_sequence(self):
        return len(self.contigs) > 1

    @property
    def missing_genes(self):
        gene_name_gen = (g.sseqid for g in list(self.genes.values()) + self.broken_genes)
        return set(gene_name_gen) ^ set(FLANK_ORDER)


class RegionData():

    def __init__(self, name, database_fp):
        self.name = name
        self.database_fp = database_fp

        self.blast_database_fp = pathlib.Path()
        self.blast_results = list()


def get_arguments():
    parser = argparse.ArgumentParser()
    # Inputs
    # TODO: need to work on generalising structure and naming of database
    parser.add_argument('--database_dir', required=True, type=pathlib.Path,
            help='Directory containing locus database')
    parser.add_argument('--query_fp', required=True, type=pathlib.Path,
            help='Input FASTA query')

    # Options
    parser.add_argument('--gene_coverage', default=0.80, type=float,
            help='Minimum percentage coverage to consider a single gene complete. [default: 0.80]')
    parser.add_argument('--gene_identity', default=0.75, type=float,
            help='Minimum percentage identity to consider a single gene complete. [default: 0.75]')
    parser.add_argument('--type_coverage', default=0.90, type=float,
            help='Minimum percentage coverage to consider a locus type. [default: 0.90]')
    parser.add_argument('--type_identity', default=0.90, type=float,
            help='Minimum percentage identity to consider a locus type. [default: 0.90]')
    parser.add_argument('--debug', action='store_const', dest='log_level', const=logging.DEBUG,
            default=logging.INFO, help='Print debug messages')
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
    check_dependencies()
    check_arguments(args)
    genes_data, region_data = initialise_database(args.database_fps)

    # Run
    search_flanking_genes(genes_data, args.query_fp, args.gene_coverage)
    # TODO: CRIT: check for complete genes dropped here
    loci_data = get_best_loci(genes_data)
    consolidate_loci_contigs(loci_data)
    type_region_two(loci_data, region_data, args.query_fp, args.type_coverage, args.type_identity)
    broken_genes = find_broken_genes(loci_data, genes_data)
    write_results(loci_data, broken_genes)


def initialise_database(database_fps):
    # Associate databases with appropriate class and region
    # TODO: must standardise naming and aggregate genes by region in a more robust way. this
    # could be done by compiling into a multi FASTA ~ region
    genes_data = dict()
    region_data = dict()
    for region_genes_data, region in zip((REGION_ONE, REGION_THREE), ('one', 'two')):
        for gene in region_genes_data:
            genes_data[gene] = GeneData(gene, region)

    for database_fp in database_fps:
        database_name = database_fp.stem
        if database_name in REGION_TWO:
            region_data[database_name] = RegionData(database_name, database_fp)
            continue
        try:
            genes_data[database_name].database_fp = database_fp
        except KeyError:
            logging.error('Can\'t match database %s to flanking gene', database_fp)
            sys.exit(1)
    return genes_data, region_data


def read_query_fasta(query_fp):
    logging.info('Collecting query nucleotide sequence')
    with query_fp.open('r') as fh:
        # TODO: cleaner way to do this?
        fasta = {desc.split(' ')[0]: seq for desc, seq in SimpleFastaParser(fh)}
    if not fasta:
        logging.error('Could not parse any valid FASTA records from %s', query_fp)
        sys.exit(1)
    return fasta


def search_flanking_genes(genes_data, query_fp, gene_coverage):
    # Blast query against database and collect complete flanking gene hits
    # TODO: q parallelise, use asyncio
    logging.info('Searching for flanking genes')
    with tempfile.TemporaryDirectory() as dh:
        for gene in genes_data.values():
            gene.blast_database_fp = create_blast_database(gene.database_fp, dh)
            blast_stdout = blast_query(query_fp, gene.blast_database_fp)
            gene.blast_results = parse_blast_stdout(blast_stdout)
            # Set hits as complete or otherwise
            complete_hits, incomplete_hits = sort_flanking_hits(gene.blast_results, gene_coverage)
            gene.complete_hits = complete_hits
            gene.incomplete_hits = incomplete_hits


def type_region_two(loci_data, region_data, query_fp, type_coverage, type_identity):
    # Blast inferred region II against references
    logging.info('Extracting region II sequence for typing')
    query_fasta = read_query_fasta(query_fp)
    logging.info('Aligning region II to known serotype sequences')
    with tempfile.TemporaryDirectory() as dh:
        for i, locus_data in enumerate(loci_data, 1):
            # Write region II, if region II is unobtainable then skip
            # TODO: set flag if region II cannot be extracted
            seq = get_region_two_sequence(locus_data, query_fasta)
            if not seq:
                continue
            seq_fp = pathlib.Path(dh, '%s.fasta' % i)
            with seq_fp.open('w') as fh:
                print('>', i, sep='', file=fh)
                for line in [seq[i:i+80] for i in range(0, len(seq), 80)]:
                    print(line, file=fh)

            # Blast region II and get a filtered hit
            for region in region_data.values():
                region.blast_database_fp = create_blast_database(region.database_fp, dh)
                blast_stdout = blast_query(seq_fp, region.blast_database_fp)
                region.blast_results = parse_blast_stdout(blast_stdout)

                type_hit = filter_type_hits(region.blast_results, type_coverage, type_identity)
                if type_hit:
                    msg = 'Found good hit for %s in locus %s'
                    logging.info(msg, region.name, locus_data.identifier)
                    locus_data.type_hits.append((region.name, type_hit))


def find_broken_genes(loci_data, genes_data):
    # TODO: this strategy needs some attention
    # Attempt to find broken genes for incomplete loci on the same contig
    logging.info('Searching for broken genes near existing loci')
    for incomplete_locus in (l for l in loci_data if not l.is_complete):
        for missing_gene in incomplete_locus.missing_genes:
            broken_gene = find_broken_gene(genes_data[missing_gene].incomplete_hits, incomplete_locus)
            incomplete_locus.broken_genes.append(broken_gene)

    # Next try to find broken genes on different contigs
    logging.info('Collecting remaining broken genes from anywhere possible')
    missing_gene_counts = dict()
    for incomplete_locus in (l for l in loci_data if not l.is_complete):
        for missing_gene in incomplete_locus.missing_genes:
            try:
                missing_gene_counts[missing_gene] += 1
            except KeyError:
                missing_gene_counts[missing_gene] = 1

    broken_genes = dict()
    for missing_gene, count in missing_gene_counts.items():
        genes = find_broken_gene(genes_data[missing_gene].incomplete_hits)
        broken_genes[missing_gene] = genes[count:]

    # Under fortunate circumstances if we can only a single incomplete loci remaining, we can
    # associate the broken genes with it
    if len([l for l in loci_data if not l.is_complete]) == 1:
        logging.info('Associating all missing genes with specific loci')
        for locus_data in loci_data:
            if locus_data.is_complete:
                locus_data.broken_genes = [g for g in broken_genes.values()]
        return
    else:
        return broken_genes


def write_results(loci_data, broken_genes):
    # TODO: clean this up
    # TODO: print much more info, decide on best format
    # TODO: write results to file. sequences? summary? annotation? png?
    logging.info('Writing results')
    for locus_data in loci_data:
        print(*locus_data.contigs, sep=',', end ='\t')
        print(*locus_data.genes, sep=',', end ='\t')
        print(*(g.sseqid for g in locus_data.broken_genes), sep=',', end ='\t')
        print(*(g.sseqid for g in locus_data.missing_genes), sep=',', end ='\t')
        print(*(r for r, h in locus_data.type_hits), sep=',', end ='\t')

        # This awful print out will go. Staying for now
        genes = [g for g in list(locus_data.genes.values()) + locus_data.broken_genes]
        positions = [int(g.qstart) for g in genes] + [int(g.qend) for g in genes]
        print(min(positions), max(positions), sep='\t')


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

    # Parameters
    if args.gene_coverage <= 0:
        logging.error('--gene_coverage must be greater than zero')
        sys.exit(1)
    if args.gene_coverage > 1:
        logging.error('--gene_coverage cannot be greater than 1.0')
        sys.exit(1)


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
    dependencies = {'blastn': {
                        'vcommand': 'blastn -version',
                        'vregex': re.compile(r'^blastn: (.+)\n'),
                        'vrequired': '2.2.28'},
                    'makeblastdb': {
                        'vcommand': 'makeblastdb -version',
                        'vregex': re.compile(r'^makeblastdb: (.+)\n'),
                        'vrequired': '2.2.28'}}
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


def create_blast_database(input_fp, output_dir):
    logging.debug('Create BLAST database for %s', input_fp)
    output_fp = pathlib.Path(output_dir, input_fp.name)
    command = 'makeblastdb -in %s -out %s -dbtype nucl' % (input_fp, output_fp)
    execute_command(command)
    return output_fp


def blast_query(query_fp, blast_db_fp):
    logging.debug('Aligning %s using BLAST to sequences in %s', query_fp, blast_db_fp)
    command = 'blastn -task blastn -db %s -query %s -outfmt "6 %s"'
    result = execute_command(command % (blast_db_fp, query_fp, ' '.join(BLAST_FORMAT)))
    return result.stdout


def parse_blast_stdout(blast_results):
    logging.debug('Parsing %s BLAST results', len(blast_results))
    line_token_gen = (line.split() for line in blast_results.split('\n'))
    return [BlastResults(*lts) for lts in line_token_gen if lts]


def sort_flanking_hits(blast_results, coverage_minimum):
    complete_genes = list()
    incomplete_genes = list()
    for result in blast_results:
        if int(result.length) / int(result.slen) >= coverage_minimum:
            complete_genes.append(result)
        else:
            incomplete_genes.append(result)
    msg = 'Found %s complete %s gene'
    msg = '%ss' % msg if len(complete_genes) > 1 else msg
    logging.info(msg, len(complete_genes), complete_genes[0].sseqid)
    return complete_genes, incomplete_genes


def get_best_loci(genes_data):
    # Pull complete flanking genes and group by contig
    logging.info('Finding best match for loci')
    genes_gen = (g for gene_data in genes_data.values() for g in gene_data.complete_hits)
    contig_genes = dict()
    for gene in genes_gen:
        try:
            contig_genes[gene.qseqid].append(gene)
        except KeyError:
            contig_genes[gene.qseqid] = [gene]

    msg = 'Region I and III genes found on %s contig'
    msg = '%ss' % msg if len(contig_genes) > 1 else msg
    logging.info(msg, len(contig_genes))

    loci_data = list()
    loci_number = 0
    for contig in contig_genes:
        genes = sorted(contig_genes[contig], key=lambda k: int(k.qstart))
        logging.info('Loci structure on contig %s: %s', contig, ' '.join(g.sseqid for g in genes))
        for start, end in find_loci_boundaries(genes):
            loci_number += 1
            logging.debug('Locus boundary found at %s %s on %s', start, end, contig)
            logging.debug('Locus contains %s', ' '.join(g.sseqid for g in genes[start:end]))
            loci_data.append(Locus(loci_number, contig, genes[start:end]))
    return loci_data


def find_loci_boundaries(genes):
    genes_names = [g.sseqid for g in genes]
    forward_alignments = get_ordered_locus_sequences(genes_names, direction='forward')
    reverse_alignments = get_ordered_locus_sequences(genes_names, direction='reverse')
    alignments = select_locus_sequences(forward_alignments + reverse_alignments)

    # Convert ranges to indices
    return [(min(a), max(a) + 1) for a in alignments]


def get_ordered_locus_sequences(genes, *, direction):
    if direction == 'forward':
        ordered = FLANK_ORDER
    elif direction == 'reverse':
        ordered = FLANK_ORDER[::-1]
    else:
        logging.error('Direction must be either forward or reverse')
        sys.exit(1)
    alignments = list()
    for qi in range(len(genes)):
        # Find 'seed'
        for si, sgene in enumerate(ordered):
            if genes[qi] == sgene:
                break
        else:
            alignments.append(range(0, 0))
            continue
        # Extend until a mismatch
        align_pairs = zip(genes[qi:], ordered[si:])
        for j, (qgene, sgene) in enumerate(align_pairs):
            if qgene != sgene:
                alignments.append(range(qi, qi + j))
                break
        else:
            alignments.append(range(qi, qi + j + 1))
    return alignments


def select_locus_sequences(locus_sequences):
    # Iteratively select best alignments
    alignment_score_groups = dict()
    for alignment in locus_sequences:
        try:
            alignment_score_groups[len(alignment)].append(alignment)
        except KeyError:
            alignment_score_groups[len(alignment)] = [alignment]
    best_alignments = list()
    for score_group in sorted(alignment_score_groups, reverse=True):
        # This is not efficient
        for alignment in alignment_score_groups[score_group]:
            # Only record alignment if there is no overlap with previously added alignment
            if not any(p in a for a in best_alignments for p in alignment):
                best_alignments.append(alignment)
    return best_alignments


def consolidate_loci_contigs(loci_data):
    # If the genome contains a full complement of cap genes but are on different contigs, we
    # attempt to resolve this here
    incomplete_loci = [l for l in loci_data if not l.is_complete]
    contigs = {l.contigs[0] for l in incomplete_loci}
    if len(contigs) < 2:
        return
    else:
        logging.info('%s contigs contain incomplete loci, attemping to consolidate', len(contigs))

    # Count the number of genes
    # If there is more than one of any gene then the solution is ambiguous
    incomplete_loci_genes = dict()
    for locus_data in loci_data:
        for gene in locus_data.genes:
            try:
                incomplete_loci_genes[gene] += 1
            except KeyError:
                incomplete_loci_genes[gene] = 1
    if any(count > 1 for count in incomplete_loci_genes.values()):
        logging.warning('Unable to consolidate genes, ambiguous solution with more than one locus')
        return

    # When we have only a single loci, consolidate all complete genes into a single Locus instance
    locus_data = Locus(len(loci_data) + 1)
    locus_data.contigs = contigs
    locus_data.genes = {gene.sseqid: gene for l in incomplete_loci for gene in l.genes.values()}
    msg = 'Consolidated loci %s to new locus %s'
    logging.info(msg, ' '.join(str(l.identifier) for l in incomplete_loci), locus_data.identifier)
    msg = 'Locus %s contains %s but is discontiguous'
    logging.info(msg , locus_data.identifier, ' '.join(locus_data.genes))

    # Remove consolidated loci from loci_data and add new instance
    loci_data.append(locus_data)
    for locus_data in incomplete_loci:
        loci_data.remove(locus_data)


def find_broken_gene(gene_hits, locus_data=None):
    min_length = 50
    if locus_data is not None:
        hit_gen = (h for h in gene_hits if h.qseqid in locus_data.contigs)
    else:
        hit_gen = (h for h in gene_hits)

    broken_genes = list()
    for gene_hit in hit_gen:
        if broken_hit_identity(gene_hit, locus_data):
            continue
        if broken_hit_distance(gene_hit, locus_data):
            continue
        # Remove hit from gene_hits
        gene_hits.remove(gene_hit)
        if locus_data is not None:
            return gene_hit
        else:
            broken_genes.append(gene_hit)
    return broken_genes


def broken_hit_identity(gene_hit, locus_data):
    # TODO: expose to command line
    max_identity_diff = 10
    average_identity = statistics.mean(float(h.pident) for h in locus_data.genes.values())
    lower_identity = average_identity - max_identity_diff
    upper_identity = average_identity + max_identity_diff
    return float(gene_hit.pident) < lower_identity or float(gene_hit.pident) > upper_identity


def broken_hit_distance(gene_hit, locus_data):
    # TODO: expose to command line
    max_dist = 5000
    gene_positions = [int(p) for h in locus_data.genes.values() for p in (h.qstart, h.qend)]
    lower_position = max(gene_positions) - max_dist
    upper_position = max(gene_positions) + max_dist
    return int(gene_hit.qstart) < lower_position or int(gene_hit.qstart) > upper_position


def get_region_two_sequence(locus_data, query_fasta):
    if 'hcsA' not in locus_data.genes or 'bexD' not in locus_data.genes:
        logging.warning('Can\'t extract region II for locus %s', locus_data.identifier)
        return

    bexD = locus_data.genes['bexD']
    hcsA = locus_data.genes['hcsA']
    if bexD.qseqid != hcsA.qseqid:
        msg = 'Can\'t extract region II for locus %s, sequences appear to be discontiguous'
        logging.warning(msg, locus_data.identifier)
        return

    # TODO: sort is likely not required, need to think through this more
    bexD_lower, bexD_upper = sorted((int(bexD.qstart), int(bexD.qend)))
    hcsA_lower, hcsA_upper = sorted((int(hcsA.qstart), int(hcsA.qend)))

    # TODO: sorted here is absolutely required
    if abs(bexD_upper - hcsA_lower) < abs(hcsA_upper - bexD_lower):
        bounds_lower, bounds_upper = sorted((bexD_upper, hcsA_lower))
    else:
        bounds_lower, bounds_upper = sorted((hcsA_upper, bexD_lower))
    return query_fasta[locus_data.contigs[0]][bounds_lower:bounds_upper]


def filter_type_hits(blast_results, coverage_minimum, identity_minimum):
    for result in blast_results:
        if int(result.length) / int(result.slen) < coverage_minimum:
            continue
        if float(result.pident) < identity_minimum:
            continue
        return result


if __name__ == '__main__':
    main()
