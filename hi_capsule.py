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
import statistics
import sys
import tempfile


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

    def __init__(self, contig, genes):
        self.contigs = [contig]
        self.genes = genes

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
        gene_name_gen = (g.sseqid for g in self.genes)
        return set(gene_name_gen) ^ set(FLANK_ORDER)


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
    genes_data = dict()
    for region_genes_data, region in zip((REGION_ONE, REGION_THREE), ('one', 'two')):
        for gene in region_genes_data:
            genes_data[gene] = GeneData(gene, region)

    for database_fp in args.database_fps:
        database_gene = database_fp.stem
        if database_gene in REGION_TWO:
            continue
        try:
            genes_data[database_gene].database_fp = database_fp
        except KeyError:
            logging.error('Can\'t match database %s to flanking gene', database_fp)
            sys.exit(1)

    # Run
    with tempfile.TemporaryDirectory() as dh:
        # Blast query against database and collect complete flanking gene hits
        # TODO: q parallelise, use asyncio
        for gene in genes_data.values():
            gene.blast_database_fp = create_blast_database(gene.database_fp, dh)
            blast_stdout = blast_query(args.query_fp, gene.blast_database_fp)
            gene.blast_results = parse_blast_stdout(blast_stdout)
            # Set hits as complete or otherwise
            complete_hits, incomplete_hits = sort_hits(gene.blast_results, args.gene_coverage)
            gene.complete_hits = complete_hits
            gene.incomplete_hits = incomplete_hits

        # Find the best matching locus sequence for each
        # TODO: CRIT: check for complete genes dropped here
        loci_data = get_best_loci(genes_data)

        # Attempt to find broken genes for incomplete loci on the same contig
        for incomplete_locus in (l for l in loci_data if not l.is_complete):
            for missing_gene in incomplete_locus.missing_genes:
                broken_gene = find_broken_gene(genes_data[missing_gene].incomplete_hits, locus_data)
                locus_data.broken_genes.append(broken_gene)

        # Next try to find broken genes on different contigs
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
        if len(l for l in loci_data if not l.is_complete) == 1:
            for locus_data in loci_data:
                if locus_data.is_complete:
                    locus_data.broken_genes = [g for g in broken_genes.values()]

        # TODO: blast region II

        # TODO: software version checking
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


def sort_hits(blast_results, coverage_minimum):
    complete_genes = list()
    incomplete_genes = list()
    for result in blast_results:
        if int(result.length) / int(result.slen) > coverage_minimum:
            complete_genes.append(result)
        else:
            incomplete_genes.append(result)
    return complete_genes, incomplete_genes


def get_best_loci(genes_data):
    # Pull complete flanking genes and group by contig
    genes_gen = (g for gene_data in genes_data.values() for g in gene_data.complete_hits)
    contig_genes = dict()
    for gene in genes_gen:
        try:
            contig_genes[gene.qseqid].append(gene)
        except KeyError:
            contig_genes[gene.qseqid] = [gene]

    loci_data = list()
    for contig in contig_genes:
        genes = sorted(contig_genes[contig], key=lambda k: int(k.qstart))
        loci = [Locus(contig, genes[s:e]) for s, e in find_loci_boundaries(genes)]
        loci_data.extend(loci)
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
            # Exclude alignments shorter than 2 genes
            if len(alignment) < 3:
                continue
            # Only record alignment if there is no overlap with previously added alignment
            if not any(p in a for a in best_alignments for p in alignment):
                best_alignments.append(alignment)
    return best_alignments


def find_broken_gene(gene_hits, locus_data=None):
    min_length = 50
    if same_contig:
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
        if same_contig:
            return gene_hit
        else:
            broken_genes.append(gene_hit)
    return broken_genes


def broken_hit_identity(gene_hit, locus_data):
    # TODO: expose to command line
    max_identity_diff = 10
    average_identity = statistics.mean(float(h.pident) for h in locus_data.genes)
    lower_identity = average_identity - max_identity_diff
    upper_identity = average_identity + max_identity_diff
    return float(gene_hit.pident) < lower_identity or float(gene_hit.pident) > upper_identity


def broken_hit_distance(gene_hit, locus_data):
    # TODO: expose to command line
    max_dist = 5000
    gene_positions = [int(p) for h in locus_data.genes for p in (h.qstart, h.qend)]
    lower_position = max(gene_positions) - max_dist
    upper_position = max(gene_positions) + max_dist
    return int(gene_hit.qstart) < lower_position or int(gene_hit.qstart) > upper_position



if __name__ == '__main__':
    main()
