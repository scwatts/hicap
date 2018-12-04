import logging
import pathlib
import tempfile
import sys


from . import annotation
from . import arguments
from . import database
from . import locus
from . import report
from . import utility


def main():
    # Init
    args = arguments.get_args()
    utility.initialise_logging(args.log_level, args.log_fp)
    utility.check_dependencies()
    arguments.check_args(args)

    # Check FASTA input contig names aren't too long for the genbank format
    if any(len(desc) > 20 for desc in utility.read_fasta(args.query_fp)):
        msg = ('One or more contig names exceed the genbank spec limit of 20 characters.'
               ' These will be truncated in the genbank output file')
        logging.warning(msg)

    # Collect ORFs from input assembly then align ORFs to database and assign ORFs to hits
    logging.info('Searching database for matches')
    orfs_all = annotation.collect_orfs(args.query_fp, args.model_fp)
    hits = database.search(orfs_all, args.database_fps, args.threads)
    hits = database.assign_hit_orfs(hits, orfs_all)

    # Find complete hits
    hits_complete = database.filter_hits(hits, coverage_min=args.gene_coverage, identity_min=args.gene_identity)
    hits_remaining = hits - hits_complete
    if not hits_complete:
        logging.info('No hits to any cap locus gene found, exiting')
        sys.exit(0)
    utility.log_search_hits_found(hits_complete)

    # Selected best complete hits and search for hits of broken/ truncated genes
    region_groups = dict()
    all_region_hits = locus.sort_hits_by_region(hits_complete)
    filter_params = {'identity_min': args.broken_gene_identity, 'length_min': args.broken_gene_length}
    logging.info('Discovering loci region clusters')
    for region, region_hits in all_region_hits.items():
        group = locus.discover_region_clusters(region_hits, hits_remaining, region, filter_params)
        region_groups[region] = group

    # If no completed hits were found for region two, attempt to find fragmented ORFs
    if not region_groups['two'].hits:
        group = locus.locate_fragmented_region_two(region_groups, hits_remaining, filter_params)
        region_groups['two'] = group

    # For any gene, attempt to find fragments proximal to previously dsicovered ORFs
    contigs_fasta = utility.read_fasta(args.query_fp)
    locus.find_proximal_fragments(region_groups, hits_remaining, contigs_fasta)

    # Collect ORFs that are not apart of the Hi cap loci in the surrounding areas
    # TODO: contig min/max(orf.start, orf.end, boundary)
    nearby_orfs = locus.collect_nearby_orfs(region_groups, orfs_all)

    # Generate output data and files
    report.write_outputs(region_groups, nearby_orfs, args)
