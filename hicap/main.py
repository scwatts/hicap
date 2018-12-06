import logging
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
    hits = database.search(orfs_all, args.gene_database_fps, args.threads)
    hits = database.assign_hit_orfs(hits, orfs_all)

    # Find complete hits
    hits_complete = database.filter_hits(hits, coverage_min=args.gene_coverage, identity_min=args.gene_identity)
    hits_remaining = hits - hits_complete
    if not hits_complete:
        logging.info('No hits to any cap locus gene found, exiting')
        sys.exit(0)
    utility.log_search_hits_found(hits_complete)

    # Selected best complete hits and search for hits of broken/ truncated genes
    locus_data = locus.LocusData()
    all_region_hits = locus.sort_hits_by_region(hits_complete)
    filter_params = {'identity_min': args.broken_gene_identity, 'length_min': args.broken_gene_length}
    logging.info('Discovering loci region clusters')
    for region, region_hits in all_region_hits.items():
        region_data = locus.discover_region_clusters(region_hits, hits_remaining, region, filter_params)
        locus_data.regions[region] = region_data

    # If no completed hits were found for region two, attempt to find fragmented ORFs
    if not locus_data.regions['two'].orf_hits:
        region_data = locus.locate_fragmented_region_two(locus_data.regions, hits_remaining, filter_params)
        locus_data.regions['two'] = region_data

    # For any gene, attempt to find fragments proximal to previously discovered ORFs
    contig_fastas = utility.read_fasta(args.query_fp)
    locus.find_proximal_fragments(locus_data.regions, hits_remaining, contig_fastas)

    # For any genes which are missing, attempt to find via basic BLAST
    locus.blast_missing_genes(locus_data.regions, contig_fastas, args.gene_database_fps)

    # Collect ORFs not apart of the Hi cap loci in surrounding areas and search for IS1016
    locus_data.is_hits = locus.discover_is1016(locus_data.regions, contig_fastas, args.is_database_fp)
    locus_data.nearby_orfs = locus.collect_nearby_orfs(locus_data, orfs_all)

    # Generate output data and files
    report.write_outputs(locus_data, args)
