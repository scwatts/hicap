import logging
import pathlib
import tempfile
import sys


from . import annotation
from . import arguments
from . import database
from . import locus
from . import utility


def main():
    # Init
    args = arguments.get_args()
    utility.initialise_logging(args.log_level, args.log_fp)
    utility.check_dependencies()
    arguments.check_args(args)

    # Collect ORFs from input assembly
    orfs_all = annotation.collect_orfs(args.query_fp)

    # Align ORFs to databases using BLAST
    with tempfile.TemporaryDirectory() as dh:
        orfs_fp = pathlib.Path(dh, 'orfs.fasta')
        with orfs_fp.open('w') as fh:
            for i, orf in enumerate(orfs_all):
                print('>%s' % i, orf.sequence, sep='\n', file=fh)
        hits = database.search(orfs_fp, args.database_fps)

    # Set the respective ORF for each hit and get contig sizes
    hits = database.assign_hit_orfs(hits, orfs_all)

    # Find complete hits
    hits_complete = database.filter_hits(hits, coverage_min=args.gene_coverage, identity_min=args.gene_identity)
    hits_remaining = hits - hits_complete
    if not hits_complete:
        logging.info('No hits to any cap locus gene found, exiting')
        sys.exit(0)

    # TODO: print some info here wrt complete hits
    region_groups = dict()
    all_region_hits = locus.region_sort_hits(hits_complete)
    filter_params = {'identity_min': args.broken_gene_identity, 'length_min': args.broken_gene_length}
    for region, region_hits in all_region_hits.items():
        group = locus.discover_region_clusters(region_hits, hits_remaining, region, filter_params)
        region_groups[region] = group

    # If no completed hits were found for region two, attempt to find fragmented ORFs
    if not region_groups['two'].hits:
        group = locus.locate_fragmented_region_two(region_groups, hits_remaining, filter_params)
        region_groups['two'] = group

    # Create genbank record, report, and summary
    for group in region_groups.values():
        print(group)
        hits_sorted = sorted(group.hits, key=lambda k: k.orf.start)
        for hit in hits_sorted:
            print(hit.sseqid, hit.orf.start, hit.orf.end, hit.orf.contig, hit.broken, sep='\t')
