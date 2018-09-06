import logging
import pathlib
import subprocess
import tempfile
import sys


from . import annotation
from . import arguments
from . import database
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

    # Find complete genes
    hits.complete = database.filter_hits(hits.all, coverage_min=args.gene_coverage, identity_min=args.gene_identity)
    if not hits.complete:
        logging.info('No hits to any cap locus gene found, exiting')
        sys.exit(0)

    # Find missing genes
    missing_genes = database.discover_missing_genes(hits.complete)
    if missing_genes:
        hits.broken = database.filter_hits(hits.remaining, identity_min=args.broken_gene_identity, length_min=args.broken_gene_length)

    # Assign hits to ORFs
    orfs_assigned = database.match_orfs_and_hits(hits.passed, orfs_all)

    # Annotate loci and serotype
    for loci in database.characterise_loci(orfs_assigned):
        # TODO: smallest of start and end
        # TODO: print region two gene names
        print(loci.contig, loci.orfs[0].start, loci.orfs[-1].end, end='\t')
        print(*(hit for orf in loci.orfs for hit in orf.hits), sep=',')


if __name__ == '__main__':
    main()
