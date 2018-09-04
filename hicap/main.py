import logging
import pathlib
import subprocess
import tempfile


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
        hits_all = database.search(orfs_fp, args.database_fps)

    # Assign hits to ORFs and record those with filtered hits
    hits_filtered = database.filter_hits(hits_all, coverage_min=args.gene_coverage,
                                            identity_min=args.gene_identity)
    orf_indices = set()
    for region, hits in hits_filtered.items():
        for hit in hits:
            orf_index = int(hit.qseqid)
            orf_indices.add(orf_index)
            try:
                orfs_all[orf_index].hits[region].append(hit)
            except KeyError:
                orfs_all[orf_index].hits[region] = [hit]

    # Extract ORFs and try to find missing genes (truncated or broken)
    orfs_complete = [orfs_all[index] for index in orf_indices]
    genes_broken = database.discover_missing_genes(orfs_complete)
    orfs_broken = database.collect_missing_orfs(genes_broken, orfs_all, hits_all, hits_filtered,
                                                    args.broken_gene_identity, args.broken_gene_length)

    # Annotate loci and serotype
    orfs = orfs_complete + orfs_broken
    for loci in database.characterise_loci(orfs):
        # TODO: smallest of start and end
        # TODO: print region two gene names
        print(loci.contig, loci.orfs[0].start, loci.orfs[-1].end, end='\t')
        print(*(hit for orf in loci.orfs for hit in orf.hits), sep=',')


if __name__ == '__main__':
    main()
