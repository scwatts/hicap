import logging
import pathlib
import subprocess
import tempfile
import sys


import Bio.SeqIO


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
    else:
        for database_name, database_hits in hits.complete.items():
            message = 'Found %s %s for %s'
            hit_quant = 'hits' if len(database_hits) > 1 else 'hit'
            logging.info(message, len(database_hits), hit_quant, database_name)

    # Find missing genes - select the best hit for each missing gene
    missing_genes = database.discover_missing_genes(hits.complete)
    if missing_genes:
        logging.info('Searching for %s missing genes', sum(count for count in missing_genes.values()))
        hits.broken = database.filter_hits(hits.remaining, identity_min=args.broken_gene_identity, length_min=args.broken_gene_length)
        hits.broken = database.select_best_hits(hits.broken, missing_genes)
        logging.info('Found %s missing genes', len([count for counts in hits.broken.values() for count in counts]))

    # Assign hits to ORFs
    logging.info('Assigning hits to ORFs')
    orfs_assigned = database.match_orfs_and_hits(hits.passed, orfs_all)

    # Annotate loci and serotype
    logging.info('Performing loci gene collation and seroptying')
    loci_blocks = database.characterise_loci(orfs_assigned)

    # Create and write out genbank record
    output_gbk_fp = pathlib.Path(args.output_dir, '%s.gbk' % args.query_fp.stem)
    genbank_data = utility.create_genbank_record(loci_blocks, args.query_fp)
    with output_gbk_fp.open('w') as fh:
        Bio.SeqIO.write(genbank_data, fh, 'genbank')

    # Print results
    output_summary_fp = pathlib.Path(args.output_dir, '%s.tsv' % args.query_fp.stem)
    serotypes = {s for lb in loci_blocks for s in lb.serotypes}
    with output_summary_fp.open('w') as fh:
        print('#', ','.join(serotypes), sep='', file=fh)
        for loci_block in loci_blocks:
            genes = list()
            for orf in loci_block.orfs:
                # There should only ever be one hit per ORF here
                database_name, [orf_hit] = list(orf.hits.items())[0]
                # TODO: is there an appreciable difference if we construct a reverse hash map
                for region, databases in database.SCHEME.items():
                    if database_name in databases:
                        orf_region = region
                        break
                # Get appropriate representation of gene name
                if region == 'two':
                    genes.append(orf_hit.sseqid)
                else:
                    genes.append(database_name)
            start = loci_block.orfs[0].start
            end = loci_block.orfs[-1].end
            print(loci_block.contig, start, end, ','.join(genes), sep='\t', file=fh)


if __name__ == '__main__':
    main()
