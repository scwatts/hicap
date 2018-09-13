import logging
import pathlib
import tempfile
import sys


import Bio.SeqIO


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
    contig_sizes = utility.get_contig_sizes(args.query_fp)

    # Identify with ORFs are near contig bounds. Explicit assignment for clarity
    orfs_all = locus.identify_orfs_near_boundaries(orfs_all, contig_sizes, 1000)

    # Find complete hits
    hits_complete = database.filter_hits(hits, coverage_min=args.gene_coverage, identity_min=args.gene_identity)
    hits_remaining = hits - hits_complete
    if not hits_complete:
        logging.info('No hits to any cap locus gene found, exiting')
        sys.exit(0)

    # TODO: print some info here wrt complete hits

    # Find region clusters. Select most likely hits. Find broken neighbouring genes
    # for region, region_hits in all_region_hits.items():
    #   clusters = discover_clusters(region_hits, hits_remaining, contig_sizes, region, filter_params)
    clusters = list()
    all_region_hits = locus.region_sort_hits(hits_complete)
    filter_params = {'identity_min': args.broken_gene_identity, 'length_min': args.broken_gene_length}
    for region, region_hits in all_region_hits.items():
        region_clusters = locus.discover_region_clusters(region_hits, hits_remaining, region, contig_sizes, filter_params)
        clusters.extend(region_clusters)


    for group in clusters:
        print(group)
        for hit in group.hits:
            print(hit.sseqid, hit.orf.start, hit.orf.end, hit.orf.contig)

    ## Find region clusters. Select most likely hits. Find broken neighbouring genes
    #hits_remaining = hits - hits_complete
    #filter_params = {'identity_min': args.broken_gene_identity, 'length_min': args.broken_gene_length}
    #region_one_clusters = database.discover_region_clusters(hits_complete, hits_remaining, 'one', contig_sizes, filter_params)
    ##region_two_clusters = database.discover_region_two_clusters(hits_complete, hits_remaining, contig_sizes, filter_params)
    #for group in region_one_clusters:
    #    print(group)
    #    for hit in group.hits:
    #        print(hit.sseqid, hit.orf.contig, hit.orf.start)

    ## Find missing genes in region one and three
    #missing_genes = database.count_missing_genes_flank(hits_complete)
    #if missing_genes:
    #    hits_remaining = hits - hits_complete
    #    filter_params = {'identity_min': args.broken_gene_identity, 'length_min': args.broken_gene_length}
    #    hits_broken = database.discover_missing_genes_flank(hits_remaining, missing_genes, filter_params)

    ## Assign hits to ORFs
    #hits_passed = hits_complete | hits_broken
    #missing_genes = database.count_missing_genes_region_two(hits_passed)
    ##orfs_assigned = database.match_orfs_and_hits(hits_passed, orfs_all)



    ## Find complete genes
    #hits.complete = database.filter_hits(hits.all, coverage_min=args.gene_coverage, identity_min=args.gene_identity)
    #if not sum(len(database_hits) for database_hits in hits.complete.values()):
    #    logging.info('No hits to any cap locus gene found, exiting')
    #    sys.exit(0)
    #else:
    #    for database_name, database_hits in hits.complete.items():
    #        message = 'Found %s %s for %s'
    #        hit_quant = 'hits' if len(database_hits) > 1 else 'hit'
    #        logging.info(message, len(database_hits), hit_quant, database_name)

    ## TODO: continuity check - accept gene groups within physical proximity. this may involve
    ## inspect of gene distance from contig boundaries

    ## Find missing genes in region one and three - select the best hit for each missing gene
    ## Places results in hits.broken
    #database.discover_missing_genes_flank(hits, args.broken_gene_identity, args.broken_gene_length)

    ## Assign hits to ORFs. The hits.pass property returns all complete and broken hits
    #logging.info('Assigning hits to ORFs')
    #orfs_assigned = database.match_orfs_and_hits(hits.passed, orfs_all)

    ## Collect region two ORFs and find any broken, missing genes and select best hit for each
    #region_two_groups = database.characterise_loci_region_two(orfs_assigned)
    #database.discover_missing_genes_region_two(hits, args.broken_gene_identity, args.broken_gene_length)

    ## Annotate loci and serotype
    #logging.info('Performing loci gene collation and seroptying')
    #contig_sizes = utility.get_contig_sizes(args.query_fp)
    #loci_blocks = database.characterise_loci(orfs_assigned, contig_sizes)

    ## Write out summary table
    #logging.info('Writing summary')
    #output_summary_fp = pathlib.Path(args.output_dir, '%s.tsv' % args.query_fp.stem)
    #utility.write_summary_data(loci_blocks, output_summary_fp)

    ## Create and write out genbank record
    #logging.info('Writing genbank')
    #output_gbk_fp = pathlib.Path(args.output_dir, '%s.gbk' % args.query_fp.stem)
    #genbank_data = utility.create_genbank_record(loci_blocks, args.query_fp)
    #with output_gbk_fp.open('w') as fh:
    #    Bio.SeqIO.write(genbank_data, fh, 'genbank')


if __name__ == '__main__':
    main()
