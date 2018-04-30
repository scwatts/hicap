import logging
import pathlib
import sys


import Bio.SeqIO


from . import annotation
from . import arguments
from . import region
from . import utility


# TODO: improve fuzzy range matching with BLAST htis and predicted ORFs


def main():
    # Init
    args = arguments.get_args()
    utility.initialise_logging(args.log_level, args.log_fp)
    utility.check_dependencies()
    arguments.check_args(args)
    databases = region.init_databases(args)

    # Find complete good hits
    for locus_region in databases:
        region.search_database(args.query_fp, databases[locus_region])

    # Resolve overlaps in region II, if any
    region.resolve_region_two_overlaps(databases['two'])

    # Attempt to locate missing, truncated genes
    region.count_and_set_missing(databases)
    region.find_missing(databases)
    if not region.total_items_found(databases):
        logging.warning('No cap locus genes found, exiting')
        sys.exit(0)

    ## Sort into contigs and extract sequence
    query_fastas = utility.read_fasta(args.query_fp)
    loci_data = region.aggregate_hits_and_sort(databases)
    region.extract_sequences(loci_data, query_fastas)

    ## Annotate loci data and write out data
    annotation.discover_orfs(loci_data)
    annotation.generate_genbank(loci_data, args.query_fp.stem)
    logging.info('Writing genbank records')
    output_gbk_fp = pathlib.Path(args.output_dir, '%s.gbk' % args.query_fp.stem)
    output_summary_fp = pathlib.Path(args.output_dir, '%s.tsv' % args.query_fp.stem)

    with output_gbk_fp.open('w') as fh:
        for locus_data in loci_data.values():
            Bio.SeqIO.write(locus_data.genbank, fh, 'genbank')

    with output_summary_fp.open('w') as fh:
        print('#', end='', file=fh)
        print(*{serotype for serotype in databases['two'].hits}, file=fh)
        for locus_data in loci_data.values():
            print(locus_data.contig, end='\t', file=fh)
            print(locus_data.qstart, end='\t', file=fh)
            print(locus_data.qend, end='\t', file=fh)
            genes = list()
            for hit in locus_data.hits:
                gene = hit.sseqid
                if hit.broken:
                    gene = '%s*' % gene
                genes.append(gene)
            print(*genes, sep=',', file=fh)


if __name__ == '__main__':
    main()
