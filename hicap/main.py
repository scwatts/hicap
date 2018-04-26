import logging
import pathlib
import sys


import Bio.SeqIO


from . import annotation
from . import arguments
from . import region
from . import utility


# TODO: split type references into genes


def main():
    # Init
    args = arguments.get_args()
    utility.initialise_logging(args.log_level, args.log_fp)
    utility.check_dependencies()
    arguments.check_args(args)
    databases = region.init_databases(args)

    # Find complete good hits
    region_hits = dict()
    for locus_region in databases:
        region_hits[locus_region] = region.search_region(args.query_fp, databases[locus_region])

    # Attempt to locate missing, truncated genes
    region.count_and_set_missing(region_hits)
    region.find_missing(region_hits, databases)
    items_found = sum(len(h.complete_hits + h.broken_hits) for r in region_hits.values() for h in r.values())
    if not items_found:
        logging.warning('No cap locus genes found, exiting')
        sys.exit(0)

    # Sort into contigs and extract sequence
    query_fastas = utility.read_fasta(args.query_fp)
    loci_data = region.aggregate_hits_and_sort(region_hits)
    region.extract_sequences(loci_data, query_fastas)

    # Annotate loci data and write out data
    annotation.discover_orfs(loci_data)
    annotation.generate_genbank(loci_data, args.query_fp.stem)
    # TODO: output some text info as well (serotype, complete genes, broken genes, etc)
    logging.info('Writing genbank records')
    for i, locus_data in enumerate(loci_data.values(), 1):
        output_fp = pathlib.Path(args.output_dir, '%s_%s.gbk' %(args.query_fp.name, i))
        with output_fp.open('w') as fh:
            Bio.SeqIO.write(locus_data.genbank, fh, 'genbank')


if __name__ == '__main__':
    main()
