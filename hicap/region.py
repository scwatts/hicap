import logging
import tempfile


from . import blast


REGION_ONE = ('bexA', 'bexB', 'bexC', 'bexD')
REGION_TWO = ('type_a', 'type_b', 'type_c', 'type_d', 'type_e', 'type_f')
REGION_THREE = ('hcsA', 'hcsB')


class Database():

    def __init__(self, region):
        self.region = region

        self.fps = list()
        self.parameters = dict()


class Hits():

    def __init__(self, name, complete_hits, partial_hits):
        self.name = name
        self.complete_hits = complete_hits
        self.partial_hits = partial_hits

        self.missing = 0
        self.broken_hits = list()


class LocusData():

    def __init__(self, contig, hits):
        self.contig = contig
        self.hits = hits
        self._positions = [p for h in self.hits for p in (h.qstart, h.qend)]
        self.qstart = min(self._positions)
        self.qend = max(self._positions)

        self.sequence = None
        self.sequence_offset = None
        self.orfs = None
        self.genbank = None


# TODO: use more robust database structure
def init_databases(args):
    '''Aggregate input database filepaths into respective region'''
    # Create gene->region map
    database_map = dict()
    for region, genes in zip(('one', 'two', 'three'), (REGION_ONE, REGION_TWO, REGION_THREE)):
        for gene in genes:
            database_map[gene] = region

    # Set filepaths and then parameters
    databases = {region: Database(region) for region in database_map.values()}
    for database_fp in args.database_fps:
        database_region = database_map[database_fp.stem]
        databases[database_region].fps.append(database_fp)

    for database in databases.values():
        if database.region in ('one', 'three'):
            database.parameters = {'coverage_min': args.gene_coverage,
                                   'identity_min': args.gene_identity,
                                   'broken_identity_min': args.broken_gene_identity,
                                   'broken_length_min': args.broken_gene_length}
        elif database.region == 'two':
            database.parameters = {'coverage_min': args.type_coverage,
                                   'identity_min': args.type_identity,
                                   'broken_identity_min': args.broken_type_identity,
                                   'broken_length_min': args.broken_type_length}
    return databases


def align_region(query_fp, database_fps):
    '''Find all alignments between query and each gene database'''
    results = dict()
    with tempfile.TemporaryDirectory() as dh:
        for database_fp in database_fps:
            blast_database_fp = blast.create_blast_database(database_fp, dh)
            blast_stdout = blast.align(query_fp, blast_database_fp)
            results[database_fp.stem] = blast.parse_blast_stdout(blast_stdout)
    return results


def search_region(query_fp, database):
    logging.info('Searching for region %s hits', database.region)
    results = dict()
    for gene, hits in align_region(query_fp, database.fps).items():
        results[gene] = filter_blast_hits(gene, hits, database.parameters)
    return results


def filter_blast_hits(name, blast_results, params):
    complete_hits = list()
    partial_hits = list()
    for result in blast_results:
        if result.length / result.slen < params['coverage_min']:
            partial_hits.append(result)
        elif result.pident < params['identity_min']:
            partial_hits.append(result)
        else:
            complete_hits.append(result)
    if complete_hits:
        logging.info('Complete %s hits: %s', name, len(complete_hits))
    else:
        logging.info('No complete %s hits found', name)
    if partial_hits:
        logging.debug('Partial %s hits: %s', name, len(partial_hits))
    return Hits(name, complete_hits, partial_hits)


def count_and_set_missing(region_hits):
    counts = count_units(region_hits)
    total_loci = max(counts.values())
    complete_loci = min(counts.values())
    incomplete_loci = total_loci - complete_loci
    msg = 'Found %s full and %s partial complement(s) of the cap locus'
    logging.info(msg, complete_loci, incomplete_loci)

    missing_count = set_missing(region_hits, total_loci)
    if incomplete_loci:
        logging.info('Missing %s gene(s) to complete partial complements', len(missing_count))


def count_units(region_hits):
    '''Get count of good hits for genes and types'''
    counts = dict()
    for gene, gene_hits in region_hits['one'].items():
        counts[gene] = len(gene_hits.complete_hits)
    for gene, gene_hits in region_hits['two'].items():
        try:
            counts['region_two'] += len(gene_hits.complete_hits)
        except KeyError:
            counts['region_two'] = len(gene_hits.complete_hits)
    for gene, gene_hits in region_hits['three'].items():
        counts[gene] = len(gene_hits.complete_hits)
    return counts


def set_missing(region_hits, total_loci):
    '''For each Hits instance, set the number of missing genes/types'''
    missing_counts = dict()
    for region_name, region in region_hits.items():
        for hits in region.values():
            missing_count = total_loci - len(hits.complete_hits)
            if region_name == 'two':
                if len(hits.complete_hits) > 0 and missing_count:
                    logging.debug('Missing %s region two', missing_count)
            elif missing_count > 0:
                logging.debug('Missing %s %s', missing_count, hits.name)
                hits.missing = missing_count
                missing_counts[hits.name] = missing_count
    return missing_counts


def find_missing(region_hits, databases):
    '''For gene Hit instance, attempt to find broken/ truncated genes'''
    for (region_name, region), database in zip(region_hits.items(), databases.values()):
        for hits in region.values():
            for hit in hits.partial_hits:
                if hits.missing == 0:
                    break
                if hit.length < database.parameters['broken_length_min']:
                    continue
                elif hit.pident < database.parameters['broken_identity_min']:
                    continue
                else:
                    hits.missing -= 1
                    hits.broken_hits.append(hit)


def aggregate_hits_and_sort(region_hits):
    '''Sort complete and broken hits into contig and then order by position'''
    contig_hits = dict()
    for region in region_hits.values():
        for region_hits in region.values():
            for hit in region_hits.complete_hits + region_hits.broken_hits:
                try:
                    contig_hits[hit.qseqid].append(hit)
                except KeyError:
                    contig_hits[hit.qseqid] = [hit]
    # TODO: expand logical here. for now we consider all contigs hits as the same locus
    for contig, hits in contig_hits.items():
        sorted_hits = sorted(hits, key=lambda k: k.qstart)
        contig_hits[contig] = LocusData(contig, sorted_hits)
    return contig_hits


def extract_sequences(contig_hits, query_fastas):
    flanking = 2000
    for contig, locus_data in contig_hits.items():
        # TODO: break large gaps
        # Python takes care of slicing out of bounds but must define behaviour for -ve numbers
        bound_lower = locus_data.qstart - flanking if locus_data.qstart > flanking else 0
        bound_upper = locus_data.qend + flanking
        locus_data.sequence = query_fastas[locus_data.contig][bound_lower:bound_upper]
        locus_data.sequence_offset = bound_lower
