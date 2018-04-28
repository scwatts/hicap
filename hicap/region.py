import logging
import tempfile


from . import blast
from . import utility


REGION_ONE = ('bexA', 'bexB', 'bexC', 'bexD')
REGION_TWO = ('type_a', 'type_b', 'type_c', 'type_d', 'type_e', 'type_f')
REGION_THREE = ('hcsA', 'hcsB')


class Database():

    def __init__(self, region):
        self.region = region

        self.fps = list()
        self.parameters = dict()
        self.hits = dict()


class Hits():

    def __init__(self, name, complete_hits, partial_hits):
        self.name = name
        self.complete_hits = complete_hits
        self.partial_hits = partial_hits

        self.missing = 0
        self.broken_hits = list()

    @property
    def all_hits(self):
        return self.complete_hits + self.broken_hits


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
        database_name = database_map[database_fp.stem]
        databases[database_name].fps.append(database_fp)

    for database in databases.values():
        database.parameters = {'coverage_min': args.gene_coverage,
                                'identity_min': args.gene_identity,
                                'broken_identity_min': args.broken_gene_identity,
                                'broken_length_min': args.broken_gene_length}
        if database.region in ('one', 'three'):
            database.filter_func = filter_flat_blast_hits
        elif database.region == 'two':
            database.filter_func = filter_multi_blast_hits
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


def search_database(query_fp, database):
    logging.info('Searching for region %s hits', database.region)
    for gene, hits in align_region(query_fp, database.fps).items():
        database.hits[gene] = database.filter_func(gene, hits, database.parameters)


def gene_hit_generator(database):
    for data in database.hits.values():
        if database.region == 'two':
            for type_hits in data.values():
                yield type_hits
        else:
            yield data


def filter_flat_blast_hits(name, blast_results, params):
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


def filter_multi_blast_hits(name, blast_results, params):
    # Group by sseqid
    blast_results_genes = dict()
    for br in blast_results:
        try:
            blast_results_genes[br.sseqid].append(br)
        except KeyError:
            blast_results_genes[br.sseqid] = [br]
    genes_hits = dict()
    for gene, br in blast_results_genes.items():
        genes_hits[gene] = filter_flat_blast_hits(gene, br, params)
    return genes_hits


def resolve_region_two_overlaps(databases):
    most_probable = list()
    sort_func = lambda k: sum(bool(v.complete_hits) for k, v in databases.hits[k].items())
    region_type_priority = sorted(databases.hits, key=sort_func, reverse=True)
    for region_type in region_type_priority:
        region_type_hits = databases.hits[region_type]
        # Most inefficient
        for gene, gene_hits in region_type_hits.items():
            for hit in gene_hits.complete_hits:
                start, end = sorted((hit.qstart, hit.qend))
                if any(utility.range_overlaps(range(start, end), r) for r in most_probable):
                    gene_hits.complete_hits.remove(hit)
                else:
                    most_probable.append(range(start, end))
    # Remove types without hits
    empty_types = list()
    for region_type, region_type_hits in databases.hits.items():
        if not all(bool(gh.complete_hits) for gh in region_type_hits.values()):
            empty_types.append(region_type)
    for empty_type in empty_types:
        del databases.hits[empty_type]


def count_and_set_missing(databases):
    counts = count_units(databases)
    total_loci = max(counts.values())
    complete_loci = min(counts.values())
    incomplete_loci = total_loci - complete_loci
    msg = 'Found %s full and %s partial complement(s) of the cap locus'
    logging.info(msg, complete_loci, incomplete_loci)

    missing_count = set_missing(databases, total_loci)
    if incomplete_loci:
        logging.info('Missing %s gene(s) to complete partial complements', len(missing_count))


def count_units(databases):
    '''Get count of good hits for genes and types'''
    counts = dict()
    for gene, gene_hits in databases['one'].hits.items():
        counts[gene] = len(gene_hits.complete_hits)
    region_two_hit_gen = (gh for th in databases['two'].hits.values() for gh in th.values())
    region_two_counts = [len(gh.complete_hits) for gh in region_two_hit_gen]
    counts['region_two'] = max(region_two_counts) if region_two_counts else 0
    for gene, gene_hits in databases['three'].hits.items():
        counts[gene] = len(gene_hits.complete_hits)
    return counts


def set_missing(databases, total_loci):
    '''For each Hits instance, set the number of missing genes/types'''
    missing_counts = dict()
    for database in databases.values():
        for hits in gene_hit_generator(database):
            missing_count = total_loci - len(hits.complete_hits)
            if missing_count > 0:
                logging.debug('Missing %s %s', missing_count, hits.name)
                hits.missing = missing_count
                missing_counts[hits.name] = missing_count
    return missing_counts


def find_missing(databases):
    '''For gene Hit instance, attempt to find broken/ truncated genes'''
    for region_name, database in databases.items():
        for hits in gene_hit_generator(database):
            for hit in hits.partial_hits:
                if hits.missing == 0:
                    break
                if hit.length < database.parameters['broken_length_min']:
                    continue
                elif hit.pident < database.parameters['broken_identity_min']:
                    continue
                else:
                    hits.missing -= 1
                    logging.info('Found truncated gene %s', hits.name)
                    hit.broken = True
                    hits.broken_hits.append(hit)


def total_items_found(databases):
    count = 0
    for database in databases.values():
        for hits in gene_hit_generator(database):
            count += len(hits.complete_hits + hits.broken_hits)
    return count


def aggregate_hits_and_sort(databases):
    '''Sort complete and broken hits into contig and then order by position'''
    contig_hits = dict()
    for database in databases.values():
        for gene_hits in gene_hit_generator(database):
            for hit in gene_hits.all_hits:
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
