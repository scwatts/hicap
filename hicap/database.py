import logging
import tempfile


from . import alignment


SCHEME = {
        'one': {'bexA', 'bexB', 'bexC', 'bexD'},
        'two': {'type_a', 'type_b', 'type_c', 'type_d', 'type_e', 'type_f'},
        'three': {'hcsA', 'hcsB'}
        }


SEROTYPES = {
        'type_a': {'acs1', 'acs2', 'acs3', 'acs4'},
        'type_b': {'bcs1', 'bcs2', 'bcs3', 'bcs4'},
        'type_c': {'ccs1', 'ccs2', 'ccs3', 'ccs4'},
        'type_d': {'dcs1', 'dcs2', 'dcs3', 'dcs4', 'dcs5'},
        'type_e': {'ecs1', 'ecs2', 'ecs3', 'ecs4', 'ecs5', 'ecs6', 'ecs7', 'ecs8'},
        'type_f': {'fcs1', 'fcs2', 'fcs3'},
        }


def search(query_fp, database_fps):
    '''Perform search via alignment of query sequences in provided database files'''
    logging.info('Searching database for matches')
    hits_all = set()
    for database_fp in database_fps:
        with tempfile.TemporaryDirectory() as dh:
            blast_database_fp = alignment.create_blast_database(database_fp, dh)
            blast_stdout = alignment.align(query_fp, blast_database_fp)
            hits_all.update(alignment.parse_blast_stdout(blast_stdout))
    return hits_all


def filter_hits(hits, coverage_min=None, identity_min=None, length_min=None):
    '''Filter hits using provided thresholds'''
    hits_filtered = set()
    for hit in hits:
        if identity_min and hit.pident < (identity_min * 100):
            continue
        if coverage_min and hit.length / hit.slen < coverage_min:
            continue
        if length_min and hit.length < length_min:
            continue
        hits_filtered.add(hit)
    return hits_filtered


def assign_hit_orfs(hits, orfs):
    '''Assign each hit their respective ORF'''
    for hit in hits:
        hit.orf = orfs[int(hit.qseqid)]
    return hits


def get_serotype_group(gene_name):
    '''Return the serotype group of a given region two gene'''
    for serotype, gene_names in SEROTYPES.items():
        if gene_name in gene_names:
            return serotype
    raise ValueError('Could not find %s gene in serotypes' % gene_name)
