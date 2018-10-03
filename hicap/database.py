import logging
import tempfile


from . import alignment


SCHEME = {
        'one': {'bexA', 'bexB', 'bexC', 'bexD'},
        'two': {'type_a', 'type_b', 'type_c', 'type_d', 'type_e', 'type_f'},
        'three': {'hcsA', 'hcsB'}
        }


SEROTYPES = {
        'type_a': {'ac1', 'ac2', 'ac3', 'ac4'},
        'type_b': {'bc1', 'bc2', 'bc3', 'bc4'},
        'type_c': {'cc1', 'cc2', 'cc3', 'cc4'},
        'type_d': {'dc1', 'dc2', 'dc3', 'dc4', 'dc5'},
        'type_e': {'ec1', 'ec2', 'ec3', 'ec4', 'ec5', 'ec6', 'ec7', 'ec8'},
        'type_f': {'fc1', 'fc2', 'fc3'},
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
        if identity_min and hit.pident < identity_min:
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
