from . import database
from . import locus


def discover_clusters(hits_complete, hits_remaining, region, contig_sizes, filter_params):
    hits_selected = select_best_genes(hits_complete)
    missing = database.SCHEME[region] - {hit.sseqid for hit in hits}
    if missing:
        all_contig_hits = locus.sort_hits_by_contig(hits_complete)
        for contig_hits in contig_hits:
            hits_broken = discover_missing_genes(hits_selected, hits_remaining, region, filter_params)
            hits_selected |= hits_broken
    return hits_selected


def select_best_genes(hits):
    '''Retaining hits which are a part of region one and three

    Priortise hit with best evalue on the same contig'''
    hits_selected = set()
    for orf_hits in locus.sort_hits_by_orf(hits).values():
        best_hit = max(orf_hits, key=lambda h: (h.evalue, h.orf.near_boundary))
        hits_selected.add(best_hit)
    return hits_selected


def discover_missing_genes(hits, hits_remaining, region, filter_params):
    '''Count missing genes and return any broken hits for region one and three'''
    if not missing:
        return set()
    hits_filtered = database.filter_hits(hits_remaining, **filter_params)
    hits_proximal = locus.find_proximal_hits(hits, hits_filtered, 10000)
    hits_broken_all = {hit for hit in hits_proximal if hit.sseqid in missing}
    broken_hits = select_best_genes(hits_broken_all)
    for hit in broken_hits:
        hit.broken = True
    return broken_hits
