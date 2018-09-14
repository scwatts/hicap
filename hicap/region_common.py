from . import database
from . import locus


def discover_clusters(hits_complete, hits_remaining, region, filter_params):
    hits_selected = select_best_hits(hits_complete)
    missing_genes = database.SCHEME[region] - {hit.sseqid for hit in hits_complete}
    filtered_hits = database.filter_hits(hits_remaining, **filter_params)
    missing_hits = {hit for hit in filtered_hits if hit.sseqid in missing_genes}

    if missing_hits:
        broken_hits = select_best_hits(missing_hits)
        for hit in broken_hits:
            hit.broken = True
        hits_selected |= broken_hits
    hits_remaining -= missing_hits
    return hits_selected


def select_best_hits(hits):
    '''Retaining hits which are a part of region one and three

    Priortise hit with best evalue on the same contig'''
    hits_selected = set()
    for orf_hits in locus.sort_hits_by_orf(hits).values():
        best_hit = max(orf_hits, key=lambda h: (h.evalue, h.orf.near_boundary))
        hits_selected.add(best_hit)
    return hits_selected
