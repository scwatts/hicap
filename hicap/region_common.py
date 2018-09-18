from . import database
from . import locus


def discover_clusters(hits_complete, hits_remaining, region, filter_params):
    hits_selected = select_best_hits(hits_complete)
    genes_missing = database.SCHEME[region] - {hit.sseqid for hit in hits_complete}
    hits_filtered = database.filter_hits(hits_remaining, **filter_params)
    hits_missing = {hit for hit in hits_filtered if hit.sseqid in genes_missing}
    # Select best hits for each discovered missing ORF
    if hits_missing:
        hits_broken = select_best_hits(hits_missing)
        for hit in hits_broken:
            hit.broken = True
        hits_selected |= hits_broken

    # Create and return locus.Group
    hits_remaining -= hits_selected
    contigs = {hit.orf.contig for hit in hits_selected}
    return locus.Group(hits_selected, contigs=contigs)


def select_best_hits(hits):
    '''Retaining hits which are a part of region one and three'''
    hits_selected = set()
    for orf_hits in locus.sort_hits_by_orf(hits).values():
        best_hit = max(orf_hits, key=lambda h: (h.evalue, h.bitscore))
        hits_selected.add(best_hit)
    return hits_selected
