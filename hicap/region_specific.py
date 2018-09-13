def select_best_region_two_genes(hits, serotype):
    '''Retaining hits which are a part of the most frequent serotype

    Fallback to other genes if none are present for a given ORF'''
    # Sort by ORF and keep hits which are most likely
    hits_selected = hits.copy()
    for orf_hits in sort_hits_by_orf(hits).values():
        # Delete none
        if len(orf_hits) <= 1:
            continue
        # Delete all but first
        if not any(hit.sseqid in SEROTYPES[serotype] for hit in orf_hits):
            for hit in list(hits)[1:]:
                hits_selected.remove(hit)
            continue
        # Delete all but first gene of most frequeny serotype
        found_hit = False
        for hit in orf_hits:
            if found_hit or hit.sseqid not in SEROTYPES[serotype]:
                hits_selected.remove(hit)
                continue
            found_hit = True
    return hits_selected
