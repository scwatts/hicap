from . import database
from . import locus


def discover_clusters(hits_complete, hits_remaining, contig_sizes, filter_params):
    serotype = most_frequent_type(hits_complete)
    hits_selected = select_best_genes(hits_complete, serotype)
    hits_broken = discover_missing_genes(hits_selected, hits_remaining, serotype, filter_params)
    return hits_selected | hits_broken, {serotype}


def select_best_genes(hits, serotype):
    '''Retaining hits which are a part of the most frequent serotype

    Fallback to other genes if none are present for a given ORF'''
    # Sort by ORF and keep hits which are most likely
    hits_selected = hits.copy()
    for orf_hits in locus.sort_hits_by_orf(hits).values():
        hits_sorted = sorted(orf_hits, key=lambda h: (h.evalue, h.orf.near_boundary))
        # Delete none
        if len(orf_hits) <= 1:
            continue
        # Delete all but first
        if not any(hit.sseqid in database.SEROTYPES[serotype] for hit in hits_sorted):
            for hit in list(hits)[1:]:
                hits_selected.remove(hit)
            continue
        # Delete all but first gene of most frequent serotype
        found_hit = False
        for hit in hits_sorted:
            if found_hit or hit.sseqid not in database.SEROTYPES[serotype]:
                hits_selected.remove(hit)
                continue
            found_hit = True
    return hits_selected


def most_frequent_type(hits):
    '''Return the most serotype which is most complete'''
    counts = {stype: list() for stype in database.SEROTYPES}
    for hit in hits:
        serotype = database.get_serotype_group(hit.sseqid)
        counts[serotype].append(hit)
    most_frequent = max(counts, key=lambda k: len(counts[k]))

    # Check for ties
    most_frequent_serotypes = {st for st in counts if counts[st] == counts[most_frequent]}
    if len(most_frequent_serotypes) == 1:
        return most_frequent
    else:
        return break_most_frequent_type_tie(hits, counts, most_frequent_serotypes)


def break_most_frequent_type_tie(hits, counts, most_frequent_serotypes):
    # Break by highest accumlative bitscore
    # This assumes that at least one ORF has n hits (where n is the number of serotypes)
    orf_hits = dict()
    for serotype in most_frequent_serotypes:
        for hit in counts[serotype]:
            try:
                org_hits[hit.orf].append(hit)
            except KeyError:
                org_hits[hit.orf] [hit]
    serotype_bitscores = {serotype: 0 for serotype in most_frequent_serotypes}
    for orf, orf_hits in orf_hits.items():
        if len(orf_hits) <= 1:
            continue
        best_hit = max(orf_hits, key=lambda k: k.bitscore)
        best_hit_serotype = database.get_serotype_group(hit.qseqid)
        serotype_bitscores[best_hit_serotype] += 1
    # Ignore any breakable ties
    return max(serotype_bitscores, key=lambda k: len(serotype_bitscores[k]))


def discover_missing_genes(hits_complete, hits_remaining, serotype, filter_params):
    '''Count missing genes and return any broken hits for region two

    Filter, locate nearby hits, and select best for each respective ORF'''
    missing = database.SEROTYPES[serotype] - {hit.sseqid for hit in hits_complete}
    if not missing:
        return set()
    hits_filtered = database.filter_hits(hits_remaining, **filter_params)
    hits_proximal = locus.find_proximal_hits(hits_complete, hits_filtered, 1000)
    hits_broken_all = {hit for hit in hits_proximal if hit.sseqid in missing}
    broken_hits = select_best_genes(hits_broken_all, serotype)
    for hit in broken_hits:
        hit.broken = True
    return broken_hits
