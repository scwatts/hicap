from . import database
from . import locus


NEIGHBOUR_DIST = 5000


def discover_clusters(hits_complete, hits_remaining, filter_params):
    # Find best hits, determine missing genes, find missing genes, and select best
    hits_selected, serotypes = select_best_genes(hits_complete, NEIGHBOUR_DIST)
    genes_missing = dict()
    for serotype in serotypes:
        missing_serotype_genes = locus.count_missing_genes(hits_selected, database.SEROTYPES[serotype])
        genes_missing.update(missing_serotype_genes)
    hits_filtered = database.filter_hits(hits_remaining, **filter_params)
    hits_missing = locus.collect_missing_genes(hits_filtered, genes_missing)
    if hits_missing:
        for hit in hits_missing:
            hit.broken = True
        hits_candidate = hits_selected | hits_missing
        hits_selected, serotypes = select_best_genes(hits_candidate, NEIGHBOUR_DIST)

    # Create and return locus.Group
    hits_remaining -= hits_selected
    contigs = {hit.orf.contig for hit in hits_selected}
    return locus.Group(hits_selected, serotypes=serotypes, contigs=contigs)


def select_best_genes(hits, distance):
    hits_selected = set()
    serotypes = set()
    orfs_hits = locus.sort_hits_by_orf(hits)
    for orf, orf_hits in orfs_hits.items():
        # TODO: efficiency could be improved
        # Get most frequent serotype hit within +/-5 kb of this orf
        start = orf.start - distance
        end = orf.end + distance
        neighbourhood_hits = collect_neighbourhood_hits(start, end, orf.contig, orfs_hits)
        serotype = most_frequent_serotype(neighbourhood_hits)
        serotypes.add(serotype)

        hits_sorted = sorted(orf_hits, key=lambda h: h.evalue)
        # Keep singular hit
        if len(hits_sorted) <= 1:
            hits_selected.update(hits_sorted)
            continue

        # Select the hits first of the serotype
        for hit in hits_sorted:
            if hit.sseqid in database.SEROTYPES[serotype]:
                hits_selected.add(hit)
                break
        else:
            # Retain the best hit when no hits of serotype found
            hits_selected.add(list(hits_sorted)[0])
    return hits_selected, serotypes


def collect_neighbourhood_hits(start, end, contig, orfs_hits):
    orfs_neighbouring = set()
    for orf, orf_hits in orfs_hits.items():
        if contig != orf.contig:
            continue
        if orf.start < start:
            continue
        if orf.end > end:
            continue
        orfs_neighbouring.update(orf_hits)
    return orfs_neighbouring


def most_frequent_serotype(hits):
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
        return break_most_frequent_type_tie(counts, most_frequent_serotypes)


def break_most_frequent_type_tie(counts, most_frequent_serotypes):
    # Break by highest accumulative bitscore
    # This assumes that at least one ORF has n hits (where n is the number of serotypes)
    orf_hits = dict()
    for serotype in most_frequent_serotypes:
        for hit in counts[serotype]:
            try:
                orf_hits[hit.orf].append(hit)
            except KeyError:
                orf_hits[hit.orf] = [hit]
    serotype_bitscores = {serotype: 0 for serotype in most_frequent_serotypes}
    for orf, orf_hits in orf_hits.items():
        if len(orf_hits) <= 1:
            continue
        best_hit = max(orf_hits, key=lambda k: k.bitscore)
        best_hit_serotype = database.get_serotype_group(best_hit.qseqid)
        serotype_bitscores[best_hit_serotype] += 1
    # Ignore any breakable ties
    return max(serotype_bitscores, key=lambda k: len(serotype_bitscores[k]))
