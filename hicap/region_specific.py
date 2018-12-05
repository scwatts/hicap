from . import database
from . import locus


NEIGHBOUR_DIST = 5000


def discover_clusters(hits_complete, hits_remaining, filter_params):
    # Find best hits, determine missing genes, find missing genes, and select best
    hits_selected, serotypes = select_best_genes(hits_complete, NEIGHBOUR_DIST)
    genes_missing = dict()
    for serotype in serotypes:
        hits_serotype = {hit for hit in hits_selected if hit.sseqid in database.SEROTYPES[serotype]}
        missing_serotype_genes = locus.count_missing_genes(hits_serotype, database.SEROTYPES[serotype])
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
        # Determine serotype current orf belongs to
        serotype = determine_serotype(orf, orf_hits, distance, orfs_hits)
        serotypes.add(serotype)

        # Get the best hit
        hit_best = perform_selection(orf_hits, serotype)
        hits_selected.add(hit_best)
    return hits_selected, serotypes


def perform_selection(hits, serotype):
    # Keep singular hit
    if len(hits) <= 1:
        return list(hits)[0]
    # Select the orfs first hit of the serotype
    hits_sorted = sorted(hits, key=lambda h: h.evalue)
    for hit in hits_sorted:
        if hit.sseqid in database.SEROTYPES[serotype]:
            return hit
    else:
        # Retain the best hit when no hits of serotype found
        return list(hits_sorted)[0]


def collect_neighbourhood_hits(start, end, contig, orfs_hits):
    orfs_neighbouring = list()
    for orf, orf_hits in orfs_hits.items():
        if contig != orf.contig:
            continue
        if orf.start < start:
            continue
        if orf.end > end:
            continue
        orfs_neighbouring.append(orf_hits)
    return orfs_neighbouring


def determine_serotype(orf, orf_hits, distance, orfs_hits):
    # Check if we have unambiguous hits
    gene_hits = {orf_hit.sseqid for orf_hit in orf_hits}
    if len(gene_hits) == 1:
        return database.get_serotype_group(*gene_hits)

    # Search for any unambiguous hits in the neighbourhood
    start = orf.start - distance
    end = orf.end + distance
    neighbour_orfs_hits = collect_neighbourhood_hits(start, end, orf.contig, orfs_hits)
    for neighbour_orf_hits in neighbour_orfs_hits:
        neighbour_gene_hits = {orf_hit.sseqid for orf_hit in neighbour_orf_hits}
        if len(neighbour_gene_hits) == 1:
            return database.get_serotype_group(*neighbour_gene_hits)

    # See if there are hits anywhere which are unambiguous
    all_unambiguous_st = set()
    for all_orf_hits in orfs_hits.values():
        all_gene_hits = {orf_hit.sseqid for orf_hit in all_orf_hits}
        if len(all_gene_hits) == 1:
            serotype = database.get_serotype_group(*all_gene_hits)
            all_unambiguous_st.add(serotype)
    # Unambiguous hits must only be present for one of the ORF hits
    unambiguous_st = all_unambiguous_st & {database.get_serotype_group(gene) for gene in gene_hits}
    if len(unambiguous_st) == 1:
        return list(unambiguous_st)[0]

    # Make best guess
    neighbourhood_hits = {hit for hits in neighbour_orfs_hits for hit in hits}
    return most_frequent_serotype(neighbourhood_hits)


def most_frequent_serotype(hits):
    counts = {stype: list() for stype in database.SEROTYPES}
    for hit in hits:
        serotype = database.get_serotype_group(hit.sseqid)
        counts[serotype].append(hit)
    most_frequent = max(counts, key=lambda k: len(counts[k]))

    # Check for ties
    most_frequent_serotypes = {st for st in counts if len(counts[st]) == len(counts[most_frequent])}
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
        best_hit = max(orf_hits, key=lambda k: k.bitscore / k.length)
        best_hit_serotype = database.get_serotype_group(best_hit.sseqid)
        serotype_bitscores[best_hit_serotype] += best_hit.bitscore / best_hit.length
    return max(serotype_bitscores, key=lambda k: serotype_bitscores[k])
