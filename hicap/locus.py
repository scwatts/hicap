from . import database
from . import region_common
from . import region_specific


RTWO_FLANK_DIST = 5000
NEARBY_FLANK_DIST = 1000


class Group:

    def __init__(self, hits, *, serotypes=None, contigs=None):
        self.hits = hits
        self.serotypes = serotypes
        self.contigs = contigs

        if self.contigs and len(self.contigs) <= 1:
            hits_sorted = sorted(self.hits, key=lambda h: h.orf.start)
            self.start = hits_sorted[0].orf.start
            self.end = hits_sorted[-1].orf.end
        else:
            self.start = None
            self.end = None


def get_gene_region(gene_name):
    '''Determine region a gene belongs to given the gene name'''
    for region, region_genes in database.SCHEME.items():
        if gene_name in region_genes:
            return region
    for serotype, region_two_genes in database.SEROTYPES.items():
        if gene_name in region_two_genes:
            return 'two'
    else:
        raise ValueError('Could not find %s gene in database scheme' % gene_name)


def discover_region_clusters(hits_complete, hits_remaining, region, filter_params):
    if region in {'one', 'three'}:
        return region_common.discover_clusters(hits_complete, hits_remaining, region, filter_params)
    else:
        return region_specific.discover_clusters(hits_complete, hits_remaining, filter_params)


def count_missing_genes(hits, expected_genes):
    hit_counts = dict.fromkeys(expected_genes, 0)
    for hit in hits:
        hit_counts[hit.sseqid] += 1
    expected_count = max(hit_counts.values())

    # In the absence of complete hits, attempt to find nonetheless
    if not expected_count:
        expected_count = 1

    missing_count = dict()
    for hit, count in hit_counts.items():
        missing_count[hit] = expected_count - count
    return missing_count


def collect_missing_genes(hits, genes_missing):
    hits_missing = set()
    gene_hits = sort_hits_by_gene(hits)
    for gene, count in genes_missing.items():
        if gene not in gene_hits:
            continue
        hits_sorted = sorted(gene_hits[gene], key=lambda k: (1-k.evalue, k.bitscore), reverse=True)
        hits_missing.update(hits_sorted[:count])
    return hits_missing


def locate_fragmented_region_two(groups, hits_remaining, filter_params):
    # Collect all possible hits for region two
    genes_rtwo_all = {gene for genes in database.SEROTYPES.values() for gene in genes}
    hits_rtwo_all = {hit for hit in hits_remaining if hit.sseqid in genes_rtwo_all}
    hits_rtwo_filtered = database.filter_hits(hits_rtwo_all, **filter_params)
    if not hits_rtwo_filtered:
        return Group({})

    # Hits upstream and downstream of region one and three
    hits_candidate = set()
    for region in ('one', 'three'):
        for contig, contig_hits in sort_hits_by_contig(groups[region].hits).items():
            getter = lambda h: h.orf
            hits_start, hits_end = get_elements_bounds(contig_hits, getter)
            range_start = hits_start - RTWO_FLANK_DIST
            range_end = hits_end + RTWO_FLANK_DIST
            hits_candidate |= collect_elements_in_bounds(range_start, range_end, contig, hits_rtwo_filtered, getter)

    # Select best hits and set them to broken
    hits_remaining -= hits_candidate
    group = region_specific.discover_clusters(hits_candidate, hits_remaining, filter_params)
    for hit in group.hits:
        hit.broken = True
    return group


def find_adjacent_fragments(hits, region, hits_remaining):
    hits_orfs = {hit.orf for hit in hits}
    hits_adjacent = set()
    for hit in hits:
        # Only search around broken hits
        if not hit.broken:
            continue

        # Find all adjacent hits with the same subject gene
        if hit.sstart < hit.send:
            start = hit.orf.start - hit.sstart + 1
            end = hit.orf.end + (hit.slen - hit.send) + 1
        else:
            start = hit.orf.start - (hit.slen - hit.send) + 1
            end = hit.orf.end + hit.sstart + 1
        getter = lambda h: h.orf
        hits_collected = collect_elements_in_bounds(start, end, hit.orf.contig, hits_remaining, getter)
        hits_collected = {h for h in hits_collected if h.sseqid == hit.sseqid and h.orf not in hits_orfs}

        # Apply some sanity filtering here - not exposed to user
        hits_filtered = set()
        for hit in hits_collected:
            if hit.evalue >= 0.01:
                continue
            if hit.length <= 60:
                continue
            hits_filtered.add(hit)

        # Select best hits - Just being careful here; in most cases this will be unnecessary
        if region in {'one', 'three'}:
            hits_selected = region_common.select_best_hits(hits_filtered)
        else:
            hits_selected = set()
            serotype = database.get_serotype_group(hit.sseqid)
            for orf_hits in sort_hits_by_orf(hits_filtered).values():
                hit_best = region_specific.perform_selection(orf_hits, serotype)
                hits_selected.add(hit_best)
        hits_adjacent |= hits_selected

    # Update hits_remaining and return
    hits_remaining -= hits_adjacent
    for hit in hits_adjacent:
        hit.broken = True
    return hits_adjacent


def collect_nearby_orfs(region_groups, orfs_all):
    hits_selected = {hit for group in region_groups.values() for hit in group.hits}
    orfs_selected = sort_hits_by_orf(hits_selected)
    orfs_remaining = set(orfs_all) - set(orfs_selected)

    nearby_orfs = set()
    for contig, contig_hits in sort_hits_by_contig(hits_selected).items():
        getter_hit = lambda h: h.orf
        hits_start, hits_end = get_elements_bounds(contig_hits, getter_hit)
        range_start = hits_start - NEARBY_FLANK_DIST
        range_end = hits_end + NEARBY_FLANK_DIST
        # TODO: we should check if these nearby ORFs are _somehow_ cap locus specific genes
        orfs = collect_elements_in_bounds(range_start, range_end, contig, orfs_remaining)
        # Apply some sanity filtering here - not exposed to user
        orfs_filtered = set()
        for orf in orfs:
            if (orf.end - orf.start) <= 100:
                continue
            orfs_filtered.add(orf)
        nearby_orfs |= orfs_filtered
    return nearby_orfs


def get_elements_bounds(elements, getter):
    # This getter confusion saves ~60 lines of repeating code. worth it? that's highly doubtful
    elements_sorted = sorted(elements, key=lambda k: getter(k).start)
    start = min(getter(elements_sorted[0]).start, getter(elements_sorted[0]).end)
    end = max(getter(elements_sorted[-1]).start, getter(elements_sorted[-1]).end)
    return start, end


def collect_elements_in_bounds(start, end, contig, elements, getter=lambda e: e):
    # TODO: optimise if required
    elements_selected = set()
    for element in elements:
        if getter(element).contig != contig:
            continue
        element_start = min(getter(element).start, getter(element).end)
        element_end = max(getter(element).start, getter(element).end)
        if start <= element_start <= end:
            elements_selected.add(element)
        elif start <= element_end <= end:
            elements_selected.add(element)
    return elements_selected


def sort_hits_by_orf(hits):
    orfs_hits = dict()
    for hit in hits:
        try:
            orfs_hits[hit.orf].add(hit)
        except KeyError:
            orfs_hits[hit.orf] = {hit}
    return orfs_hits


def sort_hits_by_gene(hits):
    gene_hits = dict()
    for hit in hits:
        try:
            gene_hits[hit.sseqid].add(hit)
        except KeyError:
            gene_hits[hit.sseqid] = {hit}
    return gene_hits


def sort_hits_by_contig(hits):
    contigs_hits = dict()
    for hit in hits:
        try:
            contigs_hits[hit.orf.contig].add(hit)
        except KeyError:
            contigs_hits[hit.orf.contig] = {hit}
    return contigs_hits


def sort_hits_by_region(hits):
    region_hits = {region: list() for region in database.SCHEME}
    for hit in hits:
        if not hit.region:
            hit.region = get_gene_region(hit.sseqid)
        region_hits[hit.region].append(hit)
    return region_hits


def sort_orfs_by_contig(orfs):
    contigs_orfs = dict()
    for orf in orfs:
        try:
            contigs_orfs[orf.contig].add(orf)
        except KeyError:
            contigs_orfs[orf.contig] = {orf}
    return contigs_orfs
