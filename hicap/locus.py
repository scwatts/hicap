from . import database
from . import region_common
from . import region_specific


class Group:

    def __init__(self, hits, serotypes=set(), contigs=set()):
        self.hits = hits
        self.serotypes = serotypes
        self.contigs = contigs

        hits_sorted = sorted(self.hits, key=lambda h: h.orf.start)
        self.start = hits_sorted[0].orf.start
        self.end = hits_sorted[-1].orf.end

        self.near_boundary = False


def identify_orfs_near_boundaries(orfs, contig_sizes, distance):
    for orf in orfs:
        if orf.start <= distance:
            orf.near_boundary = True
        elif orf.end >= (contig_sizes[orf.contig] - distance):
            orf.near_boundary = True
    return orfs


def near_contig_boundary(start, end, contig, contig_sizes, distance):
    return (start - distance) <= 0 or (contig_sizes[contig] - end) <= distance


def region_sort_hits(hits):
    '''Sort hits into a dictionary by region'''
    region_hits = {region: list() for region in database.SCHEME}
    for hit in hits:
        if not hit.region:
            hit.region = get_gene_region(hit.sseqid)
        region_hits[hit.region].append(hit)
    return region_hits


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


def discover_region_clusters(hits_complete, hits_remaining, region, contig_sizes, filter_params):
    groups = list()
    neighbour_hits = find_neighbours(hits_complete, 1000)
    for contig, hits_list in sorted(neighbour_hits, key=lambda k: len(k[1])):
        hits = set(hits_list)
        if region in {'one', 'three'}:
            hits_final = region_common.discover_clusters(hits_complete, hits_remaining, region, contig_sizes, filter_params)
        else:
            # TODO: temp
            continue
            hits_final = region_specific.discover_clusters(hits_complete, hits_remaining, contig_sizes, filter_params)
        hits_remaining ^= hits_final
        # TODO: generalise argument passing
        contigs = {hit.orf.contig for hit in hits_final}
        groups.append(Group(hits_final, contigs=contigs))
    # Link overlapping regions together
    groups = harmonise_groups(groups)
    return remove_singletons(groups, contig_sizes)


def find_neighbours(hits, distance):
    '''Find neighbouring hits within a given threshold on the same contig'''
    contig_hits = dict()
    for hit in hits:
        try:
            contig_hits[hit.orf.contig].append(hit)
        except KeyError:
            contig_hits[hit.orf.contig] = [hit]

    # TODO: use group class here?
    groups = list()
    for contig, hits in contig_hits.items():
        hits_sorted = sorted(hits, key=lambda k: k.orf.start)
        group = [hits_sorted.pop(0)]
        for hit in hits_sorted:
            if (hit.orf.start - group[-1].orf.end) <= 500:
                group.append(hit)
            else:
                groups.append((contig, group))
                group = [hit]
        groups.append((contig, group))
    return groups


def sort_hits_by_orf(hits):
    orfs_hits = dict()
    for hit in hits:
        try:
            orfs_hits[hit.orf].add(hit)
        except KeyError:
            orfs_hits[hit.orf] = {hit}
    return orfs_hits


def find_proximal_hits(hits, hits_filtered, distance):
    '''Find hits which are nearby or at contig boundaries'''
    hits_sorted = sorted(hits, key=lambda h: h.orf.start)
    hits_start = hits_sorted[0].orf.start
    hits_end = hits_sorted[-1].orf.end

    hits_proximal = set()
    near_boundary = any(hit.orf.near_boundary for hit in hits)
    for hit in hits_filtered:
        if hit.orf.end >= (hits_start - distance):
            hits_proximal.add(hit)
        elif hit.orf.start <= (hits_end + distance):
            hits_proximal.add(hit)
        elif near_boundary and hit.orf.near_boundary:
            hits_proximal.add(hit)
    return hits_proximal


def harmonise_groups(groups):
    # Aggregate groups by ORF
    orf_map = dict()
    for group in groups:
        for hit in group.hits:
            try:
                orf_map[hit.orf].add(group)
            except KeyError:
                orf_map[hit.orf] = {group}

    # Create mapping of all linked groups
    group_map = {group: {group} for group in groups}
    for orf, groups in orf_map.items():
        if len(groups) <= 1:
            continue
        # This could be far more efficient but gains would be ultimately minor
        # Lazily get all elements, merge and reassign
        current_groups = {g for group in groups for g in group_map[group]}
        for group in current_groups:
            group_map[group] = current_groups

    # Create new groups
    groups_merged = set()
    group_sets = {frozenset(group_set) for group_set in group_map.values()}
    for group_set in group_sets:
        hits = {hit for group in group_set for hit in group.hits}
        contigs = {contig for group in group_set for contig in group.contigs}
        serotypes = {stype for group in group_set for stype in group.serotypes}

        group_merged = Group(hits, serotypes=serotypes, contigs=contigs)
        group_merged.near_boundary = any(group.near_boundary for group in group_set)
        groups_merged.add(group_merged)
    return groups_merged


def remove_singletons(groups, contig_sizes):
    # TODO: we won't have to check for contigs, they'll be merged into a single group
    # Check if loci may be split across contigs
    for group in groups:
        for contig in group.contigs:
            # We have code that is so W I D E
            near_boundary = near_contig_boundary(group.start, group.end, group.contigs.pop(), contig_sizes, 1000)
            if near_boundary:
                group.near_boundary = True
                break

    # Allow singletons if there's at least two groups near contig boundaries and be one itself
    allow_singletons = sum(g.near_boundary for g in groups) > 1
    groups_final = set()
    for group in groups:
        if len(group.hits) > 1:
            groups_final.add(group)
        elif allow_singletons and group.near_boundary:
            groups_final.add(group)
    return groups_final
