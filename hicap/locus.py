from . import database
from . import region_common
from . import region_specific


class Group:

    def __init__(self, hits, serotypes=set(), contigs=set()):
        self.hits = hits
        self.serotypes = serotypes
        self.contigs = contigs

        if len(self.contigs) <= 1:
            hits_sorted = sorted(self.hits, key=lambda h: h.orf.start)
            self.start = hits_sorted[0].orf.start
            self.end = hits_sorted[-1].orf.end
        else:
            self.start = None
            self.end = None

        self.near_boundary = False


def merge_group(group_1, group_2):
    hits = group_1.hits | group_2.hits
    contigs = group_1.contigs | group_2.contigs
    serotypes = group_1.serotypes | group_2.serotypes

    group_merged = Group(hits, serotypes=serotypes, contigs=contigs)
    group_merged.near_boundary = group_1.near_boundary | group_2.near_boundary
    return group_merged


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
    region_groups = list()
    neighbour_groups = find_neighbours(hits_complete, contig_sizes, 5000)
    for group in sorted(neighbour_groups, key=lambda k: len(k.hits), reverse=True):
        if region in {'one', 'three'}:
            hits_final = region_common.discover_clusters(group.hits, hits_remaining, region, contig_sizes, filter_params)
            serotype = set()
        else:
            #hits_final, serotype = region_specific.discover_clusters(group.hits, hits_remaining, contig_sizes, filter_params)
            continue
        hits_remaining -= hits_final
        contigs = {hit.orf.contig for hit in hits_final}
        region_groups.append(Group(hits_final, serotypes=serotype, contigs=contigs))
    # Link overlapping regions together
    return harmonise_groups(region_groups)


def find_neighbours(hits, contig_sizes, distance):
    '''Find neighbouring hits within a given threshold on the same contig'''
    contig_hits = dict()
    for hit in hits:
        try:
            contig_hits[hit.orf.contig].append(hit)
        except KeyError:
            contig_hits[hit.orf.contig] = [hit]

    # Find neighbours on same contig
    groups = set()
    for contig, hits in contig_hits.items():
        hits_sorted = sorted(hits, key=lambda k: k.orf.start)
        group_hits = [hits_sorted.pop(0)]
        for hit in hits_sorted:
            if (hit.orf.start - group_hits[-1].orf.end) <= distance:
                group_hits.append(hit)
            else:
                groups.add(Group(set(group_hits), contig=contig))
                group = {hit}
        groups.add(Group(set(group_hits), contigs={contig}))

    # If more than one group, check if they're at contig bounds and merge if so
    if len(groups) <= 1:
        return groups
    groups_near_boundary = set()
    for group in groups:
        # Should only ever be only contig at this stage
        if len(group.contigs) <= 1:
            contig = group.contigs.pop()
        else:
            raise ValueError('Invalid number of contigs found during neighbourhood grouping')
        near_boundary = near_contig_boundary(group.start, group.end, contig, contig_sizes, 5000)
        if not near_boundary:
            continue
        groups_near_boundary.add(group)

    groups_merged = groups_near_boundary.pop()
    for group in groups_near_boundary:
        groups_merged = merge_group(groups_merged, group)

    groups.add(groups_merged)
    return groups - groups_near_boundary


def sort_hits_by_orf(hits):
    orfs_hits = dict()
    for hit in hits:
        try:
            orfs_hits[hit.orf].add(hit)
        except KeyError:
            orfs_hits[hit.orf] = {hit}
    return orfs_hits


def sort_hits_by_contig(hits):
    contigs_hits = dict()
    for hit in hits:
        try:
            contigs_hits[hit.orf.contig].add(hit)
        except KeyError:
            contigs_hits[hit.orf.contig] = {hit}
    return contigs_hits


def find_proximal_hits(hits, hits_filtered, distance):
    '''Find hits which are nearby or at contig boundaries'''
    hits_sorted = sorted(hits, key=lambda h: h.orf.start)
    hits_start = hits_sorted[0].orf.start
    hits_end = hits_sorted[-1].orf.end

    hits_proximal = set()
    near_boundary = any(hit.orf.near_boundary for hit in hits)
    contigs = {hit.orf.contig for hit in hits}
    for hit in hits_filtered:
        if hit.orf.contig in contigs:
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
    for group_set_frozen in group_sets:
        group_set = set(group_set_frozen)
        merged = group_set.pop()
        for unmerged in group_set:
            merged = merge_group(merged, unmerged)
        groups_merged.add(merged)
    return groups_merged
