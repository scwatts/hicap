from . import database
from . import region_common
from . import region_specific


FLANK_DIST = 5000


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


def discover_region_clusters(hits_complete, hits_remaining, region, filter_params):
    if region in {'one', 'three'}:
        return region_common.discover_clusters(hits_complete, hits_remaining, region, filter_params)
    else:
        return region_specific.discover_clusters(hits_complete, hits_remaining, filter_params)


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
            hits_sorted = sorted(contig_hits, key=lambda k: k.orf.start)
            start = min(hits_sorted[0].orf.start, hits_sorted[0].orf.end)
            end = max(hits_sorted[-1].orf.start, hits_sorted[-1].orf.end)
            hits_candidate |= collect_hits_in_bounds(start, end, contig, FLANK_DIST, hits_rtwo_filtered)

    # Select best hits and set them to broken
    hits_remaining -= hits_candidate
    group = region_specific.discover_clusters(hits_candidate, hits_remaining, filter_params)
    for hit in group.hits:
        hit.broken = True
    return group


def collect_hits_in_bounds(start, end, contig, distance, hits):
    # TODO: optimise if required
    hits_selected = set()
    bounds = ((start - distance, start), (end, end + distance))
    for hit in hits:
        if hit.orf.contig != contig:
            continue
        hit_start = min(hit.orf.start, hit.orf.end)
        hit_end = max(hit.orf.start, hit.orf.end)
        for bound_start, bound_end in bounds:
            if  hit_start >= bound_start and hit_start <= bound_end:
                hits_selected.add(hit)
            elif  hit_end >= bound_start and hit_end <= bound_end:
                hits_selected.add(hit)
    return hits_selected


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
