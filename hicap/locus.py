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


def discover_region_clusters(hits_complete, hits_remaining, region, filter_params):
    if region in {'one', 'three'}:
        hits_final = region_common.discover_clusters(hits_complete, hits_remaining, region, filter_params)
        serotype = set()
    else:
        hits_final, serotype = region_specific.discover_clusters(group.hits, hits_remaining, contig_sizes, filter_params)
    hits_remaining -= hits_final
    contigs = {hit.orf.contig for hit in hits_final}
    return Group(hits_final, serotypes=serotype, contigs=contigs)


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
