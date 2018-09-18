from . import database
from . import region_common
from . import region_specific


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


def sort_hits_by_orf(hits):
    orfs_hits = dict()
    for hit in hits:
        try:
            orfs_hits[hit.orf].add(hit)
        except KeyError:
            orfs_hits[hit.orf] = {hit}
    return orfs_hits
