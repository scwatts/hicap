import logging
import tempfile


from . import alignment


SCHEME = {
        'one': ('bexA', 'bexB', 'bexC', 'bexD'),
        'two': ('type_a', 'type_b', 'type_c', 'type_d', 'type_e', 'type_f'),
        'three': ('hcsA', 'hcsB')
        }


class Hits():

    def __init__(self, hits):
        self.all = hits
        self.complete = dict()
        self.broken = dict()


    @property
    def remaining(self):
        '''Create dict which contains hits which are unassigned'''
        _remaining = dict()
        for database, hits in self.all.items():
            assigned_hits = set()
            if database in self.complete:
                assigned_hits.update(self.complete[database])
            if database in self.broken:
                assigned_hits.update(self.broken[database])
            _remaining[database] = [hit for hit in hits if hit not in assigned_hits]
        return _remaining


    @property
    def passed(self):
        '''Create a dict for hits that pass filtering for complete or broken genes'''
        _assigned = dict()
        for database in (*self.complete, *self.broken):
            _assigned[database] = list()
            if database in self.complete:
                _assigned[database] += self.complete[database]
            if database in self.broken:
                _assigned[database] += self.broken[database]
        return _assigned


class Locus():

    def __init__(self, contig, orfs, serotype):
        self.contig = contig
        self.orfs = orfs
        self.serotype = serotype


def search(query_fp, database_fps):
    '''Perform search via alignment of query sequences in provided database files'''
    logging.info('Searching database for matches')
    hits_all = {database_fp.stem: list() for database_fp in database_fps}
    for database_fp in database_fps:
        with tempfile.TemporaryDirectory() as dh:
            blast_database_fp = alignment.create_blast_database(database_fp, dh)
            blast_stdout = alignment.align(query_fp, blast_database_fp)
            hits_all[database_fp.stem] = alignment.parse_blast_stdout(blast_stdout)
    return Hits(hits_all)


def filter_hits(hits, coverage_min=None, identity_min=None, length_min=None):
    '''Filter hits using provided thresholds'''
    hits_filtered = {database: list() for database in hits}
    for database, database_hits in hits.items():
        for hit in database_hits:
            if identity_min and hit.pident < identity_min:
                continue
            if coverage_min and hit.length / hit.slen < coverage_min:
                continue
            if length_min and hit.length < length_min:
                continue
            hits_filtered[database].append(hit)
    return hits_filtered


def discover_missing_genes(hits):
    '''Find the names of missing genes'''
    # Count hit names
    region_hits = hit_region_sort(hits)
    counts = dict()
    for hit in (*region_hits['one'], *region_hits['three']):
        try:
            counts[hit.sseqid] += 1
        except KeyError:
            counts[hit.sseqid] = 1

    # Find missing
    expected_count = int(round(sum(counts.values()) / len(counts.values()), 0))
    missing = dict()
    for gene in (*SCHEME['one'], *SCHEME['three']):
        if gene not in counts:
            missing[gene] = expected_count
        elif counts[gene] < expected_count:
            missing[gene] = expected_count  - counts[gene]
    return missing


def select_best_hits(hits, gene_counts):
    '''For a given dict of database hits, select n best hits

    gene_counts is a dict of databases mapping to a count'''
    hits_best = dict()
    for database, database_hits in hits.items():
        if database not in gene_counts:
            continue
        if len(database_hits) <= gene_counts[database]:
            hits_best[database] = database_hits
        else:
            hits_sorted = sorted(database_hits, key=lambda h: h.evalue)
            count = gene_counts[database]
            hits_best[database] = hits_sorted[:count]
    return hits_best


def hit_region_sort(hits):
    '''Sort hits into a dictionary by region'''
    region_hits = {region: list() for region in SCHEME}
    for database, database_hits in hits.items():
        for region, names in SCHEME.items():
            # TODO: does a set have an appreciable improvement?
            if database in names:
                try:
                    region_hits[region].extend(database_hits)
                except KeyError:
                    region_hits[region] = database_hits
    return region_hits


def orf_region_sort(orfs):
    '''Sort ORFs into a dictionary by region'''
    region_orfs = {region: list() for region in SCHEME}
    for orf in orfs:
        for region, names in SCHEME.items():
            if any(name in orf.hits for name in names):
                try:
                    region_orfs[region].append(orf)
                except KeyError:
                    region_orfs[region] = [orf]
    return region_orfs


def match_orfs_and_hits(hits, orfs):
    '''Assign hits to respective ORFs'''
    orf_indices = set()
    for region, hits in hits.items():
        for hit in hits:
            orf_index = int(hit.qseqid)
            orf_indices.add(orf_index)
            try:
                orfs[orf_index].hits[region].append(hit)
            except KeyError:
                orfs[orf_index].hits[region] = [hit]
    return [orfs[orf_index] for orf_index in orf_indices]


def characterise_loci(orfs):
    '''Given a list of ORFs define loci predicated on distance'''
    # Process from largest group to smallest
    groups = find_neighbours(orfs)
    loci = list()
    for contig, group_orfs in sorted(groups, key=lambda k: len(k[1])):
        # Sort ORFs by region and predict region two serotype
        region_orfs = orf_region_sort(group_orfs)
        region_two_types = list()
        for neighbourhood in find_neighbours(region_orfs['two']):
            rtwo_type = predict_region_two_type(region_orfs['two'])
            region_two_types.append(rtwo_type)
        loci.append(Locus(contig, group_orfs, region_two_types))
    return loci


def predict_region_two_type(orfs):
    '''Determine the serotype for each region two group

    Some region two genes from different serotypes share homology. We decide
    region two type by representation comparison'''
    # Count hit regions two types
    counts = dict()
    hit_gen = (hit for orf in orfs for hit in orf.hits)
    for hit in hit_gen:
        try:
            counts[hit] += 1
        except KeyError:
            counts[hit] = 1

    # Select the best represented and remove competing hits from ORFs
    rtwo_type = max(counts, key=lambda k: counts[k])
    for orf in orfs:
        if len(orf.hits) > 1 and rtwo_type in orf.hits:
            orf.hits = {rtwo_type: orf.hits[rtwo_type]}
        else:
            for other_type in sorted(counts, reverse=True, key=lambda k: counts[k]):
                if other_type in orf.hits:
                    orf.hits = {other_type: orf.hits[other_type]}
                    break
    return rtwo_type


def find_neighbours(orfs):
    '''Find neighbouring ORFs within a given threshold on the same contig'''
    contig_orfs = dict()
    for orf in orfs:
        try:
            contig_orfs[orf.contig].append(orf)
        except KeyError:
            contig_orfs[orf.contig] = [orf]

    groups = list()
    for contig, orfs in contig_orfs.items():
        orfs_sorted = sorted(orfs, key=lambda orf: orf.start)
        group = [orfs_sorted.pop(0)]
        for orf in orfs_sorted:
            if (orf.start - group[-1].end) <= 500:
                group.append(orf)
            else:
                groups.append((contig, group))
                group = [orf]
        groups.append((contig, group))
    return groups
