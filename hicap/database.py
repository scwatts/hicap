import logging
import tempfile


from . import alignment


SCHEME = {
        'one': {'bexA', 'bexB', 'bexC', 'bexD'},
        'two': {'type_a', 'type_b', 'type_c', 'type_d', 'type_e', 'type_f'},
        'three': {'hcsA', 'hcsB'}
        }


SEROTYPES = {
        'type_a': {'ac1', 'ac2', 'ac3', 'ac4'},
        'type_b': {'bc1', 'bc2', 'bc3', 'bc4'},
        'type_c': {'cc1', 'cc2', 'cc3'},
        'type_d': {'dc1', 'dc2', 'dc3', 'dc4', 'dc5'},
        'type_e': {'ec1', 'ec2', 'ec3', 'ec4', 'ec5', 'ec6', 'ec7', 'ec8'},
        'type_f': {'fc1', 'fc2', 'fc3'},
        }


def search(query_fp, database_fps):
    '''Perform search via alignment of query sequences in provided database files'''
    logging.info('Searching database for matches')
    hits_all = set()
    for database_fp in database_fps:
        with tempfile.TemporaryDirectory() as dh:
            blast_database_fp = alignment.create_blast_database(database_fp, dh)
            blast_stdout = alignment.align(query_fp, blast_database_fp)
            hits_all.update(alignment.parse_blast_stdout(blast_stdout))
    return hits_all


def filter_hits(hits, coverage_min=None, identity_min=None, length_min=None):
    '''Filter hits using provided thresholds'''
    hits_filtered = set()
    for hit in hits:
        if identity_min and hit.pident < identity_min:
            continue
        if coverage_min and hit.length / hit.slen < coverage_min:
            continue
        if length_min and hit.length < length_min:
            continue
        hits_filtered.add(hit)
    return hits_filtered


def assign_hit_orfs(hits, orfs):
    '''Assign each hit their respective ORF'''
    for hit in hits:
        hit.orf = orfs[int(hit.qseqid)]
    return hits





def discover_region_two_clusters(hits, hits_remaining, contig_sizes, filter_params):
    # Select most likely genes and find nearby broken genes
    # TODO: use group class instead?
    region_hits = region_sort_hits(hits)
    neighbour_hits = find_neighbours(region_hits['two'], 1000)
    groups = list()
    for contig, hits_list in sorted(neighbour_hits, key=lambda k: len(k[1])):
        hits = set(hits_list)
        serotype = most_frequent_type(hits)
        hits_selected = select_best_region_two_genes(hits, serotype)
        hits_broken = discover_missing_region_two(hits_selected, hits_remaining, serotype, filter_params)
        hits_final = hits_selected | hits_broken
        hits_remaining ^= hits_final
        groups.append(Group(hits_final, serotypes=serotypes, contig=contig))
    return remove_singletons(groups, contig_sizes)







def most_frequent_type(hits):
    '''Return the most serotype which is most complete'''
    counts = {stype: list() for stype in SEROTYPES}
    for hit in hits:
        serotype = get_serotype_group(hit.sseqid)
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
        best_hit_serotype = get_serotype_group(hit.qseqid)
        serotype_bitscores[best_hit_serotype] += 1
    # Ignore any breakable ties
    return max(serotype_bitscores, key=lambda k: len(serotype_bitscores[k]))


def get_serotype_group(gene_name):
    '''Return the serotype group of a given region two gene'''
    for serotype, gene_names in SEROTYPES.items():
        if gene_name in gene_names:
            return serotype
    raise ValueError('Could not find %s gene in serotypes' % gene_name)


def discover_missing_region_two(hits, hits_remaining, serotype, filter_params):
    '''Count missing genes and return any broken hits for region two

    Filter, locate nearby hits, and select best for each respective ORF'''
    missing = SEROTYPES[serotype] - {hit.sseqid for hit in hits}
    if not missing:
        return
    hits_filtered = filter_hits(hits_remaining, **filter_params)
    hits_proximal = find_proximal_hits(hits, hits_filtered, 1000)
    hits_broken_all = {hit for hit in hits_proximal if hit.sseqid if SEROTYPES[serotype]}
    broken_hits = select_best_region_two_genes(hits_broken_all, serotype)
    for hit in broken_hits:
        hit.broken = True
    return broken_hits



#def discover_missing_genes_flank(hits_remaining, missing_genes, filter_params):
#    '''Find the names of missing genes in flanking regions (region one and three)'''
#    hits_filtered = filter_hits(hits_remaining, **filter_params)
#    broken_hits = select_best_hits(hits_filtered, missing_genes)
#    for hit in broken_hits:
#        hit.broken = True
#    return broken_hits
#
#
#def discover_missing_region_two(hits, broken_gene_identity, broken_gene_length):
#    '''Find the names of missing genes in region two'''
#    pass
#
#
#def count_missing_genes_flank(hits):
#    # Count hit names
#    region_hits = region_sort_hits(hits)
#    counts = dict()
#    for hit in (*region_hits['one'], *region_hits['three']):
#        try:
#            counts[hit.sseqid] += 1
#        except KeyError:
#            counts[hit.sseqid] = 1
#
#    if len(counts.values()) < 1:
#        return dict()
#
#    # Find missing
#    expected_count = int(round(sum(counts.values()) / len(counts.values()), 0))
#    missing = dict()
#    for gene in (*SCHEME['one'], *SCHEME['three']):
#        if gene not in counts:
#            missing[gene] = expected_count
#        elif counts[gene] < expected_count:
#            missing[gene] = expected_count - counts[gene]
#    return missing
#
#
#def select_best_hits(hits, gene_counts):
#    '''For a given dict of database hits, select n best hits
#
#    gene_counts is a dict of databases mapping to a count'''
#    # Aggregate hits by sseqid
#    gene_hits = dict()
#    for hit in hits:
#        try:
#            gene_hits[hit.sseqid].append(hit)
#        except KeyError:
#            gene_hits[hit.sseqid] = [hit]
#
#    # Select best hit for each gene
#    hits_best = set()
#    for gene, hits in gene_hits.items():
#        if gene not in gene_counts:
#            continue
#        if len(hits) <= gene_counts[gene]:
#            hits_best.update(hits)
#        else:
#            hits_sorted = sorted(hits, key=lambda h: h.evalue)
#            count = gene_counts[gene]
#            hits_best.update(hits_sorted[:count])
#    return hits_best
#
#


#def orf_region_sort(orfs):
#    '''Sort ORFs into a dictionary by region'''
#    region_orfs = {region: list() for region in SCHEME}
#    for orf in orfs:
#        for region, names in SCHEME.items():
#            if any(name in orf.hits for name in names):
#                try:
#                    region_orfs[region].append(orf)
#                except KeyError:
#                    region_orfs[region] = [orf]
#    return region_orfs
#
#
#def match_orfs_and_hits(hits, orfs):
#    '''Assign hits to respective ORFs'''
#    orfs_assigned = set()
#    for hit in hits:
#        orf = orfs[int(hit.qseqid)]
#        orf.hits.add(hit)
#        # If not in set, add this ORF. Future hits update ref in all places
#        orfs_assigned.add(orf)
#    return orfs_assigned
#
#
#def characterise_loci(orfs, contig_sizes):
#    '''Given a list of ORFs define loci predicated on distance'''
#    # Process from largest group to smallest
#    groups = find_neighbours(orfs)
#    loci_blocks = list()
#    for contig, group_orfs in sorted(groups, key=lambda k: len(k[1])):
#        # Sort ORFs by region and predict region two serotype
#        region_orfs = orf_region_sort(group_orfs)
#        region_two_types = list()
#        for _, neighbourhood in find_neighbours(region_orfs['two']):
#            # Attempt to find missing or broken genes for region two
#            rtwo_type, rtwo_complete = predict_region_two_type(neighbourhood)
#            region_two_types.append(rtwo_type)
#        loci_blocks.append(Block(contig, group_orfs, region_two_types))
#    return loci_blocks
#
#
#def characterise_loci_region_two(orfs):
#    region_orfs = orf_region_sort(group_orfs)
#    region_two_groups = find_neighbours(region_orfs['two'])
#    results = list()
#    for group in region_two_groups:
#        rtwo_type, missing = predict_region_two_type(group)
#        results.append((rtwo_type, missing))
#    return results
#
#
#def predict_region_two_type(orfs):
#    '''Determine the serotype for each region two group
#
#    Some region two genes from different serotypes share homology. We decide
#    region two type by representation comparison'''
#    # Count hit regions two types
#    database_hits = dict()
#    #hit_gen = (hit for orf in orfs for hit in orf.hits)
#    hits_gen = (kv for orf in orfs for kv in orf.hits.items())
#    for (database, hits) in hits_gen:
#        for hit in hits:
#            try:
#                database_hits[database].append(hit.sseqid)
#            except KeyError:
#                database_hits[database] = [hit.sseqid]
#
#    # Select the best represented and remove competing hits from ORFs
#    # TODO: provide additional information - parital, w and wo borked genes, complete
#    rtwo_type = max(database_hits, key=lambda k: len(database_hits[k]))
#    for orf in orfs:
#        if len(orf.hits) > 1 and rtwo_type in orf.hits:
#            orf.hits = {rtwo_type: orf.hits[rtwo_type]}
#        else:
#            for other_type in sorted(database_hits, reverse=True, key=lambda k: database_hits[k]):
#                if other_type in orf.hits:
#                    orf.hits = {other_type: orf.hits[other_type]}
#                    break
#    return rtwo_type, SEROTYPES[rtwo_type] - set(database_hits[rtwo_type])
#
#
