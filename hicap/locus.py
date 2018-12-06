import pathlib
import tempfile


from . import annotation
from . import database
from . import region_common
from . import region_specific


RTWO_FLANK_DIST = 5000
NEARBY_FLANK_DIST = 1000


class LocusData:

    def __init__(self):
        self.regions = dict()
        self.is_hits = None
        self.nearby_orfs = None


class Region:

    def __init__(self, orf_hits, *, serotypes=None, contigs=None):
        self.orf_hits = orf_hits
        self.serotypes = serotypes
        self.contigs = contigs
        self.blast_hits = set()

        if self.contigs and len(self.contigs) <= 1:
            hits_sorted = sorted(self.orf_hits, key=lambda h: h.orf.start)
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
        return Region({})

    # Hits upstream and downstream of region one and three
    hits_candidate = set()
    for region in ('one', 'three'):
        for contig, contig_hits in sort_hits_by_contig(groups[region].orf_hits).items():
            hits_start, hits_end = get_elements_bounds(contig_hits)
            range_start = hits_start - RTWO_FLANK_DIST
            range_end = hits_end + RTWO_FLANK_DIST
            hits_candidate |= collect_elements_in_bounds(range_start, range_end, contig, hits_rtwo_filtered)

    # Select best hits and set them to broken
    hits_remaining -= hits_candidate
    group = region_specific.discover_clusters(hits_candidate, hits_remaining, filter_params)
    for hit in group.orf_hits:
        hit.broken = True
    return group


def find_proximal_fragments(region_groups, hits_remaining, contig_fasta):
    # Get proximal distances
    hits = {hit for region_data in region_groups.values() for hit in region_data.orf_hits}
    contig_ranges = get_proximal_ranges(hits, contig_fasta)

    # Find hits within the locus ranges
    hit_orfs = {hit.orf for hit in hits}
    hits_fragmented = set()
    for contig, contig_hits in sort_hits_by_contig(hits_remaining).items():
        if contig not in contig_ranges:
            continue
        for hit in contig_hits:
            # Skip hit if ORF is already assigned a hit or is not in any locus range
            if hit.orf in hit_orfs:
                continue
            if not any(hit.orf.start in r or hit.orf.end in r for r in contig_ranges[contig]):
                continue
            # Apply some sanity filtering here - not exposed to user
            if hit.bitscore < 200:
                continue
            hits_fragmented.add(hit)

    # Group by ORF and select best hit
    hits_selected = {'one': set(), 'two': set(), 'three': set()}
    for orf, orf_hits in sort_hits_by_orf(hits_fragmented).items():
        [region] = {database.get_region(hit.sseqid) for hit in orf_hits}
        if region in {'one', 'three'}:
            [hit_best] = region_common.select_best_hits(orf_hits)
        else:
            serotype = region_specific.determine_serotype(orf, orf_hits, region_specific.NEIGHBOUR_DIST, hits)
            hit_best = region_specific.perform_selection(orf_hits, serotype)
        hits_selected[region].add(hit_best)

    # Update region_group and hit.broken status
    for region, hits_fragmented in hits_selected.items():
        contigs_new = {hit.orf.contig for hit in hits_fragmented}
        region_groups[region].orf_hits.update(hits_fragmented)
        region_groups[region].contigs.update(contigs_new)
        hits_remaining -= hits_fragmented
        for hit in hits_fragmented:
            hit.broken = True


def get_proximal_ranges(hits, contig_fastas):
    contig_ranges = dict()
    for contig, hits in sort_hits_by_contig(hits).items():
        first_hit, *hits_sorted = sorted(hits, key=lambda h: h.orf.start)
        start = last_position = first_hit.orf.end
        contig_ranges[contig] = list()
        for hit in hits_sorted:
            if hit.orf.start - last_position > 5000:
                # Record
                end_bound = last_position + 5000
                start_bound = start - 5000
                contig_ranges[contig].append(range(start_bound, end_bound))
                # Restart
                start = hit.orf.start
            # Update
            last_position = hit.orf.end
        # Catch dangling
        end_bound = last_position + 5000
        start_bound = start - 5000
        contig_ranges[contig].append(range(start_bound, end_bound))

    # If are hits near a contig boundary, allow fragments to be found near any contig boundary
    # TODO: is there a better way to do this without a switch?
    allow_near_boundary = False
    for contig, ranges in contig_ranges.items():
        contig_start_pad = 2000
        contig_end_pad = len(contig_fastas[contig]) - 2000
        if any(min(r) <= contig_start_pad or max(r) >= contig_end_pad for r in ranges):
            allow_near_boundary = True
            break
    if allow_near_boundary:
        # Blindly add additional ranges, may overlap
        for contig in contig_fastas:
            if contig not in contig_ranges:
                contig_ranges[contig] = list()
            right_start = len(contig_fastas[contig]) - 2000
            contig_ranges[contig].append(range(0, 2000))
            contig_ranges[contig].append(range(right_start, len(contig_fastas[contig])))

    return contig_ranges


def collect_nearby_orfs(region_groups, orfs_all):
    hits_selected = {hit for group in region_groups.values() for hit in group.orf_hits}
    orfs_selected = sort_hits_by_orf(hits_selected)
    orfs_remaining = set(orfs_all) - set(orfs_selected)

    # If there are no ORFs without hits, make an early exit
    if not orfs_remaining:
        return set()

    nearby_orfs = set()
    for contig, contig_hits in sort_hits_by_contig(hits_selected).items():
        orfs = run_nearby_orf_collection(contig, contig_hits, orfs_remaining)

        # Apply some sanity filtering here - not exposed to user
        orfs_filtered = set()
        for orf in orfs:
            if (orf.end - orf.start) <= 200:
                continue
            orfs_filtered.add(orf)
        nearby_orfs |= orfs_filtered
    return nearby_orfs


def run_nearby_orf_collection(contig, contig_hits, orfs_remaining):
    hits_start, hits_end = get_elements_bounds(contig_hits)
    range_start = hits_start - NEARBY_FLANK_DIST
    range_end = hits_end + NEARBY_FLANK_DIST

    # The cap locus can be split and found at either end of a large contig, check for this
    if range_end - range_start > 60000:
        last_position = 0
        contig_hits_sorted = sorted(contig_hits, key=lambda h: h.orf.start)
        for i, hit in enumerate(contig_hits_sorted):
            if hit.orf.start - last_position > 5000:
                break
            last_position = hit.orf.end
        # Recursing on split elements
        orfs_lower = run_nearby_orf_collection(contig, contig_hits_sorted[:i], orfs_remaining)
        orfs_upper = run_nearby_orf_collection(contig, contig_hits_sorted[i:], orfs_remaining)
        return orfs_lower | orfs_upper
    else:
        return collect_elements_in_bounds(range_start, range_end, contig, orfs_remaining)


def get_elements_bounds(elements):
    elements_sorted = sorted(elements, key=lambda k: k.orf.start)
    start = min(elements_sorted[0].orf.start, elements_sorted[0].orf.end)
    end = max(elements_sorted[-1].orf.start, elements_sorted[-1].orf.end)
    return start, end


def collect_elements_in_bounds(start, end, contig, elements):
    # Apply some lambda madness
    test_element = list(elements)[0]
    if hasattr(test_element, 'start'):
        getter = lambda e: e
    elif hasattr(test_element, 'orf'):
        getter = lambda e: e.orf

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


def blast_missing_genes(region_groups, contig_fastas, database_fps):
    # Find missing genes - currently only operating as missing or not
    # TODO: uses counts - then select best n hits
    missing_genes = set()
    for region, region_data in region_groups.items():
        genes = {hit.sseqid for hit in region_data.orf_hits}
        if region in {'one', 'three'}:
            missing_genes |= genes ^ database.SCHEME[region]
        elif region == 'two':
            for serotype in region_data.serotypes:
                missing_genes |= genes ^ database.SEROTYPES[serotype]

    if not missing_genes:
        return

    # Get nucleotide sequence including and around locus
    locus_sequences, query_contigs, query_offsets = collect_proximal_locus_sequence(region_groups, contig_fastas)

    # Write out query sequences and run alignment
    with tempfile.TemporaryDirectory() as dh:
        query_fp = pathlib.Path(dh, 'locus_seq.fasta')
        with query_fp.open('w') as fh:
            for i, sequence in enumerate(locus_sequences):
                print('>%s' % i, sequence, sep='\n', file=fh)
        hits_missing = database.run_search(query_fp, database_fps)

    # Filter for missing genes and apply some sanity filtering here - not exposed to user
    hits_filtered = {hit for hit in hits_missing if hit.sseqid in missing_genes}
    hits_filtered = {hit for hit in hits_filtered if hit.bitscore >= 200}

    # Create and add annotation.SeqSection to each blast hit
    for hit in hits_filtered:
        hit.seq_section = create_seq_section(hit, query_offsets, query_contigs)

    # Ensure there are no overlaps - select if there are
    # TODO: allow small amount of overlap with ORFs
    hits_selected = set()
    hits = get_all_orf_hits(region_groups)
    orfs = {hit.orf for hit in hits}
    orf_ranges = [range(orf.start, orf.end) for orf in orfs]
    for hit in hits_filtered:
        if any(hit.seq_section.start in r or hit.seq_section.end in r for r in orf_ranges):
            continue
        hits_selected.add(hit)

    # Add hits to region group
    for hit in hits_selected:
        region = database.get_region(hit.sseqid)
        region_groups[region].blast_hits.add(hit)


def collect_proximal_locus_sequence(region_groups, contig_fastas):
    hits = get_all_orf_hits(region_groups)
    locus_ranges = get_proximal_ranges(hits, contig_fastas)
    locus_sequences = list()
    query_offsets = list()
    query_contigs = list()
    for contig, locus_ranges in locus_ranges.items():
        for locus_range in locus_ranges:
            start = locus_range.start if locus_range.start > 0 else 0
            sequence = contig_fastas[contig][start:locus_range.stop]
            locus_sequences.append(sequence)
            query_offsets.append(start)
            query_contigs.append(contig)
    return locus_sequences, query_contigs, query_offsets


def create_seq_section(hit, query_offsets, query_contigs):
    query_index = int(hit.qseqid)
    contig = query_contigs[query_index]
    start = hit.qstart + query_offsets[query_index]
    end = hit.qend + query_offsets[query_index]
    # Determine strand
    strand = -1 if hit.sstart > hit.send else 1
    # Ensure start is always lower than end
    if start > end:
        start, end = end, start
    return annotation.SeqSection(contig, start, end, strand)


def discover_is1016(region_groups, contig_fastas, database_fp):
    # Get nucleotide sequence including and around locus
    locus_sequences, query_contigs, query_offsets = collect_proximal_locus_sequence(region_groups, contig_fastas)

    # Write out query sequences and run search
    with tempfile.TemporaryDirectory() as dh:
        query_fp = pathlib.Path(dh, 'locus_seq.fasta')
        with query_fp.open('w') as fh:
            for i, sequence in enumerate(locus_sequences):
                print('>%s' % i, sequence, sep='\n', file=fh)
        hits = database.run_search(query_fp, database_fp)

    # Apply some sanity filtering here - not exposed to user
    hits_filtered = set()
    for hit in hits:
        if hit.bitscore < 200:
            continue
        if hit.length < 200:
            continue
        if hit.evalue > 0.5:
            continue
        if hit.pident < 60:
            continue
        hits_filtered.add(hit)

    # TODO: set threshold to define a hit has truncated/ broken

    # Create and add annotation.SeqSection to each blast hit
    for hit in hits_filtered:
        hit.seq_section = create_seq_section(hit, query_offsets, query_contigs)
    return hits_filtered


def get_all_hits(locus_data):
    hits_all = set()
    for region, region_data in locus_data.regions.items():
        hits_all |= region_data.orf_hits | region_data.blast_hits
    return hits_all | locus_data.is_hits


def get_all_orf_hits(region_groups):
    return {hit for region_data in region_groups.values() for hit in region_data.orf_hits}


def get_all_blast_hits(locus_data):
    return {hit for region in locus_data.regions.values() for hit in region.blast_hits}


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
        # Get contig
        if getattr(hit, 'orf'):
            contig = hit.orf.contig
        elif getattr(hit, 'seq_section'):
            contig = hit.seq_section.contig
        # Add
        try:
            contigs_hits[contig].add(hit)
        except KeyError:
            contigs_hits[contig] = {hit}
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


def get_hit_start(hit):
    '''Get the starting position in the query sequence of a hit'''
    if getattr(hit, 'orf'):
        return hit.orf.start
    elif getattr(hit, 'seq_section'):
        return hit.seq_section.start
