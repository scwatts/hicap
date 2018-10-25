import logging
import pathlib
import re


import Bio.Alphabet
import Bio.Seq
import Bio.SeqRecord
import Bio.SeqFeature
import Bio.SeqIO
import Bio.Graphics.GenomeDiagram


from . import database
from . import locus
from . import utility
from . import graphic


SEQ_PADDING = 1000


class SummaryData:

    def __init__(self):
        self.completeness = dict.fromkeys(database.SCHEME, None)
        self.truncated_genes = dict.fromkeys(database.SCHEME, None)
        self.duplicated = None
        self.multiple_contigs = None
        self.serotypes = None
        self.hits_by_contig = None


def create_report(region_groups, nearby_orfs, fasta_fp, prefix, output_dir):
    # Report
    output_report_fp = pathlib.Path(output_dir, '%s.tsv' % prefix)
    summary_data = create_summary(region_groups)
    with output_report_fp.open('w') as fh:
        write_summary(summary_data, fh)

    # Genbank - create
    # Genbank spec requires contig names/ locus names of less than 20 characters but
    # we want full contigs names in the graphic. So we'll truncate the names after
    # generating the graphic
    output_gbk_fp = pathlib.Path(output_dir, '%s.gbk' % prefix)
    genbank_data = create_genbank_record(region_groups, nearby_orfs, fasta_fp)

    # Graphic
    # TODO: legend?
    output_svg_fp = pathlib.Path(output_dir, '%s.svg' % prefix)
    graphic_data = graphic.create_graphic(genbank_data, prefix)
    svg_data = graphic.patch_graphic(graphic_data)
    with output_svg_fp.open('w') as fh:
        fh.write(svg_data)

    # Genbank - write
    for record in genbank_data:
        record.name = record.name.split()[0][:15]
        record.id = record.id.split()[0][:15]
    with output_gbk_fp.open('w') as fh:
        Bio.SeqIO.write(genbank_data, fh, 'genbank')


def create_summary(region_groups):
    summary_data = SummaryData()
    for region in ('one', 'two', 'three'):
        group = region_groups[region]
        # Completeness and truncated genes
        genes_found = {hit.sseqid for hit in group.hits}
        if region in ('one', 'three'):
            genes_expected = database.SCHEME[region]
        else:
            genes_expected = set()
            for serotype in group.serotypes:
                genes_expected.update(database.SEROTYPES[serotype])
        genes_missing = genes_expected - genes_found
        summary_data.completeness[region] = (genes_missing, len(genes_found), len(genes_expected))
        summary_data.truncated_genes[region] = {hit for hit in group.hits if hit.broken}
        # Duplication. Verbose for clarity
        if group.hits and is_duplicated(group.hits):
            summary_data.duplicated = True

    # Serotype
    summary_data.serotypes = region_groups['two'].serotypes

    # On multiple contigs
    hits_all_gen = (hit for group in region_groups.values() for hit in group.hits)
    contig_hits = locus.sort_hits_by_contig(hits_all_gen)
    if len(contig_hits) > 1:
        summary_data.multiple_contigs = True

    # Store hits sorted by contig
    summary_data.hits_by_contig = contig_hits
    return summary_data


def write_summary(data, fh):
    attributes = list()
    print('#', ','.join(data.serotypes), sep='', file=fh)
    if any(region_complete[0] for region_complete in data.completeness.values()):
        attributes.append('missing_genes')
    else:
        attributes.append('full_gene_complement')
    if any(data.truncated_genes.values()):
        attributes.append('truncated_genes')
    if data.multiple_contigs:
        attributes.append('fragmented_locus')
    if data.duplicated:
        attributes.append('duplicated')
    print('#', end='', file=fh)
    print(*attributes, sep=',', file=fh)
    print('#', end='', file=fh)
    print('contig', 'start', 'end', 'genes', sep='\t', file=fh)
    for contig, contig_hits in data.hits_by_contig.items():
        hits_sorted = sorted(contig_hits, key=lambda k: k.orf.start)
        hits_bounds = [b for hit in hits_sorted for b in (hit.orf.start, hit.orf.end)]
        region_start = min(hits_bounds)
        region_end = max(hits_bounds)
        gene_names = get_gene_names(hits_sorted)
        print(contig, region_start, region_end, ','.join(gene_names), sep='\t', file=fh)


def get_gene_names(hits):
    gene_names = list()
    for hit in hits:
        if hit.broken:
            gene_names.append(hit.sseqid + '*')
        else:
            gene_names.append(hit.sseqid)
    return gene_names


def is_duplicated(hits):
    genes_hits = locus.sort_hits_by_gene(hits)
    gene_counts = dict.fromkeys(genes_hits)
    # Sometimes a single gene is broken into multiple ORFs but wrt to duplication should
    # be considered as a single unit. We do this by setting a bound ~ to expected gene len
    for gene, gene_hits in genes_hits.items():
        hit_first, *hits_sorted = sorted(gene_hits, key=lambda hit: min(hit.orf.start, hit.orf.end))
        upper_bound = min(hit_first.orf.start, hit_first.orf.end) + hit_first.slen * 1.5
        gene_counts[gene] = 1
        for hit in hits_sorted:
            if max(hit.orf.start, hit.orf.end) >= upper_bound:
                gene_counts[gene] += 1
    return any(gene_count > 1 for gene_count in gene_counts.values())


def create_genbank_record(region_groups, nearby_orfs, fasta_fp):
    # TODO: handle sequence boundaries correctly - max(1000 padding or up to the most distant ORF)
    logging.info('Creating genbank records')
    fasta = utility.read_fasta(fasta_fp)
    hits_all = [hit for group in region_groups.values() for hit in group.hits]

    # Create base records
    contigs = {hit.orf.contig for hit in hits_all}
    position_delta = 0
    gb_records = dict()
    for contig in contigs:
        gb_records[contig] = Bio.SeqRecord.SeqRecord(seq='atgc', name=contig, id=fasta_fp.stem)


    # TODO: TEMP
    pd = dict()

    # Add hits
    for contig, contig_hits in locus.sort_hits_by_contig(hits_all).items():
        for hit in sorted(contig_hits, key=lambda k: k.orf.start):
            # TODO TEMP
            position_delta, block_sequence = get_block_sequence(contig_hits, fasta[contig], SEQ_PADDING)
            sequence = Bio.Seq.Seq(block_sequence, Bio.Alphabet.IUPAC.unambiguous_dna)
            gb_records[contig].seq = sequence
            pd[contig] = position_delta

            # Get appropriate representation of gene name
            region = hit.region if hit.region else locus.get_gene_region(hit.sseqid)
            qualifiers = {'gene': hit.sseqid, 'region': region}
            if hit.broken:
                qualifiers['note'] = 'fragment'
            feature_start = hit.orf.start - position_delta if (hit.orf.start - position_delta) > 1 else 1
            feature_end = hit.orf.end - position_delta
            feature_loc = Bio.SeqFeature.FeatureLocation(start=feature_start, end=feature_end, strand=hit.orf.strand)
            feature = Bio.SeqFeature.SeqFeature(location=feature_loc, type='CDS', qualifiers=qualifiers)
            gb_records[contig].features.append(feature)

    # Add ORFs
    orf_counter = 0
    for contig, orfs in locus.sort_orfs_by_contig(nearby_orfs).items():
        # Get current bounds - FeatureLocation.start always smallest regardless of strand
        start = min(gb_records[contig].features, key=lambda k: k.location.start)
        end = max(gb_records[contig].features, key=lambda k: k.location.end)

        # TODO TEMP
        position_delta = pd[contig]

        for orf in sorted(orfs, key=lambda o: o.start):
            orf_counter += 1
            qualifiers = {'gene': 'orf_%s' % orf_counter, 'region': 'none'}
            feature_start = orf.start - position_delta if (orf.start - position_delta) > 1 else 1
            feature_end = orf.end - position_delta
            feature_loc = Bio.SeqFeature.FeatureLocation(start=feature_start, end=feature_end, strand=orf.strand)
            feature = Bio.SeqFeature.SeqFeature(location=feature_loc, type='CDS', qualifiers=qualifiers)
            gb_records[contig].features.append(feature)

    # Add appropriately sizes sequence to records
    for contig, gb_record in gb_records.items():
        #sequence = Bio.Seq.Seq(fasta[contig], Bio.Alphabet.IUPAC.unambiguous_dna)
        pass

    return [record for record in gb_records.values()]


def get_block_sequence(hits, block_contig, margin):
    '''Extract sequence for loci with margin on either side'''
    hits_sorted = sorted(hits, key=lambda k: k.orf.start)
    orf_start = min(hits_sorted[0].orf.start, hits_sorted[0].orf.end)
    orf_end = min(hits_sorted[-1].orf.start, hits_sorted[-1].orf.end)

    # Calculate start end position for sequence with margin
    block_start = orf_start - margin if orf_start >= margin else 0
    block_end = orf_end + margin if (orf_end + margin) <= len(block_contig) else len(block_contig)

    # Return delta in position and sequence with margin
    return block_start, block_contig[block_start:block_end]
