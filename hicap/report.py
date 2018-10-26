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
    # Get contig sequences
    logging.info('Creating genbank records')
    fasta = utility.read_fasta(fasta_fp)
    hits_all = [hit for group in region_groups.values() for hit in group.hits]
    contig_sequences = collect_contig_sequences(fasta, hits_all, nearby_orfs)

    # Create base records
    position_deltas = dict()
    gb_records = dict()
    for contig, (position_delta, sequence) in contig_sequences.items():
        sequence_record=Bio.Seq.Seq(sequence, Bio.Alphabet.IUPAC.unambiguous_dna)
        gb_records[contig] = Bio.SeqRecord.SeqRecord(seq=sequence_record, name=contig, id=fasta_fp.stem)
        position_deltas[contig] = position_delta

    # Add hits
    for contig, contig_hits in locus.sort_hits_by_contig(hits_all).items():
        position_delta = position_deltas[contig]
        for hit in sorted(contig_hits, key=lambda k: k.orf.start):
            # Get appropriate representation of gene name
            region = hit.region if hit.region else locus.get_gene_region(hit.sseqid)
            qualifiers = {'gene': hit.sseqid, 'region': region}
            if hit.broken:
                qualifiers['note'] = 'fragment'
            # Create feature record
            feature = create_feature(hit.orf.start, hit.orf.end, position_delta, hit.orf.strand, qualifiers)
            gb_records[contig].features.append(feature)

    # Add ORFs
    orf_counter = 0
    for contig, orfs in locus.sort_orfs_by_contig(nearby_orfs).items():
        position_delta = position_deltas[contig]
        for orf in sorted(orfs, key=lambda o: o.start):
            orf_counter += 1
            qualifiers = {'gene': 'orf_%s' % orf_counter, 'region': 'none'}
            # Create feature record
            feature = create_feature(orf.start, orf.end, position_delta, orf.strand, qualifiers)
            gb_records[contig].features.append(feature)

    # Sort features by location
    for contig in gb_records.keys():
        gb_records[contig].features = sorted(gb_records[contig].features, key=lambda f: f.location.start)

    return [record for record in gb_records.values()]


def collect_contig_sequences(fasta, hits, nearby_orfs):
    # Sort all ORFs by contig
    hit_orfs = {hit.orf for hit in hits}
    contig_orfs = locus.sort_orfs_by_contig(hit_orfs | nearby_orfs)
    contig_sequences = dict()
    for contig, orfs in contig_orfs.items():
        # Get the most left and most right ORF associated with a hit
        orfs_sorted = sorted(orfs, key=lambda o: o.start)
        orf_start = None
        orf_end = None
        for orf in orfs_sorted:
            if orf in nearby_orfs:
                continue
            orf_end = orf
            if not orf_start:
                orf_start = orf

        # Apply sequencing padding - extend if we do not extend beyond nearby orfs
        start = orf_start.start - SEQ_PADDING
        if start > orfs_sorted[0].start:
            start = orfs_sorted[0].start
        start = max(start, 0)
        end = orf_end.end + SEQ_PADDING
        if end < orfs_sorted[-1].end:
            end = orfs_sorted[-1].end
        end = min(end, len(fasta[contig]))
        contig_sequences[contig] = (start, fasta[contig][start:end])
    return contig_sequences


def create_feature(start, end, delta, strand, qualifiers):
    feature_start = start - delta if (start - delta) > 1 else 1
    feature_end = end - delta
    feature_loc = Bio.SeqFeature.FeatureLocation(start=feature_start, end=feature_end, strand=strand)
    return Bio.SeqFeature.SeqFeature(location=feature_loc, type='CDS', qualifiers=qualifiers)
