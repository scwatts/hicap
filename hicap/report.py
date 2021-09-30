import pathlib


import Bio.SeqIO


from . import database
from . import genbank
from . import graphic
from . import locus
from . import utility


class SummaryData:

    def __init__(self):
        self.completeness = dict.fromkeys(database.SCHEME, None)
        self.truncated_genes = dict.fromkeys(database.SCHEME, None)
        self.hits_without_orfs = dict.fromkeys(database.SCHEME, None)
        self.duplicated = None
        self.multiple_contigs = None
        self.serotypes = None
        self.hits_by_contig = None
        self.is_hits = None


def write_outputs(locus_data, args):
    # Report
    prefix = pathlib.Path(args.query_fp.stem).stem if args.query_fp.name.endswith('gz') else args.query_fp.stem
    output_report_fp = pathlib.Path(args.output_dir, '%s.tsv' % prefix)
    fasta = utility.read_fasta(args.query_fp)
    summary_data = create_summary(locus_data, fasta)
    with output_report_fp.open('w') as fh:
        write_summary(summary_data, prefix, fh)

    # Genbank - create
    # We allow the user to request the full input sequence to be present in the output. The most
    # simple approach is therefore to generate a genbank record set for the image (which will
    # always have a substring of the full FASTA sequence) and another genbank record set for the
    # genbank output file if required
    # Genbank spec requires contig names/ locus names of less than 20 characters but
    # we want full contigs names in the graphic. So we'll truncate the names after
    # generating the graphic
    contig_sequences = genbank.collect_contig_sequences(fasta, locus_data)
    genbank_data = genbank.create_genbank_record(locus_data, contig_sequences)

    # Graphic
    # TODO: legend?
    # Loci on a single contig can be split at either end - rotate for graphic if required
    genbank_data = graphic.prepare_genbank(genbank_data)
    output_svg_fp = pathlib.Path(args.output_dir, '%s.svg' % prefix)
    graphic_data = graphic.create_graphic(genbank_data, prefix)
    svg_data = graphic.patch_graphic(graphic_data)
    with output_svg_fp.open('w') as fh:
        fh.write(svg_data)

    # Genbank - finalise and write
    output_gbk_fp = pathlib.Path(args.output_dir, '%s.gbk' % prefix)
    # Use full input sequence if requested
    if args.full_sequence:
        contig_sequences = {contig: (0, fasta[contig]) for contig in fasta}
        genbank_data = genbank.create_genbank_record(locus_data, contig_sequences)
    # Truncate names for legal genbank format
    for record in genbank_data:
        record.name = record.name.split()[0][:15]
        record.id = record.id.split()[0][:15]
    # Add misc_feature for cap locus (or locus fragments) - done in place
    genbank.add_locus_feature(genbank_data)
    for record in genbank_data:
        record.features = sorted(record.features, key=lambda f: f.location.start)
    with output_gbk_fp.open('w') as fh:
        Bio.SeqIO.write(genbank_data, fh, 'genbank')


def create_summary(locus_data, contig_sequences):
    summary_data = SummaryData()
    for region in ('one', 'two', 'three'):
        group = locus_data.regions[region]
        # Completeness and truncated genes
        genes_found = {hit.sseqid for hit in group.orf_hits}
        if region in ('one', 'three'):
            genes_expected = database.SCHEME[region]
        else:
            genes_expected = set()
            for serotype in group.serotypes:
                genes_expected.update(database.SEROTYPES[serotype])
        genes_missing = genes_expected - genes_found
        summary_data.completeness[region] = (genes_missing, len(genes_found), len(genes_expected))
        summary_data.truncated_genes[region] = {hit for hit in group.orf_hits if hit.broken}
        summary_data.hits_without_orfs[region] = {hit for hit in group.blast_hits}
        # Duplication. Verbose for clarity
        if group.orf_hits and is_duplicated(group.orf_hits, contig_sequences):
            summary_data.duplicated = True

    # Serotype and count of IS1016 hits
    summary_data.serotypes = locus_data.regions['two'].serotypes
    summary_data.is_hits = len(locus_data.is_hits)

    # Sort all hits by contig - use manaul sort here rather than generic function
    contig_hits = dict()
    for region_data in locus_data.regions.values():
        for hit in region_data.orf_hits | region_data.blast_hits:
            # Get contig - check which attribute is not None
            if getattr(hit, 'orf'):
                contig = hit.orf.contig
            elif getattr(hit, 'seq_section'):
                contig = hit.seq_section.contig
            # Add
            try:
                contig_hits[contig].add(hit)
            except:
                contig_hits[contig] = {hit}

    # Locus on multiple contigs
    if len(contig_hits) > 1:
        summary_data.multiple_contigs = True

    # Store hits sorted by contig
    summary_data.hits_by_contig = contig_hits
    return summary_data


def write_summary(data, prefix, fh):
    # Header
    header = ('isolate', 'predicted_serotype', 'attributes', 'genes_identified', 'locus_location',
              'region_I_genes', 'region_II_genes', 'region_III_genes', 'IS1016_hits')
    print('#', end='', file=fh)
    print(*header, sep='\t', file=fh)

    # Isolate and predict serotypes
    print(prefix, end='\t', file=fh)
    print(','.join(data.serotypes), end='\t', file=fh)

    # Locus attributes
    attributes = list()
    if any(region_complete[0] for region_complete in data.completeness.values()):
        attributes.append('missing_genes')
    else:
        attributes.append('full_gene_complement')
    if any(data.truncated_genes.values()):
        attributes.append('truncated_genes')
    if any(data.hits_without_orfs.values()):
        attributes.append('matches_without_orfs')
    if data.multiple_contigs:
        attributes.append('fragmented_locus')
    if data.duplicated:
        attributes.append('duplicated')
    # missing orfs
    print(','.join(attributes), end='\t', file=fh)

    # Genes found and contigs with location
    contig_genes = dict()
    contig_bounds = dict()
    for contig, contig_hits in data.hits_by_contig.items():
        hits_sorted = sorted(contig_hits, key=lambda h: locus.get_hit_start(h))
        contig_genes[contig] = ','.join(get_gene_names(hits_sorted))
        start = locus.get_hit_start(hits_sorted[0])
        end = locus.get_hit_end(hits_sorted[-1])
        contig_bounds[contig] = '%s:%s-%s' % (contig, start, end)
    print(*contig_genes.values(), sep=';', end='\t', file=fh)
    print(*contig_bounds.values(), sep=';', end='\t', file=fh)

    # Genes missing
    missing_genes = dict()
    for region, (genes, found, expected) in data.completeness.items():
        text = '%s/%s' % (found, expected)
        if genes:
            text = '%s (missing: %s)' % (text, ','.join(genes))
        missing_genes[region] = text
    print(*missing_genes.values(), sep='\t', end='\t', file=fh)

    # IS1016 hits
    print(data.is_hits, file=fh)


def get_gene_names(hits):
    gene_names = list()
    for hit in hits:
        if hit.broken:
            # Truncated gene
            gene_names.append(hit.sseqid + '*')
        elif getattr(hit, 'seq_section'):
            # Blast hit only - no ORF assocaited
            gene_names.append(hit.sseqid + '^')
        else:
            gene_names.append(hit.sseqid)

    return gene_names


def is_duplicated(hits, contig_sequences):
    genes_hits = locus.sort_hits_by_gene(hits)
    gene_counts = dict.fromkeys(genes_hits)
    # Sometimes a single gene is broken into multiple ORFs but wrt to duplication should
    # be considered as a single unit. We do this by setting a bound ~ to expected gene len
    for gene, gene_hits in genes_hits.items():
        hit_first, *hits_sorted = sorted(gene_hits, key=lambda h: min(h.orf.start, h.orf.end))
        upper_bound = min(hit_first.orf.start, hit_first.orf.end) + hit_first.slen * 1.5
        gene_counts[gene] = 1

        # If gene is close to a contig boundary, check nearby other boundaries for fragments
        contig_padding = hit_first.slen * 1.5
        near_bounds = near_contig_bounds(hit_first, contig_sequences, contig_padding)

        # Check for duplication
        for hit in hits_sorted:
            if max(hit.orf.start, hit.orf.end) < upper_bound:
                continue
            if near_bounds and near_contig_bounds(hit, contig_sequences, contig_padding):
                continue
            gene_counts[gene] += 1
    return any(gene_count > 1 for gene_count in gene_counts.values())


def near_contig_bounds(hit, contig_sequences, padding):
    # TODO: similar code is used elsewhere, refactor and simplify
    contig = hit.orf.contig
    contig_start_pad = padding
    contig_end_pad = len(contig_sequences[contig]) - padding
    return hit.orf.start <= contig_start_pad or hit.orf.end >= contig_end_pad
