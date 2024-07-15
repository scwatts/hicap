# hicap

Identify _cap_ locus serotype and structure in your _Haemophilus influenzae_ assemblies.

![Example locus](https://image.ibb.co/kGLC1e/type_b_example.png)

## Table of contents

- [Table of contents](#table-of-contents)
- [Introduction](#introduction)
- [Citation](#citation)
- [Requirements](#requirements)
- [Installation](#installation)
- [Outputs](#outputs)
- [Usage](#usage)
- [Example](#example)
- [License](#license)

## Introduction

The _cap_ locus of _H. influenzae_ are categorised into 6 different groups based on serology (a-f). There are three
functionally distinct regions of the _cap_ locus, designated `region I`, `region II`, and `region III`. Genes within `region
I` (`bexABCD`) and `region III` (`hcsAB`) are associated with transport and post-translation modification. The `region II`
genes encode serotype specific proteins, with each serotype (a-f) having a distinct set of genes. _cap_ loci are often
subject to structural changes (e.g. duplication, deletion) making the process of _in silico_ typing and characterisation of
loci difficult.

`hicap` automates identification of the _cap_ locus, describes the structural layout, and performs _in silico_ serotyping.

## Citation

If you use this tool, please cite the `hicap` paper:

- Watts, S. C., & Holt, K. E. (2019). hicap: <em>in silico</em> serotyping of the <em>Haemophilus influenzae</em>
  capsule locus. Journal of Clinical Microbiology, JCM.00190-19. https://doi.org/10.1128/JCM.00190-19

## Requirements

There are a couple of software dependencies that are required.

- `Python`, version 3.6 or above with the following packages:
  - `Biopython`, version 1.79 or above
  - `ReportLab`, version 3.5.0 or above
- `Prodigal`, version 2.6.1 or above
- `BLAST+`, version 2.2.28 or above. Commands used are:
  - `makeblastdb`
  - `blastn`

## Installation

The recommended method of installation is bioconda:

```bash
conda install -c bioconda -c conda-forge hicap
```

Otherwise you can install `hicap` with pip:

```bash
# Install into user directory via pip
pip3 install --user git+https://github.com/scwatts/hicap.git

# Check install
hicap --help
```

Or clone the git repo and use the `hicap-runner.py` script:

```bash
# Install into current directory by cloning
git clone https://github.com/scwatts/hicap.git

# Check install
./hicap/hicap-runner.py --help
```

If installing `hicap` by the `pip` or `clone` method, you'll need to satisfy all software dependencies manually.

## Usage

Basic usage requires an input genome assembly and an output directory for result files:

```bash
# Create output directory and run serotype prediction
mkdir -p output/
hicap --query_fp input_genome.fasta --output_dir output/
```

There are several parameters that can be set manually. The default settings have been empirically determined to be optimal on
a test dataset so different settings may be more appropriate for other input data. To list this parameters use, the extended
help for `hicap`:

```text
hicap --help_all
usage: hicap -q QUERY_FP -o OUTPUT_DIR [-d DATABASE_DIR] [--gene_coverage GENE_COVERAGE] [--gene_identity GENE_IDENTITY]
             [--broken_gene_length BROKEN_GENE_LENGTH] [--broken_gene_identity BROKEN_GENE_IDENTITY] [--log_fp LOG_FP] [--debug] [-v] [-h]
             [--help_all]

File input and output:
  -q QUERY_FP, --query_fp QUERY_FP              Input FASTA query
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR        Output directory
  -d DATABASE_DIR, --database_dir DATABASE_DIR  Directory containing locus database. [default: /home/stephen/.local/lib/python3.6/site-
                                                packages/hicap/database]
  -s, --full_sequence                           Write the full input sequence out to the genbank file rather than just the region
                                                surrounding and including the locus

Search parameters:
  --gene_coverage GENE_COVERAGE                 Minimum percentage coverage to consider a single gene complete. [default: 0.80]
  --gene_identity GENE_IDENTITY                 Minimum percentage identity to consider a single gene complete. [default: 0.70]
  --broken_gene_length BROKEN_GENE_LENGTH       Minimum length to consider a broken gene. [default: 60]
  --broken_gene_identity BROKEN_GENE_IDENTITY   Minimum percentage identity to consider a broken gene. [default: 0.80]

Other:
  --log_fp LOG_FP                               Record logging messages to file
  --debug                                       Print debug messages
  -v, --version                                 Show version number and exit
  -h, --help                                    Show this help message and exit
  --help_all                                    Display extended help
```

## Outputs

Analysis using `hicap` will generate three result files in the specified output directory:

- `Summary`: a somewhat machine parsable file with detailed summary information
- `Genbank`: a genbank file with sequence marked up with _cap_ locus annotations
- `Graphic`: a visual representation of the annotated _cap_ locus

### Summary

The summary output file contains information relating to the annotation of the _cap_ locus. Below is a description of each
column:

| Column               | Description                                                                                          |
| -------------------- | ---------------------------------------------------------------------------------------------------- |
| `isolate`            | input isolate name                                                                                   |
| `predicted_serotype` | `hicap` serotype prediction                                                                          |
| `attributes`         | attributes of the locus (e.g. missing genes, truncated genes, duplication, etc)                      |
| `genes_identified`   | _cap_ genes identified. genes on different contigs delimited by`;`. truncation shown by trailing `*` |
| `locus_location`     | location of _cap_ genes. contigs delimited by `;` and matches `gene_identified` positions            |
| `region_n_genes`     | expected and identified genes in region `n`. missing names are shown when applicable                 |
| `IS1016_hits`        | count of IS*1016* hits found                                                                         |

### Genbank

The genbank output file provides marked up input sequence data with annotations of the _cap_ locus. This file will contain a
single record for each contig which the capsular locus is identified on. The default setting is to only include the sequence
of the _cap_ locus and surrounding region. This behaviour can be overridden to include all input sequence data by specific
`--full-sequence` on the command line.

_cap_ locus regions are annotated by the `misc_feature` feature with a `/note` qualifier specifying an identifier for that
specific region. _cap_ genes are given the `CDS` feature with a `/gene` qualifier set to the appropriate gene name. ORFs
found nearby that are not typically part of the _cap_ locus are also marked as a `CDS` feature but `/gene` is set to `orf`
with an incremental counter suffix for differentiation.

All `CDS` features have a `/note` qualifier which contains various details, each separated by `;`. Possible values are
`region_one`, `region_two`, `region_three`, `fragment`, `no_orf`, `insertion_sequence`, `misc_orf`.

### Graphic

The graphical output provides a visualisation of the annotated _cap_ locus. The output format is `svg` and can be viewed in
any modern browers or imagine viewer with `svg` support. Genes of each region are coloured differently; green, red, and
yellow for `region I`, `region II`, and `region III` respectively. Truncated ORFs are coloured using a darker, desaturated
corresponding region colour. ORFs which are not part of the _cap_ locus are left grey. Regions with homology to IS*1016* are
annotated using small blue arrows.

![single track](https://i.ibb.co/WWPFFCd/Hi76.png)

Each track of the visualisation represents a contig. The _cap_ locus of the example below is located on three different
contigs.

![multi track](https://i.ibb.co/fGzyxL8/Hi83.png)

## Example

To provide an example of `hicap` usage, we will annotate the _cap_ locus in the following isolates:

| Isolate     | BioSample Accession | Notes                                                                  |
| ----------- | ------------------- | ---------------------------------------------------------------------- |
| `Hi75`      | `SAMEA33515668`     | Downloaded SRA reads (`ERX1834398`; ion torrent) and assembled locally |
| `Hi84`      | `SAMEA33522418`     | Downloaded SRA reads (`ERX1834407`; ion torrent) and assembled locally |
| `M21328`    | `SAMN09704914`      | Downloaded GenBank assembly (`GCA_003497005.1`)                        |
| `M26329`    | `SAMN09704930`      | Downloaded GenBank assembly (`GCA_003492045.1`)                        |
| `PTHi-1539` | `SAMEA4643429`      | Downloaded GenBank assembly (`GCA_900407865.1`)                        |

For `Hi75` and `Hi84` reads were assembled using SPAdes 3.13.0. All assemblies have been placed into the `example/`
directory. `hicap` can be run in series or in parallel (using `gnu parallel`).

**Serial** execution:

```bash
mkdir -p output/
for assembly_fp in example/*fasta; do
  hicap --query_fp "${assembly_fp}" --output_dir output/;
done
```

**Parallel** execution:

```bash
mkdir -p output/
parallel hicap --query_fp {} --output_dir output/ ::: example/*fasta
```

The summary output for these assemblies can be combined using `sed`:

```bash
summary_fps=(output/*tsv);
{ head -n1 "${summary_fps[0]}"; sed '/^#/d' "${summary_fps[@]}"; } > example_summary.tsv
```

This file can be viewed on the command line or imported into excel.

```bash
column -t -s$'\t' example_summary.tsv | less -S
```

To view an isolate's genbank output in `artemis`, you first must combine the multiple genbank records. This can be done
using `seqret` from `EMBOSS`. Using `M21328` as an example here:

```bash
union -sequence output/M21328.gbk -outseq output/M21328_combined.gbk -osformat genbank -feature
```

Note that while combining multiple records enables viewing in `artemis`, contig boundaries are masked.

## License

[GNU General Public License v3.0](https://www.gnu.org/licenses/gpl-3.0.en.html)
