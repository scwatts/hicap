# hicap
[![Build Status](https://travis-ci.org/scwatts/hicap.svg?branch=master)](https://travis-ci.org/scwatts/hicap)
[![Code Coverage](https://codecov.io/gh/scwatts/hicap/branch/master/graph/badge.svg)](https://codecov.io/gh/scwatts/hicap)
[![License](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)

Identify *cap* loci serotype and structure in your *Haemophilus influenzae* assemblies

**This tools remains under development and is not quite ready to use**


## Table of contents
* [Table of contents](#table-of-contents)
* [Introduction](#introduction)
* [Requirements](#requirements)
* [Installation](#installation)
* [Usage](#usage)
* [License](#license)


## Introduction
The *cap* loci of *H. influenzae* are categorised into 6 different groups based on serology (a-f). There are three
functionally distinct regions of the *cap* locus, designated `region I`, `region II`, and `region III`. Genes within `region
I` (`bexABCD`) and `region III` (`hcsAB`) are associated with membrane transport and post-translation modification. The
`region II` genes encode serotype specific proteins, with each serotype (a-f) having a distinct set of genes. *cap* loci are
often subject to structural changes (e.g. duplication, deletion) making the process of *in silico* typing and characterisation
of loci difficult.

`hicap` automates identification of the *cap* locus, describes the structural layout, and performs *in silico* serotyping.


## Requirements
There are a couple of software dependencies that are required:
* `Python`, version 3.6 or above
* `Biopython`, version 1.63 or above
* `BLAST+`, version 2.2.28 or above. Commands used are:
    * `makeblastdb`
    * `blastn`


## Installation
Recommended method of installation is via `pip`:
```bash
pip3 install git+https://github.com/scwatts/hicap
```


## Usage
```bash
./hicap.py --database_dir data/fasta/ --query_fp data/GCA_000210875.1_ASM21087v1_genomic.fasta
```


## License
[GNU General Public License v3.0](https://www.gnu.org/licenses/gpl-3.0.en.html)
