# Haemophilus Capsule Typing
[![Build Status](https://travis-ci.org/scwatts/hi_capsule.svg?branch=master)](https://travis-ci.org/scwatts/hi_capsule)
[![Code Coverage](https://codecov.io/gh/scwatts/hi_capsule/branch/master/graph/badge.svg)](https://codecov.io/gh/scwatts/hi_capsule)
[![License](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)

Identify *cap* loci serotype and structure in your *Haemophilus influenzae* assemblies.

**This tools remains under developemnt and is not quite ready to use**


## Table of contents
* [Table of contents](#table-of-contents)
* [Introduction](#introduction)
* [Requirements](#requirements)
* [Installation](#installation)
* [Usage](#usage)
* [License](#license)


## Introduction
The *cap* locus of *H. influenzae* are categorised into 6 different groups based on serology (a-f). The locus has three
functionally distinct regions, designated `region I`, `region II`, and `region III`. Genes within `region I` (`bexABCD`) and
`region III` (`hcsAB`) are associated with membrane transport and post-translation modification. The `region II` genes encode
serotype specific proteins, each having a distinct set of genes. The difficulty in typing the *H. influenzae* *cap* results
from the large structual diversity seen in the locus. Partial or whole loci subject to duplication or deletion is not
uncommon.

This tools automates identification of the *cap* locus, describes the structual layout of loci, and performs *in silico* serotyping.


## Requirements
There are a couple of software dependencies that are required:
* `Python`, version 3.6 or above
    * `biopython`
* `BLAST+`, version 2.2.28 or above
    * `makeblastdb`
    * `blastn`

## Installation
Recommended method of installation is via `pip`:
```bash
pip3 install git+https://github.com/scwatts/hi_capsule
```

## Usage
```bash
./hi_capsule.py --database_dir data/fasta/ --query_fp data/GCA_000210875.1_ASM21087v1_genomic.fasta
```


## License
[GNU General Public License v3.0](https://www.gnu.org/licenses/gpl-3.0.en.html)
