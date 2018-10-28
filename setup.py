#!/usr/bin/env python3
from __future__ import print_function # allow us to detect python2 w/o syntax error
import setuptools
import sys


import hicap


if sys.version_info < (3,5):
    msg = 'error: hicap requires Python 3.5+. Python %d.%d detected'
    print(msg % sys.version_info[:2])
    sys.exit(1)

setuptools.setup(
        name='hicap',
        version=hicap.__version__,
        description='in silico typing of H. influenzae cap loci',
        author='Stephen Watts',
        license='GPLv3',
        test_suite='tests',
        packages=setuptools.find_packages(),
        package_data={'hicap': ['database/*fasta']},
        entry_points={
                'console_scripts': ['hicap=hicap.main:main'],
            }
)
