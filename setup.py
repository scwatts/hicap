#!/usr/bin/env python3
import setuptools


import hicap


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
