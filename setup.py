#!/usr/bin/env python3
import setuptools


# Set package details
package_name = 'hicap'
package_description = 'in silico typing of H. influenzae cap locus'
package_version = '0.0.1'
author = 'Stephen Watts'
license = 'GPLv3'


# Call setup
setuptools.setup(
        name=package_name,
        version=package_version,
        license=license,
        test_suite='tests',
        packages=setuptools.find_packages(),
        package_data={'hicap': ['database/*fasta']},
        entry_points={
                'console_scripts': ['hicap=hicap.main:main'],
            }
)
