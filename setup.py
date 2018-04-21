#!/usr/bin/env python3
import setuptools


# Set package details
package_name = 'hi_capsule'
package_description = 'Hi capsule typing'
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
        scripts=['hi_capsule.py']
)
