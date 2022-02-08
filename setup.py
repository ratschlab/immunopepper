#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages
from Cython.Build import cythonize

with open('README.md') as readme_file:
    readme = readme_file.read()

setup_requirements = ['pytest-runner']

with open('requirements.txt') as f:
    requirements = list(f.readlines())

test_requirements = ['pytest']

setup(
    author="Andre Kahles",
    author_email='andre.kahles@inf.ethz.ch',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
    ],
    description="Software to translate splicing graphs into peptides",
    entry_points = {
        'console_scripts': ['immunopepper=immunopepper.immunopepper:cmd_entry'],
    },
    install_requires=requirements,
    license="MIT license",
    long_description=readme,
    include_package_data=True,
    include_dirs='.',
    keywords='immunopepper',
    name='immunopepper',
    packages=find_packages(include=['immunopepper']),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/ratschlab/immunopepper',
    version='2.0.0',
    zip_safe=False,
    ext_modules = cythonize('./immunopepper/dna_to_peptide.pyx', include_path=['.'], annotate=True)
)
