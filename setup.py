#!/usr/bin/env python
# -*- coding: utf-8 -*-


try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


readme = open('README.rst').read()

requirements = [
    'pysam',
    'python-levenshtein',
    'sexpdata',
]

test_requirements = []

setup(
    name='atactk',
    version='0.1.2',
    description="A toolkit for working with ATAC-seq data.",
    long_description=readme + '\n\n',
    author="The Parker Lab",
    author_email='parkerlab-software@umich.edu',
    url='https://github.com/ParkerLab/atactk',
    packages=['atactk', ],
    scripts=[
        'scripts/make_cut_matrix',
        'scripts/trim_adapters'
    ],
    include_package_data=True,
    install_requires=requirements,
    license="GPLv3+",
    zip_safe=False,
    keywords='atactk',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Natural Language :: English',
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
    ],
    test_suite='tests',
    tests_require=test_requirements
)
