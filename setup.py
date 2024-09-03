#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import fnmatch
from setuptools import setup, find_packages


# Path to the directory that contains this setup.py file.
base_dir = os.path.abspath(os.path.dirname(__file__))


def find_files(directory):
    matches = []

    for root, _, filenames in os.walk(directory):
        for filename in fnmatch.filter(filenames, "*"):
            matches.append(os.path.join(root, filename))

    return matches

setup(
    name="ringtail",
    version="2.0.0",
    author="Forli Lab",
    author_email="forli@scripps.edu",
    url="https://github.com/forlilab/Ringtail",
    description="Package for creating database from virtual screening files and performing filtering on results.",
    long_description=open(os.path.join(base_dir, "README.md")).read(),
    long_description_content_type="text/markdown",
    packages=find_packages(),
    include_package_data=True,
    zip_safe=False,
    license="L-GPL-v2.1",
    keywords=[
        "virtual screening",
        "molecular modeling",
        "drug discovery",
        "drug design",
        "docking",
        "autodock",
    ],
    classifiers=[
        "Environment :: Console",
        "Environment :: Other Environment",
        "Intended Audience :: Education",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU Lesser General Public License v2 or later (LGPLv2+)",
        "Natural Language :: English",
        "Operating System :: MacOS",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3.9",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Software Development :: Libraries",
    ],
    entry_points={
        'console_scripts': [
            'rt_process_vs=ringtail.cli.rt_process_vs:main',
            'rt_compare=ringtail.cli.rt_compare:main',
            'rt_db_v100_to_v110=ringtail.cli.rt_db_v100_to_v110:main',
            'rt_db_to_v200=ringtail.cli.rt_db_to_v200:main',
            'rt_generate_config_file=ringtail.cli.rt_generate_config_file:main'
        ]
    }
)
