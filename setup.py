#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import fnmatch
from setuptools import setup, find_packages


# Path to the directory that contains this setup.py file.
base_dir = os.path.abspath(os.path.dirname(__file__))

def find_files(directory):
    matches = []

    for root, dirnames, filenames in os.walk(directory):
        for filename in fnmatch.filter(filenames, '*'):
            matches.append(os.path.join(root, filename))

    return matches


setup(
    name="ringtail",
    author="Forli Lab",
    author_email="forli@scripps.edu",
    url="https://github.com/forlilab/Ringtail",
    description='Package for creating SQLite database from virtual screening DLGs and performing filtering on results.',
    long_description=open(os.path.join(base_dir, 'README.md')).read(),
    long_description_content_type="text/markdown",
    packages=find_packages(),
    scripts=["scripts/run_ringtail.py"],
    package_data={"ringtail" : ["data/*"]},
    data_files=[("", ["README.md", "LICENSE"]),
                ("scripts", find_files("scripts"))],
    include_package_data=True,
    zip_safe=False,
    install_requires=['numpy>=1.18'],
    python_requires='>=3.8.*',
    license="GPL-v3",
)