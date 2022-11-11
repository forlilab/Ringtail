#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import fnmatch
import platform
from setuptools import setup, find_packages


# Path to the directory that contains this setup.py file.
base_dir = os.path.abspath(os.path.dirname(__file__))


def find_files(directory):
    matches = []

    for root, dirnames, filenames in os.walk(directory):
        for filename in fnmatch.filter(filenames, "*"):
            matches.append(os.path.join(root, filename))

    return matches

required_modules = ["rdkit>=2022.03.2", "scipy>=1.8.0",  "meeko>=0.4", "matplotlib", "pandas"]
if platform.system() == "Darwin":  # mac
    required_modules.append("multiprocess>=0.70.13")

setup(
    name="ringtail",
    author="Forli Lab",
    author_email="forli@scripps.edu",
    url="https://github.com/forlilab/Ringtail",
    description="Package for creating database from virtual screening files and performing filtering on results.",
    long_description=open(os.path.join(base_dir, "README.md")).read(),
    long_description_content_type="text/markdown",
    packages=find_packages(),
    scripts=["scripts/rt_process_vs.py", "scripts/rt_compare.py"],
    data_files=[("", ["README.md", "LICENSE"]), ("scripts", find_files("scripts"))],
    include_package_data=True,
    zip_safe=False,
    install_requires=required_modules,
    python_requires=">=3.9.*",
    license="L-GPL-v2.1",
)
