![ringtail logo final](https://user-images.githubusercontent.com/41704502/170797800-53a9d94a-932e-4936-9bea-e2d292b0c62b.png)

*(Original artwork by Althea Hansel-Harris)*

# Ringtail
*A tool for handling results from virtual screening of molecules*

[![License: L-GPL v2.1](https://img.shields.io/badge/License-LGPLv2.1-blue.svg)](https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html)
![Build Status](https://github.com/forlilab/Ringtail/actions/workflows/python-package.yml/badge.svg?event=push)
[![Documentation Status](https://readthedocs.org/projects/ringtail/badge/?version=latest)](https://ringtail.readthedocs.io)

Ringtail is an open-source, lightweight, and highly customizable Python package for organizing, filtering, and exploring the results of molecular docking and virtual screening, from a handful of ligands to **tens of millions**. It reads collections of docking results such as SDFs from [AutoDock-6](https://github.com/forlilab/AutoDock), Docking Log Files (DLGs) from [AutoDock-GPU](https://github.com/ccsb-scripps/AutoDock-GPU), and PDBQTs from [AutoDock-Vina](https://github.com/ccsb-scripps/AutoDock-Vina) — into a compact database that stays fast to query as it grows, backed by either DuckDB (default) or SQLite. Results file parsing is parallelized across your CPUs for fast database writing.

Once your docking results are in a database, Ringtail provides a wealth of ways to apply your chemical intuition and narrow down the results to likely pharmacological hits: filter by docking score, ligand efficiency, receptor interactions, or ligand chemistry; cluster for diversity; compare hits across targets; and export exactly the molecules and data you want.

Use Ringtail however suits you, as a user-friendly [command line tool](https://ringtail.readthedocs.io/en/latest/cmdline.html) for straightforward screening, or an extensive [Python API](https://ringtail.readthedocs.io/en/latest/api.html) for building it into your own pipelines and scripts. In-depth documentation for both can be found on [ReadTheDocs](https://ringtail.readthedocs.io).

## Fast, even at scale

With the DuckDB backend, filtering stays in the seconds range as a library grows into the millions of ligands (Intel i9, 18 cores, 64 GB RAM, SSD; docking score alone, then combined with one, and two interaction filters):

| Ligands | Poses | Database size | Score | Score + 1 interaction | Score + 2 interactions |
|--:|--:|--:|--:|--:|--:|
| 100,000 | 277,048 | 0.20 GB | 1.2 s | 1.4 s | 1.4 s |
| 2,000,000 | 5,448,313 | 3.3 GB | 1.4 s | 5.4 s | 4.0 s |
| 9,039,451 | 24,801,508 | 15 GB | 2.0 s | 16.5 s | 13.8 s |


Database size scales roughly linearly with the number of stored poses, see [the changelog](https://ringtail.readthedocs.io/en/latest/changes.html) for full benchmarks.

## Installation

Ringtail requires Python ≥3.9 and is tested on Linux, macOS, and Windows. It's recommended to install Ringtail in a dedicated environment such as conda or micromamba.

```bash
$ conda create -n ringtail python=3.11
$ conda activate ringtail
```

**From conda-forge** (handles all dependencies):

```bash
$ conda install -c conda-forge ringtail
```

**From PyPI:**

```bash
$ pip install ringtail
```

When installing from PyPI you may need one or more dependencies, including `rdkit`, `scipy`, `pandas`, `prody`, and `meeko` (a Forli lab tool). DuckDB is the default storage backend for Ringtail, install if not already present (SQLite ships with Python):

```bash
$ pip install <dependency>
```

## Quick start

After installing, a virtual screen can be as simple as two commands:

```bash
# write a folder of docking results into a database
$ rt_process_vs write --docking_results results_folder/ --recursive

# filter for docking score and a specific interaction, and write the list of results to text log
$ rt_process_vs read --input_db output.db --eworst -6 --vdw_interactions A:VAL:279: --output_log hits.txt
```

Same example but using the API:

```python
from ringtail import RingtailCore

rtc = RingtailCore()
rtc.add_results_from_files(
    file_path="results_folder/",
    recursive=True,
)

rtc.filter(
    eworst=-6,
    vdw_interactions=[("A:VAL:279:", True)],
    output_log="hits.txt",
)
```


### Upgrading older databases

A database written with an earlier Ringtail will need to be upgraded to work with the current version. Use `rt_upgrade_db` with the target version; see [upgrading a database](https://ringtail.readthedocs.io/en/latest/upgrade_database.html) for details.

## Citing Ringtail

Ringtail is developed by the [Forli lab](https://forlilab.org/) at the [Center for Computational Structural Biology (CCSB)](https://ccsb.scripps.edu) at [Scripps Research](https://www.scripps.edu/).

This publication in JCIM describes the original design, implementation, and features of Ringtail:

[*Ringtail: A Python Tool for Efficient Management and Storage of Virtual Screening Results.*
Althea T. Hansel-Harris, Diogo Santos-Martins, Niccolò Bruciaferri, Andreas F. Tillack, Matthew Holcomb, and Stefano Forli.
*Journal of Chemical Information and Modeling* **2023** 63 (7), 1858-1864.
DOI: 10.1021/acs.jcim.3c00166](https://pubs.acs.org/doi/full/10.1021/acs.jcim.3c00166)

If using Ringtail in your work, please cite this publication.

## Contributing

Found a bug or have a feature request? Open an issue on [GitHub](https://github.com/forlilab/Ringtail/issues).

## License

Ringtail is released under the GNU LGPL-2.1-or-later license.
