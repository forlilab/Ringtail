![ringtail logo final](https://user-images.githubusercontent.com/41704502/170797800-53a9d94a-932e-4936-9bea-e2d292b0c62b.png)

*(Original artwork by Althea Hansel-Harris)*

# Ringtail
*A tool for handling results from virtual screening of molecules*

[![License: L-GPL v2.1](https://img.shields.io/badge/License-LGPLv2.1-blue.svg)](https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html)
![Build Status](https://github.com/forlilab/Ringtail/actions/workflows/python-package.yml/badge.svg?event=push)
[![Documentation Status](https://readthedocs.org/projects/ringtail/badge/?version=latest)](https://ringtail.readthedocs.io)

Ringtail is an open-source Python package for organizing, filtering, and exploring molecular docking results at the scale of millions of ligands. It reads [AutoDock-6](https://github.com/forlilab/AutoDock) (AD6) SDFs, [AutoDock-GPU](https://github.com/ccsb-scripps/AutoDock-GPU) DLGs, and [AutoDock-Vina](https://github.com/ccsb-scripps/AutoDock-Vina) PDBQTs into compact DuckDB (default) or SQLite databases, with parallell input file processing. 

Filter by docking score, ligand efficiency, receptor interactions, or ligand chemistry, cluster for diversity, compare hits across different targets, and export selected results as SDF or CSV files. 

Use the easy-to-use [command line tool](https://ringtail.readthedocs.io/en/latest/cmdline.html) or the extensive [Python API](https://ringtail.readthedocs.io/en/latest/api.html) for building it into your own pipelines and scripts. See the [documentation](https://ringtail.readthedocs.io) for detailed instructions.

Some of the strengths of the Python API include accepting docked RDKit molecules from AD6 and Vina result strings directly, allowing docking pipelines to write results without intermediate files. Through the API you can also screen the poses by marking them as either accepted, maybe, or rejected, and you can attach comments to track screening decisions.


## Fast at scale

With the DuckDB backend, filtering stays in the seconds range as a library grows into the millions of ligands (Intel i9, 18 cores, 64 GB RAM, SSD; docking score alone, then combined with one, and two interaction filters):

| Ligands | Poses | Database size | Score | Score + 1 interaction | Score + 2 interactions |
|--:|--:|--:|--:|--:|--:|
| 100,000 | 277,048 | 0.20 GB | 1.2 s | 1.4 s | 1.4 s |
| 2,000,000 | 5,448,313 | 3.3 GB | 1.4 s | 5.4 s | 4.0 s |
| 9,039,451 | 24,801,508 | 15 GB | 2.0 s | 16.5 s | 13.8 s |


See the [changelog](https://ringtail.readthedocs.io/en/latest/changes.html) for full benchmarks and the [compression guide](https://ringtail.readthedocs.io/en/latest/compress.html) for how to easily transfer large databases.


## Installation

Ringtail requires Python ≥3.10 and is tested on Linux, macOS, and Windows. It's recommended to install Ringtail in a dedicated environment such as conda or micromamba.

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

When installing with pip, you must install all dependencies separately, including `duckdb`, `rdkit`, `numpy`, `scipy`, `pandas`, `packaging`, and `meeko >= 0.8.0` (a Forli lab tool), and optionally `prody`. DuckDB is the default storage backend for Ringtail >=3, install if not already present (or switch backend to SQLite, which ships with Python):

```bash
$ pip install <dependency>
```

## Quick start

After installing, import and filter docking results with two commands:

Ringtail defaults to AutoDock-6 SDF input. For AutoDock-GPU DLGs, specify `--docking_mode adgpu`; for Vina PDBQTs, use `--docking_mode vina`. In the Python API, pass `docking_mode="adgpu"` or `"vina"` to `add_results_from_files()`.

Provide the receptor used for docking as a Meeko Polymer JSON or PDBQT file. Ringtail also stores up to three poses per ligand as a default, if you wish to retain all poses use `--store_all_poses`. 

```bash
# write a folder of docking results into a database
$ rt_process_vs write --docking_results results_folder/ --recursive --receptor_file receptor.json --save_receptor

# filter for docking score and a specific interaction, and write the list of results to text log
$ rt_process_vs read --input_db output.db --eworst -6 --vdw_interactions A:VAL:243: --output_log hits.txt
```

Same example but using the API:

```python
from ringtail import RingtailCore

rtc = RingtailCore()
rtc.add_results_from_files(
    docking_results="results_folder/",
    recursive=True,
    receptor_file="receptor.json",
    save_receptor=True,
)

rtc.filter(
    eworst=-6,
    vdw_interactions=[("A:VAL:243:", True)],
    output_log="hits.txt",
)
```


### Upgrading older databases

Upgrade databases from Ringtail v2 or earlier with `rt_upgrade_db -d old_database.db` (note that this will remove existing filters and bookmarks). See the [upgrade guide](https://ringtail.readthedocs.io/en/latest/upgrade_database.html) and [v3 changelog](https://ringtail.readthedocs.io/en/latest/changes.html) for database migration, renamed options, and removal of built-in plotting and PyMOL integration. 


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
