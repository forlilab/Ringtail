![ringtail logo final](https://user-images.githubusercontent.com/41704502/170797800-53a9d94a-932e-4936-9bea-e2d292b0c62b.png)

(Original artwork by Althea Hansel-Harris)


# Ringtail
Package for creating SQLite database from virtual screening results, performing filtering, and exporting results. Compatible with [AutoDock-GPU](https://github.com/ccsb-scripps/AutoDock-GPU) and [AutoDock-Vina](https://github.com/ccsb-scripps/AutoDock-Vina).

[![License: L-GPL v2.1](https://img.shields.io/badge/License-LGPLv2.1-blue.svg)](https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html)
![Build Status](https://github.com/forlilab/Ringtail/actions/workflows/python-package.yml/badge.svg?event=push)

Ringtail reads collections of Docking Log File (DLG) from virtual screenings performed with [AutoDock-GPU](https://github.com/ccsb-scripps/AutoDock-GPU)
or PDBQT results from [AutoDock-Vina](https://github.com/ccsb-scripps/AutoDock-Vina), and deposits them into
a SQLite database. It then allows for the filtering of results with numerous pre-defined options, generation of a simple result scatterplot, export of 
molecule SDFs, and export of CSVs of result data. Parsing of output files from docking is parallelized across the user's CPU.

The publication describing the design, implementation, and features of Ringtail may be found in the JCIM paper:

[_Ringtail: A Python Tool for Efficient Management and Storage of Virtual Screening Results._
Althea T. Hansel-Harris, Diogo Santos-Martins, Niccolò Bruciaferri, Andreas F. Tillack, Matthew Holcomb, and Stefano Forli.
_Journal of Chemical Information and Modeling_ **2023** 63 (7), 1858-1864.
DOI: 10.1021/acs.jcim.3c00166](https://pubs.acs.org/doi/full/10.1021/acs.jcim.3c00166)

If using Ringtail in your work, please cite this publication.

Ringtail is developed by the [Forli lab](https://forlilab.org/)
at [Scripps Research](https://www.scripps.edu/).

### Documentation
In-depth documentation can be found on [ReadTheDocs](https://ringtail.readthedocs.io).

### New in version 2.0.0 and 2.1.0

- changes in keywords used for the command line tool
- fully developed API can add docking results without using file system (for Vina)
- much faster filtering (v2.1.0)
- bug fixes
- see detailed list of changes on [ReadTheDocs](/https://ringtail.readthedocs.io)

#### Updating database to work with v2.0
If you have previously written a database with Ringtail < v2, it will need to be updated to be compatible with filtering with v2. We have included a new script `rt_db_to_v200` to perform this updated. Please note that all existing bookmarks will be removed during the update. The usage is as follows:

```
$ rt_db_to_v200 -d <v2 database 1 (required)> <v2 database 2+ (optional)>
```

Multiple databases may be specified at once. The update may take a few minutes per database.

##### Example Filtering Timings (M1Pro MacBook, ~2 million ligands)
![rt_v11_timings](https://github.com/forlilab/Ringtail/assets/41704502/eac373fc-1324-45df-b845-6697dc9d1465)

#### Updating database written with v1.0.0 to work with v1.1.0
If you have previously written a database with Ringtail v1.0.0, it will need to be updated to be compatible with filtering with v1.1.0. We have included a new script `rt_db_v100_to_v110` to perform this updated. Please note that all existing bookmarks will be removed during the update. The usage is as follows:

```
$ rt_db_v100_to_v110 -d <v1.0.0 database 1 (required)> <v1.0.0 database 2+ (optional)>
```

### Installation 
#### Create a Ringtail environment
It is necessary to create a Ringtail python environment for managing the external dependencies, conda will be used in the following examples but other environment managers such as the lightweight micromamba will also work. Please note that Ringtail requires Python 3.9, 3.10, or 3.11.

```bash
$ conda create -n Ringtail python=3.11
$ conda activate ringtail
```

#### From PyPi
Make sure your Ringtail environment is active, then install via pip

```bash
$ pip install ringtail
```

You may need to install one or more of the listed dependencies, please note that multiprocess is only necessary for MacOS. 

```bash
$ pip install <dependency>
```

Chemicalite is required and only available on conda-forge:

```bash
$ conda install -c conda-forge chemicalite
```

#### From conda-forge
Ringtail 2 is now available on conda-forge, and installation from conda-forge will handle all of the dependencies. 

```bash
$ conda install -c conda-forge ringtail
```

## Quick-start for the command line interface 
#### Create and populate a database
To write to the database we need to specify a few things:
- that we are operating in `write` mode
- source of docking results files. Docking results can be added either by providing one or more single files, a .txt file containing files, or by providing a directory containing docking results files.
- optional database name: ringtail will default to creating a database of name `output.db`
- optional docking mode: ringtail will default to assuming the files were produced by Autodock-GPU, if they are from vina specify `--mode vina`

Start by navigating to the directory containing the data, in our case test\_data/adgpu. Let us add all docking files within the path test\_data (specified by `.` meaning current directory), whose folders we can traverse recursively by specifying `--recursive`

```bash
$ cd test/test_data/adpgu/

$ rt_process_vs write --file_path . --recursive
```

We can print a summary of the contents of the database by using the optional tag `-su` or `--summary` and specifying the database database from which to `read`:

```bash
$ rt_process_vs read --input_db output.db -su

Total Stored Poses: 645
Total Unique Interactions: 183

Energy statistics:
min_docking_score: -7.93 kcal/mol
max_docking_score: -2.03 kcal/mol
1%_docking_score: -7.43 kcal/mol
10%_docking_score: -6.46 kcal/mol
min_leff: -0.62 kcal/mol
max_leff: -0.13 kcal/mol
1%_leff: -0.58 kcal/mol
10%_leff: -0.47 kcal/mol
```

#### Filtering and visualizing the data in the database
Multiple groups of filters are available in Ringtail, including simpler filters based on docking score and ligand efficiency, filters based on interactions with receptor residues, and filters based on the chemical characteristics of the ligand. Below is an example filtering with a basic docking score cutoff of -6 kcal/mol and finding any ligand that has a van der Waals interaction with a receptor valine that is residue no. 279 on the receptor:

```bash
$ rt_process_vs read --input_db output.db --eworst -6 --vdw_interactions A:VAL:279:
```

This produces an output log `output_log.txt` with the names of ligands passing the filter, as well as their binding energies. Each round of filtering is also stored in the database as a SQLite view, which we refer to as a "bookmark" (default value is `passing_results`). In addition to the easy-to-read text file, filtered ligands can be written to SDF files and visualized in PyMol (other output options like plotting are described on [ReadTheDocs](/https://ringtail.readthedocs.io)).

```bash
$ rt_process_vs read --input_db output.db --bookmark_name passing_results --export_sdf_path . --pymol
```

For quick access to the options in the command line interface, you can use the keyword `--help`:
```bash
$ rt_process_vs --help

$ rt_process_vs write --help

$ rt_process_vs read --help
```

### Scripts
The Ringtail package includes two command line oriented scripts: `rt_process_vs` and `rt_compare`. Both may be run with options specified in the command line and/or using options specified in a JSON-formatted file given with `--config`. Command line options override any conflicting options in the config file.

`rt_process_vs` serves as the primary script for the package and is used to both write docking files to a SQLite database and to perform filtering and export tasks on the database. It is designed to handle docking output files associated with a single virtual screening in a single database.

`rt_compare` is used to combine information across multiple virtual screenings (in separate databases) to allow or exclude the selection of ligands passing filters across multiple targets/models. This can be useful for filtering out promiscuous ligands, a technique commonly used in exerimental high-throughput screening. It may also be used if selection of ligands binding multiple protein structures/conformations/homologs are desired. For detailed description of usage, please see [the readthedocs.org site for ringtail](https://ringtail.readthedocs.io/en/latest/compare.html).

`rt_generate_config_file` can be ran to create a config file template

[`rt_db_to_v200`](https://github.com/forlilab/Ringtail#Updating-database-to-work-with-v200) is used to update older databases to version 2.1. 

[`rt_db_v100_to_v110`](https://github.com/forlilab/Ringtail#Updating-database-written-with-v100-to-work-with-v110) is used to update db v1.0.0 to 1.1.0. 

## Advanced usage: scripting with Ringtail API 
Ringtail has been re-designed to allow for direct use of its API for e.g., scripting purposes. This circumvents the use of the command line tools, and allows for more advanced usage.
The available operations and keywords are the same as for the command line interface, but the methods can now be accessed at a more granular level if desired. For docking engines that provides direct string output such as Vina, it is also possible to save the docking results output directly to the database as a string and thus circumventing use of the computer file system (some link to vina scripting, probably in readthedocs).

#### Create and populate a database
A ringtail core is created by instantiating a `RingtailCore` object with a database. Currently, a database can only be added upon instantiation. To add results to the database, use the `add_results_from_files` method that takes as input files and directories to upload, as well as a receptor path and database properties and how to handle the results (how many poses to save, how to deal with interactions if having vina results), and whether or not to print a summary after writing the results to the database.

```python
rtc = RingtailCore("output.db")
rtc.add_results_from_files( file_path = "test_data/", 
                            recursive = True, 
                            save_receptor = False,
                            max_poses = 3)
```

If at any point you wish to print a summary of the database, the method can be called directly:

```python
rtc.produce_summary()
```

#### Filtering and visualizing the data in the database

To filter, simply access the API method `filter` and provide desired filter values. Names of bookmark and output log for containing filtered results can be specified in the method.

```python
rtc.filter(eworst=-6, 
           vdw_interactions=[('A:VAL:279:', True)]
           bookmark_name = "passing_results",
           log_file = "filtered_results.txt")
```

To export filtered molecules in a specific bookmark to SDF files use the following method, where the `sdf_path` directory will be created if it does not already exist. Visualizing molecules in pymol is similarly accomplished by calling the `pymol` method.

```python
rtc.write_molecule_sdfs(sdf_path = ".", bookmark_name = "passing_results")

rtc.pymol(bookmark_name = "passing_results")
```

### Arguments used for API vs command line
All of the arguments used for the command line tool applies to the Ringtail API in some form. For example, bookmark names and filter values are provided when an API method is called, while the log level can be sat at instantiation or at any time during the scripting process. When using the API, instead of differentiating between an `--input_db` and `--output_db`, only one database file is operated on in a given instantiated `RingtailCore` object. A subset of the command line arguments are actual API methods (such as `--plot` or `--find_similar_ligands`) that will be called directly as methods, with optional input arguments (typically a `bookmark_name` or `ligand_name`). Each API method comes with type hints and extensive documentation. Additionally, extensive example of the use of both can be found on [readthedocs](https://ringtail.readthedocs.io/). 
