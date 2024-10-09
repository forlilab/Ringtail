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

### New in version 2.0.0 and 2.0.1

- changes in keywords used for the command line tool
- fully developed API can add docking results without using file system (for Vina)
- much faster filtering (v2.1.0)
- bug fixes
- see detailed list of changes on [ReadTheDocs](/https://ringtail.readthedocs.io)

#### Updating database to work with v2.0
If you have previously written a database with Ringtail < v2.0, it will need to be updated to be compatible with filtering with v2.0. We have included a new script `rt_db_to_v200` to perform this updated. Please note that all existing bookmarks will be removed during the update. The usage is as follows:

```
$ rt_db_to_v200 -d <v2.0 database 1 (required)> <v2.0 database 2+ (optional)>
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
Ringtail 2.0 is now available on conda-forge, and installation from conda-forge will handle all of the dependencies. 

```bash
$ conda install -c conda-forge ringtail
```

## Getting started with the command line interface 

#### Create and populate a database
Navigate to the directory containing the data, in our case test\_data/adgpu:

```bash
$ cd test/test_data/adpgu/
```
To write to the database we need to specify a few things:
- that we are operating in `write` mode
- source of docking results files. Docking results can be added either by providing one or more single files, a .txt file containing files, or by providing a directory containing docking results files.
- optional database name: ringtail will default to creating a database of name `output.db`
- optional docking mode: ringtail will default to assuming the files were produced by Autodock-GPU, if they are from vina specify `--mode vina`

Let us add all docking files within the path test\_data (specified by `.` meaning current directory), whose folders we can traverse recursively by specifying `--recursive`

```bash
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
Let us start filtering with a basic docking score cutoff of -6 kcal/mol:

```bash
$ rt_process_vs read --input_db output.db --eworst -6
```

This produces an output log `output_log.txt` with the names of ligands passing the filter, as well as their binding energies. Each round of filtering is also stored in the database as a SQLite view, which we refer to as a "bookmark" (default value is `passing_results`). 

We can also save a round of filtering with a specific bookmark name, and perform more filtering on this bookmark.
For example, start out with filtering out the compounds that are within the 5th percentile in terms of docking score and save the bookmark as `ep5`:

```bash
$ rt_process_vs read --input_db output.db --score_percentile 5 --log_file ep5_log.txt --bookmark_name ep5
```

Let's then further refine the set of molecules by applying an interaction filter for van der Waals interactions with V279 on the receptor:

```bash
$ rt_process_vs read --input_db output.db --filter_bookmark ep5 --vdw_interactions A:VAL:279: --log_file ep5_vdwV279_log.txt --bookmark_name ep5_vdwV279
```

The filtered molecules can then be exported as an e.g., SDF file which can be used for visual inspection in molecular graphics programs. At the same time, if pymol is installed, we can kick off a pymol session of the ligands

```bash
$ rt_process_vs read --input_db output.db --bookmark_name ep5_vdwV279 --export_sdf_path ep5_vdwV279_sdfs --pymol
```

#### Access help message for rt\_process_vs
```bash
$ rt_process_vs --help

$ rt_process_vs write --help

$ rt_process_vs read --help
```

### Scripts
The Ringtail package includes two command line oriented scripts: `rt_process_vs` and `rt_compare`. Both may be run with options specified in the command line and/or using options specified in a JSON-formatted file given with `--config`. Command line options override any conflicting options in the config file.

[rt_process_vs](https://github.com/forlilab/Ringtail#rt_process_vspy-documentation) serves as the primary script for the package and is used to both write docking files to a SQLite database and to perform filtering and export tasks on the database. It is designed to handle docking output files associated with a single virtual screening in a single database.

[rt_compare](https://github.com/forlilab/Ringtail#rt_comparepy-documentation) is used to combine information across multiple virtual screenings (in separate databases) to allow or exclude the selection of ligands passing filters across multiple targets/models. This can be useful for filtering out promiscuous ligands, a technique commonly used in exerimental high-throughput screening. It may also be used if selection of ligands binding multiple protein structures/conformations/homologs are desired.

[rt_generate_config_file](https://github.com/forlilab/Ringtail#rt_generate_config_filepy-documentation) can be ran to create a config file template

[rt_db_to_v200](https://github.com/forlilab/Ringtail#Updating-database-to-work-with-v200) is used to update older databases to the latest version. 

[rt_db_v100_to_v110](https://github.com/forlilab/Ringtail#Updating-database-written-with-v100-to-work-with-v110) is used to update db v1.0.0 to 1.1.0. 

#### rt_compare Documentation
The `rt_compare` script is designed to be used with databases already made and filtered. The script is used to select ligands which are shared between the given filter bookmark(s) of some virtual screenings (wanted) or exclusive to some screenings and not others (unwanted). The script uses a subset of commands similar to `rt_process_vs`.

An example of use: select ligands found in "filter_bookmark" bookmarks of database1 but not database2 (they must both contain a bookmark named "filter1"):

```bash
rt_compare --wanted database1.db --unwanted database2.db --bookmark_name filter_bookmark
```

For more detailed description of usage, please see [the readthedocs.org site for ringtail](https://ringtail.readthedocs.io/en/latest/compare.html).

## Advanced usage: scripting with Ringtail API 
Ringtail has been re-designed to allow for direct use of its API for e.g., scripting purposes. This circumvents the use of the command line tools, and allows for more advanced usage.
The available operations and keywords are the same as for the command line interface, but the methods can now be accessed at a more granular level if desired. For docking engines that provides direct string output such as Vina, it is also possible to save the docking results output directly to the database as a string and thus circumventing use of the computer file system (some link to vina scripting, probably in readthedocs).

#### Instantiating the Ringtail object
A ringtail core is created by instantiating a `RingtailCore` object with a database. Currently, a database can only be added upon instantiation.

```bash
rtc = RingtailCore("output.db")
```

Default logging level is "WARNING", and a different logger level can be set at the time of object instantiation, or later by the log level change API:

```bash
rtc = RingtailCore(db_file="output.db", logging_level="DEBUG)
# or
rtc.logger.set_level("INFO")
```

#### Populate the database
To add results to the database, use the `add_results_from_files` method that takes as input files and directories to upload,
as well as a receptor path and database properties and how to handle the resutls (how many poses to save, how to deal with interactions if having vina results),
and whether or not to print a summary after writing the results to the database.

```python
rtc.add_results_from_files( file_path = "test_data/", 
                            recursive = True, 
                            save_receptor = False,
                            max_poses = 3)
```

If at any point you wish to print a summary of the database, the method can be called directly:

```python
rtc.produce_summary()
```

The default docking mode is "dlg", and can be changed to "vina" by accessing the ringtail core property `docking_mode`. 

```python
rtc_vina = RingtailCore("output_vina.db")
rtc_vina.docking_mode = "vina"
```

Since vina does not automatically write docking results to the file system, these can be added to the database by associating them with a ligand name in a dictionary and using this dictionary as the source of results when adding to the database:

```python

vina_results = {
    "ligand1": vina_docking_ligand1_result,
    "ligand2": vina_docking_ligand2_result
}

rtc_vina.add_results_from_vina_string(results_strings = vina_results,
                                 max_poses = 2)
```

#### Filtering and visualizing the data in the database

To filter, simply access the API method `filter` and provide desired filter values. Names of bookmark and output log for containing filtered results can be specified in the method.

```python
rtc.filter(eworst=-6, 
           bookmark_name = "e6",
           log_file = "filtered_results.txt")
```

Just like with the command line tool, you can choose to filter over a bookmark that has already been created:

```python
rtc.filter(vdw_interactions=[('A:VAL:279:', True), ('A:LYS:162:', True)],
           bookmark_name = "e6vdw279162",
           filter_bookmark = "e6",
           log_file = "filtered_results_2.txt")
```

To export filtered molecules in a specific bookmark to SDF files use the following method, where the `sdf_path` directory will be created if it does not already exist:

```python
rtc.write_molecule_sdfs(sdf_path = "sdf_files", bookmark_name = "e6vdw279162")
```

One or more of the filtered ligands can be visualized in PyMol:

```python
rtc.pymol(bookmark_name = "e6vdw279162")
```

### Arguments used for API vs command line
All of the arguments used for the command line tool also applies to the Ringtail API in some form. For example, bookmark names and filter values are provided when an API method is called, while the log level can be sat at instantiation or at any time during the scripting process. Instead of differentiating between an `--input_db` and `--output_db`, only one database file is operated on in a given instantiated `RingtailCore` object. A subset of the command line arguments are actual API methods (such as `--plot` or `--find_similar_ligands`) that will be called directly, with optional input arguments (typically a `bookmark_name` or `ligand_name`). Each API method comes with type hints and extensive documentation. Additionally, a more extensive example of its use can be found on [readthedocs](https://ringtail.readthedocs.io/en/latest/). 

