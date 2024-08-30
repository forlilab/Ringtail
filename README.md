![ringtail logo final](https://user-images.githubusercontent.com/41704502/170797800-53a9d94a-932e-4936-9bea-e2d292b0c62b.png)

(Original artwork by Althea Hansel-Harris)


# Ringtail
Package for creating SQLite database from virtual screening results, performing filtering, and exporting results. Compatible with [AutoDock-GPU](https://github.com/ccsb-scripps/AutoDock-GPU) and [AutoDock-Vina](https://github.com/ccsb-scripps/AutoDock-Vina).

[![AD compat](https://img.shields.io/badge/AutoDock_Compatibility-ADGPU|Vina-brightgreen)](https://shields.io/)
[![License: L-GPL v2.1](https://img.shields.io/badge/License-LGPLv2.1-blue.svg)](https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html)
[![made-with-python](https://img.shields.io/badge/Made%20with-Python-1f425f.svg)](https://www.python.org/)
![Build Status](https://github.com/forlilab/Ringtail/actions/workflows/python-package.yml/badge.svg?event=push)

Ringtail reads collections of Docking Log File (DLG) or PDBQT results from virtual screenings performed with [AutoDock-GPU](https://github.com/ccsb-scripps/AutoDock-GPU) and [AutoDock-Vina](https://github.com/ccsb-scripps/AutoDock-Vina), respectively, and deposits them into
a SQLite database. It then allows for the filtering of results with numerous pre-defined filtering options, generation of a simple result scatterplot, export of 
molecule SDFs, and export of CSVs of result data. Result file parsing is parallelized across the user's CPU.

The publication describing the design, implementation, and features of Ringtail may be found in the JCIM paper:

[_Ringtail: A Python Tool for Efficient Management and Storage of Virtual Screening Results._
Althea T. Hansel-Harris, Diogo Santos-Martins, Niccolò Bruciaferri, Andreas F. Tillack, Matthew Holcomb, and Stefano Forli.
_Journal of Chemical Information and Modeling_ **2023** 63 (7), 1858-1864.
DOI: 10.1021/acs.jcim.3c00166](https://pubs.acs.org/doi/full/10.1021/acs.jcim.3c00166)

If using Ringtail in your work, please cite this publication.

Ringtail is developed by the [Forli lab](https://forlilab.org/) at the
[Center for Computational Structural Biology (CCSB)](https://ccsb.scripps.edu)
at [Scripps Research](https://www.scripps.edu/).

### New in version 2.0
##### Changes in keywords used for the command line tool

- `--mode` is now `--docking_mode`
- `--summary` is now `--print_summary`
- `--pattern` is now `--file_pattern`
- `--name` is now `--ligand_name`
- `--max_nr_atoms` is now `--ligand_max_atoms`
- `--smarts` is now `--ligand_substruct`
- `--smarts_idxyz` is now `--ligand_substruct_pos`
- `--smarts_join` is now `--ligand_operator`
- `--van_der_waals` is now `--vdw_interactions`
- `--hydrogen_bond` is now `--hb_interactions`
- `--reactive_res` is now `--reactive_interactions`

##### Enhancements to the codebase
- Fully developed API can use python for scripting exclusively 
- Can add docking results directly without using file system (for vina only as output comes as a string). 
- The Ringtail log is now written to a logging file in addition to STDOUT

##### Changes to code behavior
- Interaction tables: one new table has been added (`Interactions`) which references the interaction id from `Interaction_indices`, while the table `Interaction_bitvectors` has been discontinued.
- A new method to update an existing database 1.1.0 (or 1.0.0) to 2.0.0 is included. However, if the existing database was created with the duplicate handling option, there is a chance of inconsistent behavior of anything involving interactions as the Pose_ID was not used as an explicit foreign key in db v1.0.0 and v1.1.0 (see Bug fixes below).

##### Bug fixes
- The option `duplicate_handling` could previously only be applied during database creation and produced inconsistent table behavior. Option can now be applied at any time results are added to a database, and will create internally consistent tables. **Please note: if you have created tables in the past and invoking the keyword `duplicate_handling` you may have errors in the "Interaction_bitvectors" table. These errors cannot be recovered, and we recommend you re-make the database with Ringtail 2.0.0.**
- Writing SDFs from filtering bookmarks: will check that bookmark exists and has data before writing, and will now produce SDFs for any bookmarks existing bookmarks. If the bookmark results from a filtering where `max_miss` &lt; 0 it will note if the non-union bookmark is used, and if the base name for such bookmarks is provided it will default to the `basename_union` bookmark for writing the SDFs.
- Output from filtering using `max_miss` and `output_all_poses=False`(default) now producing expected behavior of outputting only one pose per ligand. Filtering for interactions `max_miss` allows any given pose for a ligand to miss `max_miss` interactions and still be considered to pass the filter. Previously, in the resulting `union` bookmark and `output_log` text file some ligands would present with more than one pose, although the option to `output_all_poses` was `False` (and thus the expectation would be one pose outputted per ligand). This would give the wrong count for how many ligands passed a filter, as some were counted more than once. 

#### Updating database to work with v2.0.0
If you have previously written a database with Ringtail < v2.0.0, it will need to be updated to be compatible with filtering with v2.0.0. We have included a new script `rt_db_to_v200.py` to perform this updated. Please note that all existing bookmarks will be removed during the update. The usage is as follows:

```
$ rt_db_to_v200.py -d <v2.0.0 database 1 (required)> <v2.0.0 database 2+ (optional)>
```

Multiple databases may be specified at once. The update may take a few minutes per database.

### New in version 1.1:
Code base and database schema version update
- [Significant filtering runtime improvements vs v1.0](https://github.com/forlilab/Ringtail/#example-filtering-timings-m1pro-macbook-2-million-ligands)
- `--summary` option for getting quick overview of data across entire dataset
- Selection of dissimilar output ligands with Morgan fingerprint or interaction fingerprint clustering
- Select similar ligands from query ligand name in previous Morgan fingerprint or interaction finger clustering groups
- Option for exporting stored receptor PDBQTs
- Filter by ligand substructure
- Filter by ligand substructure location in cartesian space
- `--max_miss` option now outputs union of interaction combinations by default, with `--enumerate_interaction_combs` option to log passing ligands/poses for individual interaction combination


##### Example Filtering Timings (M1Pro MacBook, ~2 million ligands)
![rt_v11_timings](https://github.com/forlilab/Ringtail/assets/41704502/eac373fc-1324-45df-b845-6697dc9d1465)

#### Updating database written with v1.0.0 to work with v1.1.0
If you have previously written a database with Ringtail v1.0.0, it will need to be updated to be compatible with filtering with v1.1.0. We have included a new script `rt_db_v100_to_v110.py` to perform this updated. Please note that all existing bookmarks will be removed during the update. The usage is as follows:

```
$ rt_db_v100_to_v110.py -d <v1.0.0 database 1 (required)> <v1.0.0 database 2+ (optional)>
```

Multiple databases may be specified at once. The update may take a few minutes per database.

### Dependencies 
- python (> 3.9, tested up to 3.12)
- RDKit
- SciPy
- Matplotlib
- Pandas
- chemicalite
- [Meeko](https://github.com/forlilab/Meeko) (from the Forli Lab)
- [Multiprocess](https://pypi.org/project/multiprocess/)

## README Outline
- [Installation](https://github.com/forlilab/Ringtail#installation)
- [Definitions](https://github.com/forlilab/Ringtail#definitions)
- [Getting Started Tutorial](https://github.com/forlilab/Ringtail#getting-started)
- [Scripts](https://github.com/forlilab/Ringtail#scripts)
- [rt_process_vs.py Documentation](https://github.com/forlilab/Ringtail#rt_process_vspy-documentation)
- [rt_compare.py Documentation](https://github.com/forlilab/Ringtail#rt_comparepy-documentation)
- [Python tutorials](https://github.com/forlilab/Ringtail#brief-python-tutorials)

### Installation (from PyPI)
Please note that Ringtail requires Python 3.9 or 3.10.

```bash
$ pip install ringtail
```
If using conda, `pip` installs the package in the active environment.

Also note that if using MacOS, you may need to install Multiprocess separately:
```bash
$ pip install multiprocess
```

### Installation (from source code)
```
$ conda create -n ringtail python=3.10
$ conda activate ringtail
```
After this, navigate to the desired directory for installing Ringtail and do the following:
```
$ git clone git@github.com:forlilab/Ringtail.git
$ cd Ringtail
$ pip install .
```
This will automatically fetch the required modules and install them into the current conda environment.

If you wish to make the code for Ringtail editable without having to re-run `pip install .`, instead use
```
$ pip install --editable .
```
### Test installation
If you would like to test your installation of Ringtail, a set of automated tests are included with the source code. To begin, you must install pytest in the Ringtail conda environment:
```
$ pip install -U pytest
```
Next, navigate to the `test` subdirectory within the cloned Ringtail directory and run pytest by simply calling
```
$ pytest
```
The compounds used for the testing dataset were taken from the [NCI Diversity Set V](https://wiki.nci.nih.gov/display/NCIDTPdata/Compound+Sets). The receptor used was [PDB: 4J8M](https://www.rcsb.org/structure/4J8M).

## Definitions
- __DLG__: Docking Log File, output from AutoDock-GPU.
- __PDBQT__: Modified PDB format, used for receptors (input to AutoDock-GPU and Vina) and output ligand poses from AutoDock-Vina.
- __Cluster__: Each docking result contains a number of independent runs, usually 20-50. These independent poses are then clustered by RMSD, giving groups of similar poses called clusters.
- __Pose__: The predicted ligand shape and position for single run of a single ligand in a single receptor.
- __Docking score__: The predicited binding energy from AutoDock-GPU or Vina.
- __Bookmark__: The set of ligands or ligand poses from a virtual screening passing a given set of filters. Stored within a virtual screening database as a view.
- __Ringtail__: 
> Drat, I'm not a cat!  Even though this eye-catching omnivore sports a few vaguely feline characteristics such as pointy ears, a sleek body, and a fluffy tail, the ringtail is really a member of the raccoon family. https://animals.sandiegozoo.org/animals/ringtail

## Getting started with the command line interface 
The Ringtail command line interface is orchestrated through the script `rt_process_vs.py`.
#### Create and populate a database
Navigate to the directory containing the data, in our case test_data:
```
$ cd test/test_data/
```
To write to the database we need to specify a few things:
- that we are using `write` mode
- source of docking results files. Docking results can be added either by providing one or more single files, a .txt file containing files, or by providing a directory containing docking results files.
- optional database name: ringtail will default to creating a database of name `output.db`
- optional docking mode: ringtail will default to assuming the files were produced by Autodock-GPU, if they are from vina specify `--mode vina`

Let us add all docking files within the path test_data (specified by `.` meaning current directory), whose folders we can traverse recursively by specifying `--recursive`
```
$ rt_process_vs.py write --file_path . --recursive
```
We can print a summary of the contents of the database by using the optional tag `-su` or `--summary` and specifying the database database from which to `read`:
```
$ rt_process_vs.py read --input_db output.db -su

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
```
$ rt_process_vs.py read --input_db output.db --eworst -6
```
This produces an output log `output_log.txt` with the names of ligands passing the filter, as well as their binding energies. Each round of filtering is also stored in the database as a SQLite view, which we refer to as a "bookmark" (default value is `passing_results`). 

We can also save a round of filtering with a specific bookmark name, and perform more filtering on this bookmark.
For example, start out with filtering out the compounds that are within the 5th percentile in terms of docking score and save the bookmark as `ep5`:
```
$ rt_process_vs.py read --input_db output.db --score_percentile 5 --log_file ep5_log.txt --bookmark_name ep5
```
Let's then further refine the set of molecules by applying an interaction filter for van der Waals interactions with V279 on the receptor:

```
$ rt_process_vs.py read --input_db output.db --filter_bookmark ep5 --vdw_interactions A:VAL:279: --log_file ep5_vdwV279_log.txt --bookmark_name ep5_vdwV279
```
The filtered molecules can then be exported as an e.g., SDF file which can be used for visual inspection in molecular graphics programs. At the same time, if pymol is installed, we can kick off a pymol session of the ligands

```
$ rt_process_vs.py read --input_db output.db --bookmark_name ep5_vdwV279 --export_sdf_path ep5_vdwV279_sdfs --pymol
```
#### Access help message for rt_process_vs.py
```
$ rt_process_vs.py --help
```
#### Access help message for rt_process_vs.py write mode
```
$ rt_process_vs.py write --help
```
#### Access help message for rt_process_vs.py read mode
```
$ rt_process_vs.py read --help
```

#### Ringtail arguments

| Argument         || Description                                           | Default value   | Requires interactions |
|:------------------------|:-----|:-------------------------------------------------|:----------------|----:|
|--config           | -c| Configuration JSON file to specify new default options. Overridded by command line | no default       |<tr><td colspan="5"></td></tr>
|--input_db         | -i| Database file to use instead of creating new database | no default       ||
|--bookmark_name      |-s| Name for bookmark view in database                      | passing_results  ||
|--mode          |-m| specify AutoDock program used to generate results. Available options are "dlg" and "vina". Vina mode will automatically change --pattern to \*.pdbqt   | dlg         ||
|--summary          |-su| Print summary information about virtual screening data to STDOUT. | FALSE        ||
|--verbose          |-v| Flag indicating that passing results should be printed to STDOUT. Will also include information about runtime progress. | FALSE        ||
|--debug            |-d| Flag indicating that additional debugging information (e.g. error traceback) should be printed to STDOUT. | FALSE |<tr><td colspan="5">**Write Mode**</td></tr>
|--file             |-f| DLG/Vina PDBQT file(s) to be read into database                  | no default       ||
|--file_path        |-fp| Path(s) to files to read into database            | no default       ||
|--file_list        |-fl| File(s) with list of files to read into database  | no default       ||
|--pattern          |-p| Specify pattern to search for when finding files   | \*.dlg\* / \*.pdbqt\* (vina mode)        ||
|--recursive        |-r| Flag to perform recursive subdirectory search on --file_path directory(s)  | FALSE      ||
|--append_results      |-a| Add new docking files to existing database given with --input_db  | FALSE       ||
|--duplicate_handling|-dh| Specify how dulicate results should be handled. May specify "ignore" or "replace". Unique results determined from ligand and target names and ligand pose. *NB: use of duplicate handling causes increase in database writing time*| None |
|--save_receptor    |-sr| Flag to specify that receptor file should be imported to database. Receptor file must also be specified with --receptor_file| FALSE   ||
|--output_db        |-o| Name for output database                              | output.db        ||
|--overwrite        |-ov| Flag to overwrite existing database           | FALSE       ||
|--max_poses        |-mp| Number of clusters for which to store top-scoring pose (dlg) or number of poses (vina) to save in database| 3     ||
|--store_all_poses  |-sa| Flag to indicate that all poses should be stored in database| FALSE      ||
|--interaction_tolerance|-it| Adds the interactions for poses within some tolerance RMSD range of the top pose in a cluster to that top pose. Can use as flag with default tolerance of 0.8, or give other value as desired | FALSE -> 0.8 (Å)  | Yes |
|--add_interactions  |-ai| Find interactions between ligands and receptor. Requires receptor PDBQT to be written. | FALSE      ||
|--interaction_cutoffs  |-ic| Specify distance cutoffs for measuring interactions between ligand and receptor in angstroms. Give as string, separating cutoffs for hydrogen bonds and VDW with comma (in that order). E.g. '-ic 3.7,4.0' will set the cutoff for hydrogen bonds to 3.7 angstroms and for VDW to 4.0. | 3.7,4.0     ||
|--receptor_file |-rn| Use with --save_receptor and/or --add_interactions. Give receptor PDBQT. | None      ||
|--max_proc |-mpr| Maximum number of subprocesses to spawn during database writing. | [# available CPUs]      |<tr><td colspan="5">**Read Mode**</td></tr>
|--log_file              |-l| Name for log of filtered results                      | output_log.txt   ||
|--outfields       |-of| Data fields to be written in output (log file and STDOUT). Ligand name always included. | e        ||
|--order_results    |-ord| String for field by which the passing results should be ordered in log file. | no default ||
|--output_all_poses        |-ap| Flag that if mutiple poses for same ligand pass filters, log all poses | (OFF)        ||
|--export_bookmark_csv |-xs| Name of database result bookmark or table to be exported as CSV. Output as <table_name>.csv | no default      ||
|--export_query_csv |-xq| Create csv of the requested SQL query. Output as query.csv. MUST BE PRE-FORMATTED IN SQL SYNTAX e.g. SELECT [columns] FROM [table] WHERE [conditions] | no default      ||
|--export_sdf_path|-sdf| Path for saving exported SDF files of ligand poses passing given filtering criteria | no default       ||
|--export_bookmark_db |-xdb| Export a database containing only the results found in the bookmark specified by --bookmark_name. Will save as <input_db>_<bookmark_name>.db| FALSE      ||
|--data_from_bookmark |-nd| Flag that out_fields data should be written to log for results in given --bookmark_name. Requires no filters. | FALSE       ||
|--filter_bookmark |-fb| Filter over specified bookmark, not whole Results table. | FALSE       ||
|--find_similar_ligands |-fsl| Given query ligand name, find ligands previously clustered with that ligand. User prompted at runtime to choose cluster group of interest. | no default       ||
|--plot             |-p| Flag to create scatterplot of ligand efficiency vs docking score for best pose of each ligand. Saves as [filters_file].png or out.png. | FALSE        ||
|--pymol             |-py| Flag to launch interactive LE vs Docking Score plot and PyMol session. Ligands in the bookmark specified with --bookmark_name will be ploted and displayed in PyMol when clicked on.| FALSE        |<tr><td colspan="5">PROPERTY FILTERS</td></tr>
|--eworst           |-e| Worst energy value accepted (kcal/mol)                | no default  ||
|--ebest            |-eb| Best energy value accepted (kcal/mol)                 | no default  ||
|--leworst          |-le| Worst ligand efficiency value accepted                | no default  ||
|--lebest           |-leb| Best ligand efficiency value accepted                 | no default  ||
|--score_percentile      |-pe| Worst energy percentile accepted. Give as percentage (1 for top 1%, 0.1 for top 0.1%) | 1.0  ||
|--le_percentile   |-ple| Worst ligand efficiency percentile accepted. Give as percentage (1 for top 1%, 0.1 for top 0.1%) | no default |  <tr><td colspan="5">LIGAND FILTERS</td></tr>
|--ligand_name             |-n| Search for specific ligand name. Multiple names joined by "OR". Multiple filters should be separated by commas | no default  ||
|--ligand_max_atoms     |-mna| Specify maximum number of heavy atoms a ligand may have | no default  ||
|--ligand_substruct           || SMARTS pattern(s) for substructur matching | no default  ||
|--ligand_substruct_pos     || SMARTS pattern, index of atom in SMARTS, cutoff distance, and target xyz coordinates. Finds poses in which the specified substructure atom is within the distance cutoff from the target location | no default  ||
|--ligand_operator     |-n| logical operator for multiple SMARTS | OR  | <tr><td colspan="5">INTERACTION FILTERS</td></tr>
|--vdw_interactions    |-vdw| Filter for van der Waals interaction with given receptor information.  | no default  | Yes|
|--hb_interactions    |-hb| Filter with hydrogen bonding interaction with given information. Does not distinguish between donating or accepting | no default  | Yes|
|--reactive_interactions     |-r| Filter for reation with residue containing specified information | no default  |Yes |
|--hb_count         |-hc| Filter for poses with at least this many hydrogen bonds. Does not distinguish between donating and accepting | no default  | Yes|
|--react_any        |-ra| Filter for poses with reaction with any residue       | FALSE     | Yes|
|--max_miss         |-mm| Will filter given interaction filters excluding up to max_miss interactions. Results in ![equation](https://latex.codecogs.com/svg.image?\sum_{m=0}^{m}\frac{n!}{(n-m)!*m!}) combinations for *n* interaction filters and *m* max_miss. Will log and output union of combinations unless used with `--enumerate_interaction_combs`. | 0  | Yes |
|--enumerate_interactions_combs  |-eic| When used with `--max_miss` > 0, will log ligands/poses passing each separate interaction filter combination as well as union of combinations. Can significantly increase runtime. | FALSE  | Yes <tr><td colspan="5">PASSING RESULT CLUSTERING</td></tr>
|--mfpt_cluster     |-mfpc| Cluster ligands passing given filters based on the Tanimoto distances of the Morgan fingerprints. Will output ligand with best (lowest) ligand efficiency from each cluster. Uses Butina clustering algorithm | 0.5  ||
|--interaction_cluster     |-ifpc| Cluster ligands passing given filters based on the Tanimoto distances of the interaction fingerprints. Will output ligand with best (lowest) ligand efficiency from each cluster. Uses Butina clustering algorithm | 0.5  | Yes |

---

### Scripts
The Ringtail package includes two command line oriented scripts: `rt_process_vs.py` and `rt_compare.py`. Both may be run with options specified in the command line and/or using options specified in a JSON-formatted file given with `--config`. Command line options override any conflicting options in the config file.

[rt_process_vs.py](https://github.com/forlilab/Ringtail#rt_process_vspy-documentation) serves as the primary script for the package and is used to both write docking files to a SQLite database and to perform filtering and export tasks on the database. It is designed to handle docking output files associated with a single virtual screening in a single database.

[rt_compare.py](https://github.com/forlilab/Ringtail#rt_comparepy-documentation) is used to combine information across multiple virtual screenings (in separate databases) to allow or exclude the selection of ligands passing filters across multiple targets/models. This can be useful for filtering out promiscuous ligands, a technique commonly used in exerimental high-throughput screening. It may also be used if selection of ligands binding multiple protein structures/conformations/homologs are desired.

[rt_generate_config_file.py](https://github.com/forlilab/Ringtail#rt_generate_config_filepy-documentation) can be ran to create a config file template

[rt_db_to_v200.py](https://github.com/forlilab/Ringtail#Updating-database-to-work-with-v200) is used to update older databases to the latest version. 

[rt_db_v100_to_v110.py](https://github.com/forlilab/Ringtail#Updating-database-written-with-v100-to-work-with-v110) is used to update db v1.0.0 to 1.1.0. 

#### rt_compare.py Documentation
The `rt_compare.py` script is designed to be used with databases already made and filtered. The script is used to select ligands which are shared between the given filter bookmark(s) of some virtual screenings (wanted) or exclusive to some screenings and not others (unwanted). The script uses a subset of commands similar to `rt_process_vs.py`.

An example of use: select ligands found in "filter_bookmark" bookmarks of database1 but not database2 (they must both contain a bookmark named "filter1"):
```
rt_compare.py --wanted database1.db --unwanted database2.db --bookmark_name filter_bookmark
```

For more detailed description of usage, please see [the readthedocs.org site for ringtail](https://ringtail.readthedocs.io/en/latest/compare.html).

#### rt_generate_config_file.py Documentation


## Advanced usage: scripting with Ringtail API 
Ringtail has been re-designed to allow for direct use of its API for e.g., scripting purposes. This circumvents the use of the command line tools, and allows for more advanced usage.
The available operations and keywords are the same as for the command line interface, but the methods can now be accessed at a more granular level if desired. For docking engines that provides direct string output such as Vina, it is also possible to save the docking results output directly to the database as a string and thus circumventing use of the computer file system (some link to vina scripting, probably in readthedocs).

#### Instantiating the Ringtail object
A ringtail core is created by instantiating a `RingtailCore` object with a database. Currently, a database can only be added upon instantiation.
```
rtc = RingtailCore("output.db")
```

Default logging level is "WARNING", and a different logger level can be set at the time of object instantiation, or later by the log level change API:
```
rtc = RingtailCore(db_file="output.db", logging_level="DEBUG)
# or
rtc.logger.set_level("INFO")
```

#### Populate the database
To add results to the database, use the `add_results_from_files` method that takes as input files and directories to upload,
as well as a receptor path and database properties and how to handle the resutls (how many poses to save, how to deal with interactions if having vina results),
and whether or not to print a summary after writing the results to the database.

```
rtc.add_results_from_files( file_path = "test_data/", 
                            recursive = True, 
                            save_receptor = False,
                            max_poses = 3)
```
Both files (`filesources_dict`) and processing options (`optionsdict`) can be provided as dictionaries as well or instead of the the individual options. Any provided individual options will overwrite the options provided through dictionaries. The use and prioritization of dictionaries and method attributes is true for most of the available API methods.

```
file_sources = {
    "file_path": "test_data/",
    "recursive": True,
}

writeoptions = {
    "store_all_poses": True,
    "max_proc": 4
}

rtc.add_results_from_files( filesources_dict = file_sources,
                            optionsdict = writeoptions,)
```

If at any point you wish to print a summary of the database, the method can be called directly:
```
rtc.produce_summary()
```

The default docking mode is "dlg", and can be changed to "vina" by accessing the ringtail core property `docking_mode`. 
```
rtc_vina = RingtailCore("output_vina.db")
rtc_vina.docking_mode = "vina"
```
Since vina does not automatically write docking results to the file system, these can be added to the database by associating them with a ligand name in a dictionary and using this dictionary as the source of results when adding to the database:
```
vina_docking_result1 = "long string of results"
vina_docking_result2 = "different string of results"

vina_results = {
    "ligand1": vina_docking_result1,
    "ligand2": vina_docking_result2
}

rtc_vina.add_results_from_vina_string(results_strings = vina_results,
                                 max_poses = 2)
```

#### Filtering and visualizing the data in the database

To filter, simply access the API method `filter` and provide desired filter values. Names of bookmark and output log for containing filtered results can be specified in the method.
```
rtc.filter(eworst=-6, 
           bookmark_name = "e6",
           log_file = "filtered_results.txt")
```
Just like with the command line tool, you can choose to filter over a bookmark that has already been created:
```
rtc.filter(vdw_interactions=[('A:VAL:279:', True), ('A:LYS:162:', True)],
           bookmark_name = "e6vdw279162",
           filter_bookmark = "e6",
           log_file = "filtered_results_2.txt")
```
To export filtered molecules in a specific bookmark to SDF files use the following method, where the `sdf_path` directory will be created if it does not already exist:

```
rtc.write_molecule_sdfs(sdf_path = "sdf_files", bookmark_name = "e6vdw279162")
```

One or more of the filtered ligands can be visualized in PyMol:
```
rtc.pymol(bookmark_name = "e6vdw279162")
```

### Arguments used for API vs command line
All of the arguments used for the command line tool also applies to the Ringtail API in some form. For example, bookmark names and filter values are provided when an API method is called, while the log level can be sat at instantiation or at any time during the scripting process. Instead of differentiating between an `--input_db` and `--output_db`, only one database file is operated on in a given instantiated `RingtailCore` object. A subset of the command line arguments are actual API methods (such as `--plot` or `--find_similar_ligands`) that will be called directly, with optional input arguments (typically a `bookmark_name` or `ligand_name`). Each API method comes with type hints and extensive documentation. Additionally, a more extensive example of its use can be found on [readthedocs](https://ringtail.readthedocs.io/en/latest/). 

