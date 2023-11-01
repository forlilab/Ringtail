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

### New in version 1.1:
- [Significant filtering runtime improvements vs v1.0](https://github.com/forlilab/Ringtail/#example-filtering-timings-m1pro-macbook-2-million-ligands)
- `--summary` option for getting quick overview of data across entire dataset
- Selection of dissimilar output ligands with Morgan fingerprint or interaction fingerprint clustering
- Select similar ligands from query ligand name in previous Morgan fingerprint or interaction finger clustering groups
- Option for exporting stored receptor PDBQTs
- Filter by ligand substructure
- Filter by ligand substructure location in cartesian space
- `--max_miss` option now outputs union of interaction combinations by default, with `--enumerate_interaction_combs` option to log passing ligands/poses for individual interaction combination

#### Updating database written with v1.0.0 to work with v1.1.0
If you have previously written a database with Ringtail v1.0.0, it will need to be updated to be compatible with filtering with v1.1.0. We have included a new script `rt_db_v100_to_v110.py` to perform this updated. Please note that all existing bookmarks will be removed during the update. The usage is as follows:

```
$ rt_db_v100_to_v110.py -d <v1.0.0 database 1 (required)> <v1.0.0 database 2+ (optional)>
```

Multiple databases may be specified at once. The update may take a few minutes per database.

## README Outline
- [Installation](https://github.com/forlilab/Ringtail#installation)
- [Definitions](https://github.com/forlilab/Ringtail#definitions)
- [Getting Started Tutorial](https://github.com/forlilab/Ringtail#getting-started)
- [Scripts](https://github.com/forlilab/Ringtail#scripts)
- [rt_process_vs.py Documentation](https://github.com/forlilab/Ringtail#rt_process_vspy-documentation)
- [rt_compare.py Documentation](https://github.com/forlilab/Ringtail#rt_comparepy-documentation)
- [Python tutorials](https://github.com/forlilab/Ringtail#brief-python-tutorials)

## Installation
It is recommended that you create a new Conda environment for installing Ringtail. Ringtail requires the following non-standard dependencies:
- RDKit
- SciPy
- Matplotlib
- Pandas
- chemicalite
- [Meeko](https://github.com/forlilab/Meeko) (from the Forli Lab)
- [Multiprocess](https://pypi.org/project/multiprocess/) (MacOS only)

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
#### Test installation
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
- __Cluster__: Each DLG contains a number of independent runs, usually 20-50. These independent poses are then clustered by RMSD, giving groups of similar poses called clusters.
- __Pose__: The predicted ligand shape and position for single run of a single ligand in a single receptor.
- __Docking score__: The predicited binding energy from AutoDock-GPU or Vina.
- __Bookmark__: The set of ligands or ligand poses from a virtual screening passing a given set of filters. Stored within a virtual screening database as a view.
- __Ringtail__: 
> Drat, I'm not a cat!  Even though this eye-catching omnivore sports a few vaguely feline characteristics such as pointy ears, a sleek body, and a fluffy tail, the ringtail is really a member of the raccoon family. https://animals.sandiegozoo.org/animals/ringtail

## Getting Started

Ringtail offers a wealth of database creation and filtering options. They are detailed at length below. This section will provide a quick overview of the basic usage of Ringtail from the command line. We will you the provided test data to create a database with default storage options and perform basic filtering of it.

Let us begin in the Ringtail directory. First, we must navigate to the test data directory:
```
$ cd test/test_data/
```
Now, let us create a database containing the results from only group 1. Note that these files are DLGs. If we were using Vina PDBQTs, we would need to add `--mode vina`.
```
$ rt_process_vs.py write --file_path group1
```
By default, the database we have created is called `output.db`. Let us now make a second database named `all_groups.db` containing all three groups:
```
$ rt_process_vs.py write --file_path . --recursive --output_db all_groups.db
```
The `--recursive` option tells Ringtail to scan the directories specified with `--file_path` for subdirectories containing output files (in this case, DLGs). This allowed all three group directories to be added to the database with a single --file_path option.

Now that we have created the databases, we can filter them to pull out compounds of interest. Before we do that, let's find out a little more about the data contained within the database. For this, we can use the `-s/--summary` option:
```
$rt_process_vs.py read --input_db all_groups.db -s

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
We could also have used the `--summary` option when writing the database to display this info at that time.

Now, let us start filtering with a basic docking score cutoff of -6 kcal/mol:
```
$ rt_process_vs.py read --input_db all_groups.db --eworst -6
```

This produces an output log `output_log.txt` with the names of ligands passing the filter, as well as their binding energies. Let's now do another round of filtering, this time with an energy percentile filter of 5 percent (top 5% of coumpounds by docking score). Each round of filtering is also stored in the database as a SQLite view, which we refer to as a "bookmark". We will also save this round of filtering with the bookmark name "ep5".

```
$ rt_process_vs.py read --input_db all_groups.db --score_percentile 5 --log ep5_log.txt --bookmark_name ep5
```
Now, let us further refine the set of molecules we just filtered. We will use an interaction filter for van der Waals interactions with V279 on the receptor:

```
$ rt_process_vs.py read --input_db all_groups.db --filter_bookmark ep5 --van_der_waals A:VAL:279: --log ep5_vdwV279_log.txt --bookmark_name ep5_vdwV279
```

We are now ready to export these molecules for visual inspection in your favorite molecular graphics program. We will create a new directory `ep5_vdwV279_sdfs` and store the exported molecule files there.

```
$ mkdir ep5_vdwV279_sdfs
$ rt_process_vs.py read --input_db all_groups.db --bookmark_name ep5_vdwV279 --export_sdf_path ep5_vdwV279_sdfs
```

Now we have our filtered molecules as SDF files ready for visual inspection!

## Example Filtering Timings (M1Pro MacBook, ~2 million ligands)
![rt_v11_timings](https://github.com/forlilab/Ringtail/assets/41704502/eac373fc-1324-45df-b845-6697dc9d1465)


# Extended documentation

## Scripts
The Ringtail package includes two command line oriented scripts: `rt_process_vs.py` and `rt_compare.py`. Both may be run with options specified in the command line and/or using options specified in a JSON-formatted file given with `--config`. Command line options override any conflicting options in the config file.

[rt_process_vs.py](https://github.com/forlilab/Ringtail#rt_process_vspy-documentation) serves as the primary script for the package and is used to both write docking files to a SQLite database and to perform filtering and export tasks on the database. It is designed to handle docking output files associated with a single virtual screening in a single database.

[rt_compare.py](https://github.com/forlilab/Ringtail#rt_comparepy-documentation) is used to combine information across multiple virtual screenings (in separate databases) to allow or exclude the selection of ligands passing filters across multiple targets/models. This can be useful for filtering out promiscuous ligands, a technique commonly used in exerimental high-throughput screening. It may also be used if selection of ligands binding multiple protein structures/conformations/homologs are desired.

## rt_process_vs.py Documentation
### Usage examples
#### Access help message for rt_process_vs.py
```
rt_process_vs.py --help
```
#### Access help message for rt_process_vs.py write mode
```
rt_process_vs.py write --help
```
#### Access help message for rt_process_vs.py read mode
```
rt_process_vs.py read --help
```
#### Create database named example.db from all input options
```
rt_process_vs.py write --file lig1.dlg lig2.dlg --file_path path1/ path2 --file_list filelist1.txt filelist2.txt --output_db example.db

```
Example file list
```
lig3.dlg
lig4.dlg.gz
rec1.pdbqt
```
#### Write and filter using a config file
```
rt_process_vs.py -c config_w.json write
rt_process_vs.py -c config_r.json read
```
config_w.json:

```
{
"file_path": "path1/",
"output_db": "example.db"
}
```

config_r.json:

```
{
"score_percentile": "0.1"
}
```

#### Export results from a previous filtering as a CSV
```
rt_process_vs.py write --file_path Files/
rt_process_vs.py read --input_db output.db --score_percentile 0.1 --bookmark_name filter1
rt_process_vs.py read --input_db output.db --export_bookmark_csv filter1
```
#### Create scatterplot highlighting ligands passing filters
```
rt_process_vs.py write --file_path Files/
rt_process_vs.py read --input_db output.db --score_percentile 0.1 --bookmark_name filter1
rt_process_vs.py read --input_db output.db --bookmark_name filter1 --plot
```
`all_ligands_scatter.png`

![all_ligands_scatter](https://user-images.githubusercontent.com/41704502/215909808-2edc29e9-ebdb-4f0e-a87a-a1c293687b2e.png)

### Usage Details
The script for writing a database and filtering is `rt_process_vs.py`. __This is intended to be used for a set of DLGs/Vina PDBQTs pertaining to a single target and binding site. This may include multiple ligand libraries as long as the target and binding site is the same. Be cautious when adding results from multiple screening runs, since some target information is checked and some is not.__ One receptor PDBQT may also be saved to the database.

The rt_process_vs.py script has two modes: `write` and `read`. The desired mode must be specified in the command line before any other options are given (except `-c [CONFIG]` which is given first). The `write` mode is used to create a database for a virtual screening from ADGPU DLGs or Vina PDBQTs. After this initial run, a database is created and may be read directly by rt_process_vs.py in `read` mode for subsequent filtering and export operations.

#### Write Mode
Upon calling rt_process_vs.py in `write` mode for the first time, the user must specify where the program can find files to write to the newly-created database. This is done using the
`--file`, `--file_path`, and/or `--file_list` options. Any combination of these options can be used, and multiple arguments for each are accepted. Compressed `.gz` files
are also accepted.

When searching for result files in the directory specified with `--file_path`, rt_process_vs.py will search for files with the pattern `*.dlg*` by default. This may be changed with the
`--pattern` option. Note also that, by default, Ringtail will only search the directory provided in `--file_path` and not subdirectories. Subdirectory searching
is enabled with the `--recursive` flag. If you are trying to read Vina PDBQTs, specify this with `--mode vina`. This will automatically change the file search pattern to `*.pdbqt*`. If the receptor PDBQT file is present in a directory being searched, it **must** be specified with `--receptor_file`.

To add new files to an existing database, the `--append_results` flag can be used in conjuction with `--input_db` and `--file`, `--file_path`, and/or `--file_list` options. If one is concerned about adding duplicate results, the `--duplicate_handling` option can be used to specify how duplicate entries should be handled. However, this option makes database writing significantly slower.

To overwrite an existing database, use the `--overwrite` flag.

One receptor PDBQT, corresponding to that in the DLGs, may be saved to the database using the `--save_receptor` flag. This will store the receptor file itself in a binary format in the database. The user must specify the path to the receptor file with the `--receptor_file` option. Ringtail will also throw an exception if this flag is given but no receptor is found, if the name of the receptor in any DLG does not match the receptor file, or if this flag is used with a database that already has a receptor. `--save_receptor` can be used to add a receptor to an existing database given with `--input_db`. `--save_receptor` may not be used with the `--append_results` option.

By default, the newly-created database will be named `output.db`. This name may be changed with the `--output_db` option.

By default (for DLGs), Ringtail will store the best-scored (lowest energy) binding pose from the first 3 pose clusters in the DLG. For Vina, Ringtail will store the 3 best poses. The number of clusters/poses stored may be
changed with the `--max_poses` option. The `--store_all_poses` flag may also be used to override `--max_poses` and store every pose from every file.

ADGPU is capable of performing interaction analysis at runtime, with these results being stored in the database if present. If interaction analysis is not present in the input file (including Vina PDBQTs), it may be added by Ringtail with the `--add_interactions` option. **This adds a signifcant increase to the total database write time.** Distance cutoffs for the interactions are specified with the `--interaction_cutoffs` option. Adding interactions requires that the receptor PDBQT be provided as an input by the user with the `--receptor_file` option.

The `--interaction_tolerance` option also allows the user to give more leeway for poses to pass given interaction filters. With this option, the interactions from poses within *c* angstrom RMSD of a cluster's top pose will be appended to the interactions for that top pose. The theory behind this is that this gives some sense of the "fuzziness" of a given binding pose, allowing the user to filter for interactions that may not be present for the top pose specifically, but could be easily accessible to it. When used as a flag, the `interaction_tolerance` default is 0.8 angstroms. The user may also specify their own cutoff. This option is intended for use with DLGs from AD-GPU, which clusters output poses based on RMSD.

#### Read mode
In `read` mode, an existing database is used to filter or export results.

When filtering, a text log file will be created containing the results passing the given filter(s). The default log name is `output_log.txt` and by default will include the ligand name and docking score of every pose passing filtering criteria. The log name
may be changed with the `--log` option and the information written to the log can be specified with `--outfields`. The full list of available output fields may be seen by using the `--help` option with `read` mode (see example above).
By default, only the information for the top-scoring binding pose will be written to the log. If desired, each individual passing pose can be written by using the `--output_all_poses` flag. The passing results may also be ordered in the log file using the `--order_results` option.

No filtering is performed if no filters are given. If both `--eworst` and `--score_percentile` are used together, the `--eworst` cutoff alone is used. The same is true of `--leworst` and `--le_percentile`.

In addition to the filtering options outlined in the table below, ligands passing given filters can be clustered to provide a reduced set of dissimilar ligands based on Morgan fingerprints (`--mfpt_cluster`) or interaction (`--interaction_cluster`) fingerprints. Dissimilarity is measured by Tanimoto distance and clustering is performed with the Butina clustering algorithm.

When filtering, the passing results are saved as a view in the database. This view is named `passing_results` by default. The user can specify a name for the view using the `--bookmark_name` option. Data for poses in a view may be accessed later using the `--data_from_bookmark` option. When `max_miss` > 0 is used, a view is created for each combination of interaction filters and is named `<bookmark_name>_<n>` where n is the index of the filter combination in the log file (indexing from 0).

Filtering may take from seconds to minutes, depending on the size of the database, roughly scaling as O(n) for n database Results rows (i.e. stored poses). One may also filter over a previous bookmark specified with the `--filter_bookmark` option. If using this option, the bookmarks specified by `--filter_bookmark` and `--bookmark_name` must be different.

While not quite a filtering option, the user can provide a ligand name from a previously-run clustering and re-output other ligands that were clustered with that query ligand with `--find_similar_ligands`. The user is prompted at runtime to choose a specific clustering group from which to re-output ligands. Filtering/clustering will be performed from the same command-line call prior to this similarity search, but all subsequent output tasks will be performed on the group of similar ligands obtained with this option unless otherwise specified. 

##### Other available outputs
The primary outputs from `rt_process_vs.py` are the database itself (`write` mode) and the filtering log file (`read` mode). There are several other output options as well, intended to allow the user to further explore the data from a virtual screening.

The `--plot` flag generates a scatterplot of ligand efficiency vs docking score for the top-scoring pose from each ligand. Ligands passing the given filters or in the bookmark given with `--bookmark_name` will be highlighted in red. The plot also includes histograms of the ligand efficiencies and binding energies. The plot is saved as `[filters_file].png` if a `--filters_file` is used, otherwise it is saved as `out.png`.

The `--pymol` flag also generates a scatterplot of ligand efficiency vs docking score, but only for the ligands contained in the bookmark specified with `--bookmark_name`. It also launches a PyMol session and will display the ligands in PyMol when clicked on the scatterplot. N.B.: Some users may encounter a `ConnectionRefusedError`. If this happens, try manually launching PyMol (`pymol -R`) in a separate terminal window.

Using the `--export_sdf_path` option allows the user to specify a directory to save SDF files for ligands passing the given filters or in the bookmark given with `--bookmark_name`. The SDF will contain poses passing the filter/in the bookmark ordered by increasing docking score. Each ligand is written to its own SDF. This option enables the visualization of docking results, and includes any flexible/covalent ligands from the docking. The binding energies, ligand efficiencies, and interactions are also written as properties within the SDF file, with the order corresponding to the order of the pose order.

If the user wishes to explore the data in CSV format, Ringtail provides two options for exporting CSVs. The first is `--export_bookmark_csv`, which takes a string for the name of a table or result bookmark in the database and returns the CSV of the data in that table. The file will be saved as `<table_name>.csv`.
The second option is `--export_query_csv`. This takes a string of a properly-formatted SQL query to run on the database, returning the results of that query as `query.csv`. This option allows the user full, unobstructed access to all data in the database.

As noted above, a bookmark may also be exported as a separate SQLite dabase with the `--export_bookmark_db` flag.

Finally, a receptor stored in the database may be re-exported as a PDBQT with the `--export_receptor` option.

### Interaction filter formatting and options

**Interaction filtering requires interactions to be present in database.**

The `--vdw`, `--hb`, and `--react_res` interaction filters must be specified in the order `CHAIN:RES:NUM:ATOM_NAME`. Any combination of that information may be used, as long as 3 colons are present and the information ordering between the colons is correct. All desired interactions of a given type (e.g. `--vdw`) may be specified with a single option tag (`--vdw=B:THR:276:,B:HIS:226:`) or separate tags (`--vdw=B:THR:276: --vdw=B:HIS:226:`).

The `--max_miss` option allows the user to filter by given interactions excluding up to `max_miss` interactions. This gives ![equation](https://latex.codecogs.com/svg.image?\sum_{m=0}^{m}\frac{n!}{(n-m)!*m!}) combinations for *n* interaction filters and *m* max_miss. By default, results will be given for the union of the interaction conbinations. Use with `--enumerate_interaction_combs` to log ligands/poses passing each separate interaction combination (can significantly increase runtime).

The `--smarts_idxyz` option may be used to filter for a specific ligand substructure (specified with a SMARTS string) to be placed within some distance of a given cartesian coordinate. The format for this option is `"<SMARTS pattern: str>" <index of atom in SMARTS: int> <cutoff distance: float> <target x coord: float> <target y coord: float> <target z coord: float>`.

### Exploring the database in the Command Line
View the data contained within the database using a terminal, we recommend using the [VisiData tool](https://www.visidata.org/). In addition to command line visualization, this tool has a number of other feature, including ploting. Upon opening the database with `vd`, the terminal should look like this:

![Screenshot from 2022-05-18 14-57-22](https://user-images.githubusercontent.com/41704502/169162632-3a71d338-faa1-4109-8f04-40a96ee6d24e.png)

In this example (made with DLGs), the database contains ~3 poses for 9999 discrete ligands. Each of the rows here is a separate table or view within the database. From this screen, you can easily perform the sanity checks outline below. One should note that the number of column displayed on the first screen is 1 greater than the actual number of columns in a table (the number is correct for views). To more fully explore a given table, one may use the arrow keys or mouse to navigate to it, then press `Enter/Return` to access that table/view. The user may then scroll horizontally with the arrow keys, or press `q` to return up a level.

Using `vd` is particularly helpful to examine possible interactions of interest, stored within the `Interaction_indices` table.

To exit, return to the screen shown in the image above by pressing `q`, then press `q` to exit.

### Data integrity sanity checks
There are a few quick checks the user can make to ensure that the data has been properly written from the input files to the database. Discrepancies may indicate an error occurred while writing the database or the input file format did not match that which Ringtail expected.
- The number of rows in the `Ligands` table should match the number of input ligand files
- The number of rows in the `Results` and `Interaction_bitvectors` tables should match
- Number of columns in the `Interactions_bitvectors` table should match the number of rows in the `Interaction_indices` table + 1 (+2 if using `vd`)
- The number of rows in the `Results` table should be ~`max_poses`\* `number of files` and should be less than or equal to that number. For DLGs not every ligand may have up to `max_poses`, which is why the number of rows is typically smaller than `max_poses`\* `number of DLGs`.
- No ligand should have more than `max_poses` rows in the `Results` table.
- If storing all poses, the number of rows in the Results table should match the `number of ligands` * `number of output poses`.

### Potential pitfalls
Any PDBQT files specified through any of the input options in ADGPU mode will be read by `rt_process_vs.py` as receptor files, even if the files actually represent ligands. Therefore, ligand PDBQT files should not be present in any directories given with `--file_path`.

When writing from Vina PDBQTs, ensure there are no other PDBQTs (input or receptor) in directories specified with `--file_path` UNLESS the receptor PDBQT is specified with the `--receptor_file` option.

Occassionally, errors may occur during database reading/writing that corrupt the database. This may result in the database becoming locked. If this occurs it is recommended to delete the existing database and re-write it from scratch.

If trying to read a database created with Ringtail v1.0.0 with a newer version of Ringtail, you may encounter errors related to changes to the internal database structure. If you encounter this, run the follow commands (example of database named `output.db`:
```
$ sqlite3 output.db
> ALTER TABLE Results RENAME energies_binding TO docking_score;
> ALTER TABLE Bookmarks ADD COLUMN filters;
```
If you encounter further errors related to views/bookmarks, please contact the ForliLab.

### rt_process_vs.py supported arguments

| Argument          || Description                                           | Default value   | Requires interactions |
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
|--log              |-l| Name for log of filtered results                      | output_log.txt   ||
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
|--name             |-n| Search for specific ligand name. Multiple names joined by "OR". Multiple filters should be separated by commas | no default  ||
|--max_nr_atoms     |-mna| Specify maximum number of heavy atoms a ligand may have | no default  ||
|--smarts           || SMARTS pattern(s) for substructur matching | no default  ||
|--smarts_idxyz     || SMARTS pattern, index of atom in SMARTS, cutoff distance, and target xyz coordinates. Finds poses in which the specified substructure atom is within the distance cutoff from the target location | no default  ||
|--smarts_join     |-n| logical operator for multiple SMARTS | OR  | <tr><td colspan="5">INTERACTION FILTERS</td></tr>
|--van_der_waals    |-vdw| Filter for van der Waals interaction with given receptor information.  | no default  | Yes|
|--hydrogen_bond    |-hb| Filter with hydrogen bonding interaction with given information. Does not distinguish between donating or accepting | no default  | Yes|
|--reactive_res     |-r| Filter for reation with residue containing specified information | no default  |Yes |
|--hb_count         |-hc| Filter for poses with at least this many hydrogen bonds. Does not distinguish between donating and accepting | no default  | Yes|
|--react_any        |-ra| Filter for poses with reaction with any residue       | FALSE     | Yes|
|--max_miss         |-mm| Will filter given interaction filters excluding up to max_miss interactions. Results in ![equation](https://latex.codecogs.com/svg.image?\sum_{m=0}^{m}\frac{n!}{(n-m)!*m!}) combinations for *n* interaction filters and *m* max_miss. Will log and output union of combinations unless used with `--enumerate_interaction_combs`. | 0  | Yes |
|--enumerate_interactions_combs  |-eic| When used with `--max_miss` > 0, will log ligands/poses passing each separate interaction filter combination as well as union of combinations. Can significantly increase runtime. | FALSE  | Yes <tr><td colspan="5">PASSING RESULT CLUSTERING</td></tr>
|--mfpt_cluster     |-mfpc| Cluster ligands passing given filters based on the Tanimoto distances of the Morgan fingerprints. Will output ligand with best (lowest) ligand efficiency from each cluster. Uses Butina clustering algorithm | 0.5  ||
|--interaction_cluster     |-ifpc| Cluster ligands passing given filters based on the Tanimoto distances of the interaction fingerprints. Will output ligand with best (lowest) ligand efficiency from each cluster. Uses Butina clustering algorithm | 0.5  | Yes |

---

## rt_compare.py Documentation
The `rt_compare.py` script is designed to be used with databases already made and filtered with the `rt_process_vs.py` script. The script is used to select ligands which are shared between the given filter bookmark(s) of some virtual screenings (wanted) or exclusive to some screenings and not others (unwanted). The basic process of preparing to use this script and the concept behind it is thus:

Let us assume that kinase1 is our target of interest. It has related proteins kinase1a and kinase1b. protein2 is an unrelated protein.
1. Create a database for each virtual screening on each target (kinase1.db, kinase1a.db, kinase1b.db, protein2.db)
2. Filter each database separately to get a set of virtual hits for each target. Each set of filters may be different as desired (e.g. change interaction filters for analogous residues). The bookmark within each database may be given as a single string (same bookmark name in every database) or multiple bookmark names (one per database) with the `--bookmark_name` option. If specifying multiple names, the order should match the order that the databases were provided in, beginning with wanted, then unwanted databases. The default name is `passing_results`.
3. Use `rt_compare.py` to find ligands that pass the filters for kinase1 but not kinase1a or kinase1b. This will create a log file of the same format as that output from `rt_process_vs.py`.
```
rt_compare.py --wanted kinase1.db --unwanted kinase1a.db kinase1b.db
```
4. Other usage examples and output options given below. For example, one can also select for potential dual-target ligands with
```
rt_compare.py --wanted kinase1.db protein2.db --unwanted kinase1a.db kinase1b.db
```

### Usage examples
#### Access help message for rt_compare.py
```
rt_compare.py --help
```
#### Select ligands found in "passing_results" bookmarks of vs1 but not vs2 or vs3
```
rt_compare.py --wanted vs1.db --unwanted vs2.db vs3.db
```
#### Select ligands found in "passing_results" bookmarks of vs1 and vs2 but not vs3 or vs4
```
rt_compare.py --wanted vs1.db vs2.db --unwanted vs3.db vs4.db
```
#### Select ligands found in "passing_results" bookmarks of every vs except vs4
```
rt_compare.py --wanted vs1.db vs2.db vs3.db --unwanted vs4.db
```
#### Select ligands found in "filter1" bookmarks of vs1 but not vs2
```
rt_compare.py --wanted vs1.db --unwanted vs2.db --bookmark_name filter1
```
#### Save bookmark of ligands found in "filter1" bookmarks of vs1 and vs2 but not vs3 or vs4 as "selective_bookmark" in vs1.db
```
rt_compare.py --wanted vs1.db vs2.db --unwanted vs3.db vs4.db --save_bookmark selective_bookmark
```
#### Export bookmark set of ligands found in "filter1" bookmarks of vs1 and vs2 but not vs3 or vs4 as CSV
```
rt_compare.py --wanted vs1.db vs2.db --unwanted vs3.db vs4.db --export_csv
```
### rt_compare.py supported arguments

| Argument          || Description                                           | Default value   |
|:------------------------|:-----|:-------------------------------------------------|----:|
|--config           | -c| Configuration JSON file to specify new default options. Overridded by command line | no default <tr><td colspan="4"></td></tr>
|--wanted |-w| Database files for which to include the intersection of ligands in bookmark_name(s) for all databases specified with this option.| no default|
|--unwanted |-n| Database files for which to exclude any ligands found in bookmark_name of any of the databases specified with this option. | no default|
|--bookmark_name |-sn| Name of bookmark to select ligands within. Must be present in all databases given.| passing_results|
|--log |-l| Name for log file| selective_log.txt |
|--save_bookmark| -s| Save the final selective bookmark as a view with given name in the first database specified with --wanted. | no default|
|--export_csv| -x| Save final selective bookmark as csv. Saved as [save_bookmark].csv or 'crossref.csv' if --save_bookmark not used.| FALSE|

---
## Brief python tutorials
#### Make sqlite database from current directory
```python
from ringtail import RingtailCore

opts = RingtailCore.get_defaults()
RingtailCore.set_opts(opts, ["storage_type", "db_file", "path", "recursive"], ["sqlite", "example.db", [["."]], True])

with RingtailCore(**opts) as rt_core:
    rt_core.add_results()
```
#### Convert database tables to pandas dataframes
```python
from ringtail import StorageManagerSQLite

# make database manager with connection to SQLite file vs.db
with StorageManagerSQLite("vs.db") as dbman:

    # fetch entire Results table as pandas dataframe
    results_df = dbman.to_dataframe("Results")

    # fetch entire Ligands table as pandas dataframe
    ligands_df = dbman.to_dataframe("Ligands")

    # fetch entire Receptors table as pandas dataframe
    rec_df = dbman.to_dataframe("Receptors")

    # fetch entire Interaction Indices table as pandas dataframe
    interaction_idx_df = dbman.to_dataframe("Interaction_indices")

    # fetch entire Interaction bitvectors table as pandas dataframe
    interaction_bv_df = dbman.to_dataframe("Interaction_bitvectors")

```
#### Make an ROC plot and calculate its AUC for a virtual screening with a file containing a list of known binders
```python
#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
from sklearn.metrics import auc
import matplotlib.pyplot as plt
from ringtail import RingtailCore, StorageManagerSQLite

# make database manager with connection to SQLite file output.db
with StorageManagerSQLite("output.db") as dbman:

    # fetch entire Results table as pandas dataframe
    results_df = dbman.to_dataframe("Select LigName FROM Results GROUP BY LigName ORDER BY energies_binding", table=False)

print(results_df)

with open("binders.txt", 'r') as f:
    binders = [l.strip() for l in f.readlines() if l != ""]

num_actives = len(binders)
num_decoys = len(results_df) - num_actives

roc_tpr = []
roc_fpr = []
tp = 0
fp = 0
i = 1
for l in results_df["LigName"]:
    print("Testing cutoff:", i)
    i += 1
    if l in binders:
        tp += 1
    else:
        fp += 1
    roc_tpr.append(tp / num_actives)
    roc_fpr.append(fp / num_decoys)
auc = auc(roc_fpr, roc_tpr)

with open("auc.txt", 'w') as f:
    f.write(str(auc))

plt.plot(roc_fpr, roc_tpr)
x = np.linspace(0,1)
plt.plot(x, x, linestyle="--")
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.show()
```

