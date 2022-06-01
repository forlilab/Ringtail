![ringtail logo final](https://user-images.githubusercontent.com/41704502/170797800-53a9d94a-932e-4936-9bea-e2d292b0c62b.png)



# Ringtail
Package for creating SQLite database from virtual screening results and performing filtering on results. Compatible with [AutoDock-GPU](https://github.com/ccsb-scripps/AutoDock-GPU) and [AutoDock-Vina](https://github.com/ccsb-scripps/AutoDock-Vina).

[![AD compat](https://img.shields.io/badge/AutoDock_Compatibility-ADGPU|Vina-brightgreen)](https://shields.io/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![API stability](https://img.shields.io/badge/stable%20API-no-orange)](https://shields.io/)
[![made-with-python](https://img.shields.io/badge/Made%20with-Python-1f425f.svg)](https://www.python.org/)

Ringtail reads collections of Docking Log File (DLG) or PDBQT results from virtual screenings performed with [AutoDock-GPU](https://github.com/ccsb-scripps/AutoDock-GPU) and [AutoDock-Vina](https://github.com/ccsb-scripps/AutoDock-Vina), respectively, and deposits them into
a SQLite database. It then allows for the filtering of results with numerous pre-defined filtering options, generation of simple result plots, export of 
molecule SDFs, and export of CSVs of result data. Result file parsing is parallelized across the user's CPU.

Ringtail is developed by the [Forli lab](https://forlilab.org/) at the
[Center for Computational Structural Biology (CCSB)](https://ccsb.scripps.edu)
at [Scripps Research](https://www.scripps.edu/).

## Installation
It is recommended that you create a new Conda environment for installing Ringtail. Ringtail requires the following:
- RDKit
- SciPy
- [Meeko](https://github.com/forlilab/Meeko) (from the Forli Lab)

Installation is outlined below:
```
conda create -n ringtail
conda activate ringtail
conda install -c conda-forge rdkit scipy
```
Now, navigate to the desired directory for installing Meeko and do the following:
```
$ git clone git@github.com:forlilab/Meeko.git
$ cd Meeko
$ pip install .
```
After this, navigate to the desired directory for installing Ringtail and do the following:
```
$ git clone git@github.com:forlilab/Ringtail.git
$ cd Ringtail
$ pip install .
```
If you wish to make the code for either Meeko or Ringtail editable without re-running `pip install .`, instead use
```
pip install --editable .
```
## Definitions
- __DLG__: Docking Log File, output from AutoDock-GPU.
- __PDBQT__: Modified PDB format, used for receptors (input to AutoDock-GPU and Vina) and output ligand poses from AutoDock-Vina.
- __Cluster__: Each DLG contains a number of independent runs, usually 20-50. These independent poses are then clustered by RMSD, giving groups of similar poses called clusters.
- __Pose__: The predicted ligand shape and position for single run of a single ligand in a single receptor.
- __Binding score/ binding energy__: The predicited binding energy from AutoDock.
- __Ringtail__: 
> Drat, I'm not a cat!  Even though this eye-catching omnivore sports a few vaguely feline characteristics such as pointy ears, a sleek body, and a fluffy tail, the ringtail is really a member of the raccoon family. https://animals.sandiegozoo.org/animals/ringtail


## Usage examples
#### Access help message for run_ringtail.py
```
run_ringtail.py --help
```
#### Access help message for run_ringtail.py write mode
```
run_ringtail.py write --help
```
#### Access help message for run_ringtail.py read mode
```
run_ringtail.py read --help
```
#### Create database named example.db from all input options
```
run_ringtail.py write --file lig1.dlg lig2.dlg --file_path path1/ path2 --file_list filelist1.txt filelist2.txt --output_db example.db"

```
Example file list
```
lig3.dlg
lig4.dlg.gz
rec1.pdbqt
```
#### Write and filter using a config file
```
run_ringtail.py -c config_w.json write
run_ringtail.py -c config_r.json read
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
"energy_percentile": "0.1"
}
```

#### Export results from a previous filtering as a CSV
```
run_ringtail.py write --file_path Files/
run_ringtail.py read --input_db output.db --epercentile 0.1 --subset_name filter1
run_ringtail.py read --input_db output.db --export_table_csv filter1
```
#### Create scatterplot highlighting ligands passing filters
```
run_ringtail.py write --file_path Files/
run_ringtail.py read --input_db output.db --epercentile 0.1 --subset_name filter1
run_ringtail.py read --input_db output.db --subset_name filter1 --plot
```
`all_ligands_scatter.png`

![all_ligands_scatter](https://user-images.githubusercontent.com/41704502/171295726-7315f929-edfa-49a0-b984-dacadf1a4327.png)

## Usage Details
The script for writing a database and filtering is `run_ringtail.py`. __This is intended to be used for a set of DLGs/Vina PDBQTs pertaining to a single target. This may include multiple ligand libraries as long as the target is the same. Be cautious when adding results from multiple screening runs, since some target information is checked and some is not.__ One receptor PDBQT may also be included if using with DLGs. 

The run_ringtail.py script has two modes: `write` and `read`. The desired mode must be specified in the command line before any other options are given. The `write` mode is used to create a database for a virtual screening from ADGPU DLGs or Vina PDBQTs. After this initial run, a database is created and may be read directly by run_ringtail.py in `read` mode for subsequent filtering operations.

#### Write Mode
Upon calling run_ringtail.py in `write` mode for the first time, the user must specify where the program can find files to write to the newly-created database. This is done using the
`--file`, `--file_path`, and/or `--file_list` options. Any combination of these options can be used, and multiple arguments for each are accepted. Compressed `.gz` files
are also accepted.

When searching for result files in the directory specified with `--file_path`, run_ringtail.py will search for files with the pattern `*.dlg*` by default. This may be changed with the
`--pattern` option. Note also that, by default, Ringtail will only search the directory provided in `--file_path` and not subdirectories. Subdirectory searching
is enabled with the `--recursive` flag.

To add new files to an existing database, the `--add_results` flag can be used in conjuction with `--input_db` and `--file`, `--file_path`, and/or `--file_list` options. If one is concerned about adding duplicate results, the `--duplicate_handling` option can be used to specify how duplicate entries should be handled. However, this option makes database writing significantly slower.

To overwrite an existing database, use the `--overwrite` flag.

One receptor PDBQT, corresponding to that in the DLGs, may be saved to the database using the `--save_receptor` flag. This will store the receptor file itself in a binary format in the database and is only enabled for use with DLGs. Ringtail will throw an exception if this flag is given but no receptor is found, if the name of the receptor in any DLG does not match the receptor file, or if this flag is used with a database that already has a receptor. `--save_receptor` can be used to add a receptor to an existing database given with `--input_db`. `--save_receptor` may not be used with the `--add_results` option and will not be able to be used with a database created from Vina PDBQTs.

By default, the newly-created database will be named `output.db`. This name may be changed with the `--output_db` option.

By default (for DLGs), Ringtail will store the best-scored (lowest energy) binding pose from the first 3 pose clusters in the DLG. The number of clusters stored may be
changed with the `--max_poses` option. The `--store_all_poses` flag may also be used to override `--max_poses` and store every pose from every DLG. All poses from a Vina PDBQT are saved to the database.

The `--interaction_tolerance` option also allows the user to give more leeway for poses to pass given interaction filters. With this option, the interactions from poses within *c* angstrom RMSD of a cluster's top pose will be appended to the interactions for that top pose. The theory behind this is that this gives some sense of the "fuzziness" of a given binding pose, allowing the user to filter for interactions that may not be present for the top pose specifically, but could be easily accessible to it. When used as a flag, the `interaction_tolerance` default is 0.8 angstroms. The user may also specify their own cutoff.

#### Read mode
In `read` mode, an existing database is read in and used to filter or export results.

When filtering, a text log file will be created containing the results passing the given filter(s). The default log name is `output_log.txt` and by default will include the ligand name and binding energy of every pose passing filtering criteria. The log name
may be changed with the `--log` option and the information written to the log can be specified with `--out_fields`. The full list of available output fields may be seen by using the `--help`` option with `read` mode (see example above).
By default, only the information for the top-scoring binding pose will be written to the log. If desired, each individual passing pose can be written by using the `--all_poses` flag. The passing results may also be ordered in the log file using the `--order_results` option.

No filtering is performed if no filters are given. If both `--eworst` and `--epercentile` are used together, the `--eworst` cutoff alone is used. The same is true of `--leworst` and `--leffpercentile`.

When filtering, the passing results are saved as a view in the database. This view is named `passing_results` by default. The user can specify a name for the view using the `--subset_name` option. Other data for poses in a view may be accessed later using the `--new_data_from_subset` option. When `max_miss` > 0 is used, a view is created for each combination of interaction filters and is named `<subset_name>_<n>` where n is the index of the filter combination in the log file (indexing from 0).

##### Other available outputs
The primary outputs from run_ringtail are the database itself (`write` mode) and the filtering log file (`read` mode). There are several other output options as well, intended to allow the user to further explore the data from a virtual screening.

The `--plot` flag generates a scatterplot of ligand efficiency vs binding energy for the top-scoring pose from each ligand. Ligands passing the given filters or in the subset given with `--subset_name` will be highlighted in red. The plot also includes histograms of the ligand efficiencies and binding energies. The plot is saved as `[filters_file].png` if a `--filters_file` is used, otherwise it is saved as `out.png`.

Using the `--export_sdf_path` option allows the user to specify a directory to save SDF files for ligands passing the given filters or in the subset given with `--subset_name`. The SDF will contain poses passing the filter/in the subset ordered by increasing binding energy. Each ligand is written to its own SDF. This option enables the visualization of docking results, and includes any flexible/covalent ligands from the docking. The binding energies, ligand efficiencies, and interactions are also written as properties within the SDF file, with the order corresponding to the order of the pose order.

If the user wishes to explore the data in CSV format, Ringtail provides two options for exporting CSVs. The first is `--export_subset_csv`, which takes a string for the name of a table or result subset in the database and returns the CSV of the data in that table. The file will be saved as `<table_name>.csv`.
The second option is `--export_query_csv`. This takes a string of a properly-formatted SQL query to run on the database, returning the results of that query as `query.csv`. This option allows the user full, unobstructed access to all data in the database.

## Interaction filter formatting and options

**Interaction filtering is not available for databases created with Vina PDBQTs.**

The `--vdw`, `--hb`, and `--react_res` interaction filters must be specified in the order `CHAIN:RES:NUM:ATOM_NAME`. Any combination of that information may be used, as long as 3 colons are present and the information ordering between the colons is correct. All desired interactions of a given type (e.g. `--vdw`) may be specified with a single option tag (`--vdw=B:THR:276:,B:HIS:226:`) or separate tags (`--vdw=B:THR:276: --vdw=B:HIS:226:`).

The `--max_miss` option allows the user to separately filter each combination of the given interaction filters excluding up to `max_miss` interactions. This gives ![equation](https://latex.codecogs.com/svg.image?\sum_{m=0}^{m}\frac{n!}{(n-m)!*m!}) combinations for *n* interaction filters and *m* max_miss. Results for each combination of interaction filters will be written separately in the log file. This option cannot be used with `--plot` or `--export_poses_path`.

## Exploring the database in the Command Line
View the data contained within the database using a terminal, we recommend using the [VisiData tool](https://www.visidata.org/). In addition to command line visualization, this tool has a number of other feature, including ploting. Upon opening the database with `vd`, the terminal should look like this:

![Screenshot from 2022-05-18 14-57-22](https://user-images.githubusercontent.com/41704502/169162632-3a71d338-faa1-4109-8f04-40a96ee6d24e.png)

In this example (made with DLGs), the database contains ~3 poses for 9999 discrete ligands. Each of the rows here is a separate table or view within the database. From this screen, you can easily perform the sanity checks outline below. One should note that the number of column displayed on the first screen is 1 greater than the actual number of columns in a table (the number is correct for views). To more fully explore a given table, one may use the arrow keys or mouse to navigate to it, then press `Enter/Return` to access that table/view. The user may then scroll horizontally with the arrow keys, or press `q` to return up a level.

Using `vd` is particularly helpful to examine possible interactions of interest, stored within the `Interaction_indices` table.

To exit, return to the screen shown in the image above by pressing `q`, then press `q` to exit.

## Data integrity sanity checks
There are a few quick checks the user can make to ensure that the data has been properly written from the DLGs to the database. Discrepancies may indicate an error occurred while writting the database or the DLG format did not that which Ringtail expected.
- The number of rows in the `Ligands` table should match the number of input ligand files
- The number of rows in the `Results` and `Interaction_bitvectors` tables should match (DLGs only)
- Number of columns in the `Interactions_bitvectors` table should match the number of rows in the `Interaction_indices` table + 1 (+2 if using `vd`)
- The number of rows in the `Results` table should be ~`max_poses`\* `number of DLGs` and should be less than or equal to that number. Not every ligand may have up to `max_poses`, which is why the number of rows is typically smaller than `max_poses`\* `number of DLGs`. (DLGs only)
- No ligand should have more than `max_poses` rows in the `Results` table (unless storing results from multiple virtual screenings in the same database; DLGs only).
- If reading from Vina PDBQTs/using storing all poses, the number of rows in the Results table should match the `number of ligands` * `number of output poses`.

## Potential pitfalls
When importing DLG files into a database with Ringtail, the files must have interaction analysis already performed. This is specified with `-C 1` when running [AutoDock-GPU](https://github.com/ccsb-scripps/AutoDock-GPU).

Any PDBQT files specified through any of the input options will be read by `run_ringtail.py` as receptor files, even if the files actually represent ligands. Therefore, ligand PDBQT files should not be present in any directories given with `--file_path`.

Occassionally, errors may occur during database reading/writing that corrupt the database. If this occurs and you start running into unclear errors related to the SQLite3 package, it is recommended to delete the existing database and re-write it from scratch.

## run_ringtail.py supported arguments

| Argument          || Description                                           | Default value   | Vina-compatible |
|:------------------------|:-----|:-------------------------------------------------|:----------------|----:|
|--config           | -c| Configuration JSON file to specify new default options. Overridded by command line | no default       |<tr><td colspan="5"></td></tr>
|--input_db         | -i| Database file to use instead of creating new database | no default       ||
|--subset_name      |-s| Name for subset view in database                      | passing_results  ||
|--mode          |-m| specify AutoDock program used to generate results. Available options are "ADGPU" and "Vina". Vina mode will automatically change --pattern to \*.pdbqt   | ADGPU         ||
|--verbose          |-v| Flag indicating that passing results should be printed to STDOUT | FALSE        | <tr><td colspan="5">**Write Mode**</td></tr>
|--file             |-f| DLG/Vina PDBQT/receptor file(s) to be read into database                  | no default       ||
|--file_path        |-fp| Path(s) to files to read into database            | no default       ||
|--file_list        |-fl| File(s) with list of files to read into database  | no default       ||
|--pattern          |-p| Specify pattern to search for when finding DLG files   | \*.dlg\*         ||
|--recursive        |-r| Flag to perform recursive subdirectory search on --file_path directory(s)  | FALSE      ||
|--add_results      |-a| Add new DLG files to existing database given with --input_db  | FALSE       ||
|--duplicate_handling|-dh| Specify how dulicate results should be handled. May specify "ignore" or "replace". Unique results determined from ligand and target names and ligand pose. *NB: use of duplicate handling causes increase in database writing time*| None |
|--save_receptor    |-sr| Flag to specify that receptor file should be imported to database. Receptor file must also be in a location specified with --file, --file_path, or --file_list| FALSE   |No |
|--output_db        |-o| Name for output database                              | output.db        ||
|--overwrite        |-ov| Flag to overwrite existing log and database           | FALSE       ||
|--max_poses        |-mp| Number of cluster for which to store top-scoring pose in database| 3     |No|
|--store_all_poses  |-sa| Flag to indicate that all poses should be stored in database| FALSE      |No|
|--interaction_tolerance|-it| Adds the interactions for poses within some tolerance RMSD range of the top pose in a cluster to that top pose. Can use as flag with default tolerance of 0.8, or give other value as desired | FALSE -> 0.8 (Å)  | No <tr><td colspan="5">**Read Mode**</td></tr>
|--log              |-l| Name for log of filtered results                      | output_log.txt   ||
|--out_fields       |-of| Data fields to be written in output (log file and STDOUT). Ligand name always included. | e        ||
|--order_results    |-ord| String for field by which the passing results should be ordered in log file. | no default ||
|--all_poses        |-ap| Flag that if mutiple poses for same ligand pass filters, log all poses | (OFF)        ||
|--export_subset_csv |-xs| Name of database result subset or table to be exported as CSV. Output as <table_name>.csv | no default      ||
|--export_query_csv |-xq| Create csv of the requested SQL query. Output as query.csv. MUST BE PRE-FORMATTED IN SQL SYNTAX e.g. SELECT [columns] FROM [table] WHERE [conditions] | no default      ||
|--export_sdf_path|-sdf| Path for saving exported SDF files of ligand poses passing given filtering criteria | no default       |No|
|--new_data_from_subset |-nd| Flag that out_fields data should be written to log for results in given --subset_name. Requires no filters. | FALSE       ||
|--plot             |-p| Flag to create scatterplot of ligand efficiency vs binding energy for best pose of each ligand. Saves as [filters_file].png or out.png. | FALSE        | <tr><td colspan="5">PROPERTY FILTERS</td></tr>
|--eworst           |-e| Worst energy value accepted (kcal/mol)                | no_default  ||
|--ebest            |-eb| Best energy value accepted (kcal/mol)                 | no default  ||
|--leworst          |-le| Worst ligand efficiency value accepted                | no default  ||
|--lebest           |-leb| Best ligand efficiency value accepted                 | no default  ||
|--energy_percentile      |-pe| Worst energy percentile accepted. Give as percentage (1 for top 1%, 0.1 for top 0.1%) | 1.0  ||
|--le_percentile   |-ple| Worst ligand efficiency percentile accepted. Give as percentage (1 for top 1%, 0.1 for top 0.1%) | no default |  <tr><td colspan="5">LIGAND FILTERS</td></tr>
|--name             |-n| Search for specific ligand name. Joined by "OR" with substructure search and multiple names. Multiple filters should be separated by commas | no default  ||
|--substructure     |-st| SMILES substring to search for. *Performs substring search, will not find equivalent chemical structures with different SMILES denotations.* Multiple filters should be separated by commas | no default  | No|
|--substruct_join   |-sj| Specify whether to join separate substructure searchs with "AND" or "OR". | "OR" | No  <tr><td colspan="5">INTERACTION FILTERS</td></tr>
|--van_der_waals    |-vdw| Filter for van der Waals interaction with given receptor information.  | no default  | No|
|--hydrogen_bond    |-hb| Filter with hydrogen bonding interaction with given information. Does not distinguish between donating or accepting | no default  | No|
|--reactive_res     |-r| Filter for reation with residue containing specified information | no default  | No|
|--hb_count         |-hc| Filter for poses with at least this many hydrogen bonds. Does not distinguish between donating and accepting | no default  | No|
|--react_any        |-ra| Filter for poses with reaction with any residue       | FALSE     | No|
|--max_miss         |-mm| Will separately filter each combination of given interaction filters excluding up to max_miss interactions. Results in ![equation](https://latex.codecogs.com/svg.image?\sum_{m=0}^{m}\frac{n!}{(n-m)!*m!}) combinations for *n* interaction filters and *m* max_miss. Results for each combination written separately in log file. Cannot be used with --plot or --export_poses_path. | 0  | No|
