![ringtail logo with text](https://user-images.githubusercontent.com/41704502/168933741-8b4a78f7-6e23-4d74-9c4c-4b47633e42c7.png)
# Ringtail
Package for creating SQLite database from virtual screening DLGs and performing filtering on results.

Ringtail reads collections of Docking Log File (DLG) results from virtual screenings performed with [AutoDock-GPU](https://github.com/ccsb-scripps/AutoDock-GPU) and deposits them into
a SQLite database. It then allows for the filtering of results with numerous pre-defined filtering options, generation of simple result plots, export of resulting
molecule poses, and export of CSVs of result data. DLG parsing is parallelized across the user's CPU.

Ringtail is developed by the [Forli lab](https://forlilab.org/) at the
[Center for Computational Structural Biology (CCSB)](https://ccsb.scripps.edu)
at [Scripps Research](https://www.scripps.edu/).

## Basic Usage
Upon calling Ringtail for the first time, the user must specify where the program can find DLGs to write to the newly-created database. This is done using the
`--file`, `--file_path`, and/or `--file_list` options. Any combination of these options can be used, and multiple arguments for each are accepted. dlg.gz files
are also accepted. After this initial run, a database is created and may be read directory for subsequent filtering operations.
#### Inputs
When searching for DLG files in the directory specified with `--file_path`, Ringtail will search for files with the pattern `*.dlg*`. This may be changed with the
`--pattern` option. Note also that, by default, Ringtail will only search the directory provided in `--file_path` and not subdirectories. Subdirectory searching
is enabled with the `--recursive` flag.

Once a database is written, this database can be read in directly without re-writting using the `--input_db` option. To add new DLGs to an existing database, the `--add_results` flag can be used in conjuction with `--input_db` and `--file`, `--file_path`, and/or `--file_list` options. To overwrite an existing database, use the `--overwrite` flag in combination with `--file`, `--file_path`, and/or `--file_list` options.
#### Outputs
By default, the newly-created database will be named `output.db`. This name may be changed with the `--output_db` option.

By default, Ringtail will store the best-scored (lowest energy) binding pose from the first 3 pose clusters in the DLG. The number of clusters stored may be
changed with the `--max_poses` option. The `--store_all_poses` flag may also be used to override `--max_poses` and store every pose from every DLG.

The default log name is `output_log.txt` and by default will include the ligand name and binding energy of every pose passing filtering criteria. The log name
may be changed with the `--log` option and the information written to the log can be specified with `--out_fields`. By default, only the information for the
top-scoring binding pose will be written to the log. If desire, each individual passing pose can be written by using the `--all_poses` flag.
#### Filters
When running with default settings (no user-specified filters), the only filter used is `--epercentile 1.0`. This gives the top 1% of poses by overall binding energy score.

## Interaction filter formatting
The `--vdw`, `--hb`, and `--react_res` interaction filters must be specified in the order `CHAIN:RES:NUM:ATOM_NAME`. Any combination of that information may be used, as long as 3 colons are present and the information ordering between the colons is correct. All desired interactions of a given type (e.g. `--vdw`) may be specified with a single option tag (`--vdw=B:THR:276:,B:HIS:226:`) or separate tags (`--vdw=B:THR:276: --vdw=B:HIS:226:`).

## Using filters file
Filters may also be read from a text file given with the `--filters_file` tag. Below is an example of a filters text file:
```
#     this is a comment
eworst=-0.4
ebest=-100
leworst=-0.2
#    vdw=A:THR:276:
#    vdw=B:HIS:36:
react_res=A:LYS:307:
react_count=1
hb_count=10
max_miss=1
```

## Potential pitfalls
When importing DLG files into a database with Ringtail, the files must have interaction analysis already performed. This is specified with `-C 1` when running [AutoDock-GPU](https://github.com/ccsb-scripps/AutoDock-GPU).

## Supported arguments

| Argument          | Description                                           | Default value    <tr><td colspan="3">**INPUT**</td></tr>
|:---------------------|:---------------------------------------------------|-----------------:|
|--file             | DLG file(s) to be read into database                  | no default       |
|--file_path        | Path(s) to DLG files to read into database            | no default       |
|--file_list        | File(s) with list of DLG files to read into database  | no default       |
|--save_receptors   | Flag to specify that receptor files should be imported to database. Receptor files must also be in locations specified with --file, --file_path, and/or --file_list| FALSE   |
|--recursive        | Flag to perform recursive subdirectory search on --file_path directory(s)  | FALSE      |
|--pattern          | Specify patter to serach for when finding DLG files   | \*.dlg\*         |
|--filters_file     | Text file specifying filters. Override command line filters  | no default|
|--input_db         | Database file to use instead of creating new database | no default       |
|--add_results      | Add new DLG files to existing database given with --input_db  | FALSE       |
|--conflict_handling| Specify how conflicting results should be handled. May specify "ignore" or "replace". Unique results determined from ligand and target names and ligand pose. *NB: use of conflict handling causes increase in database writing time*| None      |
|--one_receptor     | Flag to indicate that all results being added share the same receptor. Decreased runtime when using --save_receptors  | FALSE <tr><td colspan="3">**OUTPUT**</td></tr>
|--output_db        | Name for output database                              | output.db        |
|--export_table_csv | Name of database table to be exported as CSV. Output as <table_name>.csv | no default      |
|--export_query_csv | Create csv of the requested SQL query. Output as query.csv. MUST BE PRE-FORMATTED IN SQL SYNTAX e.g. SELECT [columns] FROM [table] WHERE [conditions] | no default      |
|--interaction_tolerance | Adds the interactions for poses within some tolerance RMSD range of the top pose in a cluster to that top pose. Can use as flag with default tolerance of 0.8, or give other value as desired | FALSE -> 0.8 (Å)  |
|--max_poses        | Number of cluster for which to store top-scoring pose in database| 3     |
|--store_all_poses  | Flag to indicate that all poses should be stored in database| FALSE      |
|--log              | Name for log of filtered results                      | output_log.txt   |
|--subset_name      | Name for subset view in database                      | passing_results  |
|--out_fields       | Data fields to be written in output (log file and STDOUT). Ligand name always included. | e        |
|--data_from_subset | Flag that out_fields data should be written to log for results in given --subset_name. Requires --no_filter. | FALSE       |
|--export_poses_path| Path for saving exported SDF files of ligand poses passing given filtering criteria | no default       |
|--no_print         | Flag indicating that passing results should not be printed to STDOUT | FALSE        |
|--no_filter        | Flag that no filtering should be performed            | FALSE       |
|--plot             | Flag to create scatterplot of ligand efficiency vs binding energy for best pose of each ligand. Saves as [filters_file].png or out.png. | FALSE        |
|--all_poses        | Flag that if mutiple poses for same ligand pass filters, log all poses | (OFF)        |
|--overwrite        | Flag to overwrite existing log and database           | FALSE       |
|--order_results    | String for field by which the passing results should be ordered in log file. | no default  <tr><td colspan="3">**PROPERTY FILTERS**</td></tr>
|--eworst           | Worst energy value accepted (kcal/mol)                | no_default  |
|--ebest            | Best energy value accepted (kcal/mol)                 | no default  |
|--leworst          | Worst ligand efficiency value accepted                | no default  |
|--lebest           | Best ligand efficiency value accepted                 | no default  |
|--epercentile      | Worst energy percentile accepted. Give as percentage (1 for top 1%, 0.1 for top 0.1%) | 1.0  |
|--leffpercentile   | Worst ligand efficiency percentile accepted. Give as percentage (1 for top 1%, 0.1 for top 0.1%) | no default  <tr><td colspan="3">**LIGAND FILTERS**</td></tr>
|--name             | Search for specific ligand name. Joined by "OR" with substructure search | no default  |
|--substructure     | SMILES substring to search for. *Performs substring search, will not find equivalent chemical structures with different SMILES denotations.* | no default  |
|--substruct_join   | Specify whether to join separate substructure searchs with "AND" or "OR". | "OR"  <tr><td colspan="3">**INTERACTION FILTERS**</td></tr>
|--vdw              | Filter for van der Waals interaction with given receptor information.  | no default  |
|--hb               | Filter with hydrogen bonding interaction with given information. Does not distinguish between donating or accepting | no default  |
|--max_miss         | Will separately filter each combination of given interaction filters excluding up to max_miss interactions. Results in ![equation](https://latex.codecogs.com/svg.image?\sum_{m=0}^{m}\frac{n!}{(n-m)!*m!}) combinations for *n* interaction filters and *m* max_miss. Results for each combination written separately in log file. Cannot be used with --plot or --export_poses_path. | 0  |
|--hb_count          | Filter for poses with at least this many hydrogen bonds. Does not distinguish between donating and accepting | no default  |
|--react_any         | Filter for poses with reaction with any residue       | FALSE     |
|--react_res         | Filter for reation with residue containing specified information | no default  |
