# API procedures

### Command line tools for Ringtail 
Should work just as before, but please test them too to make sure I have not made any bugs in my re-write process. 
The one option that is not completely tested yet, is using the config file to add options. This method cannot be trusted to work at the moment. 

### Creating a new ringtail core

A ringtail core is created by instantiating the object with a database. Currently, a database can only be added upon instantiation.
```
rtc = RingtailCore("output.db")
```

Default logging level is "WARNING", and a different logger level can be set at the time of object instantiation.
```
rtc = RingtailCore(db_file="output.db", logging_level="DEBUG)
```

The default docking mode is "dlg", and can be changed to "vina" by accessing the ringtail core attribute `docking_mode`:
```
rtc.docking_mode = "vina"
```

### Adding results
To add results to the database, use the `add_results_from_files` method that takes as input files and directories to upload,
as well as a receptor path and database properties and how to handle the resutls (how many poses to save, how to deal with interactions if having vina results),
and whether or not to print a summary after writing the results to the database.

```
rtc.add_results_from_files( file_path = "test_data/", 
                            recursive = True, 
                            save_receptor = False,
                            max_poses = 3,
                            summary = True)
```

Both files (`filesources_dict`) and processing options (`optionsdict`) can be provided as dictionaries as well or instead of the the individual options. Any provided individual options will overwrite the options provided through dictionaries. 

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
                            optionsdict = writeoptions,
                            summary = True)
```

##### Available options for writing to the database include:

| File options        || Description                                           | Default value   | Requires interactions |
|:------------------------|:-----|:-------------------------------------------------|:----------------|----:|
|--file             |-f| DLG/Vina PDBQT file(s) to be read into database                  | no default       ||
|--file_path        |-fp| Path(s) to files to read into database            | no default       ||
|--file_list        |-fl| File(s) with list of files to read into database  | no default       ||
|--pattern          |-p| Specify pattern to search for when finding files   | \*.dlg\* / \*.pdbqt\* (vina mode)        ||
|--recursive        |-r| Flag to perform recursive subdirectory search on --file_path directory(s)  | FALSE      ||
|--receptor_file |-rn| Use with --save_receptor and/or --add_interactions. Give receptor PDBQT. | None      ||
|--save_receptor    |-sr| Flag to specify that receptor file should be imported to database. Receptor file must also be specified with --receptor_file| FALSE   ||
| Result processing         || Description                                           | Default value   | Requires interactions |
|:------------------------|:-----|:-------------------------------------------------|:----------------|----:|
|--max_poses        |-mp| Number of clusters for which to store top-scoring pose (dlg) or number of poses (vina) to save in database| 3     ||
|--store_all_poses  |-sa| Flag to indicate that all poses should be stored in database| FALSE      ||
|--interaction_tolerance|-it| Adds the interactions for poses within some tolerance RMSD range of the top pose in a cluster to that top pose. Can use as flag with default tolerance of 0.8, or give other value as desired | FALSE -> 0.8 (Å)  | Yes |
|--add_interactions  |-ai| Find interactions between ligands and receptor. Requires receptor PDBQT to be written. | FALSE      ||
|--interaction_cutoffs  |-ic| Specify distance cutoffs for measuring interactions between ligand and receptor in angstroms. Give as string, separating cutoffs for hydrogen bonds and VDW with comma (in that order). E.g. '-ic 3.7,4.0' will set the cutoff for hydrogen bonds to 3.7 angstroms and for VDW to 4.0. | 3.7,4.0     ||
|--max_proc |-mpr| Maximum number of subprocesses to spawn during database writing. | [# available CPUs]      |<tr><td colspan="5">**Read Mode**</td></tr>
|--append_results      |-a| Add new docking files to existing database given with --input_db  | FALSE       ||
|--duplicate_handling|-dh| Specify how dulicate results should be handled. May specify "ignore" or "replace". Unique results determined from ligand and target names and ligand pose. *NB: use of duplicate handling causes increase in database writing time*| None |
|--overwrite        |-ov| Flag to overwrite existing database           | FALSE       ||


If at any point you wish to print a summary of the database, the method can be called directly:
```
rtc.print_summary()
```

### Filtering on a database 
To filter results in a database the method `filter` is used. 

```
rtc.filter( eworst=-2, 
            vdw_interactions=[('A:VAL:279:', True), ('A:LYS:162:', True)])
```

```
filters = {
    "ebest":-5.5,
    "hb_interactions": [('A:VAL:279:', True), ('A:LYS:162:', True)],
    "max_miss": 1
}

rtc.filter( filters_dict = filters)
```

A "log file" containing details from the filtering exercise will be saved in the working directory to the default filename `output_log.txt`.
If you would like to change the name of this file, this can be done by running the method `set_read_options`:
```
rtc.set_read_options(log_file = "custom_logfile_name.txt")
```




#TODO the available options are: 

#TODO what to do with the outputs

### Output options

### Options
#TODO other methods
### Using the config file



| Argument          || Description                                           | Default value   | Requires interactions |
|:------------------------|:-----|:-------------------------------------------------|:----------------|----:|
|--config           | -c| Configuration JSON file to specify new default options. Overridded by command line | no default       |<tr><td colspan="5"></td></tr>
|--bookmark_name      |-s| Name for bookmark view in database                      | passing_results  ||
|--mode          |-m| specify AutoDock program used to generate results. Available options are "dlg" and "vina". Vina mode will automatically change --pattern to \*.pdbqt   | dlg         ||
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