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

A receptor can be added to a database after the fact, using the method `save_receptor`:
```
rtc.save_receptor(receptor_file = "realistic_receptor_name.pdbqt")
```

If at any point you wish to print a summary of the database, the method can be called directly:
```
rtc.print_summary()
```

##### Available options for writing to the database include:

| File options        |Description                                           | Default value   | Requires interactions |
|:------------------------|:-------------------------------------------------|:----------------|----:|
|file             | DLG/Vina PDBQT file(s) to be read into database                  | no default       ||
|file_path        | Path(s) to files to read into database            | no default       ||
|file_list        | File(s) with list of files to read into database  | no default       ||
|pattern          | Specify pattern to search for when finding files   | \*.dlg\* / \*.pdbqt\* (vina mode)        ||
|recursive        | Flag to perform recursive subdirectory search on file_path directory(s)  | FALSE      ||
|receptor_file | Use with save_receptor and/or add_interactions. Give receptor PDBQT. | None      ||
|save_receptor    | Flag to specify that receptor file should be imported to database. Receptor file must also be specified with receptor_file| FALSE     |<tr><td colspan="5">***Result processing options***</td></tr>
|max_poses        | Number of clusters for which to store top-scoring pose (dlg) or number of poses (vina) to save in database| 3     ||
|store_all_poses  | Flag to indicate that all poses should be stored in database| FALSE      ||
|interaction_tolerance| Adds the interactions for poses within some tolerance RMSD range of the top pose in a cluster to that top pose. Can use as flag with default tolerance of 0.8, or give other value as desired | FALSE -> 0.8 (Å)  | Yes |
|add_interactions  | Find interactions between ligands and receptor. Requires receptor PDBQT to be written. | FALSE      ||
|interaction_cutoffs  | Specify distance cutoffs for measuring interactions between ligand and receptor in angstroms. Give as string, separating cutoffs for hydrogen bonds and VDW with comma (in that order). E.g. '3.7,4.0' will set the cutoff for hydrogen bonds to 3.7 angstroms and for VDW to 4.0. | 3.7,4.0     ||
|max_proc | Maximum number of subprocesses to spawn during database writing. | [# available CPUs]      
|append_results      | Add new docking files to existing database given with input_db  | FALSE       ||
|duplicate_handling| Specify how dulicate results should be handled. May specify "ignore" or "replace". Unique results determined from ligand and target names and ligand pose. *NB: use of duplicate handling causes increase in database writing time*| None |
|overwrite        | Flag to overwrite existing database           | FALSE       ||


### Filtering on a database 
To filter results in a database the method `filter` is called on the ringtail core. Filter values can be set directly in the method call:

```
num_ligands_passing_filters = rtc.filter( eworst=-2, 
                                        vdw_interactions=[('A:VAL:279:', True), ('A:LYS:162:', True)])
```

You can also create a dictionary of filters, and pass this to the `filter` method:
```
filters = {
    "ebest":-5.5,
    "hb_interactions": [('A:VAL:279:', True), ('A:LYS:162:', True)],
    "max_miss": 1
}

num_ligands_passing_filters = rtc.filter( filters_dict = filters)
```

Storage and read settings can also be set directly in the method call, for example:
```
num_ligands_passing_filters = rtc.filter( eworst=-2, 
                                        vdw_interactions=[('A:VAL:279:', True), ('A:LYS:162:', True)],
                                        log_file = "experiment424_log.txt",
                                        overwrite = False,
                                        output_all_poses = True)

```

Such settings can be set using a dictionary as well:
```
filters = {
    "ebest":-5.5,
    "hb_interactions": [('A:VAL:279:', True), ('A:LYS:162:', True)],
    "max_miss": 1
}

options = {
    "log_file": "experiment424_log.txt",
    "overwrite": False,
    "output_all_poses": True
}

num_ligands_passing_filters = rtc.filter( filters_dict = filters, options_dict = options)
```

Available filter and options are:

| Filters          || Description                                           | Default value   | Requires interactions |
|:------------------------|:-----|:-------------------------------------------------|:----------------|----:|
|eworst           | Worst energy value accepted (kcal/mol)                | no default  ||
|ebest            | Best energy value accepted (kcal/mol)                 | no default  ||
|leworst          | Worst ligand efficiency value accepted                | no default  ||
|lebest           | Best ligand efficiency value accepted                 | no default  ||
|score_percentile      | Worst energy percentile accepted. Give as percentage (1 for top 1%, 0.1 for top 0.1%) | 1.0  ||
|le_percentile    | Worst ligand efficiency percentile accepted. Give as percentage (1 for top 1%, 0.1 for top 0.1%) | no default |  <tr><td colspan="5">LIGAND FILTERS</td></tr>
|ligand_name             | Search for specific ligand name. Multiple names joined by "OR". Multiple filters should be separated by commas | no default  ||
|ligand_max_atoms     | Specify maximum number of heavy atoms a ligand may have | no default  ||
|ligand_substruct           | SMARTS pattern(s) for substructur matching | no default  ||
|ligand_substruct_pos     | SMARTS pattern, index of atom in SMARTS, cutoff distance, and target xyz coordinates. Finds poses in which the specified substructure atom is within the distance cutoff from the target location | no default  ||
|ligand_operator     | logical operator for multiple SMARTS | OR  | <tr><td colspan="5">INTERACTION FILTERS</td></tr>
|vdw_interactions    | Filter for van der Waals interaction with given receptor information.  | no default  | Yes|
|hb_interactions    | Filter with hydrogen bonding interaction with given information. Does not distinguish between donating or accepting | no default  | Yes|
|reactive_interactions     | Filter for reation with residue containing specified information | no default  |Yes |
|interactions_count         | Filter for poses with at least this many hydrogen bonds. Does not distinguish between donating and accepting | no default  | Yes|
|react_any        | Filter for poses with reaction with any residue       | FALSE     | Yes|
|max_miss         | Will filter given interaction filters excluding up to max_miss interactions. Results in ![equation](https://latex.codecogs.com/svg.image?\sum_{m=0}^{m}\frac{n!}{(n-m)!*m!}) combinations for *n* interaction filters and *m* max_miss. Will log and output union of combinations unless used with `enumerate_interaction_combs`. | 0  | <tr><td colspan="5">***Storage and read options***</td></tr>Yes |
|log_file              | Name for log of filtered results                      | output_log.txt   ||
|overwrite        | Flag to overwrite existing logfile of same name           | FALSE       ||
|bookmark_name      | Name for bookmark view in database                      | passing_results  ||
|outfields       | Data fields to be written in output (log file and STDOUT). Ligand name always included. | e        ||
|order_results    | String for field by which the passing results should be ordered in log file. | no default ||
|output_all_poses        | Flag that if mutiple poses for same ligand pass filters, log all poses | (OFF)        ||
|mfpt_cluster     | Cluster ligands passing given filters based on the Tanimoto distances of the Morgan fingerprints. Will output ligand with best (lowest) ligand efficiency from each cluster. Uses Butina clustering algorithm | 0.5  ||
|interaction_cluster     | Cluster ligands passing given filters based on the Tanimoto distances of the interaction fingerprints. Will output ligand with best (lowest) ligand efficiency from each cluster. Uses Butina clustering algorithm | 0.5  | Yes |
|enumerate_interactions_combs  | When used with `max_miss` > 0, will log ligands/poses passing each separate interaction filter combination as well as union of combinations. Can significantly increase runtime. | FALSE  | Yes|

### Output options
There are a number of output methods available to filter, view, and store the results. 

| Availble output methods          | Description                                           |  Arguments   | 
|:------------------------|:-------------------------------------------------|:----------------|
|`export_csv`| Name of database result bookmark or table to be exported as CSV. Output as <table_name>.csv | requested_data= bookmark_name, csv_name, table=True |
|`export_csv`| Create csv of the requested SQL query. Output as query.csv. MUST BE PRE-FORMATTED IN SQL SYNTAX e.g. SELECT [columns] FROM [table] WHERE [conditions] |requested_data = query string, csv_name, table=False|
|`export_bookmark_db` | Export a database containing only the results found in the specified bookmark name. Will save as <core_db_file>_<bookmark_name>.db| bookmark_name |
|`export_receptors`| Export receptor to pdbqt | None |
|`write_molecule_sdfs`| Write molecule sdfs from a given bookmark to specified path  |  sdf_path, bookmark_name   |
|`find_similar_ligands`|  Given query ligand name, find ligands previously clustered with that ligand. User prompted at runtime to choose cluster group of interest. | query_ligname |
|`get_previous_filter_data`| Get data requested in `outfields` from the bookmark of a previous filtering | outfields: str, bookmark_name" str |
|`find_similar_ligands`| Find ligands in cluster with query_ligname |query_ligname|
|`plot`| Freate scatterplot of ligand efficiency vs docking score for best pose of each ligand. Saves as "scatter.png". | save: bool |
|`pymol`| Launch interactive LE vs Docking Score plot and PyMol session. Ligands in the bookmark specified with bookmark_name will be ploted and displayed in PyMol when clicked on.  | bookmark_name |

### Using the config file
Both the command line tool and the API can make use of a configuration file, but this option is not fully vetted yet. 