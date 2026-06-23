.. _changes:

Changes in Ringtail
######################
Changes in 3.0: Faster and smaller, and lots of new functionality
******************************************************************

Ringtail 3 comes with support for the new ADNG docking engine <link>, now uses the DuckDB database engine by default, and includes an overhauled database schema leading to faster screening and filtering. There are new output options including CSV export with custom column choices without needing to know SQL, and new or improved command line tools to handle database merging and compression/decompression. Overall a Ringtail database now takes up less space, and is significantly faster filter and screen. 

Enhancements to the codebase
==============================
* Works with ADNG SDF file docking output, and allows receptor to be provided as a `json`
* DuckDB backend offered as an alternative to SQLite, with overall similar database creating times and significantly enhanced filtering times and smaller file size
* Abandoned using views to store filtered poses in favor of a long-and-skinny `Filtered_poses` and `Filters` tables, significantly speeding up filtering and especially progressive filtering
* Multiprocess now uses 'fork' start method for compatibility with multithreaded processes and enables Ringtail on Windows machines #TODO 
* A larger default chunk size of docking data is parsed before writing to the database 
* Package handling modernized to use pyproject.toml 
* Additional filters allow specification of minimun and maximum ligand molecular weight, `ligand_min_molweight` and `ligand_max_molweight`
* The method `export_bookmark_db` uses enhanced logic which speeds up the creation of a new subset database
* It's now possible to assign status/flags to poses via the API using `update_pose_status`, such as Accepted/1, Maybe/2, Rejected/3 or 0 to remove status
* New method `merge_databases` which will safely merge one or more secondary database with the database currently initialized as a Ringtail object. 
* New method to (re-)calculate interactions for vina and ADNG results. This will delete all interactions present and calculate them all anew based on current receptor data in the database and given vdw and hb cutoffs. Useful if interactions were not calculated during database creation, or if user wants to re-calculate them with new interaction cutoff distances. 
* Interaction calculations uses a k-d tree of the receptor atoms and batched lookup for all provided poses, built once instead of per pose as previously done, and all atoms in the pose are checked in batch instead of one by one, significantly reducing time to calculate interactions 
* Ringtail database version is now tracked in a `ringtail_schema_version` table (SQLite `PRAGMA user_version` has been depreceated)
* A more flexible Pytests harness for advanced users and dvelopers


Changes in command line tools
==================================================
* New field `--docking_results`, `--dr` accepts any file input including one or more folders/file paths, one or more single files, and one or more file lists (.txt). The depreceated fields `--file`, `--file_list`, and `--file_path` will remain for compatibility with existing scripts.
* `--add_interactions` has been replaced with `--no_interactions`/`-ni`, making calculating interactions the default
* New filter options for molecular weight `--ligand_min_molweight` and `--ligand_max_molweight`
* New filter input `--ligand_name_file` allows providing a .csv file of ligand names, which will be applied as a filter (e.g., for exporting SDFs for select ligands). Works the same as `--ligand_name` but allows for a larger number of provided names
* There is now one upgrade CLI `rt_upgrade_db` script with version as input
* New CLI scripts `rt_compress_db` and `rt_decompress_db` to compress and decompress a database, with optional `--eworst` and `--leworst` filters applied
* `--storage_type` can be used to specify database engine (defaults to duckdb)
* `--docking_mode` now only required during `write` 
* Writing a filter results "log" file has a new command line keyword `--output_log` (still uses shorthand `-l`), is now optional, and will only be done if `--output_log` is specified
* The command `--logfile` will write standard logger output to a file
* The previous `--filter_bookmark` has been changed to `--input_bookmark` for consistency
* New CLI flag `--print_bookmarks` prints all current bookmarks in the database
* `--file_pattern` has been discontinued, and is inferred from docking mode (which will attempt to accept file pattern as input)
* `--export_bookmark_csv` can now be used as a flag (True/False) to export bookmark given in `--bookmark_name`, or with string input as the resulting csv file name. If used in conjunction with `--bookmark_name` and `--outfields` it will produce a csv with the desired columns. Full tables (like `Interactions`) can be exported using the `--bookmar_name` tag, but this will not work with `--outfields`
* `--outfields` now uses the names of columns as they are in the database from the Results, Ligands, and Interaction_indices tables (and not a mix of column names and aliases, as before). This inculdes new fields from the interaction table, and that some old fields have changed:
    ============ ===============
      Old          New
    ============ ===============
    ligand_name   ligname 
    e             docking_score
    le            leff
    ref_rmsd      reference_rmsd
    e_inter       energies_inter
    e_vdw         energies_vdw
    e_elec        energies_electro
    e_intra       energies_intra
    rank          pose_rank
    run           run_number
    hb            num_hb
    ============ ===============

Changes to API and code behavior
================================
* Ringtail auto-detects `storage_type` for existing databases
* Methods that create bookmarks, such as filter() and cluster() now uses `output_bookmark` to name the new, resulting bookmark, instead of `bookmark_name`.
* Simplified `add_results_from_files` API has only one file input field `docking_results` which will accept a single or a list of, and a mix of files, folders, and lists of file paths. `file`, `file_list`, and `file_path` have been depreceated from the API. 
* `docking_mode` specification is more flexible, e.g., 'adgpu', 'gpu', and 'dlg' are all valid for AutoDock-GPU docking mode
* Database schema is now fully defined in `schema.py`, and any database table and column info is derived from this single source of truth
* For the methods `filter()` and `cluster()` the keys `bookmark_name` and `filter_bookmark` have been changed to `output_bookmark` and `input_bookmark`, respectively, for clarity
* The column `nr_interactions` in the Results table is now called `num_interactions` for consistency with `num_hb`
* The column `ligand_coordinates` in the Results table is now called `pose_coordinates`, and these coordinates are no longer stored as a string (sqlite: float32 BLOB and duckdb: float array)
* The column `deltas` in the Results table is now called `delta`
* The following columns have been removed from the Ligands table (information now stored in the binary rdkit Mol): `atom_index_map`, `hydrogen_parents`, and `input_model`.
* `create_rdkit_mol` now `create_rdkit_mols` and supports more efficient database fetching of binary rdkit write_molecule_sdfs
* The method `drop_bookmark` is now `delete_bookmark`
* The API `find_similar_ligads` has been replaced by `fetch_cluster_options` and `fetch_clustered_similars`, respectively
* `docking_mode` is no longer a property of the Ringtail object, only an argument for writing to the database (i.e., `add_results_from_files`)
* `access_mode` is an optional initialization argument which alters the behavior of some API methods 
* The `write_flexres_pdb` method now allows more than one ligand input, for example by providing a bookmark name all ligands in that bookmark will be used to write the same number of PDBs (there will be a warning of attempting to write more than 10 files)
* The method `export_receptor` is now `export_receptor_pdbqt`, and a modernized version to `write_flexres_pdb` has been added which is compatible with the new receptor Polymer object (this new method works without flexible residues as well)
* `write_molecule_sdfs` method argument `write_nonpassing` has been discontinued, and `ligname` (string or list of strings) has been added. It is assumed that if a `bookmark_name` is provided, only passing poses of each ligand will be written to an SD file. If no `bookmark_name` is provided, each pose of each ligand is written to the SDF.  
* The method `export_csv` has been broken into three distinct methods, `export_columns_as_csv` where one or more columns (from Results and Ligands tables + modified interaction columns) are specified and exported, `export_table_as_csv` where an entire table is exported, and `export_sql_as_csv` where the user specifies a properly formatted SQL prompt
* The options for creating plots and opening PyMol sessions (and associated methods) through the CLI have been discontinued 
* The class ResultsManager, designed to handle different multi processor options, has been removed

Bug fixes
===========
* For vina results, special docking atoms (for macrocycles and waters) may have been contributing to calculated van der Waals interactions in the database. This is no longer the case, so if e.g., a database is recreated in v3.0.0 from the original docking .PDBQTs the new database may have fewer interactions. 
* Ligand efficiency, a calculated value, is rounded to two decimal points reflect the accuracy of the numbers used to calculate it (docking score and number of atoms)
* Will only write a filter log file (e.g., `output_log.txt`) if specified
* Exporting poses with flexible receptor residues will now export all poses of a given ligand, not just the best scoring one 

#TODO add 2.3.1
Changes in 2.1.1: bug fixes and result plot enhancements
********************************************************
Enhancements
============
* Data used for plotting, and a handle to the matplotlib.figure object is available through the API for personalized plotting
* The appearance of the standard Ringtail docking results plot has been improved with hard cutoffs for docking scores and ligand efficiency above 0, and it has been made clearer the difference between plotted data that has been binned and plotted data that are single data points. The number of bins used to bin data as well as marker size is scaled to the amount of data for databases of less than 10,000 ligands for enhanced visibility. 
* Plotting faster for large databases. 
* Possible to choose printing all filtered ligands to one large, or multiple individual SDF files for the CLI using `--invidual_sdf_files`.

Bug fixes
=========
* There was a bug when using "overwrite" a database with the CLI where Ringtail would overwrite the database at the wrong time, leading to issues with data in the receptor table and functions like "add_interactions" (as of v3 please note add_interactions is no longer an option). This has now been fixed.
* Multiple bugs related to the ligand filters have been fixed, where the main issue was if the `ligand_operator` was set to "OR", "OR" would be used not only between ligand substructures but also between other parts of the query (leading to more ligands passing filtering). 


Changes in 2.1.0: enhanced filtering speed
******************************************
Enhancements to the code base
==============================
* The format of the queries produced to filter the database have been completely rewritten, reducing filtering time by at least a factor of 10 compared to 1.1.0. Extra indices were added to three of the tables to support the faster filtering speeds. 

Bug fixes
===========
* The use of the keywords `--ligand_name`, `--ligand_substruct`, and `--ligand_substruct_pos` had ambiguous behavior where if they were invoked more than once, only the last filter value would be used (as opposed to concatenating the values). They now will work by supplying multiple values to one keyword, as well as one or more values to two or more keywords. Further, `ligand_substruct_pos` now takes input as one string (`"[C][Oh] 1 1.5 -20 42 -7.1"`)as opposed to one string and five numbers (`"[C][Oh]"" 1 1.5 -20 42 -7.1`).
* `--ligand_max_atoms` counted all atoms in the ligand, including hydrogens. With bug fix it counts only heavy atoms(not hydrogens). 

Changes in 2.x: fully developed API
***************************************

Changes in keywords used for the command line tool
==================================================
* `--mode` is now `--docking_mode`
* `--summary` is now `--print_summary`
* `--pattern` is now `--file_pattern`
* `--name` is now `--ligand_name`
* `--max_nr_atoms` is now `--ligand_max_atoms`
* `--smarts` is now `--ligand_substruct`
* `--smarts_idxyz` is now `--ligand_substruct_pos`
* `--smarts_join` is now `--ligand_operator`
* `--van_der_waals` is now `--vdw_interactions`
* `--hydrogen_bond` is now `--hb_interactions`
* `--reactive_res` is now `--reactive_interactions`

Enhancements to the codebase
==============================
* Fully developed API can use python for scripting exclusively (see :ref:`API <api>` page for full description)
* Can add docking results directly without using file system (for vina only as output comes as a string). 
* The Ringtail log is now written to a logging file in addition to STDOUT if log level is det to "DEBUG". 

Changes to code behavior
=========================
* Interaction tables: one new table has been added (`Interactions`) which references the interaction id from `Interaction_indices`, while the table `Interaction_bitvectors` has been discontinued.
* A new method to update an existing database 1.1.0 (or 1.0.0) to 2.0 is included. However, if the existing database was created with the duplicate handling option, there is a chance of inconsistent behavior of anything involving interactions as the pose_id was not used as an explicit foreign key in db v1.0.0 and v1.1.0 (see Bug fixes below).

Bug fixes
===========
* The option `duplicate_handling` could previously only be applied during database creation and produced inconsistent table behavior. Option can now be applied at any time results are added to a database, and will create internally consistent tables. **Please note: if you have created tables in the past and invoking the keyword `duplicate_handling` you may have errors in the "Interaction_bitvectors" table (<2.0). These errors cannot be recovered, and we recommend you re-make the database with Ringtail 2.0.**
* Writing SDFs from filtering bookmarks: will check that bookmark exists and has data before writing, and will now produce SDFs for any bookmarks existing bookmarks. If the bookmark results from a filtering where `max_miss` &lt; 0 it will note if the non-union bookmark is used, and if the base name for such bookmarks is provided it will default to the `basename_union` bookmark for writing the SDFs.
* Output from filtering using `max_miss` and `output_all_poses=False`(default) now producing expected behavior of outputting only one pose per ligand. Filtering for interactions `max_miss` allows any given pose for a ligand to miss `max_miss` interactions and still be considered to pass the filter. Previously, in the resulting `union` bookmark and `output_log` text file some ligands would present with more than one pose, although the option to `output_all_poses` was `False` (and thus the expectation would be one pose outputted per ligand). This would give the wrong count for how many ligands passed a filter, as some were counted more than once. 
* The use of the keywords `--ligand_name`, `--ligand_substruct`, and `--ligand_substruct_pos` had ambiguous behavior where if they were invoked more than once, only the last filter value would be used (as opposed to concatenating the values). They now will work by supplying multiple values to one keyword, as well as one or more values to two or more keywords. Further, `ligand_substruct_pos` now takes input as one string (`"[C][Oh] 1 1.5 -20 42 -7.1"`)as opposed to one string and five numbers (`"[C][Oh]"" 1 1.5 -20 42 -7.1`).

Changes in 1.1.0: enhanced database performance
***********************************************

Database operations
====================
* Significant filtering runtime improvements vs v1.0.0 by using multithreaded processing
* Added the ability to print a `summary` to stdout for getting quick overview of data across entire dataset
* Added dability to export receptors stored in the database as a receptor PDBQT

Filtering and querying
=======================
* Can now select of dissimilar output ligands with Morgan fingerprint or interaction fingerprint clustering
* Can now select similar ligands from querying a ligand name used in previous Morgan fingerprint or interaction finger clustering groups
* Can filter by substructures present in the ligand 
* Can filter by ligand substructure location in cartesian space
* The option to specify how many interaction filter combinations is OK to be missed (`max_miss`) now defaults to outputting the union of interaction combinations, and when used in conjunction with the `enumerate_interaction_combs` option will log passing ligands/poses for individual interaction combination