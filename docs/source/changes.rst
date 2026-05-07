.. _changes:

Changes in Ringtail
######################
Changes in 3.x: New graphical user interface and support for DuckDB
********************************************************************

Changes in command line tools
==================================================
* New filter options for molecular weight `--ligand_min_molweight` and `--ligand_max_molweight`
* There is now one upgrade script with version as input, where there previously was one db upgrade script per version
* `--storage_type` can be used to specify database engine, with option to use 'duckdb' (defualts to sqlite)
* `--docking_mode` now only an option for `write` as it is not relevant for the `read` processes
* Writing a filter results "log" file has a new command line keyword `--output_log` (still uses shorthand `-l`), is now optional, and will only be done if `--output_log` is specified
* The command `--logfile` will write logging output to a file
* The previous `--filter_bookmark` has been changed to `--input_bookmark` for consistency
* New CLI flag `--print_bookmarks` prints all current bookmarks in the database
* `--file_pattern` has been discontinued, and is inferred from docking mode (which will attempt to accept file pattern as input)
* `--export_bookmark_csv` can now be used as a flag (True/False) to export bookmark given in `--bookmark_name`, or with string input as the resulting csv file name. If used in conjunction with `--bookmark_name` and `--outfields` it will produce a csv with the desired columns. Full tables can be exported using the `--bookmar_name` tag, but this will not work with `--outfields`
* `--outfields` now uses names of columns as they are in the database from the Results, Ligands, and Interaction_indices tables. This inculdes new fields from the interaction table, and that some old fields have changed:
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


Enhancements to the codebase
==============================
***** Fully developed API can use python for scripting exclusively (see :ref:`API <api>` page for full description)
* DuckDB offered as an alternative to SQLite, with overall similar database creating times and significantly enhanced filtering times
* `storage_type` needs only be specified when a database is first created, storage_type will automatically be detected for an existing database (detected storage_type takes precedence over specified storage_type)
* New dataclass `RingtailDefaults` which are used to a larger extend in method signatures where appropriate
* Using toml 
* Validation of docking mode adds some flexibility in how a docking engine is reference to, eg ADGPU vs GPU
* Multiprocess now uses 'fork' start method for compatibility with multithreaded processes such as the GUI
* Ligand efficiency, a calculated value, is rounded to two decimal points reflect the accuracy of the numbers used to calculate it (docking score and number of atoms)
* Database write is faster by parsing a larger number of docking results to memory before committing to the database (previously 1 at the time, now 10,000)
* `docking_mode` is no longer a property of the Ringtail object, only an argument for writing to the database (i.e., `add_results_from_files`)
* `access_mode` is an optional initialization argument which alters the behavior of some API methods (defaults to `api`)
* Additional filters for minimun and maximum ligand molecular weight, `ligand_min_molweight` and `ligand_max_molweight`
* Will only write a filter log file (e.g., `output_log.txt`) if specified, which significantly enhances filtering speeds when not used
* New API method for clustering across a pre-filtered bookmark (can technically also cluster all Results)
* Bookmarks are no longer saved as views (i.e., unrealized tables), instead the criteria of a bookmark (e.g., a dict of filters) is saved in the new `Filters` table with a unique filter_id, and all passing pose_id`s with associated filter_id`s are stored in `Filtered_poses`
* The `write_flexres_pdb` method now allows more than one ligand input, for example by providing a bookmark name all ligands in that bookmark will be used to write the same number of PDBs (there will be a warning of attempting to write more than 10 files)
* The method `export_receptor` is now `export_receptor_pdbqt`, and a modernized version to `write_flexres_pdb` has been added which is compatible with the new receptor Polymer object (this new method works without flexible residues as well)
* `write_molecule_sdfs` method argument `write_nonpassing` has been discontinued, and `ligname` (string or list of strings) has been added. It is assumed that if a `bookmark_name` is provided, only passing poses of each ligand will be written to an SD file. If no `bookmark_name` is provided, each pose of each ligand is written to the SDF.  
* The method `export_csv` has been broken into three distinct methods, `export_columns_as_csv` where one or more columns (from Results and Ligands tables + modified interaction columns) are specified and exported, `export_table_as_csv` where an entire table is exported, and `export_sql_as_csv` where the user specifies a properly formatted SQL prompt
* The method `export_bookmark_db` 
* Status selection enabled using `set_ligand_status` API, three new tables included: Accepted, Maybe, Rejected
* The method `get_plot_data` now has more input parameters including specifying x and y axis, and returning the status of a pose/point 
* The method `drop_bookmark` is now `delete_bookmark`
* New method `get_bookmark_interactions` to get interaction data from a bookmark
* Several new APIs to support the GUI, generally not useful outside the GUI
* New method `merge_databases` which will safely merge a secondary database with the database currently initialized as a Ringtail object. 
* New method to assign status (accepted, rejected, maybe) to one or more ligands, poses, and/or ligand poses in a bookmark
* New method to (re-)calculate interactions for vina and ng results. This will delete all interactions present and calculate them all anew based on current receptor data in the database and given vdw and hb cutoffs.

Changes to code behavior
=========================
* For the methods `filter()` and `cluster()` the keys `bookmark_name` and `filter_bookmark` have been changed to `output_bookmark` and `input_bookmark`, respectively, for clarity
* The column `nr_interactions` in the Results table is now called `num_interactions` for consistency with `num_hb`
* The column `ligand_coordinates` in the Results table is now called `pose_coordinates`
* The column `deltas` in the Results table is now called `delta`
* Ringtail bookmarks from e.g., filtering clustering were previously created as database views, which appear as tables that are unrealized until viewing them. This has been replaced by a `Filters` table which holds filter information (previous equivalent was `Bookmarks` table) and the poses passing a given filter are stored in a tall-skinny table `Filtered_poses`. A significant speed increase was enabled by this move, and any other behavior related to bookmarks is the same. 
* The following columns have been removed from the Ligands table (information now stored in the binary rdkit Mol): `atom_index_map`, `hydrogen_parents`, and `input_model`.

Bug fixes
===========
* For vina results, special docking atoms (for macrocycles and waters) may have been contributing to calculated van der Waals interactions in the database. This is no longer the case, so if e.g., a database is recreated in v3.0.0 from the original docking .PDBQTs the new database may have fewer interactions. 


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