.. _changes:

Changes in Ringtail
######################

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
* A new method to update an existing database 1.1.0 (or 1.0.0) to 2.0 is included. However, if the existing database was created with the duplicate handling option, there is a chance of inconsistent behavior of anything involving interactions as the Pose_ID was not used as an explicit foreign key in db v1.0.0 and v1.1.0 (see Bug fixes below).

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