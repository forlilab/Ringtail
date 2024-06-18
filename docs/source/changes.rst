.. _changes:

Changes in Ringtail
######################

Changes in 2.0.0: fully developed API
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
* Can add docking results directly without using file system (for vina only as output comes as a string). See this scripting example: (ref to vina meeko ringtail from Diogo)
* The Ringtail log is now written to a logging file in addition to STDOUT

Changes to code behavior
=========================
* Interaction tables: one new table has been added (`Interactions`). The existing `Interaction_indices` table and the table `Interaction_bitvectors` are remade every time the database is written to as opposed to being made on the go as results are added in previous Ringtail version. 
* A new method to update an existing database 1.1.0 (or 1.0.0) to 2.0.0 is included. However, if the existing database was created with the duplicate handling option, there is a chance of inconsistent behavior of anything involving interactions as the Pose_ID was not used as an explicit foreign key in db v1.0.0 and v1.1.0 (see Bug fixes below).

Bug fixes
===========
* The option `duplicate_handling` could previously only be applied during database creation and produced inconsistent table behavior. Option can now be applied at any time results are added to a database, and will create internally consistent tables. **Please note: if you have created tables in the past and invoking the keyword `duplicate_handling` you may have errors in the "Interaction_bitvectors" table. These errors cannot be recovered, and we recommend you re-make the database with Ringtail 2.0.0.**


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
