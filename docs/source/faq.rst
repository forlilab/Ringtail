.. _faq:

Frequently asked questions
#############################

Potential pitfalls
**********************
Using the command line tool: any PDBQT files specified through any of the input options in ADGPU mode will be read by `rt_process_vs` as receptor files, even if the files actually represent ligands. Therefore, ligand PDBQT files should not be present in any directories given with `--file_path`.

When writing from Vina PDBQTs, ensure there are no other PDBQTs (input or receptor) in directories specified with `file_path` UNLESS the receptor PDBQT is specified with the `receptor_file` option in the same command line/method call.

Occassionally, errors may occur during database reading/writing that corrupt the database. This may result in the database becoming locked. If this occurs it is recommended to delete the existing database and re-write it from scratch.

SQLite will make a new, empty database file whenever a creation is connected, i.e., whenever Ringtail is used. In Ringtail write-mode, if a database does not already exist, this is not a problem, and a database will be created and populated as expected. If Ringtail is used in read-mode (e.g., filtering) and there is a spelling error in the database file name Ringtail will throw an error since the database is empty. However, since a connection was made to the misspelled db to check if the db was empty or not, a new empty db file will be created with the misspelled name. Since the file is empty, Ringtail will keep throwing errors if you try to write to or read from it, and if you wish to proceed e.g., writing to the database with the formerly misspelled name, you will have to delete the empty file first. 
