.. _faq:

Frequently asked questions
#############################

Potential pitfalls
**********************
Using the command line tool: any PDBQT files specified through any of the input options in ADGPU mode will be read by `rt_process_vs` as receptor files, even if the files actually represent ligands. Therefore, ligand PDBQT files should not be present in any directories given with `--file_path`.

When writing from Vina PDBQTs, ensure there are no other PDBQTs (input or receptor) in directories specified with `file_path` UNLESS the receptor PDBQT is specified with the `receptor_file` option in the same command line/method call.

Occassionally, errors may occur during database reading/writing that corrupt the database. This may result in the database becoming locked. If this occurs it is recommended to delete the existing database and re-write it from scratch.

