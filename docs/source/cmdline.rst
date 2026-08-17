.. _cmdline:

Command line interface
#######################

The Ringtail command line interface is the easiest method to use to for exploring virtual screening results in a database. 
If this is your first time learning about Ringtail, take a look at the page :ref:`Get started <get_started>`. The current page uses the knowledge already presented on the 'Get started' page as we continue exploring the wealth of options that Ringtail offers. 

The script for writing a database and filtering is ``rt_process_vs``. This is intended to be used for a set of DLGs/Vina PDBQTs pertaining to a single target and binding site. This may include multiple ligand libraries as long as the target and binding site is the same. Be cautious when adding results from multiple screening runs. Ringtail only allows one receptor per database, and will do this by checking that the receptor name is the same. However, if you do multiple screens with e.g., different binding pockets for the same receptor, this information will not be checked, and Ringtail will allow any data to be appended to a database as long as the receptor name is the same. One receptor PDBQT may also be saved to the database.

The rt_process_vs script has two modes: ``write`` and ``read``. The desired mode must be specified in the command line before any other options are given (except ``-c [CONFIG]`` which is given first). The ``write`` mode is used to create a database for a virtual screening from ADGPU DLGs or Vina PDBQTs. After this initial run, a database is created and may be read directly by rt_process_vs in ``read`` mode for subsequent filtering and export operations.

Please note that Ringtail does not automatically have permission to perform changes outside of the working directory, so be advised that any folders or documents that Ringtail outputs will be saved in the current working directory. 

Ringtail inputs
*****************

Navigate to the data repository and chose one of several paths of adding results. 

.. code-block:: bash

  $ cd test/test_data/adpgu

Input file_sources
===================
By default, the newly-created database will be named ``output.db``. This name may be changed with the ``--output_db`` or ``-o`` option.
Ringtail allows multiple formats to provide docking results using ``--docking_results`` or ``-dr``, including paths to one or more single files, path to one or more .txt files containing docking file paths, and paths to one or more directories containing files (use ``--recursive``/``-r`` to search recursively traverse folder). You can specify one or more arguments at once:

.. code-block:: bash

    $ rt_process_vs write --docking_results lig1.dlg --docking_results lig2.dlg  folder1/ folder filelist1.txt filelist2.txt --output_db example.db

Example file list:

.. code-block:: python

    lig3.dlg
    lig4.dlg.gz
    rec1.pdbqt

Input options
===============
To include the details of a receptor in the database, it is necessary to provide explicitly state that the receptor should be saved. If the ``--save_receptor`` argument is invoked, but no PDBQT ``--receptor_file`` is provided, Ringtail will raise an error. During results processing Ringtail checks to make sure the provided receptor file matches the receptor name used for the docking. 

.. code-block:: bash

    $ rt_process_vs write --docking_results filelist1.txt --receptor_file test_data/4j8m.pdbqt.gz --save_receptor

It is possible to add docking results *or* a receptor file to a database that already exists. For this it is necessary to use the keyword ``--append_results``.
You can also specify what to do if you are adding duplicate results for a ligand, by invoking the ``--duplicate_handling`` keyword with the value ``IGNORE`` (will not add the newest duplicate) or ``REPLACE`` (will overwrite the newest duplicate). Please note that the ``--duplicate_handling`` option makes database writing significantly slower.

.. code-block:: bash

    $ rt_process_vs write --input_db output.db --docking_results test_data/group2 --append_results --duplicate_handling REPLACE

By default (for DLGs), Ringtail will store the best-scored (lowest energy) binding pose from the first 3 pose clusters in the DLG. For Vina, Ringtail will store the 3 best poses. Additional settings for writing to the database include how to handle the number of poses docked (``--max_poses``, or ``--store_all_poses`` which will overwrite the former).

ADGPU is capable of performing interaction analysis at runtime, with these results being stored in the database if present. If interaction analysis is not present in the input file (including Vina PDBQTs), it will by default be added by Ringtail during data parsing, unless the user specifies not to by using the ``--no_interactions`` flag. **This will signifcantly decrease the total database write time.** By default, the cutoff distance for being considered an interaction is 3.7 å for hydrogen bonds and 4.0å for van der Waals interactions. If the user would like to calculate interactions for vina results with other distance cutoffs, the option ``--interaction_cutoffs`` can be used. To be able to add interactions it is important that the ``--receptor_file`` is provided during database write (or that the receptor has already been saved explicitly in the database). 

The ``--interaction_tolerance`` option also allows the user to give more leeway for poses to pass given interaction filters. With this option, the interactions from poses within *c* angstrom RMSD of a cluster's top pose will be appended to the interactions for that top pose. The theory behind this is that this gives some sense of the "fuzziness" of a given binding pose, allowing the user to filter for interactions that may not be present for the top pose specifically, but could be easily accessible to it. When used as a flag, the ``--interaction_tolerance`` default is 0.8 angstroms. The user may also specify their own cutoff. This option is intended for use with DLGs from AD-GPU, which clusters output poses based on RMSD.

It is further possible to overwrite a database by use of the argument ``--overwrite``.

.. code-block:: bash

    #AD-GPU
    $ rt_process_vs write --input_db output.db --docking_results test_data/group1 --max_poses 2 --interaction_tolerance 0.8

    #vina
    $ rt_process_vs write --input_db output.db --docking_results test_data/vina --overwrite --receptor_file receptor.pdbqt --save_receptor --interaction_cutoffs 3.5,4.5

Printing a database summary
***************************
During both ``write`` and ``read`` it is possible to add the tag ``-su`` or ``--print_summary`` which will print a summary of the database to stdout.

.. code-block:: bash

    $ rt_process_vs read --input_db output.db -su

    Total Stored Poses: 645
    Total Unique Interactions: 183

    Energy statistics:
    min_docking_score: -7.93 kcal/mol
    max_docking_score: -2.03 kcal/mol
    1%_docking_score: -7.43 kcal/mol
    10%_docking_score: -6.46 kcal/mol
    min_leff: -0.62 kcal/mol
    max_leff: -0.13 kcal/mol
    1%_leff: -0.58 kcal/mol
    10%_leff: -0.47 kcal/mol

Filtering
*********
In ``read`` mode, an existing database is used to filter or export results.

When filtering, a text results log file will be created containing the results passing the given filter(s). The default log name is ``output_log.txt`` and by default will include the ligand name and docking score of every pose passing filtering criteria. The log name may be changed with the ``--output_log`` option. There are six scoring filters that include best and worst docking score/energy, best and worst ligand efficieny, and results above worst docking score or ligand efficiency percentile. Some of these are internally inconsistent: if both ``--eworst`` and ``--score_percentile`` are used together, the ``--eworst`` cutoff alone is used. The same is true of ``--leworst`` and ``--le_percentile``.

Scoring filters
=================

.. code-block:: bash

    $ rt_process_vs read --input_db output.db --score_percentile 0.1 --output_log output_log_01percent.txt

The information written to the log can be specified with ``--outfields``. The full list of available output fields may be seen by using the ``--help`` option with ``read`` mode.
By default, only the information for the top-scoring binding pose will be written to the log. If desired, each individual passing pose can be written by using the ``--output_all_poses`` flag. The passing results may also be ordered in the log file using the ``--order_results`` option.

.. code-block:: bash

    $ rt_process_vs read --input_db output.db --eworst -6 --outfields ligname,docking_score,pose_rank --order_results ref_rmsd --bookmark_name eworst6

When filtering, the passing results are also saved as a view in the database. This view is named ``passing_results`` by default. The user can specify a name for the view using the ``--bookmark_name`` option. No filtering is performed if no filters are given (see full list of filters :ref:`here <filter_kw_table>`). 
Filtering may take from seconds to minutes, depending on the size of the database, roughly scaling as O(n) for n database Results rows (i.e. stored poses). Data for poses in a view may be accessed later using the ``--data_from_bookmark`` option.

Interaction filters
=====================
It is possible to filter the docking results based on different types of interactions (hydrogen bonds onr van der waals) with specific residues. It is further possible to have ligands pass the filters while only fulfilling some of the interaction combinations in union (max number of interactions combinations missed, ``--max_miss``).
The available interaction filters are ``--hb_interactions``, ``--vdw_interactions``, and ``--reactive_interactions``. Interaction filters must be specified in the order ``CHAIN:RES:NUM:ATOM_NAME``. Any combination of that information may be used, as long as 3 colons are present and the information ordering between the colons is correct. All desired interactions of a given type (e.g. ``-vdw``) may be specified with a single option tag (``-vdw B:THR:276:,B:HIS:226:``) or separate tags (``-vdw B:THR:276: -vdw B:HIS:226:``).

The ``--max_miss`` option allows the user to filter by given interactions excluding up to ``max_miss`` interactions. This gives :math:`\sum_{m=0}^{m}\frac{n!}{(n-m)!*m!}` combinations for *n* interaction filters and *m* max_miss. By default, results will be given for the union of the interaction conbinations. Use with ``--enumerate_interaction_combs`` to log ligands/poses passing each separate interaction combination (can significantly increase runtime). ßIf ``max_miss > 0`` is used during filtering, a view is created for each combination of interaction filters and is named ``<bookmark_name>_<n>`` where n is the index of the filter combination in the log file (indexing from 0).
``--react_any`` offers an option to filtering for poses that have reactions with any residue.

.. code-block:: bash

    $ rt_process_vs read --input_db output.db --eworst -6 --hb_interactions A:VAL:279: A:LYS:162: --vdw_interactions A:VAL:279: A:LYS:162: --max_miss 1 --react_any

Ligand filters 
=================
The docked ligands can be filtered for presence of certain substrctures specified by their SMARTS string using ``--ligand_substruct``, as well as their ligand name contaning a specific phrase ``--ligand_name``. The ligand name search will include any ligand names that contain the specified phrase, and does not look for exact matches only. 
Use the keyword ``--ligand_operator`` to determine if the ligand filters should be evaluated as this ``OR`` that (default), or combined with ``AND``. ``--ligand_max_atoms`` can be used to specify maximum number of heavy atoms a ligand may have.

.. code-block:: bash

    $ rt_process_vs read --input_db output.db --ligand_substruct 'C=O' 'CC(C)(C)' --ligand_operator AND --ligand_max_atoms 5

The ``--ligand_substruct_pos`` option may be used to filter for a specific ligand substructure to be placed within some distance of a given cartesian coordinate. The format for this option is the six elements inside quotes and separated by spaces: ``"<SMARTS pattern: str> <index of atom in SMARTS: int> <cutoff distance: float> <target x coord: float> <target y coord: float> <target z coord: float>""``. 

.. code-block:: bash

    $ rt_process_vs read --input_db output.db --ligand_name cool_ligand --ligand_substruct_pos "[C][Oh] 1 1.5 -20.3 42 -7.1"

Clustering
============
In addition to the filtering options outlined in the table below, ligands passing given filters can be clustered to provide a reduced set of dissimilar ligands based on Morgan fingerprints (``--mfpt_cluster``) or interaction (``--interaction_cluster``) fingerprints. Dissimilarity is measured by Tanimoto distance and clustering is performed with the Butina clustering algorithm. Clustering can be also be performed on a bookmark that has already been saved to the database, without providing any extra filter values. In this case, the bookmark over which to cluster (or additional filtering) on is specified by ``--input_bookmark`` (must be different from ``--bookmark_name``).

.. code-block:: bash

    $ rt_process_vs read --input_db output.db --input_bookmark eworst6 --mfpt_cluster

While not quite a filtering option, the user can provide a ligand name from a previously-run clustering and re-output other ligands that were clustered with that query ligand with ``--find_similar_ligands``. The user is prompted at runtime to choose a specific clustering group from which to re-output ligands. Filtering/clustering will be performed from the same command-line call prior to this similarity search, but all subsequent output tasks will be performed on the group of similar ligands obtained with this option unless otherwise specified. 

Outputs
*********
The primary outputs from ``rt_process_vs`` are the database itself (``write`` mode) and the filtering log file (``read`` mode). There are several other output options as well, intended to allow the user to further explore the data from a virtual screening.

Using the ``--export_sdf_path`` option allows the user to specify a directory to save SDF files for ligands passing the given filters or in the bookmark given with ``--bookmark_name``. The SDF will contain poses passing the filter/in the bookmark ordered by increasing docking score, optionally the argument ``--individual_sdf_files`` can be used to output one SDF per ligand. Each ligand is written to its own SDF. This option enables the visualization of docking results, and includes any flexible/covalent ligands from the docking. The binding energies, ligand efficiencies, and interactions are also written as properties within the SDF file, with the order corresponding to the order of the pose order.

If the user wishes to explore the data in CSV format, Ringtail provides two options for exporting CSVs. The first is ``--export_bookmark_csv``, which takes a string for the name of a table or result bookmark in the database and returns the CSV of the data in that table. The file will be saved as ``<table_name>.csv``.
The second option is ``--export_query_csv``. This takes a string of a properly-formatted SQL query to run on the database, returning the results of that query as ``query.csv``. This option allows the user full, unobstructed access to all data in the database.

As noted above, a bookmark may also be exported as a separate SQLite dabase with the ``--export_bookmark_db`` flag.

Finally, a receptor stored in the database may be re-exported as a PDBQT with the ``--export_receptor_pdbqt`` option.

Export results from a previous filtering as a CSV
==================================================
Filtered poses and their select information can be exported to a csv file, where the user can specify select columns to include in the csv using ``--outfields``:
.. code-block:: bash

    $ rt_process_vs write --docking_results Files/
    $ rt_process_vs read --input_db output.db --score_percentile 0.1 --bookmark_name filter1
    $ rt_process_vs read --input_db output.db --export_bookmark_csv filter1 --outfields ligname,pose_rank,docking_score,ligand_smile


Using a config file
*********************
It is possible to populate the argument list using a config file, which needs to be in a json format. The keywords needs to correspond exactly to an argument option, and the value given can be provided as a string as you would type it using the command line interface.

.. code-block:: bash

    $ rt_process_vs -c config_w.json write
    $ rt_process_vs -c config_r.json read

.. code-block:: python 

    config_w.json:
        {
        "docking_results": "path1/",
        "output_db": "example.db"
        }

    config_r.json:
        {
        "score_percentile": "0.1"
        }

The Ringtail API can provide a config file template by running the following script. The file will be saved as ``config.json``.

.. code-block:: bash

    $ rt_generate_config_file

Logging
********
Ringtail comes with a global logger object that will write to a new text file for each time ``rt_process_vs`` is called. Any log messages will also be displayed in stdout. and the default logger level is "WARNING". It is possible to change the logger level by adding ``--debug`` for lowest level of logging (will make the process take longer) or ``--verbose`` for some additional, but not very deep, logging. To write logger output to a log file, use the tag with chosen file name, e.g., ``--logfile ringtail.log`` 

.. code-block:: bash

    $ rt_process_vs write --verbose --docking_results filelist1.txt 

Access help message
********************

.. code-block:: bash

    $ rt_process_vs --help

    $ rt_process_vs write --help

    $ rt_process_vs read --help

Available command line arguments
**********************************


Keywords pertaining to database write and file handling
========================================================
.. _input_kw_table:
.. csv-table:: Ringtail input options
    :header: "Keyword","Description","Default"
    :widths: 30, 70, 10
    
    "docking_mode", "Docking engine used to perform the molecular docking","|default_docking_mode|"
    "output_db","Name of the database to which to write the docking output","|default_output_db|"
    "docking_results", "docking file(s), path(s) to files to read into database, file(s) with list of files to read into database", "|default_docking_results|"
    "recursive", "Flag to perform recursive subdirectory search on provided directory(s)", "|default_recursive|"
    "receptor_file", "Use with save_receptor and/or for vina data if interaction calculations are wanted. Give receptor JSON or PDBQT.", "|default_receptor_file|"
    "save_receptor", "Flag to specify that receptor file should be imported to database. Receptor file must also be specified with receptor_file", "|default_save_receptor|"
    "max_poses", "Number of top-scoring poses to save in database", "|default_max_poses|"
    "store_all_poses", "Flag to indicate that all poses should be stored in database", "|default_store_all_poses|"
    "interaction_tolerance", "Adds the interactions for poses within some tolerance RMSD range of the top pose in a cluster to that top pose. Can use as flag with default tolerance of 0.8, or give other value as desired [note]_ ", "|default_interaction_tolerance|"
    "no_interactions", "If interactions for AD6 or vina results should not be calculated and stored", "False"
    "interaction_cutoffs", "Use values other than defaults for distance cutoffs for measuring interactions between ligand and receptor in angstroms. Give as string, separating cutoffs for hydrogen bonds and VDW with comma (in that order). E.g. '3.7,4.0' will set the cutoff for hydrogen bonds to 3.7 angstroms and for VDW to 4.0.", "|default_interaction_cutoffs|"
    "max_proc", "Maximum number of subprocesses to spawn during database writing.", "Num available CPUs"
    "append_results", "Add new docking files to existing database given with input_db", "|default_append_results|"
    "duplicate_handling", "Specify how dulicate results should be handled. May specify 'ignore' or 'replace'. Unique results determined from ligand and target names and ligand pose. *NB: use of duplicate handling causes increase in database writing time*", "|default_duplicate_handling|"
    "overwrite", "Flag to overwrite existing database", "|default_overwrite|"
    "storage_type", "Database engine/backend to use", "|default_storage_type|"


Keywords pertaining to filtering 
=================================
.. _filter_kw_table:
.. csv-table:: Ringtail filters
    :header: "Keyword","Description","Default"
    :widths: 30, 70, 10

    "eworst","Worst energy value accepted (kcal/mol)","|default_eworst|"
    "ebest","Best energy value accepted (kcal/mol)","|default_ebest|"
    "leworst","Worst ligand efficiency value accepted","|default_leworst|"
    "lebest","Best ligand efficiency value accepted","|default_lebest|"
    "score_percentile","Worst energy percentile accepted. Giveas percentage (1 for top 1%, 0.1 for top 0.1%)","|default_score_percentile|"
    "le_percentile","Worst ligand efficiency percentile accepted. Give as percentage (1 for top 1%, 0.1 for top 0.1%)","|default_le_percentile|"
    "ligand_name","Search for specific ligand name. Multiple names joined by 'OR'. Multiple filters should be separated by commas","|default_ligand_name|"
    "ligand_name_file","Text file with ligand names, can provide thousands instead, alternative to ligand_name","|default_ligand_name_file|"
    "ligand_max_atoms","Specify maximum number of heavy atoms a ligand may have","|default_ligand_max_atoms|"
    "ligand_substruct","SMARTS pattern(s) for substructur matching","|default_ligand_substruct|"
    "ligand_substruct_pos","SMARTS pattern, index of atom in SMARTS, cutoff distance, and target xyz coordinates. Finds poses in which the specified substructure atom is within the distance cutoff from the target location","|default_ligand_substruct_pos|"
    "ligand_operator","logical operator for multiple SMARTS","|default_ligand_operator|"
    "ligand_min_molweight","Minimum molecular weight of ligands", "|default_ligand_min_molweight|"
    "ligand_max_molweight","Maximum molecular weight of ligands", "|default_ligand_max_molweight|"
    "vdw_interactions","Filter for van der Waals interaction with given receptor information. [note]_ ","|default_vdw_interactions|"
    "hb_interactions","Filter with hydrogen bonding interaction with given information. Does not distinguish between donating or accepting. [note]_ ","|default_hb_interactions|"
    "reactive_interactions","Filter for reation with residue containing specified information. [note]_ ","|default_reactive_interactions|"
    "hb_count","Filter for poses with at least this many hydrogen bonds, inclusive (5 keeps poses with 5 or more). A negative value filters for no more than that many, also inclusive (-5 keeps poses with 5 or fewer); 0 keeps only poses with no hydrogen bonds. Does not distinguish between donating and accepting. [note]_ ","|default_hb_count|"
    "react_any","Filter for poses with reaction with any residue. [note]_ ","|default_react_any|"
    "max_miss","Will filter given interaction filters excluding up to max_miss interactions. Will log and output union of combinations unless used with `enumerate_interaction_combs`. See section for reference. [note]_", "|default_max_miss|"

.. [note] Requires interactions are calculated and present in the database.


Keywords pertaining to output of data
======================================
.. _output_kw_table:
.. csv-table:: Ringtail output options
    :header: "Keyword","Description","Default"
    :widths: 30, 70, 10

    "output_log","Name for log of filtered results","|default_output_log|"
    "overwrite","Flag to overwrite existing logfile of same name","|default_overwrite|"
    "bookmark_name","Name for bookmark view in database","|default_bookmark_name|"
    "input_bookmark","Name for bookmark to use as basis for further filtering (as opposed to all results","none" 
    "outfields","Data fields to be written in output (log file and STDOUT). Ligand name always included.","|default_outfields|"
    "order_results","String for field by which the passing results should be ordered in log file.","|default_order_results|"
    "output_all_poses","Flag that if mutiple poses for same ligand pass filters, log all poses","|default_output_all_poses|"
    "mfpt_cluster","Cluster ligands passing given filters based on the Tanimoto distances of the Morgan fingerprints. Will output ligand with best (lowest) ligand efficiency from each cluster. Uses Butina clustering algorithm","|default_mfpt_cluster|"
    "interaction_cluster","Cluster ligands passing given filters based on the Tanimoto distances of the interaction fingerprints. Will output ligand with best (lowest) ligand efficiency from each cluster. Uses Butina clustering algorithm (*)","|default_interaction_cluster|"
    "enumerate_interaction_combs","When used with `max_miss` > 0, will log ligands/poses passing each separate interaction filter combination as well as union of combinations. Can significantly increase runtime. (*)","|default_enumerate_interaction_combs|"
    "individual_sdf_files","Whether or not to output one large SDF file with all ligands, or one per ligand","|default_individual_sdf_files|"

Keywords pertaining to output methods
======================================
.. _read_kw_table:
.. csv-table:: Ringtail read/output methods
    :header: "Keyword","Description","Input options"
    :widths: 10, 30, 10

    "print_summary","Prints a summary of the database, incl. number of ligands, poses, interactions, and energy percentiles",None
    "print_bookmarks","Method that prints name of all screening/filter bookmarks to stdout",None
    "data_from_bookmark","Method that makes an output log file for an existing bookmark (specified with ``--bookmark_name``)", None
    "export_bookmark_csv", "Name of database result bookmark or table to be exported as CSV. Output as <table_name>.csv. Use with ``--outfields`` to chose columns to export", "requested_data= bookmark_name OR csv_name, table (bool)"
    "export_query_csv", "Create csv of the requested SQL query. Output as query.csv. MUST BE PRE-FORMATTED IN SQL SYNTAX e.g. SELECT [columns] FROM [table] WHERE [conditions]", "requested_data (str), csv_name (str), table (bool)"
    "export_bookmark_db", "Export a database containing only the results found in the specified bookmark name. Will save as <core_db_file>_<bookmark_name>.db", "bookmark_name (str)"
    "export_receptor_pdbqt", "Export receptor to pdbqt", None
    "export_sdf_path", "Write molecule sdfs from a given bookmark to specified path", "sdf_path (str), bookmark_name (str)"
    "find_similar_ligands", "Given query ligand name, find ligands previously clustered with that ligand. User prompted at runtime to choose cluster group of interest.", "query_ligname (str)"
    "get_previous_filter_data", "Get data requested in `outfields` from the bookmark of a previous filtering", "outfields (str), bookmark_name (str)"
    "find_similar_ligands", "Find ligands in cluster with query_ligname", "query_ligname (str)"
    "logfile","File in which to write debug logging to",None


Other command-line tools
****************************
Besides ``rt_process_vs``, Ringtail installs several focused command-line utilities, each documented on its topic page:

* ``rt_compare`` — select ligands shared between, or exclusive to, the filter bookmarks of multiple screenings (cross-target comparison). See :ref:`compare`.
* ``rt_merge`` — merge two or more Ringtail databases of the same target into one. See :ref:`big_data`.
* ``rt_compress_db`` / ``rt_decompress_db`` — optionally filter, then compress a database for transfer off an HPC, and unpack it again. See :ref:`compress`.
* ``rt_upgrade_db`` — upgrade databases made with older Ringtail versions (e.g. 1.1.0, 2.0.0) to the current schema version. See :ref:`upgrade_database`.
* ``rt_generate_config_file`` — write a template JSON configuration file of ``rt_process_vs`` options, which can be passed back with ``-c`` / ``--config``.

Run any of them with ``--help`` for the full list of options.
