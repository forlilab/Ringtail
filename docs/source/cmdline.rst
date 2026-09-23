.. _cmdline:

Command line interface
#######################

The Ringtail command line interface is the easiest method to use to for exploring virtual screening results in a database. 
If this is your first time learning about Ringtail, take a look at the page :ref:`Get started <get_started>`. The current page uses the knowledge already presented on the 'Get started' page as we continue exploring the wealth of options that Ringtail offers. 

The script for writing a database and filtering is ``rt_process_vs``. It accepts
AD6 SDFs, AutoDock-GPU DLGs, or Vina PDBQTs for one target and binding site. A
database may contain multiple ligand libraries, but only one receptor. Ringtail
checks receptor names when appending results, but cannot detect different binding
sites on the same named receptor. A receptor PDBQT or Meeko Polymer JSON may be
saved in the database.

The ``rt_process_vs`` script has ``write`` and ``read`` modes. Specify the mode before
its options; a config file, when used, precedes the mode. ``write`` creates or appends
to a database, while ``read`` filters and exports an existing database.

Please note that Ringtail does not automatically have permission to perform changes outside of the working directory, so be advised that any folders or documents that Ringtail outputs will be saved in the current working directory. 

Ringtail inputs
*****************

The write examples in this section use the AutoDock-GPU and Vina test data, and are run from the ``test`` directory of the Ringtail repository, because the paths inside the example file lists are relative to it:

.. code-block:: bash

  $ cd test

Input file_sources
===================
By default, the newly-created database will be named ``output.db``. This name may be changed with the ``--output_db`` or ``-o`` option.
Ringtail allows multiple formats to provide docking results using ``--docking_results`` or ``-dr``, including paths to one or more single files, path to one or more .txt files containing docking file paths, and paths to one or more directories containing files (use ``--recursive``/``-r`` to search recursively traverse folder). You can specify one or more arguments at once:

.. code-block:: bash

    $ rt_process_vs write --docking_mode adgpu --docking_results test_data/adgpu/group3/60239.dlg.gz test_data/adgpu/group3/60373.dlg.gz test_data/adgpu/group2/ test_data/adgpu/filelist1.txt test_data/adgpu/filelist2.txt --output_db example.db

Example file list (``test_data/adgpu/filelist1.txt``). Relative paths in a file list are resolved against the current working directory, not against the location of the list:

.. code-block:: text

    test_data/adgpu/group1/127458.dlg.gz
    test_data/adgpu/group1/173101.dlg.gz
    test_data/adgpu/group1/100729.dlg.gz

Input options
===============
To include the details of a receptor in the database, it is necessary to explicitly state that the receptor should be saved. If the ``--save_receptor`` argument is invoked, but no ``--receptor_file`` (PDBQT or Meeko Polymer JSON) is provided, Ringtail will raise an error. During results processing Ringtail checks to make sure the provided receptor file matches the receptor name used for the docking.

.. code-block:: bash

    $ rt_process_vs write --docking_mode adgpu --docking_results test_data/adgpu/filelist1.txt --receptor_file test_data/adgpu/4j8m.pdbqt.gz --save_receptor

It is possible to add docking results *or* a receptor file to a database that already exists. For this it is necessary to use the keyword ``--append_results``.
Bookmarks store the poses that passed filtering, so adding results to a database that has bookmarks would leave every bookmark out of date. In that case Ringtail asks before writing: answering ``yes`` deletes all bookmarks, filters and clusterings and adds the results, any other answer writes nothing. Use ``--yes`` to skip the prompt, e.g. in batch jobs, where there is no terminal to answer it. Status assignments (Accepted/Maybe/Rejected) and pose comments are kept.
You can also specify what to do if you are adding duplicate results for a ligand, by invoking the ``--duplicate_handling`` keyword with the value ``IGNORE`` (keeps the existing entry and does not add the new duplicate) or ``REPLACE`` (overwrites the existing entry with the new duplicate; the status assignment and comment of a replaced pose are deleted with it). Please note that the ``--duplicate_handling`` option makes database writing significantly slower.

.. code-block:: bash

    $ rt_process_vs write --input_db output.db --docking_mode adgpu --docking_results test_data/adgpu/group2 --append_results --duplicate_handling REPLACE

By default, Ringtail stores up to |default_max_poses| poses per ligand: for DLGs, the best-scored (lowest energy) pose from each of the first |default_max_poses| pose clusters, and for Vina and AD6, the |default_max_poses| best poses. Use ``--max_poses`` to change this number, or ``--store_all_poses`` (which overrides ``--max_poses``) to store every pose. The number of poses is fixed when a database is created, so results appended later must use the same setting.

For AutoDock-GPU results, Ringtail reads interactions already present in the DLG files. For AD6 and Vina results, Ringtail calculates interactions by default when a receptor is supplied. Use ``--no_interactions`` to disable calculation and reduce database write time. The default distance cutoffs are 3.7 Å for hydrogen bonds and 4.0 Å for van der Waals interactions. Use ``--interaction_cutoffs`` to specify different cutoffs. Provide the receptor using ``--receptor_file`` during database writing, or use a receptor already saved in the database.

The ``--interaction_tolerance`` option also allows the user to give more leeway for poses to pass given interaction filters. With this option, the interactions from poses within *c* angstrom RMSD of a cluster's top pose will be appended to the interactions for that top pose. The theory behind this is that this gives some sense of the "fuzziness" of a given binding pose, allowing the user to filter for interactions that may not be present for the top pose specifically, but could be easily accessible to it. When used as a flag, the ``--interaction_tolerance`` default is 0.8 angstroms. The user may also specify their own cutoff. This option is intended for use with DLGs from AD-GPU, which clusters output poses based on RMSD.

It is further possible to overwrite a database by use of the argument ``--overwrite``.

.. code-block:: bash

    #AD-GPU
    $ rt_process_vs write --output_db adgpu_2poses.db --docking_mode adgpu --docking_results test_data/adgpu/group1 --max_poses 2 --interaction_tolerance 0.8

    #vina
    $ rt_process_vs write --output_db vina.db --docking_results test_data/vina --docking_mode vina --overwrite --receptor_file test_data/vina/receptor.pdbqt --save_receptor --interaction_cutoffs 3.5,4.5

Printing a database summary
***************************
During both ``write`` and ``read`` it is possible to add the tag ``-su`` or ``--print_summary`` which will print a summary of the database to stdout.

.. code-block:: bash

    $ rt_process_vs read --input_db output.db -su

Filtering
*********
In ``read`` mode, an existing database is used to filter or export results.

The read examples from here on use the AD6 database ``output.db`` created in :ref:`Get started <get_started>` (in ``test/test_data/ad6``), which contains the ligands ``first_mol``, ``second_mol``, ``third_mol``, and ``fourth_mol``.

A text results log is written only when ``--output_log`` is supplied. Used alone, it writes ``output_log.txt``, or provide a filename to choose another name. By default, it includes ligand names and docking scores for the best-scoring passing pose of each ligand.

There are six scoring filters: best and worst docking score, best and worst ligand
efficiency, and docking-score or ligand-efficiency percentiles. Do not combine
``--eworst`` with ``--score_percentile``, or ``--leworst`` with
``--le_percentile``; Ringtail rejects conflicting cutoffs.

Scoring filters
=================

.. code-block:: bash

    $ rt_process_vs read --input_db output.db --eworst -14 --output_log output_log_e14.txt

The information written to the log can be specified with ``--outfields``. The full list of available output fields may be seen by using the ``--help`` option with ``read`` mode.
By default, only the information for the top-scoring binding pose will be written to the log. If desired, each individual passing pose can be written by using the ``--output_all_poses`` flag. The passing results may also be ordered in the log file using the ``--order_results`` option.

.. code-block:: bash

    $ rt_process_vs read --input_db output.db --eworst -6 --outfields ligname,docking_score,pose_rank --order_results reference_rmsd --bookmark_name eworst6

When filtering, the passing results are saved as a bookmark in the database. This bookmark is named ``passing_results`` by default. The user can specify a name using the ``--bookmark_name`` option. No filtering is performed if no filters are given (see full list of filters :ref:`here <filter_kw_table>`).

Before Ringtail v3, bookmarks were stored as database views that evaluated the filter query when accessed. In v3, bookmarks instead store the selection of poses that passed filtering. Because a stored selection cannot include results added later, adding docking results deletes all bookmarks after asking for consent (see *Input options* above); rerun filtering afterwards.

Filtering may take from seconds to minutes, depending on the size of the database, roughly scaling as O(n) for n database Results rows (i.e. stored poses). Data for poses in a bookmark may be written to a log later using the ``--data_from_bookmark`` option together with ``--bookmark_name`` and ``--output_log``.

Interaction filters
=====================
It is possible to filter the docking results based on different types of interactions (hydrogen bonds, van der Waals, or reactive interactions) with specific residues. It is further possible to have ligands pass the filters while only fulfilling some of the interaction combinations in union (max number of interactions combinations missed, ``--max_miss``).
The available interaction filters are ``--hb_interactions``, ``--vdw_interactions``, and ``--reactive_interactions``. Interaction filters must be specified in the order ``CHAIN:RES:NUM:ATOM_NAME``. Any combination of that information may be used, as long as 3 colons are present and the information ordering between the colons is correct. All desired interactions of a given type (e.g. ``-vdw``) may be specified with a single option tag (``-vdw A:VAL:243:,A:GLU:246:``) or separate tags (``-vdw A:VAL:243: -vdw A:GLU:246:``).

The ``--max_miss`` option allows up to the specified number of interaction filters to
be absent. By default, Ringtail stores the union of the passing combinations. Use
``--enumerate_interaction_combs`` to create a separate bookmark for each combination;
this can significantly increase runtime.
``--react_any`` offers an option to filtering for poses that have reactions with any residue.

.. code-block:: bash

    $ rt_process_vs read --input_db output.db --eworst -6 --hb_interactions A:VAL:245: --hb_interactions A:GLU:246: --vdw_interactions A:VAL:243: --vdw_interactions A:GLU:246: --max_miss 1 --react_any

Ligand filters 
=================
The docked ligands can be filtered for presence of certain substructures specified by their SMARTS string using ``--ligand_substruct``, as well as their ligand name containing a specific phrase ``--ligand_name``. The ligand name search will include any ligand names that contain the specified phrase, and does not look for exact matches only.
Use the keyword ``--ligand_operator`` to determine if the ligand filters should be evaluated as this ``OR`` that (default), or combined with ``AND``. ``--ligand_max_atoms`` can be used to specify maximum number of heavy atoms a ligand may have.

.. code-block:: bash

This example selects ligands containing both a carbonyl and a pyridine ring, with at most 15 heavy atoms (``fourth_mol``):

.. code-block:: bash

    $ rt_process_vs read --input_db output.db --ligand_substruct 'C=O' 'c1ccncc1' --ligand_operator AND --ligand_max_atoms 15

The ``--ligand_substruct_pos`` option may be used to filter for a specific ligand substructure to be placed within some distance of a given cartesian coordinate. The format for this option is the six elements inside quotes and separated by spaces: ``"<SMARTS pattern: str> <index of atom in SMARTS: int> <cutoff distance: float> <target x coord: float> <target y coord: float> <target z coord: float>"``. This example selects poses with a carbonyl oxygen (atom index 1 in ``[C]=O``) within 1.5 Å of the given coordinate (``first_mol``):

.. code-block:: bash

    $ rt_process_vs read --input_db output.db --ligand_substruct_pos "[C]=O 1 1.5 7.0 -20.6 16.9"

Clustering
============
In addition to the filtering options outlined in the table below, ligands passing given filters can be clustered to provide a reduced set of dissimilar ligands based on Morgan fingerprints (``--mfpt_cluster``) or interaction (``--interaction_cluster``) fingerprints. Dissimilarity is measured by Tanimoto distance and clustering is performed with the Butina clustering algorithm. Clustering can also be performed on a bookmark that has already been saved to the database, without providing any extra filter values. In this case, the bookmark over which to cluster (or additional filtering) on is specified by ``--input_bookmark`` (must be different from ``--bookmark_name``).

.. code-block:: bash

    $ rt_process_vs read --input_db output.db --input_bookmark eworst6 --mfpt_cluster

While not quite a filtering option, the user can provide a ligand name from a previously-run clustering and re-output other ligands that were clustered with that query ligand with ``--find_similar_ligands``. The user is prompted at runtime to choose a specific clustering group from which to re-output ligands. Filtering/clustering will be performed from the same command-line call prior to this similarity search, but all subsequent output tasks will be performed on the group of similar ligands obtained with this option unless otherwise specified. 

Outputs
*********
The primary output of ``write`` mode is the database. ``read`` mode stores filtering
results as bookmarks and can optionally write logs, SDFs, CSVs, subset databases, and
receptor files.

Using ``--export_sdf_path`` saves poses from the current filtering or selected bookmark to a directory. By default, Ringtail writes one SDF containing all ligands; ``--individual_sdf_files`` writes one SDF per ligand. Poses are ordered by increasing docking score. Binding energies, ligand efficiencies, and interactions are stored as SDF properties.

If the user wishes to explore the data in CSV format, Ringtail provides two options for exporting CSVs. The first is ``--export_bookmark_csv``. Select the bookmark or table with ``--bookmark_name``. Use ``--export_bookmark_csv`` alone for ``<bookmark_name>.csv``, or supply an output filename. Use ``--outfields`` to choose columns when exporting a bookmark.
The second option is ``--export_query_csv``. This takes a string of a properly-formatted SQL query to run on the database, returning the results of that query as ``query.csv``. This option allows the user full, unobstructed access to all data in the database.

Export a bookmark as a separate database with ``--export_bookmark_db``. The exported database uses the same backend as the source database.

Finally, a receptor stored in the database may be re-exported as a PDBQT with the ``--export_receptor_pdbqt`` option.

Export results from a previous filtering as a CSV
==================================================
Filtered poses and their select information can be exported to a csv file, where the user can specify select columns to include in the csv using ``--outfields``:

.. code-block:: bash

    $ rt_process_vs read --input_db output.db --eworst -14 --bookmark_name filter1
    $ rt_process_vs read --input_db output.db --bookmark_name filter1 --export_bookmark_csv --outfields ligname,pose_rank,docking_score,ligand_smile


Using a config file
*********************
It is possible to populate the argument list using a config file, which needs to be in a json format. The keywords need to correspond exactly to an argument option, and the value given can be provided as a string as you would type it using the command line interface. ``docking_results`` accepts a single path or a list of paths. Options given on the command line take precedence over the config file.

.. code-block:: bash

    $ rt_process_vs -c config_w.json write
    $ rt_process_vs -c config_r.json read

``config_w.json``:

.. code-block:: json

    {
      "docking_results": "path1/",
      "output_db": "example.db"
    }

``config_r.json``:

.. code-block:: json

    {
      "input_db": "example.db",
      "score_percentile": 0.1
    }

The Ringtail API can provide a config file template by running the following script. The file will be saved as ``config.json``.

.. code-block:: bash

    $ rt_generate_config_file

Logging
********
Ringtail logs warnings and errors to the terminal by default. Use ``--verbose`` for
INFO messages or ``--debug`` for DEBUG messages. Supply ``--logfile ringtail.log``
to write the same logger output to a file.

.. code-block:: bash

    $ rt_process_vs write --verbose --docking_mode adgpu --docking_results test_data/adgpu/filelist1.txt --output_db verbose.db --logfile ringtail.log

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
    "recursive", "Flag to perform recursive subdirectory search on provided directory(s)", "False"
    "receptor_file", "Receptor JSON or PDBQT used to calculate AD6/Vina interactions and/or save the receptor", "|default_receptor_file|"
    "save_receptor", "Flag to specify that receptor file should be imported to database. Receptor file must also be specified with receptor_file", "False"
    "max_poses", "Number of top-scoring poses to save in database", "|default_max_poses|"
    "store_all_poses", "Flag to indicate that all poses should be stored in database", "|default_store_all_poses|"
    "interaction_tolerance", "Adds the interactions for poses within some tolerance RMSD range of the top pose in a cluster to that top pose. Can use as flag with default tolerance of 0.8, or give other value as desired [note]_ ", "|default_interaction_tolerance|"
    "no_interactions", "If interactions for AD6 or vina results should not be calculated and stored", "False"
    "interaction_cutoffs", "Use values other than defaults for distance cutoffs for measuring interactions between ligand and receptor in angstroms. Give as string, separating cutoffs for hydrogen bonds and VDW with comma (in that order). E.g. '3.7,4.0' will set the cutoff for hydrogen bonds to 3.7 angstroms and for VDW to 4.0.", "|default_interaction_cutoffs|"
    "max_proc", "Maximum number of subprocesses to spawn during database writing.", "Num available CPUs"
    "append_results", "Add new docking files to existing database given with input_db", "False"
    "yes", "If the database has bookmarks/filters, delete them (and any clusterings) without asking so results can be added", "False"
    "duplicate_handling", "Specify how duplicate results should be handled. May specify 'ignore' or 'replace'. Unique results determined from ligand and target names and ligand pose. *NB: use of duplicate handling causes increase in database writing time*", "|default_duplicate_handling|"
    "overwrite", "Flag to overwrite existing database", "False"
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
    "score_percentile","Worst energy percentile accepted. Give as percentage (1 for top 1%, 0.1 for top 0.1%)","|default_score_percentile|"
    "le_percentile","Worst ligand efficiency percentile accepted. Give as percentage (1 for top 1%, 0.1 for top 0.1%)","|default_le_percentile|"
    "ligand_name","Search for specific ligand name. Multiple names joined by 'OR'. Multiple filters should be separated by commas","|default_ligand_name|"
    "ligand_name_file","Text file with ligand names, can provide thousands instead, alternative to ligand_name","|default_ligand_name_file|"
    "ligand_max_atoms","Specify maximum number of heavy atoms a ligand may have","|default_ligand_max_atoms|"
    "ligand_substruct","SMARTS pattern(s) for substructure matching","|default_ligand_substruct|"
    "ligand_substruct_pos","SMARTS pattern, index of atom in SMARTS, cutoff distance, and target xyz coordinates. Finds poses in which the specified substructure atom is within the distance cutoff from the target location","|default_ligand_substruct_pos|"
    "ligand_operator","logical operator for multiple SMARTS","OR"
    "ligand_min_molweight","Minimum molecular weight of ligands", "|default_ligand_min_molweight|"
    "ligand_max_molweight","Maximum molecular weight of ligands", "|default_ligand_max_molweight|"
    "vdw_interactions","Filter for van der Waals interaction with given receptor information. [note]_ ","|default_vdw_interactions|"
    "hb_interactions","Filter with hydrogen bonding interaction with given information. Does not distinguish between donating or accepting. [note]_ ","|default_hb_interactions|"
    "reactive_interactions","Filter for reaction with residue containing specified information. [note]_ ","|default_reactive_interactions|"
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
    "bookmark_name","Name for a filtering bookmark or bookmark selected for export","|default_bookmark_name|"
    "input_bookmark","Name for bookmark to use as basis for further filtering (as opposed to all results)","none"
    "outfields","Data fields to be written in filtered output. Ligand name is always included. If omitted, logs contain ligand name and docking score; full-table CSV exports retain all columns.","none"
    "order_results","String for field by which the passing results should be ordered in log file.","|default_order_results|"
    "output_all_poses","Include every passing pose, rather than only the best passing pose per ligand, in requested outputs","|default_output_all_poses|"
    "mfpt_cluster","Cluster ligands passing given filters based on the Tanimoto distances of the Morgan fingerprints. Will output ligand with best (lowest) ligand efficiency from each cluster. Uses Butina clustering algorithm","|default_mfpt_cluster|"
    "interaction_cluster","Cluster ligands passing given filters based on the Tanimoto distances of the interaction fingerprints. Will output ligand with best (lowest) ligand efficiency from each cluster. Uses Butina clustering algorithm (*)","|default_interaction_cluster|"
    "enumerate_interaction_combs","When used with `max_miss` > 0, will log ligands/poses passing each separate interaction filter combination as well as union of combinations. Can significantly increase runtime. (*)","|default_enumerate_interaction_combs|"
    "individual_sdf_files","Whether or not to output one large SDF file with all ligands, or one per ligand","False"

Keywords pertaining to output methods
======================================
.. _read_kw_table:
.. csv-table:: Ringtail read/output methods
    :header: "Keyword","Description","Input options"
    :widths: 10, 30, 10

    "print_summary","Prints a summary of the database, incl. number of ligands, poses, interactions, and energy percentiles",None
    "print_bookmarks","Method that prints name of all screening/filter bookmarks to stdout",None
    "data_from_bookmark","Method that makes an output log file for an existing bookmark (specified with ``--bookmark_name``). Requires ``--output_log``", None
    "export_bookmark_csv", "Export the bookmark or table selected by ``--bookmark_name`` as CSV. Optionally provide an output filename; use ``--outfields`` to choose bookmark columns.", "optional output filename"
    "export_query_csv", "Run the supplied SQL query and write its results to query.csv", "SQL query"
    "export_bookmark_db", "Export a database containing only the results found in the specified bookmark name. Will save as <core_db_file>_<bookmark_name>.db", "flag (uses ``--bookmark_name``)"
    "export_receptor_pdbqt", "Export receptor to pdbqt", None
    "export_sdf_path", "Write molecule sdfs from a given bookmark to specified path", "sdf_path (str), bookmark_name (str)"
    "find_similar_ligands", "Given query ligand name, find ligands previously clustered with that ligand. User prompted at runtime to choose cluster group of interest.", "query_ligname (str)"
    "logfile","File in which to write debug logging to",None


Other command-line tools
****************************
Besides ``rt_process_vs``, Ringtail installs several focused command-line utilities, each documented on its topic page:

* ``rt_compare`` — select ligands shared between, or exclusive to, the filter bookmarks of multiple screenings (cross-target comparison). See :ref:`compare`.
* ``rt_merge`` — merge two or more Ringtail databases of the same target into one. See :ref:`big_data`.
* ``rt_compress_db`` / ``rt_decompress_db`` — optionally filter, then compress a database for transfer off an HPC, and unpack it again. See :ref:`compress`.
* ``rt_recalc_interactions`` — recalculate stored pose-receptor interactions, optionally with new distance cutoffs. See :ref:`recalc_interactions_cli` below.
* ``rt_upgrade_db`` — upgrade databases made with older Ringtail versions (e.g. 1.1.0, 2.0.0) to the current schema version. See :ref:`upgrade_database`.
* ``rt_generate_config_file`` — write a template JSON configuration file of ``rt_process_vs`` options, which can be passed back with ``-c`` / ``--config``.

Run any of them with ``--help`` for the full list of options.

.. _recalc_interactions_cli:

Recalculating interactions
==========================
``rt_recalc_interactions`` deletes the interactions stored in one or more databases and calculates them again from the stored poses and receptor, so no docking files are needed. The database must contain a receptor. Work is committed in batches, and a run that is interrupted resumes where it stopped when the same command is run again; resuming with different cutoffs is refused. Recalculating changes ``num_hb`` and ``num_interactions``, so bookmarks filtered on interactions are worth re-running afterwards. See :ref:`api` for the equivalent ``add_interactions`` method.

.. code-block:: bash

    $ rt_recalc_interactions -d output.db
    $ rt_recalc_interactions -d vs1.db vs2.db --hb_cutoff 3.5 --vdw_cutoff 4.5 --yes

.. csv-table:: ``rt_recalc_interactions`` options
   :header: "Argument", "Description", "Default"
   :widths: 25, 60, 15

   "``-d``, ``--database``", "One or more database files to recalculate interactions for (required)", "none"
   "``--hb_cutoff``", "Hydrogen bond distance cutoff in angstroms", "3.7"
   "``--vdw_cutoff``", "Van der Waals distance cutoff in angstroms", "4.0"
   "``--chunk_size``", "Poses per commit. Smaller means more frequent checkpoints to resume from, and less memory", "500"
   "``-y``, ``--yes``", "Skip the confirmation prompt; required for unattended and batch-scheduler runs", "False"
   "``--debug``", "Log at DEBUG level (default is INFO)", "False"
   "``--logfile``", "Write log output to this file", "none"
