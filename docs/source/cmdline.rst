.. _cmdline:

Ringtail command line interface
###############################

The Ringtail command line interface is the easiest method to use to for exploring virtual screening results in a database. 
If this is your first time learning about Ringtail, take a look at the page :ref:`Get started <get_started>`. The current page uses the knowledge already presented on the 'Get started' page as we continue exploring the wealth of options that Ringtail offers. 

The script for writing a database and filtering is ``rt_process_vs.py``. This is intended to be used for a set of DLGs/Vina PDBQTs pertaining to a single target and binding site. This may include multiple ligand libraries as long as the target and binding site is the same. Be cautious when adding results from multiple screening runs, since some target information is checked and some is not. One receptor PDBQT may also be saved to the database.

The rt_process_vs.py script has two modes: ``write`` and ``read``. The desired mode must be specified in the command line before any other options are given (except ``-c [CONFIG]`` which is given first). The ``write`` mode is used to create a database for a virtual screening from ADGPU DLGs or Vina PDBQTs. After this initial run, a database is created and may be read directly by rt_process_vs.py in ``read`` mode for subsequent filtering and export operations.

Please note that Ringtail does not automatically have permission to perform changes outside of the working directory, so be advised that any folders or documents that Ringtail outputs will be saved in the current working directory. 

Ringtail inputs
**************************

Navigate to the data repository and chose one of several paths of adding results. 

.. code-block:: bash

  $ cd test/test_data/

Input file_sources
=======================================
By default, the newly-created database will be named ``output.db``. This name may be changed with the ``--output_db``` option.
Ringtail allows referencing docking results files in multiple ways:
``--file/-f``: path to a single file such as ``group1/1451.dlg`` (files compressed with gzip such as ``group1/1451.dlg.gz`` are also allowed)
``--file_path/-fp``: path to a folder containing results. If this folder contains additional folders, use the ``--recursive`` option to traverse all the folders within the path. Adding "/" after the path will also make it recursive.
``--file_list/-fl``: a text file containing a list of paths to docking results files. A receptor file can also be part of this list. 
For each of these options you can specify one or more arguments, and we can create a database using all the allowed input options:

.. code-block:: bash

    $ rt_process_vs.py write --file lig1.dlg lig2.dlg --file_path path1/ path2 --file_list filelist1.txt filelist2.txt --output_db example.db

Example file list:

.. code-block:: python

    lig3.dlg
    lig4.dlg.gz
    rec1.pdbqt

Input options
=======================================
To include the details of a receptor in the database, it is necessary to provide explicitly state that the receptor should be saved. If the ``--save_receptor`` argument is invoked, but no PDBQT ``--receptor_file`` is provided, Ringtail will raise an error. During results processing Ringtail checks to make sure the provided receptor file matches the receptor name used for the docking. 

.. code-block:: bash

    $ python ../scripts/rt_process_vs.py write --file_list filelist1.txt --receptor_file test_data/4j8m.pdbqt.gz --save_receptor

It is possible to add docking results *or* a receptor file to a database that already exists. For this it is necessary to use the keyword ``--append_results``.
You can also specify what to do if you are adding duplicate results for a ligand, by invoking the ``--duplicate_handling`` keyword with the value ``IGNORE`` (will not add the newest duplicate) or ``REPLACE`` (will overwrite the newest duplicate). Please note that the ``--duplicate_handling`` option makes database writing significantly slower.

.. code-block:: bash

    $ python ../scripts/rt_process_vs.py write --input_db output.db --file_path test_data/group2 --append_results --duplicate_handling REPLACE

By default (for DLGs), Ringtail will store the best-scored (lowest energy) binding pose from the first 3 pose clusters in the DLG. For Vina, Ringtail will store the 3 best poses. Additional settings for writing to the database include how to handle the number of poses docked (``--max_poses``, or ``--store_all_poses`` which will overwrite the former).

ADGPU is capable of performing interaction analysis at runtime, with these results being stored in the database if present. If interaction analysis is not present in the input file (including Vina PDBQTs), it may be added by Ringtail with the ``--add_interactions`` option. **This adds a signifcant increase to the total database write time.** Distance cutoffs for the interactions are specified with the ``--interaction_cutoffs`` option. Adding interactions requires that the receptor PDBQT be provided as an input by the user with the ``--receptor_file`` option.

The ``--interaction_tolerance`` option also allows the user to give more leeway for poses to pass given interaction filters. With this option, the interactions from poses within *c* angstrom RMSD of a cluster's top pose will be appended to the interactions for that top pose. The theory behind this is that this gives some sense of the "fuzziness" of a given binding pose, allowing the user to filter for interactions that may not be present for the top pose specifically, but could be easily accessible to it. When used as a flag, the ``--interaction_tolerance`` default is 0.8 angstroms. The user may also specify their own cutoff. This option is intended for use with DLGs from AD-GPU, which clusters output poses based on RMSD.

It is further possible to overwrite a database by use of the argument ``--overwrite``.

.. code-block:: bash

    #AD-GPU
    $ python ../scripts/rt_process_vs.py write --input_db output.db --file_path test_data/group1 --max_poses 2 --interaction_tolerance 0.8

    #vina
    $ python ../scripts/rt_process_vs.py write --input_db output.db --file_path test_data/vina --overwrite --receptor_file receptor.pdbqt --save_receptor --add_interactions --interaction_cutoffs 3.7,4.0

Printing a database summary
*****************************
During both ``write`` and ``read`` it is possible to add the tag ``-su`` or ``--summary`` which will print a summary of the database to stdout.

.. code-block:: bash

    $ rt_process_vs.py read --input_db output.db -su

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
************************
In ``read`` mode, an existing database is used to filter or export results.

When filtering, a text log file will be created containing the results passing the given filter(s). The default log name is ``output_log.txt`` and by default will include the ligand name and docking score of every pose passing filtering criteria. The log name may be changed with the ``--log_file`` option. There are six scoring filters that include best and worst docking score/energy, best and worst ligand efficieny, and results above worst docking score or ligand efficiency percentile. Some of these are internally inconsistent: if both ``--eworst`` and ``--score_percentile`` are used together, the ``--eworst`` cutoff alone is used. The same is true of ``--leworst`` and ``--le_percentile``.

Scoring filters
===================

.. code-block:: bash

    $ python ../scripts/rt_process_vs.py read --input_db output.db --score_percentile 0.1 --log_file output_log_01percent.txt

The information written to the log can be specified with ``--outfields``. The full list of available output fields may be seen by using the ``--help`` option with ``read`` mode.
By default, only the information for the top-scoring binding pose will be written to the log. If desired, each individual passing pose can be written by using the ``--output_all_poses`` flag. The passing results may also be ordered in the log file using the ``--order_results`` option.

.. code-block:: bash

    $ python ../scripts/rt_process_vs.py read --input_db output.db --eworst -6 --outfields Ligand_Name,e,rank,receptor --order_results ref_rmsd --bookmark_name eworst6

When filtering, the passing results are also saved as a view in the database. This view is named ``passing_results`` by default. The user can specify a name for the view using the ``--bookmark_name`` option. No filtering is performed if no filters are given (see full list of filters #TODO). 
Filtering may take from seconds to minutes, depending on the size of the database, roughly scaling as O(n) for n database Results rows (i.e. stored poses). Data for poses in a view may be accessed later using the ``--data_from_bookmark`` option.

Interaction filters
=====================
It is possible to filter the docking results based on different types of interactions (hydrogen bonds onr van der waals) with specific residues. It is further possible to have ligands pass the filters while only fulfilling some of the interaction combinations in union (max number of interactions combinations missed, ``--max_miss``).
The available interaction filters are ``--hb_interactions``, ``--vdw_interactions``, and ``--reactive_interactions``. Interaction filters must be specified in the order ``CHAIN:RES:NUM:ATOM_NAME``. Any combination of that information may be used, as long as 3 colons are present and the information ordering between the colons is correct. All desired interactions of a given type (e.g. ``-vdw``) may be specified with a single option tag (``-vdw B:THR:276:,B:HIS:226:``) or separate tags (``-vdw B:THR:276: -vdw B:HIS:226:``).

The ``--max_miss`` option allows the user to filter by given interactions excluding up to ``max_miss`` interactions. This gives :math:`\sum_{m=0}^{m}\frac{n!}{(n-m)!*m!}` combinations for *n* interaction filters and *m* max_miss. By default, results will be given for the union of the interaction conbinations. Use with ``--enumerate_interaction_combs`` to log ligands/poses passing each separate interaction combination (can significantly increase runtime). If ``max_miss > 0`` is used during filtering, a view is created for each combination of interaction filters and is named ``<bookmark_name>_<n>`` where n is the index of the filter combination in the log file (indexing from 0).
``--react_any`` offers an option to filtering for poses that have reactions with any residue.

.. code-block:: bash

    $ python ../scripts/rt_process_vs.py read --input_db output.db --eworst -6 --hb_interactions A:VAL:279: A:LYS:162: --vdw_interactions A:VAL:279: A:LYS:162: --max_miss 1 --react_any)

Ligand filters #TODO copy from cmdline docu
=====================
The ``--smarts_idxyz`` option may be used to filter for a specific ligand substructure (specified with a SMARTS string) to be placed within some distance of a given cartesian coordinate. The format for this option is ``"<SMARTS pattern: str>" <index of atom in SMARTS: int> <cutoff distance: float> <target x coord: float> <target y coord: float> <target z coord: float>``.

.. code-block:: bash

    $ python ../scripts/rt_process_vs.py read --input_db output.db --eworst -6 --hb_interactions A:VAL:279: A:LYS:162: --vdw_interactions 'A:VAL:279: A:LYS:162: --max_miss 1)


Clustering
==============
In addition to the filtering options outlined in the table below #TODO, ligands passing given filters can be clustered to provide a reduced set of dissimilar ligands based on Morgan fingerprints (``--mfpt_cluster``) or interaction (``--interaction_cluster``) fingerprints. Dissimilarity is measured by Tanimoto distance and clustering is performed with the Butina clustering algorithm. Clustering can be also be performed on a bookmark that has already been saved to the database, without providing any extra filter values. In this case, the bookmark over which to cluster (or additional filtering) on is specified by ``--filter_bookmark`` (must be different from ``--bookmark_name``).

.. code-block:: bash

    $ python ../scripts/rt_process_vs.py read --input_db output.db --filter_bookmark eworst6 --mfpt_cluster

While not quite a filtering option, the user can provide a ligand name from a previously-run clustering and re-output other ligands that were clustered with that query ligand with ``--find_similar_ligands``. The user is prompted at runtime to choose a specific clustering group from which to re-output ligands. Filtering/clustering will be performed from the same command-line call prior to this similarity search, but all subsequent output tasks will be performed on the group of similar ligands obtained with this option unless otherwise specified. 

Outputs
********************
The primary outputs from ``rt_process_vs.py`` are the database itself (``write`` mode) and the filtering log file (``read`` mode). There are several other output options as well, intended to allow the user to further explore the data from a virtual screening.

The ``--plot`` flag generates a scatterplot of ligand efficiency vs docking score for the top-scoring pose from each ligand. Ligands passing the given filters or in the bookmark given with ``--bookmark_name`` will be highlighted in red. The plot also includes histograms of the ligand efficiencies and binding energies. The plot is saved as ``scatter.png``.

The ``--pymol`` flag also generates a scatterplot of ligand efficiency vs docking score, but only for the ligands contained in the bookmark specified with ``--bookmark_name``. It also launches a PyMol session and will display the ligands in PyMol when clicked on the scatterplot. N.B.: Some users may encounter a ``ConnectionRefusedError``. If this happens, try manually launching PyMol (``pymol -R``) in a separate terminal window.

Using the ``--export_sdf_path`` option allows the user to specify a directory to save SDF files for ligands passing the given filters or in the bookmark given with ``--bookmark_name``. The SDF will contain poses passing the filter/in the bookmark ordered by increasing docking score. Each ligand is written to its own SDF. This option enables the visualization of docking results, and includes any flexible/covalent ligands from the docking. The binding energies, ligand efficiencies, and interactions are also written as properties within the SDF file, with the order corresponding to the order of the pose order.

If the user wishes to explore the data in CSV format, Ringtail provides two options for exporting CSVs. The first is ``--export_bookmark_csv``, which takes a string for the name of a table or result bookmark in the database and returns the CSV of the data in that table. The file will be saved as ``<table_name>.csv``.
The second option is ``--export_query_csv``. This takes a string of a properly-formatted SQL query to run on the database, returning the results of that query as ``query.csv``. This option allows the user full, unobstructed access to all data in the database.

As noted above, a bookmark may also be exported as a separate SQLite dabase with the ``--export_bookmark_db`` flag.

Finally, a receptor stored in the database may be re-exported as a PDBQT with the ``--export_receptor`` option.

Export results from a previous filtering as a CSV
==============================================

.. code-block:: bash

    $ rt_process_vs.py write --file_path Files/
    $ rt_process_vs.py read --input_db output.db --score_percentile 0.1 --bookmark_name filter1
    $ rt_process_vs.py read --input_db output.db --export_bookmark_csv filter1

Create scatterplot highlighting ligands passing filters
==============================================

.. code-block:: bash

    $ rt_process_vs.py write --file_path Files/
    $ rt_process_vs.py read --input_db output.db --score_percentile 0.1 --bookmark_name filter1
    $ rt_process_vs.py read --input_db output.db --bookmark_name filter1 --plot

    `all_ligands_scatter.png`

.. image:: https://user-images.githubusercontent.com/41704502/215909808-2edc29e9-ebdb-4f0e-a87a-a1c293687b2e.png

Using a config file
**************************
It is possible to populate the argument list using a config file, which needs to be in a json format. The keywords needs to correspond exactly to an argument option, and the value given can be provided as a string as you would type it using the command line interface.

.. code-block:: bash

    $ rt_process_vs.py -c config_w.json write
    $ rt_process_vs.py -c config_r.json read

.. code-block:: python 

    config_w.json:
        {
        "file_path": "path1/",
        "output_db": "example.db"
        }

    config_r.json:
        {
        "score_percentile": "0.1"
        }

Logging
********************
Ringtail comes with a global logger object that will write to a new text file for each time ``rt_process_vs.py`` is called. Any log messages will also be displayed in stdout. and the default logger level is "WARNING". It is possible to change the logger level by adding ``--debug`` for lowest level of logging (will make the process take longer) or ``--verbose`` for some additional, but not very deep, logging. 

.. code-block:: bash

    $ python ../scripts/rt_process_vs.py write --verbose --file_list filelist1.txt 

Access help message
**********************

.. code-block:: bash

    $ rt_process_vs.py --help

    $ rt_process_vs.py write --help

    $ rt_process_vs.py read --help

Available command line arguments
**********************************

#TODO table showing all arguments, keywords, defaults, and info
#TODO should be a separate file containing this table that is referenced, maybe even incliude code so that it renders correctly for api vs cmd? like two columns that are conditional or something
