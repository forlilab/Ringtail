.. _compare:

Compare docking results from different virtual screenings
##########################################################

The script ``rt_compare`` is designed to be used with databases already made and filtered. It is used to combine information across multiple virtual screenings to allow or exclude the selection of ligands passing filters across multiple targets/models. This can be useful for filtering out promiscuous ligands, a technique commonly used in experimental high-throughput screening. It may also be used if selection of ligands binding multiple protein structures/conformations/homologs are desired.

Programmatically, the ``rt_compare`` script is used to select ligands which are shared between the given filter bookmark(s) of some virtual screenings (``--wanted``) or exclusive to some screenings and not others (``--unwanted``). The script uses a subset of commands similar to ``rt_process_vs``.

The basic process is illustrated by a target, ``kinase1``, its related proteins
``kinase1a`` and ``kinase1b``, and an unrelated protein, ``protein2``:

#. Create one database for each target.
#. Filter each database and note the bookmark used for each selection. Filters may
   differ between targets, for example to select analogous receptor interactions.
#. Use ``rt_compare`` to find ligands that pass the wanted bookmarks and exclude
   those that pass unwanted bookmarks. Supply every database/bookmark pair with a
   separate ``--wanted`` or ``--unwanted`` option.
#. Choose any required outputs:

   * ``-l`` or ``-l comparison.txt`` writes text results in the same format as
     ``rt_process_vs`` filter output, one file per database, written next to that
     database and named ``<database name>_<log name>``
     (``<database name>_crossref_results.txt`` when ``-l`` is used without a filename).
   * ``-xs`` / ``--export_sdf`` writes the selected ligands and poses to one SDF per
     database, ``<database name>_<bookmark>.sdf``, in the current working directory.
   * ``-xd`` / ``--export_db`` exports each cross-referenced bookmark as a database,
     ``<database name>_<bookmark>.db``, next to the source database.

The comparison creates a bookmark in each participating database. Its name combines
the ``--save_bookmark`` prefix (``crossref`` by default) with the source bookmark name.

.. code-block:: bash

    $ rt_compare --wanted kinase1.db best_interactions --unwanted kinase1a.db bad_interaction_filter --unwanted kinase1b.db bad_interaction_filter

The same command can select potential dual-target ligands:


.. code-block:: bash

    $ rt_compare --wanted kinase1.db passing_results --wanted protein2.db top99percentile --unwanted kinase1a.db bad_interaction_filter --unwanted kinase1b.db bad_interaction_filter


Usage examples
****************

Select ligands found in "passing_results" bookmarks of vs1 but not vs2 or vs3
===============================================================================

.. code-block:: bash

    $ rt_compare --wanted vs1.db passing_results --unwanted vs2.db passing_results --unwanted vs3.db passing_results

Select ligands found in "passing_results" bookmarks of vs1 and vs2 but not vs3 or vs4
======================================================================================

.. code-block:: bash

    $ rt_compare --wanted vs1.db passing_results -w vs2.db passing_results --unwanted vs3.db passing_results -uw vs4.db passing_results

Select ligands found in "passing_results" bookmarks of every vs except vs4
============================================================================

.. code-block:: bash

    $ rt_compare -w vs1.db passing_results -w vs2.db passing_results -w vs3.db passing_results -uw vs4.db passing_results

Select ligands found in "filter1" bookmarks of vs1 but not "passing_results" of vs2
===================================================================================

.. code-block:: bash

    $ rt_compare -w vs1.db filter1 -uw vs2.db passing_results

Save bookmark of ligands found in "filter1" bookmarks of vs1 and vs2 but not vs3 or vs4 as "selective_bookmark_filter1" in all databases
========================================================================================================================================

.. code-block:: bash

    $ rt_compare -w vs1.db filter1 -w vs2.db filter1 -uw vs3.db filter1 -uw vs4.db filter1 --save_bookmark selective_bookmark

Export bookmark set of ligands found in "filter1" bookmarks of vs1 and vs2 but not vs3 or vs4 as SDFs per database
==================================================================================================================

.. code-block:: bash

    $ rt_compare -w vs1.db filter1 -w vs2.db filter1 -uw vs3.db filter1 -uw vs4.db filter1 --export_sdf

Access help message for rt_compare
**********************************

.. code-block:: bash

    $ rt_compare --help


Supported arguments for the comparison script
***********************************************

.. csv-table:: ``rt_compare`` options
   :header: "Argument", "Short", "Description", "Default"
   :widths: 20, 8, 52, 20

   ``--config``, ``-c``, "JSON configuration file; command-line values take precedence", "none"
   ``--wanted``, ``-w``, "Database and bookmark to include", "none"
   ``--unwanted``, ``-uw``, "Database and bookmark to exclude", "none"
   ``--store_best_pose``, ``-bp``, "Write only the best-ranked pose per ligand to the text results log", ``False``
   ``--output_log``, ``-l``, "Write text results; optionally supply a filename", "none"
   ``--save_bookmark``, ``-s``, "Prefix for bookmarks created by the comparison", ``crossref``
   ``--export_db``, ``-xd``, "Export each compared bookmark as a database", "disabled"
   ``--export_sdf``, ``-xs``, "Export each compared bookmark as SDF (a value given after the flag is ignored)", "disabled"
   ``--verbose``, ``-v``, "Set log level to INFO", "disabled"
   ``--debug``, ``-d``, "Set log level to DEBUG", "disabled"
   ``--logfile``, "", "Write logger output to a file", "none"



The ``--store_best_pose`` flag does not restrict the saved bookmarks or SDF exports.
