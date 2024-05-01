.. _compare:

Compare docking results from different virtual screenings
##########################################################


The script ``rt_compare.py`` is designed to be used with databases already made and filtered. It is used to combine information across multiple virtual screenings to allow or exclude the selection of ligands passing filters across multiple targets/models. This can be useful for filtering out promiscuous ligands, a technique commonly used in exerimental high-throughput screening. It may also be used if selection of ligands binding multiple protein structures/conformations/homologs are desired.

Programmatically, the ``rt_compare.py`` script is used to select ligands which are shared between the given filter bookmark(s) of some virtual screenings (``--wanted``) or exclusive to some screenings and not others (``--unwanted``). The script uses a subset of commands similar to ``rt_process_vs.py``.

The basic process of preparing to use this script and the concept behind it is thus:

Let us assume that kinase1 is our target of interest. It has related proteins kinase1a and kinase1b. protein2 is an unrelated protein.
1. Create a database for each virtual screening on each target (kinase1.db, kinase1a.db, kinase1b.db, protein2.db)
2. Filter each database separately to get a set of virtual hits for each target. Each set of filters may be different as desired (e.g. change interaction filters for analogous residues). The bookmark within each database may be given as a single string (same bookmark name in every database) or multiple bookmark names (one per database) with the ``--bookmark_name`` option. If specifying multiple names, the order should match the order that the databases were provided in, beginning with wanted, then unwanted databases. The default name is ``passing_results``.
3. Use ``rt_compare.py`` to find ligands that pass the filters for kinase1 but not kinase1a or kinase1b. This will create a log file of the same format as that output from ``rt_process_vs.py``.

.. code-block:: bash

    $ rt_compare.py --wanted kinase1.db --unwanted kinase1a.db kinase1b.db

4. Other usage examples and output options given below. For example, one can also select for potential dual-target ligands with


.. code-block:: bash

    $ rt_compare.py --wanted kinase1.db protein2.db --unwanted kinase1a.db kinase1b.db


Usage examples
**********************

Select ligands found in "passing_results" bookmarks of vs1 but not vs2 or vs3
============================

.. code-block:: bash

    $ rt_compare.py --wanted vs1.db --unwanted vs2.db vs3.db

Select ligands found in "passing_results" bookmarks of vs1 and vs2 but not vs3 or vs4
============================

.. code-block:: bash

    $ rt_compare.py --wanted vs1.db vs2.db --unwanted vs3.db vs4.db

Select ligands found in "passing_results" bookmarks of every vs except vs4
============================

.. code-block:: bash

    $ rt_compare.py --wanted vs1.db vs2.db vs3.db --unwanted vs4.db

Select ligands found in "filter1" bookmarks of vs1 but not vs2
============================

.. code-block:: bash

    $ rt_compare.py --wanted vs1.db --unwanted vs2.db --bookmark_name filter1

Save bookmark of ligands found in "filter1" bookmarks of vs1 and vs2 but not vs3 or vs4 as "selective_bookmark" in vs1.db
============================

.. code-block:: bash

    $ rt_compare.py --wanted vs1.db vs2.db --unwanted vs3.db vs4.db --save_bookmark selective_bookmark

Export bookmark set of ligands found in "filter1" bookmarks of vs1 and vs2 but not vs3 or vs4 as CSV
============================
.. code-block:: bash

    $ rt_compare.py --wanted vs1.db vs2.db --unwanted vs3.db vs4.db --export_csv

Access help message for rt_compare.py
****************************

.. code-block:: bash

    $ rt_compare.py --help


Supported arguments for the comparison script
**************************

+----------------+---+----------------------------------------------------------------------------------------------------------------------------------+------------------+
| Argument           | Description                                                                                                                      | Default value    |
+================+===+==================================================================================================================================+==================+
|--config        | -c| Configuration JSON file to specify new default options. Overridded by command line                                               | no default       |
+----------------+---+----------------------------------------------------------------------------------------------------------------------------------+------------------+
|--wanted        | -w| Database files for which to include the intersection of ligands in bookmark_name(s) for all databases specified with this option.| no default       |
+----------------+---+----------------------------------------------------------------------------------------------------------------------------------+------------------+
|--unwanted      | -n| Database files for which to exclude any ligands found in bookmark_name of any of the databases specified with this option.       | no default       |
+----------------+---+----------------------------------------------------------------------------------------------------------------------------------+------------------+
|--bookmark_name |-sn| Name of bookmark to select ligands within. Must be present in all databases given.                                               | passing_results  |
+----------------+---+----------------------------------------------------------------------------------------------------------------------------------+------------------+
|--log           | -l| Name for log file to which results are written                                                                                   | selective_log.txt|
+----------------+---+----------------------------------------------------------------------------------------------------------------------------------+------------------+
|--save_bookmark | -s| Save the final selective bookmark as a view with given name in the first database specified with ``--wanted``.                   | no default       |
+----------------+---+----------------------------------------------------------------------------------------------------------------------------------+------------------+
|--export_csv    | -x| Save final selective bookmark as csv. Saved as [save_bookmark].csv or 'crossref.csv' if ``--save_bookmark`` not used.            | FALSE            |
+----------------+---+----------------------------------------------------------------------------------------------------------------------------------+------------------+



