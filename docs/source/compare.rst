.. _compare:

Compare docking results from different virtual screenings
##########################################################


The script ``rt_compare`` is designed to be used with databases already made and filtered. It is used to combine information across multiple virtual screenings to allow or exclude the selection of ligands passing filters across multiple targets/models. This can be useful for filtering out promiscuous ligands, a technique commonly used in exerimental high-throughput screening. It may also be used if selection of ligands binding multiple protein structures/conformations/homologs are desired.

Programmatically, the ``rt_compare`` script is used to select ligands which are shared between the given filter bookmark(s) of some virtual screenings (``--wanted``) or exclusive to some screenings and not others (``--unwanted``). The script uses a subset of commands similar to ``rt_process_vs``.

The basic process of preparing to use this script and the concept behind it is thus:

Let us assume that kinase1 is our target of interest. It has related proteins kinase1a and kinase1b. protein2 is an unrelated protein.
1. Create a database for each virtual screening on each target (kinase1.db, kinase1a.db, kinase1b.db, protein2.db)
2. Filter each database separately to get a set of virtual hits for each target. Each set of filters may be different as desired (e.g. change interaction filters for analogous residues), and make not of the bookmark name used in each database, here let's say kinase1.db filtering was stored in bookmark ``best_interactions``, while the filters in kinase1a.db and kinase1b.db was the same ``bad_interaction_filter``. (Each databaser and bookmark needs a separate ``--(un)wanted`` keyword.)
3. Use ``rt_compare`` to find ligands that pass the filters for kinase1 but not kinase1a or kinase1b. This will produce a new bookmark in each database with the original specified bookmark name prefixed the option provuded from ``--save_bookmark`` (defaults to ``crossref``).
4. There are several output options, each will produce the same output per compared database as poses, energies, et c., will be different. 
    a. By specifying ``--log`` or ``--log comparison.txt`` the script will create a log file of the same format as the filter output from ``rt_process_vs``. 
    b. Specifying ``-xs`` or ``--export_sdf`` will create one file with all ligands and poses that passed
    c. Specifying ``-xd`` or ``--export_database`` will create a new database from the crossreferenced bookmark
    d. Other export options maybe used after the cross referencing has been performed, either by connecting to one of the Ringtail databases, and exporting based on the crossreferenced bookmark, or by performing Ringtail actions on an exported filtered database.  

.. code-block:: bash

    $ rt_compare --wanted kinase1.db best_interactions --unwanted kinase1a.db bad_interaction_filter --unwanted kinase1b.db bad_interaction_filter

4. Other usage examples and output options given below. For example, one can also select for potential dual-target ligands with


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

+----------------+----+----------------------------------------------------------------------------------------------------------------------------------+------------------+
| Argument            | Description                                                                                                                      | Default value    |
+================+====+==================================================================================================================================+==================+
|--config        | -c | Configuration JSON file to specify new default options. Overridded by command line                                               | no default       |
+----------------+----+----------------------------------------------------------------------------------------------------------------------------------+------------------+
|--wanted        | -w | Database files and associated bookmark names for which to include the intersection of ligands                                    | no default       |
+----------------+----+----------------------------------------------------------------------------------------------------------------------------------+------------------+
|--unwanted      | -uw| Database files and associated bookmark names for which to exclude the intersection of ligands                                    | no default       |
+----------------+----+----------------------------------------------------------------------------------------------------------------------------------+------------------+
|--log           | -l | Name for log file to which results are written                                                                                   | no default       |
+----------------+----+----------------------------------------------------------------------------------------------------------------------------------+------------------+
|--save_bookmark | -s | Prefix used with original bookmark name, as the bookmark that will be stored in each database after comparison                   | no default       |
+----------------+----+----------------------------------------------------------------------------------------------------------------------------------+------------------+
|--export_db     | -xd| Export the new compared bookmark in each database as a new database                                                              | no default       |
+----------------+----+----------------------------------------------------------------------------------------------------------------------------------+------------------+
|--export_sdf    | -xs| Export the new compared bookmark in each database as an SD file of all ligands and poses in that database                        | no default       |
+----------------+----+----------------------------------------------------------------------------------------------------------------------------------+------------------+



