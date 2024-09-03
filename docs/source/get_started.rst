.. _get_started:

Getting started with Ringtail using the command line interface
###############################################################

The Ringtail command line interface is orchestrated through the script ``rt_process_vs.py``.

Create and populate a database
*********************************
Navigate to the directory containing the data, in our case test_data:

.. code-block:: bash

    $ cd test/test_data/

To write to the database we need to specify a few things:
- that we are using ``write`` mode
- source of docking results files. Docking results can be added either by providing one or more single files, a .txt file containing files, or by providing a directory containing docking results files.
- optional database name: ringtail will default to creating a database of name ``output.db``
- optional docking mode: ringtail will default to assuming the files were produced by Autodock-GPU, if they are from vina specify ``--mode vina``

Let us add all docking files within the path test_data (specified by ``.`` meaning current directory), whose folders we can traverse recursively by specifying ``--recursive``

.. code-block:: bash

    $ rt_process_vs.py write --file_path . --recursive

We can print a summary of the contents of the database by using the optional tag ``-su`` or ``--summary`` and specifying the database database from which to ``read``:

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

Filtering and visualizing the data in the database
***************************************************

Let us start filtering with a basic docking score cutoff of -6 kcal/mol:

.. code-block:: bash

    $ rt_process_vs.py read --input_db output.db --eworst -6

This produces an output log ``output_log.txt`` with the names of ligands passing the filter, as well as their binding energies. Each round of filtering is also stored in the database as a SQLite view, which we refer to as a "bookmark" (default value is ``passing_results``). 

We can also save a round of filtering with a specific bookmark name, and perform more filtering on this bookmark.
For example, start out with filtering out the compounds that are within the 5th percentile in terms of docking score and save the bookmark as `ep5`:

.. code-block:: bash

    $ rt_process_vs.py read --input_db output.db --score_percentile 5 --log ep5_log.txt --bookmark_name ep5

Let's then further refine the set of molecules by applying an interaction filter for van der Waals interactions with V279 on the receptor:

.. code-block:: bash

    $ rt_process_vs.py read --input_db output.db --filter_bookmark ep5 --vdw_interactions A:VAL:279: --log ep5_vdwV279_log.txt --bookmark_name ep5_vdwV279

The filtered molecules can then be exported as an e.g., SDF file which can be used for visual inspection in molecular graphics programs. At the same time, if pymol is installed, we can kick off a pymol session of the ligands

.. code-block:: bash

    $ rt_process_vs.py read --input_db output.db --bookmark_name ep5_vdwV279 --export_sdf_path ep5_vdwV279_sdfs --pymol

Access help message for rt_process_vs.py
*****************************************

.. code-block:: bash

    $ rt_process_vs.py --help

Access help message for rt_process_vs.py write mode
***************************************************

.. code-block:: bash

    $ rt_process_vs.py write --help

Access help message for rt_process_vs.py read mode
**************************************************

.. code-block:: bash

    $ rt_process_vs.py read --help

