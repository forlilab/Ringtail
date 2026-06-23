.. _get_started:

Getting started
###############

The Ringtail command line interface is orchestrated through the script ``rt_process_vs``.

Create and populate a database
*********************************
Navigate to the directory containing the data, in our case test_data from Autodock-GPU (make sure your Ringtail environment is active):

.. code-block:: bash

    $ cd test/test_data/adgpu

To write to the database we need to specify a few things:
- that we are using ``write`` mode
- source of docking results files. Docking results can be added either by providing one or more single files, a .txt file containing files, or by providing a directory containing docking results files (or a combination of al of them).
- optional database name: ringtail will default to creating a database of name |default_output_db|
- optional docking mode: ringtail will default to assuming the files were produced by |default_docking_mode|, otherwise specify with ``--docking_mode <mode>``

Let us add all docking files within the path test_data (specified by ``.`` meaning current directory), whose folders we can traverse recursively by specifying ``--recursive``

.. code-block:: bash

    $ rt_process_vs write --docking_results . --recursive

We can print a summary of the contents of the database by using the optional tag ``-su`` or ``--summary`` and specifying the database database from which to ``read``:

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

Filtering and visualizing the data in the database
***************************************************

Let us start filtering with a basic docking score cutoff of -6 kcal/mol:

.. code-block:: bash

    $ rt_process_vs read --input_db output.db --eworst -6 --output_log

By using the argument ``-l`` / ``--output_log`` ringtail will write an output log with the default name ``output_log.txt`` with the names of ligands passing the filter, as well as their binding energies. It is not necessary to use this argument unless you want a written log of the filter results (you can also specify a different name for the file). Each round of filtering is also stored in the database as a ringtail bookmark (default value is ``passing_results``). The bookmark name is the basis for further progressive filtering as well as making use of export options such as creating SDF files. 

We can choose a specific bookmark name for any round of filtering. For example, start out with filtering out the compounds that are within the 5th percentile in terms of docking score and save the bookmark as `ep5`:

.. code-block:: bash

    $ rt_process_vs read --input_db output.db --score_percentile 5 --bookmark_name ep5

Let's then further refine the set of molecules by applying an interaction filter for van der Waals interactions with valine 279 on the receptor:

.. code-block:: bash

    $ rt_process_vs read --input_db output.db --input_bookmark ep5 --vdw_interactions A:VAL:279: --bookmark_name ep5_vdwV279

Once you are happy with the number of passing ligands, you can for example produce an SDF file with all passing ligands and their poses:

.. code-block:: bash

    $ rt_process_vs read --input_db output.db --bookmark_name ep5_vdwV279 -sdf /path/to/sdf_folder/

Ringtail offers a large number of properties on which to filter and screen the docking data, and many formats for exports suitable for collaboration and further screening. Explore e.g., the command line page for more options (:ref:`cmdline`). 


Access help message for rt_process_vs
**************************************

.. code-block:: bash

    $ rt_process_vs --help

Access options message for rt_process_vs write mode
***************************************************

.. code-block:: bash

    $ rt_process_vs write --help

Access options message for rt_process_vs read mode
**************************************************

.. code-block:: bash

    $ rt_process_vs read --help

