Getting Started
================

Ringtail offers a wealth of database creation and filtering options. They are detailed at length below. This section will provide a quick overview of the basic usage of Ringtail from the command line. We will you the provided test data to create a database with default storage options and perform basic filtering of it.

Let us begin in the Ringtail directory. First, we must navigate to the test data directory:

.. code-block:: bash

  $ cd test/test_data/

Now, let us create a database containing the results from only group 1. Note that these files are DLGs. If we were using Vina PDBQTs, we would need to add :code:`--mode vina`.

.. code-block:: bash

    $ rt_process_vs.py write --file_path group1

By default, the database we have created is called :code:`output.db`. Let us now make a second database named :code:`all_groups.db` containing all three groups:

.. code-block:: bash

    $ rt_process_vs.py write --file_path . --recursive --output_db all_groups.db

The :code:`--recursive` option tells Ringtail to scan the directories specified with :code:`--file_path` for subdirectories containing output files (in this case, DLGs). This allowed all three group directories to be added to the database with a single --file_path option.

Now that we have created the databases, we can filter them to pull out compounds of interest. Before we do that, let's find out a little more about the data contained within the database. For this, we can use the :code:`-s/--summary` option:

.. code-block:: bash

    $ rt_process_vs.py read --input_db all_groups.db -s

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

We could also have used the :code:`--summary` option when writing the database to display this info at that time.

Now, let us start filtering with a basic docking score cutoff of -6 kcal/mol:

.. code-block:: bash

    $ rt_process_vs.py read --input_db all_groups.db --eworst -6


This produces an output log :code:`output_log.txt` with the names of ligands passing the filter, as well as their binding energies. Let's now do another round of filtering, this time with an energy percentile filter of 5 percent (top 5% of coumpounds by docking score). Each round of filtering is also stored in the database as a SQLite view, which we refer to as a "bookmark". We will also save this round of filtering with the bookmark name "ep5".

.. code-block:: bash

    $ rt_process_vs.py read --input_db all_groups.db --score_percentile 5 --log ep5_log.txt --bookmark_name ep5

Now, let us further refine the set of molecules we just filtered. We will use an interaction filter for van der Waals interactions with V279 on the receptor:

.. code-block:: bash

    $ rt_process_vs.py read --input_db all_groups.db --filter_bookmark ep5 --van_der_waals A:VAL:279: --log ep5_vdwV279_log.txt --bookmark_name ep5_vdwV279


We are now ready to export these molecules for visual inspection in your favorite molecular graphics program. The :code:`export_sdf_path` will create the directory :code:`ep5_vdwV279_sdfs` in cd if it does not already existy, and store the exported molecule files there.

.. code-block:: bash
    
    $ rt_process_vs.py read --input_db all_groups.db --bookmark_name ep5_vdwV279 --export_sdf_path ep5_vdwV279_sdfs


Now we have our filtered molecules as SDF files ready for visual inspection!

## Example Filtering Timings (M1Pro MacBook, ~2 million ligands)

.. image:: https://github.com/forlilab/Ringtail/assets/41704502/eac373fc-1324-45df-b845-6697dc9d1465





# Extended documentation

## Scripts
The Ringtail package includes two command line oriented scripts: :code:`rt_process_vs.py` and :code:`rt_compare.py`. Both may be run with options specified in the command line and/or using options specified in a JSON-formatted file given with :code:`--config`. Command line options override any conflicting options in the config file.

:code:`rt_process_vs.py` serves as the primary script for the package and is used to both write docking files to a SQLite database and to perform filtering and export tasks on the database. It is designed to handle docking output files associated with a single virtual screening in a single database.

:doc:`rt_compare.py <compare>` is used to combine information across multiple virtual screenings (in separate databases) to allow or exclude the selection of ligands passing filters across multiple targets/models. This can be useful for filtering out promiscuous ligands, a technique commonly used in exerimental high-throughput screening. It may also be used if selection of ligands binding multiple protein structures/conformations/homologs are desired.

## rt_process_vs.py Documentation
### Usage examples
#### Access help message for rt_process_vs.py

.. code-block:: bash

    $ rt_process_vs.py --help

#### Access help message for rt_process_vs.py write mode

.. code-block:: bash

    $ rt_process_vs.py write --help

#### Access help message for rt_process_vs.py read mode

.. code-block:: bash

    $ rt_process_vs.py read --help

#### Create database named example.db from all input options

.. code-block:: bash

    $ rt_process_vs.py write --file lig1.dlg lig2.dlg --file_path path1/ path2 --file_list filelist1.txt filelist2.txt --output_db example.db

Example file list

.. code-block:: python

    lig3.dlg
    lig4.dlg.gz
    rec1.pdbqt

#### Write and filter using a config file

.. code-block:: bash

    $ rt_process_vs.py -c config_w.json write
    $ rt_process_vs.py -c config_r.json read

config_w.json:

.. code-block:: python 

    {
    "file_path": "path1/",
    "output_db": "example.db"
    }


config_r.json:

.. code-block:: python

    {
    "score_percentile": "0.1"
    }


#### Export results from a previous filtering as a CSV

.. code-block:: bash

    $ rt_process_vs.py write --file_path Files/
    $ rt_process_vs.py read --input_db output.db --score_percentile 0.1 --bookmark_name filter1
    $ rt_process_vs.py read --input_db output.db --export_bookmark_csv filter1

#### Create scatterplot highlighting ligands passing filters

.. code-block:: bash

    $ rt_process_vs.py write --file_path Files/
    $ rt_process_vs.py read --input_db output.db --score_percentile 0.1 --bookmark_name filter1
    $ rt_process_vs.py read --input_db output.db --bookmark_name filter1 --plot

    `all_ligands_scatter.png`

![all_ligands_scatter](https://user-images.githubusercontent.com/41704502/215909808-2edc29e9-ebdb-4f0e-a87a-a1c293687b2e.png)

### Usage Details
The script for writing a database and filtering is :code:`rt_process_vs.py`. __This is intended to be used for a set of DLGs/Vina PDBQTs pertaining to a single target and binding site. This may include multiple ligand libraries as long as the target and binding site is the same. Be cautious when adding results from multiple screening runs, since some target information is checked and some is not.__ One receptor PDBQT may also be saved to the database.

The rt_process_vs.py script has two modes: :code:`write` and :code:`read`. The desired mode must be specified in the command line before any other options are given (except :code:`-c [CONFIG]` which is given first). The :code:`write` mode is used to create a database for a virtual screening from ADGPU DLGs or Vina PDBQTs. After this initial run, a database is created and may be read directly by rt_process_vs.py in :code:`read` mode for subsequent filtering and export operations.

#### Write Mode
Upon calling rt_process_vs.py in :code:`write` mode for the first time, the user must specify where the program can find files to write to the newly-created database. This is done using the
:code:`--file`, :code:`--file_path`, and/or :code:`--file_list` options. Any combination of these options can be used, and multiple arguments for each are accepted. Compressed :code:`.gz` files
are also accepted.

When searching for result files in the directory specified with :code:`--file_path`, rt_process_vs.py will search for files with the pattern :code:`*.dlg*` by default. This may be changed with the :code:`--pattern` option. Note also that, by default, Ringtail will only search the directory provided in :code:`--file_path` and not subdirectories. Subdirectory searching
is enabled with the :code:`--recursive` flag. If you are trying to read Vina PDBQTs, specify this with :code:`--mode vina`. This will automatically change the file search pattern to :code:`*.pdbqt*`. If the receptor PDBQT file is present in a directory being searched, it **must** be specified with :code:`--receptor_file`.

To add new files to an existing database, the :code:`--append_results` flag can be used in conjuction with :code:`--input_db` and :code:`--file`, :code:`--file_path`, and/or :code:`--file_list` options. If one is concerned about adding duplicate results, the :code:`--duplicate_handling` option can be used to specify how duplicate entries should be handled. However, this option makes database writing significantly slower.

To overwrite an existing database, use the :code:`--overwrite` flag.

One receptor PDBQT, corresponding to that in the DLGs, may be saved to the database using the :code:`--save_receptor` flag. This will store the receptor file itself in a binary format in the database. The user must specify the path to the receptor file with the :code:`--receptor_file` option. Ringtail will also throw an exception if this flag is given but no receptor is found, if the name of the receptor in any DLG does not match the receptor file, or if this flag is used with a database that already has a receptor. :code:`--save_receptor` can be used to add a receptor to an existing database given with :code:`--input_db`. :code:`--save_receptor` may not be used with the :code:`--append_results` option.

By default, the newly-created database will be named :code:`output.db`. This name may be changed with the :code:`--output_db` option.

By default (for DLGs), Ringtail will store the best-scored (lowest energy) binding pose from the first 3 pose clusters in the DLG. For Vina, Ringtail will store the 3 best poses. The number of clusters/poses stored may be
changed with the :code:`--max_poses` option. The :code:`--store_all_poses` flag may also be used to override :code:`--max_poses` and store every pose from every file.

ADGPU is capable of performing interaction analysis at runtime, with these results being stored in the database if present. If interaction analysis is not present in the input file (including Vina PDBQTs), it may be added by Ringtail with the :code:`--add_interactions` option. **This adds a signifcant increase to the total database write time.** Distance cutoffs for the interactions are specified with the :code:`--interaction_cutoffs` option. Adding interactions requires that the receptor PDBQT be provided as an input by the user with the :code:`--receptor_file` option.

The :code:`--interaction_tolerance` option also allows the user to give more leeway for poses to pass given interaction filters. With this option, the interactions from poses within *c* angstrom RMSD of a cluster's top pose will be appended to the interactions for that top pose. The theory behind this is that this gives some sense of the "fuzziness" of a given binding pose, allowing the user to filter for interactions that may not be present for the top pose specifically, but could be easily accessible to it. When used as a flag, the :code:`interaction_tolerance` default is 0.8 angstroms. The user may also specify their own cutoff. This option is intended for use with DLGs from AD-GPU, which clusters output poses based on RMSD.

#### Read mode
In :code:`read` mode, an existing database is used to filter or export results.

When filtering, a text log file will be created containing the results passing the given filter(s). The default log name is :code:`output_log.txt` and by default will include the ligand name and docking score of every pose passing filtering criteria. The log name
may be changed with the :code:`--log` option and the information written to the log can be specified with :code:`--outfields`. The full list of available output fields may be seen by using the :code:`--help` option with :code:`read` mode (see example above).
By default, only the information for the top-scoring binding pose will be written to the log. If desired, each individual passing pose can be written by using the :code:`--output_all_poses` flag. The passing results may also be ordered in the log file using the :code:`--order_results` option.

No filtering is performed if no filters are given. If both :code:`--eworst` and :code:`--score_percentile` are used together, the :code:`--eworst` cutoff alone is used. The same is true of :code:`--leworst` and :code:`--le_percentile`.

In addition to the filtering options outlined in the table below, ligands passing given filters can be clustered to provide a reduced set of dissimilar ligands based on Morgan fingerprints (`--mfpt_cluster`) or interaction (`--interaction_cluster`) fingerprints. Dissimilarity is measured by Tanimoto distance and clustering is performed with the Butina clustering algorithm.

When filtering, the passing results are saved as a view in the database. This view is named :code:`passing_results` by default. The user can specify a name for the view using the :code:`--bookmark_name` option. Data for poses in a view may be accessed later using the :code:`--data_from_bookmark` option. When :code:`max_miss > 0` is used, a view is created for each combination of interaction filters and is named :code:`<bookmark_name>_<n>` where n is the index of the filter combination in the log file (indexing from 0).

Filtering may take from seconds to minutes, depending on the size of the database, roughly scaling as O(n) for n database Results rows (i.e. stored poses). One may also filter over a previous bookmark specified with the :code:`--filter_bookmark` option. If using this option, the bookmarks specified by :code:`--filter_bookmark` and :code:`--bookmark_name` must be different.

While not quite a filtering option, the user can provide a ligand name from a previously-run clustering and re-output other ligands that were clustered with that query ligand with :code:`--find_similar_ligands`. The user is prompted at runtime to choose a specific clustering group from which to re-output ligands. Filtering/clustering will be performed from the same command-line call prior to this similarity search, but all subsequent output tasks will be performed on the group of similar ligands obtained with this option unless otherwise specified. 

##### Other available outputs
The primary outputs from :code:`rt_process_vs.py` are the database itself (:code:`write` mode) and the filtering log file (:code:`read` mode). There are several other output options as well, intended to allow the user to further explore the data from a virtual screening.

The :code:`--plot` flag generates a scatterplot of ligand efficiency vs docking score for the top-scoring pose from each ligand. Ligands passing the given filters or in the bookmark given with :code:`--bookmark_name` will be highlighted in red. The plot also includes histograms of the ligand efficiencies and binding energies. The plot is saved as :code:`[filters_file].png` if a :code:`--filters_file` is used, otherwise it is saved as :code:`out.png`.

The :code:`--pymol` flag also generates a scatterplot of ligand efficiency vs docking score, but only for the ligands contained in the bookmark specified with :code:`--bookmark_name`. It also launches a PyMol session and will display the ligands in PyMol when clicked on the scatterplot. N.B.: Some users may encounter a :code:`ConnectionRefusedError`. If this happens, try manually launching PyMol (`pymol -R`) in a separate terminal window.

Using the :code:`--export_sdf_path` option allows the user to specify a directory to save SDF files for ligands passing the given filters or in the bookmark given with :code:`--bookmark_name`. The SDF will contain poses passing the filter/in the bookmark ordered by increasing docking score. Each ligand is written to its own SDF. This option enables the visualization of docking results, and includes any flexible/covalent ligands from the docking. The binding energies, ligand efficiencies, and interactions are also written as properties within the SDF file, with the order corresponding to the order of the pose order.

If the user wishes to explore the data in CSV format, Ringtail provides two options for exporting CSVs. The first is :code:`--export_bookmark_csv`, which takes a string for the name of a table or result bookmark in the database and returns the CSV of the data in that table. The file will be saved as :code:`<table_name>.csv`.
The second option is :code:`--export_query_csv`. This takes a string of a properly-formatted SQL query to run on the database, returning the results of that query as :code:`query.csv`. This option allows the user full, unobstructed access to all data in the database.

As noted above, a bookmark may also be exported as a separate SQLite dabase with the :code:`--export_bookmark_db` flag.

Finally, a receptor stored in the database may be re-exported as a PDBQT with the :code:`--export_receptor` option.

### Interaction filter formatting and options

**Interaction filtering requires interactions to be present in database.**

The :code:`--vdw`, :code:`--hb`, and :code:`--react_res` interaction filters must be specified in the order :code:`CHAIN:RES:NUM:ATOM_NAME`. Any combination of that information may be used, as long as 3 colons are present and the information ordering between the colons is correct. All desired interactions of a given type (e.g. :code:`--vdw`) may be specified with a single option tag (`--vdw=B:THR:276:,B:HIS:226:`) or separate tags (`--vdw=B:THR:276: --vdw=B:HIS:226:`).

The :code:`--max_miss` option allows the user to filter by given interactions excluding up to :code:`max_miss` interactions. This gives ![equation](https://latex.codecogs.com/svg.image?\sum_{m=0}^{m}\frac{n!}{(n-m)!*m!}) combinations for *n* interaction filters and *m* max_miss. By default, results will be given for the union of the interaction conbinations. Use with :code:`--enumerate_interaction_combs` to log ligands/poses passing each separate interaction combination (can significantly increase runtime).

The :code:`--smarts_idxyz` option may be used to filter for a specific ligand substructure (specified with a SMARTS string) to be placed within some distance of a given cartesian coordinate. The format for this option is :code:`"<SMARTS pattern: str>" <index of atom in SMARTS: int> <cutoff distance: float> <target x coord: float> <target y coord: float> <target z coord: float>`.
