
.. _api:

API procedures
###############

The Ringtail API allows for a more advanced, flexible use of ringtail where the user can create their own scripts. In the case of the docking engine AutoDock-Vina that outputs docking results as strings, Ringtail can be directly integrated in the virtual screen pipeline and add docking results to the database as soon as the docking engine outputs the results. 

The Ringtail API can thus process output files from both AutoDock-GPU (DLGs) and vina (PDBQTs), as well as vina strings. Ringtail is intended to be used for a set of docking results for a single target and binding site. This may include multiple ligand libraries as long as the target and binding site is the same. Be cautious when adding results from multiple screening runs, since some target information is checked and some is not. One receptor PDBQT may also be saved to the database.

Unlike the command line interace (:ref:`cmdline`) the API does not need to be specified for a write or read mode. It works by instantiating a RingtailCore object, and performing actions on that object. The API also offers some extra flexibility compared to the command line interface, for example it is possible to produce an internally inconsistent database that includes docking results from both AutoDock-GPU and vina, or saving docking results with different number of poses. Due to the nature of how the Ringtail API can be used, as long as a RingtailCore object has been instantiated you can keep adding results without e.g., specifying that you are appending to an existing database. 

Please note that Ringtail does not automatically have permission to perform changes outside of the working directory, so be advised that any folders or documents that Ringtail outputs will be saved in the current working directory. 

Ringtail inputs
****************

Start by creating an instance of the RingtailCore class. The object will be created with the default database file name of ``output.db`` unless otherwise specified. The logger level may also be set at this time, and defaults to ``"WARNING"``. The ``docking_mode`` defaults to ``"dlg"`` (AutoDock-GPU) and can be set at instantiation or changed at any time by directly setting the docking_mode property of the RingtailCore object. 

.. code-block:: python

    rtc = RingtailCore()
    #or
    rtc = RingtailCore(db_file = "output.db", logging_level = "DEBUG", docking_mode = "vina")
    rtc.docking_mode = "dlg"

Writing docking results to the database
=======================================
To add results from docking results files, the method ``add_results_from_files`` is used. It allows one or multiple sources of results, and a range of options pertinent to the storage handling and the results processing can be set at this time. Please note that both sources of results and settings for database writing can be provided as single options or a dictionary of allowed keywords and values. If both are provided, any individual options will overwrite that given in the dictionary. The path to the receptor file is considered part of the results options. 

.. code-block:: python

    rtc.add_results_from_files( file = ["lig1.dlg", "lig2.dlg"]
                                file_path = "path1/path2", 
                                file_list = "filelist.txt",
                                recursive = True, 
                                receptor_file = "receptor.pdbqt",
                                save_receptor = True,)
    
Example file list:

.. code-block:: python

    lig3.dlg
    lig4.dlg.gz
    rec1.pdbqt

Adding a receptor
-------------------
The receptor can be added to a database by itself, either to a populated database without a receptor, or to an empty database that you anticipate to fill with docking results later. 

.. code-block:: python
    
    rtc.save_receptor(receptor_file = "receptor.pdbqt")

Printing a database summary
---------------------------
If at any point you wish to print a summary of the contents of the database, the method can be called directly. 

.. code-block:: python

    rtc.produce_summary()

Input options
==============
The Ringtail API uses the same options that are used in the command line interface. Relevant to adding results to the database, including how many poses of a docked ligand to save, and how to handle any duplicated ligands. 

Handling of duplicate and existing results
-------------------------------------------
With the Ringtial API you can keep adding results using the same object without specifying whether or not to ``append_results``, which is contrary to the command line interface where one command line call corresponds to one ringtail core object and one connection to the database.
You can specify what to do if you are adding duplicate results for a ligand, by invoking the ``duplicate_handling`` keyword with the value ``IGNORE`` (will not add the newest duplicate) or ``REPLACE`` (will overwrite the newest duplicate). Please note that the ``duplicate_handling`` option makes database writing significantly slower.

.. code-block:: python

    rtc.add_results_from_files( file_path = "path1/",
                                duplicate_handling = "REPLACE")

Handling interaction parameters
----------------------------------
ADGPU is capable of performing interaction analysis at runtime, with these results being stored in the database if present. If interaction analysis is not present in the input file (including Vina PDBQTs), it may be added by Ringtail with the ``add_interactions`` option. **This adds a signifcant increase to the total database write time.** Distance cutoffs for the interactions are specified with the ``interaction_cutoffs`` option. Adding interactions requires that the receptor has already been added to the database, or by supplying the receptor PDBQT as one of the inputs.

.. code-block:: python

    rtc.docking_mode = "vina"
    rtc.add_results_from_files( file = "lig1.pdbqt"
                                add_interactions = True,
                                receptor_file = "receptor.pdbqt",
                                save_receptor = True,
                                interaction_cutoffs = [3.7, 4.0])

The ``interaction_tolerance`` option also allows the user to give more leeway for poses to pass given interaction filters. With this option, the interactions from poses within *c* angstrom RMSD of a cluster's top pose will be appended to the interactions for that top pose. The theory behind this is that this gives some sense of the "fuzziness" of a given binding pose, allowing the user to filter for interactions that may not be present for the top pose specifically, but could be easily accessible to it. When used as a flag, the ``interaction_tolerance`` default is 0.8 angstroms. The user may also specify their own cutoff. This option is intended for use with DLGs from AD-GPU, which clusters output poses based on RMSD.

.. code-block:: python

    rtc.docking_mode = "dlg"
    rtc.add_results_from_files( file_path = "path1/",
                                duplicate_handling = "REPLACE",
                                interaction_tolerance = 0.6)

Number of poses to save
-------------------------
By default (for DLGs), Ringtail will store the best-scored (lowest energy) binding pose from the first 3 pose clusters in the DLG. For Vina, Ringtail will store the 3 best poses. Additional settings for writing to the database include how to handle the number of poses docked (``max_poses``, or ``store_all_poses`` which will overwrite the former).

.. code-block:: python

    rtc.add_results_from_files( file_path = "path2"
                                max_poses = 5)

Filtering
**********

Docking results stored in the Ringtail database can be filtered using the ``filter`` method. When filtering, a text log file will be created containing the results passing the given filter(s). The default log name is ``output_log.txt`` and by default will include the ligand name (``Ligand_Name``) and docking score (``e``) of every pose passing filtering criteria. The name of the filter log name may be changed using the ``log_file`` keyword. 

Scoring filters
=================
There are six scoring filters that include best (``ebest``) and worst docking score/energy (``eworst``), best and worst ligand efficieny (``lebest`` and ``leworst``), and results above worst docking score or ligand efficiency percentile (``score_percentile`` and ``le_percentile``, respecitvely). Some of these are internally inconsistent: if both ``eworst`` and ``score_percentile`` are used together, the ``eworst`` cutoff alone is used. The same is true of ``leworst`` and ``le_percentile``.

.. code-block:: python

    rtc.filter(score_percentile = 0.1, log_file = "output_log_01percent.txt")

The information written to the log file can be specified with ``outfields``. The full list of available output fields may be seen in the documentation/"hover-over" over the method.
By default, only the information for the top-scoring binding pose will be written to the log. If desired, each individual passing pose can be written by using ``output_all_poses = True``. The passing results may also be ordered in the log file using the ``order_results`` option.

.. code-block:: python

    rtc.filter(eworst = -6, outfields = "Ligand_Name,e,rank,receptor", order_results = "ref_rmsd", bookmark_name = "eworst6")

When filtering, the passing results are also saved as a view (or bookmark) in the database. This view is named ``passing_results`` by default. The user can specify a name for the view with the ``bookmark_name`` keyword. No filtering is performed if no filters are given (see full list of filters #REF). 
Filtering may take from seconds to minutes, depending on the size of the database, roughly scaling as O(n) for n database Results rows (i.e. stored poses). Data for poses in a view may be accessed later using the ``get_previous_filter_data`` method.

.. code-block:: python

    rtc.get_previous_filter_data(outfields = "Ligand_Name,e,rank", bookmark_name = "eworst6", log_file = "previously_filtered_results.txt")

Interaction filters
=====================
It is possible to filter the docking results based on different types of interactions (hydrogen bonds and van der waals interactions) with specific residues. It is further possible to have ligands pass the filters while only fulfilling some of the interaction combinations in union (max number of interactions combinations missed, ``max_miss``).
The available interaction filters are ``hb_interactions``, ``vdw_interactions``, and ``reactive_interactions``. Interaction filters must be specified as the interaction specifications in the order ``CHAIN:RES:NUM:ATOM_NAME``. Any combination of that information may be used, as long as 3 colons are present and the information ordering between the colons is correct. All desired interactions of a given type is specified as a list of one or more tuples of specified reactions and weather to show results that includes ``(":::", True)`` or exclude ``(":::", False)`` them as shown below for ``vdw_interactions``:

.. code-block:: python

    rtc.filter( eworst=-2,
                vdw_interactions=[('A:VAL:279:', True), ('A:LYS:162:', True)])

The ``max_miss`` keywords allows the user to filter by given interactions excluding up to ``max_miss`` interactions. This gives :math:`\sum_{m=0}^{m}\frac{n!}{(n-m)!*m!}` combinations for *n* interaction filters and *m* max_miss. By default, results will be given for the union of the interaction conbinations. Use with ``enumerate_interaction_combs = True`` to log ligands/poses passing each separate interaction combination (can significantly increase runtime). If ``max_miss > 0`` is used during filtering, a view is created for each combination of interaction filters and is named ``<bookmark_name>_<n>`` where n is the index of the filter combination in the log file (indexing from 0).
``react_any`` offers an option to filtering for poses that have reactions with any residue.

.. code-block:: python

    rtc.filter( eworst=-6,
                vdw_interactions=[('A:VAL:279:', True), ('A:LYS:162:', True)],
                hb_interactions = [("A:VAL:279:", True), ("A:LYS:162:)", True)],
                max_miss = 1,
                react_any = True)


Ligand filters 
===============
Several filters pertaining to the SMARTS structure of the ligand can be used. For example, the ``ligand_substruct_pos`` keyword may be used to filter for a specific ligand substructure (specified with a SMARTS string) to be placed within some distance of a given cartesian coordinate. The format for this option is ``"<SMARTS pattern: str>" <index of atom in SMARTS: int> <cutoff distance: float> <target x coord: float> <target y coord: float> <target z coord: float>``.
ligand_name: Specify ligand name(s). Will combine name filters with 'OR'.
ligand_substruct: SMARTS pattern(s) for substructure matching.
ligand_substruct_pos: SMARTS pattern(s) for substructure matching, e.g., ['[Oh]C', 0, 1.2, -5.5, 10.0, 15.5] -> [smart_string, index_of_positioned_atom, cutoff_distance, x, y, z].
ligand_max_atoms: Maximum number of heavy atoms a ligand may have.
ligand_operator: Logical join operator for multiple SMARTS.

.. code-block:: python

    rtc.filter(ligand_substruct=["[Oh]C"], ligand_substruct_pos=["[Oh]C", 0, 1.2, -5.5, 10.0, 15.5])


Clustering
============
In addition to the filtering options outlined in the table below, ligands passing given filters can be clustered to provide a reduced set of dissimilar ligands based on Morgan fingerprints (``mfpt_cluster``) or interaction (``interaction_cluster``) fingerprints. Dissimilarity is measured by Tanimoto distance (float input to the cluster keyword) and clustering is performed with the Butina clustering algorithm. Clustering can be also be performed on a bookmark that has already been saved to the database, without providing any extra filter values. In this case, the bookmark over which to cluster (or additional filtering) on is specified by ``filter_bookmark`` (must be different from ``bookmark_name`` that contains previously filtered results).

.. code-block:: python

    rtc.filter( filter_bookmark = "eworst6",
                mfpt_cluster = 0.6)

While not quite a filtering option, the user can provide a ligand name from a previously-run clustering and re-output other ligands that were clustered with that query ligand with the method ``find_similar_ligands``. The user is prompted at runtime to choose a specific clustering group from which to re-output ligands. Filtering/clustering will be performed from the same command-line call prior to this similarity search, but all subsequent output tasks will be performed on the group of similar ligands obtained with this option unless otherwise specified. 

.. code-block:: python

    rtc.find_similar_ligands("ligand_name")


Output options
***************
There are multiple options to output and visualize data in Ringtail.

Creating scatter plot of ligands
=================================
The method ``plot`` generates a scatterplot of ligand efficiency vs docking score for the top-scoring pose from each ligand. Ligands passing the given filters or in the bookmark given with the keyword ``bookmark_name`` will be highlighted in red. The plot also includes histograms of the ligand efficiencies and binding energies. The plot is saved as ``scatter.png``.

.. code-block:: python

    rtc.plot()

Evaluating molecules in 3d space with PyMol
=============================================
The ``pymol`` flag generates a scatterplot of ligand efficiency vs docking score as well, but only for the ligands contained in the bookmark specified with ``bookmark_name``. It also launches a PyMol session and will display the ligands in PyMol when clicked on the scatterplot. N.B.: Some users may encounter a ``ConnectionRefusedError``. If this happens, try manually launching PyMol (``pymol -R``) in a separate terminal window.

.. code-block:: python

    rtc.pymol(bookmark_name = "eworst6")

Export molecule SDF files
==========================
The method ``write_molecule_sdfs`` will write SDF files for each ligand passing the filter and saved in a specified bookmark (can also include those who don't pass by invoking the ``write_nonpassing = True`` option). The files will be saved to the path specified in the method call. If none is specified, the files will be saved in the current working directory. The SDF will contain poses passing the filter/in the bookmark ordered by increasing docking score. Each ligand is written to its own SDF. This option enables the visualization of docking results, and includes any flexible/covalent ligands from the docking. The binding energies, ligand efficiencies, and interactions are also written as properties within the SDF file, with the order corresponding to the order of the pose order.

.. code-block:: python

    rtc.write_molecule_sdfs(sdf_path = "sdf_files", bookmark_name = "eworst6")

Exporting tables or query results as CSV files
==============================================
If the user wishes to explore the data in CSV format, Ringtail provides two options for exporting CSVs. First, you can export a database table or bookmark (``requested_data``) to a csv file with a name (``csv_name``) specified in the method call. In this case one must specify that the type of the ``requested_data`` is of database type table. 

.. code-block:: python
    
    rtc.export_csv(requested_data = "Ligands", csv_name = "Ligand_table.csv", table = True)

It is also possible to write a database query and export the results of the query to a csv file. In this case, the requested data must be a properly formatted SQL query string. User needs to specify that the ``requested_data`` is not provided directly as a table. 

.. code-block:: python

    query_string = "SELECT docking_score, leff, Pose_ID, LigName FROM Results"
    rtc.export_csv(requested_data = query_string, csv_name = "query_results.csv", table = False)

Creating a new database from a bookmark
=======================================
A bookmark may also be exported as a separate SQLite dabase with the ``export_bookmark_db`` method. This will create a database of name ``<current_db_name>_<bookmark_name>.db``. This is currently only possible if using SQLite.

.. code-block:: python 

    rtc.export_bookmark_db(bookmark_name = "eworst6")

    #results in creation of output_eworst6.db

Finally, a receptor stored in the database may be re-exported as a PDBQT with the ``export_receptor`` method. This will save the receptor PDBQT in the current working directory. 

.. code-block:: python 

    rtc.export_bookmark_db()

Some usage notes
******************

For many of these operations, if you do not specify a bookmark name Ringtail will simply use the bookmark that was last used for operations in the object. If it is a newly instantiated object without a "recent" bookmark, it will look for a bookmark of the default name ``passing_results``. 
Most methods accept both individual options as well as grouped options in a dictionary format. In each of these cases, for arguments that are duplicated between the two formats individual options will overwrite that given in the dictionary. 

Logging module
===============
Ringtail comes with a global logger object that will write to a new text file for each time a new ``RingtailCore()`` object is instantiated. Any log messages at or above specified level will write to the log file as well as be displayed in stdout . The default logger level is "WARNING". While logger level can be set at time of instantiation, it is also possible to change it later by accessing the logger object directly:

.. code-block:: python

    rtc = RingtailCore("database.db")
    # log level defaults to "WARNING"
    logger.setLevel("INFO")

Available Ringtail arguments
******************************

Keywords pertaining to database write and file handling
========================================================

.. csv-table:: Ringtail input options
    :header: "Keyword","Description","Default value"
    :widths: 30, 70, 10

    "file", "DLG/Vina PDBQT file(s) to be read into database", None
    "file_path", "Path(s) to files to read into database", None
    "file_list", "File(s) with list of files to read into database", None
    "pattern", "Specify pattern to search for when finding files", "'dlg' or 'pdbqt'"
    "recursive", "Flag to perform recursive subdirectory search on file_path directory(s)", FALSE
    "receptor_file", "Use with save_receptor and/or add_interactions. Give receptor PDBQT.", None
    "save_receptor", "Flag to specify that receptor file should be imported to database. Receptor file must also be specified with receptor_file", FALSE
    "max_poses", "Number of clusters for which to store top-scoring pose (dlg) or number of poses (vina) to save in database", 3
    "store_all_poses", "Flag to indicate that all poses should be stored in database", FALSE
    "interaction_tolerance", "Adds the interactions for poses within some tolerance RMSD range of the top pose in a cluster to that top pose. Can use as flag with default tolerance of 0.8, or give other value as desired [note]_ ", "0.8 Å if used"
    "add_interactions", "Find interactions between ligands and receptor. Requires receptor PDBQT to be written.", FALSE
    "interaction_cutoffs", "Specify distance cutoffs for measuring interactions between ligand and receptor in angstroms. Give as string, separating cutoffs for hydrogen bonds and VDW with comma (in that order). E.g. '3.7,4.0' will set the cutoff for hydrogen bonds to 3.7 angstroms and for VDW to 4.0.", "3.7,4.0"
    "max_proc", "Maximum number of subprocesses to spawn during database writing.", "number of available CPUs or fewer"
    "append_results", "Add new docking files to existing database given with input_db", FALSE
    "duplicate_handling", "Specify how dulicate results should be handled. May specify 'ignore' or 'replace'. Unique results determined from ligand and target names and ligand pose. *NB: use of duplicate handling causes increase in database writing time*", None
    "overwrite", "Flag to overwrite existing database", FALSE


Keywords pertaining to filtering and read/output
================================================

.. csv-table:: Ringtail filters
    :header: "Keyword","Description","Default value"
    :widths: 30, 70, 10

    "eworst","Worst energy value accepted (kcal/mol)",None
    "ebest","Best energy value accepted (kcal/mol)",None
    "leworst","Worst ligand efficiency value accepted",None
    "lebest","Best ligand efficiency value accepted",None
    "score_percentile","Worst energy percentile accepted. Giveas percentage (1 for top 1%, 0.1 for top 0.1%)",1.0
    "le_percentile","Worst ligand efficiency percentile accepted. Give as percentage (1 for top 1%, 0.1 for top 0.1%)",None
    "ligand_name","Search for specific ligand name. Multiple names joined by 'OR'. Multiple filters should be separated by commas",None
    "ligand_max_atoms","Specify maximum number of heavy atoms a ligand may have",None
    "ligand_substruct","SMARTS pattern(s) for substructur matching",None
    "ligand_substruct_pos","SMARTS pattern, index of atom in SMARTS, cutoff distance, and target xyz coordinates. Finds poses in which the specified substructure atom is within the distance cutoff from the target location",None
    "ligand_operator","logical operator for multiple SMARTS","OR"
    "vdw_interactions","Filter for van der Waals interaction with given receptor information. [note]_ ",None
    "hb_interactions","Filter with hydrogen bonding interaction with given information. Does not distinguish between donating or accepting. [note]_ ",None
    "reactive_interactions","Filter for reation with residue containing specified information. [note]_ ",None
    "hb_count","Filter for poses with at least this many hydrogen bonds. Does not distinguish between donating and accepting. [note]_ ",None
    "react_any","Filter for poses with reaction with any residue. [note]_ ",FALSE
    "max_miss","Will filter given interaction filters excluding up to max_miss interactions. Will log and output union of combinations unless used with `enumerate_interaction_combs`. See section for reference. [note]_", 0

.. [note] Requires interactions are calculated and present in the database.

Keywords pertaining to output of data
======================================
.. csv-table:: Ringtail output options
    :header: "Keyword","Description","Default value"
    :widths: 30, 70, 10

    "log_file","Name for log of filtered results","output_log.txt"
    "overwrite","Flag to overwrite existing logfile of same name",FALSE
    "bookmark_name","Name for bookmark view in database","passing_results"
    "outfields","Data fields to be written in output (log file and STDOUT). Ligand name always included.","e"
    "order_results","String for field by which the passing results should be ordered in log file.",None
    "output_all_poses","Flag that if mutiple poses for same ligand pass filters, log all poses",FALSE
    "mfpt_cluster","Cluster ligands passing given filters based on the Tanimoto distances of the Morgan fingerprints. Will output ligand with best (lowest) ligand efficiency from each cluster. Uses Butina clustering algorithm",0.5
    "interaction_cluster","Cluster ligands passing given filters based on the Tanimoto distances of the interaction fingerprints. Will output ligand with best (lowest) ligand efficiency from each cluster. Uses Butina clustering algorithm (*)",0.5
    "enumerate_interactions_combs","When used with ``max_miss`` > 0, will log ligands/poses passing each separate interaction filter combination as well as union of combinations. Can significantly increase runtime. (*)",FALSE

