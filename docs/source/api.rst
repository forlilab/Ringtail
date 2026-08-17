
.. _api:

API procedures
###############

The Ringtail API allows for a more advanced and flexible use of Ringtail where the user can create their own scripts. In the case of the docking engines that can output docking results directly as a variable, instead of writing to the file system, including vina and AD6, Ringtail can be directly integrated in the virtual screen pipeline and write docking results to the database directly. 

The Ringtail API can thus process output files from both AD6 (SDFs), AutoDock-GPU (DLGs), and vina (PDBQTs), as well as AD6 dicts and vina strings. Ringtail is intended to be used for a set of docking results for a single target and binding site. This may include multiple ligand libraries as long as the target and binding site is the same as Ringtail only allows one receptor per database (checks are performed based on provided receptor name). However, if you do multiple screens with e.g., different binding pockets for the same receptor, this information will not be checked, and Ringtail will allow any data to be appended to a database as long as the receptor name is the same. One receptor JSON (or PDBQT) may also be saved to the database. Ringtail has tools for performing comparisons of ligands with different targets (:ref:`compare`).

Unlike the command line interface (:ref:`cmdline`) the API does not need to be specified for a write or read mode. It works by instantiating a RingtailCore object, and performing actions on that object. The API also offers some extra flexibility compared to the command line interface, for example it is possible to produce an internally inconsistent database that includes saving docking results with different number of poses. Due to the nature of how the Ringtail API can be used, as long as a RingtailCore object has been instantiated you can keep adding results without e.g., specifying that you are appending to an existing database. 

Ringtail inputs
****************

Start by creating an instance of the RingtailCore class. The object will be created with the default database file name of |default_output_db| unless otherwise specified. This is also the time to specify what database backend to use, if different from the default (|default_storage_type|).

.. code-block:: python

    rtc = RingtailCore()
    # equivalent to 
    rtc = RingtailCore(db_file = "output.db",storage_type="duckdb")

The Ringtail logger can be accessed through ``setup_logging``, to e.g., change log level or add an output log file:

.. code-block:: python

    from ringtail import setup_logging
    setup_logging(level="DEBUG", logfile="rt.log").


Writing docking results to the database
=======================================
To add results from docking results files, the method ``add_results_from_files`` is used. It allows one or multiple sources of results, and a range of options pertinent to the storage handling and the results processing can be set at this time. Please note that both sources of results and settings for database writing can be provided as single options or a dictionary of allowed keywords and values. If both are provided, any individual options will overwrite that given in the dictionary. The path to the receptor file is considered part of the results options. 

.. code-block:: python

    rtc.add_results_from_files( docking_results = ["lig1.sdf", "lig2.sdf", "path1/path2","filelist.txt"],
                                recursive = True, 
                                receptor_file = "receptor.pdbqt",
                                save_receptor = True,)
    
Example file list:

.. code-block:: python

    lig3.sdf
    lig4.sdf.gz

Adding a receptor
-------------------
The receptor can be added to a database by itself, either to a populated database without a receptor, or to an empty database that you anticipate to fill with docking results later. 
As of version 3.0.0, you can add either a pdbqt file or a receptor meeko.Polymer object as a .json file, or even the Polymer object itself.

.. code-block:: python
    
    rtc.save_receptor(receptor = "receptor_polymer.json")
    # equivalent to
    polymer = meeko.Polymer()
    rtc.save_receptor(receptor = polymer)

Writing docked Mols directly to the database
--------------------------------------------
The latest docking engine AD6 can produce json strings representing the ligand Mol after docking. These can be written directly to a Ringtail database without touching the file system in between, using the method ``add_mol``. ``chunk_size`` is used to determine how many Mols are processed before writing to the database, higher numbers usually meands faster write times, but uses more RAM, smaller number less RAM but longer write times. 

.. code-block:: python

    # pre process output to rdkit Mol
    def docked_mols():
        for result in my_docking_pipeline():  # AD6 docking output generator
            yield result.to_rdkit_mol()

    rtc.add_mol(
        mols=docked_mols(),
        docking_mode="ad6",
        chunk_size=10000,
        calculate_interactions=True,
    )

If using the SQLite backend, it is necessary to call ``rtc.finalize_write()`` at the end.

Writing vina results directly to the database
---------------------------------------------
Similarly to writing AD6 mols directly to a database, vina docking results can be streamed directly to a database as well using ``add_results_from_vina_string``.

.. code-block:: python
    def vina_results():
            for name, pose_pdbqt in my_vina_pipeline():  # vina docking output generator
                yield {name: pose_pdbqt}

    rtc.add_results_from_vina_string(
        results=vina_results(), 
        max_poses=3,
        chunk_size=10000,
        calculate_interactions=True,
    )

If using the SQLite backend, it is necessary to call ``rtc.finalize_write()`` at the end.

Printing a database summary
---------------------------
If at any point you wish to receive summary data from the database (eg, number of ligands, docking top percentiles) a dictionary of summary data is produced by the method db_summary_data. 

.. code-block:: python

    summary_data, requested_fields = rtc.db_summary_data()

Input options
==============
The Ringtail API uses the same options that are used in the command line interface. Relevant to adding results to the database, including how many poses of a docked ligand to save, and how to handle any duplicated ligands. 

Handling of duplicate and existing results
-------------------------------------------
With the Ringtial API you can keep adding results using the same object without specifying whether or not to ``append_results`` as is needed in the command line tool. You can specify what to do if you are adding duplicate results for a ligand, by invoking the ``duplicate_handling`` keyword with the value ``IGNORE`` (will not add the newest duplicate) or ``REPLACE`` (will overwrite the newest duplicate). Please note that the ``duplicate_handling`` option makes database writing significantly slower.

.. code-block:: python

    rtc.add_results_from_files( docking_results = "path1/",
                                duplicate_handling = "REPLACE")

Handling interaction parameters
----------------------------------
Only AutoDock-GPU is capable of performing interaction analysis at runtime, with these results being stored in the database if present. If interaction analysis is not present in the input file (including AD6 SDFs and Vina PDBQTs), it will by default be added by Ringtail during data parsing, unless the user specifies ``calculate_interactions=False`` option. **This will signifcantly decrease the total database write time.** By default, the cutoff distance for being considered an interaction is 3.7 å for hydrogen bonds and 4.0 å for van der Waals interactions. If the user would like to calculate interactions with other distance cutoffs, the option ``interaction_cutoffs`` can be used. To be able to add interactions it is important that the receptor file is provided during database write (or that the receptor has already been saved explicitly in the database). 

.. code-block:: python

    rtc.add_results_from_files( docking_results = "ligands.sdf",
                                docking_mode = "ad6",
                                receptor_file = "receptor.json",
                                save_receptor = True,
                                interaction_cutoffs = [3.7, 4.0])
    # or
    rtc.add_results_from_files( docking_results = "ligands.sdf",
                                docking_mode="ad6",
                                calculate_interactions = False)

The ``interaction_tolerance`` option for AD-GPU, which clusters output poses based on RMSD, allows the user to give more leeway for poses to pass given interaction filters. With this option, the interactions from poses within *c* angstrom RMSD of a cluster's top pose will be appended to the interactions for that top pose. The theory behind this is that this gives some sense of the "fuzziness" of a given binding pose, allowing the user to filter for interactions that may not be present for the top pose specifically, but could be easily accessible to it. 

.. code-block:: python

    rtc = RingtailCore()
    rtc.add_results_from_files( docking_results = "path1/",
                                docking_mode = "adgpu",
                                duplicate_handling = "REPLACE",
                                interaction_tolerance = 0.6)

Number of poses to save
-------------------------
By default Ringtail will store the best-scored (lowest energy) binding pose from the first |default_max_poses| poses (pose clusters for AD-GPU). Additional settings for writing to the database include how to handle the number of poses docked (``max_poses``, or ``store_all_poses`` which will overwrite the former).

.. code-block:: python

    rtc.add_results_from_files( docking_results = "path2",
                                max_poses = 5)

Iteratively appending to an SQLite database
-------------------------------------------
The DuckDB backend requires no extra processing once data has been added to the database. However, if you opt to use the SQLite backend, there is a final step where tables are indexed. If you are adding data iteratively through e.g., a for-loop and adding some number at files at once, it is time-consuming (and not necessary) to do this every iteration. Instead, you can invoke the keyword ``finalize=False``, and run the finalization method separately at the end:

.. code-block:: python

    for folder in enumerate("path_with_many_folders"):
        rtc.add_results_from_files( docking_results = folder,
                                    finalize = False)
    
    rtc.finalize_write()

Filtering
**********

Docking results stored in the Ringtail database can be screened/filtered using the ``filter`` method. The filtered poses will be stored under the chosen ``bookmark_name`` (defaults to |default_bookmark_name|), which can be accessed and used for further filtering, or output options such as writing SDF files. Optionally, after filtering a text log file can be created containing the results passing the given filter(s) if specified with the key word ``output_log``. The output log will include the ligand name (``ligname``) and docking score (``docking_score``) of every pose passing filtering criteria, and additional or other output fields can be specified with ``outfields``.

Scoring filters
=================
There are six scoring filters that include best (``ebest``) and worst docking score/energy (``eworst``), best and worst ligand efficieny (``lebest`` and ``leworst``), and results above worst docking score or ligand efficiency percentile (``score_percentile`` and ``le_percentile``, respecitvely). Some of these are internally inconsistent: if both ``eworst`` and ``score_percentile`` are used together, the ``eworst`` cutoff alone is used. The same is true of ``leworst`` and ``le_percentile``. The method can output the number of passing ligands and bookmark name, or an iterable of the database results. 

.. code-block:: python

    num_passing_ligands, bookmark_name = rtc.filter(score_percentile = 0.1)
    print(f"{num_passing_ligands} ligands passing in bookmark {bookmark_name}")
    # 3 ligands passing in bookmark passing_results

    passing_items = rtc.filter(score_percentile = 0.1, return_iter=True, outfields="ligname,pose_rank")
    for name,rank in passing_items:
        print(f"Ligand {name} passed, ranked {rank}")
        # Ligand lig1 passed, ranked 1
        # Ligand lig1 passed, ranked 2
        # Ligand lig5 passed, ranked 1


By default, only the information for the top-scoring binding pose will be written to the log. If desired, each individual passing pose can be written by using ``output_all_poses = True``. The passing results may also be ordered in the log file using the ``order_results`` option.

.. code-block:: python

    rtc.filter(eworst = -6, outfields = "ligname,docking_score,pose_rank", order_results = "reference_rmsd", output_bookmark = "eworst6")

Filtering may take from seconds to minutes, depending on the size of the database, roughly scaling as O(n) for n database Results rows (i.e. stored poses).

Interaction filters
=====================
It is possible to filter the docking results based on different types of interactions (hydrogen bonds, van der waals interactions, and reactive interactions) with specific residues. It is further possible to have ligands pass the filters while only fulfilling some of the interaction combinations in union (max number of interactions combinations missed, ``max_miss``) which gives some flexibility e.g., when you have a set of interactions of importance, but can accept that not all ligands interacts for every chosen interaction.
The available interaction filters where you can specify to lesser or greater level of detail are ``hb_interactions``, ``vdw_interactions``, and ``reactive_interactions``. Interaction filters must be specified as the interaction specifications in the order ``CHAIN:RES:NUM:ATOM_NAME``. Any combination of that information may be used, as long as 3 colons are present and the information ordering between the colons is correct. All desired interactions of a given type is specified as a list of one or more tuples of specified reactions and weather to show results that includes ``(":::", True)`` or exclude ``(":::", False)`` them as shown below for ``vdw_interactions``:

.. code-block:: python

    rtc.filter( eworst=-2,
                vdw_interactions=[('A:VAL:279:', True), ('A:LYS::', True)])

The ``max_miss`` keywords allows the user to filter by given interactions excluding up to ``max_miss`` interactions. This gives :math:`\sum_{m=0}^{m}\frac{n!}{(n-m)!*m!}` combinations for *n* interaction filters and *m* max_miss. By default, results will be given for the union of the interaction conbinations. Use with ``enumerate_interaction_combs = True`` to log ligands/poses passing each separate interaction combination. If ``max_miss > 0`` is used during filtering, a bookmark is created for each combination of interaction filters and is named ``<bookmark_name>_<n>`` where n is the index of the filter combination in the log file (indexing from 0). A union of all the enumerated bookmarks will be returned as the bookmark name ``<bookmark_name>_union``.

.. code-block:: python

    rtc.filter( eworst=-6,
                vdw_interactions=[('A:VAL:279:', True), (':LYS:162:', True)],
                hb_interactions = [("A:VAL:279:", True), ("A:LYS::)", True)],
                max_miss = 1,
                react_any = True)

``react_any`` offers an option to filtering for poses that have reactions with any residue, and it's possible to filter based on number of hydrogen bonds.

.. code-block:: python

    rtc.filter( hb_count=10,
                react_any = True)

``hb_count`` takes a plain integer and both directions are inclusive: ``hb_count=10`` keeps
poses with ten or more hydrogen bonds. A negative value inverts the direction, so
``hb_count=-10`` keeps poses with ten or fewer, and ``hb_count=0`` keeps only poses with no
hydrogen bonds at all.


Ligand filters 
===============
Several filters pertaining to the SMARTS structure of the ligand can be used. For example, ligands can be filtered for presence of certain substrctures specified by a SMARTS string using ``ligand_substruct``. The ligand name search will include any ligand names that contain the specified phrase, and does not look for exact matches only. Use the keyword ``ligand_operator`` to determine if the ligand filters should be evaluated as this ``OR`` that (default), or combined with ``AND``. 

.. code-block:: python

    rtc.filter(ligand_substruct=["[Oh]C", "C=O"], ligand_operator="AND")

The ``ligand_substruct_pos`` option may be used to filter for a specific ligand substructure to be placed within some distance of a given cartesian coordinate. The format for this option using the API is as a list of the six elements: ``["<SMARTS pattern: str>"," <index of atom in SMARTS: int>, <cutoff distance: float>, <target x coord: float>, <target y coord: float>, <target z coord: float>]``. If seachring for more than one ``ligand_substruct_pos`` make the value a list of lists. 

.. code-block:: python

    rtc.filter(ligand_substruct_pos=[["C=O", 1, 10, 102, 106, 154], ['[C][Oh]', 1, 10, 102, 106, 154]])

It's possible to filter based on either ``ligand_max_atoms``, which specifies maximum number of heavy atoms a ligand may have (anything but hydrogens), or giving a min and/or max molecular weight (g/mol).

.. code-block:: python

    rtc.filter(ligand_max_atoms=5)
    # or
    rtc.filter(ligand_min_molweight=50,ligand_max_molweight=350)

It may be desirable to filter and thus export a set of ligands based on their names, for example if a collaborator provides a list of compounds they've identified from other forms of screening. This can be done by either supplying ``ligand_name`` which is one or more strings that will be treated with wildcards ``*<ligand_name>*`` when searching the database (limited  to 50 names), or by supplying a .csv ``ligand_name_file`` which can have unlimited number of ligand names, but will not apply wildcards in its search. 

.. code-block:: python

    rtc.filter(ligand_name="my_favorite_library_chem")
    # or
    rtc.filter(ligand_name_file="my_collaborators_gave_me_this_ligand_list.csv")


Clustering
============
In addition to the filtering options outlined in the table below, ligands passing given filters can be clustered to provide a reduced set of dissimilar ligands based on Morgan fingerprints (``mfpt_cluster``) or interaction (``interaction_cluster``) fingerprints. Dissimilarity is measured by Tanimoto distance (float input to the cluster keyword) and clustering is performed with the Butina clustering algorithm. Clustering can be performed while filtering (will cluster the resulting poess), or on an existing bookmark by using the ``cluster`` method. In this case, the bookmark over which to cluster (or additional filtering) on is specified by ``input_bookmark`` (must be different from ``bookmark_name`` that contains previously filtered results).

.. code-block:: python

    num_clustered_poses, _ = rtc.filter(eworst = -6,
                                        mfpt_cluster = 0.6
                                        output_bookmark="clustered")
    print(num_clustered_poses)
    # 5

    # equivalent to 

    _, bookmark_name = rtc.filter(eworst = -6)
    num_clustered_poses, _ = rtc.cluster(cluster_type="mfp",
                                        cutoff=0.6,
                                        output_bookmark="clustered",
                                        input_bookmark=bookmark_name)
    print(num_clustered_poses)
    # 5

The user can provide a ligand name from a previously-run clustering and re-output other ligands that were clustered with that query ligand, see section below. 


Output options
***************
There are multiple options to output and visualize data in Ringtail.

Export molecule SDF files
==========================
The method ``write_molecule_sdfs`` will write SDF files for each ligand passing the filter and saved in a specified bookmark, and any additional ligands listed by name. By default it will write all ligands to one SDF file, though this can be changed by setting ``all_in_one`` to ``False``. The files will be saved to the path specified by ``sdf_path``. If none is specified, the files will be saved in the current working directory. The binding energies, ligand efficiencies, and interactions are written as properties within the SDF file, with the order corresponding to the order of the pose order.

.. code-block:: python

    rtc.write_molecule_sdfs(sdf_path = "sdf_files", bookmark_name = "eworst6", all_in_one=False)

Exporting select columns, tables, or query results as CSV files
===============================================================
If the user wishes to explore the data in CSV format, Ringtail provides three options. 
First, you can export chosen columns for any given bookmark including interaction data:

.. code-block:: python

    rtc.export_columns_as_csv(columns=["ligname","docking_score","interaction_type","rec_resname"],bookmark_name="eworst6")

For a more flexible option, you can write your own SQL and produce an output CSV file from that:

.. code-block:: python
    
    rtc.export_sql_as_csv(sql = "SELECT ligname, docking_score, leff, pose_id, ligand_smile FROM Results JOIN Ligands ON Results.ligand_id=Ligands.ligand_id", csv_name = "my_query.csv")

Lastly, you can simply dump an entire database table to CSV, though please note this can be slow for large databases. A useful option may be for the table ``Interaction_indices`` which describes each unique interaction in the database. 

.. code-block:: python

    rtc.export_table_as_csv(table="Interaction_indices", csv_name = "unique_interactions.csv")

Creating a new database from a bookmark
=======================================
A bookmark may be exported as a separate database with the ``export_bookmark_db`` method. This is particularly useful if you start with a large database, have filtered down to a reasonable subset, and wish to continue working with a smaller and faster Ringtial database (especially for e.g., downloading from an HPC to a personal computer). This will create a database of name ``<current_db_name>_<bookmark_name>.db``, unless a full path is specified with db_filepath.

.. code-block:: python 

    rtc.export_bookmark_db(bookmark_name = "eworst6",db_filepath="/Users/mydata/smol.db")

Exporting receptor information
==============================

A receptor stored in the database may be retrieved and eithed used as an object using the ``get_receptor_object`` method which returns a ReceptorData object. This dataclass holds receptor name, pdbqt blob_str if any, and the polymer_json string if any. 

.. code-block:: python 
    
    from meeko import Polymer
    rec_obj = rtc.get_receptor_object()
    rec_json = rec_obj.polymer_json
    # create a Meeko polymer object
    polymer = Polymer.from_json(rec_json)


It's possible to export the original PDBQT using ``export_receptor_pdbqt``, which will write the .pdbqt file to the working directory.

.. code-block:: python 

    rtc.export_receptor_pdbqt()

Exporting a flexible residue receptor PDBQT
===========================================
If docking was performed with a receptor with flexible residues, it's possible to export a PDB with the receptor confirmation given one or more ligands, or given a filter bookmark. For this it's necessary to either have a receptor polymer stored in the database, or it can be provided at runtime for databases created prior to the use of meeko Polymers. 

.. code-block:: python

    rtc.write_flexres_pdb(bookmark_name="best_filters", filename="flexres_of_best_filters.pdb")


Writing an output log for an existing bookmark
===============================================
Data for poses in a bookmark may be written to an output log using the ``get_previous_filter_data`` method.

.. code-block:: python

    rtc.get_previous_filter_data(outfields = "ligname,docking_score,pose_rank", bookmark_name = "eworst6", output_log = "previously_filtered_results.txt")

Find similar ligands to a clustered ligand
==========================================
If you'd like to write to an output log file all ligands clustered to a ligand of choice, a string of methods are used together. First, fetch the available clustering groups for the select ligand, then pass the chosen cluster ID to retrieve similar ligands. 

.. code-block:: python

    cluster_options = rtc.fetch_cluster_options("my_best_ligand")
    # options is a list of (cluster_id, cluster_window, name) tuples
    cluster_id = options[0][0]  # choose the desired cluster
    ligands, bookmark_name, cluster_name = rtc.fetch_clustered_similars(
        "my_best_ligand", cluster_id, output_log="cluster_to_my_best_ligand.txt"
    )

Advanced API usage
*******************

Merge databases that share the target receptor
==============================================
It's possible to perform a virtual screen on *one* target receptor where results are split across multiple Ringtail databases, for example when separating work across multiple nodes on an HPC. These can be merged with ``merge_databases``. A RingtailCore object is instantiated by one database (recommended to choose the largest in the set). This primary database is backed up by default before the merging databases are merged into it. Please note that the merging/secondary database files will not be altered, only the primary database will be written to. 
Explicit paths and/or wildcard glob file patterns are accepted.

.. code-block:: python

    rtc = RingtailCore("largest.db") # primary db
    rtc.merge_databases(["batch1.db", "batch2.db", "more_batches/*.db"])

The primary database itself and any duplicate paths are excluded automatically, and all
databases must be schema-compatible (v3). Use ``backup=False`` to skip bacup of the primary database. 

Crossreference instances of ligands betweeh different target databases
======================================================================
If screens have been performed with the same molecule library but with *different* target receptors it's possible to compare, or cross reference, what ligands occur in selections/filter bookmarks of the different screens using ``cross_reference_databases``. Each database is provided together with a bookmark name (screened data) from where the comparison will happen. The comparison relies on unique ligand names, and assumes that any ligand that shares a name across databases are the same. You can provide both "wanted" databases and "unwanted", for example if one of the screens produced a set of binding modes or pharmacophores you want to exclude from your final selection. A list of wanted and unwanted databases/bookmarks is provided as a tuple (``(database_path, bookmark_name)``). It's possible to set bookmark_name to None if you wish to compare the enture database. A bookmark prefixed with ``bookmark_prefix`` is written into each participating database for later use.
The selection (second tuple element) may also be a **status table** name (``"accepted"``, ``"maybe"`` or ``"rejected"``) instead of a bookmark. This lets you cross reference by acceptance status assigned via ``update_pose_status`` — for example, finding the ligands marked ``Accepted`` in two different target screens:

.. code-block:: python

    rtc = RingtailCore("target_A.db")
    num_shared, new_bookmarks, _ = rtc.cross_reference_databases(
        wanted_dbs=[("target_A.db", "accepted"),
                    ("target_B.db", "accepted")],
    )
    # new_bookmarks -> {database_path: "crossref_accepted"}

.. code-block:: python

    rtc = RingtailCore("target_A.db")
    num_shared, new_bookmarks, _ = rtc.cross_reference_databases(
        wanted_dbs=[("target_A.db", "passing_results"),
                    ("target_B.db", "passing_results")],
        unwanted_dbs=[("offtarget.db", "passing_results")],
        bookmark_prefix="selective",
    )
    print(f"{num_shared} ligands pass in A and B but not the off-target")
    # new_bookmarks -> {database_path: created_bookmark_name}

This method is used in the command line tool :ref:`rt_compare <compare>`.

Retrieve ligand RDKit Mols 
===========================
You can create RDKit Mols to be used elsewhere, using the methods ``fetch_select_ligands_poses`` and ``create_rdkit_mols_by_ligand``. The first method builds a list of poses for which to build Mols based on a. ligand name(s), b. already known pose_ids, and either of these can be bracketed by c. bookmark name which can also be used in isolation (i.e., all poses in this bookmark). The resulting dict is passed as input to the create mol method.

.. code-block:: python

    ligand_pose_dict = rtc.fetch_select_ligands_poses(bookmark_name = "eworst6")

    mols = rtc.create_rdkit_mols_by_ligand( ligand_poses=ligand_pose_dict, 
                                            include_interactions=False, 
                                            include_comment=True)
    for ligmol in rtc.create_rdkit_mols_by_ligand(  ligand_pose_dict
                                                    include_interactions=False, 
                                                    include_comment=True):
        mol = ligmol.mol
        print(ligmol.ligname, mol.GetNumConformers(), "poses")
        # per-pose property lists (order matches the conformers/poses):
        print(ligmol.properties["Binding energies"])
        print(ligmol.properties["Ligand efficiencies"])
        print(ligmol.properties["Comment"])

Each yielded ``LigandMol`` has the fields ``ligname``, ``mol``, ``properties``, ``flexres_per_pose``, and ``flexres_residues``. Interaction strings are embedded by default ``include_comment=True`` adds per-pose comments, and ``flexres_data=rtc.make_receptor_flexres_mols()`` also builds flexible-residue mols.

Retrieve interactions for selection
====================================


Calculating interactions after building a Ringtail database
============================================================
If you have created a Ringtail database with ``calculate_interactions=False``, or you would like to re-calculate pose-receptor interactions with different interaction distance limtis, this is possible using the method ``add_interactions``. Please note, though, this calculation can take a long time for large databases. If calculation fails part way through, it is possible to restart it simply by running the method again, as Ringtail keeps track of the progress with a transaction tracking table until the process finished. 

.. code-block:: python

    rtc.add_interactions(hb_cutoff=3.9,vdw_cutoff = 4.2)


Flagging or commenting a poses
=====================================
Ringtail 3 comes with the ability to flag or assign status to poses, meaning the pose will be flagged as "Accepted"/1, "Maybe"/2, or "Rejected"/3 (0 is no status), giving more reversible tracking abilities in Ringtail while screening is ongoing. The flagging tables "Accepted", "Maybe", and "Rejected" are completely useable the way bookmarks are, and can be used as basis for further filtering, or export options. It's also possible to add a comment to a pose, stored in the schema table "pose_comments". Please note you need to know the ``pose_id`` to use these methods, which can be retrieved e.g., by ``return_iter=True`` and including ``pose_id`` in the ``outfields`` in filter(). 

.. code-block:: python

    # accept poses
    rtc.update_pose_status([1,5,10,17],1)

    # reject poses
    rtc.update_pose_status([2,4,7],3)

    # remove flags
    rtc.update_pose_status([1,2,4,5],0)

    # add comments
    rtc.set_pose_comment(10,"This pose has a really neat ring")

Bookmark management
===================
You can get a list of all bookmarks, as well as delete bookmarks based on name.
This will not delete any docking data, simply remove the association between poses and any given filter bookmark (and free its name up for re-use). 

.. code-block:: python

    bookmark_list = rtc.get_bookmark_names()

    for bookmark in bookmark_list:
        rtc.delete_bookmark(bookmark)

Writing raw SQL
===============
It's possible to perform SQL queries directly using the method ``db_query``, this assumes prior knowledge of SQL and backend dialect. The method accepts parameters, as well as having the ability to write, by committing changes, so use with caution. Ringtail does not currently have guards against writes performed usign this method (or any other direct sql access). 

.. code-block:: python

    query = "SELECT docking_score, leff, pose_id, ligand_smile FROM Results JOIN Ligands ON Results.ligand_id=Ligands.ligand_id WHERE Ligands.ligname = ?;"
    ligand_poes_list = rtc.db_query(query, params=("myligand",))