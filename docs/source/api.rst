
.. _api:

API procedures
###############

The Ringtail API allows for a more advanced and flexible use of Ringtail where the user can create their own scripts. In the case of the docking engines that can output docking results directly as a variable, instead of writing to the file system, including vina and AD6, Ringtail can be directly integrated in the virtual screen pipeline and write docking results to the database directly. 

The Ringtail API can thus process output files from AD6 (SDFs), AutoDock-GPU (DLGs), and vina (PDBQTs), as well as docked AD6 RDKit molecules and vina result strings. Ringtail is intended to be used for a set of docking results for a single target and binding site. This may include multiple ligand libraries as long as the target and binding site is the same as Ringtail only allows one receptor per database (checks are performed based on provided receptor name). However, if you do multiple screens with e.g., different binding pockets for the same receptor, this information will not be checked, and Ringtail will allow any data to be appended to a database as long as the receptor name is the same. One receptor JSON (or PDBQT) may also be saved to the database. Ringtail has tools for performing comparisons of ligands with different targets (:ref:`compare`).

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
    setup_logging(level="DEBUG", logfile="rt.log")


Writing docking results to the database
=======================================
To add results from docking results files, the method ``add_results_from_files`` is used. It allows one or multiple sources of results, and a range of options pertinent to the storage handling and the results processing can be set at this time. ``docking_results`` accepts a single string or a list of strings, each of which may be a docking results file, a directory, or a text file listing result files.

.. code-block:: python

    rtc.add_results_from_files( docking_results = ["lig1.sdf", "lig2.sdf", "path1/path2","filelist.txt"],
                                recursive = True,
                                receptor_file = "receptor.json",
                                save_receptor = True,)

Example file list. Relative paths in a file list are resolved against the current working directory, not against the location of the list:

.. code-block:: text

    lig3.sdf
    lig4.sdf.gz

Adding a receptor
-------------------
The receptor can be added to a database by itself, either to a populated database without a receptor, or to an empty database that you anticipate to fill with docking results later. 
As of version 3.0.0, you can add either a PDBQT file, a Meeko Polymer JSON file,
or a dictionary mapping a receptor name to a ``Polymer`` object.

.. code-block:: python

    from pathlib import Path
    from meeko import Polymer

    rtc.save_receptor(receptor = "helix--scofu01.json")
    # equivalent object form
    polymer = Polymer.from_json(Path("helix--scofu01.json").read_text())
    rtc.save_receptor(receptor={"helix--scofu01": polymer})

When a file is given, the receptor name is taken from the file name (without extension); in the dictionary form, the key is stored as the receptor name. Use the same name in both forms, since results added later with a receptor file of a different name are rejected as belonging to another receptor.

Writing docked Mols directly to the database
--------------------------------------------
AD6 docking results can be converted to RDKit molecules and written directly to a
Ringtail database with ``add_mol``, without intermediate files. ``chunk_size`` controls
how many molecules are processed per database write; larger chunks generally improve
throughput while using more memory. To calculate interactions, save a receptor to the
database first (see above), or pass it using ``receptor_string``. No receptor is needed
when ``calculate_interactions=False``.

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
Vina PDBQT result strings can likewise be streamed directly with
``add_results_from_vina_string``. To calculate interactions, save a receptor first, or
pass it using ``receptor_string``. No receptor is needed when
``calculate_interactions=False``.

.. code-block:: python

    def vina_results():
        for name, pose_pdbqt in my_vina_pipeline():  # Vina docking output generator
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
With the Ringtail API you can keep adding results using the same object without specifying whether or not to ``append_results`` as is needed in the command line tool. You can specify what to do if you are adding duplicate results for a ligand, by invoking the ``duplicate_handling`` keyword with the value ``IGNORE`` (keeps the existing entry and does not add the new duplicate) or ``REPLACE`` (overwrites the existing entry with the new duplicate; the status assignment and comment of a replaced pose are deleted with it). Please note that the ``duplicate_handling`` option makes database writing significantly slower.

.. code-block:: python

    rtc.add_results_from_files( docking_results = "path1/",
                                duplicate_handling = "REPLACE")

Bookmarks store the poses that passed filtering, so adding results to a database that has bookmarks would leave every bookmark out of date. ``add_results_from_files``, ``add_mol``, and ``add_results_from_vina_string`` therefore raise ``OptionError`` in that case, and write nothing, unless ``consent=True`` is passed, which deletes all bookmarks, filters and clusterings before adding the results. Status assignments and pose comments are kept. ``has_filter_data()`` tells you in advance whether consent will be needed.

.. code-block:: python

    if rtc.has_filter_data():
        print("Adding results will delete these bookmarks:", rtc.get_bookmark_names())
    rtc.add_results_from_files(docking_results="new_batch/", consent=True)

Handling interaction parameters
----------------------------------
For AutoDock-GPU results, Ringtail reads interactions already present in the DLG files. For AD6 and Vina results, Ringtail calculates interactions by default when a receptor is supplied. Use ``calculate_interactions=False`` to disable calculation and reduce database write time. The default distance cutoffs are 3.7 Å for hydrogen bonds and 4.0 Å for van der Waals interactions. Use ``interaction_cutoffs`` to specify different cutoffs. Provide the receptor using ``receptor_file`` during database writing, or use a receptor already saved in the database.

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
By default Ringtail stores up to |default_max_poses| top-ranked poses per ligand (for AD-GPU, the best-scored pose from each of the first |default_max_poses| pose clusters). Use ``max_poses`` to change this number, or ``store_all_poses`` (which overrides ``max_poses``) to store every pose.

.. code-block:: python

    rtc.add_results_from_files( docking_results = "path2",
                                max_poses = 5)

Iteratively appending to an SQLite database
-------------------------------------------
The DuckDB backend requires no extra processing once data has been added to the database. However, if you opt to use the SQLite backend, there is a final step where tables are indexed. If you are adding data iteratively through e.g., a for-loop and adding some number at files at once, it is time-consuming (and not necessary) to do this every iteration. Instead, you can invoke the keyword ``finalize=False``, and run the finalization method separately at the end:

.. code-block:: python

    from pathlib import Path

    for folder in Path("path_with_many_folders").iterdir():
        rtc.add_results_from_files( docking_results = str(folder),
                                    finalize = False)
    
    rtc.finalize_write()

Filtering
**********

Docking results stored in the Ringtail database can be screened with ``filter``. Passing
poses are stored under ``output_bookmark`` (|default_bookmark_name| by default) for
further filtering or export. Supplying ``output_log`` writes ligand names and docking
scores for the best-scoring passing pose of each ligand; set ``output_all_poses=True``
to include every passing pose, and use ``outfields`` to select other columns.

Before Ringtail v3, bookmarks were stored as database views that evaluated the filter query when accessed. In v3, bookmarks instead store the selection of poses that passed filtering. Because a stored selection cannot include results added later, adding docking results requires ``consent=True`` and deletes all bookmarks (see *Handling of duplicate and existing results* above); rerun filtering afterwards.

Combining filter groups
=======================
Use ``filters=`` to combine nested AND/OR groups. Supply either a ``Filters`` object or a dictionary; this cannot be combined with individual filter arguments such as ``eworst=``. In the dictionary, supply the operator using the ``"op"`` key (``"and"`` or ``"or"``), and list the groups under ``"children"``. Criteria within each leaf dictionary are combined with AND.

.. code-block:: python

    rtc.filter(
        filters={
            "op": "or",
            "children": [
                {"eworst": -9, "hb_interactions": [("A:GLU:246:", True)]},
                {"eworst": -11},
            ],
        },
        output_bookmark="combined_hits",
    )

This selects poses scoring ≤ −9 with a hydrogen bond to ``A:GLU:246:``, or poses scoring ≤ −11.

Scoring filters
=================
There are six scoring filters: best (``ebest``) and worst (``eworst``) docking
score, best (``lebest``) and worst (``leworst``) ligand efficiency, and docking-score
or ligand-efficiency percentiles (``score_percentile`` and ``le_percentile``).
``eworst`` cannot be combined with ``score_percentile``; likewise, ``leworst`` cannot
be combined with ``le_percentile``; either combination raises ``OptionError``. The
percentiles are given in percent, so ``score_percentile=0.1`` keeps the top 0.1% of
docking scores. The method returns the number of passing ligands and bookmark name,
or result rows when ``return_iter=True`` (one row per ligand, or one per passing pose
with ``output_all_poses=True``).

.. code-block:: python

    num_passing_ligands, bookmark_name = rtc.filter(score_percentile = 0.1)
    print(f"{num_passing_ligands} ligands passing in bookmark {bookmark_name}")
    # 2 ligands passing in bookmark passing_results

    passing_items = rtc.filter(score_percentile = 0.1, return_iter=True, output_all_poses=True, outfields="ligname,pose_rank")
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
Results can be filtered by hydrogen-bond, van der Waals, or reactive interactions
with specific receptor residues. Use ``hb_interactions``, ``vdw_interactions``, or
``reactive_interactions`` with specifications in ``CHAIN:RES:NUM:ATOM_NAME`` order.
Any field may be blank, but all three colons are required. Each API value is a tuple
of the specification and a Boolean indicating whether the interaction is required
(``True``) or excluded (``False``).

.. code-block:: python

    rtc.filter( eworst=-2,
                vdw_interactions=[('A:VAL:243:', True), ('A:GLU::', True)])

``max_miss`` allows up to the specified number of interaction filters to be absent.
By default, Ringtail stores the union of the passing combinations in
``output_bookmark``. With ``enumerate_interaction_combs=True``, the union is instead
stored as ``<output_bookmark>_union``, and a numbered bookmark
``<output_bookmark>_<n>`` is created for each combination that has passing poses.

.. code-block:: python

    rtc.filter( eworst=-6,
                vdw_interactions=[('A:VAL:243:', True), (':GLU:246:', True)],
                hb_interactions = [("A:VAL:245:", True), ("A:GLU:246:", True)],
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
Several filters pertaining to the SMARTS structure of the ligand can be used. For example, ligands can be filtered for presence of certain substructures specified by a SMARTS string using ``ligand_substruct``. The ligand name search will include any ligand names that contain the specified phrase, and does not look for exact matches only. Use the keyword ``ligand_operator`` to determine if the ligand filters should be evaluated as this ``OR`` that (default), or combined with ``AND``.

.. code-block:: python

    rtc.filter(ligand_substruct=["[Oh]C", "C=O"], ligand_operator="AND")

The ``ligand_substruct_pos`` option may be used to filter for a specific ligand substructure to be placed within some distance of a given cartesian coordinate. The format for this option using the API is as a list of the six elements: ``[<SMARTS pattern: str>, <index of atom in SMARTS: int>, <cutoff distance: float>, <target x coord: float>, <target y coord: float>, <target z coord: float>]``. If searching for more than one ``ligand_substruct_pos`` make the value a list of lists.

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
Passing ligands can be clustered into dissimilar representatives using Morgan
fingerprints (``mfpt_cluster``) or interaction fingerprints
(``interaction_cluster``). Ringtail uses Tanimoto distance and Butina clustering.
Clustering can run during filtering, or afterward with ``cluster`` using a distinct
``input_bookmark`` and ``output_bookmark``.

.. code-block:: python

    num_clustered_poses, _ = rtc.filter(eworst = -6,
                                        mfpt_cluster = 0.6,
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
The method ``write_molecule_sdfs`` will write SDF files for each ligand passing the filter and saved in a specified bookmark. Use ``ligname`` to export only the listed ligands; if a bookmark is also given, only the listed ligands' poses in that bookmark are written. By default it will write all ligands to one SDF file, though this can be changed by setting ``all_in_one`` to ``False``. The files will be saved to the path specified by ``sdf_path``. If none is specified, the files will be saved in the current working directory. The binding energies, ligand efficiencies, and interactions are written as SDF properties in the corresponding pose order.

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
A bookmark may be exported as a separate database with the ``export_bookmark_db`` method. This is particularly useful if you start with a large database, have filtered down to a reasonable subset, and wish to continue working with a smaller and faster Ringtail database (especially for e.g., downloading from an HPC to a personal computer). This will create a database of name ``<current_db_name>_<bookmark_name>.db``, unless a full path is specified with db_filepath.

.. code-block:: python 

    rtc.export_bookmark_db(bookmark_name = "eworst6",db_filepath="/Users/mydata/smol.db")

Exporting receptor information
==============================

A receptor stored in the database may be retrieved as an object using the ``get_receptor_object`` method, which returns a ``ReceptorData`` object. This dataclass holds the receptor name, PDBQT blob string if present, and Polymer JSON string if present.

.. code-block:: python 
    
    from meeko import Polymer
    rec_obj = rtc.get_receptor_object()
    rec_json = rec_obj.polymer_json
    # create a Meeko polymer object
    polymer = Polymer.from_json(rec_json)


It's possible to export the receptor as a PDBQT using ``export_receptor_pdbqt``, which will write ``<receptor_name>.pdbqt`` to the working directory. If the receptor was saved as a PDBQT, the original file is written; if it was saved as a Meeko Polymer, the PDBQT is generated from the Polymer.

.. code-block:: python 

    rtc.export_receptor_pdbqt()

Exporting a receptor PDB for selected poses
===========================================
If docking was performed with a receptor with flexible residues, it's possible to export a PDB with the receptor conformation given one or more ligands, or given a filter bookmark. For this it's necessary to either have a receptor polymer stored in the database, or it can be provided at runtime (``receptor_polymer``) for databases created prior to the use of meeko Polymers. One PDB is written per ligand, named ``<filename>_<ligname>.pdb``, with one MODEL per pose. Writing more than 10 ligands requires ``consent=True``; without it, nothing is written.

.. code-block:: python

    rtc.write_flexres_pdb(bookmark_name="eworst6", filename="flexres_eworst6.pdb")
    # writes flexres_eworst6_<ligname>.pdb for each ligand in the bookmark


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
    # cluster_options is a list of (cluster_id, clustered_bookmark, cluster_name) tuples,
    # e.g. (1, "clustered_preclust", "mfp_0.6")
    cluster_id = cluster_options[0][0]  # choose the desired cluster
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
databases must be schema-compatible (v3). Use ``backup=False`` to skip backup of the primary database.

Crossreference instances of ligands between different target databases
======================================================================
If screens have been performed with the same molecule library but with *different* target receptors it's possible to compare, or cross reference, what ligands occur in selections/filter bookmarks of the different screens using ``cross_reference_databases``. Each database is provided together with a bookmark name (screened data) from where the comparison will happen. The comparison relies on unique ligand names, and assumes that any ligand that shares a name across databases are the same. You can provide both "wanted" databases and "unwanted", for example if one of the screens produced a set of binding modes or pharmacophores you want to exclude from your final selection. A list of wanted and unwanted databases/bookmarks is provided as a tuple (``(database_path, bookmark_name)``). It's possible to set bookmark_name to None if you wish to compare the entire database. A bookmark prefixed with ``bookmark_prefix`` is written into each participating database for later use.
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

    for ligmol in rtc.create_rdkit_mols_by_ligand(
        ligand_pose_dict,
        include_interactions=False,
        include_comment=True,
    ):
        mol = ligmol.mol
        print(ligmol.ligname, mol.GetNumConformers(), "poses")
        # per-pose property lists (order matches the conformers/poses):
        print(ligmol.properties["Binding energies"])
        print(ligmol.properties["Ligand efficiencies"])
        # only the poses that have a comment contribute an entry:
        print(ligmol.properties["Comment"])

Each yielded ``LigandMol`` is a named tuple with the fields ``ligname``, ``mol``, ``flexres_per_pose``, ``properties``, and ``flexres_residues``, in that order. ``Binding energies`` and ``Ligand efficiencies`` hold one value per pose. Interaction strings are embedded by default (``include_interactions=True``). ``include_comment=True`` adds pose comments. The ``Interactions`` and ``Comment`` lists only contain entries for poses that have them, so they do not necessarily line up with the conformers. ``flexres_data=rtc.make_receptor_flexres_mols()`` also builds flexible-residue mols.

Calculating interactions after building a Ringtail database
============================================================
If you have created a Ringtail database with ``calculate_interactions=False``, or you would like to re-calculate pose-receptor interactions with different interaction distance limits, this is possible using the method ``add_interactions``. Please note, though, this calculation can take a long time for large databases. If calculation fails part way through, it is possible to restart it simply by running the method again, as Ringtail keeps track of the progress with a transaction tracking table until the process finished.

.. code-block:: python

    # database written with calculate_interactions=False
    rtc.add_interactions(hb_cutoff=3.9, vdw_cutoff=4.2)

    # database that already has interactions (including every AutoDock-GPU database)
    rtc.add_interactions(hb_cutoff=3.9, vdw_cutoff=4.2, consent=True)

Recalculating over a database that already has interactions deletes them, so it requires ``consent=True``; without it, nothing is changed and the returned dict reports ``completed=False``. Finishing an interrupted run does not, since it only computes the poses that were never reached. Pass ``backup=True`` to have Ringtail clone the database before it deletes anything; the copy is written next to the original with ``.bk`` appended. There is also a command line equivalent, ``rt_recalc_interactions``, which takes one or more databases.

The interactions are computed against the receptor stored in the database, so one has to be there: a database built without ``save_receptor=True`` has no receptor, which is easy to end up with since ``receptor_file`` is optional for vina. ``add_interactions`` raises ``RTCoreError`` in that case, and does so before deleting anything, so nothing is lost — add the receptor with ``save_receptor`` and call it again.

Because the calculation is long, it can report progress and be stopped. ``progress_callback`` is called with ``(poses_done, poses_total)`` after every committed batch, and ``should_cancel`` is checked between batches; returning ``True`` from it stops the run on a committed boundary and leaves the database resumable. The method returns a dict saying what happened.

.. code-block:: python

    import threading

    stop_flag = threading.Event()  # call stop_flag.set() from elsewhere to stop the run

    result = rtc.add_interactions(
        consent=True,
        backup=True,
        chunk_size=200,
        progress_callback=lambda done, total: print(f"{done}/{total}"),
        should_cancel=stop_flag.is_set,
    )
    if not result["completed"]:
        print(f"stopped after {result['poses_done']} of {result['poses_total']} poses")

A run that was cancelled or interrupted leaves its tracking table behind, which is the only record that the database is half recomputed — the poses that were not reached have no interactions at all until it is finished. ``interaction_recalc_status()`` reports whether that is the case, and at which cutoffs the unfinished run was working. Resuming at any other cutoffs raises ``OptionError``, since the database would otherwise hold two different calculations with no record of which pose got which.

.. code-block:: python

    status = rtc.interaction_recalc_status()
    if status["pending"]:
        hb, vdw = status["cutoffs"]
        rtc.add_interactions(hb_cutoff=hb, vdw_cutoff=vdw)

Recalculating changes ``Results.num_hb`` and ``Results.num_interactions``, so bookmarks that were filtered on interactions no longer describe what their filters would now select. Their poses are untouched, but ``bookmarks_with_interaction_filters()`` lists the ones worth re-running afterwards.


Flagging or commenting on poses
===============================
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
It's possible to perform SQL queries directly using the method ``db_query``, this assumes prior knowledge of SQL and backend dialect. The method accepts parameters, as well as having the ability to write, by committing changes, so use with caution.

.. warning::
    Ringtail has no guards against writes performed using this method, or any other direct SQL access. On DuckDB, statements run in autocommit unless a transaction is open, so ``commit=False`` does not prevent a write from persisting.

.. code-block:: python

    query = "SELECT docking_score, leff, pose_id, ligand_smile FROM Results JOIN Ligands ON Results.ligand_id=Ligands.ligand_id WHERE Ligands.ligname = ?;"
    ligand_poses_list = rtc.db_query(query, params=("myligand",))
