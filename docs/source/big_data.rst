.. _big_data:

Strategies for large virtual screens
####################################

Screens of millions of ligands can produce large databases, fortunately Ringtail has several features that keep this manageable. The recommended pattern is **rough cut filter early, carry only what you need, and reassemble smaller pieces later if needed**.

Choose the backend and write settings
======================================

* **Use the DuckDB backend** (``--storage_type duckdb`` or ``RingtailCore(storage_type="duckdb")``) for large screens: filtering is substantially faster, and database-creation throughput stays roughly flat as the database grows.
* In Ringtail v3 pose coordinates are stored in a compact native format, so databases are considerably smaller than in earlier versions for the same data.
* When creating a database programmatically for AD6 or vina docking results, stream results with ``add_mol`` / ``add_results_from_vina_string`` and tune ``chunk_size`` to trade memory for write speed (see :ref:`api`).

Filter, export, compress
========================

The most common way to shrink a finished screen database for analysis or transfer:

#. **Rough filter** applied to the full database to eliminate poor docking results, stored in a bookmark. Alternatively, run ``rt_process_vs read -i <database>.db -su`` to get a rough idea of the docking score distribution in the database, and use this to inform your rough cut filter. 
#. **Export** that bookmark to a small, standalone database with ``export_bookmark_db`` resulting in a compact database of only the hits.
#. **Compress** it for transfer with ``rt_compress_db`` / ``compress_file`` (see :ref:`compress`).

``rt_compress_db`` can do all three in one call with ``-e`` / ``-le``:

.. code-block:: bash

    $ rt_compress_db -i screen.db -e -9.0 -le -0.4
    # -> screen_filtered.db.zst : filtered, exported, and compressed in one step

Batched screens: per-batch filter, then merge
==============================================

Very large screens are often docked in batches (e.g. one HPC job per chunk of the library), giving one database per batch rather than a single giant one. To recombine them into one compact database of hits:

#. **Filter each batch database** with docking score or ligand efficiency (use **absolute** cutoffs — see the warning).
#. **Export** each batch's bookmark to a small database with ``export_bookmark_db``.
#. **Merge** the small per-batch databases into one with ``merge_databases`` (all batches must share the same target receptor, which they do for a single screen — see :ref:`api`).

.. code-block:: python

    from ringtail import RingtailCore

    # per batch: filter to hits and export a small database
    for batch in ["batch1.db", "batch2.db", "batch3.db"]:
        rtc = RingtailCore(batch)
        rtc.filter(eworst=-9.0, output_bookmark="hits")
        rtc.export_bookmark_db(bookmark_name="hits",
                               db_filepath=batch.replace(".db", "_hits.db"))

    # combine the small hit databases into one
    combined = RingtailCore("batch1_hits.db")
    combined.merge_databases(["batch2_hits.db", "batch3_hits.db"])

The same merge is available from the command line with ``rt_merge``, which merges one or more secondary databases into a primary one (the primary is backed up first unless ``--dont_backup_db1`` is given):

.. code-block:: bash

    $ rt_merge --primary_db batch1_hits.db --secondary_db batch2_hits.db batch3_hits.db

.. warning::
   **Relative filters do not compose across batched databases.** Percentile filters
   (``score_percentile`` / ``le_percentile``) are computed *per database*, so filtering each
   batch by percentile and then merging is **not** the same as ranking the whole screen. When
   you plan to merge batched databases, filter with **absolute** cutoffs (``eworst`` /
   ``leworst``). If you need true screen-wide percentiles, append everything into one database
   first and filter once.

Is "filter then merge" the right approach?
==========================================

Yes — *when the screen is already split into separate databases*. It keeps every step small (no single multi-GB database to filter or move) and yields one tidy database of hits. If you can instead combine all batches into **one** database as you dock, the simpler path is to filter that once and ``export_bookmark_db`` the result; merging is only needed when the data started out in batches.
