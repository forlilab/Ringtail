.. _compress:

Compressing databases for transfer
##################################

Ringtail databases from large screens (>>millions of ligands) can be tens of gigabytes. If transfer is desired or required, for example off an HPC unto a personal computer, Ringtail can shrink the file either by simply compressing, or by applying a simple filter and then compressing. The methods will be described both for the command line as well as scripting with the API. The original database is never modified destructively.


From the command line
=====================

``rt_compress_db`` optionally filters a database by docking score and/or ligand efficiency into a minimal database (using ``export_bookmark_db`` if filtering), then compresses it. Without any filters, it simply compresses the database as-is. The compressed database is unpacked with ``rt_decompress_db``.

.. code-block:: bash
    
    # just compress
    $ rt_compress_db -i screen.db
    # -> screen.db.zst

    # filter, then compress
    $ rt_compress_db -i screen.db -e -9.0 -le -0.4 --compressor zstd --level 18
    # -> screen_filtered.db.zst  (only the passing ligands/poses)

    # unpack/decompress
    $ rt_decompress_db -i screen_filtered.db.zst
    # -> screen_filtered.db

``rt_compress_db`` arguments:

* ``-i`` / ``--input`` — input Ringtail database (required)
* ``-o`` / ``--output`` — compressed database path (default: input name + the method's extension)
* ``-e`` / ``--eworst`` — highest/worst docking score to keep 
* ``-le`` / ``--leworst`` — worst ligand efficiency to keep
* ``--compressor`` — ``zstd`` (default), ``gzip``, or ``xz`` (``zstd`` falls back to ``gzip`` if the binary is missing).
* ``--level`` — compression level (default 18; ``zstd`` 1–22, ``gzip`` 1–9, ``xz`` 0–9).
* ``--keep-db`` — keep the intermediate uncompressed filtered database.

``rt_decompress_db`` only needs the path to the compressed database (``-i``); the method is read from the
``.zst`` / ``.gz`` / ``.xz`` extension, and ``-o`` optionally sets the output path.

From the API
============

The same compression is available as ``ringtail.util.compress_file`` and ``decompress_file``, which by itself will simply compress the database file. To produce a filtered, minimal database first, filter and export bookmark as a subset database, then compress it.

.. code-block:: python

    from ringtail import RingtailCore
    from ringtail.util import compress_file, decompress_file

    # (optional) filter down to a small standalone database first
    rtc = RingtailCore("screen.db")
    rtc.filter(eworst=-9.0, leworst=-0.4, output_bookmark="hits")
    rtc.export_bookmark_db(bookmark_name="hits", db_filepath="hits.db")

    # compress (original left intact); returns the compressed db path actually written
    artifact = compress_file("hits.db", method="zstd", level=18)   # -> "hits.db.zst"

    # ...transfer the compressed db, then on the other side:
    db_path = decompress_file("hits.db.zst")                       # -> "hits.db"

.. note::
   ``zstd`` (the default) is a fast, multithreaded compressor; if the ``zstd`` binary is not
   on the ``PATH``, Ringtail automatically falls back to ``gzip``. The source file is never
   modified, moved, or deleted.

See :ref:`big_data` for how compression fits into a broader strategy for large screens.
