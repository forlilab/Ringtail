.. _faq:

Frequently asked questions
#############################


Locked or corrupt database
--------------------------
Occassionally, errors may occur during database reading/writing that corrupt the database. This may result in the database becoming locked. First, find any processes that still have the database file open:

.. code-block:: bash

    lsof /path/to/output.db

This lists each process holding the file open; the PID is shown in the second column. Kill it by PID (use ``-9`` only if it does not stop):

.. code-block:: bash

    kill <PID>
    kill -9 <PID>

Alternatively, kill every process using the file in one step:

.. code-block:: bash

    fuser -k /path/to/output.db

If the database is still locked or corrupted after this, it is recommended to delete the existing database and re-write it from scratch.

PDBQT file parsing issues
-------------------------
When writing from Vina PDBQTs, ensure there are no other PDBQTs (input or receptor) in directories specified within `docking_results` UNLESS the receptor PDBQT is specified with the `receptor_file` option in the same command line/method call.
