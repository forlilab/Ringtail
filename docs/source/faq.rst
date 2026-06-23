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

SQLite auto creation of db file
--------------------------------
If using the SQLite backend, please note that SQLite will always make a new, empty database file whenever a creation is connected, i.e., whenever Ringtail is used, even if no data is written to it (e.g., if database creation fails). In Ringtail write-mode, if a database does not already exist, this is not a problem, and a database will be created and populated as expected. If Ringtail is used in read-mode (e.g., filtering) and there is a spelling error in the database file name Ringtail will throw an error since the database is empty. However, since a connection was made to the misspelled db to check if the db was empty or not, a new empty db file will be created with the misspelled name. Since the file is empty, Ringtail will keep throwing errors if you try to write to or read from it, and if you wish to proceed e.g., writing to the database with the formerly misspelled name, you will have to delete the empty file first. 

PDBQT file parsing issues
-------------------------
When writing from Vina PDBQTs, ensure there are no other PDBQTs (input or receptor) in directories specified within `docking_results` UNLESS the receptor PDBQT is specified with the `receptor_file` option in the same command line/method call.
