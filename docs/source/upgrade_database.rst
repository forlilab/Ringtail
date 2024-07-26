.. _upgrade_database:

Updating database written with v1.0.0/v1.1.0 to work with v2.0.0
#################################################################

If you have previously written a database with Ringtail v<2.0.0, it will need to be updated to be compatible with the newest v2.0.0 Ringtail package. We have included a script ``rt_db_to_v200.py`` to perform this updated. Please note that all existing bookmarks will be removed during the update. The usage is as follows:

.. code-block:: bash

    $ rt_db_to_v200.py -d <v2.0.0 database 1 (required)> <v2.0.0 database 2+ (optional)>


Multiple databases may be specified at once. The update may take a few minutes per database.

Updating database written with v1.0.0 to work with v1.1.0
##########################################################

If you have previously written a database with Ringtail v1.0.0, it will need to be updated to be compatible with filtering with v1.1.0. We have included a script ``rt_db_v100_to_v110.py`` to perform this updated. Please note that all existing bookmarks will be removed during the update. The usage is as follows:

.. code-block:: bash

    $ rt_db_v100_to_v110.py -d <v1.0.0 database 1 (required)> <v1.0.0 database 2+ (optional)>


Multiple databases may be specified at once. The update may take a few minutes per database.

If trying to read a database created with Ringtail v1.0.0 with a newer version of Ringtail, you may encounter errors related to changes to the internal database structure. If you encounter this, run the follow commands (example of database named ``output.db``):

.. code-block:: bash

    $ sqlite3 output.db
    # opens database in sqlite
    > ALTER TABLE Results RENAME energies_binding TO docking_score;
    > ALTER TABLE Bookmarks ADD COLUMN filters;

If you encounter further errors related to views/bookmarks, please contact the ForliLab.