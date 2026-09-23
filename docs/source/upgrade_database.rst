.. _upgrade_database:

Upgrade any Ringtail database to v3
###################################

If you have previously written a database with Ringtail v<3.0, it will need to be updated to be compatible with the newest v3 Ringtail package. The CLI ``rt_upgrade_db`` performs this upgrade. All existing filters and bookmarks are removed during the update.

.. code-block:: bash

    $ rt_upgrade_db -d old_database_1.db old_database_2.db

Before upgrading, ``rt_upgrade_db`` asks for confirmation: type ``yes`` to proceed. Any other answer leaves the databases unchanged and exits with status 1. Databases are upgraded in place and no backup is made, so copy any database you want to keep in its original version first.

Multiple databases may be specified at once. The update may take a while depending on the size of the database.

If you need to upgrade an older database to work with any version other than the latest, simply specify the schema version you need to upgrade to (``1.1.0``, ``2.0.0``, or ``3.0.0``, the default):

.. code-block:: bash

    $ rt_upgrade_db -d older_database_1.db --version 2.0.0
