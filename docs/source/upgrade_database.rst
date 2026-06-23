.. _upgrade_database:

Upgrade any Ringtail database to v3
###################################

If you have previously written a database with Ringtail v<3.0, it will need to be updated to be compatible with the newest v3 Ringtail package. The CLI ``rt_upgrade_db.py`` will perform this upgrade, and please note that all existing bookmarks and screening tables will be removed during the update. The usage is as follows:

.. code-block:: bash

    $ rt_upgrade_db -d old_database_1.db (required) old_database_2+.db (optional)


Multiple databases may be specified at once. The update may take a while depending on the size of the database.

If you need to upgrade an older database to work with any version other than the latest, simply specify the schema version you need to upgrade to:
.. code-block:: bash

    $ rt_upgrade_db -d older_database_1.db --version 2.0.0

