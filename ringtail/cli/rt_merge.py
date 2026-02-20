# Script to merge two databases

"""
This script merges two (or more) ringtail databases of db schema version 2.1 or newer, and it aims to maintain all relationships, primary keys, foreign keys, etc.
The script will take a primary and a secondary database, and merge the secondary database into the primary database.
By default a backup will be made of the primary database before merging, but the the user has the option to suppress this behavior and proceed without a backup.
If the merge is not successful, the primary database may have the new merge table, and one or more corrupt data tables. To restore simply delete the corrupt file,
and remove '.backup' from the filename of the backup file.

It is recommended that the user chooses as the primary database, whichever file size is larger, as this will minimize the data transfer.

In the "new" merged database there will be two tables describing the merge history:
    1. name of db file merged into it, datetime merge started, datetime merge ended
    2. table of primary keys: dbfile, table (bc only one PK per table), original value of PK, new value of PK in merge table
    The new PK value will in most instances be created by adding the max PK value from primary database to all PKs from secondary database.

If the user wishes to merge databases of filtered results only they need to run the ringtail command make db from filtered results, and use that database in the merge process.

The script does not look for or handle duplicates in the Results or Interactions table. The Ligands table will by default not allow duplicates.
Since Ringtail only allows one receptor per database, the receptor names in both databases will be checked, and the merge will only proceed if the receptor name matches
(as well as PRAGMA user_version). The script will delete existing bookmarks and ligand_cluster tables, as these can be dependent on primary keys which may have changed.
A log file of the process will be written by default, and the name of this log file can be specified by user.

--------- Example usage ---------
# merge with backup
python rt_merge.py -db1 db1.db -db2 db2.db

# merge without backup
python rt_merge.py -db1 db1.db -db2 db2.db --dont_backup_db1

# merge with user specified log file
python rt_merge.py -db1 db1.db -db2 db2.db -l mergetest.log

# merge with multiple databases with a common name, using wildcard (*) in the filename. Assuming you e.g., have databases db1.db, db2.db, db3.db, ..., dbn.db.
python rt_merge.py -db1 db1.db -db2 db*.db

"""
import argparse
import sys
import logging
from ringtail import RingtailCore


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-db1",
        "--primary_db",
        help="full path to primary database, of which a copy will be made, and the other database merged with",
        required=True,
        default=argparse.SUPPRESS,
    )
    parser.add_argument(
        "-db2",
        "--secondary_db",
        help="full path to secondary database that will be merged into the primary database",
        default=argparse.SUPPRESS,
        nargs="+",
    )
    parser.add_argument(
        "--dont_backup_db1",
        help="option to not back up primary database before merging with database2",
        action="store_true",
    )
    parser.add_argument(
        "--debug",
        help="if logging at debug level",
        action="store_true",
    )

    args = parser.parse_args()
    if args.debug:
        logging_level = "DEBUG"
    else:
        logging_level = "INFO"
    rtc = RingtailCore(args.primary_db, logging_level=logging_level)
    logger = logging.getLogger("RingtailLogger")

    logger.debug(str(args))
    logger.info(f"Database {args.primary_db} ready to merge.")

    # pass raw patterns — ringtailcore handles glob expansion and self-filtering
    merging_dbs = (
        args.secondary_db
        if isinstance(args.secondary_db, list)
        else [args.secondary_db]
    )

    failed = rtc.merge_databases(
        merging_dbs=merging_dbs, backup=not args.dont_backup_db1
    )
    for db, error in failed:
        print(f"FAILED TO MERGE: {db}: {error}")
        with open("failed_databases.txt", "a") as f:
            f.write(f"{db}\n")

    num_total = len(merging_dbs)
    unsuccessful = len(failed)
    num_successful = num_total - unsuccessful
    if num_successful:
        print(f"Merging of {num_successful} out of {num_total} databases complete.")
    if unsuccessful:
        print(f"{unsuccessful} databases failed to merge.")


if __name__ == "__main__":
    sys.exit(main())
