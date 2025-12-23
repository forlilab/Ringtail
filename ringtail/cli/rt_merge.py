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
    )
    parser.add_argument(
        "--dont_backup_db1",
        help="option to not back up primary database before merging with database2",
        action="store_true",
    )
    parser.add_argument(
        "-l",
        "--log_filename",
        help="filename to write log",
        default="dbmerge.log",
    )

    args = parser.parse_args()
    if args.log_filename:
        log_level = "DEBUG"
    else:
        "INFO"
    rtc = RingtailCore(args.primary_db, logging_level=log_level)
    logger = logging.getLogger("RingtailLogger")
    if args.log_filename:
        file_handler = logging.FileHandler(args.log_filename)
        formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    logger.debug(str(args))
    db1 = args.primary_db
    logger.info(f"Database {db1} ready to merge. ")
    db2 = args.secondary_db
    # determine if single or multiple merging databases
    if "*" in db2:
        import glob

        merging_db = glob.glob(db2)
        # remove primary db
        if db1 in merging_db:
            merging_db.remove(db1)
            logger.info(
                "Removed the primary database from the list of merging databases."
            )
        logger.info(f"Multifile merging: {merging_db}")
    else:
        if db1 == db2:
            logger.error("Cannot merge a database with itself.")
            sys.exit(1)
        merging_db = [db2]
        logger.info(f"Single file merging: {db2}")

    # merge each database
    for index, db in enumerate(merging_db):
        try:
            if index == 0 and not args.dont_backup_db1:
                rtc.merge_databases(db, backup=True)
            else:
                rtc.merge_databases(db, backup=False)
        except Exception as e:
            logger.error(f"Database {db} failed to merge: {str(e)}.")
            with open("failed_databases.txt", "a") as f:
                f.write(f"{db}\n")
        else:
            logger.info(f"merging of {db} complete.")


if __name__ == "__main__":
    sys.exit(main())
