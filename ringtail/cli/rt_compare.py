#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#

import time
import argparse
import json
import sys
import os
from ringtail import StorageManager, StorageManagerSQLite
from ringtail import OutputManager, OptionError
from ringtail import OptionError
from ringtail import logutils
import traceback


def cmdline_parser(defaults={}):

    conf_parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        add_help=False,
    )
    conf_parser.add_argument(
        "-c",
        "--config",
        help="specify a JSON-format file containing the option definitions. NOTE: options defined here will be overridden by command line options!",
    )
    confargs, remaining_argv = conf_parser.parse_known_args()

    defaults = {
        "wanted": None,
        "unwanted": None,
        "bookmark_name": ["passing_results"],
        "log": "selective_log.txt",
        "save_bookmark": None,
        "export_csv": None,
    }

    config = json.loads(
        json.dumps(defaults)
    )  # using dict -> str -> dict as a safe copy method

    if confargs.config is not None:
        logger.info("Reading options from config file")
        with open(confargs.config) as f:
            c = json.load(f)
            config.update(c)

    parser = argparse.ArgumentParser(
        usage="Please see GitHub for full usage details.",
        description="Script for filtering unique passing ligands across multiple virtual screenings. Takes databases created and filtered with rt_process_vs.py.",
        epilog="""

        REQUIRED PACKAGES
                Requires RDkit, SciPy, Meeko.\n

        AUTHOR
                Written by Althea Hansel-Harris. Based on code by Stefano Forli, PhD, Andreas Tillack, PhD, and Diogo Santos-Martins, PhD.\n

        REPORTING BUGS
                Please report bugs to:
                AutoDock mailing list   http://autodock.scripps.edu/mailing_list\n

        COPYRIGHT
                Copyright (C) 2022 Stefano Forli Laboratory, Center for Computational Structural Biology,
                             The Scripps Research Institute.
                GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>
        """,
        exit_on_error=False,
    )

    parser.add_argument(
        "--wanted",
        "-w",
        help="Database for which ligands MUST be included in bookmark",
        nargs="+",
        type=str,
        metavar="[DATABASE_FILE].db",
    )
    parser.add_argument(
        "--unwanted",
        "-n",
        help="Database for which ligands MUST NOT be included in bookmark",
        nargs="+",
        type=str,
        metavar="[DATABASE_FILE].db",
        action="store",
    )
    parser.add_argument(
        "--bookmark_name",
        "-sn",
        help="Name of filter bookmark that ligands should be compared accross. Must be present in all databases",
        type=str,
        metavar="STRING",
        action="store",
        nargs="+",
    )
    parser.add_argument(
        "--log",
        "-l",
        help="Name for log file of passing ligands and data",
        type=str,
        metavar="[LOG FILE].txt",
        action="store",
    )
    parser.add_argument(
        "--save_bookmark",
        "-s",
        help="Name for bookmark of passing cross-reference ligands to be saved as in first database given with --wanted",
        type=str,
        metavar="STRING",
        action="store",
    )
    parser.add_argument(
        "--export_csv",
        "-x",
        help="Save final cross-referenced bookmark as csv. Saved as [save_bookmark].csv or 'crossref.csv' if --save_bookmark not used.",
        action="store_true",
    )
    parser.add_argument(
        "--verbose", "-v", help="Verbose output while running", action="store_true"
    )
    parser.add_argument(
        "--database_type",
        "-dbt",
        help="Specify what database type was used to construct the databases. NOTE: All databases have to be created with the same type.",
        action="store",
    )

    parser.set_defaults(**config)
    args = parser.parse_args(remaining_argv)

    # check that name was given with save_bookmark
    if args.save_bookmark == "":
        raise IOError(
            "--save_bookmark option used but no bookmark name given. Must specify bookmark name as string."
        )

    return args

def main():
    time0 = time.perf_counter()
    logger = logutils.LOGGER
    logger.info("Starting a ringtail database compare process")
    try:
        args = cmdline_parser()

        # make sure we have a positive db and at least one other database
        if args.wanted is None:
            raise IOError(
                "No wanted database found. Must specify an included database."
            )
        else:
            db_count = len(args.wanted)
            if args.unwanted is not None:
                db_count += len(args.unwanted)
            if db_count < 2:
                raise IOError("Must specify at least two databases for comparison.")

        # set logging level
        debug = True
        if debug:
            level = logger.set_level("DEBUG")
        elif args.verbose:
            level = logger.set_level("INFO")
        else:
            level = logger.set_level("WARNING")
        logger.add_filehandler("ringtail", level)
        wanted_dbs = args.wanted
        unwanted_dbs = args.unwanted

        # check that ref database exists
        if not os.path.exists(wanted_dbs[0]):
            logger.critical("Wanted database {0} not found!".format(wanted_dbs[0]))

        ref_db = wanted_dbs[0]
        wanted_dbs = wanted_dbs[1:]

        previous_bookmarkname = args.bookmark_name[0]
        original_bookmark_name = args.bookmark_name[0]
        # make list of bookmark names for compared databases
        bookmark_list = []
        if len(args.bookmark_name) > 1:
            num_dbs = len(args.wanted)
            if args.unwanted is not None:
                num_dbs += len(args.unwanted)
            if len(args.bookmark_name) != num_dbs:
                raise OptionError(
                    "Not enough bookmark names given for the number of databases given. If specifying more than one bookmark name, must give bookmark name for every database to be compared."
                )
            bookmark_list = args.bookmark_name[1:]
        else:
            for db in wanted_dbs:
                bookmark_list.append(original_bookmark_name)
            if unwanted_dbs is not None:
                for db in unwanted_dbs:
                    bookmark_list.append(original_bookmark_name)

        logger.info("Starting cross-reference process")

        # set database type to default sqlite if not specified in cmdline args
        if args.database_type is None:
            args.database_type = "sqlite"
        # ensure the specified database types are supported
        storagemanager = StorageManager.check_storage_compatibility(args.database_type)

        dbman = storagemanager(ref_db)

        last_db = None
        num_wanted_dbs = len(wanted_dbs)
        # storageman is a context manager, and keeps connection to the database open within the `with` statement
        with dbman:
            for idx, db in enumerate(wanted_dbs):
                logger.info(f"cross-referencing {db}")
                if not os.path.exists(db):
                    logger.critical("Wanted database {0} not found!".format(db))
                print(
                    f"passing in these parameters: {db}, {previous_bookmarkname}, {bookmark_list[idx]}, {last_db}"
                )
                previous_bookmarkname, number_passing_ligands = dbman.crossref_filter(
                    db,
                    previous_bookmarkname,
                    bookmark_list[idx],
                    selection_type="+",
                    old_db=last_db,
                )
                last_db = db

            if unwanted_dbs is not None:
                for idx, db in enumerate(unwanted_dbs):
                    logger.info(f"cross-referencing {db}")
                    if not os.path.exists(db):
                        logger.critical("Unwanted database {0} not found!".format(db))
                    previous_bookmarkname, number_passing_ligands = (
                        dbman.crossref_filter(
                            db,
                            previous_bookmarkname,
                            bookmark_list[idx + num_wanted_dbs],
                            selection_type="-",
                            old_db=last_db,
                        )
                    )
                    last_db = db

            logger.info("Writing log")
            output_manager = OutputManager(log_file=args.log)
            with output_manager:
                if args.save_bookmark is not None:
                    output_manager.write_results_bookmark_to_log(args.save_bookmark)
                output_manager.log_num_passing_ligands(number_passing_ligands)
                final_bookmark = dbman.fetch_bookmark(previous_bookmarkname)
                output_manager.write_filter_log(final_bookmark)

            if args.save_bookmark is not None:
                dbman.create_bookmark_from_temp_table(
                    previous_bookmarkname,
                    args.save_bookmark,
                    original_bookmark_name,
                    args.wanted,
                    args.unwanted,
                )

            if args.export_csv:
                if args.save_bookmark is not None:
                    csv_name = args.save_bookmark + ".csv"
                else:
                    csv_name = "crossref.csv"
                dataframe = dbman.to_dataframe(previous_bookmarkname, table=True)
                dataframe.to_csv(csv_name)
                logger.info("Exported bookmark to csv")

        logger.info(
            "Time to cross-reference: {0:.0f} seconds".format(
                time.perf_counter() - time0
            )
        )
        logger.info("Database comparison process complete.")

    except Exception as e:
        tb = traceback.format_exc()
        logger.debug(tb)
        logger.critical(str(e))
        logger.error(
            "Error encountered while cross-referencing. If error states 'Error while getting number of passing ligands', please confirm that given bookmark names are correct."
        )
        sys.exit(1)
    return

if __name__ == "__main__":
    sys.exit(main())