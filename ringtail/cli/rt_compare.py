#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#

import time
import argparse
import json
import sys
import os
from ringtail import RingtailCore
from ringtail import OutputManager
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
        "log": None,
        "save_bookmark": "crossref",
        "export_sdf": None,
        "export_db": None,
    }

    config = json.loads(
        json.dumps(defaults)
    )  # using dict -> str -> dict as a safe copy method

    if confargs.config is not None:
        logutils.LOGGER.info("Reading options from config file")
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
        help="Database and bookmark name for which ligands MUST be included in",
        nargs=2,
        type=str,
        action="append",
        metavar="<DATABASE_FILE>.db <bookmark_name>",
    )
    parser.add_argument(
        "--unwanted",
        "-uw",
        help="Database and bookmark name for which ligands MUST NOT be included in",
        nargs=2,
        type=str,
        metavar="<DATABASE_FILE>.db <bookmark_name>",
        action="append",
    )
    parser.add_argument(
        "--verbose", "-v", help="Verbose output while running", action="store_true"
    )
    parser.add_argument(
        "--debug",
        "-d",
        help="Debug logging (includes log file) output while running",
        action="store_true",
    )
    parser.add_argument(
        "--log",
        "-l",
        help="Name for crossref log file of passing ligands and data",
        type=str,
        metavar="<LOG FILE>.txt",
        action="store",
    )
    parser.add_argument(
        "--save_bookmark",
        "-s",
        help="Prefix for bookmark of passing cross-reference ligands to be saved in each database, defaults to crossref_<compared_bookmark>",
        type=str,
        metavar="STRING",
        action="store",
    )
    parser.add_argument(
        "--export_sdf",
        "-xs",
        help="Exports all crossreferenced ligands to one SDF per database. File name will be the compared bookmark name prefixed with chosen phrase, defaults to crossref_<bookmark>.sdf",
        const="crossref",
        nargs="?",
        action="store",
    )
    parser.add_argument(
        "--store_best_pose",
        "-bp",
        help="Will only store in bookmark and export the best ranked pose for each ligand.",
        action="store_true",
    )
    parser.add_argument(
        "--export_db",
        "-xd",
        help="Exports all crossreferenced ligands to one new database per database. File name will be the compared bookmark name prefixed with chosen phrase, defaults to crossref_<bookmark>.db",
        action="store_true",
    )

    parser.set_defaults(**config)
    args = parser.parse_args(remaining_argv)

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
        if args.verbose:
            logger.set_level("INFO")
        elif args.debug:
            logger.set_level("DEBUG")
        else:
            logger.set_level("WARNING")

        wanted_dbs = args.wanted
        unwanted_dbs = args.unwanted
        if not unwanted_dbs:
            unwanted_dbs = []
        filter_dict = {"wanted": wanted_dbs, "unwanted": unwanted_dbs}
        ### Ensure we have (enough) valid database files before proceeding
        # check that each wanted database exists
        for index, database_bookmark in enumerate(reversed(wanted_dbs)):
            if not os.path.exists(database_bookmark[0]):
                logger.warning(
                    f"Wanted database {database_bookmark[0]} not found, removing from list!"
                )
                del wanted_dbs[index]
        # needs to be at least one wanted database
        if not wanted_dbs:
            logger.critical(
                "There are no valid databases in the --wanted list, cannot run cross comparison!"
            )
        # check that each unwanted database exists
        for index, database_bookmark in enumerate(reversed(unwanted_dbs)):
            if not os.path.exists(database_bookmark[0]):
                logger.warning(
                    f"Unwanted database {database_bookmark[0]} not found, removing from list!"
                )
                del unwanted_dbs[index]
        # check that we have at least two databases for cross comparison
        if len(wanted_dbs) + len(unwanted_dbs) < 2:
            logger.critical(
                "There are less than two valid databases specified in the script, cannot run cross comparison!"
            )
        if len(wanted_dbs) + len(unwanted_dbs) > 10:
            logger.critical(
                "There are more than 10 specified databases, please choose 10 or less databases total. Cannot run cross comparison!"
            )
        # use first database as main connetion for cross referencing
        print("requested wanted databases: ", wanted_dbs)
        rtc = RingtailCore(wanted_dbs[0][0])
        num_shared_ligands, db_new_bookmarks = rtc.cross_reference_databases(
            wanted_dbs=wanted_dbs,
            unwanted_dbs=unwanted_dbs,
            bookmark_prefix=args.save_bookmark,
            store_best_pose=args.store_best_pose,
        )
        print("Number of ligands passing wanted minus unwanted: ", num_shared_ligands)
        if num_shared_ligands == 0:
            logger.warning("No ligands found passing cross comparison.")
        else:
            # output options
            for db, bookmark in db_new_bookmarks.items():
                rtc = RingtailCore(db)

                if args.log:
                    log_file = str(db.split(".")[0]) + "_" + args.log
                    logger.info(f"Writing log text file for {db} bookmark {bookmark}")
                    rtc._write_filter_output(
                        bookmark,
                        "wanted, unwanted",
                        num_shared_ligands,
                        log_file,
                        output_all_poses=True,
                    )
                    logger.info(f"Wrote {db} bookmark {bookmark} results to log file.")

                if args.export_sdf:
                    rtc.write_molecule_sdfs(bookmark_name=bookmark)
                    logger.info(f"Exported {db} bookmark {bookmark} to SDF.")

                if args.export_db:
                    rtc.export_bookmark_db(bookmark)
                    logger.info(f"Exported {db} bookmark {bookmark} to a new database.")

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
        logger.error("Error encountered while cross-referencing, please check inputs!")
        sys.exit(1)
    return


if __name__ == "__main__":
    sys.exit(main())
