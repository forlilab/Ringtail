#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#

import time
import argparse
import json
import sys
import os
from ringtail import RingtailCore, setup_logging, get_logger

logger = get_logger(__name__)
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
                Written by Althea Hansel-Harris and May-Linn Paulsen. Based on code by Stefano Forli, PhD, Andreas Tillack, PhD, and Diogo Santos-Martins, PhD.\n

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
        "--logfile",
        help="Write log output to this file (useful with --debug).",
        type=str,
        metavar="FILE.log",
        default=None,
    )
    parser.add_argument(
        "--output_log",
        "-l",
        help="Name for crossref results file of passing ligands and data",
        type=str,
        metavar="<RESULTS FILE>.txt",
        const="crossref_results.txt",
        nargs="?",
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
    try:
        args = cmdline_parser()

        # set logging level
        loglvl = "DEBUG" if args.debug else "INFO" if args.verbose else "WARNING"
        setup_logging(level=loglvl, logfile=args.logfile)
        logger.info("Starting a ringtail database compare process")

        # make sure we have a positive db and at least one other database
        if args.wanted is None:
            raise IOError(
                "No wanted database found. Must specify an included database."
            )

        db_count = len(args.wanted) + (len(args.unwanted) if args.unwanted else 0)
        if db_count < 2:
            raise IOError("Must specify at least two databases for comparison.")

        wanted_dbs = []
        for db in args.wanted:
            if os.path.exists(db[0]):
                wanted_dbs.append(db)
            else:
                logger.warning(
                    f"Wanted database {db[0]} not found, removing from list."
                )

        unwanted_dbs = []
        for db in args.unwanted or []:
            if os.path.exists(db[0]):
                unwanted_dbs.append(db)
            else:
                logger.warning(
                    f"Unwanted database {db[0]} not found, removing from list."
                )

        # needs to be at least one wanted database
        if not wanted_dbs:
            raise IOError(
                "There are no valid databases in the --wanted list, cannot run cross comparison!"
            )
        # check that we have at least two databases for cross comparison
        if len(wanted_dbs) + len(unwanted_dbs) < 2:
            raise IOError(
                "There are less than two valid databases specified in the script, cannot run cross comparison!"
            )
        if len(wanted_dbs) + len(unwanted_dbs) > 10:
            raise IOError(
                "There are more than 10 specified databases, please choose 10 or less databases total. Cannot run cross comparison!"
            )
        # use first database as main connetion for cross referencing
        rtc = RingtailCore(wanted_dbs[0][0])
        num_shared_ligands, db_new_bookmarks, filter_dict = (
            rtc.cross_reference_databases(
                wanted_dbs=wanted_dbs,
                unwanted_dbs=unwanted_dbs,
                bookmark_prefix=args.save_bookmark,
            )
        )

        if num_shared_ligands == 0:
            print("No ligands found passing cross comparison.")
            return 0

        print("Number of ligands passing wanted minus unwanted: ", num_shared_ligands)

        for db, bookmark in db_new_bookmarks.items():
            rtc = RingtailCore(db)

            if args.output_log:
                output_log = str(db.split(".")[0]) + "_" + args.output_log
                logger.info(f"Writing log text file for {db} bookmark {bookmark}")
                rtc.write_filter_output(
                    bookmark,
                    filter_dict,
                    num_shared_ligands,
                    output_log,
                    output_all_poses=not args.store_best_pose,
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
        return 0

    except Exception as e:
        tb = traceback.format_exc()
        logger.debug(tb)
        logger.critical(str(e))
        sys.exit(1)


if __name__ == "__main__":
    sys.exit(main())
