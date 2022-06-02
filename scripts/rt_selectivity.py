#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#

import time
import argparse
import json
import sys
from ringtail import DBManagerSQLite
from ringtail import Outputter
from ringtail import DatabaseError, OutputError
import logging
import traceback


def cmdline_parser(defaults={}):

    conf_parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        add_help=False)
    conf_parser.add_argument(
        '-c',
        '--config',
        help=
        'specify a JSON-format file containing the option definitions. NOTE: options defined here will be overridden by command line options!'
    )
    confargs, remaining_argv = conf_parser.parse_known_args()

    defaults = {
        "positive_selection": None,
        "negative_selection": None,
        "subset_name": 'passing_results',
        "log": "selective_log.txt",
        "save_subset": None,
        "export_csv": None
    }

    config = json.loads(json.dumps(
        defaults))  # using dict -> str -> dict as a safe copy method

    if confargs.config is not None:
        logging.info("Reading options from config file")
        with open(confargs.config) as f:
            c = json.load(f)
            config.update(c)

    parser = argparse.ArgumentParser(
        usage="Please see GitHub for full usage details.",
        description=
        "Script for filtering unique passing ligands across multiple virtual screenings. Takes databases created and filtered with run_ringtail.py.",
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
        exit_on_error=False)

    parser.add_argument(
        "--positive_selection",
        '-p',
        help="Database for which ligands MUST be included in subset",
        nargs='+',
        type=str,
        metavar="[DATABASE_FILE].db",
        required=True)
    parser.add_argument(
        "--negative_selection",
        '-n',
        help="Database for which ligands MUST NOT be included in subset",
        nargs='+',
        type=str,
        metavar="[DATABASE_FILE].db",
        action='store')
    parser.add_argument(
        "--subset_name",
        '-sn',
        help=
        "Name of filter subset that ligands should be compared accross. Must be present in all databases",
        type=str,
        metavar="STRING",
        action='store')
    parser.add_argument("--log",
                        '-l',
                        help="Name for log file of passing ligands and data",
                        type=str,
                        metavar="[LOG FILE].txt",
                        action='store')
    parser.add_argument(
        "--save_subset",
        '-s',
        help=
        "Name for subset of passing cross-reference ligands to be saved as in first database given with --positive_selection",
        type=str,
        metavar="STRING",
        action='store')
    parser.add_argument(
        "--export_csv",
        '-x',
        help=
        "Save final cross-referenced subset as csv. Saved as [save_subset].csv or 'crossref.csv' if --save_subset not used.",
        action='store_true')
    parser.add_argument("--verbose",
                        '-v',
                        help="Verbose output while running",
                        action='store_true')

    parser.set_defaults(**config)
    args = parser.parse_args(remaining_argv)

    # check that name was given with save_subset
    if args.save_subset == "":
        raise IOError("--save_subset option used but no subset name given. Must specify subset name as string.")

    return args


if __name__ == '__main__':
    time0 = time.perf_counter()

    try:
        args = cmdline_parser()

        # set logging level
        debug = False
        if debug:
            level = logging.DEBUG
        elif args.verbose:
            level = logging.INFO
        else:
            level = logging.WARNING
        logging.basicConfig(level=level)

        db_opts = {
            "write_db_flag": False,
            "add_results": False,
            'num_clusters': None,
            "order_results": None,
            "log_distinct_ligands": None,
            "interaction_tolerance": None,
            "results_view_name": args.subset_name,
            "store_all_poses": None,
            "overwrite": None,
            "conflict_opt": None,
            "mode": None
        }

        positive_dbs = args.positive_selection
        negative_dbs = args.negative_selection

        ref_db = positive_dbs[0]
        positive_dbs = positive_dbs[1:]

        previous_subsetname = args.subset_name

        logging.info("Starting cross-reference process")

        dbman = DBManagerSQLite(db_opts)

        for db in positive_dbs:
            logging.info(f"cross-referencing {db}")
            previous_subsetname = dbman.crossref_filter(db,
                                                        previous_subsetname,
                                                        selection_type="+")

        for db in negative_dbs:
            logging.info(f"cross-referencing {db}")
            previous_subsetname = dbman.crossref_filter(db,
                                                        previous_subsetname,
                                                        selection_type="-")

        logging.info("Writing log")
        output_manager = Outputter(args.log)
        if args.save_subset is not None:
            output_manager.write_results_subset_to_log(args.save_subset)
        number_passing_ligands = dbman.get_number_passing_ligands()
        output_manager.log_num_passing_ligands(number_passing_ligands)
        final_subset = dbman.fetch_view(previous_subsetname)
        output_manager.write_log(final_subset)

        if args.save_subset is not None:
            dbman.save_temp_subset(args.save_subset)

        if args.export_csv:
            if args.save_subset is not None:
                csv_name = args.save_subset + ".csv"
            else:
                csv_name = "crossref.csv"
            output_manager.export_csv(previous_subsetname, csv_name, True)

    except Exception as e:
        tb = traceback.format_exc()
        logging.debug(tb)
        logging.critical(e)
        sys.exit(1)

    finally:
        dbman.close_db_crossref()
