#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail script for updating older databases to v3.0.0
#

import argparse
from ringtail import RingtailCore, setup_logging, get_logger

logger = get_logger(__name__)
import sys
import traceback
import time


def cmdline_parser():

    parser = argparse.ArgumentParser(
        prog="rt_upgrade_db",
        description="Given one or multiple Ringtail databases made with older versions (e.g., 1.1.0, 2.0.0), will upgrade them to the latest version 3.0.0, unless given a specific version (e.g., intermediate 2.0.0). Can only upgrade, not downgrade.",
    )
    parser.add_argument(
        "-d",
        "--database",
        help="Database file(s) made with older Ringtail versions (1.0.0, 1.1.0, 2.0.0)",
        nargs="+",
        type=str,
        action="store",
        required=True,
    )
    parser.add_argument(
        "--version",
        help="Version to upgrade to, eg 3.0.0 (default)",
        default="3.0.0",
        type=str,
        action="store",
    )
    parser.add_argument(
        "--debug",
        help="if logging at debug level",
        action="store_true",
    )
    parser.add_argument(
        "--logfile",
        help="Write log output to this file (useful with --debug).",
        type=str,
        metavar="FILE.log",
        default=None,
    )

    args = parser.parse_args()

    return args


def main():
    time0 = time.perf_counter()
    try:
        args = cmdline_parser()
        setup_logging(level="DEBUG" if args.debug else "INFO", logfile=args.logfile)

        version = args.version
        # validate version
        if version not in ["1.1.0", "2.0.0", "3.0.0"]:
            raise ValueError(
                f"Requested upgrade version {version} not recognized. Please choose 1.1.0, 2.0.0, or 3.0.0."
            )

        logger.warning(
            "WARNING: All existing filters and bookmarks in database will be dropped during database update!"
        )
        consent = input("Type 'yes' if you wish to continue: ") == "yes"

        for db in args.database:
            rtcore = RingtailCore(db)
            rtcore.update_database_version(consent=consent, new_version=version)

        logger.info(
            "Time to upgrade: {0:.0f} seconds".format(time.perf_counter() - time0)
        )

        return 0

    except Exception as e:
        tb = traceback.format_exc()
        logger.debug(tb)
        logger.critical(str(e))
        return 1


if __name__ == "__main__":
    sys.exit(main())
