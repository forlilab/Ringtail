#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail script for updating v1.0.0 databases to v1.1.0
#

import argparse
from ringtail import RingtailCore
import logging
import sys


def main():
    logging.basicConfig(
        level=logging.INFO, stream=sys.stdout, filemode="w", format="%(message)s"
    )
    # get name(s) of dbs to update from command line
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
    )
    parser.add_argument(
        "-v",
        "--version",
        help="Version to upgrade to, either 1.1.0 or 2.0.0",
        default="3.0.0",
        type=str,
        action="store",
    )
    args = parser.parse_args()

    consent = False
    version = args.version
    # validate version
    if version not in ["1.1.0", "2.0.0", "3.0.0"]:
        print(
            f"Requested upgrade version (-{version}-) not recognized. Please choose either -1.1.0-, -2.0.0-, or the default -3.0.0-."
        )

    for db in args.database:
        rtcore = RingtailCore(db)
        consent = rtcore.update_database_version(consent, version)
    return


if __name__ == "__main__":
    sys.exit(main())
