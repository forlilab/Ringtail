#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail script for updating v1.1.0 databases to v2.0.0
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
        prog="rt_db_to_v200",
        description="Given one or multiple Ringtail databases made with v1.1.0, will update them to be compatible with v2.0.0",
    )

    parser.add_argument(
        "-d",
        "--database",
        help="Database file(s) made with Ringtail v1.0.0 or v1.1.0 to update to v2.0.0",
        nargs="+",
        type=str,
        action="store",
    )
    args = parser.parse_args()

    consent = False

    for db in args.database:
        # create a new ringtailcore for each db file to be converted
        rtcore = RingtailCore(db)
        consent = rtcore.update_database_version(consent, new_version="2.0.0")
    return


if __name__ == "__main__":
    sys.exit(main())
