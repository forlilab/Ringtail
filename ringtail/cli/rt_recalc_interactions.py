#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail script for recalculating interactions in an existing database
#

import argparse
from ringtail import RingtailCore, setup_logging, get_logger
from ringtail.ringtailoptions import RingtailDefaults

logger = get_logger(__name__)
import sys
import traceback
import time


def cmdline_parser():

    parser = argparse.ArgumentParser(
        prog="rt_recalc_interactions",
        description=(
            "Recalculate the interactions in one or more existing Ringtail databases, "
            "optionally at new cutoffs. Existing interactions are deleted and computed "
            "again from the stored poses and receptor, so no docking files are needed. "
            "Work is committed in batches and a database that is interrupted resumes "
            "where it stopped, which is what makes this usable on large databases and "
            "under HPC job time limits."
        ),
    )
    parser.add_argument(
        "-d",
        "--database",
        help="Ringtail database file(s) to recalculate interactions for",
        nargs="+",
        type=str,
        action="store",
        required=True,
    )
    parser.add_argument(
        "--hb_cutoff",
        help=f"Hydrogen bond distance cutoff in angstroms (default {RingtailDefaults.interaction_cutoffs[0]})",
        default=RingtailDefaults.interaction_cutoffs[0],
        type=float,
        action="store",
    )
    parser.add_argument(
        "--vdw_cutoff",
        help=f"Van der Waals distance cutoff in angstroms (default {RingtailDefaults.interaction_cutoffs[1]})",
        default=RingtailDefaults.interaction_cutoffs[1],
        type=float,
        action="store",
    )
    parser.add_argument(
        "--chunk_size",
        help=(
            "Poses per commit. Smaller means more frequent checkpoints to resume from "
            f"and less memory (default {RingtailDefaults.chunk_size // 10})"
        ),
        default=RingtailDefaults.chunk_size // 10,
        type=int,
        action="store",
    )
    parser.add_argument(
        "-y",
        "--yes",
        help=(
            "Skip the confirmation prompt. Required for unattended and batch-scheduler "
            "runs, where there is no terminal to answer it"
        ),
        action="store_true",
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

        if args.chunk_size < 1:
            raise ValueError(f"--chunk_size must be at least 1, got {args.chunk_size}")

        logger.warning(
            "WARNING: existing interactions in the database will be deleted and "
            "recalculated. Results.num_hb and Results.num_interactions change with "
            "them, so anything filtered or bookmarked on those is worth revisiting."
        )
        consent = args.yes or input("Type 'yes' if you wish to continue: ") == "yes"
        if not consent:
            logger.info("Consent not given, nothing was changed.")
            return 0

        failed = []
        for db in args.database:
            logger.info(f"Recalculating interactions for {db}")
            try:
                rtcore = RingtailCore(db)
                rtcore.add_interactions(
                    hb_cutoff=args.hb_cutoff,
                    vdw_cutoff=args.vdw_cutoff,
                    consent=True,
                    chunk_size=args.chunk_size,
                )
            except Exception as e:
                # One unreadable database should not abandon the rest of the batch;
                # an interrupted one resumes on the next run either way.
                logger.error(f"Could not recalculate interactions for {db}: {e}")
                logger.debug(traceback.format_exc())
                failed.append(db)

        logger.info(
            "Time to recalculate interactions: {0:.0f} seconds".format(
                time.perf_counter() - time0
            )
        )
        if failed:
            logger.critical(f"Failed for {len(failed)} database(s): {', '.join(failed)}")
            return 1

        return 0

    except Exception as e:
        tb = traceback.format_exc()
        logger.debug(tb)
        logger.critical(str(e))
        return 1


if __name__ == "__main__":
    sys.exit(main())
