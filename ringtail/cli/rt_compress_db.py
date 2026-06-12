#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Prepare a Ringtail database for transfer off an HPC: optionally filter it down
# by docking score and/or ligand efficiency to a minimal, self-contained database,
# then compress it as small as practical.
#
#   * With -e/--eworst and/or -le/--leworst: writes a NEW database containing only
#     the passing ligands/poses (via export_bookmark_db), then compresses it.
#   * Without thresholds: compresses the input database as-is.
#
# The input/original database is NEVER deleted, moved, or overwritten. With a
# filter, a temporary bookmark is added to the source and removed again (so the
# source must be writable); its data is left unchanged. Unpack with
# rt_decompress_db.

import argparse
import os

from ringtail import RingtailCore, setup_logging, get_logger
from ringtail.util import compress_file, detect_db_type

logger = get_logger(__name__)

_TMP_BOOKMARK = "rt_compress_tmp"


def cmdline_parser():
    parser = argparse.ArgumentParser(
        prog="rt_compress_db",
        description=(
            "Filter a Ringtail database by docking score / ligand efficiency (optional) "
            "into a minimal self-contained database, then compress it for transfer. "
            "The original database is never modified destructively."
        ),
    )
    parser.add_argument("-i", "--input", required=True, help="input Ringtail database")
    parser.add_argument(
        "-o",
        "--output",
        default=None,
        help="output artifact path (default derived from input + method extension)",
    )
    parser.add_argument(
        "-e", "--eworst", type=float, default=None, help="worst docking score to keep"
    )
    parser.add_argument(
        "-le",
        "--leworst",
        type=float,
        default=None,
        help="worst ligand efficiency to keep",
    )
    parser.add_argument(
        "--compressor",
        choices=["zstd", "gzip", "xz"],
        default="zstd",
        help="compression method (default zstd; falls back to gzip if zstd missing)",
    )
    parser.add_argument(
        "--level",
        type=int,
        default=18,
        help="compression level (default 18; ~3x smaller, good speed/ratio balance)",
    )
    parser.add_argument(
        "--keep-db",
        action="store_true",
        help="keep the intermediate uncompressed filtered database",
    )
    parser.add_argument("--debug", action="store_true", help="debug-level logging")
    parser.add_argument("--logfile", type=str, default=None, help="write log to FILE")
    return parser.parse_args()


def main():
    args = cmdline_parser()
    setup_logging(level="DEBUG" if args.debug else "INFO", logfile=args.logfile)

    in_size = os.path.getsize(args.input) / 1e9
    logger.info(f"{args.input}: {detect_db_type(args.input)}, {in_size:.2f} GB")

    filtering = args.eworst is not None or args.leworst is not None
    created_subset = None  # path of the intermediate db we made (to clean up)

    if filtering:
        # default subset path: <input stem>_filtered.db (must not already exist)
        if args.output is not None:
            base = args.output
            for ext in (".zst", ".gz", ".xz"):
                if base.endswith(ext):
                    base = base[: -len(ext)]
                    break
            subset_db = base
        else:
            subset_db = args.input.removesuffix(".db") + "_filtered.db"
        if os.path.exists(subset_db):
            raise SystemExit(f"Refusing to overwrite existing file: {subset_db}")

        rtc = RingtailCore(args.input)
        filt = {}
        if args.eworst is not None:
            filt["eworst"] = args.eworst
        if args.leworst is not None:
            filt["leworst"] = args.leworst
        logger.info(f"Filtering {args.input} with {filt} ...")
        num_passing, _ = rtc.filter(output_bookmark=_TMP_BOOKMARK, **filt)

        if not num_passing:
            rtc.delete_bookmark(_TMP_BOOKMARK)  # clean up the empty bookmark
            raise SystemExit("No ligands passed the filter; nothing to export.")
        logger.info(f"{num_passing} ligands passed; writing filtered database ...")

        # build the minimal subset DB, then drop the temp bookmark from the source
        out_db = rtc.export_bookmark_db(_TMP_BOOKMARK, subset_db)
        rtc.delete_bookmark(_TMP_BOOKMARK)
        if not out_db:
            raise SystemExit(f"Export failed (target may exist): {subset_db}")
        created_subset = out_db
        db_to_compress = out_db
    else:
        logger.info("No filter given; compressing the database as-is (read-only).")
        db_to_compress = args.input

    artifact = compress_file(
        db_to_compress, args.output, method=args.compressor, level=args.level
    )

    if created_subset and not args.keep_db:
        os.remove(created_subset)

    out_size = os.path.getsize(artifact) / 1e9
    logger.info(
        f"Wrote {artifact}: {out_size:.2f} GB "
        f"({in_size / out_size:.1f}x smaller than the {in_size:.2f} GB input). "
        f"Unpack with: rt_decompress_db -i {artifact}"
    )


if __name__ == "__main__":
    main()
