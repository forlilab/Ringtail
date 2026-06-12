#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Unpack a Ringtail database artifact produced by rt_compress_db. The compression
# method is inferred from the extension (.zst / .gz / .xz). The compressed input
# artifact is left intact.

import argparse
import os

from ringtail import setup_logging, get_logger
from ringtail.util import decompress_file, detect_db_type

logger = get_logger(__name__)


def cmdline_parser():
    parser = argparse.ArgumentParser(
        prog="rt_decompress_db",
        description=(
            "Decompress a Ringtail database artifact (.zst/.gz/.xz) produced by "
            "rt_compress_db back into a usable database file."
        ),
    )
    parser.add_argument(
        "-i", "--input", required=True, help="compressed artifact (.zst/.gz/.xz)"
    )
    parser.add_argument(
        "-o",
        "--output",
        default=None,
        help="output database path (default: input with the compression suffix removed)",
    )
    parser.add_argument("--debug", action="store_true", help="debug-level logging")
    parser.add_argument("--logfile", type=str, default=None, help="write log to FILE")
    return parser.parse_args()


def main():
    args = cmdline_parser()
    setup_logging(level="DEBUG" if args.debug else "INFO", logfile=args.logfile)

    in_size = os.path.getsize(args.input) / 1e9
    logger.info(f"{args.input}: {in_size:.2f} GB")
    db_path = decompress_file(args.input, args.output)

    out_size = os.path.getsize(db_path) / 1e9
    try:
        dialect = detect_db_type(db_path)
    except Exception as e:
        logger.warning(f"Decompressed file did not look like a Ringtail database: {e}")
        dialect = "unknown"
    logger.info(f"Wrote {db_path}: {dialect}, {out_size:.2f} GB")


if __name__ == "__main__":
    main()
