#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Convert a Ringtail v3 database's pose_coordinates from the legacy JSON-text
# format to the compact native format:
#   * DuckDB : FLOAT[][]            (ALP-compressed, ~8x smaller than the text)
#   * SQLite : packed float32 BLOB  (~5x smaller, ~65x faster to read)
#
# DuckDB cannot ALTER/DROP a column on a table referenced by foreign keys
# (Interactions -> Results.pose_id), so the DuckDB path RECREATES the database
# into a fresh file, which also compacts it ~3-5x (ideal before downloading a
# database off an HPC). SQLite converts in place, then VACUUMs to reclaim space.

import argparse
import os
import sqlite3
import time

from ringtail import setup_logging, get_logger
from ringtail.util import detect_db_type

logger = get_logger(__name__)


def _convert_sqlite(db_file: str):
    from ringtail import RingtailCore

    t0 = time.perf_counter()
    RingtailCore(db_file).convert_pose_coordinates_to_native()
    logger.info(f"convert: {time.perf_counter() - t0:.1f}s")
    t1 = time.perf_counter()
    conn = sqlite3.connect(db_file)
    conn.execute("VACUUM")
    conn.close()
    logger.info(f"VACUUM: {time.perf_counter() - t1:.1f}s")


def _convert_duckdb(db_file: str, keep_backup: bool):
    import duckdb
    from ringtail import RingtailCore
    from ringtail.schema import TABLE_SCHEMAS, RESULTS_SCHEMA, build_create_table

    probe = duckdb.connect(db_file, read_only=True)
    dtype = probe.execute(
        "SELECT data_type FROM information_schema.columns "
        "WHERE table_name = 'Results' AND column_name = 'pose_coordinates'"
    ).fetchone()
    probe.close()
    if dtype and "[]" in dtype[0]:
        logger.info("pose_coordinates is already native (FLOAT[][]); nothing to do.")
        return

    new_file = db_file + ".native.tmp"
    if os.path.exists(new_file):
        os.remove(new_file)

    t0 = time.perf_counter()
    rtc = RingtailCore(new_file)
    with rtc.storageman as sm:
        con = sm.conn
        sm._create_tables()
        try:
            sm.ensure_gui_tables()
        except Exception:
            pass
        con.execute(f"ATTACH '{db_file}' AS old (READ_ONLY)")

        src_tables = {
            r[0]
            for r in con.execute(
                "SELECT table_name FROM duckdb_tables() WHERE database_name = 'old'"
            ).fetchall()
        }
        dst_tables = {
            r[0]
            for r in con.execute(
                "SELECT table_name FROM duckdb_tables() WHERE database_name != 'old'"
            ).fetchall()
        }

        # create any source tables the standard schema didn't (clusters, merge, etc.)
        for name in src_tables - dst_tables:
            schema = TABLE_SCHEMAS.get(name.lower())
            if schema is not None:
                for stmt in build_create_table(name, schema, "duckdb"):
                    con.execute(stmt)
                dst_tables.add(name)

        rcols = list(RESULTS_SCHEMA.columns.keys())
        results_select = ", ".join(
            "CAST(pose_coordinates AS FLOAT[][]) AS pose_coordinates"
            if c == "pose_coordinates"
            else c
            for c in rcols
        )
        for name in sorted(src_tables):
            if name not in dst_tables:
                logger.warning(f"skipping {name} (no destination schema)")
                continue
            if name == "Results":
                con.execute(
                    f"INSERT INTO Results ({', '.join(rcols)}) "
                    f"SELECT {results_select} FROM old.Results"
                )
            else:
                cols = [
                    r[0]
                    for r in con.execute(
                        "SELECT column_name FROM duckdb_columns() "
                        f"WHERE database_name = 'old' AND table_name = '{name}'"
                    ).fetchall()
                ]
                collist = ", ".join(cols)
                con.execute(
                    f"INSERT INTO {name} ({collist}) SELECT {collist} FROM old.{name}"
                )
            logger.info(f"copied {name} ({time.perf_counter() - t0:.1f}s)")
        con.execute("DETACH old")
        con.execute("CHECKPOINT")

    backup = db_file + ".bak"
    os.replace(db_file, backup)
    os.replace(new_file, db_file)
    logger.info(f"recreate total: {time.perf_counter() - t0:.1f}s")
    if keep_backup:
        logger.info(f"original kept at {backup}")
    else:
        os.remove(backup)


def cmdline_parser():
    parser = argparse.ArgumentParser(
        prog="rt_convert_coordinates",
        description=(
            "Convert a Ringtail v3 database's pose_coordinates to the compact native "
            "format (DuckDB FLOAT[][]; SQLite packed float32 BLOB). DuckDB recreates "
            "and compacts the file; SQLite converts in place + VACUUM."
        ),
    )
    parser.add_argument(
        "-d", "--database", nargs="+", required=True, help="Ringtail v3 database file(s)"
    )
    parser.add_argument(
        "--keep-backup",
        action="store_true",
        help="DuckDB only: keep the original as <db>.bak instead of deleting it",
    )
    parser.add_argument("--debug", action="store_true", help="debug-level logging")
    parser.add_argument("--logfile", type=str, default=None, help="write log to FILE")
    return parser.parse_args()


def main():
    args = cmdline_parser()
    setup_logging(level="DEBUG" if args.debug else "INFO", logfile=args.logfile)
    for db_file in args.database:
        before = os.path.getsize(db_file) / 1e9
        dialect = detect_db_type(db_file)
        logger.info(f"{db_file}: {dialect}, {before:.1f} GB")
        if dialect == "duckdb":
            _convert_duckdb(db_file, args.keep_backup)
        elif dialect == "sqlite":
            _convert_sqlite(db_file)
        else:
            logger.error(f"unsupported storage type: {dialect}")
            continue
        after = os.path.getsize(db_file) / 1e9
        logger.info(f"done: {before:.1f} GB -> {after:.1f} GB")


if __name__ == "__main__":
    main()
