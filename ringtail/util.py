#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail static utility methods
#

import os
import shutil
import subprocess
from typing import Union
from .logutils import get_logger

logger = get_logger(__name__)

# file extension per compression method
_COMPRESS_EXT = {"zstd": ".zst", "gzip": ".gz", "xz": ".xz"}


def numlist2str(list: list, separator: str) -> str:
    """
    Joines item in a list by specified string separator

    Args:
        list (list): list to be joined
        separator (str): string item to separate the items in the list

    Returns:
        str: list as a string separated by separator
    """
    return separator.join([str(x) for x in list])


def iterate_nested(obj):
    """
    File inputs can come in multiple levels of nested lists, this method unpacks them

    Args:
        obj (list[list[list[etc]]]): None or nested lists

    Returns:
        None: if input is None

    Yields:
        str: should be unpacked paths to docking results
    """
    if obj is None:
        return None
    elif isinstance(obj, list):
        for item in obj:
            yield from iterate_nested(item)
    else:
        yield obj


def valid_bookmark_name(name: str) -> Union[str, None]:
    """Checks that bookmark name adheres to sqlite naming conventions of alphanumerical and limited symbols.

    Args:
        name (str): bookmark name

    Returns:
        str, None: bookmark name if valid else None

    """
    import re

    if any(c.isupper() for c in name):
        logger.warning(
            f"Bookmark name '{name}' contains uppercase letters, converting to lowercase."
        )
        name = name.lower()

    return name if re.match(r"^[a-z0-9_]*$", name) else None


def detect_db_type(filepath: str) -> str:
    """
    Attempts to detect storage type of an existing database file

    Args:
        filepath (str): _description_

    Raises:
        TypeError: _description_

    Returns:
        str: _description_
    """
    with open(filepath, "rb") as f:
        header = f.read(16)
    if header.startswith(b"SQLite"):
        return "sqlite"
    elif header.find(b"DUCK") != -1:
        return "duckdb"
    else:
        try:
            # check if empty sqlite file
            import sqlite3

            sqlite3.connect(filepath)
            return "sqlite"
        except Exception:
            raise TypeError(f"Database type not recognized in file {filepath}")


def compress_file(
    src: str, dst: str = None, method: str = "zstd", level: int = 18
) -> str:
    """Compress a file to a single artifact, leaving ``src`` untouched.

    Args:
        src: path to the file to compress (never modified, moved, or deleted).
        dst: output artifact path. If None, ``src`` + the method's extension.
        method: "zstd" (default), "gzip", or "xz".
        level: compression level (zstd 1-22, gzip 1-9, xz 0-9).

    Returns:
        str: the artifact path actually written (extension reflects the method
        actually used — note that "zstd" falls back to "gzip" if the zstd binary
        is unavailable).
    """
    if method not in _COMPRESS_EXT:
        raise ValueError(f"Unknown compression method '{method}'")

    # zstd is a multithreaded shell-out; fall back to stdlib gzip if absent
    if method == "zstd" and shutil.which("zstd") is None:
        logger.warning("zstd binary not found on PATH — falling back to gzip.")
        method = "gzip"

    # normalize dst: strip any known compression suffix, then add the one for the
    # method actually used (so an explicit dst + gzip fallback can't double-suffix)
    if dst is None:
        dst = src
    for _ext in _COMPRESS_EXT.values():
        if dst.endswith(_ext):
            dst = dst[: -len(_ext)]
            break
    dst += _COMPRESS_EXT[method]

    logger.info(f"Compressing {src} -> {dst} ({method} -{level})")
    if method == "zstd":
        # -T0 = all cores; -f overwrite; --long for better ratio on large files
        subprocess.run(
            ["zstd", f"-{level}", "-T0", "--long=27", "-f", "-q", "-o", dst, src],
            check=True,
        )
    elif method == "gzip":
        import gzip

        with open(src, "rb") as fin, gzip.open(dst, "wb", compresslevel=min(level, 9)) as fout:
            shutil.copyfileobj(fin, fout, length=1024 * 1024)
    elif method == "xz":
        import lzma

        with open(src, "rb") as fin, lzma.open(dst, "wb", preset=min(level, 9)) as fout:
            shutil.copyfileobj(fin, fout, length=1024 * 1024)

    return dst


def decompress_file(src: str, dst: str = None) -> str:
    """Decompress an artifact produced by :func:`compress_file`.

    Method is inferred from the extension (.zst/.gz/.xz). ``src`` is left intact.

    Args:
        src: compressed artifact path.
        dst: output path. If None, ``src`` with the compression suffix removed.

    Returns:
        str: the decompressed file path.
    """
    ext = os.path.splitext(src)[1].lower()
    if dst is None:
        dst = src[: -len(ext)] if ext in (".zst", ".gz", ".xz") else src + ".out"

    logger.info(f"Decompressing {src} -> {dst}")
    if ext == ".zst":
        if shutil.which("zstd") is None:
            raise RuntimeError(
                "zstd binary not found on PATH; cannot decompress a .zst artifact. "
                "Install zstd or decompress it on a machine that has it."
            )
        subprocess.run(["zstd", "-d", "-f", "-q", "-o", dst, src], check=True)
    elif ext == ".gz":
        import gzip

        with gzip.open(src, "rb") as fin, open(dst, "wb") as fout:
            shutil.copyfileobj(fin, fout, length=1024 * 1024)
    elif ext == ".xz":
        import lzma

        with lzma.open(src, "rb") as fin, open(dst, "wb") as fout:
            shutil.copyfileobj(fin, fout, length=1024 * 1024)
    else:
        raise ValueError(f"Unrecognized compression extension '{ext}' for {src}")

    return dst


def db_alias_from_path(db_path: str) -> str:
    import re
    import pathlib

    db_name = pathlib.Path(db_path).name.split(".")[0]
    # can only have underscpres
    sanitized = re.sub(r"[^a-zA-Z0-9_]", "_", db_name)
    # cannot start with a number
    if sanitized[0].isdigit():
        sanitized = "_" + sanitized

    return sanitized
