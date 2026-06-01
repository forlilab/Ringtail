#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail static utility methods
#

from typing import Union
from .logutils import get_logger

logger = get_logger(__name__)


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
