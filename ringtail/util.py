#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail static utility methods
#

from typing import Union
from .logutils import get_logger

logger = get_logger(__name__)


def caller_info(skip=2):
    import inspect

    """Get the name of a caller in the format module.class.method.

    https://gist.github.com/lee-pai-long/d3004225e1847b84acb4fbba0c2aea91
    Copied from: https://gist.github.com/techtonik/2151727

    :arguments:
        - skip (integer): Specifies how many levels of stack
                          to skip while getting caller name.
                          skip=1 means "who calls me",
                          skip=2 "who calls my caller" etc.

    :returns:
        - package (string): caller package.
        - module (string): caller module.
        - klass (string): caller classname if one otherwise None.
        - caller (string): caller function or method (if a class exist).
        - line (int): the line of the call.
        - An empty string is returned if skipped levels exceed stack height.
    """
    stack = inspect.stack()
    start = 0 + skip
    if len(stack) < start + 1:
        return ""
    parentframe = stack[start][0]

    # module and packagename.
    module_info = inspect.getmodule(parentframe)
    if module_info:
        mod = module_info.__name__.split(".")
        package = mod[0]
        try:
            module = mod[1]
        except IndexError:
            module = ""

    # class name.
    klass = None
    if "self" in parentframe.f_locals:
        klass = parentframe.f_locals["self"].__class__.__name__

    # method or function name.
    caller = None
    if parentframe.f_code.co_name != "<module>":  # top level usually
        caller = parentframe.f_code.co_name

    # call line.
    line = parentframe.f_lineno

    # Remove reference to frame
    # See: https://docs.python.org/3/library/inspect.html#the-interpreter-stack
    del parentframe

    return package, module, klass, caller, line


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
        bool: true if bookmark name is valid

    """
    import re

    if any(c.isupper() for c in name):
        logger.warning(
            f"Bookmark name '{name}' contains uppercase letters, converting to lowercase."
        )
        name = name.lower()

    return name if re.match(r"^[a-z0-9_]*$", name) else None


def ligand_sdf_to_pdb(sdf_file: str):

    import os
    from rdkit import Chem
    from rdkit.Chem import AllChem

    suppl = Chem.SDMolSupplier(sdf_file)
    for mol in suppl:
        name = mol.GetProp("_Name")

    if mol.GetNumConformers() == 0 or not mol.GetConformer(0).Is3D():
        mol = Chem.AddHs(mol, addCoords=True)
        AllChem.EmbedMolecule(mol)
    else:
        mol = Chem.AddHs(mol, addCoords=True)
    Chem.MolToPDBFile(mol, os.path.join(name + ".pdb"), flavor=0)


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
