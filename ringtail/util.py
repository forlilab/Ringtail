#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail static utility methods
#


def split_dict(dict: dict, items: list) -> tuple:
    """Utility method that takes one dictionary and splits it into two based on the listed keys

    Args:
        dict (dict): original dictionary
        items (list): ist of keys to use for separation

    Returns:
        tuple: original dict minus the removed items and new dict containing the items removed from the original dict
    """

    new_dict = {}

    for key in items:
        new_dict[key] = dict.pop(key)

    return dict, new_dict


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
        except:
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


def valid_bookmark_name(name) -> bool:
    """Checks that bookmark name adheres to sqlite naming conventions of alphanumerical and limited symbols.

    Args:
        name (str): bookmark name

    Returns:
        bool: true if bookmark name is valid

    """
    import re

    regex = "^[A-Za-z0-9_]*$"
    return re.match(regex, name)


docking_mode_file_ext = {"dlg": "dlg", "vina": "pdbqt"}


docking_mode_aliases = {"dlg": ["gpu", "adgpu", "dlg"], "vina": ["vina"]}


statuses = ["accepted", "maybe", "rejected"]


def generate_not_implemented_message():
    import inspect

    frame = inspect.currentframe()
    outer_frame = frame.f_back
    method_name = outer_frame.f_code.co_name
    class_name = None

    # Try to find the class name from 'self'
    if "self" in outer_frame.f_locals:
        class_name = outer_frame.f_locals["self"].__class__.__name__

    if class_name:
        return (
            f"Method '{method_name}' must be implemented in subclass of '{class_name}'."
        )
    else:
        return f"Method '{method_name}' must be implemented in subclass."


def raise_not_implemented():
    raise NotImplementedError(generate_not_implemented_message())


def ligand_sdf_to_pdb(sdf_file: str):

    import os
    from rdkit import Chem
    from rdkit.Chem import AllChem

    suppl = Chem.SDMolSupplier(sdf_file)
    for mol in suppl:
        name = mol.GetProp("_Name")

    if mol.GetNumConformers() == 0 or not mol.GetConformer(0).Is3D():
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
    else:
        mol = Chem.AddHs(mol, addCoords=True)
    Chem.MolToPDBFile(mol, os.path.join(name + ".pdb"), flavor=0)
