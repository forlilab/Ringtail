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
