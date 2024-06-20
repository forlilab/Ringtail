#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail static utility methods
#

# region #-#-#- Util method -#-#-#


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


# endregion
