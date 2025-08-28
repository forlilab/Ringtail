#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail options handler
#

import os
from .exceptions import OptionError
from .logutils import LOGGER as logger
from .util import iterate_nested
from dataclasses import dataclass, asdict
import copy


class TypeSafe:
    """Class that handles safe typesetting of values of a specified built-in type.
    Any attribute can be set as a TypeSafe object, this ensures its type is checked whenever it is changed.
    This makes the attribute of type 'object' as opposed to its actual type. To return the value of an attribute as
    a native type value, you can create a '__getattribute__' method in the class that holds the attribute (see e.g., RTOptions).

    It is the hope to extend this to work with custom types, such as "percentage" (float with a max and min value),
    and direcotry (string that must end with '/').

    Args:
        object_name (str): name of type safe instance
        type (type): any of the native types in python that the instance must adhere to
        default (any): default value of the object, can be any including None
        value (any): value of type type assigned to instance, can be same or different than default

    Raises:
        OptionError: if wrong type is attempted.
    """

    def __init__(self, default, type, object_name):
        self.object_name = object_name
        self.type = type
        self.default = default
        self.value = self.default

    def __setattr__(self, name, value):
        """set attribute method that does the type checking, using native data types in python.
        The only 'exception' is allows float numbers to be written as a float or as an integer (but integers must always be integers).
        If a value of the wrong type is attempted set, the attribute value will be reset to the default value.
        Args:
            name (str): name of the attribute
            value (any): value to assign to the attribute
        """

        if name == "value":
            if type(value) in self.type:
                self.__dict__["value"] = value
            elif float in self.type and type(value) in [float, int]:
                self.__dict__["value"] = float(value)
            else:
                self.__dict__["value"] = self.default
                if value is not None:
                    raise OptionError(
                        f"Object {self.object_name} was assigned a value of type {type(value)} but is only allowed as type {self.type}."
                    )
        else:
            self.__dict__[name] = value


@dataclass
class RingtailDefaults:
    # maybe reconsider
    docking_mode: str = "adgpu"
    output_db: str = "output.db"
    storage_type: str = "duckdb"
    max_proc: int = None
    max_poses: int = 3
    store_all_poses: bool = False
    bookmark_name: str = "passing_results"
    enumerate_interaction_combs: bool = False
    log_file: str = None
    output_all_poses: bool = False
    input_db: str = None
    print_summary: bool = None
    verbose: bool = None
    debug: bool = None
    file: str = None
    file_path: str = None
    file_list: str = None
    file_pattern: str = None
    recursive: bool = None
    save_receptor: bool = None
    receptor_file: str = None
    append_results: bool = None
    duplicate_handling: str = None
    overwrite: bool = None
    interaction_tolerance: float = None
    add_interactions: bool = None
    interaction_cutoffs: str = "3.7,4.0"
    outfields: str = "LigName,docking_score"
    order_results: str = None
    mfpt_cluster: float = None
    interaction_cluster: float = None
    export_bookmark_csv: str = None
    export_query_csv: str = None
    export_bookmark_db: str = None
    export_receptor: bool = None
    export_sdf_path: str = None
    individual_sdf_files: bool = None
    data_from_bookmark: bool = None
    filter_bookmark: str = None
    find_similar_ligands: bool = None
    plot: bool = None
    pymol: bool = None


def default_dict(dataclass):
    return asdict(dataclass)


class ResultsObject:
    """Class that handles sources of data to be written including ligand data paths and how
    to traverse them, and options to store receptor.
    """

    def __init__(self):
        self.file = None
        self.file_path = None
        self.file_list = None
        self.recursive_path_traverse: bool = None
        self.receptor_file_path: str = None
        self.save_receptor: bool = None
        self.receptor_string: str = None
        self.strings: dict = None

    @property
    def target_name(self):
        if self.receptor_file_path and os.path.exists(self.receptor_file_path):
            return os.path.basename(self.receptor_file_path).split(".")[0]
        else:
            return None

    @property
    def has_results(self):
        results = bool(
            any(iterate_nested(self.file))
            or any(iterate_nested(self.file_path))
            or any(iterate_nested(self.file_list))
            or self.strings
        )
        return results

    @property
    def has_file_results(self):
        results = bool(
            any(iterate_nested(self.file))
            or any(iterate_nested(self.file_path))
            or any(iterate_nested(self.file_list))
        )
        return results


class Filters:
    """Object that holds all optional filters."""

    def __init__(self, filters: dict = {}):
        self.eworst: float = None
        self.ebest: float = None
        self.lebest: float = None
        self.leworst: float = None
        self.score_percentile: float = None
        self.le_percentile: float = None

        self.vdw_interactions: list = []
        self.hb_interactions: list = []
        self.reactive_interactions: list = []
        self.hb_count: int = None
        self.react_any: bool = None
        self.max_miss: int = 0

        self.ligand_name: str = None
        self.ligand_operator: str = None
        self.ligand_substruct: str = None
        self.ligand_substruct_pos: list = None
        self.ligand_max_atoms: int = None
        self.ligand_min_molweight: float = None
        self.ligand_max_molweight: float = None
        if filters:
            for key, value in filters.items():
                if hasattr(self, key):
                    if value is not None:
                        setattr(self, key, value)

        self.checks()

    def asdict(self):
        return copy.deepcopy(vars(self))

    def checks(self):
        """Ensures all values are internally consistent and valid. Runs once after all values are set initially,
        then every time a value is changed."""
        if self.eworst is not None and self.score_percentile is not None:
            logger.warning(
                "Cannot use 'eworst' cutoff with 'score_percentile'. Overiding 'score_percentile' with 'eworst'."
            )
            self.score_percentile = None

        if self.leworst is not None and self.le_percentile is not None:
            logger.warning(
                "Cannot use 'eworst' cutoff with 'le_percentile'. Overiding 'le_percentile' with 'leworst'."
            )
            self.le_percentile = None

        if self.score_percentile is not None and (
            self.score_percentile < 0 or self.score_percentile > 100
        ):
            raise OptionError(
                f"Given 'score_percentile' {self.score_percentile} not allowed. Should be within percentile range of 0-100."
            )

        if self.le_percentile is not None and (
            self.le_percentile < 0 or self.le_percentile > 100
        ):
            raise OptionError(
                f"Given 'score_percentile' {self.le_percentile} not allowed. Should be within percentile range of 0-100."
            )

        if self.ligand_operator not in ["OR", "AND"] and (
            self.ligand_substruct or self.ligand_substruct_pos
        ):
            logger.debug(f"'ligand_operator' set to default 'OR'.")
            self.ligand_operator = "OR"

        if self.max_miss < 0:
            raise OptionError("'max_miss' must be greater than or equal to 0.")

        if self.ligand_max_atoms and (
            self.ligand_min_molweight or self.ligand_max_molweight
        ):
            raise OptionError(
                "Cannot filter based on both max heavy atoms and mol weight restrictions."
            )

    @classmethod
    def get_filter_keys(self, group) -> list:
        """Provide keys associated with each of the filter groups.
        Args:
            group (str): includese property filters, interaction filters, ligand filters, or all filters
        Returns:
            list of filter keywords associated with the specified group(s)
        """

        if group.lower() not in ["property", "interaction", "ligand", "all"]:
            raise OptionError(
                f'{group.lower()} is not a valid filter group. Please use "property", "interactions", "ligand", or "all'
            )

        filter_groups = {
            "property": [
                "eworst",
                "ebest",
                "leworst",
                "lebest",
                "score_percentile",
                "le_percentile",
            ],
            "interaction": [
                "vdw_interactions",
                "hb_interactions",
                "reactive_interactions",
            ],
            "ligand": [
                "ligand_name",
                "ligand_substruct",
                "ligand_substruct_pos",
                "ligand_max_atoms",
                "ligand_operator",
                "ligand_min_molweight",
                "ligand_max_molweight",
            ],
        }
        if group.lower() == "all":
            list = []
            for i in filter_groups:
                list.extend(filter_groups[i])
            return list
        else:
            list = filter_groups[group.lower()]
        return list
