#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail options handler
#

import os
from .exceptions import OptionError
from .logutils import LOGGER
from .util import iterate_nested
from dataclasses import dataclass, asdict
import copy


docking_modes = {"adgpu": {"adgpu", "dlg", "gpu"}, "vina": {"vina", "pdbqt"}}

docking_mode_file_ext = {"adgpu": "dlg", "vina": "pdbqt"}

docking_alias_to_mode = {}
for canonical, aliases in docking_modes.items():
    for alias in aliases:
        docking_alias_to_mode[alias] = canonical

statuses = ["accepted", "maybe", "rejected"]


def validate_docking_mode(docking_mode: str):
    """Method that validates specified AutoDock program used to generate results.

    Args:
        docking_mode (str): string that describes docking mode

    Raises:
        RTCoreError: if docking_mode is not supported
    """
    if type(docking_mode) is not str:
        LOGGER.warning(
            f'The given docking mode was not given as a string, it will be set to default value "{RingtailDefaults.docking_mode}".'
        )
        return RingtailDefaults.docking_mode
    elif docking_mode.lower() not in docking_alias_to_mode:
        raise NotImplementedError(
            f"Docking mode {docking_mode} is not supported. Please choose between {docking_modes}."
        )
    else:
        LOGGER.debug(
            f"Docking mode {docking_mode} is valid, will be used for results parsing."
        )
        # return canonical docking mode
        return docking_alias_to_mode[docking_mode.lower()]


def validate_file_pattern(docking_mode: str, file_pattern: str = None) -> str:
    """
    Assumes docking_mode has been validated

    Args:
        docking_mode (str): _description_
        file_pattern (str): _description_

    Returns:
        str: _description_
    """
    if not file_pattern:
        return f"*.{docking_mode_file_ext[docking_mode]}*"
    if docking_mode_file_ext[docking_mode] not in file_pattern:
        raise OptionError(
            f"Requested file pattern {file_pattern} does not contain the necessary extension for docking mode {docking_mode}."
        )
    else:
        return file_pattern


@dataclass
class RingtailDefaults:
    # maybe reconsider
    docking_mode: str = "adgpu"
    output_db: str = "output.db"
    storage_type: str = "sqlite"
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
    no_interactions: bool = None
    interaction_cutoffs: tuple[float] = (
        3.7,
        4.0,
    )  # HB CUTOFF,VDW CUTOFF
    outfields: tuple[str] = (
        "LigName",
        "docking_score",
    )
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
        self.file_pattern: str = None
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
            LOGGER.warning(
                "Cannot use 'eworst' cutoff with 'score_percentile'. Overiding 'score_percentile' with 'eworst'."
            )
            self.score_percentile = None

        if self.leworst is not None and self.le_percentile is not None:
            LOGGER.warning(
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
            LOGGER.debug(f"'ligand_operator' set to default 'OR'.")
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
