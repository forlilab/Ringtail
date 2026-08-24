#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail options handler
#

import os
from dataclasses import dataclass, asdict, fields

from .exceptions import OptionError, RTCoreError
from .logutils import get_logger

logger = get_logger(__name__)

docking_modes = {
    "adgpu": {"adgpu", "dlg", "gpu"},
    "vina": {"vina", "pdbqt"},
    "ad6": {"ad6", "adng", "ng", "zeta"},
}

docking_mode_file_ext = {"adgpu": "dlg", "vina": "pdbqt", "ad6": "sdf"}

docking_alias_to_mode = {
    alias: mode for mode, aliases in docking_modes.items() for alias in aliases
}


statuses = {1: "accepted", 2: "maybe", 3: "rejected", 0: None}


def validate_docking_mode(docking_mode: str):
    """Method that validates specified AutoDock program used to generate results.

    Args:
        docking_mode (str): string that describes docking mode

    Raises:
        RTCoreError: if docking_mode is not supported
    """
    if not isinstance(docking_mode, str):
        raise RTCoreError(
            f"The given docking mode was not given as a string: {type(docking_mode)}."
        )
    elif docking_mode.lower() not in docking_alias_to_mode:
        raise RTCoreError(
            f"Docking mode {docking_mode} is not supported. Please choose between {docking_modes}."
        )
    else:
        logger.debug(
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
    docking_mode: str = "ad6"
    output_db: str = "output.db"
    storage_type: str = "duckdb"
    max_proc: int = None
    max_poses: int = 3
    store_all_poses: bool = False
    bookmark_name: str = "passing_results"
    enumerate_interaction_combs: bool = False
    output_log: str = None
    output_all_poses: bool = False
    input_db: str = None
    print_summary: bool = None
    verbose: bool = None
    debug: bool = None
    docking_results: str = None
    recursive: bool = None
    save_receptor: bool = None
    receptor_file: str = None
    append_results: bool = None
    duplicate_handling: str = None
    overwrite: bool = None
    interaction_tolerance: float = None
    calculate_interactions: bool = True
    interaction_cutoffs: tuple[float] = (3.7, 4.0)  # HB CUTOFF,VDW CUTOFF
    outfields: tuple[str] = None
    order_results: str = None
    mfpt_cluster: float = None
    interaction_cluster: float = None
    export_bookmark_csv: str = None
    export_query_csv: str = None
    export_bookmark_db: str = None
    export_receptor_pdbqt: bool = None
    export_sdf_path: str = None
    individual_sdf_files: bool = None
    data_from_bookmark: bool = None
    input_bookmark: str = None
    find_similar_ligands: bool = None
    chunk_size: int = 5000

    @classmethod
    def all_fields(cls):
        return [f.name for f in fields(cls)]


def default_dict(dataclass):
    return asdict(dataclass)


def ringtail_defaults() -> dict:
    """Default values for the non-filter options.

    Filter options are deliberately absent: every one of them is either genuinely unset
    (None) or declares its default on the CLI argument itself, so there is nothing
    meaningful for an option-defaults dict to say about them.
    """
    return default_dict(RingtailDefaults())


class ResultsObject:
    """Class that handles sources of data to be written including ligand data paths and how
    to traverse them, and options to store receptor.
    """

    def __init__(self):
        self.docking_results = None
        self.recursive_path_traverse: bool = None
        self.receptor_file_path: str = None
        self.save_receptor: bool = None
        self.receptor_string: str = None

    @property
    def target_name(self):
        if self.receptor_file_path and os.path.exists(self.receptor_file_path):
            return os.path.basename(self.receptor_file_path).split(".")[0]
        else:
            return None
