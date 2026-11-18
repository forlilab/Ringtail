#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail options handler
#

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


def validate_docking_mode(docking_mode: str) -> str:
    """Method that validates specified AutoDock program used to generate results.

    Args:
        docking_mode (str): string that describes docking mode

    Returns:
        str: lower case of canonical docking mode

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


def check_file_pattern(docking_mode: str, file_pattern: str = None) -> str:
    """
    Checks if the file pattern to look for is the correct file pattern for the
    give ndocking mode

    Args:
        docking_mode (str): _description_
        file_pattern (str): _description_

    Returns:
        str: file pattern inferred from docking mode if not given, else given file pattern
    """
    if not file_pattern:
        return f"*.{docking_mode_file_ext[docking_mode]}*"
    if docking_mode_file_ext[docking_mode] not in file_pattern:
        raise OptionError(
            f"Requested file pattern {file_pattern} does not contain the necessary extension for docking mode {docking_mode}."
        )
    else:
        return file_pattern


def validate_bookmark_name(name: str) -> str | None:
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

    from .schema import ALL_TABLE_NAMES, CANDIDATES_NAME

    if name in ALL_TABLE_NAMES or name == CANDIDATES_NAME:
        logger.warning(
            f"Bookmark name '{name}' is reserved for a table and cannot be used."
        )
        return None
    return name if re.match(r"^[a-z0-9_]+$", name) else None


@dataclass
class RingtailDefaults:
    docking_mode: str = "ad6"
    output_db: str = "output.db"
    storage_type: str = "duckdb"
    max_proc: int | None = None
    max_poses: int = 3
    store_all_poses: bool | None = False
    bookmark_name: str = "passing_results"
    enumerate_interaction_combs: bool = False
    output_log: str | None = None
    output_all_poses: bool = False
    input_db: str | None = None
    print_summary: bool = None
    verbose: bool | None = None
    debug: bool | None = None
    docking_results: str | None = None
    recursive: bool | None = None
    save_receptor: bool | None = None
    receptor_file: str | None = None
    append_results: bool | None = None
    duplicate_handling: str | None = None
    overwrite: bool | None = None
    interaction_tolerance: float | None = None
    calculate_interactions: bool = True
    interaction_cutoffs: tuple[float, float] = (3.7, 4.0)  # HB cutoff,VDW cutoff
    outfields: tuple[str, ...] | None = None
    order_results: str | None = None
    mfpt_cluster: float | None = None
    interaction_cluster: float | None = None
    export_bookmark_csv: str | None = None
    export_query_csv: str | None = None
    export_bookmark_db: str | None = None
    export_receptor_pdbqt: bool | None = None
    export_sdf_path: str | None = None
    individual_sdf_files: bool | None = None
    data_from_bookmark: bool | None = None
    input_bookmark: str | None = None
    find_similar_ligands: bool | None = None
    chunk_size: int = 5000

    @classmethod
    def all_fields(cls):
        return [f.name for f in fields(cls)]


def ringtail_defaults() -> dict:
    """Default values for the non-filter options.

    Filter options are deliberately absent: every one of them is either genuinely unset
    (None) or declares its default on the CLI argument itself, so there is nothing
    meaningful for an option-defaults dict to say about them.
    """
    return asdict(RingtailDefaults())


def num_poses_to_store(max_poses: int, store_all_poses: bool) -> int:
    """
    Determine number of poses to store based on Ringtail user arguments

    Args:
        max_poses (int): number of poses given as integer
        store_all_poses (bool): yes/no to store all poses

    Returns:
        int: number of poses to store, -1 for all, positive int for anything else
    """
    if store_all_poses:
        return -1
    else:
        return max_poses
