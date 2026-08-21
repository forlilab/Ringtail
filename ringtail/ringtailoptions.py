#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail options handler
#

import os
from pathlib import Path
from typing import Callable
from .exceptions import OptionError, RTCoreError
from .logutils import get_logger
from .schema import (
    HB_COUNT_COLUMN,
    INTERACTION_TYPE_LETTERS,
    NUMERIC_FILTER_SCHEMA,
)

logger = get_logger(__name__)
from dataclasses import dataclass, asdict, field, fields

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
    defaults = {}
    defaults.update(default_dict(RingtailDefaults()))
    defaults.update(Filters().asdict())
    return defaults


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


@dataclass
class Filters:
    """Object that holds all optional filters."""

    eworst: float = None
    ebest: float = None
    lebest: float = None
    leworst: float = None
    score_percentile: float = None
    le_percentile: float = None
    vdw_interactions: list = field(default_factory=list)
    hb_interactions: list = field(default_factory=list)
    reactive_interactions: list = field(default_factory=list)
    hb_count: int = None
    react_any: bool = None
    max_miss: int = 0
    ligand_name: str = None
    ligand_name_file: str = None
    ligand_operator: str = None
    ligand_substruct: str = None
    ligand_substruct_pos: list = None
    ligand_max_atoms: int = None
    ligand_min_molweight: float = None
    ligand_max_molweight: float = None

    def __post_init__(self):
        if self.vdw_interactions is None:
            self.vdw_interactions = []
        if self.hb_interactions is None:
            self.hb_interactions = []
        if self.reactive_interactions is None:
            self.reactive_interactions = []
        self.checks()

    def asdict(self):
        return asdict(self)

    def checks(self):
        """Ensures all values are internally consistent and valid. Runs once after all values are
        set initially"""
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

        # normalize ligand names to one flat list so the count limit is meaningful
        if self.ligand_name:
            from .util import iterate_nested

            self.ligand_name = [
                name
                for item in iterate_nested(self.ligand_name)
                for name in item.split(",")
            ]
            if len(self.ligand_name) > 50:
                raise OptionError(
                    "The number of provided ligand names is too large, please prepare as a csv and use 'ligand_name_file' instead."
                )

        if (
            self.ligand_name_file
            and Path(self.ligand_name_file).suffix.lower() != ".csv"
        ):
            raise OptionError(
                f"The file of ligand names needs to be a csv file, cannot proceed with {Path(self.ligand_name_file).suffix.lower()}."
            )

    def to_criteria(self, percentile_cutoff: Callable = None) -> dict:
        """Group these filters into the structured criteria a query builder renders.

        Pure option-domain work: maps each filter onto the column it constrains (per
        ``schema.NUMERIC_FILTER_SCHEMA``), buckets criteria by kind, and normalizes values.
        Emits no SQL — the storage manager turns the ``(column, operator, value)`` tuples
        into dialect-specific predicates.

        Percentile filters are the exception: resolving a percentile to a cutoff requires
        reading the data, so the caller injects a resolver rather than this class holding a
        database handle.

        Args:
            percentile_cutoff (Callable[[float, str], float], optional): resolves
                (percentile, column) to a concrete cutoff. Required only when
                'score_percentile' or 'le_percentile' is set.

        Returns:
            dict: any of ``numeric`` (list of ``(column, operator, value)``), ``interactions``
            (list of ``[type_letter, chain, resname, resid, atom, wanted]``), ``max_miss``
            (present only alongside ``interactions``), and ``ligand`` (in-memory ligand
            filters, applied by RDKit rather than SQL). Empty buckets are omitted, so an
            empty dict means no filters were set.
        """
        numeric = []
        interactions = []
        ligand = {}
        max_miss = 0
        interaction_keys = self.get_filter_keys("interaction")
        ligand_keys = self.get_filter_keys("ligand")

        for key, value in self.asdict().items():
            if value is None:
                continue

            if key in NUMERIC_FILTER_SCHEMA:
                column, operator = NUMERIC_FILTER_SCHEMA[key]
                if key.endswith("_percentile"):
                    if percentile_cutoff is None:
                        raise OptionError(
                            f"Filtering on '{key}' needs a percentile resolver; "
                            "pass 'percentile_cutoff' to to_criteria()."
                        )
                    value = percentile_cutoff(value, column)
                numeric.append((column, operator, value))

            elif key == "hb_count":
                # positive means "at least this many"; negative means "no more than",
                # and 0 lands here so that "no hydrogen bonds" stays expressible
                if value > 0:
                    numeric.append((HB_COUNT_COLUMN, ">=", value))
                else:
                    numeric.append((HB_COUNT_COLUMN, "<=", -value))

            elif key in interaction_keys:
                for interact in value:
                    # interact is ["chain:res:resno:resatom", wanted(bool)]
                    interaction_string = INTERACTION_TYPE_LETTERS[key] + ":" + interact[0]
                    interactions.append(interaction_string.split(":") + [interact[1]])

            elif key == "react_any" and value:
                interactions.append(["R", "", "", "", "", True])

            elif key == "max_miss":
                max_miss = value

            elif key in ligand_keys or key == "ligand_name_file":
                if key == "ligand_substruct_pos":
                    for f in value:
                        f[1] = int(f[1])
                        for index in range(2, 6):
                            f[index] = float(f[index])
                ligand[key] = value

        criteria = {}
        if numeric:
            criteria["numeric"] = numeric
        if interactions:
            criteria["interactions"] = interactions
            criteria["max_miss"] = max_miss
        if ligand:
            criteria["ligand"] = ligand
        return criteria


    @classmethod
    def get_filter_keys(cls, group: str) -> list:
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
            lst = []
            for i in filter_groups:
                lst.extend(filter_groups[i])
            return lst
        else:
            lst = filter_groups[group.lower()]
        return lst
