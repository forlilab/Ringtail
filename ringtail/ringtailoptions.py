#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail options handler
#


import os
from .exceptions import OptionError, RTCoreError
from .util import iterate_nested
from dataclasses import dataclass, asdict, fields
import copy

from .logutils import LOGGER as logger

# Make class for options

# Make class for holding results data


docking_modes = {
    "adgpu": {"adgpu", "dlg", "gpu"},
    "vina": {"vina", "pdbqt"},
}

docking_mode_file_ext = {"adgpu": "dlg", "vina": "pdbqt"}

docking_alias_to_mode = {}
for canonical, aliases in docking_modes.items():
    for alias in aliases:
        docking_alias_to_mode[alias] = canonical


def validate_docking_mode(docking_mode: str):
    """Method that validates specified AutoDock program used to generate results.

    Args:
        docking_mode (str): string that describes docking mode

    Raises:
        RTCoreError: if docking_mode is not supported
    """
    if type(docking_mode) is not str:
        logger.warning(
            f'The given docking mode was not given as a string, it will be set to default value "{RingtailDefaults.docking_mode}".'
        )
        return RingtailDefaults.docking_mode
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
        mode = docking_alias_to_mode[docking_mode.lower()]
        return f"*.{docking_mode_file_ext[mode]}*"
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
    output_log: str = None
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
    interaction_cutoffs: tuple[float] = (3.7, 4.0)  # HB CUTOFF,VDW CUTOFF
    add_interactions: bool = True
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
    export_receptor_pdbqt: bool = None
    export_sdf_path: str = None
    individual_sdf_files: bool = None
    data_from_bookmark: bool = None
    filter_bookmark: str = None
    find_similar_ligands: bool = None
    plot: bool = None
    pymol: bool = None

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
        self.ligand_name_file: str = None
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


@dataclass
class RingtailOptions:
    """Class that handles options for the storage (database) manager class, including
    conflict handling, and results clustering and ordering."""

    duplicate_handling: str = None
    overwrite: bool = None
    bookmark_name: str = "passing_results"
    filter_bookmark: str = None
    order_results: str = None
    outfields: str = "Ligand_name, e"
    output_all_poses: bool = None
    mfpt_cluster: float = None
    interaction_cluster: float = None
    store_all_poses: bool = False
    max_poses: int = 3
    add_interactions: bool = True
    interaction_tolerance: float = None
    interaction_cutoffs: tuple[float, float] = (3.7, 4.0)
    max_proc: int = 2

    order_options = {
        "e",
        "le",
        "delta",
        "ref_rmsd",
        "e_inter",
        "e_vdw",
        "e_elec",
        "e_intra",
        "n_interact",
        "rank",
        "run",
        "hb",
    }


class OutputOptions:
    """Class that holds options related to reading and output from the database, including format for
    result export and alternate ways of displaying the data (plotting)."""

    options = {
        "log_file": {
            "default": None,
            "type": str,
            "description": "By default, read and filtering results are saved in 'output_log.txt'; if this option is used, ligands and requested info passing the filters will be written to specified file.",
        },
        "export_sdf_path": {
            "default": "",
            "type": str,
            "description": "Specify the path where to save poses of ligands passing the filters (SDF format); if the directory does not exist, it will be created; if it already exist, it will throw an error, unless the 'overwrite' is used  NOTE: the log file will be automatically saved in this path. Ligands will be stored as SDF files in the order specified.",
        },
        "individual_sdf_files": {
            "default": False,
            "type": bool,
            "description": "Use if you like to print chosen molecules to individual SDF files, as opposed to one big SDF.",
        },
        "enumerate_interaction_combs": {
            "default": None,
            "type": bool,
            "description": "When used with 'max_miss' > 0, will log ligands/poses passing each separate interaction filter combination as well as union of combinations. Can significantly increase runtime.",
        },
    }

    def __init__(self):
        super().initialize_from_dict(self.options, self.__class__.__name__)

    def checks(self):
        """Ensures all values are internally consistent and valid. Runs once after all values are set initially,
        then every time a value is changed."""
        if hasattr(self, "enumerate_interaction_combs"):
            if (
                self.export_sdf_path is not None
                and not self.export_sdf_path == ""
                and not self.export_sdf_path.endswith("/")
            ):
                self.export_sdf_path += "/"


class ReadOptions:
    """Object that holds choices and default values for read and export modes, mostly used for the command line interface."""

    options = {
        "plot": {
            "default": None,
            "type": bool,
            "description": "Makes scatterplot of LE vs Best Energy, saves as scatter.png.",
        },
        "pymol": {
            "default": None,
            "type": bool,
            "description": "Lauch PyMOL session and plot of ligand efficiency vs docking score for molecules in bookmark specified with --bookmark_name. Will display molecule in PyMOL when clicked on plot. Will also open receptor if given.",
        },
        "export_bookmark_db": {
            "default": None,
            "type": bool,
            "description": "Export a database containing only the results found in the bookmark specified by --bookmark_name. Will save as <input_db>_<bookmark_name>.db",
        },
        "export_receptor": {
            "default": None,
            "type": bool,
            "description": "Export stored receptor pdbqt. Will write to current directory.",
        },
        "data_from_bookmark": {
            "default": None,
            "type": bool,
            "description": "Write log of --outfields data for bookmark specified by --bookmark_name. Must use without any filters.",
        },
        "find_similar_ligands": {
            "default": None,
            "type": str,
            "description": "Allows user to find similar ligands to given ligand name based on previously performed morgan fingerprint or interaction clustering.",
        },
        "export_bookmark_csv": {
            "default": None,
            "type": str,
            "description": "Create csv of the bookmark given with bookmark_name. Output as <bookmark_name>.csv. Can also export full database tables.",
        },
        "export_query_csv": {
            "default": None,
            "type": str,
            "description": "Create csv of the requested SQL query. Output as query.csv. MUST BE PRE-FORMATTED IN SQL SYNTAX e.g. SELECT [columns] FROM [table] WHERE [conditions]",
        },
    }

    def __init__(self):
        super().initialize_from_dict(self.options, self.__class__.__name__)

    def checks(self):
        """Ensures all values are internally consistent and valid. Runs once after all values are set initially,
        then every time a value is changed."""
        pass


class GeneralOptions:
    """Object that holds choices and default values for miscellaneous arguments used for the command line interface only."""

    options = {
        "docking_mode": {
            "default": "dlg",
            "type": str,
            "description": "specify AutoDock program used to generate results. Available options are 'DLG' and 'vina'. Will automatically change --file_pattern to *.dlg* for DLG and *.pdbqt* for vina.",
        },
        "db_file": {
            "default": "output.db",
            "type": str,
            "description": "DB file for which to use for all Ringtail activities.",
        },
        "verbose": {
            "default": None,
            "type": bool,
            "description": "Print results passing filtering criteria to STDOUT and to log. NOTE: runtime may be slower option used.",
        },
        "debug": {
            "default": None,
            "type": bool,
            "description": "Print additional error information to STDOUT and to log.",
        },
        "print_summary": {
            "default": None,
            "type": bool,
            "description": "prints summary information about stored data to STDOUT.",
        },
    }

    def __init__(self):
        super().initialize_from_dict(self.options, self.__class__.__name__)

    def checks(self):
        """Ensures all values are internally consistent and valid. Runs once after all values are set initially,
        then every time a value is changed."""
        pass
