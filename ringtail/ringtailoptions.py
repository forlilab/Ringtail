#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail options handler
#

import os
from .exceptions import OptionError
from .logutils import LOGGER as logger
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
            if type(value) == self.type:
                self.__dict__["value"] = value
            elif self.type == float and type(value) in [float, int]:
                self.__dict__["value"] = float(value)
            else:
                self.__dict__["value"] = self.default
                if value is not None:
                    raise OptionError(
                        f"Object {self.object_name} was assigned a value of type {type(value)} but is only allowed as type {self.type}."
                    )
        else:
            self.__dict__[name] = value


class RTOptions:
    """Holds standard methods for the ringtail option child classes.
    Options can be added using this format:
    options = {
        "":{
            "default":'',
            "type":'',
            "description": ""
        },
    }

    """

    @classmethod
    def is_valid_path(self, path):
        """Checks if path exist in current directory.

        Args:
            path (str)

        Returns:
            bool: if path exist

        """
        if os.path.exists(path):
            return True
        else:
            return False

    @staticmethod
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

    def initialize_from_dict(self, dict: dict, name):
        """Initializes a child objects using the values available in their option dictionary.

        Args:
            dict (dict): of attributes to be initialized to the object
            name (str): name of the childclass/object
        """
        for item, info in dict.items():
            setattr(
                self,
                item,
                TypeSafe(default=info["default"], type=info["type"], object_name=item),
            )
        logger.info(
            f"A {name} object was created with default values for all attributes."
        )

    def todict(self):
        """Return class and its attributes as a dict of native types and not as objects (which they are if they are type checked using TypeSafe)."""
        dict = {}
        for item in self.options.keys():
            dict[item] = getattr(self, item)
        returndict = copy.deepcopy(dict)
        return returndict

    def __getattribute__(self, attribute):
        """Ensures attribute is returned as valeu and not an object.
        Args:
            attribute (str): name of attribute
        """
        dataobject = object.__getattribute__(self, attribute)
        if isinstance(dataobject, TypeSafe):
            value = dataobject.value
        else:
            value = dataobject
        return value

    def __setattr__(self, attribute, value):
        """Overloaded set attribute method that ensures attribute is type checked during setting.

        Args:
            attribute (str): name of attribute
            value (any): value given to attribute
        """
        if not hasattr(self, attribute):
            object.__setattr__(self, attribute, value)
        else:
            dataobject = object.__getattribute__(self, attribute)
            dataobject.__setattr__("value", value)
        # runs internal (consistency) checks specified in each child class
        self.checks()


class InputStrings(RTOptions):
    """Class that handles docking results strings from vina docking, with options to store receptor.
    Takes docking results string as a dictionary of:
    {ligand_name: docking_result}"""

    options = {
        "results_strings": {
            "default": None,
            "type": dict,
            "description": "A dictionary of ligand names and ligand docking output results. Currently only valid for vina docking",
        },
        "receptor_file": {
            "default": None,
            "type": str,
            "description": "Use with Vina mode. Give file for receptor PDBQT.",
        },
        "save_receptor": {
            "default": None,
            "type": bool,
            "description": "Saves receptor PDBQT to database. Receptor location must be specied with in 'receptor_file'.",
        },
        "target": {
            "default": None,
            "type": str,
            "description": "Name of receptor. This field is autopopulated if 'receptor_file' is supplied.",
        },
    }

    def __init__(self):
        super().initialize_from_dict(self.options, self.__class__.__name__)

    def checks(self):
        """Ensures all values are internally consistent and valid. Runs once after all values are set initially,
        then every time a value is changed."""
        if hasattr(
            self, "target"
        ):  # ensures last item in the option dictionary has been
            if type(self.results_strings) == str:
                pass
                self.results_strings = list(self.results_strings)
            if type(self.target) != str:
                if self.receptor_file is None:
                    pass
                elif self.receptor_file is not None and self.is_valid_path(
                    self.receptor_file
                ):
                    self.target = os.path.basename(self.receptor_file).split(".")[0]
                else:
                    raise OptionError(
                        "The receptor PDBQT file path is not valid. Please check location of receptor file and 'receptor_file' option."
                    )


class InputFiles(RTOptions):
    """Class that handles sources of data to be written including ligand data paths and how
    to traverse them, and options to store receptor.
    """

    options = {
        "file": {
            "default": None,
            "type": list,
            "description": "Ligand docking output file to save. Compressed (.gz) files allowed. Only results files associated the same receptor allowed.",
        },
        "file_path": {
            "default": None,
            "type": list,
            "description": "Directory(s) containing docking output files to save. Compressed (.gz) files allowed",
        },
        "file_list": {
            "default": None,
            "type": list,
            "description": "Text file(s) containing the list of docking output files to save; relative or absolute paths are allowed. Compressed (.gz) files allowed.",
        },
        "file_pattern": {
            "default": None,
            "type": str,
            "description": "Specify which pattern to use when searching for result files to process (only with 'file_path').",
        },
        "recursive": {
            "default": None,
            "type": bool,
            "description": "Enable recursive directory scan when 'file_path' is used.",
        },
        "receptor_file": {
            "default": None,
            "type": str,
            "description": "Use with Vina mode. Give file for receptor PDBQT.",
        },
        "save_receptor": {
            "default": None,
            "type": bool,
            "description": "Saves receptor PDBQT to database. Receptor location must be specied with in 'receptor_file'.",
        },
        "target": {
            "default": None,
            "type": str,
            "description": "Name of receptor. This field is autopopulated if 'receptor_file' is supplied.",
        },
    }

    def __init__(self):
        super().initialize_from_dict(self.options, self.__class__.__name__)

    def checks(self):
        """Ensures all values are internally consistent and valid. Runs once after all values are set initially,
        then every time a value is changed."""
        if hasattr(
            self, "target"
        ):  # ensures last item in the option dictionary has been
            if type(self.target) != str:
                if self.receptor_file is None:
                    pass
                elif self.receptor_file is not None and self.is_valid_path(
                    self.receptor_file
                ):
                    self.target = os.path.basename(self.receptor_file).split(".")[0]
                else:
                    raise OptionError(
                        "The receptor PDBQT file path is not valid. Please check location of receptor file and 'receptor_file' option."
                    )


class ResultsProcessingOptions(RTOptions):
    """Class that holds database write options that affects write time, such as how to
    break up data files, number of computer processes to use, and and how many poses to store.
    """

    options = {
        "store_all_poses": {
            "default": False,
            "type": bool,
            "description": "Store all poses from input files. Overrides 'max_poses'.",
        },
        "max_poses": {
            "default": 3,
            "type": int,
            "description": "Store top pose for top n clusters.",
        },
        "add_interactions": {
            "default": False,
            "type": bool,
            "description": "Find interactions between ligand poses and receptor and save to database. Requires receptor PDBQT to be given with input files (all modes) and 'receptor_file' to be specified with Vina mode. SIGNIFICANTLY INCREASES DATBASE WRITE TIME.",
        },
        "interaction_tolerance": {
            "default": None,
            "type": float,
            "description": "Will add the interactions for poses within some tolerance RMSD range of the top pose in a cluster to that top pose. Can use as flag with default tolerance of 0.8 for cmd line tool, or give other value as desired (cmd line and api). Only compatible with ADGPU mode.",
        },
        "interaction_cutoffs": {
            "default": [3.7, 4.0],
            "type": list,
            "description": "Use with 'add_interactions', specify distance cutoffs for measuring interactions between ligand and receptor in angstroms. Give as string, separating cutoffs for hydrogen bonds and VDW with comma (in that order). E.g. '-ic 3.7,4.0' will set the cutoff for hydrogen bonds to 3.7 angstroms and for VDW to 4.0. These are the default cutoffs.",
        },
        "max_proc": {
            "default": None,
            "type": int,
            "description": "Maximum number of processes to create during parallel file parsing. Defaults to number of CPU processors.",
        },
    }

    def __init__(self):
        super().initialize_from_dict(self.options, self.__class__.__name__)

    def checks(self):
        """Ensures all values are internally consistent and valid. Runs once after all values are set initially,
        then every time a value is changed."""
        if hasattr(self, "max_proc"):
            if self.store_all_poses is True and self.interaction_tolerance is not None:
                logger.warning(
                    "Cannot use 'interaction_tolerance' with 'store_all_poses'. Removing 'interaction_tolerance'."
                )
                self.interaction_tolerance = None


class StorageOptions(RTOptions):
    """Class that handles options for the storage (database) manager class, including
    conflict handling, and results clustering and ordering."""

    options = {
        "duplicate_handling": {
            "default": None,
            "type": str,
            "description": "Specify how duplicate Results rows should be handled when inserting into database. Options are 'ignore' or 'replace'. Default behavior (no option provided) will allow duplicate entries.",
        },
        "filter_bookmark": {
            "default": None,
            "type": str,
            "description": "Perform filtering over specified bookmark.",
        },
        "overwrite": {
            "default": None,
            "type": bool,
            "description": "This option will allow overwriting of the database (in 'write'/add files-mode) and filtering log_file (in 'read'/filtering mode).",
        },
        "order_results": {
            "default": None,
            "type": str,
            "description": """Stipulates how to order the results when written to the log file. By default will be ordered by order results were added to the database. ONLY TAKES ONE OPTION.
                            available fields are:  
                            'e' (docking_score), 
                            'le' (ligand efficiency), 
                            'delta' (delta energy from best pose), 
                            'ref_rmsd' (RMSD to reference pose), 
                            'e_inter' (intermolecular energy), 
                            'e_vdw' (van der waals energy), 
                            'e_elec' (electrostatic energy), 
                            'e_intra' (intermolecular energy), 
                            'n_interact' (number of interactions), 
                            'rank' (rank of ligand pose), 
                            'run' (run number for ligand pose), 
                            'hb' (hydrogen bonds); """,
        },
        "outfields": {
            "default": "Ligand_name,e",
            "type": str,
            "description": """Defines which fields are used when reporting the results (to stdout and to the log file). Fields are specified as comma-separated values, e.g. 'outfields=e,le,hb'; by default, docking_score (energy) and ligand name are reported. Ligand always reported in first column available fields are: \n
                            'Ligand_name' (Ligand name), 
                            'e' (docking_score), 
                            'le' (ligand efficiency), 
                            'delta' (delta energy from best pose), 
                            'ref_rmsd' (RMSD to reference pose), 
                            'e_inter' (intermolecular energy), 
                            'e_vdw' (van der waals energy), 
                            'e_elec' (electrostatic energy), 
                            'e_intra' (intermolecular energy), 
                            'n_interact' (number of iteractions), 
                            'ligand_smile' , 
                            'rank' (rank of ligand pose), 
                            'run' (run number for ligand pose), 
                            'hb' (hydrogen bonds), 
                            'receptor' (receptor name);""",
        },
        "output_all_poses": {
            "default": None,
            "type": bool,
            "description": "By default, will output only top-scoring pose passing filters per ligand. This flag will cause each pose passing the filters to be logged.",
        },
        "mfpt_cluster": {
            "default": None,
            "type": float,
            "description": "Cluster filtered ligands by Tanimoto distance of Morgan fingerprints with Butina clustering and output ligand with lowest ligand efficiency from each cluster. Useful for selecting chemically dissimilar ligands.",
        },
        "interaction_cluster": {
            "default": None,
            "type": float,
            "description": "Cluster filtered ligands by Tanimoto distance of interaction fingerprints with Butina clustering and output ligand with lowest ligand efficiency from each cluster. Useful for enhancing selection of ligands with diverse interactions.",
        },
        "bookmark_name": {
            "default": "passing_results",
            "type": str,
            "description": "name for resulting book mark file. Default value is 'passing_results'",
        },
    }

    def __init__(self):
        super().initialize_from_dict(self.options, self.__class__.__name__)

    def checks(self):
        """Ensures all values are internally consistent and valid. Runs once after all values are set initially,
        then every time a value is changed."""
        if hasattr(self, "bookmark_name"):
            # Makes sure all default values have been set once for internal comparisons
            if self.bookmark_name is not None and not RTOptions.valid_bookmark_name(
                self.bookmark_name
            ):
                raise OptionError(
                    "The chosen bookmark name {0} is not valid, as it contains symbols other than letters, numbers, and underscore (_)".format(
                        self.bookmark_name
                    )
                )
            if self.duplicate_handling is not None:
                if self.duplicate_handling.upper() not in ["IGNORE", "REPLACE"]:
                    logger.warning(
                        f"--duplicate_handling option {self.duplicate_handling} not allowed. Reverting to default behavior."
                    )
                    self.duplicate_handling = None
            if (
                self.order_results is not None
                and self.order_results not in self.order_options
            ):
                raise OptionError(
                    "Requested ording option that is not available. Please see --help for available options."
                )
            # Make sure we include ligand name in output columns
            if self.outfields is not None and "Ligand_name" not in self.outfields:
                self.outfields = "Ligand_name," + self.outfields

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


class OutputOptions(RTOptions):
    """Class that holds options related to reading and output from the database, including format for
    result export and alternate ways of displaying the data (plotting)."""

    options = {
        "log_file": {
            "default": "output_log.txt",
            "type": str,
            "description": "By default, read and filtering results are saved in 'output_log.txt'; if this option is used, ligands and requested info passing the filters will be written to specified file.",
        },
        "export_sdf_path": {
            "default": "",
            "type": str,
            "description": "Specify the path where to save poses of ligands passing the filters (SDF format); if the directory does not exist, it will be created; if it already exist, it will throw an error, unless the 'overwrite' is used  NOTE: the log file will be automatically saved in this path. Ligands will be stored as SDF files in the order specified.",
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


class Filters(RTOptions):
    """Object that holds all optional filters."""

    options = {
        "eworst": {
            "default": None,
            "type": float,
            "description": "Specify the worst energy value accepted.",
        },
        "ebest": {
            "default": None,
            "type": float,
            "description": "Specify the best energy value accepted.",
        },
        "leworst": {
            "default": None,
            "type": float,
            "description": "Specify the worst ligand efficiency value accepted.",
        },
        "lebest": {
            "default": None,
            "type": float,
            "description": "Specify the best ligand efficiency value accepted.",
        },
        "score_percentile": {
            "default": None,
            "type": float,
            "description": "Specify the worst energy percentile accepted. Express as percentage e.g. 1 for top 1 percent.",
        },
        "le_percentile": {
            "default": None,
            "type": float,
            "description": "Specify the worst ligand efficiency percentile accepted. Express as percentage e.g. 1 for top 1 percent.",
        },
        "vdw_interactions": {
            "default": [],
            "type": list,
            "description": "Define van der Waals interactions with residue as [-][CHAIN]:[RES]:[NUM]:[ATOM_NAME]. E.g., [('A:VAL:279:', True), ('A:LYS:162:', True)] -> [('chain:resname:resid:atomname', <wanted (bool)>), ('chain:resname:resid:atomname', <wanted (bool)>)].",
        },
        "hb_interactions": {
            "default": [],
            "type": list,
            "description": "Define HB (ligand acceptor or donor) interaction as [-][CHAIN]:[RES]:[NUM]:[ATOM_NAME]. E.g., [('A:VAL:279:', True), ('A:LYS:162:', True)] -> [('chain:resname:resid:atomname', <wanted (bool)>), ('chain:resname:resid:atomname', <wanted (bool)>)].",
        },
        "reactive_interactions": {
            "default": [],
            "type": list,
            "description": "Check if ligand reacted with specified residue as [-][CHAIN]:[RES]:[NUM]:[ATOM_NAME]. E.g., [('A:VAL:279:', True), ('A:LYS:162:', True)] -> [('chain:resname:resid:atomname', <wanted (bool)>), ('chain:resname:resid:atomname', <wanted (bool)>)].",
        },
        "hb_count": {
            "default": None,
            "type": list,
            "description": "Accept ligands with at least the requested number of HB interactions. If a negative number is provided, then accept ligands with no more than the requested number of interactions. E.g., [('hb_count', 5)].",
        },
        "react_any": {
            "default": None,
            "type": bool,
            "description": "Check if ligand reacted with any residue.",
        },
        "max_miss": {
            "default": 0,
            "type": int,
            "description": "Will compute all possible combinations of interaction filters excluding up to 'max_miss' number of interactions from given set. Default will only return union of poses interaction filter combinations. Use with 'enumerate_interaction_combs' for enumeration of poses passing each individual combination of interaction filters.",
        },
        "ligand_name": {
            "default": [],
            "type": list,
            "description": "Specify ligand name(s). Will combine name filters with 'OR'.",
        },
        "ligand_substruct": {
            "default": [],
            "type": list,
            "description": "SMARTS pattern(s) for substructure matching.",
        },
        "ligand_substruct_pos": {
            "default": [],
            "type": list,
            "description": "SMARTS pattern(s) for substructure matching, e.g., [''[Oh]C' 0 1.2 -5.5 10.0 15.5'] -> ['smart_string index_of_positioned_atom cutoff_distance x y z'].",
        },
        "ligand_max_atoms": {
            "default": None,
            "type": int,
            "description": "Maximum number of heavy atoms a ligand may have.",
        },
        "ligand_operator": {
            "default": "OR",
            "type": str,
            "description": "Logical join operator for multiple SMARTS.",
        },
    }

    def __init__(self):
        super().initialize_from_dict(self.options, self.__class__.__name__)

    def checks(self):
        """Ensures all values are internally consistent and valid. Runs once after all values are set initially,
        then every time a value is changed."""
        if hasattr(self, "ligand_operator"):
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

            if self.ligand_operator not in ["OR", "AND"]:
                raise OptionError(
                    f"Given 'ligand_operator' {self.ligand_operator} not allowed. Must be 'OR' or 'AND'."
                )

            if self.max_miss < 0:
                raise OptionError("'max_miss' must be greater than or equal to 0.")

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


class ReadOptions(RTOptions):
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


class GeneralOptions(RTOptions):
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
