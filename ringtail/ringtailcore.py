#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail virtual screening manager
#
import itertools
import os
from typing import Callable, Union
from collections import defaultdict
import json
from pathlib import Path
import matplotlib.pyplot as plt
from meeko import RDKitMolCreate
from rdkit import Chem

from .util import *
from .logutils import LOGGER
from .ringtailoptions import (
    RingtailDefaults,
    Filters,
    ResultsObject,
    default_dict,
    statuses,
    validate_docking_mode,
)
from .exceptions import (
    RTCoreError,
    OutputError,
    StorageError,
    ResultsProcessingError,
    OptionError,
)
from .resultsmanager import ResultsManager
from .receptormanager import ReceptorManager
from .outputmanager import OutputManager
from .storagemanager import StorageManager
from .storagemanager_duckdb import StorageManagerDuckDB, HAS_DUCK
from .storagemanager_sqlite import StorageManagerSQLite, HAS_SQLITE

storage_types = {}
if HAS_SQLITE:
    storage_types.update({"sqlite": StorageManagerSQLite})
if HAS_DUCK:
    storage_types.update({"duckdb": StorageManagerDuckDB})


def get_valid_storageclass(storage_type) -> StorageManager:
    """Checks if chosen storage type has been implemented

    Args:
        storage_type (str): name of the storage type

    Raises:
        NotImplementedError: raised if seelected storage type has not been implemented

    Returns:
        class: of implemented storage type
    """

    if storage_type in storage_types:
        return storage_types[storage_type]
    else:
        raise NotImplementedError(
            f"Given storage type {storage_type} is not implemented."
        )


class RingtailCore:
    """Core class for coordinating different actions on virtual screening
    including adding results to storage, filtering and clusteirng, and outputting data as
    rdkit molecules, plotting docking results, and visualizing select ligands in pymol.

    Attributes:
        db_file (str): name of database file being operated on
        docking_mode (str): specifies what docking mode has been used for the results in the database
        storageman (StorageManager): Interface module with database
        resultsman (ResultsManager): Module to deal with results processing before adding to database
        storagetype (str): database type
        LOGGER (RaccoonLogger): handles logging
        _run_mode (str): refers to whether ringtail is ran from the command line or through direct API use, where the former is more restrictive
    """

    def __init__(
        self,
        db_file: str = RingtailDefaults.output_db,
        storage_type: str = RingtailDefaults.storage_type,
        logging_level: str = "WARNING",
        access_mode: str = "api",
    ):
        """Initialize ringtail core, and create a storageman object with the db file.
        Can set logger level here, otherwise change it by LOGGER.setLevel("level")

        Args:
            db_file (str): Database file to initialize core with. Defaults to "output.db".
            storage_type (str, optional): Database setup to use for interacting with the database.
                        If db file exist, it will attempt to detect storage type. Defaults to "sqlite".
            logging_level (str, optional): Global logger level. Defaults to "DEBUG".
            access_mode (str, optional): through what process Ringtail is accessed. Defaults to "api".
        """

        # Initiate logging
        LOGGER.set_level(logging_level)
        # create log file handler if log level is debug
        if LOGGER.level() == "DEBUG":
            LOGGER.add_filehandler()

        # try to find existing storage type, will overwrite input
        if Path(db_file).is_file():
            storage_type = detect_db_type(db_file)

        # Check if storage type is implemented
        try:
            storageman = get_valid_storageclass(storage_type)
        except NotImplementedError as e:
            LOGGER.error(e)
            raise e

        # initialize remaining variables and storagemanager
        self.storagetype = storage_type
        self.db_file = db_file
        self.storageman: StorageManager = storageman(db_file)

        LOGGER.info(
            f"[     New RingtailCore object initialized with database file {db_file}    ]"
        )

        self._run_mode = access_mode

    def get_previous_docking_mode(self) -> Union[None, str]:
        """
        Will return the last used docking mode/engine of the database, normalized to the standardized
        names for docking modes used in ringtail (e.g., see RingtailOptions)

        Returns:
            Union[None, str]: None or previous docking mode
        """
        with self.storageman as sm:
            return validate_docking_mode(sm.get_previous_docking_mode())

    def add_pose(
        self,
        pose_data: Union[dict, list[dict]],
        duplicate_handling: str = RingtailDefaults.duplicate_handling,
        docking_mode: str = RingtailDefaults.docking_mode,
        store_all_poses: bool = RingtailDefaults.store_all_poses,
        max_poses: int = RingtailDefaults.max_poses,
        no_interactions: bool = RingtailDefaults.no_interactions,
        interaction_tolerance: float = RingtailDefaults.interaction_tolerance,
        interaction_cutoffs: list = RingtailDefaults.interaction_cutoffs,
    ):
        interaction_vars = {
            "no_interactions": no_interactions,
            "interaction_tolerance": interaction_tolerance,
            "interaction_cutoffs": interaction_cutoffs,
        }
        write_options = {
            "duplicate_handling": duplicate_handling,
            "docking_mode": docking_mode,
            "num_poses": self._num_poses_to_store(max_poses, store_all_poses),
        }
        pass

    # region write to database
    def add_results_from_files(
        self,
        file: str = RingtailDefaults.file,
        file_path: str = RingtailDefaults.file_path,
        file_list: str = RingtailDefaults.file_list,
        file_pattern: str = None,
        recursive: bool = RingtailDefaults.recursive,
        receptor_file: str = None,
        save_receptor: bool = RingtailDefaults.save_receptor,
        duplicate_handling: str = RingtailDefaults.duplicate_handling,
        overwrite: bool = RingtailDefaults.overwrite,
        docking_mode: str = RingtailDefaults.docking_mode,
        store_all_poses: bool = RingtailDefaults.store_all_poses,
        max_poses: int = RingtailDefaults.max_poses,
        no_interactions: bool = RingtailDefaults.no_interactions,
        interaction_tolerance: float = RingtailDefaults.interaction_tolerance,
        interaction_cutoffs: list = RingtailDefaults.interaction_cutoffs,
        max_proc: int = RingtailDefaults.max_proc,
        finalize: bool = True,
    ):
        """
        Call storage manager to process result files and add to database. Creates or adds to an existing a database.
        Options can be provided as a dict or as individual options. If both are provided, individual options will overwrite those from the dictionary.

        Args:
            file (str, optional: list(str)): ligand result file
            file_path (str, optional: list(str)): list of folders containing one or more result files
            file_list (str, optional: list(str)): list of ligand result file(s)
            file_pattern (str, optional): optional custom file pattern, such as MY_FAV_PROT_*.pdbqt*
            recursive (bool): used to recursively search file_path for folders inside folders
            receptor_file (str): string containing the receptor .pdbqt
            save_receptor (bool): whether or not to store the full receptor details in the database (needed for some things)
            duplicate_handling (str, options): specify how duplicate Results rows should be handled when inserting into database. Options are "ignore" or "replace". Default behavior will allow duplicate entries.
            overwrite (bool): whether or not to overwrite existing storage
            docking_mode (str): what docking engine was used to perform the docking
            store_all_poses (bool): store all ligand poses, does it take precedence over max poses?
            max_poses (int): how many poses to save (ordered by soem score?)
            no_interactions (bool): use if not wanting Ringtail to calculate and store interactions, only in vina mode
            interaction_tolerance (float): longest ångström distance that is considered interaction?
            interaction_cutoffs (list): ångström distance cutoffs for x and y interaction
            max_proc (int): max number of computer processors to use for file reading
            finalize (bool, optional): whether to finalize database creation and write indices, or wait, for example when you write to the db in batches

        Raises:
            OptionError

        """
        results_data = ResultsObject()
        results_data.file = file
        results_data.file_path = file_path
        results_data.file_list = file_list
        results_data.file_pattern = file_pattern
        results_data.recursive_path_traverse = recursive
        results_data.receptor_file_path = receptor_file
        results_data.save_receptor = save_receptor

        if not results_data.has_results and not results_data.save_receptor:
            raise OptionError(
                "At least one input option needs to be used: file, file_path, file_list, or save_receptor"
            )

        # If there are ligand files present, process ligand data
        if results_data.has_results or results_data.save_receptor:
            self._add_results(
                results_data,
                duplicate_handling,
                overwrite,
                docking_mode,
                self._num_poses_to_store(max_poses, store_all_poses),
                no_interactions,
                interaction_tolerance,
                interaction_cutoffs,
                max_proc,
                finalize,
            )

    def add_results_from_vina_string(
        self,
        results: dict = None,
        receptor_file: str = RingtailDefaults.save_receptor,
        save_receptor: bool = RingtailDefaults.save_receptor,
        duplicate_handling: str = RingtailDefaults.duplicate_handling,
        overwrite: bool = RingtailDefaults.overwrite,
        docking_mode: str = "vina",
        store_all_poses: bool = RingtailDefaults.store_all_poses,
        max_poses: int = RingtailDefaults.max_poses,
        no_interactions: bool = RingtailDefaults.no_interactions,
        interaction_cutoffs: list = RingtailDefaults.interaction_cutoffs,
        max_proc: int = RingtailDefaults.max_proc,
        finalize: bool = True,
    ):
        """
        Call storage manager to process the given vina output string and add to database.
        Options can be provided as a dict or as individual options.
        Creates or adds to an existing a database.

        Args:
            results (dict): string containing the ligand identified and docking results as a dictionary
            receptor_file (str): string containing the receptor .pdbqt
            save_receptor (bool): whether or not to store the full receptor details in the database (needed for some things)
            duplicate_handling (str, options): specify how duplicate Results rows should be handled when inserting into database. Options are "ignore" or "replace". Default behavior will allow duplicate entries.
            overwrite (bool): whether or not to overwrite existing storage
            docking_mode (str): what docking engine was used to perform the docking
            store_all_poses (bool): store all ligand poses, does it take precedence over max poses?
            max_poses (int): how many poses to save (ordered by soem score?)
            no_interactions (bool): use if not wanting Ringtail to calculate and store interactions, only in vina mode
            interaction_cutoffs (list): ångström distance cutoffs for x and y interaction
            max_proc (int): max number of computer processors to use for file reading
            finalize (bool, optional): whether to finalize database creation and write indices, or wait, for example when you write to the db in batches

        Raises:
            OptionError

        """
        # Method currently only works with vina output, set automatically
        if validate_docking_mode(docking_mode) != "vina":
            docking_mode = "vina"

        # create results string object
        results_data = ResultsObject()
        results_data.strings = results
        results_data.receptor_file_path = receptor_file
        results_data.save_receptor = save_receptor

        if not results_data.strings and not results_data.save_receptor:
            raise OptionError(
                "At least one input option needs to be used: 'results', or 'save_receptor'"
            )

        # If there are any docking data strings, process docking results
        if results_data.strings or results_data.save_receptor:
            self._add_results(
                results_data,
                duplicate_handling,
                overwrite,
                docking_mode,
                self._num_poses_to_store(max_poses, store_all_poses),
                no_interactions,
                None,
                interaction_cutoffs,
                max_proc,
                finalize,
            )

    def finalize_write(self):
        """
        Finalize database write by creating indices
        """
        with self.storageman:
            self.storageman.finalize_database_write()

    def save_receptor(self, receptor: Union[str, dict]):
        """
        Add receptor to database. Accets .pdbqt file as well as a polymer .json or a polymer object

        Args:
            receptor_file (str): path to receptor file

        """

        receptor_blob = None
        receptor_jsons = None
        if type(receptor) == str:
            # get all extensions, lower case, account for zipped files
            extensions = Path(receptor).suffixes
            extensions = [e.lower() for e in extensions]
            if ".json" in extensions:
                receptor_name, receptor_jsons = self._process_receptor_polymer(receptor)
            elif ".pdbqt" in extensions:
                receptor_name, receptor_blob = self._process_receptor_pdbqt(receptor)
        else:
            from meeko import Polymer

            try:
                receptor_name = receptor_name.keys()[0]
                receptor_polymer: Polymer = receptor[receptor_name]
                receptor_jsons = receptor_polymer.to_json()
            except:
                raise OptionError(
                    f"Unable to serialize provided receptor object, cannot proceed to save receptor."
                )

        # check if db has receptor added
        db_name, db_blob, db_jsons = self.get_receptor_object()
        if db_name:
            # make all checks if it is possible to add the requested receptor data
            if receptor_name != db_name:
                raise OptionError(
                    "You attempted to add a receptor to a database that already has a receptor in it, and they do not appear to be the same receptor."
                )
            # if trying to add blob and there is blob, return
            if receptor_blob and db_blob:
                LOGGER.warning(
                    "The Receptors table already has a .pdbqt blob present, will not add another."
                )
                return
            # if trying to add polymer json and there is polymer json, return
            if receptor_jsons and db_jsons:
                LOGGER.warning(
                    "The Receptors table already has a polymer json string present, will not add another."
                )
                return
        with self.storageman:
            # if you made it this far, it is OK to add whatever you brought
            if receptor_blob:
                self.storageman.insert_receptor_blob(receptor_blob, receptor_name)
                LOGGER.info("Receptor pdbqt blob data was added to the database.")
            elif receptor_jsons:
                self.storageman.insert_receptor_polymer(receptor_jsons, receptor_name)
                LOGGER.info("Receptor polymer json string was added to the database.")

    # endregion

    # region filter, export

    def produce_summary(
        self, columns=["docking_score", "leff"], percentiles=[1, 10]
    ) -> None:
        """Print summary of data in storage to sdout

        Args:
            columns (list(str)): data columns used to prepare summary
            percentiles (list(int)): cutoff percentiles for the summary

        """
        with self.storageman:
            summary_data = self.storageman.fetch_summary_data(columns, percentiles)
        if "summary_data" not in locals():
            raise RTCoreError(f"Summary data is empty, please check the database.")
        print("Total Stored Ligands          :", summary_data.pop("num_ligands"))
        print("Total Stored Poses            :", summary_data.pop("num_poses"))
        print(
            "Total Unique Interactions     :",
            summary_data.pop("num_unique_interactions"),
        )
        print(
            "Number Interacting Residues   :",
            summary_data.pop("num_interacting_residues"),
        )

        colon_col = 18
        print("\nEnergy statistics:")
        print("=======================================")
        for col in columns:
            if col == "docking_score":
                min_e = summary_data["min_docking_score"]
                max_e = summary_data["max_docking_score"]
                print(f"Energy (min)      : {min_e:.2f} kcal/mol")
                print(f"Energy (max)      : {max_e:.2f} kcal/mol")
            elif col == "leff":
                min_le = summary_data["min_leff"]
                max_le = summary_data["max_leff"]
                print(f"LE     (min)      : {min_le:.2f} kcal/mol/heavyatom")
                print(f"LE     (max)      : {max_le:.2f} kcal/mol/heavyatom")
            else:
                min_col = summary_data[f"min_{col}"]
                max_col = summary_data[f"max_{col}"]
                print(f"{col} (min) : {min_col}")
                print(f"{col} (max) : {max_col}")
        if percentiles != [] and percentiles is not None:
            print("----------------------------------------")

            for col in columns:
                for p in percentiles:
                    if col == "docking_score":
                        p_string = f"Energy (top {p}% )"
                        p_string += " " * (colon_col - len(p_string))
                        print(
                            f"{p_string}: {summary_data[f'{p}%_docking_score']:.2f} kcal/mol"
                        )
                    elif col == "leff":
                        p_string = f"LE     (top {p}% )"
                        p_string += " " * (colon_col - len(p_string))
                        print(
                            f"{p_string}: {summary_data[f'{p}%_leff']:.2f} kcal/mol/heavyatom"
                        )
                    else:
                        p_string = f"{col} (top {p}%)"
                        p_string += " " * (colon_col - len(p_string))
                        print(f"{p_string}: {summary_data[f'{p}%_{col}']:.2f}")

    def filter(
        self,
        eworst=None,
        ebest=None,
        leworst=None,
        lebest=None,
        score_percentile=None,
        le_percentile=None,
        vdw_interactions=None,
        hb_interactions=None,
        reactive_interactions=None,
        hb_count=None,
        react_any=None,
        max_miss=0,
        ligand_name=None,
        ligand_operator=None,
        ligand_substruct=None,
        ligand_substruct_pos=None,
        ligand_max_atoms=None,
        ligand_min_molweight=None,
        ligand_max_molweight=None,
        # other processing options:
        enumerate_interaction_combs: bool = None,
        output_all_poses: bool = None,
        mfpt_cluster: float = None,
        interaction_cluster: float = None,
        log_file: str = None,
        outfields: str = list(RingtailDefaults.outfields),
        order_results: str = None,
        bookmark_name: str = RingtailDefaults.bookmark_name,
        filter_bookmark: str = None,
        return_iter: bool = False,
    ) -> Union[int, iter]:
        """Prepare list of filters, then filters and writes all passing pose_ids to a bookmark of given name. Creates an output log text file of all ligand (or all poses, if requested) docking results that passes.
        If clustering is requested, it will first perform the filtering, then cluster, and create a new bookmark of name <bookmark_<cluster_type>_cluster> for each requested cluster type.
        In the case of clustering, the representative ligands will be the ones written to the output log.

        Args:
            Filters:
                eworst (float): specify the worst energy value accepted
                ebest (float): specify the best energy value accepted
                leworst (float): specify the worst ligand efficiency value accepted
                lebest (float): specify the best ligand efficiency value accepted
                score_percentile (float): specify the worst energy percentile accepted. Express as percentage e.g. 1 for top 1 percent.
                le_percentile (float): specify the worst ligand efficiency percentile accepted. Express as percentage e.g. 1 for top 1 percent.
                vdw_interactions (list[tuple]): define van der Waals interactions with residue as [-][CHAIN]:[RES]:[NUM]:[ATOM_NAME]. E.g., [('A:VAL:279:', True), ('A:LYS:162:', True)] -> [('chain:resname:resid:atomname', <wanted (bool)>), ('chain:resname:resid:atomname', <wanted (bool)>)]
                hb_interactions (list[tuple]): define HB (ligand acceptor or donor) interaction as [-][CHAIN]:[RES]:[NUM]:[ATOM_NAME]. E.g., [('A:VAL:279:', True), ('A:LYS:162:', True)] -> [('chain:resname:resid:atomname', <wanted (bool)>), ('chain:resname:resid:atomname', <wanted (bool)>)]
                reactive_interactions (list[tuple]): check if ligand reacted with specified residue as [-][CHAIN]:[RES]:[NUM]:[ATOM_NAME]. E.g., [('A:VAL:279:', True), ('A:LYS:162:', True)] -> [('chain:resname:resid:atomname', <wanted (bool)>), ('chain:resname:resid:atomname', <wanted (bool)>)]
                hb_count (list[tuple]): accept ligands with at least the requested number of HB interactions. If a negative number is provided, then accept ligands with no more than the requested number of interactions. E.g., [('hb_count', 5)]
                react_any (bool): check if ligand reacted with any residue
                max_miss (int): Will compute all possible combinations of interaction filters excluding up to max_miss numer of interactions from given set. Default will only return union of poses interaction filter combinations. Use with 'enumerate_interaction_combs' for enumeration of poses passing each individual combination of interaction filters.
                ligand_name (list[str]): specify ligand name(s). Will combine name filters with OR, e.g., [["lig1", "lig2"]]
                ligand_substruct (list[str]): SMARTS, index of atom in SMARTS, cutoff dist, and target XYZ coords, e.g., [["ccc", "CN"]]
                ligand_substruct_pos (list[list[type]]): SMARTS pattern(s) for substructure matching, e.g., [["[Oh]C", 0, 1.2, -5.5, 10.0, 15.5]] -> [["smart_string", index_of_positioned_atom, cutoff_distance, x, y, z]]
                ligand_max_atoms (int): Maximum number of heavy atoms a ligand may have
                ligand_min_molweight (float): min molweight (inclusive, g/mol) for the ligand
                ligand_max_molweight (float): max molweight (inclusive, g/mol) for the ligand
                ligand_operator (str): logical join operator for multiple SMARTS (default: OR), either AND or OR
            Ligand results options:
                enumerate_interaction_combs (bool): When used with `max_miss` > 0, will log ligands/poses passing each separate interaction filter combination as well as union of combinations. Can significantly increase runtime.
                output_all_poses (bool): By default, will output only top-scoring pose passing filters per ligand. This flag will cause each pose passing the filters to be logged.
                mfpt_cluster (float): Cluster filtered ligands by Tanimoto distance of Morgan fingerprints with Butina clustering and output ligand with lowest ligand efficiency from each cluster. Default clustering cutoff is 0.5. Useful for selecting chemically dissimilar ligands.
                interaction_cluster (float): Cluster filtered ligands by Tanimoto distance of interaction fingerprints with Butina clustering and output ligand with lowest ligand efficiency from each cluster. Default clustering cutoff is 0.5. Useful for enhancing selection of ligands with diverse interactions.
                log_file (str): by default, results are saved in `output_log.txt`; if this option is used, ligands and requested info passing the filters will be written to specified file
                outfields (str): defines which fields are used when reporting the results (to stdout and to the log file); fields are specified as comma-separated values, e.g. `--outfields=e,le,hb`; by default, docking_score (energy) and ligand name are reported; ligand always reported in first column. Any column in Results and Ligands tables are available, although one of the following is recommended: \n
                    "Ligand_name" (Ligand name),
                    "e" (docking_score),
                    "le" (ligand efficiency),
                    "delta" (delta energy from best pose),
                    "ref_rmsd" (RMSD to reference pose),
                    "e_inter" (intermolecular energy),
                    "e_vdw" (van der waals energy),
                    "e_elec" (electrostatic energy),
                    "e_intra" (intermolecular energy),
                    "n_interact" (number of iteractions),
                    "ligand_smile" ,
                    "rank" (rank of ligand pose),
                    "run" (run number for ligand pose),
                    "hb" (hydrogen bonds),
                    "receptor" (receptor name)
                order_results (str): Stipulates how to order the results when written to the log file. By default will be ordered by order results were added to the database. ONLY TAKES ONE OPTION. Some available fields are shown, although any sortable column in Ligands or Results tables are OK:\n
                    "e" (docking_score),
                    "le" (ligand efficiency),
                    "delta" (delta energy from best pose),
                    "ref_rmsd" (RMSD to reference pose),
                    "e_inter" (intermolecular energy),
                    "e_vdw" (van der waals energy),
                    "e_elec" (electrostatic energy),
                    "e_intra" (intermolecular energy),
                    "n_interact" (number of interactions),
                    "rank" (rank of ligand pose),
                    "run" (run number for ligand pose),
                    "hb" (hydrogen bonds);
                bookmark_name (str): name for resulting book mark file. Default value is 'passing_results'
                filter_bookmark (str): name of bookmark to perform filtering over
                return_iter (bool): return an iterable of all of the filtering results

        Returns:
            int: number of ligands passing filter
            iter (optional): an iterable of all of the filtering results

        """

        filters = Filters(
            {
                "eworst": eworst,
                "ebest": ebest,
                "leworst": leworst,
                "lebest": lebest,
                "score_percentile": score_percentile,
                "le_percentile": le_percentile,
                "vdw_interactions": vdw_interactions,
                "hb_interactions": hb_interactions,
                "reactive_interactions": reactive_interactions,
                "hb_count": hb_count,
                "react_any": react_any,
                "max_miss": max_miss,
                "ligand_name": ligand_name,
                "ligand_operator": ligand_operator,
                "ligand_substruct": ligand_substruct,
                "ligand_substruct_pos": ligand_substruct_pos,
                "ligand_max_atoms": ligand_max_atoms,
                "ligand_min_molweight": ligand_min_molweight,
                "ligand_max_molweight": ligand_max_molweight,
            }
        )

        # Compatibility check with docking mode
        if react_any:
            with self.storageman as sm:
                if sm.get_previous_docking_mode() == "vina":
                    LOGGER.warning(
                        "Cannot use reaction filters with Vina mode. Removing react_any filter."
                    )
                react_any = False
        if not valid_bookmark_name(bookmark_name):
            raise OptionError(
                "Bookmark name contains illegal symbols, please only use letters, numbers, and underscore."
            )
        # make sure enumerate_interaction_combs always true if max_miss = 0, since we don't ever worry about the union in this case
        if max_miss == 0:
            enumerate_interaction_combs = True

        # guard against unsing percentile filter with all_poses
        if output_all_poses and not (score_percentile is None or le_percentile is None):
            LOGGER.warning(
                "Cannot return all passing poses with percentile filter. Will only log best pose."
            )
            output_all_poses = False
        # check that all filter values are valid
        filters.checks()

        LOGGER.info("Starting to filter results")
        # get possible permutations of interaction with max_miss excluded
        if enumerate_interaction_combs and max_miss > 0:
            interaction_combs = self._generate_interaction_combinations(
                filters.asdict(), max_miss
            )
            if len(interaction_combs) > 1:
                write_one_bookmark = False
            else:
                write_one_bookmark = True
                filter_dict = self._prepare_interaction_combo_filters(
                    filters.asdict(), combination
                )
        else:
            filter_dict = filters.asdict()
            write_one_bookmark = True

        # parse clustering
        clustering = {}
        if mfpt_cluster or interaction_cluster:
            if mfpt_cluster:
                clustering["mfp"] = mfpt_cluster
            if interaction_cluster:
                clustering["ifp"] = interaction_cluster

        if write_one_bookmark:
            with self.storageman:
                num_passing_ligands = self.storageman.filter_results(
                    all_filters=filter_dict,
                    bookmark_name=bookmark_name,
                    filtering_bookmark=filter_bookmark,
                )
            print_string = ""
        # else produce a bookmark for each interaction combination
        elif not write_one_bookmark:
            for ic_idx, combination in enumerate(interaction_combs):

                filters_copy = filters.asdict()
                # prepare Filter object with only desired interaction combination for storageManager
                temp_filters = self._prepare_interaction_combo_filters(
                    filters_copy, combination
                )
                iterated_bookmark_name = bookmark_name + "_" + str(ic_idx)
                # process each interaction combination filter
                with self.storageman:
                    num_passing_ligands = self.storageman.filter_results(
                        temp_filters,
                        iterated_bookmark_name,
                        filter_bookmark,
                    )
                if num_passing_ligands:
                    print(
                        f"\nNumber passing ligands filter combination {str(combination)}:",
                        num_passing_ligands,
                    )
                    if log_file:
                        self.write_filter_output(
                            iterated_bookmark_name,
                            temp_filters,
                            num_passing_ligands,
                            log_file,
                            outfields,
                            output_all_poses,
                            order_results,
                        )
                else:
                    print("\nNo ligands passing filters")
                    LOGGER.warning(
                        f"WARNING: No ligands found passing filter combination {str(combination)}."
                    )
                    with self.storageman:
                        self.storageman.delete_bookmark(iterated_bookmark_name)

            # Process the union of the max miss combinations of interactions
            with self.storageman:
                num_passing_ligands, bookmark_name = self.storageman.get_maxmiss_union(
                    len(interaction_combs), bookmark_name, filters_copy
                )
            # rename so can be reused in common method below
            filter_dict = filters_copy
            print_string = " in max_miss union"

        if num_passing_ligands:
            print(
                f"\nNumber passing ligands{print_string}:",
                num_passing_ligands,
            )
            if clustering:
                bookmark_name, num_passing_ligands = self._parse_clustering(
                    clustering, bookmark_name
                )
            # use original filters for final output log write
            if log_file:
                self.write_filter_output(
                    bookmark_name,
                    filters.asdict(),
                    num_passing_ligands,
                    log_file,
                    outfields,
                    output_all_poses,
                    order_results,
                    clustering,
                )
        else:
            LOGGER.warning(f"WARNING: No ligands found passing filters{print_string}.")
            print(f"\nNo ligands passing filters{print_string}")

        if return_iter:
            with self.storageman:
                formatted_query = self.storageman.get_bookmark_selection(
                    bookmark_name, outfields, not output_all_poses, order_results
                )
                # have to make it into list of tuples since storageman uses row factory
                iter = [
                    tuple(row)
                    for row in self.storageman.db_query(formatted_query).fetchall()
                ]
            return iter
        else:
            return num_passing_ligands

    def cluster(
        self, bookmark_name: str, type: str = "mfp", cutoff: float = 0.5
    ) -> str:
        """
        Clusters data in a bookmark using given cluster algorithm and with given cutoff distance.

        Args:
            bookmark_name (str): bookmark name that describes what data to cluster
            type (str, optional): cluster algirithm. Defaults to "mfp".
            cutoff (float, optional): cutoff distance for each cluster. Defaults to 0.5.

        Returns:
            str: bookmark name describing the representative poses for the clustered data
        """
        with self.storageman as sm:
            bookmark_name, num_clusters = sm.cluster_data(
                bookmark_name,
                type,
                cutoff,
            )
        LOGGER.info(
            f"Number of clusters: {num_clusters}.\nPassing poses saved to {bookmark_name}."
        )
        return bookmark_name

    def write_filter_output(
        self,
        bookmark_name: str,
        filters: dict,
        num_passing_ligands: int,
        log_file: str = RingtailDefaults.log_file,
        output_fields: Union[str, list] = RingtailDefaults.outfields,
        output_all_poses: bool = RingtailDefaults.output_all_poses,
        order_results: str = RingtailDefaults.order_results,
        clustering: dict = {},
        max_miss_combs: int = None,
        append_to_log: bool = None,
    ):
        """
        Args:
            bookmark_name (str): describing the passing poses
            filters (dict): filters used to pass the poses
            num_passing_ligands (int): number passing ligands
            log_file (str): name of file to write to
            output_fields (Union[str, list]): output columns to include
            output_all_poses (bool): whether or not to output one or all passing poses
            order_results (str): if ordering the results by a specific column
            clustering (dict, optional): dict describing any clustering done. Defaults to {}.
            max_miss_combs (int, optional): how many max_miss combinations were used. Defaults to None.
            append_to_log (bool, optional): whether or not to append to existing logfile. Defaults to None.

        Raises:
            OptionError
        """
        if not order_results:
            order_results = "ligname"
        with self.storageman:
            formatted_query = self.storageman.get_bookmark_selection(
                bookmark_name, output_fields, not output_all_poses, order_results
            )

        # format output_log string
        if clustering:
            cluster_string = f"Morgan Fingerprints butina clustering cutoff: {clustering.get('mfp')}\nInteraction Fingerprints clustering cutoff: {clustering.get('ifp')}"
        else:
            cluster_string = "No clustering performed\n."
        with OutputManager(log_file, append_to_log) as opm:
            if max_miss_combs:
                opm.write_maxmiss_union_header()
                opm.write_bookmarkname_in_log(bookmark_name)
            else:
                opm.write_filtervalues_in_log(
                    filters, [], bookmark_name, cluster_string
                )
            with self.storageman:
                opm.write_filter_results_in_log(
                    self.storageman.db_query(formatted_query).fetchall()
                )
            opm.log_num_passing_ligands(num_passing_ligands)

    def write_flexres_pdb(
        self,
        receptor_polymer: Union[object, str] = None,
        ligname: Union[str, list] = None,
        bookmark_name: str = None,
        filename: str = "receptor.pdb",
        get_consent: Callable = None,
    ) -> tuple[Chem.rdchem.Mol, dict]:
        """
        Writes a receptor pdb with flexible residues based on the ligand provided

        Args:
            receptor_polymer (Polymer, json string, .json file): version of receptor produced by meeko,
                                        or a json representation that can be rebuilt to a polymer, including
                                        valid json string and valid json file
            ligname (str | list): ligand name for which the receptor flexible residue info should be collected
            filename (str): name of the output pdb, extension is optional, will default to '.pdb'
            bookmark_name (str, optional): if provided, it will only export flex res for ligand poses passing bookmark filters

        Returns
            Chem.rdchem.Mol: the Mol object for the given ligand
            dict: dictionary describing the Mols for the flexible residues
        """
        # check database if polymer not provided, else error
        if not receptor_polymer:
            _, _, polymer_string = self.get_receptor_object()
            if not polymer_string:
                raise OptionError(
                    "No receptor polymer provided and no polymer found in the database. Cannot create flexres PDBs."
                )
            else:
                receptor_polymer = polymer_string

        from meeko.export_flexres import pdb_updated_flexres_from_rdkit

        ligands_poses = self._fetch_select_ligands_poses(
            ligand_names=ligname, bookmark_name=bookmark_name
        )
        if len(ligands_poses) > 10:
            length = len(ligands_poses)

            def _api_cmd_consent(length: int):
                consent = input(
                    f"WARNING: You are requesting to prepare pdbs for {length} ligands/flexible residue combinations.\n Are you sure you wish to proceed? If so, type 'yes': "
                )
                return consent.strip().lower() == "yes"

            if get_consent is None:
                get_consent = _api_cmd_consent

            if not get_consent(length):
                LOGGER.critical(
                    "Consent not given for large number of pdbs to write, exiting."
                )
                return

        # ensure polymer is correct type, else check if string or path
        if type(receptor_polymer) == str:
            from meeko import Polymer

            _, polymer_string = self._process_receptor_polymer(receptor_polymer)
            polymer = Polymer.from_json(polymer_string)
        else:
            polymer = receptor_polymer

        flexres_data = self.make_receptor_flexres_mols()
        lig_flex_mol = {}
        for ligand, poses in ligands_poses.items():
            ligand_mol, flexres_mols, _, flexible_residues = self.create_rdkit_mol(
                ligand, poses, flexres_data
            )

            if filename:
                # if providing filename, make sure it has .pdb extension
                # and add ligand name
                root, ext = os.path.splitext(filename)
                if not ext:
                    ext = ".pdb"
                path = root + ligand + ext
            else:
                receptor_name, _, _ = self.get_receptor_object()
                path = receptor_name + ligand + ".pdb"

            flexmoldict = {}
            # string in list of strings
            for index, flexres in enumerate(flexible_residues):
                # res id is chain:resnum
                flexmoldict[f"{flexres[4]}:{flexres[-3:]}"] = flexres_mols[index]

            pdb_str = pdb_updated_flexres_from_rdkit(polymer, flexmoldict)
            # write pdb string to file
            with open(path, "w") as file:
                file.write(pdb_str)

            lig_flex_mol[ligand] = {
                "ligand_mol": ligand_mol,
                "flexmoldict": flexmoldict,
            }

        return lig_flex_mol

    def write_molecule_sdfs(
        self,
        bookmark_name: str = None,
        sdf_path: str = ".",
        all_in_one: bool = True,
        ligname: Union[str, list] = None,
    ):
        """
        Have output manager write molecule sdf files for passing results in given results bookmark

        Args:
            bookmark_name (str): bookmark name from which to export Mols
            sdf_path (str, optional): Optional path existing or to be created in cd where SDF files will be saved
            all_in_one (bool, optional): If True will write all molecules to one SDF (separated by $$$$), if False will write one molecule pre SDF
            ligname (str | list): ligand name for which the receptor flexible residue info should be collected
        Raises:
            StorageError: if bookmark or data not found
        """

        try:
            ligands_poses = self._fetch_select_ligands_poses(
                ligand_names=ligname, bookmark_name=bookmark_name
            )
            flexres_data = self.make_receptor_flexres_mols()
            all_mols = {}
            for ligname, poses in ligands_poses.items():

                mol, flexres_mols, properties, _ = self.create_rdkit_mol(
                    ligname,
                    poses,
                    flexres_data,
                )
                all_mols[ligname] = {
                    "ligand": mol,
                    "flex_residues": flexres_mols,
                    "properties": properties,
                }

        except StorageError as e:
            LOGGER.error(str(e))
            return

        if all_mols is None:
            LOGGER.error(
                f"Selected bookmark {bookmark_name} does not exist or does not have any data, cannot write molecule SDFs."
            )
            return

        for ligname, info in all_mols.items():
            # determine filename
            if all_in_one:
                # will write one SDF file for all molecules in bookmark (_None if no bookmark present)
                db_file_name = os.path.splitext(os.path.basename(self.db_file))[0]
                sdf_file_name = f"{db_file_name}_{bookmark_name}.sdf"
                LOGGER.info(f"Writing {ligname} to {sdf_file_name}")
            else:
                # filename is name of ligand
                sdf_file_name = ligname + ".sdf"
                LOGGER.info("Writing " + ligname + ".sdf")

            if (
                sdf_path is not None
                and not sdf_path == ""
                and not os.path.isdir(sdf_path)
            ):
                os.makedirs(sdf_path)
                LOGGER.info(
                    "Specified directory for SDF files was created in current working directory."
                )
            with OutputManager() as opm:
                opm.write_out_mol(
                    sdf_file_name,
                    info["ligand"],
                    info["flex_residues"],
                    info["properties"],
                    sdf_path,
                )

    def find_similar_ligands(
        self, query_ligname: str, log_file: str = "cluster_log.txt"
    ) -> int:
        """
        Find ligands in cluster with query_ligname

        Args:
            query_ligname (str): name of the ligand in the ligand table to look for similars to

        Returns:
            int: number of ligands that are similar
        """

        number_similar = 0
        with self.storageman:
            similar_ligands, bookmark_name, cluster_name = (
                self.storageman.fetch_clustered_similars(query_ligname)
            )
            if similar_ligands is not None:
                number_similar = len(similar_ligands)
                with OutputManager(log_file) as opm:
                    opm.write_find_similar_header(query_ligname, cluster_name)
                    opm.write_bookmarkname_in_log(bookmark_name)
                    opm.write_filter_results_in_log(similar_ligands)
                    opm.log_num_passing_ligands(number_similar)
            print(
                "\nNumber of similar ligands:",
                number_similar,
            )
        return number_similar

    def export_columns_as_csv(
        self,
        columns: list[str],
        bookmark_name: str,
        csv_name: str = "columns_output.csv",
    ):
        """
        Exports specified columns from Results or bookmark (or status table) as csv

        Args:
            columns (list[str]): list of string columns to be exported
            bookmark_name (str): data limiting table/bookmark to export for
            csv_name (str, optional): name of csv file. Defaults to "columns_output.csv".
        """
        with self.storageman as sm:
            sql = sm.prepare_column_export_query(columns, bookmark_name)

        self._export_csv(sql, csv_name, False)

    def export_table_as_csv(self, table: str, csv_name: str = "table_output.csv"):
        """
        Exports specified table (all of it) as a csv file

        Args:
            table (str): table to be exported
            csv_name (str, optional): name of csv file. Defaults to "table_output.csv".
        """
        self._export_csv(table, csv_name, True)

    def export_sql_as_csv(self, sql: str, csv_name: str = "sql_output.csv"):
        """
        Exports sql-formatted query as csv, assumes user has used appropriate sql terminology for their database engine.

        Args:
            sql (str): query as a string
            csv_name (str, optional): name of csv file. Defaults to "sql_output.csv".
        """
        self._export_csv(sql, csv_name, False)

    def export_bookmark_db(self, bookmark_name: str, db_filepath: str = None) -> str:
        """Export database containing data from bookmark

        Args:
            bookmark_name (str): name for bookmark_db

        Returns:
            str: name of the new, exported database
        """
        if db_filepath is None:
            db_filepath = self.db_file.rstrip(".db") + "_" + bookmark_name + ".db"

        LOGGER.info("Exporting bookmark database")
        if os.path.exists(db_filepath):
            LOGGER.warning(
                "Requested export DB name already exists. Please rename or remove existing database. New database not exported."
            )
            return
        with self.storageman:
            self.storageman.create_subset_database(bookmark_name, db_filepath)

        return db_filepath

    def export_receptor(self):
        """
        Export receptor in database to pdbqt
        """
        # TODO export pdbqt, need a differnet method for json
        recname, recstr, recjson = self.get_receptor_object()
        if recstr is None:
            LOGGER.warning(f"No receptor pdbqt stored for {recname}. Export failed.")
            return
        output_manager = OutputManager()
        output_manager.write_receptor_pdbqt(recname, recstr)

    def get_previous_filter_data(
        self,
        bookmark_name: str,
        outfields: str = "Ligname,docking_score",
        order_results: str = "ligname",
        log_file: str = None,
    ):
        """Get data requested in outfields from the
        results bookmark of a previous filtering

        Args:
            bookmark_name (str): bookmark for which the filters were used
            outfields (str, optional): columns to write to the output. Defaults to "Ligname,docking_score".
            order_results (str, optional): how to order the data. Defaults to "ligname".
            log_file (str, optional): log file name to write to. Defaults to None.
        """
        if log_file is None:
            log_file = "output_log.txt"
        if bookmark_name is None:
            raise OptionError("A bookmark name has to be provided")
        with self.storageman:
            new_data = self.storageman.fetch_data_for_passing_results(
                bookmark_name, outfields, order_results
            )
            with OutputManager(log_file) as opm:
                opm.write_filter_results_in_log(new_data)

    def cross_reference_databases(
        self,
        wanted_dbs: list[tuple[str, Union[str, None]]] = None,
        unwanted_dbs: list[tuple[str, Union[str, None]]] = None,
        bookmark_prefix: str = "crossref",
    ) -> tuple[int, dict, list, dict]:
        """
        Method to cross reference two or more databases. Will attach all other databases
        to the "current" database, which is always the first wanted database (will have the
        active ringtail object).

        Will create an intersect of ligand names in the wanted databases+bookmarks,
        and delete any ligands found in the union of the unwanted databases+bookmarks.

        A bookmark with specified prefix to the original bookmark name will be stored in each database,
        making the crossreferencing data easily available for later use.

        Args:
            wanted_dbs (list[tuple[str, Union[str, None]]], optional): _description_. Defaults to None.
            unwanted_dbs (list[tuple[str, Union[str, None]]], optional): _description_. Defaults to None.
            bookmark_prefix (str, optional): _description_. Defaults to "crossref".

        Returns:
            tuple[int, dict, list, dict]: Number of "passing" ligands, dict of database {database
            file path: new bookmark name from cross ref}, list of passing ligand names, dict of dbs and how they were compared
        """
        with self.storageman as sm:
            return sm.crossref_databases(wanted_dbs, unwanted_dbs, bookmark_prefix)

    def get_gui_plot_data(
        self,
        bookmark_name: str = None,
        include_status: bool = False,
        x_axis: str = "docking_score",
        y_axis: str = "leff",
        limit: int = None,
    ):
        """
        Data from Ringtail formatted to work with the GUI plotting tool. Will most likely depreceate the get_plot_data and functionality
        outside of GUI, and change name of this method accordingly

        Args:
            bookmark_name (str, optional): _description_. Defaults to None.
            include_status (bool, optional): _description_. Defaults to False.
            x_axis (str, optional): _description_. Defaults to "docking_score".
            y_axis (str, optional): _description_. Defaults to "leff".
            limit (int, optional): _description_. Defaults to None.

        Returns:
            _type_: _description_
        """

        with self.storageman:
            return self.storageman.get_gui_plot_data(
                bookmark_name, include_status, x_axis, y_axis, limit
            )

    def get_plot_data(
        self,
        bookmark_name: str = None,
        only_passing: bool = False,
        include_status: bool = False,
        x_axis: str = "docking_score",
        y_axis: str = "leff",
        limit: int = None,
    ) -> tuple[
        list[tuple[float, float]], list[tuple[float, float, int, str, Union[None, str]]]
    ]:
        """
        Get ligand efficiency and energy for all docking data and for ligands that passed
        filtering in specified bookmark. Each tuple in the respective lists contains
        docking_score, leff, pose_id, and ligand name. If requested, will also include
        pose_id status if previously assigned

        Args:
            bookmark_name (str): name for which to fetch passing poses
            include_status (bool, optional): Whether to include pose_id status in the data (a 5th column). Defaults to False

        Returns:
            list(tuple), list(tuple): [all_data <only docking score and ligand efficiency], [filtered_data <also includes pose id and ligand name>]
        """
        with self.storageman:
            all_data, passing_data = self.storageman.get_plot_data(
                bookmark_name, only_passing, include_status, x_axis, y_axis, limit
            )

        return all_data, passing_data

    def get_receptor_object(self, pdbqt_str_only: bool = False) -> tuple:
        """
        Gets the receptor object from the database

        Args:
            pdbqt_str_only (bool, optional): if json is in db, return it as pdbqt str instead of json string. Defaults to False.

        Returns:
            tuple[str,str, str]: receptor name and receptor blob and receptor serialized polymer as strings
        """
        with self.storageman as sm:
            rec_data = sm.fetch_receptor_object()
            if pdbqt_str_only and (polymer_json := rec_data.get("polymer")):
                polymer = ReceptorManager.polymer_json2pdbqt_str(polymer_json)
            else:
                polymer = rec_data.get("polymer")
            if rec_data is not None:
                return (
                    rec_data.get("RecName"),
                    ReceptorManager.blob2str(rec_data.get("receptor_object")),
                    polymer,
                )
            else:
                return None, None, None

    # endregion

    # region visualize ringtail

    def plot(
        self, bookmark_name: str, save: bool = True, return_fig_handle: bool = False
    ) -> Union[None, plt.Figure]:
        """
        Get data needed for creating Ligand Efficiency vs
        Energy scatter plot from storageManager. Call OutputManager to create plot.
        Option to save the plot and close it immediately, or keep it open and save it manually later.

        Args:
            bookmark_name (str): bookmark from which to fetch filtered data to plot
            save (bool): whether to save plot to cd. Will save and close figure
            return_fig_handle (bool): use to return a handle to the matplotlib figure instead of saving or showing figure

        Returns:
            matplotlib.pyplot.figure (optional): will not show figure if returning figure handle
        """
        with self.storageman:
            bookmark_filters = self.storageman.fetch_filters_from_bookmark(
                bookmark_name
            )  # fetches the filters used to produce the bookmark

        if bookmark_filters:
            max_miss = bookmark_filters["max_miss"]
            if max_miss > 0:
                raise OptionError(
                    "Cannot use --plot with --max_miss > 0. Can plot for desired bookmark with --bookmark_name."
                )

        LOGGER.info("Creating plot of results")
        all_data, passing_data = self.get_plot_data(bookmark_name)
        xdata = []
        ydata = []
        # add to list as docking_score/energy and ligand_efficiency
        for line in all_data:
            # handle empty db rows
            if None in line:
                continue
            xdata.append(line[0])
            ydata.append(line[1])

        # base number of bins on data size
        datalength = len(xdata)
        # scale some plot parameters to size of dataset
        if datalength > 1000:
            num_of_bins = 100
            markersize = 20
        # for smaller dataset, scale num of bins and markersize to size of dataset
        else:
            num_of_bins = max(1, round(datalength / 10))
            markersize = 60 - (datalength / 25)

        # plot the data
        output_manager = OutputManager()
        fig = output_manager.plot_all_data(xdata, ydata, num_of_bins)

        if passing_data != []:  # handle if no passing ligands
            xaxis = []
            yaxis = []
            for line in passing_data:
                # energy
                xaxis.append(line[0])
                # leff
                yaxis.append(line[1])

            output_manager.plot_single_points(xaxis, yaxis, markersize)

        if save:
            output_manager.save_scatterplot()

        if return_fig_handle:
            return fig
        else:
            plt.show()

    def display_pymol(self, bookmark_name, viewer_name: str = "pymol"):
        # TODO update to new gui viewer paradigm
        """
        Launch pymol session and plot of LE vs docking score. Displays molecules when clicked.

        Args:
            bookmark_name (str): bookmark name to use in pymol. 'None' uses the whole db?
        """
        available_viewers = {
            "pymol": self.launch_pymol,
            "molview": self.launch_molviewer,
        }
        if viewer_name.lower() not in available_viewers.keys():
            raise OptionError(
                f"Selected viewer type {viewer_name} not in list of available viewers: {available_viewers.keys()}. Please pick a valid molecular viewer."
            )

        viewer = available_viewers[viewer_name]()
        self.make_clickable_plot(viewer, bookmark_name, viewer_name)

    def launch_molviewer(self):
        # placeholder for forli lab mol viewer
        return True

    def launch_pymol(self):
        """doc string"""
        import subprocess
        from rdkit.Chem import PyMol
        import socket
        import time

        # launch pymol session
        process = subprocess.Popen(
            ["pymol", "-R"], stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )

        # ensure pymol was opened
        def _is_pymol_running(host="localhost", port=9123):
            """Check if PyMOL's remote server is responding."""
            with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
                return sock.connect_ex((host, port)) == 0

        # Wait until PyMOL is fully loaded
        timeout = 10  # max wait time in seconds
        start_time = time.time()

        while not _is_pymol_running():
            if time.time() - start_time > timeout:
                break
            time.sleep(0.5)

        if not _is_pymol_running():
            LOGGER.error("PyMOL failed to start.")

        try:
            pymol = PyMol.MolViewer()
        except ConnectionRefusedError as e:
            raise RTCoreError(
                "Error establishing connection with PyMol. Try manually launching PyMol with `pymol -R` in another terminal window."
            ) from e
        return pymol, process

    def make_clickable_plot(
        self,
        viewer,
        bookmark_name: str,
        viewer_type: str,
        canvas=None,
    ):
        """should be handed a molecular viewer, plot requested data, and
        have the plotted data be clickable which click will display them
        in the molecular viewer"""
        # NOTE not being maintained as of 2025/10

        if hasattr(self, "cid"):
            # disconnect any old connections
            canvas.mpl_disconnect(self.cid)

        # in case plot exist
        if canvas is None:
            axes = plt.axes()
        else:
            axes = canvas.axes

        poseIDs = {}
        with self.storageman:
            # fetch data for passing ligands
            _, passing_data = self.storageman.get_plot_data(
                bookmark_name, only_passing=True
            )

        for line in passing_data:
            axes.plot(line[0], line[1], ".r", mfc="None", picker=5)
            poseIDs[(line[0], line[1])] = (
                line[2],
                line[3],
            )  # line[0] is LE, line[1] is docking score, line[2] is pose_ID, line[3] is LigName
        axes.set_ylabel("Ligand Efficiency (kcal/mol/heavy atom)")
        axes.set_xlabel("Docking Score (kcal/mol)")
        axes.set_title("Passing Docking Poses")

        # check if receptor in db
        receptor = self.get_receptor_object()
        if receptor[1]:

            # load receptor if it exist in database
            rec_name = receptor[0]
            rec_string = receptor[1]

            rec_string_list = rec_string.split("\n")
            if viewer_type.lower() == "pymol":
                import tempfile

                rec_string = ReceptorManager.blob2str(receptor[1])
                # with the rdkit pymol api, easiest to read receptor from file
                with tempfile.NamedTemporaryFile(suffix=".pdbqt") as temp_file:
                    temp_file.write(rec_string.encode("utf-8"))
                    temp_file.flush()
                    temp_file_path = temp_file.name
                    viewer.LoadFile(temp_file_path, rec_name)
                    # center view on receptor
                    viewer.server.do("zoom")
            else:
                viewer.addpdbqt(
                    rec_string_list,
                    name=rec_name,
                    parent=None,
                    metadata=None,
                    assign_radii_and_bonds=True,
                )
                o = viewer.objects
                o.summary
                viewer.showspheres(o[0])
                viewer.autozoom()
        else:
            LOGGER.debug(
                "No receptor information in the database, receptor will not be displayed."
            )

        def _onpick(event):
            # get info about the point
            line = event.artist
            # coordinates is x (leff) and y (e) axis
            coords = tuple([c[0] for c in line.get_data()])
            chosen_pose = poseIDs[coords]
            ligname = chosen_pose[1]
            pose_id = chosen_pose[0]
            LOGGER.info(f"LigName: {ligname}; Pose_ID: {pose_id}")

            # make rdkit mol for poseid
            mol, flexres_mols, _, flexible_residues = self.create_rdkit_mol(
                ligname, [pose_id]
            )

            # update viewer
            if viewer_type.lower() == "pymol":
                viewer.ShowMol(mol, name=ligname, showOnly=False)
                for idx, resmol in enumerate(flexres_mols):
                    viewer.ShowMol(
                        resmol,
                        name=ligname + "_" + flexible_residues[idx],
                        showOnly=False,
                    )
            else:
                o = viewer.objects
                # this prop is required by _load_rdkit
                mol.SetProp("_Name", f"{Chem.MolToSmiles(mol)}")
                metadata = {"source": f"name::{Chem.MolToSmiles(mol)}"}

                viewer.hidesticks(o[-1])
                viewer.addrdkit(mol, "ligand smiles:" + Chem.MolToSmiles(mol))
                viewer.showsticks(o[-1])
                viewer.autozoom(o[-1])
                viewer.colorbyelement(o[-1], carbon_color="grey")

        if not canvas:
            fig = plt.gcf()
            self.cid = fig.canvas.mpl_connect("pick_event", _onpick)
            plt.show()
        else:
            # the connection ID terminates once a new plot is made, need to keep track of
            self.cid = canvas.mpl_connect("pick_event", _onpick)
            # need to return updated mol
            return canvas

    # endregion

    # region utilities

    def db_query(self, query: str, params=(), commit=False) -> Union[None, iter]:
        """
        Run a db query and return iterable if applicable

        Args:
            query (str): sql formatted string
            params (tuple, optional): parameters, assumes placehodlers in query

        Returns:
            iter: if any
        """
        with self.storageman as sm:
            return sm.db_query(query, params, commit).fetchall()

    def set_ligand_status(
        self,
        status: str,
        ligands: Union[str, list[str]] = None,
        pose_ids: Union[int, list[int]] = None,
        bookmark_name: str = None,
    ):
        """
        Will update the status of specified ligand-pose-bookmark combination. E.g., if only a ligand is provided, all its poses will be given
        the desired status. If just poses are provided, they will be given chosen status. If only a bookmark is provided, all poses in that
        bookmark will be given the status. If one or more ligands are given + a bookmark, the poses of those ligands in that bookmark will
        be given the status. If pose_ids are also specified, they will be given the set status even if not in the bookmark.

        Args:
            status (str): _description_
            ligand_name (Union[str, list[str]], optional): _description_. Defaults to None.
            pose_id (Union[int, list[int]], optional): _description_. Defaults to None.
            bookmark_name (str, optional): _description_. Defaults to None.
        """
        status = status.lower()

        ligands_poses = self._fetch_select_ligands_poses(
            ligands, pose_ids, bookmark_name
        )
        # parse to just poses, I don't need ligand names for this
        parsed_poses = [
            pose for pose_list in ligands_poses.values() for pose in pose_list
        ]
        self.update_pose_status(pose_id=parsed_poses, status=status)

    def update_pose_status(self, pose_id: Union[int, list[int]], status: str):
        """
        Updates the status of a given pose so that it is only ever in one status table

        Args:
            pose_id (Union[int, list[int]]): pose id for which a status is assigned, can be a single pose or a list
            status (str): new status

        Raises:
            OptionError
        """
        status = status.lower()
        if status not in statuses:
            raise OptionError(
                f"""The selected status {status} is not a valid option for Ringtail. Please choose between {','.join(statuses)}"""
            )

        with self.storageman as sm:
            if status == "accepted":
                sm.accept_pose(pose_id)
            elif status == "maybe":
                sm.maybe_pose(pose_id)
            elif status == "rejected":
                sm.reject_pose(pose_id)

    def get_data_table_pointer(
        self, table: str, length: int = 100, starting_pose_id: int = 0
    ) -> iter:
        """
        Returns a pointer or cursor (iterable) to the data in the Results table,
        if table is a bookmark it will limit the Results data to the poses represented
        in the bookmark.

        Args:
            table (str): name of table or bookmark
            length (int, optional): number of rows to collect. Defaults to 100.
            starting_pose_id (int, optional): pose id to start with. Defaults to 0.

        Returns:
            iter: iterable/cursor pointing to data in given table or bookmark
        """
        with self.storageman as sm:
            return sm.fetch_viewable_data_columns_from(table, length, starting_pose_id)

    def get_exportable_columns(self) -> dict[(str, list)]:
        """
        Creates a dict of useful export columns from Results, Ligands, and Interactions

        Returns:
            dict: table: [columns]
        """
        with self.storageman as sm:
            return sm.get_useful_columns()

    def get_row_from_pose(self, table: str, pose_id: int) -> int:
        """
        Method to get table row for a a pose id in given table, used to support gui

        Args:
            table (str): table name
            pose_id (int):

        Returns:
            int: row in given table
        """
        with self.storageman as sm:
            return sm.pose_row_in_table(table, pose_id)

    def get_limited_table_data(
        self, table: str, length: int = 100, starting_row_id: int = 1, reverse=False
    ) -> dict[list[str], list]:
        """
        Returns a pointer or cursor (iterable) to the data in the Results table,
        if table is a bookmark it will limit the Results data to the poses represented
        in the bookmark.

        Args:
            table (str): name of table or bookmark
            length (int, optional): number of rows to collect. Defaults to 100.
            starting_row_id (int, optional): rowid to start fetching from. Defaults to 0.
            reverse (bool, optional): whether to collect data in descending/reverse order. Defaults to False

        Returns:
            dict[list[str], list]: dict of headers and data
        """
        with self.storageman as sm:
            return sm.fetch_viewable_data_columns_from(
                table, length, starting_row_id, reverse
            )

    def get_table_columns(
        self,
        table: str,
        columns: list = ["*"],
        length: int = None,
        starting_rowid: int = 0,
    ) -> tuple[list[str], list[dict]]:
        """
        Will get requested table data for a table given one or more columns.
        Data will be limited by a certain length, and can be retrieved from a desired
        rowid.

        Args:
            table (str): name of table or bookmark
            columns (list, optional): list of columns to retrieve. Defaults to ["*"].
            length (int, optional): number of rows to collect. Defaults to 100.
            starting_rowid (int, optional): rowid to start with. Defaults to 0.

        Returns:
            tuple[list[str], list[dict]]: list of column names, and list of dicts where each dict is one row,
                                            and column is the key, value is the row-col cell value
        """
        with self.storageman as sm:
            return sm.fetch_columns_from_table_as_dicts(
                table, columns, length, starting_rowid
            )

    def get_query_data(self, query: str) -> tuple[list[str], list[dict]]:
        """
        Will return data requested in an sql formatted query. User must ensure the sql
        string is formatted correctly for the type of database used.

        Args:
            query (str): sql query formatted to sqlite database

        Returns:
            tuple[list[str], list[dict]]: list of column names, and list of dicts where each dict is one row,
                                            and column is the key, value is the row-col cell value
        """

        with self.storageman as sm:
            return sm.get_query_data_as_dicts(query)

    def table_length(self, table: str) -> int:
        """
        Get length of table or bookmark

        Args:
            table (str): name of table or bookmark

        Returns:
            int: number of poses in table/bookmark
        """

        with self.storageman:
            return self.storageman.table_length(table)

    def get_starting_rowid(self, table: str):
        """
        Returns the starting row id for a table. For a normal table like Results or status table, this will be
        1, but for a bookmark, this will be the rowid it start with in the Filtered_poses table

        Args:
            table (str): table or bookmark name

        Returns:
            int: rowid
        """
        with self.storageman:
            return self.storageman.get_starting_rowid(table)

    def database_is_empty(self) -> bool:
        """
        Determines whether or not there is data in the database

        Returns:
            bool: True if db is empty
        """
        with self.storageman as sm:
            return sm.db_empty()

    def get_range_of_e_le(self, table: str = "Results") -> tuple[float, float]:
        """
        Get min and max of e/docking_score and ligand efficiency/le/leff/

        Args:
            table (str): table limit data, e.g., either Results or a bookmark name

        Returns:
            tuple: e_min, e_max, le_min, le_max
        """
        with self.storageman as sm:
            return sm.get_range_of_e_le(table)

    def get_percentiles(
        self, column: str, num_bins: int = 20, table: str = "Results"
    ) -> tuple[list[int], list[float]]:
        """
        Will calculate percentiles for a given column and given number of bins to divide the data in.
        Will group the data by ligand_id, so it will be per ligand and not per pose id.

        Args:
            column (str): what column to calculate percentile for. must be numeric
            num_bins (int, optional): how many percentile bins data should be divided into. Defaults to 20.
            table (str, optional): whether the column is in Results or filtered results (i.e., bookmark). Defaults to "Results".

        Returns:
            tuple[list[int],list[float]]: list of percentiles as bins, and list of edge of each bin
        """
        with self.storageman as sm:
            return sm.calculate_percentiles(column, num_bins, table)

    def all_database_tables(self) -> list:
        """
        Returns a list of all table names in the database

        Returns:
            list: list of table names
        """
        with self.storageman as sm:
            return sm.tables_in_db()

    def delete_bookmark(self, bookmark_name: str):
        """Drops specified bookmark from the database

        Args:
            bookmark_name (str): name of bookmark to be dropped.
        """

        with self.storageman:
            self.storageman.delete_bookmark(bookmark_name)
        LOGGER.info(
            f"Bookmark {bookmark_name} was dropped from the database {self.db_file}"
        )

    def get_bookmark_names(self) -> list:
        """
        Method to retrieve all bookmark names in a database

        Returns:
            list: of all bookmarks in a database
        """
        with self.storageman:
            return self.storageman.get_all_bookmark_names()

    def get_bookmark_interactions(self, bookmark_name: str):
        """
        Get all interactions represented by the poses in a bookmark
        (could also be status table)

        Args:
            bookmark_name (str): _description_

        Returns:
            dict: key is column name,
        """
        with self.storageman as sm:
            return sm.fetch_bookmark_interactions(bookmark_name)

    def create_rdkit_mol(
        self,
        ligname: str,
        pose_ids: list = None,
        flexres_data: tuple = None,
    ) -> tuple:
        """Creates rdkit molecule for given ligand, for all included pose_ids

        Args:
            ligname (string): ligand name
            pose_ids (list[int], optional): list of poses to consider for the ligand
            flexres_data (tuple, optional): flexible residue data if any. Defaults to None.

        Raises:
            OutputError: raises error if there is an issue with determining a flexible residue identity

        Returns:
            tuple: (ligand rdkit mol, [flexres rdkit mols], {properties for ligand})


        Note: needs to be ran inside a storageman context manager, will not be able to access the temporary table otherwise.
        """
        # get rdkit data
        with self.storageman as sm:
            rdmol, atom_index_map, hydrogen_parents = (
                sm.fetch_ligand_rdkit_relevant_info(ligname)
            )

        mol = Chem.Mol(rdmol)
        atom_indices = json.loads(atom_index_map)
        ligand_saved_coords = []

        if flexres_data:
            # unpack flexres data
            flexres_mols, flexres_info, flexres_saved_coords, flexres_residues = (
                flexres_data
            )
        else:
            flexres_mols, flexres_info, flexres_saved_coords, flexres_residues = (
                [],
                [],
                [],
                [],
            )

        properties = {
            "Binding energies": [],
            "Ligand effiencies": [],
            "Interactions": [],
        }

        if type(pose_ids) == int:
            pose_ids = [pose_ids]

        # update mols with additional data from select poses
        (
            mol,
            flexres_mols,
            ligand_saved_coords,
            flexres_saved_coords,
            properties,
        ) = self._update_mols_with_poses(
            pose_ids,
            atom_indices,
            mol,
            flexres_mols,
            flexres_info,
            ligand_saved_coords,
            flexres_saved_coords,
            properties,
        )

        # add ligand name to properties
        properties["_Name"] = ligname
        # add hydrogens to mols
        lig_h_parents = [int(idx) for idx in json.loads(hydrogen_parents)]
        mol = RDKitMolCreate.add_hydrogens(mol, ligand_saved_coords, lig_h_parents)
        flexres_hparents = []
        for idx, res in enumerate(flexres_mols):
            flexres_hparents = flexres_info[idx][2]
            flexres_mols[idx] = RDKitMolCreate.add_hydrogens(
                res, flexres_saved_coords[idx], flexres_hparents
            )

        return mol, flexres_mols, properties, flexres_residues

    def make_receptor_flexres_mols(
        self,
    ) -> tuple[list, list, list, list]:
        """
        Makes rdkit.Chem.Mols for the receptor based on flexible residues

        Raises:
            OutputError: _description_

        Returns:
            tuple[list, list, list, list]: _description_
        """

        mols = []
        info = []
        saved_coords = []
        residues = []

        with self.storageman as sm:
            flexible_residues, flexres_atomnames = sm.fetch_flexres_info(1)
        if flexible_residues is None:
            flexible_residues, flexres_atomnames = [], []
        elif flexible_residues != []:  # converts string to list
            flexible_residues = json.loads(flexible_residues)
            flexres_atomnames = json.loads(flexres_atomnames)
        # make flexible residue molecules
        for res, res_ats in zip(flexible_residues, flexres_atomnames):
            saved_coords.append([])
            resname = res[:3]
            res_ats = [
                at.strip() for at in res_ats
            ]  # strip out whitespace around atom names
            (
                res_smiles,
                res_index_map,
                res_h_parents,
            ) = RDKitMolCreate.guess_flexres_smiles(resname, res_ats)
            if res_smiles is None:  # catch error in guessing smiles
                raise OutputError(
                    f"Error while creating Mol for flexible residue {res}: unrecognized residue or incorrect atomtypes"
                )
            frm = Chem.MolFromSmiles(res_smiles)
            frm.SetProp("resinfo", res)
            mols.append(frm)
            info.append((res_smiles, res_index_map, res_h_parents))
            residues.append(res)

            return mols, info, saved_coords, residues

    def merge_databases(self, merging_db: str, backup: bool = True):
        """
        Method that merges two databases, ensuring integrity of primary and foreign keys.
        The merging will create a new table if needed, that keeps track of the primary key
        in the original and the merged database on a per-table basis. Another table will also
        keep track of how many databases has been merged into the primary database.
        The merging will ensure the two databases are -compatible based on the receptor only-.
        PLEASE NOTE: If two databases has been docked with dlg and vina respectively,
            these will be allowed to merge.

        Args:
            merging_db (str): path to database being merged into current db
            backup (bool, optional): whether or not to back up the current db before merging. Defaults to True.
        """

        LOGGER.warning(
            "If you have performed clustering, this data will be lost in the new, merged database. "
        )
        # create merge tables
        with self.storageman as sm:
            sm.merge_databases(**{"merging_db": merging_db, "backup": backup})

    def update_database_version(
        self, consent: bool = False, new_version: str = "3.0.0", backup: bool = False
    ):
        """
        Updates/upgrades sqlite database schema 1.0.0 through 3.0.0. The upgrade may take a while to
        complete, as it has to upgrade via each major upgrade, e.g., it will not upgrade straight
        from 1.0.0 to 3.0.0, but rather 1.0.0 -> 1.1.0 -> 2.0.0 -> 3.0.0 etc.

        Please #NOTE that currently no backup is made of the database, so please do this prior to
            upgrading.


        Args:
            consent (bool, optional): You have to give explicit consent to perform upgrade. Defaults to False.
            new_version (str, optional): What version you want to upgrade to. Defaults to "3.0.0".

        Returns:
            bool: consent/success
        """

        return self.storageman.update_database_version(new_version, consent, backup)

    @staticmethod
    def defaults(mode: str = None) -> dict:
        """
        Creates a dict of all Ringtail options.

        Return:
            str: json string with options
        """
        defaults = {}
        ringtail_defs = default_dict(RingtailDefaults())
        if mode == "cli":
            for key, value in ringtail_defs.items():
                if isinstance(value, list):
                    ringtail_defs[key] = ",".join(map(str, ringtail_defs[key]))
        defaults.update(default_dict(RingtailDefaults()))
        defaults.update(Filters().asdict())

        return defaults

    @staticmethod
    def generate_config_file_template():
        """Outputs to "config.json in current working directory

        Returns:
            str: file name of config file with template including default values
        """

        filename = "config.json"
        with open(filename, "w") as f:
            f.write(json.dumps(RingtailCore.defaults(), indent=4))
        return filename

    # endregion

    # region private methods

    def _num_poses_to_store(self, max_poses: int, store_all_poses: bool) -> int:
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

    def _fetch_select_ligands_poses(
        self,
        ligand_names: Union[str, list] = None,
        pose_ids: Union[int, list] = None,
        bookmark_name: str = None,
    ) -> dict[str, list[int]]:
        """
        Basically give this method anything, and get back a dict with ligand names and
        the pose ids selected by all the inputs. Can be a combo of ligand names, pose ids,
        and bookmark/status table. PS The combo only pose_ids and bookmark will work as if just
        the pose_ids were given.

        Args:
            ligand_names (Union[str, list], optional): _description_. Defaults to None.
            pose_ids (Union[int, list], optional): _description_. Defaults to None.
            bookmark_name (str, optional): _description_. Defaults to None.

        Returns:
            dict[str, list[int]]: _description_
        """
        select_ligand_poses = defaultdict(list)

        # if just given bookmark, get all ligands in that bookmark
        if bookmark_name:
            if not self._validate_bookmark_name(bookmark_name):
                raise RTCoreError(
                    f"Requested data, -{bookmark_name}-, is not a valid bookmark."
                )
            if not ligand_names and not pose_ids:
                with self.storageman as sm:
                    select_ligand_poses = sm.fetch_lignames_and_poses_for_selection(
                        bookmark_name
                    )

        # if just pose_id(s), get their ligand names and add to dict
        # assumes that the specific IDs were given for a reason, bookmark will not be considered
        if pose_ids:
            # ensure list
            if type(pose_ids) == int:
                pose_ids = [pose_ids]
            with self.storageman as sm:
                # get ligand name and add to dict
                for pose_id in pose_ids:
                    select_ligand_poses[sm.get_ligname_from_pose(pose_id)].append(
                        pose_id
                    )

        if ligand_names:
            # ensure list
            if type(ligand_names) == str:
                ligand_names = [ligand_names]
            for ligand_name in ligand_names:
                # get all poses for that ligand, possibly bracketed by bookmark or status table
                with self.storageman as sm:
                    poses = sm.fetch_selected_ligand_poses(ligand_name, bookmark_name)
                    select_ligand_poses[ligand_name].extend(poses)
        return select_ligand_poses

    def _update_mols_with_poses(
        self,
        pose_ids,
        atom_indices,
        mol,
        flexres_mols,
        flexres_info,
        ligand_saved_coords,
        flexres_saved_coords,
        properties,
    ) -> tuple:
        """Add requested poses and their information to a given ligand Mol

        Args:
            pose_ids (list[int]): list of pose ids to be evaluated
            atom_indices (list): List of ints indicating mapping of coordinate indices to smiles indices
            mol (RDKit.Chem.Mol): RDKit molecule for ligand
            flexres_mols (list): list of rdkit molecules for flexible residues
            flexres_info (list): list of tuples containing info for each flexible residue (res_smiles, res_index_map, res_h_parents)
            ligand_saved_coords (list): list of coordinates to save for adding hydrogens later
            flexres_saved_coords (list): list of lists of flexres coords to save for adding hydrogens later
            properties (dict): Dictionary of lists of properties, with each element corresponding to that conformer in the rdkit mol

            Returns:
                tuple: mol, flexres_mols, ligand_saved_coords, flexres_saved_coords, properties
        """
        with self.storageman as sm:
            for (
                Pose_ID,
                docking_score,
                leff,
                ligand_pose,
                flexres_pose,
            ) in sm.fetch_rdkit_relevant_pose_properties(pose_ids):
                # fetch info about pose interactions and format into string with format <type>-<chain>:<resname>:<resnum>:<atomname>:<atomnumber>, joined by commas
                interactions = self.storageman.fetch_pose_interactions(Pose_ID)
                # if that pose id has interactions
                if interactions is not None:
                    # make a list of all of them
                    interactions_list = []
                    # for each interaction row, make into a string according to format above
                    for interaction_info in interactions:
                        interaction = (
                            interaction_info[0] + "-" + ":".join(interaction_info[1:])
                        )
                        interactions_list.append(interaction)

                    interactions_str = ", ".join(interactions_list)
                    properties["Interactions"].append(interactions_str)

                # add properties to dictionary lists
                properties["Binding energies"].append(docking_score)
                properties["Ligand effiencies"].append(leff)
                # get pose coordinate info
                ligand_pose = json.loads(ligand_pose)
                flexres_pose = json.loads(flexres_pose)
                mol = RDKitMolCreate.add_pose_to_mol(mol, ligand_pose, atom_indices)
                for fr_idx, fr_mol in enumerate(flexres_mols):
                    flexres_mols[fr_idx] = RDKitMolCreate.add_pose_to_mol(
                        fr_mol, flexres_pose[fr_idx], flexres_info[fr_idx][1]
                    )
                    flexres_saved_coords[fr_idx].append(flexres_pose[fr_idx])
                ligand_saved_coords.append(ligand_pose)
        return mol, flexres_mols, ligand_saved_coords, flexres_saved_coords, properties

    def _process_receptor_polymer(self, receptor_file: str) -> tuple[str, str]:
        """
        Returns the json string of the contents of a receptor polymer file.
        Will try even if extension is not .json

        Args:
            receptor_file (str): (preferably) .json file with a serialized meeko.Polymer object

        Raises:
            OptionError

        Returns:
            str: string representation of the meeko.Polymer
        """

        if ".json" in receptor_file:
            receptor_name = Path(receptor_file).stem
            polypath = Path(receptor_file)
            if polypath.is_file():
                with open(polypath) as f:
                    polymer_string = f.read()
            else:
                raise OptionError(
                    "It appears a path to a receptor polymer json file was provided, but the file does not exist: ",
                    receptor_file,
                )
        else:
            try:
                polymer_string = json.loads(receptor_file)
            except json.JSONDecodeError:
                raise OptionError(
                    "The provided input for receptor_polymer was provided is not a valid Polymer object, not a valid path or valid json, cannot proceed."
                )
        return receptor_name, polymer_string

    def _process_receptor_pdbqt(self, receptor_file: str) -> tuple[str, bytes]:
        """
        Makes dict of rec name and receptor blob from pdbqt file (can be zipped)

        Args:
            receptor_file (str): pdbqt file

        Returns:
            tuple[str, bytes]: receptor name and bytes/blob object
        """
        receptor_name, receptor_blob = ReceptorManager.make_receptor_blob(receptor_file)
        return receptor_name, receptor_blob

    def _add_results(
        self,
        results: ResultsObject,
        duplicate_handling: str,
        overwrite: bool,
        docking_mode: str,
        num_poses: int,
        no_interactions: bool,
        interaction_tolerance: float,
        interaction_cutoffs: list,
        max_proc: int,
        finalize: bool,
    ):
        """Method that is agnostic of results type, and will do the actual call to storage manager to process result files and add to database.

        Args:
            results (ResultsObject): type checked and validated results object
            duplicate_handling (str): specify how duplicate Results rows should be handled when inserting into database. Options are "ignore" or "replace". Default behavior will allow duplicate entries.
            overwrite (bool): whether or not to overwrite existing storage
            docking_mode (str): docking engine used for the docking, describes file type/results style
            num_poses (int): max how many poses to store for a ligand, -1 is all
            no_interactions (bool): use if not wanting Ringtail to calculate and store interactions, only in vina mode
            interaction_tolerance (float): longest ångström distance that is considered interaction?
            interaction_cutoffs (list): ångström distance cutoffs for x and y interaction
            max_proc (int): max number of computer processors to use for file reading
            finalize (bool): whether or not this is last db write, will add indices

        Raises:
            OptionError
        """
        # validate docking mode
        docking_mode = validate_docking_mode(docking_mode)
        # Docking mode compatibility check
        if docking_mode == "vina" and interaction_tolerance:
            LOGGER.warning(
                "Cannot use interaction_tolerance with Vina mode. Removing interaction_tolerance."
            )
            interaction_tolerance = None
        # make sure we don't try to calculate interactions for adgpu
        if docking_mode == "adgpu":
            no_interactions = True

        processing_options = {
            "num_poses": num_poses,
            "no_interactions": no_interactions,
            "interaction_tolerance": interaction_tolerance,
            "interaction_cutoffs": interaction_cutoffs,
            "max_proc": max_proc,
        }

        # Process results files and handle database versioning
        with self.storageman:
            if overwrite:
                self.storageman.overwrite_storage()

            self.storageman.check_storage_ready(self._run_mode, docking_mode, num_poses)
        if results.save_receptor:
            self.save_receptor(results.receptor_file_path)

        # Prepare the results manager with the provided docking results sources
        self.resultsman = ResultsManager(
            self.db_file,
            docking_mode,
            get_valid_storageclass(self.storagetype),
            duplicate_handling,
        )

        # need receptor file contents if adding interaction
        if not no_interactions:
            # grab receptor info from database, this assumes there is only one receptor in the database
            # try pdbqt first
            results.receptor_string = self.get_receptor_object()[1]
            if not results.receptor_string:
                # if fail, try json
                try:
                    results.receptor_string = self.get_receptor_object()[2]
                except:
                    raise ResultsProcessingError(
                        "Trying to calculate interactions, but cannot find the receptor in the database. Please ensure to include the receptor_file and save_receptor if the receptor has not already been added to the database."
                    )

        LOGGER.info("Adding results...")
        self.resultsman.process_docking_data(results, processing_options)
        with self.storageman:
            if finalize:
                self.storageman.finalize_database_write()

    def _generate_interaction_combinations(self, filters: dict, max_miss: int) -> list:
        """Recursive function to list of tuples of possible interaction filter combinations, excluding up to max_miss interactions per filtering round

        Args:
            filters (dict): filter object with all filters
            max_miss (int): Maximum number of interactions to be excluded

        Returns:
            list: of interaction combinations for that specific combo
        """

        all_interactions = []
        for _type in Filters.get_filter_keys("interaction"):
            interactions = filters[_type]
            for interact in interactions:
                all_interactions.append(_type + "-" + interact[0])
        # warn if max_miss greater than number of interactions
        if max_miss > len(all_interactions):
            LOGGER.warning(
                "Requested max_miss options greater than number of interaction filters given. Defaulting to max_miss = number interaction filters"
            )
            max_miss = len(all_interactions)

        # BASE CASE:
        if max_miss == 0:
            return [tuple(all_interactions)]
        else:
            combinations = list(
                itertools.combinations(
                    all_interactions, len(all_interactions) - max_miss
                )
            )
            return combinations + self._generate_interaction_combinations(
                filters, max_miss=max_miss - 1
            )

    def _prepare_interaction_combo_filters(
        self, filters: dict, interaction_combination
    ):
        """Takes desired interaction combination, formats Filter object to dict, removes interactions not in given interaction_combination

        Args:
            filters (dict): original filters to be used as basis for creating specific filter setup given an interaction combo
            interaction_combination (list): list of interactions to be included in this round of filtering

        Returns:
            dict: filters for storageman
        """
        for itype in Filters.get_filter_keys("interaction"):
            itype_interactions = filters[itype]
            for interaction in itype_interactions:
                if itype + "-" + interaction[0] not in interaction_combination:
                    filters[itype].remove(interaction)
        # we want a match for all interactions when enumerating combinations
        filters["max_miss"] = 0

        return filters

    def _validate_bookmark_name(self, bookmark_name: str) -> Union[str, None]:
        """
        Ensures a bookmark exist, and runs checks on it. Particularly useful for dealing with union and interaction combinations,
        as user may provide the base bookmark_name but the filtering auto-generates derived bookmark_names for the combos.

        Args:
            bookmark_name (str)

        Raises:
            StorageError

        Returns:
            str: bookmark name if valid
        """
        with self.storageman:
            if self.storageman.is_bookmark(bookmark_name):
                # check if has max_miss filter
                bookmark_filters = self.storageman.fetch_filters_from_bookmark(
                    bookmark_name
                )
                try:
                    max_miss_present = bool(
                        bookmark_filters["max_miss"] > 0
                        and not bookmark_filters["enumerate_interaction_combs"]
                        and not "_union" in bookmark_name
                    )
                except:
                    #  in case bookmark query string does not contain the phrase 'max_miss', carry on
                    pass
                else:
                    if max_miss_present:
                        LOGGER.warning(
                            "'max_miss' used in filtering, but the bookmark used for sdfs writing is not the union of the search"
                        )
            # if bookmark name is not in the database
            else:
                # does bookmark name + _union resolve the issue
                if self.storageman.is_bookmark(bookmark_name + "_union"):
                    bookmark_name = bookmark_name + "_union"
                    LOGGER.warning(
                        "Requested 'export_sdf_path' with 'max_miss' and 'enumerate_interaction_combs' used in the filtering process. Exported SDFs will be for union of interaction combinations."
                    )
                elif bookmark_name.lower() in statuses:
                    # that is valid, will be handled properly
                    pass
                # if not, raise error
                else:
                    LOGGER.error(
                        f"Filtering bookmark {bookmark_name} does not exist in database. "
                    )
                    return None

            if not self.storageman.get_passing_poses_count(bookmark_name) and not (
                bookmark_name.lower() in statuses
            ):
                raise StorageError(
                    "Given results bookmark exists but does not have any data. Cannot write passing molecule SDFs"
                )
        return bookmark_name

    def _parse_clustering(
        self, cluster_data: dict, bookmark_name: str
    ) -> tuple[str, int]:
        """
        Takes dictionary with one or more cluster constraints and performs clustering one after another

        Args:
            cluster_data (dict): containing cluster types to be performed and their cutoffs
            bookmark_name (str): bookmark over which to perform the clustering

        Returns:
            tuple[str,int]: bookmark name of representative poses as well as count of clusters
        """
        with self.storageman:
            count_poses = self.storageman.get_passing_poses_count(bookmark_name, False)

        LOGGER.info(f"Preparing to cluster {count_poses} passing poses.")
        if len(cluster_data) > 1:
            LOGGER.warning(
                "N.B.: If using both interaction and morgan fingerprint clustering, the morgan fingerprint clustering will be performed first."
            )
        if cluster_data.get("ifp"):
            bookmark_name = self.cluster(
                bookmark_name,
                "ifp",
                cluster_data.get("ifp"),
            )
        if cluster_data.get("mfp"):
            bookmark_name = self.cluster(
                bookmark_name,
                "mfp",
                cluster_data.get("mfp"),
            )
        with self.storageman:
            count_reps = self.storageman.get_passing_poses_count(bookmark_name, True)
        print(
            f"\nNumber of cluster representative ligands:",
            count_reps,
        )
        return bookmark_name, count_reps

    def _export_csv(self, requested_data: str, csv_name: str, table=False):
        """Get requested data from database, export as CSV

        Args:
            requested_data (str): Table name or SQL-formatted query
            csv_name (str): Name for exported CSV file
            table (bool): flag indicating is requested data is a table name
        """
        with self.storageman:
            df = self.storageman.to_dataframe(requested_data, table=table)
            df.to_csv(csv_name)

    # endregion
