#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail virtual screening manager
#
import itertools
import os
from typing import Callable, Union
from collections import defaultdict
import functools
from collections.abc import Iterable
import json
from pathlib import Path
from meeko import RDKitMolCreate
from rdkit import Chem

from .interactions import find_interactions, make_interaction_finder
from .util import *
from .logutils import get_logger

logger = get_logger(__name__)
from .ringtailoptions import (
    RingtailDefaults,
    Filters,
    ResultsObject,
    statuses,
    validate_docking_mode,
)
from .schema import ORDER_RESULT_SCHEMA, OUTFIELD_BY_TABLE
from .exceptions import (
    RTCoreError,
    OutputError,
    StorageError,
    ResultsProcessingError,
    OptionError,
    MergeError,
    RingtailError,
)
from .mpmanager import MPManager
from .receptormanager import ReceptorManager
from .outputmanager import OutputManager
from .storagemanager import StorageManager
from .storagemanager_duckdb import StorageManagerDuckDB, HAS_DUCK
from .storagemanager_sqlite import StorageManagerSQLite, HAS_SQLITE
from .parsers import process_docked_mol, VinaMoleculeSupplier

storage_types = {}
if HAS_SQLITE:
    storage_types.update({"sqlite": StorageManagerSQLite})
if HAS_DUCK:
    storage_types.update({"duckdb": StorageManagerDuckDB})


def _wrap_exceptions(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except RingtailError:
            raise
        except Exception as e:
            raise RTCoreError(
                f"{func.__name__} failed with unexpected error: {e}"
            ) from e

    return wrapper


@_wrap_exceptions
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
        raise OptionError(f"Given storage type {storage_type} is not implemented.")


class RingtailCore:
    """Core class for coordinating different actions on virtual screening
    including adding results to storage, filtering and clusteirng

    Attributes:
        db_file (str): name of database file being operated on
        storageman (StorageManager): Interface module with database
        storagetype (str): database type
        _run_mode (str): refers to whether ringtail is ran from the command line or through direct API use, where the former is more restrictive
    """

    def __init__(
        self,
        db_file: str = RingtailDefaults.output_db,
        storage_type: str = RingtailDefaults.storage_type,
        access_mode: str = "api",
    ):
        """Initialize ringtail core, and create a storageman object with the db file.
        Can set logger level here, otherwise change it by logger.setLevel("level")

        Args:
            db_file (str): Database file to initialize core with. Defaults to "output.db".
            storage_type (str, optional): Database setup to use for interacting with the database.
                        If db file exist, it will attempt to detect storage type. Defaults to "sqlite".
            access_mode (str, optional): through what process Ringtail is accessed. Defaults to "api".
        """

        # try to find existing storage type, will overwrite input
        if Path(db_file).is_file():
            storage_type = detect_db_type(db_file)

        # Check if storage type is implemented
        storageman = get_valid_storageclass(storage_type)

        # initialize remaining variables and storagemanager
        self.storagetype = storage_type
        self.db_file = db_file
        self.storageman: StorageManager = storageman(db_file)

        logger.info(
            f"[     New RingtailCore object initialized with database file {db_file}    ]"
        )

        self._run_mode = access_mode.lower()

    def build_gui_tables(self):
        with self.storageman as sm:
            sm.ensure_gui_tables()

    @_wrap_exceptions
    def get_previous_docking_mode(self) -> Union[None, str]:
        """
        Will return the last used docking mode/engine of the database, normalized to the standardized
        names for docking modes used in ringtail (e.g., see RingtailOptions)

        Returns:
            Union[None, str]: None or previous docking mode
        """
        with self.storageman as sm:
            return validate_docking_mode(sm.get_previous_docking_mode())

    # region write to database
    @_wrap_exceptions
    def add_results_from_files(
        self,
        docking_results: "str | list[str]" = RingtailDefaults.docking_results,
        recursive: bool = RingtailDefaults.recursive,
        receptor_file: str = None,
        save_receptor: bool = RingtailDefaults.save_receptor,
        duplicate_handling: str = RingtailDefaults.duplicate_handling,
        overwrite: bool = RingtailDefaults.overwrite,
        docking_mode: str = RingtailDefaults.docking_mode,
        store_all_poses: bool = RingtailDefaults.store_all_poses,
        max_poses: int = RingtailDefaults.max_poses,
        calculate_interactions: bool = RingtailDefaults.calculate_interactions,
        interaction_tolerance: float = RingtailDefaults.interaction_tolerance,
        interaction_cutoffs: list = RingtailDefaults.interaction_cutoffs,
        max_proc: int = RingtailDefaults.max_proc,
        finalize: bool = True,
    ):
        """
        Call storage manager to process result files and add to database. Creates or adds to an existing a database.
        Options can be provided as a dict or as individual options. If both are provided, individual options will overwrite those from the dictionary.

        Args:
            docking_results (str | list[str]): one or more docking result files, directories to scan,
                or text files listing result paths. Type is auto-detected per item.
            recursive (bool): used to recursively search directories in docking_results
            receptor_file (str): string containing the receptor .pdbqt
            save_receptor (bool): whether or not to store the full receptor details in the database (needed for some things)
            duplicate_handling (str, options): specify how duplicate Results rows should be handled when inserting into database. Options are "ignore" or "replace". Default behavior will allow duplicate entries.
            overwrite (bool): whether or not to overwrite existing storage
            docking_mode (str): what docking engine was used to perform the docking
            store_all_poses (bool): store all ligand poses, does it take precedence over max poses?
            max_poses (int): how many poses to save (ordered by soem score?)
            calculate_interactions (bool): use if not wanting Ringtail to calculate and store interactions, only in vina mode
            interaction_tolerance (float): longest ångström distance that is considered interaction?
            interaction_cutoffs (list): ångström distance cutoffs for x and y interaction
            max_proc (int): max number of computer processors to use for file reading
            finalize (bool, optional): whether to finalize database creation and write indices, or wait, for example when you write to the db in batches

        Raises:
            OptionError

        """
        results = ResultsObject()
        results.docking_results = docking_results
        results.recursive_path_traverse = recursive
        results.receptor_file_path = receptor_file
        results.save_receptor = save_receptor

        if not results.has_results and not results.save_receptor:
            raise OptionError(
                "At least one input option needs to be used: docking_results or save_receptor"
            )

        # validate docking mode
        docking_mode = validate_docking_mode(docking_mode)
        num_poses = self._num_poses_to_store(max_poses, store_all_poses)

        # Docking mode compatibility check
        if docking_mode == "vina" and interaction_tolerance:
            logger.warning(
                "Cannot use interaction_tolerance with Vina mode. Removing interaction_tolerance."
            )
            interaction_tolerance = None
        if num_poses == -1 and interaction_tolerance:
            logger.warning(
                "Cannot use 'interaction_tolerance' and 'store_all_poses'. Removing 'interaction_tolerance'."
            )
            interaction_tolerance = None

        # make sure we don't try to calculate interactions for adgpu
        if docking_mode == "adgpu":
            calculate_interactions = False

        # Process results files and handle database versioning
        with self.storageman:
            if overwrite:
                self.storageman.overwrite_storage()

            self.storageman.check_storage_ready(self._run_mode, docking_mode, num_poses)
        if results.save_receptor:
            self.save_receptor(results.receptor_file_path)

        # Prepare the results manager with the provided docking results sources
        mpmanager = MPManager(
            max_proc=max_proc, receptor_file_path=results.receptor_file_path
        )

        # need receptor file contents if adding interaction
        if calculate_interactions:
            # grab receptor info from database, this assumes there is only one receptor in the database
            if results.receptor_file_path:
                import pathlib

                # parse the receptor string for interaction calculation
                # make receptor string
                if ".pdbqt" in pathlib.Path(results.receptor_file_path).suffixes:
                    logger.debug(
                        "Building receptor string from .pdbqt file for interaction calculation"
                    )
                    results.receptor_string = ReceptorManager.receptor_str_from_file(
                        results.receptor_file_path
                    )
                elif ".json" in pathlib.Path(results.receptor_file_path).suffixes:
                    logger.debug(
                        "Building receptor string from polymer .json file for interaction calculation"
                    )
                    with open(results.receptor_file_path, "r") as f:
                        polymer_json = f.read()
                    results.receptor_string = ReceptorManager.polymer_json2pdbqt_str(
                        polymer_json
                    )
                else:
                    logger.warning(
                        f"Receptor file {results.receptor_file_path} has unrecognized extension — receptor string not set, interactions will not be calculated."
                    )
                if results.receptor_string:
                    logger.debug(
                        f"Receptor string set successfully (length={len(results.receptor_string)})"
                    )
                else:
                    logger.warning(
                        "Receptor string is empty after parsing receptor file — interactions will not be calculated."
                    )
            else:
                logger.debug(
                    "No receptor file path provided, attempting to load receptor string from database"
                )
                # try pdbqt first
                receptor_data = self.get_receptor_object()
                results.receptor_string = receptor_data[1]
                if not results.receptor_string:
                    # if fail, try json
                    try:
                        results.receptor_string = receptor_data[2]
                    except IndexError:
                        raise ResultsProcessingError(
                            "Trying to calculate interactions, but cannot find the receptor in the database. Please ensure to include the receptor_file and save_receptor if the receptor has not already been added to the database."
                        )
                if results.receptor_string:
                    logger.debug(
                        f"Receptor string loaded from database (length={len(results.receptor_string)})"
                    )
                else:
                    logger.warning(
                        "Could not load receptor string from database — interactions will not be calculated."
                    )

        logger.debug(
            f"These are the provided docking results for database storage: {results.docking_results}"
        )

        logger.info("Adding results...")
        mpmanager.process_results(
            results,
            docking_mode,
            num_poses,
            interaction_tolerance,
            calculate_interactions,
            interaction_cutoffs,
            self.storageman.db_file,
            self.storageman.__class__,
            5000,
            duplicate_handling,
        )
        if finalize:
            self.finalize_write()

    @_wrap_exceptions
    def add_results_from_vina_string(
        self,
        results: Union[dict, Iterable[dict]],
        duplicate_handling: str = RingtailDefaults.duplicate_handling,
        store_all_poses: bool = RingtailDefaults.store_all_poses,
        max_poses: int = RingtailDefaults.max_poses,
        calculate_interactions: bool = RingtailDefaults.calculate_interactions,
        interaction_cutoffs: list = RingtailDefaults.interaction_cutoffs,
        receptor_string: str = None,
        chunk_size: int = 5000,
        finalize: bool = False,
    ):
        """Add Vina docking results from a dict or iterable of dicts to the database.

        Args:
            results (dict | Iterable[dict]): ligand name → PDBQT string mapping(s)
            duplicate_handling (str): how to handle duplicate Results rows ("ignore" or "replace")
            store_all_poses (bool): store all poses regardless of max_poses
            max_poses (int): maximum number of poses to store per ligand
            calculate_interactions (bool): whether to calculate and store interactions
            interaction_cutoffs (list): [hb_cutoff, vdw_cutoff] in ångströms
            receptor_string (str): receptor PDBQT or polymer JSON string; fetched from DB if omitted
            chunk_size (int): number of ligands to buffer before flushing to DB
            finalize (bool): reserved for future use
        """
        if isinstance(results, dict):
            results = [{key: value} for key, value in results.items()]
        elif not isinstance(results, Iterable):
            results = [results]

        interaction_finder = self._build_interaction_finder(
            calculate_interactions, receptor_string, interaction_cutoffs
        )
        num_poses = self._num_poses_to_store(max_poses, store_all_poses)
        vina_parser = VinaMoleculeSupplier(num_poses)

        self._insert_non_file_results(
            "vina",
            results,
            lambda r: vina_parser.parse_vina_results(
                r,
                interaction_finder=interaction_finder,
                calculate_interactions=calculate_interactions,
                is_file=False,
            ),
            duplicate_handling,
            num_poses,
            chunk_size,
            finalize,
        )

    @_wrap_exceptions
    def add_mol(
        self,
        mols: Union[Iterable[Chem.Mol], Chem.Mol],
        docking_mode: str = "adng",
        calculate_interactions: bool = RingtailDefaults.calculate_interactions,
        interaction_cutoffs: list = RingtailDefaults.interaction_cutoffs,
        receptor_string: str = None,
        chunk_size: int = 5000,
        duplicate_handling: str = RingtailDefaults.duplicate_handling,
        finalize: bool = False,
    ):
        """Add RDKit Mol objects directly to the database.

        All poses in the provided mols are stored. To limit poses, filter before calling.

        Args:
            mols (Chem.Mol | Iterable[Chem.Mol]): docked molecule(s) with pose properties
            docking_mode (str): must be "adng"
            calculate_interactions (bool): whether to calculate and store interactions
            interaction_cutoffs (list): [hb_cutoff, vdw_cutoff] in ångströms
            receptor_string (str): receptor PDBQT or polymer JSON string; fetched from DB if omitted
            chunk_size (int): number of ligands to buffer before flushing to DB
            duplicate_handling (str): how to handle duplicate Results rows
            finalize (bool): reserved for future use
        """
        docking_mode = validate_docking_mode(docking_mode)
        if docking_mode != "adng":
            raise OptionError(
                f"Docking mode {docking_mode} is not currently valid for adding Mols directly."
            )
        if not isinstance(mols, Iterable):
            mols = [mols]

        interaction_finder = self._build_interaction_finder(
            calculate_interactions, receptor_string, interaction_cutoffs
        )

        self._insert_non_file_results(
            docking_mode,
            mols,
            lambda m: process_docked_mol(
                m,
                interaction_finder=interaction_finder,
                calculate_interactions=calculate_interactions,
            ),
            duplicate_handling,
            -1,
            chunk_size,
            finalize,
        )

    @_wrap_exceptions
    def finalize_write(self):
        """
        Finalize database write by creating indices
        """
        with self.storageman:
            self.storageman.finalize_database_write()

    @_wrap_exceptions
    def save_receptor(self, receptor: Union[str, dict]):
        """
        Add receptor to database. Accets .pdbqt file as well as a polymer .json or a polymer object

        Args:
            receptor_file (str): path to receptor file

        """
        with self.storageman:
            self.storageman.check_storage_ready(self._run_mode)

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
            except Exception as e:
                raise OptionError(
                    f"Unable to serialize provided receptor object, cannot proceed to save receptor: {str(e)}."
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
                logger.warning(
                    "The Receptors table already has a .pdbqt blob present, will not add another."
                )
                return
            # if trying to add polymer json and there is polymer json, return
            if receptor_jsons and db_jsons:
                logger.warning(
                    "The Receptors table already has a polymer json string present, will not add another."
                )
                return
        with self.storageman:
            # if you made it this far, it is OK to add whatever you brought
            if receptor_blob:
                self.storageman.insert_receptor_blob(receptor_blob, receptor_name)
                logger.info("Receptor pdbqt blob data was added to the database.")
            elif receptor_jsons:
                self.storageman.insert_receptor_polymer(receptor_jsons, receptor_name)
                logger.info("Receptor polymer json string was added to the database.")

    @_wrap_exceptions
    def add_interactions(
        self,
        hb_cutoff: float = RingtailDefaults.interaction_cutoffs[0],
        vdw_cutoff: float = RingtailDefaults.interaction_cutoffs[1],
        receptor_string: str = None,
        get_consent: Callable = None,
        chunk_size: int = 500,
    ):
        """
        ###NOTE SLOOOOOOW

        Args:
            hb_cutoff (float, optional): _description_. Defaults to RingtailDefaults.interaction_cutoffs[0].
            vdw_cutoff (float, optional): _description_. Defaults to RingtailDefaults.interaction_cutoffs[1].
            receptor_string (str, optional): _description_. Defaults to None.
            get_consent (Callable, optional): _description_. Defaults to None.
            chunk_size (int, optional): _description_. Defaults to 10.

        Raises:
            RTCoreError: _description_

        Returns:
            _type_: _description_
        """
        from .parsers import generate_interaction_tuples

        # make sure user knows risk of calculating interaction in db with existing interactions
        track_table_name = "recomputed_interactions"
        if self.table_length("Interactions") > 0:

            def _api_cmd_consent():
                consent = input(
                    "WARNING: Calculating interactions for a database with existing interactions will delete all existing interactions.\n Are you sure you wish to proceed? If so, type 'yes': "
                )
                return consent.strip().lower() == "yes"

            if get_consent is None:
                get_consent = _api_cmd_consent

            if not get_consent():
                logger.critical(
                    "Consent not given for deleting and re-calculating interactions, exiting."
                )
                return
            else:
                with self.storageman as sm:
                    # Check if track table name exist, then do not clear but use
                    if not track_table_name in sm.tables_in_db():
                        success = sm.clear_interaction_tables()
                        if not success:
                            raise RTCoreError(
                                "Trouble while clearing existing interaction tables."
                            )
                    else:
                        logger.info(
                            "A table tracking processed pose interaction exists.\nWill continue processing poses and not recompute those that have already been processed."
                        )

        # get receptor representation
        if not receptor_string:
            receptor_string = self.get_receptor_representation()

        # accounting
        db_commit_counter = 0
        interactions = []
        processed_poseids = []
        results_counts = []

        with self.storageman as sm:
            # create a "temporary" table that keeps track of what pose_ids have been processed

            sm.create_transaction_tracking_table(track_table_name)

            # process in batches
            for pose in sm._stream_query(
                f"""SELECT R.pose_id, R.pose_coordinates, L.rdmol 
                FROM Results AS R 
                JOIN Ligands AS L ON L.ligand_id=R.ligand_id 
                LEFT JOIN {track_table_name} tt ON R.pose_id = tt.pose_id
                    WHERE tt.pose_id IS NULL;""",
            ):
                pose_id, coordinates, rdbin = pose
                mol = Chem.Mol(rdbin)
                interaction_dicts, num_hb, num_int = find_interactions(
                    [
                        (
                            {"pose_id": pose_id},
                            json.loads(coordinates),
                        )
                    ],
                    mol,
                    receptor_string,
                    hb_cutoff,
                    vdw_cutoff,
                )
                interactions.extend(generate_interaction_tuples(interaction_dicts))
                results_counts.append(
                    {"pose_id": pose_id, "num_hb": num_hb[0], "num_int": num_int[0]}
                )
                processed_poseids.append((pose_id,))
                db_commit_counter += 1

                if db_commit_counter > chunk_size - 1:
                    sm.post_insert_interactions(
                        interactions,
                        results_counts,
                        processed_poseids,
                        track_table_name,
                    )
                    interactions = []
                    processed_poseids = []
                    results_counts = []
                    db_commit_counter = 0

            #  insert remaining data after for loop
            if interactions:
                sm.post_insert_interactions(
                    interactions,
                    results_counts,
                    processed_poseids,
                    track_table_name,
                )
            sm._delete_table(track_table_name)

    # endregion

    # region filter, export

    @_wrap_exceptions
    def db_summary_data(
        self, columns=["docking_score", "leff"], percentiles=[1, 10]
    ) -> tuple[dict, dict]:
        """
        Produces dictionary of docking score and ligand efficiency data for 1st and 10th percentile


        Args:
            columns (list(str)): data columns used to prepare summary
            percentiles (list(int)): cutoff percentiles for the summary

        Returns:
            tuple[dict,dict]: relevant summary data and input columns and percentiles
        """
        with self.storageman:
            return self.storageman.fetch_summary_data(columns, percentiles), {
                "columns": columns,
                "percentiles": percentiles,
            }

    @_wrap_exceptions
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
        output_log: str = None,
        outfields: str = None,
        order_results: str = None,
        output_bookmark: str = RingtailDefaults.bookmark_name,
        input_bookmark: str = None,
        return_iter: bool = False,
        ligand_name_file=None,
    ) -> Union[tuple[int, str], iter]:
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
                ligand_name_file (str); csv file containing a list of ligand names to be filtered for
            Ligand results options:
                enumerate_interaction_combs (bool): When used with `max_miss` > 0, will log ligands/poses passing each separate interaction filter combination as well as union of combinations. Can significantly increase runtime.
                output_all_poses (bool): By default, will output only top-scoring pose passing filters per ligand. This flag will cause each pose passing the filters to be logged.
                mfpt_cluster (float): Cluster filtered ligands by Tanimoto distance of Morgan fingerprints with Butina clustering and output ligand with lowest ligand efficiency from each cluster. Default clustering cutoff is 0.5. Useful for selecting chemically dissimilar ligands.
                interaction_cluster (float): Cluster filtered ligands by Tanimoto distance of interaction fingerprints with Butina clustering and output ligand with lowest ligand efficiency from each cluster. Default clustering cutoff is 0.5. Useful for enhancing selection of ligands with diverse interactions.
                output_log (str): by default, results are saved in `output_log.txt`; if this option is used, ligands and requested info passing the filters will be written to specified file
                outfields (str): Comma-separated columns to include in output, e.g.
                    ``"docking_score,leff,num_hb"``. Ligand name is always returned first.
                    Available fields are the keys of ``ringtail.schema.OUTFIELD_SCHEMA``.
                order_results (str): Single column name to sort results by. Takes one value.
                    Available fields are the keys of ``ringtail.schema.ORDER_RESULT_SCHEMA``.
                output_bookmark (str): name for resulting book mark file. Default value is 'passing_results', can only contain lower case letters, numbers, and underscore (_)
                input_bookmark (str): name of bookmark to perform filtering over
                return_iter (bool): return an iterable of all of the filtering results

        Returns:
            tuple[int, str]: number of ligands passing filter and final bookmark name (may change if eg filtering and clustering)
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
                "ligand_name_file": ligand_name_file,
            }
        )

        # Compatibility check with docking mode
        if react_any:
            with self.storageman as sm:
                if sm.get_previous_docking_mode() == "vina":
                    logger.warning(
                        "Cannot use reaction filters with Vina mode. Removing react_any filter."
                    )
                react_any = False
        output_bookmark = valid_bookmark_name(output_bookmark)
        if output_bookmark is None:
            raise OptionError(
                "Bookmark name contains illegal symbols, please only use letters, numbers, and underscore."
            )
        # make sure enumerate_interaction_combs always true if max_miss = 0, since we don't ever worry about the union in this case
        if max_miss == 0:
            enumerate_interaction_combs = True

        # guard against unsing percentile filter with all_poses
        if output_all_poses and not (score_percentile is None or le_percentile is None):
            logger.warning(
                "Cannot return all passing poses with percentile filter. Will only log best pose."
            )
            output_all_poses = False
        # check that all filter values are valid
        filters.checks()

        logger.info("Starting to filter results")
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
                    filters.asdict(), interaction_combs[0]
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

        # when clustering, filter writes to {name}_preclust to avoid collision with the final
        # cluster output, which takes the user-provided output_bookmark
        effective_bookmark = (
            f"{output_bookmark}_preclust" if clustering else output_bookmark
        )

        if write_one_bookmark:
            with self.storageman:
                num_passing_ligands = self.storageman.filter_results(
                    all_filters=filter_dict,
                    output_bookmark=effective_bookmark,
                    input_bookmark=input_bookmark,
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
                iterated_bookmark_name = output_bookmark + "_" + str(ic_idx)
                # process each interaction combination filter
                with self.storageman:
                    num_passing_ligands = self.storageman.filter_results(
                        temp_filters,
                        iterated_bookmark_name,
                        input_bookmark,
                    )
                if num_passing_ligands:
                    logger.info(
                        f"\nNumber passing ligands filter combination {str(combination)}:{num_passing_ligands}"
                    )
                    if output_log:
                        self.write_filter_output(
                            iterated_bookmark_name,
                            temp_filters,
                            num_passing_ligands,
                            output_log,
                            outfields,
                            output_all_poses,
                            order_results,
                        )
                else:
                    logger.warning(
                        f"WARNING: No ligands found passing filter combination {str(combination)}."
                    )
                    with self.storageman:
                        self.storageman.delete_bookmark(iterated_bookmark_name)

            # Process the union of the max miss combinations of interactions
            with self.storageman:
                num_passing_ligands, output_bookmark = (
                    self.storageman.get_maxmiss_union(
                        len(interaction_combs), output_bookmark, filters_copy
                    )
                )
            # union result is the pre-cluster input for the max_miss path
            effective_bookmark = output_bookmark
            # rename so can be reused in common method below
            filter_dict = filters_copy
            print_string = " in max_miss union"

        if num_passing_ligands:
            logger.info(
                f"\nNumber passing ligands{print_string}: {num_passing_ligands}"
            )
            if clustering:
                _, num_passing_ligands = self._parse_clustering(
                    clustering, output_bookmark, effective_bookmark
                )
            # use original filters for final output log write
            if output_log:
                self.write_filter_output(
                    output_bookmark,
                    filters.asdict(),
                    num_passing_ligands,
                    output_log,
                    outfields,
                    output_all_poses,
                    order_results,
                    clustering,
                )
        else:
            logger.warning(f"WARNING: No ligands found passing filters{print_string}.")

        if return_iter:
            with self.storageman:
                formatted_query = self.storageman.get_bookmark_selection(
                    output_bookmark, outfields, not output_all_poses, order_results
                )
                # have to make it into list of tuples since storageman uses row factory
                iter = [
                    tuple(row)
                    for row in self.storageman.db_query(formatted_query).fetchall()
                ]
            return iter
        else:
            return num_passing_ligands, output_bookmark

    @_wrap_exceptions
    def cluster(
        self,
        output_bookmark: str,
        type: str = "mfp",
        cutoff: float = 0.5,
        input_bookmark: Union[str, None] = None,
    ) -> tuple[str, int]:
        """
        Clusters data and saves representatives to output_bookmark.

        Args:
            output_bookmark (str): bookmark name for the cluster representatives
            type (str, optional): cluster algorithm. Defaults to "mfp".
            cutoff (float, optional): cutoff distance for each cluster. Defaults to 0.5.
            input_bookmark (str | None, optional): bookmark to cluster over; None clusters all results.

        Returns:
            tuple[str, int]: output_bookmark, number of cluster representatives
        """
        with self.storageman as sm:
            generated_name, num_clusters = sm.cluster_data(
                input_bookmark,
                type,
                cutoff,
            )
            sm.rename_bookmark(generated_name, output_bookmark)
        logger.info(
            f"Number of clusters: {num_clusters}.\nPassing poses saved to {output_bookmark}."
        )
        return output_bookmark, num_clusters

    @_wrap_exceptions
    def write_filter_output(
        self,
        bookmark_name: str,
        filters: dict,
        num_passing_ligands: int,
        output_log: str = RingtailDefaults.output_log,
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
            output_log (str): name of file to write to
            output_fields (Union[str, list]): output columns to include
            output_all_poses (bool): whether or not to output one or all passing poses
            order_results (str): if ordering the results by a specific column
            clustering (dict, optional): dict describing any clustering done. Defaults to {}.
            max_miss_combs (int, optional): how many max_miss combinations were used. Defaults to None.
            append_to_log (bool, optional): whether or not to append to existing logfile. Defaults to None.

        Raises:
            OptionError
        """
        bookmark_name = self._normalize_bookmark_name(bookmark_name)
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
        with OutputManager(output_log, append_to_log) as opm:
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

    @_wrap_exceptions
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
        bookmark_name = self._normalize_bookmark_name(bookmark_name)
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
        if length := len(ligands_poses) > 10:

            def _api_cmd_consent(length: int):
                consent = input(
                    f"WARNING: You are requesting to prepare pdbs for {length} ligands/flexible residue combinations.\n Are you sure you wish to proceed? If so, type 'yes': "
                )
                return consent.strip().lower() == "yes"

            if get_consent is None:
                get_consent = _api_cmd_consent

            if not get_consent(length):
                logger.critical(
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
        grouped = self.create_rdkit_mols_by_ligand(ligands_poses, flexres_data)
        lig_flex_mol = {}
        for ligand, ligand_mol, all_flexres_per_pose, _, flexible_residues in grouped:
            if filename:
                root, ext = os.path.splitext(filename)
                if not ext:
                    ext = ".pdb"
                path = f"{root}_{ligand}{ext}"
            else:
                receptor_name, _, _ = self.get_receptor_object()
                path = f"{receptor_name}_{ligand}.pdb"

            pdb_parts = []
            for pose_idx, flexres_mols in enumerate(all_flexres_per_pose, start=1):
                flexmoldict = {}
                for index, flexres in enumerate(flexible_residues):
                    residue = flexres.split(":")[1]
                    chain = residue[0]
                    resnum = residue[1:]
                    flexmoldict[f"{chain}:{resnum}"] = flexres_mols[index]
                pdb_str = pdb_updated_flexres_from_rdkit(polymer, flexmoldict)
                pdb_parts.append(f"MODEL{pose_idx:>9}\n{pdb_str}ENDMDL\n")

            with open(path, "w") as file:
                file.write("".join(pdb_parts))

            lig_flex_mol[ligand] = {
                "ligand_mol": ligand_mol,
                "flexmoldict": {
                    f"{flexres.split(':')[1][0]}:{flexres.split(':')[1][1:]}": all_flexres_per_pose[
                        0
                    ][
                        i
                    ]
                    for i, flexres in enumerate(flexible_residues)
                },
            }

        return lig_flex_mol

    @_wrap_exceptions
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
        bookmark_name = self._normalize_bookmark_name(bookmark_name)

        try:
            ligands_poses = self._fetch_select_ligands_poses(
                ligand_names=ligname, bookmark_name=bookmark_name
            )
            flexres_data = self.make_receptor_flexres_mols()

        except StorageError as e:
            logger.error(str(e))
            raise RTCoreError(
                f"An error occurred while fetching molecule data for SDFs: {e}"
            ) from e

        opm = OutputManager()
        db_file_name = os.path.splitext(os.path.basename(self.db_file))[0]
        if sdf_path is not None and sdf_path != "" and not os.path.isdir(sdf_path):
            os.makedirs(sdf_path)
            logger.info(
                "Specified directory for SDF files was created in current working directory."
            )

        if all_in_one:
            all_pose_ids = [pid for poses in ligands_poses.values() for pid in poses]
            mol_results = self.create_rdkit_mols(all_pose_ids, flexres_data)
            if not mol_results:
                logger.error(
                    f"Selected bookmark {bookmark_name} does not exist or does not have any data, cannot write molecule SDFs."
                )
                return
            sdf_file_name = f"{db_file_name}_{bookmark_name}.sdf"
            logger.info(f"Writing {len(mol_results)} poses to {sdf_file_name}")
            opm.write_out_mols(
                sdf_file_name,
                (
                    (mol, flexres_mols, properties)
                    for _, mol, _, flexres_mols, properties, _ in mol_results
                ),
                sdf_path,
            )
        else:
            grouped = self.create_rdkit_mols_by_ligand(ligands_poses, flexres_data)
            if not grouped:
                logger.error(
                    f"Selected bookmark {bookmark_name} does not exist or does not have any data, cannot write molecule SDFs."
                )
                return
            for mol_ligname, mol, all_flexres, properties, _ in grouped:
                sdf_file_name = f"{mol_ligname}.sdf"
                logger.info(
                    f"Writing {mol.GetNumConformers()} conformers to {sdf_file_name}"
                )
                opm.write_out_mol(sdf_file_name, mol, all_flexres, properties, sdf_path)

    @_wrap_exceptions
    def find_similar_ligands(
        self, query_ligname: str, output_log: str = "cluster_log.txt"
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
                with OutputManager(output_log) as opm:
                    opm.write_find_similar_header(query_ligname, cluster_name)
                    opm.write_bookmarkname_in_log(bookmark_name)
                    opm.write_filter_results_in_log(similar_ligands)
                    opm.log_num_passing_ligands(number_similar)
            logger.info(f"\nNumber of similar ligands: {number_similar}")
        return number_similar

    @_wrap_exceptions
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
        bookmark_name = self._normalize_bookmark_name(bookmark_name)
        with self.storageman as sm:
            sql = sm.prepare_column_export_query(columns, bookmark_name)

        self._export_csv(sql, csv_name, False)

    @_wrap_exceptions
    def export_table_as_csv(self, table: str, csv_name: str = "table_output.csv"):
        """
        Exports specified table (all of it) as a csv file

        Args:
            table (str): table to be exported
            csv_name (str, optional): name of csv file. Defaults to "table_output.csv".
        """
        self._export_csv(table, csv_name, True)

    @_wrap_exceptions
    def export_sql_as_csv(self, sql: str, csv_name: str = "sql_output.csv"):
        """
        Exports sql-formatted query as csv, assumes user has used appropriate sql terminology for their database engine.

        Args:
            sql (str): query as a string
            csv_name (str, optional): name of csv file. Defaults to "sql_output.csv".
        """
        self._export_csv(sql, csv_name, False)

    @_wrap_exceptions
    def export_bookmark_db(self, bookmark_name: str, db_filepath: str = None) -> str:
        """Export database containing data from bookmark

        Args:
            bookmark_name (str): name for bookmark_db

        Returns:
            str: name of the new, exported database
        """
        bookmark_name = self._normalize_bookmark_name(bookmark_name)
        if db_filepath is None:
            db_filepath = self.db_file.removesuffix(".db") + "_" + bookmark_name + ".db"

        logger.info("Exporting bookmark database")
        if os.path.exists(db_filepath):
            logger.warning(
                "Requested export DB name already exists. Please rename or remove existing database. New database not exported."
            )
            return
        with self.storageman:
            self.storageman.create_subset_database(bookmark_name, db_filepath)

        return db_filepath

    @_wrap_exceptions
    def export_receptor_pdbqt(self):
        """
        Export receptor in database to pdbqt
        """
        recname, recstr, recjson = self.get_receptor_object()
        if recstr is None:
            if recjson:
                recstr = ReceptorManager.polymer_json2pdbqt_str(recjson)
            else:
                raise OptionError(
                    f"No receptor data stored for {recname}. Export failed."
                )
        output_manager = OutputManager()
        output_manager.write_receptor_pdbqt(recname, recstr)

    @_wrap_exceptions
    def get_available_outfields(self) -> dict[str, str]:
        """Return all column names available as output fields with their descriptions.

        Returns:
            dict[str, str]: mapping of column name to human-readable description.
                No database connection is required.
        """
        return {
            col: desc
            for cols in self.get_exportable_columns().values()
            for col, desc in cols.items()
        }

    def get_previous_filter_data(
        self,
        bookmark_name: str,
        outfields: str = "Ligname,docking_score",
        order_results: str = "ligname",
        output_log: str = None,
    ):
        """Get data requested in outfields from the
        results bookmark of a previous filtering

        Args:
            bookmark_name (str): bookmark for which the filters were used
            outfields (str, optional): columns to write to the output. Defaults to "Ligname,docking_score".
            order_results (str, optional): how to order the data. Defaults to "ligname".
            output_log (str, optional): log file name to write to. Defaults to None.
        """
        bookmark_name = self._normalize_bookmark_name(bookmark_name)
        if output_log is None:
            return
        if bookmark_name is None:
            raise OptionError("A bookmark name has to be provided")
        with self.storageman:
            new_data = self.storageman.fetch_data_for_passing_results(
                bookmark_name, outfields, order_results
            )
            with OutputManager(output_log) as opm:
                opm.write_filter_results_in_log(new_data)

    @_wrap_exceptions
    def cross_reference_databases(
        self,
        wanted_dbs: list[tuple[str, Union[str, None]]] = None,
        unwanted_dbs: list[tuple[str, Union[str, None]]] = None,
        bookmark_prefix: str = "crossref",
        alternative_database_names: dict = None,
    ) -> tuple[int, dict, dict]:
        """
        Method to cross reference two or more databases. Will attach all other databases
        to the "current" database, which is always the first wanted database (will have the
        active ringtail object).

        Will create an intersect of ligand names in the wanted databases+bookmarks,
        and omit any ligands found in the union of the unwanted databases+bookmarks.

        A bookmark with specified prefix to the original bookmark name will be stored in each database,
        making the crossreferencing data easily available for later use.

        Args:
            wanted_dbs (list[tuple[str, Union[str, None]]], optional): _description_. Defaults to None.
            unwanted_dbs (list[tuple[str, Union[str, None]]], optional): _description_. Defaults to None.
            bookmark_prefix (str, optional): _description_. Defaults to "crossref".
            alternative_database_names (dict, optional):  {path: alt name}. Defaults to None.

        Returns:
            tuple[int, dict, list, dict]: Number of "passing" ligands, dict of database {database
            file path: new bookmark name from cross ref}, list of passing ligand names, dict of dbs and how they were compared
        """
        with self.storageman as sm:
            return sm.crossref_databases(
                wanted_dbs, unwanted_dbs, bookmark_prefix, alternative_database_names
            )

    @_wrap_exceptions
    def detach_database(self, name: str) -> None:
        """Detach a previously attached database by its alias, no-op if not attached.

        Args:
            name (str): alias of the database to detach
        """
        with self.storageman as sm:
            attached = [db[0] for db in sm._check_attached() if db]
            if name in attached:
                sm.detach_db(name)

    @_wrap_exceptions
    def sanitize_dbpath_to_alias(self, db_path: str) -> str:
        """
        Makes a db safe alias for a database based on its path (cannot start with a number, only underscore symbols)

        Args:
            db_path (str): _description_

        Returns:
            str: _description_
        """
        return db_alias_from_path(db_path)

    @_wrap_exceptions
    def get_receptor_representation(self) -> str:
        with self.storageman as sm:
            rec_data = sm.fetch_receptor_object()
            if receptor_blob := rec_data.get("receptor_object"):
                return ReceptorManager.blob2str(receptor_blob)
            elif polymer_json := rec_data.get("polymer"):
                return polymer_json

    @_wrap_exceptions
    def get_receptor_object(self, as_pdbqt_str_only: bool = False) -> tuple:
        """
        Gets the receptor object from the database

        Args:
            as_pdbqt_str_only (bool, optional): if json is in db, return it as pdbqt str instead of json string. Defaults to False.

        Returns:
            tuple[str,str, str]: receptor name and receptor blob and receptor serialized polymer as strings
        """
        with self.storageman as sm:
            rec_data = sm.fetch_receptor_object()
            polymer_json = rec_data.get("polymer")
            if as_pdbqt_str_only and polymer_json:
                polymer = ReceptorManager.polymer_json2pdbqt_str(polymer_json)
            else:
                polymer = polymer_json
            if rec_data is not None:
                return (
                    rec_data.get("RecName"),
                    ReceptorManager.blob2str(rec_data.get("receptor_object")),
                    polymer,
                )
            else:
                return None, None, None

    # endregion

    # region utilities

    @_wrap_exceptions
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

    @_wrap_exceptions
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
        bookmark_name = self._normalize_bookmark_name(bookmark_name)
        status = status.lower()
        str_to_int = {v: k for k, v in statuses.items() if v is not None}
        str_to_int.update({"none": 0, "0": 0})
        if status not in str_to_int:
            raise OptionError(
                f"The selected status '{status}' is not valid. "
                f"Choose from: {', '.join(str_to_int.keys())}"
            )
        status = str_to_int[status]

        ligands_poses = self._fetch_select_ligands_poses(
            ligands, pose_ids, bookmark_name
        )
        # parse to just poses, I don't need ligand names for this
        parsed_poses = [
            pose for pose_list in ligands_poses.values() for pose in pose_list
        ]
        self.update_pose_status(pose_id=parsed_poses, status=status)

    @_wrap_exceptions
    def update_pose_status(self, pose_id: Union[int, list[int]], status: int):
        """
        Updates the status of a given pose so that it is only ever in one status table

        Args:
            pose_id (Union[int, list[int]]): pose id for which a status is assigned, can be a single pose or a list
            status (int): new status (0=unassign, 1=accepted, 2=maybe, 3=rejected)

        Raises:
            OptionError
        """
        if status not in statuses.keys():
            raise OptionError(
                f"""The selected status {status} is not a valid option for Ringtail. Please choose between {','.join(statuses)}"""
            )

        with self.storageman as sm:
            if status == 1:
                sm.accept_pose(pose_id)
            elif status == 2:
                sm.maybe_pose(pose_id)
            elif status == 3:
                sm.reject_pose(pose_id)
            else:
                sm.remove_status(pose_id)

    @_wrap_exceptions
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
        table = self._normalize_bookmark_name(table)
        with self.storageman as sm:
            return sm.fetch_viewable_data_columns_from(table, length, starting_pose_id)

    @_wrap_exceptions
    def get_exportable_columns(self) -> dict[str, dict[str, str]]:
        """Return exportable columns grouped by table with descriptions.

        Returns:
            dict[str, dict[str, str]]: ``{table_name: {col_name: description}}``
                for Results, Ligands, and Interaction_indices.
                No database connection is required.
        """
        return {
            table: {name: col.description for name, col in cols.items()}
            for table, cols in OUTFIELD_BY_TABLE.items()
        }

    @_wrap_exceptions
    def set_pose_comment(self, pose_id: int, comment: str) -> None:
        """Write or clear a user comment for a pose."""
        with self.storageman as sm:
            sm.set_pose_comment(pose_id, comment)

    @_wrap_exceptions
    def get_pose_comment(self, pose_id: int) -> "str | None":
        """Return the comment for a pose, or None if none exists."""
        with self.storageman as sm:
            return sm.get_pose_comment(pose_id)

    @_wrap_exceptions
    def get_row_from_pose(self, table: str, pose_id: int) -> int:
        """
        Method to get table row for a a pose id in given table, used to support gui

        Args:
            table (str): table name
            pose_id (int):

        Returns:
            int: row in given table
        """
        table = self._normalize_bookmark_name(table)
        with self.storageman as sm:
            return sm.pose_row_in_table(table, pose_id)

    @_wrap_exceptions
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
        table = self._normalize_bookmark_name(table)
        with self.storageman as sm:
            return sm.fetch_viewable_data_columns_from(
                table, length, starting_row_id, reverse
            )

    @_wrap_exceptions
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
        table = self._normalize_bookmark_name(table)
        with self.storageman as sm:
            return sm.fetch_columns_from_table_as_dicts(
                table, columns, length, starting_rowid
            )

    @_wrap_exceptions
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

    @_wrap_exceptions
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

    @_wrap_exceptions
    def get_ligand_name(self, pose_id: int = None):
        return self.db_query(
            "SELECT ligname FROM Ligands WHERE ligand_id = (SELECT ligand_id FROM Results WHERE pose_id = ? LIMIT 1);",
            (pose_id,),
        )[0][0]

    @_wrap_exceptions
    def get_starting_rowid(self, table: str):
        """
        Returns the starting row id for a table. For a normal table like Results or status table, this will be
        1, but for a bookmark, this will be the rowid it start with in the Filtered_poses table

        Args:
            table (str): table or bookmark name

        Returns:
            int: rowid
        """
        table = self._normalize_bookmark_name(table)
        with self.storageman:
            return self.storageman.get_starting_rowid(table)

    @_wrap_exceptions
    def database_is_empty(self) -> bool:
        """
        Determines whether or not there is data in the database

        Returns:
            bool: True if db is empty
        """
        with self.storageman as sm:
            return sm.db_empty()

    @_wrap_exceptions
    def get_range_of_e_le(
        self, table: str = "Results"
    ) -> tuple[float, float, float, float]:
        """
        Get min and max of e/docking_score and ligand efficiency/le/leff/

        Args:
            table (str): table limit data, e.g., either Results or a bookmark name

        Returns:
            tuple: e_min, e_max, le_min, le_max
        """
        with self.storageman as sm:
            return sm.get_range_of_e_le(table)

    @_wrap_exceptions
    def get_pose_interactions(self, pose_id: int) -> list[tuple]:
        """
        Gets all interactions for a given pose as a list of tuples with the interaction specification

        Args:
            pose_id (int): _description_

        Returns:
            list'tuple: interaction specifications as interaction_type,rec_chain,rec_resname,rec_resid,rec_atom,rec_atomid
        """
        with self.storageman as sm:
            return sm.fetch_pose_interactions(pose_id)

    @_wrap_exceptions
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

    @_wrap_exceptions
    def get_histogram_data(
        self, column: str, num_bins: int = 50, bookmark_name: str = None
    ) -> tuple[list[float], list[int]]:
        """SQL-side histogram counts for a numeric column in a bookmark.

        Args:
            column (str): numeric column name ("docking_score" or "leff")
            num_bins (int): number of bins. Defaults to 50.
            bookmark_name (str, optional): bookmark or table name. Defaults to None (→ Results).

        Returns:
            tuple[list[float], list[int]]: (bin_edges, counts) where len(bin_edges) == len(counts) + 1
        """
        bookmark_name = self._normalize_bookmark_name(bookmark_name)
        with self.storageman as sm:
            return sm.fetch_histogram_counts(column, num_bins, bookmark_name)

    @_wrap_exceptions
    def all_database_tables(self) -> list:
        """
        Returns a list of all table names in the database

        Returns:
            list: list of table names
        """
        with self.storageman as sm:
            return sm.tables_in_db()

    @_wrap_exceptions
    def delete_bookmark(self, bookmark_name: str):
        """Drops specified bookmark from the database

        Args:
            bookmark_name (str): name of bookmark to be dropped.
        """
        bookmark_name = self._normalize_bookmark_name(bookmark_name)
        with self.storageman:
            self.storageman.delete_bookmark(bookmark_name)
        logger.info(
            f"Bookmark {bookmark_name} was dropped from the database {self.db_file}"
        )

    @_wrap_exceptions
    def get_bookmark_names(self, db_name: str = None, db_path: str = None) -> list:
        """
        Method to retrieve all bookmark names in a database

        Args:
            db_name (str, optional): if attaching a database. Defaults to None.
            db_path (str, optional): path to attached database. Defaults to None.

        Returns:
            list: of all bookmarks in a database
        """
        with self.storageman:
            return self.storageman.get_all_bookmark_names(db_name, db_path)

    @_wrap_exceptions
    def get_bookmark_interactions(self, bookmark_name: str):
        """
        Get all interactions represented by the poses in a bookmark
        (could also be status table)

        Args:
            bookmark_name (str): _description_

        Returns:
            dict: key is column name,
        """
        bookmark_name = self._normalize_bookmark_name(bookmark_name)
        with self.storageman as sm:
            return sm.fetch_bookmark_interactions(bookmark_name)

    @_wrap_exceptions
    def create_rdkit_mols(
        self,
        pose_ids,
        flexres_data: tuple = None,
        include_interactions: bool = True,
    ) -> list:
        """Create one RDKit mol per pose_id. 1 DB call (no interactions) or 2 (with interactions).

        Args:
            pose_ids (int | list[int]): one or more pose ids
            flexres_data (tuple, optional): output of make_receptor_flexres_mols(). Defaults to None.
            include_interactions (bool): whether to fetch and embed interaction strings. Defaults to True.

        Returns:
            list of (pose_id, mol, ligname, flexres_mols, properties, flexres_residues)
        """
        if isinstance(pose_ids, int):
            pose_ids = [pose_ids]

        if flexres_data:
            flexres_mols_tmpl, flexres_info, _, flexres_residues = flexres_data
        else:
            flexres_mols_tmpl, flexres_info, flexres_residues = [], [], []

        with self.storageman as sm:
            pose_rows = sm.fetch_rdkit_pose_properties(pose_ids)
            if include_interactions:
                interactions_by_pose = (
                    sm.fetch_pose_interactions(pose_ids) if include_interactions else {}
                )
            else:
                interactions_by_pose = []

        result = []
        for (
            pose_id,
            rdmol_binary,
            ligname,
            docking_score,
            leff,
            pose_coords_json,
            flexres_coords_json,
            pose_rank,
        ) in pose_rows:
            mol = Chem.Mol(rdmol_binary)
            pose_coordinates = json.loads(pose_coords_json)
            if mol.GetNumAtoms() != len(pose_coordinates):
                logger.error(
                    f"{mol.GetNumAtoms()=} differs from {len(pose_coordinates)=}"
                )
                continue

            mol = self._add_conformer(mol, pose_coordinates)
            if any(a.GetNumImplicitHs() > 0 for a in mol.GetAtoms()):
                mol = Chem.AddHs(mol, addCoords=True)
            mol.SetProp("_Name", ligname)

            # Build properties dict (single-element lists matching old format)
            if include_interactions:
                interactions = interactions_by_pose.get(pose_id, [])
                interactions_str = ", ".join(
                    f"{row[0]}-{':'.join(str(x) for x in row[1:])}"
                    for row in interactions
                )
            else:
                interactions_str = None
            properties = {
                "Binding energies": [docking_score],
                "Ligand effiencies": [leff],
                "Interactions": [interactions_str] if interactions_str else [],
                "pose_rank": pose_rank,
            }
            if not mol.HasProp("_Name"):
                properties["_Name"] = ligname

            # Flexres: fresh copy of templates for each pose
            pose_flexres_mols = []
            if flexres_mols_tmpl and flexres_coords_json:
                flexres_pose_coords = json.loads(flexres_coords_json)
                for fr_idx, fr_mol in enumerate(flexres_mols_tmpl):
                    fr_copy = Chem.Mol(fr_mol)
                    fr_copy = RDKitMolCreate.add_pose_to_mol(
                        fr_copy,
                        flexres_pose_coords[fr_idx],
                        flexres_info[fr_idx][1],
                    )
                    fr_copy = RDKitMolCreate.add_hydrogens(
                        fr_copy,
                        [flexres_pose_coords[fr_idx]],
                        flexres_info[fr_idx][2],
                        True,
                    )
                    pose_flexres_mols.append(fr_copy)

            result.append(
                (pose_id, mol, ligname, pose_flexres_mols, properties, flexres_residues)
            )

        return result

    def create_rdkit_mols_by_ligand(
        self,
        ligands_poses: dict,
        flexres_data: tuple = None,
        include_interactions: bool = True,
    ):
        """One multi-conformer mol per ligand; all poses become conformers.

        Parses rdmol_binary once per ligand; subsequent poses reuse an in-memory
        topology copy instead of re-deserializing the binary. Yields one ligand
        at a time so callers can process/display progressively.

        Args:
            ligands_poses (dict): {ligname: [pose_id, ...]} from _fetch_select_ligands_poses
            flexres_data (tuple, optional): output of make_receptor_flexres_mols()
            include_interactions (bool): whether to embed interaction strings

        Yields:
            (ligname, mol_N_conformers, per_pose_flexres_mols, merged_properties, flexres_residues)
            where per_pose_flexres_mols is a list of flexres_mol lists, one per pose
        """
        all_pose_ids = [pid for poses in ligands_poses.values() for pid in poses]

        with self.storageman as sm:
            pose_rows = sm.fetch_rdkit_pose_properties(all_pose_ids)
            interactions_by_pose = (
                sm.fetch_pose_interactions(all_pose_ids) if include_interactions else {}
            )

        # Group raw DB rows by ligand before any RDKit work
        by_ligand = {}
        for row in pose_rows:
            by_ligand.setdefault(row[2], []).append(row)  # row[2] = ligname

        if flexres_data:
            flexres_mols_tmpl, flexres_info, _, flexres_residues = flexres_data
        else:
            flexres_mols_tmpl, flexres_info, flexres_residues = [], [], []

        for ligname, rows in by_ligand.items():
            try:
                # Deserialize binary ONCE per ligand → topology skeleton (no H, no conformers)
                topology = Chem.Mol(rows[0][1])
            except Exception:
                logger.error(f"Failed to deserialize topology for {ligname}, skipping")
                continue

            base_mol = None
            merged_props = {
                "Binding energies": [],
                "Ligand effiencies": [],
                "Interactions": [],
            }
            all_flexres = []

            for pose_id, _, _, docking_score, leff, coords_json, fr_coords_json in rows:
                try:
                    coords = json.loads(coords_json)

                    # In-memory copy of topology (no binary parse for poses 2..N)
                    tmp = Chem.Mol(topology)
                    if tmp.GetNumAtoms() != len(coords):
                        logger.error(
                            f"atom/coord mismatch for pose {pose_id}, skipping"
                        )
                        continue
                    tmp = self._add_conformer(tmp, coords)
                    tmp = Chem.AddHs(tmp, addCoords=True)
                except Exception:
                    logger.error(f"Failed to process pose {pose_id}, skipping")
                    continue

                if base_mol is None:
                    base_mol = tmp
                    base_mol.SetProp("_Name", ligname)
                else:
                    base_mol.AddConformer(tmp.GetConformer(), assignId=True)

                interactions = interactions_by_pose.get(pose_id, [])
                interactions_str = ", ".join(
                    f"{r[0]}-{':'.join(str(x) for x in r[1:])}" for r in interactions
                )
                merged_props["Binding energies"].append(docking_score)
                merged_props["Ligand effiencies"].append(leff)
                if interactions_str:
                    merged_props["Interactions"].append(interactions_str)

                pose_flexres = []
                if flexres_mols_tmpl and fr_coords_json:
                    try:
                        flexres_pose_coords = json.loads(fr_coords_json)
                        for fr_idx, fr_mol in enumerate(flexres_mols_tmpl):
                            fr_copy = Chem.Mol(fr_mol)
                            fr_copy = RDKitMolCreate.add_pose_to_mol(
                                fr_copy,
                                flexres_pose_coords[fr_idx],
                                flexres_info[fr_idx][1],
                            )
                            fr_copy = RDKitMolCreate.add_hydrogens(
                                fr_copy,
                                [flexres_pose_coords[fr_idx]],
                                flexres_info[fr_idx][2],
                                True,
                            )
                            pose_flexres.append(fr_copy)
                    except Exception:
                        logger.error(
                            f"Failed to process flexres for pose {pose_id}, skipping flexres"
                        )
                        pose_flexres = []
                all_flexres.append(pose_flexres)

            if base_mol is not None:
                yield (ligname, base_mol, all_flexres, merged_props, flexres_residues)

    @_wrap_exceptions
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

    @_wrap_exceptions
    def merge_databases(self, merging_dbs: list[str], backup: bool = True):
        """
        Merges one or more databases into the primary database. Accepts a list
        of explicit paths, glob patterns (e.g. "path/to/*.db"), or a mix of both.
        Duplicate paths and the primary database itself are automatically excluded.

        Args:
            merging_dbs (list[str]): paths or glob patterns for databases to merge
            backup (bool, optional): back up current db before merging. Defaults to True.

        Returns:
            list[tuple[str, str]]: list of (db_path, error_message) for failed merges
        """
        import glob

        if isinstance(merging_dbs, str):
            merging_dbs = [merging_dbs]

        # Expand glob patterns and collect unique paths
        primary_abspath = os.path.abspath(self.db_file)
        seen = set()
        expanded_dbs = []
        for pattern in merging_dbs:
            matches = glob.glob(pattern)
            if not matches:
                logger.warning(f"No files matched pattern: {pattern}")
            for db in matches:
                abspath = os.path.abspath(db)
                if abspath == primary_abspath:
                    logger.info(f"Skipping {db} (same as primary database)")
                    continue
                if abspath not in seen:
                    seen.add(abspath)
                    expanded_dbs.append(db)

        if not expanded_dbs:
            raise OptionError("No valid databases to merge after filtering.")

        logger.warning(
            "If you have performed filtering or clustering, this data will be lost in the new, merged database."
        )
        logger.info(f"Merging {len(expanded_dbs)} databases into {self.db_file}")

        with self.storageman as sm:
            if backup:
                sm.clone()
            sm.prepare_for_merging()
            failed = []
            for merging_db in expanded_dbs:
                try:
                    sm.merge_database(merging_db=merging_db)
                    logger.info(f"Successfully merged {merging_db}.")
                except MergeError as e:
                    logger.error(f"Database {merging_db} failed to merge: {e}")
                    failed.append((merging_db, str(e)))

            sm.complete_merging()

        return failed

    @_wrap_exceptions
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

    @_wrap_exceptions
    def db_compatibility_check(self, database_path: str) -> bool:
        """
        Checks if the database in the provided path is compatible for operations with the
        current Ringtail database.

        Args:
            database_path (str): _description_

        Returns:
            bool: whether or not compatible
        """
        return detect_db_type(self.db_file) == detect_db_type(database_path)

    @staticmethod
    def defaults() -> dict:
        """
        Creates a dict of all Ringtail options.

        Return:
            str: json string with options
        """
        from .ringtailoptions import ringtail_defaults

        return ringtail_defaults()

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

    def _insert_non_file_results(
        self,
        docking_mode: str,
        items: Iterable,
        parse_fn: Callable,
        duplicate_handling: str,
        num_poses: int,
        chunk_size: int,
        finalize: bool,
    ):
        """
        Method to process non-file results

        Args:
            docking_mode (str): _description_
            items (Iterable): _description_
            parse_fn (Callable): _description_
            duplicate_handling (str): _description_
            num_poses (int): _description_
            chunk_size (int): _description_
            finalize (bool): _description_
        """
        with self.storageman:
            self.storageman.check_storage_ready(self._run_mode, docking_mode, num_poses)

        buffer = {"ligands": [], "poses": [], "interactions": []}
        with self.storageman as sm:
            for item in items:
                result = parse_fn(item)
                for key in buffer:
                    buffer[key].extend(result.get(key, []))
                if len(buffer["ligands"]) >= chunk_size:
                    sm.insert_data(buffer, duplicate_handling)
                    buffer = {"ligands": [], "poses": [], "interactions": []}
            if buffer["ligands"]:
                sm.insert_data(buffer, duplicate_handling)

        if finalize:
            self.finalize_write()

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

    def _build_interaction_finder(
        self,
        calculate_interactions: bool,
        receptor_string: str,
        interaction_cutoffs: list[float, float],
    ):
        """
        Create an InteractionFinder object

        Args:
            calculate_interactions (bool): _description_
            receptor_string (str): _description_
            interaction_cutoffs (list[float,float]): hb and vdw max distance for interaction

        Returns:
            _type_: _description_
        """
        if not calculate_interactions:
            return None
        if not receptor_string:
            receptor_string = self.get_receptor_representation()
        return make_interaction_finder(receptor_string, interaction_cutoffs)

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

    def _add_conformer(self, mol: Chem.Mol, coordinates: list[list]) -> Chem.Mol:
        """
        Adds a conformer with given coordinates to a mol

        Args:
            mol (Chem.Mol): _description_
            coordinates (list[list]): _description_

        Returns:
            Chem.Mol: _description_
        """
        from rdkit import Geometry

        n_atoms = mol.GetNumAtoms()
        conf = Chem.Conformer(n_atoms)
        for i, (x, y, z) in enumerate(coordinates):
            conf.SetAtomPosition(i, Geometry.Point3D(float(x), float(y), float(z)))
        mol.AddConformer(conf, assignId=True)
        return mol

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
            logger.warning(
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
                except KeyError:
                    #  in case bookmark query string does not contain the phrase 'max_miss', carry on
                    pass
                else:
                    if max_miss_present:
                        logger.warning(
                            "'max_miss' used in filtering, but the bookmark used for sdfs writing is not the union of the search"
                        )
            # if bookmark name is not in the database
            else:
                # does bookmark name + _union resolve the issue
                if self.storageman.is_bookmark(bookmark_name + "_union"):
                    bookmark_name = bookmark_name + "_union"
                    logger.warning(
                        "Requested 'export_sdf_path' with 'max_miss' and 'enumerate_interaction_combs' used in the filtering process. Exported SDFs will be for union of interaction combinations."
                    )
                elif bookmark_name.lower() in statuses.values():
                    # that is valid, will be handled properly
                    pass
                # if not, raise error
                else:
                    raise OptionError(
                        f"Filtering bookmark {bookmark_name} does not exist in database. "
                    )

            if not self.storageman.get_passing_poses_count(bookmark_name) and not (
                bookmark_name.lower() in statuses
            ):
                raise StorageError(
                    "Given results bookmark exists but does not have any data. Cannot write passing molecule SDFs"
                )
        return bookmark_name

    def _parse_clustering(
        self,
        cluster_data: dict,
        output_bookmark: str,
        input_bookmark: Union[str, None] = None,
    ) -> tuple[str, int]:
        """
        Takes dictionary with one or more cluster constraints and performs clustering one after another.

        Args:
            cluster_data (dict): containing cluster types to be performed and their cutoffs
            output_bookmark (str): bookmark name for the final cluster representatives
            input_bookmark (str | None, optional): bookmark to cluster over; None clusters all results.

        Returns:
            tuple[str,int]: output_bookmark, count of cluster representatives
        """
        if input_bookmark is not None:
            with self.storageman:
                count_poses = self.storageman.get_passing_poses_count(
                    input_bookmark, False
                )
            logger.info(f"Preparing to cluster {count_poses} passing poses.")
        else:
            logger.info("Preparing to cluster all results.")

        if sum(1 for v in cluster_data.values() if v) > 1:
            logger.warning(
                "N.B.: If using both interaction and morgan fingerprint clustering, the morgan fingerprint clustering will be performed first."
            )
        count_reps = 0
        if cluster_data.get("mfp"):
            _, count_reps = self.cluster(
                output_bookmark,
                "mfp",
                cluster_data.get("mfp"),
                input_bookmark,
            )
        if cluster_data.get("ifp"):
            # if mfp was also run, chain: cluster the mfp output
            ifp_input = output_bookmark if cluster_data.get("mfp") else input_bookmark
            _, count_reps = self.cluster(
                output_bookmark,
                "ifp",
                cluster_data.get("ifp"),
                ifp_input,
            )
        logger.info(f"\nNumber of cluster representative ligands: {count_reps}")
        return output_bookmark, count_reps

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

    def _normalize_bookmark_name(self, bookmark_name: str) -> str:
        if bookmark_name is None:
            return None
        return bookmark_name.lower()

    # endregion
