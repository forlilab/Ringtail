#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail storage adaptors
#

import sqlite3
import time
import json
import pandas as pd
from .logutils import LOGGER as logger
import sys
from signal import signal, SIGINT
from rdkit import Chem
from rdkit import DataStructs
from rdkit.ML.Cluster import Butina
import numpy as np
import time
from importlib.metadata import version
from .ringtailoptions import Filters
from .util import numlist2str
from .exceptions import (
    StorageError,
    DatabaseInsertionError,
    DatabaseConnectionError,
    DatabaseTableCreationError,
)
from .exceptions import DatabaseQueryError, DatabaseViewCreationError, OptionError

import multiprocess


class StorageManager:

    _db_schema_ver = "2.0.0"

    # "db_schema_ver":list("compatible code versions")
    _db_schema_code_compatibility = {
        "1.0.0": ["1.0.0"],
        "1.1.0": ["1.1.0"],
        "2.0.0": ["2.0.0", "2.1.0"],
    }

    """Base class for a generic virtual screening database object.
    This class holds some of the common API for StorageManager child classes. 
    Each child class will implement their own functions to write to and read from the database

    Attributes: 
        _db_schema_ver (str): current database schema version
        _db_schema_code_compatibility (dict): dictionary showing compatibility of code base versions with relational database schema versions
    """

    # region database access
    def check_storage_compatibility(storage_type):
        """Checks if chosen storage type has been implemented

        Args:
            storage_type (str): name of the storage type

        Raises:
            NotImplementedError: raised if seelected storage type has not been implemented

        Returns:
            class: of implemented storage type
        """

        storage_types = {
            "sqlite": StorageManagerSQLite,
        }
        if storage_type in storage_types:
            return storage_types[storage_type]
        else:
            raise NotImplementedError(
                f"Given storage type {storage_type} is not implemented."
            )

    def __init__(self):
        """Initialize instance variables common to all StorageManager subclasses"""
        self.logger = logger
        self.closed_connection = False

    def __enter__(self):
        """Used to access the database if using storage manager as a context manager

        Raises:
            StorageError

        Returns:
            instance: of class with open database connection
        """
        try:
            self._open_storage()
        except StorageError as e:
            raise e
        else:
            return self

    def __exit__(self, exc_type, exc_value, tb):
        """Used to close the database if using storage manager as a context manager

        Args:
            exc_type (_type_): error exit parameter requirerd when using context manager
            exc_value (_type_): error exit parameter requirerd when using context manager
            tb (_type_): error exit parameter requirerd when using context manager

        Returns:
            instance: of class with closed database connection
        """
        if not self.closed_connection:
            self.close_storage()
        if exc_type:
            if exc_type == Exception:
                self.logger.error(str(exc_value))
            else:
                raise
        return self

    def _sigint_handler(self, signal_received, frame):
        """Handles and reports if program is interrupted through the terminal"""
        self.logger.critical("Ctrl + C pressed, keyboard interupt initiated")
        self.__exit__(None, None, None)
        sys.exit(0)

    def finalize_database_write(self):
        """
        Methods to finalize when a database has been written to, and saving the current database schema to the sqlite database.
        """
        # index certain tables
        self._create_indices()
        # set version of the database
        self._set_ringtail_db_schema_version(self._db_schema_ver)
        self.logger.info("Database write session completed successfully.")

    def close_storage(self, attached_db=None, vacuum=False):
        """Close connection to database

        Args:
            attached_db (str, optional): name of attached DB (not including file extension)
            vacuum (bool, optional): indicates that database should be vacuumed before closing
        """
        if attached_db is not None:
            self._detach_db(attached_db)
        # close any open cursors
        self._close_open_cursors()
        # vacuum database
        if vacuum:
            self._vacuum()
        # close db itself
        self._close_connection()
        self.closed_connection = True

    # endregion

    # region write data
    def insert_data(
        self,
        results_array,
        ligands_array,
        interaction_list,
        receptor_array=[],
        insert_receptor=False,
    ):
        """Inserts data from all arrays returned from results manager.

        Args:
            results_array (list): list of data to be stored in Results table
            ligands_array (list): list of data to be stored in Ligands table
            interaction_list (list): list of data to be stored in interaction tables
            receptor_array (list): list of data to be stored in Receptors table
            insert_receptor (bool, optional): flag indicating that receptor info should inserted
        """
        Pose_IDs, duplicates = self._insert_results(results_array)
        self._insert_ligands(ligands_array)
        if insert_receptor and receptor_array != []:
            # first checks if there is receptor info already in the db
            receptors = self.fetch_receptor_objects()
            # insert receptor if database does not have already have a receptor entry
            if len(receptors) == 0:
                self._insert_receptors(receptor_array)
        # insert interactions if they are present
        if interaction_list != []:
            self.insert_interactions(Pose_IDs, interaction_list, duplicates)

    def insert_interactions(self, Pose_IDs: list, interactions_list, duplicates):
        """Takes list of interactions, inserts into database

        Args:
            Pose_IDs (list(int)): list of pose ids assigned while writing the current results to database
            interactions_list (list): List of tuples for interactions in form
                ("type", "chain", "residue", "resid", "recname", "recid")
            duplicates (list(Pose_ID)): any duplicates identified in "insert_results", if duplicate handling has been specified

        """

        # for each pose id, list
        interaction_rows = []
        for index, Pose_ID in enumerate(Pose_IDs):
            # add interaction if unique, returns index of interaction
            pose_interactions = [
                ((Pose_ID,) + self._insert_interaction_index_row(interaction_tuple))
                for interaction_tuple in interactions_list[index]
            ]
            # adds each pose_interaction row to list
            interaction_rows.extend(pose_interactions)
        self._insert_interaction_rows(interaction_rows, duplicates)

    # endregion

    # region read data
    def filter_results(self, all_filters: dict, suppress_output=False) -> iter:
        """Generate and execute database queries from given filters.

        Args:
            all_filters (dict): dict containing all filters. Expects format and keys corresponding to ringtail.Filters().todict()
            suppress_output (bool): prints filtering summary to sdout

        Returns:
             iter: iterable, such as an sqlite cursor, of passing results
        """
        # before we do anything, check that the DB version matches the version number of our module
        rt_version_same, db_rt_version = self.check_ringtaildb_version()
        if not rt_version_same:
            # NOTE will cause error when any version int is > 10
            # catch version 1.0.0 where returned db_rt_version will be 0
            if db_rt_version == 0:
                db_rt_version = 100
            raise StorageError(
                f"Input database was created with Ringtail v{'.'.join([i for i in db_rt_version[:2]] + [db_rt_version[2:]])}. Confirm that this matches current Ringtail version and use Ringtail update script(s) to update database if needed."
            )
        # create view of passing results
        filter_results_str, view_query = self._generate_result_filtering_query(
            all_filters
        )
        self.logger.debug(f"Query for filtering results: {filter_results_str}")

        # if max_miss> and we are enumerating interaction combinations, we want to give each passing view a new name by changing the self.bookmark_name
        if self.view_suffix is not None:
            self.current_bookmark_name = self.bookmark_name + "_" + self.view_suffix
        else:
            self.current_bookmark_name = self.bookmark_name

        # make sure we keep Pose_ID in view
        self.create_bookmark(
            name=self.current_bookmark_name, query=view_query, filters=all_filters
        )

        # perform filtering
        if suppress_output:
            return None

        self.logger.debug("Running filtering query...")
        time0 = time.perf_counter()
        filtered_results = self._run_query(filter_results_str).fetchall()
        self.logger.debug(
            f"Time to run query: {time.perf_counter() - time0:.2f} seconds"
        )
        return filtered_results

    def check_passing_bookmark_exists(self, bookmark_name: str | None = None):
        """Checks if bookmark name is in database

        Args:
            bookmark_name (str, optional): name of bookmark name to check if exist, or else will use storageman bookmark_name attribute

        Returns:
            bool: indicates if bookmark_name exists in the current database
        """
        if bookmark_name is None:
            bookmark_name = self.bookmark_name

        bookmark_exists = bookmark_name in self.get_all_bookmark_names()

        return bookmark_exists

    def get_plot_data(self, bookmark_name: str = None, only_passing=False):
        """This function is expected to return an ascii plot
        representation of the results

        Args:
            bookmark_name (str): name of bookmark for which to fetch passing data. Will use default bookmark name if None. Returns empty list if bookmark does not exist.
            only_passing (bool): Only return data for passing ligands. Will return empty list for all data.

        Returns:
            tuple: cursors as (<all data cursor>, <passing data cursor>)
        """

        # checks if we have filtered by looking for view name in list of view names
        if self.check_passing_bookmark_exists(bookmark_name):
            if only_passing:
                return [], self._fetch_passing_plot_data(bookmark_name)
            else:
                return self._fetch_all_plot_data(), self._fetch_passing_plot_data(
                    bookmark_name
                )
        else:
            return self._fetch_all_plot_data(), []

    def crossref_filter(
        self,
        new_db: str,
        bookmark1_name: str,
        bookmark2_name: str,
        selection_type="-",
        old_db=None,
    ) -> tuple:
        """Selects ligands found or not found in the given bookmark in both current db and new_db. Stores as temp view

        Args:
            new_db (str): file name for database to attach
            bookmark1_name (str): string for name of first bookmark/temp table to compare
            bookmark2_name (str): string for name of second bookmark to compare
            selection_type (str): "+" or "-" indicating if ligand names should ("+") or should not "-" be in both databases
            old_db (str, optional): file name for previous database

        Returns:
            tuple: (name of new bookmark (str), number of ligands passing new bookmark (int))
        """

        if old_db is not None:
            self._detach_db(old_db.split(".")[0])  # remove file extension

        new_db_name = new_db.split(".")[0]  # remove file extension

        self._attach_db(new_db, new_db_name)

        if selection_type == "-":
            select_str = "NOT IN"
        elif selection_type == "+":
            select_str = "IN"
        else:
            raise StorageError(f"Unrecognized selection type {selection_type}")

        temp_name = "temp_" + str(self.temptable_suffix)
        self._create_temp_table(temp_name)
        temp_insert_query = self._generate_selective_insert_query(
            bookmark1_name, bookmark2_name, select_str, new_db_name, temp_name
        )

        self._insert_into_temp_table(temp_insert_query)

        num_passing = self._get_number_passing_ligands(temp_name)

        self.temptable_suffix += 1

        return temp_name, num_passing

    def prune(self):
        """Deletes rows from results, ligands, and interactions in a bookmark
        if they do not pass filtering criteria
        """
        self._delete_from_results()
        self._delete_from_ligands()
        self._delete_from_interactions_not_in_view()

    # endregion

    # region common keywords
    @classmethod
    def _data_kw_groups(cls, group):
        """Method containing lists of keywords in specific data groups, used to associate data with database columns

        Args:
            group (str): group of whose keywords are needed, including
                stateVar_keys
                ligand_data_keys
                interaction_data_kws
                outfield_options

        Returns:
            list: of keywords belonging to a specific group
        """
        groups = {
            "stateVar_keys": ["pose_about", "pose_translations", "pose_quarternions"],
            "ligand_data_keys": [
                "cluster_rmsds",
                "ref_rmsds",
                "scores",
                "leff",
                "delta",
                "intermolecular_energy",
                "vdw_hb_desolv",
                "electrostatics",
                "flex_ligand",
                "flexLigand_flexReceptor",
                "internal_energy",
                "torsional_energy",
                "unbound_energy",
            ],
            "interaction_data_kws": [
                "type",
                "chain",
                "residue",
                "resid",
                "recname",
                "recid",
            ],
            "outfield_options": [
                "Ligand_name",
                "e",
                "le",
                "delta",
                "ref_rmsd",
                "e_inter",
                "e_vdw",
                "e_elec",
                "e_intra",
                "n_interact",
                "interactions",
                "fname",
                "ligand_smile",
                "rank",
                "run",
                "hb",
                "source_file",
            ],
        }
        return groups[group]

    field_to_column_name = {
        "Ligand_name": "LigName",
        "e": "docking_score",
        "le": "leff",
        "delta": "deltas",
        "ref_rmsd": "reference_rmsd",
        "e_inter": "energies_inter",
        "e_vdw": "energies_vdw",
        "e_elec": "energies_electro",
        "e_intra": "energies_intra",
        "n_interact": "nr_interactions",
        "interactions": "interactions",
        "ligand_smile": "ligand_smile",
        "rank": "pose_rank",
        "run": "run_number",
        "hb": "num_hb",
        "receptor": "receptor",
    }
    # endregion


class StorageManagerSQLite(StorageManager):
    """SQLite-specific StorageManager subclass

    Attributes:
        conn (SQLite.conn): Connection to database
        open_cursors (list): list of cursors that were not closed by the function that created them.
            Will be closed by close_connection method.
        db_file (str): database name
        overwrite (bool): switch to overwrite database if it exists
        order_results (str): what column name will be used to order results once read
        outfields (str): data fields/columns to include when reading and outputting data
        filter_bookmark (str): name of bookmark that filtering will be performed over
        output_all_poses (bool): whether or not to output all poses of a ligand
        mfpt_cluster (float): distance in ångströms to cluster ligands based on morgan fingerprints
        interaction_cluster (float): distance in ångströms to cluster ligands based on interactions
        bookmark_name (str): name of current bookmark being written to or read from
        duplicate_handling (str): optional attribute to deal with insertion of ligands already in the database

        current_bookmark_name (str): name of last view to have been written to in the database
        filtering_window (str): name of bookmark/view being filtered on
        index_columns (list)
        view_suffix (int): current suffix for views
        temptable_suffix (int): current suffix for temporary tables
        field_to_column_name (dict): Dictionary for converting ringtail options into DB column names
    """

    def __init__(
        self,
        db_file: str = None,
        overwrite: bool = None,
        order_results: str = None,
        outfields: str = None,
        filter_bookmark: str = None,
        output_all_poses: bool = None,
        mfpt_cluster: float = None,
        interaction_cluster: float = None,
        bookmark_name: str = None,
        duplicate_handling: str = None,
    ):
        self.db_file = db_file
        self.overwrite = overwrite
        self.order_results = order_results
        self.outfields = outfields
        self.output_all_poses = output_all_poses
        self.mfpt_cluster = mfpt_cluster
        self.interaction_cluster = interaction_cluster
        self.filter_bookmark = filter_bookmark
        self.bookmark_name = bookmark_name
        self.duplicate_handling = duplicate_handling
        super().__init__()

        self.energy_filter_sqlite_call_dict = {
            "eworst": "docking_score < {value}",
            "ebest": "docking_score > {value}",
            "leworst": "leff < {value}",
            "lebest": "leff > {value}",
        }
        self.view_suffix = None
        self.temptable_suffix = 0
        self.open_cursors = []

    # region Methods for inserting into/removing from the database
    def _create_tables(self):
        """
        Creates all tables needed for a Ringtail database of a specific version
        """
        self._create_results_table()
        self._create_ligands_table()
        self._create_receptors_table()
        self._create_interaction_index_table()
        self._create_interaction_table()
        self._create_bookmark_table()
        self._create_db_properties_table()

    @classmethod
    def format_for_storage(cls, ligand_dict: dict) -> tuple:
        """takes file dictionary from the file parser, formats required storage format

        Args:
            ligand_dict (dict): Dictionary containing data from the fileparser

        Returns:
            tuple: of lists ([result_row_1, result_row_2,...],
                    ligand_row,
                    [interaction_tuple_1, interaction_tuple_2, ...])
        """

        # initialize row holders
        result_rows = []
        interaction_dictionaries = []
        interaction_tuples = []
        saved_pose_idx = 0  # save index of last saved pose
        cluster_saved_pose_map = {}  # save mapping of cluster number to saved_pose_idx

        # do the actual result formating
        # For each run we save, we add its interaction dict to the interaction_dictionaries list and save its other data
        # We also save a mapping of the its cluster number to the index in interaction_dictionaries
        # Then, when we find a pose to tolerate interactions for, we lookup the index to append the interactions to from cluster_saved_pose_map
        # Finally, we calculate the interaction tuple lists for each pose
        for idx, run_number in enumerate(ligand_dict["sorted_runs"]):
            cluster = ligand_dict["cluster_list"][idx]
            # save everything if this is a cluster top pose
            if run_number in ligand_dict["poses_to_save"]:
                result_rows.append(
                    cls._generate_results_row(ligand_dict, idx, run_number)
                )
                cluster_saved_pose_map[cluster] = saved_pose_idx
                saved_pose_idx += 1
                if ligand_dict["interactions"] != []:
                    interaction_dictionaries.append([ligand_dict["interactions"][idx]])
            elif run_number in ligand_dict["tolerated_interaction_runs"]:
                # adds to list started by best-scoring pose in cluster
                if cluster not in cluster_saved_pose_map:
                    continue
                interaction_dictionaries[cluster_saved_pose_map[cluster]].append(
                    ligand_dict["interactions"][idx]
                )

        for idx, pose_interactions in enumerate(interaction_dictionaries):
            if not any(pose_interactions):  # skip any empty dictionaries
                continue
            interaction_tuples.append(
                cls._generate_interaction_tuples(pose_interactions)
            )
        return (
            result_rows,
            cls._generate_ligand_row(ligand_dict),
            interaction_tuples,
            cls._generate_receptor_row(ligand_dict),
        )

    def _create_ligands_table(self):
        """Create table for ligands. Columns are:
        LigName             VARCHAR NOT NULL,
        ligand_smile        VARCHAR[],
        atom_index_map      VARCHAR[],
        hydrogen_parents    VARCHAR[],
        input_model         VARCHAR[]

        Raises:
            DatabaseTableCreationError: Description

        """
        ligand_table = """CREATE TABLE IF NOT EXISTS Ligands (
            LigName             VARCHAR NOT NULL PRIMARY KEY ON CONFLICT IGNORE,
            ligand_smile        VARCHAR[],
            ligand_rdmol        MOL,
            atom_index_map      VARCHAR[],
            hydrogen_parents    VARCHAR[],
            input_model         VARCHAR[])"""

        try:
            cur = self.conn.cursor()
            cur.execute(ligand_table)
            cur.close()
        except sqlite3.OperationalError as e:
            raise DatabaseTableCreationError(
                "Error while creating ligands table. If database already exists, use --overwrite to drop existing tables"
            ) from e

    @classmethod
    def _generate_ligand_row(cls, ligand_dict):
        """writes row to be inserted into ligand table

        Args:
            ligand_dict (dict): Dictionary of ligand data from parser

        Returns:
            List: List of data to be written as row in ligand table. Format:
                [ligand_name, ligand_smile, ligand_index_map,
                ligand_h_parents, input_model]
        """
        ligand_name = ligand_dict["ligname"]
        ligand_smile = ligand_dict["ligand_smile_string"]
        ligand_index_map = json.dumps(ligand_dict["ligand_index_map"])
        ligand_h_parents = json.dumps(ligand_dict["ligand_h_parents"])
        input_model = json.dumps(ligand_dict["ligand_input_model"])

        return [
            ligand_name,
            ligand_smile,
            ligand_index_map,
            ligand_h_parents,
            input_model,
        ]

    def _insert_ligands(self, ligand_array):
        """Takes array of ligand rows, inserts into Ligands table.

        Args:
            ligand_array (np.ndarray): Numpy array of arrays
                containing formatted ligand rows

        Raises:
            DatabaseInsertionError

        """

        sql_insert = """INSERT INTO Ligands (
        LigName,
        ligand_smile,
        ligand_rdmol,
        atom_index_map,
        hydrogen_parents,
        input_model
        ) VALUES
        (?,?,mol_from_smiles(?),?,?,?)"""

        ## repeat smiles in the third position of ligand array, to create rdmol
        for ligand_entry in ligand_array:
            smiles = ligand_entry[1]
            ligand_entry.insert(2, smiles)

        try:
            cur = self.conn.cursor()
            cur.executemany(sql_insert, ligand_array)
            self.conn.commit()
            cur.close()

        except sqlite3.OperationalError as e:
            raise DatabaseInsertionError("Error while inserting ligands.") from e

    def _delete_from_ligands(self):
        """Remove rows from ligands table if they did not pass filtering

        Raises:
            StorageError
        """
        try:
            cur = self.conn.cursor()
            cur.execute(
                "DELETE FROM Ligands WHERE LigName NOT IN (SELECT LigName from Results WHERE Pose_ID IN (SELECT Pose_ID FROM {view}))".format(
                    view=self.bookmark_name
                )
            )
            self.conn.commit()
            cur.close()
        except sqlite3.OperationalError as e:
            raise StorageError(
                f"Error occured while pruning Ligands not in {self.bookmark_name}"
            ) from e

    def _create_results_table(self):
        """Creates table for results. Columns are:
        Pose_ID             INTEGER PRIMARY KEY AUTOINCREMENT,
        LigName             VARCHAR NOT NULL,
        receptor            VARCHAR[],
        pose_rank           INT[],
        run_number          INT[],
        docking_score    FLOAT(4),
        leff                FLOAT(4),
        deltas              FLOAT(4),
        cluster_rmsd        FLOAT(4),
        cluster_size        INT[],
        reference_rmsd      FLOAT(4),
        energies_inter      FLOAT(4),
        energies_vdw        FLOAT(4),
        energies_electro    FLOAT(4),
        energies_flexLig    FLOAT(4),
        energies_flexLR     FLOAT(4),
        energies_intra      FLOAT(4),
        energies_torsional  FLOAT(4),
        unbound_energy      FLOAT(4),
        nr_interactions     INT[],
        num_hb              INT[],
        about_x             FLOAT(4),
        about_y             FLOAT(4),
        about_z             FLOAT(4),
        trans_x             FLOAT(4),
        trans_y             FLOAT(4),
        trans_z             FLOAT(4),
        axisangle_x         FLOAT(4),
        axisangle_y         FLOAT(4),
        axisangle_z         FLOAT(4),
        axisangle_w         FLOAT(4),
        dihedrals           VARCHAR[],
        ligand_coordinates         VARCHAR[],
        flexible_res_coordinates   VARCHAR[]

        Raises:
            DatabaseTableCreationError: Description
        """

        sql_results_table = """CREATE TABLE IF NOT EXISTS Results (
            Pose_ID             INTEGER PRIMARY KEY AUTOINCREMENT,
            LigName             VARCHAR NOT NULL,
            receptor            VARCHAR[],
            pose_rank           INT[],
            run_number          INT[],
            docking_score    FLOAT(4),
            leff                FLOAT(4),
            deltas              FLOAT(4),
            cluster_rmsd        FLOAT(4),
            cluster_size        INT[],
            reference_rmsd      FLOAT(4),
            energies_inter      FLOAT(4),
            energies_vdw        FLOAT(4),
            energies_electro    FLOAT(4),
            energies_flexLig    FLOAT(4),
            energies_flexLR     FLOAT(4),
            energies_intra      FLOAT(4),
            energies_torsional  FLOAT(4),
            unbound_energy      FLOAT(4),
            nr_interactions     INT[],
            num_hb              INT[],
            about_x             FLOAT(4),
            about_y             FLOAT(4),
            about_z             FLOAT(4),
            trans_x             FLOAT(4),
            trans_y             FLOAT(4),
            trans_z             FLOAT(4),
            axisangle_x         FLOAT(4),
            axisangle_y         FLOAT(4),
            axisangle_z         FLOAT(4),
            axisangle_w         FLOAT(4),
            dihedrals           VARCHAR[],
            ligand_coordinates         VARCHAR[],
            flexible_res_coordinates   VARCHAR[]
            ); """

        try:
            cur = self.conn.cursor()
            cur.execute(sql_results_table)
            cur.close()
        except sqlite3.OperationalError as e:
            raise DatabaseTableCreationError(
                "Error while creating results table. If database already exists, use 'overwrite' to drop existing tables"
            ) from e

    @classmethod
    def _generate_results_row(cls, ligand_dict, pose_rank, run_number):
        """generate list of lists of ligand values to be
            inserted into sqlite database

        Args:
            ligand_dict (dict): Dictionary of ligand data from parser
            pose_rank (int): Rank of pose to generate row for
                all runs for the given ligand
            run_number (int): Run number of pose to generate row for
                all runs for the given ligand

        Returns:
            List: List of pose data to be inserted into Results table.
            In same order as expected in insert_results:
            LigName, [0]
            receptor, [2]
            pose_rank, [3]
            run_number, [4]
            cluster_rmsd, [5]
            reference_rmsd, [6]
            docking_score, [7]
            leff, [8]
            deltas, [9]
            energies_inter, [10]
            energies_vdw, [11]
            energies_electro, [12]
            energies_flexLig, [13]
            energies_flexLR, [14]
            energies_intra, [15]
            energies_torsional, [16]
            unbound_energy, [17]
            nr_interactions, [18]
            num_hb, [19]
            cluster_size, [20]
            about_x, [21]
            about_y, [22]
            about_z, [23]
            trans_x, [24]
            trans_y, [25]
            trans_z, [26]
            axisangle_x, [27]
            axisangle_y, [28]
            axisangle_z, [29]
            axisangle_w, [30]
            dihedrals, [31]
            ligand_coordinates, [32]
            flexible_res_coordinates [33]
        """

        # # # # # # get pose-specific data

        # check if run is best for a cluster.
        # We are only saving the top pose for each cluster
        ligand_data_list = [
            ligand_dict["ligname"],
            ligand_dict["receptor"],
            pose_rank + 1,
            int(run_number),
        ]
        # get energy data
        for key in cls._data_kw_groups("ligand_data_keys"):
            if ligand_dict[key] == []:  # guard against incomplete data
                ligand_data_list.append(None)
            else:
                ligand_data_list.append(ligand_dict[key][pose_rank])
        if ligand_dict["interactions"] != [] and any(
            ligand_dict["interactions"][pose_rank]
        ):  # catch lack of interaction data
            # add interaction count
            ligand_data_list.append(ligand_dict["interactions"][pose_rank]["count"][0])
            if int(ligand_dict["interactions"][pose_rank]["count"][0]) != 0:
                # count number H bonds, add to ligand data list
                ligand_data_list.append(
                    ligand_dict["interactions"][pose_rank]["type"].count("H")
                )
            else:
                ligand_data_list.append(0)
        else:
            ligand_data_list.extend(
                [
                    None,
                    None,
                ]
            )
        # Add the cluster size for the cluster this pose belongs to
        ligand_data_list.append(
            ligand_dict["cluster_sizes"][ligand_dict["cluster_list"][pose_rank]]
        )
        # add statevars
        for key in cls._data_kw_groups("stateVar_keys"):
            if ligand_dict[key] == []:
                if key == "pose_about" or key == "pose_translations":
                    ligand_data_list.extend(
                        [
                            None,
                            None,
                            None,
                        ]
                    )
                if key == "pose_quarternions":
                    ligand_data_list.extend(
                        [
                            None,
                            None,
                            None,
                            None,
                        ]
                    )
                continue
            stateVar_data = ligand_dict[key][pose_rank]
            if stateVar_data != []:
                for dim in stateVar_data:
                    ligand_data_list.append(dim)
        dihedral_string = ""
        if ligand_dict["pose_dihedrals"] != []:
            pose_dihedrals = ligand_dict["pose_dihedrals"][pose_rank]
            for dihedral in pose_dihedrals:
                dihedral_string = dihedral_string + json.dumps(dihedral) + ", "
        ligand_data_list.append(dihedral_string)

        # add coordinates
        # convert to string for storage as VARCHAR
        ligand_data_list.append(json.dumps(ligand_dict["pose_coordinates"][pose_rank]))
        ligand_data_list.append(
            json.dumps(ligand_dict["flexible_res_coordinates"][pose_rank])
        )

        return ligand_data_list

    def _check_unique_results_row(self, result_data):
        """Checks if a pose ID is uniquely represented in the result table, based on the following [index in result_data] columns:
        [0] LigName,
        [1] receptor,
        [20] about_x,
        [21] about_y,
        [22] about_z,
        [23] trans_x,
        [24] trans_y,
        [25] trans_z,
        [26] axisangle_x,
        [27] axisangle_y,
        [28] axisangle_z,
        [29] axisangle_w,
        [30] dihedrals,

        #NOTE Please note that this method will only identify one duplicate in the table. If there are more than one duplicates, it will just deal with the earliest entry

        Args:
            result_data (list): data packet coming from the results processing

        Raises:
            DatabaseQueryError

        Returns:
            Pose_ID (int): returns the Pose_ID of the duplicate if found, returns -1 of no duplicate found

        """
        # create list of the data that is to be considered unique
        unique_data_indices = [0, 1, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30]
        unique_data = [result_data[index] for index in unique_data_indices]

        try:
            cur = self.conn.cursor()
            query = """SELECT Pose_ID 
                        FROM Results 
                        WHERE 
                        LigName=?
                        AND receptor=?
                        AND about_x=?
                        AND about_y=?
                        AND about_z=?
                        AND trans_x=?
                        AND trans_y=?
                        AND trans_z=?
                        AND axisangle_x=?
                        AND axisangle_y=?
                        AND axisangle_z=?
                        AND axisangle_w=?
                        AND dihedrals=?;"""

            cur.execute(query, unique_data)
            row = cur.fetchone()
            if row is None:
                Pose_ID = -1
                self.logger.debug("Duplicate row not found.")
            else:
                Pose_ID = row[0]
                self.logger.debug(f"Duplicate row found for Pose_ID {Pose_ID}")
            cur.close()

            return Pose_ID

        except sqlite3.OperationalError as e:
            raise DatabaseQueryError(
                "Error while looking for unique result row."
            ) from e

    def _insert_results(self, results_array):
        """Takes array of database rows to insert, adds data to results table. Will handle duplicates if specified

        Args:
            results_array (np.ndAaray): numpy array of arrays containing
                formatted result rows

        Returns:
            Pose_ID (list(int)): returns the pose ids for the ligand written to results, these are used to ensure internal consistency when writing to the interaction table
            duplicates (list(int)): list of pose ids that are duplicates, if duplicate handling is specified. Filled with None if not specified or not duplicate

        Raises:
            DatabaseInsertionError
        """

        sql_insert = """INSERT INTO Results (
                        LigName,
                        receptor,
                        pose_rank,
                        run_number,
                        cluster_rmsd,
                        reference_rmsd,
                        docking_score,
                        leff,
                        deltas,
                        energies_inter,
                        energies_vdw,
                        energies_electro,
                        energies_flexLig,
                        energies_flexLR,
                        energies_intra,
                        energies_torsional,
                        unbound_energy,
                        nr_interactions,
                        num_hb,
                        cluster_size,
                        about_x,
                        about_y,
                        about_z,
                        trans_x,
                        trans_y,
                        trans_z,
                        axisangle_x,
                        axisangle_y,
                        axisangle_z,
                        axisangle_w,
                        dihedrals,
                        ligand_coordinates,
                        flexible_res_coordinates
                        ) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?);"""

        try:
            Pose_IDs = []
            duplicates = []
            cur = self.conn.cursor()
            # for each pose/docking result
            for result in results_array:
                Pose_ID = (
                    -1
                )  # nonsensical table index to initialize row index if checking for duplicates
                if self.duplicate_handling:
                    Pose_ID = self._check_unique_results_row(result)

                if Pose_ID != -1:  # row exists in table
                    duplicates.append(Pose_ID)
                    # row exist, evaluate if ignore or replace
                    if self.duplicate_handling.upper() == "IGNORE":
                        # do not add the new, duplicated row
                        pass
                    elif self.duplicate_handling.upper() == "REPLACE":
                        # update the existing row with the new results
                        # reformat sqlite query to update
                        sql_replace = sql_insert.replace(
                            "INSERT INTO Results", "UPDATE Results SET"
                        )
                        sql_replace = sql_replace.replace("VALUES", "=")
                        sql_replace = sql_replace.replace(";", " WHERE Pose_ID = ?;")
                        result.append(
                            Pose_ID
                        )  # add pose ID to the data being processed in sqlite statement
                        cur.execute(sql_replace, result)

                else:  # row does not exist
                    duplicates.append(None)
                    cur.execute(sql_insert, result)
                    Pose_ID = cur.lastrowid
                # create list of pose ids just processed
                Pose_IDs.append(Pose_ID)

            self.conn.commit()
            cur.close()

            return Pose_IDs, duplicates

        except sqlite3.OperationalError as e:
            raise DatabaseInsertionError("Error while inserting results.") from e

    def _delete_from_results(self):
        """Remove rows from results table if they did not pass filtering

        Raises:
            StorageError
        """
        try:
            cur = self.conn.cursor()
            cur.execute(
                "DELETE FROM Results WHERE Pose_ID NOT IN (SELECT Pose_ID FROM {view})".format(
                    view=self.bookmark_name
                )
            )
            self.conn.commit()
            cur.close()
        except sqlite3.OperationalError as e:
            raise StorageError(
                f"Error occured while pruning Results not in {self.bookmark_name}"
            ) from e

    def _create_receptors_table(self):
        """Create table for receptors. Columns are:
        Receptor_ID         INTEGER PRIMARY KEY AUTOINCREMENT,
        RecName             VARCHAR,
        box_dim             VARCHAR[],
        box_center          VARCHAR[],
        grid_spacing        FLOAT(4),
        flexible_residues   VARCHAR[],
        flexres_atomnames   VARCHAR[],
        receptor_object     BLOB

        Raises:
            DatabaseTableCreationError: Description
        """
        receptors_table = """CREATE TABLE IF NOT EXISTS Receptors (
            Receptor_ID         INTEGER PRIMARY KEY AUTOINCREMENT,
            RecName             VARCHAR,
            box_dim             VARCHAR[],
            box_center          VARCHAR[],
            grid_spacing        FLOAT(4),
            flexible_residues   VARCHAR[],
            flexres_atomnames   VARCHAR[],
            receptor_object     BLOB
        )"""

        try:
            cur = self.conn.cursor()
            cur.execute(receptors_table)
            cur.close()
        except sqlite3.OperationalError as e:
            raise DatabaseTableCreationError(
                "Error while creating receptor table. If database already exists, use --overwrite to drop existing tables"
            ) from e

    def _insert_receptors(self, receptor_array):
        """Takes array of receptor rows, inserts into Receptors table

        Args:
            receptor_array (list): List of lists
                containing formatted ligand rows

        Raises:
            DatabaseInsertionError
        """
        sql_insert = """INSERT INTO Receptors (
        RecName,
        box_dim,
        box_center,
        grid_spacing,
        flexible_residues,
        flexres_atomnames
        ) VALUES
        (?,?,?,?,?,?)"""

        try:
            cur = self.conn.cursor()
            cur.executemany(sql_insert, receptor_array)
            self.conn.commit()
            cur.close()

        except sqlite3.OperationalError as e:
            raise DatabaseInsertionError("Error while inserting receptor.") from e

    @classmethod
    def _generate_receptor_row(cls, ligand_dict):
        """Writes row to be inserted into receptor table

        Args:
            ligand_dict (dict): Dictionary of ligand data from parser
        """

        rec_name = ligand_dict["receptor"]
        box_dim = json.dumps(ligand_dict["grid_dim"])
        box_center = json.dumps(ligand_dict["grid_center"])
        grid_spacing = ligand_dict["grid_spacing"]
        if grid_spacing != "":
            grid_spacing = float(grid_spacing)
        flexible_residues = json.dumps(ligand_dict["flexible_residues"])
        flexres_atomnames = json.dumps(ligand_dict["flexres_atomnames"])

        return [
            rec_name,
            box_dim,
            box_center,
            grid_spacing,
            flexible_residues,
            flexres_atomnames,
        ]

    def insert_receptor_blob(self, receptor, rec_name):
        """Takes object of Receptor class, updates the column in Receptor table

        Args:
            receptor (bytes): bytes receptor object to be inserted into DB
            rec_name (string): Name of receptor. Used to insert into correct row of DB

        Raises:
            DatabaseInsertionError: Description
        """
        # Check if there is already a row for the receptor
        cur = self.conn.execute("SELECT COUNT(*) FROM Receptors")
        count = cur.fetchone()[0]
        if count == 0:
            # Insert receptor statement
            query = f"""INSERT INTO Receptors (
                      RecName,
                      receptor_object)
                      VALUES (?,?)"""

        else:
            query = """UPDATE Receptors SET RecName = ?, receptor_object = ? WHERE Receptor_ID == 1"""
        try:
            cur = self.conn.execute(query, (rec_name, receptor))
            self.conn.commit()
            cur.close()
        except sqlite3.OperationalError as e:
            raise DatabaseInsertionError(
                "Error while adding receptor blob to database"
            ) from e

    def _create_db_properties_table(self):
        """Create table of database properties used during write session to the database. Columns are:
        DB_write_session int (primary key)
        docking_mode (vina or dlg)
        num_of_poses ("all" or int)

        Raises:
            DatabaseTableCreationError
        """
        sql_str = """CREATE TABLE IF NOT EXISTS DB_properties (
        DB_write_session    INTEGER PRIMARY KEY AUTOINCREMENT,
        docking_mode        VARCHAR[],
        number_of_poses     VARCHAR[])"""

        try:
            cur = self.conn.cursor()
            cur.execute(sql_str)
            cur.close()
        except sqlite3.OperationalError as e:
            raise DatabaseTableCreationError(
                "Error while creating db properties table. If database already exists, use --overwrite to drop existing tables"
            ) from e

    def _insert_db_properties(self, docking_mode: str, number_of_poses: str):
        """Insert db properties into database properties table

        Args:
            docking_mode (str): docking mode for the current dataset being written
            number_of_poses (str): number of poses written to database in current session, either "all" or specified max_poses

        Raises:
            DatabaseInsertionError
        """
        sql_insert = """INSERT INTO DB_properties (
        docking_mode,
        number_of_poses
        ) VALUES (?,?)"""

        try:
            cur = self.conn.cursor()
            cur.execute(sql_insert, [docking_mode, number_of_poses])
            self.conn.commit()
            cur.close()

        except sqlite3.OperationalError as e:
            raise DatabaseInsertionError(
                "Error while inserting database properties info into DB_properties table"
            ) from e

    def _create_interaction_index_table(self):
        """create table of data for each unique interaction, will be remade everytime db is written to.
        Columns are:
        interaction_id      INTEGER PRIMARY KEY AUTOINCREMENT,
        interaction_type    VARCHAR[],
        rec_chain           VARCHAR[],
        rec_resname         VARCHAR[],
        rec_resid           VARCHAR[],
        rec_atom            VARCHAR[],
        rec_atomid          VARCHAR[]

        Raises:
            DatabaseTableCreationError: Description

        """
        interaction_index_table = """CREATE TABLE Interaction_indices (
                                        interaction_id      INTEGER PRIMARY KEY,
                                        interaction_type    VARCHAR[],
                                        rec_chain           VARCHAR[],
                                        rec_resname         VARCHAR[],
                                        rec_resid           VARCHAR[],
                                        rec_atom            VARCHAR[],
                                        rec_atomid          VARCHAR[],
                                        UNIQUE (interaction_type, rec_chain, rec_resname, rec_resid, rec_atom, rec_atomid) ON CONFLICT IGNORE );
                                        """

        try:
            cur = self.conn.cursor()
            cur.execute("""DROP TABLE IF EXISTS Interaction_indices""")
            cur.execute(interaction_index_table)
            cur.close()
        except sqlite3.OperationalError as e:
            raise DatabaseTableCreationError(
                f"Error while creating interaction index table: {e}"
            ) from e

    def _create_interaction_table(self):
        """Create table a "tall-skinny" table of each pose-interaction.
        This table enables proper handling of duplicates if specified.
        Columns are:
        interaction_pose_id INTERGER PRIMARY KEY AUTOINCREMENT,
        Pose_ID             INTEGER FOREIGN KEY from RESULTS,
        interaction_id      INTEGER FOREIGN KEY from Interaction_indices

        Raises:
            DatabaseTableCreationError: Description
        """

        interaction_table = """CREATE TABLE IF NOT EXISTS Interactions (
        interaction_pose_ID INTEGER PRIMARY KEY AUTOINCREMENT,
        Pose_ID   INTEGER,
        interaction_id INTEGER,
        FOREIGN KEY (Pose_ID) REFERENCES Results(Pose_ID),
        FOREIGN KEY (interaction_id) REFERENCES Interaction_indices(interaction_id))"""

        try:
            cur = self.conn.cursor()
            cur.execute(interaction_table)
            cur.close()
        except sqlite3.OperationalError as e:
            raise DatabaseTableCreationError(
                "Error while creating interactions table. If database already exists, use 'overwrite' to drop existing tables"
            ) from e

    def _insert_interaction_rows(self, interaction_rows, duplicates):
        """Inserts the interaction data into a "tall-and-skinny" table, with a primary autoincremented key and a Pose_ID that is 1-to-1 with Results table.
        Table will contain as many rows with the same Pose_ID as that pose has interactions.

        Args:
            interaction_rows (list(tuple)): list of tuples containing the interaction data
            duplicates (list(int)): list of pose_ids from results table deemed duplicates, can also contain Nones, will be treated according to self.duplicate_handling

        Raises:
            DatabaseInsertionError
        """
        sql_insert = """INSERT INTO Interactions 
                            (Pose_ID,
                            interaction_id)
                            VALUES (?,?);"""
        try:
            cur = self.conn.cursor()
            if not self.duplicate_handling:  # add all results
                cur.executemany(sql_insert, interaction_rows)
            else:
                # first, add any poses that are not duplicates
                non_duplicates = [
                    interaction_row
                    for interaction_row in interaction_rows
                    if interaction_row[0] not in duplicates
                ]
                # check if there are duplicates or if duplicates list contains only None
                duplicates_exist = bool(duplicates.count(None) != len(duplicates))
                cur.executemany(sql_insert, non_duplicates)

                # only look for values to replace if there are duplicate pose ids
                if self.duplicate_handling == "REPLACE" and duplicates_exist:
                    # delete all rows pertaining to duplicated pose_ids
                    duplicated_pose_ids = [id for id in duplicates if id is not None]
                    self._delete_interactions(duplicated_pose_ids)
                    # insert the interaction tuples for the new pose_ids
                    duplicates_only = [
                        interaction_row
                        for interaction_row in interaction_rows
                        if interaction_row[0] in duplicates
                    ]
                    cur.executemany(sql_insert, duplicates_only)

                elif self.duplicate_handling == "IGNORE":
                    # ignore and don't add any poses that are duplicates
                    pass
            self.conn.commit()
            cur.close()

        except sqlite3.OperationalError as e:
            raise DatabaseInsertionError(
                f"Error while inserting an interaction row: {e}"
            ) from e

    @classmethod
    def _generate_interaction_tuples(cls, interaction_dictionaries: list):
        """takes dictionary of file results, formats as
        list of tuples for interactions

        Args:
            interaction_dictionaries (list): List of pose interaction
                dictionaries from parser

        Returns:
            list: List of tuples of interaction data
        """
        interactions = set()
        for pose_interactions in interaction_dictionaries:
            count = pose_interactions["count"][0]
            for i in range(int(count)):
                interactions.add(
                    tuple(
                        pose_interactions[kw][i]
                        for kw in cls._data_kw_groups("interaction_data_kws")
                    )
                )

        return list(interactions)

    def _insert_interaction_index_row(self, interaction_tuple) -> tuple:
        """
        Writes unique interactions and returns the interaction_id of the given interaction

        Args:
            interaction_tuple (tuple): (rec_chain, rec_resname, rec_resid, rec_atom, rec_atomid)

        Returns:
            tuple: if interaction index (int_index,)

        Raises:
            DatabaseInsertionError
        """
        # to insert interaction if unique
        sql_insert = """INSERT OR IGNORE INTO Interaction_indices (interaction_id, interaction_type,rec_chain,rec_resname,rec_resid,rec_atom,rec_atomid) 
                        VALUES (?,?,?,?,?,?,?);"""
        # to get interaction_id from the given interaction
        sql_query = f"""SELECT interaction_id FROM Interaction_indices 
        WHERE interaction_type = ?
        AND rec_chain = ?
        AND rec_resname = ?
        AND rec_resid = ?
        AND rec_atom = ?
        AND rec_atomid = ?"""

        try:
            cur = self.conn.cursor()
            cur.execute(sql_query, interaction_tuple)
            self.conn.commit()
            interaction_index = cur.fetchall()
            if not interaction_index:
                # get table length and use that as index
                interaction_index = (self._get_length_of_table("Interaction_indices"),)
                # create and insert new interaction id
                input_tuple = interaction_index + interaction_tuple
                cur.execute(sql_insert, input_tuple)
                self.conn.commit()
            else:
                interaction_index = interaction_index[0]
            cur.close()
            return interaction_index
        except sqlite3.OperationalError as e:
            raise DatabaseInsertionError(
                f"Error inserting unique interaction tuples in index table: {e}"
            ) from e

    def _delete_interactions(self, Pose_IDs):
        """Remove rows from interactions table where pose id is represented in Pose_IDs

        Args:
            Pose_IDs (list(int)): list of pose ids to delete from the table

        Raises:
            StorageError: Description
        """
        Pose_IDs_string = ",".join(map(str, Pose_IDs))
        sql_delete = f"DELETE FROM Interactions WHERE Pose_ID IN ({Pose_IDs_string});"
        try:
            cur = self.conn.cursor()
            cur.execute(sql_delete)
            self.conn.commit()
            cur.close()

        except sqlite3.OperationalError as e:
            raise StorageError(
                "Error while deleting rows in the Interaction table"
            ) from e

    def _delete_from_interactions_not_in_view(self):
        """Remove rows from interactions table if they did not pass filtering.

        Raises:
            StorageError: Description
        """
        try:
            cur = self.conn.cursor()
            cur.execute(
                "DELETE FROM Interactions WHERE Pose_ID NOT IN (SELECT Pose_ID FROM {view})".format(
                    view=self.bookmark_name
                )
            )
            self.conn.commit()
            cur.close()
        except sqlite3.OperationalError as e:
            raise StorageError(
                f"Error occured while pruning Interactions not in {self.bookmark_name}"
            ) from e

    def _create_bookmark_table(self):
        """Create table of bookmark names and their queries. Columns are:
        Bookmark_name
        Query

        Raises:
            DatabaseTableCreationError: Description
        """
        sql_str = """CREATE TABLE IF NOT EXISTS Bookmarks (
        Bookmark_name       VARCHAR[] PRIMARY KEY,
        Query               VARCHAR[],
        filters             VARCHAR[])"""

        try:
            cur = self.conn.cursor()
            cur.execute(sql_str)
            cur.close()
        except sqlite3.OperationalError as e:
            raise DatabaseTableCreationError(
                "Error while creating bookmark table. If database already exists, use --overwrite to drop existing tables"
            ) from e

    def _insert_cluster_data(
        self, clusters: list, poseid_list: list, cluster_type: str, cluster_cutoff: str
    ):
        """Insert cluster data into ligand cluster table

        Args:
            clusters (list)
            poseid_list (list)
            cluster_type (str)
            cluster_cutoff (str)
        """
        cur = self.conn.cursor()
        cur.execute(
            "CREATE TABLE IF NOT EXISTS Ligand_clusters (pose_id  INT[] UNIQUE)"
        )
        ligand_cluster_columns = self._fetch_ligand_cluster_columns()
        column_name = (
            f"{self.bookmark_name}_{cluster_type}_{cluster_cutoff.replace('.', 'p')}"
        )
        if column_name not in ligand_cluster_columns:
            cur.execute(f"ALTER TABLE Ligand_clusters ADD COLUMN {column_name}")
        for ci, cl in enumerate(clusters):
            for i in cl:
                poseid = poseid_list[i]
                cur.execute(
                    f"INSERT INTO Ligand_clusters (pose_id, {column_name}) VALUES (?,?) ON CONFLICT (pose_id) DO UPDATE SET {column_name}=excluded.{column_name}",
                    (poseid, ci),
                )

        cur.close()
        self.conn.commit()

    def _create_indices(self):
        """Create index for specified tables and columns. 'ak' stands for 'alternate key' and is prepended to index name to avoid naming conflicts

        Raises:
            StorageError
        """
        try:
            cur = self.conn.cursor()
            self.logger.debug("Creating columns index...")

            cur.execute(
                "CREATE INDEX IF NOT EXISTS ak_results ON Results(LigName, docking_score, leff, deltas, reference_rmsd, energies_inter, energies_vdw, energies_electro, energies_intra, nr_interactions, run_number, pose_rank, num_hb)"
            )
            cur.execute("CREATE INDEX IF NOT EXISTS ak_poseid ON Results(Pose_id)")
            cur.execute(
                "CREATE INDEX IF NOT EXISTS ak_intind ON Interaction_indices(interaction_type, rec_chain, rec_resname, rec_resid, rec_atom, rec_atomid)"
            )
            cur.execute(
                "CREATE INDEX IF NOT EXISTS ak_interactions ON Interactions(Pose_id, interaction_id)"
            )
            cur.execute("CREATE INDEX IF NOT EXISTS ak_ligands ON Ligands(LigName)")
            self.conn.commit()
            cur.close()
            self.logger.info(
                "Indicies were created for specified Results and Interaction_indices columns."
            )
        except sqlite3.OperationalError as e:
            raise StorageError("Error occurred while indexing") from e

    def _delete_table(self, table_name: str):
        """
        Method to delete a table

        Args:
            table_name (str): table to be dropped

        """
        query = f"""DROP TABLE IF EXISTS {table_name};"""
        return self._run_query(query)

    # endregion

    # region Methods for dealing with bookmarks/views and temporary tables
    def get_all_bookmark_names(self):
        """Get all bookmarks in sql database as a list of names. Bookmarks are a specific type of sqlite-views
        whose information is stored in the Bookmarks table.

        Returns:
            list: of bookmark names
        """
        try:
            cur = self.conn.cursor()
            cur.execute("SELECT Bookmark_name FROM Bookmarks;")
            bookmark_names = [name[0] for name in cur.fetchall()]
            cur.close()

        except sqlite3.OperationalError as e:
            raise StorageError(
                "Error occured while fetching existing bookmark names"
            ) from e

        return bookmark_names

    def set_bookmark_suffix(self, suffix):
        """Sets internal bookmark_suffix variable

        Args:
            suffix (str): suffix to attached to bookmark-related queries or creation
        """
        if not isinstance(suffix, str):
            self.view_suffix = str(suffix)
        else:
            self.view_suffix = suffix

    def fetch_filters_from_bookmark(self, bookmark_name: str | None = None):
        """Method that will retrieve filter values used to construct bookmark

        Args:
            bookmark_name (str, optional): can get filter values for given bookmark, or filter values from currently active bookmark in storageman

            Returns:
                dict: containing the filter data
        """
        if bookmark_name is None:
            bookmark_name = self.bookmark_name
        sql_query = (
            f"SELECT filters FROM Bookmarks where Bookmark_name = '{bookmark_name}'"
        )
        filters = self._run_query(sql_query).fetchone()
        if not filters: 
            return {}

        return json.loads(filters[0])

    def bookmark_has_rows(self, bookmark_name: str) -> bool:
        """
        Method that checks if a given bookmark has any data in it

        Args:
            bookmark_name (str): view to check

        Returns:
            bool: True if more than zero rows in bookmark
        """

        count = self._run_query(
            "SELECT COUNT(*) FROM {0};".format(bookmark_name)
        ).fetchone()[0]
        return bool(count > 0)

    def get_current_bookmark_name(self):
        """returns current bookmark name

        Returns:
            str: name of last passing results bookmark used by database
        """
        return self.current_bookmark_name

    def fetch_bookmark(self, bookmark_name: str) -> sqlite3.Cursor:
        """returns SQLite cursor of all fields in bookmark

        Args:
            bookmark_name (str): name of bookmark to retrieve

        Returns:
            sqlite3.Cursor: cursor of requested view
        """
        return self._run_query(f"SELECT * FROM {bookmark_name}")

    def create_bookmark(self, name, query, temp=False, add_poseID=False, filters={}):
        """Takes name and selection query and creates a bookmark of name.
        Bookmarks are Ringtail specific views that whose information is stored in the 'Bookmark' table.
        #FIXME bug where ligand filter only results are not added as bookmarks

        Args:
            name (str): Name for bookmark which will be created
            query (str): SQLite-formated query used to create bookmark
            temp (bool, optional): Flag if bookmark should be temporary
            add_poseID (bool, optional): Add Pose_ID column to bookmark
            filters (dict, optional): a dict of filters used to construct the query
        """

        if add_poseID:
            query = query.replace("SELECT ", "SELECT Pose_ID, ", 1)
        self.logger.info("Creating bookmark...")

        if temp:
            temp_flag = "TEMP "
        else:
            temp_flag = ""

        bookmark_query = f"CREATE {temp_flag}VIEW {name} AS {query}"
        self._create_view(name, bookmark_query)
        self._insert_bookmark_info(name, bookmark_query, filters)
        self.logger.debug(
            f"Created bookmark from the following query: {bookmark_query}"
        )

    def _create_view(self, name, query):
        """takes name and selection query,
            creates view of query stored as name.

        Args:
            name (str): Name for view which will be created
            query (str): SQLite-formated query used to create view

        Raises:
            DatabaseViewCreationError
        """
        # check that name does not start with int, this causes a sqlite error
        if name[0].isdigit():
            raise DatabaseViewCreationError(
                f"The given view name {name} starts with a number, view names may not start with digit."
            )
        cur = self.conn.cursor()
        # drop old view if there is one
        try:
            cur.execute(f"DROP VIEW IF EXISTS {name}")
            cur.execute(query)
            cur.close()
        except sqlite3.OperationalError as e:
            raise DatabaseViewCreationError(
                f"Error ({e}) creating view from query \n{query}"
            ) from e

    def _insert_bookmark_info(self, name: str, sqlite_query: str, filters={}):
        """Insert bookmark info into bookmark table (bookmarks are a Ringtail specific type of views)

        Args:
            name (str): name for bookmark
            sqlite_query (str): sqlite query used to generate bookmark
            filters (dict): filters used to generate bookmark

        Raises:
            DatabaseInsertionError
        """
        sql_insert = """INSERT OR REPLACE INTO Bookmarks (
        Bookmark_name,
        Query,
        filters
        ) VALUES (?,?,?)"""

        try:
            cur = self.conn.cursor()
            cur.execute(sql_insert, [name, sqlite_query, json.dumps(filters)])
            self.conn.commit()
            cur.close()

        except sqlite3.OperationalError as e:
            raise DatabaseInsertionError(
                "Error while inserting Bookmark info into Bookmark table"
            ) from e

    def drop_bookmark(self, bookmark_name: str):
        """Drops specified bookmark from database

        Args:
            bookmark_name (str): bookmark to be dropped

        Raises:
            DatabaseInsertionError
        """

        query_drop = f"DROP VIEW IF EXISTS {bookmark_name}"
        query_delete = f"DELETE FROM Bookmarks WHERE Bookmark_name = '{bookmark_name}'"

        try:
            cur = self.conn.execute(query_drop)
            cur.execute(query_delete)
            self.conn.commit()
            cur.close()
            self.logger.info(f"Dropped bookmark {bookmark_name}.")
        except sqlite3.OperationalError as e:
            raise DatabaseInsertionError(
                f"Error while attempting to drop bookmark {bookmark_name}"
            ) from e

    def create_temp_table_from_bookmark(self):
        """Method that creates a temporary table named "passing_temp".
        Please note that this table will be dropped as soon as the database connection closes.
        """
        cur = self.conn.cursor()
        cur.execute(
            f"CREATE TEMP TABLE passing_temp AS SELECT * FROM {self.bookmark_name}"
        )
        cur.close()
        self.logger.debug(
            "Creating a temporary table of passing ligands named 'passing_temp'."
        )

    def create_bookmark_from_temp_table(
        self,
        temp_table_name,
        bookmark_name,
        original_bookmark_name,
        wanted_list,
        unwanted_list=[],
    ):
        """Resaves temp bookmark stored in self.current_bookmark_name as new permenant bookmark

        Args:
            bookmark_name (str): name of bookmark to save last temp bookmark as
            original_bookmark_name (str): name of original bookmark
            wanted_list (list): List of wanted database names
            unwanted_list (list, optional): List of unwanted database names
            temp_table_name (str): name of temporary table
        """
        self._create_view(
            bookmark_name,
            "SELECT * FROM {0} WHERE Pose_ID in (SELECT Pose_ID FROM {1})".format(
                original_bookmark_name, temp_table_name
            ),
            add_poseID=False,
        )
        compare_bookmark_str = "Comparison. Wanted: "
        compare_bookmark_str += ", ".join(wanted_list)
        if unwanted_list is not None:
            compare_bookmark_str += ". Unwanted: " + ", ".join(unwanted_list)
        self._insert_bookmark_info(bookmark_name, compare_bookmark_str)

    def _create_temp_table(self, table_name):
        """create temporary table with given name and with ligand name and pose_id information

        Args:
            table_name (str): name for temp table

        Raises:
            DatabaseTableCreationError
        """

        create_table_str = (
            f"CREATE TEMP TABLE {table_name}(Pose_ID PRIMARY KEY, LigName)"
        )
        try:
            cur = self.conn.cursor()
            cur.execute(create_table_str)
            self.conn.commit()
            cur.close()
        except sqlite3.OperationalError as e:
            raise DatabaseTableCreationError(
                f"Error while creating temporary table {table_name}"
            ) from e

    def _insert_into_temp_table(self, query):
        """Execute insertion into temporary table

        Args:
            query (str): Insertion command

        Raises:
            DatabaseInsertionError
        """
        try:
            cur = self.conn.cursor()
            cur.execute(query)
            self.conn.commit()
            cur.close()
        except sqlite3.OperationalError as e:
            raise DatabaseInsertionError(
                f"Error while inserting into temporary table"
            ) from e

    # endregion

    # region Methods for getting information from database
    def fetch_receptor_object_by_name(self, rec_name):
        """Returns Receptor object from database for given rec_name

        Args:
            rec_name (str): Name of receptor to return object for

        Returns:
        str: receptor object as a string
        """

        cursor = self._run_query(
            """SELECT receptor_object FROM Receptors WHERE RecName LIKE '{0}'""".format(
                rec_name
            )
        )
        return str(cursor.fetchone()[0])

    def fetch_receptor_objects(self):
        """Returns all Receptor objects from database

        Args:
            rec_name (str): Name of receptor to return object for

        Returns:
            iter (tuple): of receptor names and objects
        """

        cursor = self._run_query("SELECT RecName, receptor_object FROM Receptors")
        return cursor.fetchall()

    def fetch_data_for_passing_results(self) -> iter:
        """Will return SQLite cursor with requested data for outfields for poses that passed filter in self.bookmark_name

        Returns:
            iter: sqlite cursor of data from passing data

        Raises:
            OptionError
        """
        output_fields = self.outfields
        if type(output_fields) == str:
            output_fields = output_fields.replace(" ", "")
            output_fields_list = output_fields.split(",")
        elif type(output_fields) == list:
            output_fields_list = output_fields
        else:
            raise OptionError(
                f"The output fields {outfield_string} were provided in the wrong format {type(output_fields)}. Please provide a string or a list."
            )
        outfield_string = "LigName, " + ", ".join(
            [self.field_to_column_name[field] for field in output_fields_list]
        )

        query = (
            "SELECT "
            + outfield_string
            + f" FROM Results WHERE Pose_ID IN (SELECT Pose_ID FROM {self.bookmark_name})"
        )
        return self._run_query(query)

    def fetch_flexres_info(self):
        """fetch flexible residues names and atomname lists

        Returns:
            tuple: (flexible_residues, flexres_atomnames)
        """
        try:
            cur = self.conn.cursor()
            cur.execute("SELECT flexible_residues, flexres_atomnames FROM Receptors")
            info = cur.fetchone()
            cur.close()
            return info
        except sqlite3.OperationalError as e:
            raise DatabaseQueryError("Error retrieving flexible residue info") from e

    def fetch_passing_ligand_output_info(self) -> iter:
        """fetch information required by vsmanager for writing out molecules

        Returns:
            iter: contains LigName, ligand_smile,
                atom_index_map, hydrogen_parents
        """
        query = "SELECT LigName, ligand_smile, atom_index_map, hydrogen_parents FROM Ligands WHERE LigName IN (SELECT DISTINCT LigName FROM passing_temp)"
        return self._run_query(query)

    def fetch_single_ligand_output_info(self, ligname) -> str:
        """get output information for given ligand

        Args:
            ligname (str): ligand name

        Raises:
            DatabaseQueryError

        Returns:
            str: information containing smiles, atom and index mapping, and hydrogen parents
        """
        try:
            cur = self.conn.cursor()
            cur.execute(
                f"SELECT LigName, ligand_smile, atom_index_map, hydrogen_parents FROM Ligands WHERE LigName LIKE '{ligname}'"
            )
            info = cur.fetchone()
            cur.close()
            return info
        except sqlite3.OperationalError as e:
            raise DatabaseQueryError(
                f"Error retrieving ligand info for {ligname}"
            ) from e

    def fetch_single_pose_properties(self, pose_ID: int) -> iter:
        """fetch coordinates for pose given by pose_ID

        Args:
            pose_ID (int): name of ligand to fetch coordinates for

        Returns:
            iter: SQLite cursor that contains Pose_ID, docking_score, leff, ligand_coordinates,
                flexible_res_coordinates, flexible_residues
        """
        query = f"SELECT Pose_ID, docking_score, leff, ligand_coordinates, flexible_res_coordinates FROM Results WHERE Pose_ID={pose_ID}"
        return self._run_query(query)

    def fetch_interaction_info_by_index(self, interaction_idx) -> tuple:
        """Returns tuple containing interaction info for given interaction_idx

        Args:
            interaction_idx (int): interaction index to fetch info for

        Returns:
            tuple: tuple of info for requested interaction
        """
        query = "SELECT * FROM Interaction_indices WHERE interaction_id = {0}".format(
            interaction_idx
        )
        return self._run_query(query).fetchone()[1:]  # cut off interaction index

    def fetch_pose_interactions(self, Pose_ID) -> iter:
        """
        Fetch all interactions parameters belonging to a Pose_ID

        Args:
            Pose_ID (int): pose id, 1-1 with Results table

        Returns:
            iter: of interaction information for given Pose_ID
        """
        # check if table exist
        cur = self._run_query(
            """SELECT name FROM sqlite_master WHERE type='table' AND name='Interactions';"""
        )
        if len(cur.fetchall()) == 0:
            return None

        query = f"""SELECT ii.interaction_type, ii.rec_chain, ii.rec_resname, ii.rec_resid, ii.rec_atom, ii.rec_atomid 
        FROM Interaction_indices ii 
        JOIN Interactions i ON i.interaction_id = ii.interaction_id
        WHERE i.Pose_ID = {Pose_ID}"""

        return self._run_query(query).fetchall()

    def count_receptors_in_db(self):
        """returns number of rows in Receptors table where receptor_object already has blob

        Returns:
            int: number of rows in receptors table
            str: name of receptor if present in table

        Raises:
            DatabaseQueryError
        """
        try:
            cur = self.conn.execute(
                "SELECT COUNT(*) FROM Receptors WHERE receptor_object NOT NULL"
            )
            row_count = cur.fetchone()[0]
            cur.close()
            return row_count
        except sqlite3.OperationalError as e:
            raise DatabaseQueryError(
                "Error occurred while fetching number of receptor rows containing PDBQT blob"
            ) from e

    def _fetch_all_plot_data(self):
        """Fetches cursor for best energies and leff for all ligands

        Returns:
             iter: SQLite Cursor containing docking_score,
                leff for the first pose for each ligand
        """
        return self._run_query(
            "SELECT docking_score, leff FROM Results GROUP BY LigName"
        )

    def _fetch_passing_plot_data(self, bookmark_name: str | None = None):
        """Fetches cursor for best energies and leffs for
            ligands passing filtering

        Args:
            bookmark_name (str): name for bookmark for which to fetch data. None will return data for default bookmark_name

        Returns:
            iter: SQL Cursor containing docking_score,
                leff for the first pose for passing ligands
        """
        if bookmark_name is None:
            bookmark_name = self.bookmark_name

        return self._run_query(
            f"SELECT docking_score, leff, Pose_ID, LigName FROM Results WHERE LigName IN (SELECT DISTINCT LigName FROM {bookmark_name}) GROUP BY LigName"
        )

    def _fetch_ligand_cluster_columns(self):
        """fetching columns from Ligand_clusters table

        Raises:
            IndexError

        Returns:
            list: columns from ligand clusters table
        """
        try:
            return [
                c[1]
                for c in self._run_query(
                    "PRAGMA table_info(Ligand_clusters)"
                ).fetchall()
            ][1:]
        except IndexError:
            raise IndexError(
                "Error fetching columns from Ligand_clusters table. Confirm that ligand clustering has been previously performed."
            )

    def _fetch_results_column_names(self):
        """Fetches list of string for column names in results table

        Returns:
            list: List of strings of results table column names

        Raises:
            StorageError
        """
        try:
            return [
                column_tuple[1]
                for column_tuple in self.conn.execute("PRAGMA table_info(Results)")
            ]
        except sqlite3.OperationalError as e:
            raise StorageError(
                "Error while fetching column names from Results table"
            ) from e

    def to_dataframe(self, requested_data: str, table=True) -> pd.DataFrame:
        """Returns a panda dataframe of table or query given as requested_data

        Args:
            requested_data (str): String containing SQL-formatted query or table name
            table (bool): Flag indicating if requested_data is table name or not

        Returns:
            pd.DataFrame: dataframe of requested data
        """
        if table:
            return pd.read_sql_query(
                "SELECT * FROM {0}".format(requested_data), self.conn
            )
        else:
            return pd.read_sql_query(requested_data, self.conn)

    def _get_length_of_table(self, table_name: str):
        """
        Finds the rowcount/length of a table based on the rowid

        Args:
            table_name (str): name of table to count the length of

        Returns:
            int: length of the table
        """
        query = """SELECT COUNT(rowid) from {0}""".format(table_name)

        return self._run_query(query).fetchone()[0]

    def fetch_summary_data(
        self, columns=["docking_score", "leff"], percentiles=[1, 10]
    ) -> dict:
        """Collect summary data for database:
            Num Ligands
            Num stored poses
            Num unique interactions

            min, max, percentiles for columns in columns

        Args:
            columns (list (str)): columns to be displayed and used in summary
            percentiles (list(int)): percentiles to consider

        Returns:
            dict: of data summary
        """
        try:
            summary_data = {}
            cur = self.conn.cursor()
            summary_data["num_ligands"] = cur.execute(
                "SELECT COUNT(LigName) FROM Ligands"
            ).fetchone()[0]
            if summary_data["num_ligands"] == 0:
                raise StorageError("There is no ligand data in the database. ")
            summary_data["num_poses"] = cur.execute(
                "SELECT COUNT(Pose_id) FROM Results"
            ).fetchone()[0]
            summary_data["num_unique_interactions"] = cur.execute(
                "SELECT COUNT(interaction_id) FROM Interaction_indices"
            ).fetchone()[0]
            summary_data["num_interacting_residues"] = cur.execute(
                "SELECT COUNT(*) FROM (SELECT interaction_id FROM Interaction_indices GROUP BY interaction_type,rec_resid,rec_chain)"
            ).fetchone()[0]

            allowed_columns = self._fetch_results_column_names()
            for col in columns:
                if col not in allowed_columns:
                    raise StorageError(
                        f"Requested summary column {col} not found in Results table! Available columns: {allowed_columns}"
                    )
                summary_data[f"min_{col}"] = cur.execute(
                    f"SELECT MIN({col}) FROM Results"
                ).fetchone()[0]
                summary_data[f"max_{col}"] = cur.execute(
                    f"SELECT MAX({col}) FROM Results"
                ).fetchone()[0]
                for p in percentiles:
                    summary_data[f"{p}%_{col}"] = self._calc_percentile_cutoff(p, col)

            return summary_data

        except sqlite3.OperationalError as e:
            raise StorageError("Error while fetching summary data!") from e

    # endregion

    # region Methods dealing with filtered results

    def _get_number_passing_ligands(self, bookmark_name: str | None = None):
        """Returns count of the number of ligands that
            passed filtering criteria

        Args:
            bookmark_name (str): bookmark name to query

        Returns:
            int: Number of passing ligands

        Raises:
            DatabaseQueryError
        """
        if bookmark_name is None:
            bookmark_name = self.current_bookmark_name
        try:
            cur = self.conn.cursor()
            cur.execute(f"SELECT COUNT(DISTINCT LigName) FROM {bookmark_name}")
            n_ligands = int(cur.fetchone()[0])
            cur.close()
            return n_ligands
        except sqlite3.OperationalError as e:
            raise DatabaseQueryError(
                "Error while getting number of passing ligands"
            ) from e

    def get_maxmiss_union(self, total_combinations: int):
        """Get results that are in union considering max miss

        Args:
            total_combinations (int): numer of possible combinations

        Returns:
            iter: of passing results
        """
        selection_strs = []
        view_strs = []
        outfield_list = self._generate_outfield_list()
        for i in range(total_combinations):
            selection_strs.append(
                f"""SELECT {", ".join(outfield_list)} FROM {self.bookmark_name + '_' + str(i)}"""
            )
            view_strs.append(f"SELECT * FROM {self.bookmark_name + '_' + str(i)}")

        bookmark_name = f"{self.bookmark_name}_union"
        self.logger.debug("Saving union bookmark...")
        union_view_query = " UNION ".join(view_strs)
        union_select_query = " UNION ".join(selection_strs)
        if not self.output_all_poses:
            # if not outputting all poses, it is necessary to "create" the view (each of which had a grouping statement), then group by in the final view
            union_view_query = (
                "SELECT * FROM (" + union_view_query + ") GROUP BY LigName"
            )
            union_select_query = (
                "SELECT * FROM (" + union_select_query + ") GROUP BY LigName"
            )
        self.create_bookmark(bookmark_name, union_view_query)
        self.logger.debug("Running union query...")
        return self._run_query(union_select_query)

    def fetch_clustered_similars(self, ligname: str):
        """Given ligname, returns poseids for similar poses/ligands from previous clustering. User prompted at runtime to choose cluster.

        Args:
            ligname (str): ligname for ligand to find similarity with

        Raises:
            ValueError: wrong terminal input
            DatabaseQueryError
        """
        self.logger.warning(
            "N.B.: When finding similar ligands, export tasks (i.e. SDF export) will be for the selected similar ligands, NOT ligands passing given filters."
        )
        cur = self.conn.cursor()

        ligand_cluster_columns = self._fetch_ligand_cluster_columns()
        print(
            "Here are the existing clustering groups. Please ensure that you query ligand(s) is a part of the group you select."
        )
        print(
            "   Choice number   |   Underlying filter bookmark   |   Morgan or interaction fingerprint?   |   cutoff   "
        )
        print(
            "----------------------------------------------------------------------------------------------------------"
        )
        for i, col in enumerate(ligand_cluster_columns):
            col_info = col.split("_")
            option_list = (
                [str(i)]
                + ["_".join(col_info[:-2])]
                + [col_info[-2]]
                + [col_info[-1].replace("p", ".")]
            )
            print(f"{'    |    '.join(option_list)}")
        cluster_choice = input(
            "Please specify choice number for the cluster you would like to return similar ligands from: "
        )
        try:
            cluster_col_choice = ligand_cluster_columns[int(cluster_choice)]
        except ValueError:
            raise ValueError(
                f"Given cluster number {cluster_choice} cannot be converted to int. Please be sure you are specifying integer."
            )
        query_ligand_cluster = cur.execute(
            f"SELECT {cluster_col_choice} FROM Ligand_clusters WHERE pose_id IN (SELECT Pose_ID FROM Results WHERE LigName LIKE '{ligname}')"
        ).fetchone()
        if query_ligand_cluster is None:
            raise DatabaseQueryError(
                f"Requested ligand name {ligname} not found in cluster {cluster_col_choice}!"
            )
        query_ligand_cluster = query_ligand_cluster[0]  # extract from tuple
        sql_query = f"SELECT LigName FROM Results WHERE Pose_ID IN (SELECT pose_id FROM Ligand_clusters WHERE {cluster_col_choice}={query_ligand_cluster}) GROUP BY LigName"
        view_query = f"SELECT * FROM Results WHERE Pose_ID IN (SELECT pose_id FROM Ligand_clusters WHERE {cluster_col_choice}={query_ligand_cluster}) GROUP BY LigName"

        self.bookmark_name = f"similar_{ligname}_{cluster_col_choice}"
        self.create_bookmark(self.bookmark_name, view_query)

        return self._run_query(sql_query), self.bookmark_name, cluster_col_choice

    def fetch_passing_pose_properties(self, ligname):
        """fetch coordinates for poses passing filter for given ligand

        Args:
            ligname (str): name of ligand to fetch coordinates for

        Returns:
            iter: SQLite cursor that contains Pose_ID, docking_score, leff, ligand_coordinates,
                flexible_res_coordinates, flexible_residues
        """
        query = f"SELECT Pose_ID, docking_score, leff, ligand_coordinates, flexible_res_coordinates FROM Results WHERE Pose_ID IN (SELECT Pose_ID FROM passing_temp WHERE LigName LIKE '{ligname}')"
        return self._run_query(query)

    def fetch_nonpassing_pose_properties(self, ligname):
        """fetch coordinates for poses of ligname which did not pass the filter

        Args:
            ligname (str): name of ligand to fetch coordinates for

        Returns:
            iter: SQLite cursor that contains Pose_ID, docking_score, leff, ligand_coordinates,
                flexible_res_coordinates, flexible_residues
        """
        query = f"SELECT Pose_ID, docking_score, leff, ligand_coordinates, flexible_res_coordinates FROM Results WHERE LigName LIKE '{ligname}' AND Pose_ID NOT IN (SELECT Pose_ID FROM passing_temp)"
        return self._run_query(query)

    def _calc_percentile_cutoff(self, percentile: float, column="docking_score"):
        """Make query for percentile by calculating energy or leff cutoff

        Args:
            percentile (float): cutoff percentile
            column (str, optional): string indicating column for percentile to be calculated over

        Returns:
            float: effective cutoff value of results based on percentile
        """
        # get total number of ligands
        try:
            self.logger.debug(f"Generating percentile filter query for {column}")
            cur = self.conn.cursor()
            cur.execute("SELECT COUNT(LigName) FROM Ligands")
            n_ligands = int(cur.fetchone()[0])
            n_passing = int((percentile / 100) * n_ligands)
            # find energy cutoff
            counter = 0
            for i in cur.execute(
                f"SELECT {column} FROM Results GROUP BY LigName ORDER BY {column}"
            ):
                if counter == n_passing:
                    cutoff = i[0]
                    break
                counter += 1
            self.logger.debug(f"{column} percentile cutoff is {cutoff}")
            return cutoff
        except sqlite3.OperationalError as e:
            raise StorageError("Error while generating percentile query") from e

    def _generate_outfield_list(self):
        """list describing outfields to be written

        Returns:
            list: column names for which the data is to be displayed

        Raises:
            OptionError
        """
        # parse requested output fields and convert to column names in database
        outfields_list = self.outfields.split(",")
        for outfield in outfields_list:
            if outfield not in self._data_kw_groups("outfield_options"):
                raise OptionError(
                    f"{outfield} is not a valid output option. Please see rt_process_vs.py --help for allowed options"
                )
        return [self.field_to_column_name[field] for field in outfields_list]

    def _process_filters_for_query(self, filters_dict: dict):
        # NOTE this method can maybe be a main class method once we get more database types
        """
        Method that reformats the filters to the specified database columns, handles less than/more than filters, etc

        Args:
            filters_dict (dict): all Ringtail filters, okay to contain None

        Returns:
            list: list of numerical filters formatted to be inserted in a query
            list: list of interaction filters formatted to be inserted in a query
        """
        # write energy filters and compile list of interactions to search for
        numerical_filters = []
        interaction_filters = []
        ligand_filters = {}
        energy_filter_col_name = {
            "eworst": "docking_score",
            "ebest": "docking_score",
            "leworst": "leff",
            "lebest": "leff",
            "score_percentile": "docking_score",
            "le_percentile": "leff",
        }
        for filter_key, filter_value in filters_dict.items():
            # filter dict contains all possible filters, are None if not specified by user
            if filter_value is None:
                continue
            # if filter has to do with docking energies
            if filter_key in energy_filter_col_name:
                if filter_key == "score_percentile" or filter_key == "le_percentile":
                    # convert from percent to decimal
                    cutoff = self._calc_percentile_cutoff(
                        filter_value, energy_filter_col_name[filter_key]
                    )
                    numerical_filters.append(
                        f"{energy_filter_col_name[filter_key]} < {cutoff}"
                    )
                else:
                    numerical_filters.append(
                        self.energy_filter_sqlite_call_dict[filter_key].format(
                            value=filter_value
                        )
                    )

            # write hb count filter(s)
            if filter_key == "hb_count":
                for k, v in filter_value:
                    if k != "hb_count":
                        self.logger.warning(
                            f"An unrecognized interaction count filter was found: {k}, which will not be included in the filtering."
                        )
                        continue
                    if v > 0:
                        numerical_filters.append(f"num_hb > {v}")
                    else:
                        # if value is negative, it means less than specified number of hydrogen bonds
                        numerical_filters.append(f"num_hb <= {-v}")
            interaction_name_to_letter = {
                "vdw_interactions": "V",
                "hb_interactions": "H",
                "reactive_interactions": "R",
            }
            # reformat interaction filters as list
            if filter_key in Filters.get_filter_keys("interaction"):
                for interact in filter_value:
                    # interact has format ["chain:res:resno:resatom", bool(include or exclude interaction)]
                    interaction_string = (
                        interaction_name_to_letter[filter_key] + ":" + interact[0]
                    )
                    # add bool flag for included (T) or excluded (F) interaction
                    interaction_filters.append(
                        interaction_string.split(":") + [interact[1]]
                    )
            # add react_any flag as interaction filter if not None
            if filter_key == "react_any" and filter_value:
                interaction_filters.append(["R", "", "", "", "", True])

            # if filter has to do with ligands and SMARTS
            if filter_key in Filters.get_filter_keys("ligand"):
                if filter_key == "ligand_substruct_pos" and filter_value:
                    # go through each item and make sure the numbers are cast from string to numbers
                    for filter in filter_value:
                        # cast second item to int
                        filter[1] = int(filter[1])
                        # cast last four items to float
                        for index in range(2, 6):
                            filter[index] = float(filter[index])

                ligand_filters[filter_key] = filter_value

            if filter_key == "max_miss":
                max_miss = filter_value
        # put all processed filter in a dict
        processed_filters = {}
        if len(numerical_filters) > 0:
            processed_filters["num_filters"] = numerical_filters
        if len(interaction_filters) > 0:
            processed_filters["int_filters"] = interaction_filters
            processed_filters["max_miss"] = max_miss
        if len(ligand_filters) > 0:
            processed_filters["lig_filters"] = ligand_filters
        return processed_filters

    def _generate_result_filtering_query(self, filters_dict):
        """takes lists of filters, writes sql filtering string

        Args:
            filters_dict (dict): dict of filters. Keys names and value formats must match those found in the Filters class

        Returns:
            str: SQLite-formatted string for filtering query
        """
        # table to filter over
        filtering_window = "Results"

        outfield_columns = self._generate_outfield_list()
        num_query = ""
        int_query = ""
        lig_query = ""
        ligand_substruct_queries = []
        join_stmnt = ""

        # if filtering over a bookmark (i.e., already filtered results) as opposed to a whole database
        if self.filter_bookmark is not None:
            if self.filter_bookmark == self.bookmark_name:
                # cannot write data from bookmark_a to bookmark_a
                self.logger.error(
                    f"Specified 'filter_bookmark' and 'bookmark_name' are the same: {self.bookmark_name}"
                )
                raise OptionError(
                    "'filter_bookmark' and 'bookmark_name' cannot be the same! Please rename 'bookmark_name'"
                )
            # cannot use percentile for an already reduced dataset
            if (
                filters_dict["score_percentile"] is not None
                or filters_dict["le_percentile"] is not None
            ):
                raise OptionError(
                    "Cannot use 'score_percentile' or 'le_percentile' with 'filter_bookmark'."
                )
            # filtering window can be specified bookmark, as opposed to entire database using Results table
            filtering_window = self.filter_bookmark

        # process filter values to lists and dicts that are easily incorporated in sql queries
        processed_filters = self._process_filters_for_query(filters_dict)

        # raise error if no filters are present and no clusterings
        if not processed_filters and not clustering:
            raise DatabaseQueryError(
                "Ringtail query strings are empty, please check filter options."
            )

        # check if clustering
        clustering = bool(self.mfpt_cluster or self.interaction_cluster)
        # if clustering without filtering
        if clustering:
            # allows for clustering without filtering
            self.logger.info("Preparing to cluster results")
            unclustered_query = f"SELECT R.Pose_id FROM {filtering_window} R "
            if not processed_filters and filtering_window == "Results":
                self.logger.warning(
                    "If clustering is not performed on a pre-filtered bookmark, the clustering process can be slow."
                )
        else:
            # start with empty string, will prepend SELECT statement later
            unclustered_query = ""

        # create query string from filters if present
        if processed_filters:
            # start stringing together queries
            # check what filters are present, and prepare them as partial queries
            if "num_filters" in processed_filters:
                num_query = " AND ".join(
                    ["R." + filter for filter in processed_filters["num_filters"]]
                )
            # check for interactions and prepare for query
            if "int_filters" in processed_filters:
                # if interaction filters are present and valid, two lists of included and excluded interactions are returned
                # each item in the lists to be joined by "AND", and each item within the list item (if >1) to be joined by "OR"
                include_interactions, exclude_interactions = (
                    self._prepare_interaction_indices_for_filtering(
                        interaction_list=processed_filters["int_filters"]
                    )
                )
                # ensure there are interactions in the list after processing
                if bool(exclude_interactions or include_interactions):
                    # prepare partial queries for the different interaction combinations
                    int_query = self._prepare_interaction_filtering_query(
                        include_interactions,
                        exclude_interactions,
                        processed_filters["max_miss"],
                    )
            # check if ligand filters and prepare for query
            if "lig_filters" in processed_filters:
                lig_filters = processed_filters["lig_filters"]
                ligand_queries = []
                if (
                    "ligand_substruct" in lig_filters
                    or "ligand_name" in lig_filters
                    or "ligand_max_atoms" in lig_filters
                ):
                    ligand_queries.append(
                        self._generate_ligand_filtering_query(lig_filters)
                    )
                # if complex ligand filter, generate partial query
                if "ligand_substruct_pos" in lig_filters:
                    for substruct_pos in lig_filters["ligand_substruct_pos"]:
                        temp_lig_filter = lig_filters
                        temp_lig_filter["ligand_substruct_pos"] = substruct_pos
                        ligand_substruct_queries.append(
                            self._ligand_substructure_position_filter(temp_lig_filter)
                        )
                    join_stmnt = " " + lig_filters["ligand_operator"] + " "
                # join all ligand queries that are not empty
                lig_query = " AND ".join(
                    [lig_filter for lig_filter in ligand_queries if lig_filter]
                )
            # if filter queries exist for each group, string them together appropriately
            if int_query:
                # add with a join statement
                unclustered_query += (
                    "JOIN (" + int_query + ") I ON R.Pose_ID = I.Pose_ID "
                )
            if lig_query:
                # add with a join statement
                unclustered_query += (
                    "JOIN (" + lig_query + ") L ON R.LigName = L.LigName "
                )
                # these two queries are joined on the Results table, after the multiple table spanning queries
            if num_query or ligand_substruct_queries:
                # add condition
                unclustered_query += "WHERE "
                # add numerical part of query
                if num_query:
                    unclustered_query += num_query
                # if both numerical and ligand_substruct_pos handle appropriately
                if num_query and ligand_substruct_queries:
                    unclustered_query += " AND " + join_stmnt.join(
                        ligand_substruct_queries
                    )
                # if not, only the ligand_substruct_pos sets the WHERE condition
                else:
                    unclustered_query += join_stmnt.join(ligand_substruct_queries)

        # if clustering is requested, do that before saving view or filtering results for output
        if clustering:
            # add appropriate select
            try:
                query = self._prepare_cluster_query(unclustered_query)
                query = "WHERE " + query
            except OptionError as e:
                raise e
        else:
            # if not clustering, rename query
            query = unclustered_query

        # choose columns to be selected from filtering_window
        query_select_string = f"""SELECT {", ".join("R." + column for column in outfield_columns)} FROM {filtering_window} R """
        # adding if we only want to keep one pose per ligand (will keep first entry)
        if not self.output_all_poses:
            query += " GROUP BY R.LigName "
        # add how to order results
        if self.order_results:
            query += " ORDER BY " + self.field_to_column_name[self.order_results]

        output_query = query_select_string + query
        view_query = f"SELECT * FROM {filtering_window} R " + query
        return output_query, view_query

    def _prepare_cluster_query(self, unclustered_query: str) -> str | None:
        """
        These methods will take data returned from unclustered filter query, then run the cluster query and cluster the filtered data.
        This will output pose_ids that are representative of the clusters, and these pose_ids will be returned so that
        they can be added to the unclustered query in the main filtering method.

        Args:
            unclustered_query (str): query containing none or some filters that defines over which ligands the clustering should happen

        Returns:
            str: (reduced) query to include in overall filter query if clustering returned results
        """
        if self.interaction_cluster and self.mfpt_cluster:
            self.logger.warning(
                "N.B.: If using both interaction and morgan fingerprint clustering, the morgan fingerprint clustering will be performed on the results staus post interaction fingerprint clustering."
            )

        def _clusterFps(fps, cutoff):
            """
            https://macinchem.org/2023/03/05/options-for-clustering-large-datasets-of-molecules/

            Args:
                fps (): fingerprints
                cutoff distance (float)
            """

            # first generate the distance matrix:
            dists = []
            nfps = len(fps)
            inputs = []

            def gen(fps):
                for i in range(1, len(fps)):
                    yield (i, fps)

            def mp_wrapper(input_tpl):
                i, fps = input_tpl
                return DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])

            with multiprocess.Pool() as p:
                inputs = gen(fps)
                for sims in p.imap(mp_wrapper, inputs):
                    dists.extend([1 - x for x in sims])

            # now cluster the data:
            cs = Butina.ClusterData(dists, nfps, cutoff, isDistData=True)
            return cs

        cluster_query_string = None

        if self.interaction_cluster:
            cluster_query = f"SELECT Pose_ID, leff FROM Results WHERE Pose_ID IN ({unclustered_query})"
            # resulting data
            poseid_leffs = self._run_query(cluster_query).fetchall()
            cluster_poseids = (
                "("
                + ",".join(map(str, [(poseid_leff[0]) for poseid_leff in poseid_leffs]))
                + ")"
            )
            poseid_bvs = self._generate_interaction_bitvectors(cluster_poseids)

            # create a list of tuples from the query data and bitvector
            poseid_leff_bvs = [
                (poseid_leff[0], poseid_leff[1], poseid_bvs[str(poseid_leff[0])])
                for poseid_leff in poseid_leffs
            ]
            # index 2 is the bitvector string element
            bclusters = _clusterFps(
                [
                    DataStructs.CreateFromBitString(poseid_leff_bv[2])
                    for poseid_leff_bv in poseid_leff_bvs
                ],
                self.interaction_cluster,
            )
            self.logger.info(
                f"Number of interaction fingerprint butina clusters: {len(bclusters)}"
            )

            # select ligand from each cluster with best ligand efficiency
            # interaction clusters representative pose ids
            int_rep_poseids = []

            for cluster in bclusters:
                # element 1 in individual pose id item is the ligand efficiency (leff)
                c_leffs = np.array(
                    [poseid_leff_bvs[cluster_element][1] for cluster_element in cluster]
                )
                # element 0 ([0]) in each leff_poseid_ifps row is the pose_id
                best_lig_c = poseid_leff_bvs[cluster[np.argmin(c_leffs)]][0]
                int_rep_poseids.append(str(best_lig_c))

            # element 0 ([0]) in each leff_poseid_ifps row is the pose_id
            self._insert_cluster_data(
                bclusters,
                [l[0] for l in poseid_leff_bvs],
                "ifp",
                str(self.interaction_cluster),
            )

            # catch if no pose_ids returned
            if int_rep_poseids == []:
                raise OptionError(
                    "No passing results prior to clustering. Clustering not performed."
                )
            else:
                cluster_query_string = "R.Pose_ID = " + " OR R.Pose_ID = ".join(
                    int_rep_poseids
                )
                # if more clustering
                if self.mfpt_cluster is not None:
                    # carry the pose ids returned by this cluster to the MFPT clustering
                    unclustered_query = (
                        f"SELECT R.Pose_ID FROM Results WHERE {cluster_query_string}"
                    )

        if self.mfpt_cluster:
            cluster_query = f"SELECT R.Pose_ID, R.leff, mol_morgan_bfp(L.ligand_rdmol, 2, 1024) FROM Ligands L INNER JOIN Results R ON R.LigName = L.LigName WHERE R.Pose_ID IN ({unclustered_query})"
            poseid_leff_mfps = self._run_query(cluster_query).fetchall()
            bclusters = _clusterFps(
                [DataStructs.CreateFromBinaryText(mol[2]) for mol in poseid_leff_mfps],
                self.mfpt_cluster,
            )
            self.logger.info(
                f"Number of Morgan fingerprint butina clusters: {len(bclusters)}"
            )
            # select ligand from each cluster with best ligand efficiency
            fp_rep_poseids = []
            for c in bclusters:
                c_leffs = np.array([poseid_leff_mfps[i][1] for i in c])
                best_lig_c = poseid_leff_mfps[c[np.argmin(c_leffs)]][0]
                fp_rep_poseids.append(str(best_lig_c))

            self._insert_cluster_data(
                bclusters,
                [l[0] for l in poseid_leff_mfps],
                "mfp",
                str(self.mfpt_cluster),
            )

            # catch if no pose_ids returned
            if fp_rep_poseids == []:
                raise OptionError(
                    "No passing results prior to clustering. Clustering not performed."
                )
            else:
                cluster_query_string = "R.Pose_ID = " + " OR R.Pose_ID = ".join(
                    fp_rep_poseids
                )
        return cluster_query_string

    def _prepare_interaction_filtering_query(
        self, include_interactions: list, exclude_interactions: list, max_miss: int
    ) -> str:
        """
        Method that prepares a partial query for interactions

        Args:
            include_interactions (list): interactions a pose should have
            exclude_interactions (list): interactions a pose should not have
            max_miss (int): max number of the provided interactions a pose_id is allowed to miss

        Returns:
            str: partial query to include in main filter query
        """
        # nonsensical number to count an interaction if it satisfies an incomplete ("wildcard") interaction
        nonsense_counter = -10000
        num_of_interactions = (
            len(include_interactions) + len(exclude_interactions) - max_miss
        )

        def _prepare_indices_for_query(interactions: list):
            """
            Method to organize a list of indices into those to be included in an OR and an AND statement

            Args:
                interactions (list): list of indices organized by how they were queried from the db

            Returns:
                list, list: indices organized in lists appropriate for the query
            """
            and_interactions = []
            or_interactions = []
            # add the included OR indices
            for index_list in interactions:
                # one list per interaction, mode indices if interaction had a wildcard
                for index in index_list:
                    # if index not already represented
                    if index not in and_interactions:
                        # add to list
                        and_interactions.append(index)
                # if index list has more than one element, they should also be combined in an "OR" statement
                if len(index_list) > 1:
                    # adds a string element to the list
                    or_interactions.append(index_list)
            return and_interactions, or_interactions

        # prepare lists of indices ready to be cast to tuples and strings
        if include_interactions:
            and_include_interactions, or_include_interactions = (
                _prepare_indices_for_query(include_interactions)
            )
        else:
            and_include_interactions = []
            or_include_interactions = []
        if exclude_interactions:
            and_exclude_interactions, or_exclude_interactions = (
                _prepare_indices_for_query(exclude_interactions)
            )
        else:
            and_exclude_interactions = []
            or_exclude_interactions = []

        # building the query
        # 1. select pose id, call CASE, in paranthesis because grouping with different query
        query = "SELECT Pose_ID FROM (SELECT Pose_ID "
        if or_include_interactions or or_exclude_interactions:
            # add the case statements
            query += ", CASE "
            # 2. list all OR statements
            if or_include_interactions:
                for interactions in or_include_interactions:
                    # iterate the nonsense counter
                    nonsense_counter += 1
                    query += (
                        "WHEN interaction_id IN ("
                        + numlist2str(interactions, ",")
                        + f") THEN {nonsense_counter} "
                    )
            if or_exclude_interactions:
                for interactions in or_exclude_interactions:
                    # iterate the nonsense counter
                    nonsense_counter += 1
                    query += (
                        "WHEN interaction_id NOT IN ("
                        + numlist2str(interactions, ",")
                        + f") THEN {nonsense_counter} "
                    )
            query += "ELSE interaction_id END "
        else:
            query += ", interaction_id "
        # 3. proceed with all interactions
        query += "AS filtered_interactions FROM Interactions WHERE "
        if and_include_interactions:
            query += (
                "interaction_id IN ("
                + numlist2str(and_include_interactions, ",")
                + ") "
            )
            # if both include and exclude are there, need "AND"
            if and_exclude_interactions:
                query += "AND "
        if and_exclude_interactions:
            query += (
                "interaction_id NOT IN ("
                + numlist2str(and_exclude_interactions, ",")
                + ") "
            )
        # 4. add grouping and wildcard for total interactions minus max_miss, essentially
        query += f") GROUP BY Pose_ID HAVING COUNT(DISTINCT filtered_interactions) >= ({num_of_interactions}) "

        return query

    def _ligand_substructure_position_filter(self, ligand_filters_dict: dict) -> str:
        """
        Method that takes all ligand filters in the presence of a ligand_substruct_pos filter, and reduces the query to
        " IN pose_ids" based on what pose_ids passed the ligand filters

        Args:
            ligand_filters_dict (dict): all specified ligand filters

        Raises:
            OptionError

        Returns:
            str: partial query that identifies pose ids passing the ligand substructure filter
        """
        queries = []
        nr_args_per_group = 6
        nr_smarts = int(
            len(ligand_filters_dict["ligand_substruct_pos"]) / nr_args_per_group
        )
        # create temporary table with molecules that pass all smiles
        tmp_lig_filters = {"ligand_operator": ligand_filters_dict["ligand_operator"]}
        if "ligand_max_atoms" in ligand_filters_dict:
            tmp_lig_filters["ligand_max_atoms"] = ligand_filters_dict[
                "ligand_max_atoms"
            ]
        tmp_lig_filters["ligand_substruct"] = [
            ligand_filters_dict["ligand_substruct_pos"][i * nr_args_per_group]
            for i in range(nr_smarts)
        ]
        cmd = self._generate_ligand_filtering_query(tmp_lig_filters)
        cmd = cmd.replace(
            "SELECT L.LigName FROM Ligands L",
            "SELECT "
            "R.Pose_ID, "
            "L.LigName, "
            "L.ligand_smile, "
            "L.atom_index_map, "
            "R.ligand_coordinates "
            "FROM Ligands L INNER JOIN Results R ON R.LigName = L.LigName",
        )
        cmd = "CREATE TEMP TABLE passed_smarts AS " + cmd
        cur = self.conn.cursor()
        cur.execute("DROP TABLE IF EXISTS passed_smarts")
        cur.execute(cmd)
        smarts_loc_filters = []
        for i in range(nr_smarts):
            smarts = ligand_filters_dict["ligand_substruct_pos"][
                i * nr_args_per_group + 0
            ]
            index = int(
                ligand_filters_dict["ligand_substruct_pos"][i * nr_args_per_group + 1]
            )
            sqdist = (
                float(
                    ligand_filters_dict["ligand_substruct_pos"][
                        i * nr_args_per_group + 2
                    ]
                )
                ** 2
            )
            x = float(
                ligand_filters_dict["ligand_substruct_pos"][i * nr_args_per_group + 3]
            )
            y = float(
                ligand_filters_dict["ligand_substruct_pos"][i * nr_args_per_group + 4]
            )
            z = float(
                ligand_filters_dict["ligand_substruct_pos"][i * nr_args_per_group + 5]
            )
            # save filter for bookmark
            smarts_loc_filters.append((smarts, index, x, y, z))
            poses = self._run_query("SELECT * FROM passed_smarts")
            pose_id_list = []
            smartsmol = Chem.MolFromSmarts(smarts)
            for pose_id, ligname, smiles, idxmap, coords in poses:
                mol = Chem.MolFromSmiles(smiles)
                idxmap = [int(value) - 1 for value in json.loads(idxmap)]
                idxmap = {
                    idxmap[j * 2]: idxmap[j * 2 + 1]
                    for j in range(int(len(idxmap) / 2))
                }
                for hit in mol.GetSubstructMatches(smartsmol):
                    xyz = [
                        float(value) for value in json.loads(coords)[idxmap[hit[index]]]
                    ]
                    d2 = (xyz[0] - x) ** 2 + (xyz[1] - y) ** 2 + (xyz[2] - z) ** 2
                    if d2 <= sqdist:
                        pose_id_list.append(str(pose_id))
                        break  # add pose only once
            if len(pose_id_list) > 0:
                queries.append("R.Pose_ID IN ({0})".format(",".join(pose_id_list)))
        cur.close()
        if not queries:
            raise OptionError(
                "There are no ligands passing the 'ligand_substruct_pos' filter, please revise your filter query."
            )

        return "".join(queries)

    def _generate_interaction_bitvectors(self, pose_ids: str) -> dict:
        """
        Method to generate a dict of generate bitvector strings from pose_ids

        Args:
            pose_ids (str): query formatted list of pose_ids (as tuple)

        Returns:
            dict: of "pose_id":"bitvector"
        """
        # create a list of 0 items the length of interaction_indices table
        ii_length = self._get_length_of_table("Interaction_indices")
        # for each pose id, get a list of interaction_indices from joining the two tables i and ii
        poseid_intind_query = f"""SELECT Pose_ID, interaction_id
                                    FROM Interactions
                                    WHERE Pose_ID IN {pose_ids}"""
        poseid_intinds = self._run_query(poseid_intind_query).fetchall()
        # cast string to list
        pose_id_list = list((pose_ids.strip("()")).split(","))
        # make dict of pose id and bitvector
        poseid_bvlist = {(pose_id): [0] * ii_length for pose_id in pose_id_list}
        # iterate over the tuple results from the query
        for poseid_intind in poseid_intinds:
            poseid_bvlist[str(poseid_intind[0])][poseid_intind[1] - 1] = 1

        # join list as string without any delimiter
        poseid_bv = {
            key: "".join(map(str, value)) for (key, value) in poseid_bvlist.items()
        }
        # return dict of pose id as string and bitvector
        return poseid_bv

    def _prepare_interaction_indices_for_filtering(self, interaction_list: list):
        """
        Prepare lists of interaction indices where they are grouped by whether or not they should be evaluated as "AND" or "OR",
        and whether to be excluded or included in the passing filter poses

        Args:
            interaction_list (list): list of interactions

        Raises:
            OptionError

        Returns:
            list: two lists of indices for interactions to exclude and to include
        """
        # initialize variables
        exclude_interactions = []
        include_interactions = []
        interaction_not_found = []

        # figure out if each interaction is in database, make a list of list of indices for each interaction
        for interaction in interaction_list:
            # get all interaction indices matching the interaction filter (returns more than one index if filter has a "wildcard")
            interaction_index_tuples = self._get_interaction_indices(interaction[:-1])
            # make list of indices from iterable cursor tuples (should create empty list if no results)
            interaction_indices = [i[0] for i in interaction_index_tuples]
            # catch if interaction not found in database
            if interaction_indices == []:
                if interaction == ["R", "", "", "", "", True]:
                    self.logger.warning(
                        "Given 'react_any' filter, no reactive interactions found. Excluded from filtering."
                    )
                else:
                    # create string representation of ecah interaction not found
                    interaction_not_found.append(":".join(interaction[:4]))
                continue  # ends this iteration of the for loop

            # create a list of lists for interactions to either include or exclude
            if interaction[-1] is True:
                include_interactions.append(interaction_indices)
            elif interaction[-1] is False:
                exclude_interactions.append(interaction_indices)
            else:
                raise OptionError(
                    "Unrecognized flag in interaction. Please contact Forli Lab with traceback and context."
                )
        # if one or more interactions not found, raise error
        if interaction_not_found:
            raise OptionError(
                f"The following interactions do not exist in the database: {interaction_not_found} not found in the database. Please check for spelling errors or remove from filter."
            )
        else:
            return include_interactions, exclude_interactions

    def _get_interaction_indices(self, interaction_list) -> iter:
        """takes list of interaction info and looks up corresponding interaction index

        Args:
            interaction_list (list): List containing interaction info
                in format [<interaction_type>, <rec_chain>, <rec_resname>,
                <rec_resid>, <rec_atom>]

        Returns:
            iter: sqlite cursor with the interaction index/indices
        """
        interaction_info = [
            "interaction_type",
            "rec_chain",
            "rec_resname",
            "rec_resid",
            "rec_atom",
        ]
        len_interaction_info = len(interaction_info)
        sql_string = "SELECT interaction_id FROM Interaction_indices WHERE "

        sql_string += " AND ".join(
            [
                "{column} LIKE '{value}'".format(
                    column=interaction_info[i], value=interaction_list[i]
                )
                for i in range(len_interaction_info)
                if interaction_list[i] != ""
            ]
        )

        return self._run_query(sql_string).fetchall()

    def _generate_ligand_filtering_query(self, ligand_filters: dict) -> str:
        """write string to select from ligand table

        Args:
            ligand_filters (list): List of filters on ligand table

        Returns:
            str: SQLite-formatted query
        """

        sql_ligand_string = "SELECT L.LigName FROM Ligands L WHERE"
        if "ligand_operator" in ligand_filters:
            logical_operator = ligand_filters["ligand_operator"]
        else:
            self.logger.info(
                "A logical operator to combine ligand filters were not provided, will use the default value 'OR'."
            )
            logical_operator = "OR"
        if logical_operator is None:
            logical_operator = "AND"
        for kw in ligand_filters.keys():
            fils = ligand_filters[kw]
            if kw == "ligand_name":
                for name in fils:
                    if name == "":
                        continue
                    name_sql_str = " L.LigName LIKE '%{value}%' OR".format(value=name)
                    sql_ligand_string += name_sql_str
            if kw == "ligand_max_atoms" and ligand_filters[kw] is not None:
                maxatom_sql_str = " mol_num_hvyatms(ligand_rdmol) <= {} {}".format(
                    ligand_filters[kw], logical_operator
                )
                sql_ligand_string += maxatom_sql_str
            if kw == "ligand_substruct":
                for smarts in fils:
                    # check for hydrogens in smarts pattern
                    smarts_mol = Chem.MolFromSmarts(smarts)
                    for atom in smarts_mol.GetAtoms():
                        if atom.GetAtomicNum() == 1:
                            raise DatabaseQueryError(
                                f"Given ligand substructure filter {smarts} contains explicit hydrogens. Please re-run query with SMARTs without hydrogen."
                            )
                    substruct_sql_str = " mol_is_substruct(ligand_rdmol, mol_from_smarts('{smarts}')) {logical_operator}".format(
                        smarts=smarts, logical_operator=logical_operator
                    )
                    sql_ligand_string += substruct_sql_str
        if sql_ligand_string.endswith("AND"):
            sql_ligand_string = sql_ligand_string.rstrip("AND")
        if sql_ligand_string.endswith("OR"):
            sql_ligand_string = sql_ligand_string.rstrip("OR")

        return sql_ligand_string

    def _generate_selective_insert_query(
        self, bookmark1_name, bookmark2_name, select_str, new_db_name, temp_table
    ):
        """Generates string to select ligands found/not found in the given bookmark in both current db and new_db

        Args:
            bookmark1_name (str): name of bookmark to cross-reference for main db
            bookmark2_name (str): name of bookmark to cross-reference for attached db
            select_str (str): "IN" or "NOT IN" indicating if ligand names should or should not be in both databases
            new_db_name (str): name of attached db
            temp_table (str): name of temporary table to store passing results in

        Returns:
            str: sqlite formatted query string
        """
        return f"INSERT INTO {temp_table} SELECT Pose_ID, LigName FROM {bookmark1_name} WHERE LigName {select_str} (SELECT LigName FROM {new_db_name}.{bookmark2_name})"

    # endregion

    # region Database operations
    def _open_storage(self):
        """Create connection to db. Then, check if db needs to be created.
        If self.overwrite drop existing tables and initialize new tables

        Raises:
            StorageError
        """
        try:
            self.conn = self._create_connection()
            signal(
                SIGINT, self._sigint_handler
            )  # signal handler to catch keyboard interupts
            if self._db_empty() or self.overwrite:  # write and drop tables as necessary
                if not self._db_empty():
                    self._drop_existing_tables()
                self._create_tables()
                self._set_ringtail_db_schema_version(self._db_schema_ver)

            self.logger.info(f"Ringtail connected to database {self.db_file}.")
        except Exception as e:
            raise StorageError(f"Errow while creating or connecting to database: {e}.")

    def check_storage_ready(
        self, run_mode: str, docking_mode: str, store_all_poses: bool, max_poses: int
    ):
        """Check that storage is ready before proceeding.

        Args:
            run_mode (str): if ringtail is ran using cmd line interface or api
            docking_mode (str): what docking engine was used to produce results
            store_all_poses (bool): overrwrites max poses
            max_poses (int): max poses to save to db

        Raises:
            StorageError
            OptionError: if database options are not compatible
        """
        count = self.conn.execute("SELECT COUNT (*) FROM DB_properties").fetchone()[0]

        compatible = True
        if count < 1:
            self.logger.info(
                "Adding results to an existing database that is currently empty of docking results."
            )
        else:
            compatibility_string = "The following database properties do not agree with the properties last used for this database: \n"
            try:
                cur = self.conn.execute(
                    "SELECT * FROM DB_properties ORDER BY DB_write_session DESC LIMIT 1"
                )
                (_, last_docking_mode, num_of_poses) = cur.fetchone()
                if docking_mode != last_docking_mode:
                    compatible = False
                    compatibility_string += f"Current docking mode is {docking_mode} but last used docking mode of database is {last_docking_mode}.\n"
                if num_of_poses == "all" != store_all_poses:
                    compatible = False
                    compatibility_string += f"Current number of poses saved is {max_poses} but database was previously set to 'store_all_poses'.\n"
                elif int(num_of_poses) != max_poses:
                    compatible = False
                    compatibility_string += f"Current number of poses saved is {max_poses} but database was previously set to {num_of_poses}."
            except Exception as e:
                raise e
            finally:
                cur.close()

        if not compatible:
            if run_mode == "cmd":
                raise OptionError(compatibility_string)
            elif run_mode == "api":
                self.logger.warning(compatibility_string)

        # write current database properties to database
        if store_all_poses:
            number_of_poses = "all"
        else:
            number_of_poses = str(max_poses)
        self._insert_db_properties(docking_mode, number_of_poses)
        self.logger.info("Storage compatibility has been checked.")

    def clone(self, backup_name=None):
        """Creates a copy of the db

        Args:
            backup_name (str, optional): name of the cloned database
        """
        if backup_name is None:
            backup_name = self.db_file + ".bk"
        bck = sqlite3.connect(backup_name)
        with bck:
            self.conn.backup(bck, pages=1)
        bck.close()

    def _set_ringtail_db_schema_version(self, db_version: str = "2.0.0"):
        """Will check current stoarge manager db schema version and only set if it is compatible with the code base version (i.e., version(ringtail)).

        Raises:
            StorageError: if versions are incompatible
        """
        # check that code base is compatible with db schema version
        code_version = version("ringtail")
        if code_version in self._db_schema_code_compatibility[db_version]:
            rtdb_version = db_version.replace(".", "")
            # if so, proceed to set db schema version
            cur = self.conn.cursor()
            cur.execute(f"PRAGMA user_version = {rtdb_version}")
            self.conn.commit()
            cur.close()
            self.logger.info("Database version set to {0}".format(rtdb_version))
        else:
            raise StorageError(
                f"Code base version {code_version} is not compatible with database schema version {db_version}."
            )

    def check_ringtaildb_version(self):
        """
        Checks the database version and confirms whether the code base is compatible with it

        Returns:
            bool: whether or not db is compatible with the code base
            str: current database version
        """
        cur = self.conn.cursor()
        db_version = str(cur.execute("PRAGMA user_version").fetchone()[0])
        db_schema_ver = ".".join([*db_version])
        if version("ringtail") in self._db_schema_code_compatibility[db_schema_ver]:
            is_compatible = True
            self.logger.debug(
                "Database version {0} is compatible with code base version {1}".format(
                    db_schema_ver, version("ringtail")
                )
            )
        else:
            is_compatible = False
            self.logger.warning(
                "Database version {0} is NOT compatible with code base version {1}".format(
                    db_schema_ver, version("ringtail")
                )
            )
        cur.close()
        return is_compatible, db_version

    def update_database_version(self, new_version, consent=False):
        """method that updates sqlite database schema 1.0.0 or 1.1.0 to 1.1.0 or 2.0.0

        #NOTE: If you created the database with the duplicate handling option, there is a chance of inconsistent behavior of anything involving interactions as
        the Pose_ID was not used as an explicit foreign key in db v1.0.0 and v1.1.0.

        Args:
            consent (bool, optional): variable to ensure consent to update database is explicit

        Returns:
            bool
        """
        # create cursor
        cur = self.conn.cursor()

        # get consent, same for both
        if not consent:
            self.logger.warning(
                "WARNING: All existing bookmarks in database will be dropped during database update!"
            )
            consent = input("Type 'yes' if you wish to continue: ") == "yes"
        if not consent:
            self.logger.critical("Consent not given for database update. Cancelling...")
            sys.exit(1)

        # get views and drop them
        self.logger.info(f"Updating {self.db_file}...")
        views = cur.execute(
            "SELECT name FROM sqlite_master WHERE type = 'view'"
        ).fetchall()
        for v in views:
            cur.execute(f"DROP VIEW IF EXISTS {v[0]}")
        # delete all rows in bookmarks table
        cur.execute("DELETE FROM Bookmarks")

        # if current version is 1.0.0
        if self.check_ringtaildb_version()[1] == "1.0.0":
            # reformat for v1.1.0
            cur.execute(
                "ALTER TABLE Results RENAME COLUMN energies_binding TO docking_score"
            )
            cur.execute("ALTER TABLE Bookmarks ADD COLUMN filters")
            cur.execute(
                "CREATE INDEX IF NOT EXISTS ak_results ON Results(LigName, docking_score, leff, deltas, reference_rmsd, energies_inter, energies_vdw, energies_electro, energies_intra, nr_interactions, run_number, pose_rank, num_hb)"
            )
            cur.execute(
                "CREATE INDEX IF NOT EXISTS ak_intind ON Interaction_indices(interaction_type, rec_chain, rec_resname, rec_resid, rec_atom, rec_atomid)"
            )
            try:
                self.conn.commit()
                cur.close()
            except sqlite3.OperationalError as e:
                raise DatabaseConnectionError(
                    f"Error while updating database from v1.0.0 to v1.1.0: {e}"
                ) from e
        # if you only wanted to upgrade to v1.1.0, stop here
        if new_version == "1.1.0":
            self._set_ringtail_db_schema_version("1.1.0")  # set explicit version
        elif new_version == "2.0.0":
            # major table updates and sets db version inside method
            self._update_db_110_to_200()

        return consent

    def _update_db_110_to_200(self):
        """
        Method to update from database v 1.1.0 to 2.0.0, will remove bitvetor table and create Interactions table

        Raises:
            DatabaseConnectionError
            StorageError
        """
        # delete interaction table if necessary
        self._delete_table("Interactions")
        # create interaction table
        self._create_interaction_table()

        # get all interaction bitvector tuples
        cur = self.conn.cursor()
        cur.execute("SELECT * FROM Interaction_bitvectors")

        pose_indices = []
        # for each table entry
        for entry in cur:
            # pose id is firste element of tuple
            pose_id = entry[0]
            # enumerate the remaining (1:) tuple data which are all the bits
            for index, bit in enumerate(entry[1:]):
                # if column is "1" it means that (index+1) interaction was active
                if bit == 1:
                    # index will correspond to the Interaction_index table if +1
                    pose_indices.append((pose_id, index + 1))

        try:
            # just populate the Interaction table straight
            cur.executemany(
                """INSERT INTO Interactions (Pose_id, Interaction_id) VALUES (?,?)""",
                pose_indices,
            )
            # drop old bitvector table
            cur.execute("""DROP TABLE IF EXISTS Interaction_bitvectors;""")
            # create new indixes
            cur.execute("CREATE INDEX IF NOT EXISTS ak_poseid ON Results(Pose_id)")
            cur.execute(
                "CREATE INDEX IF NOT EXISTS ak_interactions ON Interactions(Pose_id, interaction_id)"
            )
            cur.execute("CREATE INDEX IF NOT EXISTS ak_ligands ON Ligands(LigName)")
            self.conn.commit()
            self._set_ringtail_db_schema_version("2.0.0")  # set explicit version
        except sqlite3.OperationalError as e:
            raise DatabaseConnectionError(
                f"Error while creating new interaction tables: {e}"
            ) from e
        except StorageError as e:
            raise StorageError(
                f"Error while setting the database schema version: {e}"
            ) from e

    def _create_connection(self):
        """Creates database connection to self.db_file

        Returns:
            SQLite.conn: Connection object to self.db_file

        Raises:
            DatabaseConnectionError
        """
        try:
            con = sqlite3.connect(self.db_file)
            try:
                con.enable_load_extension(True)
                con.load_extension("chemicalite")
                con.enable_load_extension(False)
            except sqlite3.OperationalError as e:
                self.logger.critical(
                    "Failed to load chemicalite cartridge. Please ensure chemicalite is installed with `conda install -c conda-forge chemicalite`."
                )
                raise e
            cursor = con.execute("PRAGMA synchronous = OFF;")
            cursor.execute("PRAGMA journal_mode = MEMORY;")
            con.commit()
            cursor.close()
        except sqlite3.OperationalError as e:
            raise DatabaseConnectionError(
                "Error while establishing database connection"
            ) from e
        return con

    def _close_connection(self):
        """Closes connection to database"""
        self.logger.info("Closing database")
        self.conn.close()

    def _close_open_cursors(self):
        """closes any cursors stored in self.open_cursors.
        Resets self.open_cursors to empty list
        """
        for cur in self.open_cursors:
            cur.close()

        self.open_cursors = []

    def _db_empty(self):
        """empty database, for example if overwrite

        Returns:
            bool: whether or not db is empty
        """
        cur = self.conn.execute(
            "SELECT COUNT(*) name FROM sqlite_master WHERE type='table';"
        )
        tablecount = cur.fetchone()[0]
        cur.close()
        return True if tablecount == 0 else False

    def _vacuum(self):
        """SQLite vacuum rebuilds the database file, repacking it into a minimal amount of disk space

        Raises:
            DatabaseInsertionError
        """
        try:
            cur = self.conn.cursor()
            cur.execute("VACUUM")
            self.conn.commit()
            cur.close()
        except sqlite3.OperationalError as e:
            raise DatabaseInsertionError(f"Error while vacuuming DB") from e

    def _attach_db(self, new_db, new_db_name):
        """Attaches new database file to current database

        Args:
            new_db (str): file name for database to attach
            new_db_name (str): name of new database

        Raises:
            StorageError
        """
        attach_str = f"ATTACH DATABASE '{new_db}' AS {new_db_name}"

        try:
            cur = self.conn.cursor()
            cur.execute(attach_str)
            self.conn.commit()
            cur.close()
        except sqlite3.OperationalError as e:
            raise StorageError(f"Error occurred while attaching {new_db}") from e

    def _detach_db(self, new_db_name):
        """Detaches new database file from current database

        Args:
            new_db_name (str): db name for database to detach

        Raises:
            StorageError
        """
        detach_str = f"DETACH DATABASE {new_db_name}"

        try:
            cur = self.conn.cursor()
            cur.execute(detach_str)
            self.conn.commit()
            cur.close()
        except sqlite3.OperationalError as e:
            raise StorageError(f"Error occurred while detaching {new_db_name}") from e

    def _drop_existing_tables(self):
        """drop any existing tables.
        Will only be called if self.overwrite is true

        Raises:
            StorageError
        """

        # fetch existing tables
        cur = self.conn.cursor()
        tables = self._fetch_existing_table_names()

        # drop tables
        for table in tables:
            # cannot drop this, so we catch it instead
            if table[0] == "sqlite_sequence":
                continue
            try:
                cur.execute("DROP TABLE {table_name}".format(table_name=table[0]))
            except sqlite3.OperationalError as e:
                raise StorageError(
                    "Error occurred while dropping table {0}".format(table[0])
                ) from e
        cur.close()

    def _fetch_existing_table_names(self):
        """Returns list of all tables in database

        Returns:
            list: list of table names

        Raises:
            DatabaseQueryError
        """

        try:
            cur = self.conn.cursor()
            cur.execute("SELECT name FROM sqlite_schema WHERE type='table';")
            return cur.fetchall()
        except sqlite3.OperationalError as e:
            raise DatabaseQueryError(
                "Error while getting names of existing database tables"
            ) from e

    def _run_query(self, query):
        """Executes provided SQLite query. Returns cursor for results.
            Since cursor remains open, added to list of open cursors

        Args:
            query (str): Formated SQLite query as string

        Returns:
            SQLite cursor: Contains results of query
        """
        try:
            cur = self.conn.cursor()
            cur.execute(query)
            self.open_cursors.append(cur)
        except sqlite3.OperationalError as e:
            raise DatabaseQueryError(
                "Unable to execute query {0}: {1}".format(query, e)
            ) from e
        return cur

    # endregion
