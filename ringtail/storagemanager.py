#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail storage adaptors
#

import sqlite3
import time
import json
import os.path
import pandas as pd
from .logutils import LOGGER as logger
import sys
from signal import signal, SIGINT
from rdkit import Chem
from rdkit import DataStructs
from rdkit.ML.Cluster import Butina
import numpy as np
from typing import Union
import time
from importlib.metadata import version
from .ringtailoptions import Filters
from .util import numlist2str
from .exceptions import (
    StorageError,
    DatabaseInsertionError,
    DatabaseConnectionError,
    DatabaseTableCreationError,
    DatabaseQueryError,
    DatabaseViewCreationError,
    OptionError,
    MergeError,
)
from .clustermanager import *
from .querybuilder import QueryBuilder


class StorageManager:
    # NOTE gotta be careful with schema
    _db_schema_ver = "2.0.0"

    # "db_schema_ver":list("compatible code versions")
    _db_schema_code_compatibility = {
        "1.0.0": ["1.0.0"],
        "1.1.0": ["1.1.0"],
        "2.0.0": ["2.0.0", "2.1.0", "2.1.1", "2.1.2"],
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
        self.closed_connection = False
        self.keyboard_interrupt_allowed = False

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
            logger.error(str(e))
            raise e from None
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
            logger.error(str(exc_value))
            raise exc_value from None
        return self

    def _sigint_handler(self, signal_received, frame):
        """Handles and reports if program is interrupted through the terminal"""
        logger.critical("Ctrl + C pressed, keyboard interupt initiated")
        self.__exit__(None, None, None)
        sys.exit(0)

    def finalize_database_write(self):
        """
        Methods to finalize when a database has been written to, and saving the current database schema to the sqlite database.
        """
        # index certain tables
        self._create_indices()
        logger.info("Database write session completed successfully.")

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
        insert_dict: dict,
        write_options: dict = {},
    ):
        """Inserts data from all arrays returned from results manager.

        Args:
            insert_dict
        """
        # parse ligand info form dict into list of ligands
        ligands_array = [
            docked_ligand["ligand_row"] for docked_ligand in insert_dict.values()
        ]
        # get unique ligand_id (will not add duplicate, instead return existing ligand_id)
        ligand_ids = self._insert_ligands(ligands_array)

        # add ligand ids to results array and make result array list of poses
        results_array = []
        for index, docked_ligand in enumerate(insert_dict.values()):
            for pose in docked_ligand["poses_results"]:
                results_array.append([ligand_ids[index]] + pose)
        # get unique pose ids and duplicate handling info
        Pose_IDs, duplicates = self._insert_results(results_array, write_options)
        # check if are interactions:
        if any(
            docked_ligand.get("poses_interactions") not in (None, [])
            for docked_ligand in insert_dict.values()
        ):

            self.insert_interactions(
                Pose_IDs,
                insert_dict,
                duplicates,
                write_options.get("duplicate_handling"),
            )

    def insert_receptor(self, receptor_data):
        receptors = self.fetch_receptor_objects()
        # insert receptor if database does not have already have a receptor entry
        if not receptors:
            self._insert_receptors(receptor_data)

    def insert_interactions(
        self,
        Pose_IDs: list,
        all_data: dict,
        duplicates: list,
        duplicate_handling: str = None,
    ):
        """Takes list of interactions, inserts into database

        Args:
            interactions_list (list): List of tuples for interactions in form
                ("type", "chain", "residue", "resid", "recname", "recid")
            duplicates (list(Pose_ID)): any duplicates identified in "insert_results", if duplicate handling has been specified

        """
        if duplicate_handling:
            if duplicate_handling.lower() not in ["ignore", "replace"]:
                logger.warning(
                    f"duplicate_handling option {duplicate_handling} not allowed. Reverting to default behavior."
                )
                duplicate_handling = None
        interaction_rows = self._insert_and_format_interactions(Pose_IDs, all_data)
        # NOTE here as of dev meeting
        self._insert_interaction_rows(interaction_rows, duplicates, duplicate_handling)

    # endregion

    # region read data
    def filter_results(
        self,
        all_filters: Filters,
        bookmark_name: str,
        filtering_bookmark: str = None,
        clustering: dict = {},
    ) -> iter:
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
            raise StorageError(
                f"Input database was created with Ringtail v{'.'.join([i for i in db_rt_version[:2]] + [db_rt_version[2:]])}. Confirm that this matches current Ringtail version and use Ringtail update script(s) to update database if needed."
            )

        # get the final filter query, has a {selection} place holder
        filter_query: str = self._generate_result_filtering_query(
            all_filters, bookmark_name, filtering_bookmark
        )

        filtered_poses = filter_query.format(selection=" R.Pose_ID ")
        logger.debug(f"Query for filtering results: {filter_query}")

        # perform filtering
        logger.debug("Running filtering query...")
        time0 = time.perf_counter()

        self.populate_filter_tables(
            name=bookmark_name,
            query=filtered_poses,
            filters=all_filters,
        )
        logger.debug(
            f"Time to filter results: {time.perf_counter() - time0:.2f} seconds"
        )
        count = self.get_passing_poses_count(bookmark_name, True)

        return count

    def bookmark_exists(self, bookmark_name: str):
        """Checks if bookmark name is in database

        Args:
            bookmark_name (str): name of bookmark name to check if exist

        Returns:
            bool: indicates if bookmark_name exists in the current database
        """

        return bool(bookmark_name in self.get_all_bookmark_names())

    def get_plot_data(self, bookmark_name: str, only_passing=False):
        """This function is expected to return an ascii plot
        representation of the results

        Args:
            bookmark_name (str): name of bookmark for which to fetch passing data. Returns empty list if bookmark does not exist.
            only_passing (bool): Only return data for passing ligands. Will return empty list for all data.

        Returns:
            tuple: cursors as (<all data cursor>, <passing data cursor>)
        """
        # TODO messy method
        # checks if we have filtered by looking for view name in list of view names
        if self.bookmark_exists(bookmark_name):
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
        temp_table_suffix: int,
        selection="NOT IN",
        old_db=None,
    ) -> tuple:
        """Selects ligands found or not found in the given bookmark in both current db and new_db. Stores as temp view

        Args:
            new_db (str): file name for database to attach
            bookmark1_name (str): string for name of first bookmark/temp table to compare
            bookmark2_name (str): string for name of second bookmark to compare
            selection (str): "IN" or "NOT IN" indicating if ligand names should or should not be in both databases
            old_db (str, optional): file name for previous database

        Returns:
            tuple: (name of new bookmark (str), number of ligands passing new bookmark (int))
        """
        if old_db is not None:
            self._detach_db(old_db.split(".")[0])  # remove file extension
        new_db_name = new_db.split(".")[0]  # remove file extension
        self._attach_db(new_db, new_db_name)

        selection = selection.upper().strip()
        if selection not in ["NOT IN", "IN"]:
            raise StorageError(f"Unrecognized selection type {selection}")

        temp_name = "temp_" + str(temp_table_suffix)
        self._create_crossref_temp_table(temp_name)
        temp_insert_query = self._generate_selective_insert_query(
            bookmark1_name, bookmark2_name, selection, new_db_name, temp_name
        )
        self.db_update(temp_insert_query, ())

        counting = QueryBuilder()
        count_pool = QueryBuilder()
        count_pool_string = (
            count_pool.SELECT("tt.pose_id")
            .FROM(temp_name, "tt")
            .JOIN("Results", "r", "pose_id")
            .GROUP_BY("r.ligand_id")
            .build()[0]
        )
        counting.WITH_SUBQUERY("count_pool", count_pool_string).SELECT("COUNT(*)").FROM(
            "count_pool", "cp"
        )
        num_passing = tuple(self.db_query(counting.build()[0]).fetchone())[0]
        print("\n\n Number passing the cross referenced filters: ", num_passing)

        return temp_name, num_passing

    def prune_nonpassing(self, bookmark_name: str):
        """Deletes rows from results, ligands, and interactions in a bookmark
        if they do not pass filtering criteria
        """
        self._delete_from_results(bookmark_name)
        self._delete_from_ligands(bookmark_name)
        self._delete_from_interactions_not_in_view(bookmark_name)
        # self._clear_bookmarks()
        # TODO probably needs to clea Filters and filtered_results

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
        }
        return groups[group]

    # endregion


class StorageManagerSQLite(StorageManager):
    """SQLite-specific StorageManager subclass

    Attributes:
        conn (SQLite.conn): Connection to database
        open_cursors (list): list of cursors that were not closed by the function that created them.
            Will be closed by close_connection method.
        db_file (str): database name
    """

    def __init__(
        self,
        db_file: str = None,
    ):
        self.db_file = db_file
        super().__init__()

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
        self._create_db_properties_table()
        self._create_filtering_tables()

        self._set_ringtail_db_schema_version(self._db_schema_ver)

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
        # TODO #TODO this is where I should organize it all in one dict, because this method only acts
        # on one file data packet/ligand_dict at the time
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
                # Check how things are parsed here, might not be most efficient
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
        # TODO this should have same length as above, so maybe I can zip them and evaluate together. This would read much clearer
        for idx, pose_interactions in enumerate(interaction_dictionaries):
            if not any(pose_interactions):  # skip any empty dictionaries
                continue
            interaction_tuples.append(
                cls._generate_interaction_tuples(pose_interactions)
            )
        data_dict = {
            "poses_results": result_rows,
            "poses_interactions": interaction_tuples,
            "ligand_row": cls._generate_ligand_row(ligand_dict),
            "receptor_row": cls._generate_receptor_row(ligand_dict),
        }
        return data_dict
        # return (
        #     result_rows,
        #     cls._generate_ligand_row(ligand_dict),
        #     interaction_tuples,
        #     cls._generate_receptor_row(ligand_dict),
        # )

    def _create_ligands_table(self):
        """Create table for ligands. Columns are:
        ligand_id           INTEGER PRIMARY KEY AUTOINCREMENT,
        LigName             VARCHAR NOT NULL UNIQUE ON CONFLICT IGNORE,
        ligand_smile        VARCHAR[],
        rdmol               BLOB,
        atom_index_map      VARCHAR[],
        hydrogen_parents    VARCHAR[],
        input_model         VARCHAR[]

        Raises:
            DatabaseTableCreationError: Description

        """
        ligand_table = """CREATE TABLE IF NOT EXISTS Ligands (
            ligand_id           INTEGER PRIMARY KEY AUTOINCREMENT,
            LigName             VARCHAR NOT NULL UNIQUE ON CONFLICT IGNORE,
            ligand_smile        VARCHAR[],
            rdmol        BLOB,
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
                [ligand_name, ligand_smile, ligand_rdbin, ligand_index_map,
                ligand_h_parents, input_model]
        """
        ligand_name = ligand_dict["ligname"]
        ligand_smile = ligand_dict["ligand_smile_string"]
        ligand_mol = Chem.MolFromSmiles(ligand_smile)
        # sanitize the ligand
        Chem.SanitizeMol(ligand_mol)
        ligand_rdbin = ligand_mol.ToBinary()
        ligand_index_map = json.dumps(ligand_dict["ligand_index_map"])
        ligand_h_parents = json.dumps(ligand_dict["ligand_h_parents"])
        input_model = json.dumps(ligand_dict["ligand_input_model"])

        return [
            ligand_name,
            ligand_smile,
            ligand_rdbin,
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
        rdmol,
        atom_index_map,
        hydrogen_parents,
        input_model
        ) VALUES
        (?,?,?,?,?,?)
        ON CONFLICT(LigName) DO 
            UPDATE SET LigName=excluded.LigName
        RETURNING ligand_id"""
        ligand_ids = []
        try:
            cur = self.conn.cursor()
            for ligand in ligand_array:
                ligand_id = cur.execute(sql_insert, ligand).fetchone()[0]
                ligand_ids.append(ligand_id)
            self.conn.commit()
            cur.close()
            return ligand_ids

        except sqlite3.OperationalError as e:
            raise DatabaseInsertionError("Error while inserting ligands.") from e

    def _delete_from_ligands(self, bookmark_name: str):
        """Remove rows from ligands table if they did not pass filtering

        Raises:
            StorageError
        """
        passing_poses_query = self.get_bookmark_poses_query(bookmark_name)

        try:
            self.db_update(
                f"DELETE FROM Ligands WHERE ligand_id NOT IN (SELECT ligand_id from Results WHERE Pose_ID IN ({passing_poses_query}))",
                (),
            )
        except sqlite3.OperationalError as e:
            raise StorageError(
                f"Error occured while pruning Ligands not in {bookmark_name}"
            ) from e

    def _create_results_table(self):
        """Creates table for results. Columns are:
        Pose_ID             INTEGER PRIMARY KEY AUTOINCREMENT,
        ligand_id           INTEGER FOREIGN KEY from Ligands,
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
            ligand_id           INT[],
            receptor            VARCHAR[],
            pose_rank           INT[],
            run_number          INT[],
            docking_score       FLOAT(4),
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
            flexible_res_coordinates   VARCHAR[],
            FOREIGN KEY (ligand_id) REFERENCES Ligands(ligand_id)
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
        # TODO needs to be ligand_id, needs that input
        ligand_data_list = [
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
        [0] ligand_id,
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
                        ligand_id=?
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
                logger.debug("Duplicate row not found.")
            else:
                Pose_ID = row[0]
                logger.debug(f"Duplicate row found for Pose_ID {Pose_ID}")
            cur.close()

            return Pose_ID

        except sqlite3.OperationalError as e:
            raise DatabaseQueryError(
                "Error while looking for unique result row."
            ) from e

    def _insert_results(self, results_array, options: dict):
        """Takes array of database rows to insert, adds data to results table. Will handle duplicates if specified

        Args:
            results_array (np.ndAaray): numpy array of arrays containing
                formatted result rows #TODO I am not sure this is actually numpy array, might just be a list

        Returns:
            Pose_ID (list(int)): returns the pose ids for the ligand written to results, these are used to ensure internal consistency when writing to the interaction table
            duplicates (list(int)): list of pose ids that are duplicates, if duplicate handling is specified. Filled with None if not specified or not duplicate

        Raises:
            DatabaseInsertionError
        """

        sql_insert = """INSERT INTO Results (
                        ligand_id,
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
                # TODO duplicate handling stuff
                if options.get("duplicate_handling"):
                    Pose_ID = self._check_unique_results_row(result)

                if Pose_ID != -1:  # row exists in table
                    duplicates.append(Pose_ID)
                    # row exist, evaluate if ignore or replace
                    if options.get("duplicate_handling").upper() == "IGNORE":
                        # do not add the new, duplicated row
                        pass
                    elif options.get("duplicate_handling").upper() == "REPLACE":
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

    def _delete_from_results(self, bookmark_name: str):
        """Remove rows from results table if they did not pass filtering

        Raises:
            StorageError
        """
        passing_poses_query = self.get_bookmark_poses_query(bookmark_name)
        try:
            self.db_update(
                f"DELETE FROM Results WHERE Pose_ID NOT IN ({passing_poses_query})", ()
            )
        except sqlite3.OperationalError as e:
            raise StorageError(
                f"Error occured while pruning Results not in {bookmark_name}"
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
                containing formatted receptor rows

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
            self.db_update(sql_insert, [receptor_array])

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

    def _insert_and_format_interactions(self, pose_ids, docking_data: dict):
        # insert interactions if they are present
        interaction_list = []
        pose_counter = 0
        for docked_ligand in docking_data.values():
            # for each pose of that ligand
            for pose_interactions in docked_ligand.get("poses_interactions"):
                # for each interaction of that pose
                Pose_ID = pose_ids[pose_counter]
                pose_interactions_with_poseid = [
                    (
                        (Pose_ID,)
                        + tuple(self._insert_interaction_index_row(interaction_tuple))
                    )
                    for interaction_tuple in pose_interactions
                ]

                interaction_list.extend(pose_interactions_with_poseid)
                pose_counter += 1

        return interaction_list

    def _insert_interaction_rows(
        self, interaction_rows, duplicates, duplicate_handling
    ):
        """Inserts the interaction data into a "tall-and-skinny" table, with a primary autoincremented key and a Pose_ID that is 1-to-1 with Results table.
        Table will contain as many rows with the same Pose_ID as that pose has interactions.

        Args:
            interaction_rows (list(tuple)): list of tuples containing the interaction data
            duplicates (list(int)): list of pose_ids from results table deemed duplicates, can also contain Nones, will be treated according to duplicate_handling
            duplicate_handling (str): how to handle duplicates
        Raises:
            DatabaseInsertionError
        """
        sql_insert = """INSERT INTO Interactions 
                            (Pose_ID,
                            interaction_id)
                            VALUES (?,?);"""
        try:
            cur = self.conn.cursor()
            if not duplicate_handling:  # add all results
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
                if duplicate_handling == "REPLACE" and duplicates_exist:
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

                elif duplicate_handling == "IGNORE":
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

    def _delete_from_interactions_not_in_view(self, bookmark_name: str):
        """Remove rows from interactions table if they did not pass filtering.

        Raises:
            StorageError: Description
        """
        passing_poses_query = self.get_bookmark_poses_query(bookmark_name)
        try:
            self.db_update(
                f"DELETE FROM Interactions WHERE Pose_ID NOT IN ({passing_poses_query})",
                (),
            )
            # remove unused interaction indices, if any
            self.db_update(
                """DELETE FROM Interaction_indices WHERE interaction_id IN
                            (SELECT ii.interaction_id FROM Interaction_indices ii 
                            LEFT JOIN Interactions i ON ii.interaction_id=i.interaction_id 
                            WHERE i.interaction_id IS NULL);""",
                (),
            )

        except sqlite3.OperationalError as e:
            raise StorageError(
                f"Error occured while pruning Interactions not in {bookmark_name}"
            ) from e

    def _create_filtering_tables(self):
        # Create filters table keeping track of filter id etc
        filters_sql = """CREATE TABLE IF NOT EXISTS Filters (
        filter_id           INTEGER PRIMARY KEY AUTOINCREMENT,
        name                VARCHAR[],
        query               VARCHAR[],
        filters             VARCHAR[]);"""

        filter_pose_sql = """CREATE TABLE IF NOT EXISTS Filtered_poses (
        filter_id           INTEGER,
        pose_id             INTEGER,
        FOREIGN KEY(filter_id) REFERENCES Filters(filter_id),
        FOREIGN KEY(pose_id) REFERENCES Results(pose_id));"""

        try:
            cur = self.conn.cursor()
            cur.execute(filters_sql)
            cur.execute(filter_pose_sql)
            cur.close()
        except sqlite3.OperationalError as e:
            raise DatabaseTableCreationError(
                "Error while creating bookmark table. If database already exists, use --overwrite to drop existing tables"
            ) from e

    def _insert_cluster_data(
        self,
        clusters: list,
        poseid_list: list,
        cluster_type: str,
        cluster_cutoff: str,
        bookmark_name: str,
    ):
        """Insert cluster data into ligand cluster table

        Args:
            clusters (list)
            poseid_list (list)
            cluster_type (str)
            cluster_cutoff (str)
            bookmark_name (str)
        """
        cur = self.conn.cursor()
        cur.execute(
            "CREATE TABLE IF NOT EXISTS Ligand_clusters (pose_id  INT[] UNIQUE)"
        )
        ligand_cluster_columns = self._fetch_ligand_cluster_columns()
        # TODO doesnt have to be this way
        column_name = (
            f"{bookmark_name}_{cluster_type}_{cluster_cutoff.replace('.', 'p')}"
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
            logger.debug("Creating columns index...")
            # TODO check performance without this index (would save space)
            cur.execute(
                "CREATE INDEX IF NOT EXISTS ak_results ON Results(docking_score, leff, deltas, reference_rmsd, energies_inter, energies_vdw, energies_electro, energies_intra, nr_interactions, run_number, pose_rank, num_hb)"
            )
            cur.execute("CREATE INDEX IF NOT EXISTS ak_poseid ON Results(Pose_id)")
            # TODO can I trim this one too, to the main one that filters? and shouldn't it include interaction_id
            cur.execute(
                "CREATE INDEX IF NOT EXISTS ak_intind ON Interaction_indices(interaction_type, rec_chain, rec_resname, rec_resid, rec_atom, rec_atomid)"
            )
            cur.execute(
                "CREATE INDEX IF NOT EXISTS ak_interactions ON Interactions(Pose_id, interaction_id)"
            )
            cur.execute("CREATE INDEX IF NOT EXISTS ak_ligands ON Ligands(ligand_id)")
            self.conn.commit()
            cur.close()
            logger.info(
                "Indicies were created for specified Results and Interaction_indices columns."
            )
        except sqlite3.OperationalError as e:
            raise StorageError("Error occurred while indexing") from e

    def _delete_table(self, table_name: str, db_alias: str = None):
        """
        Method to delete a table

        Args:
            table_name (str): table to be dropped

        """

        if db_alias:
            name = db_alias + "." + table_name
        else:
            name = table_name
        query = f"""DROP TABLE IF EXISTS {name};"""
        return self.db_query(query)

    def create_merge_tables(self) -> str:
        try:
            cur = self.conn.cursor()
            # create mergedata table: merge_id (PK), dbfile, timestamp, numofrows table
            mergetbl_sql = """CREATE TABLE IF NOT EXISTS merged_tables (
            merge_id                INTEGER PRIMARY KEY AUTOINCREMENT,
            dbfile                  VARCHAR[],
            merge_start             DATETIME DEFAULT CURRENT_TIMESTAMP);"""
            cur.execute(mergetbl_sql)

            # create PK table: merge_id(FK), table, original_val, merge_val
            pktable_sql = """CREATE TABLE IF NOT EXISTS PK_conversions (
            merge_id        INTEGER,
            table_name      VARCHAR[],
            original_PK     INTEGER,
            merged_PK       INTEGER,
            FOREIGN KEY(merge_id) REFERENCES merged_tables(merge_id));"""
            cur.execute(pktable_sql)
            cur.execute(
                "CREATE INDEX IF NOT EXISTS ak_merge ON PK_conversions(merge_id, original_PK)"
            )
            self.conn.commit()
        except Exception as e:
            raise StorageError(e) from e

    # endregion

    # region Methods for dealing with bookmarks/views and temporary tables
    def get_all_bookmark_names(self):
        """Get all bookmarks in sql database as a list of names. Bookmarks are a Ringtail-specific saved query (much like views)

        Returns:
            list: of bookmark names
        """
        try:
            cur = self.conn.cursor()
            # TODO good place to use row factory?
            cur.execute("SELECT name FROM Filters;")
            bookmark_names = [name[0] for name in cur.fetchall()]
            cur.close()

        except sqlite3.OperationalError as e:
            raise StorageError(
                "Error occured while fetching existing bookmark names"
            ) from e

        return bookmark_names

    def get_passing_poses_count(
        self, bookmark_name, grouped_by_ligand: bool = False
    ) -> int:
        if grouped_by_ligand:
            group_statement = "GROUP BY ligand_id"
        else:
            group_statement = ""

        query = QueryBuilder()
        query.SELECT("COUNT(r.pose_id)").FROM("results", "r").JOIN_BOOKMARK(
            bookmark_name, "bm"
        )
        if grouped_by_ligand:
            query.GROUP_BY("r.ligand_id")
        query_string = query.build(count=True)[0]

        if self.db_query(query_string).fetchone():
            return self.db_query(query_string).fetchone()[0]
        else:
            return 0

    @classmethod
    # TODO make this a method in query builder maybe?
    def get_bookmark_poses_query(self, bookmark_name: str):
        query = QueryBuilder()
        return query.SELECT("pose_id").FROM_BOOKMARK(bookmark_name).build()[0]

    def delete_bookmark(self, bookmark_name: str):
        # get filter id
        filter_id = self.db_query(
            "SELECT filter_id from Filters WHERE name = ?", (bookmark_name,)
        ).fetchone()["filter_id"]
        # delete from filters
        self.db_update(
            "DELETE FROM Filters WHERE filter_id = ?",
            (filter_id,),
        )
        # delete from filtered_poses table
        self.db_update(
            "DELETE FROM Filtered_poses WHERE filter_id = ?",
            (filter_id,),
        )
        logger.info(
            f"The bookmark {bookmark_name} and its associated filter data has been deleted."
        )

    def format_editable_filter_query(
        self, bookmark_name, results_str="R", filter_str="f"
    ) -> str:
        """
        Pre-formats a sqlite query string for retrieving {selection} columns from filtered poses.
        The string the {selection} as well as a {group_statement} when .format() after retrieval

        Args:
            bookmark_name (_type_): _description_
            results_str (str, optional): _description_. Defaults to "R".
            filter_str (str, optional): _description_. Defaults to "f".

        Returns:
            str: _description_
        """
        query = """SELECT {{selection}} FROM 
                    Results {result_alias} JOIN ( 
                            SELECT Pose_id FROM filtered_poses 
                            WHERE filter_id = 
                                (SELECT filter_id FROM Filters 
                                WHERE name = '{bookmark_name}')) {filter_alias}            
                    ON {result_alias}.Pose_ID = {filter_alias}.Pose_ID {{group_statement}}"""
        formatted_query = query.format(
            bookmark_name=bookmark_name,
            result_alias=results_str,
            filter_alias=filter_str,
        )
        return formatted_query

    def get_bookmark_selection(self, bookmark_name: str, selection: str):
        if not selection.startswith("R."):
            selection = "R." + selection

        return self.format_editable_filter_query(bookmark_name).format(
            selection=selection, group_statement=""
        )

    def fetch_filters_from_bookmark(self, bookmark_name: str) -> dict:
        """Method that will retrieve filter values used to construct bookmark

        Args:
            bookmark_name (str, optional): can get filter values for given bookmark, or filter values from currently active bookmark in storageman

            Returns:
                dict: containing the filter data
        """
        sql_query = "SELECT filters FROM Filters where name = ?"

        filters = self.db_query(sql_query, (bookmark_name,)).fetchone()
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

        count = self.db_query(f"SELECT COUNT(*) FROM {bookmark_name};").fetchone()[0]
        return bool(count > 0)

    def populate_filter_tables(self, name, query: str, filters={}) -> bool:

        # fetch filtered poses
        passing_poses_tuples = self.db_query(query).fetchall()
        passing_poses = [row[0] for row in passing_poses_tuples]
        if passing_poses:
            # check if bookmark exists
            if self.bookmark_exists(name):
                logger.warning(
                    f"The bookmark {name} already exists, and will be overwritten by the current filter."
                )
                self.delete_bookmark(name)
            filter_sql = """INSERT INTO Filters (name,query,filters) VALUES (?,?,?) RETURNING filter_id;"""
            try:
                filter_id = self.db_query(
                    filter_sql,
                    (
                        name,
                        query,
                        json.dumps(filters),
                    ),
                ).fetchone()[0]

                filter_pose_sql = f"""
                    INSERT INTO Filtered_poses (filter_id, pose_id) VALUES (?,?);"""
                params = [(filter_id, pose_id) for pose_id in passing_poses]
                self.db_update(filter_pose_sql, params)
            except Exception as e:
                raise StorageError(
                    f"Problems while writing filtered poses to database: {e}"
                ) from e
            else:
                logger.info("Successfully wrote filtered poses to database.")
        return bool(passing_poses)

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
            logger.info(f"Dropped bookmark {bookmark_name}.")
        except sqlite3.OperationalError as e:
            raise DatabaseInsertionError(
                f"Error while attempting to drop bookmark {bookmark_name}"
            ) from e

    def create_bookmark_from_temp_table(
        self,
        temp_table_name,
        bookmark_name,
        original_bookmark_name,
        wanted_list=[],
        unwanted_list=[],
    ):
        """Resaves temp bookmark stored in bookmark_name as new permenant bookmark

        Args:
            bookmark_name (str): name of bookmark to save last temp bookmark as
            original_bookmark_name (str): name of original bookmark
            wanted_list (list): List of wanted database names
            unwanted_list (list, optional): List of unwanted database names
            temp_table_name (str): name of temporary table
        """
        query = QueryBuilder()
        subq = QueryBuilder()
        subq_string = subq.SELECT("t.pose_id").FROM(temp_table_name, "t").build()[0]
        query_string = (
            query.SELECT("bm.pose_id")
            .FROM_BOOKMARK(original_bookmark_name, "bm")
            .WHERE(f"bm.pose_id IN ({subq_string})")
            .build()[0]
        )
        filters = {
            "comparison_wanted": ", ".join(wanted_list),
            "comparison_unwanted": ", ".join(unwanted_list),
        }

        self.populate_filter_tables(bookmark_name, query_string, filters)

    def _create_crossref_temp_table(self, table_name):
        """create temporary table with given name and with ligand name and pose_id information

        Args:
            table_name (str): name for temp table

        Raises:
            DatabaseTableCreationError
        """

        create_table_str = f"CREATE TEMP TABLE {table_name}(Pose_ID, ligname)"
        try:
            cur = self.conn.cursor()
            cur.execute(create_table_str)
            self.conn.commit()
            cur.close()
        except sqlite3.OperationalError as e:
            raise DatabaseTableCreationError(
                f"Error while creating temporary table {table_name}"
            ) from e

    def get_filterid_from_name(self, bookmark_name: str) -> int:
        return self.db_query(
            f"""SELECT filter_id FROM Filters WHERE name = ?;""",
            (bookmark_name,),
        ).fetchone()[0]

    # endregion

    # region Methods for getting information from database
    def fetch_receptor_object_by_name(self, rec_name):
        """Returns Receptor object from database for given rec_name

        Args:
            rec_name (str): Name of receptor to return object for

        Returns:
        str: receptor object as a string
        """

        cursor = self.db_query(
            """SELECT receptor_object FROM Receptors WHERE RecName LIKE ?""",
            (rec_name,),
        )
        return str(cursor.fetchone()[0])

    def fetch_receptor_objects(self):
        """Returns all Receptor objects from database

        Args:
            rec_name (str): Name of receptor to return object for

        Returns:
            iter (tuple): of receptor names and objects
        """

        cursor = self.db_query("SELECT RecName, receptor_object FROM Receptors")
        return cursor.fetchall()

    def fetch_data_for_passing_results(
        self, bookmark_name: str, outfields: Union[str, list], order_results: str = None
    ) -> iter:
        """Will return SQLite cursor with requested data for outfields for poses that passed filter in bookmark_name

        Returns:
            iter: sqlite cursor of data from passing data

        Raises:
            OptionError
        """
        outfields_list = self.format_output_fields(
            outfields, ligands_alias="L", results_alias="R"
        )

        bookmark_selection = self.get_bookmark_poses_query(bookmark_name)

        query = QueryBuilder()
        query.SELECT(*outfields_list).FROM("Results", "R").WHERE(
            f"R.pose_id IN ({bookmark_selection})"
        ).JOIN("ligands", "L", "ligand_id", "results").GROUP_BY("R.ligand_id")
        if order_results:
            order_by = self.format_orderby(order_results)
            query.ORDER_BY(order_by)

        return self.db_query(query.build()[0])

    def add_output_fields_to_query(self, outfields, bookmark_name):
        bookmark_selection = self.format_editable_filter_query(bookmark_name).format(
            selection="pose_id"
        )
        if outfields:
            outfield_dict = self.format_output_fields(outfields)
            outfield_string = outfield_dict.get("selection")
            join = outfield_dict.get("join")
        else:
            outfield_string = "Results.*"

        query = f"SELECT {outfield_string} FROM (SELECT * FROM Results WHERE Pose_ID IN ({bookmark_selection}))"
        if join:
            query = query + join

        return query

    def fetch_flexres_info(self, receptor):
        """fetch flexible residues names and atomname lists

        Returns:
            tuple: (flexible_residues, flexres_atomnames)
        """
        if type(receptor) == int:
            selection = "receptor_id = ?"
        elif type(receptor) == str:
            selection = "recname = ?"
        try:
            query = f"SELECT flexible_residues, flexres_atomnames FROM Receptors WHERE {selection}"
            info = self.db_query(query, (receptor,)).fetchone()
            if info is None:
                info = [], []
            return info
        except sqlite3.OperationalError as e:
            raise DatabaseQueryError("Error retrieving flexible residue info") from e

    def fetch_passing_ligands_rdkit_relevant_info(self, bookmark_name: str) -> iter:
        """fetch information required by vsmanager for writing out molecules

        Returns:
            iter: contains LigName, rdmol,
                atom_index_map, hydrogen_parents
        """
        bookmark_selection = self.get_bookmark_selection(bookmark_name, "ligand_id")

        query = f"""SELECT LigName, rdmol, atom_index_map, hydrogen_parents 
        FROM Ligands WHERE ligand_id IN (SELECT DISTINCT ligand_id FROM ({bookmark_selection}))"""
        return self.db_query(query)

    def fetch_ligand_rdkit_relevant_info(self, ligname) -> str:
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
                f"SELECT LigName, rdmol, atom_index_map, hydrogen_parents FROM Ligands WHERE LigName = '{ligname}'"
            )
            info = cur.fetchone()
            cur.close()
            return info
        except sqlite3.OperationalError as e:
            raise DatabaseQueryError(
                f"Error retrieving ligand info for {ligname}"
            ) from e

    def fetch_single_pose_properties(self, pose_ID: int) -> iter:
        # TODO where is this used
        """fetch coordinates for pose given by pose_ID

        Args:
            pose_ID (int): pose id to fetch coordinates for

        Returns:
            iter: SQLite cursor that contains Pose_ID, docking_score, leff, ligand_coordinates,
                flexible_res_coordinates, flexible_residues
        """
        query = f"SELECT Pose_ID, docking_score, leff, ligand_coordinates, flexible_res_coordinates FROM Results WHERE Pose_ID=?"
        return self.db_query(query, (pose_ID,))

    def fetch_interaction_info_by_index(self, interaction_idx) -> tuple:
        """Returns tuple containing interaction info for given interaction_idx

        Args:
            interaction_idx (int): interaction index to fetch info for

        Returns:
            tuple: tuple of info for requested interaction
        """
        query = "SELECT * FROM Interaction_indices WHERE interaction_id = ?"
        return self.db_query(query, (interaction_idx,)).fetchone()[
            1:
        ]  # cut off interaction index

    def fetch_pose_interactions(self, Pose_ID) -> iter:
        """
        Fetch all interactions parameters belonging to a Pose_ID

        Args:
            Pose_ID (int): pose id, 1-1 with Results table

        Returns:
            iter: of interaction information for given Pose_ID
        """
        # check if table exist
        cur = self.db_query(
            """SELECT name FROM sqlite_master WHERE type='table' AND name='Interactions';"""
        )
        if len(cur.fetchall()) == 0:
            return None

        query = f"""SELECT ii.interaction_type, ii.rec_chain, ii.rec_resname, ii.rec_resid, ii.rec_atom, ii.rec_atomid 
        FROM Interaction_indices ii 
        JOIN Interactions i ON i.interaction_id = ii.interaction_id
        WHERE i.Pose_ID = ?"""

        return self.db_query(query, (Pose_ID,)).fetchall()

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
        # TODO thi sshould be merged with next method, and just use columns as needed in the output
        # confusing that they give different kinds of data
        """Fetches cursor for best energies and leff for all ligands

        Returns:
             iter: SQLite Cursor containing docking_score,
                leff for the first pose for each ligand
        """
        return self.db_query(
            "SELECT docking_score, leff FROM Results GROUP BY ligand_id"
        )

    def _fetch_passing_plot_data(self, bookmark_name: str):
        """Fetches cursor for best energies and leffs for
            ligands passing filtering

        Args:
            bookmark_name (str): name for bookmark for which to fetch data. None will return data for default bookmark_name

        Returns:
            iter: SQL Cursor containing docking_score,
                leff for the first pose for passing ligands
        """
        query = QueryBuilder()
        query.SELECT("R.docking_score", "R.leff", "R.Pose_ID", "l.LigName").FROM(
            "Results", "R"
        ).JOIN("Ligands", "L", "ligand_id").WHERE(
            f"r.Pose_id IN ({self.get_bookmark_poses_query(bookmark_name)})"
        ).GROUP_BY(
            "r.ligand_id"
        )

        return self.db_query(query.build()[0])

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
                for c in self.db_query("PRAGMA table_info(Ligand_clusters)").fetchall()
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
        # TODO cna probably be optimized a bit
        # TODO this also will not work with bookmarks rn
        """Returns a panda dataframe of table or query given as requested_data

        Args:
            requested_data (str): String containing SQL-formatted query or table name
            table (bool): Flag indicating if requested_data is table name or not

        Returns:
            pd.DataFrame: dataframe of requested data
        """
        if table:
            return pd.read_sql_query(f"SELECT * FROM {requested_data}", self.conn)
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
        query = f"""SELECT COUNT(rowid) from {table_name}"""

        return self.db_query(query).fetchone()[0]

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
                "SELECT COUNT(ligand_id) FROM Ligands"
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
    def get_maxmiss_union(
        self, total_combinations: int, bookmark_name: str, all_filters={}
    ):
        """Get results that are in union considering max miss

        Args:
            total_combinations (int): numer of possible combinations

        Returns:
            iter: of passing results
        """
        enumerated_bookmark_queries = []
        existing_bookmarks = self.get_all_bookmark_names()
        for i in range(total_combinations):
            bmn = bookmark_name + "_" + str(i)
            if bmn in existing_bookmarks:
                result_alias = "R_" + str(i)
                filter_alias = "f_" + str(i)
                selection = f"{result_alias}.pose_id"
                partial_query = self.format_editable_filter_query(
                    bmn, result_alias, filter_alias
                ).format(selection=selection, group_statement="")
                enumerated_bookmark_queries.append(partial_query)

        bookmark_name = f"{bookmark_name}_union"
        logger.debug("Saving union bookmark...")
        union_view_query = " UNION ".join(enumerated_bookmark_queries)
        updated_query = union_view_query
        logger.debug("Running union query...")
        self.populate_filter_tables(
            name=bookmark_name, query=updated_query, filters=all_filters
        )

        count = self.get_passing_poses_count(bookmark_name, True)
        if not count:
            bookmark_name = None
        return count, bookmark_name

    def fetch_clustered_similars(self, ligname: str):
        """Given ligname, returns poseids for similar poses/ligands from previous clustering. User prompted at runtime to choose cluster.

        Args:
            ligname (str): ligname for ligand to find similarity with

        Raises:
            ValueError: wrong terminal input
            DatabaseQueryError
        """
        logger.warning(
            "N.B.: When finding similar ligands, export tasks (i.e. SDF export) will be for the selected similar ligands, NOT ligands passing given filters."
        )
        cur = self.conn.cursor()
        # TODO should be able to just check what clusters the ligand is a part of, and then offer them here
        # nothing somehow what filters were involved before the clustering perhaps?
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
        # TODO join with ligands on ligname
        query_ligand_cluster = cur.execute(
            f"SELECT {cluster_col_choice} FROM Ligand_clusters WHERE pose_id IN (SELECT Pose_ID FROM Results WHERE ligand_id = (SELECT ligand_id FROM Ligands WHERE LigName =  '{ligname}'))"
        ).fetchone()
        if query_ligand_cluster is None:
            raise DatabaseQueryError(
                f"Requested ligand name {ligname} not found in cluster {cluster_col_choice}!"
            )
        query_ligand_cluster = query_ligand_cluster[0]  # extract from tuple

        # TODO this is not ideal or optimized, but it does what it needs to do.
        # This is definitely an output place where group by ligname is needed
        sql_query = f"""
        SELECT LigName FROM Ligands WHERE ligand_id IN (SELECT ligand_id FROM Results WHERE Pose_ID IN 
        (SELECT pose_id FROM Ligand_clusters WHERE {cluster_col_choice}={query_ligand_cluster}))
        GROUP BY ligand_id"""

        # view_query = f"""
        # SELECT * FROM Results WHERE Pose_ID IN
        # (SELECT pose_id FROM Ligand_clusters WHERE {cluster_col_choice}={query_ligand_cluster})
        # GROUP BY LigName"""
        bookmark_name = f"similar_{ligname}_{cluster_col_choice}"
        # TODO needed?
        # probably write to passing results instead
        # self.create_bookmark(bookmark_name, view_query)

        return self.db_query(sql_query).fetchall(), bookmark_name, cluster_col_choice

    def fetch_rdkit_relevant_pose_properties(self, pose_ids: list) -> iter:
        placeholders = ",".join(["?"] * len(pose_ids))
        query = f"""
        SELECT Pose_ID, docking_score, leff, ligand_coordinates, flexible_res_coordinates 
        FROM Results WHERE Pose_ID IN ({placeholders})
        """
        return self.db_query(query, pose_ids)

    def fetch_passing_pose_properties(self, ligname, bookmark_name):
        """fetch coordinates for poses passing filter for given ligand

        Args:
            ligname (str): name of ligand to fetch coordinates for

        Returns:
            iter: SQLite cursor that contains Pose_ID, docking_score, leff, ligand_coordinates,
                flexible_res_coordinates, flexible_residues
        """
        bookmark_selection = self.get_bookmark_selection(
            bookmark_name, "R.pose_Id, R.ligand_id"
        )
        query = f"""
        SELECT Pose_ID, docking_score, leff, ligand_coordinates, flexible_res_coordinates 
        FROM Results WHERE Pose_ID IN 
        (({bookmark_selection}) WHERE ligand_id = (SELECT ligand_id from Ligands WHERE LigName = ?))"""
        return self.db_query(query, (ligname,))

    def fetch_nonpassing_pose_properties(self, ligname: str, bookmark_name):
        """fetch coordinates for poses of ligname which did not pass the filter

        Args:
            ligname (str): name of ligand to fetch coordinates for

        Returns:
            iter: SQLite cursor that contains Pose_ID, docking_score, leff, ligand_coordinates,
                flexible_res_coordinates, flexible_residues
        """
        bookmark_selection = self.get_bookmark_selection(bookmark_name, "Pose_id")
        query = f"""SELECT Pose_ID, docking_score, leff, ligand_coordinates, flexible_res_coordinates 
        FROM Results WHERE ligand_id = (SELECT ligand_id from Ligands where LigName = ?) AND Pose_ID NOT IN ({bookmark_selection})"""
        return self.db_query(query, (ligname,))

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
            logger.debug(f"Generating percentile filter query for {column}")
            cur = self.conn.cursor()
            cur.execute("SELECT COUNT(ligand_id) FROM Ligands")
            n_ligands = int(cur.fetchone()[0])
            n_passing = int((percentile / 100) * n_ligands)
            # find energy cutoff
            counter = 0
            for i in cur.execute(
                f"SELECT {column} FROM Results GROUP BY ligand_id ORDER BY {column}"
            ):
                if counter == n_passing:
                    cutoff = i[0]
                    break
                counter += 1
            logger.debug(f"{column} percentile cutoff is {cutoff}")
            return cutoff
        except sqlite3.OperationalError as e:
            raise StorageError("Error while generating percentile query") from e

    def calculate_percentile_from_value(self, docking_score_max=None, leff_max=None):
        # TODO can replace with kwargs, and check in the create table statement or something if the column
        # name is actually a numerical column or not
        if docking_score_max and leff_max:
            logger.warning(
                "Can not calculate percentil for both docking score and ligand efficiency, will proceed with just docking score"
            )
            leff_max = None
        if docking_score_max:
            column = "docking_score"
            value = docking_score_max
        elif leff_max:
            column = "leff"
            value = leff_max

        cur = self.conn.cursor()
        cur.execute("SELECT COUNT(ligand_id) FROM Ligands")
        n_ligands_total = int(cur.fetchone()[0])

        cur.execute(
            f"SELECT COUNT(*) FROM (SELECT Pose_ID FROM Results WHERE {column} < {value} GROUP BY ligand_id);"
        )
        n_ligands_passing = int(cur.fetchone()[0])

        return n_ligands_passing / n_ligands_total * 100

    def format_orderby(self, column_name: str):
        columns, aliased_columns = self.get_possible_output_columns()
        if column_name.lower() in columns:
            index = columns.index(column_name.lower())
            order_by = aliased_columns[index].format(
                Ligands_alias="L", Results_alias="R"
            )
            return order_by

    def format_output_fields(
        self, outfields: Union[str, list], results_alias="R", ligands_alias="L"
    ) -> str:
        """Handles string or list input of column names to be outputted, will make sure LigName
        is in the list, and make sure all options are valid

        Returns:
            list: column names for which the data is to be displayed that needs formatting with table alias
                for which table they belong to

        Raises:
            OptionError
        """
        if type(outfields) == str:
            outfields = outfields.replace(" ", "")
            outfields_list = outfields.split(",")
        elif type(outfields) == list:
            outfields_list = outfields
        else:
            logger.warning(
                "The provided outfields is not in a usable format (string or list). Will only use ligname"
            )
            outfields_list = []
        table_formatted_outfields = []
        if "ligname" not in [field.lower() for field in outfields_list]:
            outfields_list.insert(0, "LigName")
        possible_columns, table_formatted_columns = self.get_possible_output_columns()

        for outfield in outfields_list:
            if outfield.lower() in possible_columns:
                table_formatted_outfields.append(
                    table_formatted_columns[possible_columns.index(outfield.lower())]
                )
            else:
                logger.warning(
                    f"{outfield} is not a valid output option, and will be removed from the output columns. Please see rt_process_vs.py --help for allowed options"
                )
        formatted_outfields = [
            outfield.format(Ligands_alias=ligands_alias, Results_alias=results_alias)
            for outfield in table_formatted_outfields
        ]

        return formatted_outfields

    def get_possible_output_columns(self, tables=["Results", "Ligands"]):
        """
        _summary_

        Args:
            tables (list, optional): _description_. Defaults to ["Results", "Ligands"].

        Returns:
            columns (list[str]): _description_
            columns_with_tablename (list[str.format]): needs formatted with table_alias for use
        """
        columns = []
        columns_with_tablename = []
        for table in tables:
            columns_info = self.db_query(f"PRAGMA table_info({table})").fetchall()
            columns.extend([col[1].lower() for col in columns_info])
            columns_with_tablename.extend(
                [
                    (
                        "{{{table_alias}_alias}}.{col}".format(
                            col=col[1], table_alias=table
                        )
                    )
                    for col in columns_info
                ]
            )

        return columns, columns_with_tablename

    def get_numeric_columns(self, table_name: str) -> iter:
        """
        Method to get the names of all numeric columns, for example for allowable sorting options

        Args:
            table_name (str): table name to evaluate

        Returns:
            iter: column names that has a numeric type
        """
        return self.db_query(
            f"""SELECT
                            name
                        FROM
                            pragma_table_info('{table_name}')
                            WHERE CASE 
                                WHEN UPPER(type) LIKE '%INT%' THEN 'numerical'
                                WHEN UPPER(type) LIKE '%REAL%' THEN 'numerical'
                                WHEN UPPER(type) LIKE '%NUM%' THEN 'numerical'
                                WHEN UPPER(type) LIKE '%DEC%' THEN 'numerical'
                                WHEN UPPER(type) LIKE '%FLOAT%' THEN 'numerical'
                                WHEN UPPER(type) LIKE '%DOUBLE%' THEN 'numerical'
                            END ='numerical';"""
        ).fetchall()

    def _process_filters_for_query(self, filters_dict: dict):
        score_maxmin_to_sqlite_call = {
            "eworst": "docking_score < {value}",
            "ebest": "docking_score > {value}",
            "leworst": "leff < {value}",
            "lebest": "leff > {value}",
        }
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
                        score_maxmin_to_sqlite_call[filter_key].format(
                            value=filter_value
                        )
                    )

            # write hb count filter(s)
            if filter_key == "hb_count":
                for k, v in filter_value:
                    if k != "hb_count":
                        logger.warning(
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

    def _generate_result_filtering_query(
        self, filters_dict, bookmark_name: str, filter_bookmark: str
    ):
        """takes lists of filters, writes sql filtering string

        Args:
            filters_dict (dict): dict of filters. Keys names and value formats must match those found in the Filters class

        Returns:
            str: SQLite-formatted string for filtering query
        """
        # table to filter over
        filtering_window = "Results"
        num_query = ""
        int_query = ""
        ligname_query = ""
        partial_filter_query = ""
        rdkit_query = False

        # if filtering over a bookmark (i.e., already filtered results) as opposed to a whole database
        if filter_bookmark is not None:
            if filter_bookmark == bookmark_name:
                # cannot write data from bookmark_a to bookmark_a
                logger.error(
                    f"Specified 'filter_bookmark' and 'bookmark_name' are the same: {bookmark_name}"
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
            filtering_window = filter_bookmark

        # process filter values to lists and dicts that are easily incorporated in sql queries
        processed_filters = self._process_filters_for_query(filters_dict)

        if not processed_filters:
            raise OptionError("No filters were provided, cannot filter.")

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
            if "ligand_name" in lig_filters:
                lig_names = lig_filters["ligand_name"]
                ligname_query = " OR ".join(
                    [f"LigName LIKE '%{ligname}%' " for ligname in lig_names if ligname]
                )
                # TODO can this be made ligand_id?
                ligname_query = "SELECT ligand_id FROM Ligands WHERE " + ligname_query
            # rdkit queries need to be handled in memory separate from the main query
            if (
                "ligand_substruct" in lig_filters
                or "ligand_substruct_pos" in lig_filters
                or "ligand_max_atoms" in lig_filters
            ):
                rdkit_query = True

        ### Join each of the filter groups
        # if filter queries exist for each group, string them together appropriately
        if int_query:
            # add with a join statement
            partial_filter_query += f"JOIN ({int_query}) I ON R.Pose_ID = I.Pose_ID "
        if ligname_query:
            # add with a join statement
            partial_filter_query += (
                f"JOIN ({ligname_query}) L ON R.ligand_id = L.ligand_id "
            )
            # these two queries are joined on the Results table, after the multiple table spanning queries
        if num_query:
            # add quantitative condition
            partial_filter_query += "WHERE " + num_query

        # run the rdkit queries first
        if rdkit_query:
            lignames_poseids_with_substructs = {}
            # run the partial_filter_query with appropriate select statement, keeping track of passing pose ids
            partial_rdkit_query = f" FROM {filtering_window} R JOIN Ligands ON Ligands.ligand_id = R.ligand_id {partial_filter_query}"

            # get dict of ligands and pose ids passing all filters including substruct
            lignames_poseids_with_substructs = self._perform_rdkit_filtering(
                partial_rdkit_query, lig_filters
            )
            # make one list of all pose_ids passing so far
            passing_pose_ids = []
            for _, poseids in lignames_poseids_with_substructs.items():
                # loop through each list of pose ids and join them into a query
                passing_pose_ids.extend(poseids)

            # create new partial_filter_query with passing pose ids from the ligand queries
            partial_filter_query = " R.Pose_ID IN ({0})".format(
                ",".join(map(str, passing_pose_ids))
            )

            partial_filter_query = " WHERE " + partial_filter_query

        filter_query = (
            f"SELECT {{selection}} FROM {filtering_window} R {partial_filter_query}"
        )

        return filter_query

    def _perform_rdkit_filtering(
        self, partial_query: str, ligand_filters: dict
    ) -> dict:
        """Will run a filtering query with regular filters, then and perform in-memory filtering on
        whichever of the three RDKit filters substruct, substruct_pos, and maxheavtatoms.
        The method streams 100 ligands into memory at once from the cursor (while the cursor is
        still reading if database is big) and performs filtering on the one batch at the time.
        One row/ligand-pose_id from the stream is evaluated at once, from fastest to slowest filter:
        First, check if passing max_heavy_atoms filter if specified,
        Second, check if passing substructure filter if specified (does not account for substruct_pos substructures),
        Third, check if passing substruct_pos filter if specified
        Ligand-pose_id will only advance through the queries if it passes prior specified filter.
        If using logical operator "OR" to combine a set of substructs or substruct_pos', it will only test a ligand-pose_id
        until it passes one of the given filters, and will not test the remaining filter values for that filter.
        If a ligand passes a filter where pose is not relevant, it will add all pose_ids for that ligand
        in the resuling filtered_ligand dict.

        This method has several internal, private methods including
        _smarts_to_mol: creates mol object of substructure/SMARTS and checks for explicit hydrogens
        _stream_query: method that "streams" from the database and yields db cursor rows
        _substructure_position_calculation: checks if a substructure is in the right location
        _ligand_indexmap: calculates an index map for the ligand.

        Args:
            partial_query (str): partial query string for regular filters without a SELECT statement
            ligand_filters (dict): List of filters on ligand table

        Returns:
            dict: dict of ligand idsnames and all their pose ids passing regular+rdkit filters
        """
        maxatoms = 0
        position = False
        substruct_mols = []
        # handle single value filters
        if "ligand_operator" in ligand_filters:
            logical_operator = ligand_filters["ligand_operator"]
        if "ligand_max_atoms" in ligand_filters:
            maxatoms = ligand_filters["ligand_max_atoms"]

        def _smarts_to_mol(smarts: str) -> Chem.Mol:
            """
            creates a mol object from a smarts pattern if no explicit hydrogens

            Args:
                smarts (str): substructure smarts to match

            Raises:
                OptionError

            Returns:
                Mol: rdkit mol object
            """
            smarts_mol = Chem.MolFromSmarts(smarts)
            # make sure SMARTS is valid
            for atom in smarts_mol.GetAtoms():
                if atom.GetAtomicNum() == 1:
                    raise OptionError(
                        f"Given ligand substructure filter {smarts} contains explicit hydrogens. Please re-run query with SMARTs without hydrogen."
                    )
            return smarts_mol

        select_statement = "SELECT R.Pose_id, Ligands.ligand_id, Ligands.rdmol "
        # handle substruct
        if "ligand_substruct" in ligand_filters:
            for substruct in ligand_filters["ligand_substruct"]:
                substruct_mols.append(_smarts_to_mol(substruct))
        # whether or not doing a substruct position filter
        if "ligand_substruct_pos" in ligand_filters:
            substruct_pos = ligand_filters["ligand_substruct_pos"]
            position = True
            # we need additional info if doing position search
            select_statement += ", Ligands.atom_index_map, R.ligand_coordinates "
        # build full query
        query = select_statement + partial_query

        def _stream_query(query: str, batch_size: int = 100):
            """
            cursor stream generator

            Args:
                query (str): sql query
                batch_size (int): how many cursor hits to read into memory at once

            Yields:
                cursor batch: cursor results for the query in batch increments, where row can be read
            """
            cursor = self.conn.cursor()
            cursor.execute(query)

            while True:
                batch = cursor.fetchmany(batch_size)
                if not batch:
                    break
                for row in batch:
                    yield row

        def _substructure_position_calculation(
            idxmap, ligand_coordinates, smartsmol, filter
        ) -> bool:
            """
            Method that checks whether or not a substructure specified by smarts
            is present on a ligand in the specified location (by means of the filter values).

            Args:
                idxmap (dict): index map of the ligand
                ligand_coordinates (json): ligand coordinates
                smartsmol (mol): mol object of the smarts pattern
                filter (list[str]): filter values from user

            Returns:
                bool: whether or not the smarts in the pose qualified in the right position
            """
            # unpack filter values
            index = filter[0]
            sqdist = filter[1] ** 2
            x = filter[2]
            y = filter[3]
            z = filter[4]

            # calculate xyz space coordinates
            xyz = [
                float(value)
                for value in json.loads(ligand_coordinates)[idxmap[smartsmol[index]]]
            ]
            # calculate the sum of squares distances
            d2 = (xyz[0] - x) ** 2 + (xyz[1] - y) ** 2 + (xyz[2] - z) ** 2
            if d2 <= sqdist:
                return True
            else:
                return False

        def _ligand_indexmap(atom_index_map) -> dict:
            # TODO there is now an extra column to consider
            """
            Method that converts the atom index mapping in the database to one usable by the filter method

            Args:
                atom_index_map (json): indices of all atoms in the ligand

            Returns:
                dict: filter ready index map of the ligand
            """
            idxmap = [int(value) - 1 for value in json.loads(atom_index_map)]

            return {
                idxmap[j * 2]: idxmap[j * 2 + 1] for j in range(int(len(idxmap) / 2))
            }

        ligands_checked = 0
        filtered_ligands = {}
        for ligandrow in _stream_query(query):
            # substruct and maxatoms do not discriminate on poses, check if ligand has already been accounted for
            if not position and ligandrow["ligand_id"] in filtered_ligands:
                filtered_ligands[ligandrow["ligand_id"]].append(ligandrow["pose_id"])
            # the real workhorse
            else:
                # deserialize binary rdmol
                ligand_mol = Chem.Mol(ligandrow["rdmol"])
                # check if qualify for maxatoms
                if maxatoms > 0:
                    if not ligand_mol.GetNumHeavyAtoms() <= maxatoms:
                        # continue for ligandrow in _stream_query, ligand did not pass
                        continue

                # if there are substructures in the search
                if substruct_mols:
                    # count how many matches
                    SMARTS_match = 0
                    # check for each SMARTS if is substruct
                    for smarts_mol in substruct_mols:
                        # is_substruct = ligand_mol.GetSubstructMatch(smarts_mol)
                        if ligand_mol.GetSubstructMatch(smarts_mol):
                            SMARTS_match += 1
                            ligands_checked += 1
                            # if ligand substruct operator is OR, only one match is needed to qualify the ligand
                            if logical_operator == "OR":
                                # breaks the smarts_mol in substruct_mols loop
                                break
                    # check if ligand passed enough substruct queries
                    if (logical_operator == "OR" and SMARTS_match < 1) or (
                        logical_operator == "AND" and SMARTS_match < len(substruct_mols)
                    ):
                        # continue for ligandrow in _stream_query, ligand did not pass
                        continue

                # queries with substructure in a certain position
                if position:
                    # count how many matches
                    SMARTS_pos_match = 0
                    # check for each SMARTS if is substruct
                    for filterrow in substruct_pos:
                        # filterrow[0] should be the smarts pattern
                        smarts_mol = _smarts_to_mol(filterrow[0])
                        ligand_index_map = _ligand_indexmap(ligandrow["atom_index_map"])
                        ligand_coordinates = ligandrow["ligand_coordinates"]
                        # filterrow [1:] should be indices, distance allowance, and coordinates for smarts match
                        substruct_pos_filter = filterrow[1:]
                        for hit in ligand_mol.GetSubstructMatches(smarts_mol):
                            filter_match = _substructure_position_calculation(
                                ligand_index_map,
                                ligand_coordinates,
                                hit,
                                substruct_pos_filter,
                            )
                            if filter_match:
                                # count each passing filter
                                SMARTS_pos_match += 1
                                # break for hit in ligand_mol.GetSubstructMatches
                                break
                        if logical_operator == "OR" and SMARTS_pos_match == 1:
                            # break for filterrow in substruct_pos
                            break
                    # check if ligand passed enough substruct queries
                    if (logical_operator == "OR" and SMARTS_pos_match < 1) or (
                        logical_operator == "AND"
                        and SMARTS_pos_match < len(substruct_pos)
                    ):
                        # continue for ligandrow in _stream_query, ligand did not pass
                        continue

            # ligand only makes it here if it passed all specified rdkit filters
            # add pose id to list if ligand already in the list
            if ligandrow["ligand_id"] in filtered_ligands:

                filtered_ligands[ligandrow["ligand_id"]].append(ligandrow["pose_id"])
            # add new ligand in the list of passing ligands
            else:
                filtered_ligands[ligandrow["ligand_id"]] = [ligandrow["pose_id"]]

        return filtered_ligands

    def cluster_data(
        self,
        bookmark_name: str,
        cluster_type: str = "mfpt",
        cutoff: float = 0.5,
    ):
        logger.debug("Preparing to cluster")
        time0 = time.perf_counter()

        query = QueryBuilder()
        if cluster_type.lower() == "ifp":
            query.SELECT("r.pose_id", "r.leff").FROM_BOOKMARK(bookmark_name, "bm").JOIN(
                "Results", "R", "pose_id"
            )
            pose_ids, leffs = zip(*self.db_query(query.build()[0]).fetchall())
            leffs = list(leffs)
            pose_id_bitvectors = self._generate_interaction_bitvectors(pose_ids)
            pose_ids = list(pose_ids)
            bitvectors = list(pose_id_bitvectors.values())

            ibc = InteractionBitvectorCluster(pose_ids, leffs, bitvectors, cutoff)
            clusters, representatives = ibc.cluster()
            # the representatives needs to be written to the bookmark
            # create new bookmark
        elif cluster_type.lower() == "mfp":
            query.SELECT("r.pose_id", "r.leff", "l.rdmol").FROM("Results", "R").JOIN(
                "ligands", "l", "ligand_id", "results"
            ).WHERE(f"r.pose_id IN ({self.get_bookmark_poses_query(bookmark_name)})")

            pose_ids, leffs, rdmols = zip(*self.db_query(query.build()[0]).fetchall())
            mfpc = MorganFingerprintCluster(pose_ids, leffs, rdmols, cutoff)
            clusters, representatives = mfpc.cluster()

        # print some stats
        max_len = max(len(lst) for lst in clusters)
        min_len = min(len(lst) for lst in clusters)
        logger.info(f"Number of {cluster_type} clusters: {len(clusters)}")
        logger.info(
            f"Biggest {cluster_type} cluster contains {max_len} poses while the smallest cluster contains {min_len} poses."
        )
        # write to db
        self._insert_cluster_data(
            clusters,
            pose_ids,
            cluster_type.lower(),
            str(cutoff),
            bookmark_name,
        )

        bookmark_name = bookmark_name + "_" + cluster_type.lower() + "_clustered"
        clustered_poses = QueryBuilder()
        clustered_poses.SELECT("pose_id").FROM("results").WHERE(
            f"pose_id IN ({','.join(representatives)})"
        )

        self.populate_filter_tables(
            bookmark_name,
            clustered_poses.build()[0],
            {"cluster_type": cluster_type, "cutoff": cutoff},
        )

        logger.info(f"Time to cluster data: {time.perf_counter() - time0:.2f} seconds")

        return bookmark_name, len(clusters)

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

    def _generate_interaction_bitvectors(self, pose_ids: tuple[str]) -> dict:
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
        poseid_intind_query = """SELECT Pose_ID, interaction_id
                                    FROM Interactions
                                    WHERE Pose_ID IN ({placeholders});""".format(
            placeholders=",".join(["?"] * len(pose_ids))
        )
        poseid_intinds = self.db_query(poseid_intind_query, pose_ids).fetchall()
        # make dict of pose id and bitvector
        poseid_bvlist = {(pose_id): [0] * ii_length for pose_id in pose_ids}
        # iterate over the tuple results from the query
        for poseid_intind in poseid_intinds:
            poseid_bvlist[(poseid_intind[0])][poseid_intind[1] - 1] = 1

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
                    logger.warning(
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

        return self.db_query(sql_string).fetchall()

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
        # DUh, #TODO, they won't have the same pose ids for the same poses, so I hve to do this based on ligand name!
        # ligand id will also be different
        query = QueryBuilder()
        subq = QueryBuilder()
        subq_string = (
            subq.SELECT("l.ligname")
            .FROM_BOOKMARK(f"{bookmark2_name}", "bm2", new_db_name)
            .JOIN(
                f"{new_db_name}.results",
                "r",
                "pose_Id",
                f"bm2",
            )
            .JOIN(f"{new_db_name}.ligands", "l", "ligand_id", f"{new_db_name}.results")
            .build()[0]
        )
        query_string = (
            query.INSERT_INTO(temp_table)
            .SELECT("bm1.pose_id", "l.ligname")
            .FROM_BOOKMARK(bookmark1_name, "bm1")
            .JOIN("results", "r", "pose_id", "bm1")
            .JOIN("ligands", "l", "ligand_id", "results")
            .WHERE(f"l.ligname {select_str} ({subq_string})")
        ).build()[0]

        return query_string

    # endregion

    # region Database operations

    def overwrite_storage(self):
        """
        Will drop all tables in the database.
        """
        if not self._db_empty():
            self._drop_existing_tables()
            logger.info("Tables in existing database were dropped.")

    def _open_storage(self):
        """Create connection to db. Then, check if db needs to be created.
        If self.overwrite drop existing tables and initialize new tables

        Raises:
            StorageError
        """
        try:
            # check to see if file exist, and if it does, check that version is matching
            if os.path.isfile(self.db_file) and os.path.getsize(self.db_file) > 0:
                self.conn = self._create_connection()
                compatible, db_version = self.check_ringtaildb_version()
                # TODO this one is a problem
                if not compatible and not self.overwrite:
                    raise StorageError(
                        f"The database is of version {db_version} which is not compatible with the code base of version {version('ringtail')}"
                    )
            else:
                logger.info("Creating a new database file.")
                self.conn = self._create_connection()
            if self.keyboard_interrupt_allowed:
                signal(
                    SIGINT, self._sigint_handler
                )  # signal handler to catch keyboard interupts
            logger.debug(f"Ringtail connected to database {self.db_file}.")
        except Exception as e:
            raise StorageError(f"Error while creating or connecting to database: {e}.")

    def check_storage_ready(
        self, run_mode: str, docking_mode: str, store_all_poses: bool, max_poses: int
    ):
        """Check that storage is ready before proceeding, and creates new tables if needed

        Args:
            run_mode (str): if ringtail is ran using cmd line interface or api
            docking_mode (str): what docking engine was used to produce results
            store_all_poses (bool): overrwrites max poses
            max_poses (int): max poses to save to db

        Raises:
            StorageError
            OptionError: if database options are not compatible
        """
        if self._db_empty():
            self._create_tables()

        count = self.conn.execute("SELECT COUNT (*) FROM DB_properties").fetchone()[0]

        compatible = True
        if count < 1:
            logger.info(
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
            else:
                logger.warning(compatibility_string)

        # write current database properties to database
        if store_all_poses:
            number_of_poses = "all"
        else:
            number_of_poses = str(max_poses)
        self._insert_db_properties(docking_mode, number_of_poses)
        logger.debug("Storage compatibility has been checked and is ensured.")
        # cannot use Signal/keyboard interrupt in the GUI bc it uses multiple threads
        if run_mode != "gui":
            self.keyboard_interrupt_allowed = True

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
        logger.info(f"Database {self.db_file} was backed up to {backup_name}.")

    def _set_ringtail_db_schema_version(self, db_version: str = "2.0.0"):
        """Will check current storage manager db schema version and only set if it is compatible with the code base version (i.e., version(ringtail)).

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
            logger.info("Database version set to {0}".format(rtdb_version))
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
        if db_version == "0":
            # ringtail 1.0.0 did not have a user version, so catch if database has contents and version 0
            cur.execute(
                "SELECT EXISTS(SELECT 1 FROM sqlite_master WHERE type='table' AND name='Results');"
            )
            # if db version is 0 but has a results table, it is 1.0.0
            if cur.fetchone()[0] != 0:
                db_version = "100"
            # else empty or corrupt database
            else:
                raise StorageError(
                    f"The database requested {self.db_file} does not exist or does not have any tables. Check for spelling errors, else the database may be corrupt (delete the file before using the same name again)"
                )
        db_schema_ver = ".".join([*db_version])
        if version("ringtail") in self._db_schema_code_compatibility[db_schema_ver]:
            is_compatible = True
        else:
            is_compatible = False
            logger.warning(
                f"Database version {db_schema_ver} is NOT compatible with code base version {version('ringtail')}"
            )
        cur.close()
        return is_compatible, db_version

    def merge_databases(self, merging_db: str, backup: bool = True):
        # back up main database
        if backup:
            self.clone()
        # attach incoming database and check compatibility
        merging_db_alias = self._attach_db(merging_db, "merging")
        if not self._check_if_db_is_compatible(merging_db_alias, 200):
            raise StorageError(
                "Trying to merge two databases of incompatible or too old versions, cannot proceed."
            )
        # create new tables to keep track of merger
        self.create_merge_tables()
        # add to merging table the absolute path
        mergedb_abspath = str(os.path.abspath(merging_db))
        # insert merging db name first
        cur = self.conn.execute(
            """INSERT INTO merged_tables(dbfile) VALUES (?)""",
            (mergedb_abspath,),
        )
        # receptor compatibility check
        receptorcheck_sql = """
        SELECT CASE 
            WHEN Receptors.RecName = merging_receptors.RecName THEN 'True'
            ELSE 'False'
        END AS comparison_result 
        FROM Receptors 
        JOIN merging.Receptors AS merging_receptors 
            ON Receptors.receptor_id = merging_receptors.receptor_id"""

        try:
            assert self.db_query(receptorcheck_sql).fetchone()[0] == "True"
        except AssertionError:
            raise StorageError(
                f"The receptors in the merging databases are not the same. \nThese databases cannot be merged."
            )
        else:
            logger.info(
                "The two databases are of compatible version and receptors. Merging will proceed."
            )

        # get the active merge_id to ensure we use the correct merge data moving forward
        merge_id = self.db_query(
            "SELECT last_insert_rowid() FROM merged_tables"
        ).fetchone()[0]

        # delete incompatible tables
        # for main db
        self._drop_views()
        self._delete_table("Ligand_clusters")
        # for attached db
        self._drop_views(merging_db_alias)
        self._delete_table("Ligand_clusters", merging_db_alias)

        # merge tables
        try:
            self._merge_db_properties_table(merge_id)
            logger.info("The 'db_properties' table has been merged.")

            self._merge_ligands_table()
            logger.info("The 'Ligands' table has been merged.")

            self._merge_results_table(merge_id)
            logger.info("The 'Results' table has been merged.")

            self._merge_interaction_tables(merge_id)
            logger.info(
                "The 'Interaction_indices' and 'Interactions' tables have been merged."
            )
            # self._detach_db(merging_db_alias)
        except Exception as e:
            raise MergeError(f"Error during database merging: {e}") from e
        else:
            logger.info(
                f"The database {merging_db} has been successfully merged into {self.db_file}."
            )
        finally:
            cur.close()

    def _merge_interaction_tables(self, merge_id):
        # Interaction definitions are unique and results independent, so we only insert those that are new with updated PK,
        # and assign existing interaction_ids to those already described in primary db
        convert_ii_sql = """INSERT INTO PK_conversions (
        merge_id,
        table_name,
        original_PK,
        merged_PK) SELECT 
        ?,
        "Interaction_indices", 
        interaction_id,
            CASE 
                WHEN EXISTS (
                    SELECT 1
                    FROM Interaction_indices
                    WHERE 
                        merging.Interaction_indices.interaction_type = Interaction_indices.interaction_type
                        AND merging.Interaction_indices.rec_chain = Interaction_indices.rec_chain
                        AND merging.Interaction_indices.rec_resname = Interaction_indices.rec_resname
                        AND merging.Interaction_indices.rec_resid = Interaction_indices.rec_resid
                        AND merging.Interaction_indices.rec_atom = Interaction_indices.rec_atom
                        AND merging.Interaction_indices.rec_atomid = Interaction_indices.rec_atomid
                    ) 
                THEN (
                    SELECT main.Interaction_indices.interaction_id
                    FROM main.Interaction_indices
                    WHERE 
                        merging.Interaction_indices.interaction_type = Interaction_indices.interaction_type
                        AND merging.Interaction_indices.rec_chain = Interaction_indices.rec_chain
                        AND merging.Interaction_indices.rec_resname = Interaction_indices.rec_resname
                        AND merging.Interaction_indices.rec_resid = Interaction_indices.rec_resid
                        AND merging.Interaction_indices.rec_atom = Interaction_indices.rec_atom
                        AND merging.Interaction_indices.rec_atomid = Interaction_indices.rec_atomid
                    )
                ELSE merging.Interaction_indices.interaction_id + (SELECT MAX(interaction_id) FROM Interaction_indices)
            END AS new_interaction_id
        FROM merging.Interaction_indices;"""

        # then inserting only those that aren't already in the table
        insert_interaction_indices = """INSERT INTO Interaction_indices (
        interaction_id,
        interaction_type,
        rec_chain,
        rec_resname,
        rec_resid,
        rec_atom,
        rec_atomid)
        SELECT 
            (SELECT merged_PK FROM PK_conversions WHERE original_PK = interaction_id and merge_id = ? AND table_name = 'Interaction_indices') new_id,
            interaction_type,
            rec_chain,
            rec_resname,
            rec_resid,
            rec_atom,
            rec_atomid
        FROM merging.Interaction_indices WHERE new_id > (SELECT MAX(interaction_id) FROM Interaction_indices);
        """

        # Adding new data to Interactions table with (updated) pose_id and interaction_id
        insert_interactions = """
        INSERT INTO Interactions (
        Pose_ID,
        interaction_id
        )    SELECT P.merged_pk as pose_id, II.merged_pk as interaction_id
                FROM merging.Interactions I
                LEFT JOIN (SELECT original_PK, merged_pk
                FROM PK_conversions
                WHERE table_name = 'Results' 
                AND merge_id = ?) P ON (I.Pose_ID = P.original_PK)
            LEFT JOIN (SELECT original_PK, merged_pk
                FROM PK_conversions
                WHERE table_name = 'Interaction_indices' 
                AND merge_id = ?)  II ON (I.Interaction_ID = II.original_PK);"""

        try:
            cur = self.conn.cursor()
            cur.execute(
                convert_ii_sql,
                (merge_id,),
            )
            cur.execute(insert_interaction_indices, (merge_id,))
            cur.execute(
                insert_interactions,
                (
                    merge_id,
                    merge_id,
                ),
            )
            self.conn.commit()
        except Exception as e:
            raise Exception(
                f"Error during update and insertion of interactions: {e}"
            ) from e

    def _merge_results_table(self, merge_id):

        cur = self.conn.cursor()
        convert_poseid_sql = """INSERT INTO PK_conversions (
        merge_id,
        table_name,
        original_PK,
        merged_PK) SELECT 
        ?,
        "Results", 
        Pose_ID,
        Pose_ID + (SELECT MAX(Pose_ID) FROM Results) 
        FROM merging.Results;"""

        # insert results with updated pose_ids
        insert_Results = """INSERT INTO Results (
            Pose_ID,
            ligand_id,
            receptor,
            pose_rank,
            run_number,
            docking_score,
            leff,
            deltas,
            cluster_rmsd,
            cluster_size,
            reference_rmsd,
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
            flexible_res_coordinates) 
        SELECT 
            II.merged_PK as pose_id,
            I.ligand_id,
            I.receptor,
            I.pose_rank,
            I.run_number,
            I.docking_score,
            I.leff,
            I.deltas,
            I.cluster_rmsd,
            I.cluster_size,
            I.reference_rmsd,
            I.energies_inter,
            I.energies_vdw,
            I.energies_electro,
            I.energies_flexLig,
            I.energies_flexLR,
            I.energies_intra,
            I.energies_torsional,
            I.unbound_energy,
            I.nr_interactions,
            I.num_hb,
            I.about_x,
            I.about_y,
            I.about_z,
            I.trans_x,
            I.trans_y,
            I.trans_z,
            I.axisangle_x,
            I.axisangle_y,
            I.axisangle_z,
            I.axisangle_w,
            I.dihedrals,
            I.ligand_coordinates,
            I.flexible_res_coordinates
        FROM merging.Results I
                LEFT JOIN (SELECT original_PK, merged_PK 
                FROM PK_conversions 
                WHERE table_name = 'Results' 
                AND merge_id = ?) II ON II.original_PK = I.Pose_ID;"""

        try:
            cur.execute(
                convert_poseid_sql,
                (merge_id,),
            )
            cur.execute(insert_Results, (merge_id,))
            # commit after each large table in case there are memory constraints
            self.conn.commit()
        except Exception as e:
            raise StorageError(f"Error during insertion of results: {e}") from e

    def _merge_ligands_table(self):
        try:
            cur = self.conn.cursor()
            # Ligand table has unique constraints and primary key is ligand_id so no checks needed
            merge_ligands = """INSERT INTO Ligands 
                SELECT * FROM merging.Ligands"""
            cur.execute(merge_ligands)
            self.conn.commit()
        except Exception as e:
            raise StorageError("Error encountered while merging Ligands table") from e

    def _merge_db_properties_table(self, merge_id):
        try:
            cur = self.conn.cursor()
            convert_dbprop_sql = """INSERT INTO PK_conversions (
                merge_id,
                table_name,
                original_PK,
                merged_PK) SELECT 
                ?,
                "db_properties", 
                DB_write_session,
                DB_write_session + (SELECT MAX(DB_write_session) FROM db_properties) 
                FROM merging.db_properties;"""
            cur.execute(
                convert_dbprop_sql,
                (merge_id,),
            )

            insert_dbprops_sql = """INSERT INTO DB_properties (
                DB_write_session,
                docking_mode,
                number_of_poses)
                SELECT 
                    (SELECT merged_PK FROM PK_conversions WHERE original_PK = DB_write_session and merge_id = ? and table_name = 'DB_properties'),
                    docking_mode,
                    number_of_poses
                FROM merging.DB_properties;"""
            cur.execute(insert_dbprops_sql, (merge_id,))
            self.conn.commit()
        except Exception as e:
            raise StorageError(
                "Error encountered while merging db_properties table"
            ) from e

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
            logger.warning(
                "WARNING: All existing bookmarks in database will be dropped during database update!"
            )
            consent = input("Type 'yes' if you wish to continue: ") == "yes"
        if not consent:
            logger.critical("Consent not given for database update. Cancelling...")
            sys.exit(1)

        # get views and drop them
        logger.info(f"Updating {self.db_file}...")
        self._drop_views()

        # if current version is 1.0.0
        if self.check_ringtaildb_version()[1] == "1.0.0":
            # reformat for v1.1.0
            cur.execute(
                "ALTER TABLE Results RENAME COLUMN energies_binding TO docking_score"
            )
            cur.execute("ALTER TABLE Bookmarks ADD COLUMN filters")
            cur.execute(
                "CREATE INDEX IF NOT EXISTS ak_results ON Results(ligand_id, docking_score, leff, deltas, reference_rmsd, energies_inter, energies_vdw, energies_electro, energies_intra, nr_interactions, run_number, pose_rank, num_hb)"
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
            # TODO shouldn't this use the create_indices method
            cur.execute("""DROP TABLE IF EXISTS Interaction_bitvectors;""")
            # create new indixes
            cur.execute("CREATE INDEX IF NOT EXISTS ak_poseid ON Results(Pose_id)")
            cur.execute(
                "CREATE INDEX IF NOT EXISTS ak_interactions ON Interactions(Pose_id, interaction_id)"
            )
            cur.execute("CREATE INDEX IF NOT EXISTS ak_ligands ON Ligands(ligand_id)")
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

    def _drop_views(self, db_alias: str = None):
        if db_alias:
            alias_string = db_alias + "."
        else:
            alias_string = ""
        query = f"SELECT name FROM {alias_string}sqlite_master WHERE type = 'view'"
        cur = self.conn.execute(query)
        views = cur.fetchall()
        for v in views:
            cur.execute(f"DROP VIEW IF EXISTS {alias_string}{v[0]}")
        # delete all rows in bookmarks table
        cur.execute(f"DELETE FROM {alias_string}Bookmarks")

    def _create_connection(self):
        """Creates database connection to self.db_file

        Returns:
            SQLite.conn: Connection object to self.db_file

        Raises:
            DatabaseConnectionError
        """
        try:
            con = sqlite3.connect(self.db_file)
            con.row_factory = sqlite3.Row
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
        logger.debug("Closing database")
        self.conn.close()

    def _close_open_cursors(self):
        # TODO
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
            "SELECT COUNT(*) name FROM sqlite_master WHERE type='table' AND name <> 'sqlite_sequence';"
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

    def _attach_db(self, new_db: str, new_db_alias: str = "attached_db") -> str:
        """Attaches new database file to current database

        Args:
            new_db (str): file name for database to attach
            new_db_name (str): name of new database

        Raises:
            StorageError
        """
        attach_str = f"ATTACH DATABASE '{new_db}' AS {new_db_alias}"

        try:
            cur = self.conn.cursor()
            cur.execute(attach_str)
            self.conn.commit()
            cur.close()
        except sqlite3.OperationalError as e:
            raise StorageError(f"Error occurred while attaching {new_db}") from e
        else:
            logger.info(f"Attached database {new_db} aliased as {new_db_alias}.")
            return new_db_alias

    def _detach_db(self, new_db_alias):
        """Detaches new database file from current database

        Args:
            new_db_name (str): db name for database to detach

        Raises:
            StorageError
        """
        detach_str = f"DETACH DATABASE {new_db_alias}"
        try:
            cur = self.conn.cursor()
            cur.execute(detach_str)
            self.conn.commit()
            cur.close()
        except sqlite3.OperationalError as e:
            raise StorageError(f"Error occurred while detaching {new_db_alias}") from e
        else:
            logger.info(f"Detached database aliased as {new_db_alias}.")

    def _drop_existing_tables(self):
        """drop any existing tables.

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

    def db_query(self, query, params: tuple = ()):
        """Executes provided SQLite query. Returns cursor for results.
            Since cursor remains open, added to list of open cursors

        Args:
            query (str): Formated SQLite query as string

        Returns:
            SQLite cursor: Contains results of query
        """
        try:
            cur = self.conn.cursor()
            cur.execute(query, params)
            self.open_cursors.append(cur)
        except sqlite3.OperationalError as e:
            raise DatabaseQueryError(
                f"Unable to execute query {query} with given parameters {params}: {e}"
            ) from e
        return cur

    def db_update(self, query: str, parameters: list[tuple], commit=True) -> iter:
        if type(parameters) == tuple:
            parameters = [parameters]
        elif type(parameters) != list:
            raise OptionError(
                "Cannot use a non-list or non-tuple as insert parameters, please format appropriately."
            )
        try:
            cur = self.conn.cursor()
            cur.executemany(query, parameters)
            if commit:
                self.conn.commit()
        except sqlite3.OperationalError as e:
            raise DatabaseInsertionError(f"Error while committing insert query") from e

    # endregion
