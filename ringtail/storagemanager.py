#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail storage adaptors
#

import sqlite3
import time
import json
import pandas as pd
import logging
import typing
import sys
from signal import signal, SIGINT
from rdkit import Chem
from rdkit import DataStructs
from rdkit.ML.Cluster import Butina
import numpy as np
import time
import pkg_resources

try:
    import cPickle as pickle
except ImportError:
    import pickle

from .filters import Filters
from .exceptions import (
    StorageError,
    DatabaseInsertionError,
    DatabaseConnectionError,
    DatabaseTableCreationError,
)
from .exceptions import DatabaseQueryError, DatabaseViewCreationError, OptionError


class StorageManager:
    """Prototype class for a generic VS database object
    this class defines the API of each
    StorageManager class (including sub-classess)
    which will implement their own functions to return the data requested

    Attributes:
        conflict_opt (str): string indicating how conficting entries should be handled
        current_view_name (str): name of current results view
        db_file (string): Name of file containing database
        field_to_column_name (dictionary): Dictionary for
            converting command-line field options into DB column names
        interaction_filter_types (set): Set for types of interaction filters
        interactions_initialized_flag (boolean): Flag indicating if
            interaction tables have been created
        output_all_poses (boolean): Flag dictating, if true, all poses from passing ligands will be logged instead of just top one
        next_unique_interaction_idx (int): Index for the next unique
            interaction to be added to interaction_index table
        open_cursors (list): Storage for any DB cursors which are opened
            and not closed by the function that opened them. Will be closed by
            close_connection method.
        opts (dict): Dictionary of database options
        order_results (string): string indicating result field that passing
            results should be ordered by in log.
        overwrite_flag (boolean): Flag indictating that Storage should drop
            all existing tables and add nes data
        passing_results_view_name (string): Name for the view of passing
            results to be created after filtering
        temptable_suffix (int): suffix number for temp table
        unique_interactions (dict): Dictionary for storing unique interactions
            to be written in interaction_index table
        view_suffix (str): suffix to add to end of view. Used with max_miss
    """

    # initialize dictionary processing kw lists
    stateVar_keys = ["pose_about", "pose_translations", "pose_quarternions"]
    ligand_data_keys = [
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
    ]
    interaction_data_kws = [
        "type",
        "chain",
        "residue",
        "resid",
        "recname",
        "recid",
    ]

    def __init__(self):
        """Initialize instance variables common to all StorageManager subclasses"""

        self.unique_interactions = {}
        self.next_unique_interaction_idx = 1
        self.interactions_initialized_flag = False
        self.closed_connection = False

    @classmethod
    def get_defaults(cls, storage_type):
        storage_types = {
            "sqlite": StorageManagerSQLite,
        }
        return storage_types[storage_type].get_defaults()

    @classmethod
    def get_default_types(cls, storage_type):
        storage_types = {
            "sqlite": StorageManagerSQLite,
        }
        return typing.get_type_hints(storage_types[storage_type].__init__)

    def __enter__(self):
        self.open_storage()
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        if not self.closed_connection:
            self.close_storage()

    def _sigint_handler(self, signal_received, frame):
        logging.critical("Ctrl + C pressed, keyboard interupt initiated")
        self.__exit__(None, None, None)
        sys.exit(0)

    # # # # # # # # # # # # # # # # # # #
    # # # Common StorageManager methods # # #
    # # # # # # # # # # # # # # # # # # #

    def get_plot_data(self, only_passing=False):
        """this function is expected to return an ascii plot
        representation of the results

        Returns:
            Tuple: cursors as [<all data cursor>, <passing data cursor>]
        """

        # checks if we have filtered by looking for view name in list of view names
        if self.check_passing_view_exists():
            if only_passing:
                return [], self._fetch_passing_plot_data()
            else:
                return self._fetch_all_plot_data(), self._fetch_passing_plot_data()
        else:
            return self._fetch_all_plot_data(), []

    def prune(self):
        """Deletes rows from results, ligands, and interactions
        if they do not pass filtering criteria
        """
        self._delete_from_results()
        self._delete_from_ligands()
        self._delete_from_interactions()

    def check_passing_view_exists(self):
        """Return if self.results_view_name in database

        Returns:
            Bool: indicates if self.results_view_name exists
        """
        return self.results_view_name in [
            name[0] for name in self._fetch_view_names().fetchall()
        ]

    def close_storage(self, attached_db=None, vacuum=False):
        """close connection to database

        Args:
            attached_db (str, optional): name of attached DB (not including file extension)
            vacuum (bool, optional): indicates that database should be vacuumed before closing
        """
        if attached_db is not None:
            self._detach_db(attached_db)
        # drop indices created when filtering
        self._remove_indices()
        # close any open cursors
        self._close_open_cursors()
        # vacuum database
        if vacuum:
            self._vacuum()
        # close db itself
        self._close_connection()
        self.closed_connection = True

    def _add_unique_interactions(self, interactions_list):
        """takes list of interaction tuple lists. Examines
        self.unique_interactions, add interactions if not already inserted.

        self.unique_interactions {(interaction tuple): unique_interaction_idx}

        Args:
            interactions_list (List): List of interaction tuples
            for insertion into database
        """

        for pose in interactions_list:
            for interaction_tuple in pose:
                if interaction_tuple not in self.unique_interactions:
                    self.unique_interactions[
                        interaction_tuple
                    ] = self.next_unique_interaction_idx
                    if self.interactions_initialized_flag or self.append_results:
                        self._insert_one_interaction(interaction_tuple)
                        self._make_new_interaction_column(
                            self.next_unique_interaction_idx
                        )
                    self.next_unique_interaction_idx += 1

    @classmethod
    def format_for_storage(cls, ligand_dict: dict) -> tuple:
        """takes file dictionary from the file parser, formats required storage format

        Args:
            ligand_dict (dict): Dictionary containing data from the fileparser
        """

        raise NotImplementedError

    def insert_data(
        self,
        results_array,
        ligands_array,
        interaction_array,
        receptor_array=[],
        insert_receptor=False,
    ):
        """Summary

        Args:
            data_dictionaries (list): list of dictionaries of data to be stored
            insert_receptor (bool, optional): flag indicating that receptor info should inserted
        """
        raise NotImplementedError

    def set_view_suffix(self, suffix):
        """Sets internal view_suffix variable

        Args:
            suffix(str): suffix to attached to view-related queries or creation
        """

        self.view_suffix = suffix

    def filter_results(
        self,
        all_filters: dict,
        suppress_output=False
    ) -> iter:
        """Generate and execute database queries from given filters.

        Args:
            all_filters (dict): dict containing all filters. Expects format and keys corresponding to ringtail.Filters().to_dict()


        Returns:
            SQLite Cursor: Cursor of passing results
        """
        # create view of passing results
        filter_results_str, view_query = self._generate_result_filtering_query(
            all_filters
        )
        logging.debug(filter_results_str)
        # if max_miss is not 0, we want to give each passing view a new name by changing the self.results_view_name
        if self.view_suffix is not None:
            self.current_view_name = self.results_view_name + "_" + self.view_suffix
        else:
            self.current_view_name = self.results_view_name

        self._create_view(
            self.current_view_name, view_query
        )  # make sure we keep Pose_ID in view
        self._insert_bookmark_info(self.current_view_name, view_query, all_filters)
        # perform filtering
        if suppress_output:
            return None

        logging.debug("Running filtering query...")
        time0 = time.perf_counter()
        filtered_results = self._run_query(filter_results_str)
        logging.debug(f"Time to run query: {time.perf_counter() - time0:.2f} seconds")

        # get number of passing ligands
        return filtered_results

    def fetch_data_for_passing_results(self) -> iter:
        """Will return SQLite cursor with requested data for outfields for poses that passed filter in self.results_view_name

        Returns:
            sqlite cursor: cursor of data from passing data
        """
        return self._run_query(self._generate_results_data_query(self.outfields))

    def _fetch_view_names(self):
        """Returns DB curor with the names of all view in DB

        Returns:
            sqlite cursor: cursor of view names
        """
        return self._run_query(self._generate_view_names_query())

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
            selection_type (string): "+" or "-" indicating if ligand names should ("+") or should not "-" be in both databases
            old_db (str, optional): file name for previous database
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

        num_passing = self.get_number_passing_ligands(temp_name)

        self.temptable_suffix += 1

        return temp_name, num_passing

    # # # # # # # # # # # # # # # # #
    # # # Child-specific methods # # #
    # # # # # # # # # # # # # # # # #

    def insert_results(self, results_array):
        """takes array of database rows to insert, adds data to results table

        Args:
            results_array (list): list of lists
                containing formatted result rows

        """
        raise NotImplementedError

    def insert_ligands(self, ligand_array):
        """Takes array of ligand rows, inserts into Ligands table.

        Args:
            ligand_array (list): List of lists
                containing formatted ligand rows

        """
        raise NotImplementedError

    def insert_receptors(self, receptor_array):
        """Takes array of receptor rows, inserts into Receptors table

        Args:
            receptor_array (list): List of lists
                containing formatted ligand rows
        """
        raise NotImplementedError

    def insert_interactions(self, interactions_list):
        """generic function for inserting interactions from given
            interaction list into DB

        Args:
            interactions_list (list): List of tuples for interactions
                in form
                ("type", "chain", "residue", "resid", "recname", "recid")
        """
        raise NotImplementedError

    def save_receptor(self, receptor, rec_name):
        """Takes object of Receptor class, updates the column in Receptor table for the row with rec_name

        Args:
            receptor (Receptor): Receptor object to be inserted into DB
            rec_name (string): Name of receptor. Used to insert into correct row of DB

        Raises:
            NotImplementedError: Description
        """
        raise NotImplementedError

    def fetch_receptor_object_by_name(self, rec_name):
        """Returns Receptor object from database for given rec_name

        Args:
            rec_name (string): Name of receptor to return object for

        Raises:
            NotImplementedError: Description
        """
        raise NotImplementedError

    def get_results(self):
        """Gets all fields for filtered results

        No Longer Returned:
            DB cursor: Cursor with all fields and rows in passing results view

        Raises:
            NotImplementedError: Description
        """
        raise NotImplementedError

    def clone(self):
        """Creates a copy of the db

        Raises:
            NotImplementedError: Description
        """
        raise NotImplementedError

    def remake_bookmarks(self):
        """Reads all views from Bookmarks table and remakes them"""
        raise NotImplementedError

    def get_number_passing_ligands(self):
        """Returns count of ligands that passed filtering criteria

        No Longer Returned:
            Int: Number of passing ligands

        Raises:
            NotImplementedError: Description
        """
        raise NotImplementedError

    def fetch_flexres_info(self):
        """fetch flexres names and atomname_lists

        No Longer Returned:
            DB cursor: contains
                flexible_residues, flexres_atomnames

        Raises:
            NotImplementedError: Description
        """
        raise NotImplementedError

    def fetch_passing_ligand_output_info(self):
        """fetch information required by vsmanager for writing out molecules

        No Longer Returned:
            DB cursor: contains
                LigName, ligand_smile, atom_index_map, hydrogen_parents

        Raises:
            NotImplementedError: Description
        """
        raise NotImplementedError

    def fetch_passing_pose_properties(self, ligname):
        """fetch coordinates for poses passing filter for given ligand

        Args:
            ligname (string): name of ligand to return coordinates for

        No Longer Returned:
            DB cursor: contains
                Pose_ID, docking_score, leff, ligand_coordinates, flexible_res_coordinates, flexible_residues

        Raises:
            NotImplementedError: Description
        """
        raise NotImplementedError

    def fetch_nonpassing_pose_properties(self, ligname):
        """fetch coordinates for poses of ligname which did not pass the filter

        Args:
            ligname (string): name of ligand to fetch coordinates for

        No Longer Returned:
            DB cursor: contains
                Pose_ID, docking_score, leff, ligand_coordinates, flexible_res_coordinates, flexible_residues

        Raises:
            NotImplementedError: Description
        """
        raise NotImplementedError

    def fetch_interaction_bitvector(self, pose_id):
        """Returns tuple containing interaction bitvector line for given pose_id

        Args:
            pose_id (int): pose id to fetch interaction bitvector for

        No Longer Returned:
            tuple: tuple representing interaction bitvector

        Raises:
            NotImplementedError: Description
        """
        raise NotImplementedError

    def fetch_interaction_info_by_index(self, interaction_idx):
        """Returns tuple containing interaction info for given interaction_idx

        Args:
            interaction_idx (int): interaction index to fetch info for

        No Longer Returned:
            tuple: tuple of info for requested interaction

        Raises:
            NotImplementedError: Description
        """
        raise NotImplementedError

    def to_dataframe(self, requested_data: str, table=True) -> pd.DataFrame:
        """Returns dataframe of table or query given as requested_data

        Args:
            requested_data (str): String containing SQL-formatted query or table name
            table (bool): Flag indicating if requested_data is table name or not

        Raises:
            NotImplementedError: Description
        """
        raise NotImplementedError

    def fetch_view(self, viewname):
        """returns SQLite cursor of all fields in viewname

        Args:
            viewname (TYPE): Description

        Raises:
            NotImplementedError: Description
        """
        raise NotImplementedError

    def save_temp_bookmark(self, bookmark_name, original_bookmark_name):
        """Resaves temp bookmark stored in self.current_view_name as new permenant bookmark

        Args:
            bookmark_name (string): name of bookmark to save last temp bookmark as
            original_bookmark_name (TYPE): Description

        Deleted Parameters:
            orginal_bookmark_name (string): Name of original bookmark to pull data from

        Raises:
            NotImplementedError: Description
        """
        raise NotImplementedError

    def _create_connection(self):
        """Creates database connection to self.db_file

        No Longer Returned:
            DB connection: Connection object to self.db_file

        Raises:
            NotImplementedError: Description

        """
        raise NotImplementedError

    def _close_connection(self, attached_db=None):
        """Closes connection to database

        Args:
            attached_db (None, optional): Description

        Raises:
            NotImplementedError: Description
        """
        raise NotImplementedError

    def _close_open_cursors(self):
        """closes any cursors stored in self.open_cursors.
        Resets self.open_cursors to empty list

        Raises:
            NotImplementedError: Description
        """
        raise NotImplementedError

    def open_storage(self):
        """Create connection to db. Then, check if db needs to be written.
        If so, (if self.overwrite drop existing tables and )
        initialize the tables

        Raises:
            NotImplementedError: Description
        """
        raise NotImplementedError

    def _fetch_existing_table_names(self):
        """Returns list of all tables in database

        Raises:
            NotImplementedError: Description
        """
        raise NotImplementedError

    def _drop_existing_tables(self):
        """drop any existing tables. Will only be called
        if self.overwrite is true

        Raises:
            NotImplementedError: Description
        """
        raise NotImplementedError

    def _drop_existing_views(self):
        """drop any existing views. Will only be called
        if self.overwrite is true

        Raises:
            NotImplementedError: Description
        """
        raise NotImplementedError

    def _run_query(self, query):
        """Executes provided SQLite query. Returns cursor for results

        Args:
            query (string): Formated SQLite query as string

        No Longer Returned:
            DB cursor: Contains results of query

        Raises:
            NotImplementedError: Description
        """
        raise NotImplementedError

    def create_indices(self):
        """Create indices for columns in self.index_columns

        Args:
            index_lignames (bool, optional): flag indicating that index should be created over ligand names

        Raises:
            StorageError: Description
        """
        raise NotImplementedError

    def _remove_indices(self):
        """Removes idx_filter_cols and idx_ligname"""
        raise NotImplementedError

    def _create_view(self, name, query, temp=False):
        """takes name and selection query,
            creates view of query stored as name.

        Args:
            name (string): Name for view which will be created
            query (string): DB-formated query which will be used to create view
            temp (bool, optional): Flag if view should be temporary

        Raises:
            NotImplementedError: Description
        """
        raise NotImplementedError

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
        flexible_residues   VARCHAR[],
        flexible_res_coordinates   VARCHAR[]

        Raises:
            NotImplementedError: Description
        """
        raise NotImplementedError

    def _create_ligands_table(self):
        """Create table for ligands. Columns are:
        LigName             VARCHAR NOT NULL,
        ligand_smile        VARCHAR[],
        ligand_rdmol        MOL,
        atom_index_map      VARCHAR[],
        hydrogen_parents    VARCHAR[],
        input_model         VARCHAR[]

        Raises:
            NotImplementedError: Description

        """
        raise NotImplementedError

    def _create_receptors_table(self):
        """Create table for receptors. Columns are:
        Receptor_ID         INTEGER PRIMARY KEY AUTOINCREMENT,
        RecName                VARCHAR NOT NULL,
        box_dim             VARCHAR[],
        box_center          VARCHAR[],
        grid_spacing        INT[],
        flexible_residues   VARCHAR[],
        flexres_atomnames   VARCHAR[],
        receptor_object     BLOB

        Raises:
            NotImplementedError: Description
        """
        raise NotImplementedError

    def _create_interaction_index_table(self):
        """create table of data for each unique interaction. Columns are:
        interaction_id      INTEGER PRIMARY KEY AUTOINCREMENT,
        interaction_type    VARCHAR[],
        rec_chain           VARCHAR[],
        rec_resname         VARCHAR[],
        rec_resid           VARCHAR[],
        rec_atom            VARCHAR[],
        rec_atomid          VARCHAR[]

        Raises:
            NotImplementedError: Description

        """
        raise NotImplementedError

    def _create_interaction_bv_table(self):
        """Create table of interaction bits for each pose. Columns are:
        Pose_ID INTERGER PRIMARY KEY AUTOINCREMENT
        Interaction_1
        Interaction_2
        ...
        Interaction_n

        Raises:
            NotImplementedError: Description

        """
        raise NotImplementedError

    def _create_bookmark_table(self):
        """Create table of bookmark names and their queries. Columns are:
        Bookmark_name
        Query

        Raises:
            NotImplementedError: Description
        """
        raise NotImplementedError

    def _insert_bookmark_info(self, name: str, query: str):
        """Insert bookmark info into bookmark table

        Args:
            name (str): name for bookmark
            query (str): query used to generate bookmark

        Raises:
            NotImplementedError: Description
        """
        raise NotImplementedError

    def _insert_unique_interactions(self, unique_interactions):
        """Inserts interaction data for unique interactions
            into Interaction_index table

        Args:
            unique_interactions (list): List of tuples of
                interactions to be inserted

        Raises:
            NotImplementedError: Description
        """
        raise NotImplementedError

    def _insert_one_interaction(self, interaction):
        """Insert interaction data for a single new interaction
            into the interaction indices table

        Args:
            interaction (tuple): Tuple of interaction data
            (interaction_type, rec_chain, rec_resname,
            rec_resid, rec_atom, rec_atomid)

        Raises:
            NotImplementedError: Description

        """
        raise NotImplementedError

    def _make_new_interaction_column(self, column_number):
        """Add column for new interaction to interaction bitvector table

        Args:
            column_number (int): Index for new interaction

        Raises:
            NotImplementedError: Description
        """
        raise NotImplementedError

    def _fetch_all_plot_data(self):
        """Fetches cursor for best energies and leff for all ligands

        No Longer Returned:
            DB Cursor: Cursor containing docking_score,
                leff for the first pose for each ligand

        Raises:
            NotImplementedError: Description
        """
        raise NotImplementedError

    def _generate_plot_all_results_query(self):
        """Make DB-formatted query string to get docking_score,
            leff of first pose for each ligand

        No Longer Returned:
            String: DB-formatted query string

        Raises:
            NotImplementedError: Description
        """
        raise NotImplementedError

    def _fetch_passing_plot_data(self):
        """Fetches cursor for best energies and leffs for ligands
            passing filtering

        No Longer Returned:
            SQLite cursor: Cursor containing docking_score,
                leff for the first pose for passing ligands

        Raises:
            NotImplementedError: Description
        """
        raise NotImplementedError

    def _generate_plot_passing_results_query(self):
        """Make DB-formatted query string to get docking_score,
            leff of first pose for passing ligands

        No Longer Returned:
            String: DB-formatted query string

        Raises:
            NotImplementedError: Description
        """
        raise NotImplementedError

    def _generate_result_filtering_query(
        self, results_filters_list, ligand_filters_list, output_fields
    ):
        """takes lists of filters, writes sql filtering string

        Args:
            results_filters_list (list): list of tuples where
                (filter column/key, filtering cutoff)
            ligand_filters_list (list): list of filters on ligand information
            output_fields (list): List of result column data to for output

        No Longer Returned:
            String: DB-formatted string for filtering query

        No Longer Raises:
            KeyError: Raises KeyError if user requests result ordering by
                invalid or multiple options
        """
        raise NotImplementedError

    def _generate_interaction_index_filtering_query(self, interaction_list):
        """takes list of interaction info for a given ligand, looks up
            corresponding interaction index

        Args:
            interaction_list (List): List containing interaction info
            in format
            [<interaction_type>, <rec_chain>, <rec_resname>,
            <rec_resid>, <rec_atom>]

        No Longer Returned:
            String: DB-formated query on Interaction_indices table

        Raises:
            NotImplementedError: Description
        """
        raise NotImplementedError

    def _generate_interaction_filtering_query(self, interaction_index_list):
        """takes list of interaction indices and searches for ligand ids which
            have those interactions

        Args:
            interaction_index_list (list): List of interaction indices

        No Longer Returned:
            String: DB-formatted query

        Raises:
            NotImplementedError: Description
        """
        raise NotImplementedError

    def _generate_ligand_filtering_query(self, ligand_filters):
        """write string to select from ligand table

        Args:
            ligand_filters (list): List of filters on ligand table

        No Longer Returned:
            String: DB-formatted query

        Raises:
            NotImplementedError: Description
        """
        raise NotImplementedError

    def _generate_results_data_query(self, output_fields):
        """Generates SQLite-formatted query string to select outfields data
            for ligands in self.results_view_name

        Args:
            output_fields (List): List of result column data for output

        Raises:
            NotImplementedError: Description
        """
        raise NotImplementedError

    def _generate_interaction_bitvectors(self, interactions_list):
        """takes string of interactions and makes bitvector

        Args:
            interactions_list (list): list of list of tuples. Inner lists
            contain interaction tuples for the saved poses for a single ligand

        No Longer Returned:
            List: List of bitvectors for saved poses

        Raises:
            NotImplementedError: Description
        """
        raise NotImplementedError

    def _insert_interaction_bitvectors(self, bitvectors):
        """Insert bitvectors of interaction data into database

        Args:
            bitvectors (List): List of lists With inner list representing
                interaction bitvector for a pose

        Raises:
            NotImplementedError: Description
        """
        raise NotImplementedError

    def _generate_percentile_rank_window(self):
        """makes window with percentile ranks for percentile filtering

        No Longer Returned:
            String: DB-formatted string for creating percent ranks on energies
                binding and leff

        Raises:
            NotImplementedError: Description
        """
        raise NotImplementedError

    def _calc_percentile_cutoff(self, percentile, column="docking_score"):
        """Make query for percentile by calculating energy or leff cutoff

        Args:
            percentile (float): cutoff percentile
            column (str, optional): string indicating column for cutoff to be calculated for
        """
        raise NotImplementedError

    def _fetch_results_column_names(self):
        """Fetches list of string for column names in results table

        No Longer Returned:
            List: List of strings of results table column names

        Raises:
            NotImplementedError: Description
        """
        raise NotImplementedError

    def _delete_from_results(self):
        """Remove rows from results table if they did not pass filtering

        Raises:
            NotImplementedError: Description
        """
        raise NotImplementedError

    def _delete_from_ligands(self):
        """Remove rows from ligands table if they did not pass filtering

        Raises:
            NotImplementedError: Description
        """
        raise NotImplementedError

    def _delete_from_interactions(self):
        """Remove rows from interactions bitvector table
        if they did not pass filtering

        Raises:
            NotImplementedError: Description
        """
        raise NotImplementedError

    def _generate_view_names_query(self):
        """Generate string to return names of views in database

        Raises:
            NotImplementedError: Description
        """
        raise NotImplementedError

    def _attach_db(self, new_db, new_db_name):
        """Attaches new database file to current database

        Args:
            new_db (string): file name for database to attach
            new_db_name (TYPE): Description

        Raises:
            NotImplementedError: Description
        """
        raise NotImplementedError

    def _detach_db(self, new_db_name):
        """Detaches new database file from current database

        Args:
            new_db_name (string): db name for database to detach

        Raises:
            NotImplementedError: Description
        """
        raise NotImplementedError

    def _create_temp_table(self, table_name):
        """create temporary table with given name

        Args:
            table_name (string): name for temp table
        """
        raise NotImplementedError

    def _generate_selective_insert_query(
        self, bookmark1_name, bookmark2_name, select_str, new_db_name, temp_table
    ):
        """Generates string to select ligands found/not found in the given bookmark in both current db and new_db

        Args:
            bookmark1_name (string): name of bookmark to cross-reference for main db
            bookmark2_name (string): name of bookmark to cross-reference for attached db
            select_str (string): "IN" or "NOT IN" indicating if ligand names should or should not be in both databases
            new_db_name (str): name of attached db
            temp_table (str): name of temporary table to store passing results in
        """
        raise NotImplementedError

    def _insert_into_temp_table(self, query):
        """Execute insertion into temporary table

        Args:
            query (str): Insertion command
        """
        raise NotImplementedError

    def _vacuum(self):
        raise NotImplementedError

    def check_storage_empty(self):
        """Check that storage is empty before proceeding.

        Raises:
            StorageError: Description
        """
        raise NotImplementedError


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


class StorageManagerSQLite(StorageManager):
    """SQLite-specific StorageManager subclass

    Attributes:
        conn (SQLite connection): Connection to database
        energy_filter_sqlite_call_dict (dictionary): Dictionary for
            translating filter options
        interactions_initialized_flag (bool): flag indicating if interaction tables intitialized
        next_unique_interaction_idx (int): idx for next interaction to be inserted into Interaction_indices table
        open_cursors (list): list of cursors that were not closed by the function that created them

    """

    def __init__(
        self,
        db_file: str = "output.db",
        append_results: bool = False,
        order_results: str = None,
        outfields: str = "Ligand_name,e",
        filter_bookmark: str = None,
        output_all_poses: bool = None,
        mfpt_cluster: float = None,
        interaction_cluster: float = None,
        results_view_name: str = "passing_results",
        overwrite: bool = None,
        conflict_opt: str = None,
        _stop_at_defaults=False,
    ):
        """Initialize superclass and subclass-specific instance variables

        Args:
            db_file (str): database file name
            opts (dict, optional): Dictionary of database options
            # TODO update this
        """
        self.append_results = append_results
        self.order_results = order_results
        if "Ligand_name" not in outfields:  # make sure we are outputting the ligand name
            outfields = "Ligand_name," + outfields
        self.outfields = outfields
        self.output_all_poses = output_all_poses
        self.mfpt_cluster = mfpt_cluster
        self.interaction_cluster = interaction_cluster
        self.filter_bookmark = filter_bookmark
        self.results_view_name = results_view_name
        self.overwrite = overwrite
        self.conflict_opt = conflict_opt
        self.db_file = db_file
        if _stop_at_defaults:
            return

        super().__init__()

        self.outfield_options = [
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
        ]
        self.order_options = {
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

        self.field_to_column_name = {
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
        self.interaction_name_to_letter = {
            "vdw_interactions": "V", 
            "hb_interactions": "H", 
            "reactive_interactions": "R",
            }

        self.view_suffix = None

        self.temptable_suffix = 0

        self.filtering_window = "Results"

        self.index_columns = []

        # keep track of any open cursors
        self.open_cursors = []

        self.energy_filter_sqlite_call_dict = {
            "eworst": "docking_score < {value}",
            "ebest": "docking_score > {value}",
            "leworst": "leff < {value}",
            "lebest": "leff > {value}",
        }

        self.energy_filter_col_name = {
            "eworst": "docking_score",
            "ebest": "docking_score",
            "leworst": "leff",
            "lebest": "leff",
            "score_percentile": "docking_score",
            "le_percentile": "leff",
        }

    @classmethod
    def get_defaults(cls):
        return cls(_stop_at_defaults=True).__dict__

    @classmethod
    def get_default_types(cls):
        return typing.get_type_hints(cls.__init__)

    # # # # # # # # # # # # # # # # #
    # # # # #Public methods # # # # #
    # # # # # # # # # # # # # # # # #
    @classmethod
    def format_for_storage(cls, ligand_dict: dict) -> tuple:
        """takes file dictionary from the file parser, formats required storage format

        Args:
            ligand_dict (dict): Dictionary containing data from the fileparser

        Returns:
            tuple: Tuple of lists ([result_row_1, result_row_2,...],
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

    @classmethod
    def _generate_results_row(cls, ligand_dict, pose_rank, run_number):
        """generate list of lists of ligand values to be
            inserted into sqlite database

        Args:
            ligand_dict (Dictionary): Dictionary of ligand data from parser
            pose_rank (int): Rank of pose to generate row for
                all runs for the given ligand
            run_number (int): Run number of pose to generate row for
                all runs for the given ligand

        Returns:
            List: List of pose data to be inserted into Results table.
            In same order as expected in _insert_results:
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
        for key in cls.ligand_data_keys:
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
            # Add the cluster size for the cluster this pose belongs to
            ligand_data_list.append(
                ligand_dict["cluster_sizes"][ligand_dict["cluster_list"][pose_rank]]
            )
        else:
            ligand_data_list.extend(
                [
                    None,
                    None,
                    None,
                ]
            )
        # add statevars
        for key in cls.stateVar_keys:
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

    @classmethod
    def _generate_ligand_row(cls, ligand_dict):
        """writes row to be inserted into ligand table

        Args:
            ligand_dict (Dictionary): Dictionary of ligand data from parser

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

    @classmethod
    def _generate_receptor_row(cls, ligand_dict):
        """Writes row to be inserted into receptor table

        Args:
            ligand_dict (Dictionary): Dictionary of ligand data from parser
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

    @classmethod
    def _generate_interaction_tuples(cls, interaction_dictionaries):
        """takes dictionary of file results, formats as
        list of tuples for interactions

        Args:
            interaction_dictionaries (List): List of pose interaction
            dictionaries from parser

        Returns:
            List: List of tuples of interaction data
        """
        interactions = set()
        for pose_interactions in interaction_dictionaries:
            count = pose_interactions["count"][0]
            for i in range(int(count)):
                interactions.add(
                    tuple(pose_interactions[kw][i] for kw in cls.interaction_data_kws)
                )

        return list(interactions)

    def insert_data(
        self,
        results_array,
        ligands_array,
        interaction_array,
        receptor_array=[],
        insert_receptor=False,
    ):
        """Summary

        Args:
            results_array (list): list of data to be stored in Results table
            ligands_array (list): list of data to be stored in Ligands table
            interactions_array (list): list of data to be stored in interaction tables
            receptor_array (list): list of data to be stored in Receptors table
            insert_receptor (bool, optional): flag indicating that receptor info should inserted
        """
        self.insert_results(results_array)
        self.insert_ligands(ligands_array)
        if insert_receptor and receptor_array != []:
            self.insert_receptors(receptor_array)
        if interaction_array != []:
            self.insert_interactions(interaction_array)

    def insert_results(self, results_array):
        """takes array of database rows to insert, adds data to results table

        Args:
            results_array (numpy array): numpy array of arrays containing
                formatted result rows

        Raises:
            DatabaseInsertionError: Description

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
        ) VALUES \
        (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)"""

        try:
            cur = self.conn.cursor()
            cur.executemany(sql_insert, results_array)
            self.conn.commit()
            cur.close()

        except sqlite3.OperationalError as e:
            raise DatabaseInsertionError("Error while inserting results.") from e

    def insert_ligands(self, ligand_array):
        """Takes array of ligand rows, inserts into Ligands table.

        Args:
            ligand_array (numpy array): Numpy array of arrays
                containing formatted ligand rows

        Raises:
            DatabaseInsertionError: Description

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

    def insert_receptors(self, receptor_array):
        """Takes array of receptor rows, inserts into Receptors table

        Args:
            receptor_array (list): List of lists
                containing formatted ligand rows

        Raises:
            DatabaseInsertionError: Description
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

    def _fetch_existing_interactions(self):
        """return cursor of interactions in Interaction_indices table

        Returns:
            sqlite cursor: cursor of
        """
        query = """SELECT interaction_id, interaction_type, rec_chain, rec_resname, rec_resid, rec_atom, rec_atomid from Interaction_indices"""
        return self._run_query(query)

    def insert_interactions(self, interactions_list):
        """Takes list of interactions, inserts into database

        Args:
            interactions_list (list): List of tuples for interactions in form
            ("type", "chain", "residue", "resid", "recname", "recid")
        """
        # populate unique interactions from existing database
        if self.append_results:
            existing_unique_interactions = self._fetch_existing_interactions()
            for interaction in existing_unique_interactions:
                self.unique_interactions[interaction[1:]] = interaction[0]

            self.next_unique_interaction_idx = (
                interaction[0] + 1
            )  # sets next index for next unique interaction

        self._add_unique_interactions(interactions_list)

        # check if we need to initialize the interaction bv table and
        # insert first set of interaction
        if not self.interactions_initialized_flag and not self.append_results:
            self._create_interaction_bv_table()
            self._insert_unique_interactions(list(self.unique_interactions.keys()))
            self.interactions_initialized_flag = True

        self._insert_interaction_bitvectors(
            self._generate_interaction_bitvectors(interactions_list)
        )

    def save_receptor(self, receptor):
        """Takes object of Receptor class, updates the column in Receptor table

        Args:
            receptor (bytes): bytes receptor object to be inserted into DB

        Deleted Parameters:
            rec_name (string): Name of receptor. Used to insert into correct row of DB

        Raises:
            DatabaseInsertionError: Description
        """

        sql_update = (
            """UPDATE Receptors SET receptor_object = ? WHERE Receptor_ID == 1"""
        )

        try:
            cur = self.conn.cursor()
            cur.execute(sql_update, (receptor,))
            self.conn.commit()
            cur.close()

        except sqlite3.OperationalError as e:
            raise DatabaseInsertionError(
                "Error while adding receptor blob to database"
            ) from e

    def fetch_receptor_object_by_name(self, rec_name):
        """Returns Receptor object from database for given rec_name

        Args:
            rec_name (string): Name of receptor to return object for
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
            rec_name (string): Name of receptor to return object for
        """

        cursor = self._run_query(
            "SELECT RecName, receptor_object FROM Receptors"
        )
        return cursor.fetchall()

    def clone(self, backup_name=None):
        """Creates a copy of the db"""
        if backup_name is None:
            backup_name = self.db_file + ".bk"
        bck = sqlite3.connect(backup_name)
        with bck:
            self.conn.backup(bck, pages=1)
        bck.close()

    def set_ringtaildb_version(self):
        rt_version = pkg_resources.get_distribution("ringtail").version.replace(".", "")
        cur = self.conn.cursor()
        cur.execute(f"PRAGMA user_version = {rt_version}")
        self.conn.commit()
        cur.close()

    def check_ringtaildb_version(self):
        cur = self.conn.cursor()
        db_version = str(cur.execute("PRAGMA user_version").fetchone()[0])
        cur.close()
        return db_version == pkg_resources.get_distribution("ringtail").version.replace(".", ""), db_version

    def update_database(self, consent=False):
        cur = self.conn.cursor()
        # get views and drop them
        if not consent:
            logging.warning("WARNING: All existing bookmarks in database will be dropped during database update!")
            consent = input("Type 'yes' if you wish to continue: ") == 'yes'
        if not consent:
            logging.critical("Consent not given for database update. Cancelling...")
            sys.exit(1)

        logging.info(f"Updating {self.db_file}...")
        views = cur.execute("SELECT name FROM sqlite_master WHERE type = 'view'").fetchall()
        for v in views:
            cur.execute(f"DROP VIEW IF EXISTS {v[0]}")
        # delete all rows in bookmarks table
        cur.execute("DELETE FROM Bookmarks")

        # reformat for v1.1.0
        cur.execute("ALTER TABLE Results RENAME COLUMN energies_binding TO docking_score")
        cur.execute("ALTER TABLE Bookmarks ADD COLUMN filters")
        cur.execute("CREATE INDEX allind ON Results(LigName, docking_score, leff, deltas, reference_rmsd, energies_inter, energies_vdw, energies_electro, energies_intra, nr_interactions, run_number, pose_rank, num_hb)")
        self.conn.commit()
        cur.close()

        self.set_ringtaildb_version()

        return consent

    def remake_bookmarks(self):
        """Reads all views from Bookmarks table and remakes them"""
        try:
            bookmark_info = self._run_query("SELECT * from Bookmarks")
            for bookmark_name, query in bookmark_info:
                if "Comparision. Wanted:" in query:  # cannot remake comparison bookmark
                    continue
                self._create_view(bookmark_name, query)
        except sqlite3.OperationalError as e:
            raise StorageError("Error while remaking views") from e

    def get_results(self):
        """Gets all fields for filtered results

        Returns:
            SQLite cursor: Cursor with all fields
                and rows in passing results view
        """
        # check if we have previously filtered and saved view
        return self._run_query(
            "SELECT * FROM {passing_view}".format(passing_view=self.results_view_name)
        )

    def get_maxmiss_union(self, total_combinations: int):
        """"""
        selection_strs = []
        view_strs = []
        outfield_str = self._generate_outfield_string()
        for i in range(total_combinations):
            selection_strs.append(f"SELECT {outfield_str} FROM {self.results_view_name + '_' + str(i)}")
            view_strs.append(f"SELECT * FROM {self.results_view_name + '_' + str(i)}")

        view_name = f"{self.results_view_name}_union"
        logging.debug("Saving union bookmark...")
        self._create_view(view_name, " UNION ".join(view_strs))
        self._insert_bookmark_info(view_name, " UNION ".join(view_strs))
        logging.debug("Running union query...")
        return self._run_query(" UNION ".join(selection_strs))

    def fetch_summary_data(self, columns=["docking_score", "leff"], percentiles=[1,10]) -> dict:
        """Collect summary data for database:
            Num Ligands
            Num stored poses
            Num unique interactions
            
            min, max, percentiles for columns in columns
        """
        try:
            summary_data = {}
            cur = self.conn.cursor()
            summary_data["num_ligands"] = cur.execute("SELECT COUNT(*) FROM Ligands").fetchone()[0]
            summary_data["num_poses"] = cur.execute("SELECT COUNT(*) FROM Results").fetchone()[0]
            summary_data["num_unique_interactions"] = cur.execute("SELECT COUNT(*) FROM Interaction_indices").fetchone()[0]
            summary_data["num_interacting_residues"] = cur.execute("SELECT COUNT(*) FROM (SELECT interaction_id FROM Interaction_indices GROUP BY interaction_type,rec_resid,rec_chain)").fetchone()[0]

            allowed_columns = self._fetch_results_column_names()
            for col in columns:
                if col not in allowed_columns:
                    raise StorageError(f"Requested summary column {col} not found in Results table! Available columns: {allowed_columns}")
                summary_data[f"min_{col}"] = cur.execute(f"SELECT MIN({col}) FROM Results").fetchone()[0]
                summary_data[f"max_{col}"] = cur.execute(f"SELECT MAX({col}) FROM Results").fetchone()[0]
                for p in percentiles:
                    summary_data[f"{p}%_{col}"] = self._calc_percentile_cutoff(p, col)

            return summary_data
                
        except sqlite3.OperationalError as e:
            raise StorageError("Error while fetching summary data!") from e

    def get_number_passing_ligands(self, bookmark_name=None):
        """Returns count of the number of ligands that
            passed filtering criteria

        Returns:
            Int: Number of passing ligands

        Raises:
            DatabaseQueryError: Description
        """
        if bookmark_name is None:
            bookmark_name = self.current_view_name
        try:
            cur = self.conn.cursor()
            cur.execute(
                "SELECT COUNT(DISTINCT LigName) FROM {results_view}".format(
                    results_view=bookmark_name
                )
            )
            n_ligands = int(cur.fetchone()[0])
            cur.close()
            return n_ligands
        except sqlite3.OperationalError as e:
            raise DatabaseQueryError(
                "Error while getting number of passing ligands"
            ) from e

    def fetch_flexres_info(self):
        """fetch flexres names and atomname lists

        Returns:
            Tuple: (flexible_residues, flexres_atomnames)
        """
        try:
            cur = self.conn.cursor()
            cur.execute("SELECT flexible_residues, flexres_atomnames FROM Receptors")
            info = cur.fetchone()
            cur.close()
            return info
        except sqlite3.OperationalError as e:
            raise DatabaseQueryError("Error retrieving flexible residue info") from e

    def fetch_passing_ligand_output_info(self):
        """fetch information required by vsmanager for writing out molecules

        Returns:
            SQLite cursor: contains LigName, ligand_smile,
                atom_index_map, hydrogen_parents
        """
        query = "SELECT LigName, ligand_smile, atom_index_map, hydrogen_parents FROM Ligands WHERE LigName IN (SELECT DISTINCT LigName FROM passing_temp)"
        return self._run_query(query)

    def fetch_single_ligand_output_info(self, ligname):
        try:
            cur = self.conn.cursor()
            cur.execute(f"SELECT LigName, ligand_smile, atom_index_map, hydrogen_parents FROM Ligands WHERE LigName LIKE '{ligname}'")
            info = cur.fetchone()
            cur.close()
            return info
        except sqlite3.OperationalError as e:
            raise DatabaseQueryError(f"Error retrieving ligand info for {ligname}") from e

    def fetch_clustered_similars(self, ligname: str):
        """Given ligname, returns poseids for similar poses/ligands from previous clustering. User prompted at runtime to choose cluster.

        Args:
            ligname (str): ligname for ligand to find similarity with
        """
        logging.warning("N.B.: When finding similar ligands, export tasks (i.e. SDF export) will be for the selected similar ligands, NOT ligands passing given filters.")
        cur = self.conn.cursor()

        ligand_cluster_columns = self._fetch_ligand_cluster_columns()
        print("Here are the existing clustering groups. Please ensure that you query ligand(s) is a part of the group you select.")
        print("   Choice number   |   Underlying filter bookmark   |   Morgan or interaction fingerprint?   |   cutoff   ")
        print("----------------------------------------------------------------------------------------------------------")
        for i, col in enumerate(ligand_cluster_columns):
            col_info = col.split("_")
            option_list = [str(i)] + ["_".join(col_info[:-2])] + [col_info[-2]] + [col_info[-1].replace("p", ".")]
            print(f"{'    |    '.join(option_list)}")
        cluster_choice = input("Please specify choice number for the cluster you would like to return similar ligands from: ")
        try:
            cluster_col_choice = ligand_cluster_columns[int(cluster_choice)]
        except ValueError:
            raise ValueError(f"Given cluster number {cluster_choice} cannot be converted to int. Please be sure you are specifying integer.")
       
        query_ligand_cluster = cur.execute(f"SELECT {cluster_col_choice} FROM Ligand_clusters WHERE pose_id IN (SELECT Pose_ID FROM Results WHERE LigName LIKE '{ligname}')").fetchone()
        if query_ligand_cluster is None:
            raise DatabaseQueryError(f"Requested ligand name {ligname} not found in cluster {cluster_col_choice}!")
        query_ligand_cluster = query_ligand_cluster[0]  # extract from tuple
        sql_query = f"SELECT LigName FROM Results WHERE Pose_ID IN (SELECT pose_id FROM Ligand_clusters WHERE {cluster_col_choice}={query_ligand_cluster}) GROUP BY LigName"
        view_query = f"SELECT * FROM Results WHERE Pose_ID IN (SELECT pose_id FROM Ligand_clusters WHERE {cluster_col_choice}={query_ligand_cluster}) GROUP BY LigName"

        view_name = f"similar_{ligname}_{cluster_col_choice}"
        self._create_view(view_name, view_query)
        self._insert_bookmark_info(name=view_name, sqlite_query=sql_query)

        self.results_view_name = view_name

        return self._run_query(sql_query), view_name, cluster_col_choice


    def create_temp_passing_table(self):
        cur = self.conn.cursor()
        cur.execute(f"CREATE TEMP TABLE passing_temp AS SELECT * FROM {self.results_view_name}")
        cur.close()

    def fetch_passing_pose_properties(self, ligname):
        """fetch coordinates for poses passing filter for given ligand

        Args:
            ligname (string): name of ligand to fetch coordinates for

        Returns:
            SQLite cursor: contains Pose_ID, docking_score, leff, ligand_coordinates,
                flexible_res_coordinates, flexible_residues
        """
        query = "SELECT Pose_ID, docking_score, leff, ligand_coordinates, flexible_res_coordinates FROM Results WHERE Pose_ID IN (SELECT Pose_ID FROM passing_temp WHERE LigName LIKE '{ligand}')".format(ligand=ligname)
        return self._run_query(query)

    def fetch_nonpassing_pose_properties(self, ligname):
        """fetch coordinates for poses of ligname which did not pass the filter

        Args:
            ligname (string): name of ligand to fetch coordinates for

        Returns:
            SQLite cursor: contains Pose_ID, docking_score, leff, ligand_coordinates,
                flexible_res_coordinates, flexible_residues
        """
        query = "SELECT Pose_ID, docking_score, leff, ligand_coordinates, flexible_res_coordinates FROM Results WHERE LigName LIKE '{ligand}' AND Pose_ID NOT IN (SELECT Pose_ID FROM passing_temp)".format(
            ligand=ligname,
        )
        return self._run_query(query)

    def fetch_single_pose_properties(self, pose_ID:int):
        """fetch coordinates for pose given by pose_ID

        Args:
            pose_ID (int): name of ligand to fetch coordinates for

        Returns:
            SQLite cursor: contains Pose_ID, docking_score, leff, ligand_coordinates,
                flexible_res_coordinates, flexible_residues
        """
        query = f"SELECT Pose_ID, docking_score, leff, ligand_coordinates, flexible_res_coordinates FROM Results WHERE Pose_ID={pose_ID}"
        return self._run_query(query)

    def fetch_interaction_bitvector(self, pose_id):
        """Returns tuple containing interaction bitvector line for given pose_id

        Args:
            pose_id (int): pose id to fetch interaction bitvector for

        Returns:
            tuple: tuple representing interaction bitvector
            None: if no interactions in database
        """
        # catch if database does not have interactions
        table_names = [table[0] for table in self._fetch_existing_table_names()]
        if "Interaction_bitvectors" not in table_names:
            return None

        query = "SELECT * FROM Interaction_bitvectors WHERE Pose_ID = {0}".format(
            pose_id
        )
        return self._run_query(query).fetchone()[1:]  # cut off pose id

    def fetch_interaction_info_by_index(self, interaction_idx):
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

    def get_current_view_name(self):
        """returns current view name

        Returns:
            string: name of last passing results view used by database
        """
        return self.current_view_name

    def check_receptors_saved(self):
        """returns number of rows in Receptors table where receptor_object already has blob

        Returns:
            int: number of rows in receptors table

        Raises:
            DatabaseQueryError: Description
        """
        try:
            cur = self.conn.cursor()
            cur.execute("SELECT COUNT(*) FROM Receptors WHERE receptor_object NOT NULL")
            row_count = cur.fetchone()[0]
            cur.close()
            return row_count
        except sqlite3.OperationalError as e:
            raise DatabaseQueryError(
                "Error occurred while fetching number of receptor rows containing PDBQT blob"
            ) from e

    def to_dataframe(self, requested_data: str, table=True) -> pd.DataFrame:
        """Returns dataframe of table or query given as requested_data

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

    def fetch_view(self, viewname: str) -> sqlite3.Cursor:
        """returns SQLite cursor of all fields in viewname

        Args:
            viewname (str): name of view to retrieve

        Returns:
            sqlite3.Cursor: cursor of requested view
        """
        return self._run_query(f"SELECT * FROM {viewname}")

    def save_temp_table(
        self,
        temp_table_name,
        bookmark_name,
        original_bookmark_name,
        wanted_list,
        unwanted_list=[],
    ):
        """Resaves temp bookmark stored in self.current_view_name as new permenant bookmark

        Args:
            bookmark_name (string): name of bookmark to save last temp bookmark as
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
        compare_bookmark_str = "Comparision. Wanted: "
        compare_bookmark_str += ", ".join(wanted_list)
        if unwanted_list is not None:
            compare_bookmark_str += ". Unwanted: " + ", ".join(unwanted_list)
        self._insert_bookmark_info(bookmark_name, compare_bookmark_str)

    # # # # # # # # # # # # # # # # #
    # # # # #Private methods # # # # #
    # # # # # # # # # # # # # # # # #

    def _create_connection(self):
        """Creates database connection to self.db_file

        Returns:
            SQLite connection: Connection object to self.db_file

        Raises:
            DatabaseConnectionError: Description

        """
        try:
            con = sqlite3.connect(self.db_file)
            try:
                con.enable_load_extension(True)
                con.load_extension("chemicalite")
                con.enable_load_extension(False)
            except sqlite3.OperationalError as e:
                logging.critical("Failed to load chemicalite cartridge. Please ensure chemicalite is installed with `conda install -c conda-forge chemicalite`.")
                raise e
            cursor = con.cursor()
            cursor.execute("PRAGMA synchronous = OFF")
            cursor.execute("PRAGMA journal_mode = MEMORY")
            cursor.close()
        except sqlite3.OperationalError as e:
            raise DatabaseConnectionError(
                "Error while establishing database connection"
            ) from e
        return con

    def _close_connection(self):
        """Closes connection to database"""
        logging.info("Closing database")
        self.conn.close()

    def _close_open_cursors(self):
        """closes any cursors stored in self.open_cursors.
        Resets self.open_cursors to empty list
        """
        for cur in self.open_cursors:
            cur.close()

        self.open_cursors = []

    def open_storage(self):
        """Create connection to db. Then, check if db needs to be written.
        If so, (if self.overwrite drop existing tables and )
        initialize the tables
        """
        
        self.conn = self._create_connection()
        
        # register signal handler to catch keyboard interupts
        signal(SIGINT, self._sigint_handler)

        # if we want to overwrite old db, drop existing tables
        if self.overwrite:
            self._drop_existing_tables()
        # create tables in db
        self._create_results_table()
        self._create_ligands_table()
        self._create_receptors_table()
        self._create_interaction_index_table()
        self._create_bookmark_table()

    def _fetch_existing_table_names(self):
        """Returns list of all tables in database

        Returns:
            list: list of table names

        Raises:
            DatabaseQueryError: Description
        """

        try:
            cur = self.conn.cursor()
            cur.execute("SELECT name FROM sqlite_schema WHERE type='table';")
            return cur.fetchall()
        except sqlite3.OperationalError as e:
            raise DatabaseQueryError(
                "Error while getting names of existing database tables"
            ) from e

    def _drop_existing_tables(self):
        """drop any existing tables.
        Will only be called if self.overwrite is true

        Raises:
            StorageError: Description
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

    def _drop_existing_views(self):
        """Drop any existing views
        Will only be called if self.overwrite is true

        Raises:
            StorageError: Description
        """
        # fetch existing views
        try:
            cur = self.conn.cursor()
            cur.execute("SELECT name FROM sqlite_schema WHERE type='view';")
            views = cur.fetchall()

            # drop views
            for view in views:
                # cannot drop this, so we catch it instead
                if view[0] == "sqlite_sequence":
                    continue
                cur.execute("DROP VIEW {view_name}".format(view_name=view[0]))
            cur.close()
        except sqlite3.OperationalError as e:
            raise StorageError(
                "Error occured while dropping existing database views"
            ) from e

    def _run_query(self, query):
        """Executes provided SQLite query. Returns cursor for results.
            Since cursor remains open, added to list of open cursors

        Args:
            query (string): Formated SQLite query as string

        Returns:
            SQLite cursor: Contains results of query
        """
        try:
            cur = self.conn.cursor()
            cur.execute(query)
            self.open_cursors.append(cur)
        except sqlite3.OperationalError as e:
            raise DatabaseQueryError("Unable to execute query {0}".format(query)) from e
        return cur

    def create_indices(self):
        """Create index containing possible filter and order by columns
        """
        try:
            cur = self.conn.cursor()
            logging.debug("Creating columns index...")
            cur.execute("CREATE INDEX allind ON Results(LigName, docking_score, leff, deltas, reference_rmsd, energies_inter, energies_vdw, energies_electro, energies_intra, nr_interactions, run_number, pose_rank, num_hb)")
            self.conn.commit()
            cur.close()
        except sqlite3.OperationalError as e:
            raise StorageError("Error occurred while indexing") from e

    def _remove_indices(self):
        """Removes idx_filter_cols and idx_ligname"""
        try:
            cur = self.conn.cursor()
            cur.execute("DROP INDEX IF EXISTS idx_filter_cols")
            cur.execute("DROP INDEX IF EXISTS idx_ligname")
            cur.close()
        except sqlite3.OperationalError as e:
            raise StorageError("Error while dropping indices") from e

    def _create_view(self, name, query, temp=False, add_poseID=False):
        """takes name and selection query,
            creates view of query stored as name.

        Args:
            name (string): Name for view which will be created
            query (string): SQLite-formated query used to create view
            temp (bool, optional): Flag if view should be temporary
            add_poseID (bool, optional): Add Pose_ID column to view

        Raises:
            DatabaseViewCreationError: Description
        """
        # check that bookmark does not start with int, this causes a sqlite error
        if name[0].isdigit():
            raise DatabaseViewCreationError(f"Bookmark names may not start with digit. Given bookmark name {name}.")
        cur = self.conn.cursor()
        if add_poseID:
            query = query.replace("SELECT ", "SELECT Pose_ID, ", 1)
        logging.info("Creating bookmark...")
        # drop old view if there is one
        try:
            if temp:
                temp_flag = "TEMP "
            else:
                temp_flag = ""
            cur.execute("DROP VIEW IF EXISTS {name}".format(name=name))
            cur.execute(
                "CREATE {temp_flag}VIEW {name} AS {query}".format(
                    name=name, query=query, temp_flag=temp_flag
                )
            )
            logging.debug(
                "CREATE {temp_flag}VIEW {name} AS {query}".format(
                    name=name, query=query, temp_flag=temp_flag
                )
            )
            cur.close()
        except sqlite3.OperationalError as e:
            raise DatabaseViewCreationError(
                "Error creating view from query \n{0}".format(query)
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
        flexible_residues   VARCHAR[],
        flexible_res_coordinates   VARCHAR[]

        Raises:
            DatabaseTableCreationError: Description
        """
        unique_string = ""
        if self.conflict_opt is not None:
            unique_string = """, UNIQUE(LigName, receptor, about_x, about_y, about_z,
                   trans_x, trans_y, trans_z,
                   axisangle_x, axisangle_y, axisangle_z, axisangle_w,
                   dihedrals, flexible_residues) ON CONFLICT {0}""".format(
                self.conflict_opt
            )

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
            {0}
        );
        """.format(
            unique_string
        )

        try:
            cur = self.conn.cursor()
            cur.execute(sql_results_table)
            cur.close()
        except sqlite3.OperationalError as e:
            raise DatabaseTableCreationError(
                "Error while creating results table. If database already exists, use --overwrite to drop existing tables"
            ) from e

    def _create_receptors_table(self):
        """Create table for receptors. Columns are:
        Receptor_ID         INTEGER PRIMARY KEY AUTOINCREMENT,
        RecName             VARCHAR,
        box_dim             VARCHAR[],
        box_center          VARCHAR[],
        grid_spacing        INT[],
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
            grid_spacing        INT[],
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

    def _create_interaction_index_table(self):
        """create table of data for each unique interaction. Columns are:
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
        interaction_index_table = """CREATE TABLE IF NOT EXISTS Interaction_indices (
            interaction_id      INTEGER PRIMARY KEY AUTOINCREMENT,
            interaction_type    VARCHAR[],
            rec_chain           VARCHAR[],
            rec_resname         VARCHAR[],
            rec_resid           VARCHAR[],
            rec_atom            VARCHAR[],
            rec_atomid          VARCHAR[])"""

        try:
            cur = self.conn.cursor()
            cur.execute(interaction_index_table)
            cur.close()
        except sqlite3.OperationalError as e:
            raise DatabaseTableCreationError(
                "Error while creating interaction index table. If database already exists, use --overwrite to drop existing tables"
            ) from e

    def _create_interaction_bv_table(self):
        """Create table of interaction bits for each pose. Columns are:
        Pose_ID INTERGER PRIMARY KEY AUTOINCREMENT
        Interaction_1
        Interaction_2
        ...
        Interaction_n

        Raises:
            DatabaseTableCreationError: Description

        """
        interact_columns_str = (
            " INTEGER,\n".join(
                [
                    "Interaction_" + str(i + 1)
                    for i in range(len(self.unique_interactions))
                ]
            )
            + " INTEGER"
        )

        bv_table = """CREATE TABLE IF NOT EXISTS Interaction_bitvectors (
        Pose_ID INTEGER PRIMARY KEY AUTOINCREMENT,
        {columns})""".format(
            columns=interact_columns_str
        )

        try:
            cur = self.conn.cursor()
            cur.execute(bv_table)
            cur.close()
        except sqlite3.OperationalError as e:
            raise DatabaseTableCreationError(
                "Error while creating interaction bitvector table. If database already exists, use --overwrite to drop existing tables"
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

    def _insert_bookmark_info(self, name: str, sqlite_query: str, filters = {}):
        """Insert bookmark info into bookmark table

        Args:
            name (str): name for bookmark
            sqlite_query (str): sqlite query used to generate bookmark
            filters (dict): filters used to generate bookmark

        Raises:
            DatabaseInsertionError: Description
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

    def _insert_unique_interactions(self, unique_interactions):
        """Inserts interaction data for unique interactions
            into Interaction_index table

        Args:
            unique_interactions (list): List of tuples of interactions
                to be inserted

        Raises:
            DatabaseInsertionError: Description
        """
        sql_insert = """INSERT INTO Interaction_indices (
        interaction_type,
        rec_chain,
        rec_resname,
        rec_resid,
        rec_atom,
        rec_atomid
        ) VALUES (?,?,?,?,?,?)"""

        try:
            cur = self.conn.cursor()
            cur.executemany(sql_insert, unique_interactions)
            self.conn.commit()
            cur.close()

        except sqlite3.OperationalError as e:
            raise DatabaseInsertionError(
                "Error while inserting unique interactions into interaction index table"
            ) from e

    def _insert_one_interaction(self, interaction):
        """Insert interaction data for a single new interaction
            into the interaction indices table

        Args:
            interaction (tuple): Tuple of interaction data
                (interaction_type, rec_chain, rec_resname,
                rec_resid, rec_atom, rec_atomid)

        Raises:
            DatabaseInsertionError: Description

        """
        sql_insert = """INSERT INTO Interaction_indices (
        interaction_type,
        rec_chain,
        rec_resname,
        rec_resid,
        rec_atom,
        rec_atomid
        ) VALUES (?,?,?,?,?,?)"""

        try:
            cur = self.conn.cursor()
            cur.execute(sql_insert, interaction)
            self.conn.commit()
            cur.close()

        except sqlite3.OperationalError as e:
            raise DatabaseInsertionError(
                "Error inserting interaction {interact} into interaction index table".format(
                    interact=str(interaction)
                )
            ) from e

    def _make_new_interaction_column(self, column_number):
        """Add column for new interaction to interaction bitvector table

        Args:
            column_number (int): Index for new interaction

        Raises:
            StorageError: Description
        """
        add_column_str = """ALTER TABLE Interaction_bitvectors ADD COLUMN Interaction_{n_inter}""".format(
            n_inter=str(column_number)
        )
        try:
            cur = self.conn.cursor()
            cur.execute(add_column_str)
            self.conn.commit()
            cur.close()

        except sqlite3.OperationalError as e:
            raise StorageError(
                "Error adding column for Interaction_{num} to interaction bitvector table".format(
                    num=str(column_number)
                )
            ) from e

    #TODO combine functions here
    def _fetch_all_plot_data(self):
        """Fetches cursor for best energies and leff for all ligands

        Returns:
            SQLite Cursor: Cursor containing docking_score,
                leff for the first pose for each ligand
        """
        return self._run_query(self._generate_plot_all_results_query())

    def _generate_plot_all_results_query(self):
        """Make SQLite-formatted query string to get docking_score,
            leff of first pose of all ligands

        Returns:
            String: SQLite-formatted query string
        """
        return "SELECT docking_score, leff FROM Results GROUP BY LigName"

    def _fetch_passing_plot_data(self):
        """Fetches cursor for best energies and leffs for
            ligands passing filtering

        Returns:
            SQLite cursor: Cursor containing docking_score,
                leff for the first pose for passing ligands
        """
        return self._run_query(self._generate_plot_passing_results_query())

    def _generate_plot_passing_results_query(self):
        """Make SQLite-formatted query string to get docking_score,
            leff of first pose for passing ligands

        Returns:
            String: SQLite-formatted query string
        """
        return "SELECT docking_score, leff, Pose_ID, LigName FROM Results WHERE LigName IN (SELECT DISTINCT LigName FROM {results_view}) GROUP BY LigName".format(
            results_view=self.results_view_name
        )

    def _generate_outfield_string(self):

        # parse requested output fields and convert to column names in database
        outfields_list = self.outfields.split(",")
        for outfield in outfields_list:
            if outfield not in self.outfield_options:
                raise OptionError(
                    "{out_f} is not a valid output option. Please see rt_process_vs.py --help for allowed options".format(
                        out_f=outfield
                    )
                )
        return", ".join(
            [self.field_to_column_name[field] for field in outfields_list]
        )

    def _generate_result_filtering_query(
        self,
        filters_dict
    ):
        """takes lists of filters, writes sql filtering string

        Args:
            filters_dict (dict): dict of filters. Keys names and value formats must match those found in the Filters class

        Returns:
            String: SQLite-formatted string for filtering query
        """
        # before we do anything, check that the DB version matches the version number of our module
        rt_version_same, db_rt_version = self.check_ringtaildb_version()
        if not rt_version_same:
            # TODO: will cause error when any version int is > 10
            # catch version 1.0.0 where returned db_rt_version will be 0
            if db_rt_version == 0:
                db_rt_version = 100
            raise StorageError(f"Input database was created with Ringtail v{'.'.join([i for i in db_rt_version[:2]] + [db_rt_version[2:]])}. Confirm that this matches current Ringtail version and use Ringtail update script(s) to update database if needed.")

        outfield_string = self._generate_outfield_string()

        if self.filter_bookmark is not None:
            if self.filter_bookmark == self.results_view_name:
                logging.error(
                    f"Specified filter_bookmark and bookmark_name are the same: {self.results_view_name}"
                )
                raise OptionError(
                    "--filter_bookmark and --bookmark_name cannot be the same! Please rename --bookmark_name"
                )
            self.filtering_window = self.filter_bookmark

        # write energy filters and compile list of interactions to search for
        queries = []
        interaction_filters = []

        for filter_key, filter_value in filters_dict.items():
            if filter_value is None:
                continue
            if filter_key in self.energy_filter_col_name:
                self.index_columns.append(self.energy_filter_col_name[filter_key])
                if filter_key == "score_percentile" or filter_key == "le_percentile":
                    # convert from percent to decimal
                    cutoff = self._calc_percentile_cutoff(filter_value, self.energy_filter_col_name[filter_key])
                    queries.append( f"{self.energy_filter_col_name[filter_key]} < {cutoff}"
                    )
                else:
                    queries.append(
                        self.energy_filter_sqlite_call_dict[filter_key].format(
                            value=filter_value
                        )
                    )

            # write hb count filter(s)
            if filter_key == "interactions_count":
                for k,v in filter_value:
                    # TODO implement other interaction count filters
                    if k != "hb_count":
                        continue
                    self.index_columns.append("num_hb")
                    if v > 0:
                        queries.append("num_hb > {value}".format(value=v))
                    else:
                        queries.append("num_hb <= {value}".format(value=-1 * v))

            # reformat interaction filters as list
            if filter_key in Filters.get_interaction_filter_keys():
                for interact in filter_value:
                    interaction_string = filter_key + ":" + interact[0]
                    interaction_filters.append(
                        interaction_string.split(":") + [interact[1]]
                    )  # add bool flag for included (T) or excluded (F) interaction

            # add react_any flag as interaction filter
            # check if react_any is true
            if filter_key == "react_any" and filter_value:
                interaction_filters.append(["reactive_interactions", "", "", "", "", True])

        # for each interaction filter, get the index
        # from the interactions_indices table
        interaction_queries = []
        for interaction in interaction_filters:
            interaction = [self.interaction_name_to_letter[interaction[0]]] + interaction[1:]
            interaction_filter_indices = []
            interact_index_str = self._generate_interaction_index_filtering_query(
                interaction[:-1]
            )  # remove bool include/exclude flag
            interaction_indices = self._run_query(interact_index_str)
            for i in interaction_indices:
                interaction_filter_indices.append(i[0])

            # catch if interaction not found in results
            if interaction_filter_indices == []:
                if interaction == ["R", "", "", "", "", True]:
                    logging.warning(
                        "Given --react_any filter, no reactive interactions found. Excluded from filtering."
                    )
                else:
                    logging.warning(
                        "Interaction {i} not found in results, excluded from filtering".format(
                            i=":".join(interaction[:4])
                        )
                    )
                continue
            # determine include/exclude string
            if interaction[-1] is True:
                include_str = "IN"
            elif interaction[-1] is False:
                include_str = "NOT IN"
            else:
                raise RuntimeError(
                    "Unrecognized flag in interaction. Please contact Forli Lab with traceback and context."
                )
            # find pose ids for ligands with desired interactions
            interaction_queries.append(
                "Pose_ID {include_str} ({interaction_str})".format(
                    include_str=include_str,
                    interaction_str=self._generate_interaction_filtering_query(
                        interaction_filter_indices
                    ),
                )
            )

        # add ligand filters
        ligand_filters_dict = {k:v for k, v in filters_dict.items() if k in Filters.get_ligand_filter_keys()}
        if filters_dict["ligand_substruct"] != [] or filters_dict["ligand_name"] != []:
            ligand_query_str = self._generate_ligand_filtering_query(
                        ligand_filters_dict
                    )
            queries.append(
                "LigName IN ({ligand_str})".format(
                    ligand_str=ligand_query_str
                )
            )
        if len(ligand_filters_dict["ligand_substruct_pos"]):
            nr_args_per_group = 6
            nr_smarts = int(len(ligand_filters_dict["ligand_substruct_pos"]) / nr_args_per_group)
            # create temporary table with molecules that pass all smiles
            tmp_lig_filters = {"ligand_operator": ligand_filters_dict["ligand_operator"]}
            if "ligand_max_atoms" in ligand_filters_dict:
                tmp_lig_filters["ligand_max_atoms"] = ligand_filters_dict["ligand_max_atoms"]
            tmp_lig_filters["ligand_substruct"] = [ligand_filters_dict["ligand_substruct_pos"][i * nr_args_per_group] for i in range(nr_smarts)]
            cmd = self._generate_ligand_filtering_query(tmp_lig_filters)
            cmd = cmd.replace(
                "SELECT LigName FROM Ligands",
                "SELECT "
                    "Results.Pose_ID, "
                    "Ligands.LigName, "
                    "Ligands.ligand_smile, "
                    "Ligands.atom_index_map, "
                    "Results.ligand_coordinates "
                    "FROM Ligands INNER JOIN Results ON Results.LigName = Ligands.LigName"
            )
            cmd = "CREATE TEMP TABLE passed_smarts AS " + cmd
            cur = self.conn.cursor()
            cur.execute("DROP TABLE IF EXISTS passed_smarts")
            cur.execute(cmd)
            smarts_loc_filters = []
            for i in range(nr_smarts):
                smarts =      ligand_filters_dict["ligand_substruct_pos"][i * nr_args_per_group + 0]
                index =   int(ligand_filters_dict["ligand_substruct_pos"][i * nr_args_per_group + 1])
                sqdist = float(ligand_filters_dict["ligand_substruct_pos"][i * nr_args_per_group + 2])**2
                x =     float(ligand_filters_dict["ligand_substruct_pos"][i * nr_args_per_group + 3])
                y =     float(ligand_filters_dict["ligand_substruct_pos"][i * nr_args_per_group + 4])
                z =     float(ligand_filters_dict["ligand_substruct_pos"][i * nr_args_per_group + 5])
                # save filter for bookmark
                smarts_loc_filters.append((smarts, index, x, y, z))
                poses = self._run_query("SELECT * FROM passed_smarts")
                pose_id_list = []
                smartsmol = Chem.MolFromSmarts(smarts)
                for pose_id, ligname, smiles, idxmap, coords in poses:
                    mol = Chem.MolFromSmiles(smiles)
                    idxmap = [int(value)-1 for value in json.loads(idxmap)]
                    idxmap = {idxmap[j*2]: idxmap[j*2+1] for j in range(int(len(idxmap)/2))}
                    for hit in mol.GetSubstructMatches(smartsmol):
                        xyz = [float(value) for value in json.loads(coords)[idxmap[hit[index]]]]
                        d2 = (xyz[0] - x)**2 + (xyz[1] - y)**2 + (xyz[2] - z)**2
                        if d2 <= sqdist:
                            pose_id_list.append(str(pose_id))
                            break # add pose only once
                queries.append("Pose_ID IN ({0})".format(",".join(pose_id_list)))
            cur.close()
        # format query string
        # raise error if query string is empty
        if queries == [] and interaction_queries == []:
            raise DatabaseQueryError(
                "Query strings are empty. Please check filter options and ensure requested interactions are present."
            )

        sql_string = output_str = """SELECT {out_columns} FROM {window} WHERE """.format(
            out_columns=outfield_string, window=self.filtering_window
        )
        if interaction_queries == []:
            joined_queries = " AND ".join(queries)
            sql_string = sql_string + joined_queries
            unclustered_query = f"SELECT Pose_id FROM {self.filtering_window} WHERE " + joined_queries
        else:
            with_stmt = f"WITH subq as (SELECT Pose_id FROM {self.filtering_window}) "
            if queries != []:
                with_stmt = with_stmt[:-2] + f" WHERE {' AND '.join(queries)}) "
            joined_interact_queries = " AND ".join(interaction_queries)
            sql_string = with_stmt + sql_string + joined_interact_queries
            unclustered_query = f"SELECT Pose_id FROM {self.filtering_window} WHERE " + joined_interact_queries

        # adding if we only want to keep
        # one pose per ligand (will keep first entry)
        if not self.output_all_poses:
            sql_string += " GROUP BY LigName"

        # add how to order results
        if self.order_results is not None:
            try:
                sql_string += (
                    " ORDER BY " + self.field_to_column_name[self.order_results]
                )
            except KeyError:
                raise RuntimeError(
                    "Please ensure you are only requesting one option for --order_results and have written it correctly"
                ) from None

        # if clustering is requested, do that before saving view or filtering results for output
        #Define clustering setup
        def clusterFps(fps, cutoff):  #https://www.macinchem.org/reviews/clustering/clustering.php

            # first generate the distance matrix:
            dists = []
            nfps = len(fps)
            for i in range(1,nfps):
                sims = DataStructs.BulkTanimotoSimilarity(fps[i],fps[:i])
                dists.extend([1-x for x in sims])

            # now cluster the data:
            cs = Butina.ClusterData(dists,nfps,cutoff,isDistData=True)
            return cs

        if self.interaction_cluster is not None:
            logging.warning("WARNING: Interaction fingerprint clustering is memory-constrained. Using overly-permissive filters with clustering may cause issues.")# TODO: remove this memory bottleneck
            cluster_query = f"SELECT Results.leff, Interaction_bitvectors.* FROM Interaction_bitvectors INNER JOIN Results ON Results.Pose_ID = Interaction_bitvectors.Pose_id WHERE Results.Pose_ID IN ({unclustered_query})"
            if interaction_queries != []:
                cluster_query = with_stmt + cluster_query
            leff_poseid_ifps = self._run_query(cluster_query).fetchall()
            def make_bitstring(pose_bv):
                bs = ""
                for i in pose_bv:
                    if i is None:
                        bs += "0"
                    elif i == 1:
                        bs += "1"
                    else:
                        raise RuntimeError(f"Unrecognized character {i} in interaction bitvector.")
                return bs

            bclusters = clusterFps([DataStructs.CreateFromBitString(make_bitstring(pose[2:])) for pose in leff_poseid_ifps], self.interaction_cluster)
            logging.info(f"Number of interaction fingerprint butina clusters: {len(bclusters)}")

            # select ligand from each cluster with best ligand efficiency
            int_rep_poseids = []
            for c in bclusters:
                c_leffs = np.array([leff_poseid_ifps[i][0] for i in c])  # beware magic numbers
                best_lig_c = leff_poseid_ifps[c[np.argmin(c_leffs)]][1]
                int_rep_poseids.append(str(best_lig_c))

            self._insert_cluster_data(bclusters, [l[1] for l in leff_poseid_ifps], "ifp", str(self.interaction_cluster))

            # catch if no pose_ids returned
            if int_rep_poseids == []:
                logging.warning("No passing results prior to clustering. Clustering not performed.")
            else:
                if self.mfpt_cluster is None:
                    sql_string = output_str + "Pose_ID=" + " OR Pose_ID=".join(int_rep_poseids)
                else:
                    unclustered_query = f"SELECT Pose_ID FROM Results WHERE {'Pose_ID=' + ' OR Pose_ID='.join(int_rep_poseids)}"

        if self.mfpt_cluster is not None:
            logging.warning("WARNING: Ligand morgan fingerprint clustering is memory-constrained. Using overly-permissive filters with clustering may cause issues.")# TODO: remove this memory bottleneck
            logging.warning("N.B.: If using both interaction and morgan fingerprint clustering, the morgan fingerprint clustering will be performed on the results staus post interaction fingerprint clustering.")
            cluster_query = f"SELECT Results.Pose_ID, Results.leff, mol_morgan_bfp(Ligands.ligand_rdmol, 2, 1024) FROM Ligands INNER JOIN Results ON Results.LigName = Ligands.LigName WHERE Results.Pose_ID IN ({unclustered_query})"
            if interaction_queries != []:
                cluster_query = with_stmt + cluster_query
            poseid_leff_mfps = self._run_query(cluster_query).fetchall()
            bclusters = clusterFps([DataStructs.CreateFromBinaryText(mol[2]) for mol in poseid_leff_mfps], self.mfpt_cluster)
            logging.info(f"Number of Morgan fingerprint butina clusters: {len(bclusters)}")
            
            # select ligand from each cluster with best ligand efficiency
            fp_rep_poseids = []
            for c in bclusters:
                c_leffs = np.array([poseid_leff_mfps[i][1] for i in c])
                best_lig_c = poseid_leff_mfps[c[np.argmin(c_leffs)]][0]
                fp_rep_poseids.append(str(best_lig_c))

            self._insert_cluster_data(bclusters, [l[0] for l in poseid_leff_mfps], "mfp", str(self.mfpt_cluster))

            # catch if no pose_ids returned
            if fp_rep_poseids == []:
                logging.warning("No passing results prior to clustering. Clustering not performed.")
            else:
                sql_string = output_str + "Pose_ID=" + " OR Pose_ID=".join(fp_rep_poseids)

        return sql_string, sql_string.replace("""SELECT {out_columns} FROM {window}""".format(
            out_columns=outfield_string, window=self.filtering_window
        ), f"SELECT * FROM {self.filtering_window}")  # sql_query, view_query

    def _fetch_ligand_cluster_columns(self):
        try:
            return [c[1] for c in self._run_query("PRAGMA table_info(Ligand_clusters)").fetchall()][1:]
        except IndexError:
            raise IndexError("Error fetching columns from Ligand_clusters table. Confirm that ligand clustering has been previously performed.")

    def _insert_cluster_data(self, clusters: list, poseid_list: list, cluster_type: str, cluster_cutoff: str):
        cur = self.conn.cursor()
        cur.execute("CREATE TABLE IF NOT EXISTS Ligand_clusters (pose_id  INT[] UNIQUE)")
        ligand_cluster_columns = self._fetch_ligand_cluster_columns()
        column_name = f"{self.results_view_name}_{cluster_type}_{cluster_cutoff.replace('.', 'p')}"
        if column_name not in ligand_cluster_columns:
            cur.execute(f"ALTER TABLE Ligand_clusters ADD COLUMN {column_name}")
        for ci, cl in enumerate(clusters):
            for i in cl:
                poseid = poseid_list[i]
                cur.execute(f"INSERT INTO Ligand_clusters (pose_id, {column_name}) VALUES (?,?) ON CONFLICT (pose_id) DO UPDATE SET {column_name}=excluded.{column_name}", (poseid, ci))

        cur.close()
        self.conn.commit()
    
    def _generate_interaction_index_filtering_query(self, interaction_list):
        """takes list of interaction info for a given ligand,
            looks up corresponding interaction index

        Args:
            interaction_list (List): List containing interaction info
                in format [<interaction_type>, <rec_chain>, <rec_resname>,
                <rec_resid>, <rec_atom>]

        Returns:
            String: SQLite-formated query on Interaction_indices table
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

        return sql_string

    def _generate_interaction_filtering_query(self, interaction_index_list):
        """takes list of interaction indices and searches for ligand ids
            which have those interactions

        Args:
            interaction_index_list (list): List of interaction indices

        Returns:
            String: SQLite-formatted query
        """

        return "SELECT Pose_id FROM (SELECT * FROM Interaction_bitvectors WHERE Pose_ID IN subq) WHERE " + " OR ".join(
            [
                "Interaction_{index_n} = 1".format(index_n=index)
                for index in interaction_index_list
            ]
        )

    def _generate_ligand_filtering_query(self, ligand_filters):
        """write string to select from ligand table

        Args:
            ligand_filters (list): List of filters on ligand table

        Returns:
            String: SQLite-formatted query, Dict: dictionary of filters and values
        """

        sql_ligand_string = "SELECT LigName FROM Ligands WHERE"
        logical_operator = ligand_filters["ligand_operator"]
        if logical_operator is None:
            logical_operator = "AND"
        for kw in ligand_filters.keys():
            fils = ligand_filters[kw]
            if kw == "ligand_name":
                for name in fils:
                    if name == "":
                        continue
                    name_sql_str = " LigName LIKE '%{value}%' OR".format(value=name)
                    sql_ligand_string += name_sql_str
            if kw == "ligand_max_atoms" and ligand_filters[kw] is not None:
                maxatom_sql_str = " mol_num_atms(ligand_rdmol) <= {} {}".format(ligand_filters[kw], logical_operator)
                sql_ligand_string += maxatom_sql_str
            if kw == "ligand_substruct":
                for smarts in fils:
                    # check for hydrogens in smarts pattern
                    smarts_mol = Chem.MolFromSmarts(smarts)
                    for atom in smarts_mol.GetAtoms():
                        if atom.GetAtomicNum() == 1:
                            raise DatabaseQueryError(f"Given ligand substructure filter {smarts} contains explicit hydrogens. Please re-run query with SMARTs without hydrogen.")
                    substruct_sql_str = " mol_is_substruct(ligand_rdmol, mol_from_smarts('{smarts}')) {logical_operator}".format(
                        smarts=smarts, logical_operator=logical_operator)
                    sql_ligand_string += substruct_sql_str
        if sql_ligand_string.endswith("AND"):
            sql_ligand_string = sql_ligand_string.rstrip("AND")
        if sql_ligand_string.endswith("OR"):
            sql_ligand_string = sql_ligand_string.rstrip("OR")
        
        return sql_ligand_string

    def _generate_results_data_query(self, output_fields):
        """Generates SQLite-formatted query string to select outfields data for ligands in self.results_view_name

        Args:
            output_fields (List): List of result column data for output

        Returns:
            str: string to select data from passing results view
        """
        outfield_string = "LigName, " + ", ".join(
            [self.field_to_column_name[field] for field in output_fields]
        )

        return (
            "SELECT "
            + outfield_string
            + " FROM Results WHERE Pose_ID IN (SELECT Pose_ID FROM {0})".format(
                self.results_view_name
            )
        )

    def _generate_interaction_bitvectors(self, interactions_list):
        """takes string of interactions and makes bitvector

        Args:
            interactions_list (list): list of list of tuples.
                Inner lists contain interaction tuples for
                saved poses for a single ligand

        Returns:
            List: List of bitvectors for saved poses
        """
        bitvectors_list = []
        for pose_interactions in interactions_list:
            pose_bitvector = [None] * len(self.unique_interactions)
            for interaction_tuple in pose_interactions:
                interaction_idx = self.unique_interactions[interaction_tuple]
                # index corrected for interaction_indices starting at 1
                pose_bitvector[interaction_idx - 1] = 1

            bitvectors_list.append(pose_bitvector)
        return bitvectors_list

    def _insert_interaction_bitvectors(self, bitvectors):
        """Insert bitvectors of interaction data into database

        Args:
            bitvectors (List): List of lists With inner list representing
                interaction bitvector for a pose

        Raises:
            DatabaseInsertionError: Description
        """
        interaction_columns = range(len(self.unique_interactions))
        column_str = ""
        filler_str = ""
        for i in interaction_columns:
            column_str += "Interaction_" + str(i + 1) + ", "
            filler_str += "?,"
        column_str = column_str.rstrip(", ")
        filler_str = filler_str.rstrip(",")

        # make sure we have unique interactions, otherwise insert empty rows for those poses
        if self.unique_interactions == {}:
            sql_insert = """INSERT INTO Interaction_bitvectors DEFAULT VALUES"""
        else:
            sql_insert = """INSERT INTO Interaction_bitvectors ({columns}) VALUES ({fillers})""".format(
                columns=column_str, fillers=filler_str
            )

        try:
            cur = self.conn.cursor()
            cur.executemany(sql_insert, bitvectors)
            self.conn.commit()
            cur.close()

        except sqlite3.OperationalError as e:
            raise DatabaseInsertionError("Error while inserting bitvectors") from e

    def _generate_percentile_rank_window(self):
        """makes window with percentile ranks for percentile filtering

        Returns:
            String: SQLite-formatted string for creating
                percent ranks on docking_score and leff
        """
        column_names = ",".join(self._fetch_results_column_names())
        return "SELECT {columns}, PERCENT_RANK() OVER (ORDER BY docking_score) score_percentile_rank, PERCENT_RANK() OVER (ORDER BY leff) leff_percentile_rank FROM Results Group BY LigName".format(
            columns=column_names
        )

    def _calc_percentile_cutoff(self, percentile: float, column="docking_score"):
        """Make query for percentile by calculating energy or leff cutoff

        Args:
            percentile (float): cutoff percentile
            column (str, optional): string indicating column for percentile to be calculated over
        """
        # get total number of ligands
        try:
            logging.debug(f"Generating percentile filter query for {column}")
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
            logging.debug(f"{column} percentile cutoff is {cutoff}")
            return cutoff
        except sqlite3.OperationalError as e:
            raise StorageError("Error while generating percentile query") from e

    def _fetch_results_column_names(self):
        """Fetches list of string for column names in results table

        Returns:
            List: List of strings of results table column names

        Raises:
            StorageError: Description
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

    def _delete_from_results(self):
        """Remove rows from results table if they did not pass filtering

        Raises:
            StorageError: Description
        """
        try:
            cur = self.conn.cursor()
            cur.execute(
                "DELETE FROM Results WHERE Pose_ID NOT IN (SELECT Pose_ID FROM {view})".format(
                    view=self.results_view_name
                )
            )
            self.conn.commit()
            cur.close()
        except sqlite3.OperationalError as e:
            raise StorageError(
                f"Error occured while pruning Results not in {self.results_view_name}"
            ) from e

    def _delete_from_ligands(self):
        """Remove rows from ligands table if they did not pass filtering

        Raises:
            StorageError: Description
        """
        try:
            cur = self.conn.cursor()
            cur.execute(
                "DELETE FROM Ligands WHERE LigName NOT IN (SELECT LigName from Results WHERE Pose_ID IN (SELECT Pose_ID FROM {view}))".format(
                    view=self.results_view_name
                )
            )
            self.conn.commit()
            cur.close()
        except sqlite3.OperationalError as e:
            raise StorageError(
                f"Error occured while pruning Ligands not in {self.results_view_name}"
            ) from e

    def _delete_from_interactions(self):
        """Remove rows from interactions bitvector table
        if they did not pass filtering

        Raises:
            StorageError: Description
        """
        try:
            cur = self.conn.cursor()
            cur.execute(
                "DELETE FROM Interaction_bitvectors WHERE Pose_ID NOT IN (SELECT Pose_ID FROM {view})".format(
                    view=self.results_view_name
                )
            )
            self.conn.commit()
            cur.close()
        except sqlite3.OperationalError as e:
            raise StorageError(
                f"Error occured while pruning Interaction_bitvectors not in {self.results_view_name}"
            ) from e

    def _generate_view_names_query(self):
        """Generate string to return names of views in database"""
        return "SELECT name FROM sqlite_schema WHERE type = 'view'"

    def _attach_db(self, new_db, new_db_name):
        """Attaches new database file to current database

        Args:
            new_db (string): file name for database to attach
            new_db_name (str): name of new database

        Raises:
            StorageError: Description
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
            new_db_name (string): db name for database to detach
        """
        detach_str = f"DETACH DATABASE {new_db_name}"

        try:
            cur = self.conn.cursor()
            cur.execute(detach_str)
            self.conn.commit()
            cur.close()
        except sqlite3.OperationalError as e:
            raise StorageError(f"Error occurred while detaching {new_db_name}") from e

    def _create_temp_table(self, table_name):
        """create temporary table with given name

        Args:
            table_name (string): name for temp table
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

    def _generate_selective_insert_query(
        self, bookmark1_name, bookmark2_name, select_str, new_db_name, temp_table
    ):
        """Generates string to select ligands found/not found in the given bookmark in both current db and new_db

        Args:
            bookmark1_name (string): name of bookmark to cross-reference for main db
            bookmark2_name (string): name of bookmark to cross-reference for attached db
            select_str (string): "IN" or "NOT IN" indicating if ligand names should or should not be in both databases
            new_db_name (str): name of attached db
            temp_table (str): name of temporary table to store passing results in
        """
        return "INSERT INTO {0} SELECT Pose_ID, LigName FROM {1} WHERE LigName {2} (SELECT LigName FROM {3}.{4})".format(
            temp_table, bookmark1_name, select_str, new_db_name, bookmark2_name
        )

    def _insert_into_temp_table(self, query):
        """Execute insertion into temporary table

        Args:
            query (str): Insertion command
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

    def _vacuum(self):
        try:
            cur = self.conn.cursor()
            cur.execute("VACUUM")
            self.conn.commit()
            cur.close()
        except sqlite3.OperationalError as e:
            raise DatabaseInsertionError(f"Error while vacuuming DB") from e

    def check_storage_empty(self):
        """Check that storage is empty before proceeding.

        Raises:
            StorageError: Description
        """
        if not self.overwrite and not self.append_results:
            try:
                cur = self.conn.cursor()
                cur.execute("SELECT COUNT(*) FROM Results")
                if cur.fetchone()[0] != 0:
                    raise StorageError(
                        "Database already exists. Use --overwrite or --append_results if wanting to replace or append to existing database."
                    )
            except Exception as e:
                raise e
            finally:
                cur.close()
