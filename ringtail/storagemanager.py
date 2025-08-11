#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail storage adaptors
#

import time
import json
import os.path
import pandas as pd
from .logutils import LOGGER as logger
import sys
import re
from signal import signal, SIGINT
from rdkit import Chem
from typing import Union
import time
from importlib.metadata import version
from .ringtailoptions import Filters
from .util import numlist2str, raise_not_implemented, statuses
from .exceptions import (
    StorageError,
    DatabaseInsertionError,
    DatabaseConnectionError,
    DatabaseTableCreationError,
    DatabaseQueryError,
    OptionError,
    MergeError,
)
from .clustermanager import *
from .querybuilder import QueryBuilder, QueryBuilderSQLite
from collections import defaultdict

try:
    import sqlite3

    HAS_SQLITE = True
except ImportError:
    HAS_SQLITE = False


class StorageManager:
    # NOTE gotta be careful with schema
    _db_schema_ver = "3.0.0"
    QueryBuilder = QueryBuilder

    """Base class for a generic virtual screening database object.
    This class holds some of the common API for StorageManager child classes. 
    Each child class will implement their own functions to write to and read from the database

    Attributes: 
        _db_schema_ver (str): current database schema version
        _db_schema_code_compatibility (dict): dictionary showing compatibility of code base versions with relational database schema versions
        active_connection (bool): if there is a current open connection to database
        keyboard_interrupt_allowed (bool): if Ctrl+Z will work, for example not supported in a GUI
    """

    # region database access
    def __init__(self):
        self.keyboard_interrupt_allowed = False
        self.active_connection = True

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
        if self.active_connection:
            self.close_storage()
            self.active_connection = False
        if exc_type:
            logger.error(str(exc_value))
            raise exc_value from None
        return self

    def _sigint_handler(self, signal_received, frame):
        """Handles and reports if program is interrupted through the terminal"""
        logger.critical("Ctrl + C pressed, keyboard interupt initiated")
        self.__exit__(None, None, None)
        sys.exit(0)

    # endregion

    # region public api

    @classmethod
    def format_for_storage(cls, ligand_dict: dict) -> dict:
        # TODO the only sql specific about this method is the methods it uses
        """Takes file dictionary from the file parser, formats required storage format. Only handles docking data for one ligand at the time.
        For each run we save, we add its interaction dict to the interaction_dictionaries list and save its other data
        We also save a mapping of the its cluster number to the index in interaction_dictionaries
        Then, when we find a pose to tolerate interactions for, we lookup the index to append the interactions to from cluster_saved_pose_map
        Finally, we calculate the interaction tuple lists for each pose

        Args:
            ligand_dict (dict): Dictionary containing data from the fileparser

        Returns:
            dict: with storage formatted rows, including: results rows per pose, interactions per pose, one ligand row, and one receptor row

        """
        result_rows = []
        interaction_dictionaries = []
        interaction_tuples = []
        saved_pose_idx = 0  # save index of last saved pose
        cluster_saved_pose_map = {}  # save mapping of cluster number to saved_pose_idx

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

    def insert_data(
        self,
        docking_data: dict,
        write_options: dict = {},
    ):
        """Inserts data from all arrays returned from results manager. Can have one or more ligands in docking_data

        Args:
            docking_data (dict): docking results to be inserted, where key is ligand name and value is data to be written
            write_options (dict): options for how to write the data, primarily how to treat duplicates
        """
        # parse ligand info form dict into list of ligands
        ligands_array = [
            docked_ligand["ligand_row"] for docked_ligand in docking_data.values()
        ]
        # get unique ligand_id (will not add duplicate, instead return existing ligand_id)
        ligand_ids = self._insert_ligands(ligands_array)

        # add ligand ids to results array and make result array list of poses
        results_array = []
        for index, docked_ligand in enumerate(docking_data.values()):
            for pose in docked_ligand["poses_results"]:
                results_array.append([ligand_ids[index]] + pose)
        # get unique pose ids and duplicate handling info
        Pose_IDs, duplicates = self._insert_results(results_array, write_options)
        # check if are interactions:
        if any(
            docked_ligand.get("poses_interactions") not in (None, [])
            for docked_ligand in docking_data.values()
        ):

            self.insert_interactions(
                Pose_IDs,
                docking_data,
                duplicates,
                write_options.get("duplicate_handling"),
            )

    def insert_receptor(self, receptor_data: list):
        """

        Method to insert receptor information into the database

        Args:
            receptor_data (list): of receptor data, ordered by columns in the db
        """
        receptors = self.fetch_receptor_object()
        # insert receptor if database does not have already have a receptor entry
        if not receptors:
            self._insert_receptors(receptor_data)

    def insert_interactions(
        self,
        pose_IDs: list,
        docking_data: dict,
        duplicates: list,
        duplicate_handling: str = None,
    ):
        """Looks for interactions in the docking data, and inserts new interaction definitions
        as well as interaction rows if there are any.

        Args:
            pose_IDs (list[int]): pose ids currently being processed
            docking_data (dict): all parsed data from the docking fror those poses
            duplicates (list(Pose_ID)): any duplicates identified in "insert_results", if duplicate handling has been specified
            duplicate_handling (str): how to treat any duplicates

        """
        if duplicate_handling:
            if duplicate_handling.lower() not in ["ignore", "replace"]:
                logger.warning(
                    f"duplicate_handling option {duplicate_handling} not allowed. Reverting to default behavior."
                )
                duplicate_handling = None
        interaction_rows = self._insert_and_format_interactions(pose_IDs, docking_data)
        # NOTE here as of dev meeting
        self._insert_interaction_rows(interaction_rows, duplicates, duplicate_handling)

    def insert_receptor_blob(self, receptor: bytes, rec_name: str):
        """Takes object of Receptor class, updates the column in Receptor table

        Args:
            receptor (bytes): bytes receptor object to be inserted into DB
            rec_name (string): Name of receptor. Used to insert into correct row of DB
        """
        # TODO make virtual
        self._insert_receptor_blob(receptor, rec_name)

    def filter_results(
        self,
        all_filters: Filters,
        bookmark_name: str,
        filtering_bookmark: str = None,
    ) -> int:
        """Generate and execute database queries from given filters.

        Args:
            all_filters (dict): dict containing all filters
            bookmark_name (str): name of bookmark in which to save the data
            filtering_bookmark (str, optional): if filtering not across all data, but a pre-filtered bookmark

        Returns:
             int: number of passing ligands
        """
        # before we do anything, check that the DB version matches the version number of our module
        rt_version_same, db_rt_version = self.check_ringtaildb_version()
        if not rt_version_same:
            # NOTE will cause error when any version int is > 10
            raise StorageError(
                f"Input database was created with Ringtail v{'.'.join([i for i in db_rt_version[:2]] + [db_rt_version[2:]])}. Confirm that this matches current Ringtail version and use Ringtail upgrade script to upgrade database if needed."
            )

        # get the final filter query, has a {selection} place holder
        filter_query: str = self._generate_result_filtering_query(
            all_filters, bookmark_name, filtering_bookmark
        )

        logger.debug(f"Query for filtering results: {filter_query}")

        # perform filtering
        logger.debug("Running filtering query...")
        time0 = time.perf_counter()
        if filtering_bookmark == None:
            filtering_bookmark = "Results"

        self._populate_filter_tables(
            name=bookmark_name,
            query=filter_query,
            filters=all_filters,
            filtering_bookmark=filtering_bookmark,
        )
        logger.debug(
            f"Time to filter results: {time.perf_counter() - time0:.2f} seconds"
        )
        count = self.get_passing_ligands_count(bookmark_name)

        return count

    def get_all_bookmark_names(self) -> list[str]:
        """Get all bookmarks in sql database as a list of names. Bookmarks are a Ringtail-specific saved query (much like views)

        Returns:
            list: of bookmark names
        """
        try:
            query = self.QueryBuilder()
            query.SELECT("name").FROM("Filters")
            cur = self.db_query(*query.build())
            bookmark_names = [row["name"].lower() for row in cur.fetchall()]

        except Exception as e:
            raise StorageError(
                "Error occured while fetching existing bookmark names"
            ) from e

        return bookmark_names

    def get_passing_poses_count(
        self, bookmark_name: str, grouped_by_ligand: bool = False
    ) -> int:
        """
        Count poses in bookmark

        Args:
            bookmark_name (str): bookmark name in which to count
            grouped_by_ligand (bool, optional): if grouping by ligand, essentially returns passing ligands

        Returns:
            int: number of poses (optionally grouped by ligand) in bookmark
        """
        query = self.QueryBuilder()
        query.SELECT("r.pose_id").FROM("results", "r").IN_BOOKMARK(bookmark_name)
        if grouped_by_ligand:
            query.GROUP_BY("r.ligand_id")
        query_string = query.build(count=True)[0]
        if self.db_query(query_string).fetchone():
            return self.db_query(query_string).fetchone()[0]
        else:
            return 0

    def get_passing_ligands_count(self, bookmark_name: str) -> int:
        """
        Get number of passing ligands in bookmark name

        Args:
            bookmark_name (str): bookmark that defines passing

        Returns:
            int: number of ligands
        """
        return self.get_passing_poses_count(bookmark_name, True)

    def delete_bookmark(self, bookmark_name: str):
        """
        Deletes bookmark (i.e., Filters table) and its associated poses (filtered_poses table)

        Args:
            bookmark_name (str): bookmark to delete
        """
        # get filter id
        query = self.QueryBuilder()
        query.SELECT("filter_id").FROM("Filters").WHERE("name=?", bookmark_name)

        filter_id = self.db_query(*query.build()).fetchone()[0]
        # delete from filters
        query = self.QueryBuilder()
        self.db_query(
            *query.DELETE_FROM("Filters").WHERE("filter_id = ?", filter_id).build()
        )
        # delete from filtered_poses table
        query = self.QueryBuilder()
        self.db_query(
            *query.DELETE_FROM("Filtered_poses")
            .WHERE("filter_id = ?", filter_id)
            .build(),
            commit=True,
        )
        logger.info(
            f"The bookmark {bookmark_name} and its associated filter data has been deleted."
        )

    def get_ligname_from_pose(self, pose_id: int) -> str:
        """
        Get ligand name given a pose_id

        Args:
            pose_id (int): pose id for which to get ligand name

        Returns:
            str: ligand name
        """
        query = self.QueryBuilder()
        query.SELECT("L.LigName").FROM("Ligands", "L").JOIN(
            "Results", "R", "ligand_id"
        ).WHERE("R.pose_id =?", pose_id)
        return self.db_query(*query.build()).fetchone()[0]

    def get_maxmiss_union(
        self, total_combinations: int, bookmark_name: str, all_filters={}
    ):
        """Get results that are in union considering max miss

        Args:
            total_combinations (int): numer of possible combinations

        Returns:
            iter: of passing results
        """
        enumerated_bookmarks = []
        existing_bookmarks = self.get_all_bookmark_names()
        for i in range(total_combinations):
            bmn = bookmark_name + "_" + str(i)
            if bmn in existing_bookmarks:
                enumerated_bookmarks.append(f"'{bmn}'")

        subq = self.QueryBuilder()
        subq.SELECT("filter_id").FROM("Filters").WHERE(
            f"name IN ({', '.join(enumerated_bookmarks)})"
        )
        query = self.QueryBuilder()
        query.SELECT("DISTINCT pose_id").FROM("filtered_poses").WHERE(
            f"filter_id IN ({subq.build()[0]})"
        )
        bookmark_name = f"{bookmark_name}_union"
        logger.debug("Saving union bookmark...")
        logger.debug("Running union query...")
        self._populate_filter_tables(
            name=bookmark_name, query=query.build()[0], filters=all_filters
        )

        count = self.get_passing_ligands_count(bookmark_name)

        if not count:
            bookmark_name = None
        return count, bookmark_name

    def cluster_data(
        self,
        bookmark_name: str,
        cluster_type: str = "mfpt",
        cutoff: float = 0.5,
    ) -> tuple:
        """
        Clusters data in a given bookmark. Will create a new bookmark with the format
        <bookmark_name>_<cluster-type>_clustered containing the representative poses
        for the clusters

        Args:
            bookmark_name (str): bookmark name with poses to cluster
            cluster_type (str, optional): type of clustering. Defaults to "mfpt".
            cutoff (float, optional): cutoff cluster distance. Defaults to 0.5.

        Returns:
            tuple (str, int): clustered bookmark name, number of clusters
        """
        logger.debug("Preparing to cluster")
        time0 = time.perf_counter()

        query = self.QueryBuilder()
        # TODO handle if table is Results
        if cluster_type.lower() == "ifp":
            query.SELECT("r.pose_id", "r.leff").FROM("Results", "R")

            if self.is_bookmark(bookmark_name):
                query.IN_BOOKMARK(bookmark_name)
            elif self._is_statustable(bookmark_name):
                query.JOIN(bookmark_name, "T", "pose_id")
            # tuple, tuple
            pose_ids, leffs = zip(*self.db_query(query.build()[0]).fetchall())
            leffs = list(leffs)
            pose_id_bitvectors = self._generate_interaction_bitvectors(pose_ids)
            pose_ids = list(pose_ids)
            bitvectors = list(pose_id_bitvectors.values())

            ibc = InteractionBitvectorCluster(pose_ids, leffs, bitvectors, cutoff)
            clusters, representatives = ibc.cluster()

        elif cluster_type.lower() == "mfp":
            query.SELECT("r.pose_id", "r.leff", "l.rdmol").FROM("Results", "R").JOIN(
                "ligands", "l", "ligand_id", "results"
            )

            if self.is_bookmark(bookmark_name):
                query.IN_BOOKMARK(bookmark_name)
            elif self._is_statustable(bookmark_name):
                query.JOIN(bookmark_name, "T", "pose_id")

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
        cluster_bookmark_name = self._insert_cluster_data(
            clusters,
            pose_ids,
            cluster_type.lower(),
            str(cutoff),
            bookmark_name,
        )

        clustered_poses = self.QueryBuilder()
        clustered_poses.SELECT("pose_id").FROM("results").WHERE(
            f"pose_id IN ({','.join(representatives)})"
        )

        self._populate_filter_tables(
            cluster_bookmark_name,
            clustered_poses.build()[0],
            {"cluster_type": cluster_type, "cutoff": cutoff},
            bookmark_name,
        )

        logger.info(f"Time to cluster data: {time.perf_counter() - time0:.2f} seconds")

        return cluster_bookmark_name, len(clusters)

    # endregion

    # region virtual public api
    def close_storage(self, attached_db=None, vacuum=False):
        """Close connection to database

        Args:
            attached_db (str, optional): name of attached DB (not including file extension)
            vacuum (bool, optional): indicates that database should be vacuumed before closing
        """
        raise_not_implemented()

    def check_storage_ready(
        self, run_mode: str, docking_mode: str, store_all_poses: bool, max_poses: int
    ):
        """
        Check that storage is ready and compatible with options before proceeding, and creates new tables if needed

        Args:
            run_mode (str): if ringtail is ran using cmd line interface or api
            docking_mode (str): what docking engine was used to produce results
            store_all_poses (bool): overrwrites max poses
            max_poses (int): max poses to save to db
        """
        raise_not_implemented()

    def check_ringtaildb_version(self):
        """
        Checks the database version and confirms whether the code base is compatible with it

        Returns:
            bool: whether or not db is compatible with the code base
            str: current database version
        """
        raise_not_implemented()

    def overwrite_storage(self):
        """
        Will drop all tables in the database.
        """
        raise_not_implemented()

    def finalize_database_write(self):
        """
        Methods to finalize when a database has been written to, including creating indices
        """
        # index certain tables
        self._create_indices()
        logger.info("Database write session completed successfully.")

    def clone(self, backup_name=None):
        """Creates a copy of the db

        Args:
            backup_name (str, optional): name of the cloned database
        """
        raise_not_implemented()

    def update_database_version(self, new_version, consent=False):
        """method that updates sqlite database schema 1.0.0 or 1.1.0 to 1.1.0 or 2.0.0

        #NOTE: If you created the database with the duplicate handling option, there is a chance of inconsistent behavior of anything involving interactions as
        the Pose_ID was not used as an explicit foreign key in db v1.0.0 and v1.1.0.

        Args:
            consent (bool, optional): variable to ensure consent to update database is explicit

        Returns:
            bool: final consent
        """
        raise_not_implemented()

    def db_empty(self):
        """empty database, for example if overwrite

        Returns:
            bool: whether or not db is empty
        """
        raise_not_implemented()

    def db_query(self, query, params: iter) -> iter:
        """Executes provided sql query. Returns iter for results.

        Args:
            query (str): Formated sql query as string
            params (iter): parameters to be used in query (assumes query as appropriate place holders)

        Returns:
            iter: Contains results of query
        """
        raise_not_implemented()

    def db_update(self, query: str, parameters: list[tuple], commit=True) -> iter:
        """
        A db query that also commits if/when specified

        Args:
            query (str): sql formatted query string
            parameters (list[tuple]): assumes appropriate place holders in query
            commit (bool, optional): whether to commit the transaction in open connection. Defaults to True.

        Raises:
            OptionError
            DatabaseInsertionError

        Returns:
            iter: if requesting return value(s)
        """
        raise_not_implemented()

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
            merging_db (str): path to database being merged into current
            backup (bool, optional): whether or not to back up current database before
                merging another database into it. Defaults to True.

        Raises:
            MergeError
        """
        raise_not_implemented()

    def table_length(self, table) -> int:
        """
        Get length of table or bookmark

        Args:
            table (str): name of table or bookmark

        Returns:
            int: number of poses in table/bookmark
        """
        raise_not_implemented()

    def crossref_filter(
        self,
        new_db: str,
        bookmark1_name: str,
        bookmark2_name: str,
        temp_table_suffix: int = 0,
        selection="NOT IN",
        old_db=None,
    ) -> tuple:
        """Selects ligands found or not found in the given bookmark in both current db and new_db.
        Stores as a temporary table, only accessible within the same database connection.

        Args:
            new_db (str): file name for database to attach
            bookmark1_name (str): string for name of first bookmark/temp table to compare
            bookmark2_name (str): string for name of second bookmark to compare
            temp_table_suffix (int, optional): if comparing more than set of bookmarks in one database connection, use this to give different temp table names
            selection (str, optional): "IN" or "NOT IN" indicating if ligand names should or should not be in both databases
            old_db (str, optional): file name for previous database

        Returns:
            tuple: (name of new bookmark (str), number of ligands passing new bookmark (int))
        """
        raise_not_implemented()

    def prune_nonpassing(self, bookmark_name: str):
        """
        Used when creating a new database from filtered data, will remove the data
        that did not pass filtering

        Args:
            bookmark_name (str): bookmark name which has the only poses to save
        """

        self._delete_from_results(bookmark_name)
        self._delete_from_ligands(bookmark_name)
        self._delete_from_interactions_not_in_view(bookmark_name)
        self._clear_bookmarks()

    def get_bookmark_selection(
        self,
        bookmark_name: str,
        selection: Union[list, str],
        group_by: str = None,
        order_results: str = None,
    ) -> str:
        """
        Generates query to gather chosen columns based on passing poses in a bookmark

        Args:
            bookmark_name (str): bookmark name from which to get the passing poses
            selection (Union[list, str]): what columns to have in the output query
            group_by (str, optional): whether or not to group the output by a column
            order_results (str, optional): Whether or not to order by a column

        Raises:
            OptionError: _description_

        Returns:
            str: sql string that describes selection of data from bookmark
        """
        raise_not_implemented()

    def fetch_filters_from_bookmark(self, bookmark_name: str) -> dict:
        """Method that will retrieve filter values used to construct bookmark

        Args:
            bookmark_name (str): bookmark for which to get filters

            Returns:
                dict: containing the filter data
        """
        raise_not_implemented()

    def fetch_filters_and_filterwindow(self, bookmark_name: str) -> dict:
        """Method that will retrieve filter values used to construct bookmark
        and the filter window used as basis

        Args:
            bookmark_name (str): bookmark which was the result of the filtering

            Returns:
                tuple(dict, str): containing the filter data and filter window
        """
        raise_not_implemented()

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
        raise_not_implemented()

    def get_plot_data(
        self,
        bookmark_name: str,
        only_passing: bool = False,
        include_status: bool = False,
        x_axis: str = "docking_score",
        y_axis: str = "leff",
    ):
        """This function gathers two docking results columns (docking score and ligand efficienct) from all data,
        as well as pose_id and ligand name from given bookmark. Can request the data just for poses in the bookmark.

        Args:
            bookmark_name (str): name of bookmark for which to fetch passing data. Returns empty list if bookmark does not exist.
            only_passing (bool): Only return data for passing ligands. Will return empty list for all data.
            include_status (bool): look for status tables and include if requested
            x_axis (str, optional): Defaults to "docking_score".
            y_axis (str, optional): Defaults to "leff".

        Returns:
            tuple: cursors as (<all data cursor>, <passing data cursor>)
        """
        raise_not_implemented()

    def fetch_receptor_object(self) -> tuple:
        """Returns all Receptor objects from database

        Args:
            rec_name (str): Name of receptor to return object for

        Returns:
            tuple: of receptor name and object
        """
        raise_not_implemented()

    def fetch_clustered_similars(self, ligname: str):
        """Given ligname, returns poseids for similar poses/ligands from previous clustering. User prompted at runtime to choose cluster.

        Args:
            ligname (str): ligname for ligand to find similarity with

        Raises:
            ValueError: wrong terminal input
        """
        raise_not_implemented()

    def fetch_rdkit_relevant_pose_properties(self, pose_ids: list) -> iter:
        """
        Gets molecular data that is needed to create rdkit mols for a given list of poses

        Args:
            pose_ids (list): pose ids for which to collect molecular data

        Returns:
            iter: of the following columns Pose_ID, docking_score, leff, ligand_coordinates, flexible_res_coordinates
        """
        raise_not_implemented()

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
        raise_not_implemented()

    def fetch_data_for_passing_results(
        self, bookmark_name: str, outfields: Union[str, list], order_results: str = None
    ) -> iter:
        """Will return SQLite cursor with requested data for outfields for poses that passed filter in bookmark_name

        Returns:
            iter: sqlite cursor of data from passing data

        Raises:
            OptionError
        """
        raise_not_implemented()

    def fetch_flexres_info(self, receptor):
        """fetch flexible residues names and atomname lists

        Returns:
            tuple: (flexible_residues, flexres_atomnames)
        """
        raise_not_implemented()

    def fetch_passing_ligands_rdkit_relevant_info(self, bookmark_name: str) -> iter:
        """fetch information required by vsmanager for writing out molecules

        Returns:
            iter: contains LigName, rdmol,
                atom_index_map, hydrogen_parents
        """
        raise_not_implemented()

    def fetch_ligand_rdkit_relevant_info(self, ligname: str) -> tuple:
        """fetch information required by vsmanager for writing out molecules

        Returns:
            tuple: contains rdmol, atom_index_map, hydrogen_parents
        """
        raise_not_implemented()

    def fetch_lignames_and_poses_for_selection(
        self, selection: str
    ) -> dict[str, list[int]]:
        """
        Creates a dictionary of ligands and the selected poses that appear in a selection,
        such as a bookmark or a status table (where only poses are given).

        Args:
            selection (str): name of bookmark or status table

        Returns:
            dict[str, list[int]]: ligand name is keyword, value is list of poses in given selection
        """
        raise_not_implemented()

    def fetch_pose_interactions(self, Pose_ID) -> iter:
        """
        Fetch all interactions parameters belonging to a Pose_ID

        Args:
            Pose_ID (int): pose id, 1-1 with Results table

        Returns:
            iter: of interaction information for given Pose_ID
        """
        raise_not_implemented()

    def count_receptors_in_db(self):
        """returns number of rows in Receptors table where receptor_object already has blob

        Returns:
            int: number of rows in receptors table
            str: name of receptor if present in table
        """
        raise_not_implemented()

    def to_dataframe(self, requested_data: str, table=True) -> pd.DataFrame:
        """Returns a panda dataframe of table or query given as requested_data

        Args:
            requested_data (str): String containing SQL-formatted query or table name
            table (bool): Flag indicating if requested_data is table name or not

        Returns:
            pd.DataFrame: dataframe of requested data
        """
        raise_not_implemented()

    def get_range_of_e_le(self, table: str) -> tuple:
        """
        Get min and max of e/docking_score and ligand efficiency/le/leff/

        Args:
            table (str): table limit data, e.g., either Results or a bookmark name

        Returns:
            tuple: e_min, e_max, le_min, le_max
        """
        raise_not_implemented()

    def create_status_tables(self):
        """
        Creates status tables if needed
        """
        raise_not_implemented()

    def accept_pose(self, pose_id: int):
        """
        Will add pose_id to accepted, and delete from maybe and rejected if needed

        Args:
            pose_id (int)
        """
        raise_not_implemented()

    def maybe_pose(self, pose_id: int):
        """
        Will add pose_id to maybe, and delete from accepted and rejected if needed

        Args:
            pose_id (int)
        """
        raise_not_implemented()

    def reject_pose(self, pose_id: int):
        """
        Will add pose_id to rejected, and delete from accepted and maybe if needed

        Args:
            pose_id (int)
        """
        raise_not_implemented()

    def fetch_viewable_data_columns_from(
        self, table: str, length: int, last_pose_id: int = 0
    ) -> iter:
        """
        Makes a selection of columns and includes the status of the pose

        Returns:
            iter: iterable/cursor of the columns
        """
        raise_not_implemented()

    def fetch_columns_from_table_as_dicts(
        self, table: str, columns: list, length: int = 500, starting_rowid: int = 0
    ) -> tuple[list[str], list[dict]]:
        """
        Will get requested table data for a table given one or more columns.
        Data will be limited by a certain length, and can be retrieved from a desired
        rowid.

        Args:
            table (str): name of table or bookmark
            columns (list, optional): list of columns to retrieve. Defaults to ["*"].
            length (int, optional): number of rows to collect. Defaults to 500.
            starting_rowid (int, optional): rowid to start with. Defaults to 0.

        Returns:
            tuple[list[str], list[dict]]: list of column names, and list of dicts where each dict is one row,
                                            and column is the key, value is the row-col cell value
        """
        raise_not_implemented()

    def get_query_data_as_dicts(self, query: str) -> tuple[list[str], list[dict]]:
        """
        Will return data requested in an sql formatted query

        Args:
            query (str): sql query formatted to sqlite database

        Returns:
            tuple[list[str], list[dict]]: list of column names, and list of dicts where each dict is one row,
                                            and column is the key, value is the row-col cell value
        """
        raise_not_implemented()

    def is_bookmark(self, table: str) -> bool:
        """
        Returns True if table name is actually a bookmark

        Args:
            table (str): name of table or bookmark to check

        Returns:
            bool: if table name is a bookmark
        """
        raise_not_implemented()

    def pose_row_in_table(self, table: str, pose_id: int) -> Union[None, int]:
        """
        Find the row id of a pose in a given table

        Args:
            table (str)
            pose_id (int)

        Returns:
            Union[None, int]: rowid if any
        """
        raise_not_implemented()

    def tables_in_db(self) -> list:
        """
        Returns a list of all table names in the database

        Returns:
            list: list of table names
        """
        raise_not_implemented()

    def calculate_percentiles(
        self, column: str, num_bins: int, table: str
    ) -> tuple[list[int], list[float]]:
        """
        Will calculate percentiles for a given column and given number of bins to divide the data in.
        Will group the data by ligand_id, so it will be per ligand and not per pose id.

        Args:
            column (str): what column to calculate percentile for. must be numeric
            num_bins (int): how many percentile bins data should be divided into
            table (str): whether the column is in Results or filtered results (i.e., bookmark)

        Raises:
            OptionError: if column given is not numeric and in results

        Returns:
            tuple[list[int],list[float]]: list of percentiles as bins, and list of edge of each bin
        """
        raise_not_implemented()

    def get_previous_docking_mode(self) -> Union[None, str]:
        """
        Checks the docking_mode last used in a database write session

        Returns:
            Union[None, str]: docking_mode if any
        """
        raise_not_implemented()

    def fetch_selected_ligand_poses(self, ligand_name: str, selection: str):
        """
        Gets only the poses of a given ligand that are present in give selection (e.g., a bookmark)

        Args:
            ligand_name (str)
            selection (str): status table or bookmark name

        Returns:
            list[int]: selected poses for ligand
        """
        raise_not_implemented()

    def get_starting_rowid(self, table: str) -> int:
        """
        Starting row id for a table, will be 1 for regular tables, and 1 or non-1 for bookmarks
        (whose rows are inside Filtered_poses)

        Args:
            table (str): table or bookmark name

        Returns:
            int: first row id belonging to that selection
        """
        raise_not_implemented()

    def get_bookmark_selection(
        self,
        bookmark_name: str,
        selection: Union[list, str],
        group_by: str = None,
        order_results: str = None,
    ) -> str:
        """
        Generates query to gather chosen columns based on passing poses in a bookmark

        Args:
            bookmark_name (str): bookmark name from which to get the passing poses
            selection (Union[list, str]): what columns to have in the output query
            group_by (str, optional): whether or not to group the output by a column
            order_results (str, optional): Whether or not to order by a column

        Raises:
            OptionError: _description_

        Returns:
            str: sql string that describes selection of data from bookmark
        """

        if selection and selection != "*":
            outfields_list = self._format_output_fields(selection, "R", "L")
        elif selection == "*":
            raise OptionError(
                "Output fields/columns cannot be 'all'/'*', please select one or more specific columns, or use the default."
            )
        # start formatting write query
        query = self.QueryBuilder()
        # select stuff from results where pose id in filter poses join ligands for extra fields
        query.SELECT(*outfields_list).FROM("Results", "R").JOIN(
            "Ligands", "L", "ligand_id"
        ).WHERE(f"R.pose_id IN ({self._get_bookmark_poses_query(bookmark_name)})")

        if group_by:
            query.GROUP_BY("r.ligand_id")
        if order_results:
            order_by = self._format_orderby(order_results)
            if order_by:
                query.ORDER_BY(order_by)

        return query.build()[0]

    # endregion

    # region private methods

    def _open_storage(self):
        """Create connection to db. Then, check if db needs to be created.

        Raises:
            StorageError
        """
        try:
            # check to see if file exist, and if it does, check that version is matching
            if os.path.isfile(self.db_file) and os.path.getsize(self.db_file) > 0:
                self.conn = self._create_connection()
                compatible, db_version = self.check_ringtaildb_version()
                if not compatible:
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

    def _create_tables(self) -> None:
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
    def _generate_results_row(cls, ligand_dict, pose_rank, run_number) -> list:
        """generate list of lists of ligand values to be
            inserted into duckdb database for a given pose

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
        # # # # # # get pose-specific data

        # check if run is best for a cluster.
        # We are only saving the top pose for each cluster
        ligand_data_list = [
            ligand_dict["receptor"],
            pose_rank + 1,
            int(run_number),
        ]
        # get energy data
        for key in ligand_data_keys:
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
        state_var_keys = ["pose_about", "pose_translations", "pose_quarternions"]
        # add statevars
        for key in state_var_keys:
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
    def _generate_ligand_row(cls, ligand_dict: dict) -> list:
        """writes row to be inserted into ligand table

        Args:
            ligand_dict (dict): Dictionary of ligand data from parser

        Returns:
            list: List of data to be written as row in ligand table. Format:
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

    @classmethod
    def _generate_receptor_row(cls, ligand_dict: dict) -> list:
        """Writes row to be inserted into receptor table

        Args:
            ligand_dict (dict): Dictionary of ligand data from parser

        Returns:
            list: receptor row columns
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
    def _generate_interaction_tuples(cls, interaction_dictionaries: list) -> list:
        """takes dictionary of file results, formats as
        list of tuples for interactions

        Args:
            interaction_dictionaries (list): List of pose interaction
                dictionaries from parser

        Returns:
            list: List of tuples of interaction data
        """
        interaction_keywords = [
            "type",
            "chain",
            "residue",
            "resid",
            "recname",
            "recid",
        ]
        interactions = set()
        for pose_interactions in interaction_dictionaries:
            count = pose_interactions["count"][0]
            for i in range(int(count)):
                interactions.add(
                    tuple(pose_interactions[kw][i] for kw in interaction_keywords)
                )

        return list(interactions)

    def _check_unique_results_row(self, result_data: list) -> int:
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

        query = self.QueryBuilder()
        query.SELECT("Pose_ID").FROM("Results").WHERE(
            """ligand_id=?
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
                    AND dihedrals=?""",
            *unique_data,
        )
        try:
            data = self.db_query(*query.build())
        except Exception as e:
            raise e
        row = data.fetchone()
        if row is None:
            Pose_ID = -1
            logger.debug("Duplicate row not found.")
        else:
            Pose_ID = row[0]
            logger.debug(f"Duplicate row found for Pose_ID {Pose_ID}")

    def _create_results_table(self):
        pass

    def _create_ligands_table(self):
        pass

    def _create_receptors_table(self):
        pass

    def _create_interaction_index_table(self):
        pass

    def _create_interaction_table(self):
        pass

    def _create_db_properties_table(self):
        pass

    def _create_filtering_tables(self):
        pass

    def _set_ringtail_db_schema_version(self):
        pass

    def _insert_ligands(self, ligand_array: list) -> list:
        raise_not_implemented()

    def _create_indices(self):
        raise_not_implemented()

    def _insert_results(self, results_array, options):
        raise_not_implemented()

    def _insert_receptors(self, receptor_array):
        raise_not_implemented()

    def _insert_and_format_interactions(
        self, pose_ids: list, docking_data: dict
    ) -> list:
        """
        This method will evaluate the docking data, and determine whether or not there are
        interactions to be processed and written.


        Args:
            pose_ids (list[int]): pose ids being processed
            docking_data (dict): all docking data

        Returns:
            interactions_list (list): List of tuples for interactions in form
                ("type", "chain", "residue", "resid", "recname", "recid")
        """
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

    def _insert_interaction_index_row(self, interaction_tuple) -> int:
        """
        Writes unique interactions and returns the interaction_id of the given interaction

        Args:
            interaction_tuple (tuple): (rec_chain, rec_resname, rec_resid, rec_atom, rec_atomid)

        Returns:
            int: interaction index

        Raises:
            DatabaseInsertionError
        """
        raise_not_implemented()

    def _insert_interaction_rows(
        self, interaction_rows, duplicates, duplicate_handling
    ):
        raise_not_implemented()

    def _generate_result_filtering_query(
        self, filters_dict, bookmark_name, filter_bookmark
    ):
        raise_not_implemented()

    def _populate_filter_tables(self, name, query, filters):
        raise_not_implemented()

    def _insert_receptor_blob(self, receptor, rec_name):
        raise_not_implemented()

    def _clear_bookmarks(self):
        """
        Clears all filters and filtered poses
        """
        query = self.QueryBuilder()
        self.db_query(*(query.DROP_IF_EXISTS("Filters").build()))
        query = self.QueryBuilder()
        self.db_query(*(query.DROP_IF_EXISTS("Filtered_poses").build()))
        self._create_filtering_tables()

    def _process_filters_for_query(self, filters_dict: dict):
        score_maxmin_to_duckdb_call = {
            "eworst": "docking_score <= {value}",
            "ebest": "docking_score >= {value}",
            "leworst": "leff <= {value}",
            "lebest": "leff >= {value}",
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
                        score_maxmin_to_duckdb_call[filter_key].format(
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
        from rdkit.Chem import Descriptors

        maxatoms = 0
        position = False
        substruct_mols = []
        # handle single value filters
        if "ligand_operator" in ligand_filters:
            logical_operator = ligand_filters["ligand_operator"]
        if "ligand_max_atoms" in ligand_filters:
            maxatoms = ligand_filters["ligand_max_atoms"]
        if "ligand_min_molweight" in ligand_filters:
            min_mw = ligand_filters["ligand_min_molweight"]
        else:
            min_mw = None
        if "ligand_max_molweight" in ligand_filters:
            max_mw = ligand_filters["ligand_max_molweight"]
        else:
            max_mw = None

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
            cursor = self.db_query(query)

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
                        # continue for ligandrow in _stream_query, ligand did not pass num atoms filter
                        continue
                # check if mol weight filters are present
                if min_mw or max_mw:
                    if not (
                        (min_mw is None or Descriptors.MolWt(ligand_mol) >= min_mw)
                        and (max_mw is None or Descriptors.MolWt(ligand_mol) <= max_mw)
                    ):
                        # continue for ligandrow in _stream_query, ligand did not pass molweight filter
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

    def _get_bookmark_poses_query(self, bookmark_name: str) -> str:
        """
        Creates a query that retrieves all poses from a bookmark, that can be used in other queries

        Args:
            bookmark_name (str): bookmark for which to create the query

        Returns:
            str: query representing the poses in a bookmark
        """
        return self.QueryBuilder.bookmark_query(bookmark_name)

    # endregion


class StorageManagerSQLite(StorageManager):
    """SQLite-specific StorageManager subclass

    Attributes:
        db_file (str): database name
        conn (SQLite.conn): Connection to database
    """

    # "db_schema_ver":list("compatible code versions")
    _db_schema_code_compatibility = {
        "1.0.0": ["1.0.0"],
        "1.1.0": ["1.1.0"],
        "2.0.0": ["2.0.0", "2.1.0", "2.1.1", "2.1.2"],
        "3.0.0": ["3.0.0"],
    }
    QueryBuilder = QueryBuilderSQLite

    def __init__(
        self,
        db_file: str = None,
    ):
        self.db_file = db_file
        super().__init__()

    # region Methods for creating and inserting into tables the database

    def _create_ligands_table(self, name="Ligands") -> None:
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
        ligand_table = f"""CREATE TABLE IF NOT EXISTS {name} (
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

    def _insert_ligands(self, ligand_array: list) -> list:
        """Takes array of ligand rows, inserts into Ligands table.

        Args:
            ligand_array (list[list]): list of lists containing formatted ligand rows

        Returns:
            list: of ligand IDs just inserted

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

    def _create_results_table(self, name="Results"):
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

        sql_results_table = f"""CREATE TABLE IF NOT EXISTS {name} (
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

        except sqlite3.OperationalError as e:
            raise DatabaseQueryError(
                "Error while looking for unique result row."
            ) from e

    def _insert_results(self, results_array, options: dict) -> tuple:
        """Takes list of database rows to insert, adds data to results table. Will handle duplicates if specified

        Args:
            results_array (list): list of arrays containing formatted result rows
            options (dict): includes options on how to handle duplicates if there are any

        Returns:
            list[int]: returns the pose ids for the ligand written to results, these are used to ensure internal consistency when writing to the interaction table
            list[int]: found duplicates (knowledge may be needed in other methods)

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

    def _insert_receptor_blob(self, receptor: bytes, rec_name: str):
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

    def _insert_interaction_index_row(self, interaction_tuple) -> int:
        """
        Writes unique interactions and returns the interaction_id of the given interaction

        Args:
            interaction_tuple (tuple): (rec_chain, rec_resname, rec_resid, rec_atom, rec_atomid)

        Returns:
            int: interaction index

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

    def _create_filtering_tables(self):
        """
        Creates a Filter table which includes filter_id (PK), name (bookmark_name), sqlite formatted query,
        and dictionary of filters used, as well as Filtered_poses, which uses filter_id as FK,
        and lists all poses passing that filter_id

        Raises:
            DatabaseTableCreationError
        """
        # Create filters table keeping track of filter id etc
        filters_sql = """CREATE TABLE IF NOT EXISTS Filters (
        filter_id           INTEGER PRIMARY KEY AUTOINCREMENT,
        name                VARCHAR[],
        query               VARCHAR[],
        filters             VARCHAR[],
        filter_window       VARCHAR[]);"""

        filter_pose_sql = """CREATE TABLE IF NOT EXISTS Filtered_poses (
        filter_id           INTEGER,
        pose_id             INTEGER,
        FOREIGN KEY(filter_id) REFERENCES Filters(filter_id),
        FOREIGN KEY(pose_id) REFERENCES Results(pose_id));"""

        try:
            cur = self.conn.cursor()
            cur.execute(filters_sql)
            cur.execute(filter_pose_sql)
            self.conn.commit()
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
    ) -> str:
        """Insert cluster data into ligand cluster table

        Args:
            clusters (list[list]): list of clusters
            poseid_list (list): representative poses for each sluter
            cluster_type (str): how clustering was performed
            cluster_cutoff (str): distance to representative pose
            bookmark_name (str): bookmark name which is clustered over

        Returns:
            str: name of cluster bookmark
        """
        # TODO
        cur = self.conn.cursor()
        cur.execute(
            "CREATE TABLE IF NOT EXISTS Ligand_clusters (pose_id  INT[] UNIQUE)"
        )
        ligand_cluster_columns = self._fetch_ligand_cluster_columns()
        column_name = (
            f"{bookmark_name}_{cluster_type}_{cluster_cutoff.replace('.', 'p')}"
        )
        try:
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

            return column_name
        except sqlite3.OperationalError as e:
            raise StorageError("Error occurred while inserting cluster data") from e

    def _create_indices(self):
        """Create index for specified tables and columns. 'ak' stands for 'alternate key' and is prepended to index name to avoid naming conflicts

        Raises:
            StorageError
        """
        try:
            cur = self.conn.cursor()
            logger.debug("Creating columns indices...")
            cur.execute(
                "CREATE INDEX IF NOT EXISTS ak_results ON Results(docking_score, leff)"
            )
            cur.execute(
                "CREATE INDEX IF NOT EXISTS ak_resultids ON Results(Pose_id, ligand_id)"
            )
            cur.execute(
                "CREATE INDEX IF NOT EXISTS ak_interactions ON Interactions(Pose_id, interaction_id)"
            )
            cur.execute("CREATE INDEX IF NOT EXISTS ak_ligands ON Ligands(ligand_id)")
            self.conn.commit()
            cur.close()
            logger.info(
                "Indicies were created for specified Results, Ligands, and Interaction_indices columns."
            )
        except sqlite3.OperationalError as e:
            raise StorageError("Error occurred while indexing") from e

    # endregion

    # region merge databases

    def merge_databases(self, merging_db: str, backup: bool = True):
        """
        Method that merges two databases, ensuring integrity of primary and foreign keys.
        The merging will create a new table if needed, that keeps track of the primary key
        in the original and the merged database on a per-table basis. Another table will also
        keep track of how many databases has been merged into the primary database.
        Each merge session is given a merge_id.
        The merging will ensure the two databases are -compatible based on the receptor only-.
        PLEASE NOTE: If two databases have been docked with dlg and vina respectively,
            these will be allowed to merge.

        Args:
            merging_db (str): path to database being merged into current
            backup (bool, optional): whether or not to back up current database before
                merging another database into it. Defaults to True.

        Raises:
            MergeError
        """
        # back up main database
        if backup:
            self.clone()
        # attach incoming database and check compatibility
        merging_db_alias = self._attach_db(merging_db, "merging")
        if not self._check_if_db_compatible_for_merge(merging_db_alias):
            raise StorageError(
                "Trying to merge two databases of incompatible or too old versions, cannot proceed."
            )
        # create new tables to keep track of merger
        self._create_merge_tables()
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
        self._delete_filter_data()
        self._delete_table("Ligand_clusters")
        # for attached db
        self._delete_filter_data(merging_db_alias)
        self._delete_table("Ligand_clusters", merging_db_alias)

        # merge tables
        try:
            self._merge_db_properties_table(merge_id)
            logger.info("The 'db_properties' table has been merged.")

            self._merge_ligands_and_results_tables(merge_id)
            logger.info("The 'Ligands' and 'Results' tables have been merged.")

            self._merge_interaction_tables(merge_id)
            logger.info(
                "The 'Interaction_indices' and 'Interactions' tables have been merged."
            )
        except Exception as e:
            raise MergeError(f"Error during database merging: {e}") from e
        else:
            logger.info(
                f"The database {merging_db} has been successfully merged into {self.db_file}.\n Rebuilding indices."
            )
            self._cleanup_storage(merging_db_alias, vacuum=True, reindex=True)
            logger.info("The final database has neem cleaned up, and indices rebuilt.")
        finally:
            cur.close()

    def _create_merge_tables(self):
        """
        Creates tables necessary when merging two or more databases

        Raises:
            StorageError
        """
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

    def _merge_interaction_tables(self, merge_id: int):
        """
        Merges the interaction tables. Interaction definitions are unique and independent of the Results table, so we only
        insert those that are new with updated PK, and assign existing interaction_ids to those already described in primary db

        Args:
            merge_id (int): merge session id

        Raises:
            Exception
        """
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

    def _merge_ligands_and_results_tables(self, merge_id: int):
        """
        Merges first the Ligands table, then the Results table, maintaining ligand_id and pose_id as primary keys,
        where their relationship to the original PK is kept track of in the mering datble. Duplicate ligands will maintain
        the ligand_id from the main database.

        Args:
            merge_id (int): merge session id

        Raises:
            StorageError
        """
        # convert ligand_ids and log
        convert_ligand_ids_sql = """INSERT INTO PK_conversions (
        merge_id,
        table_name,
        original_PK,
        merged_PK) SELECT 
        ?,
        "Ligands", 
        ligand_id,
            CASE 
                WHEN EXISTS (
                    SELECT 1
                    FROM Ligands
                    WHERE 
                        merging.Ligands.LigName = Ligands.LigName
                        AND merging.Ligands.ligand_smile = Ligands.ligand_smile
                        AND merging.Ligands.rdmol = Ligands.rdmol
                        AND merging.Ligands.atom_index_map = Ligands.atom_index_map
                        AND merging.Ligands.hydrogen_parents = Ligands.hydrogen_parents
                        AND merging.Ligands.input_model = Ligands.input_model
                    ) 
                THEN (
                    SELECT main.Ligands.ligand_id
                    FROM main.Ligands
                    WHERE 
                        merging.Ligands.LigName = Ligands.LigName
                        AND merging.Ligands.ligand_smile = Ligands.ligand_smile
                        AND merging.Ligands.rdmol = Ligands.rdmol
                        AND merging.Ligands.atom_index_map = Ligands.atom_index_map
                        AND merging.Ligands.hydrogen_parents = Ligands.hydrogen_parents
                        AND merging.Ligands.input_model = Ligands.input_model
                    )
                ELSE merging.Ligands.ligand_id + (SELECT MAX(ligand_id) FROM Ligands)
            END AS new_ligand_id
        FROM merging.Ligands;"""

        # then inserting only those that aren't already in the table
        insert_new_ligands = """INSERT INTO Ligands (
        ligand_id,
        LigName,
        ligand_smile,
        rdmol,
        atom_index_map,
        hydrogen_parents,
        input_model)
        SELECT 
            (SELECT merged_PK FROM PK_conversions WHERE original_PK = ligand_id and merge_id = ? AND table_name = 'Ligands') new_id,
            LigName,
            ligand_smile,
            rdmol,
            atom_index_map,
            hydrogen_parents,
            input_model
        FROM merging.Ligands WHERE new_id > (SELECT MAX(ligand_id) FROM Ligands);
        """

        # convert pose_id
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
            pose.merged_PK as pose_id,
            ligand.merged_PK as ligand_id,
            mr.receptor,
            mr.pose_rank,
            mr.run_number,
            mr.docking_score,
            mr.leff,
            mr.deltas,
            mr.cluster_rmsd,
            mr.cluster_size,
            mr.reference_rmsd,
            mr.energies_inter,
            mr.energies_vdw,
            mr.energies_electro,
            mr.energies_flexLig,
            mr.energies_flexLR,
            mr.energies_intra,
            mr.energies_torsional,
            mr.unbound_energy,
            mr.nr_interactions,
            mr.num_hb,
            mr.about_x,
            mr.about_y,
            mr.about_z,
            mr.trans_x,
            mr.trans_y,
            mr.trans_z,
            mr.axisangle_x,
            mr.axisangle_y,
            mr.axisangle_z,
            mr.axisangle_w,
            mr.dihedrals,
            mr.ligand_coordinates,
            mr.flexible_res_coordinates
        FROM merging.Results mr
        LEFT JOIN (
            SELECT original_PK, merged_PK 
            FROM PK_conversions 
            WHERE table_name = 'Results' 
            AND merge_id = ?
            ) pose 
        ON pose.original_PK = mr.Pose_ID
        LEFT JOIN (
            SELECT original_PK, merged_PK 
            FROM PK_conversions 
            WHERE table_name = 'Ligands' 
            AND merge_id = ?
            ) ligand
        ON (mr.ligand_id = ligand.original_PK);"""

        try:
            cur = self.conn.cursor()
            # insert and get new ligand ids
            cur.execute(
                convert_ligand_ids_sql,
                (merge_id,),
            )
            cur.execute(insert_new_ligands, (merge_id,))
            # get new pose_ids
            cur.execute(
                convert_poseid_sql,
                (merge_id,),
            )
            # insert results with new pose_ids and ligand_ids
            cur.execute(
                insert_Results,
                (
                    merge_id,
                    merge_id,
                ),
            )

            self.conn.commit()
        except Exception as e:
            raise StorageError(
                f"Error encountered while merging Ligands and Results tables: \n{str(e)}"
            ) from e

    def _merge_db_properties_table(self, merge_id: int):
        """
        Merges database properties table, but importantly will not check for property compatibility

        Args:
            merge_id (int): merge session id

        Raises:
            StorageError
        """
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

    # endregion

    # region Methods for dealing with bookmarks and filtering

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
            if self.is_bookmark(filter_bookmark):
                filtering_window = f"""(SELECT * FROM Results WHERE Pose_id IN ({self.QueryBuilder.bookmark_query(filter_bookmark)}))"""
            elif self._is_statustable(filter_bookmark):
                filtering_window = f"""(SELECT * FROM Results WHERE Pose_id IN (SELECT Pose_ID FROM {filter_bookmark}))"""

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
                lig_names = lig_filters.pop("ligand_name")
                ligname_query = " OR ".join(
                    [f"LigName LIKE '%{ligname}%' " for ligname in lig_names if ligname]
                )
                ligname_query = "SELECT ligand_id FROM Ligands WHERE " + ligname_query
            # rdkit queries need to be handled in memory separate from the main query
            if lig_filters:
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

        filter_query = f"SELECT R.pose_id FROM {filtering_window} R {partial_filter_query} ORDER BY R.pose_id"

        return filter_query

    def _populate_filter_tables(
        self, name, query: str, filters={}, filtering_bookmark: str = ""
    ) -> bool:
        """
        Will run a filter query and determine if there are passing poses, in which case all relevant
        data is written to the database

        Args:
            name (str): name of new bookmark
            query (str): query that defines what poses to insert
            filters (dict, optional): filters or restrictions used
            filtering_bookmark (str, optional): If filters were performed across an existing obokmark. Defaults to None.

        Raises:
            StorageError
            OptionError

        Returns:
            bool: whether or not there are poses passing the filter
        """

        # fetch filtered poses
        passing_poses_tuples = self.db_query(query).fetchall()
        passing_poses = [row[0] for row in passing_poses_tuples]
        if passing_poses:
            # make sure bookmark name is not a table name
            if name in self.tables_in_db():
                raise OptionError(
                    f"Bookmark name {name} is the same as an existing table in the database, and cannot be used."
                )
            # check if bookmark exists
            if self.is_bookmark(name):
                logger.warning(
                    f"The bookmark {name} already exists, and will be overwritten by the current filter."
                )
                self.delete_bookmark(name)
            filter_sql = """INSERT INTO Filters (name,query,filters,filter_window) VALUES (?,?,?,?) RETURNING filter_id;"""
            try:
                filter_id = self.db_query(
                    filter_sql,
                    (name.lower(), query, json.dumps(filters), filtering_bookmark),
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

    # endregion

    # region crossreferencing filtered databases

    def crossref_filter(
        self,
        new_db: str,
        bookmark1_name: str,
        bookmark2_name: str,
        temp_table_suffix: int = 0,
        selection="NOT IN",
        old_db=None,
    ) -> tuple:
        """Selects ligands found or not found in the given bookmark in both current db and new_db.
        Stores as a temporary table, only accessible within the same database connection.

        Args:
            new_db (str): file name for database to attach
            bookmark1_name (str): string for name of first bookmark/temp table to compare
            bookmark2_name (str): string for name of second bookmark to compare
            temp_table_suffix (int, optional): if comparing more than set of bookmarks in one database connection, use this to give different temp table names
            selection (str, optional): "IN" or "NOT IN" indicating if ligand names should or should not be in both databases
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

        num_passing = self._count_ligands_in_temptable(temp_name)
        print("\n\n Number passing the cross referenced filters: ", num_passing)

        return temp_name, num_passing

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
        query = self.QueryBuilder()
        subq = self.QueryBuilder()
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

        self._populate_filter_tables(bookmark_name, query_string, filters)

    def _delete_from_ligands(self, bookmark_name: str):
        """Remove rows from ligands table if they did not pass filtering

        Raises:
            StorageError
        """
        passing_poses_query = self._get_bookmark_poses_query(bookmark_name)

        try:
            self.db_update(
                f"DELETE FROM Ligands WHERE ligand_id NOT IN (SELECT ligand_id from Results WHERE Pose_ID IN ({passing_poses_query}))",
                (),
            )
        except sqlite3.OperationalError as e:
            raise StorageError(
                f"Error occured while pruning Ligands not in {bookmark_name}"
            ) from e

    def _delete_from_results(self, bookmark_name: str):
        """Remove rows from results table if they did not pass filtering

        Raises:
            StorageError
        """
        passing_poses_query = self._get_bookmark_poses_query(bookmark_name)
        try:
            self.db_update(
                f"DELETE FROM Results WHERE Pose_ID NOT IN ({passing_poses_query})", ()
            )
        except sqlite3.OperationalError as e:
            raise StorageError(
                f"Error occured while pruning Results not in {bookmark_name}"
            ) from e

    def _delete_from_interactions_not_in_view(self, bookmark_name: str):
        """Remove rows from interactions table if they were not used for poses that passed filtering.

        Args:
            bookmark_name (str): defines which poses are passing

        Raises:
            StorageError: Description
        """
        passing_poses_query = self._get_bookmark_poses_query(bookmark_name)
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

    def _count_ligands_in_temptable(self, temp_name: str) -> int:
        """
        Counts ligands represented in the temporary table

        Args:
            temp_name (str): name of temporary table

        Returns:
            int: number of poses in temporary table
        """
        counting = self.QueryBuilder()
        count_pool = self.QueryBuilder()
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
        return tuple(self.db_query(counting.build()[0]).fetchone())[0]

    def _create_crossref_temp_table(self, table_name: str):
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
        query = self.QueryBuilder()
        subq = self.QueryBuilder()
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

    # region fetching specific columns

    def _format_orderby(self, column_name: str) -> str:
        """
        Ensures chosen order by column is a valid choice

        Args:
            column_name (str): column to order by

        Returns:
            str: column to order by with appropriate alias
        """
        columns, aliased_columns = self._get_possible_output_columns()
        if column_name.lower() in columns:
            index = columns.index(column_name.lower())
            order_by = aliased_columns[index].format(
                Ligands_alias="L", Results_alias="R"
            )
            return order_by
        else:
            return None

    def _format_output_fields(
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
        possible_columns, table_formatted_columns = self._get_possible_output_columns()

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

    def _get_possible_output_columns(self, tables=["Results", "Ligands"]):
        """
        Gets all column names from given tables

        Args:
            tables (list, optional): Defaults to ["Results", "Ligands"].

        Returns:
            columns (list[str]): list of column names for all listed tables
            columns_with_tablename (list[str.format]): needs formatted with table_alias (one per table) for use
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

    def _get_numeric_columns(self, table_name: str) -> list:
        """
        Method to get the names of all numeric columns in a table, for example for
        allowable sorting options

        Args:
            table_name (str): table name to evaluate

        Returns:
            list: column names that has a numeric type
        """
        return [
            table[0]
            for table in self.db_query(
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
        ]

    def get_range_of_e_le(self, table: str) -> tuple:
        """
        Get min and max of e/docking_score and ligand efficiency/le/leff/

        Args:
            table (str): table limit data, e.g., either Results or a bookmark name

        Returns:
            tuple: e_min, e_max, le_min, le_max
        """

        query = self.QueryBuilder()
        query.SELECT(
            "MIN(R.docking_score)",
            "MAX(R.docking_score)",
            "MIN(R.leff)",
            "MAX(R.leff)",
        ).FROM("Results", "R")
        if self.is_bookmark(table):
            query.JOIN("Filtered_poses", "fp", "pose_id").JOIN(
                "Filters", "f", "filter_id", to="Filtered_poses"
            ).WHERE("f.name = ?", table)
        elif self._is_statustable(table):
            query.JOIN(table, "T", "pose_id")

        return self.db_query(*query.build()).fetchall()[0]

    def fetch_receptor_object(self) -> Union[None, tuple]:
        """Returns all Receptor objects from database

        Returns:
            tuple: of receptor name and object
        """

        cursor = self.db_query(
            "SELECT RecName, receptor_object FROM Receptors"
        ).fetchone()
        if cursor:
            return tuple(cursor)
        else:
            return None

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

    def fetch_data_for_passing_results(
        self, bookmark_name: str, outfields: Union[str, list], order_results: str = None
    ) -> iter:
        """Will return SQLite cursor with requested data for outfields for poses that passed filter in bookmark_name

        Returns:
            iter: sqlite cursor of data from passing data

        Raises:
            OptionError
        """
        outfields_list = self._format_output_fields(
            outfields, ligands_alias="L", results_alias="R"
        )

        bookmark_selection = self._get_bookmark_poses_query(bookmark_name)

        query = self.QueryBuilder()
        query.SELECT(*outfields_list).FROM("Results", "R").WHERE(
            f"R.pose_id IN ({bookmark_selection})"
        ).JOIN("ligands", "L", "ligand_id", "results").GROUP_BY("R.ligand_id")
        if order_results:
            order_by = self._format_orderby(order_results)
            if order_by:
                query.ORDER_BY(order_by)

        return self.db_query(query.build()[0])

    def fetch_filters_from_bookmark(self, bookmark_name: str) -> dict:
        """Method that will retrieve filter values used to construct bookmark

        Args:
            bookmark_name (str): bookmark for which to get filters

            Returns:
                dict: containing the filter data
        """
        query = self.QueryBuilder()
        query.SELECT("filters").FROM("Filters").WHERE("name = ?", bookmark_name)

        filters = self.db_query(*query.build()).fetchone()
        if not filters:
            return {}

        return json.loads(filters[0])

    def fetch_filters_and_filterwindow(self, bookmark_name: str) -> tuple[dict, str]:
        """Method that will retrieve filter values used to construct bookmark
        and the filter window used as basis

        Args:
            bookmark_name (str): bookmark which was the result of the filtering

            Returns:
                tuple(dict, str): containing the filter data and filter window
        """
        if not self.is_bookmark(bookmark_name):
            return {}, ""

        query = self.QueryBuilder()
        query.SELECT("filters", "filter_window").FROM("Filters").WHERE(
            "name = ?", bookmark_name
        )
        filters, filter_window = self.db_query(*query.build()).fetchone()

        return json.loads(filters), filter_window

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
        query = self.QueryBuilder()
        query.SELECT("ligname", "rdmol", "atom_index_map", "hydrogen_parents").FROM(
            "Ligands", "L"
        ).WHERE(
            f"""L.ligand_id IN (SELECT DISTINCT ligand_id FROM ({self.get_bookmark_selection(bookmark_name, "ligand_id")}))"""
        )
        return self.db_query(query.build()[0])

    def fetch_ligand_rdkit_relevant_info(self, ligname: str) -> tuple:
        """fetch information required by vsmanager for writing out molecules

        Returns:
            tuple: contains rdmol, atom_index_map, hydrogen_parents
        """
        query = self.QueryBuilder()
        query.SELECT("rdmol", "atom_index_map", "hydrogen_parents").FROM(
            "Ligands"
        ).WHERE(f"ligname = ?", ligname)
        return self.db_query(*query.build()).fetchone()

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

    def _fetch_ligand_cluster_columns(self) -> list:
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

    def _fetch_results_column_names(self) -> list:
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
        query_ligand_cluster = cur.execute(
            f"SELECT {cluster_col_choice} FROM Ligand_clusters WHERE pose_id IN (SELECT Pose_ID FROM Results WHERE ligand_id = (SELECT ligand_id FROM Ligands WHERE LigName =  '{ligname}'))"
        ).fetchone()
        if query_ligand_cluster is None:
            raise DatabaseQueryError(
                f"Requested ligand name {ligname} not found in cluster {cluster_col_choice}!"
            )
        query_ligand_cluster = query_ligand_cluster[0]  # extract from tuple

        # TODO this is not ideal or optimized, but it does what it needs to do.
        sql_query = f"""
        SELECT LigName FROM Ligands WHERE ligand_id IN (SELECT ligand_id FROM Results WHERE Pose_ID IN 
        (SELECT pose_id FROM Ligand_clusters WHERE {cluster_col_choice}={query_ligand_cluster}))
        GROUP BY ligand_id"""

        bookmark_name = f"similar_{ligname}_{cluster_col_choice}"

        return self.db_query(sql_query).fetchall(), bookmark_name, cluster_col_choice

    def fetch_rdkit_relevant_pose_properties(self, pose_ids: list) -> iter:
        """
        Gets molecular data that is needed to create rdkit mols for a given list of poses

        Args:
            pose_ids (list): pose ids for which to collect molecular data

        Returns:
            iter: of the following columns pose_id, docking_score, leff, ligand_coordinates, flexible_res_coordinates
        """
        placeholders = ",".join(["?"] * len(pose_ids))
        query = f"""
        SELECT pose_id, docking_score, leff, ligand_coordinates, flexible_res_coordinates 
        FROM Results WHERE Pose_ID IN ({placeholders})
        """
        return self.db_query(query, pose_ids).fetchall()

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

    def calculate_percentiles(
        self, column: str, num_bins: int, table: str
    ) -> tuple[list[int], list[float]]:
        """
        Will calculate percentiles for a given column and given number of bins to divide the data in.
        Will group the data by ligand_id, so it will be per ligand and not per pose id.

        Args:
            column (str): what column to calculate percentile for. must be numeric
            num_bins (int): how many percentile bins data should be divided into
            table (str): whether the column is in Results or filtered results (i.e., bookmark)

        Raises:
            OptionError: if column given is not numeric and in results

        Returns:
            tuple[list[int],list[float]]: list of percentiles as bins, and list of edge of each bin
        """

        if not column in self._get_numeric_columns("Results"):
            raise OptionError(
                f"Requested column {column} in not numeric, percentiles cannot be calcualted."
            )
        query = self.QueryBuilder()
        query.SELECT(f"{column}").FROM("Results")
        if self.is_bookmark(table):
            query.IN_BOOKMARK(table)
        elif self._is_statustable(table):
            query.JOIN(table, "T", "pose_id")
        query.GROUP_BY("ligand_id")
        values = [val[0] for val in self.db_query(query.build()[0]).fetchall()]

        bins = np.linspace(0, 100, num_bins + 1)
        bin_edges = np.percentile(values, bins)
        return bins, bin_edges

    def get_plot_data(
        self,
        bookmark_name: str,
        only_passing: bool = False,
        include_status: bool = False,
        x_axis: str = "docking_score",
        y_axis: str = "leff",
        limit: int = None,
    ):
        """This function gathers two docking results columns (docking score and ligand efficienct) from all data,
        as well as pose_id and ligand name from given bookmark. Can request the data just for poses in the bookmark.

        Args:
            bookmark_name (str): name of bookmark for which to fetch passing data. Returns empty list if bookmark does not exist.
            only_passing (bool): Only return data for passing ligands. Will return empty list for all data.
            include_status (bool): look for status tables and include if requested
            x_axis (str, optional): Defaults to "docking_score".
            y_axis (str, optional): Defaults to "leff".

        Returns:
            tuple: cursors as (<all data cursor>, <passing data cursor>)
        """
        all_data_query = self.QueryBuilder()
        all_data_query.SELECT("docking_score", "leff").FROM("Results")
        bookmark_query = self.QueryBuilder()
        bookmark_query.SELECT(
            "R." + x_axis, "R." + y_axis, "R." + "Pose_ID", "L." + "LigName"
        )
        if limit:
            bookmark_query.LIMIT(limit)

        if self.is_bookmark(bookmark_name):
            if include_status:
                bookmark_query.SELECT_STATUS()
            bookmark_query.FROM("Results", "R").IN_BOOKMARK(bookmark_name).JOIN(
                "Ligands", "L", "ligand_id"
            )

            if only_passing:
                all_data = []
            else:
                all_data = self.db_query(all_data_query.build()[0]).fetchall()
            passing_data = self.db_query(bookmark_query.build()[0]).fetchall()

        elif self._is_table(bookmark_name) and bookmark_name.lower() != "results":
            # will assume it is a status table
            if include_status:
                bookmark_query.SELECT(f"""'{bookmark_name.lower()}' as status""")
            bookmark_query.FROM("Results", "R").JOIN(
                bookmark_name, "T", "pose_id"
            ).JOIN("Ligands", "L", "ligand_id")

            if only_passing:
                all_data = []
            else:
                all_data = self.db_query(all_data_query.build()[0]).fetchall()
            passing_data = self.db_query(bookmark_query.build()[0]).fetchall()

        else:
            all_data = self.db_query(all_data_query.build()[0]).fetchall()

            if include_status:
                bookmark_query.SELECT_STATUS()
            bookmark_query.FROM("Results", "R").JOIN("Ligands", "L", "ligand_id")

            passing_data = self.db_query(bookmark_query.build()[0]).fetchall()

        return all_data, passing_data

    # endregion

    # region general database operations

    def to_dataframe(self, requested_data: str, table=True) -> pd.DataFrame:
        """Returns a panda dataframe of table or query given as requested_data

        Args:
            requested_data (str): String containing SQL-formatted query or table name
            table (bool): Flag indicating if requested_data is table name or not

        Returns:
            pd.DataFrame: dataframe of requested data
        """
        if table:
            query = self.QueryBuilder()
            if requested_data in self.get_all_bookmark_names():
                query.SELECT("*").FROM("Results").IN_BOOKMARK(requested_data)
            else:
                # assume it is a table
                query.SELECT("*").FROM("Results")
            return pd.read_sql_query(query.build()[0], self.conn)
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

    def overwrite_storage(self):
        """
        Will drop all tables in the database.
        """
        if not self.db_empty():
            self._drop_existing_tables()
            logger.info("Tables in existing database were dropped.")

    def get_previous_docking_mode(self) -> Union[None, str]:
        """
        Checks the docking_mode last used in a database write session

        Returns:
            Union[None, str]: docking_mode if any
        """
        if self.db_empty():
            return None
        docking_mode = self.conn.execute(
            "SELECT docking_mode FROM DB_properties ORDER BY DB_write_session DESC LIMIT 1"
        ).fetchone()
        return docking_mode[0].lower() if docking_mode else None

    def check_storage_ready(
        self, run_mode: str, docking_mode: str, store_all_poses: bool, max_poses: int
    ):
        """Check that storage is ready

        Raises:
            OptionError: if database options are not compatible
        """
        if self.db_empty():
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

    def _set_ringtail_db_schema_version(self, db_version: str = "3.0.0"):
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

    def check_ringtaildb_version(self) -> tuple[bool, str]:
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
        return is_compatible, db_schema_ver

    def _check_if_db_compatible_for_merge(self, merging_db_alias: str) -> bool:
        """
        Method that checks if the database merging into main is compatible with main,
        and checks if both databases are of appropriately high version where merge has
        been implemented

        Args:
            merging_db_alias (str): alias for the database being merged into main

        Returns:
            bool: if the two databases are compatible
        """
        main_version = self.db_query("PRAGMA main.user_version").fetchone()[0]
        merging_version = self.db_query(
            f"PRAGMA {merging_db_alias}.user_version"
        ).fetchone()[0]

        if main_version != merging_version:
            return False
        if main_version < 300:
            if main_version == 200:
                logger.error(
                    "The database is enabled for merging, but only with an older version of the code. Please contact code maintainer."
                )
            return False

        return True

    def tables_in_db(self) -> list:
        """
        Returns a list of all table names in the database

        Returns:
            list: list of table names
        """
        return [
            name[0].lower()
            for name in self.db_query(
                "SELECT name FROM sqlite_master WHERE type='table';"
            ).fetchall()
        ]

    def update_database_version(self, new_version, consent=False):
        """method that updates sqlite database schema 1.0.0 through 3.0.0.
        The way it currently works, it has to upgrade via each major upgrade, e.g., it will not upgrade straight
        from 1.0.0 to 3.0.0, but rather 1.0.0 -> 1.1.0 -> 2.0.0 -> 3.0.0

        Args:
            consent (bool, optional): variable to ensure consent to update database is explicit

        Returns:
            bool: final consent
        """
        self.conn = self._create_connection()
        # get consent, same for both
        if not consent:
            logger.warning(
                "WARNING: All existing filters and bookmarks in database will be dropped during database update!"
            )
            consent = input("Type 'yes' if you wish to continue: ") == "yes"
        if not consent:
            logger.critical("Consent not given for database update. Cancelling...")
            sys.exit(1)

        original_version = self.check_ringtaildb_version()[1]
        print(
            f"Upgrading {self.db_file} of version {original_version} to version {new_version}:"
        )

        # upgrade to 1.1.0
        if original_version in ["1.0.0", "1.1.0"]:
            logger.warning(
                "If you created the database with the duplicate handling option, there is a chance of inconsistent behavior of anything involving interactions as the Pose_ID was not used as an explicit foreign key in db v1.0.0 and v1.1.0."
            )
            if original_version == "1.0.0":
                self._update_db_100_to_110()
                print("\n\nSuccessfully upgraded to 1.1.0!\n\n")

            # upgrade to 2.0.0
            if new_version in ["2.0.0", "3.0.0"]:
                self._update_db_110_to_200()
                print("\n\nSuccessfully upgraded to 2.0.0!\n\n")

        # upgrade to 3.0.0
        if new_version == "3.0.0" and original_version == "2.0.0":
            self._update_db_200_to_300()
            print("\n\nSuccessfully upgraded to 3.0.0!\n\n")

        return consent

    def _update_db_100_to_110(self):
        """
        Will update a database of version 1.0.0 to 1.1.0, which renames energies_binding to docking_score in Results, adds a column in Bookmarks to store filters dict,
        and add indices to Results and Interaction_indices

        Raises:
            DatabaseConnectionError
        """
        self._drop_views()
        # create cursor
        cur = self.conn.cursor()
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
            self._set_ringtail_db_schema_version("1.1.0")
        except sqlite3.OperationalError as e:
            raise DatabaseConnectionError(
                f"Error while updating database from v1.0.0 to v1.1.0: {e}"
            ) from e

    def _update_db_110_to_200(self):
        """
        Method to update from database v 1.1.0 to 2.0.0,mainly removes the bitvetor table and creates Interactions table
        where interaction just lists Pose_id and interaction_id in a long-skinny table

        Raises:
            DatabaseConnectionError
            StorageError
        """
        self._drop_views()
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
            self.conn.commit()
            # index certain tables
            self._create_indices()
            self._set_ringtail_db_schema_version("2.0.0")  # set explicit version
        except sqlite3.OperationalError as e:
            raise DatabaseConnectionError(
                f"Error while creating new interaction tables: {e}"
            ) from e
        except StorageError as e:
            raise StorageError(
                f"Error while setting the database schema version: {e}"
            ) from e

    def _update_db_200_to_300(self):
        """
        Upgrades database from 2.0.0 to 3.0.0.
        This includes
        - converting the a chemicalite object in the Ligands table to a serialized blob (removing chemicalite dependency)
        - giving ligands a ligand_id which is used in Results instead of LigName
        - removes the use of views for storing filtered data, instead adds a Filtered_poses table to store all passing poses
        - keeps bookmark table but gives each bookmark an id which is used in the Filtered_poses table
        - removes some of the rarely used indices and adds a few others for minimizing db file size

        Raises:
            StorageError
        """
        # drop views first, because they depend on tables that will be altered
        self._drop_views()
        self._delete_table("Bookmarks")
        # create Filter and filtered_poses tables
        self._create_filtering_tables()
        # drop indices
        indices = self.db_query(
            "SELECT name FROM sqlite_master WHERE type == 'index'"
        ).fetchall()
        for index in indices:
            index_name = index[0]
            try:
                self.db_query(f"DROP INDEX {index_name}")
            except:
                pass
        # create new, empty ligands table
        self._create_ligands_table("Ligands_new")

        # create a temp connection function
        def _smile_to_rdbin(smile):
            """
            Temporary db connection method that will use rdkit to convert smiles to Mol
            inline in the sql query

            Args:
                smile (str): smiles describing ligand

            Returns:
                blob: binary Chem.rdchem.Mol ready to insert in db
            """
            try:
                mol = Chem.MolFromSmiles(smile)
                if mol is None:
                    return None
                Chem.SanitizeMol(mol)
                return mol.ToBinary()
            except Exception:
                return None

        self.conn.create_function("smile_to_rdbin", 1, _smile_to_rdbin)

        # populate with data from original ligands table, will autogenerate ligand_id PK
        self.db_query(
            """INSERT INTO Ligands_new (
                LigName,
                ligand_smile,
                rdmol,
                atom_index_map,
                hydrogen_parents,
                input_model) 
            SELECT 
                LigName,
                ligand_smile,
                smile_to_rdbin(ligand_smile),
                atom_index_map,
                hydrogen_parents,
                input_model FROM Ligands;""",
            commit=True,
        )
        # ensure row numbers are the same
        original_length = self.table_length("Ligands")
        new_length = self.table_length("Ligands_new")
        if original_length != new_length:
            raise StorageError(
                "Problems while upgrading database, Ligands table did not copy properly."
            )
        # delete old table
        self._delete_table("Ligands")
        # rename new table
        self.db_query("ALTER TABLE Ligands_new RENAME TO Ligands;", commit=True)

        # update results table to use ligand_id from Ligands
        self.db_query("ALTER TABLE Results ADD COLUMN ligand_id INTEGER;")
        # populate ligand_id in Results
        self.db_query(
            """UPDATE Results
                        SET ligand_id = (
                            SELECT ligand_id FROM Ligands
                            WHERE Ligands.LigName = Results.LigName);"""
        )

        # create new Results table without LigName column
        self._create_results_table("Results_new")
        # insert data from original to new Results
        self.db_query(
            """INSERT INTO Results_new (
                        pose_id,
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
                        flexible_res_coordinates)
                    SELECT
                        pose_id,
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
                      FROM Results""",
            commit=True,
        )
        # ensure row numbers are the same
        original_length = self.table_length("Results")
        new_length = self.table_length("Results_new")
        if original_length != new_length:
            raise StorageError(
                "Problems while upgrading database, Results table did not copy properly."
            )
        # delete old table
        self._delete_table("Results")
        # rename new table
        self.db_query("ALTER TABLE Results_new RENAME TO Results;", commit=True)

        # build new indices if you got this far successfully
        self._create_indices()
        self.db_query("REINDEX", commit=True)
        self._set_ringtail_db_schema_version("3.0.0")

    def _delete_filter_data(self, db_alias: str = None):
        """
        Empties all data in the filter and filtered_poses tables

        Args:
            db_alias (str, optional): if needing to empty tables from a connected, aliased database. Defaults to None.
        """
        if db_alias:
            alias_string = db_alias + "."
        else:
            alias_string = ""
        self.db_query(f"DELETE FROM {alias_string}Filters")
        # delete all rows in bookmarks table
        self.db_query(f"DELETE FROM {alias_string}filtered_poses", commit=True)

    def _drop_views(self, db_alias: str = None):
        """
        Will drop views and clear bookmark table

        Args:
            db_alias (str, optional): if needing to drop views from a connected, aliased database. Defaults to None.
        """
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

    def close_storage(self, attached_db=None, vacuum=None):
        """
        Closes storage

        Args:
            attached_db (str, optional): alias of attached database. Defaults to None.
            vacuum (bool, optional): whether or not to vacuum file to save space. Defaults to None.
        """
        if attached_db or vacuum:
            self._cleanup_storage(attached_db, vacuum)

        # close db itself
        logger.debug("Closing database")
        self.conn.close()

    def _cleanup_storage(
        self, attached_db_alias: str = None, vacuum: bool = None, reindex=False
    ):
        """
        Cleans up storage, especially useful for situations where two databases
        hvae been connected/attached. Will detach the database, reindex main database,
        and vacuum the file if requested.

        Args:
            attached_db (str, optional): Name of attached database, if any. Defaults to None.
            vacuum (bool, optional): rebuilds db file to minimize space. Defaults to None.
            reindex (bool, optional): deletes and reruns all indixes. Defaults to False.

        """
        if attached_db_alias is not None:
            self._detach_db(attached_db_alias)
        if reindex:
            self.db_query("REINDEX", commit=True)
        # vacuum database
        if vacuum:
            self._vacuum()

    def db_empty(self):
        """empty database, for example if overwrite

        Returns:
            bool: whether or not db is empty
        """
        cur = self.conn.execute(
            "SELECT COUNT(*) name FROM sqlite_master WHERE type='table' AND name <> 'sqlite_sequence';"
        )
        tablecount = cur.fetchone()[0]
        if tablecount > 0:
            cur.execute("SELECT COUNT(*) FROM Results")
            datacount = cur.fetchone()[0]
        else:
            datacount = 0
        cur.close()
        return True if datacount == 0 else False

    def table_length(self, table: str) -> int:
        """
        Get length of table or bookmark

        Args:
            table (str): name of table or bookmark

        Returns:
            int: number of poses in table/bookmark
        """
        if self._is_table(table):
            query = f"""SELECT COUNT(*) FROM {table};"""
            params = ()
        elif self.is_bookmark(table):
            query = """SELECT COUNT(*) FROM Filtered_poses WHERE filter_id = (SELECT filter_id FROM Filters WHERE name = ?);"""
            params = (table,)
        else:
            logger.error(f"Table -{table}- does not exist in the database.")
            return None

        return self.db_query(query, params).fetchone()[0]

    def pose_row_in_table(self, table: str, pose_id: int) -> Union[None, int]:
        """
        Find the row id of a pose in a given table

        Args:
            table (str)
            pose_id (int)

        Returns:
            Union[None, int]: rowid if any
        """
        query = self.QueryBuilder()
        query.SELECT("rowid")
        if self.is_bookmark(table):
            query.FROM("Filtered_poses").WHERE(
                "filter_id = (SELECT filter_id from Filters WHERE name = ?)", table
            )
        else:
            query.FROM(table)
        query.WHERE("pose_id = ?", pose_id)
        row = self.db_query(*query.build()).fetchone()
        if row:
            return row[0]
        else:
            return None

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
        return self.db_query(query, commit=True)

    def db_query(self, query, params: tuple = (), commit=False) -> sqlite3.Cursor:
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
            if commit:
                self.conn.commit()
        except sqlite3.OperationalError as e:
            raise DatabaseQueryError(
                f"Unable to execute query -{query}- with given parameters -{params}-: -{e}-"
            ) from e
        return cur

    def db_update(self, query: str, parameters: list[tuple], commit=True) -> iter:
        """
        A db query that also commits if/when specified

        Args:
            query (str): sqlite formatted query string
            parameters (list[tuple]): assumes appropriate place holders in query
            commit (bool, optional): whether to commit the transaction in open connection. Defaults to True.

        Raises:
            OptionError
            DatabaseInsertionError

        Returns:
            iter: if requesting return value(s)
        """
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

    def _is_table(self, table: str) -> bool:
        """
        Returns True if table name is actually a bookmark

        Args:
            table (str): name of table or bookmark to check

        Returns:
            bool: if table name is a bookmark
        """

        if table.lower() in self.tables_in_db():
            return True
        else:
            return False

    def is_bookmark(self, table: str) -> bool:
        """
        Returns True if table name is actually a bookmark

        Args:
            table (str): name of table or bookmark to check

        Returns:
            bool: if table name is a bookmark
        """

        if table and table.lower() in self.get_all_bookmark_names():
            return True
        else:
            return False

    def _is_statustable(self, table: str) -> bool:
        """
        Returns True if table name is actually a status table (table with poses who have been assigned a status like accept, reject, maybe)

        Args:
            table (str): name of table or bookmark to check

        Returns:
            bool: if table name is a status table
        """
        # TODO hardcoded, need to find better way
        if table.lower() in ["accepted", "maybe", "rejected"]:
            return True
        else:
            return False

    # endregion

    # region GUI specific API
    def fetch_viewable_data_columns_from(
        self, table: str, length: int, starting_rowid: int = 0, reverse=False
    ) -> dict[list[str], list]:
        """
        Makes a selection of columns and includes the status of the pose

        Returns:
            dict[list[str], list]: dict of headers and data
        """

        if reverse:
            where_operator = "<="
        else:
            where_operator = ">="
        query = self.QueryBuilder()
        query.FROM("Results", "R")
        status_assignement = """CASE
            WHEN EXISTS (SELECT 1 FROM Accepted s WHERE s.pose_id = R.pose_ID) THEN 'accepted'
            WHEN EXISTS (SELECT 1 FROM Rejected s WHERE s.pose_id = R.pose_ID) THEN 'rejected'
            WHEN EXISTS (SELECT 1 FROM Maybe s WHERE s.pose_id = R.pose_ID) THEN 'maybe'
            ELSE 'not evaluated'
        END AS status,"""

        if table.lower() == "results":
            rowid = "R.rowid"

        elif self._is_statustable(table):
            query.JOIN(table, "T", "pose_id")
            rowid = "T.rowid"
            # status assignement doesn't make sense for status tables
            status_assignement = f"""'{table.lower()}',"""

        elif self.is_bookmark(table):
            query.JOIN("filtered_poses", "fp", "pose_id").JOIN(
                "filters", "f", "filter_id", "filtered_poses"
            ).WHERE("f.name = ?", table)
            rowid = "fp.rowid"

        ordered_columns = f"""
        {status_assignement}
        R.Pose_ID, L.LigName, R.docking_score, 
        R.leff, R.cluster_size, R.cluster_rmsd, 
        R.pose_rank, R.num_hb, R.receptor, R.run_number, 
        R.deltas, R.nr_interactions, R.unbound_energy, 
        R.reference_RMSD, R.energies_inter, R.energies_vdw, 
        R.energies_electro, R.energies_flexLig, R.energies_flexLR, 
        R.energies_intra, R.energies_torsional, R.about_x, R.about_y, 
        R.about_z, R.trans_x, R.trans_y, R.trans_z, R.axisangle_x, 
        R.axisangle_y, R.axisangle_z, R.axisangle_w, R.dihedrals, {rowid}"""

        query.SELECT(ordered_columns)

        query.JOIN("Ligands", "L", "ligand_id", "results").WHERE(
            f"{rowid} {where_operator} ?", starting_rowid
        ).ORDER_BY(rowid).LIMIT(length).DESC(reverse)

        cursor = self.db_query(*query.build())
        headers = [desc[0] for desc in cursor.description]
        data = cursor.fetchall()
        return {"headers": headers, "data": data}

    def get_starting_rowid(self, table: str) -> int:
        """
        Starting row id for a table, will be 1 for regular tables, and 1 or non-1 for bookmarks
        (whose rows are inside Filtered_poses)

        Args:
            table (str): table or bookmark name

        Returns:
            int: first row id belonging to that selection
        """
        query = self.QueryBuilder()
        query.SELECT("MIN(rowid)")

        if self._is_table(table):
            query.FROM(table)
        elif self.is_bookmark(table):
            query.FROM("Filtered_poses").WHERE(
                "filter_id = (SELECT filter_id FROM Filters WHERE name = ?)", table
            )

        else:
            logger.error(f"Table -{table}- does not exist in the database.")
            return None
        return self.db_query(*query.build()).fetchone()[0]

    def fetch_columns_from_table_as_dicts(
        self, table: str, columns: list, length: int = 500, starting_rowid: int = 0
    ) -> tuple[list[str], list[dict]]:
        """
        Will get requested table data for a table given one or more columns.
        Data will be limited by a certain length, and can be retrieved from a desired
        rowid.

        Args:
            table (str): name of table or bookmark
            columns (list, optional): list of columns to retrieve. Defaults to ["*"].
            length (int, optional): number of rows to collect. Defaults to 500.
            starting_rowid (int, optional): rowid to start with. Defaults to 0.

        Returns:
            tuple[list[str], list[dict]]: list of column names, and list of dicts where each dict is one row,
                                            and column is the key, value is the row-col cell value
        """
        query = self.QueryBuilder()
        query.SELECT(",".join(columns)).FROM(table)

        if length:
            query.LIMIT(length)
        if starting_rowid:
            query.WHERE(f"{table}.rowid = {starting_rowid}")

        return self.get_query_data_as_dicts(query.build()[0])

    def fetch_lignames_and_poses_for_selection(
        self, selection: str
    ) -> dict[str, list[int]]:
        """
        Creates a dictionary of ligands and the selected poses that appear in a selection,
        such as a bookmark or a status table (where only poses are given).

        Args:
            selection (str): name of bookmark or status table

        Returns:
            dict[str, list[int]]: ligand name is keyword, value is list of poses in given selection
        """
        query = self.QueryBuilder()
        query.SELECT("L.Ligname", "r.pose_id").FROM("Ligands", "L").JOIN(
            "Results", "R", "ligand_id"
        )
        if self.is_bookmark(selection):
            query.IN_BOOKMARK(selection)
        elif selection.lower() in statuses:
            query.WHERE(f"R.pose_id IN {selection}")
        else:
            logger.error(
                f"-{selection}- is not a valid selection for this method. Please provide a bookmark_name or a status table."
            )
            return
        rows = self.db_query(query.build()[0]).fetchall()
        ligand_poses = defaultdict(list)
        for name, id in rows:
            ligand_poses[name].append(id)
        # convert to normal dict
        return dict(ligand_poses)

    def fetch_selected_ligand_poses(self, ligand_name: str, selection: str):
        """
        Gets only the poses of a given ligand that are present in give selection (e.g., a bookmark)

        Args:
            ligand_name (str)
            selection (str): status table or bookmark name

        Returns:
            list[int]: selected poses for ligand
        """
        query = self.QueryBuilder()
        query.SELECT("R.Pose_id").FROM("Results", "R").WHERE(
            "L.LigName = ?", ligand_name
        ).JOIN("Ligands", "L", "ligand_id")

        if self.is_bookmark(selection):
            query.IN_BOOKMARK(selection)
        elif selection == None:
            # get all poses for ligand, no WHERE clause for pose_id
            pass
        elif selection.lower() in statuses:
            query.WHERE(f"R.Pose_ID IN (SELECT Pose_ID FROM {selection})")
        else:
            logger.error(
                f"-{selection}- is not a valid selection for this method. Please provide a bookmark_name or a status table."
            )
            return
        return [row[0] for row in self.db_query(*query.build()).fetchall()]

    def get_query_data_as_dicts(self, query: str) -> tuple[list[str], list[dict]]:
        """
        Will return data requested in an sqlite formatted query

        Args:
            query (str): sql query formatted to sqlite database

        Returns:
            tuple[list[str], list[dict]]: list of column names, and list of dicts where each dict is one row,
                                            and column is the key, value is the row-col cell value
        """
        rows = self.db_query(query).fetchall()
        column_names = rows[0].keys() if rows else []
        return list(column_names), [dict(row) for row in rows]

    def create_status_tables(self) -> None:
        """
        Creates pose status tables if needed
        """
        self.db_update(
            f"""
                CREATE TABLE IF NOT EXISTS Accepted 
                (Pose_ID INTEGER UNIQUE,
                FOREIGN KEY (Pose_ID) REFERENCES Results(Pose_ID)
                );""",
            (),
        )

        # create maybe table
        self.db_update(
            f"""
                CREATE TABLE IF NOT EXISTS Maybe 
                (Pose_ID INTEGER UNIQUE,
                FOREIGN KEY (Pose_ID) REFERENCES Results(Pose_ID)
                );""",
            (),
        )

        # create rejected table
        self.db_update(
            f"""
                CREATE TABLE IF NOT EXISTS Rejected 
                (Pose_ID INTEGER UNIQUE,
                FOREIGN KEY (Pose_ID) REFERENCES Results(Pose_ID)
                );""",
            (),
        )

        return None

    def accept_pose(self, pose_id: int):
        """
        Will add pose_id to accepted, and delete from maybe and rejected if needed

        Args:
            pose_id (int)
        """
        self.db_update(
            """INSERT OR IGNORE INTO Accepted (pose_id) VALUES (?);""",
            (pose_id,),
            commit=False,
        )
        self.db_update(
            """DELETE FROM Maybe WHERE Pose_id = ?;""", (pose_id,), commit=False
        )
        self.db_update(
            """DELETE FROM Rejected WHERE Pose_id = ?;""", (pose_id,), commit=True
        )

    def maybe_pose(self, pose_id: int):
        """
        Will add pose_id to maybe, and delete from accepted and rejected if needed

        Args:
            pose_id (int)
        """
        self.db_update(
            """INSERT OR IGNORE INTO Maybe (pose_id) VALUES (?);""",
            (pose_id,),
            commit=False,
        )
        self.db_update(
            """DELETE FROM Accepted WHERE Pose_id = ?;""", (pose_id,), commit=False
        )
        self.db_update(
            """DELETE FROM Rejected WHERE Pose_id = ?;""", (pose_id,), commit=True
        )

    def reject_pose(self, pose_id: int):
        """
        Will add pose_id to rejected, and delete from accepted and maybe if needed

        Args:
            pose_id (int)
        """
        self.db_update(
            """INSERT OR IGNORE INTO Rejected (pose_id) VALUES (?);""",
            (pose_id,),
            commit=False,
        )
        self.db_update(
            """DELETE FROM Accepted WHERE Pose_id = ?;""", (pose_id,), commit=False
        )
        self.db_update(
            """DELETE FROM Maybe WHERE Pose_id = ?;""", (pose_id,), commit=True
        )

    # endregion
