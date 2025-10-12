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
from signal import signal, SIGINT
from rdkit import Chem
from typing import Union
import time
from importlib.metadata import version
from .ringtailoptions import Filters, statuses
from .exceptions import StorageError, OptionError, VersionError
from .clustermanager import *
from .querybuilder import QueryBuilder
from collections import defaultdict


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

    # region setup
    def __init__(self):
        self.keyboard_interrupt_allowed = False
        self.db_file: str

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

    # endregion

    # region public api

    @classmethod
    def format_for_storage(cls, ligand_dict: dict) -> dict:
        """Takes file dictionary from the file parser, formats required storage format. Only handles docking data for one ligand at the time.
        For each run we save, we add its interaction dict to the interaction_dictionaries list and save its other data
        We also save a mapping of the its cluster number to the index in interaction_dictionaries
        Then, when we find a pose to tolerate interactions for, we lookup the index to append the interactions to from cluster_map_to_run_pose
        Finally, we calculate the interaction tuple lists for each pose

        Args:
            ligand_dict (dict): Dictionary containing data from the fileparser

        Returns:
            dict: with storage formatted rows, including: results rows per pose, interactions per pose, one ligand row, and one receptor row

        """
        result_rows = []
        interaction_dictionaries = []
        cluster_map_to_run_pose = {}
        ligand_row = cls._generate_ligand_row(ligand_dict)

        # iterates and essentially creates pose rows
        for idx, run_number in enumerate(ligand_dict["sorted_runs"]):
            cluster = ligand_dict["cluster_list"][idx]
            pose_rank = idx + 1
            # save everything if this is a cluster top pose
            current_interactions = {
                "ligand_name": ligand_row[0],
            }
            if run_number in ligand_dict["poses_to_save"]:
                cluster_map_to_run_pose[cluster] = {
                    "run_number": run_number,
                    "pose_rank": pose_rank,
                }
                current_interactions.update(cluster_map_to_run_pose[cluster])
                # Check how things are parsed here, might not be most efficient
                result_rows.append(
                    cls._generate_results_row(ligand_dict, idx, run_number)
                )
                if ligand_dict["interactions"] != []:
                    # here is where I have to add run number ligname and rank
                    current_interactions.update(ligand_dict["interactions"][idx])
                    interaction_dictionaries.append(current_interactions)
            elif run_number in ligand_dict["tolerated_interaction_runs"]:
                # adds to list started by best-scoring pose in cluster
                if cluster not in cluster_map_to_run_pose.keys():
                    continue
                else:
                    current_interactions.update(
                        {
                            "run_number": cluster_map_to_run_pose[cluster][
                                "run_number"
                            ],
                            "pose_rank": cluster_map_to_run_pose[cluster]["pose_rank"],
                        }
                    )
                current_interactions.update(ligand_dict["interactions"][idx])
                interaction_dictionaries.append(current_interactions)

        interactions_tuples = cls._generate_interaction_tuples(interaction_dictionaries)

        data_dict = {
            "ligand": ligand_row,
            "poses": result_rows,
            "interactions": interactions_tuples,
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
        self._insert_ligands(docking_data["ligands"])
        interaction_data = docking_data["interactions"]
        # deduplicate by using a set comprehension, then convert to list
        just_interactions = list({interaction[3:] for interaction in interaction_data})

        self._insert_interaction_index_rows(just_interactions)

        self._insert_docking_data(
            docking_data["poses"], interaction_data, write_options
        )
        self.conn.commit()

    def insert_receptor_basic_info(self, receptor_data: list):
        """

        Method to insert receptor information into the database

        Args:
            receptor_data (list): of receptor data, ordered by columns in the db
        """
        receptors = self.fetch_receptor_object()
        # insert receptor if database does not have already have a receptor entry
        if not receptors.get("RecName"):
            self._insert_receptors(receptor_data)

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
        """Get all bookmarks in sql database as a list of names.

        Returns:
            list: of bookmark names
        """
        query = self.QueryBuilder()
        query.SELECT("name").FROM("Filters")
        cur = self.db_query(query.build()[0])
        bookmark_names = [row[0].lower() for row in cur.fetchall()]

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
        # make sure bookmark needs deleting
        if not self.is_bookmark(bookmark_name):
            return
        # get filter id
        query = self.QueryBuilder()
        query.SELECT("filter_id").FROM("Filters").WHERE("name=?", bookmark_name)

        filter_id = self.conn.execute(*query.build()).fetchone()[0]
        # delete from filtered_poses table
        query = self.QueryBuilder()
        self.conn.execute(
            *query.DELETE_FROM("Filtered_poses")
            .WHERE("filter_id = ?", filter_id)
            .build()
        )
        self.conn.commit()
        # delete from filters
        query = self.QueryBuilder()
        self.conn.execute(
            *query.DELETE_FROM("Filters").WHERE("filter_id = ?", filter_id).build()
        )
        self.conn.commit()
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
        logger.debug("Running interaction union query")
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
        num_clusters = len(clusters)
        logger.info(f"Number of {cluster_type} clusters: {num_clusters}")
        logger.info(
            f"Biggest {cluster_type} cluster contains {max_len} poses while the smallest cluster contains {min_len} poses."
        )
        # write to db
        cluster_bookmark_name = self._insert_cluster_data(
            clusters,
            [int(item) for item in representatives],
            cluster_type.lower(),
            str(cutoff),
            bookmark_name,
        )
        if type(cluster_bookmark_name) == tuple:
            logger.info("Clustering has been ran before, old bookmark will be used.")
            num_clusters = cluster_bookmark_name[1]
            cluster_bookmark_name = cluster_bookmark_name[0]
        else:

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

        return cluster_bookmark_name, num_clusters

    def finalize_database_write(self):
        """
        Methods to finalize when a database has been written to, including creating indices
        """
        # index certain tables
        self._create_indices()
        logger.info("Database write session completed successfully.")

    def create_subset_database(self, bookmark_name: str, database_name: str):
        """
        Creates an empty database with all tables, attaches to main database,
        and populates tables based on the poses present in bookmark

        Args:
            bookmark_name (str): filter used to determine poses
            database_name (str): name of new database

        Raises:
            StorageError:
        """
        # create the database
        new_db = self.__class__(database_name)
        # make tables
        with new_db:
            new_db._create_tables(commit=False)
        logger.debug(f"Created a new database with empty tables: {database_name}.")
        # attach the incoming database
        alias = "subset"
        self._attach_db(database_name, alias)
        logger.debug(
            f"Database {database_name} attached to main database {self.db_file}."
        )
        # ligands
        ligands = f"""
        INSERT INTO {alias}.Ligands
        SELECT * FROM main.Ligands
        WHERE main.Ligands.ligand_id IN (
            SELECT ligand_id from main.Results
            WHERE pose_id IN (
                SELECT Pose_id FROM main.filtered_poses
                WHERE filter_id =
                    (SELECT filter_id FROM main.Filters
                    WHERE name = '{bookmark_name}'))
                    );"""
        self.db_query(ligands)
        logger.debug("Ligands have been copied into the new subset database.")
        # receptor
        receptor = f"""
        INSERT INTO {alias}.Receptors
        SELECT * FROM main.Receptors;
        """
        self.db_query(receptor)
        logger.debug("The receptor have been copied into the new subset database.")
        # results
        poses = f"""
        INSERT INTO {alias}.Results
        SELECT * FROM main.Results
        WHERE main.Results.pose_id IN (
            SELECT Pose_id FROM main.filtered_poses
            WHERE filter_id =
                (SELECT filter_id FROM main.Filters
                WHERE name = '{bookmark_name}')
                );"""
        self.db_query(poses)
        logger.debug("Results have been copied into the new subset database.")
        # interaction_indices
        interaction_indices = f"""
        INSERT INTO {alias}.Interaction_indices
        SELECT * FROM main.Interaction_indices
        WHERE main.Interaction_indices.interaction_id IN (
            SELECT interaction_id FROM main.Interactions
            WHERE pose_id IN (
                SELECT Pose_id FROM main.filtered_poses
                WHERE filter_id =
                    (SELECT filter_id FROM main.Filters
                    WHERE name = '{bookmark_name}'))
            );"""
        self.db_query(interaction_indices)
        # interactions
        interactions = f"""
        INSERT INTO {alias}.Interactions
        SELECT * FROM main.Interactions
        WHERE main.Interactions.pose_id IN (
            SELECT Pose_id FROM main.filtered_poses
            WHERE filter_id =
                (SELECT filter_id FROM main.Filters
                WHERE name = '{bookmark_name}')
                );"""
        self.db_query(interactions)
        logger.debug("Interactions have been copied into the new subset database.")
        try:
            self.conn.commit()
        except Exception as e:
            raise StorageError("Problems while creating a subset database: ", str(e))

        self._detach_db(alias)
        logger.info(f"Subset database {database_name} has been successfully created.")

    def fetch_pose_interactions(self, Pose_ID) -> iter:
        """
        Fetch all interactions parameters belonging to a Pose_ID

        Args:
            Pose_ID (int): pose id, 1-1 with Results table

        Returns:
            iter: of interaction information for given Pose_ID
        """
        # check if table exist
        if not "Interactions" in self.tables_in_db():
            return

        query = self.QueryBuilder()
        query.SELECT(
            "ii.interaction_type",
            "ii.rec_chain",
            "ii.rec_resname",
            "ii.rec_resid",
            "ii.rec_atom",
            "ii.rec_atomid",
        ).FROM("Interaction_indices", "ii").JOIN(
            "Interactions", "i", "interaction_id"
        ).WHERE(
            "i.pose_id = ?", Pose_ID
        )

        return self.db_query(*query.build()).fetchall()

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
            query.GROUP_BY("l.ligname")
        if order_results:
            order_by = self._format_orderby(order_results)
            if order_by:
                query.ORDER_BY(order_by)
        return query.build()[0]

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
            query.SELECT("*").FROM("Results")
            if self.is_bookmark(requested_data):
                query.IN_BOOKMARK(requested_data)
            elif self._is_statustable(requested_data):
                query.WHERE(f"pose_id in {requested_data}")
            return pd.read_sql_query(query.build()[0], self.conn)
        else:
            return pd.read_sql_query(requested_data, self.conn)

    def get_query_data_as_dicts(self, query: str) -> tuple[list[str], list[dict]]:
        """
        Will return data requested in an duckdb formatted query

        Args:
            query (str): sql query formatted to duckdb database

        Returns:
            tuple[list[str], list[dict]]: list of column names, and list of dicts where each dict is one row,
                                            and column is the key, value is the row-col cell value
        """
        cur = self.db_query(query)
        rows = cur.fetchall()
        column_names = [desc[0] for desc in cur.description] if cur.description else []
        dict_rows = [dict(zip(column_names, row)) for row in rows]
        return column_names, dict_rows

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
        query = self.QueryBuilder()
        query.SELECT("docking_mode").FROM("DB_properties").ORDER_BY(
            "DB_write_session"
        ).DESC("DB_write_session").LIMIT(1)
        docking_mode = self.db_query(query.build()[0]).fetchone()
        return docking_mode[0].lower() if docking_mode else None

    def check_storage_ready(self, run_mode: str, docking_mode: str, num_poses: int):
        """
        Check that storage is ready

        Args:
            run_mode (str): _description_
            docking_mode (str): _description_
            num_poses (int): _description_

        Raises:
            OptionError: if database options are not compatible
        """
        if self.db_empty():
            self._create_tables()
        query = self.QueryBuilder()
        query.SELECT("COUNT(*)").FROM("DB_properties")
        count = self.db_query(query.build()[0]).fetchone()[0]

        compatible = True
        if count < 1:
            logger.info(
                "Adding results to an existing database that is currently empty of docking results."
            )
        else:
            compatibility_string = "The following database properties do not agree with the properties last used for this database: \n"
            try:
                query = self.QueryBuilder()
                query.SELECT("*").FROM("DB_properties").ORDER_BY(
                    "DB_write_session"
                ).DESC("DB_write_session").LIMIT(1)
                cur = self.db_query(query.build()[0])

                (_, last_docking_mode, last_num_poses) = cur.fetchone()
                if docking_mode != last_docking_mode:
                    compatible = False
                    compatibility_string += f"Current docking mode is {docking_mode} but last used docking mode of database is {last_docking_mode}.\n"

                if last_num_poses == -1 != num_poses:
                    compatible = False
                    compatibility_string += f"Current number of poses saved is {num_poses} but database was previously set to 'store_all_poses'.\n"
                elif last_num_poses != num_poses:
                    compatible = False
                    compatibility_string += f"Current number of poses saved is {num_poses} but database was previously set to {last_num_poses}."
            except Exception as e:
                raise e

        # command line does not allow incompatibility, but API (and I suppose GUI) does
        if not compatible:
            if run_mode == "cli":
                raise OptionError(compatibility_string)
            else:
                logger.warning(compatibility_string)

        # write current database properties to database
        self._insert_db_properties(docking_mode, num_poses)
        logger.debug("Storage compatibility has been checked and is ensured.")
        # cannot use Signal/keyboard interrupt in the GUI bc it uses threading
        if run_mode != "gui":
            self.keyboard_interrupt_allowed = True

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

    def db_empty(self):
        """Checks if database is empty (has rows in Results)

        Returns:
            bool: whether or not db is empty
        """
        if self.tables_in_db():
            datacount = self.table_length("Results")
        else:
            datacount = 0
        return True if datacount == 0 else False

    def count_receptors_in_db(self):
        """returns number of rows in Receptors table where receptor_object already has blob

        Returns:
            int: number of rows in receptors table
            str: name of receptor if present in table

        Raises:
            DatabaseQueryError
        """
        query = self.QueryBuilder()
        row_count = self.db_query(
            *query.SELECT("COUNT(*)")
            .FROM("Receptors")
            .WHERE("receptor_object NOT NULL")
            .build()
        ).fetchone()[0]
        return row_count

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
        bin_edges = np.percentile(values, 100 - bins)
        return bins, bin_edges

    def fetch_data_for_passing_results(
        self, bookmark_name: str, outfields: Union[str, list], order_results: str = None
    ) -> iter:
        """Will return duckdb cursor with requested data for outfields for poses that passed filter in bookmark_name

        Returns:
            iter: duckdb cursor of data from passing data

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
        ).JOIN("ligands", "L", "ligand_id", "results").GROUP_BY("L.ligname")
        if order_results:
            order_by = self._format_orderby(order_results)
            if order_by:
                query.ORDER_BY(order_by)

        return self.db_query(query.build()[0]).fetchall()

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

    def fetch_receptor_object(self) -> Union[None, dict]:
        """Returns all Receptor objects from database

        Returns:
            dict: of receptor name and object and/or polymer json if column exist
        """

        columns = ["RecName", "receptor_object"]

        if self._check_if_column_in_table("Receptors", "polymer"):
            columns.append("polymer")

        query = self.QueryBuilder()
        query.SELECT(*columns).FROM("Receptors")
        # check if polymer column exist
        row = self.db_query(query.build()[0]).fetchone()
        if not row:
            return {}
        else:
            return dict(zip(columns, row))

    def fetch_flexres_info(self, receptor: Union[str, int]):
        """fetch flexible residues names and atomname lists

        Args:
            receptor (Union[str, int]): receptor descriptor, either receptor_id or receptor name

        Returns:
            tuple: (flexible_residues, flexres_atomnames)
        """
        if type(receptor) == int:
            selection = "receptor_id = ?"
        elif type(receptor) == str:
            selection = "recname = ?"
        query = self.QueryBuilder()
        query.SELECT("flexible_residues", "flexres_atomnames").FROM("Receptors").WHERE(
            selection, receptor
        )
        info = self.db_query(*query.build()).fetchone()
        if info is None:
            info = [], []
        return info

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

    def fetch_rdkit_relevant_pose_properties(self, pose_ids: list) -> iter:
        """
        Gets molecular data that is needed to create rdkit mols for a given list of poses

        Args:
            pose_ids (list): pose ids for which to collect molecular data

        Returns:
            iter: of the following columns pose_id, docking_score, leff, ligand_coordinates, flexible_res_coordinates
        """
        placeholders = ",".join(["?"] * len(pose_ids))

        query = self.QueryBuilder()
        query.SELECT(
            "pose_id",
            "docking_score",
            "leff",
            "ligand_coordinates",
            "flexible_res_coordinates",
        ).FROM("Results").WHERE(f"Pose_ID IN ({placeholders})", *pose_ids)
        return self.db_query(*query.build()).fetchall()

    def get_gui_plot_data(
        self,
        bookmark_name: str,
        include_status: bool = False,
        x_axis: str = "docking_score",
        y_axis: str = "leff",
        limit: int = None,
    ):
        """This function gathers two docking results columns (docking score and ligand efficienct) from all data,
        as well as pose_id and ligand name from given bookmark. Can request the data just for poses in the bookmark.

        Args:
            bookmark_name (str): name of bookmark for which to fetch passing data. Returns empty list if bookmark does not exist.
            include_status (bool): look for status tables and include if requested
            x_axis (str, optional): Defaults to "docking_score".
            y_axis (str, optional): Defaults to "leff".

        Returns:
            tuple: cursors as (<all data cursor>, <passing data cursor>)
        """

        bookmark_query = self.QueryBuilder()
        bookmark_query.SELECT(
            "R." + x_axis,
            "R." + y_axis,
            "R." + "Pose_ID",
            "L." + "LigName",
            "L." + "ligand_smile",
        )
        if limit:
            bookmark_query.LIMIT(limit)

        if self.is_bookmark(bookmark_name):
            if include_status:
                bookmark_query.SELECT_STATUS()
            bookmark_query.FROM("Results", "R").IN_BOOKMARK(bookmark_name).JOIN(
                "Ligands", "L", "ligand_id"
            )

            data = self.db_query(bookmark_query.build()[0]).fetchall()

        elif self._is_table(bookmark_name) and bookmark_name.lower() != "results":
            # will assume it is a status table
            if include_status:
                bookmark_query.SELECT(f"""'{bookmark_name.lower()}' as status""")
            bookmark_query.FROM("Results", "R").JOIN(
                bookmark_name, "T", "pose_id"
            ).JOIN("Ligands", "L", "ligand_id")

            data = self.db_query(bookmark_query.build()[0]).fetchall()

        else:
            if include_status:
                bookmark_query.SELECT_STATUS()
            bookmark_query.FROM("Results", "R").JOIN("Ligands", "L", "ligand_id")

            data = self.db_query(bookmark_query.build()[0]).fetchall()

        return data

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
            ELSE ''
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

    def fetch_bookmark_interactions(self, bookmark_name: str) -> pd.DataFrame:
        # TODO
        query = f"""
        SELECT DISTINCT 
            II.interaction_type,
            II.rec_chain,
            II.rec_resname,
            II.rec_resid,
            II.rec_atom
        FROM Interaction_indices AS II
        JOIN Interactions AS I ON I.interaction_id=II.interaction_id
        WHERE I.pose_id IN (
            SELECT Pose_id FROM filtered_poses 
                WHERE filter_id = 
                    (SELECT filter_id FROM Filters 
                    WHERE name = '{bookmark_name}')
                    );"""
        # iterable of tuples
        interactions = self.db_query(query).fetchall()
        interactions_df = pd.DataFrame(
            interactions,
            columns=[
                "interaction_type",
                "rec_chain",
                "rec_resname",
                "rec_resid",
                "rec_atom",
            ],
        )

        return interactions_df

    def accept_pose(self, pose_ids: Union[int, list[int]]):
        """
        _summary_

        Args:
            pose_ids (Union[int, list[int]]): _description_
        """
        """
        Will add pose_id to accepted, and delete from maybe and rejected if needed

        Args:
            pose_ids (Union[int, list[int]])
        """
        if isinstance(pose_ids, int):
            pose_ids = [pose_ids]
        pose_ids = [(pose,) for pose in pose_ids]

        self.db_update(
            """INSERT OR IGNORE INTO Accepted (pose_id) VALUES (?);""",
            pose_ids,
            commit=False,
        )
        self.db_update(
            """DELETE FROM Maybe WHERE Pose_id = ?;""", pose_ids, commit=False
        )
        self.db_update(
            """DELETE FROM Rejected WHERE Pose_id = ?;""", pose_ids, commit=True
        )

    def maybe_pose(self, pose_ids: Union[int, list[int]]):
        """
        Will add pose_id to maybe, and delete from accepted and rejected if needed

        Args:
            pose_ids (Union[int, list[int]])
        """
        if isinstance(pose_ids, int):
            pose_ids = [pose_ids]
        pose_ids = [(pose,) for pose in pose_ids]

        self.db_update(
            """INSERT OR IGNORE INTO Maybe (pose_id) VALUES (?);""",
            pose_ids,
            commit=False,
        )
        self.db_update(
            """DELETE FROM Accepted WHERE Pose_id = ?;""", pose_ids, commit=False
        )
        self.db_update(
            """DELETE FROM Rejected WHERE Pose_id = ?;""", pose_ids, commit=True
        )

    def reject_pose(self, pose_ids: Union[int, list[int]]):
        """
        Will add pose_id to rejected, and delete from accepted and maybe if needed

        Args:
            pose_ids (Union[int, list[int]])
        """
        if isinstance(pose_ids, int):
            pose_ids = [pose_ids]
        pose_ids = [(pose,) for pose in pose_ids]

        self.db_update(
            """INSERT OR IGNORE INTO Rejected (pose_id) VALUES (?);""",
            pose_ids,
            commit=False,
        )
        self.db_update(
            """DELETE FROM Accepted WHERE Pose_id = ?;""", pose_ids, commit=False
        )
        self.db_update(
            """DELETE FROM Maybe WHERE Pose_id = ?;""", pose_ids, commit=True
        )

    def prepare_column_export_query(self, columns: dict, bookmark: str) -> str:
        query = self.QueryBuilder()
        for _, column in columns.items():
            if column[0]:
                if column[1]:
                    # alias
                    query.SELECT(f"{column[0]} AS {column[1]}")
                else:
                    query.SELECT(column[0])
        # not the most flexible way to allow column selection
        query.FROM("Results")
        query.JOIN("Ligands", "Ligands", on="ligand_id", to="Results")

        if self.is_bookmark(bookmark):
            query.IN_BOOKMARK(bookmark)
        elif self._is_statustable(bookmark):
            query.JOIN(bookmark, bookmark, "pose_id")

        return query.build()[0]

    # endregion

    # region virtual public api
    def close_storage(self, attached_db=None, vacuum=False):
        """Close connection to database

        Args:
            attached_db (str, optional): name of attached DB (not including file extension)
            vacuum (bool, optional): indicates that database should be vacuumed before closing
        """
        raise NotImplementedError

    def check_ringtaildb_version(self):
        """
        Checks the database version and confirms whether the code base is compatible with it

        Returns:
            bool: whether or not db is compatible with the code base
            str: current database version
        """
        raise NotImplementedError

    def insert_receptor_blob(self, receptor: bytes, rec_name: str):
        """Takes object of Receptor class, updates the column in Receptor table

        Args:
            receptor (bytes): bytes receptor object to be inserted into DB
            rec_name (string): Name of receptor. Used to insert into correct row of DB
        """
        raise NotImplementedError

    def insert_receptor_polymer(self, receptor: str, rec_name: str):
        """Takes object of Receptor class, updates the column in Receptor table

        Args:
            receptor (str): json string representation of a receptor meeko.Polymer oobject to be inserted into DB
            rec_name (str): Name of receptor. Used to insert into correct row of DB
        """
        raise NotImplementedError

    def clone(self, backup_name: str = None) -> str:
        """Creates a copy of the db

        Args:
            backup_name (str, optional): name of the cloned database

        Returns:
            str: path of backed up database
        """
        raise NotImplementedError

    def update_database_version(self, new_version: str, consent=False, backup=False):
        """method that updates sqlite database schema 1.0.0 or 1.1.0 to 1.1.0 or 2.0.0

        #NOTE: If you created the database with the duplicate handling option,
        # there is a chance of inconsistent behavior of anything involving interactions as
        the Pose_ID was not used as an explicit foreign key in db v1.0.0 and v1.1.0.

        Args:
            new_version (str): _description_
            consent (bool, optional): variable to ensure consent to update database is explicit
            backup (bool, optional): _description_. Defaults to False.

        Returns:
            bool: final consent
        """
        raise NotImplementedError

    def db_query(self, query: str, params: iter) -> iter:
        """Executes provided sql query. Returns iter for results.

        Args:
            query (str): Formated sql query as string
            params (iter): parameters to be used in query (assumes query as appropriate place holders)

        Returns:
            iter: Contains results of query
        """
        raise NotImplementedError

    def db_update(self, query: str, parameters: list[tuple], commit=True) -> iter:
        """
        A db query that uses executemany

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
        raise NotImplementedError

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
        raise NotImplementedError

    def table_length(self, table: str) -> int:
        """
        Get length of table or bookmark

        Args:
            table (str): name of table or bookmark

        Returns:
            int: number of poses in table/bookmark
        """
        raise NotImplementedError

    def fetch_clustered_similars(self, ligname: str):
        """Given ligname, returns poseids for similar poses/ligands from previous clustering. User prompted at runtime to choose cluster.

        Args:
            ligname (str): ligname for ligand to find similarity with

        Raises:
            ValueError: wrong terminal input
        """
        raise NotImplementedError

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
        raise NotImplementedError

    def _create_status_tables(self):
        """
        Creates status tables if needed
        """
        raise NotImplementedError

    def pose_row_in_table(self, table: str, pose_id: int) -> Union[None, int]:
        """
        Find the row id of a pose in a given table

        Args:
            table (str)
            pose_id (int)

        Returns:
            Union[None, int]: rowid if any
        """
        raise NotImplementedError

    def tables_in_db(self) -> list:
        """
        Returns a list of all table names in the database

        Returns:
            list: list of table names
        """
        raise NotImplementedError

    def get_starting_rowid(self, table: str) -> int:
        """
        Starting row id for a table, will be 1 for regular tables, and 1 or non-1 for bookmarks
        (whose rows are inside Filtered_poses)

        Args:
            table (str): table or bookmark name

        Returns:
            int: first row id belonging to that selection
        """
        raise NotImplementedError

    # endregion

    # region private api

    def _create_tables(self, commit=True) -> None:
        """
        Creates all tables needed for a Ringtail database
        """
        self._create_ligands_table()
        self._create_results_table()
        self._create_receptors_table()
        self._create_interaction_index_table()
        self._create_interaction_table()
        self._create_db_properties_table()
        self._create_filtering_tables()
        self._create_status_tables()

        self._set_ringtail_db_schema_version(self._db_schema_ver)
        if commit:
            self.conn.commit()

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
            ligname [business key]
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
        # append ligand name as business key to Ligands table
        ligand_data_list.append(ligand_dict["ligname"])

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
        """Takes dictionary of file results, formats as list of tuples for interactions.
        To each interaction description is added business keys/columns that identifies
        which results row/pose each interaction belongs to

        Args:
            interaction_dictionaries (list): List of pose interaction
                dictionaries from parser

        Returns:
            list: of tuples of interaction data
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
                # include the three values that uniquely identify a pose within a docking run
                interactions.add(
                    (
                        pose_interactions["ligand_name"],
                        pose_interactions["run_number"],
                        pose_interactions["pose_rank"],
                    )
                    + tuple(pose_interactions[kw][i] for kw in interaction_keywords)
                )

        return list(interactions)

    def _insert_docking_data(
        self, results: list[list], interactions: list[list], options: dict
    ):
        """Takes list of database rows to insert, adds data to results table.
        First stages data in temporary tables, then handles duplicates (if requested),
        and finally transfers data from temporary tables into permanent storage and
        commits once at the end.

        Args:
            results (list): list of arrays containing formatted result rows
            interactions (list): list of interactions
            options (dict): includes options on how to handle duplicates if there are any

        """

        self._create_temporary_results_tables()
        self._insert_results_in_temp_tables(results, interactions)
        dupl_handl = options.get("duplicate_handling")
        # handle duplicates if requested
        if dupl_handl and dupl_handl.lower() == "replace":
            # first delete duplicate results and interactions
            # then insert all the new ones indiscrimenately
            self._delete_old_duplicate_results()
        elif dupl_handl and dupl_handl.lower() == "ignore":
            # delete from incoming data any duplicates
            # then insert all new data indiscrimenately
            self._delete_new_duplicate_results()
        # then, move results from temp tables to database
        self._move_tempresults_to_database()
        logger.info(
            f"Results ({len(results)} rows) and interactions ({len(interactions)} rows) have been added to the database"
        )

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
            poseid_list (list): representative poses for each cluster
            cluster_type (str): how clustering was performed
            cluster_cutoff (str): distance to representative pose
            bookmark_name (str): bookmark name which is clustered over

        Returns:
            str: name of cluster bookmark
        """
        self._create_cluster_tables()

        cluster_name = f"{cluster_type}_{str(cluster_cutoff)}"

        cluster_bookmark = f"{bookmark_name}_{cluster_name}"
        # check if clusters already exist
        clusters_exist = self._cluster_exists(cluster_name, bookmark_name)
        if clusters_exist:
            logger.warning(
                f"This cluster has been ran before, will reuse bookmark {cluster_bookmark}."
            )
            return (cluster_bookmark, clusters_exist[0])

        cluster_id = self._insert_new_cluster_info(
            cluster_name, "", bookmark_name, len(clusters)
        )
        print(
            f"Length of clusters coming in: {len(clusters)}, and length of representative\n poses: {len(poseid_list)}"
        )
        cluster_groups = []
        pose_rows = []
        for group_index, cluster in enumerate(clusters):
            cluster = list(cluster)
            representative_pose = poseid_list[group_index]
            cluster_groups.append([cluster_id, group_index, representative_pose])
            # make sure we add the representative pose
            cluster.append(representative_pose)
            for pose in cluster:
                pose_rows.append([cluster_id, group_index, pose])

        self._insert_clusters(cluster_groups, pose_rows)
        self.conn.commit()

        return cluster_bookmark

    def _open_storage(self):
        """Create connection to db. Then, check if db needs to be created.

        Raises:
            StorageError
            VersionError
        """
        try:
            # check to see if file exist, and if it does, check that version is matching
            if os.path.isfile(self.db_file) and os.path.getsize(self.db_file) > 0:
                self.conn = self._create_connection()
                compatible, db_version = self.check_ringtaildb_version()
                if not compatible:
                    raise VersionError(
                        f"The database is of version {db_version} which is not compatible with the code base of version {version('ringtail')}"
                    )
            else:
                logger.info("Creating a new database file.")
                self.conn = self._create_connection()
            if self.keyboard_interrupt_allowed:
                signal(
                    SIGINT, self._sigint_handler
                )  # signal handler to catch keyboard interupts
            logger.debug(
                f"Ringtail connected to database {self.db_file} with connection: {self.conn}"
            )
        except VersionError as e:
            raise
        except Exception as e:
            raise StorageError(f"Error while creating or connecting to database: {e}.")

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

    def _drop_existing_tables(self):
        """Drops existing tables, in order of foreign key dependency"""
        # first, delete tables with foreign keys
        self._delete_table("merged_tables")
        self._delete_table("PK_conversions")
        self._delete_table("Pose_clusters")
        self._delete_table("Cluster_groups")
        self._delete_table("Clusters")
        self._delete_table("Filtered_poses")
        self._delete_table("Filters")
        self._delete_table("Interactions")
        self._delete_table("Interaction_indices")
        for table in statuses:
            self._delete_table(table.capitalize())
        self._delete_table("Results")
        # then, fetch remaining tables
        tables = self.tables_in_db()
        for table in tables:
            self._delete_table(table)
        # check if has sequences
        self._delete_nontables()
        self.conn.commit()

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
        query = self.QueryBuilder()
        query.DROP_IF_EXISTS(name)
        return self.db_query(query.build()[0])

    def _is_statustable(self, table: str) -> bool:
        """
        Returns True if table name is actually a status table (table with poses who have been assigned a status like accept, reject, maybe)

        Args:
            table (str): name of table or bookmark to check

        Returns:
            bool: if table name is a status table
        """
        if table.lower() in statuses:
            return True
        else:
            return False

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

    def _process_filters_for_query(self, filters_dict: dict) -> dict:
        """
        Method that reformats the filters to the specified database columns, handles less than/more than filters, etc

        Args:
            filters_dict (dict): all Ringtail filters, okay to contain None

        Returns:
            dict: of lists of numerical, interaction, and ligand filters + maxmiss
        """
        score_maxmin_to_sql = {
            "eworst": "docking_score <= {value}",
            "ebest": "docking_score >= {value}",
            "leworst": "leff <= {value}",
            "lebest": "leff >= {value}",
        }

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
                        score_maxmin_to_sql[filter_key].format(value=filter_value)
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

        select_statement = "SELECT R.pose_id, Ligands.ligand_id, Ligands.rdmol "
        headers = [
            "pose_id",
            "ligand_id",
            "rdmol",
        ]
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
            headers.extend(["atom_index_map", "ligand_coordinates"])
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
            ligandict = dict(zip(headers, ligandrow))
            # substruct and maxatoms do not discriminate on poses, check if ligand has already been accounted for
            if not position and ligandict["ligand_id"] in filtered_ligands:
                filtered_ligands[ligandict["ligand_id"]].append(ligandict["pose_id"])
            # the real workhorse
            else:
                # deserialize binary rdmol
                ligand_mol = Chem.Mol(ligandict["rdmol"])
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
                        ligand_index_map = _ligand_indexmap(ligandict["atom_index_map"])
                        ligand_coordinates = ligandict["ligand_coordinates"]
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
            if ligandict["ligand_id"] in filtered_ligands:

                filtered_ligands[ligandict["ligand_id"]].append(ligandict["pose_id"])
            # add new ligand in the list of passing ligands
            else:
                filtered_ligands[ligandict["ligand_id"]] = [ligandict["pose_id"]]

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
        query = self.QueryBuilder()
        query.SELECT("Pose_ID", "interaction_id").FROM("Interactions").WHERE(
            f"""pose_id IN ({",".join(["?"] * len(pose_ids))})""", *pose_ids
        )

        poseid_intinds = self.db_query(*query.build()).fetchall()
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

    def _get_interaction_indices(self, interaction_list: list) -> iter:
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

        where_clause = " AND ".join(
            [
                f"""{interaction_info[i]} = '{interaction_list[i]}'"""
                for i in range(len_interaction_info)
                if interaction_list[i] != ""
            ]
        )
        query = self.QueryBuilder()
        query.SELECT("interaction_id").FROM("Interaction_indices").WHERE(where_clause)

        return self.db_query(*query.build()).fetchall()

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
            int: number of passing poses
        """
        count_query = f"""SELECT COUNT(pose_id) FROM ({query});"""
        num_passing_poses = self.conn.execute(count_query).fetchone()[0]

        if not num_passing_poses:
            return 0

        # make sure bookmark name is not a table name
        if name in self.tables_in_db():
            raise OptionError(
                f"Bookmark name {name} is the same as an existing table in the database, and cannot be used."
            )
        # overwrite bookmark if it exists
        if self.is_bookmark(name):
            logger.warning(
                f"The bookmark {name} already exists, and will be overwritten by the current filter."
            )
            self.delete_bookmark(name)

        self._begin_transaction()
        insert_query = self.QueryBuilder()
        insert_query.INSERT_INTO(
            "Filters", "name", "query", "filters", "filter_window"
        ).RETURNING("filter_id")

        # insert filter info, return filter_id
        filter_id = self.conn.execute(
            insert_query.build()[0],
            (name.lower(), query, json.dumps(filters), filtering_bookmark),
        ).fetchone()[0]

        # insert filtered poses with filter_id
        insert_query = f"""
            WITH results_poses AS (
                {query})
            INSERT INTO Filtered_poses(filter_id, pose_id) 
            SELECT {filter_id}, pose_id FROM results_poses;
            """
        self.conn.execute(insert_query)
        self.conn.commit()

        return num_passing_poses

    def get_useful_columns(self):
        return {
            "Results": [
                "docking_score",
                "pose_rank",
            ],
            "Ligands": ["ligname", "smiles"],
        }

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
            columns_info = self._fetch_table_column_names(table)
            columns.extend(columns_info)
            columns_with_tablename.extend(
                [
                    ("{{{table_alias}_alias}}.{col}".format(col=col, table_alias=table))
                    for col in columns_info
                ]
            )

        return columns, columns_with_tablename

    # endregion

    # region virtual private methods
    def _create_connection(self):
        """Creates database connection to self.db_file

        Returns:
            <db type>conn: Connection object to self.db_file

        Raises:
            DatabaseConnectionError
        """
        raise NotImplementedError

    def _begin_transaction(self):
        raise NotImplementedError

    def _rollback(self):
        raise NotImplementedError

    def _create_results_table(self, name="Results"):
        """Creates table for results. Columns are:
        Pose_ID             INTEGER PRIMARY KEY AUTOINCREMENT,
        ligand_id           INTEGER FOREIGN KEY from Ligands,
        receptor            VARCHAR,
        pose_rank           INTEGER,
        run_number          INTEGER,
        docking_score       FLOAT,
        leff                FLOAT,
        deltas              FLOAT,
        cluster_rmsd        FLOAT,
        cluster_size        INTEGER,
        reference_rmsd      FLOAT,
        energies_inter      FLOAT,
        energies_vdw        FLOAT,
        energies_electro    FLOAT,
        energies_flexLig    FLOAT,
        energies_flexLR     FLOAT,
        energies_intra      FLOAT,
        energies_torsional  FLOAT,
        unbound_energy      FLOAT,
        nr_interactions     INTEGER,
        num_hb              INTEGER,
        about_x             FLOAT,
        about_y             FLOAT,
        about_z             FLOAT,
        trans_x             FLOAT,
        trans_y             FLOAT,
        trans_z             FLOAT,
        axisangle_x         FLOAT,
        axisangle_y         FLOAT,
        axisangle_z         FLOAT,
        axisangle_w         FLOAT,
        dihedrals           VARCHAR,
        ligand_coordinates         VARCHAR,
        flexible_res_coordinates   VARCHAR

        """
        raise NotImplementedError

    def _create_ligands_table(self, name="Ligands") -> None:
        """Create table for ligands. Columns are:
        ligand_id           INTEGER PRIMARY KEY AUTOINCREMENT,
        LigName             VARCHAR NOT NULL UNIQUE ON CONFLICT IGNORE,
        ligand_smile        VARCHAR,
        rdmol               BLOB,
        atom_index_map      VARCHAR,
        hydrogen_parents    VARCHAR,
        input_model         VARCHAR
        """
        raise NotImplementedError

    def _create_receptors_table(self):
        """Create table for receptors. Has primary key although only one receptor allowed"""
        raise NotImplementedError

    def _create_interaction_index_table(self):
        """Creates a table describing unique interactions in the database"""
        raise NotImplementedError

    def _create_interaction_table(self):
        """Creates a table of each pose-interaction combination."""
        raise NotImplementedError

    def _create_db_properties_table(self):
        """Create table of database properties used during write session to the database. Columns are:
        DB_write_session int (primary key)
        docking_mode (vina or adgpu)
        num_of_poses ("all" or str(int))
        """
        raise NotImplementedError

    def _insert_db_properties(self, docking_mode: str, number_of_poses: str):
        """Insert db properties into database properties table

        Args:
            docking_mode (str): docking mode for the current dataset being written
            number_of_poses (str): number of poses written to database in current session, either "all" or specified max_poses
        """
        raise NotImplementedError

    def _create_temporary_results_tables(self):
        """
        Creates temporary tables for results and interactions, which will be
        used for staging incoming data.
        """
        raise NotImplementedError

    def _insert_results_in_temp_tables(
        self, results_array: list, interactions_array: list
    ):
        """
        Inserts docking results and interactions into their respective
        temporary tables

        Args:
            results_array (list): list of result rows
            interactions_array (list): list of interaction rows
        """
        raise NotImplementedError

    def _move_tempresults_to_database(self):
        """Inserts data from the temporary results tables to their permanent
        database equivalents"""
        raise NotImplementedError

    def _delete_new_duplicate_results(self):
        """Checks if a pose is uniquely represented in the Results table,
        and deletes it from the staged incoming data if duplicated.
        Based on the following columns:
        ligname,
        receptor,
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
        dihedrals
        """
        raise NotImplementedError

    def _delete_old_duplicate_results(self):
        """Checks if a pose is uniquely represented in the Results table,
        and deletes it from Results if duplicated.
        Based on the following columns:
        ligname,
        receptor,
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
        dihedrals
        """
        raise NotImplementedError

    def _create_filtering_tables(self):
        """
        Creates a Filter table which includes filter_id (PK), name (bookmark_name), duckdb formatted query,
        and dictionary of filters used, as well as Filtered_poses, which uses filter_id as FK,
        and lists all poses passing that filter_id
        """
        raise NotImplementedError

    def _set_ringtail_db_schema_version(self, db_version: str = "3.0.0"):
        """
        Will check current storage manager db schema version and only set if it
        is compatible with the code base version (i.e., version(ringtail)).

        Args:
            db_version (str, optional): _description_. Defaults to "3.0.0".

        Raises:
            StorageError: _description_
        """
        pass

    def _insert_ligands(self, ligands: list):
        """Takes array of ligand rows, inserts into Ligands table.

        Args:
            ligand_array (list[list]): list of lists containing formatted ligand rows

        """
        raise NotImplementedError

    def _create_indices(self):
        """
        Creates specific indices on tables for those databases that use indices
        """
        pass

    def _insert_receptors(self, receptor_array: list):
        """Takes array of receptor rows, inserts into Receptors table

        Args:
            receptor_array (list): List of lists
                containing formatted receptor rows
        """
        raise NotImplementedError

    def _insert_interaction_index_rows(self, interactions: list[tuple]):
        """
        Writes unique interactions to database

        Args:
            interaction_tuple (list[tuple]): [(interaction_type, rec_chain, rec_resname, rec_resid, rec_atom, rec_atomid)]
        """
        raise NotImplementedError

    def _delete_nontables(self):
        """
        Deletes objects in the database that are not tables
        """
        pass

    def _generate_result_filtering_query(
        self, filters_dict, bookmark_name, filter_bookmark
    ):
        raise NotImplementedError

    def _create_cluster_tables(self):
        """
        Creates cluster tables if they don't already exist
        """
        raise NotImplementedError

    def _cluster_exists(
        self, cluster_name: str, cluster_window: str
    ) -> Union[int, None]:
        """
        Checks if a cluster already exists, based on what window was clustered over,
        and the standardized cluster name. Method will not work if cluster name starts
        to be non-standardized

        Args:
            cluster_name (str):
            cluster_window (str):

        Returns:
            int: number of clusters in that cluster if any, else None
        """
        raise NotImplementedError

    def _insert_new_cluster_info(
        self, name: str, description: str, cluster_window: str, length: int
    ) -> int:
        """
        Inserts the basic info about a clustering exercise, but not the clustered data itself

        Args:
            name (str): name of cluster (standardized)
            description (str): placeholder for if/when more info is needed about clusters
            cluster_window (str): what was clustered over, bookmark or all results
            length (int): number of clusters

        Returns:
            int: cluster id of the new inserted cluster
        """
        raise NotImplementedError

    def _insert_clusters(self, cluster_groups: list, pose_rows: list):
        """
        Inserts all cluster data, including each grouping and its representative
        pose, and all poses involved and which clusters they belong to

        Args:
            cluster_groups (list): each cluster from the clustering exercise
            pose_rows (list): pose and cluster id and group id
        """
        raise NotImplementedError

    def _get_numeric_columns(self, table_name: str) -> list:
        """
        Method to get the names of all numeric columns in a table, for example for
        allowable sorting options

        Args:
            table_name (str): table name to evaluate

        Returns:
            list: column names that has a numeric type
        """
        raise NotImplementedError

    def _check_if_column_in_table(self, table_name: str, column_name: str) -> bool:
        """
        Checks if a column exist in given table

        Args:
            table_name (str):
            column_name (str):

        Returns:
            bool: True if column there, False if not
        """
        raise NotImplementedError

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
        raise NotImplementedError

    def _detach_db(self, new_db_alias):
        """Detaches new database file from current database

        Args:
            new_db_name (str): db name for database to detach

        Raises:
            StorageError
        """
        raise NotImplementedError

    def _attach_db(self, new_db: str, new_db_alias: str = "attached_db") -> str:
        """Attaches new database file to current database

        Args:
            new_db (str): file name for database to attach
            new_db_name (str): name of new database

        """
        raise NotImplementedError

    def _create_crossref_temp_table(self, table_name: str):
        """create temporary table with given name and with ligand name and pose_id information

        Args:
            table_name (str): name for temp table
        """
        raise NotImplementedError

    # endregion

    # region cross referencing filtered databases

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
        self.db_query(temp_insert_query, commit=True)

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
        return self.db_query(counting.build()[0]).fetchone()[0]

    # endregion
