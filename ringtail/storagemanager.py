#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail storage adaptors
#

import time
import json
import os.path
from abc import ABC, abstractmethod
from pathlib import Path
import pandas as pd
from .logutils import get_logger

logger = get_logger(__name__)
from .util import db_alias_from_path, numlist2str
import sys
from signal import signal, SIGINT
from rdkit import Chem
from typing import Any, Union, ClassVar
from importlib.metadata import version
from .ringtailoptions import Filters, statuses
from .exceptions import StorageError, OptionError, VersionError, MergeError
from .clustermanager import *
from .querybuilder import QueryBuilder
from .schema import (
    NUMERIC_TYPES,
    OUTFIELD_BY_TABLE,
    CANDIDATES_SUBQ,
    CANDIDATES_NAME,
    LIGANDS_ONLY_COLS,
    TABLE_SCHEMAS,
    LIGANDS_SCHEMA,
    RESULTS_SCHEMA,
    DB_PROPERTIES_SCHEMA,
    INTERACTION_INDICES_SCHEMA,
    INTERACTIONS_SCHEMA,
    FILTERS_SCHEMA,
    FILTERED_POSES_SCHEMA,
    CLUSTERS_SCHEMA,
    CLUSTER_GROUPS_SCHEMA,
    POSE_CLUSTERS_SCHEMA,
    MERGED_TABLES_SCHEMA,
    PK_CONVERSIONS_SCHEMA,
    RECEPTORS_SCHEMA,
    STATUS_TABLE_SCHEMA,
    POSE_COMMENTS_SCHEMA,
    build_create_table,
)
from collections import defaultdict


class StorageManager(ABC):
    _db_schema_ver = "3.0.0"
    QueryBuilder = QueryBuilder

    # Whether Filtered_poses carries a denormalized ligand_id column. SQLite sets
    # this True so passing-ligand counts are a single-table COUNT(DISTINCT) instead
    # of a costly JOIN into the row-store Results table. DuckDB leaves it False
    # (its columnar count-JOIN is cheap) so existing DuckDB databases stay valid.
    _filtered_poses_has_ligand_id: ClassVar[bool] = False

    """Base class for a generic virtual screening database object.
    This class holds some of the common API for StorageManager child classes. 
    Each child class will implement their own functions to write to and read from the database

    Attributes: 
        _db_schema_ver (str): current database schema version
        _db_schema_code_compatibility (dict): dictionary showing compatibility of code base versions with relational database schema versions
        active_connection (bool): if there is a current open connection to database
        keyboard_interrupt_allowed (bool): if Ctrl+Z will work, for example not supported in a GUI
    """

    dialect: ClassVar[str]

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
            raise
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
            raise exc_value
        return self

    def _sigint_handler(self, signal_received, frame):
        """Handles and reports if program is interrupted through the terminal"""
        logger.critical("Ctrl + C pressed, keyboard interupt initiated")
        self.__exit__(None, None, None)
        sys.exit(0)

    # endregion

    # region public api

    def insert_data(self, docking_data: dict, duplicate_handling: str = None):
        """Inserts data from all arrays returned from results manager. Can have one or more ligands in docking_data

        Args:
            docking_data (dict): docking results to be inserted, where key is ligand name and value is data to be written
            duplicate_handling (str, optional): Describes how to handle duplicates, if any. Defaults to None.
        """
        ligands = docking_data.get("ligands", None)
        if ligands:
            self._insert_ligands(ligands)
        interaction_data = docking_data.get("interactions", None)
        if interaction_data:
            self._insert_interaction_index_rows(interaction_data)
        poses = docking_data.get("poses", None)
        if poses:
            # interactions has ligname and pose rank
            self._insert_docking_data(poses, interaction_data, duplicate_handling)
        self.conn.commit()

    def post_insert_interactions(
        self,
        interactions: list[tuple],
        interaction_counts: list[dict],
        processed_pose_ids: list[tuple],
        tracking_table_name: str,
    ):
        """
        Method to insert interactions into an already populated database, for example when
        recalculating interactions with new interaction distances

        Args:
            interactions (list[tuple]): list of interaction tuples
            interaction_counts (list[dict]): interaction counts for the Results table
            processed_pose_ids (list[tuple]): psoe_ids whose interactions have been (re-)calculated
            tracking_table_name (str): table that tracks (re-)processed pose interactions
        """
        self._begin_transaction()

        self._insert_interaction_index_rows(interactions)
        self._insert_interactions(interactions)
        self._update_interaction_counts(interaction_counts)
        self._insert_completed_poses(processed_pose_ids, tracking_table_name)

        self.conn.commit()

    def insert_receptor_basic_info(self, receptor_data: list):
        """

        Method to insert receptor information into the database

        Args:
            receptor_data (list): of receptor data, ordered by columns in the db
        """
        receptors = self.fetch_receptor_object()
        # insert receptor if database does not have already have a receptor entry
        if not receptors.get("recname"):
            self._insert_receptors(receptor_data)

    def filter_results(
        self,
        all_filters: Filters,
        output_bookmark: str,
        input_bookmark: str = None,
    ) -> int:
        """Generate and execute database queries from given filters.

        Args:
            all_filters (dict): dict containing all filters
            output_bookmark (str): name of bookmark in which to save the data
            input_bookmark (str, optional): if filtering not across all data, but a pre-filtered bookmark

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
            all_filters, output_bookmark, input_bookmark
        )

        logger.debug(f"Query for filtering results: {filter_query}")

        # perform filtering
        logger.debug("Running filtering query...")
        time0 = time.perf_counter()
        if input_bookmark is None:
            input_bookmark = "Results"

        self._populate_filter_tables(
            name=output_bookmark,
            query=filter_query,
            filters=all_filters,
            input_bookmark=input_bookmark,
        )
        logger.debug(
            f"Time to filter results: {time.perf_counter() - time0:.2f} seconds"
        )
        count = self.get_passing_ligands_count(output_bookmark)

        return count

    def get_all_bookmark_names(
        self, db_name: str = None, db_path: str = None
    ) -> list[str]:
        """Get all bookmarks in sql database as a list of names.

        Args:
            db_name (str, optional): if attaching a database. Defaults to None.
            db_path (str, optional): path to attached database. Defaults to None.

        Returns:
            list: of bookmark names
        """
        if db_path:
            self._attach_db(db_path, db_name)
        query = self.QueryBuilder()
        query.SELECT("name").FROM("Filters", db_name=db_name)
        cur = self.db_query(query.build()[0])
        bookmark_names = [row[0].lower() for row in cur.fetchall()]

        if db_path:
            self.detach_db(db_name)

        return bookmark_names

    def get_passing_poses_count(
        self, bookmark_name: str, grouped_by_ligand: bool = False
    ) -> int:
        """
        Count poses (or distinct ligands) in a bookmark.

        Args:
            bookmark_name (str): bookmark name in which to count
            grouped_by_ligand (bool, optional): if True, returns number of distinct ligands

        Returns:
            int: number of poses or ligands in the bookmark
        """
        filter_row = self.conn.execute(
            "SELECT filter_id FROM Filters WHERE name = ?", (bookmark_name.lower(),)
        ).fetchone()
        if filter_row is None:
            return 0
        filter_id = filter_row[0]
        if grouped_by_ligand:
            if self._filtered_poses_has_ligand_id:
                # SQLite: ligand_id is denormalized into Filtered_poses — single-table count
                return self.conn.execute(
                    "SELECT COUNT(DISTINCT ligand_id) FROM Filtered_poses WHERE filter_id = ?",
                    (filter_id,),
                ).fetchone()[0]
            # DuckDB: columnar JOIN to Results reads only the ligand_id column — cheap
            return self.conn.execute(
                "SELECT COUNT(DISTINCT R.ligand_id) FROM Filtered_poses FP "
                "JOIN Results R ON R.pose_id = FP.pose_id WHERE FP.filter_id = ?",
                (filter_id,),
            ).fetchone()[0]
        return self.conn.execute(
            "SELECT COUNT(*) FROM Filtered_poses WHERE filter_id = ?",
            (filter_id,),
        ).fetchone()[0]

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

    def rename_bookmark(self, old_name: str, new_name: str) -> None:
        """Rename a bookmark. Clusters.cluster_window references remain valid because
        intermediate bookmarks are never renamed — only the final cluster output is.

        Args:
            old_name (str): current bookmark name
            new_name (str): desired bookmark name
        """
        if self.is_bookmark(new_name):
            logger.warning(
                f"The bookmark {new_name} already exists, and will be overwritten."
            )
            self.delete_bookmark(new_name)
        self.conn.execute(
            "UPDATE Filters SET name = ? WHERE name = ?", (new_name, old_name)
        )
        self.conn.commit()
        logger.info(f"Bookmark '{old_name}' renamed to '{new_name}'.")

    def get_maxmiss_union(
        self, total_combinations: int, bookmark_name: str, all_filters=None
    ) -> tuple[int, str]:
        """
        Get results that are in union considering max miss

        Args:
            total_combinations (int): numer of possible combinations
            bookmark_name (str): name of bookmark to store
            all_filters (dict, optional): All filters used. Defaults to None.

        Returns:
            int: num passing poses
            str: union bookmark name
        """
        all_filters = all_filters or {}
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
        union_cols = (
            "DISTINCT pose_id, ligand_id"
            if self._filtered_poses_has_ligand_id
            else "DISTINCT pose_id"
        )
        query = self.QueryBuilder()
        query.SELECT(union_cols).FROM("filtered_poses").WHERE(
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

    def crossref_databases(
        self,
        wanted_dbs: list[tuple[str, Union[str, None]]] = None,
        unwanted_dbs: list[tuple[str, Union[str, None]]] = None,
        bookmark_prefix: str = "crossref",
        store_best_pose: bool = False,
        alternative_database_names: dict = None,
    ) -> tuple[int, dict, dict]:
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
            store_best_pose (bool, optional): To store just the best pose for a ligand vs all poses in the database. Defaults to False/storing all poses.
            alternative_database_names (dict, optional):  {path: alt name}. Defaults to None.

        Returns:
            int: number of passing crossref ligands
            dict: {db_file: new_bookmark_name,}
            dict: "filter" {"wanted":[list of wanted dbs],"unwanted":[list of unwanted dbs]
        """
        alternative_database_names = alternative_database_names or {}
        processed_wanted = []
        for index, database in enumerate(wanted_dbs):
            if index == 0:
                db_name = "main"
                alternative_database_names[database[0]] = db_name
            else:
                if database[0] in alternative_database_names.keys():
                    db_name = alternative_database_names[database[0]]
                else:
                    db_name = db_alias_from_path(database[0])
                    alternative_database_names[database[0]] = db_name
            processed_wanted.append(
                (
                    db_name,
                    database[0],
                    database[1],
                )
            )
            if index > 0:
                self._attach_db(database[0], db_name)

        processed_unwanted = []
        for index, database in enumerate(unwanted_dbs or []):
            if database[0] in alternative_database_names.keys():
                db_name = alternative_database_names[database[0]]
            else:
                db_name = db_alias_from_path(database[0])
                alternative_database_names[database[0]] = db_name
            processed_unwanted.append(
                (
                    db_name,
                    database[0],
                    database[1],
                )
            )
            self._attach_db(database[0], db_name)
        # make subqueries/common table expressions for each set of ligands
        selection_query = """
        {db_name} AS 
            (SELECT ligname FROM {db_name}.Ligands 
            INNER JOIN {db_name}.Results ON 
                {db_name}.Results.ligand_id={db_name}.Ligands.ligand_id 
                WHERE pose_id IN ({bookmark_poses}))"""
        wanted_ctes = []
        for db_name, _, bookmark in processed_wanted:
            wanted_ctes.append(
                selection_query.format(
                    db_name=db_name,
                    bookmark_poses=self._get_bookmark_poses_query(
                        bookmark.lower(), db_name
                    ),
                )
            )

        unwanted_ctes = []
        for db_name, _, bookmark in processed_unwanted:
            unwanted_ctes.append(
                selection_query.format(
                    db_name=db_name,
                    bookmark_poses=self._get_bookmark_poses_query(
                        bookmark.lower(), db_name
                    ),
                )
            )

        all_ctes = "WITH\n" + ", ".join(wanted_ctes + unwanted_ctes)

        wanted_joins = " ".join(
            f"INNER JOIN {db[0]} ON {processed_wanted[0][0]}.ligname = {db[0]}.ligname"
            for db in processed_wanted[1:]
        )

        # build exclusion union
        exclude_union = " UNION ".join(
            f"SELECT ligname FROM {db[0]}" for db in processed_unwanted
        )

        if unwanted_dbs:
            unwanted_ctes = f""",
                unwanted AS (
                    {exclude_union}
                )"""
            unwanted_where = """WHERE i.ligname NOT IN (SELECT ligname FROM unwanted)"""
        else:
            unwanted_ctes = ""
            unwanted_where = ""

        query = f"""
        {all_ctes},
        wanted AS (
            SELECT {processed_wanted[0][0]}.ligname
            FROM {processed_wanted[0][0]} {wanted_joins}
        ){unwanted_ctes}
        SELECT DISTINCT i.ligname
        FROM wanted i
        {unwanted_where}
        """

        # make a bookmark in each database with the results
        approved_ligand_names = [lig[0] for lig in self.db_query(query).fetchall()]
        dbs_new_bookmark_names = {}

        # don't proceed if no approved ligands
        if len(approved_ligand_names) < 1:
            return 0, {}, {}

        # write bookmarks in each of the databases
        filter_dict = {"wanted": wanted_dbs, "unwanted": unwanted_dbs}

        for db_alias, path, bookmark in processed_wanted + processed_unwanted:
            if alternative_database_names[path].lower() in ["", "main"]:
                self._populate_filter_tables(
                    name=f"{bookmark_prefix}_{bookmark}",
                    query=self._crossref_bookmark_builder(
                        approved_ligand_names, store_best_pose
                    ),
                    filters=filter_dict,
                    input_bookmark=bookmark,
                )
            else:
                # detach
                self.detach_db(db_alias)
                # make storageman object
                with type(self)(path) as db:
                    db._populate_filter_tables(
                        name=f"{bookmark_prefix}_{bookmark}",
                        query=self._crossref_bookmark_builder(
                            approved_ligand_names, store_best_pose
                        ),
                        filters=filter_dict,
                        input_bookmark=bookmark,
                    )

        for _, file, bookmark in processed_wanted + processed_unwanted:
            dbs_new_bookmark_names.update({file: f"{bookmark_prefix}_{bookmark}"})

        return (len(approved_ligand_names), dbs_new_bookmark_names, filter_dict)

    def cluster_data(
        self,
        bookmark_name: Union[str, None],
        cluster_type: str = "mfpt",
        cutoff: float = 0.5,
    ) -> tuple[str, int]:
        """
        Clusters data in a given bookmark. Will create a new bookmark with the format
        <bookmark_name>_<cluster-type>_clustered containing the representative poses
        for the clusters

        Args:
            bookmark_name (str | None): bookmark name with poses to cluster; None clusters all results
            cluster_type (str, optional): type of clustering. Defaults to "mfpt".
            cutoff (float, optional): cutoff cluster distance. Defaults to 0.5.

        Returns:
            tuple (str, int): clustered bookmark name, number of clusters
        """
        logger.debug("Preparing to cluster")
        time0 = time.perf_counter()

        internal_name = bookmark_name if bookmark_name is not None else "results"
        query = self.QueryBuilder()

        if cluster_type.lower() == "ifp":
            query.SELECT("r.pose_id", "r.leff").FROM("Results", "R")

            if bookmark_name is not None:
                self._apply_table_filter(query, bookmark_name)
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

            if bookmark_name is not None:
                self._apply_table_filter(query, bookmark_name)

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
            representatives,
            cluster_type.lower(),
            str(cutoff),
            internal_name,
        )
        if isinstance(cluster_bookmark_name, tuple):
            logger.info("Clustering has been ran before, old bookmark will be used.")
            num_clusters = cluster_bookmark_name[1]
            cluster_bookmark_name = cluster_bookmark_name[0]
        else:

            clustered_cols = (
                ["pose_id", "ligand_id"]
                if self._filtered_poses_has_ligand_id
                else ["pose_id"]
            )
            clustered_poses = self.QueryBuilder()
            clustered_poses.SELECT(*clustered_cols).FROM("results").WHERE(
                f"pose_id IN ({','.join(str(r) for r in representatives)})"
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

    def merge_database(self, merging_db: str):
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

        Raises:
            MergeError
        """
        # attach incoming database and check compatibility
        merging_db_alias = self._attach_db(merging_db, "merging")
        if not self._db_compatible_for_merge(merging_db_alias):
            raise MergeError(
                "Trying to merge two databases of incompatible or too old versions, cannot proceed."
            )

        # add to merging table the absolute path
        mergedb_abspath = str(os.path.abspath(merging_db))
        merge_id = self._get_merge_id(mergedb_abspath)

        # receptor compatibility check
        if self._merging_receptors_compatible():
            logger.info(
                "The two databases have compatible receptors. Merging will proceed."
            )
        else:
            raise MergeError(
                f"The receptors in the merging databases are not the same. \nThese databases cannot be merged."
            )

        # merge tables
        try:
            # merge db_properties table
            self._merge_db_properties_table(merge_id)
            self.conn.commit()
            logger.info("The 'db_properties' table has been merged.")

            # merge Ligands and Results tables
            self._merge_ligands_and_results_tables(merge_id)
            self.conn.commit()
            logger.info("The 'Ligands' and 'Results' tables have been merged.")

            # merge Interactions and Interaction_indices tables
            self._merge_interaction_tables(merge_id)
            self.conn.commit()
            logger.info(
                "The 'Interaction_indices' and 'Interactions' tables have been merged."
            )
        except Exception as e:
            self._rollback_merge(merge_id)
            raise MergeError(f"Merging failed and was rolled back: {e}") from e
        else:
            self._sync_auto_increment_state()
            logger.info(
                f"The database {merging_db} has been successfully merged into {self.db_file}."
            )
        finally:
            self._cleanup_storage(merging_db_alias, vacuum=True, reindex=True)
            logger.info("The final database has been cleaned up, and indices rebuilt.")

    def complete_merging(self):
        """
        Completes the merging process by creating empty filter tables, as well
        as resetting any auto incremented values for safety
        """
        # Rebuild indices dropped in prepare_for_merging()
        self._create_indices()
        self._create_filtering_tables()
        self._sync_auto_increment_state()

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
            WHERE pose_id IN ({self.QueryBuilder.bookmark_query(bookmark_name, "main")}));"""
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
        WHERE main.Results.pose_id IN ({self.QueryBuilder.bookmark_query(bookmark_name, "main")});"""
        self.db_query(poses)
        logger.debug("Results have been copied into the new subset database.")
        # interaction_indices
        interaction_indices = f"""
        INSERT INTO {alias}.Interaction_indices
        SELECT * FROM main.Interaction_indices
        WHERE main.Interaction_indices.interaction_id IN (
            SELECT interaction_id FROM main.Interactions
            WHERE pose_id IN ({self.QueryBuilder.bookmark_query(bookmark_name, "main")}));"""
        self.db_query(interaction_indices)
        # interactions
        interactions = f"""
        INSERT INTO {alias}.Interactions
        SELECT * FROM main.Interactions
        WHERE main.Interactions.pose_id IN ({self.QueryBuilder.bookmark_query(bookmark_name, "main")});"""
        self.db_query(interactions)
        logger.debug("Interactions have been copied into the new subset database.")
        # receptor
        dbprop = f"""
        INSERT INTO {alias}.db_properties
        SELECT * FROM main.db_properties;
        """
        self.db_query(dbprop)
        logger.debug("The db_properties have been copied into the new subset database.")
        try:
            self.conn.commit()
        except Exception as e:
            raise StorageError(f"Problems while creating a subset database: {e}") from e

        self.detach_db(alias)
        logger.info(f"Subset database {database_name} has been successfully created.")

    def fetch_pose_interactions(self, pose_ids) -> Union[list, dict]:
        """
        Fetch interaction parameters for one pose or a batch of poses.

        Args:
            pose_ids (int | list[int]): single pose id or list of pose ids

        Returns:
            list: interaction rows for a single pose (backward-compatible)
            dict[int, list]: {pose_id: [interaction_rows]} for a list input
        """
        batch = isinstance(pose_ids, list)
        if not batch:
            pose_ids = [pose_ids]
        if "interactions" not in self.tables_in_db():
            return {} if batch else None
        placeholders = ",".join(["?"] * len(pose_ids))
        query = self.QueryBuilder()
        query.SELECT(
            "i.pose_id",
            "ii.interaction_type",
            "ii.rec_chain",
            "ii.rec_resname",
            "ii.rec_resid",
            "ii.rec_atom",
            "ii.rec_atomid",
        ).FROM("Interaction_indices", "ii").JOIN(
            "Interactions", "i", "interaction_id"
        ).WHERE(
            f"i.pose_id IN ({placeholders})", *pose_ids
        )
        rows = self.db_query(*query.build()).fetchall()
        if not batch:
            return [
                row[1:] for row in rows
            ]  # strip pose_id column, backward-compatible
        result = {}
        for pose_id, *interaction in rows:
            result.setdefault(pose_id, []).append(interaction)
        return result

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
        if selection == "*":
            raise OptionError(
                "Output fields/columns cannot be 'all'/'*', please select one or more specific columns, or use the default."
            )
        outfields_list = self._format_output_fields(selection, "R", "L")
        # start formatting write query
        query = self.QueryBuilder()
        # select stuff from results where pose id in filter poses join ligands for extra fields
        query.SELECT(*outfields_list).FROM("Results", "R").JOIN(
            "Ligands", "L", "ligand_id"
        ).WHERE(f"R.pose_id IN ({self._get_bookmark_poses_query(bookmark_name)})")
        if group_by:
            query.GROUP_BY("L.ligname")
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
            if self.is_bookmark(requested_data):
                query.SELECT("*").FROM("Results")
                query.IN_BOOKMARK(requested_data)
            elif self._is_statustable(requested_data):
                query.SELECT("*").FROM("Results")
                query.WHERE(f"pose_id in {requested_data}")
            else:
                query.SELECT("*").FROM(requested_data)
            return self._execute_to_df(query.build()[0])
        else:
            return self._execute_to_df(requested_data)

    def get_query_data_as_dicts(self, query: str) -> tuple[list[str], list[dict]]:
        """
        Will return data requested in an SQL-formatted query

        Args:
            query (str): SQL-formatted query string

        Returns:
            tuple[list[str], list[dict]]: list of column names, and list of dicts where each dict is one row,
                                            and column is the key, value is the row-col cell value
        """
        cur = self.db_query(query)
        rows = cur.fetchall()
        column_names = (
            [desc[0].lower() for desc in cur.description] if cur.description else []
        )
        dict_rows = [dict(zip(column_names, row)) for row in rows]
        return column_names, dict_rows

    def reset_screening_tables(self):
        """
        Deletes and recreates (clears) non results tables for filtering, clustering, gui tables
        """
        self._remove_screening_tables()
        self._create_screening_tables()

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
        query.SELECT("docking_mode").FROM("db_properties").ORDER_BY(
            "DB_write_session"
        ).DESC("DB_write_session").LIMIT(1)
        docking_mode = self.db_query(query.build()[0]).fetchone()
        return docking_mode[0].lower() if docking_mode else None

    def check_storage_ready(
        self, run_mode: str, docking_mode: str = None, num_poses: int = None
    ):
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
        if docking_mode:
            query = self.QueryBuilder()
            query.SELECT("COUNT(*)").FROM("db_properties")
            count = self.db_query(query.build()[0]).fetchone()[0]

            compatible = True
            if count < 1:
                logger.info(
                    "Adding results to an existing database that is currently empty of docking results."
                )
            else:
                compatibility_string = "The following database properties do not agree with the properties last used for this database: \n"
                query = self.QueryBuilder()
                query.SELECT("*").FROM("db_properties").ORDER_BY(
                    "DB_write_session"
                ).DESC("DB_write_session").LIMIT(1)
                cur = self.db_query(query.build()[0])

                _, last_docking_mode, last_num_poses = cur.fetchone()
                if docking_mode != last_docking_mode:
                    compatible = False
                    compatibility_string += f"Current docking mode is {docking_mode} but last used docking mode of database is {last_docking_mode}.\n"

                if last_num_poses == -1 != num_poses:
                    compatible = False
                    compatibility_string += f"Current number of poses saved is {num_poses} but database was previously set to 'store_all_poses'.\n"
                elif last_num_poses != num_poses:
                    compatible = False
                    compatibility_string += f"Current number of poses saved is {num_poses} but database was previously set to {last_num_poses}."

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
        elif self._is_candidates_table(table):
            query.JOIN(CANDIDATES_SUBQ, "T", "pose_id")

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

        col_spec = RESULTS_SCHEMA.columns.get(column)
        if col_spec is None or col_spec.sql_type not in NUMERIC_TYPES:
            raise OptionError(
                f"Requested column {column} is not a numeric Results column, percentiles cannot be calculated."
            )
        query = self.QueryBuilder()
        query.SELECT(f"{column}").FROM("Results")
        self._apply_table_filter(query, table)
        query.GROUP_BY("ligand_id")
        values = [val[0] for val in self.db_query(query.build()[0]).fetchall()]

        bins = np.linspace(0, 100, num_bins + 1)
        bin_edges = np.percentile(values, bins)
        return bins, bin_edges

    def fetch_histogram_counts(
        self, column: str, num_bins: int, bookmark_name: str
    ) -> tuple[list[float], list[int]]:
        """SQL-side histogram aggregation. Returns (bin_edges, counts) where
        len(bin_edges) == num_bins + 1 and len(counts) == num_bins.

        Raw SQL is used intentionally to avoid QueryBuilderDuck.GROUP_BY() wrapping
        all SELECT items in ANY_VALUE(), which breaks computed-column aggregates.

        Args:
            column (str): numeric Results column (e.g. "docking_score", "leff")
            num_bins (int): number of histogram bins
            bookmark_name (str): bookmark, status table, or "Results"

        Returns:
            tuple[list[float], list[int]]: bin_edges and per-bin counts
        """
        col_spec = RESULTS_SCHEMA.columns.get(column)
        if col_spec is None or col_spec.sql_type not in NUMERIC_TYPES:
            raise OptionError(f"Column {column} is not a numeric Results column")

        # Build filtered base subquery via QueryBuilder + _apply_table_filter.
        # base_params is always [] — all filter conditions embed values as SQL literals.
        base_q = self.QueryBuilder()
        base_q.SELECT(column).FROM("Results")
        if bookmark_name and bookmark_name.lower() != "results":
            self._apply_table_filter(base_q, bookmark_name)
        base_q.WHERE(f"{column} IS NOT NULL")
        base_sql, base_params = base_q.build()

        row = self.db_query(
            f"SELECT MIN({column}), MAX({column}) FROM ({base_sql}) sub", base_params
        ).fetchone()
        if not row or row[0] is None:
            return [], []
        min_val, max_val = float(row[0]), float(row[1])
        if max_val == min_val:
            total = self.db_query(
                f"SELECT COUNT(*) FROM ({base_sql}) sub", base_params
            ).fetchone()[0]
            return [min_val - 0.5, min_val + 0.5], [total]
        bin_width = (max_val - min_val) / num_bins

        # Raw SQL kept for GROUP BY on computed column — avoids QueryBuilderDuck
        # wrapping bin_idx in ANY_VALUE(), which breaks the aggregation.
        bin_sql = f"""
            SELECT
                CAST(({column} - {min_val}) / {bin_width} AS INTEGER) AS bin_idx,
                COUNT(*) AS cnt
            FROM ({base_sql}) sub
            GROUP BY bin_idx
            ORDER BY bin_idx
        """
        rows = self.db_query(bin_sql, base_params).fetchall()

        counts = [0] * num_bins
        for bin_idx, cnt in rows:
            idx = min(int(bin_idx), num_bins - 1)
            if 0 <= idx < num_bins:
                counts[idx] = cnt
        bin_edges = [min_val + i * bin_width for i in range(num_bins + 1)]
        return bin_edges, counts

    def fetch_data_for_passing_results(
        self, bookmark_name: str, outfields: Union[str, list], order_results: str = None
    ) -> Any:
        """Will return cursor with requested data for outfields for poses that passed filter in bookmark_name

        Args:
            bookmark_name (str): bookmark for which to retrieve data
            outfields (Union[str, list]): columns to include
            order_results (str, optional): if ordering by a column. Defaults to None.

        Returns:
            Any: cursor of data from passing data

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
        query = self.QueryBuilder()
        try:
            row = self.db_query(
                *query.SELECT("recname", "receptor_object", "polymer")
                .FROM("Receptors")
                .build()
            ).fetchone()
            columns = ["recname", "receptor_object", "polymer"]
        except Exception:
            row = self.db_query(
                *query.SELECT("recname", "receptor_object").FROM("Receptors").build()
            ).fetchone()
            columns = ["recname", "receptor_object"]
        return dict(zip(columns, row)) if row else {}

    def fetch_flexres_info(self, receptor: Union[str, int]):
        """fetch flexible residues names and atomname lists

        Args:
            receptor (Union[str, int]): receptor descriptor, either receptor_id or receptor name

        Returns:
            tuple: (flexible_residues, flexres_atomnames)
        """
        if isinstance(receptor, int):
            selection = "receptor_id = ?"
        elif isinstance(receptor, str):
            selection = "recname = ?"
        query = self.QueryBuilder()
        query.SELECT("flexible_residues", "flexres_atomnames").FROM("Receptors").WHERE(
            selection, receptor
        )
        info = self.db_query(*query.build()).fetchone()
        if info is None:
            info = [], []
        return info

    @staticmethod
    def _serialize_pose_coordinates(coords):
        """Serialize a list of [x, y, z] floats for storage in the pose_coordinates
        column. Base/DuckDB store the native nested list (DuckDB FLOAT[][]); SQLite
        overrides this to pack a float32 BLOB. Returns the value unchanged here.
        """
        return coords

    @staticmethod
    def _deserialize_pose_coordinates(stored):
        """Inverse of _serialize_pose_coordinates: return a list of [x, y, z] floats.
        Tolerates legacy JSON-text values from databases written before the native
        coordinate format (so un-migrated DuckDB databases still read).
        """
        if stored is None:
            return None
        if isinstance(stored, str):
            return json.loads(stored)
        return [list(atom) for atom in stored]

    def fetch_rdkit_pose_properties(self, pose_ids) -> list:
        """
        Gets molecular data needed to create rdkit mols for the given pose(s), including
        ligand mol binary and name via a JOIN with Ligands.

        Args:
            pose_ids (int | list[int]): one or more pose ids

        Returns:
            list of (pose_id, rdmol_binary, ligname, docking_score, leff, pose_coordinates, flexible_res_coordinates, pose_rank)
        """
        if isinstance(pose_ids, int):
            pose_ids = [pose_ids]
        placeholders = ",".join(["?"] * len(pose_ids))
        query = self.QueryBuilder()
        query.SELECT(
            "R.pose_id",
            "L.rdmol",
            "L.ligname",
            "R.docking_score",
            "R.leff",
            "R.pose_coordinates",
            "R.flexible_res_coordinates",
            "R.pose_rank",
        ).FROM("Results", "R").JOIN("Ligands", "L", "ligand_id").WHERE(
            f"R.pose_id IN ({placeholders})", *pose_ids
        )
        rows = self.db_query(*query.build()).fetchall()
        # pose_coordinates (index 5) is stored natively; return it as a list of [x,y,z]
        return [
            (*row[:5], self._deserialize_pose_coordinates(row[5]), *row[6:])
            for row in rows
        ]

    def fetch_columns_from_table_as_dicts(
        self, table: str, columns: list, length: int = None, starting_rowid: int = 0
    ) -> tuple[list[str], list[dict]]:
        """Get requested columns from Results, a bookmark, or a status table.

        Automatically JOINs Ligands when Ligands-only columns (ligname, ligand_smile,
        rdmol) are requested, and applies the appropriate WHERE/JOIN for bookmarks and
        status tables. Pass "status" as a column name to include pose status.

        Args:
            table (str): Results table name, bookmark name, or status table name
            columns (list): column names to select; "status" is a special sentinel
            length (int, optional): max rows to return. Defaults to None (all rows).
            starting_rowid (int, optional): R.rowid lower bound. Defaults to 0.

        Returns:
            tuple[list[str], list[dict]]: column names and list of row dicts
        """
        # Direct physical table that is not Results/bookmark/status — naive query
        if (
            self._is_table(table)
            and not self.is_bookmark(table)
            and not self._is_statustable(table)
            and table.lower() != "results"
        ):
            naive = self.QueryBuilder()
            naive.SELECT(",".join(columns)).FROM(table)
            if length:
                naive.LIMIT(length)
            if starting_rowid:
                naive.WHERE(f"{table}.rowid >= {starting_rowid}")
            return self.get_query_data_as_dicts(naive.build()[0])

        include_status = "status" in columns
        data_cols = [c for c in columns if c != "status"]
        needs_ligands = any(c in LIGANDS_ONLY_COLS for c in data_cols)

        def _qualify(col):
            return f"L.{col}" if col in LIGANDS_ONLY_COLS else f"R.{col}"

        query = self.QueryBuilder()
        for col in data_cols:
            query.SELECT(_qualify(col))
        query.FROM(RESULTS_SCHEMA.name, "R")
        if needs_ligands:
            query.JOIN(LIGANDS_SCHEMA.name, "L", "ligand_id")

        if include_status:
            if self._is_statustable(table):
                query.SELECT(f"'{table.lower()}' as status")
            else:
                query.SELECT_STATUS()
        self._apply_table_filter(query, table)

        if length:
            query.LIMIT(length)
        if starting_rowid:
            query.WHERE(f"R.rowid >= {starting_rowid}")

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
            WHEN EXISTS (SELECT 1 FROM Accepted s WHERE s.pose_id = R.pose_id) THEN 'accepted'
            WHEN EXISTS (SELECT 1 FROM Rejected s WHERE s.pose_id = R.pose_id) THEN 'rejected'
            WHEN EXISTS (SELECT 1 FROM Maybe s WHERE s.pose_id = R.pose_id) THEN 'maybe'
            ELSE ''
        END AS status,"""

        if table.lower() == "results":
            rowid = "R.rowid"

        elif self._is_statustable(table):
            query.JOIN(table, "T", "pose_id")
            rowid = "T.rowid"
            # status assignement doesn't make sense for status tables
            status_assignement = f"""'{table.lower()}' AS status,"""

        elif self._is_candidates_table(table):
            query.JOIN(CANDIDATES_SUBQ, "T", "pose_id")
            rowid = "R.rowid"
            # keep CASE status_assignement — shows real accepted/maybe per row

        elif self.is_bookmark(table):
            query.JOIN("filtered_poses", "fp", "pose_id").JOIN(
                "filters", "f", "filter_id", "filtered_poses"
            ).WHERE("f.name = ?", table)
            rowid = "fp.rowid"

        else:
            raise ValueError(
                f"Table '{table}' is not a results table, status table, or bookmark."
            )

        query.JOIN("Pose_comments", "Pc", "pose_id", kind="LEFT")

        ordered_columns = f"""
        {status_assignement}
        L.ligname, R.pose_id,  R.pose_rank,
        R.docking_score, R.leff, COALESCE(Pc.comment, '') AS comment, R.cluster_size,
        R.cluster_rmsd, R.num_hb, R.receptor, R.run_number,
        R.delta, R.num_interactions, R.unbound_energy,
        R.reference_RMSD, R.energies_inter, R.energies_vdw,
        R.energies_electro, R.energies_flexLig, R.energies_flexLR,
        R.energies_intra, R.energies_torsional, {rowid}"""

        query.SELECT(ordered_columns)

        query.JOIN("Ligands", "L", "ligand_id", "results").WHERE(
            f"{rowid} {where_operator} ?", starting_rowid
        ).ORDER_BY(rowid).LIMIT(length).DESC(reverse)

        cursor = self.db_query(*query.build())
        headers = [desc[0] for desc in cursor.description]
        data = cursor.fetchall()
        return {"headers": headers, "data": data}

    def fetch_lignames_and_poses_for_selection(
        self,
        selection: str = None,
        ligand_names: list[str] = None,
        pose_ids: list[int] = None,
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
        # general selection
        query.SELECT("L.Ligname", "r.pose_id").FROM("Ligands", "L").JOIN(
            "Results", "R", "ligand_id"
        )
        # narrow it by selection table/bookmark
        if selection is not None:
            if self.is_bookmark(selection):
                query.IN_BOOKMARK(selection)
            elif selection.lower() in statuses.values():
                query.WHERE(f"R.pose_id IN {selection}")
            elif self._is_candidates_table(selection):
                query.JOIN(CANDIDATES_SUBQ, "T", "pose_id")
            else:
                logger.error(
                    f"-{selection}- is not a valid selection for this method. Please provide a bookmark_name or a status table."
                )
                return
        # narrow to pose id(s)
        if pose_ids is not None:
            placeholders = ", ".join("?" * len(pose_ids))
            query.WHERE(f"R.pose_id IN ({placeholders})", *pose_ids)
        # narrow to ligand name(s)
        if ligand_names is not None:
            placeholders = ", ".join("?" * len(ligand_names))
            query.WHERE(f"L.ligname IN ({placeholders})", *ligand_names)

        rows = self.db_query(*query.build()).fetchall()
        ligand_poses = defaultdict(list)
        for name, id in rows:
            ligand_poses[name].append(id)
        return dict(ligand_poses)

    def fetch_bookmark_interactions(self, bookmark_name: str) -> pd.DataFrame:
        """
        Get all interactions represented in a bookmark (as opposed to the whole database)
        Mostly used for GUI applications

        Args:
            bookmark_name (str): from which to get pose interactions from

        Returns:
            pd.DataFrame: interactions presented as a dataframe
        """
        query = f"""
        SELECT DISTINCT 
            II.interaction_type,
            II.rec_chain,
            II.rec_resname,
            II.rec_resid,
            II.rec_atom
        FROM Interaction_indices AS II
        JOIN Interactions AS I ON I.interaction_id=II.interaction_id
        WHERE I.pose_id IN ({self.QueryBuilder.bookmark_query(bookmark_name)});"""
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
            """DELETE FROM Maybe WHERE pose_id = ?;""", pose_ids, commit=False
        )
        self.db_update(
            """DELETE FROM Rejected WHERE pose_id = ?;""", pose_ids, commit=True
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
            """DELETE FROM Accepted WHERE pose_id = ?;""", pose_ids, commit=False
        )
        self.db_update(
            """DELETE FROM Rejected WHERE pose_id = ?;""", pose_ids, commit=True
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
            """DELETE FROM Accepted WHERE pose_id = ?;""", pose_ids, commit=False
        )
        self.db_update(
            """DELETE FROM Maybe WHERE pose_id = ?;""", pose_ids, commit=True
        )

    def remove_status(self, pose_ids: Union[int, list[int]]):
        """
        Will remove status associated with a pose

        Args:
            pose_ids (Union[int, list[int]])
        """
        if isinstance(pose_ids, int):
            pose_ids = [pose_ids]
        pose_ids = [(pose,) for pose in pose_ids]
        self.db_update(
            "DELETE FROM Accepted WHERE pose_id = ?;", pose_ids, commit=False
        )
        self.db_update("DELETE FROM Maybe WHERE pose_id = ?;", pose_ids, commit=False)
        self.db_update("DELETE FROM Rejected WHERE pose_id = ?;", pose_ids, commit=True)

    def get_candidate_lignames(self) -> list[str]:
        """
        #TODO doc string

        Returns:
            list[str]: _description_
        """
        sql = """SELECT DISTINCT ligname FROM Ligands JOIN Results ON Ligands.ligand_id = Results.ligand_id WHERE Results.pose_id IN ({0});""".format(
            CANDIDATES_SUBQ
        )
        return [row[0] for row in self.db_query(sql).fetchall()]

    def prepare_column_export_query(
        self, columns: dict[str, str], bookmark: str
    ) -> str:
        """
        Writes a query based on what columns are requested. Resolves each column
        to its source table automatically.

        Args:
            columns (list): list of column name strings to export
            bookmark (str): bookmark or status table to filter results by

        Returns:
            str: SQL query string
        """
        # Build a flat map of column_name -> table from all queryable tables
        column_to_table = {
            col: table for table, cols in OUTFIELD_BY_TABLE.items() for col in cols
        }

        columns_by_table = {}
        for col in columns:
            table = column_to_table.get(col.lower())
            if table is None:
                raise OptionError(f"Column '{col}' not found in any queryable table.")
            columns_by_table.setdefault(table, []).append(col)

        included_tables = [t.lower() for t in columns_by_table]

        query = self.QueryBuilder()
        for col in columns:
            query.SELECT(col)
        query.FROM("Results")

        # join to tables that are represented in the columns
        if "ligands" in included_tables:
            query.JOIN("Ligands", "Ligands", on="ligand_id", to="Results")
        if "interaction_indices" in included_tables:
            query.JOIN("Interactions", "Interactions", "pose_id", "Results")
            query.JOIN(
                "Interaction_indices",
                "Interaction_indices",
                "interaction_id",
                "Interactions",
            )

        if self.is_bookmark(bookmark):
            query.IN_BOOKMARK(bookmark)
        elif self._is_statustable(bookmark):
            query.JOIN(bookmark, bookmark, "pose_id")

        return query.build()[0]

    def clear_interaction_tables(self) -> bool:
        """
        Clears (deletes and reuilds) the two interaction tables

        Raises:
            StorageError: _description_

        Returns:
            bool: True if clearing was successful
        """
        try:
            self._delete_table("Interaction_indices")
            self._delete_table("Interactions")
            self._create_interaction_table()
            self._create_interaction_index_table()
            self.conn.commit()
        except Exception as e:
            raise StorageError(f"Error during interaction tables clearing: {e}") from e
        if (
            self.table_length("Interactions") > 0
            or self.table_length("Interaction_indices") > 0
        ):
            return False
        else:
            return True

    def ensure_gui_tables(self) -> None:
        """
        Ensures gui-specific tables exist in database
        """

        for sql in build_create_table(
            POSE_COMMENTS_SCHEMA.name, POSE_COMMENTS_SCHEMA, self.dialect
        ):
            self.db_query(sql)
        self.conn.commit()

    # endregion

    # region virtual public api
    @abstractmethod
    def close_storage(self, attached_db=None, vacuum=False):
        """Close connection to database

        Args:
            attached_db (str, optional): name of attached DB (not including file extension)
            vacuum (bool, optional): indicates that database should be vacuumed before closing
        """
        ...

    @abstractmethod
    def check_ringtaildb_version(self):
        """
        Checks the database version and confirms whether the code base is compatible with it

        Returns:
            bool: whether or not db is compatible with the code base
            str: current database version
        """
        ...

    @abstractmethod
    def insert_receptor_blob(self, receptor: bytes, rec_name: str):
        """Takes object of Receptor class, updates the column in Receptor table

        Args:
            receptor (bytes): bytes receptor object to be inserted into DB
            rec_name (string): Name of receptor. Used to insert into correct row of DB
        """
        ...

    @abstractmethod
    def insert_receptor_polymer(self, receptor: str, rec_name: str):
        """Takes object of Receptor class, updates the column in Receptor table

        Args:
            receptor (str): json string representation of a receptor meeko.Polymer oobject to be inserted into DB
            rec_name (str): Name of receptor. Used to insert into correct row of DB
        """
        ...

    @abstractmethod
    def clone(self, backup_name: str = None) -> str:
        """Creates a copy of the db

        Args:
            backup_name (str, optional): name of the cloned database

        Returns:
            str: path of backed up database
        """
        ...

    @abstractmethod
    def create_transaction_tracking_table(self, table_name: str = "tracking_table"):
        """
        Creates a table to track pose ids to which are part of a multiple transaction
        database operation, which may occur over multiple commits and database connections.

        Args:
            table_name (str, optional): Name of the table. Defaults to "tracking_table".
        """
        ...

    def update_database_version(self, new_version: str, backup=False):
        """Updates database schema from older versions to new_version.

        Args:
            new_version (str): target version string, e.g. "3.0.0"
            backup (bool, optional): clone database before upgrading. Defaults to False.
        """
        raise NotImplementedError

    @abstractmethod
    def db_query(self, query: str, params: tuple = ()) -> Any:
        """Executes provided sql query. Returns cursor for results.

        Args:
            query (str): Formated sql query as string
            params (tuple): parameters to be used in query (assumes query as appropriate place holders)

        Returns:
            Any: cursor containing results of query
        """
        ...

    @abstractmethod
    def db_update(self, query: str, parameters: list[tuple], commit=True) -> None:
        """
        A db query that uses executemany

        Args:
            query (str): sql formatted query string
            parameters (list[tuple]): assumes appropriate place holders in query
            commit (bool, optional): whether to commit the transaction in open connection. Defaults to True.

        Raises:
            OptionError
            DatabaseInsertionError
        """
        ...

    @abstractmethod
    def prepare_for_merging(self):
        """
        Prepares the database for merging
        """
        ...

    @abstractmethod
    def table_length(self, table: str) -> int:
        """
        Get length of table or bookmark

        Args:
            table (str): name of table or bookmark

        Returns:
            int: number of poses in table/bookmark
        """
        ...

    @abstractmethod
    def fetch_cluster_options(self, ligname: str) -> list[tuple]:
        """Return available clustering groups for selection.

        Args:
            ligname (str): ligname for ligand to find similarity with

        Returns:
            list[tuple]: list of (cluster_id, cluster_window, name) tuples
        """
        ...

    @abstractmethod
    def fetch_clustered_similars(
        self, ligname: str, cluster_id: int
    ) -> tuple[list, str, str]:
        """Given ligname and a chosen cluster_id, return similar ligands.

        Args:
            ligname (str): ligname for ligand to find similarity with
            cluster_id (int): cluster to search within

        Returns:
            tuple[list, str, str]: (ligands, bookmark_name, cluster_name)
        """
        ...

    @abstractmethod
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
        ...

    @abstractmethod
    def pose_row_in_table(self, table: str, pose_id: int) -> Union[None, int]:
        """
        Find the row id of a pose in a given table

        Args:
            table (str)
            pose_id (int)

        Returns:
            Union[None, int]: rowid if any
        """
        ...

    @abstractmethod
    def tables_in_db(self) -> list:
        """
        Returns a list of all table names in the database

        Returns:
            list: list of table names
        """
        ...

    @abstractmethod
    def get_starting_rowid(self, table: str) -> int:
        """
        Starting row id for a table, will be 1 for regular tables, and 1 or non-1 for bookmarks
        (whose rows are inside Filtered_poses)

        Args:
            table (str): table or bookmark name

        Returns:
            int: first row id belonging to that selection
        """
        ...

    @abstractmethod
    def set_pose_comment(self, pose_id: int, comment: str) -> None:
        """Write (or clear) a user comment for a pose."""
        ...

    @abstractmethod
    def get_pose_comment(self, pose_id: int) -> "str | None":
        """Return the comment for a pose, or None if none exists."""
        ...

    @abstractmethod
    def detach_db(self, new_db_alias):
        """Detaches new database file from current database

        Args:
            new_db_name (str): db name for database to detach

        Raises:
            StorageError
        """
        ...

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

    def _create_ligands_table(self, name=LIGANDS_SCHEMA.name):
        """Create table for ligands

        Args:
            name (str, optional): Defaults to "Ligands".
        """
        for sql in build_create_table(name, LIGANDS_SCHEMA, self.dialect):
            self.db_query(sql)

    def _create_results_table(self, name=RESULTS_SCHEMA.name):
        """Creates table for results."""
        for sql in build_create_table(name, RESULTS_SCHEMA, self.dialect):
            self.db_query(sql)

    def _create_receptors_table(self):
        """Create table for receptors."""
        for sql in build_create_table(
            RECEPTORS_SCHEMA.name, RECEPTORS_SCHEMA, self.dialect
        ):
            self.db_query(sql)

    def _create_interaction_index_table(self):
        """Creates a table describing unique interactions in the database."""
        for sql in build_create_table(
            INTERACTION_INDICES_SCHEMA.name, INTERACTION_INDICES_SCHEMA, self.dialect
        ):
            self.db_query(sql)

    def _create_interaction_table(self):
        """Creates a table of each pose-interaction combination."""
        for sql in build_create_table(
            INTERACTIONS_SCHEMA.name, INTERACTIONS_SCHEMA, self.dialect
        ):
            self.db_query(sql)

    def _create_db_properties_table(self):
        """Create table of database properties used during write session to the database."""
        for sql in build_create_table(
            DB_PROPERTIES_SCHEMA.name, DB_PROPERTIES_SCHEMA, self.dialect
        ):
            self.db_query(sql)

    def _create_filtering_tables(self):
        """
        Creates Filters (bookmark metadata) and Filtered_poses (bookmark members) tables.
        """
        for sql in build_create_table(
            FILTERS_SCHEMA.name, FILTERS_SCHEMA, self.dialect
        ):
            self.db_query(sql)
        for sql in build_create_table(
            FILTERED_POSES_SCHEMA.name, FILTERED_POSES_SCHEMA, self.dialect
        ):
            self.db_query(sql)

    def _create_cluster_tables(self):
        """
        Creates cluster tables if they don't already exist.
        """
        for sql in build_create_table(
            CLUSTERS_SCHEMA.name, CLUSTERS_SCHEMA, self.dialect
        ):
            self.db_query(sql)
        for sql in build_create_table(
            CLUSTER_GROUPS_SCHEMA.name, CLUSTER_GROUPS_SCHEMA, self.dialect
        ):
            self.db_query(sql)
        for sql in build_create_table(
            POSE_CLUSTERS_SCHEMA.name, POSE_CLUSTERS_SCHEMA, self.dialect
        ):
            self.db_query(sql)

    def _create_status_tables(self) -> None:
        """
        Creates pose status tables (Accepted, Maybe, Rejected) if needed.
        """
        for name in ("Accepted", "Maybe", "Rejected"):
            for sql in build_create_table(name, STATUS_TABLE_SCHEMA, self.dialect):
                self.db_query(sql)
        self.conn.commit()

    def _insert_docking_data(
        self, results: list[list], interactions: list[list], duplicate_handling: str
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
        # handle duplicates if requested
        if duplicate_handling and duplicate_handling.lower() == "replace":
            # first delete duplicate results and interactions
            # then insert all the new ones indiscrimenately
            self._delete_old_duplicate_results()
        elif duplicate_handling and duplicate_handling.lower() == "ignore":
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
        existing_cluster_id_tuple = self._cluster_exists(cluster_name, bookmark_name)
        if existing_cluster_id_tuple:
            logger.warning(
                f"A previous cluster bookmark under the same name exists ({cluster_bookmark}), and will be deleted."
            )
            existing_cluster_id = existing_cluster_id_tuple[0]
            self._delete_cluster(existing_cluster_id)
            self.delete_bookmark(cluster_bookmark)

        cluster_id = self._insert_new_cluster_info(
            cluster_name, "", bookmark_name, len(clusters)
        )
        logger.debug(
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
        self._delete_table(MERGED_TABLES_SCHEMA.name)
        self._delete_table(PK_CONVERSIONS_SCHEMA.name)
        self._delete_table(POSE_CLUSTERS_SCHEMA.name)
        self._delete_table(CLUSTER_GROUPS_SCHEMA.name)
        self._delete_table(CLUSTERS_SCHEMA.name)
        self._delete_table(FILTERED_POSES_SCHEMA.name)
        self._delete_table(FILTERS_SCHEMA.name)
        self._delete_table(INTERACTIONS_SCHEMA.name)
        self._delete_table(INTERACTION_INDICES_SCHEMA.name)
        for status in statuses.values():
            if status:
                self._delete_table(status.capitalize())
        self._delete_table(RESULTS_SCHEMA.name)
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
            db_alias (str): Optional alias for attached databases

        """

        if db_alias:
            name = db_alias + "." + table_name
        else:
            name = table_name
        query = self.QueryBuilder()
        query.DROP_IF_EXISTS(name)
        return self.db_query(query.build()[0])

    def _get_numeric_columns(self, table: str) -> list:
        """
        Get all numeric (float, int) columns from a table

        Args:
            table (str)

        Returns:
            list: column names that are numeric
        """
        return [
            name
            for name, col in TABLE_SCHEMAS[table.lower()].columns.items()
            if col.sql_type in NUMERIC_TYPES
        ]

    def _fetch_table_column_names(self, table: str) -> list:
        """
        Get all column names in a table

        Args:
            table (str)

        Returns:
            list: column names in table schema
        """
        return list(TABLE_SCHEMAS[table.lower()].columns)

    def _apply_table_filter(self, query: QueryBuilder, table: str) -> QueryBuilder:
        """Apply the appropriate pose filter for bookmark, status table, or candidates table.

        Args:
            query (QueryBuilder): query to modify in-place
            table (str): bookmark name, status table name, or candidates table name

        Returns:
            QueryBuilder: the same query, modified
        """
        if self.is_bookmark(table):
            query.IN_BOOKMARK(table)
        elif self._is_statustable(table):
            query.JOIN(table, "T", "pose_id")
        elif self._is_candidates_table(table):
            query.JOIN(CANDIDATES_SUBQ, "T", "pose_id")
        return query

    def _is_statustable(self, table: str) -> bool:
        """
        Returns True if table name is actually a status table (table with poses who have been assigned a status like accept, reject, maybe)

        Args:
            table (str): name of table or bookmark to check

        Returns:
            bool: if table name is a status table
        """
        if table.lower() in statuses.values():
            return True
        else:
            return False

    def _is_candidates_table(self, table: str) -> bool:
        """
        Returns True if table name is the virtual Candidates table (Accepted ∪ Maybe).

        Args:
            table (str): tabale name to check

            Returns:
                bool: if this is candidates table or not
        """
        return table.lower() == CANDIDATES_NAME

    def _format_orderby(self, column_name: str) -> Union[str, None]:
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
        max_miss = 0
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
                if filter_key == "ligand_name":
                    from .util import iterate_nested

                    filter_value = [
                        name
                        for item in iterate_nested(filter_value)
                        for name in item.split(",")
                    ]
                    if len(filter_value) > 50:
                        raise OptionError(
                            "The number of provided ligand names is too large, please prepare as a csv and use 'ligand_name_file' instead."
                        )
                ligand_filters[filter_key] = filter_value

            if filter_key == "max_miss":
                max_miss = filter_value
            # parse ligand name file if present
            if filter_key == "ligand_name_file":
                if not Path(filter_value).suffix.lower() == ".csv":
                    raise OptionError(
                        f"The file of ligand names needs to be a csv file, cannot proceed with {Path(filter_value).suffix.lower()}."
                    )
                else:
                    ligand_filters[filter_key] = filter_value
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
        self, filters_dict: dict, bookmark_name: str, input_bookmark: str = None
    ) -> str:
        """Takes dict of filters, writes SQL filtering string.

        Each filter type becomes an independent WHERE condition combined with AND:
        - Numerical: direct column conditions on Results
        - Include interactions: pose_id IN (GROUP BY/HAVING subquery)
        - Exclude interactions: pose_id NOT IN (simple subquery)
        - Ligand name: ligand_id IN (Ligands subquery)
        - RDKit (substruct/position/atoms): in-memory post-filter

        Args:
            filters_dict (dict): Keys and value formats must match Filters class
            bookmark_name (str): Name of bookmark to which passing poses are saved
            input_bookmark (str, optional): Filter over pre-filtered data instead of entire database

        Returns:
            str: SQL query returning pose_ids of passing results
        """
        filtering_window = "Results"

        if input_bookmark is not None:
            if input_bookmark == bookmark_name:
                logger.error(
                    f"Specified 'input_bookmark' and 'bookmark_name' are the same: {bookmark_name}"
                )
                raise OptionError(
                    "'input_bookmark' and 'bookmark_name' cannot be the same! Please rename 'bookmark_name'"
                )
            if (
                filters_dict["score_percentile"] is not None
                or filters_dict["le_percentile"] is not None
            ):
                raise OptionError(
                    "Cannot use 'score_percentile' or 'le_percentile' with 'input_bookmark'."
                )
            if self.is_bookmark(input_bookmark):
                filtering_window = f"(SELECT * FROM Results WHERE pose_id IN ({self.QueryBuilder.bookmark_query(input_bookmark)}))"
            elif self._is_statustable(input_bookmark):
                filtering_window = f"(SELECT * FROM Results WHERE pose_id IN (SELECT pose_id FROM {input_bookmark}))"

        processed_filters = self._process_filters_for_query(filters_dict)
        if not processed_filters:
            raise OptionError("No filters were provided, cannot filter.")

        where_conditions = []
        lig_filters = {}
        rdkit_query = False

        if "num_filters" in processed_filters:
            where_conditions.append(
                " AND ".join(f"R.{f}" for f in processed_filters["num_filters"])
            )

        if "int_filters" in processed_filters:
            include_interactions, exclude_interactions = (
                self._prepare_interaction_indices_for_filtering(
                    interaction_list=processed_filters["int_filters"]
                )
            )
            if include_interactions:
                # EXISTS over the (interaction_id, pose_id) covering index. SQLite does a
                # fast per-pose point-lookup; DuckDB's optimizer de-correlates this into a
                # hash semi-join, so the same form is optimal for both dialects.
                where_conditions.append(
                    self._build_include_interaction_exists_conditions(
                        include_interactions, processed_filters["max_miss"]
                    )
                )
            if exclude_interactions:
                where_conditions.append(
                    f"R.pose_id NOT IN ({self._build_exclude_interaction_subquery(exclude_interactions)})"
                )

        if "lig_filters" in processed_filters:
            lig_filters = processed_filters["lig_filters"]
            ligname_query = ""
            if "ligand_name" in lig_filters:
                lig_names = lig_filters.pop("ligand_name")
                ligname_query = "SELECT ligand_id FROM Ligands WHERE " + " OR ".join(
                    f"ligname LIKE '%{ligname}%'" for ligname in lig_names if ligname
                )
            elif "ligand_name_file" in lig_filters:
                csv_path = lig_filters.pop("ligand_name_file")
                self._create_ligname_temp_table(csv_path)
                ligname_query = "SELECT ligand_id FROM Ligands JOIN tmp_lignames ON ligname = tmp_lignames.ligandname"
            if ligname_query:
                where_conditions.append(f"R.ligand_id IN ({ligname_query})")
            if lig_filters:
                rdkit_query = True

        if rdkit_query:
            rdkit_partial = f" FROM {filtering_window} R JOIN Ligands ON Ligands.ligand_id = R.ligand_id"
            if where_conditions:
                rdkit_partial += f" WHERE {' AND '.join(where_conditions)}"
            passing_pose_ids = [
                pid
                for poseids in self._perform_rdkit_filtering(
                    rdkit_partial, lig_filters
                ).values()
                for pid in poseids
            ]
            where_conditions = [
                f"R.pose_id IN ({','.join(map(str, passing_pose_ids))})"
            ]

        where_sql = (
            f"WHERE {' AND '.join(where_conditions)}" if where_conditions else ""
        )
        select_cols = (
            "R.pose_id, R.ligand_id"
            if self._filtered_poses_has_ligand_id
            else "R.pose_id"
        )
        return f"SELECT {select_cols} FROM {filtering_window} R {where_sql}"

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
        _substructure_position_calculation: checks if a substructure is in the right location

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

        select_statement = (
            "SELECT R.pose_id, Ligands.ligand_id, Ligands.rdmol, Ligands.ligand_smile"
        )
        headers = ["pose_id", "ligand_id", "rdmol", "ligand_smile"]
        # handle substruct
        if "ligand_substruct" in ligand_filters:
            for substruct in ligand_filters["ligand_substruct"]:
                substruct_mols.append(_smarts_to_mol(substruct))
        # whether or not doing a substruct position filter
        if "ligand_substruct_pos" in ligand_filters:
            substruct_pos = ligand_filters["ligand_substruct_pos"]
            position = True
            # we need additional info if doing position search
            select_statement += ", R.pose_coordinates "
            headers.append("pose_coordinates")
        # build full query
        query = select_statement + partial_query

        def _substructure_position_calculation(
            pose_coordinates, hit_atom_indices, filter
        ) -> bool:
            """
            Method that checks whether or not a substructure specified by smarts
            is present on a ligand in the specified location (by means of the filter values).

            Args:
                pose_coordinates (json): ligand coordinates matching rdkit smiles indices
                hit_atom_indices (tuple[int]): indices for each of the atoms in the smarts pattern
                filter (list[str]): filter values from user

            Returns:
                bool: whether or not the smarts in the pose qualified in the right position
            """
            # unpack filter values
            # index in the smarts pattern that should be within coordinates
            index = filter[0]
            sqdist = filter[1] ** 2
            x = filter[2]
            y = filter[3]
            z = filter[4]
            # calculate xyz space coordinates (pose_coordinates already a list of [x,y,z])
            xyz = [
                float(value)
                for value in pose_coordinates[hit_atom_indices[index]]
            ]
            # calculate the sum of squares distances
            d2 = (xyz[0] - x) ** 2 + (xyz[1] - y) ** 2 + (xyz[2] - z) ** 2
            if d2 <= sqdist:
                return True
            else:
                return False

        ligands_checked = 0
        filtered_ligands = {}
        for ligandrow in self._stream_query(query):
            ligandict = dict(zip(headers, ligandrow))
            # substruct and maxatoms do not discriminate on poses, check if ligand has already been accounted for
            if not position and ligandict["ligand_id"] in filtered_ligands:
                filtered_ligands[ligandict["ligand_id"]].append(ligandict["pose_id"])
            # the real workhorse
            else:
                # deserialize binary rdmol
                ligand_mol = Chem.Mol(ligandict["rdmol"])
                ligand_mol = Chem.RemoveHs(ligand_mol)
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
                        pose_coordinates = self._deserialize_pose_coordinates(
                            ligandict["pose_coordinates"]
                        )
                        # filterrow [1:] should be indices, distance allowance, and coordinates for smarts match
                        substruct_pos_filter = filterrow[1:]
                        for hit_indices in ligand_mol.GetSubstructMatches(smarts_mol):
                            filter_match = _substructure_position_calculation(
                                pose_coordinates,
                                hit_indices,
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

    def _get_bookmark_poses_query(self, bookmark_name: str, alias: str = "") -> str:
        """
        Creates a query that retrieves all poses from a bookmark, that can be used in other queries

        Args:
            bookmark_name (str): bookmark for which to create the query

        Returns:
            str: query representing the poses in a bookmark
        """
        return self.QueryBuilder.bookmark_query(bookmark_name, alias)

    def _generate_interaction_bitvectors(self, pose_ids: tuple[str]) -> dict:
        """
        Method to generate a dict of generate bitvector strings from pose_ids

        Args:
            pose_ids (str): query formatted list of pose_ids (as tuple)

        Returns:
            dict: of "pose_id":"bitvector"
        """
        # create a list of 0 items the length of interaction_indices table
        query = self.QueryBuilder()
        query.SELECT("pose_id", "interaction_id").FROM("Interactions").WHERE(
            f"""pose_id IN ({",".join(["?"] * len(pose_ids))})""", *pose_ids
        )

        iis = [
            row[0]
            for row in self.db_query(
                "SELECT interaction_id FROM Interaction_indices;"
            ).fetchall()
        ]
        ii_to_index = {ii: idx for idx, ii in enumerate(iis)}
        ii_length = len(iis)

        poseid_intinds = self.db_query(*query.build()).fetchall()

        # build bitvectors with NumPy
        poseid_bv = {}
        for pose_id in pose_ids:
            bv = np.zeros(ii_length, dtype=np.uint8)
            poseid_bv[pose_id] = bv

        for pose_id, interaction_id in poseid_intinds:
            if interaction_id in ii_to_index:
                poseid_bv[pose_id][ii_to_index[interaction_id]] = 1

        # convert to string
        poseid_bv = {key: "".join(map(str, value)) for key, value in poseid_bv.items()}
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

    def _build_include_interaction_exists_conditions(
        self, include_groups: list, max_miss: int
    ) -> str:
        """Build correlated EXISTS conditions checking per-pose interaction requirements.

        Uses the (interaction_id, pose_id) covering index: for each candidate pose SQLite
        does a fast point-lookup instead of materializing the full interaction result set.
        For selective energy filters this is 10–100x faster than the IN (subquery) approach.

        Args:
            include_groups (list[list[int]]): each inner list is a group of interaction_ids
                that satisfy one required interaction
            max_miss (int): number of groups a pose is allowed to miss

        Returns:
            str: SQL WHERE fragment using R.pose_id (caller must alias the Results table as R)
        """

        def _group_exists(ids: list) -> str:
            return (
                f"EXISTS (SELECT 1 FROM Interactions I "
                f"WHERE I.interaction_id IN ({numlist2str(ids, ',')}) AND I.pose_id = R.pose_id)"
            )

        if max_miss == 0:
            return " AND ".join(_group_exists(g) for g in include_groups)

        # General case: pose must satisfy at least (n_groups - max_miss) groups
        cases = " + ".join(
            f"(CASE WHEN {_group_exists(g)} THEN 1 ELSE 0 END)" for g in include_groups
        )
        n_required = len(include_groups) - max_miss
        return f"({cases}) >= {n_required}"

    def _build_include_interaction_subquery(
        self, include_groups: list, max_miss: int
    ) -> str:
        """Build a subquery returning pose_ids that satisfy at least (n_groups - max_miss)
        of the required interaction groups. Wildcards within a group are collapsed to a
        single count via CASE group numbering.

        Args:
            include_groups (list[list[int]]): each inner list is a group of interaction_ids
                that satisfy one required interaction (multiple IDs arise from wildcard matches)
            max_miss (int): number of groups a pose is allowed to miss

        Returns:
            str: SQL subquery returning qualifying pose_ids
        """
        # Simple case: one group, no misses — SELECT DISTINCT avoids GROUP BY + COUNT(DISTINCT) temp B-trees
        if len(include_groups) == 1 and max_miss == 0:
            ids_str = numlist2str(include_groups[0], ",")
            return f"SELECT DISTINCT pose_id FROM Interactions WHERE interaction_id IN ({ids_str})"

        all_ids = [iid for group in include_groups for iid in group]
        case_clauses = " ".join(
            f"WHEN interaction_id IN ({numlist2str(group, ',')}) THEN {i + 1}"
            for i, group in enumerate(include_groups)
        )
        n_required = len(include_groups) - max_miss
        return (
            f"SELECT pose_id FROM Interactions "
            f"WHERE interaction_id IN ({numlist2str(all_ids, ',')}) "
            f"GROUP BY pose_id "
            f"HAVING COUNT(DISTINCT CASE {case_clauses} END) >= {n_required}"
        )

    def _build_exclude_interaction_subquery(self, exclude_groups: list) -> str:
        """Build a subquery returning pose_ids that have any of the excluded interactions.
        Used as NOT IN to reject poses with forbidden interactions.

        Args:
            exclude_groups (list[list[int]]): each inner list is a group of interaction_ids
                for one forbidden interaction (multiple IDs from wildcard matches)

        Returns:
            str: SQL subquery returning pose_ids to exclude
        """
        all_ids = [iid for group in exclude_groups for iid in group]
        return (
            f"SELECT DISTINCT pose_id FROM Interactions "
            f"WHERE interaction_id IN ({numlist2str(all_ids, ',')})"
        )

    def _get_interaction_indices(self, interaction_list: list) -> list[tuple]:
        """takes list of interaction info and looks up corresponding interaction index

        Args:
            interaction_list (list): List containing interaction info
                in format [<interaction_type>, <rec_chain>, <rec_resname>,
                <rec_resid>, <rec_atom>]

        Returns:
            list[tuple]: list of tuples with the interaction index/indices
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
        self, name, query: str, filters={}, input_bookmark: str = ""
    ) -> bool:
        """
        Will run a filter query and determine if there are passing poses, in which case all relevant
        data is written to the database

        Args:
            name (str): name of new bookmark
            query (str): query that defines what poses to insert
            filters (dict, optional): filters or restrictions used
            input_bookmark (str, optional): If filters were performed across an existing obokmark. Defaults to None.

        Raises:
            StorageError
            OptionError

        Returns:
            int: number of passing poses
        """
        # make sure bookmark name is not a table name
        if name in self.tables_in_db():
            raise OptionError(
                f"Bookmark name {name} is the same as an existing table in the database, and cannot be used."
            )

        # Capture the old filter_id now so we can delete it only AFTER confirming
        # the new filter has results — prevents losing the old bookmark on empty filters.
        old_filter_id = None
        if self.is_bookmark(name):
            logger.warning(
                f"The bookmark {name} already exists, and will be overwritten by the current filter."
            )
            old_filter_id = self.conn.execute(
                "SELECT filter_id FROM Filters WHERE name = ?", (name.lower(),)
            ).fetchone()[0]

        self._begin_transaction()
        insert_query = self.QueryBuilder()
        insert_query.INSERT_INTO(
            "Filters", "name", "query", "filters", "filter_window"
        ).RETURNING("filter_id")

        # insert filter info, return filter_id
        filter_id = self.conn.execute(
            insert_query.build()[0],
            (name.lower(), query, json.dumps(filters), input_bookmark),
        ).fetchone()[0]

        # single-pass INSERT — no COUNT pre-flight. SQLite also stores ligand_id
        # (carried for free by the filter query's Results scan) so passing-ligand
        # counts avoid a JOIN back into the row-store Results table.
        if self._filtered_poses_has_ligand_id:
            insert_poses = f"""
                WITH results_poses AS (
                    {query})
                INSERT INTO Filtered_poses(filter_id, pose_id, ligand_id)
                SELECT {filter_id}, pose_id, ligand_id FROM results_poses;
                """
        else:
            insert_poses = f"""
                WITH results_poses AS (
                    {query})
                INSERT INTO Filtered_poses(filter_id, pose_id)
                SELECT {filter_id}, pose_id FROM results_poses;
                """
        self.conn.execute(insert_poses)

        # O(1) existence check — works for both SQLite and DuckDB
        row = self.conn.execute(
            "SELECT 1 FROM Filtered_poses WHERE filter_id = ? LIMIT 1", (filter_id,)
        ).fetchone()

        if row is None:
            # Nothing matched — clean up the empty filter entry; old bookmark untouched
            self.conn.execute(
                "DELETE FROM Filters WHERE filter_id = ?", (filter_id,)
            )
            self.conn.commit()
            return 0

        # Results exist — now safe to remove the old bookmark atomically
        if old_filter_id is not None:
            self.conn.execute(
                "DELETE FROM Filtered_poses WHERE filter_id = ?", (old_filter_id,)
            )
            self.conn.execute(
                "DELETE FROM Filters WHERE filter_id = ?", (old_filter_id,)
            )

        self.conn.commit()
        return 1

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

    def _stream_query(self, query: str, batch_size: int = 1000):
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

    def _execute_to_df(self, sql: str) -> pd.DataFrame:
        """
        Helper function for compatibility for eg sqlite queries that have a different
        implementation than duckdb of pandas

        Args:
            sql (str): query

        Returns:
            pd.DataFrame: data in data frame
        """
        return pd.read_sql_query(sql, self.conn)

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
        """
        Begin a transaction
        """
        raise NotImplementedError

    def _rollback(self):
        """
        Roll back transaction
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
        """
        raise NotImplementedError

    def _delete_old_duplicate_results(self):
        """Checks if a pose is uniquely represented in the Results table,
        and deletes it from Results if duplicated.
        Based on the following columns:
        ligname,
        receptor,
        pose_coordinates,
        flexible_res_coordinates
        """
        raise NotImplementedError

    def _insert_interactions(self, interactions: list[tuple]):
        """
        Inserts interaction tuples with pose_id s (IMPORTANT) into the interaction table.

        Args:
            interactions (list[tuple]): _description_
        """
        raise NotImplementedError

    def _update_interaction_counts(self, data: list[dict]):
        """
        Data is a dict that is expected to contain pose_id, num_int, and num_hb as keys.

        Args:
            data (list[dict]): _description_
        """

    def _insert_completed_poses(self, pose_ids: list[tuple], tracking_table: str):
        """
        Inserts processed poses into process tracking table

        Args:
            pose_ids (list[int]): _description_
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

    @staticmethod
    def _interaction_index_fields(r) -> tuple:
        """Extract the 6 interaction-description fields from either a dict or an InteractionRecord."""
        if isinstance(r, dict):
            return (
                r["type"],
                r["chain"],
                r["residue"],
                r["resid"],
                r["recname"],
                r["recid"],
            )
        return (r.type, r.chain, r.residue, r.resid, r.recname, r.recid)

    def _insert_interaction_index_rows(self, interactions: list):
        """
        Writes unique interactions to database

        Args:
            interactions (list): list of InteractionRecord or dicts with interaction description fields
        """
        raise NotImplementedError

    def _delete_nontables(self):
        """
        Deletes objects in the database that are not tables
        """
        pass

    def _calc_percentile_cutoff(self, percentile: float, column="docking_score"):
        """Make query for percentile by calculating energy or leff cutoff

        Args:
            percentile (float): cutoff percentile
            column (str, optional): string indicating column for percentile to be calculated over

        Returns:
            float: effective cutoff value of results based on percentile
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

    def _delete_cluster(self, cluster_id: int):
        """
        Deletes all data associated with a cluster: Pose_clusters,
        Cluster_groups, and Clusters rows for the given cluster_id.

        Args:
            cluster_id (int): cluster id to delete
        """
        raise NotImplementedError

    def _format_output_fields(
        self, outfields: Union[str, list], results_alias="R", ligands_alias="L"
    ) -> list:
        """Handles string or list input of column names to be outputted, will make sure ligname
        is in the list, and make sure all options are valid

        Returns:
            list: column names formatted with table aliases

        Raises:
            OptionError
        """
        if isinstance(outfields, str):
            outfields = outfields.replace(" ", "")
            outfields_list = outfields.split(",")
        elif isinstance(outfields, (list, tuple)):
            outfields_list = list(outfields)
        else:
            logger.warning(
                "The provided outfields is not in a usable format (string or list). Will only use ligname"
            )
            outfields_list = []
        if "ligname" not in [field.lower() for field in outfields_list]:
            outfields_list.insert(0, "ligname")
        possible_columns, table_formatted_columns = self._get_possible_output_columns()
        table_formatted_outfields = []
        for outfield in outfields_list:
            if outfield.lower() in possible_columns:
                table_formatted_outfields.append(
                    table_formatted_columns[possible_columns.index(outfield.lower())]
                )
            else:
                logger.warning(
                    f"{outfield} is not a valid output option, and will be removed from the output columns. Please see rt_process_vs.py --help for allowed options"
                )
        return [
            outfield.format(Ligands_alias=ligands_alias, Results_alias=results_alias)
            for outfield in table_formatted_outfields
        ]

    def _attach_db(self, new_db: str, new_db_alias: str = "attached_db") -> str:
        """Attaches new database file to current database

        Args:
            new_db (str): file name for database to attach
            new_db_name (str): name of new database

        """
        raise NotImplementedError

    def _db_compatible_for_merge(self, merging_db_alias: str) -> bool:
        """
        Method that checks if the database merging into main is compatible with main,
        and checks if both databases are of appropriately high version where merge has
        been implemented
        #NOTE the ringtail duckdb currently does not track schema version

        Args:
            merging_db_alias (str): alias for the database being merged into main

        Returns:
            bool: if the two databases are compatible
        """
        raise NotImplementedError

    def _get_merge_id(self, mergingdb_path: str) -> int:
        """
        Gets the merge id for the databse in given path

        Args:
            mergingdb_path (str): path to merging (secondary) database

        Returns:
            int: merge id returend by database
        """
        raise NotImplementedError

    def _merging_receptors_compatible(self) -> str:
        """
        Checks if the receptor names in the two databases are a mathc

        Returns:
            str: returns True or False string
        """
        raise NotImplementedError

    def _merge_db_properties_table(self, merge_id: int):
        """
        Merges database properties table, but importantly will not check for property compatibility

        Args:
            merge_id (int): merge session id

        Raises:
            MergeError
        """
        raise NotImplementedError

    def _merge_ligands_and_results_tables(self, merge_id: int):
        """
        Merges first the Ligands table, then the Results table, maintaining ligand_id and pose_id as primary keys,
        where their relationship to the original PK is kept track of in the mering datble. Duplicate ligands will maintain
        the ligand_id from the main database.

        Args:
            merge_id (int): merge session id

        Raises:
            MergeError
        """
        raise NotImplementedError

    def _merge_interaction_tables(self, merge_id: int):
        """
        Merges the interaction tables. Interaction definitions are unique and independent of the Results table, so we only
        insert those that are new with updated PK, and assign existing interaction_ids to those already described in primary db

        Args:
            merge_id (int): merge session id

        Raises:
            MergeError
        """
        raise NotImplementedError

    def _sync_auto_increment_state(self):
        """Reset any auto incremented sequences or counts according to last successful merge"""
        raise NotImplementedError

    def _rollback_merge(self, merge_id: int):
        """
        Remove all data inserted by a specific merge session

        Args:
            merge_id (int): Merge session for which to delete associated data
        """
        raise NotImplementedError

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
        raise NotImplementedError

    def _crossref_bookmark_builder(
        self, ligand_list: list[str], store_best_pose: bool
    ) -> str:
        """
        creates formatted sql for building the crossreferencing bookmark

        Args:
            ligand_list (list[str]): list of ligand names to be considered
            store_best_pose (bool): whether to get best pose or all poses

        Returns:
            str: formatted sql
        """
        raise NotImplementedError

    @abstractmethod
    def _remove_screening_tables(self): ...

    @abstractmethod
    def _create_screening_tables(self): ...

    # endregion
