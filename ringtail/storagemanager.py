#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail storage adaptor
#

import time
import json
import os.path
import re
from abc import ABC, abstractmethod
import numpy as np
import pandas as pd
from .logutils import get_logger

logger = get_logger(__name__)

# how to upgrade a database this version of Ringtail cannot open
UPGRADE_HINT = (
    "Upgrade Ringtail 1.x/2.x databases with rt_upgrade_db. Pre-release 3.0 databases "
    "need rt_upgrade_from_v3alpha.py, ask your package administrator for it."
)

from .util import db_alias_from_path, numlist2str
import sys
from signal import signal, SIGINT
from rdkit import Chem
from typing import Any, ClassVar
from ._version import __version__
from packaging.version import Version
from .filters import Filters
from .ringtailoptions import statuses
from .exceptions import (
    StorageError,
    OptionError,
    VersionError,
    MergeError,
    DatabaseQueryError,
)
from .clustermanager import *
from .querybuilder import QueryBuilder
from .types import PoseData
from .schema import (
    NUMERIC_TYPES,
    NO_FLEXRES_JSON,
    OUTFIELD_BY_TABLE,
    CANDIDATES_SUBQ,
    CANDIDATES_NAME,
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
    RECEPTORS_SCHEMA,
    STATUS_TABLE_SCHEMA,
    POSE_COMMENTS_SCHEMA,
    CUSTOM_POSE_DATA_SCHEMA,
    CUSTOM_LIGAND_DATA_SCHEMA,
    MINIMIZED_COLUMNS,
    ROW_ID_OWNERS,
    ROW_ID_TABLES,
    SCHEMA_VERSION_SCHEMA,
    SCHEMA_VERSION,
    build_create_table,
)
from collections import defaultdict


class StorageManager(ABC):
    """Base class for a generic virtual screening database object.
    This class holds some of the common API for StorageManager child classes.
    Each child class will implement their own functions to write to and read from the database

    The current schema version lives in schema.SCHEMA_VERSION and is stamped into
    the ringtail_schema_version table, see check_ringtaildb_version.

    Attributes:
        active_connection (bool): if there is a current open connection to database
        keyboard_interrupt_allowed (bool): if Ctrl+Z will work, for example not supported in a GUI
        QueryBuilder: dialect of querybuilder
    """

    QueryBuilder = QueryBuilder

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
            logger.debug(str(exc_value))
            raise exc_value
        return self

    def _sigint_handler(self, signal_received, frame):
        """Handles and reports if program is interrupted through the terminal"""
        logger.critical("Ctrl + C pressed, keyboard interupt initiated")
        self.__exit__(None, None, None)
        sys.exit(0)

    # endregion

    # region public api

    def get_pose_comment(self, pose_id: int) -> str | None:
        """
        Gets the comment for a pose.

        Args:
            pose_id (int): pose to get the comment for

        Returns:
            str | None: comment, or None if no comment or comment table exists
        """
        if not self._is_table("Pose_comments"):
            return None

        query = self.QueryBuilder()
        query.SELECT("comment").FROM("Pose_comments")
        query.WHERE("pose_id = ?", pose_id)
        row = self.db_query(*query.build()).fetchone()
        return row[0] if row else None

    def set_pose_comment(self, pose_id: int, comment: str) -> None:
        """
        Sets the comment for a pose, or removes it if empty.

        Args:
            pose_id (int): pose to update
            comment (str): comment to store
        """
        self.ensure_gui_tables()
        query = self.QueryBuilder()

        if comment:
            query.INSERT_INTO("Pose_comments", "pose_id", "comment")
            sql = (
                query.build()[0]
                + " ON CONFLICT(pose_id) DO UPDATE SET comment = excluded.comment"
            )
            self.db_query(sql, (pose_id, comment))
        else:
            query.DELETE_FROM("Pose_comments").WHERE("pose_id = ?", pose_id)
            self.db_query(*query.build())

        self.conn.commit()

    def create_transaction_tracking_table(self, table_name: str = "tracking_table"):
        """
        Creates a table to track pose ids to which are part of a multiple transaction
        database operation, which may occur over multiple commits and database connections.

        Args:
            table_name (str, optional): Name of the table. Defaults to "tracking_table".
        """
        if self._is_table(table_name):
            logger.info(
                "Interaction calculations have been started in a previous "
                "database connection, and not yet finished."
            )

        table_sql = f"""
            CREATE TABLE IF NOT EXISTS {self._safe_quoting(table_name)} (
                pose_id INTEGER,
                hb_cutoff DOUBLE,
                vdw_cutoff DOUBLE,
                FOREIGN KEY (pose_id) REFERENCES Results(pose_id)
            );
        """
        self.db_query(table_sql, commit=True)

    def pose_row_in_table(self, table: str, pose_id: int) -> int | None:
        """
        Gets the paging key of a pose in a table or selection.

        Args:
            table (str): table, bookmark, or Candidates
            pose_id (int): pose to locate

        Returns:
            int | None: paging key, or None if the pose is not in the selection
        """
        query = self.QueryBuilder()
        if self.is_bookmark(table):
            key = "pose_id" if self.dialect == "sqlite" else "rowid"
            query.SELECT(key).FROM("Filtered_poses").WHERE(
                "filter_id = (SELECT filter_id FROM Filters WHERE name = ?)",
                table,
            )
        elif self._is_candidates_table(table):
            query.SELECT("rowid").FROM("Results")
            self._apply_table_filter(query, table)
        else:
            query.SELECT("rowid").FROM(self._safe_quoting(table))

        query.WHERE("pose_id = ?", pose_id)
        row = self.db_query(*query.build()).fetchone()
        return row[0] if row else None

    def get_starting_rowid(self, table: str) -> int | None:
        """
        Gets the lowest paging key in a table or selection.

        Args:
            table (str): table, bookmark, or Candidates

        Returns:
            int | None: lowest paging key, or None if empty or not found
        """
        query = self.QueryBuilder()

        if self._is_table(table):
            query.SELECT("MIN(rowid)").FROM(self._safe_quoting(table))
        elif self.is_bookmark(table):
            key = "pose_id" if self.dialect == "sqlite" else "rowid"
            query.SELECT(f"MIN({key})").FROM("Filtered_poses").WHERE(
                "filter_id = (SELECT filter_id FROM Filters WHERE name = ?)",
                table,
            )
        elif self._is_candidates_table(table):
            query.SELECT("MIN(rowid)").FROM("Results")
            self._apply_table_filter(query, table)
        else:
            logger.error(f"Table -{table}- does not exist in the database.")
            return None

        return self.db_query(*query.build()).fetchone()[0]


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
        custom_data = docking_data.get("custom_data", None)
        if poses:
            # interactions has ligname and pose rank
            self._insert_docking_data(
                poses, interaction_data, custom_data, duplicate_handling
            )
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
            return
        # a receptor saved from its file has no box or flexible residue info yet
        _, box_dim, box_center, grid_spacing, flexres, flexres_atomnames = receptor_data
        self.db_query(
            f"""UPDATE Receptors SET
                box_dim = COALESCE(box_dim, ?),
                box_center = COALESCE(box_center, ?),
                grid_spacing = COALESCE(grid_spacing, ?),
                flexible_residues = CASE WHEN COALESCE(flexible_residues, '{NO_FLEXRES_JSON}') = '{NO_FLEXRES_JSON}'
                    THEN ? ELSE flexible_residues END,
                flexres_atomnames = CASE WHEN COALESCE(flexres_atomnames, '{NO_FLEXRES_JSON}') = '{NO_FLEXRES_JSON}'
                    THEN ? ELSE flexres_atomnames END
            WHERE receptor_id = 1;""",
            (box_dim, box_center, grid_spacing, flexres, flexres_atomnames),
            commit=True,
        )

    def filter_results(
        self,
        filters: Filters,
        output_bookmark: str,
        input_bookmark: str = None,
    ) -> int:
        """Generate and execute database queries from given filters.

        Args:
            filters (Filters): the filter specification. Flat filtering is a one-leaf tree
                and nested filtering a deeper one, so both render through the same path.
            output_bookmark (str): name of bookmark in which to save the data
            input_bookmark (str, optional): if filtering not across all data, but a pre-filtered bookmark

        Returns:
            int: number of passing ligands
        """
        filters = Filters.from_dict(filters)
        # before we do anything, check that the DB schema version is usable
        rt_version_same, db_rt_version = self.check_ringtaildb_version()
        if not rt_version_same:
            raise StorageError(
                f"Database schema version '{db_rt_version}' is not compatible with this "
                f"Ringtail (schema {SCHEMA_VERSION}). {UPGRADE_HINT}"
            )

        # get the final filter query, has a {selection} place holder
        filter_query: str = self._generate_filtering_query(
            filters, output_bookmark, input_bookmark
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
            filters=filters.to_dict(),
            input_bookmark=input_bookmark,
        )
        logger.debug(
            f"Time to filter results: {time.perf_counter() - time0:.2f} seconds"
        )
        count = self.get_passing_ligands_count(output_bookmark)

        return count

    def fetch_cluster_options(self, ligname: str) -> list[tuple]:
        """
        Gets clustering sessions that contain one or more poses of the ligand.

        Args:
            ligname (str): ligand name

        Returns:
            list[tuple]: list of (cluster_id, cluster_window, name)
        """
        membership = self.QueryBuilder()
        membership.SELECT("PC.cluster_id").FROM("Pose_clusters", "PC")
        membership.JOIN("Results", "R", "pose_id", "Pose_clusters")
        membership.JOIN("Ligands", "L", "ligand_id", "Results")
        membership.WHERE("L.ligname = ?", ligname)
        membership_sql, params = membership.build()

        query = self.QueryBuilder()
        query.SELECT(
            "C.cluster_id", "C.cluster_window", "C.name"
        ).FROM("Clusters", "C")
        query.WHERE(f"C.cluster_id IN ({membership_sql})", *params)
        query.ORDER_BY("C.cluster_id")

        return self.db_query(*query.build()).fetchall()

    def fetch_clustered_similars(
        self, ligname: str, cluster_id: int
    ) -> tuple[list, str, str]:
        """
        Gets ligands sharing a cluster group with the given ligand.

        Args:
            ligname (str): ligand to find similar ligands for
            cluster_id (int): clustering session to search within

        Returns:
            tuple[list, str, str]: (ligands, bookmark_name, cluster_name)

        Raises:
            OptionError: if the clustering session does not exist or the ligand
                has no poses in it
        """
        query = self.QueryBuilder()
        query.SELECT("name").FROM("Clusters").WHERE(
            "cluster_id = ?", cluster_id
        )
        row = self.db_query(*query.build()).fetchone()
        if row is None:
            raise OptionError(f"Clustering session {cluster_id} does not exist.")
        cluster_name = row[0].replace(".", "p")

        query = self.QueryBuilder()
        query.SELECT("DISTINCT PC.cluster_group").FROM("Pose_clusters", "PC")
        query.JOIN("Results", "R", "pose_id", "Pose_clusters")
        query.JOIN("Ligands", "L", "ligand_id", "Results")
        query.WHERE("PC.cluster_id = ?", cluster_id)
        query.WHERE("L.ligname = ?", ligname)
        groups = [row[0] for row in self.db_query(*query.build()).fetchall()]

        if not groups:
            raise OptionError(
                f"Ligand '{ligname}' has no poses in clustering session {cluster_id}."
            )

        placeholders = ",".join("?" for _ in groups)
        query = self.QueryBuilder()
        query.SELECT("DISTINCT L.ligname").FROM("Pose_clusters", "PC")
        query.JOIN("Results", "R", "pose_id", "Pose_clusters")
        query.JOIN("Ligands", "L", "ligand_id", "Results")
        query.WHERE("PC.cluster_id = ?", cluster_id)
        query.WHERE(f"PC.cluster_group IN ({placeholders})", *groups)
        ligands = self.db_query(*query.build()).fetchall()

        bookmark_name = f"similar_{ligname}_{cluster_name}"
        return ligands, bookmark_name, cluster_name

    def table_length(self, table: str, distinct_column: str = None) -> int:
        """
        Gets the number of rows or distinct values in a table or pose selection.

        Args:
            table (str): table name, bookmark, status table, or Candidates
            distinct_column (str, optional): column whose distinct values to count.
                Pose selections support queryable columns, including custom properties.

        Returns:
            int: number of rows, or distinct non-null values

        Raises:
            OptionError: if the selection or column is not valid
        """
        physical_table = (
            self._is_table(table)
            and table.lower() != "results"
            and not self._is_statustable(table)
            and not self._is_candidates_table(table)
        )

        field = None
        if physical_table:
            query = self.QueryBuilder().FROM(self._safe_quoting(table))
            if distinct_column is not None:
                columns = {
                    name.lower(): name
                    for name in self.fetch_dynamic_table_column_names(table)
                }
                column = columns.get(distinct_column.lower())
                if column is None:
                    raise OptionError(
                        f"Column '{distinct_column}' not found in table '{table}'."
                    )
                field = self._safe_quoting(column)
        else:
            sources = {}
            if distinct_column is not None:
                column_to_table = self.get_queryable_columns()
                source = column_to_table.get(distinct_column.lower())
                if source is None:
                    raise OptionError(
                        f"Column '{distinct_column}' not found in any queryable table."
                    )
                sources[distinct_column] = source

            query, fields = self._prepare_column_query(sources, table)
            if distinct_column is not None:
                field = fields[distinct_column]

        query.SELECT(
            "COUNT(*)" if field is None else f"COUNT(DISTINCT {field})"
        )
        return self.db_query(*query.build()).fetchone()[0]
    
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
            db_name = db_name or "attached_db"
            self._attach_db(db_path, db_name)
        query = self.QueryBuilder()
        query.SELECT("name").FROM("Filters", db_name=db_name)
        cur = self.db_query(query.build()[0])
        bookmark_names = [row[0].lower() for row in cur.fetchall()]

        if db_path:
            self.detach_db(db_name)

        return bookmark_names

    def get_bookmark_dependents(self, bookmark_name: str) -> list[str]:
        """
        Get bookmarks that were produced by filtering over the given bookmark.

        The parent name is recorded in Filters.filter_window when filtering is run
        with an input_bookmark. It is provenance only: there is no foreign key, and
        a dependent's poses live in its own Filtered_poses rows. Deleting the parent
        therefore leaves dependents fully queryable, it only leaves their recorded
        lineage pointing at a bookmark that no longer exists.

        Args:
            bookmark_name (str): bookmark to find the dependents of

        Returns:
            list: of dependent bookmark names
        """
        query = self.QueryBuilder()
        query.SELECT("name").FROM("Filters")
        query.WHERE("LOWER(filter_window) = ?", bookmark_name.lower())

        rows = self.db_query(*query.build()).fetchall()
        return [row[0].lower() for row in rows]

    def fetch_bookmarks_with_interaction_filters(self) -> list[str]:
        """
        Bookmarks whose stored filters select on interactions.

        Recalculating interactions changes what those filters would now select, so
        these are the bookmarks worth revisiting afterwards. Provenance only, like
        get_bookmark_dependents: their poses stay exactly where they are.

        Filters.filters holds the filter set as JSON. It is searched recursively
        rather than by key lookup because a tiered filter stores its criteria nested
        inside a boolean tree, and a flat lookup would report those bookmarks as
        interaction-free.

        Also counted: interaction-fingerprint clusterings, and every bookmark derived
        from an affected one (filtered or clustered with it as input_bookmark), since
        its selection was made from the old interactions too.

        Returns:
            list: names of bookmarks filtered on any interaction criterion
        """
        _interaction_filters = (
            "vdw_interactions",
            "hb_interactions",
            "reactive_interactions",
            "react_any",
        )

        def uses_interactions(node) -> bool:
            if isinstance(node, dict):
                for key, value in node.items():
                    if key == "hb_count" and value is not None:
                        return True
                    if key in _interaction_filters and value:
                        return True
                    if key == "cluster_type" and str(value).lower() == "ifp":
                        return True
                    if uses_interactions(value):
                        return True
            elif isinstance(node, list):
                return any(uses_interactions(item) for item in node)
            return False

        query = self.QueryBuilder()
        query.SELECT("name", "filters", "filter_window").FROM("Filters")
        rows = self.db_query(*query.build()).fetchall()
        affected = set()
        for name, filters, _ in rows:
            if not filters:
                continue
            try:
                decoded = json.loads(filters)
            except (ValueError, TypeError):
                # A filter set we cannot read is not a reason to fail the caller,
                # which is only assembling a warning.
                logger.debug(f"Could not parse stored filters for bookmark {name}.")
                continue
            if uses_interactions(decoded):
                affected.add(name.lower())

        # follow the input_bookmark lineage down to everything derived from them
        grew = True
        while grew:
            grew = False
            for name, _, parent in rows:
                if (
                    name.lower() not in affected
                    and parent
                    and parent.lower() in affected
                ):
                    affected.add(name.lower())
                    grew = True
        return [name.lower() for name, _, _ in rows if name.lower() in affected]

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
            if self.dialect == "sqlite":
                # SQLite: ligand_id is denormalized into Filtered_poses
                return self.conn.execute(
                    "SELECT COUNT(DISTINCT ligand_id) FROM Filtered_poses WHERE filter_id = ?",
                    (filter_id,),
                ).fetchone()[0]
            # DuckDB: columnar JOIN to Results reads only the ligand_id column
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

    def has_filter_data(self) -> bool:
        """Whether the database holds any bookmarks/filters or clusterings.

        Returns:
            bool: True if Filters or Clusters has rows
        """
        return any(
            self._is_table(table)
            and self.db_query(f"SELECT COUNT(*) FROM {table}").fetchone()[0] > 0
            for table in (FILTERS_SCHEMA.name, CLUSTERS_SCHEMA.name)
        )

    def delete_all_filter_data(self) -> None:
        """Delete every bookmark, filter and clustering, keeping the (empty) tables.

        Status tables and pose comments are screening decisions, not filter output,
        and are kept.
        """
        # children before parents, for the backends that enforce foreign keys
        for table in (
            POSE_CLUSTERS_SCHEMA.name,
            CLUSTER_GROUPS_SCHEMA.name,
            CLUSTERS_SCHEMA.name,
            FILTERED_POSES_SCHEMA.name,
            FILTERS_SCHEMA.name,
        ):
            if self._is_table(table):
                self.db_query(f"DELETE FROM {table}")
        self.conn.commit()
        logger.warning("All bookmarks, filters and clusterings have been deleted.")

    def rename_bookmark(self, old_name: str, new_name: str) -> None:
        """Rename a bookmark, and repoint everything that refers to it by name.

        Args:
            old_name (str): current bookmark name
            new_name (str): desired bookmark name
        """
        old_name = old_name.lower()
        new_name = new_name.lower()
        if old_name == new_name:
            return
        if self.is_bookmark(new_name):
            logger.warning(
                f"The bookmark {new_name} already exists, and will be overwritten."
            )
            self.delete_bookmark(new_name)
        self.conn.execute(
            "UPDATE Filters SET name = ? WHERE LOWER(name) = ?", (new_name, old_name)
        )
        # dependents record their parent by name, not by filter_id
        self.conn.execute(
            "UPDATE Filters SET filter_window = ? WHERE LOWER(filter_window) = ?",
            (new_name, old_name),
        )
        # Clusters.name is the standardized cluster descriptor, not a bookmark;
        # only cluster_window names a bookmark
        if "clusters" in self.tables_in_db():
            self.conn.execute(
                "UPDATE Clusters SET cluster_window = ? WHERE LOWER(cluster_window) = ?",
                (new_name, old_name),
            )
        self.conn.commit()
        logger.info(f"Bookmark '{old_name}' renamed to '{new_name}'.")

    def get_maxmiss_union(
        self,
        total_combinations: int,
        bookmark_name: str,
        all_filters=None,
        input_bookmark: str = None,
    ) -> tuple[int, str]:
        """
        Get results that are in union considering max miss

        Args:
            total_combinations (int): numer of possible combinations
            bookmark_name (str): name of bookmark to store
            all_filters (dict, optional): All filters used. Defaults to None.
            input_bookmark (str, optional): bookmark the combinations were filtered over. Defaults to None.

        Returns:
            int: num passing ligands
            str: union bookmark name, None if nothing passed
        """
        all_filters = all_filters or {}
        enumerated_bookmarks = []
        existing_bookmarks = self.get_all_bookmark_names()
        for i in range(total_combinations):
            bmn = bookmark_name + "_" + str(i)
            if bmn in existing_bookmarks:
                enumerated_bookmarks.append(f"'{bmn}'")
        if not enumerated_bookmarks:
            return 0, None

        subq = self.QueryBuilder()
        subq.SELECT("filter_id").FROM("Filters").WHERE(
            f"name IN ({', '.join(enumerated_bookmarks)})"
        )
        union_cols = (
            "DISTINCT pose_id, ligand_id"
            if self.dialect == "sqlite"
            else "DISTINCT pose_id"
        )
        query = self.QueryBuilder()
        query.SELECT(union_cols).FROM("filtered_poses").WHERE(
            f"filter_id IN ({subq.build()[0]})"
        )
        bookmark_name = f"{bookmark_name}_union"
        logger.debug("Running interaction union query")
        self._populate_filter_tables(
            name=bookmark_name,
            query=query.build()[0],
            filters=all_filters,
            input_bookmark=input_bookmark,
        )

        count = self.get_passing_ligands_count(bookmark_name)

        if not count:
            bookmark_name = None
        return count, bookmark_name

    def crossref_databases(
        self,
        wanted_dbs: list[tuple[str, str | None]] = None,
        unwanted_dbs: list[tuple[str, str | None]] = None,
        bookmark_prefix: str = "crossref",
        alternative_database_names: dict = None,
    ) -> tuple[int, dict, dict]:
        """
        Method to cross reference two or more databases. Will attach all other databases
        to the "current" database, which is always the first wanted database (will have the
        active ringtail object).

        Will create an intersect of ligand names in the wanted databases+bookmarks,
        and delete any ligands found in the union of the unwanted databases+bookmarks.

        A bookmark with specified prefix to the original bookmark name will be stored in each wanted database,
        making the crossreferencing data easily available for later use.

        The second element of each (database, scope) tuple may be a bookmark name
        OR a status table name ("accepted"/"maybe"/"rejected", any case), letting
        the intersect/exclude be scoped by acceptance status as well as by bookmark.

        Args:
            wanted_dbs (list[tuple[str, str | None]], optional): (database_path, scope) tuples, where scope is a bookmark, a status table, or None/"Results" for all results. Defaults to None.
            unwanted_dbs (list[tuple[str, str | None]], optional): (database_path, scope) tuples to exclude, scoped the same way. Defaults to None.
            bookmark_prefix (str, optional): prefix of the bookmark written in each wanted database. Defaults to "crossref".
            alternative_database_names (dict, optional):  {path: alt name}. Defaults to None.

        Raises:
            OptionError: if a scope is not a bookmark or status table of its database

        Returns:
            int: number of passing crossref ligands
            dict: {db_file: new_bookmark_name,}
            dict: "filter" {"wanted":[list of wanted dbs],"unwanted":[list of unwanted dbs]
        """
        alternative_database_names = alternative_database_names or {}
        # Aliases must be unique per attached database, the file stem alone is not
        used_aliases = {"main", db_alias_from_path(self.db_file).lower()}
        used_aliases.update(
            name.lower() for name in alternative_database_names.values()
        )

        def _unique_alias(path: str) -> str:
            base = f"xref_{db_alias_from_path(path)}"
            alias, n = base, 1
            while alias.lower() in used_aliases:
                n += 1
                alias = f"{base}_{n}"
            used_aliases.add(alias.lower())
            return alias

        # one alias per database, the first wanted database is this one
        aliases = {}
        entries = []
        for kind, databases in (
            ("wanted", wanted_dbs),
            ("unwanted", unwanted_dbs or []),
        ):
            for path, scope in databases:
                if not aliases:
                    aliases[path] = "main"
                elif path not in aliases:
                    aliases[path] = alternative_database_names.get(
                        path
                    ) or _unique_alias(path)
                entries.append((kind, aliases[path], path, scope))

        attached = []
        try:
            for path, alias in aliases.items():
                if alias != "main":
                    self._attach_db(path, alias)
                    attached.append(alias)
            for _, alias, path, scope in entries:
                if (
                    scope is None
                    or scope.lower() == "results"
                    or self._is_statustable(scope)
                    or self._is_candidates_table(scope)
                ):
                    continue
                bookmarks = self.get_all_bookmark_names(
                    None if alias == "main" else alias
                )
                if scope.lower() not in bookmarks:
                    raise OptionError(
                        f"{path or self.db_file} has no bookmark or status table named {scope!r}."
                    )

            # one common table expression of ligand names per (database, scope)
            selection_query = """
            {name} AS
                (SELECT ligname FROM {db_name}.Ligands
                INNER JOIN {db_name}.Results ON
                    {db_name}.Results.ligand_id={db_name}.Ligands.ligand_id
                    WHERE pose_id IN ({scope_poses}))"""
            names = [f"set_{i}" for i in range(len(entries))]
            ctes = [
                selection_query.format(
                    name=name,
                    db_name=alias,
                    scope_poses=self._get_scope_poses_query(scope, alias),
                )
                for name, (_, alias, _, scope) in zip(names, entries)
            ]
            wanted = [
                name for name, entry in zip(names, entries) if entry[0] == "wanted"
            ]
            unwanted = [
                name for name, entry in zip(names, entries) if entry[0] == "unwanted"
            ]
            wanted_joins = " ".join(
                f"INNER JOIN {name} ON {wanted[0]}.ligname = {name}.ligname"
                for name in wanted[1:]
            )
            unwanted_cte = ""
            unwanted_where = ""
            if unwanted:
                exclude_union = " UNION ".join(
                    f"SELECT ligname FROM {name}" for name in unwanted
                )
                unwanted_cte = f", unwanted AS ({exclude_union})"
                unwanted_where = "WHERE i.ligname NOT IN (SELECT ligname FROM unwanted)"
            query = f"""
            WITH {", ".join(ctes)},
            wanted AS (
                SELECT {wanted[0]}.ligname
                FROM {wanted[0]} {wanted_joins}
            ){unwanted_cte}
            SELECT DISTINCT i.ligname
            FROM wanted i
            {unwanted_where}
            """
            approved_ligand_names = [lig[0] for lig in self.db_query(query).fetchall()]
        finally:
            for alias in attached:
                self.detach_db(alias)

        # don't proceed if no approved ligands
        if not approved_ligand_names:
            return 0, {}, {}

        # write bookmarks in each of the databases
        filter_dict = {"wanted": wanted_dbs, "unwanted": unwanted_dbs}
        bookmark_query = self._crossref_bookmark_builder(approved_ligand_names)
        dbs_new_bookmark_names = {}
        for kind, alias, path, scope in entries:
            # unwanted databases hold none of the selected ligands
            if kind == "unwanted":
                continue
            # lowercase the suffix so the stored bookmark name stays SQLite-safe
            # and re-queryable when the scope is a (capitalized) status table.
            new_bookmark_name = f"{bookmark_prefix}_{(scope or 'results').lower()}"
            if alias == "main":
                self._populate_filter_tables(
                    name=new_bookmark_name,
                    query=bookmark_query,
                    filters=filter_dict,
                    input_bookmark=scope,
                )
            else:
                with type(self)(path) as db:
                    db._populate_filter_tables(
                        name=new_bookmark_name,
                        query=bookmark_query,
                        filters=filter_dict,
                        input_bookmark=scope,
                    )
            dbs_new_bookmark_names[path] = new_bookmark_name

        return (len(approved_ligand_names), dbs_new_bookmark_names, filter_dict)

    def cluster_data(
        self,
        bookmark_name: str | None,
        cluster_type: str = "mfp",
        cutoff: float = 0.5,
    ) -> tuple[str, int]:
        """
        Clusters data in a given bookmark. Will create a new bookmark named
        <bookmark_name>_<cluster_type>_<cutoff> containing the representative poses
        for the clusters

        Args:
            bookmark_name (str | None): bookmark name with poses to cluster, None clusters all results
            cluster_type (str, optional): "mfp" or "ifp". Defaults to "mfp".
            cutoff (float, optional): cutoff cluster distance. Defaults to 0.5.

        Returns:
            tuple (str, int): clustered bookmark name, number of clusters

        Raises:
            OptionError: if cluster_type is not "mfp" or "ifp"
        """
        if cluster_type.lower() not in ("mfp", "ifp"):
            raise OptionError(
                f"cluster_type must be 'mfp' or 'ifp', not {cluster_type!r}."
            )
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
        clustered_cols = (
            ["pose_id", "ligand_id"] if self.dialect == "sqlite" else ["pose_id"]
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
    
    def fetch_summary_data(
        self,
        columns: list[str] = ["docking_score", "leff"],
        percentiles: list[int] = [1, 10],
    ) -> dict:
        """
        Will collect summary data for the database: number of ligands, poses, unique interactions and interacting residues,
        and min, max and percentile cutoffs for each numeric column in columns, including custom columns

        Args:
            columns (list[str]): columns to summarize, non-numeric columns are skipped with a warning
            percentiles (list[int]): percentiles to calculate cutoffs for

        Raises:
            StorageError: if the database is empty

        Returns:
            dict: summary data
        """
        if self.db_empty():
            raise StorageError("There is no data in the database.")

        query = """
        SELECT
            (SELECT COUNT(ligand_id) FROM Ligands) AS num_ligands,
            (SELECT COUNT(pose_id) FROM Results) AS num_poses,
            (SELECT COUNT(interaction_id) FROM Interaction_indices) AS num_unique_interactions,
            (SELECT COUNT(*)
            FROM (SELECT DISTINCT interaction_type, rec_resid, rec_chain
                    FROM Interaction_indices) residues
            ) AS num_interacting_residues;"""
        data = self.db_query(query).fetchone()
        summary_data = {
            "num_ligands": data[0],
            "num_poses": data[1],
            "num_unique_interactions": data[2],
            "num_interacting_residues": data[3],
        }

        ranges = self.get_column_ranges(columns, None)
        for col in columns:
            if col not in ranges:
                logger.warning(
                    f"Requested column {col} is not a numeric Results or custom column and will not be used for the summary."
                )
                continue
            summary_data[f"min_{col}"], summary_data[f"max_{col}"] = ranges[col]
            for p in percentiles:
                summary_data[f"{p}%_{col}"] = self._calc_percentile_cutoff(p, col)

        return summary_data

    def merge_database(self, merging_db: str):
        """
        Method that merges two databases, ensuring integrity of primary and foreign keys.
        The merging will create a new table if needed, that keeps track of the primary key
        in the original and the merged database on a per-table basis. Another table will also
        keep track of how many databases has been merged into the primary database.
        The merging checks that the two databases have the same receptor and a
        compatible schema version. Databases produced from different docking formats
        may be merged when both checks pass.

        Args:
            merging_db (str): path to database being merged into current

        Raises:
            MergeError
        """
        # attach incoming database and check compatibility
        merging_db_alias = self._attach_db(merging_db, "merging")
        if not merging_db_alias:
            raise StorageError(
                f"{type(self).__name__}._attach_db did not return the alias it attached "
                f"{merging_db} as, so the database could not be detached again."
            )
        try:
            if not self._db_compatible_for_merge(merging_db_alias):
                raise MergeError(
                    "Trying to merge two databases of incompatible or too old versions, cannot proceed."
                )
            if not self._merging_receptors_compatible():
                raise MergeError(
                    f"The receptors in the merging databases are not the same. \nThese databases cannot be merged."
                )
            logger.info(
                "The two databases have compatible receptors. Merging will proceed."
            )
            if not self.db_query(
                f"SELECT COUNT(*) FROM {merging_db_alias}.Results"
            ).fetchone()[0]:
                raise MergeError(f"{merging_db} has no results to merge.")

            # add to merging table the absolute path
            mergedb_abspath = str(os.path.abspath(merging_db))
            merge_id = self._get_merge_id(mergedb_abspath)
            # ids above these were added by this merge, the rest were already here
            last_ids = self.last_row_ids()

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
                self._rollback_merge(merge_id, last_ids)
                raise MergeError(f"Merging failed and was rolled back: {e}") from e
            else:
                self._sync_auto_increment_state()
                logger.info(
                    f"The database {merging_db} has been successfully merged into {self.db_file}."
                )
        finally:
            self._cleanup_storage(merging_db_alias, vacuum=True, reindex=True)
            logger.info("The final database has been cleaned up, and indices rebuilt.")

    def last_row_ids(self) -> dict:
        """Largest id of each growing table, so delete_rows_after can undo a failed write.

        Returns:
            dict: _description_"""
        return {
            column: self._max_id(table, column) if self._is_table(table) else 0
            for column, table in ROW_ID_OWNERS.items()
        }

    def delete_rows_after(self, last_ids: dict):
        """Rolls back, then deletes every row added after last_row_ids.

        Args:
            last_ids (dict): output of last_row_ids from before the write
        """
        self._rollback()
        for table, column in ROW_ID_TABLES:
            if self._is_table(table):
                self.db_query(
                    f"DELETE FROM {table} WHERE {column} > ?", (last_ids[column],)
                )
        self.conn.commit()

    def complete_merging(self):
        """
        Completes the merging process by creating empty filter, status and comment tables, as well
        as resetting any auto incremented values for safety
        """
        # Rebuild indices dropped in prepare_for_merging()
        self._create_indices()
        self._create_screening_tables()
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
        # a bookmark, status table, candidates or all results, checked before the new file exists
        scope_poses = self._get_scope_poses_query(bookmark_name, "main")
        # create the database
        new_db = self.__class__(database_name)
        # make tables
        has_comments = self._is_table(POSE_COMMENTS_SCHEMA.name)
        with new_db:
            new_db._create_tables(commit=False)
            if has_comments:
                # Pose_comments is created lazily
                new_db.ensure_gui_tables()
        logger.debug(f"Created a new database with empty tables: {database_name}.")
        # attach the incoming database
        alias = "subset"
        self._attach_db(database_name, alias)
        logger.debug(
            f"Database {database_name} attached to main database {self.db_file}."
        )

        # columns by name, since an upgraded database can hold them in another order
        def cols(schema) -> str:
            return ", ".join(
                name
                for name, col in schema.columns.items()
                if not (col.sqlite_only and self.dialect != "sqlite")
            )

        # ligands
        ligands = f"""
        INSERT INTO {alias}.Ligands ({cols(LIGANDS_SCHEMA)})
        SELECT {cols(LIGANDS_SCHEMA)} FROM main.Ligands
        WHERE main.Ligands.ligand_id IN (
            SELECT ligand_id from main.Results
            WHERE pose_id IN ({scope_poses}));"""
        self.db_query(ligands)
        logger.debug("Ligands have been copied into the new subset database.")
        # receptor
        receptor = f"""
        INSERT INTO {alias}.Receptors ({cols(RECEPTORS_SCHEMA)})
        SELECT {cols(RECEPTORS_SCHEMA)} FROM main.Receptors;
        """
        self.db_query(receptor)
        logger.debug("The receptor have been copied into the new subset database.")
        # results
        poses = f"""
        INSERT INTO {alias}.Results ({cols(RESULTS_SCHEMA)})
        SELECT {cols(RESULTS_SCHEMA)} FROM main.Results
        WHERE main.Results.pose_id IN ({scope_poses});"""
        self.db_query(poses)
        logger.debug("Results have been copied into the new subset database.")
        # interaction_indices
        interaction_indices = f"""
        INSERT INTO {alias}.Interaction_indices ({cols(INTERACTION_INDICES_SCHEMA)})
        SELECT {cols(INTERACTION_INDICES_SCHEMA)} FROM main.Interaction_indices
        WHERE main.Interaction_indices.interaction_id IN (
            SELECT interaction_id FROM main.Interactions
            WHERE pose_id IN ({scope_poses}));"""
        self.db_query(interaction_indices)
        # interactions
        interactions = f"""
        INSERT INTO {alias}.Interactions ({cols(INTERACTIONS_SCHEMA)})
        SELECT {cols(INTERACTIONS_SCHEMA)} FROM main.Interactions
        WHERE main.Interactions.pose_id IN ({scope_poses});"""
        self.db_query(interactions)
        logger.debug("Interactions have been copied into the new subset database.")
        # pose statuses and comments — hand curation from visual inspection, and the
        # most expensive data in the database to recreate, so it travels with the subset
        for status_table in ("Accepted", "Maybe", "Rejected"):
            if not self._is_table(status_table):
                continue  # older database predating the status tables
            self.db_query(f"""
            INSERT INTO {alias}.{status_table} ({cols(STATUS_TABLE_SCHEMA)})
            SELECT {cols(STATUS_TABLE_SCHEMA)} FROM main.{status_table}
            WHERE main.{status_table}.pose_id IN ({scope_poses});""")
        logger.debug("Pose statuses have been copied into the new subset database.")
        if has_comments:
            comments = f"""
            INSERT INTO {alias}.{POSE_COMMENTS_SCHEMA.name} ({cols(POSE_COMMENTS_SCHEMA)})
            SELECT {cols(POSE_COMMENTS_SCHEMA)} FROM main.{POSE_COMMENTS_SCHEMA.name}
            WHERE main.{POSE_COMMENTS_SCHEMA.name}.pose_id IN ({scope_poses});"""
            self.db_query(comments)
            logger.debug("Pose comments have been copied into the new subset database.")
        # receptor
        dbprop = f"""
        INSERT INTO {alias}.db_properties ({cols(DB_PROPERTIES_SCHEMA)})
        SELECT {cols(DB_PROPERTIES_SCHEMA)} FROM main.db_properties;
        """
        self.db_query(dbprop)
        logger.debug("The db_properties have been copied into the new subset database.")
        try:
            self.conn.commit()
        except Exception as e:
            raise StorageError(f"Problems while creating a subset database: {e}") from e

        self.detach_db(alias)
        # rows were copied with their ids, so id sequences must catch up
        with new_db:
            new_db._sync_auto_increment_state()
        logger.info(f"Subset database {database_name} has been successfully created.")

    def fetch_pose_interactions(self, pose_ids) -> list | dict:
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
        selection: list | str,
        group_by: str = None,
        order_results: str = None,
    ) -> str:
        """
        Generates query to gather chosen columns based on passing poses in a bookmark

        Args:
            bookmark_name (str): bookmark name from which to get the passing poses
            selection (list | str): what columns to have in the output query
            group_by (str, optional): if set, gets the best passing pose per ligand
            order_results (str, optional): Whether or not to order by a column

        Raises:
            OptionError: if selection is '*' or the bookmark is not valid

        Returns:
            str: sql string that describes selection of data from bookmark
        """
        if selection == "*":
            raise OptionError(
                "Output fields/columns cannot be 'all'/'*', please select one or more specific columns, or use the default."
            )
        column_to_table = self.get_queryable_columns()
        outfields = self._format_output_fields(selection)
        sources = {col: column_to_table[col.lower()] for col in outfields}

        # Include the ordering column's table even if the column is not selected.
        order_column = None
        if order_results and order_results.lower() in column_to_table:
            order_column = order_results
            sources[order_column] = column_to_table[order_column.lower()]

        query, fields = self._prepare_column_query(sources, bookmark_name)
        query.SELECT(*(fields[col] for col in outfields))
        if group_by:
            query.WHERE(
                self._get_best_pose_per_ligand_condition(
                    bookmark_name, results_alias="Results"
                )
            )
        if order_column:
            query.ORDER_BY(fields[order_column])
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
            if (
                self.is_bookmark(requested_data)
                or self._is_statustable(requested_data)
                or self._is_candidates_table(requested_data)
            ):
                query.SELECT("*").FROM("Results")
                self._apply_table_filter(query, requested_data)
            else:
                query.SELECT("*").FROM(requested_data)
            return self._execute_to_df(query.build()[0])
        else:
            return self._execute_to_df(requested_data)

    # TODO this is a method i think may have some duplication floating around
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

    def get_previous_docking_mode(self) -> str | None:
        """
        Checks the docking_mode last used in a database write session

        Returns:
            str | None: docking_mode if any
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

    def get_column_ranges(
        self, columns: list[str], table: str
    ) -> dict[str, tuple[float | int | None, float | int | None]]:
        """
        Gets MIN and MAX of given column(s), including custom properties if they were inferred as numeric
        at time of database write

        Args:
            columns (list[str]): column names, non-numeric columns are left out
            table (str): bookmark, status table, candidates, or None/"Results" for all results

        Returns:
            dict[str, tuple[float | int | None, float | int | None]]: column:(min, max,)
        """
        query, fields = self._numeric_column_query(columns, table)
        if not fields:
            return {}
        query.SELECT(
            *(f"{agg}({field})" for field in fields.values() for agg in ("MIN", "MAX"))
        )
        row = self.db_query(*query.build()).fetchone()
        return {col: (row[2 * i], row[2 * i + 1]) for i, col in enumerate(fields)}

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
        query, fields = self._numeric_column_query([column], table)
        if not fields:
            raise OptionError(
                f"Requested column {column} is not a numeric Results column, percentiles cannot be calculated."
            )
        aggregate = "MIN" if column.lower() in MINIMIZED_COLUMNS else "MAX"
        query.SELECT(f"{aggregate}({fields[column]})").GROUP_BY("Results.ligand_id")

        values = [
            val[0]
            for val in self.db_query(query.build()[0]).fetchall()
            if val[0] is not None
        ]
        if not values:
            raise OptionError(
                f"No numeric values for {column!r} in selection {table!r}."
            )

        bins = np.linspace(0, 100, num_bins + 1)
        bin_edges = np.percentile(values, bins)
        return bins, bin_edges

    def fetch_histogram_counts(
        self, column: str, num_bins: int, table: str
    ) -> tuple[list[float], list[int]]:
        """
        Will count the poses in each of num_bins equal-width bins of a numeric column, including custom columns.

        Args:
            column (str): numeric Results or custom column
            num_bins (int): number of bins
            table (str): bookmark, status table, candidates, or None/"Results" for all results

        Raises:
            OptionError: if column is not numeric, or table is not a valid selection

        Returns:
            tuple[list[float], list[int]]: num_bins + 1 bin edges and num_bins counts, both empty if there are no values
        """
        values, fields = self._numeric_column_query([column], table)
        if not fields:
            raise OptionError(
                f"Requested column {column} is not a numeric Results column, cannot make a histogram."
            )
        field = fields[column]
        values.SELECT(f"{field} AS val").WHERE(f"{field} IS NOT NULL")
        values_sql, params = values.build()

        min_val, max_val, total = self.db_query(
            f"SELECT MIN(val), MAX(val), COUNT(*) FROM ({values_sql}) sub", params
        ).fetchone()
        if not total:
            return [], []
        min_val, max_val = float(min_val), float(max_val)
        if min_val == max_val:
            return [min_val - 0.5, min_val + 0.5], [total]
        bin_width = (max_val - min_val) / num_bins

        # the minimum gives exactly 0 in double precision, a float32 column on DuckDB would not
        position = f"(CAST(val AS DOUBLE) - {min_val}) / {bin_width}"
        # position is never negative, so SQLite's truncating CAST floors it, DuckDB's CAST rounds
        if self.dialect == "duckdb":
            raw_bin = f"CAST(FLOOR({position}) AS INTEGER)"
        else:
            raw_bin = f"CAST({position} AS INTEGER)"
        last_bin = num_bins - 1
        # raw SQL, QueryBuilderDuck.GROUP_BY would wrap the bin expression in ANY_VALUE
        # the maximum lands at num_bins, the CASE puts it in the last bin (closed right edge, like numpy)
        rows = self.db_query(
            f"""SELECT CASE WHEN {raw_bin} > {last_bin} THEN {last_bin} ELSE {raw_bin} END AS bin_idx,
                COUNT(*)
            FROM ({values_sql}) sub
            GROUP BY bin_idx""",
            params,
        ).fetchall()

        counts = [0] * num_bins
        for bin_idx, count in rows:
            counts[bin_idx] = count
        bin_edges = [min_val + i * bin_width for i in range(num_bins + 1)]
        return bin_edges, counts

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

    def fetch_receptor_object(self) -> dict | None:
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
            query = self.QueryBuilder()
            row = self.db_query(
                *query.SELECT("recname", "receptor_object").FROM("Receptors").build()
            ).fetchone()
            columns = ["recname", "receptor_object"]
        return dict(zip(columns, row)) if row else {}

    def fetch_flexres_info(self, receptor: str| int) -> tuple[list, list]:
        """fetch flexible residues names and atomname lists

        Always returns lists. The columns are decoded separately because databases from
        previous versions mix the spellings of "none" across them.

        Args:
            receptor (str | int): receptor descriptor, either receptor_id or receptor name

        Returns:
            tuple[list, list]: (flexible_residues, flexres_atomnames)
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
            return [], []
        residues = self._decode_flexres_column(info[0])
        atomnames = self._decode_flexres_column(info[1])
        if bool(residues) != bool(atomnames):
            logger.warning(
                f"Receptor {receptor} stores flexible residues without their atom "
                "names, or the reverse. No flexible residue mols can be built from it."
            )
        return residues, atomnames

    def fetch_custom_pose_properties(
        self, pose_ids: int | list[int]
    ) -> dict[int, dict[str, object]]:
        """
        Gets custom properties for the given poses, excluding missing values.

        Args:
            pose_ids (int | list[int]): one or more pose ids

        Returns:
            dict[int, dict[str, object]]: {pose_id: {property: value}}
        """
        if isinstance(pose_ids, int):
            pose_ids = [pose_ids]
        if not pose_ids or not self._is_table("custom_pose_data"):
            return {}

        columns = [
            name
            for name in self.fetch_dynamic_table_column_names("custom_pose_data")
            if name.lower() != "pose_id"
        ]
        if not columns:
            return {}

        placeholders = ",".join("?" for _ in pose_ids)
        query = self.QueryBuilder()
        query.SELECT(
            "pose_id",
            *(self._safe_quoting(col) for col in columns),
        ).FROM("custom_pose_data").WHERE(
            f"pose_id IN ({placeholders})", *pose_ids
        )

        rows = self.db_query(*query.build()).fetchall()
        return {
            row[0]: {
                col: value
                for col, value in zip(columns, row[1:])
                if value is not None
            }
            for row in rows
        }

    def fetch_rdkit_pose_properties(
        self,
        pose_ids,
        include_flexres: bool = False,
        include_comment: bool = False,
    ) -> list[PoseData]:
        """
        Gets molecular data needed to create rdkit mols for the given pose(s), including
        ligand mol binary and name via a JOIN with Ligands.

        Only the columns a caller asks for are fetched: the flexible-residue
        coordinate blob and the pose comment are both optional, so callers that
        don't need them (e.g. the GUI viewer) don't pay to read them.

        Args:
            pose_ids (int | list[int]): one or more pose ids
            include_flexres (bool): also fetch flexible_res_coordinates (needed for
                SDF/PDB export). When False, PoseData.flexres_coordinates is None.
            include_comment (bool): also LEFT JOIN Pose_comments for the user comment.
                When False (or the table is absent), PoseData.comment is ''.

        Returns:
            list[PoseData]
        """
        if isinstance(pose_ids, int):
            pose_ids = [pose_ids]
        # Pose_comments is a GUI-only table that may not exist in CLI/older databases
        has_comments = include_comment and self._is_table(POSE_COMMENTS_SCHEMA.name)
        placeholders = ",".join(["?"] * len(pose_ids))
        query = self.QueryBuilder()
        query.SELECT(
            "R.pose_id",
            "L.rdmol",
            "L.ligname",
            "R.docking_score",
            "R.leff",
            "R.pose_coordinates",
            "R.pose_rank",
        )
        if include_flexres:
            query.SELECT("R.flexible_res_coordinates")
        if has_comments:
            query.SELECT("COALESCE(Pc.comment, '') AS comment")
        query.FROM("Results", "R").JOIN("Ligands", "L", "ligand_id")
        if has_comments:
            query.JOIN("Pose_comments", "Pc", "pose_id", kind="LEFT")
        query.WHERE(f"R.pose_id IN ({placeholders})", *pose_ids)
        rows = self.db_query(*query.build()).fetchall()

        result = []
        for row in rows:
            # core columns (pose_coordinates stored natively -> list of [x,y,z])
            pose_id, rdmol, ligname, docking_score, leff, coords, pose_rank = row[:7]
            extra = row[7:]
            flexres = extra[0] if include_flexres else None
            comment = (
                (extra[1] if include_flexres else extra[0]) if has_comments else ""
            )
            result.append(
                PoseData(
                    pose_id,
                    rdmol,
                    ligname,
                    docking_score,
                    leff,
                    self._deserialize_pose_coordinates(coords),
                    pose_rank,
                    flexres,
                    comment,
                )
            )
        return result
    
    def fetch_columns_from_table_as_dicts(
        self, table: str, columns: list, length: int = None, starting_rowid: int = 0
    ) -> tuple[list[str], list[dict]]:
        """
        Gets requested columns as dictionaries from a table or a selection of poses.
        For Results, bookmarks, status tables and Candidates, joins tables as needed
        to include custom properties and other requested columns.

        Args:
            table (str): table name, bookmark name, or Candidates
            columns (list): column names to select. For pose selections, "*" includes
                Results and custom columns, and "status" includes pose status.
            length (int, optional): max rows to return. Defaults to None (all rows).
            starting_rowid (int, optional): rowid lower bound. Defaults to 0.

        Returns:
            tuple[list[str], list[dict]]: column names and list of row dicts
        """
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

        column_to_table = self.get_queryable_columns()
        include_status = any(col.lower() == "status" for col in columns)

        data_cols = []
        for col in columns:
            if col.lower() == "status":
                continue
            if col == "*":
                data_cols.extend(RESULTS_SCHEMA.columns)
                data_cols.extend(
                    name
                    for name, source in column_to_table.items()
                    if source == "custom_pose_data"
                )
            else:
                data_cols.append(col)

        sources = {}
        for col in data_cols:
            source = column_to_table.get(col.lower())
            if source is None:
                raise OptionError(f"Column '{col}' not found in any queryable table.")
            sources[col] = source

        query, fields = self._prepare_column_query(sources, table)
        query.SELECT(*fields.values())

        if include_status:
            if self._is_statustable(table):
                query.SELECT(f"'{table.lower()}' AS status")
            else:
                query.SELECT_STATUS()
        if length:
            query.LIMIT(length)
        if starting_rowid:
            query.WHERE(f"Results.rowid >= {starting_rowid}")
        if length or starting_rowid:
            query.ORDER_BY("Results.rowid")

        return self.get_query_data_as_dicts(query.build()[0])

    def fetch_pose_page(
        self, table: str, length: int, starting_rowid: int = 0, reverse=False
    ) -> dict[list[str], list]:
        """
        Gets a page of poses with status, comments, and a paging key.
        Custom columns are included after leff and before the comment.

        Args:
            table (str): Results, bookmark, status table, or Candidates
            length (int): maximum number of rows to return
            starting_rowid (int): paging key to start from, inclusive
            reverse (bool): whether to get rows in descending order

        Returns:
            dict: headers and data, with the paging key last in each row
        """

        if reverse:
            where_operator = "<="
        else:
            where_operator = ">="
        query = self.QueryBuilder()
        query.FROM("Results", "R")

        if table.lower() == "results":
            rowid = "R.rowid"

        elif self._is_statustable(table):
            query.JOIN(table, "T", "pose_id")
            rowid = "T.rowid"

        elif self._is_candidates_table(table):
            query.JOIN(CANDIDATES_SUBQ, "T", "pose_id")
            rowid = "R.rowid"

        elif self.is_bookmark(table):
            query.JOIN("filtered_poses", "fp", "pose_id").JOIN(
                "filters", "f", "filter_id", "filtered_poses"
            ).WHERE("f.name = ?", table)
            # SQLite's Filtered_poses is WITHOUT ROWID, pose_id is unique within a bookmark
            rowid = "fp.rowid" if self.dialect == "duckdb" else "fp.pose_id"

        else:
            raise ValueError(
                f"Table '{table}' is not a results table, status table, or bookmark."
            )

        if self._is_statustable(table):
            query.SELECT(f"'{table.lower()}' AS status")
        else:
            query.SELECT_STATUS()

        query.JOIN("Pose_comments", "Pc", "pose_id", kind="LEFT")
        query.JOIN("Ligands", "L", "ligand_id", "Results")

        column_to_table = self.get_queryable_columns()
        custom_columns = [
            name for name, source in column_to_table.items()
            if source == "custom_pose_data" and name not in {"status", "comment"}
        ]
        if custom_columns:
            query.JOIN("custom_pose_data", "CT", "pose_id", "Results", kind="LEFT")

        ordered_columns = [
            "ligname", "pose_id", "pose_rank", "docking_score", "leff",
            *custom_columns,
            "comment", "cluster_size", "cluster_rmsd", "num_hb", "receptor",
            "run_number", "delta", "num_interactions", "unbound_energy",
            "reference_RMSD", "energies_inter", "energies_vdw", "energies_electro",
            "energies_flexlig", "energies_flexlr", "energies_intra", "energies_torsional",
        ]
        selected_columns = set()
        for column in ordered_columns:
            name = column.lower()
            if name in selected_columns:
                continue
            selected_columns.add(name)
            if name == "comment":
                query.SELECT("COALESCE(Pc.comment, '') AS comment")
            else:
                source = query.aliased(column_to_table[name])
                query.SELECT(f"{source}.{self._safe_quoting(column)}")

        query.SELECT(rowid)
        query.WHERE(
            f"{rowid} {where_operator} ?", starting_rowid
        ).ORDER_BY(rowid).LIMIT(length).DESC(reverse)

        cursor = self.db_query(*query.build())
        headers = [desc[0] for desc in cursor.description]
        data = cursor.fetchall()
        return {"headers": headers, "data": data}

    def fetch_tracked_cutoffs(self, track_table_name: str):
        """
        The interaction cutoffs a partially finished run was using, if any.

        Every row carries the same pair, so one row answers it.

        Args:
            track_table_name (str): the transaction-tracking table.

        Returns:
            tuple[float, float] | None: (hb_cutoff, vdw_cutoff), or None if no pose
                has been processed yet.
        """
        query = self.QueryBuilder()
        query.SELECT("hb_cutoff", "vdw_cutoff").FROM(track_table_name).LIMIT(1)
        sql, params = query.build()
        row = self.db_query(sql, params).fetchone()
        if row is None or row[0] is None:
            return None
        return (float(row[0]), float(row[1]))

    def fetch_poses_pending_interactions(self, track_table_name: str, limit: int):
        """
        Fetch up to limit poses that still need interaction calculation: poses in
        Results whose pose_id is not yet recorded in the tracking table.

        Returns one batch and closes instead of streaming

        Ordered by ligand so a batch holds each ligand's poses together: the caller
        prepares one meeko molsetup per ligand and reuses it across that ligand's
        poses

        Args:
            track_table_name (str): name of the transaction-tracking table created
                by create_transaction_tracking_table.
            limit (int): maximum number of poses to return.

        Returns:
            list[tuple]: (ligand_id, pose_id, pose_coordinates, rdmol) per pending pose.
        """
        query = self.QueryBuilder()
        query.SELECT("L.ligand_id", "R.pose_id", "R.pose_coordinates", "L.rdmol")
        query.FROM("Results", "R").JOIN("Ligands", "L", "ligand_id")
        query.JOIN(track_table_name, "tt", "pose_id", kind="LEFT")
        query.WHERE("tt.pose_id IS NULL")
        query.ORDER_BY("L.ligand_id, R.pose_id")
        query.LIMIT(limit)
        sql, params = query.build()
        return self.db_query(sql, params).fetchall()

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
            self._apply_table_filter(query, selection)
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
        query = self.QueryBuilder()
        query.SELECT(
            "DISTINCT II.interaction_type",
            "II.rec_chain",
            "II.rec_resname",
            "II.rec_resid",
            "II.rec_atom",
        ).FROM("Interaction_indices", "II")
        query.JOIN("Interactions", "I", "interaction_id", "Interaction_indices")
        query.WHERE(
            f"I.pose_id IN ({self._get_scope_poses_query(bookmark_name)})"
        )
        # iterable of tuples
        interactions = self.db_query(*query.build()).fetchall()
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

    def accept_pose(self, pose_ids: int | list[int]):
        """
        Will add pose_id to accepted, and delete from maybe and rejected if needed

        Args:
            pose_ids (int | list[int])
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

    def maybe_pose(self, pose_ids: int | list[int]):
        """
        Will add pose_id to maybe, and delete from accepted and rejected if needed

        Args:
            pose_ids (int | list[int])
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

    def reject_pose(self, pose_ids: int | list[int]):
        """
        Will add pose_id to rejected, and delete from accepted and maybe if needed

        Args:
            pose_ids int | list[int])
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

    def remove_status(self, pose_ids: int | list[int]):
        """
        Will remove status associated with a pose

        Args:
            pose_ids (int | list[int])
        """
        if isinstance(pose_ids, int):
            pose_ids = [pose_ids]
        pose_ids = [(pose,) for pose in pose_ids]
        self.db_update(
            "DELETE FROM Accepted WHERE pose_id = ?;", pose_ids, commit=False
        )
        self.db_update("DELETE FROM Maybe WHERE pose_id = ?;", pose_ids, commit=False)
        self.db_update("DELETE FROM Rejected WHERE pose_id = ?;", pose_ids, commit=True)
    
    def get_column_metadata(self) -> dict[str, dict[str, str]]:
        '''
        Describe queryable columns, including any custom pose properties

        Returns:
            dict[str, dict[str, str]]: colum:{table,sql type, description if any}
        '''
        columns = {
            name.lower(): {
                "table": table,
                "sql_type": column.sql_type,
                "description": column.description or "",
            }
            for table, table_columns in OUTFIELD_BY_TABLE.items()
            for name, column in table_columns.items()
        }

        if self._is_table("custom_pose_data"):
            for name, sql_type in self.fetch_dynamic_table_column_names(
                "custom_pose_data"
            ).items():
                if name.lower() != "pose_id":
                    columns[name.lower()] = {
                        "table": "custom_pose_data",
                        "sql_type": sql_type,
                        "description": "",
                    }

        # Resolve shared identifiers to their primary tables.
        for name, table, schema in (
            ("pose_id", "Results", RESULTS_SCHEMA),
            ("ligand_id", "Ligands", LIGANDS_SCHEMA),
        ):
            column = schema.columns[name]
            columns[name] = {
                "table": table,
                "sql_type": column.sql_type,
                "description": column.description or "",
            }

        return columns

    def get_queryable_columns(self) -> dict[str, str]:
        '''
        Return column-to-table mappings, including custom pose properties

        Returns:
            dict[str, str]: {col:name}
        '''
        return {
            name: metadata["table"]
            for name, metadata in self.get_column_metadata().items()
        }

    def prepare_column_export_query(
        self,
        columns: dict[str, str] | list[str],
        bookmark: str| None = None,
    ) -> str:
        """
        Writes a query based on what columns are requested. Resolves each column
        to its source table automatically.

        Args:
            columns (dict | list): columns to export. A mapping of
                ``{column: header}`` renames the column in the output; a header of
                None, or a plain list of column names, keeps the column's own name.
            bookmark (str | None): bookmark, status table, or Candidates selection
                used to filter results. ``None`` or ``Results`` exports all results.

        Returns:
            str: SQL query string
        """
        headers = (
            dict(columns)
            if isinstance(columns, dict)
            else {col: None for col in columns}
        )

        column_to_table = self.get_queryable_columns()
        sources = {}
        for col in headers:
            table = column_to_table.get(col.lower())
            if table is None:
                raise OptionError(f"Column '{col}' not found in any queryable table.")
            sources[col] = table

        query, fields = self._prepare_column_query(sources, bookmark)

        for col, header in headers.items():
            field = fields[col]
            if header is None:
                query.SELECT(field)
                continue

            if not re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", str(header)):
                raise OptionError(
                    f"'{header}' is not usable as a column name. Use letters, "
                    "digits and underscores, and do not start with a digit."
                )
            query.SELECT(f'{field} AS "{header}"')

        return query.build()[0]

    def fetch_dynamic_table_column_names(self, table: str) -> dict[str, str]:
        '''
        Get columns from tables that may have user defined and dynamic columns

        Args:
            table (str): table

        Returns:
            dict[str, str]: col_name: SQL type
        '''
        return dict(
            self.db_query(
                "SELECT name, type FROM pragma_table_info(?)",
                [table],
            ).fetchall()
        )

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

    def ensure_gui_tables(self):
        """
        Ensures gui-specific tables exist in database
        """

        for sql in build_create_table(
            POSE_COMMENTS_SCHEMA.name, POSE_COMMENTS_SCHEMA, self.dialect
        ):
            self.db_query(sql)
        self.conn.commit()

    def check_ringtaildb_version(self) -> tuple[bool, str]:
        """Check the database schema version and whether this code can use it.

        The version is read from the ringtail_schema_version table. If that table
        is missing the lookup raises DatabaseQueryError: the database predates
        versioning (a pre-release lab build) and must be upgraded — reported as
        incompatible so callers prompt for the upgrade script. Otherwise
        compatibility is by major.minor (e.g. a future 3.1 schema is incompatible).

        Returns:
            tuple[bool, str]: (compatible, database schema version)
        """
        try:
            db_version = self._read_schema_version()
        except DatabaseQueryError:
            logger.warning(
                f"Database has no schema version, it predates versioning and must be "
                f"upgraded. {UPGRADE_HINT}"
            )
            return False, "unversioned"

        if not db_version:
            logger.warning(f"Database schema version is unset. {UPGRADE_HINT}")
            return False, "unversioned"

        db_v = Version(db_version)
        code_v = Version(SCHEMA_VERSION)
        compatible = (db_v.major, db_v.minor) == (code_v.major, code_v.minor)
        if not compatible:
            logger.warning(
                f"Database schema version {db_version} is NOT compatible with this "
                f"Ringtail (schema {SCHEMA_VERSION})."
            )
        return compatible, db_version

    def insert_receptor_blob(self, receptor: bytes, rec_name: str):
        """Takes object of Receptor class, updates the column in Receptor table

        Args:
            receptor (bytes): bytes receptor object to be inserted into DB
            rec_name (string): Name of receptor. Used to insert into correct row of DB
        """
        self._upsert_receptor_column("receptor_object", receptor, rec_name)

    def insert_receptor_polymer(self, receptor: str, rec_name: str):
        """Takes object of Receptor class, updates the column in Receptor table

        Args:
            receptor (str): json string representation of a receptor meeko.Polymer oobject to be inserted into DB
            rec_name (str): Name of receptor. Used to insert into correct row of DB
        """
        self._upsert_receptor_column("polymer", receptor, rec_name)

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
    def clone(self, backup_name: str = None) -> str:
        """Creates a copy of the db

        Args:
            backup_name (str, optional): name of the cloned database

        Returns:
            str: path of backed up database
        """
        ...

    @abstractmethod
    def update_database_version(self, new_version: str, backup=False):
        """Only SQLite databases from Ringtail 1.x/2.x can be upgraded, others are only checked to be current.

        Args:
            new_version (str): target version string, e.g. "3.0.0"
            backup (bool, optional): unused

        Raises:
            VersionError: if the schema is not current, e.g. a pre-release 3.0 database
        """
        ...

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
    def tables_in_db(self, temp: bool = False) -> list:
        """
        Returns a list of all table names in the database

        Args:
            temp (bool): whether it's a temporary table. Defaults to False

        Returns:
            list: list of table names
        """
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

    @abstractmethod
    def convert_pose_coordinates_to_native(self):
        """Convert pose_coordinates stored as JSON text to the backend's native format, a no-op if already converted."""
        ...

    @abstractmethod
    def check_attached(self) -> list[tuple]:
        """
        Checks what, if any, databases are attached

        Returns:
            list: attached databases
        """
        ...

    # endregion

    # region private api

    def _cluster_exists(
        self, cluster_name: str, cluster_window: str
    ) -> tuple[int] | None:
        """
        Gets the id of an existing clustering session with the given name and window.

        Args:
            cluster_name (str): standardized cluster name
            cluster_window (str): selection that was clustered

        Returns:
            tuple[int] | None: row containing cluster_id, or None if not found
        """
        query = self.QueryBuilder()
        query.SELECT("cluster_id").FROM("Clusters")
        query.WHERE("name = ?", cluster_name)
        query.WHERE("cluster_window = ?", cluster_window)
        return self.db_query(*query.build()).fetchone()

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
        query = self.QueryBuilder()
        query.INSERT_INTO(
            "Clusters",
            "name",
            "description",
            "cluster_window",
            "num_clusters",
        ).RETURNING("cluster_id")
        return self.db_query(
            query.build()[0],
            (name, description, cluster_window, length),
        ).fetchone()[0]

    def _insert_clusters(self, cluster_groups: list, pose_rows: list):
        """
        Inserts all cluster data, including each grouping and its representative
        pose, and all poses involved and which clusters they belong to

        Args:
            cluster_groups (list): each cluster from the clustering exercise
            pose_rows (list): pose and cluster id and group id
        """
        query = self.QueryBuilder()
        query.INSERT_INTO(
            "Cluster_groups", "cluster_id", "cluster_group", "representative"
        )
        self.db_update(query.build()[0], cluster_groups, commit=False)

        query = self.QueryBuilder()
        query.INSERT_INTO(
            "Pose_clusters", "cluster_id", "cluster_group", "pose_id"
        )
        self.db_update(query.build()[0], pose_rows, commit=False)

    def _delete_cluster(self, cluster_id: int):
        """
        Deletes all data associated with a clustering session.

        Args:
            cluster_id (int): clustering session to delete
        """
        for table in ("Pose_clusters", "Cluster_groups", "Clusters"):
            query = self.QueryBuilder()
            query.DELETE_FROM(table).WHERE("cluster_id = ?", cluster_id)
            self.db_query(*query.build())

    def _insert_completed_poses(
        self, pose_ids: list[tuple], tracking_table: str
    ):
        """
        Inserts processed poses and interaction cutoffs into the tracking table.

        Args:
            pose_ids (list[tuple]): rows of (pose_id, hb_cutoff, vdw_cutoff)
            tracking_table (str): table tracking poses with calculated interactions
        """
        query = self.QueryBuilder()
        query.INSERT_INTO(
            self._safe_quoting(tracking_table),
            "pose_id",
            "hb_cutoff",
            "vdw_cutoff",
        )
        self.db_update(query.build()[0], pose_ids, commit=False)

    def _insert_receptors(self, receptor_array: list):
        """Takes array of receptor rows, inserts into Receptors table

        Args:
            receptor_array (list): List of lists
                containing formatted receptor rows
        """
        query = self.QueryBuilder()
        query.INSERT_INTO(
            "Receptors",
            "recname",
            "box_dim",
            "box_center",
            "grid_spacing",
            "flexible_residues",
            "flexres_atomnames",
        )
        self.db_query(query.build()[0], receptor_array, commit=True)

    def _insert_db_properties(self, docking_mode: str, number_of_poses: str):
        """Insert db properties into database properties table

        Args:
            docking_mode (str): docking mode for the current dataset being written
            number_of_poses (str): number of poses written to database in current session, either "all" or specified max_poses
        """
        query = self.QueryBuilder()
        query.INSERT_INTO(
            "db_properties", "docking_mode", "number_of_poses"
        )
        self.db_query(
            query.build()[0],
            [docking_mode, number_of_poses],
            commit=True,
        )

    def _merging_receptors_compatible(self) -> bool:
        """
        Checks if the receptor names in the two databases match.

        Returns:
            bool: True if matching receptor rows have the same non-null name
        """
        query = self.QueryBuilder()
        query.SELECT("R.recname = M.recname").FROM("Receptors", "R")
        query.JOIN("merging.Receptors", "M", "receptor_id", "Receptors")
        row = self.db_query(*query.build()).fetchone()
        return bool(row and row[0])

    def _get_merge_id(self, mergingdb_path: str) -> int:
        """
        Creates a merge record for the given database and gets its id.

        Args:
            mergingdb_path (str): path to merging database

        Returns:
            int: id of the new merge record
        """
        query = self.QueryBuilder()
        query.INSERT_INTO("merged_tables", "dbfile").RETURNING("merge_id")
        return self.db_query(query.build()[0], [mergingdb_path]).fetchone()[0]

    def _delete_old_duplicate_results(self):
        """Checks if a pose is uniquely represented in the Results table,
        and deletes it from Results if duplicated.
        Based on the following columns:
        ligname,
        receptor,
        pose_coordinates,
        flexible_res_coordinates
        """
        query = self.QueryBuilder()
        query.SELECT("DISTINCT R.pose_id").FROM("Results", "R")
        query.JOIN(
            "Results_temp", "RT",
            condition="""
                (RT.receptor = R.receptor
                 OR (RT.receptor IS NULL AND R.receptor IS NULL))
                AND RT.pose_coordinates = R.pose_coordinates
                AND (RT.flexible_res_coordinates = R.flexible_res_coordinates
                     OR (RT.flexible_res_coordinates IS NULL
                         AND R.flexible_res_coordinates IS NULL))
            """,
        )
        query.JOIN(
            "Ligands", "L",
            condition="RT.ligname = L.ligname AND R.ligand_id = L.ligand_id",
        )
        delete_pose_ids_list = {
            row[0] for row in self.db_query(*query.build()).fetchall()
        }
        placeholders = ",".join("?" for _ in delete_pose_ids_list)
        if delete_pose_ids_list:
            self._delete_pose_screening_data(delete_pose_ids_list)

            tables = ["Interactions", "Results"]
            if self._is_table("custom_pose_data"):
                tables.append("custom_pose_data")

            for table in tables:
                query = self.QueryBuilder()
                query.DELETE_FROM(table).WHERE(
                    f"pose_id IN ({placeholders})", *delete_pose_ids_list
                )
                self.db_query(*query.build())

    def _delete_new_duplicate_results(self):
        """Checks if a pose is uniquely represented in the Results table,
        and deletes it from the staged incoming data if duplicated.
        Based on the following columns:
        ligname,
        receptor,
        pose_coordinates,
        flexible_res_coordinates
        """
        tables = ["Interactions_temp"]
        if self._is_table("Custom_temp", temp=True):
            tables.append("Custom_temp")

        for table in tables:
            matching = self.QueryBuilder()
            matching.SELECT("1").FROM("Results_temp", "RT")
            matching.JOIN(
                "Results", "R",
                condition="""
                    (RT.receptor = R.receptor
                     OR (RT.receptor IS NULL AND R.receptor IS NULL))
                    AND RT.pose_coordinates = R.pose_coordinates
                    AND (RT.flexible_res_coordinates = R.flexible_res_coordinates
                         OR (RT.flexible_res_coordinates IS NULL
                             AND R.flexible_res_coordinates IS NULL))
                """,
            )
            matching.JOIN(
                "Ligands", "L",
                condition="RT.ligname = L.ligname AND R.ligand_id = L.ligand_id",
            )
            matching.WHERE(f"RT.ligname = {table}.ligname")
            matching.WHERE(f"RT.pose_rank = {table}.pose_rank")
            matching.WHERE(f"RT.run_number = {table}.run_number")
            matching_sql, params = matching.build()

            query = self.QueryBuilder()
            query.DELETE_FROM(table).WHERE(f"EXISTS ({matching_sql})", *params)
            self.db_query(*query.build())

        matching = self.QueryBuilder()
        matching.SELECT("1").FROM("Results", "R")
        matching.JOIN(
            "Ligands", "L",
            condition="RT.ligname = L.ligname",
        )
        matching.WHERE("""
            (RT.receptor = R.receptor
             OR (RT.receptor IS NULL AND R.receptor IS NULL))
            AND RT.pose_coordinates = R.pose_coordinates
            AND (RT.flexible_res_coordinates = R.flexible_res_coordinates
                 OR (RT.flexible_res_coordinates IS NULL
                     AND R.flexible_res_coordinates IS NULL))
            AND L.ligand_id = R.ligand_id
        """)
        matching_sql, params = matching.build()

        query = self.QueryBuilder()
        query.DELETE_FROM("Results_temp AS RT").WHERE(
            f"EXISTS ({matching_sql})", *params
        )
        self.db_query(*query.build())

    def _crossref_bookmark_builder(self, ligand_list: list[str]) -> str:
        """
        creates formatted sql selecting all poses of the given ligands, for the crossreferencing bookmark

        Args:
            ligand_list (list[str]): list of ligand names to be considered

        Returns:
            str: formatted sql
        """
        query = self.QueryBuilder()
        query.SELECT("R.pose_id").FROM("Results", "R")
        if self.dialect == "sqlite":
            query.SELECT("R.ligand_id")

        query.JOIN("Ligands", "L", "ligand_id", "Results")

        if ligand_list:
            ligand_names = ", ".join(
                "'" + name.replace("'", "''") + "'"
                for name in ligand_list
            )
            query.WHERE(f"L.ligname IN ({ligand_names})")
        else:
            query.WHERE("1 = 0")

        return query.build()[0]

    def _merge_db_properties_table(self, merge_id: int):
        """
        Merges database properties table, but importantly will not check for property compatibility

        Args:
            merge_id (int): merge session id

        Raises:
            MergeError
        """
        try:
            convert_dbprop_sql = """INSERT INTO PK_conversions (
                merge_id,
                table_name,
                original_pk,
                merged_pk) SELECT 
                ?,
                'db_properties', 
                DB_write_session,
                DB_write_session + (SELECT MAX(DB_write_session) FROM db_properties) 
                FROM merging.db_properties;"""
            self.db_query(
                convert_dbprop_sql,
                [merge_id],
            )
            insert_dbprops_sql = """INSERT INTO db_properties (
                DB_write_session,
                docking_mode,
                number_of_poses)
                SELECT 
                    (SELECT merged_pk FROM PK_conversions WHERE original_pk = DB_write_session and merge_id = ? and table_name = 'db_properties'),
                    docking_mode,
                    number_of_poses
                FROM merging.db_properties;"""
            self.db_query(insert_dbprops_sql, [merge_id])

        except Exception as e:
            raise MergeError(
                f"Error encountered while merging db_properties table: {str(e)}"
            ) from e

    def _merge_ligands_and_results_tables(self, merge_id: int):
        """
        Merges first the Ligands table, then the Results table, maintaining ligand_id and pose_id as primary keys,
        where their relationship to the original PK is kept track of in the mering datble. Duplicate ligands will maintain
        the ligand_id from the main database.

        Args:
            merge_id (int): merge session id

        Raises:
            StorageError: if merging fails in SQLite
            MergeError: if merging fails in DuckDB
        """
        # convert ligand_ids and log
        convert_ligand_ids_sql = """INSERT INTO PK_conversions (
        merge_id,
        table_name,
        original_pk,
        merged_pk) SELECT
        ?,
        'Ligands',
        ligand_id,
            CASE
                WHEN EXISTS (
                    SELECT 1
                    FROM Ligands
                    WHERE
                        merging.Ligands.ligname = Ligands.ligname
                    ) 
                THEN (
                    SELECT main.Ligands.ligand_id
                    FROM main.Ligands
                    WHERE
                        merging.Ligands.ligname = Ligands.ligname
                    )
                ELSE merging.Ligands.ligand_id + (SELECT MAX(ligand_id) FROM Ligands)
            END AS new_ligand_id
        FROM merging.Ligands;"""

        # then inserting only those that aren't already in the table
        insert_new_ligands = """INSERT INTO Ligands (
        ligand_id,
        ligname,
        ligand_smile,
        rdmol)
        SELECT 
            (SELECT merged_pk FROM PK_conversions WHERE original_pk = ligand_id and merge_id = ? AND table_name = 'Ligands') new_id,
            ligname,
            ligand_smile,
            rdmol
        FROM merging.Ligands WHERE new_id > (SELECT MAX(ligand_id) FROM Ligands);
        """

        # convert pose_id
        convert_poseid_sql = """INSERT INTO PK_conversions (
        merge_id,
        table_name,
        original_pk,
        merged_pk) SELECT 
        ?,
        'Results', 
        pose_id,
        pose_id + (SELECT MAX(pose_id) FROM Results) 
        FROM merging.Results;"""

        # insert results with updated pose_ids
        insert_Results = """INSERT INTO Results (
            pose_id,
            ligand_id,
            receptor,
            pose_rank,
            run_number,
            docking_score,
            leff,
            delta,
            cluster_rmsd,
            cluster_size,
            reference_rmsd,
            energies_inter,
            energies_vdw,
            energies_electro,
            energies_flexlig,
            energies_flexlr,
            energies_intra,
            energies_torsional,
            unbound_energy,
            num_interactions,
            num_hb,
            pose_coordinates,
            flexible_res_coordinates) 
        SELECT 
            pose.merged_pk as pose_id,
            ligand.merged_pk as ligand_id,
            mr.receptor,
            mr.pose_rank,
            mr.run_number,
            mr.docking_score,
            mr.leff,
            mr.delta,
            mr.cluster_rmsd,
            mr.cluster_size,
            mr.reference_rmsd,
            mr.energies_inter,
            mr.energies_vdw,
            mr.energies_electro,
            mr.energies_flexlig,
            mr.energies_flexlr,
            mr.energies_intra,
            mr.energies_torsional,
            mr.unbound_energy,
            mr.num_interactions,
            mr.num_hb,
            mr.pose_coordinates,
            mr.flexible_res_coordinates
        FROM merging.Results mr
        LEFT JOIN (
            SELECT original_pk, merged_pk 
            FROM PK_conversions 
            WHERE table_name = 'Results' 
            AND merge_id = ?
            ) pose 
        ON pose.original_pk = mr.pose_id
        LEFT JOIN (
            SELECT original_pk, merged_pk 
            FROM PK_conversions 
            WHERE table_name = 'Ligands' 
            AND merge_id = ?
            ) ligand
        ON (mr.ligand_id = ligand.original_pk);"""

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

        except Exception as e:
            error_type = StorageError if self.dialect == "sqlite" else MergeError
            raise error_type(
                f"Error encountered while merging Ligands and Results tables: \n{str(e)}"
            ) from e

    def _merge_interaction_tables(self, merge_id: int):
        """
        Merges the interaction tables. Interaction definitions are unique and independent of the Results table, so we only
        insert those that are new with updated PK, and assign existing interaction_ids to those already described in primary db

        Args:
            merge_id (int): merge session id

        Raises:
            MergeError
        """
        convert_ii_sql = """INSERT INTO PK_conversions (
        merge_id,
        table_name,
        original_pk,
        merged_pk) SELECT 
        ?,
        'Interaction_indices', 
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
                ELSE merging.Interaction_indices.interaction_id + (SELECT COALESCE(MAX(interaction_id), 0) FROM Interaction_indices)
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
            (SELECT merged_pk FROM PK_conversions WHERE original_pk = interaction_id and merge_id = ? AND table_name = 'Interaction_indices') new_id,
            interaction_type,
            rec_chain,
            rec_resname,
            rec_resid,
            rec_atom,
            rec_atomid
        FROM merging.Interaction_indices WHERE new_id > (SELECT COALESCE(MAX(interaction_id), 0) FROM Interaction_indices);
        """

        # Adding new data to Interactions table with (updated) pose_id and interaction_id
        insert_interactions = """
        INSERT INTO Interactions (
        pose_id,
        interaction_id
        )    SELECT P.merged_pk as pose_id, II.merged_pk as interaction_id
                FROM merging.Interactions I
                LEFT JOIN (SELECT original_pk, merged_pk
                FROM PK_conversions
                WHERE table_name = 'Results' 
                AND merge_id = ?) P ON (I.pose_id = P.original_pk)
            LEFT JOIN (SELECT original_pk, merged_pk
                FROM PK_conversions
                WHERE table_name = 'Interaction_indices' 
                AND merge_id = ?) II ON (I.interaction_id = II.original_pk);"""

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

        except Exception as e:
            raise MergeError(
                f"Error during update and insertion of interactions: {e}"
            ) from e
        
    def _create_tables(self, commit:bool=True):
        """
        Creates all tables needed for a Ringtail database

        Args:
            commit (bool): defaults to True
        """
        self._create_ligands_table()
        self._create_results_table()
        self._create_receptors_table()
        self._create_interaction_index_table()
        self._create_interaction_table()
        self._create_db_properties_table()
        self._create_filtering_tables()
        self._create_status_tables()
        self._create_schema_version_table()

        self._stamp_schema_version(SCHEMA_VERSION)
        if commit:
            self.conn.commit()

    def _create_schema_version_table(self):
        """Create the ringtail_schema_version table (authoritative schema version)."""
        for sql in build_create_table(
            SCHEMA_VERSION_SCHEMA.name, SCHEMA_VERSION_SCHEMA, self.dialect
        ):
            self.db_query(sql)

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

    def _create_custom_pose_table(self):
        """Creates a table of each pose-interaction combination."""
        for sql in build_create_table(
            CUSTOM_POSE_DATA_SCHEMA.name, CUSTOM_POSE_DATA_SCHEMA, self.dialect
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

    def _create_screening_tables(self):
        """
        Create the filtering, status and GUI tables a screening session needs.
        """
        self._create_filtering_tables()
        self._create_status_tables()
        self.ensure_gui_tables()

    def _insert_docking_data(
        self,
        results: list[list],
        interactions: list[list],
        custom_data: list[dict] | None,
        duplicate_handling: str,
    ):
        """Takes list of database rows to insert, adds data to results table.
        First stages data in temporary tables, then handles duplicates (if requested),
        and finally transfers data from temporary tables into permanent storage and
        commits once at the end.

        Args:
            results (list): list of arrays containing formatted result rows
            interactions (list): list of interactions
            custom_data (list | None): custom data supplied as properties in the RDKit Mol
            options (dict): includes options on how to handle duplicates if there are any

        """
        custom_cols = {}
        if custom_data:
            custom_cols = {
                key: float if isinstance(value, (int, float)) else type(value)
                for entry in custom_data
                for key, value in entry.items()
            }
            # ensure the identifying columns use integers

            custom_cols.update(
                {
                    "ligname": str,
                    "run_number": int,
                    "pose_rank": int,
                }
            )

        self._create_temporary_results_tables(custom_columns=custom_cols)
        self._insert_results_in_temp_tables(results, interactions, custom_data)
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
            # the representative is picked from the cluster's own members
            for pose in cluster:
                pose_rows.append([cluster_id, group_index, pose])

        self._insert_clusters(cluster_groups, pose_rows)
        self.conn.commit()

        return cluster_bookmark

    def _check_for_custom_data(self) -> bool:
        """
        Checks if there is custom data in the current active insert session
        """
        if self._is_table(
            "Custom_temp", temp=True
        ):  # does temp tables show up in table list?
            # ensure custom table exists
            self._create_custom_pose_table()
            self._ensure_custom_columns()
            return True
        return False

    def _ensure_custom_columns(self):
        """
        Adds columns if necessary to permanent custom table from temporary custom table
        """
        identifiers = {"ligname", "run_number", "pose_rank"}

        permanent_columns = self.fetch_dynamic_table_column_names("custom_pose_data")
        existing = {name.lower() for name in permanent_columns}

        tmeporary_columns = self.fetch_dynamic_table_column_names("Custom_temp")

        missing_columns = {
            name: sql_type
            for name, sql_type in tmeporary_columns.items()
            if name.lower() not in existing | identifiers
        }
        if missing_columns:
            self._update_custom_columns("custom_pose_data", missing_columns)

    def _update_custom_columns(self, table: str, columns: dict[str, str]):
        for col, sql_type in columns.items():
            quoted_col = self._safe_quoting(col)
            self.db_query(f"""
                ALTER TABLE {table} 
                ADD COLUMN {quoted_col} {sql_type};
                """)

    def _is_table(self, table: str, temp: bool = False) -> bool:
        """
        Returns True if table name is actually a bookmark

        Args:
            table (str): name of table or bookmark to check

        Returns:
            bool: if table name is a bookmark
        """
        if table.lower() in self.tables_in_db(temp):
            return True
        else:
            return False

    def _drop_existing_tables(self):
        """Drops all tables, retrying those still referenced by a foreign key"""
        tables = self.tables_in_db()
        while tables:
            remaining = []
            for table in tables:
                try:
                    self._delete_table(table)
                except DatabaseQueryError:
                    remaining.append(table)
            if len(remaining) == len(tables):
                raise StorageError(f"Could not drop tables {remaining}.")
            tables = remaining
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

    def _get_numeric_columns(self) -> dict[str, str]:
        '''
        Return numeric Results and custom columns, excluding identifiers

        Returns:
            dict[str, str]: {col:table}
        '''
        return {
            name: metadata["table"]
            for name, metadata in self.get_column_metadata().items()
            if metadata["table"] in {"Results", "custom_pose_data"}
            and metadata["sql_type"] in NUMERIC_TYPES
            and name not in {"pose_id", "ligand_id"}
        }
    def _apply_table_filter(self, query: QueryBuilder, table: str) -> QueryBuilder:
        """Restricts a query on Results to the poses of a bookmark, status table or candidates.

        Args:
            query (QueryBuilder): query to modify in-place
            table (str): bookmark, status table, candidates, or None/"Results" for all results

        Raises:
            OptionError: if table is none of these

        Returns:
            QueryBuilder: the same query, modified
        """
        if table is None or table.lower() == "results":
            return query
        return query.WHERE(
            f"{query.aliased('Results')}.pose_id IN ({self._get_scope_poses_query(table)})"
        )

    def _prepare_column_query(
        self, sources: dict[str, str], table: str | None = None
    ) -> tuple[QueryBuilder, dict[str, str]]:
        """
        Makes a Results query, joining tables as needed for the requested columns
        and restricting poses to given selection

        Args:
            sources (dict[str, str]): {column: source table}
            table (str | None): bookmark, status table, Candidates, or
                None/"Results" for all results

        Returns:
            tuple[QueryBuilder, dict[str, str]]: query without SELECT, and
                {column: quoted source.column}
        """
        query = self.QueryBuilder().FROM("Results")
        included_tables = {source.lower() for source in sources.values()}

        if "ligands" in included_tables:
            query.JOIN("Ligands", "Ligands", "ligand_id", "Results")

        if "interaction_indices" in included_tables:
            query.JOIN("Interactions", "Interactions", "pose_id", "Results")
            query.JOIN(
                "Interaction_indices",
                "Interaction_indices",
                "interaction_id",
                "Interactions",
            )

        if "custom_pose_data" in included_tables:
            query.JOIN(
                "custom_pose_data", "custom_pose_data", "pose_id", "Results",
                kind="LEFT",
            )

        self._apply_table_filter(query, table)
        fields = {
            col: f"{source}.{self._safe_quoting(col)}"
            for col, source in sources.items()
        }
        return query, fields

    def _numeric_column_query(
        self, columns: list[str], table: str
    ) -> tuple[QueryBuilder, dict[str, str]]:
        """
        Will make a query in Results for numeric Results or custom columns, including custom data if any

        Args:
            columns (list[str]): column names
            table (str): bookmark, status table, candidates, or None/"Results" for all results

        Raises:
            OptionError: if table is not a valid selection

        Returns:
            tuple[QueryBuilder, dict[str, str]]: sql query without SELECT: quoted source.column field
        """
        numeric_columns = self._get_numeric_columns()
        sources = {
            col: numeric_columns[col.lower()]
            for col in columns
            if col.lower() in numeric_columns
        }

        return self._prepare_column_query(sources, table)

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

    #TODO future: allow custom column filtering
    def _generate_filtering_query(
        self,
        filters: Filters,
        bookmark_name: str,
        input_bookmark: str = None,
    ) -> str:
        """Takes filters, writes SQL filtering string.

        One path for both kinds of filtering, since a flat filter set is a one-leaf tree.
        Each leaf becomes an AND-group of WHERE conditions (numerical column conditions,
        include interactions via EXISTS, exclude interactions via NOT IN, ligand name via
        ligand_id IN, RDKit criteria as an in-memory scan rendered to pose_id IN), and
        leaves combine by their parent node's operator. RDKit criteria are supported
        anywhere in the tree, but each group holding one costs its own ligand scan.

        Args:
            filters (Filters): the filter specification to render
            bookmark_name (str): Name of bookmark to which passing poses are saved
            input_bookmark (str, optional): Filter over pre-filtered data instead of entire database

        Returns:
            str: SQL query returning pose_ids of passing results
        """
        filtering_window = "Results"

        if input_bookmark is not None and input_bookmark.lower() == "results":
            input_bookmark = None

        if input_bookmark is not None:
            if input_bookmark == bookmark_name:
                logger.error(
                    f"Specified 'input_bookmark' and 'bookmark_name' are the same: {bookmark_name}"
                )
                raise OptionError(
                    "'input_bookmark' and 'bookmark_name' cannot be the same! Please rename 'bookmark_name'"
                )
            if filters.uses_percentile():
                raise OptionError(
                    "Cannot use 'score_percentile' or 'le_percentile' with 'input_bookmark'."
                )
            filtering_window = f"(SELECT * FROM Results WHERE pose_id IN ({self._get_scope_poses_query(input_bookmark)}))"

        if filters.is_flat():
            # a lone leaf is the user's whole query, so a missing interaction is an error
            # rather than something that quietly makes one branch false
            predicate, lig_filters = self._render_filter_group(filters.leaf)
            where_conditions = [predicate] if predicate else []

            # RDKit substruct/position/max-atoms filters can't be expressed in SQL —
            # run the SQL part first, then post-filter the resulting poses in memory.
            if lig_filters:
                where_conditions = [
                    self._rdkit_pose_id_predicate(
                        predicate, lig_filters, filtering_window
                    )
                ]
        else:
            boolean_predicate = self._render_filter_expression(
                filters, filtering_window
            )
            if not boolean_predicate:
                raise OptionError("No filters were provided, cannot filter.")
            where_conditions = [boolean_predicate]

        where_sql = (
            f"WHERE {' AND '.join(where_conditions)}" if where_conditions else ""
        )
        select_cols = (
            "R.pose_id, R.ligand_id" if self.dialect == "sqlite" else "R.pose_id"
        )
        return f"SELECT {select_cols} FROM {filtering_window} R {where_sql}"

    def _render_filter_group(self, filters: Filters, strict: bool = True) -> tuple:
        """Render one AND-group of processed filters into a single SQL predicate.

        This is the per-group core shared by both the flat and nested filter paths:
        numerical + include/exclude interaction + ligand-name conditions joined by AND.
        RDKit ligand filters (substruct/position/max-atoms) can't be expressed in SQL,
        so they're returned to the caller to apply as an in-memory post-filter.

        Args:
            filters (Filters): the filter group to render. Grouping and column mapping are
                delegated to ``Filters.to_criteria``. Percentile cutoffs are resolved here,
                since that needs to read the data.
            strict (bool): True for flat filtering, where this group is the user's whole
                query: a group that sets no filters at all, or names an interaction absent
                from the database, is an error. False for a nested filter leaf, where an
                empty group is simply a no-op node and a missing interaction only makes that
                leaf false, so the rest of the boolean tree still evaluates.

        Raises:
            OptionError: if strict and the group sets no filters at all

        Returns:
            tuple[str | None, dict]: (AND-joined WHERE predicate — no surrounding
            parentheses or WHERE keyword, or None if the group renders nothing, remaining
            RDKit ligand filters)
        """
        criteria = filters.to_criteria(percentile_cutoff=self._calc_percentile_filter)
        if not criteria:
            if strict:
                raise OptionError("No filters were provided, cannot filter.")
            return None, {}

        where_conditions = []
        rdkit_lig_filters = {}

        if "numeric" in criteria:
            where_conditions.append(
                " AND ".join(
                    f"R.{column} {operator} {value}"
                    for column, operator, value in criteria["numeric"]
                )
            )

        if "interactions" in criteria:
            include_interactions, exclude_interactions = (
                self._prepare_interaction_indices_for_filtering(
                    interaction_list=criteria["interactions"], strict=strict
                )
            )
            if include_interactions:
                # EXISTS over the (interaction_id, pose_id) covering index. SQLite does a
                # fast per-pose point-lookup, DuckDB's optimizer de-correlates this into a
                # hash semi-join, so the same form is optimal for both dialects.
                where_conditions.append(
                    self._build_include_interaction_exists_conditions(
                        include_interactions, criteria["max_miss"]
                    )
                )
            if exclude_interactions:
                where_conditions.append(
                    f"R.pose_id NOT IN ({self._build_exclude_interaction_subquery(exclude_interactions)})"
                )

        if "ligand" in criteria:
            # copy so popping ligand_name/ligand_name_file doesn't mutate the caller's dict
            lig_filters = dict(criteria["ligand"])
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
                rdkit_lig_filters = lig_filters

        predicate = " AND ".join(where_conditions) if where_conditions else None
        return predicate, rdkit_lig_filters

    def _rdkit_pose_id_predicate(
        self, predicate: str, rdkit_filters: dict, filtering_window: str
    ) -> str:
        """Resolve RDKit (in-memory) ligand criteria to a ``pose_id IN (...)`` predicate.

        RDKit criteria (SMARTS/substructure, molweight, max heavy atoms) can't be expressed
        as SQL, so they're evaluated with one scan over ``filtering_window``, pre-narrowed by
        ``predicate`` so only poses that already passed the SQL part are streamed.

        The returned predicate *subsumes* ``predicate`` - the pose_id set was computed with it
        applied - so callers must use it in place of their SQL predicate, not in addition to it.
        Because the conjunction is baked into the id set, the result also composes correctly
        under OR in a nested filter expression.

        Args:
            predicate (str): SQL predicate used to narrow the scan, may be empty or None
            rdkit_filters (dict): RDKit ligand filters to apply in memory
            filtering_window (str): FROM target for the scan

        Returns:
            str: ``R.pose_id IN (...)``, with a NULL id list when nothing passes. NULL rather
            than an empty list because DuckDB rejects ``IN ()`` as a syntax error (SQLite
            accepts it), both dialects evaluate ``IN (NULL)`` as false.
        """
        rdkit_partial = (
            f" FROM {filtering_window} R "
            f"JOIN Ligands ON Ligands.ligand_id = R.ligand_id"
        )
        if predicate:
            rdkit_partial += f" WHERE {predicate}"
        passing = [
            pid
            for poseids in self._perform_rdkit_filtering(
                rdkit_partial, rdkit_filters
            ).values()
            for pid in poseids
        ]
        ids = ",".join(map(str, passing)) if passing else "NULL"
        return f"R.pose_id IN ({ids})"

    def _render_filter_expression(
        self, node, filtering_window: str = "Results", narrow_with: str = None
    ) -> str:
        """Recursively render a nested boolean filter expression into one SQL predicate.

        A node is either a Filters (a boolean node) or a Filter (a leaf). Boolean nodes
        wrap their children in parentheses, leaves render via ``_render_filter_group``.
        Arbitrary nesting depth is supported - the same recursion renders
        ``(a AND (b OR (c AND d)))``.

        RDKit leaves (SMARTS/substructure, molweight, max heavy atoms) can't be SQL, so
        each is evaluated up front to the set of pose_ids that pass it (one in-memory
        scan over ``filtering_window``) and rendered as ``pose_id IN (...)`` - which then
        composes into AND/OR like any other predicate. This is why an RDKit criterion is
        cheap as SQL to combine but costs one scan per leaf.

        Under an AND node the non-RDKit children are rendered first and their predicate is
        used to narrow the RDKit children's scans, so only poses that already passed the
        SQL siblings are streamed. Under OR no such narrowing is valid, since a sibling
        failing doesn't disqualify a pose.

        Args:
            node (Filters | Filter): a boolean node or a leaf
            filtering_window (str): FROM target for evaluating RDKit leaves
            narrow_with (str, optional): predicate from AND siblings, used only to
                pre-narrow RDKit scans

        Returns:
            str: parenthesized SQL predicate (empty string if the node contributes nothing)
        """
        if isinstance(node, Filters):
            op = node.op.upper()
            children = list(node.children)
            if op == "AND":
                # split so the cheap SQL siblings can narrow the expensive scans
                scans = [c for c in children if self._node_needs_rdkit(c)]
                plain = [c for c in children if c not in scans]
            else:
                scans, plain = [], children

            rendered = [
                r
                for r in (
                    self._render_filter_expression(child, filtering_window, narrow_with)
                    for child in plain
                )
                if r
            ]
            if scans:
                sibling_narrowing = " AND ".join(
                    p for p in ([narrow_with] if narrow_with else []) + rendered if p
                )
                rendered += [
                    r
                    for r in (
                        self._render_filter_expression(
                            child, filtering_window, sibling_narrowing or None
                        )
                        for child in scans
                    )
                    if r
                ]
            if not rendered:
                return ""
            if len(rendered) == 1:
                return rendered[0]
            return "(" + f" {op} ".join(rendered) + ")"

        # leaf: rendered by the same per-group core the flat path uses
        predicate, rdkit_lig_filters = self._render_filter_group(node, strict=False)
        if rdkit_lig_filters:
            # evaluate this leaf's RDKit criterion to a pose_id set (one scan), pre-narrowed
            # by this leaf's own SQL part and by any AND siblings
            pre = " AND ".join(p for p in (predicate, narrow_with) if p)
            return (
                "("
                + self._rdkit_pose_id_predicate(
                    pre, rdkit_lig_filters, filtering_window
                )
                + ")"
            )
        if not predicate:
            return ""
        return "(" + predicate + ")"

    @staticmethod
    def _node_needs_rdkit(node) -> bool:
        """True if this node, or anything under it, needs an in-memory ligand scan."""
        if isinstance(node, Filters):
            return any(leaf.has_rdkit_criteria() for leaf in node.leaves())
        return node.has_rdkit_criteria()

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
            xyz = [float(value) for value in pose_coordinates[hit_atom_indices[index]]]
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
    
    def _get_best_pose_per_ligand_condition(
        self, bookmark_name: str, results_alias: str = "R"
    ) -> str:
        """Select the best-scoring pose in a bookmark for each ligand.

        ``pose_id`` provides deterministic ordering when two poses have the same
        docking score. The correlated subquery is supported by both storage
        backends and avoids selecting arbitrary rows with ``GROUP BY``.
        """
        bookmark_poses = self._get_scope_poses_query(bookmark_name)
        return f"""
            {results_alias}.pose_id = (
                SELECT best.pose_id
                FROM Results AS best
                WHERE best.ligand_id = {results_alias}.ligand_id
                    AND best.pose_id IN ({bookmark_poses})
                ORDER BY best.docking_score ASC, best.pose_id ASC
                LIMIT 1
            )
        """

    def _get_scope_poses_query(self, scope_name: str, alias: str = "") -> str:
        """
        Returns the pose_id-selecting subquery for a selection name, in any case:
        a bookmark, a status table (Accepted/Maybe/Rejected), candidates, or Results.

        Args:
            scope_name (str | None): bookmark name, status table name, or None/"Results" for all results
            alias (str): attached-database alias to qualify the table with

        Raises:
            OptionError: if scope_name is none of these in that database

        Returns:
            str: query selecting the pose_ids in that scope
        """
        if scope_name is None or scope_name.lower() == "results":
            prefix = f"{alias}." if alias else ""
            return f"SELECT pose_id FROM {prefix}Results"
        if self._is_candidates_table(scope_name):
            prefix = f"{alias}." if alias else ""
            return f"SELECT pose_id FROM {prefix}Accepted UNION SELECT pose_id FROM {prefix}Maybe"
        if self._is_statustable(scope_name):
            prefix = f"{alias}." if alias else ""
            return f"SELECT pose_id FROM {prefix}{scope_name.capitalize()}"
        db_name = None if alias in ("", "main") else alias
        if scope_name.lower() not in self.get_all_bookmark_names(db_name):
            raise OptionError(
                f"{scope_name!r} is not a bookmark, status table, {CANDIDATES_NAME} or Results in this database."
            )
        return self.QueryBuilder.bookmark_query(scope_name.lower(), alias)

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

    def _prepare_interaction_indices_for_filtering(
        self, interaction_list: list, strict: bool = True
    ):
        """
        Prepare lists of interaction indices where they are grouped by whether or not they should be evaluated as "AND" or "OR",
        and whether to be excluded or included in the passing filter poses

        An *included* interaction that has no matching indices in the database is kept as an
        empty group rather than dropped: an unsatisfiable requirement must fail every pose
        (it renders as a false predicate), and keeping it preserves the group count that
        ``max_miss`` is measured against. An *excluded* interaction that isn't in the
        database is trivially satisfied, so it is dropped.

        Args:
            interaction_list (list): list of interactions
            strict (bool): if True, an interaction not present in the database is an error.
                If False (nested filter leaves), warn instead so that one unsatisfiable
                group does not abort the whole query.

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
                        "Given 'react_any' filter, but the database has no reactive interactions, so no pose passes it."
                    )
                    include_interactions.append([])
                    continue  # ends this iteration of the for loop
                # create string representation of ecah interaction not found
                interaction_not_found.append(":".join(interaction[:4]))
                if interaction[-1] is True:
                    # keep the unsatisfiable requirement as an empty group so it renders
                    # false instead of vanishing, and still counts toward max_miss
                    include_interactions.append([])
                # an excluded interaction absent from the database is trivially satisfied
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
        # if one or more interactions not found, raise error (flat filtering) or warn and
        # let the unsatisfiable groups render false (nested filter leaves)
        if interaction_not_found:
            if strict:
                raise OptionError(
                    f"The following interactions do not exist in the database: {interaction_not_found} not found in the database. Please check for spelling errors or remove from filter."
                )
            logger.warning(
                f"The following interactions do not exist in the database: {interaction_not_found}. "
                "Any filter group requiring them cannot match and is evaluated as false."
            )
        return include_interactions, exclude_interactions

    def _build_include_interaction_exists_conditions(
        self, include_groups: list, max_miss: int
    ) -> str:
        """Build correlated EXISTS conditions checking per-pose interaction requirements.

        Uses the (interaction_id, pose_id) covering index: for each candidate pose SQLite
        does a fast point-lookup instead of materializing the full interaction result set.
        For selective energy filters this is 10-100x faster than the IN (subquery) approach.

        Args:
            include_groups (list[list[int]]): each inner list is a group of interaction_ids
                that satisfy one required interaction. An empty group is a requirement with
                no matching interaction in the database and renders as false.
            max_miss (int): number of groups a pose is allowed to miss

        Returns:
            str: SQL WHERE fragment using R.pose_id (caller must alias the Results table as R)
        """

        def _group_exists(ids: list) -> str:
            if not ids:
                # no pose can satisfy this requirement, also avoids emitting "IN ()"
                return "0=1"
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
            input_bookmark (str, optional): If filters were performed across an existing bookmark. Defaults to None.

        Raises:
            StorageError
            OptionError

        Returns:
            int: 1 if the bookmark was written, 0 if no poses passed (and no bookmark is kept)
        """
        # make sure bookmark name is not a table name
        if name.lower() in self.tables_in_db():
            raise OptionError(
                f"Bookmark name {name} is the same as an existing table in the database, and cannot be used."
            )

        # the old bookmark is replaced, also when nothing passes the new filter
        old_filter_id = None
        if self.is_bookmark(name):
            logger.warning(
                f"The bookmark {name} already exists, and will be overwritten by the current filter."
            )
            lookup = self.QueryBuilder()
            lookup.SELECT("filter_id").FROM("Filters")
            lookup.WHERE("name = ?", name.lower())
            old_filter_id = self.conn.execute(*lookup.build()).fetchone()[0]

        last_ids = self.last_row_ids()
        self._begin_transaction()
        lookup = self.QueryBuilder()
        lookup.SELECT("COALESCE(MAX(filter_id), 0) + 1").FROM("Filters")
        new_filter_id = self.conn.execute(*lookup.build()).fetchone()[0]

        insert_query = self.QueryBuilder()
        insert_query.INSERT_INTO(
            "Filters",
            "filter_id",
            "name",
            "query",
            "filters",
            "filter_window",
        ).RETURNING("filter_id")

        # insert filter info, return filter_id
        filter_id = self.conn.execute(
            insert_query.build()[0],
            (
                new_filter_id,
                name.lower(),
                query,
                json.dumps(filters),
                # None when clustering all results; stored as filter_results does
                (input_bookmark or "Results").lower(),
            ),
        ).fetchone()[0]

        # single-pass INSERT, no COUNT pre-flight
        if self.dialect == "sqlite":
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
        try:
            self.conn.execute(insert_poses)
        except Exception:
            # DuckDB has already committed the new bookmark's Filters row
            self.delete_rows_after(last_ids)
            raise

        lookup = self.QueryBuilder()
        lookup.SELECT("1").FROM("Filtered_poses")
        lookup.WHERE("filter_id = ?", filter_id).LIMIT(1)
        row = self.conn.execute(*lookup.build()).fetchone()

        if row is None:
            # reusing a bookmark name still means replacement
            delete = self.QueryBuilder()
            delete.DELETE_FROM("Filters").WHERE("filter_id = ?", filter_id)
            self.conn.execute(*delete.build())

        # remove the old bookmark atomically
        if old_filter_id is not None:
            for table in ("Filtered_poses", "Filters"):
                delete = self.QueryBuilder()
                delete.DELETE_FROM(table).WHERE("filter_id = ?", old_filter_id)
                self.conn.execute(*delete.build())

        self.conn.commit()
        return 0 if row is None else 1
    
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

    def _delete_pose_screening_data(self, pose_ids: list[int]):
        """Delete status assignments and comments for the given poses, e.g. before
        those poses are replaced by duplicate_handling="replace".

        Args:
            pose_ids (Iterable[int]): poses whose status and comment to delete
        """
        pose_ids = tuple(pose_ids)
        if not pose_ids:
            return
        placeholders = ",".join("?" for _ in pose_ids)
        for table in ("Accepted", "Maybe", "Rejected", POSE_COMMENTS_SCHEMA.name):
            if self._is_table(table):
                self.db_query(
                    f"DELETE FROM {table} WHERE pose_id IN ({placeholders})", pose_ids
                )

    def _max_id(self, table: str, column: str) -> int:
        '''
        _summary_


        
        '''
        """Returns the largest id in a table column, 0 if the table is empty.

        Args:
            table (str): _description_
            column (str): which column to get max id for

        Returns:
                int
        """
        return self.db_query(
            f"SELECT COALESCE(MAX({column}), 0) FROM {table}"
        ).fetchone()[0]

    @staticmethod
    def _decode_flexres_column(value) -> list:
        """One stored flexres column as a list."""
        if not value:  # SQL NULL (row written by insert_receptor_blob/_polymer), ""
            return []
        if isinstance(value, (str, bytes, bytearray)):
            value = json.loads(value)  # "null" -> None, "[]" -> []
        return list(value) if value else []

    def _stamp_schema_version(self, schema_version: str = SCHEMA_VERSION):
        """Append a row to ringtail_schema_version recording the schema version and
        the ringtail package version that wrote it."""
        query = self.QueryBuilder()
        query.INSERT_INTO(
            SCHEMA_VERSION_SCHEMA.name, "schema_version", "ringtail_version"
        )
        self.db_query(
            query.build()[0],
            (schema_version, __version__),
            commit=True,
        )

    def _read_schema_version(self, db_alias: str = None) -> str:
        """Latest stamped schema version of the current database (db_alias=None) or
        an attached one. Returns the version string, or None if the table exists but
        is empty. Raises DatabaseQueryError if the table is absent (unversioned db).
        Built via QueryBuilder so the optional db-qualification and ORDER BY/LIMIT are
        rendered for the active dialect.
        """
        query = self.QueryBuilder()
        query.SELECT("schema_version").FROM(
            SCHEMA_VERSION_SCHEMA.name, db_name=db_alias
        ).ORDER_BY("id").DESC(True).LIMIT(1)
        row = self.db_query(query.build()[0]).fetchone()
        return row[0] if row and row[0] else None

    def _has_results_table(self) -> bool:
        """Whether the database holds Ringtail content (a Results table). Used to
        tell a brand-new/empty database from a pre-versioning one when the schema
        version table is absent."""
        query = self.QueryBuilder()
        query.SELECT("1").FROM(RESULTS_SCHEMA.name).LIMIT(1)
        try:
            self.db_query(query.build()[0])
        except DatabaseQueryError:
            return False
        return True

    def _upsert_receptor_column(self, column: str, value, rec_name: str) -> None:
        """Write one column of the single receptor row, creating the row if absent.

        A new row gets empty flexres arrays explicitly rather than leaving them to the
        column default, so databases created before that default also stay decodable.

        Args:
            column (str): receptor column to write, e.g. "polymer"
            value: value for that column
            rec_name (str): receptor name
        """
        if self.table_length("Receptors") == 0:
            query = (
                f"INSERT INTO Receptors (recname, {column}, "
                "flexible_residues, flexres_atomnames) VALUES (?,?,?,?);"
            )
            params = (rec_name, value, NO_FLEXRES_JSON, NO_FLEXRES_JSON)
        else:
            query = (
                f"UPDATE Receptors SET recname = ?, {column} = ? WHERE receptor_id = 1;"
            )
            params = (rec_name, value)
        self.db_query(query, params, commit=True)

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

    def _calc_percentile_filter(
        self, percentile: float, column: str = "docking_score"
    ) -> tuple[float, str]:
        """Return a percentile cutoff and its comparison operator.

        Positive percentiles that retain at least one ligand include the cutoff.
        Percentiles whose floor-rounded count is zero use an exclusive comparison
        against the best value and therefore return no ligands.
        """
        ligand_count = int(
            self.db_query(
                "SELECT COUNT(*) FROM (SELECT ligand_id FROM Results GROUP BY ligand_id)"
            ).fetchone()[0]
        )
        cutoff = self._calc_percentile_cutoff(percentile, column)
        accepted_count = int(percentile * ligand_count / 100)
        return cutoff, "<=" if accepted_count > 0 else "<"
   
    def _calc_percentile_cutoff(self, percentile: float, column="docking_score"):
        """
        Will return the per-ligand cutoff for a percentile of a numeric column, including custom columns.
        Keeps floor(percentile * ligand count) ligands, so the cutoff is the best value of the last ligand kept.
        If that count is zero it returns the best value, and the caller compares exclusively to keep no ligands.

        Args:
            percentile (float): cutoff percentile
            column (str, optional): numeric Results or custom column. Defaults to "docking_score".

        Raises:
            OptionError: if column is not numeric
            StorageError: if no ligand has a value for column

        Returns:
            float: cutoff value for the percentile
        """
        per_ligand, fields = self._numeric_column_query([column], None)
        if not fields:
            raise OptionError(
                f"Requested column {column} is not a numeric Results column, percentiles cannot be calculated."
            )
        field = fields[column]
        aggregate = "MIN" if column.lower() in MINIMIZED_COLUMNS else "MAX"
        direction = "ASC" if aggregate == "MIN" else "DESC"
        per_ligand.SELECT(f"{aggregate}({field}) AS best_value").WHERE(
            f"{field} IS NOT NULL"
        ).GROUP_BY("Results.ligand_id")
        per_ligand_sql, params = per_ligand.build()

        ligand_count = self.db_query(
            f"SELECT COUNT(*) FROM ({per_ligand_sql}) per_ligand", params
        ).fetchone()[0]
        if ligand_count == 0:
            raise StorageError(
                f"Cannot calculate a percentile, no ligand has a value for {column}."
            )

        accepted_count = int(percentile * ligand_count / 100)
        cutoff, value_type = self.db_query(
            f"""SELECT best_value, typeof(best_value) FROM ({per_ligand_sql}) per_ligand
            ORDER BY best_value {direction} LIMIT 1 OFFSET ?""",
            params + [max(accepted_count - 1, 0)],
        ).fetchone()
        if self.dialect == "duckdb" and value_type == "FLOAT":
            # DuckDB FLOAT is 4 bytes, -6.66 reads back as -6.659999847412109 and compares wrongly, its shortest form compares exactly
            cutoff = float(str(np.float32(cutoff)))
        logger.debug(f"{column} percentile cutoff is {cutoff}")
        return cutoff

    def _format_output_fields(self, outfields: str | list) -> list[str]:
        """
        Handles string or list input of column names to be outputted, will make sure ligname
        is in the list, and make sure all options are valid

        Args:
            outfields (str | list): column names, or None for ligname and docking_score

        Returns:
            list[str]: valid column names, including custom pose columns
        """
        if outfields is None:
            outfields_list = ["ligname", "docking_score"]
        elif isinstance(outfields, str):
            outfields_list = [col.strip() for col in outfields.split(",")]
        elif isinstance(outfields, (list, tuple)):
            outfields_list = list(outfields)
        else:
            logger.warning(
                "The provided outfields is not in a usable format (string or list). "
                "Using ligname and docking_score."
            )
            outfields_list = ["ligname", "docking_score"]
        if "ligname" not in [field.lower() for field in outfields_list]:
            outfields_list.insert(0, "ligname")
        possible_columns = self.get_queryable_columns()
        valid_outfields = []
        for outfield in outfields_list:
            if outfield.lower() in possible_columns:
                valid_outfields.append(outfield)
            else:
                logger.warning(
                    f"{outfield} is not a valid output option, and will be removed from the output columns. Please see rt_process_vs.py --help for allowed options"
                )
        return valid_outfields

    def _db_compatible_for_merge(self, merging_db_alias: str) -> bool:
        """Whether the database being merged in is compatible with this one.

        Both must carry the same major.minor schema version, read from each
        database's ringtail_schema_version table. A database with no such table
        predates versioning and must be upgraded before merging (see
        check_ringtaildb_version). Dialect-agnostic via QueryBuilder: the current
        database is read unqualified and the incoming one via its attach alias.

        Args:
            merging_db_alias (str): alias under which the merging database is attached

        Returns:
            bool: if the two databases are compatible for merging
        """
        try:
            main_version = self._read_schema_version()
            merging_version = self._read_schema_version(merging_db_alias)
        except DatabaseQueryError:
            logger.error(
                f"One of the databases has no schema version (it predates versioning) "
                f"and must be upgraded before merging. {UPGRADE_HINT}"
            )
            return False

        if not (main_version and merging_version):
            return False

        main_v, merging_v = Version(main_version), Version(merging_version)
        if (main_v.major, main_v.minor) != (merging_v.major, merging_v.minor):
            logger.error(
                f"Cannot merge: schema versions differ (this database "
                f"{main_version}, merging database {merging_version})."
            )
            return False
        return True

    def _rollback_merge(self, merge_id: int, last_ids: dict):
        """
        Remove all data inserted by a specific merge session

        Args:
            merge_id (int): Merge session for which to delete associated data
            last_ids (dict): output of last_row_ids from before the merge
        """
        self.delete_rows_after(last_ids)
        self.db_query("DELETE FROM PK_conversions WHERE merge_id = ?", (merge_id,))
        self.db_query("DELETE FROM merged_tables WHERE merge_id = ?", (merge_id,))
        self.conn.commit()
        self._sync_auto_increment_state()

    def _safe_quoting(self, text: str) -> str:
        """
        Applies safe quotation marks for database strings that may have space or quotation marks
        in places that make them sql unsafde

        Args:
            text (str): string

        Returns:
            str: db safe string
        """
        return '"' + text.replace('"', '""') + '"'

    @staticmethod
    def _serialize_pose_coordinates(coords):
        """Serialize a list of [x, y, z] floats for storage in the pose_coordinates
        column. Base/DuckDB store the native nested list (DuckDB FLOAT[][]), SQLite
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
                # Only an initialized database (one with a Results table) needs a
                # version check. An empty/uninitialized file — e.g. DuckDB writes a
                # header on connect, so it is >0 bytes but has no tables — is treated
                # like a brand-new database and created into.
                if self._has_results_table():
                    compatible, db_version = self.check_ringtaildb_version()
                    if not compatible:
                        raise VersionError(
                            f"Database schema version '{db_version}' is not compatible "
                            f"with this Ringtail (schema {SCHEMA_VERSION}). {UPGRADE_HINT}"
                        )
                else:
                    logger.info("Creating a new database file.")
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

    # endregion

    # region virtual private methods
    @abstractmethod
    def _create_connection(self):
        """Creates database connection to self.db_file

        Returns:
            <db type>conn: Connection object to self.db_file

        Raises:
            DatabaseConnectionError
        """
        ...

    @abstractmethod
    def _begin_transaction(self):
        """
        Begin a transaction
        """
        ...

    @abstractmethod
    def _rollback(self):
        """
        Roll back transaction
        """
        ...

    @abstractmethod
    def _create_temporary_results_tables(self):
        """
        Creates temporary tables for results and interactions, which will be
        used for staging incoming data.
        """
        ...

    @abstractmethod
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
        ...

    @abstractmethod
    def _move_tempresults_to_database(self):
        """Inserts data from the temporary results tables to their permanent
        database equivalents"""
        ...

    @abstractmethod
    def _insert_interactions(self, interactions: list[tuple]):
        """
        Inserts interaction tuples with pose_id s (IMPORTANT) into the interaction table.

        Args:
            interactions (list[tuple]): _description_
        """
        ...

    @abstractmethod
    def _update_interaction_counts(self, data: list[dict]):
        """
        Data is a dict that is expected to contain pose_id, num_int, and num_hb as keys.

        Args:
            data (list[dict]): _description_
        """
        ...

    @abstractmethod
    def _insert_ligands(self, ligands: list):
        """Takes array of ligand rows, inserts into Ligands table.

        Args:
            ligand_array (list[list]): list of lists containing formatted ligand rows

        """
        ...

    def _create_indices(self):
        """
        Creates specific indices on tables for those databases that use indices
        """
        pass

    @abstractmethod
    def _insert_interaction_index_rows(self, interactions: list):
        """
        Writes unique interactions to database

        Args:
            interactions (list): list of InteractionRecord or dicts with interaction description fields
        """
        ...

    def _delete_nontables(self):
        """
        Deletes objects in the database that are not tables
        """
        pass

    @abstractmethod
    def _attach_db(self, new_db: str, new_db_alias: str = "attached_db") -> str:
        """Attaches new database file to current database

        Args:
            new_db (str): file name for database to attach
            new_db_name (str): name of new database

        """
        ...

    @abstractmethod
    def _sync_auto_increment_state(self):
        """Reset any auto incremented sequences or counts according to last successful merge"""
        ...

    @abstractmethod
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
        ...

    # endregion
