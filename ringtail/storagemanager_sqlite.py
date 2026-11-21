#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail storage adaptor for SQLite3
#

from .logutils import get_logger

logger = get_logger(__name__)
import os
import shutil
import json
import numpy as np
from packaging.version import Version
from rdkit import Chem
from .exceptions import (
    StorageError,
    DatabaseInsertionError,
    DatabaseConnectionError,
    DatabaseQueryError,
    MergeError,
)
from .clustermanager import *
from .storagemanager import StorageManager
from .querybuilder import QueryBuilderSQLite
from .schema import (
    build_create_table,
    build_create_indices,
    CANDIDATES_SUBQ,
    RESULTS_SCHEMA,
    INTERACTION_INDICES_SCHEMA,
    INTERACTIONS_SCHEMA,
    MERGED_TABLES_SCHEMA,
    PK_CONVERSIONS_SCHEMA,
    LIGANDS_SCHEMA,
    FILTERS_SCHEMA,
    FILTERED_POSES_SCHEMA,
    CLUSTERS_SCHEMA,
    CLUSTER_GROUPS_SCHEMA,
    POSE_CLUSTERS_SCHEMA,
)

try:
    import sqlite3

    HAS_SQLITE = True
except ImportError:
    HAS_SQLITE = False


class StorageManagerSQLite(StorageManager):
    """SQLite-specific StorageManager subclass

    Attributes:
        db_file (str): database name
        conn (SQLite.conn): Connection to database
    """

    QueryBuilder = QueryBuilderSQLite

    def __init__(
        self,
        db_file: str = None,
    ):
        self.db_file = db_file
        super().__init__()
        self.conn: sqlite3.Connection
        self.dialect = "sqlite"

    # region Methods for creating and inserting into tables the database

    def _insert_ligands(self, ligands: list):
        """Takes list of ligand rows, inserts into Ligands table using executemany.

        Args:
            ligand_array (list[list]): list of lists containing formatted ligand rows

        Raises:
            DatabaseInsertionError

        """

        sql_insert = """INSERT INTO Ligands (
        ligname,
        ligand_smile,
        rdmol
        ) VALUES
        (?,?,?)
        ON CONFLICT(ligname) DO NOTHING"""

        self.db_update(sql_insert, ligands, commit=False)

    def _create_temporary_results_tables(
        self, custom_columns: dict[str, type] | None = None
    ):
        """
        Creates temporary tables for results and interactions, which will be
        used for staging incoming data, uncluding custom data if provided in the docking files.
        Args:
            custom_columns (dict[str, type] | None, optional): col: type if any custom columns. Defaults to None.
        """

        create_temp_results = """
        CREATE TEMP TABLE IF NOT EXISTS
        Results_temp(
            receptor            VARCHAR,
            pose_rank           INTEGER,
            run_number          INTEGER,
            docking_score       FLOAT,
            leff                FLOAT,
            delta               FLOAT,
            cluster_rmsd        FLOAT,
            cluster_size        INTEGER,
            reference_rmsd      FLOAT,
            energies_inter      FLOAT,
            energies_vdw        FLOAT,
            energies_electro    FLOAT,
            energies_flexlig    FLOAT,
            energies_flexlr     FLOAT,
            energies_intra      FLOAT,
            energies_torsional  FLOAT,
            unbound_energy      FLOAT,
            num_interactions     INTEGER,
            num_hb              INTEGER,
            pose_coordinates  BLOB,
            flexible_res_coordinates   VARCHAR,
            ligname             VARCHAR);
            """
        create_temp_interactions = """
        CREATE TEMP TABLE IF NOT EXISTS Interactions_temp(
            ligname             VARCHAR,
            run_number          INTEGER,
            pose_rank           INTEGER,
            interaction_type    VARCHAR,
            rec_chain           VARCHAR,
            rec_resname         VARCHAR,
            rec_resid           VARCHAR,
            rec_atom            VARCHAR,
            rec_atomid          VARCHAR);
        """
        create_temp_mapping_table = """
        CREATE TEMP TABLE IF NOT EXISTS pose_map(
            pose_id             INTEGER,
            ligand_id           INTEGER,
            run_number          INTEGER,
            pose_rank           INTEGER);
        """
        # create temporary tables
        self.db_query(create_temp_results)
        self.db_query(create_temp_interactions)
        self.db_query(create_temp_mapping_table)

        self.db_query("DROP TABLE IF EXISTS Custom_temp")
        if custom_columns:
            sql_types = {
                int: "INTEGER",
                float: "DOUBLE",
                bool: "DOUBLE",
                str: "VARCHAR",
            }
            columns = []
            for name, python_type in custom_columns.items():
                quoted_name = self._safe_quoting(name)
                columns.append(f"{quoted_name} {sql_types[python_type]}")

            column_definitions = ",\n".join(columns)
            create_temp_custom = f"""
            CREATE TEMP TABLE Custom_temp({column_definitions});
            """
            self.db_query(create_temp_custom)

    def _insert_results_in_temp_tables(
        self,
        results_array: list[dict],
        interactions_array: list,
        custom_data: list | None = None,
    ):
        """
        Inserts docking results and interactions into their respective
        temporary tables

        Args:
            results_array (list): list of result rows
            interactions_array (list): list of interaction rows
            custom_data (list | None, optional): custom data if any. Defaults to None.
        """
        temp_results_insert = """
            INSERT INTO Results_temp(
                receptor,
                pose_rank,
                run_number,
                cluster_rmsd,
                reference_rmsd,
                docking_score,
                leff,
                delta,
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
                cluster_size,
                pose_coordinates,
                flexible_res_coordinates,
                ligname)
            VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
            """
        # pack pose_coordinates (a list of [x,y,z] floats) into a float32 BLOB
        packed = [
            tuple(
                r._replace(
                    pose_coordinates=self._serialize_pose_coordinates(
                        r.pose_coordinates
                    )
                )
            )
            for r in results_array
        ]
        self.db_update(temp_results_insert, packed, commit=False)

        temp_int_insert = """
            INSERT INTO Interactions_temp (
                ligname,
                run_number,
                pose_rank,
                interaction_type,
                rec_chain,
                rec_resname,
                rec_resid,
                rec_atom,
                rec_atomid)
            VALUES (?,?,?,?,?,?,?,?,?)
            """

        self.db_update(
            temp_int_insert, [tuple(r) for r in interactions_array], commit=False
        )

        if custom_data:
            columns = list(dict.fromkeys(key for entry in custom_data for key in entry))
            quoted_columns = [self._safe_quoting(col) for col in columns]
            rows = [tuple(entry.get(col) for col in columns) for entry in custom_data]
            placeholders = ", ".join("?" for _ in columns)
            temp_custom_insert = f"""
                INSERT INTO Custom_temp ({",".join(quoted_columns)})
                VALUES ({placeholders});"""
            self.db_update(temp_custom_insert, rows, commit=False)

    def _move_tempresults_to_database(self):
        """Inserts data from the temporary results tables to their permanent
        database equivalents"""

        temp_to_results = """
            INSERT INTO Results (
                ligand_id,
                receptor,
                pose_rank,
                run_number,
                cluster_rmsd,
                reference_rmsd,
                docking_score,
                leff,
                delta,
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
                cluster_size,
                pose_coordinates,
                flexible_res_coordinates
                ) 
            SELECT 
                L.ligand_id, 
                T.receptor,
                T.pose_rank,
                T.run_number,
                T.cluster_rmsd,
                T.reference_rmsd,
                T.docking_score,
                T.leff,
                T.delta,
                T.energies_inter,
                T.energies_vdw,
                T.energies_electro,
                T.energies_flexlig,
                T.energies_flexlr,
                T.energies_intra,
                T.energies_torsional,
                T.unbound_energy,
                T.num_interactions,
                T.num_hb,
                T.cluster_size,
                T.pose_coordinates,
                T.flexible_res_coordinates
            FROM Results_temp AS T
            JOIN Ligands AS L ON L.ligname = T.ligname
            RETURNING pose_id, ligand_id, run_number, pose_rank;"""

        temp_to_interaction = """
            INSERT INTO Interactions(pose_id, interaction_id)
            SELECT M.pose_id, II.interaction_id
            FROM Interactions_temp IT
            JOIN pose_map M
                ON IT.pose_rank = M.pose_rank
                AND IT.run_number = M.run_number
            JOIN Ligands L
                ON M.ligand_id = L.ligand_id
                AND IT.ligname = L.ligname
            JOIN Interaction_indices II
                ON II.interaction_type = IT.interaction_type
                AND II.rec_chain = IT.rec_chain
                AND II.rec_resname = IT.rec_resname
                AND II.rec_resid = IT.rec_resid
                AND II.rec_atom = IT.rec_atom
                AND II.rec_atomid = IT.rec_atomid;"""
        mapping = self.db_query(temp_to_results)
        self.db_update(
            "INSERT INTO pose_map (pose_id, ligand_id, run_number, pose_rank) VALUES (?, ?, ?, ?)",
            mapping,
        )
        self.db_query(temp_to_interaction)
        # first check if we need to add new columns to

        has_custom = self._check_for_custom_data()
        if has_custom:
            identifiers = {"ligname", "run_number", "pose_rank"}
            custom_columns = [
                name
                for name in self.fetch_dynamic_table_column_names("Custom_temp")
                if name.lower() not in identifiers
            ]
            quoted_columns = [self._safe_quoting(name) for name in custom_columns]
            insert_columns = ", ".join(quoted_columns)
            select_columns = ", ".join(f"CT.{col}" for col in quoted_columns)
            temp_to_custom_pose_data = f"""
                INSERT INTO custom_pose_data(pose_id, {insert_columns})
                SELECT M.pose_id, {select_columns}
                FROM Custom_temp AS CT
                JOIN pose_map AS M
                    ON CT.pose_rank = M.pose_rank
                    AND CT.run_number = M.run_number
                JOIN Ligands AS L
                    ON M.ligand_id = L.ligand_id
                    AND CT.ligname = L.ligname;"""
            self.db_query(temp_to_custom_pose_data)
        # clear temp tables
        self.db_query("DELETE FROM pose_map")
        self.db_query("DELETE FROM Interactions_temp")
        self.db_query("DELETE FROM Results_temp")
        if has_custom:
            self.db_query("DELETE FROM Custom_temp")



    def _insert_interactions(self, interactions: list[tuple]):
        """
        Inserts interaction tuples with pose_id s (IMPORTANT) into the interaction table.

        Args:
            interactions (list[tuple]): _description_
        """
        # create staging table
        create_temp_interactions = """
        CREATE TEMP TABLE IF NOT EXISTS Interactions_temp(
            pose_id             INTEGER,
            interaction_type    VARCHAR,
            rec_chain           VARCHAR,
            rec_resname         VARCHAR,
            rec_resid           VARCHAR,
            rec_atom            VARCHAR,
            rec_atomid          VARCHAR);
        """
        self.db_query(create_temp_interactions)

        # insert into staging table
        insert_sql = """
            INSERT INTO Interactions_temp (
                pose_id, 
                interaction_type, 
                rec_chain, 
                rec_resname, 
                rec_resid, 
                rec_atom, 
                rec_atomid)
            VALUES (
                :pose_id,
                :type,
                :chain,
                :residue,
                :resid,
                :recname,
                :recid);
            """
        self.db_update(insert_sql, interactions, commit=False)

        # move to permanent table
        temp_to_interaction = """
            INSERT INTO Interactions(pose_id, interaction_id)
            SELECT IT.pose_id, II.interaction_id
            FROM Interactions_temp IT
            JOIN Interaction_indices II
                ON II.interaction_type = IT.interaction_type
                AND II.rec_chain = IT.rec_chain
                AND II.rec_resname = IT.rec_resname
                AND II.rec_resid = IT.rec_resid
                AND II.rec_atom = IT.rec_atom
                AND II.rec_atomid = IT.rec_atomid;"""
        self.db_query(temp_to_interaction)

        # delete staging table
        self._delete_table("Interactions_temp")

    def _update_interaction_counts(self, data: list[dict]):
        """
        Data is a dict that is expected to contain pose_id, num_int, and num_hb as keys.

        Args:
            data (list[dict]): _description_
        """
        insert_sql = """
        UPDATE RESULTS SET 
        num_interactions = :num_int,
        num_hb = :num_hb
        WHERE
        pose_id = :pose_id;
        """
        self.db_update(insert_sql, data, commit=False)




    def _insert_interaction_index_rows(self, interactions: list[dict]):
        """
        Writes unique interactions to database

        Args:
            interactions (list[dict]): [(interaction_type, rec_chain, rec_resname, rec_resid, rec_atom, rec_atomid)]
        """
        # to insert interaction if unique
        sql_insert = """INSERT OR IGNORE INTO Interaction_indices (interaction_type,rec_chain,rec_resname,rec_resid,rec_atom,rec_atomid)
                        VALUES (?,?,?,?,?,?);"""
        self.db_update(
            sql_insert,
            [self._interaction_index_fields(r) for r in interactions],
            commit=False,
        )

    def create_transaction_tracking_table(self, table_name: str = "tracking_table"):
        """
        Creates a table to track pose ids to which are part of a multiple transaction
        database operation, which may occur over multiple commits and database connections.

        Args:
            table_name (str, optional): Name of the table. Defaults to "tracking_table".
        """
        # check if interaction calcs have been started (and not completed)
        exists = bool(
            self.db_query(
                f"""SELECT name FROM sqlite_master WHERE type='table' AND name='{table_name}';"""
            ).fetchone()
        )
        if exists:
            logger.info(
                "Interaction calculations have been started in a previous database connection, and not yet finished."
            )

        table_sql = f"""CREATE TABLE IF NOT EXISTS {table_name} (
        pose_id INTEGER,
        hb_cutoff DOUBLE,
        vdw_cutoff DOUBLE,
        FOREIGN KEY (pose_id) REFERENCES Results(pose_id)
        );"""
        self.db_query(table_sql, commit=True)





    def _create_indices(self):
        """Create alternate-key indices ('ak_*') for queryable tables."""
        logger.debug("Creating columns indices...")
        for table, schema in [
            (RESULTS_SCHEMA.name, RESULTS_SCHEMA),
            (INTERACTIONS_SCHEMA.name, INTERACTIONS_SCHEMA),
            (INTERACTION_INDICES_SCHEMA.name, INTERACTION_INDICES_SCHEMA),
            (FILTERED_POSES_SCHEMA.name, FILTERED_POSES_SCHEMA),
        ]:
            for sql in build_create_indices(table, schema):
                self.db_query(sql)
        self.db_query("ANALYZE;")
        self.conn.commit()
        logger.info(
            "Indices created for Results, Interactions, and Interaction_indices."
        )

    def _remove_screening_tables(self):
        '''
        Deletes all data included in screening such as clusters, filtered data, status assignments, and comments
        '''
        self._delete_table(POSE_CLUSTERS_SCHEMA.name)
        self._delete_table(CLUSTER_GROUPS_SCHEMA.name)
        self._delete_table(CLUSTERS_SCHEMA.name)
        self._delete_table(FILTERED_POSES_SCHEMA.name)
        self._delete_table(FILTERS_SCHEMA.name)
        self._delete_table("Accepted")
        self._delete_table("Maybe")
        self._delete_table("Rejected")
        self._delete_table("Pose_comments")

    # endregion

    # region merge databases
    def prepare_for_merging(self):
        """
        Prepares for merging by dropping tables incompatible with merging
        """
        # create new tables to keep track of merger
        self._create_merge_tables()
        # delete incompatible tables and sequences for main db
        self._remove_screening_tables()
        # Drop user-created indices on tables that grow during merge
        for table in ("Results", "Interactions"):
            indices = self.db_query(
                f"SELECT name FROM sqlite_master WHERE type='index' AND tbl_name='{table}' AND name NOT LIKE 'sqlite_%'"
            ).fetchall()
            for (index_name,) in indices:
                self.db_query(f"DROP INDEX IF EXISTS {index_name}")


    def _get_merge_id(self, mergingdb_path: str) -> int:
        """
        Gets the merge id for the databse in given path

        Args:
            mergingdb_path (str): path to merging (secondary) database

        Returns:
            int: merge id returend by database
        """
        cur = self.conn.cursor()
        cur.execute(
            "INSERT INTO merged_tables(dbfile) VALUES (?)",
            [mergingdb_path],
        )
        return cur.lastrowid

    def _create_merge_tables(self):
        """
        Creates tables necessary when merging two or more databases

        Raises:
            StorageError
        """
        try:
            cur = self.conn.cursor()
            for sql in build_create_table(
                MERGED_TABLES_SCHEMA.name, MERGED_TABLES_SCHEMA, self.dialect
            ):
                cur.execute(sql)
            for sql in build_create_table(
                PK_CONVERSIONS_SCHEMA.name, PK_CONVERSIONS_SCHEMA, self.dialect
            ):
                cur.execute(sql)
            cur.execute(
                "CREATE INDEX IF NOT EXISTS ak_merge ON PK_conversions(merge_id, original_pk)"
            )
        except Exception as e:
            raise StorageError(e) from e




    def _sync_auto_increment_state(self):
        """Reset all sequences to current MAX values in their tables"""
        pass

    # endregion

    # region Methods for dealing with bookmarks and filtering

    def _create_ligname_temp_table(self, csv_path: str):
        """
        Reads ligand names from a CSV and loads them into a temporary table for joining.

        Args:
            csv_path (str)

        """
        import csv

        with open(csv_path, newline="") as f:
            names = [row[0].strip() for row in csv.reader(f) if row]
            self.db_query("DROP TABLE IF EXISTS tmp_lignames")
            self.db_query("CREATE TEMP TABLE tmp_lignames (ligandname TEXT)")
            self.db_update("INSERT INTO tmp_lignames VALUES (?)", [(n,) for n in names])


    # endregion

    # region data
   
    @staticmethod
    def _serialize_pose_coordinates(coords):
        """SQLite has no array type, so pack the [x, y, z] floats into a float32 BLOB
        (~5x smaller and ~65x faster to read than the old JSON text)."""
        if coords is None:
            return None
        return np.asarray(coords, dtype=np.float32).tobytes()

    @staticmethod
    def _deserialize_pose_coordinates(stored):
        """Unpack a float32 BLOB into a list of [x, y, z] floats. Tolerates legacy
        JSON-text values from databases written before the BLOB format."""
        if stored is None:
            return None
        if isinstance(stored, (bytes, bytearray)):
            return np.frombuffer(stored, dtype=np.float32).reshape(-1, 3).tolist()
        return json.loads(stored)

    # endregion

    # region general database operations

    def _begin_transaction(self):
        """
        Begin a transaction
        """
        self.conn.execute("BEGIN TRANSACTION;")

    def _rollback(self):
        """
        Roll back transaction, if one is open
        """
        self.conn.rollback()

    def tables_in_db(self, temp: bool = False) -> list:
        """
        Returns a list of all table names in the database

        Returns:
            list: list of table names
        """
        if temp:
            return [
                name[0].lower()
                for name in self.db_query(
                    "SELECT name FROM sqlite_temp_master WHERE type='table';"
                ).fetchall()
            ]

        tables = [
            name[0].lower()
            for name in self.db_query(
                "SELECT name FROM sqlite_master WHERE type='table';"
            ).fetchall()
        ]
        if "sqlite_sequence" in tables:
            tables.remove("sqlite_sequence")
        return tables

    def clone(self, backup_name: str = None) -> str:
        """Creates a copy of the db

        Args:
            backup_name (str, optional): name of the cloned database

        Returns:
            str: path of backed up database
        """
        if backup_name is None:
            backup_name = self.db_file + ".bk"
        bck = sqlite3.connect(backup_name)
        with bck:
            self.conn.backup(bck, pages=1)
        bck.close()
        logger.info(f"Database {self.db_file} was backed up to {backup_name}.")
        return backup_name

    def _set_ringtail_db_schema_version(self, db_version: str = "3.0.0"):
        """Legacy stamp PRAGMA user_version. Only the intermediate steps of the v1/v2
        upgrade chain still do this; ringtail_schema_version is authoritative.
        """
        self.conn.execute(f"PRAGMA user_version = {db_version.replace('.', '')}")

    def _legacy_pragma_version(self) -> str:
        """Dotted version of a Ringtail 1.x/2.x database from PRAGMA user_version, e.g. '1.0.0'.

        Raises:
            StorageError: if the database is not a Ringtail 1.x/2.x database
        """
        cur = self.conn.cursor()
        db_version = str(cur.execute("PRAGMA user_version").fetchone()[0])
        if db_version == "0":
            # 1.0.0 never set user_version
            results_cols = {
                row[1].lower()
                for row in cur.execute("PRAGMA table_info(Results)").fetchall()
            }
            if "energies_binding" in results_cols:
                db_version = "100"
            else:
                cur.close()
                raise StorageError(
                    f"{self.db_file} is not a Ringtail 1.x or 2.x database, so it cannot be upgraded."
                )
        cur.close()
        return ".".join([*db_version])

    def convert_pose_coordinates_to_native(self, chunk_size: int = 10000):
        """Upgrade an existing database whose pose_coordinates is stored as JSON text
        (VARCHAR) to a packed float32 BLOB (~5x smaller, ~65x faster to read).
        Idempotent: a no-op if the column is already a BLOB. Processed in chunks so a
        multi-GB column never loads fully into memory. Run once per DB.
        """
        sample = self.conn.execute(
            "SELECT pose_coordinates FROM Results WHERE pose_coordinates IS NOT NULL LIMIT 1"
        ).fetchone()
        if sample is not None and isinstance(sample[0], (bytes, bytearray)):
            logger.info("pose_coordinates is already a packed BLOB; nothing to do.")
            return
        logger.info("Converting pose_coordinates VARCHAR -> packed float32 BLOB ...")
        self.conn.execute(
            "ALTER TABLE Results ADD COLUMN pose_coordinates__native BLOB"
        )
        # snapshot ids first so we never read and write the same statement concurrently
        pose_ids = [
            r[0] for r in self.conn.execute("SELECT pose_id FROM Results").fetchall()
        ]
        for i in range(0, len(pose_ids), chunk_size):
            chunk = pose_ids[i : i + chunk_size]
            ph = ",".join("?" * len(chunk))
            rows = self.conn.execute(
                f"SELECT pose_id, pose_coordinates FROM Results WHERE pose_id IN ({ph})",
                chunk,
            ).fetchall()
            updates = [
                (
                    (
                        self._serialize_pose_coordinates(json.loads(txt))
                        if txt is not None
                        else None
                    ),
                    pid,
                )
                for pid, txt in rows
            ]
            self.conn.executemany(
                "UPDATE Results SET pose_coordinates__native = ? WHERE pose_id = ?",
                updates,
            )
        self.conn.execute("ALTER TABLE Results DROP COLUMN pose_coordinates")
        self.conn.execute(
            "ALTER TABLE Results RENAME COLUMN pose_coordinates__native TO pose_coordinates"
        )
        self.conn.commit()
        logger.info("pose_coordinates converted. Run VACUUM to reclaim space.")

    def update_database_version(self, new_version, backup=False):
        """Upgrades a Ringtail 1.x/2.x database step by step to new_version, restoring it if a step fails.

        Args:
            new_version (str): target version, "1.1.0", "2.0.0" or "3.0.0"
            backup (bool, optional): keep the copy taken before upgrading. Defaults to False.
        """
        steps = [
            ("1.0.0", "1.1.0", self._update_db_100_to_110),
            ("1.1.0", "2.0.0", self._update_db_110_to_200),
            ("2.0.0", "3.0.0", self._update_db_200_to_300),
        ]
        self.conn = self._create_connection()
        try:
            try:
                current = self._read_schema_version()
            except DatabaseQueryError:
                current = None
            if current:
                logger.info(
                    f"{self.db_file} is already at schema version {current}, nothing to upgrade."
                )
                return

            version = self._legacy_pragma_version()
            logger.info(
                f"Upgrading {self.db_file} of version {version} to version {new_version}:"
            )
            if version in ["1.0.0", "1.1.0"]:
                logger.warning(
                    "If you created the database with the duplicate handling option, there is a chance of inconsistent behavior of anything involving interactions as the pose_id was not used as an explicit foreign key in db v1.0.0 and v1.1.0."
                )
            backup_file = self.clone()
            try:
                for start, end, step in steps:
                    if version == start and Version(end) <= Version(new_version):
                        step()
                        version = end
                        logger.info(f"Successfully upgraded to {end}.")
            except Exception as e:
                self.conn.close()
                shutil.copy(backup_file, self.db_file)
                raise StorageError(
                    f"Upgrading {self.db_file} failed, it was restored to its original version: {e}"
                ) from e
            finally:
                if not backup:
                    os.remove(backup_file)
        finally:
            self.conn.close()

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
            "CREATE INDEX IF NOT EXISTS ak_results ON Results(LigName, docking_score, leff, deltas, reference_rmsd, energies_inter, energies_vdw, energies_electro, energies_intra, nr_interactions, run_number, pose_rank, num_hb)"
        )
        cur.execute(
            "CREATE INDEX IF NOT EXISTS ak_intind ON Interaction_indices(interaction_type, rec_chain, rec_resname, rec_resid, rec_atom, rec_atomid)"
        )
        try:
            cur.close()
            self._set_ringtail_db_schema_version("1.1.0")
            self.conn.commit()
        except sqlite3.OperationalError as e:
            raise DatabaseConnectionError(
                f"Error while updating database from v1.0.0 to v1.1.0: {e}"
            ) from e

    def _update_db_110_to_200(self):
        """
        Method to update from database v 1.1.0 to 2.0.0,mainly removes the bitvetor table and creates Interactions table
        where interaction just lists pose_id and interaction_id in a long-skinny table

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
                """INSERT INTO Interactions (pose_id, Interaction_id) VALUES (?,?)""",
                pose_indices,
            )
            # drop old bitvector table
            cur.execute("""DROP TABLE IF EXISTS Interaction_bitvectors;""")
            # only Interactions is indexable yet, the 3.0.0 step indexes the rest
            for sql in build_create_indices(
                INTERACTIONS_SCHEMA.name, INTERACTIONS_SCHEMA
            ):
                cur.execute(sql)
            self._set_ringtail_db_schema_version("2.0.0")  # set explicit version
            self.conn.commit()
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
        - giving ligands a ligand_id which is used in Results instead of ligname
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
        self._create_status_tables()
        # drop indices
        indices = self.db_query(
            "SELECT name FROM sqlite_master WHERE type == 'index'"
        ).fetchall()
        for index in indices:
            index_name = index[0]
            try:
                self.db_query(f"DROP INDEX {index_name}")
            except Exception:
                pass
        # create new, empty ligands table
        self._create_ligands_table("Ligands_new")

        # create a temp connection function
        def _smiles_to_rdbin(
            smiles: str,
            h_parents: list[str | int],
            coordinates: list,
            index_map: list[str],
            ligname: str,
        ):
            """
            Temporary db connection method that will use rdkit to convert smiles to Mol
            inline in the sql query. Produces a mol with atoms ordered to match the
            PDBQT coordinate list, so that sequential coordinate assignment works.

            Args:
                smiles (str): smiles describing ligand
                h_parents: JSON list of 1-indexed pairs [parent_mol_idx, h_pdbqt_idx, ...]
                coordinates: JSON list of [x, y, z] coordinate lists (one per PDBQT atom)
                index_map: JSON list of 1-indexed pairs [mol_idx, pdbqt_idx, ...]
                ligname: ligand name for error reporting

            Returns:
                blob: binary Chem.rdchem.Mol ready to insert in db
            """
            try:
                import json
                import traceback

                h_parents_typed = [int(hp) for hp in json.loads(h_parents)]
                coordinates = [
                    [float(pt) for pt in coord] for coord in json.loads(coordinates)
                ]
                original_mol = Chem.MolFromSmiles(smiles)
                if original_mol is None:
                    return None

                index_map = [int(idx) for idx in json.loads(index_map)]
                n_heavy = original_mol.GetNumAtoms()
                n_coords = len(coordinates)

                # Build heavy atom mapping: smiles_mol_idx (0-indexed) -> pdbqt_pos (0-indexed)
                heavy_mol_to_pdbqt = {}
                for i in range(0, len(index_map), 2):
                    mol_idx = index_map[i] - 1
                    pdbqt_idx = index_map[i + 1] - 1
                    heavy_mol_to_pdbqt[mol_idx] = pdbqt_idx

                if n_coords <= n_heavy or not h_parents_typed:
                    # No H in coordinates, just reorder heavy atoms by PDBQT position
                    mapped = sorted(heavy_mol_to_pdbqt.items(), key=lambda x: x[1])
                    newOrder = [mol_idx for mol_idx, _ in mapped]
                    mol = Chem.RenumberAtoms(original_mol, newOrder)
                    mol.RemoveAllConformers()
                    return mol.ToBinary(
                        Chem.PropertyPickleOptions.MolProps
                        | Chem.PropertyPickleOptions.AtomProps
                        | Chem.PropertyPickleOptions.BondProps
                    )

                # H atoms present in coordinates — need to include them in PDBQT order
                # Step 1: Add ALL H to the SMILES mol (preserves heavy atom indices)
                mol_with_h = Chem.AddHs(original_mol)

                # Step 2: Find which mol indices the polar H got, using same
                # logic as meeko's add_hydrogens (iterate parent's H neighbors,
                # pick first unused)
                used_h = set()
                h_mol_to_pdbqt = {}
                num_hydrogens = len(h_parents_typed) // 2
                for i in range(num_hydrogens):
                    parent_mol_idx = h_parents_typed[2 * i] - 1  # 0-indexed
                    h_pdbqt_idx = h_parents_typed[2 * i + 1] - 1  # 0-indexed
                    parent_atom = mol_with_h.GetAtomWithIdx(parent_mol_idx)
                    candidate_hydrogens = [
                        atom.GetIdx()
                        for atom in parent_atom.GetNeighbors()
                        if atom.GetAtomicNum() == 1
                    ]
                    for h_idx in candidate_hydrogens:
                        if h_idx not in used_h:
                            used_h.add(h_idx)
                            h_mol_to_pdbqt[h_idx] = h_pdbqt_idx
                            break

                # Step 3: Combine mappings and renumber so ALL atoms are in PDBQT order
                all_mol_to_pdbqt = {**heavy_mol_to_pdbqt, **h_mol_to_pdbqt}
                # Mapped atoms sorted by PDBQT position (these go first)
                mapped_atoms = sorted(all_mol_to_pdbqt.items(), key=lambda x: x[1])
                # Unmapped atoms (non-polar H without PDBQT coords) go last
                unmapped_atoms = [
                    i
                    for i in range(mol_with_h.GetNumAtoms())
                    if i not in all_mol_to_pdbqt
                ]
                newOrder = [mol_idx for mol_idx, _ in mapped_atoms] + unmapped_atoms
                mol_reordered = Chem.RenumberAtoms(mol_with_h, newOrder)

                # Step 4: Remove unmapped H (no PDBQT coordinates) from the end
                n_mapped = len(mapped_atoms)
                n_total = mol_reordered.GetNumAtoms()
                if n_total > n_mapped:
                    rwmol = Chem.RWMol(mol_reordered)
                    for idx in range(n_total - 1, n_mapped - 1, -1):
                        rwmol.RemoveAtom(idx)
                    Chem.SanitizeMol(rwmol)
                    mol_final = rwmol.GetMol()
                else:
                    mol_final = mol_reordered

                mol_final.RemoveAllConformers()
                mol_final.SetProp("_Name", ligname)
                return mol_final.ToBinary(
                    Chem.PropertyPickleOptions.MolProps
                    | Chem.PropertyPickleOptions.AtomProps
                    | Chem.PropertyPickleOptions.BondProps
                )
            except Exception as e:
                traceback.print_exc()
                raise

        self.conn.create_function("smile_to_rdbin", 5, _smiles_to_rdbin)
        # populate with data from original ligands table, will autogenerate ligand_id PK
        self.db_query(
            """INSERT INTO Ligands_new (
                ligname,
                ligand_smile,
                rdmol) 
            SELECT 
                ligname,
                ligand_smile,
                smile_to_rdbin(
                    ligand_smile,
                    hydrogen_parents,
                    (SELECT ligand_coordinates FROM Results 
                        WHERE Results.ligname = Ligands.ligname LIMIT 1),
                    atom_index_map,
                    ligname
                    )
                FROM Ligands;""",
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
        self._delete_table(LIGANDS_SCHEMA.name)
        # rename new table
        self.db_query("ALTER TABLE Ligands_new RENAME TO Ligands;", commit=True)

        # update results table to use ligand_id from Ligands
        try:
            self.db_query("ALTER TABLE Results ADD COLUMN ligand_id INTEGER;")
        except DatabaseQueryError:
            pass
        # populate ligand_id in Results
        self.db_query("""UPDATE Results
                        SET ligand_id = (
                            SELECT ligand_id FROM Ligands
                            WHERE Ligands.ligname = Results.ligname);""")

        # create new Results table without ligname column
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
                        delta,
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
                        cluster_size,
                        pose_coordinates,
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
                        energies_flexlig,
                        energies_flexlr,
                        energies_intra,
                        energies_torsional,
                        unbound_energy,
                        nr_interactions,
                        num_hb,
                        cluster_size,
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
        self.db_query("ALTER TABLE Results_new RENAME TO Results;")
        # coordinates were copied as JSON text
        self.convert_pose_coordinates_to_native()

        # v1 has no db_properties, v2 stored number_of_poses as text ('all' for every pose)
        # and, before 2.3.0, the docking mode alias the user typed
        if "db_properties" in self.tables_in_db():
            self.db_query("ALTER TABLE db_properties RENAME TO db_properties_v2")
            self._create_db_properties_table()
            self.db_query(
                "INSERT INTO db_properties (db_write_session, docking_mode, number_of_poses) "
                "SELECT DB_write_session, "
                "CASE LOWER(docking_mode) WHEN 'dlg' THEN 'adgpu' WHEN 'gpu' THEN 'adgpu' "
                "WHEN 'pdbqt' THEN 'vina' ELSE LOWER(docking_mode) END, "
                "CASE WHEN number_of_poses = 'all' THEN -1 ELSE CAST(number_of_poses AS INTEGER) END "
                "FROM db_properties_v2"
            )
            self._delete_table("db_properties_v2")
        else:
            self._create_db_properties_table()
        # columns that v1/v2 Receptors lack
        receptor_cols = {
            row[1].lower()
            for row in self.db_query("PRAGMA table_info(Receptors)").fetchall()
        }
        for col in ("flexres_atomnames", "polymer"):
            if col not in receptor_cols:
                self.db_query(f"ALTER TABLE Receptors ADD COLUMN {col} VARCHAR")

        # build new indices if you got this far successfully
        self._create_indices()
        self.db_query("REINDEX")
        self._create_schema_version_table()
        self._stamp_schema_version("3.0.0")
        self.conn.commit()

    def _drop_views(self, db_alias: str = None):
        """
        Legacy method for schema  =< 2. Will drop views and clear bookmark table

        Args:
            db_alias (str, optional): if needing to drop views from a connected, aliased database. Defaults to None.
        """
        if "bookmarks" not in self.tables_in_db():
            return
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
            cursor = con.execute("PRAGMA synchronous = OFF;")
            cursor.execute("PRAGMA journal_mode = MEMORY;")
            cursor.execute("PRAGMA cache_size = -65536;")  # 64 MB page cache
            cursor.execute("PRAGMA temp_store = MEMORY;")
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
        logger.debug(f"Closing database connection {self.conn}")
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
            self.detach_db(attached_db_alias)
        if reindex:
            self.db_query("REINDEX", commit=True)
        # vacuum database
        if vacuum:
            self._vacuum()

    def _vacuum(self):
        """SQLite vacuum rebuilds the database file, repacking it into a minimal amount of disk space

        Raises:
            DatabaseInsertionError
        """
        try:
            # sqlite refuses to VACUUM inside an open transaction
            self.conn.commit()
            cur = self.conn.cursor()
            cur.execute("VACUUM;")
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

    def check_attached(self) -> list[tuple]:
        """
        Check what databases are attached

        Returns:
            list[tuple]: attached database names
        """
        return self.db_query("PRAGMA database_list;").fetchall()

    def detach_db(self, new_db_alias):
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

    def db_update(self, query: str, parameters: list[tuple], commit=True):
        """
        A db query that also commits if/when specified

        Args:
            query (str): sqlite formatted query string
            parameters (list[tuple]): assumes appropriate place holders in query
            commit (bool, optional): whether to commit the transaction in open connection. Defaults to True.

        Raises:
            DatabaseInsertionError

        """
        try:
            self.conn.executemany(query, parameters)
            if commit:
                self.conn.commit()
        except sqlite3.OperationalError as e:
            raise DatabaseInsertionError(f"Error while committing insert query") from e

    # endregion

    # region GUI specific API

    def pose_row_in_table(self, table: str, pose_id: int) -> int |None:
        """
        Find the row id of a pose in a given table, for a bookmark its pose_id (Filtered_poses has no rowid)

        Args:
            table (str)
            pose_id (int)

        Returns:
            int | None: rowid if any
        """
        query = self.QueryBuilder()
        if self.is_bookmark(table):
            # Filtered_poses is WITHOUT ROWID, pose_id is unique within a bookmark
            query.SELECT("pose_id").FROM("Filtered_poses").WHERE(
                "filter_id = (SELECT filter_id from Filters WHERE name = ?)", table
            )
        elif self._is_candidates_table(table):
            query.SELECT("rowid").FROM("Results").WHERE(f"pose_id IN {CANDIDATES_SUBQ}")
        else:
            query.SELECT("rowid").FROM(table)
        query.WHERE("pose_id = ?", pose_id)
        row = self.db_query(*query.build()).fetchone()
        if row:
            return row[0]
        else:
            return None




    # endregion
