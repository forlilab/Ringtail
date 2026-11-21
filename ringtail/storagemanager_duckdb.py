#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail storage adaptor for DuckDB
#

import re
import pandas as pd
from .logutils import get_logger

logger = get_logger(__name__)
from .exceptions import (
    StorageError,
    DatabaseInsertionError,
    DatabaseConnectionError,
    DatabaseQueryError,
    MergeError,
)
from .clustermanager import *
from .storagemanager import StorageManager
from .querybuilder import QueryBuilderDuck
from .schema import (
    build_create_table,
    CANDIDATES_SUBQ,
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
)

try:
    import duckdb

    HAS_DUCK = True
except ImportError:
    HAS_DUCK = False


class StorageManagerDuckDB(StorageManager):
    """DuckDB-specific StorageManager subclass

    Attributes:
        db_file (str): database name
        conn (duckdb.DuckDBPyConnection): Connection to database
    """

    QueryBuilder = QueryBuilderDuck

    def __init__(
        self,
        db_file: str = None,
    ):
        self.db_file = db_file
        super().__init__()
        self.conn: duckdb.DuckDBPyConnection
        self.dialect = "duckdb"

    # region Methods for creating and inserting into tables the database

    def _insert_ligands(self, ligands: list):
        """
        Writes a list of ligands to the Ligands table by first writing
        to a pandas dataframe, which is then registered in the database
        connection, and used to write to the main table.

        Args:
            ligands (list): list of list of expected ligand info
                in expected order
        """

        df = pd.DataFrame(
            ligands,
            columns=[
                "ligname",
                "ligand_smile",
                "rdmol",
            ],
        )
        self.conn.register("df_view", df)
        # to insert interaction if unique
        sql_insert = """
        INSERT INTO Ligands (
            ligname,
            ligand_smile,
            rdmol
            ) 
        SELECT * FROM df_view
        ON CONFLICT(ligname) DO NOTHING;"""

        self.conn.execute(sql_insert)

    def _create_temporary_results_tables(
        self, custom_columns: dict[str, type] | None = None
    ):
        """
        Creates temporary tables for results and interactions, which will be
        used for staging incoming data.
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
            pose_coordinates         FLOAT[][],
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
        # create temporary tables
        self.db_query(create_temp_results)
        self.db_query(create_temp_interactions)

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
        results_array: list,
        interactions_array: list,
        custom_data: list | None = None,
    ):
        """
        Inserts docking results and interactions into their respective temporary
        tables via a pandas dataframe

        Args:
            results_array (list): list of result rows
            interactions_array (list): list of interaction rows
        """
        # insert results rows
        res_df = pd.DataFrame(
            results_array,
        )
        self.conn.register("incoming_poses", res_df)
        temp_insert = f"""
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
                SELECT receptor,
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
                    CAST(pose_coordinates AS FLOAT[][]),
                    flexible_res_coordinates,
                    ligname FROM incoming_poses;"""

        self.db_query(temp_insert)

        int_df = pd.DataFrame(
            interactions_array,
            columns=[
                "ligname",
                "run_number",
                "pose_rank",
                "type",
                "chain",
                "residue",
                "resid",
                "recname",
                "recid",
            ],
        )
        self.conn.register("incoming_interaction", int_df)
        temp_interaction_insert = """
            INSERT INTO Interactions_temp (ligname, run_number, pose_rank, interaction_type, rec_chain, rec_resname, rec_resid, rec_atom, rec_atomid)
            SELECT df.ligname, df.run_number, df.pose_rank, df.type, df.chain, df.residue, df.resid, df.recname, df.recid
            FROM incoming_interaction AS df;"""
        self.db_query(temp_interaction_insert)

        if custom_data:
            # dynamically create columns
            custom_df = pd.DataFrame(custom_data)
            columns = [self._safe_quoting(col) for col in custom_df.columns]
            self.conn.register("incoming_custom_data", custom_df)
            temp_custom_insert = f"""
                INSERT INTO Custom_temp ({",".join(columns)})
                SELECT {",".join(["df." + col for col in columns])}
                FROM incoming_custom_data AS df;"""
            self.db_query(temp_custom_insert)

    def _move_tempresults_to_database(self):
        """
        Inserts data from the temporary results tables to their permanent
        database equivalents
        """
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
            JOIN Ligands AS L ON L.ligname = T.ligname;"""

        temp_to_interaction = """
            INSERT INTO Interactions(pose_id, interaction_id)
            SELECT M.pose_id, II.interaction_id
            FROM Interactions_temp AS IT
            JOIN pose_map AS M
                ON IT.pose_rank = M.pose_rank
                AND IT.run_number = M.run_number
            JOIN Ligands AS L
                ON M.ligand_id = L.ligand_id
                AND IT.ligname = L.ligname
            JOIN Interaction_indices AS II
                ON II.interaction_type = IT.interaction_type
                AND II.rec_chain = IT.rec_chain
                AND II.rec_resname = IT.rec_resname
                AND II.rec_resid = IT.rec_resid
                AND II.rec_atom = IT.rec_atom
                AND II.rec_atomid = IT.rec_atomid;"""

        create_temp_mapping_table = """
        CREATE TEMP TABLE IF NOT EXISTS pose_map(
            pose_id             INTEGER,
            ligand_id           INTEGER,
            run_number          INTEGER,
            pose_rank           INTEGER);
        """
        self.db_query(create_temp_mapping_table)
        cur_max = self.conn.execute(
            "SELECT COALESCE(MAX(pose_id), 0) FROM Results"
        ).fetchone()[0]
        self.db_query(temp_to_results)
        self.db_query(
            "INSERT INTO pose_map (pose_id, ligand_id, run_number, pose_rank) "
            "SELECT pose_id, ligand_id, run_number, pose_rank FROM Results "
            "WHERE pose_id > ?",
            [cur_max],
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
        int_df = pd.DataFrame(
            interactions,
            columns=[
                "pose_id",
                "type",
                "chain",
                "residue",
                "resid",
                "recname",
                "recid",
            ],
        )
        self.conn.register("incoming_interaction", int_df)
        temp_interaction_insert = """
            INSERT INTO Interactions_temp (pose_id, interaction_type, rec_chain, rec_resname, rec_resid, rec_atom, rec_atomid)
            SELECT df.pose_id, df.type, df.chain, df.residue, df.resid, df.recname, df.recid
            FROM incoming_interaction AS df;"""

        self.db_query(temp_interaction_insert)

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
        updates_df = pd.DataFrame(data, columns=["pose_id", "num_int", "num_hb"])
        self.conn.register("updates", updates_df)

        self.db_query("""
            UPDATE Results
            SET num_interactions = u.num_int,
                num_hb = u.num_hb
            FROM updates u
            WHERE Results.pose_id = u.pose_id;
        """)




    def _insert_interaction_index_rows(self, interactions: list[dict]):
        """
        Writes unique interactions to database via a pandas dataframe

        Args:
            interactions (list[dict]): [(interaction_type, rec_chain, rec_resname, rec_resid, rec_atom, rec_atomid)]
        """
        df = pd.DataFrame(
            [self._interaction_index_fields(r) for r in interactions],
            columns=["type", "chain", "residue", "resid", "recname", "recid"],
        ).drop_duplicates()
        self.conn.register("df_view", df)
        # to insert interaction if unique
        alt_sql_insert = """
                    INSERT INTO Interaction_indices(interaction_type,rec_chain,rec_resname,rec_resid,rec_atom,rec_atomid) 
                    SELECT * FROM df_view
                    ON CONFLICT DO NOTHING;
                    """
        self.db_query(alt_sql_insert)






    def _delete_nontables(self):
        """
        Deletes objects in the database that are not tables (handled separately)
        """
        sequences = self.db_query("""SELECT sequence_name FROM duckdb_sequences();""")
        if sequences:
            for sequence in sequences.fetchall():
                self.db_query(f"DROP SEQUENCE {sequence[0]};")

    def _remove_screening_tables(self):
        '''
        Deletes all data included in screening such as clusters, filtered data, status assignments, and comments
        '''
        self._delete_table(POSE_CLUSTERS_SCHEMA.name)
        self._delete_table(CLUSTER_GROUPS_SCHEMA.name)
        self._delete_table(CLUSTERS_SCHEMA.name)
        self.db_query("DROP SEQUENCE IF EXISTS seq_clusters_cluster_id;")
        self._delete_table(FILTERED_POSES_SCHEMA.name)
        self._delete_table(FILTERS_SCHEMA.name)
        self.db_query("DROP SEQUENCE IF EXISTS seq_filters_filter_id;")
        self._delete_table("Accepted")
        self._delete_table("Maybe")
        self._delete_table("Rejected")
        self._delete_table("Pose_comments")

    # endregion

    # region merge databases
    def prepare_for_merging(self):
        """
        Prepares for merging by dropping sequences and tables incompatible with merging
        """
        # create new tables to keep track of merger
        self._create_merge_tables()
        # delete incompatible tables and sequences for main db
        self._remove_screening_tables()
        # sync sequences to current max values, critical for databases
        # created via export_bookmark_db where data was inserted with
        # explicit IDs but sequences were initialized starting at 1
        self._sync_auto_increment_state()


    def _create_merge_tables(self):
        """
        Creates tables necessary when merging two or more databases

        Raises:
            StorageError
        """
        try:
            for sql in build_create_table(
                MERGED_TABLES_SCHEMA.name, MERGED_TABLES_SCHEMA, self.dialect
            ):
                self.conn.execute(sql)
            for sql in build_create_table(
                PK_CONVERSIONS_SCHEMA.name, PK_CONVERSIONS_SCHEMA, self.dialect
            ):
                self.conn.execute(sql)

            # If table already has data, reset sequence to continue from max
            max_id = self.conn.execute(
                "SELECT COALESCE(MAX(merge_id), 0) FROM merged_tables"
            ).fetchone()[0]
            if max_id > 0:
                self.conn.execute("DROP SEQUENCE seq_merged_tables_merge_id;")
                self.conn.execute(
                    f"CREATE SEQUENCE seq_merged_tables_merge_id START {max_id + 1};"
                )

        except Exception as e:
            raise MergeError(e) from e




    def _sync_auto_increment_state(self):
        """Advance sequences that lag behind the MAX id in their table."""
        for table, column in (
            (DB_PROPERTIES_SCHEMA.name, "DB_write_session"),
            (LIGANDS_SCHEMA.name, "ligand_id"),
            (RESULTS_SCHEMA.name, "pose_id"),
            (INTERACTION_INDICES_SCHEMA.name, "interaction_id"),
            (INTERACTIONS_SCHEMA.name, "interaction_pose_id"),
        ):
            default = self.conn.execute(
                "SELECT column_default FROM information_schema.columns "
                "WHERE lower(table_name) = lower(?) AND lower(column_name) = lower(?)",
                (table, column),
            ).fetchone()
            match = re.search(r"nextval\('([^']+)'\)", (default or [""])[0] or "")
            if not match:
                continue
            seq_name = match.group(1)
            max_val = self.conn.execute(
                f"SELECT COALESCE(MAX({column}), 0) FROM {table}"
            ).fetchone()[0]
            self.conn.execute("BEGIN TRANSACTION")
            next_val = self.conn.execute(f"SELECT nextval('{seq_name}')").fetchone()[0]
            # a sequence in use by a column default cannot be restarted, only advanced
            if next_val < max_val:
                self.conn.execute(
                    f"SELECT MAX(nextval('{seq_name}')) FROM range({max_val - next_val})"
                )
            # nextval is only persisted by a transaction that also writes
            self.conn.execute(f"COMMENT ON SEQUENCE {seq_name} IS NULL")
            self.conn.execute("COMMIT")

    # endregion

    # region Methods for dealing with bookmarks and filtering

    def _create_ligname_temp_table(self, csv_path: str):
        """
        Loads ligand names from CSV into a temporary table using DuckDB's native CSV reader.

        Args:
            csv_path (str)

        """
        self.conn.execute("DROP TABLE IF EXISTS tmp_lignames")
        self.conn.execute(
            f"CREATE TEMP TABLE tmp_lignames AS "
            f"SELECT CAST(column0 AS VARCHAR) AS ligandname FROM read_csv('{csv_path}', header=false)"
        )


    # endregion

    # region general database operations

    def _begin_transaction(self):
        """
        Currently not implemented in duckdb
        """
        pass

    def _rollback(self):
        """
        No-op, since _begin_transaction does not open a transaction
        """
        pass

    def tables_in_db(self, temp: bool = False) -> list:
        """
        Returns a list of all table names in the database

        Returns:
            list: list of table names
        """
        return [row[0].lower() for row in self.db_query("SHOW TABLES").fetchall()]

    def clone(self, backup_name: str = None) -> str:
        """Creates a copy of the db

        Args:
            backup_name (str, optional): name of the cloned database

        Returns:
            str: path of backed up database
        """
        import shutil

        if backup_name is None:
            backup_name = self.db_file + ".bk"
        # Windows refuses to read an open DuckDB file, and closing writes the WAL into it
        self.conn.close()
        try:
            shutil.copy(self.db_file, backup_name)
        finally:
            self.conn = self._create_connection()

        logger.info(f"Database {self.db_file} was backed up to {backup_name}.")
        return backup_name

    def convert_pose_coordinates_to_native(self):
        """No-op for DuckDB databases written by this version: pose_coordinates is
        already a native FLOAT[][]. Only pre-release v3-alpha files hold JSON text,
        and those need a full recreate, since DuckDB cannot ALTER a column on a
        table other tables reference by foreign key.
        """
        dtype = self.conn.execute(
            "SELECT data_type FROM information_schema.columns "
            "WHERE table_name = 'Results' AND column_name = 'pose_coordinates'"
        ).fetchone()
        if dtype and "[]" in dtype[0]:
            logger.info("pose_coordinates is already a native array; nothing to do.")
            return
        raise NotImplementedError(
            "This DuckDB database stores pose_coordinates as JSON text, which only "
            "pre-release v3-alpha files do. Converting requires recreating the "
            "database; ask your package administrator for rt_upgrade_from_v3alpha.py."
        )
    
    def update_database_version(self, new_version, backup=False):
        '''
        duckdb currently up to date with code version major 3

        Args:
            new_version (_type_): _description_
            backup (bool, optional): _description_. Defaults to False.
        '''
        # opening the database checks its schema version
        with self:
            pass
        logger.info(
            f"{self.db_file} database is compatible and no migration was performed."
        )

    def _create_connection(self):
        """Creates database connection to self.db_file

        Returns:
            duckdb.DuckDBPyConnection: Connection object to self.db_file

        Raises:
            DatabaseConnectionError
        """
        try:
            conn = duckdb.connect(self.db_file)
            conn.execute("SET enable_progress_bar=false")
        except duckdb.OperationalError as e:
            raise DatabaseConnectionError(
                "Error while establishing database connection"
            ) from e
        return conn

    def close_storage(self, attached_db=None, vacuum=None):
        """
        Closes storage

        Args:
            attached_db (str, optional): alias of attached database. Defaults to None.
            vacuum (bool, optional): whether or not to vacuum file to save space. Defaults to None.
        """

        # close db itself
        logger.debug(f"Closing database connection {self.conn}")
        self.conn.close()

    def _cleanup_storage(
        self, attached_db_alias: str = None, vacuum: bool = None, reindex=False
    ):
        """
        Cleans up storage, especially useful for situations where two databases
        have been connected/attached. Will detach the database.

        Args:
            attached_db (str, optional): Name of attached database, if any. Defaults to None.
            vacuum (bool, optional): for sqlite compatibility. duckdb VACUUM only
                rebuilds statistics; compacting a file needs COPY FROM DATABASE.
            reindex (bool, optional): for sqlite compatibility, duckdb manages its own
        """
        if attached_db_alias is not None:
            self.detach_db(attached_db_alias)

    def _attach_db(self, new_db: str, new_db_alias: str = "attached_db") -> str:
        """Attaches new database file to current database

        Args:
            new_db (str): file name for database to attach
            new_db_name (str): name of new database

        Returns:
            str: alias the database was attached as

        Raises:
            StorageError
        """
        attach = f"""ATTACH '{new_db}' AS {new_db_alias};"""
        self.db_query(attach)
        return new_db_alias

    def check_attached(self) -> list:
        return self.db_query("SHOW DATABASES;").fetchall()

    def detach_db(self, new_db_alias):
        """Detaches new database file from current database

        Args:
            new_db_name (str): db name for database to detach

        Raises:
            StorageError
        """
        self.db_query(f"""DETACH {new_db_alias};""")

    def db_query(
        self, query, params: tuple = (), commit=False
    ) -> duckdb.DuckDBPyConnection:
        """Executes provided duckdb query. Returns cursor for results.
            Since cursor remains open, added to list of open cursors

        Args:
            query (str): Formated duckdb query as string

        Returns:
            duckdb cursor: Contains results of query
        """

        try:
            cur = self.conn.execute(query, params)
            if commit:
                self.conn.commit()
        except duckdb.Error as e:
            raise DatabaseQueryError(
                f"Unable to execute query -{query}- with given parameters -{params}-: -{e}-"
            ) from e
        return cur

    def db_update(self, query: str, parameters: list[list], commit=True):
        """
        A db query that also commits if/when specified

        Args:
            query (str): duckdb formatted query string
            parameters (list[list]): assumes appropriate place holders in query
            commit (bool, optional): whether to commit the transaction in open connection. Defaults to True.

        Raises:
            DatabaseInsertionError
        """
        try:
            self.conn.executemany(query, parameters)
            if commit:
                self.conn.commit()
        except duckdb.OperationalError as e:
            raise DatabaseInsertionError(f"Error while committing insert query") from e

    # endregion

    # region GUI specific API





    # endregion
