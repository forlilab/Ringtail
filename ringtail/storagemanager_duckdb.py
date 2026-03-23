#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail storage adaptor for DuckDB
#

import os.path
import json
import pandas as pd
from .logutils import LOGGER as logger
from typing import Union
from importlib.metadata import version
from .util import numlist2str
from .exceptions import (
    StorageError,
    DatabaseInsertionError,
    DatabaseConnectionError,
    DatabaseQueryError,
    OptionError,
    MergeError,
)
from .clustermanager import *
from .storagemanager import StorageManager
from .querybuilder import QueryBuilderDuck

try:
    import duckdb

    HAS_DUCK = True
except:
    HAS_DUCK = False


class StorageManagerDuckDB(StorageManager):
    """DuckDB-specific StorageManager subclass

    Attributes:
        db_file (str): database name
        conn (duckdb.DuckDBPyConnection): Connection to database
    """

    # "db_schema_ver":list("compatible code versions")
    _db_schema_code_compatibility = {
        "3.0.0": ["3.0.0"],
    }
    QueryBuilder = QueryBuilderDuck

    def __init__(
        self,
        db_file: str = None,
    ):
        self.db_file = db_file
        super().__init__()
        self.conn: duckdb.DuckDBPyConnection

    # region Methods for creating and inserting into tables the database

    def _create_ligands_table(self, name="Ligands"):
        """
        Creates ligands table with a primary key defined by a
        sequence

        Args:
            name (str, optional): Defaults to "Ligands".
        """
        ligand_table = f"""
            CREATE SEQUENCE IF NOT EXISTS seq_ligandid START 1;
            CREATE TABLE IF NOT EXISTS {name} (
            ligand_id           INTEGER DEFAULT nextval('seq_ligandid') PRIMARY KEY,
            LigName             VARCHAR NOT NULL UNIQUE,
            ligand_smile        VARCHAR,
            rdmol               BLOB)"""

        self.db_query(ligand_table)

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
                "LigName",
                "ligand_smile",
                "rdmol",
            ],
        )
        self.conn.register("df_view", df)
        # to insert interaction if unique
        sql_insert = """
        INSERT INTO Ligands (
            LigName,
            ligand_smile,
            rdmol
            ) 
        SELECT * FROM df_view
        ON CONFLICT(LigName) DO NOTHING;"""

        self.conn.execute(sql_insert)

    def _create_results_table(self, name="Results"):
        """
        Creates table for results, including a sequence creating the primary key
        pose_id as well as referencing foreign key ligand_id

        Args:
            name (str, optional): _description_. Defaults to "Results".
        """

        sql_results_table = f"""
            CREATE SEQUENCE IF NOT EXISTS seq_poseid START 1;
            CREATE TABLE IF NOT EXISTS {name} (
            Pose_ID             INTEGER DEFAULT nextval('seq_poseid') PRIMARY KEY,
            ligand_id           INTEGER REFERENCES Ligands(ligand_id),
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
            energies_flexLig    FLOAT,
            energies_flexLR     FLOAT,
            energies_intra      FLOAT,
            energies_torsional  FLOAT,
            unbound_energy      FLOAT,
            num_interactions     INTEGER,
            num_hb              INTEGER,
            pose_coordinates         VARCHAR,
            flexible_res_coordinates   VARCHAR,
            ); """

        self.db_query(sql_results_table)

    def _create_temporary_results_tables(self):
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
            energies_flexLig    FLOAT,
            energies_flexLR     FLOAT,
            energies_intra      FLOAT,
            energies_torsional  FLOAT,
            unbound_energy      FLOAT,
            num_interactions     INTEGER,
            num_hb              INTEGER,
            pose_coordinates         VARCHAR,
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

    def _insert_results_in_temp_tables(
        self, results_array: list, interactions_array: list
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
                    energies_flexLig,
                    energies_flexLR,
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
                    energies_flexLig,
                    energies_flexLR,
                    energies_intra,
                    energies_torsional,
                    unbound_energy,
                    num_interactions,
                    num_hb,
                    cluster_size,
                    pose_coordinates,
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
                energies_flexLig,
                energies_flexLR,
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
                T.energies_flexLig,
                T.energies_flexLR,
                T.energies_intra,
                T.energies_torsional,
                T.unbound_energy,
                T.num_interactions,
                T.num_hb,
                T.cluster_size,
                T.pose_coordinates,
                T.flexible_res_coordinates
            FROM Results_temp AS T
            JOIN Ligands AS L ON L.LigName = T.LigName
            RETURNING pose_id, ligand_id, run_number, pose_rank;"""

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

        mapping = self.db_query(temp_to_results).fetchall()

        create_temp_mapping_table = """
        CREATE TEMP TABLE IF NOT EXISTS pose_map(
            pose_id             INTEGER,
            ligand_id           INTEGER,
            run_number          INTEGER,
            pose_rank           INTEGER);
        """

        self.db_query(create_temp_mapping_table)
        self.db_update(
            "INSERT INTO pose_map (pose_id, ligand_id, run_number, pose_rank) VALUES (?, ?, ?, ?);",
            mapping,
            commit=False,
        )
        self.db_query(temp_to_interaction)

        # clear temp tables
        self.db_query("DELETE FROM pose_map")
        self.db_query("DELETE FROM Interactions_temp")
        self.db_query("DELETE FROM Results_temp")

    def _delete_new_duplicate_results(self):
        """Checks if a pose is uniquely represented in the Results table,
        and deletes it from the staged incoming data if duplicated.
        Based on the following columns:
        ligname,
        receptor,
        pose_coordinates
        flexible_res_coordinates
        """
        delete_int_sql = """
        DELETE FROM Interactions_temp
        WHERE EXISTS (
            SELECT 1
            FROM Results_temp AS RT
            JOIN Results AS R
                ON RT.receptor    = R.receptor
                AND RT.pose_coordinates = R.pose_coordinates
                AND RT.flexible_res_coordinates = R.flexible_res_coordinates
            JOIN Ligands AS L
                ON RT.ligname = L.ligname
            WHERE 
                RT.ligname   = Interactions_temp.ligname
                AND RT.pose_rank = Interactions_temp.pose_rank
                AND RT.run_number = Interactions_temp.run_number
        );
        """
        delete_res_sql = """
        DELETE FROM Results_temp AS RT
        WHERE EXISTS (
            SELECT 1
            FROM Results AS R
            JOIN Ligands AS L
                ON RT.ligname = L.ligname
            WHERE 
                RT.receptor  = R.receptor
                AND RT.pose_coordinates = R.pose_coordinates
                AND RT.flexible_res_coordinates = R.flexible_res_coordinates
                AND L.ligand_id    = R.ligand_id);
        """
        self.db_query(delete_int_sql)
        self.db_query(delete_res_sql)

    def _delete_old_duplicate_results(self):
        """Checks if a pose is uniquely represented in the Results table,
        and deletes it from Results if duplicated.
        Based on the following columns:
        ligname,
        receptor,
        pose_coordinates
        flexible_res_coordinates
        """
        delete_sql = """
        WITH target_poseid AS (
            SELECT R.pose_id
            FROM Results AS R
            JOIN Results_temp AS RT
                ON RT.receptor = R.receptor
                AND RT.pose_coordinates = R.pose_coordinates
                AND RT.flexible_res_coordinates = R.flexible_res_coordinates
            JOIN Ligands AS L
                ON RT.ligname = L.ligname
            )
        DELETE FROM Interactions
        WHERE pose_id IN (SELECT pose_id from target_poseid)
        returning pose_id;
        """
        delete_pose_ids = self.db_query(delete_sql).fetchall()
        delete_pose_ids_list = {row[0] for row in delete_pose_ids}
        placeholders = ",".join("?" for _ in delete_pose_ids_list)
        if delete_pose_ids:
            self.db_query(
                f"""DELETE FROM Results WHERE pose_id IN ({placeholders});""",
                delete_pose_ids_list,
            )

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

        self.db_query(
            """
            UPDATE Results
            SET num_interactions = u.num_int,
                num_hb = u.num_hb
            FROM updates u
            WHERE Results.pose_id = u.pose_id;
        """
        )

    def _insert_completed_poses(self, pose_ids: list[tuple], tracking_table: str):
        """
        Inserts processed poses into process tracking table

        Args:
            pose_ids (list[int]): _description_
        """
        query = f"""
        INSERT INTO {tracking_table} (
        pose_id) VALUES (?);
        """
        self.db_update(query, pose_ids, commit=False)

    def _create_receptors_table(self):
        """Create table for receptors. Has primary key although only one receptor allowed"""
        receptors_table = """
            CREATE SEQUENCE IF NOT EXISTS seq_receptorid START 1;
            CREATE TABLE IF NOT EXISTS Receptors (
            Receptor_ID         INTEGER DEFAULT nextval('seq_receptorid') PRIMARY KEY,
            RecName             VARCHAR,
            box_dim             VARCHAR,
            box_center          VARCHAR,
            grid_spacing        FLOAT,
            flexible_residues   VARCHAR,
            flexres_atomnames   VARCHAR,
            receptor_object     BLOB,
            polymer             JSON
        );"""
        self.db_query(receptors_table)

    def _insert_receptors(self, receptor_array: list):
        """Takes array of receptor rows, inserts into Receptors table

        Args:
            receptor_array (list): List of lists
                containing formatted receptor rows
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

        self.db_query(sql_insert, receptor_array, commit=True)

    def insert_receptor_blob(self, receptor: bytes, rec_name: str):
        """Takes object of Receptor class, updates the column in Receptor table

        Args:
            receptor (bytes): bytes receptor object to be inserted into DB
            rec_name (string): Name of receptor. Used to insert into correct row of DB
        """
        # Check if there is already a row for the receptor
        # Check if there is already a row for the receptor
        count = self.table_length("Receptors")

        if count == 0:
            # Insert receptor statement
            query = f"""INSERT INTO Receptors (
                      RecName,
                      receptor_object)
                      VALUES (?,?);"""

        else:
            query = """UPDATE Receptors SET RecName = ?, receptor_object = ? WHERE Receptor_ID == 1"""
        self.db_query(query, (rec_name, receptor), commit=True)

    def insert_receptor_polymer(self, receptor: str, rec_name: str):
        """Takes object of Receptor class, updates the column in Receptor table

        Args:
            receptor (str): json string representation of a receptor meeko.Polymer oobject to be inserted into DB
            rec_name (str): Name of receptor. Used to insert into correct row of DB
        """
        # Check if there is already a row for the receptor
        count = self.table_length("Receptors")

        if count == 0:
            # Insert receptor statement
            query = f"""INSERT INTO Receptors (
                      RecName,
                      polymer)
                      VALUES (?,?);"""

        else:
            query = """UPDATE Receptors SET RecName = ?, polymer = ? WHERE Receptor_ID == 1;"""
        self.db_query(query, (rec_name, receptor), commit=True)

    def _create_db_properties_table(self):
        """Create table of database properties used during write session to the database. Columns are:
        DB_write_session int (primary key)
        docking_mode (vina or adgpu)
        num_of_poses (int, -1 for all, otherwise value)
        """

        sql_str = """
        CREATE SEQUENCE IF NOT EXISTS seq_dbwriteid START 1;
        CREATE TABLE IF NOT EXISTS DB_properties (
        DB_write_session    INTEGER DEFAULT nextval('seq_dbwriteid') PRIMARY KEY,
        docking_mode        VARCHAR,
        number_of_poses     INTEGER);"""

        self.db_query(sql_str)

    def _insert_db_properties(self, docking_mode: str, number_of_poses: str):
        """Insert db properties into database properties table

        Args:
            docking_mode (str): docking mode for the current dataset being written
            number_of_poses (str): number of poses written to database in current session, either "all" or specified max_poses
        """
        sql_insert = """INSERT INTO DB_properties (
        docking_mode,
        number_of_poses
        ) VALUES (?,?);"""
        self.db_query(sql_insert, [docking_mode, number_of_poses], commit=True)

    def _create_interaction_index_table(self):
        """Creates a table describing unique interactions in the database"""
        interaction_index_table = """
        CREATE SEQUENCE IF NOT EXISTS seq_interactionid START 1;
        CREATE TABLE IF NOT EXISTS Interaction_indices (
        interaction_id      INTEGER DEFAULT nextval('seq_interactionid') PRIMARY KEY,
        interaction_type    VARCHAR,
        rec_chain           VARCHAR,
        rec_resname         VARCHAR,
        rec_resid           VARCHAR,
        rec_atom            VARCHAR,
        rec_atomid          VARCHAR,
        UNIQUE (interaction_type, rec_chain, rec_resname, rec_resid, rec_atom, rec_atomid));
        """
        self.db_query(interaction_index_table)

    def _create_interaction_table(self):
        """Creates a table of each pose-interaction combination."""

        interaction_table = """
        CREATE SEQUENCE IF NOT EXISTS seq_interactionposeid START 1;
        CREATE TABLE IF NOT EXISTS Interactions (
        interaction_pose_ID INTEGER DEFAULT nextval('seq_interactionposeid') PRIMARY KEY,
        Pose_ID   INTEGER REFERENCES Results(pose_id),
        interaction_id INTEGER REFERENCES Interaction_indices(interaction_id))"""

        self.db_query(interaction_table)

    def _insert_interaction_index_rows(self, interactions: list[dict]):
        """
        Writes unique interactions to database via a pandas dataframe

        Args:
            interactions (list[dict]): [(interaction_type, rec_chain, rec_resname, rec_resid, rec_atom, rec_atomid)]
        """
        df = pd.DataFrame(
            interactions,
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
        FOREIGN KEY (pose_id) REFERENCES Results(pose_id)
        );"""
        self.db_query(table_sql, commit=True)

    def _create_filtering_tables(self):
        """
        Creates a Filter table which includes filter_id (PK), name (bookmark_name), duckdb formatted query,
        and dictionary of filters used, as well as Filtered_poses, which uses filter_id as FK,
        and lists all poses passing that filter_id
        """
        filters_sql = """
        CREATE SEQUENCE IF NOT EXISTS seq_filterid START 1;
        CREATE TABLE IF NOT EXISTS Filters (
        filter_id           INTEGER DEFAULT nextval('seq_filterid') PRIMARY KEY,
        name                VARCHAR,
        query               VARCHAR,
        filters             VARCHAR,
        filter_window       VARCHAR);"""

        filter_pose_sql = """
        CREATE TABLE IF NOT EXISTS Filtered_poses (
        filter_id           INTEGER REFERENCES Filters(filter_id),
        pose_id             INTEGER REFERENCES Results(pose_id));"""

        self.db_query(filters_sql)
        self.db_query(filter_pose_sql)

    def _create_indices(self):
        """
        Duckdb does not use indexing like other databases
        """
        pass

    def _create_cluster_tables(self):
        """
        Creates cluster tables if they don't already exist
        """
        # create a cluster description table
        self.db_query(
            """

        CREATE SEQUENCE IF NOT EXISTS seq_clusterid START 1;
            CREATE TABLE IF NOT EXISTS
            Clusters (
                cluster_id INTEGER DEFAULT nextval('seq_clusterid') PRIMARY KEY,
                name VARCHAR,
                description VARCHAR,
                cluster_window VARCHAR,
                num_clusters INTEGER
                );
            """
        )
        # create a table with cluster_id and cluster_group and representative pose
        self.db_query(
            """
            CREATE TABLE IF NOT EXISTS
            Cluster_groups (
                cluster_id INTEGER REFERENCES Clusters(cluster_id),
                cluster_group INTEGER,
                representative INTEGER REFERENCES Results(pose_id)
                )
            """
        )
        # create a table of pose and cluster_id and cluster_group
        self.db_query(
            """
            CREATE TABLE IF NOT EXISTS 
            Pose_clusters (
                cluster_id INTEGER REFERENCES Clusters(cluster_id),
                cluster_group INTEGER,
                pose_id INTEGER REFERENCES Results(pose_id)
                );
            """
        )

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
        return self.db_query(
            "SELECT cluster_id FROM Clusters WHERE name=? and cluster_window=?",
            (
                cluster_name,
                cluster_window,
            ),
        ).fetchone()

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
        return self.db_query(
            """
            INSERT INTO Clusters (name, description, cluster_window, num_clusters)
                       VALUES (?,?,?,?) RETURNING cluster_id;""",
            (
                name,
                description,
                cluster_window,
                length,
            ),
        ).fetchone()[0]

    def _insert_clusters(self, cluster_groups: list, pose_rows: list):
        """
        Inserts all cluster data, including each grouping and its representative
        pose, and all poses involved and which clusters they belong to

        Args:
            cluster_groups (list): each cluster from the clustering exercise
            pose_rows (list): cluster_id, group_id, pose_od
        """
        self.db_update(
            """INSERT INTO Cluster_groups VALUES (?,?,?);""",
            cluster_groups,
            commit=False,
        )

        self.db_update(
            """
                INSERT INTO Pose_clusters VALUES (?,?,?);
                """,
            pose_rows,
            commit=False,
        )

    def _delete_cluster(self, cluster_id):
        # delete all poses and cluster groups associated with cluster id
        self.conn.execute(
            """
            DELETE FROM Pose_clusters WHERE cluster_id = ?;""",
            (cluster_id,),
        )
        # delete representative pose id for each cluster group
        self.conn.execute(
            """
            DELETE FROM Cluster_groups WHERE cluster_id = ?;""",
            (cluster_id,),
        )
        # delete the cluster information
        self.conn.execute(
            """
            DELETE FROM Clusters WHERE cluster_id = ?;""",
            (cluster_id,),
        )

    # endregion

    # region merge databases
    def prepare_for_merging(self):
        """
        Prepares for merging by dropping sequences and tables incompatible with merging
        """
        # create new tables to keep track of merger
        self._create_merge_tables()
        # delete incompatible tables and sequences for main db
        self._delete_table("Pose_clusters")
        self._delete_table("Cluster_groups")
        self.db_query("DROP SEQUENCE IF EXISTS seq_clusterid;")
        self._delete_table("Clusters")
        self._delete_table("Filtered_poses")
        self._delete_table("Filters")
        self.db_query("DROP SEQUENCE IF EXISTS seq_filterid;")
        # sync sequences to current max values, critical for databases
        # created via export_bookmark_db where data was inserted with
        # explicit IDs but sequences were initialized starting at 1
        self._sync_auto_increment_state()

    def _get_merge_id(self, mergingdb_path: str) -> int:
        """
        Gets the merge id for the databse in given path

        Args:
            mergingdb_path (str): path to merging (secondary) database

        Returns:
            int: merge id returend by database
        """
        return self.conn.execute(
            """INSERT INTO merged_tables(dbfile) VALUES (?) RETURNING merge_id""",
            [mergingdb_path],
        ).fetchone()[0]

    def _merging_receptors_compatible(self) -> bool:
        """
        Checks if the receptor names in the two databases are a mathc

        Returns:
            bool: True if receptor allows merging
        """
        receptorcheck_sql = """
        SELECT CASE 
            WHEN Receptors.RecName = merging_receptors.RecName THEN 'True'
            ELSE 'False'
        END AS comparison_result 
        FROM Receptors 
        JOIN merging.Receptors AS merging_receptors 
            ON Receptors.receptor_id = merging_receptors.receptor_id;"""
        receptor_comp = self.db_query(receptorcheck_sql).fetchone()
        if receptor_comp:
            if receptor_comp[0] == "True":
                return True
        return False

    def _create_merge_tables(self):
        """
        Creates tables necessary when merging two or more databases

        Raises:
            StorageError
        """
        try:
            # create mergedata table: merge_id (PK), dbfile, timestamp, numofrows table
            mergetbl_sql = """
            CREATE SEQUENCE IF NOT EXISTS seq_mergeid START 1;
            CREATE TABLE IF NOT EXISTS merged_tables (
            merge_id                INTEGER DEFAULT nextval('seq_mergeid') PRIMARY KEY,
            dbfile                  VARCHAR,
            merge_start             TIMESTAMP DEFAULT CURRENT_TIMESTAMP);"""
            self.conn.execute(mergetbl_sql)

            # create PK table: merge_id(FK), table, original_val, merge_val
            pktable_sql = """CREATE TABLE IF NOT EXISTS PK_conversions (
            merge_id        INTEGER REFERENCES merged_tables(merge_id),
            table_name      VARCHAR,
            original_PK     INTEGER,
            merged_PK       INTEGER);"""
            self.conn.execute(pktable_sql)

            # If table already has data, reset sequence to continue from max
            max_id = self.conn.execute(
                "SELECT COALESCE(MAX(merge_id), 0) FROM merged_tables"
            ).fetchone()[0]
            if max_id > 0:
                self.conn.execute("DROP SEQUENCE seq_mergeid;")
                self.conn.execute(f"CREATE SEQUENCE seq_mergeid START {max_id + 1};")

        except Exception as e:
            raise MergeError(e) from e

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
                original_PK,
                merged_PK) SELECT 
                ?,
                'DB_properties', 
                DB_write_session,
                DB_write_session + (SELECT MAX(DB_write_session) FROM db_properties) 
                FROM merging.db_properties;"""
            self.conn.execute(
                convert_dbprop_sql,
                [merge_id],
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
            self.conn.execute(insert_dbprops_sql, [merge_id])

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
            MergeError
        """
        # convert ligand_ids and log
        convert_ligand_ids_sql = """INSERT INTO PK_conversions (
        merge_id,
        table_name,
        original_PK,
        merged_PK) SELECT
        ?,
        'Ligands',
        ligand_id,
            CASE
                WHEN EXISTS (
                    SELECT 1
                    FROM Ligands
                    WHERE
                        merging.Ligands.LigName = Ligands.LigName
                    )
                THEN (
                    SELECT main.Ligands.ligand_id
                    FROM main.Ligands
                    WHERE
                        merging.Ligands.LigName = Ligands.LigName
                    )
                ELSE merging.Ligands.ligand_id + (SELECT MAX(ligand_id) FROM Ligands)
            END AS new_ligand_id
        FROM merging.Ligands;"""

        # then inserting only those that aren't already in the table
        insert_new_ligands = """INSERT INTO Ligands (
        ligand_id,
        LigName,
        ligand_smile,
        rdmol)
        SELECT 
            (SELECT merged_PK FROM PK_conversions WHERE original_PK = ligand_id and merge_id = ? AND table_name = 'Ligands') new_id,
            LigName,
            ligand_smile,
            rdmol
        FROM merging.Ligands WHERE new_id > (SELECT MAX(ligand_id) FROM Ligands);
        """

        # convert pose_id
        convert_poseid_sql = """INSERT INTO PK_conversions (
        merge_id,
        table_name,
        original_PK,
        merged_PK) SELECT 
        ?,
        'Results', 
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
            delta,
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
            num_interactions,
            num_hb,
            pose_coordinates,
            flexible_res_coordinates) 
        SELECT 
            pose.merged_PK as pose_id,
            ligand.merged_PK as ligand_id,
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
            mr.energies_flexLig,
            mr.energies_flexLR,
            mr.energies_intra,
            mr.energies_torsional,
            mr.unbound_energy,
            mr.num_interactions,
            mr.num_hb,
            mr.pose_coordinates,
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

        except Exception as e:
            raise MergeError(
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
        # If main database does not have interactions, it will be incompatible to
        # merge in interactions from other databases
        if not self.table_length("Interactions") > 0:
            return
        convert_ii_sql = """INSERT INTO PK_conversions (
        merge_id,
        table_name,
        original_PK,
        merged_PK) SELECT 
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
        # TODO #BUG I think the issue is here, because it has interacton_pose_id
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
                AND merge_id = ?) II ON (I.Interaction_ID = II.original_PK);"""

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

        return True

    def _sync_auto_increment_state(self):
        """Reset all sequences to current MAX values in their tables"""
        sequence_map = {
            "seq_dbwriteid": ("DB_properties", "DB_write_session"),
            "seq_ligandid": ("Ligands", "ligand_id"),
            "seq_poseid": ("Results", "Pose_ID"),
            "seq_interactionid": ("Interaction_indices", "interaction_id"),
            "seq_interactionposeid": ("Interactions", "interaction_pose_ID"),
        }
        for seq_name, (table, column) in sequence_map.items():
            max_val = self.conn.execute(
                f"SELECT COALESCE(MAX({column}), 0) FROM {table}"
            ).fetchone()[0]
            self.conn.execute(f"DROP SEQUENCE IF EXISTS {seq_name};")
            self.conn.execute(f"CREATE SEQUENCE {seq_name} START {max_val + 1};")

    def _rollback_merge(self, merge_id: int):
        """
        Remove all data inserted by a specific merge session

        Args:
            merge_id (int): Merge session for which to delete associated data
        """

        # delete from Interactions where pose_id was added by this merge
        self.conn.execute(
            """
            DELETE FROM Interactions WHERE Pose_ID IN (
                SELECT merged_PK FROM PK_conversions 
                WHERE merge_id = ? AND table_name = 'Results')""",
            (merge_id,),
        )

        # delete from Interaction_indices
        self.conn.execute(
            """
            DELETE FROM Interaction_indices WHERE interaction_id IN (
                SELECT merged_PK FROM PK_conversions 
                WHERE merge_id = ? AND table_name = 'Interaction_indices'
                AND merged_PK != original_PK)""",
            (merge_id,),
        )

        # delete from Results
        self.conn.execute(
            """
            DELETE FROM Results WHERE Pose_ID IN (
                SELECT merged_PK FROM PK_conversions 
                WHERE merge_id = ? AND table_name = 'Results')""",
            (merge_id,),
        )

        # delete from Ligands (only newly added, not deduplicated ones)
        self.conn.execute(
            """
            DELETE FROM Ligands WHERE ligand_id IN (
                SELECT merged_PK FROM PK_conversions 
                WHERE merge_id = ? AND table_name = 'Ligands'
                AND merged_PK != original_PK)""",
            (merge_id,),
        )

        # delete from DB_properties
        self.conn.execute(
            """
            DELETE FROM DB_properties WHERE DB_write_session IN (
                SELECT merged_PK FROM PK_conversions 
                WHERE merge_id = ? AND table_name = 'DB_properties')""",
            (merge_id,),
        )

        # clean up tracking
        self.conn.execute("DELETE FROM PK_conversions WHERE merge_id = ?", (merge_id,))
        self.conn.execute("DELETE FROM merged_tables WHERE merge_id = ?", (merge_id,))

        # commit
        self.conn.commit()

        # reset sequences to current max values
        self._sync_auto_increment_state()

    # endregion

    # region Methods for dealing with bookmarks and filtering

    def _generate_result_filtering_query(
        self, filters_dict: dict, bookmark_name: str, filter_bookmark: str = None
    ) -> str:
        """
        Takes dict of filters, writes sql filtering string

        Args:
            filters_dict (dict): Keys names and value formats must match those found in the Filters class
            bookmark_name (str): Name of resulting bookmark to which passing poses are saved
            filter_bookmark (str, optional): Option to filter across pre-filtered data instead of entire database. Defaults to None

        Returns:
            str: database engine formatted query string
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
        elif isinstance(outfields, Union[list, tuple]):
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
        # TODO Inaccurate patch, need to specify how to aggregate if grouing the results
        # this is added here assuming it will be grouped, and is for testing duckdb performance
        # only. Will need to be made more rigorous if including duckdb in prod.
        duck_formatted_outfields = [f"{field}" for field in formatted_outfields]

        return duck_formatted_outfields

    def crossref_databases(
        self,
        wanted_dbs: list[tuple[str, Union[str, None]]] = None,
        unwanted_dbs: list[tuple[str, Union[str, None]]] = None,
        bookmark_prefix: str = "crossref",
    ) -> tuple[int, dict, list, dict]:
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

        Returns:
            tuple[int, dict, list, dict]: Number of "passing" ligands, dict of database {database
            file path: new bookmark name from cross ref}, list of passing ligand names, dict of dbs and how they were compared
        """
        processed_wanted = []
        for index, database in enumerate(wanted_dbs):
            if index == 0:
                db_name = "main"
            else:
                db_name = f"db_{index + 1}"
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
            db_name = f"db_{index + 1 + len(processed_wanted)}"
            processed_unwanted.append(
                (
                    db_name,
                    database[0],
                    database[1],
                )
            )
            self._attach_db(database[0], db_name)

        # make subqueries/common table expressions for each set of ligands
        selection_query = """{db_name} AS (SELECT LigName FROM {db_name}.Ligands INNER JOIN {db_name}.Results on {db_name}.Results.ligand_id={db_name}.Ligands.ligand_id WHERE Pose_id IN ({bookmark_poses}))"""
        wanted_ctes = []
        for db_name, _, bookmark in processed_wanted:
            wanted_ctes.append(
                selection_query.format(
                    db_name=db_name,
                    bookmark_poses=self._get_bookmark_poses_query(bookmark, db_name),
                )
            )

        unwanted_ctes = []
        for db_name, _, bookmark in processed_unwanted:
            unwanted_ctes.append(
                selection_query.format(
                    db_name=db_name,
                    bookmark_poses=self._get_bookmark_poses_query(bookmark, db_name),
                )
            )

        all_ctes = "WITH\n" + ", ".join(wanted_ctes + unwanted_ctes)

        wanted_joins = " ".join(
            f"INNER JOIN {db[0]} ON {processed_wanted[0][0]}.LigName = {db[0]}.LigName"
            for db in processed_wanted[1:]
        )

        # build exclusion union
        exclude_union = " UNION ".join(
            f"SELECT LigName FROM {db[0]}" for db in processed_unwanted
        )

        if unwanted_dbs:
            unwanted_ctes = f""",
                unwanted AS (
                    {exclude_union}
                )"""
            unwanted_where = """WHERE i.LigName NOT IN (SELECT LigName FROM unwanted)"""
        else:
            unwanted_ctes = ""
            unwanted_where = ""

        query = f"""
        {all_ctes},
        wanted AS (
            SELECT {processed_wanted[0][0]}.LigName
            FROM {processed_wanted[0][0]} {wanted_joins}
        ){unwanted_ctes}
        SELECT DISTINCT i.LigName
        FROM wanted i
        {unwanted_where}
        """

        # the query now works, time to make a bookmark in each database with the results
        approved_ligand_names = [lig[0] for lig in self.db_query(query).fetchall()]
        dbs_new_bookmark_names = {}

        # don't proceed if no approved ligands
        if len(approved_ligand_names) < 1:
            return 0, {}, [], {}
        else:
            for _, file, bookmark in processed_wanted + processed_unwanted:
                dbs_new_bookmark_names.update({file: f"{bookmark_prefix}_{bookmark}"})
            return (
                len(approved_ligand_names),
                dbs_new_bookmark_names,
                approved_ligand_names,
                {
                    "wanted": processed_wanted,
                    "unwanted": processed_unwanted,
                },
            )

    # endregion

    # region data

    def _check_if_column_in_table(self, table_name: str, column_name: str) -> bool:
        """
        Checks if a column exist in given table

        Args:
            table_name (str):
            column_name (str):

        Returns:
            bool: True if column there, False if not
        """
        column_tuples = self.db_query(
            f"""SELECT column_name
                    FROM duckdb_columns()
                    WHERE table_name = '{table_name}';"""
        ).fetchall()
        column_names = [column[0] for column in column_tuples]
        return bool(column_name in column_names)

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

    def _fetch_table_column_names(self, table: str) -> list:
        """Fetches list of string for column names in results table

        Returns:
            list: List of strings of results table column names

        Raises:
            StorageError
        """
        data = self.db_query(f"PRAGMA table_info({table})").fetchall()
        return [column_tuple[1].lower() for column_tuple in data]

    def fetch_summary_data(
        self,
        columns: list[str] = ["docking_score", "leff"],
        percentiles: list[int] = [1, 10],
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
        if self.db_empty():
            raise StorageError("There is no data in the database.")

        query = """
        SELECT 
            (SELECT COUNT(ligand_id) FROM Ligands) AS num_ligands,
            (SELECT COUNT(Pose_id) FROM Results) AS num_poses,
            (SELECT COUNT(interaction_id) FROM Interaction_indices) AS num_unique_interactions,
            (SELECT COUNT(*) 
            FROM (SELECT ANY_VALUE(interaction_id) 
                    FROM Interaction_indices 
                    GROUP BY interaction_type, rec_resid, rec_chain)
            ) AS num_interacting_residues;"""
        data = self.db_query(query).fetchone()
        summary_data = {
            "num_ligands": data[0],
            "num_poses": data[1],
            "num_unique_interactions": data[2],
            "num_interacting_residues": data[3],
        }

        allowed_columns = self._fetch_table_column_names("Results")
        for col in columns:
            if col not in allowed_columns:
                logger.warning(
                    f"Requested column {col} not found in Results table and will not be used for the summary. Allowed columns: {allowed_columns}"
                )
                columns.pop(col)
                continue
            data = self.db_query(
                f"SELECT MIN({col}), MAX({col}) FROM Results"
            ).fetchone()
            summary_data[f"min_{col}"] = data[0]
            summary_data[f"max_{col}"] = data[1]

            for p in percentiles:
                summary_data[f"{p}%_{col}"] = self._calc_percentile_cutoff(p, col)

        return summary_data

    def fetch_clustered_similars(self, ligname: str) -> tuple[list, str, str]:
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

        cluster_info = self.db_query(
            "SELECT cluster_id, cluster_window, name FROM Clusters;"
        ).fetchall()

        print(
            "Here are the existing clustering groups. Please ensure that you query ligand(s) is a part of the group you select."
        )
        print(
            "   Choice number   |   Underlying filter bookmark   |   Morgan or interaction fingerprint_cutoff   "
        )
        print(
            "----------------------------------------------------------------------------------------------------------"
        )
        cluster_options = []
        for cluster_id, filter_window, name in cluster_info:
            cluster_options.append(cluster_id)
            option_list = [str(cluster_id), filter_window, name]
            print(f"{'             |    '.join(option_list)}")

        cluster_choice = input(
            "Please specify choice number for the cluster you would like to return similar ligands from: "
        )
        if not int(cluster_choice) in cluster_options:
            raise ValueError(
                f"Given cluster number {cluster_choice} does not exist in the database. Please be sure you are specifying an integer in the given cluster options."
            )
        # get group(s) that poses of that ligand belongs to, make set so unique groups only
        groups = {
            row[0]
            for row in self.db_query(
                """
                SELECT cluster_group FROM pose_clusters
                    WHERE cluster_id = ? 
                    AND pose_id IN (
                        SELECT pose_id FROM Results
                        WHERE ligand_id = (
                            SELECT ligand_id FROM Ligands 
                                WHERE ligname = ?));
                """,
                (
                    cluster_id,
                    ligname,
                ),
            ).fetchall()
        }
        # group with those poses
        input_params = [cluster_id, *groups]
        placeholders = ",".join(["?"] * len(groups))
        ligands = self.db_query(
            f"""
            SELECT L.LigName FROM Results AS R
            JOIN Ligands AS L
                ON R.ligand_id = L.ligand_id
            WHERE R.pose_id IN (
                SELECT pose_id FROM Pose_clusters 
                WHERE cluster_id = ?
                AND cluster_group IN ({placeholders}))
            GROUP BY L.LigName;
            """,
            input_params,
        ).fetchall()

        cluster_name = (
            self.db_query(
                "SELECT name FROM Clusters WHERE cluster_id = ?;", [cluster_id]
            )
            .fetchone()[0]
            .replace(".", "p")
        )

        bookmark_name = f"similar_{ligname}_{cluster_name}"

        return ligands, bookmark_name, cluster_name

    def _calc_percentile_cutoff(self, percentile: float, column="docking_score"):
        """Make query for percentile by calculating energy or leff cutoff

        Args:
            percentile (float): cutoff percentile
            column (str, optional): string indicating column for percentile to be calculated over

        Returns:
            float: effective cutoff value of results based on percentile
        """
        # get total number of ligands
        MIN_columns = [
            "docking_score",
            "leff",
            "energies_inter",
            "energies_vdw",
            "energies_electro",
            "energies_flexLig",
            "energies_flexLR",
            "energies_intra",
            "energies_torsional",
            "unbound_energy",
        ]
        # whether the best value in a column is lowest or highest value
        if column in MIN_columns:
            lim_kw = "MIN"
        else:
            lim_kw = "MAX"
        try:
            logger.debug(f"Generating percentile filter query for {column}")
            # use duckdb internal method to calculate percentile, may differ slightly from sqlite manual calc
            percentile_fraction = percentile / 100
            query = f"""
            SELECT quantile_disc(min_val, {percentile_fraction}) AS cutoff
            FROM (
                SELECT {lim_kw}({column}) AS min_val
                FROM Results
                GROUP BY ligand_id
            ) percentile_cutoff
            """
            cutoff = self.db_query(query).fetchone()[0]
            logger.debug(f"{column} percentile cutoff is {cutoff}")
            return cutoff
        except duckdb.OperationalError as e:
            raise StorageError("Error while generating percentile query") from e

    # endregion

    # region general database operations

    def _begin_transaction(self):
        """
        Begin a transaction
        """
        pass
        # begin transaction can block with big chunks during results insert
        # self.conn.execute("BEGIN TRANSACTION;")

    def _rollback(self):
        """
        Roll back transaction
        """
        self.conn.execute("ROLLBACK;")

    def tables_in_db(self) -> list:
        """
        Returns a list of all table names in the database

        Returns:
            list: list of table names
        """
        return [row[0].lower() for row in self.db_query("SHOW TABLES").fetchall()]

    def _delete_nontables(self):
        """
        Deletes objects in the database that are not tables (handled separately)
        """
        sequences = self.db_query("""SELECT sequence_name FROM duckdb_sequences();""")
        if sequences:
            for sequence in sequences.fetchall():
                self.db_query(f"DROP SEQUENCE {sequence[0]};")

    def _get_length_of_table(self, table_name: str):
        """
        Finds the rowcount/length of a table based on the rowid

        Args:
            table_name (str): name of table to count the length of

        Returns:
            int: length of the table
        """
        query = f"""SELECT COUNT(*) from {table_name};"""

        return self.db_query(query).fetchone()[0]

    def clone(self, backup_name: str = None) -> str:
        import shutil

        """Creates a copy of the db

        Args:
            backup_name (str, optional): name of the cloned database
        
        Returns:
            str: path of backed up database
        """
        if backup_name is None:
            backup_name = self.db_file + ".bk"
        shutil.copy(self.db_file, backup_name)

        logger.info(f"Database {self.db_file} was backed up to {backup_name}.")
        return backup_name

    def check_ringtaildb_version(self) -> tuple[bool, str]:
        # TODO rejigger for duckdb, metadata table?
        """
        Checks the database version and confirms whether the code base is compatible with it

        Returns:
            bool: whether or not db is compatible with the code base
            str: current database version
        """
        db_schema_ver = "3.0.0"
        if version("ringtail") in self._db_schema_code_compatibility[db_schema_ver]:
            is_compatible = True
        else:
            is_compatible = False
            logger.warning(
                f"Database version {db_schema_ver} is NOT compatible with code base version {version('ringtail')}"
            )
        return is_compatible, db_schema_ver

    def _create_connection(self):
        """Creates database connection to self.db_file

        Returns:
            duckdb.DuckDBPyConnection: Connection object to self.db_file

        Raises:
            DatabaseConnectionError
        """
        try:
            conn = duckdb.connect(self.db_file)
            # TODO database settings?
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
        # TODO maybe, check attached + vacuuming

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
            vacuum (bool, optional): not used for duckdb, present for signature compatibility
            reindex (bool, optional): not used for duckdb, present for signature compatibility
        """
        if attached_db_alias is not None:
            self._detach_db(attached_db_alias)

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

    def _vacuum(self):
        """duckdb doesn't really have this, can VACUUM ANALYZE <table>"""
        pass

    def _attach_db(self, new_db: str, new_db_alias: str = "attached_db") -> str:
        """Attaches new database file to current database

        Args:
            new_db (str): file name for database to attach
            new_db_name (str): name of new database

        Raises:
            StorageError
        """
        attach = f"""ATTACH '{new_db}' AS {new_db_alias};"""
        self.db_query(attach)

    def _detach_db(self, new_db_alias):
        """Detaches new database file from current database

        Args:
            new_db_name (str): db name for database to detach

        Raises:
            StorageError
        """
        self.db_query(f"""DETACH {new_db_alias};""")

    def db_query(self, query, params: tuple = (), commit=False) -> iter:
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

    def _create_status_tables(self) -> None:
        """
        Creates pose status tables if needed
        """
        self.db_query(
            f"""
                CREATE TABLE IF NOT EXISTS Accepted 
                (Pose_ID INTEGER UNIQUE,
                FOREIGN KEY (Pose_ID) REFERENCES Results(Pose_ID)
                );""",
        )

        # create maybe table
        self.db_query(
            f"""
                CREATE TABLE IF NOT EXISTS Maybe 
                (Pose_ID INTEGER UNIQUE,
                FOREIGN KEY (Pose_ID) REFERENCES Results(Pose_ID)
                );""",
        )

        # create rejected table
        self.db_query(
            f"""
                CREATE TABLE IF NOT EXISTS Rejected 
                (Pose_ID INTEGER UNIQUE,
                FOREIGN KEY (Pose_ID) REFERENCES Results(Pose_ID)
                );""",
            commit=True,
        )

        return None

    # endregion
