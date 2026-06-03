#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail storage adaptors
#

from .logutils import get_logger

logger = get_logger(__name__)
import sys
from rdkit import Chem
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
            pose_coordinates  VARCHAR,
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

    def _insert_results_in_temp_tables(
        self, results_array: list[dict], interactions_array
    ):
        """
        Inserts docking results and interactions into their respective
        temporary tables

        Args:
            results_array (list): list of result rows
            interactions_array (list): list of interaction rows
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
            VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
            """
        self.db_update(
            temp_results_insert, [tuple(r) for r in results_array], commit=False
        )

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
            JOIN Ligands AS L ON L.ligname = T.ligname
            RETURNING pose_id, ligand_id, run_number, pose_rank;"""

        temp_to_interaction = """
            INSERT INTO Interactions(pose_id, interaction_id)
            SELECT M.pose_id, II.interaction_id
            FROM Interactions_temp IT
            JOIN pose_map M
                ON IT.ligname = L.ligname
                AND IT.pose_rank = M.pose_rank
                AND IT.run_number = M.run_number
            JOIN Ligands L
                ON M.ligand_id = L.ligand_id
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
        pose_coordinates,
        flexible_res_coordinates
        """
        delete_int_sql = """
        DELETE FROM Interactions_temp
        WHERE EXISTS (
            SELECT 1
            FROM Results_temp AS RT
            JOIN Results AS R
                ON RT.receptor    = R.receptor
                AND RT.pose_coordinates     = R.pose_coordinates
                AND RT.flexible_res_coordinates     = R.flexible_res_coordinates
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
                AND RT.pose_coordinates     = R.pose_coordinates
                AND RT.flexible_res_coordinates     = R.flexible_res_coordinates
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
        pose_coordinates,
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
                tuple(delete_pose_ids_list),
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

    def _insert_receptors(self, receptor_array):
        """Takes array of receptor rows, inserts into Receptors table

        Args:
            receptor_array (list): List of lists
                containing formatted receptor rows

        Raises:
            DatabaseInsertionError
        """
        sql_insert = """INSERT INTO Receptors (
        recname,
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
        count = self.table_length("Receptors")

        if count == 0:
            # Insert receptor statement
            query = f"""INSERT INTO Receptors (
                      recname,
                      receptor_object)
                      VALUES (?,?);"""

        else:
            query = """UPDATE Receptors SET recname = ?, receptor_object = ? WHERE receptor_id == 1"""
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
                      recname,
                      polymer)
                      VALUES (?,?);"""

        else:
            query = """UPDATE Receptors SET recname = ?, polymer = ? WHERE receptor_id == 1;"""
        self.db_query(query, (rec_name, receptor), commit=True)

    def _insert_db_properties(self, docking_mode: str, number_of_poses: str):
        """Insert db properties into database properties table

        Args:
            docking_mode (str): docking mode for the current dataset being written
            number_of_poses (str): number of poses written to database in current session, either "all" or specified max_poses
        """
        sql_insert = """INSERT INTO db_properties (
        docking_mode,
        number_of_poses
        ) VALUES (?,?);"""
        self.db_query(sql_insert, [docking_mode, number_of_poses], commit=True)

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
        FOREIGN KEY (pose_id) REFERENCES Results(pose_id)
        );"""
        self.db_query(table_sql, commit=True)

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
            pose_rows (list): pose and cluster id and group id
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
        Prepares for merging by dropping tables incompatible with merging
        """
        # create new tables to keep track of merger
        self._create_merge_tables()
        # delete incompatible tables and sequences for main db
        self._delete_table("Pose_clusters")
        self._delete_table("Cluster_groups")
        self._delete_table("Clusters")
        self._delete_table("Filtered_poses")
        self._delete_table("Filters")
        # Drop user-created indices on tables that grow during merge
        for table in ("Results", "Interactions"):
            indices = self.db_query(
                f"SELECT name FROM sqlite_master WHERE type='index' AND tbl_name='{table}' AND name NOT LIKE 'sqlite_%'"
            ).fetchall()
            for (index_name,) in indices:
                self.db_query(f"DROP INDEX IF EXISTS {index_name}")

    def _merging_receptors_compatible(self) -> bool:
        """
        Checks if the receptor names in the two databases are a mathc

        Returns:
            bool: True if receptor allows merging
        """
        receptorcheck_sql = """
        SELECT CASE 
            WHEN Receptors.recname = merging_receptors.recname THEN 'True'
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
                "CREATE INDEX IF NOT EXISTS ak_merge ON PK_conversions(merge_id, original_PK)"
            )
        except Exception as e:
            raise StorageError(e) from e

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
                'db_properties', 
                DB_write_session,
                DB_write_session + (SELECT MAX(DB_write_session) FROM db_properties) 
                FROM merging.db_properties;"""
            cur.execute(
                convert_dbprop_sql,
                (merge_id,),
            )

            insert_dbprops_sql = """INSERT INTO db_properties (
                DB_write_session,
                docking_mode,
                number_of_poses)
                SELECT 
                    (SELECT merged_PK FROM PK_conversions WHERE original_PK = DB_write_session and merge_id = ? and table_name = 'db_properties'),
                    docking_mode,
                    number_of_poses
                FROM merging.db_properties;"""
            cur.execute(insert_dbprops_sql, (merge_id,))
        except Exception as e:
            raise StorageError(
                "Error encountered while merging db_properties table"
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
            (SELECT merged_PK FROM PK_conversions WHERE original_PK = ligand_id and merge_id = ? AND table_name = 'Ligands') new_id,
            ligname,
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
        ON pose.original_PK = mr.pose_id
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
            raise StorageError(
                f"Error encountered while merging Ligands and Results tables: \n{str(e)}"
            ) from e

    def _merge_interaction_tables(self, merge_id: int):
        """
        Merges the interaction tables. Interaction definitions are unique and independent of the Results table, so we only
        insert those that are new with updated PK, and assign existing interaction_ids to those already described in primary db

        Args:
            merge_id (int): merge session id

        Raises:
            Exception
        """
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
        insert_interactions = """
        INSERT INTO Interactions (
        pose_id,
        interaction_id
        )    SELECT P.merged_pk as pose_id, II.merged_pk as interaction_id
                FROM merging.Interactions I
                LEFT JOIN (SELECT original_PK, merged_pk
                FROM PK_conversions
                WHERE table_name = 'Results' 
                AND merge_id = ?) P ON (I.pose_id = P.original_PK)
            LEFT JOIN (SELECT original_PK, merged_pk
                FROM PK_conversions
                WHERE table_name = 'Interaction_indices' 
                AND merge_id = ?)  II ON (I.interaction_id = II.original_PK);"""

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

    def _sync_auto_increment_state(self):
        """Reset all sequences to current MAX values in their tables"""
        pass

    def _rollback_merge(self, merge_id: int):
        """
        Remove all data inserted by a specific merge session

        Args:
            merge_id (int): Merge session for which to delete associated data
        """

        # delete from Interactions where pose_id was added by this merge
        self.conn.execute(
            """
            DELETE FROM Interactions WHERE pose_id IN (
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
            DELETE FROM Results WHERE pose_id IN (
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

        # delete from db_properties
        self.conn.execute(
            """
            DELETE FROM db_properties WHERE DB_write_session IN (
                SELECT merged_PK FROM PK_conversions 
                WHERE merge_id = ? AND table_name = 'db_properties')""",
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
        self, filters_dict, bookmark_name: str, input_bookmark: str
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
        if input_bookmark is not None:
            if input_bookmark == bookmark_name:
                # cannot write data from bookmark_a to bookmark_a
                logger.error(
                    f"Specified 'input_bookmark' and 'bookmark_name' are the same: {bookmark_name}"
                )
                raise OptionError(
                    "'input_bookmark' and 'bookmark_name' cannot be the same! Please rename 'bookmark_name'"
                )
            # cannot use percentile for an already reduced dataset
            if (
                filters_dict["score_percentile"] is not None
                or filters_dict["le_percentile"] is not None
            ):
                raise OptionError(
                    "Cannot use 'score_percentile' or 'le_percentile' with 'input_bookmark'."
                )
            # filtering window can be specified bookmark, as opposed to entire database using Results table
            if self.is_bookmark(input_bookmark):
                filtering_window = f"""(SELECT * FROM Results WHERE pose_id IN ({self.QueryBuilder.bookmark_query(input_bookmark)}))"""
            elif self._is_statustable(input_bookmark):
                filtering_window = f"""(SELECT * FROM Results WHERE pose_id IN (SELECT pose_id FROM {input_bookmark}))"""

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
                    [f"ligname LIKE '%{ligname}%' " for ligname in lig_names if ligname]
                )
                ligname_query = "SELECT ligand_id FROM Ligands WHERE " + ligname_query
            elif "ligand_name_file" in lig_filters:
                csv_path = lig_filters.pop("ligand_name_file")
                self._create_ligname_temp_table(csv_path)
                ligname_query = "SELECT ligand_id FROM Ligands JOIN tmp_lignames ON ligname = tmp_lignames.ligandname"
            # rdkit queries need to be handled in memory separate from the main query
            if lig_filters:
                rdkit_query = True

        ### Join each of the filter groups
        # if filter queries exist for each group, string them together appropriately
        if int_query:
            # add with a join statement
            partial_filter_query += f"JOIN ({int_query}) I ON R.pose_id = I.pose_id "
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
            partial_filter_query = " R.pose_id IN ({0})".format(
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
        query = "SELECT pose_id FROM (SELECT pose_id "
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
        query += f") GROUP BY pose_id HAVING COUNT(DISTINCT filtered_interactions) >= ({num_of_interactions}) "

        return query

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

    def _format_output_fields(
        self, outfields: Union[str, list], results_alias="R", ligands_alias="L"
    ) -> str:
        """Handles string or list input of column names to be outputted, will make sure ligname
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
            outfields_list.insert(0, "ligname")
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
        ligand_sql_string = ", ".join(f"'{l}'" for l in ligand_list)
        # construct the query
        if store_best_pose:
            return f"""
            SELECT R.pose_id FROM Results AS R
            JOIN (
                SELECT ligand_id, MIN(pose_rank) as best_pose
                FROM Results
                GROUP BY ligand_id
                ) AS sel
            ON R.ligand_id = sel.ligand_id
            AND R.pose_rank = sel.best_pose
            JOIN Ligands AS L
            ON R.ligand_id = L.ligand_id
            WHERE L.ligname
            IN ({ligand_sql_string})
            """
        else:
            return f"""
            SELECT pose_id FROM Results
            WHERE ligand_id IN (
                SELECT ligand_id FROM Ligands 
                WHERE ligname IN ({ligand_sql_string})
                )
                """

    # endregion

    # region data

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
            (SELECT COUNT(pose_id) FROM Results) AS num_poses,
            (SELECT COUNT(interaction_id) FROM Interaction_indices) AS num_unique_interactions,
            (SELECT COUNT(*) 
            FROM (SELECT interaction_id
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

    def fetch_cluster_options(self, ligname: str) -> list[tuple]:
        """Return available clustering groups as (cluster_id, cluster_window, name) tuples.

        Args:
            ligname (str): ligname for ligand to find similarity with

        Returns:
            list[tuple]: list of (cluster_id, cluster_window, name)
        """
        return self.db_query(
            "SELECT cluster_id, cluster_window, name FROM Clusters;"
        ).fetchall()

    def fetch_clustered_similars(
        self, ligname: str, cluster_id: int
    ) -> tuple[list, str, str]:
        """Given ligname and a chosen cluster_id, return similar ligands.

        Args:
            ligname (str): ligname for ligand to find similarity with
            cluster_id (int): cluster to search within

        Returns:
            tuple[list, str, str]: (ligands, bookmark_name, cluster_name)

        Raises:
            DatabaseQueryError
        """
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
            SELECT L.ligname FROM Results AS R
            JOIN Ligands AS L
                ON R.ligand_id = L.ligand_id
            WHERE R.pose_id IN (
                SELECT pose_id FROM Pose_clusters
                WHERE cluster_id = ?
                AND cluster_group IN ({placeholders}))
            GROUP BY L.ligname;
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

    # endregion

    # region general database operations

    def _begin_transaction(self):
        """
        Begin a transaction
        """
        self.conn.execute("BEGIN TRANSACTION;")

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
        tables = [
            name[0].lower()
            for name in self.db_query(
                "SELECT name FROM sqlite_master WHERE type='table';"
            ).fetchall()
        ]
        if "sqlite_sequence" in tables:
            tables.remove("sqlite_sequence")
        return tables

    def _create_indices(self):
        """Create alternate-key indices ('ak_*') for queryable tables."""
        logger.debug("Creating columns indices...")
        for table, schema in [
            (RESULTS_SCHEMA.name, RESULTS_SCHEMA),
            (INTERACTIONS_SCHEMA.name, INTERACTIONS_SCHEMA),
            (INTERACTION_INDICES_SCHEMA.name, INTERACTION_INDICES_SCHEMA),
        ]:
            for sql in build_create_indices(table, schema):
                self.db_query(sql)
        self.conn.commit()
        logger.info(
            "Indices created for Results, Interactions, and Interaction_indices."
        )

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

    def _set_ringtail_db_schema_version(self, db_version: str = "3.0.0"):
        """
        Will check current storage manager db schema version and only set if it
        is compatible with the code base version (i.e., version(ringtail)).

        Args:
            db_version (str, optional): _description_. Defaults to "3.0.0".

        Raises:
            StorageError: _description_
        """

        # check that code base is compatible with db schema version
        code_version = version("ringtail")
        if code_version in self._db_schema_code_compatibility[db_version]:
            rtdb_version = db_version.replace(".", "")
            # if so, proceed to set db schema version
            self.conn.execute(f"PRAGMA user_version = {rtdb_version}")
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

    def _db_compatible_for_merge(self, merging_db_alias: str) -> bool:
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

    def update_database_version(self, new_version, backup=False):
        """Updates sqlite database schema 1.0.0 through 3.0.0.
        Upgrades step-by-step through each major version, e.g. 1.0.0 -> 1.1.0 -> 2.0.0 -> 3.0.0.

        Args:
            new_version (str): target version string
            backup (bool, optional): clone database before upgrading. Defaults to False.
        """
        self.conn = self._create_connection()
        if backup:
            self.clone()

        original_version = self.check_ringtaildb_version()[1]
        logger.info(
            f"Upgrading {self.db_file} of version {original_version} to version {new_version}:"
        )

        # upgrade to 1.1.0
        if original_version in ["1.0.0", "1.1.0"]:
            logger.warning(
                "If you created the database with the duplicate handling option, there is a chance of inconsistent behavior of anything involving interactions as the pose_id was not used as an explicit foreign key in db v1.0.0 and v1.1.0."
            )
            if original_version == "1.0.0":
                self._update_db_100_to_110()
                logger.info("\n\nSuccessfully upgraded to 1.1.0!\n\n")

            # upgrade to 2.0.0
            if new_version in ["2.0.0", "3.0.0"]:
                self._update_db_110_to_200()
                logger.info("\n\nSuccessfully upgraded to 2.0.0!\n\n")

        # upgrade to 3.0.0
        if new_version == "3.0.0" and original_version == "2.0.0":
            self._update_db_200_to_300()
            logger.info("\n\nSuccessfully upgraded to 3.0.0!\n\n")

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
            "CREATE INDEX IF NOT EXISTS ak_results ON Results(ligand_id, docking_score, leff, delta, reference_rmsd, energies_inter, energies_vdw, energies_electro, energies_intra, num_interactions, run_number, pose_rank, num_hb)"
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
            # index certain tables
            self._create_indices()
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
        #TODO
            - receptor column?
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
            h_parents: list[Union[str, int]],
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
                    return mol.ToBinary(Chem.PropertyPickleOptions.AllProps)

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
                return mol_final.ToBinary(Chem.PropertyPickleOptions.AllProps)
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
                        energies_flexLig,
                        energies_flexLR,
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
                        energies_flexLig,
                        energies_flexLR,
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

        # build new indices if you got this far successfully
        self._create_indices()
        self.db_query("REINDEX")
        self._set_ringtail_db_schema_version("3.0.0")
        self.conn.commit()

    def _drop_views(self, db_alias: str = None):
        """
        Will drop views and clear bookmark table

        Args:
            db_alias (str, optional): if needing to drop views from a connected, aliased database. Defaults to None.
        """
        if "Bookmarks" not in self.tables_in_db():
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
        elif self._is_candidates_table(table):
            query = f"SELECT COUNT(*) FROM {CANDIDATES_SUBQ} AS _c;"
            params = ()
        else:
            logger.error(f"Table -{table}- does not exist in the database.")
            return None

        return self.db_query(query, params).fetchone()[0]

    def _vacuum(self):
        """SQLite vacuum rebuilds the database file, repacking it into a minimal amount of disk space

        Raises:
            DatabaseInsertionError
        """
        try:
            cur = self.conn.cursor()
            cur.execute("VACUUM;")
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

    def _check_attached(self) -> iter:
        """
        Check what databases are attached

        Returns:
            iter: attached database names
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
        elif self._is_candidates_table(table):
            query.FROM("Results").WHERE(f"pose_id IN {CANDIDATES_SUBQ}")
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
        elif self._is_candidates_table(table):
            query.FROM("Results").WHERE(f"pose_id IN {CANDIDATES_SUBQ}")
        else:
            logger.error(f"Table -{table}- does not exist in the database.")
            return None
        return self.db_query(*query.build()).fetchone()[0]

    def set_pose_comment(self, pose_id: int, comment: str) -> None:
        if comment:
            self.db_query(
                "INSERT INTO Pose_comments (pose_id, comment) VALUES (?, ?)"
                " ON CONFLICT(pose_id) DO UPDATE SET comment = excluded.comment",
                (pose_id, comment),
            )
        else:
            self.db_query("DELETE FROM Pose_comments WHERE pose_id = ?", (pose_id,))
        self.conn.commit()

    def get_pose_comment(self, pose_id: int) -> "str | None":
        row = self.db_query(
            "SELECT comment FROM Pose_comments WHERE pose_id = ?", (pose_id,)
        ).fetchone()
        return row[0] if row else None

    # endregion
