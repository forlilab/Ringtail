#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail storage adaptor for DuckDB
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
from .ringtailoptions import Filters
from .util import numlist2str, statuses
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
from .storagemanager import StorageManager
from .querybuilder import QueryBuilderDuck
from collections import defaultdict

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

    # region Methods for creating and inserting into tables the database

    def _create_ligands_table(self, name="Ligands"):
        """
        Creates ligands table

        Args:
            name (str, optional): _description_. Defaults to "Ligands".
        """
        ligand_table = f"""
            CREATE SEQUENCE seq_ligandid START 1;
            CREATE TABLE IF NOT EXISTS {name} (
            ligand_id           INTEGER DEFAULT nextval('seq_ligandid') PRIMARY KEY,
            LigName             VARCHAR NOT NULL UNIQUE,
            ligand_smile        VARCHAR,
            rdmol               BLOB,
            atom_index_map      VARCHAR,
            hydrogen_parents    VARCHAR,
            input_model         VARCHAR)"""

        self.db_query(ligand_table, commit=True)

    def _insert_ligands(self, ligand_array: list = []) -> list:
        """
        duck db implementation of parent method

        """
        df = pd.DataFrame(
            ligand_array,
            columns=[
                "LigName",
                "ligand_smile",
                "rdmol",
                "atom_index_map",
                "hydrogen_parents",
                "input_model",
            ],
        )
        self.conn.register("df_view", df)
        # to insert interaction if unique
        sql_insert = """
        INSERT INTO Ligands (
            LigName,
            ligand_smile,
            rdmol,
            atom_index_map,
            hydrogen_parents,
            input_model
            ) 
        SELECT * FROM df_view
        ON CONFLICT(LigName) DO NOTHING;"""

        self.conn.execute(sql_insert)

    def _create_results_table(self, name="Results"):
        """Creates table for results."""

        sql_results_table = f"""
            CREATE SEQUENCE seq_poseid START 1;
            CREATE TABLE IF NOT EXISTS {name} (
            Pose_ID             INTEGER DEFAULT nextval('seq_poseid') PRIMARY KEY,
            ligand_id           INTEGER REFERENCES Ligands(ligand_id),
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
            flexible_res_coordinates   VARCHAR,
            ); """

        self.db_query(sql_results_table, commit=True)

    def _create_temporary_results_tables(self):
        try:
            create_temp_results = """
            DROP SEQUENCE IF EXISTS set_temp_poseid;
            CREATE SEQUENCE set_temp_poseid START 1;
            CREATE TEMP TABLE 
            Results_temp(
                temp_poseid         INTEGER DEFAULT nextval('set_temp_poseid') PRIMARY KEY,
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
                flexible_res_coordinates   VARCHAR,
                ligname             VARCHAR);
                """
            create_temp_interactions = """
            CREATE TEMP TABLE Interactions_temp(
                temp_poseid        INTEGER References Results_temp(temp_poseid),
                interaction_type    VARCHAR,
                rec_chain           VARCHAR,
                rec_resname         VARCHAR,
                rec_resid           VARCHAR,
                rec_atom            VARCHAR,
                rec_atomid          VARCHAR);
            """
            # create temporary tables
            self.conn.execute(create_temp_results)
            self.conn.execute(create_temp_interactions)
        except duckdb.OperationalError as e:
            raise DatabaseInsertionError(
                "Error while creating temporary results tables."
            ) from e

    def _insert_results_in_temp_tables(
        self, results_array: list, interactions_array: list
    ):
        # insert results rows
        res_df = pd.DataFrame(
            results_array,
            columns=[
                "receptor",
                "pose_rank",
                "run_number",
                "cluster_rmsd",
                "reference_rmsd",
                "docking_score",
                "leff",
                "deltas",
                "energies_inter",
                "energies_vdw",
                "energies_electro",
                "energies_flexLig",
                "energies_flexLR",
                "energies_intra",
                "energies_torsional",
                "unbound_energy",
                "nr_interactions",
                "num_hb",
                "cluster_size",
                "about_x",
                "about_y",
                "about_z",
                "trans_x",
                "trans_y",
                "trans_z",
                "axisangle_x",
                "axisangle_y",
                "axisangle_z",
                "axisangle_w",
                "dihedrals",
                "ligand_coordinates",
                "flexible_res_coordinates",
                "ligname",
            ],
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
                    flexible_res_coordinates,
                    ligname) 
                SELECT * FROM incoming_poses;"""

        self.conn.execute(temp_insert)
        int_df = pd.DataFrame(
            interactions_array,
            columns=[
                "ligand_name",
                "pose_rank",
                "run_number",
                "interaction_type",
                "rec_chain",
                "rec_resname",
                "rec_resid",
                "rec_atom",
                "rec_atomid",
            ],
        )
        self.conn.register("incoming_interaction", int_df)
        temp_interaction_insert = """
            INSERT INTO Interactions_temp (temp_poseid, interaction_type, rec_chain, rec_resname, rec_resid, rec_atom, rec_atomid)
            SELECT RT.temp_poseid, df.interaction_type, df.rec_chain, df.rec_resname, df.rec_resid, df.rec_atom, df.rec_atomid
            FROM incoming_interaction AS df
            JOIN Results_temp RT
                ON RT.ligname = df.ligand_name
                AND RT.pose_rank = df.pose_rank
                AND RT.run_number = df.run_number;"""
        self.conn.execute(temp_interaction_insert)

    def _move_tempresults_to_database(self, commit: bool = True):
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
                T.deltas,
                T.energies_inter,
                T.energies_vdw,
                T.energies_electro,
                T.energies_flexLig,
                T.energies_flexLR,
                T.energies_intra,
                T.energies_torsional,
                T.unbound_energy,
                T.nr_interactions,
                T.num_hb,
                T.cluster_size,
                T.about_x,
                T.about_y,
                T.about_z,
                T.trans_x,
                T.trans_y,
                T.trans_z,
                T.axisangle_x,
                T.axisangle_y,
                T.axisangle_z,
                T.axisangle_w,
                T.dihedrals,
                T.ligand_coordinates,
                T.flexible_res_coordinates
            FROM Results_temp AS T
            JOIN Ligands AS L ON L.LigName = T.LigName;"""

        temp_to_interaction = """
            INSERT INTO Interactions(pose_id, interaction_id)
            SELECT R.pose_id, II.interaction_id
            FROM Results AS R
            JOIN Ligands AS L
                ON R.ligand_id    = L.ligand_id
            JOIN Results_temp AS RT 
                ON RT.receptor=R.receptor
                AND RT.about_x=R.about_x
                AND RT.about_y=R.about_y
                AND RT.about_z=R.about_z
                AND RT.trans_x=R.trans_x
                AND RT.trans_y=R.trans_y
                AND RT.trans_z=R.trans_z
                AND RT.axisangle_x=R.axisangle_x
                AND RT.axisangle_y=R.axisangle_y
                AND RT.axisangle_z=R.axisangle_z
                AND RT.axisangle_w=R.axisangle_w
                AND RT.dihedrals=R.dihedrals
            JOIN Interactions_temp AS IT
                ON RT.temp_poseid      = IT.temp_poseid
            JOIN Interaction_indices AS II
                ON II.interaction_type = IT.interaction_type
                AND II.rec_chain        = IT.rec_chain
                AND II.rec_resname      = IT.rec_resname
                AND II.rec_resid        = IT.rec_resid
                AND II.rec_atom         = IT.rec_atom
                AND II.rec_atomid       = IT.rec_atomid;
        """
        try:
            self.conn.execute(temp_to_results)
            self.conn.execute(temp_to_interaction)
            if commit:
                self.conn.commit()
        except duckdb.OperationalError as e:
            raise DatabaseInsertionError(
                "Error while moving results from temp table to the database."
            ) from e

    def _create_receptors_table(self):
        """Create table for receptors. Columns are:
        Receptor_ID         INTEGER PRIMARY KEY AUTOINCREMENT,
        RecName             VARCHAR,
        box_dim             VARCHAR,
        box_center          VARCHAR,
        grid_spacing        FLOAT,
        flexible_residues   VARCHAR,
        flexres_atomnames   VARCHAR,
        receptor_object     BLOB

        """
        receptors_table = """
            CREATE SEQUENCE seq_receptorid START 1;
            CREATE TABLE IF NOT EXISTS Receptors (
            Receptor_ID         INTEGER DEFAULT nextval('seq_receptorid') PRIMARY KEY,
            RecName             VARCHAR,
            box_dim             VARCHAR,
            box_center          VARCHAR,
            grid_spacing        FLOAT,
            flexible_residues   VARCHAR,
            flexres_atomnames   VARCHAR,
            receptor_object     BLOB
        )"""
        self.db_query(receptors_table, commit=True)

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

        except duckdb.OperationalError as e:
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
        except duckdb.OperationalError as e:
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

        sql_str = """
        CREATE SEQUENCE seq_dbwriteid START 1;
        CREATE TABLE IF NOT EXISTS DB_properties (
        DB_write_session    INTEGER DEFAULT nextval('seq_dbwriteid') PRIMARY KEY,
        docking_mode        VARCHAR,
        number_of_poses     VARCHAR)"""

        self.db_query(sql_str, commit=True)

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

        except duckdb.OperationalError as e:
            raise DatabaseInsertionError(
                "Error while inserting database properties info into DB_properties table"
            ) from e

    def _create_interaction_index_table(self):
        """create table of data for each unique interaction, will be remade everytime db is written to.
        Columns are:
        interaction_id      INTEGER PRIMARY KEY AUTOINCREMENT,
        interaction_type    VARCHAR,
        rec_chain           VARCHAR,
        rec_resname         VARCHAR,
        rec_resid           VARCHAR,
        rec_atom            VARCHAR,
        rec_atomid          VARCHAR

        Raises:
            DatabaseTableCreationError: Description

        """
        interaction_index_table = """
        CREATE SEQUENCE seq_interactionid START 1;
        CREATE TABLE Interaction_indices (
        interaction_id      INTEGER DEFAULT nextval('seq_interactionid') PRIMARY KEY,
        interaction_type    VARCHAR,
        rec_chain           VARCHAR,
        rec_resname         VARCHAR,
        rec_resid           VARCHAR,
        rec_atom            VARCHAR,
        rec_atomid          VARCHAR,
        UNIQUE (interaction_type, rec_chain, rec_resname, rec_resid, rec_atom, rec_atomid));
        """

        self.db_query("""DROP TABLE IF EXISTS Interaction_indices""")
        self.db_query(interaction_index_table, commit=True)

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

        interaction_table = """
        CREATE SEQUENCE seq_interactionposeid START 1;
        CREATE TABLE IF NOT EXISTS Interactions (
        interaction_pose_ID INTEGER DEFAULT nextval('seq_interactionposeid') PRIMARY KEY,
        Pose_ID   INTEGER REFERENCES Results(pose_id),
        interaction_id INTEGER REFERENCES Interaction_indices(interaction_id))"""

        self.db_query(interaction_table, commit=True)

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

        except duckdb.OperationalError as e:
            raise DatabaseInsertionError(
                f"Error while inserting an interaction row: {e}"
            ) from e

    def _insert_interaction_index_rows(self, interactions: list[list[str]]) -> int:
        """
        Writes unique interactions and returns the interaction_id of the given interaction

        Args:
            interaction_tuple (tuple): (rec_chain, rec_resname, rec_resid, rec_atom, rec_atomid)

        Returns:
            int: interaction index

        Raises:
            DatabaseInsertionError
        """
        df = pd.DataFrame(
            interactions,
            columns=[
                "interaction_type",
                "rec_chain",
                "rec_resname",
                "rec_resid",
                "rec_atom",
                "rec_atomid",
            ],
        )
        self.conn.register("df_view", df)
        # to insert interaction if unique
        alt_sql_insert = """
                    INSERT INTO Interaction_indices(interaction_type,rec_chain,rec_resname,rec_resid,rec_atom,rec_atomid) 
                    SELECT * FROM df_view
                    ON CONFLICT DO NOTHING;
                    """
        try:
            self.conn.execute(alt_sql_insert)
        except duckdb.OperationalError as e:
            raise DatabaseInsertionError(
                f"Error inserting unique interactions in index table: {e}"
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

        except duckdb.OperationalError as e:
            raise StorageError(
                "Error while deleting rows in the Interaction table"
            ) from e

    def _create_filtering_tables(self):
        """
        Creates a Filter table which includes filter_id (PK), name (bookmark_name), duckdb formatted query,
        and dictionary of filters used, as well as Filtered_poses, which uses filter_id as FK,
        and lists all poses passing that filter_id

        Raises:
            DatabaseTableCreationError
        """
        filters_sql = """
        CREATE SEQUENCE seq_filterid START 1;
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
        self.db_query(filter_pose_sql, commit=True)

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
        # TODO rethink this later
        cur = self.conn.cursor()
        cur.execute(
            "CREATE TABLE IF NOT EXISTS Ligand_clusters (pose_id  INTEGER UNIQUE)"
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
        except duckdb.OperationalError as e:
            raise StorageError("Error occurred while inserting cluster data") from e

    def _create_indices(self):
        """Create index for specified tables and columns. 'ak' stands for 'alternate key' and is prepended to index name to avoid naming conflicts

        Raises:
            StorageError
        """
        self.db_query(
            "CREATE INDEX IF NOT EXISTS ak_results ON Results(docking_score, leff)"
        )
        self.db_query(
            "CREATE INDEX IF NOT EXISTS ak_resultids ON Results(Pose_id, ligand_id)"
        )
        self.db_query(
            "CREATE INDEX IF NOT EXISTS ak_interactions ON Interactions(Pose_id, interaction_id)"
        )
        self.db_query(
            "CREATE INDEX IF NOT EXISTS ak_ligands ON Ligands(ligand_id)", commit=True
        )

        logger.info(
            "Indicies were created for specified Results, Ligands, and Interaction_indices columns."
        )

    # endregion

    # region merge databases
    # TODO this is a biggie, start with the basic functionality so I can test speed
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
            dbfile                  VARCHAR,
            merge_start             DATETIME DEFAULT CURRENT_TIMESTAMP);"""
            cur.execute(mergetbl_sql)

            # create PK table: merge_id(FK), table, original_val, merge_val
            pktable_sql = """CREATE TABLE IF NOT EXISTS PK_conversions (
            merge_id        INTEGER,
            table_name      VARCHAR,
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
        # TODO save for later
        """takes lists of filters, writes sql filtering string

        Args:
            filters_dict (dict): dict of filters. Keys names and value formats must match those found in the Filters class

        Returns:
            str: duckdb-formatted string for filtering query
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
        # TODO this method is huge, great potential for query builder
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
        # NOTE Inaccurate patch, need to specify how to aggregate if grouing the results
        # this is added here assuming it will be grouped, and is for testing duckdb performance
        # only. Will need to be made more rigorous if including duckdb in prod.
        duck_formatted_outfields = [
            f"ANY_VALUE({field})" for field in formatted_outfields
        ]

        return duck_formatted_outfields

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
        query.SELECT("ANY_VALUE(r.pose_id)").FROM("results", "r").IN_BOOKMARK(
            bookmark_name
        )
        if grouped_by_ligand:
            query.GROUP_BY("r.ligand_id")
        query_string = query.build(count=True)[0]
        # TODO problematic query: ", query_string
        if self.db_query(query_string).fetchone():
            return self.db_query(query_string).fetchone()[0]
        else:
            return 0

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

    # endregion

    # region crossreferencing filtered databases

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
        except duckdb.OperationalError as e:
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
        except duckdb.OperationalError as e:
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

        except duckdb.OperationalError as e:
            raise StorageError(
                f"Error occured while pruning Interactions not in {bookmark_name}"
            ) from e

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
        except duckdb.OperationalError as e:
            raise DatabaseTableCreationError(
                f"Error while creating temporary table {table_name}"
            ) from e

    # endregion

    # region fetching specific columns

    def _get_possible_output_columns(self, tables=["Results", "Ligands"]):
        # TODO database specific proably
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
            """SELECT name FROM duckdb_master WHERE type='table' AND name='Interactions';"""
        )
        if len(cur.fetchall()) == 0:
            return None

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

    def tables_in_db(self) -> list:
        """
        Returns a list of all table names in the database

        Returns:
            list: list of table names
        """
        return [row[0].lower() for row in self.db_query("SHOW TABLES").fetchall()]

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
        except duckdb.OperationalError as e:
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

        except duckdb.OperationalError as e:
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
        except duckdb.OperationalError as e:
            raise StorageError("Error while generating percentile query") from e

    # endregion

    # region general database operations

    def _get_length_of_table(self, table_name: str):
        """
        Finds the rowcount/length of a table based on the rowid

        Args:
            table_name (str): name of table to count the length of

        Returns:
            int: length of the table
        """
        query = f"""SELECT COUNT(*) from {table_name}"""

        return self.db_query(query).fetchone()[0]

    def clone(self, backup_name=None):
        """Creates a copy of the db

        Args:
            backup_name (str, optional): name of the cloned database
        """
        # TODO
        if backup_name is None:
            backup_name = self.db_file + ".bk"
        bck = duckdb.connect(backup_name)
        with bck:
            self.conn.backup(bck, pages=1)
        bck.close()
        logger.info(f"Database {self.db_file} was backed up to {backup_name}.")

    def _set_ringtail_db_schema_version(self, db_version: str = "3.0.0"):
        # TODO rejigger for duckdb, metadata table?
        # TODO this might be the only onw I need working for testing
        """Will check current storage manager db schema version and only set if it is compatible with the code base version (i.e., version(ringtail)).

        Raises:
            StorageError: if versions are incompatible
        """
        return
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
        # TODO rejigger for duckdb, metadata table?
        """
        Checks the database version and confirms whether the code base is compatible with it

        Returns:
            bool: whether or not db is compatible with the code base
            str: current database version
        """
        # cur = self.conn.cursor()
        # db_version = str(cur.execute("PRAGMA user_version").fetchone()[0])
        # if db_version == "0":
        #     # ringtail 1.0.0 did not have a user version, so catch if database has contents and version 0
        #     cur.execute(
        #         "SELECT EXISTS(SELECT 1 FROM duckdb_master WHERE type='table' AND name='Results');"
        #     )
        #     # if db version is 0 but has a results table, it is 1.0.0
        #     if cur.fetchone()[0] != 0:
        #         db_version = "100"
        #     # else empty or corrupt database
        #     else:
        #         raise StorageError(
        #             f"The database requested {self.db_file} does not exist or does not have any tables. Check for spelling errors, else the database may be corrupt (delete the file before using the same name again)"
        #         )
        # db_schema_ver = ".".join([*db_version])
        db_schema_ver = "3.0.0"
        if version("ringtail") in self._db_schema_code_compatibility[db_schema_ver]:
            is_compatible = True
        else:
            is_compatible = False
            logger.warning(
                f"Database version {db_schema_ver} is NOT compatible with code base version {version('ringtail')}"
            )
        # cur.close()
        return is_compatible, db_schema_ver

    def _check_if_db_compatible_for_merge(self, merging_db_alias: str) -> bool:
        # TODO rejigger for duckdb, metadata table?
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

    def _create_connection(self) -> duckdb.DuckDBPyConnection:
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
        """#TODO duckdb doesn't really have this, can VACUUM ANALYZE <table>"""
        pass

    def _attach_db(self, new_db: str, new_db_alias: str = "attached_db") -> str:
        """Attaches new database file to current database

        Args:
            new_db (str): file name for database to attach
            new_db_name (str): name of new database

        Raises:
            StorageError
        """
        pass

    def _detach_db(self, new_db_alias):
        """Detaches new database file from current database

        Args:
            new_db_name (str): db name for database to detach

        Raises:
            StorageError
        """
        pass

    def db_query(self, query, params: tuple = (), commit=False) -> iter:
        """Executes provided duckdb query. Returns cursor for results.
            Since cursor remains open, added to list of open cursors

        Args:
            query (str): Formated duckdb query as string

        Returns:
            duckdb cursor: Contains results of query
        """

        try:
            cur = self.conn.cursor()
            cur.execute(query, params)
            if commit:
                self.conn.commit()
        except duckdb.OperationalError as e:
            raise DatabaseQueryError(
                f"Unable to execute query -{query}- with given parameters -{params}-: -{e}-"
            ) from e
        return cur

    def db_update(self, query: str, parameters: list[list], commit=True) -> iter:
        """
        A db query that also commits if/when specified

        Args:
            query (str): duckdb formatted query string
            parameters (list[list]): assumes appropriate place holders in query
            commit (bool, optional): whether to commit the transaction in open connection. Defaults to True.

        Raises:
            OptionError
            DatabaseInsertionError

        Returns:
            iter: if requesting return value(s)
        """
        if not parameters:
            return
        if type(parameters) != list:
            raise OptionError(
                "Input for duckdb execute many needs to be list[list], wrong type used:",
                parameters,
            )
        if type(parameters[0]) != list:
            if type(parameters[0]) == tuple:
                # just convert tuple to list
                parameters = [list(row) for row in parameters]
            else:
                raise OptionError(
                    "Input for duckdb execute many needs to be list[list], wrong type used:",
                    parameters,
                )
        try:
            cur = self.conn.cursor()
            cur.executemany(query, parameters)
            if commit:
                self.conn.commit()
        except duckdb.OperationalError as e:
            raise DatabaseInsertionError(f"Error while committing insert query") from e

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

    def create_status_tables(self) -> None:
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
