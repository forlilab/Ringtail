from dbmanager import DBManager

import fileparsers


class DBManager():
    """ DOCUMENTATION GOES HERE """
    def __init__(self, fname, opts={}):
        self.fname = fname
        self.opts = {}
        self.write_flag = self.opts["write_db_flag"]
        self.conn = self.__create_connection()
        self.bv_table_flag = False

        if self.write_flag:
            #create tables in db
            self.__create_results_table()
            self.__create_ligands_table()
            self.__create_interaction_index_table()


    ##########################
    ##### Public methods #####
    ##########################

    def add(self, dict_of_results):
        """ DOCUMENTATION GOES HERE """
        e_data = dict_of_results['energy']
        self._add_energy(e_data)
        pass

    def get_results(self, dict_of_filters):
        """

        dict_of_filters = {
                'energy' : {'min': None, 'max':None},
                'leff' : {'min': None, 'max':None},
                'interactions' : {},
        }


        """
        # filter by energiy, if requested
        # filter by LE if requested
        pass

    def insert_interactions(self, unique_interactions = None, unique_interactions_split = None, new_interactions = None, interaction_idx_counter = None):
        #check if we need to initialize the interaction bv table and insert first set of interactions
        if new_interactions == None:
            self.__create_interaction_bv_table(len(unique_interactions))
            self.__insert_unique_interaction(np.array(unique_interactions_split))
            self.bv_table_flag = True

        #otherwise, insert individual new interactions
        else:
            for interaction in new_interactions:
                self.__insert_one_interaction(interaction)
                self.__make_new_interaction_column(interaction_idx_counter)
                interaction_idx_counter += 1

    def insert_interaction_BVs(self, bitvectors, num_discrete_interactions):
        #print("Inserting interaction bitvectors...")
        interaction_columns = range(num_discrete_interactions)
        column_str = ""
        filler_str = "?,"
        for i in interaction_columns:
            column_str += "Interaction_"+ str(i+1) + ", "
            filler_str += "?,"
        column_str = column_str.rstrip(", ")
        filler_str = filler_str.rstrip(",")

        sql_insert = '''INSERT INTO Interaction_bitvectors (Pose_ID, {columns}) VALUES ({fillers})'''.format(columns = column_str, fillers = filler_str)

        try:
            cur = self.conn.cursor()
            cur.executemany(sql_insert, bitvectors)
            conn.commit()
            
        except Exception as e:
            print(e)

        cur.close()

    def __create_interaction_bv_table(self, num_discrete_interactions):
        interaction_columns = range(num_discrete_interactions)
        interact_columns_str = ""
        for i in interaction_columns:
            interact_columns_str += "Interaction_"+ str(i+1) + " INTEGER,\n"

        interact_columns_str = interact_columns_str.rstrip(",\n")

        bv_table = """CREATE TABLE Interaction_bitvectors (
        Pose_ID INTEGER PRIMARY KEY,
        {columns})""".format(columns = interact_columns_str)

        try:
            cur = self.conn.cursor()
            cur.execute(bv_table)
        except Exception as e:
            print(e)

        cur.close()

    def __insert_unique_interactions(self, unique_interactions):
        #rint("Inserting interactions...")
        sql_insert = '''INSERT INTO Interaction_indices (
        interaction_type,
        rec_chain,
        rec_resname,
        rec_resid,
        rec_atom,
        rec_atomid
        ) VALUES (?,?,?,?,?,?)'''

        try:
            cur = self.conn.cursor()
            cur.executemany(sql_insert, unique_interactions)
            conn.commit()
        
        except Exception as e:
            print(e)

        cur.close()

    def __insert_one_interaction(self, interaction):
        #print("Inserting interactions...")
        sql_insert = '''INSERT INTO Interaction_indices (
        interaction_type,
        rec_chain,
        rec_resname,
        rec_resid,
        rec_atom,
        rec_atomid
        ) VALUES (?,?,?,?,?,?)'''

        try:
            cur = self.conn.cursor()
            cur.execute(sql_insert, interaction)
            conn.commit()
    
        except Exception as e:
            print(e)

        cur.close()

    def __make_new_interaction_column(self, column_number):
        add_column_str = '''ALTER TABLE Interaction_bitvectors ADD COLUMN Interaction_{n_inter}'''.format(n_inter = str(column_number))
        try:
            cur = self.conn.cursor()
            cur.execute(add_column_str)
            conn.commit()
            
        except Exception as e:
            print(e)

        cur.close()


    ###########################
    ##### Private methods #####
    ###########################

    def __create_connection(self):
        con = None
        try:
            con = sqlite3.connect(self.db_file)
        except Exception as e:
            print(e)
        return con

    def __create_results_table(self):
        sql_results_table = """CREATE TABLE Results (
            Pose_ID             INTEGER PRIMARY KEY AUTOINCREMENT,
            LigName             VARCHAR NOT NULL,
            ligand_smile        VARCHAR[],
            pose_rank           INT[],
            run_number          INT[],
            energies_binding    FLOAT(4),          -- 8 x 4 Bytes = 32 Bytes
            leff                FLOAT(4),
            deltas              FLOAT(4),
            cluster_rmsd        FLOAT(4),
            reference_rmsd      FLOAT(4),
            energies_inter      FLOAT(4),
            energies_vdw        FLOAT(4),
            energies_electro    FLOAT(4),
            energies_flexLig    FLOAT(4),
            energies_flexLR     FLOAT(4),
            energies_intra      FLOAT(4),
            energies_torsional  FLOAT(4),
            unbound_energy      FLOAT(4),
            nr_interactions     INT[],          -- 2 Bytes                                      (-> 46 Bytes)
            num_hb              INT[],
            about_x             FLOAT(4),          -- 3 x 4 Bytes = 12 Bytes
            about_y             FLOAT(4),
            about_z             FLOAT(4),
            trans_x             FLOAT(4),          -- 3 x 4 Bytes = 12 Bytes
            trans_y             FLOAT(4),
            trans_z             FLOAT(4),
            axisangle_x         FLOAT(4),          -- 4 x 4 Bytes = 16 Bytes                       (-> 30 Bytes)
            axisangle_y         FLOAT(4),
            axisangle_z         FLOAT(4),
            axisangle_w         FLOAT(4),
            dihedrals           VARCHAR[]          -- 8-ish Bytes per dihedral value
        );
        -- total number of Bytes estimated per row for average 50 interactions and 10 torsions:
        --     (64 + 76 + 50 * 21 + 10 * 8) Bytes = 1,270 Bytes
        -- for 50 ga_runs per result: (1,270 * 50) Bytes = 63,500 Bytes (62 KB)
        -- for 10,000 results (typical package size) this amounts to about 62 KB * 10,000 = 591 MB"""

        create_index = "CREATE INDEX index_ligand_run on Results (LigName, run_number)"

        try:
            cur = self.conn.cursor()
            cur.execute(sql_results_table)
            cur.execute(create_index)
        except Exception as e:
            print(e)

        cur.close()

    def __create_ligands_table(self):
        ligand_table = """CREATE TABLE Ligands (
            LigName             VARCHAR NOT NULL PRIMARY KEY,
            ligand_smile        VARCHAR NOT NULL,
            input_pdbqt         VARCHAR[],
            best_binding        FLOAT(4),
            best_run            SMALLINT[])"""

        try:
            cur = self.conn.cursor()
            cur.execute(ligand_table)
        except Exception as e:
            print(e)

        cur.close()

    def __create_interaction_index_table(self):
        """create table of data about each unique interaction"""
        interaction_index_table = """CREATE TABLE Interaction_indices (
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
        except Exception as e:
            print(e)

        cur.close()

    def _add_energy(self, e_data):
        """  """
        # sql code
        pass

    def _add_leff(self, le_data):
        pass


    def _get_energy(self, e_data):
        """  """
        # sql code
        pass

    def _get_leff(self, le_data):
        pass

