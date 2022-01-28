from dbmanager import DBManager

class DBManager():
    """ DOCUMENTATION GOES HERE """
    def __init__(self, opts={}):
        self.fname = opts['sqlFile']
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

    def get_top_energies_leffs(self):
        self.__fetch_best_energies_leff()
        return self.energies, self.leffs, self.plot_data

    def filter_results(self, results_filters_list, output_fields):
        filter_results_str = self.__write_result_filtering_str_sqlite(results_filters_list, output_fields)
        filtered_results = self.__select_from_db(filter_results_str)
        return filtered_results

    def filter_ligands(self, ligand_filters):
        filter_ligands_str = self.__write_ligand_filtering_sql(ligand_filters)
        filtered_ligands = self.__select_from_db(filter_ligands_str)
        return filtered_ligands

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

    def __select_from_db(self, filter_sql_string):
        cur = self.conn.cursor()
        cur.execute(filter_sql_string)

        rows = cur.fetchall()

        cur.close()

        return rows

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

    def __fetch_best_energies_leff(self):
        conn = self.conn
        plot_data_str = self.__write_plot_results_sql()
        plot_data = self.select_from_db(plot_data_str)

        energies = []
        leffs = []
        for row in plot_data:
            if int(row[3]) == 1: #only one point per ligand, for best binding energy
                energies.append(row[1])
                leffs.append(row[2])

        self.energies = energies
        self.leffs = leffs
        self.plot_data = plot_data
        self.fetch_energies_le_flag = True

    def __write_plot_results_sql(self):
        sql_string = "SELECT LigName, energies_binding, leff, pose_rank FROM Results"

        return sql_string

    def __write_result_filtering_str_sqlite(self, results_filters_list, output_fields):
        """ takes list of filters, writes sql filtering string"""

        #parse requested output fields and convert to column names in database
        outfield_dict = {"e":"energies_binding",
                        "le":"leff",
                        "delta":"deltas",
                        "ref_rmsd":"reference_rmsd",
                        "e_inter":"energies_inter",
                        "e_vdw":"energies_vdw",
                        "e_elec":"energies_electro",
                        "e_intra":"energies_intra",
                        "n_interact":"nr_interactions",
                        "interactions":"interactions",
                        "fname":"LigName",
                        "ligand_smile":"ligand_smile",
                        "rank":"pose_rank",
                        "run":"run_number",
                        "hb":"num_hb"}
        outfield_string = ""
        for field in output_fields:
            if field in outfield_dict.keys():
                outfield_string += outfield_dict[field] + ", "
        outfield_string = outfield_string.rstrip(", ")

        interaction_filters = []

        sql_string = """SELECT {out_columns} FROM Results WHERE """.format(out_columns = outfield_string)

        energy_filter_sql_call_dict = {"eworst":"energies_binding < {value} AND ",
                        "ebest":"energies_binding > {value} AND ",
                        "leworst":"leff < {value} AND ",
                        "lebest":"leff > {value} AND "
                        }

        for filter_tuple in results_filters_list:
            filter_key = filter_tuple[0]
            filter_value = filter_tuple[1]

            #write energy filters
            if filter_key in energy_filter_sql_call_dict:
                energy_fil_str = energy_filter_sql_call_dict[filter_key].format(value = filter_value)
                sql_string += energy_fil_str
            """if filter_key == 'eworst':
                eworst_str = "energies_binding < {value} AND ".format(value = filter_value)
                sql_string += eworst_str

            if filter_key == 'ebest':
                ebest_str = "energies_binding > {value} AND ".format(value = filter_value)
                sql_string += ebest_str

            if filter_key == 'leworst':
                leworst_str = "leff < {value} AND ".format(value = filter_value)
                sql_string += leworst_str

            if filter_key == 'lebest':
                lebest_str = "leff > {value} AND ".format(value = filter_value)
                sql_string += lebest_str

            if filter_key == 'epercentile':
                requested_percentile = filter_value
                self.calculate_energy_percentile(requested_percentile)
                epercentile_string = "energies_binding < {value} AND ".format(value = self.energy_percentile)
                sql_string += epercentile_string

            if filter_key == 'leffpercentile':
                requested_percentile = filter_value
                self.calculate_leff_percentile(requested_percentile)
                leff_percentile_string = "leff < {value} AND ".format(value = self.leff_percentile)
                sql_string += leff_percentile_string"""

            #write interaction filters
            if filter_key == 'hb_count':
                if filter_value > 0:
                    hb_count_str = "num_hb > {value} AND ".format(value = filter_value)
                    sql_string += hb_count_str
                else:
                    hb_count_str = "num_hb < {value} AND ".format(value = -1*filter_value)
                    sql_string += hb_count_str

            interaction_filter_keys = ["V",
            "H",
            "R"]
            #compile list of interactions to search for
            for key in interaction_filter_keys:
                if filter_key == key:
                    for interact in filter_value:
                        interaction_string = key + ":" + interact[0]
                        interaction_filters.append(interaction_string.split(":"))

        interaction_filter_indices = []
        #for each interaction, get the index from the interactions_indices table
        for interaction in interaction_filters:
            interact_index_str = self.__write_interaction_index_filtering_str(interaction)
            interaction_indices = self.__select_from_db(interact_index_str)
            for i in interaction_indices:
                interaction_filter_indices.append(i[0])

        #find pose ids for ligands with desired interactions
        if interaction_filter_indices != []:
            interaction_filter_str = self.__write_interaction_filtering_str(interaction_filter_indices)
            sql_string += "Pose_ID IN (" + interaction_filter_str + ")"

        if sql_string.endswith("AND "):
            sql_string = sql_string.rstrip("AND ")
        if sql_string.endswith("OR "):
            sql_string = sql_string.rstrip("OR ")

        return sql_string

    def __write_interaction_index_filtering_str(self, interaction_list):
        """takes list of interaction info for a given ligand, looks up corresponding interaction index"""
        interaction_info = ["interaction_type", "rec_chain", "rec_resname", "rec_resid", "rec_atom"]
        sql_string = """SELECT interaction_id FROM Interaction_indices WHERE """

        for i in range(4):
            item = interaction_list[i]
            column_name = interaction_info[i]
            if item != "":
                sql_string += "{column} LIKE '%{value}%' AND ".format(column = column_name, value = item)

        sql_string = sql_string.rstrip("AND ")
        return sql_string

    def __write_interaction_filtering_str(self, interaction_index_list):
        """takes list of interaction indices and searches for ligand ids which have those interactions"""

        sql_string = """SELECT Pose_id FROM Interaction_bitvectors WHERE """

        for index in interaction_index_list:
            add_str = "Interaction_{index_n} = 1 OR ".format(index_n = index)
            sql_string += add_str

        sql_string = sql_string.rstrip("OR ")

        return sql_string

    def __write_ligand_filtering_sql(self, ligand_filters):
        """write string to select from ligand table"""

        sql_ligand_string = "SELECT LigName, best_binding FROM Ligands WHERE "

        substruct_flag = ligand_filters['F'][0].upper()
        for kw in ligand_filters.keys():
            fils = ligand_filters[kw]
            if kw == 'N':
                for name in fils:
                    name_sql_str = "LigName LIKE '%{value}%' OR ".format(value = name)
                    sql_ligand_string += name_sql_str
            if kw == "S":
                for substruct in fils:
                    substruct_sql_str = "ligand_smile LIKE '%{value}%' {flag} ".format(value = substruct, flag = substruct_flag)
                    sql_ligand_string += substruct_sql_str

        if sql_ligand_string.endswith("AND "):
            sql_ligand_string = sql_ligand_string.rstrip("AND ")
        if sql_ligand_string.endswith("OR "):
            sql_ligand_string = sql_ligand_string.rstrip("OR ")

        return sql_ligand_string

