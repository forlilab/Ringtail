import sqlite3
import numpy as np

class DBManager():
    """ Prototype class for a generic VS database object
        this class defines the API of each DBManager object (including sub-classess)
        which will implement their own functions to return the data requested
    """
    def __init__(self, opts={}):
        self.db_file = opts['sqlFile']
        self.opts = opts
        self.order_results = self.opts['order_results']
        self.log_distinct_ligands = self.opts["log_distinct_ligands"]
        self.write_db_flag = self.opts["write_db_flag"]
        self.num_clusters = self.opts["num_clusters"]
        #initialize dictionary processing kw lists
        self.interaction_data_kws = ["type", "chain", "residue", "resid", "recname", "recid"]
        self.ligand_data_keys = ["cluster_rmsds",
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
        "unbound_energy"]
        self.ligand_interaction_keys = ["type",
        "chain",
        "residue",
        "resid",
        "recname",
        "recid"]
        self.stateVar_keys = ["pose_about",
        "pose_translations",
        "pose_quarternions"]

    def get_results(self):
        """ generic function for retrieving results"""
        raise NotImplemented


    def get_top_scores(self):
        """ this function is expected to return [data in which format?] """
        raise NotImplemented

    def get_plot_data(self):
        """ this function is expected to return an ascii plot representation of the results"""
        # TODO this function could be actually implemented here, if the
        # plotting mechanism will be common to all the child classes, too
        # for example, this function could contain the actuall call to the ASCII library,
        # and call a _fetch_plot_data() function that will be implemented in each child class
        return self._fetch_plot_data()

    def prune(self):
        """ do we want a method for deleting rows not satisfying a given requirement? """
        raise NotImplemented

    def clone(self):
        """ this function would be useful for creating a copy of the db that can be then pruned?
            (i.e., a lightweight version of the DB containing only the distilled information)
        """
        pass

    def close_connection(self):
        """close connection to database"""
        self._close_connection()

class DBManagerSQLite(DBManager):
    """ DOCUMENTATION GOES HERE """
    def __init__(self, opts = {}):
        super().__init__(opts)
        self.conn = self._create_connection()

        self.unique_interactions = {}
        self.next_unique_interaction_idx = 1
        self.interactions_initialized_flag = False

        self.field_to_column_name = {"e":"energies_binding",
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

        self.energy_filter_sqlite_call_dict = {"eworst":"energies_binding < {value}",
                        "ebest":"energies_binding > {value}",
                        "leworst":"leff < {value}",
                        "lebest":"leff > {value}",
                        "epercentile":"energy_percentile_rank < {value}",
                        "leffpercentile":"leff_percentile_rank < {value}"
                        }

        self.interaction_filter_types = ["V",
            "H",
            "R"]

        self._initialize_db()

    ##########################
    ##### Public methods #####
    ##########################
    def insert_results(self, results_array):
        """takes array of database rows to insert, adds data to results table
        Inputs:
        -conn: database connection
        -ligand_array: numpy array of arrays"""

        #print("Inserting results...")
        sql_insert = """INSERT INTO Results (
        LigName,
        ligand_smile,
        pose_rank,
        run_number,
        cluster_rmsd,
        reference_rmsd,
        energies_binding,
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
        ) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)"""

        try:
            cur = self.conn.cursor()
            cur.executemany(sql_insert, results_array)
            self.conn.commit()
            cur.close()
            
        except Exception as e:
            print(e)
            raise e

    def insert_ligands(self, ligand_array):
        #print("Inserting ligand data...")
        sql_insert = '''INSERT INTO Ligands (
        LigName,
        ligand_smile,
        input_pdbqt,
        best_binding,
        best_run
        ) VALUES
        (?,?,?,?,?)'''

        try:
            cur = self.conn.cursor()
            cur.executemany(sql_insert, ligand_array)
            self.conn.commit()
            cur.close()
            
        except Exception as e:
            print(e)
            raise e
        

    def insert_interactions(self, pose_id_list, interactions_list):

        self._add_unique_interactions(interactions_list)

        #check if we need to initialize the interaction bv table and insert first set of interaction
        if not self.interactions_initialized_flag:
            self._create_interaction_bv_table()
            self._insert_unique_interactions(np.array(list(self.unique_interactions.keys())))
            self.interactions_initialized_flag = True

        self._insert_interaction_bitvectors(self._generate_interaction_bitvectors(pose_id_list, interactions_list))

    def get_top_energies_leffs(self):
        self._fetch_best_energies_leff()
        return self.energies, self.leffs, self.plot_data

    def filter_results(self, results_filters_list, output_fields):
        filter_results_str = self._generate_result_filtering_str_sqlite(results_filters_list, output_fields)
        print(filter_results_str)
        filtered_results = self._run_query(filter_results_str)
        print(filtered_results)
        return filtered_results

    def filter_ligands(self, ligand_filters):
        filter_ligands_str = self._write_ligand_filtering_sql(ligand_filters)
        filtered_ligands = self._run_query(filter_ligands_str)
        return filtered_ligands

    def format_rows_from_dict(self, ligand_dict):
        """takes file dictionary from the mpreader, formats into rows for the database insertion"""

        #initialize row holders
        result_rows = []
        interaction_tuples = []

        for i in range(len(ligand_dict["sorted_runs"])):
            cluster_top_pose_runs = self._find_cluster_top_pose_runs(ligand_dict)
            run_number = int(ligand_dict["sorted_runs"][i])
            if run_number in cluster_top_pose_runs:
                result_rows.append(self._generate_results_row(ligand_dict, i, run_number))
                interaction_tuples.append(self._generate_interaction_tuples(ligand_dict["interactions"][i]))

        return (result_rows, self._generate_ligand_row(ligand_dict), interaction_tuples)

    ###########################
    ##### Private methods #####
    ###########################

    def _create_connection(self):
        con = None
        try:
            con = sqlite3.connect(self.db_file)
        except Exception as e:
            print(e)
            raise e
        return con

    def _close_connection(self):
        print("Closing database")
        self.conn.close()

    def _initialize_db(self):
        """check if db needs to be written. If so, initialize the tables"""
        if not self.write_db_flag:
            return
        #create tables in db
        self._create_results_table()
        self._create_ligands_table()
        self._create_interaction_index_table()


    def _run_query(self, query):
        """self.conn = self._create_connection()
        cur = self.conn.cursor()
        cur.execute(query+";")

        rows = cur.fetchall()
        print(rows)

        cur.close()
        self.conn.close()

        return rows"""
        cur = self.conn.cursor()
        cur.execute(query)
        return cur

        
        """conn = sqlite3.connect('test2_3.db')
        cursor=conn.cursor()
        cursor.execute(query)

        for row in cursor:
            print(row)"""


    def _create_results_table(self):
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
            #cur.execute(create_index)
            cur.close()
        except Exception as e:
            print(e)
            raise e

        

    def _create_ligands_table(self):
        ligand_table = """CREATE TABLE Ligands (
            LigName             VARCHAR NOT NULL PRIMARY KEY,
            ligand_smile        VARCHAR NOT NULL,
            input_pdbqt         VARCHAR[],
            best_binding        FLOAT(4),
            best_run            SMALLINT[])"""

        try:
            cur = self.conn.cursor()
            cur.execute(ligand_table)
            cur.close()
        except Exception as e:
            print(e)
            raise e

        

    def _create_interaction_index_table(self):
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
            cur.close()
        except Exception as e:
            print(e)
            raise e

        

    def _create_interaction_bv_table(self):
        interact_columns_str = " INTEGER,\n".join(["Interaction_"+ str(i+1) for i in range(len(self.unique_interactions))]) + " INTEGER"

        bv_table = """CREATE TABLE Interaction_bitvectors (
        Pose_ID INTEGER PRIMARY KEY AUTOINCREMENT,
        {columns})""".format(columns = interact_columns_str)

        try:
            cur = self.conn.cursor()
            cur.execute(bv_table)
            cur.close()
        except Exception as e:
            print(e)
            raise e

        

    def _insert_unique_interactions(self, unique_interactions):
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
            self.conn.commit()
            cur.close()
        
        except Exception as e:
            print(e)
            raise e

        

    def _insert_one_interaction(self, interaction):
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
            self.conn.commit()
            cur.close()
    
        except Exception as e:
            print(e)
            raise e

        

    def _make_new_interaction_column(self, column_number):
        add_column_str = '''ALTER TABLE Interaction_bitvectors ADD COLUMN Interaction_{n_inter}'''.format(n_inter = str(column_number))
        try:
            cur = self.conn.cursor()
            cur.execute(add_column_str)
            self.conn.commit()
            cur.close()
            
        except Exception as e:
            print(e)
            raise e

        

    def _fetch_plot_data(self):
        conn = self.conn
        plot_data_query = self._generate_plot_results_query()
        return self._run_query(plot_data_query)

    def _generate_plot_results_query(self):

        return "SELECT energies_binding, leff FROM Results GROUP BY LigName"

    def _generate_result_filtering_str_sqlite(self, results_filters_list, output_fields):
        """ takes list of filters, writes sql filtering string"""

        #parse requested output fields and convert to column names in database
        
        outfield_string = ", ".join([self.field_to_column_name[field] for field in output_fields.split(",") if field in self.field_to_column_name])

        interaction_filters = []

        #write energy filters
        energy_filter_sql_query = []
        filtering_window = "Results"
        for filter_key, filter_value in results_filters_list:
            if filter_key in self.energy_filter_sqlite_call_dict:
                if filter_key == "epercentile" or filter_key == "leffpercentile":
                    #convert from percent to decimal
                    filter_value = str(float(filter_value)/100)
                    #reset filtering window to include generated percentile_ranks
                    filtering_window = self._generate_percentile_rank_window()
                energy_filter_sql_query.append(self.energy_filter_sqlite_call_dict[filter_key].format(value = filter_value))
                    
            
            #write hb count filter(s)
            if filter_key == 'hb_count':
                if filter_value > 0:
                    energy_filter_sql_query.append("num_hb > {value}".format(value = filter_value))
                else:
                    energy_filter_sql_query.append("num_hb < {value}".format(value = -1*filter_value))

        sql_string = """SELECT {out_columns} FROM {window} WHERE """.format(out_columns = outfield_string, window = filtering_window)
        sql_string += " AND ".join(energy_filter_sql_query)

        #compile list of interactions to search for
        for key in self.interaction_filter_types:
            if filter_key == key:
                for interact in filter_value:
                    interaction_string = key + ":" + interact[0]
                    interaction_filters.append(interaction_string.split(":"))

        interaction_filter_indices = []
        #for each interaction, get the index from the interactions_indices table
        for interaction in interaction_filters:
            interact_index_str = self._write_interaction_index_filtering_str(interaction)
            interaction_indices = self._run_query(interact_index_str)
            for i in interaction_indices:
                interaction_filter_indices.append(i[0])

        #find pose ids for ligands with desired interactions
        if interaction_filter_indices != []:
            interaction_filter_str = self._write_interaction_filtering_str(interaction_filter_indices)
            sql_string += " AND Pose_ID IN (" + interaction_filter_str + ")"

       #adding if we only want to keep one pose per ligand (will keep first entry)
        if self.log_distinct_ligands:
            sql_string += " GROUP BY LigName"

        #add how to order results
        if self.order_results != None:
            try:
                sql_string += " ORDER BY " + self.field_to_column_name[self.order_results]
            except KeyError:
                print("Please ensure you are only requesting one option for --order_results and have written it correctly")
                raise KeyError

        return sql_string

    def _write_interaction_index_filtering_str(self, interaction_list):
        """takes list of interaction info for a given ligand, looks up corresponding interaction index"""
        interaction_info = ["interaction_type", "rec_chain", "rec_resname", "rec_resid", "rec_atom"]
        sql_string = """SELECT interaction_id FROM Interaction_indices WHERE """

        sql_string += " AND ".join(["{column} LIKE '%{value}%'".format(column = interaction_info[i], value = interaction_list[i]) for i in range(4)]) 

        return sql_string

    def _write_interaction_filtering_str(self, interaction_index_list):
        """takes list of interaction indices and searches for ligand ids which have those interactions"""

        return """SELECT Pose_id FROM Interaction_bitvectors WHERE """ + " OR ".join(["Interaction_{index_n} = 1".format(index_n = index) for index in interaction_index_list])

    def _write_ligand_filtering_sql(self, ligand_filters):
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

    def _find_cluster_top_pose_runs(self, ligand_dict):
        """returns list of the run numbers for the top run in the top self.num_clusters clusters"""
        try:
            cluster_top_pose_runs = ligand_dict["cluster_top_poses"][:self.num_clusters] #will only select top n clusters. Default 3
        except IndexError:
            cluster_top_pose_runs = ligand_dict["cluster_top_poses"] #catch indexerror if not enough clusters for given ligand

        return cluster_top_pose_runs

    def _generate_results_row(self, ligand_dict, pose_rank, run_number):
        """generate list of lists of ligand values to be inserted into sqlite database"""

        ######get pose-specific data
        
        #check if run is best for a cluster. We are only saving the top pose for each cluster
        ligand_data_list = [ligand_dict["ligname"], ligand_dict["ligand_smile_string"], pose_rank+1, run_number]
        #get energy data
        for key in self.ligand_data_keys:
            ligand_data_list.append(ligand_dict[key][pose_rank])

        #add interaction count
        ligand_data_list.append(ligand_dict["interactions"][pose_rank]["count"][0]) 
        #count number H bonds, add to ligand data list
        ligand_data_list.append(ligand_dict["interactions"][pose_rank]["type"].count("H"))

        for key in self.stateVar_keys:
            stateVar_data = ligand_dict[key][pose_rank]
            for dim in stateVar_data:
                ligand_data_list.append(dim)
        pose_dihedrals = ligand_dict["pose_dihedrals"][pose_rank]
        dihedral_string = ""
        for dihedral in pose_dihedrals:
            dihedral_string = dihedral_string + str(dihedral) + ", "
        ligand_data_list.append(dihedral_string)

        return ligand_data_list

    def _generate_ligand_row(self, ligand_dict):
        """writes row to be inserted into ligand table"""
        ligand_name = ligand_dict["ligname"]
        ligand_smile = ligand_dict["ligand_smile_string"]
        input_pdbqt = "\n".join(ligand_dict["ligand_input_pdbqt"])
        best_binding = ligand_dict["scores"][0]
        best_run = ligand_dict["sorted_runs"][0]

        return [ligand_name, ligand_smile, input_pdbqt, best_binding, best_run]

    def _generate_interaction_tuples(self, interaction_dictionary):
        """takes dictionary of file results, formats list of tuples for interactions"""
        self.count = interaction_dictionary["count"][0]
        interactions = []
        for i in range(int(self.count)):
            interactions.append(tuple(interaction_dictionary[kw][i] for kw in self.interaction_data_kws))

        return interactions

    def _add_unique_interactions(self, interactions_list):
        """takes list of interaction tuple lists. Examines self.unique_interactions, add interactions if not already inserted.

        self.unique_interactions {(interaction tuple): unique_interaction_idx}"""

        for pose in interactions_list:
            for interaction_tuple in pose:
                if interaction_tuple not in self.unique_interactions:
                    self.unique_interactions[interaction_tuple] = self.next_unique_interaction_idx
                    if self.interactions_initialized_flag:
                        self._insert_one_interaction(interaction_tuple)
                        self._make_new_interaction_column(self.next_unique_interaction_idx)
                    self.next_unique_interaction_idx += 1

    def _generate_interaction_bitvectors(self, pose_id_list, interactions_list):
        """takes string of interactions and all unique interactions and makes bitvector"""
        bitvectors_list = []
        for i in range(len(pose_id_list)):
            pose_bitvector = [None]*len(self.unique_interactions)
            for interaction_tuple in interactions_list[i]:
                interaction_idx = self.unique_interactions[interaction_tuple]
                pose_bitvector[interaction_idx - 1] = 1 #index corrected for interaction indices starting at 1 in sqlite table

            bitvectors_list.append(pose_bitvector)

        return bitvectors_list

    def _insert_interaction_bitvectors(self, bitvectors):
        #print("Inserting interaction bitvectors...")
        interaction_columns = range(len(self.unique_interactions))
        column_str = ""
        filler_str = ""
        for i in interaction_columns:
            column_str += "Interaction_"+ str(i+1) + ", "
            filler_str += "?,"
        column_str = column_str.rstrip(", ")
        filler_str = filler_str.rstrip(",")

        sql_insert = '''INSERT INTO Interaction_bitvectors ({columns}) VALUES ({fillers})'''.format(columns = column_str, fillers = filler_str)

        try:
            cur = self.conn.cursor()
            cur.executemany(sql_insert, bitvectors)
            self.conn.commit()
            cur.close()
            
        except Exception as e:
            print(e)
            raise e

    def _generate_percentile_rank_window(self):
        """makes window with percentile ranks for percentile filtering"""
        column_names = ",".join(self._fetch_results_column_names())
        return "SELECT {columns}, PERCENT_RANK() OVER (ORDER BY energies_binding) energy_percentile_rank, PERCENT_RANK() OVER (ORDER BY leff) leff_percentile_rank FROM Results".format(columns = column_names)

    def _fetch_results_column_names(self):
        return [column_tuple[1] for column_tuple in self.conn.execute("PRAGMA table_info(Results)")]




