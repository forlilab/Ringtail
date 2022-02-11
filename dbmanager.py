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
        self.interaction_tolerance_cutoff = self.opts["interaction_tolerance"]
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

        #keep track of any open cursors
        self.open_cursors = []

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
        return self._fetch_all_plot_data(), self._fetch_passing_plot_data()

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
        #close any open cursors
        for cur in self.open_cursors:
            cur.close()
        #close db itself
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
        

    def insert_interactions(self, interactions_list):

        self._add_unique_interactions(interactions_list)

        #check if we need to initialize the interaction bv table and insert first set of interaction
        if not self.interactions_initialized_flag:
            self._create_interaction_bv_table()
            self._insert_unique_interactions(np.array(list(self.unique_interactions.keys())))
            self.interactions_initialized_flag = True

        self._insert_interaction_bitvectors(self._generate_interaction_bitvectors(interactions_list))

    def get_top_energies_leffs(self):
        self._fetch_best_energies_leff()
        return self.energies, self.leffs, self.plot_data

    def filter_results(self, results_filters_list, ligand_filters_list, output_fields):
        #create view of passing results
        filter_results_str = self._generate_result_filtering_str_sqlite(results_filters_list, ligand_filters_list, output_fields)
        print(filter_results_str)
        self._create_view("passing_results", filter_results_str)
        #perform filtering
        filtered_results = self._run_query(filter_results_str)
        #get number of passing ligands
        return filtered_results

    def get_number_passing_ligands(self):
        cur = self.conn.cursor()
        cur.execute("SELECT COUNT(DISTINCT LigName) FROM passing_results")
        return cur.fetchone()

    def format_rows_from_dict(self, ligand_dict):
        """takes file dictionary from the mpreader, formats into rows for the database insertion"""

        #initialize row holders
        result_rows = []
        interaction_dictionaries = []
        interaction_tuples = []

        cluster_top_pose_runs = self._find_cluster_top_pose_runs(ligand_dict)
        if self.interaction_tolerance_cutoff != None:
            tolerated_interaction_runs = self._find_tolerated_interactions(ligand_dict)
        else:
            tolerated_interaction_runs = []

        for i in range(len(ligand_dict["sorted_runs"])):
            run_number = int(ligand_dict["sorted_runs"][i])
            if run_number in cluster_top_pose_runs: #save everything if this is a cluster top pose
                if interaction_tuples != []: #don't save interaction data from previous cluster for first cluster
                    pose_interactions = self._generate_interaction_tuples(interaction_dictionaries) #will generate tuples across all dictionaries for last cluster
                    interaction_tuples.append(pose_interactions)
                    if self.interaction_tolerance_cutoff != None and result_rows != []: #only update previous entry if we added tolerated interactions
                        result_rows[-1][17] = len(pose_interactions) + int(result_rows[-1][17]) #update number of interactions
                        result_rows[-1][18] = sum(1 for interaction in pose_interactions if interaction[0] == "H") + result_rows[-1][18] #count and update number of hydrogen bonds
                    interaction_dictionaries = [] #clear the list for the new cluster
                result_rows.append(self._generate_results_row(ligand_dict, i, run_number))
                interaction_dictionaries.append(ligand_dict["interactions"][i])
            elif run_number in tolerated_interaction_runs:
                interaction_dictionaries.append(ligand_dict["interactions"][i]) #adds to list started by best-scoring pose in cluster
        
        pose_interactions = self._generate_interaction_tuples(interaction_dictionaries) #will generate tuples across all dictionaries for last cluster
        interaction_tuples.append(pose_interactions)
        if self.interaction_tolerance_cutoff != None: #only update if we added interactions
            result_rows[-1][17] = len(pose_interactions) + int(result_rows[-1][17]) #update number of interactions
            result_rows[-1][18] = sum(1 for interaction in pose_interactions if interaction[0] == "H") + int(result_rows[-1][18]) #count and update number of hydrogen bonds

        if len(interaction_tuples) > 3:
            print("HELP", len(interaction_tuples))

        return (result_rows, self._generate_ligand_row(ligand_dict), interaction_tuples)

    ###########################
    ##### Private methods #####
    ###########################

    def _create_connection(self):
        con = None
        try:
            con = sqlite3.connect(self.db_file)
            cursor = con.cursor()
            cursor.execute('PRAGMA synchronous = OFF')
            cursor.execute('PRAGMA journal_mode = MEMORY')
            cursor.close() 
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
        self.open_cursors.append(cur)
        return cur

        
        """conn = sqlite3.connect('test2_3.db')
        cursor=conn.cursor()
        cursor.execute(query)

        for row in cursor:
            print(row)"""

    def _create_view(self, name, query):
        """takes name and selection query, creates view of query."""
        cur = self.conn.cursor()
        #drop old view if there is one
        cur.execute("DROP VIEW IF EXISTS {name}".format(name = name))
        cur.execute("CREATE VIEW {name} AS {query}".format(name = name, query = query))
        cur.close()



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

        try:
            cur = self.conn.cursor()
            cur.execute(sql_results_table)
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

    def _fetch_all_plot_data(self):
        return self._run_query(self._generate_plot_all_results_query())

    def _generate_plot_all_results_query(self):

        return "SELECT energies_binding, leff FROM Results GROUP BY LigName"

    def _fetch_passing_plot_data(self):

        return self._run_query(self._generate_plot_passing_results_query())

    def _generate_plot_passing_results_query(self):

        return "SELECT energies_binding, leff FROM Results WHERE LigName IN (SELECT DISTINCT LigName FROM passing_results) GROUP BY LigName"

    def _generate_result_filtering_str_sqlite(self, results_filters_list, ligand_filters_list, output_fields):
        """ takes list of filters, writes sql filtering string"""

        #parse requested output fields and convert to column names in database
        
        outfield_string = "LigName, " + ", ".join([self.field_to_column_name[field] for field in output_fields.split(",") if field in self.field_to_column_name])

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
            sql_string += " AND Pose_ID IN ({interaction_str})".format(interaction_str = self._write_interaction_filtering_str(interaction_filter_indices))

        #add ligand filters
        if ligand_filters_list != []:
            sql_string += " AND LigName IN ({ligand_str})".format(ligand_str = self._write_ligand_filtering_sql(ligand_filters_list))

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

        sql_ligand_string = "SELECT LigName FROM Ligands WHERE "

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

    def _generate_interaction_tuples(self, interaction_dictionaries):
        """takes dictionary of file results, formats list of tuples for interactions"""
        interactions = set()
        for pose_interactions in interaction_dictionaries:
            count = pose_interactions["count"][0]
            for i in range(int(count)):
                interactions.add(tuple(pose_interactions[kw][i] for kw in self.interaction_data_kws))

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

    def _generate_interaction_bitvectors(self, interactions_list):
        """takes string of interactions and all unique interactions and makes bitvector"""
        bitvectors_list = []
        for pose_interactions in interactions_list:
            pose_bitvector = [None]*len(self.unique_interactions)
            for interaction_tuple in pose_interactions:
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

    def _find_tolerated_interactions(self, ligand_dict):
        """take ligand dict and finds which poses we should save the interactions for as tolerated interactions"""
        tolerated_runs = []
        for i in range(len(ligand_dict["sorted_runs"])):
            if float(ligand_dict["cluster_rmsds"][i]) <= self.interaction_tolerance_cutoff:
                tolerated_runs.append(ligand_dict["sorted_runs"][i])
        return tolerated_runs
