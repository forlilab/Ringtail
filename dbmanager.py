"""Summary
"""
import sqlite3
import json


class DBManager():
    """Prototype class for a generic VS database object
    this class defines the API of each
    DBManager class (including sub-classess)
    which will implement their own functions to return the data requested

    Attributes:
        db_file (string): Name of file containing database
        field_to_column_name (dictionary): Dictionary for
            converting command-line field options into DB column names
        interaction_data_kws (list): List of keywords for different
            pieces of interaction data
        interaction_filter_types (set): Set for types of interaction filters
        interaction_tolerance_cutoff (float): RMSD cutoff for interactions
            to be added to the top pose for stored clusters
        interactions_initialized_flag (boolean): Flag indicating if
            interaction tables have been created
        ligand_data_keys (list): List of keywords used to look up ligand
            data in ligand dictionaries from parser
        log_distinct_ligands (boolean): Flag dictating, if true, only top
            pose for given ligand will be stored
        next_unique_interaction_idx (int): Index for the next unique
            interaction to be added to interaction_index table
        num_clusters (int): Number of ligand clusters that top pose
            should be stored for.
        open_cursors (list): Storage for any DB cursors which are opened
            and not closed by the function that opened them. Will be closed by
            close_connection method.
        opts (dict): Dictionary of database options
        order_results (string): string indicating result field that passing
            results should be ordered by in log.
        overwrite_flag (boolean): Flag indictating that DBMan should drop
            all existing tables and add nes data
        passing_results_view_name (string): Name for the view of passing
            results to be created after filtering
        stateVar_keys (list): List of strings for the different state variables
        store_all_poses_flag (boolean): Flag indicating that all poses should
            be stored, not just the top pose from top N clusters
            where N=self.num_clusters
        unique_interactions (dict): Dictionary for storing unique interactions
            to be written in interaction_index table
        write_db_flag (boolean): Flag indicating that DBMan will be
            writing new data to DB
    """

    def __init__(self, opts={}):
        """Initialize instance variables common to all DBMan subclasses

        Args:
            opts (dict): Dictionary of database options
        """
        self.db_file = opts['sqlFile']
        self.opts = opts
        self.order_results = self.opts['order_results']
        self.log_distinct_ligands = self.opts["log_distinct_ligands"]
        self.write_db_flag = self.opts["write_db_flag"]
        self.num_clusters = self.opts["num_clusters"]
        self.interaction_tolerance_cutoff = self.opts["interaction_tolerance"]
        self.passing_results_view_name = self.opts["results_view_name"]
        self.store_all_poses_flag = self.opts["store_all_poses"]
        self.overwrite_flag = self.opts["overwrite"]
        # initialize dictionary processing kw lists
        self.interaction_data_kws = [
            "type", "chain", "residue", "resid", "recname", "recid"
        ]
        self.ligand_data_keys = [
            "cluster_rmsds", "ref_rmsds", "scores", "leff", "delta",
            "intermolecular_energy", "vdw_hb_desolv", "electrostatics",
            "flex_ligand", "flexLigand_flexReceptor", "internal_energy",
            "torsional_energy", "unbound_energy"
        ]
        """self.ligand_interaction_keys = ["type",
        "chain",
        "residue",
        "resid",
        "recname",
        "recid"]"""
        self.stateVar_keys = [
            "pose_about", "pose_translations", "pose_quarternions"
        ]

        self.unique_interactions = {}
        self.next_unique_interaction_idx = 1
        self.interactions_initialized_flag = False

        self.field_to_column_name = {
            "e": "energies_binding",
            "le": "leff",
            "delta": "deltas",
            "ref_rmsd": "reference_rmsd",
            "e_inter": "energies_inter",
            "e_vdw": "energies_vdw",
            "e_elec": "energies_electro",
            "e_intra": "energies_intra",
            "n_interact": "nr_interactions",
            "interactions": "interactions",
            "ligand_smile": "ligand_smile",
            "rank": "pose_rank",
            "run": "run_number",
            "hb": "num_hb",
            "source_file": "source_file"
        }
        self.interaction_filter_types = {"V", "H", "R"}

        # keep track of any open cursors
        self.open_cursors = []

        self._initialize_db()

    # # # # # # # # # # # # # # # # # # #
    # # # Common DBManager methods # # #
    # # # # # # # # # # # # # # # # # # #

    def get_plot_data(self):
        """this function is expected to return an ascii plot
        representation of the results

        Returns:
            Tuple: cursors as [<all data cursor>, <passing data cursor>]
        """

        # checks if we have filtered by looking for view name in list of view names
        if self.check_passing_view_exists():
            return self._fetch_all_plot_data(), self._fetch_passing_plot_data()
        else:
            return self._fetch_all_plot_data(), []

    def prune(self):
        """Deletes rows from results, ligands, and interactions
        if they do not pass filtering criteria
        """
        self._delete_from_results()
        self._delete_from_ligands()
        self._delete_from_interactions()

    def check_passing_view_exists(self):
        """Return if self.passing_results_view_name in database
        """
        return self.passing_results_view_name in [name[0] for name in self._fetch_view_names().fetchall()]

    def close_connection(self):
        """close connection to database
        """
        # close any open cursors
        self._close_open_cursors()
        # close db itself
        self._close_connection()

    def _find_cluster_top_pose_runs(self, ligand_dict):
        """returns list of the run numbers for the top run in the
        top self.num_clusters clusters

        Args:
            ligand_dict (Dictionary): Dictionary of ligand data from parser

        Returns:
            List: List of run numbers to save, which are the top runs from
            the first self.num_clusters clusters
        """
        try:
            # will only select top n clusters. Default 3
            cluster_top_pose_runs = ligand_dict[
                "cluster_top_poses"][:self.num_clusters]
        except IndexError:
            # catch indexerror if not enough clusters for given ligand
            cluster_top_pose_runs = ligand_dict[
                "cluster_top_poses"]

        return cluster_top_pose_runs

    def _generate_results_row(self, ligand_dict, pose_rank, run_number):
        """generate list of lists of ligand values to be
            inserted into sqlite database

        Args:
            ligand_dict (Dictionary): Dictionary of ligand data from parser
            pose_rank (int): Rank of pose to generate row for
                all runs for the given ligand
            run_number (int): Run number of pose to generate row for
                all runs for the given ligand

        Returns:
            List: List of pose data to be inserted into Results table.
            In same order as expected in _insert_results:
            LigName,
            ligand_smile,
            source_file,
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
            dihedrals,
            ligand_coordinates,
            flexible_residues,
            flexible_res_coordinates
        """

        # # # # # # get pose-specific data

        # check if run is best for a cluster.
        # We are only saving the top pose for each cluster
        ligand_data_list = [
            ligand_dict["ligname"], ligand_dict["ligand_smile_string"],
            ligand_dict["source_file"], pose_rank + 1, run_number
        ]
        # get energy data
        for key in self.ligand_data_keys:
            ligand_data_list.append(ligand_dict[key][pose_rank])

        # add interaction count
        ligand_data_list.append(
            ligand_dict["interactions"][pose_rank]["count"][0])
        # count number H bonds, add to ligand data list
        ligand_data_list.append(
            ligand_dict["interactions"][pose_rank]["type"].count("H"))

        # add statevars
        for key in self.stateVar_keys:
            stateVar_data = ligand_dict[key][pose_rank]
            for dim in stateVar_data:
                ligand_data_list.append(dim)
        pose_dihedrals = ligand_dict["pose_dihedrals"][pose_rank]
        dihedral_string = ""
        for dihedral in pose_dihedrals:
            dihedral_string = dihedral_string + json.dumps(dihedral) + ", "
        ligand_data_list.append(dihedral_string)

        # add coordinates
        # convert to string for storage as VARCHAR
        ligand_data_list.append(
            json.dumps(ligand_dict["pose_coordinates"]
                       [pose_rank]))
        ligand_data_list.append(json.dumps(ligand_dict["flexible_residues"]))
        ligand_data_list.append(
            json.dumps(ligand_dict["flexible_res_coordinates"][pose_rank]))

        return ligand_data_list

    def _generate_ligand_row(self, ligand_dict):
        """writes row to be inserted into ligand table

        Args:
            ligand_dict (Dictionary): Dictionary of ligand data from parser

        Returns:
            List: List of data to be written as row in ligand table. Format:
            [ligand_name, ligand_smile, ligand_index_map,
            ligand_h_parents, input_pdbqt]
        """
        ligand_name = ligand_dict["ligname"]
        ligand_smile = ligand_dict["ligand_smile_string"]
        ligand_index_map = json.dumps(ligand_dict["ligand_index_map"])
        ligand_h_parents = json.dumps(ligand_dict["ligand_h_parents"])
        input_pdbqt = json.dumps(ligand_dict["ligand_input_pdbqt"])

        return [
            ligand_name, ligand_smile, ligand_index_map, ligand_h_parents,
            input_pdbqt
        ]

    def _generate_interaction_tuples(self, interaction_dictionaries):
        """takes dictionary of file results, formats as
        list of tuples for interactions

        Args:
            interaction_dictionaries (List): List of pose interaction
            dictionaries from parser

        Returns:
            List: List of tuples of interaction data
        """
        interactions = set()
        for pose_interactions in interaction_dictionaries:
            count = pose_interactions["count"][0]
            for i in range(int(count)):
                interactions.add(
                    tuple(pose_interactions[kw][i]
                          for kw in self.interaction_data_kws))

        return list(interactions)

    def _add_unique_interactions(self, interactions_list):
        """takes list of interaction tuple lists. Examines
        self.unique_interactions, add interactions if not already inserted.

        self.unique_interactions {(interaction tuple): unique_interaction_idx}

        Args:
            interactions_list (List): List of interaction tuples
            for insertion into database
        """

        for pose in interactions_list:
            for interaction_tuple in pose:
                if interaction_tuple not in self.unique_interactions:
                    self.unique_interactions[
                        interaction_tuple] = self.next_unique_interaction_idx
                    if self.interactions_initialized_flag:
                        self._insert_one_interaction(interaction_tuple)
                        self._make_new_interaction_column(
                            self.next_unique_interaction_idx)
                    self.next_unique_interaction_idx += 1

    def _find_tolerated_interactions(self, ligand_dict):
        """take ligand dict and finds which poses we should save the
        interactions for as tolerated interactions for the top pose
        of the cluster. These runs are within the
        <self.interaction_tolerance_cutoff> angstroms RMSD of the top pose
        for a given cluster. All data for the cluster's top pose is saved.

        Args:
            ligand_dict (Dictionary): Dictionary of ligand data from parser

        Returns:
            List: List of run numbers of tolerated runs
        """
        tolerated_runs = []
        for i in range(len(ligand_dict["sorted_runs"])):
            if float(ligand_dict["cluster_rmsds"]
                     [i]) <= self.interaction_tolerance_cutoff:
                tolerated_runs.append(ligand_dict["sorted_runs"][i])
        return tolerated_runs

    def format_rows_from_dict(self, ligand_dict):
        """takes file dictionary from the file parser, formats into rows for
            the database insertion

        Args:
            ligand_dict (dict): Dictionary containing data from the fileparser

        Returns:
            Tuple: Tuple of lists ([result_row_1, result_row_2,...],
                            ligand_row,
                            [interaction_tuple_1, interaction_tuple_2, ...])
        """

        # initialize row holders
        result_rows = []
        interaction_dictionaries = []
        interaction_tuples = []

        # find run numbers for poses we want to save
        if not self.store_all_poses_flag:
            poses_to_save = self._find_cluster_top_pose_runs(ligand_dict)
        else:
            poses_to_save = ligand_dict["sorted_runs"]

        # find poses we want to save tolerated interactions for
        if self.interaction_tolerance_cutoff is not None:
            tolerated_interaction_runs = self._find_tolerated_interactions(
                ligand_dict)
        else:
            tolerated_interaction_runs = []

        # do the actual result formating
        for i in range(len(ligand_dict["sorted_runs"])):
            run_number = int(ligand_dict["sorted_runs"][i])
            # save everything if this is a cluster top pose
            if run_number in poses_to_save:
                # don't save interaction data from previous
                # cluster for first cluster
                if result_rows != []:
                    pose_interactions = self._generate_interaction_tuples(
                        interaction_dictionaries
                    )
                    # generate tuples across all dictionaries for last cluster
                    interaction_tuples.append(pose_interactions)
                    # update previous entry if tolerated interactions added
                    if self.interaction_tolerance_cutoff is not None and result_rows != []:
                        result_rows[-1][17] = len(pose_interactions) + int(
                            result_rows[-1]
                            [17])  # update number of interactions
                        result_rows[-1][18] = sum(
                            1 for interaction in pose_interactions
                            if interaction[0] == "H") + result_rows[-1][
                                18]  # update number of hydrogen bonds
                    interaction_dictionaries = [
                    ]  # clear the list for the new cluster
                result_rows.append(
                    self._generate_results_row(ligand_dict, i, run_number))
                interaction_dictionaries.append(ligand_dict["interactions"][i])
            elif run_number in tolerated_interaction_runs:
                # adds to list started by best-scoring pose in cluster
                interaction_dictionaries.append(
                    ligand_dict["interactions"]
                    [i])

        pose_interactions = self._generate_interaction_tuples(
            interaction_dictionaries
        )  # will generate tuples across all dictionaries for last cluster
        interaction_tuples.append(pose_interactions)
        # only update if we added interactions
        if self.interaction_tolerance_cutoff is not None:
            result_rows[-1][17] = len(pose_interactions) + int(
                result_rows[-1][17])  # update number of interactions
            result_rows[-1][18] = sum(
                1 for interaction in pose_interactions
                if interaction[0] == "H") + int(
                    result_rows[-1]
                    [18])  # count and update number of hydrogen bonds

        return (result_rows, self._generate_ligand_row(ligand_dict),
                interaction_tuples)

    def filter_results(self, results_filters_list, ligand_filters_list,
                       output_fields):
        """Generate and execute database queries from given filters.

        Args:
            results_filters_list (list): list of tuples with first element
                indicating column to filter and second element
                indicating passing value
            ligand_filters_list (TYPE): list of tuples with first element
                indicating column to filter and second element
                indicating passing value
            output_fields (list): List of fields (columns) to be
                included in log

        Returns:
            SQLite Cursor: Cursor of passing results
        """
        # create view of passing results
        filter_results_str = self._generate_result_filtering_query(
            results_filters_list, ligand_filters_list, output_fields)
        print(filter_results_str)
        self._create_view(
            self.passing_results_view_name,
            filter_results_str)  # make sure we keep Pose_ID in view
        # perform filtering
        filtered_results = self._run_query(filter_results_str)
        # get number of passing ligands
        return filtered_results

    def _fetch_view_names(self):
        """Returns DB curor with the names of all view in DB
        """
        return self._run_query(self._generate_view_names_query())

    # # # # # # # # # # # # # # # # #
    # # # Child-specific methods # # #
    # # # # # # # # # # # # # # # # #

    def insert_results(self, results_array):
        """takes array of database rows to insert, adds data to results table

        Args:
            results_array (numpy array): numpy array of arrays
                containing formatted result rows

        """
        raise NotImplementedError

    def insert_ligands(self, ligand_array):
        """Takes array of ligand rows, inserts into Ligands table.

        Args:
            ligand_array (numpy array): Numpy array of arrays
                containing formatted ligand rows

        """

    def insert_interactions(self, interactions_list):
        """generic function for inserting interactions from given
            interaction list into DB

        Args:
            interactions_list (list): List of tuples for interactions
                in form
                ("type", "chain", "residue", "resid", "recname", "recid")
        """
        raise NotImplementedError

    def get_results(self):
        """Gets all fields for filtered results

        Returns:
            DB cursor: Cursor with all fields and rows in passing results view
        """
        raise NotImplementedError

    def clone(self):
        """Creates a copy of the db
        """
        raise NotImplementedError

    def get_number_passing_ligands(self):
        """Returns count of ligands that passed filtering criteria

        Returns:
            Int: Number of passing ligands
        """
        raise NotImplementedError

    def fetch_passing_ligand_output_info(self):
        """fetch information required by vsmanager for writing out molecules

        Returns:
            DB cursor: contains
                LigName, ligand_smile, atom_index_map, hydrogen_parents
        """
        raise NotImplementedError

    def fetch_passing_pose_coordinates(self, ligname):
        """fetch coordinates for poses passing filter for given ligand

        Args:
            ligname (string): name of ligand to return coordinates for

        Returns:
            DB cursor: contains
                ligand_coordinates, flexible_res_coordinates, flexible_residues
        """
        raise NotImplementedError

    def fetch_nonpassing_pose_coordinates(self, ligname):
        """fetch coordinates for poses of ligname which did not pass the filter

        Args:
            ligname (string): name of ligand to fetch coordinates for

        Returns:
            DB cursor: contains
                ligand_coordinates, flexible_res_coordinates, flexible_residues
        """
        raise NotImplementedError

    def _create_connection(self):
        """Creates database connection to self.db_file

        Returns:
            DB connection: Connection object to self.db_file

        """
        raise NotImplementedError

    def _close_connection(self):
        """Closes connection to database
        """
        raise NotImplementedError

    def _close_open_cursors(self):
        """closes any cursors stored in self.open_cursors.
            Resets self.open_cursors to empty list
        """
        raise NotImplementedError

    def _initialize_db(self):
        """Create connection to db. Then, check if db needs to be written.
            If so, (if self.overwrite_flag drop existing tables and )
            initialize the tables
        """
        raise NotImplementedError

    def _drop_existing_tables(self):
        """drop any existing tables. Will only be called
            if self.overwrite_flag is true
        """
        raise NotImplementedError

    def _drop_existing_views(self):
        """drop any existing views. Will only be called
            if self.overwrite_flag is true
        """
        raise NotImplementedError

    def _run_query(self, query):
        """Executes provided SQLite query. Returns cursor for results

        Args:
            query (string): Formated SQLite query as string

        Returns:
            DB cursor: Contains results of query
        """
        raise NotImplementedError

    def _create_view(self, name, query):
        """takes name and selection query,
            creates view of query stored as name.

        Args:
            name (string): Name for view which will be created
            query (string): DB-formated query which will be used to create view
        """
        raise NotImplementedError

    def _create_results_table(self):
        """Creates table for results. Columns are:
            Pose_ID             INTEGER PRIMARY KEY AUTOINCREMENT,
            LigName             VARCHAR NOT NULL,
            ligand_smile        VARCHAR[],
            source_file         VARCHAR[],
            pose_rank           INT[],
            run_number          INT[],
            energies_binding    FLOAT(4),
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
            nr_interactions     INT[],
            num_hb              INT[],
            about_x             FLOAT(4),
            about_y             FLOAT(4),
            about_z             FLOAT(4),
            trans_x             FLOAT(4),
            trans_y             FLOAT(4),
            trans_z             FLOAT(4),
            axisangle_x         FLOAT(4),
            axisangle_y         FLOAT(4),
            axisangle_z         FLOAT(4),
            axisangle_w         FLOAT(4),
            dihedrals           VARCHAR[],
            ligand_coordinates         VARCHAR[],
            flexible_residues   VARCHAR[],
            flexible_res_coordinates   VARCHAR[]
        """
        raise NotImplementedError

    def _create_ligands_table(self):
        """Create table for ligands. Columns are:
            LigName             VARCHAR NOT NULL,
            ligand_smile        VARCHAR[],
            atom_index_map      VARCHAR[],
            hydrogen_parents    VARCHAR[],
            input_pdbqt         VARCHAR[]

        """
        raise NotImplementedError

    def _create_interaction_index_table(self):
        """create table of data for each unique interaction. Columns are:
            interaction_id      INTEGER PRIMARY KEY AUTOINCREMENT,
            interaction_type    VARCHAR[],
            rec_chain           VARCHAR[],
            rec_resname         VARCHAR[],
            rec_resid           VARCHAR[],
            rec_atom            VARCHAR[],
            rec_atomid          VARCHAR[]

        """
        raise NotImplementedError

    def _create_interaction_bv_table(self):
        """Create table of interaction bits for each pose. Columns are:
            Pose_ID INTERGER PRIMARY KEY AUTOINCREMENT
            Interaction_1
            Interaction_2
            ...
            Interaction_n

        """
        raise NotImplementedError

    def _insert_unique_interactions(self, unique_interactions):
        """Inserts interaction data for unique interactions
            into Interaction_index table

        Args:
            unique_interactions (list): List of tuples of
                interactions to be inserted
        """
        raise NotImplementedError

    def _insert_one_interaction(self, interaction):
        """Insert interaction data for a single new interaction
            into the interaction indices table

        Args:
            interaction (tuple): Tuple of interaction data
            (interaction_type, rec_chain, rec_resname,
            rec_resid, rec_atom, rec_atomid)

        """
        raise NotImplementedError

    def _make_new_interaction_column(self, column_number):
        """Add column for new interaction to interaction bitvector table

        Args:
            column_number (int): Index for new interaction
        """
        raise NotImplementedError

    def _fetch_all_plot_data(self):
        """Fetches cursor for best energies and leff for all ligands

        Returns:
            DB Cursor: Cursor containing energies_binding,
                leff for the first pose for each ligand
        """
        raise NotImplementedError

    def _generate_plot_all_results_query(self):
        """Make DB-formatted query string to get energies_binding,
            leff of first pose for each ligand

        Returns:
            String: DB-formatted query string
        """
        raise NotImplementedError

    def _fetch_passing_plot_data(self):
        """Fetches cursor for best energies and leffs for ligands
            passing filtering

        Returns:
            SQLite cursor: Cursor containing energies_binding,
                leff for the first pose for passing ligands
        """
        raise NotImplementedError

    def _generate_plot_passing_results_query(self):
        """Make DB-formatted query string to get energies_binding,
            leff of first pose for passing ligands

        Returns:
            String: DB-formatted query string
        """
        raise NotImplementedError

    def _generate_result_filtering_query(self, results_filters_list,
                                         ligand_filters_list, output_fields):
        """takes lists of filters, writes sql filtering string

        Args:
            results_filters_list (list): list of tuples where
                (filter column/key, filtering cutoff)
            ligand_filters_list (list): list of filters on ligand information
            output_fields (list): List of result column data to for output

        Returns:
            String: DB-formatted string for filtering query

        Raises:
            KeyError: Raises KeyError if user requests result ordering by
                invalid or multiple options
        """
        raise NotImplementedError

    def _generate_interaction_index_filtering_query(self, interaction_list):
        """takes list of interaction info for a given ligand, looks up
            corresponding interaction index

        Args:
            interaction_list (List): List containing interaction info
            in format
            [<interaction_type>, <rec_chain>, <rec_resname>,
            <rec_resid>, <rec_atom>]

        Returns:
            String: DB-formated query on Interaction_indices table
        """
        raise NotImplementedError

    def _generate_interaction_filtering_query(self, interaction_index_list):
        """takes list of interaction indices and searches for ligand ids which
            have those interactions

        Args:
            interaction_index_list (list): List of interaction indices

        Returns:
            String: DB-formatted query
        """
        raise NotImplementedError

    def _generate_ligand_filtering_query(self, ligand_filters):
        """write string to select from ligand table

        Args:
            ligand_filters (list): List of filters on ligand table

        Returns:
            String: DB-formatted query
        """
        raise NotImplementedError

    def _generate_interaction_bitvectors(self, interactions_list):
        """takes string of interactions and makes bitvector

        Args:
            interactions_list (list): list of list of tuples. Inner lists
            contain interaction tuples for the saved poses for a single ligand

        Returns:
            List: List of bitvectors for saved poses
        """
        raise NotImplementedError

    def _insert_interaction_bitvectors(self, bitvectors):
        """Insert bitvectors of interaction data into database

        Args:
            bitvectors (List): List of lists With inner list representing
                interaction bitvector for a pose
        """
        raise NotImplementedError

    def _generate_percentile_rank_window(self):
        """makes window with percentile ranks for percentile filtering

        Returns:
            String: DB-formatted string for creating percent ranks on energies
                binding and leff
        """
        raise NotImplementedError

    def _fetch_results_column_names(self):
        """Fetches list of string for column names in results table

        Returns:
            List: List of strings of results table column names
        """
        raise NotImplementedError

    def _delete_from_results(self):
        """Remove rows from results table if they did not pass filtering
        """
        raise NotImplementedError

    def _delete_from_ligands(self):
        """Remove rows from ligands table if they did not pass filtering
        """
        raise NotImplementedError

    def _delete_from_interactions(self):
        """Remove rows from interactions bitvector table
            if they did not pass filtering
        """
        raise NotImplementedError

    def _generate_view_names_query(self):
        """Generate string to return names of views in database
        """
        raise NotImplementedError

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


class DBManagerSQLite(DBManager):
    """SQLite-specific DBManager subclass

    Attributes:
        conn (SQLite connection): Connection to database
        energy_filter_sqlite_call_dict (dictionary): Dictionary for
            translating filter options

    """

    def __init__(self, opts={}):
        """Initialize superclass and subclass-specific instance variables

        Args:
            opts (dict, optional): Dictionary of database options
        """
        super().__init__(opts)

        self.energy_filter_sqlite_call_dict = {
            "eworst": "energies_binding < {value}",
            "ebest": "energies_binding > {value}",
            "leworst": "leff < {value}",
            "lebest": "leff > {value}",
            "epercentile": "energy_percentile_rank < {value}",
            "leffpercentile": "leff_percentile_rank < {value}"
        }

    # # # # # # # # # # # # # # # # #
    # # # # #Public methods # # # # #
    # # # # # # # # # # # # # # # # #
    def insert_results(self, results_array):
        """takes array of database rows to insert, adds data to results table

        Args:
            results_array (numpy array): numpy array of arrays containing
                formatted result rows

        """

        sql_insert = """INSERT INTO Results (
        LigName,
        ligand_smile,
        source_file,
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
        dihedrals,
        ligand_coordinates,
        flexible_residues,
        flexible_res_coordinates
        ) VALUES \
        (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)"""

        try:
            cur = self.conn.cursor()
            cur.executemany(sql_insert, results_array)
            self.conn.commit()
            cur.close()

        except Exception as e:
            print(e)
            raise e

    def insert_ligands(self, ligand_array):
        """Takes array of ligand rows, inserts into Ligands table.

        Args:
            ligand_array (numpy array): Numpy array of arrays
                containing formatted ligand rows

        """
        sql_insert = '''INSERT INTO Ligands (
        LigName,
        ligand_smile,
        atom_index_map,
        hydrogen_parents,
        input_pdbqt
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
        """Takes list of interactions, inserts into database

        Args:
            interactions_list (list): List of tuples for interactions in form
            ("type", "chain", "residue", "resid", "recname", "recid")
        """
        self._add_unique_interactions(interactions_list)

        # check if we need to initialize the interaction bv table and
        # insert first set of interaction
        if not self.interactions_initialized_flag:
            self._create_interaction_bv_table()
            self._insert_unique_interactions(
                list(self.unique_interactions.keys()))
            self.interactions_initialized_flag = True

        self._insert_interaction_bitvectors(
            self._generate_interaction_bitvectors(interactions_list))

    def clone(self):
        """Creates a copy of the db
        """
        bck = sqlite3.connect(self.db_file + ".bk")
        with bck:
            self.conn.backup(bck, pages=1)
        bck.close()

    def get_results(self):
        """Gets all fields for filtered results

        Returns:
            SQLite cursor: Cursor with all fields
                and rows in passing results view
        """
        # check if we have previously filtered and saved view
        return self._run_query("SELECT * FROM {passing_view}".format(
            passing_view=self.passing_results_view_name))

    def get_number_passing_ligands(self):
        """Returns count of the number of ligands that
            passed filtering criteria

        Returns:
            Int: Number of passing ligands
        """
        cur = self.conn.cursor()
        cur.execute(
            "SELECT COUNT(DISTINCT LigName) FROM {results_view}".format(
                results_view=self.passing_results_view_name))
        n_ligands = int(cur.fetchone()[0])
        cur.close()
        return n_ligands

    def fetch_passing_ligand_output_info(self):
        """fetch information required by vsmanager for writing out molecules

        Returns:
            SQLite cursor: contains LigName, ligand_smile,
                atom_index_map, hydrogen_parents
        """
        query = "SELECT LigName, ligand_smile, atom_index_map, hydrogen_parents FROM Ligands WHERE LigName IN (SELECT DISTINCT LigName FROM {results_view})".format(
            results_view=self.passing_results_view_name)
        return self._run_query(query)

    def fetch_passing_pose_coordinates(self, ligname):
        """fetch coordinates for poses passing filter for given ligand

        Args:
            ligname (string): name of ligand to fetch coordinates for

        Returns:
            SQLite cursor: contains ligand_coordinates,
                flexible_res_coordinates, flexible_residues
        """
        query = "SELECT ligand_coordinates, flexible_res_coordinates, flexible_residues FROM Results WHERE Pose_ID IN (SELECT Pose_ID FROM {results_view} WHERE LigName LIKE '%{ligand}%')".format(
            results_view=self.passing_results_view_name, ligand=ligname)
        return self._run_query(query)

    def fetch_nonpassing_pose_coordinates(self, ligname):
        """fetch coordinates for poses of ligname which did not pass the filter

        Args:
            ligname (string): name of ligand to fetch coordinates for

        Returns:
            SQLite cursor: contains ligand_coordinates,
                flexible_res_coordinates, flexible_residues
        """
        query = "SELECT ligand_coordinates, flexible_res_coordinates, flexible_residues FROM Results WHERE LigName LIKE '%{ligand}%' AND Pose_ID NOT IN (SELECT Pose_ID FROM {results_view})".format(
            ligand=ligname, results_view=self.passing_results_view_name)
        return self._run_query(query)

    # # # # # # # # # # # # # # # # #
    # # # # #Private methods # # # # #
    # # # # # # # # # # # # # # # # #

    def _create_connection(self):
        """Creates database connection to self.db_file

        Returns:
            SQLite connection: Connection object to self.db_file

        """
        con = None
        try:
            con = sqlite3.connect(self.db_file)
            cursor = con.cursor()
            cursor.execute('PRAGMA synchronous = OFF')
            cursor.execute('PRAGMA journal_mode = MEMORY')
            cursor.close()
        except Exception as e:
            print("Error while creating database connection")
            raise e
        return con

    def _close_connection(self):
        """Closes connection to database
        """
        print("Closing database")
        self.conn.close()

    def _close_open_cursors(self):
        """closes any cursors stored in self.open_cursors.
            Resets self.open_cursors to empty list
        """
        for cur in self.open_cursors:
            cur.close()

        self.open_cursors = []

    def _initialize_db(self):
        """Create connection to db. Then, check if db needs to be written.
            If so, (if self.overwrite_flag drop existing tables and )
            initialize the tables
        """
        self.conn = self._create_connection()

        if not self.write_db_flag:
            return
        # if we want to overwrite old db, drop existing tables
        if self.overwrite_flag:
            self._drop_existing_tables()
        # create tables in db
        self._create_results_table()
        self._create_ligands_table()
        self._create_interaction_index_table()

    def _drop_existing_tables(self):
        """drop any existing tables.
        Will only be called if self.overwrite_flag is true
        """

        # fetch existing tables
        cur = self.conn.cursor()
        cur.execute("SELECT name FROM sqlite_schema WHERE type='table';")
        tables = cur.fetchall()

        # drop tables
        for table in tables:
            # cannot drop this, so we catch it instead
            if table[0] == "sqlite_sequence":
                continue
            cur.execute("DROP TABLE {table_name}".format(table_name=table[0]))
        cur.close()

    def _drop_existing_views(self):
        """Drop any existing views
        Will only be called if self.overwrite_flag is true
        """
        # fetch existing views
        cur = self.conn.cursor()
        cur.execute("SELECT name FROM sqlite_schema WHERE type='view';")
        views = cur.fetchall()

        # drop views
        for view in views:
            # cannot drop this, so we catch it instead
            if view[0] == "sqlite_sequence":
                continue
            cur.execute("DROP VIEW {view_name}".format(view_name=view[0]))
        cur.close()

    def _run_query(self, query):
        """Executes provided SQLite query. Returns cursor for results.
            Since cursor remains open, added to list of open cursors

        Args:
            query (string): Formated SQLite query as string

        Returns:
            SQLite cursor: Contains results of query
        """
        cur = self.conn.cursor()
        cur.execute(query)
        self.open_cursors.append(cur)
        return cur

    def _create_view(self, name, query):
        """takes name and selection query,
            creates view of query stored as name.

        Args:
            name (string): Name for view which will be created
            query (string): SQLite-formated query used to create view
        """
        cur = self.conn.cursor()
        query = query.replace("SELECT ", "SELECT Pose_ID, ", 1)
        # drop old view if there is one
        cur.execute("DROP VIEW IF EXISTS {name}".format(name=name))
        cur.execute("CREATE VIEW {name} AS {query}".format(name=name,
                                                           query=query))
        cur.close()

    def _create_results_table(self):
        """Creates table for results. Columns are:
            Pose_ID             INTEGER PRIMARY KEY AUTOINCREMENT,
            LigName             VARCHAR NOT NULL,
            ligand_smile        VARCHAR[],
            source_file         VARCHAR[],
            pose_rank           INT[],
            run_number          INT[],
            energies_binding    FLOAT(4),
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
            nr_interactions     INT[],
            num_hb              INT[],
            about_x             FLOAT(4),
            about_y             FLOAT(4),
            about_z             FLOAT(4),
            trans_x             FLOAT(4),
            trans_y             FLOAT(4),
            trans_z             FLOAT(4),
            axisangle_x         FLOAT(4),
            axisangle_y         FLOAT(4),
            axisangle_z         FLOAT(4),
            axisangle_w         FLOAT(4),
            dihedrals           VARCHAR[],
            ligand_coordinates         VARCHAR[],
            flexible_residues   VARCHAR[],
            flexible_res_coordinates   VARCHAR[]
        """
        sql_results_table = """CREATE TABLE Results (
            Pose_ID             INTEGER PRIMARY KEY AUTOINCREMENT,
            LigName             VARCHAR NOT NULL,
            ligand_smile        VARCHAR[],
            source_file         VARCHAR[],
            pose_rank           INT[],
            run_number          INT[],
            energies_binding    FLOAT(4),
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
            nr_interactions     INT[],
            num_hb              INT[],
            about_x             FLOAT(4),
            about_y             FLOAT(4),
            about_z             FLOAT(4),
            trans_x             FLOAT(4),
            trans_y             FLOAT(4),
            trans_z             FLOAT(4),
            axisangle_x         FLOAT(4),
            axisangle_y         FLOAT(4),
            axisangle_z         FLOAT(4),
            axisangle_w         FLOAT(4),
            dihedrals           VARCHAR[],
            ligand_coordinates         VARCHAR[],
            flexible_residues   VARCHAR[],
            flexible_res_coordinates   VARCHAR[]
        );
        """

        try:
            cur = self.conn.cursor()
            cur.execute(sql_results_table)
            cur.close()
        except Exception as e:
            print(
                "Error while creating results table. If database already exists, use --overwrite to drop existing tables"
            )
            raise e

    def _create_ligands_table(self):
        """Create table for ligands. Columns are:
            LigName             VARCHAR NOT NULL,
            ligand_smile        VARCHAR[],
            atom_index_map      VARCHAR[],
            hydrogen_parents    VARCHAR[],
            input_pdbqt         VARCHAR[]

        """
        ligand_table = """CREATE TABLE Ligands (
            LigName             VARCHAR NOT NULL,
            ligand_smile        VARCHAR[],
            atom_index_map      VARCHAR[],
            hydrogen_parents    VARCHAR[],
            input_pdbqt         VARCHAR[])"""

        try:
            cur = self.conn.cursor()
            cur.execute(ligand_table)
            cur.close()
        except Exception as e:
            print(
                "Error while creating ligands table. If database already exists, use --overwrite to drop existing tables"
            )
            raise e

    def _create_interaction_index_table(self):
        """create table of data for each unique interaction. Columns are:
            interaction_id      INTEGER PRIMARY KEY AUTOINCREMENT,
            interaction_type    VARCHAR[],
            rec_chain           VARCHAR[],
            rec_resname         VARCHAR[],
            rec_resid           VARCHAR[],
            rec_atom            VARCHAR[],
            rec_atomid          VARCHAR[]

        """
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
            print(
                "Error while creating interaction index table. If database already exists, use --overwrite to drop existing tables"
            )
            raise e

    def _create_interaction_bv_table(self):
        """Create table of interaction bits for each pose. Columns are:
            Pose_ID INTERGER PRIMARY KEY AUTOINCREMENT
            Interaction_1
            Interaction_2
            ...
            Interaction_n

        """
        interact_columns_str = " INTEGER,\n".join([
            "Interaction_" + str(i + 1)
            for i in range(len(self.unique_interactions))
        ]) + " INTEGER"

        bv_table = """CREATE TABLE Interaction_bitvectors (
        Pose_ID INTEGER PRIMARY KEY AUTOINCREMENT,
        {columns})""".format(columns=interact_columns_str)

        try:
            cur = self.conn.cursor()
            cur.execute(bv_table)
            cur.close()
        except Exception as e:
            print(
                "Error while creating interaction bitvector table. If database already exists, use --overwrite to drop existing tables"
            )
            raise e

    def _insert_unique_interactions(self, unique_interactions):
        """Inserts interaction data for unique interactions
            into Interaction_index table

        Args:
            unique_interactions (list): List of tuples of interactions
                to be inserted
        """
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
            print(
                "Error while inserting unique interactions into interaction index table"
            )
            raise e

    def _insert_one_interaction(self, interaction):
        """Insert interaction data for a single new interaction
            into the interaction indices table

        Args:
            interaction (tuple): Tuple of interaction data
                (interaction_type, rec_chain, rec_resname,
                rec_resid, rec_atom, rec_atomid)

        """
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
            print(
                "Error inserting interaction {interact} into interaction index table"
                .format(interact=str(interaction)))
            raise e

    def _make_new_interaction_column(self, column_number):
        """Add column for new interaction to interaction bitvector table

        Args:
            column_number (int): Index for new interaction
        """
        add_column_str = '''ALTER TABLE Interaction_bitvectors ADD COLUMN Interaction_{n_inter}'''.format(
            n_inter=str(column_number))
        try:
            cur = self.conn.cursor()
            cur.execute(add_column_str)
            self.conn.commit()
            cur.close()

        except Exception as e:
            print(
                "Error adding column for Interaction_{num} to interaction bitvector table"
                .format(num=str(column_number)))
            raise e

    def _fetch_all_plot_data(self):
        """Fetches cursor for best energies and leff for all ligands

        Returns:
            SQLite Cursor: Cursor containing energies_binding,
                leff for the first pose for each ligand
        """
        return self._run_query(self._generate_plot_all_results_query())

    def _generate_plot_all_results_query(self):
        """Make SQLite-formatted query string to get energies_binding,
            leff of first pose of all ligands

        Returns:
            String: SQLite-formatted query string
        """
        return "SELECT energies_binding, leff FROM Results GROUP BY LigName"

    def _fetch_passing_plot_data(self):
        """Fetches cursor for best energies and leffs for
            ligands passing filtering

        Returns:
            SQLite cursor: Cursor containing energies_binding,
                leff for the first pose for passing ligands
        """
        return self._run_query(self._generate_plot_passing_results_query())

    def _generate_plot_passing_results_query(self):
        """Make SQLite-formatted query string to get energies_binding,
            leff of first pose for passing ligands

        Returns:
            String: SQLite-formatted query string
        """
        return "SELECT energies_binding, leff FROM Results WHERE LigName IN (SELECT DISTINCT LigName FROM {results_view}) GROUP BY LigName".format(
            results_view=self.passing_results_view_name)

    def _generate_result_filtering_query(self, results_filters_list,
                                         ligand_filters_list, output_fields):
        """takes lists of filters, writes sql filtering string

        Args:
            results_filters_list (list): list of tuples where
                (filter column/key, filtering cutoff)
            ligand_filters_list (list): list of filters on ligand information
            output_fields (list): List of result column data for output

        Returns:
            String: SQLite-formatted string for filtering query

        Raises:
            KeyError: Raises KeyError if user requests
                result ordering by invalid or multiple options
        """

        # parse requested output fields and convert to column names in database

        outfield_string = "LigName, " + ", ".join(
            [self.field_to_column_name[field] for field in output_fields])
        filtering_window = "Results"

        # write energy filters and compile list of interactions to search for
        energy_filter_sql_query = []
        interaction_filters = []

        for filter_key, filter_value in results_filters_list:
            if filter_key in self.energy_filter_sqlite_call_dict:
                if filter_key == "epercentile" or filter_key == "leffpercentile":
                    # convert from percent to decimal
                    filter_value = str(float(filter_value) / 100)
                    # reset filtering window to include percentile_ranks
                    filtering_window = "({percentile_window})".format(
                        percentile_window=self.
                        _generate_percentile_rank_window())
                energy_filter_sql_query.append(
                    self.energy_filter_sqlite_call_dict[filter_key].format(
                        value=filter_value))

            # write hb count filter(s)
            if filter_key == 'hb_count':
                if filter_value > 0:
                    energy_filter_sql_query.append(
                        "num_hb > {value}".format(value=filter_value))
                else:
                    energy_filter_sql_query.append(
                        "num_hb < {value}".format(value=-1 * filter_value))

            # add interaction filters to list
            if filter_key in self.interaction_filter_types:
                for interact in filter_value:
                    interaction_string = filter_key + ":" + interact[0]
                    interaction_filters.append(interaction_string.split(":"))

            # add react_any flag as interaction filter
            # check if react_any is true
            if filter_key == "react_any" and filter_value:
                interaction_filters.append(["R", "", "", "", ""])

        # initialize query string
        sql_string = """SELECT {out_columns} FROM {window} WHERE """.format(
            out_columns=outfield_string, window=filtering_window)

        # add energy filters to our query string
        sql_string += " AND ".join(energy_filter_sql_query)

        # for each interaction filter, get the index
        # from the interactions_indices table
        for interaction in interaction_filters:
            interaction_filter_indices = []
            interact_index_str = self._generate_interaction_index_filtering_query(
                interaction)
            interaction_indices = self._run_query(interact_index_str)
            for i in interaction_indices:
                interaction_filter_indices.append(i[0])

            # catch if interaction not found in results
            if interaction_filter_indices == []:
                print(
                    "Interaction {i} not found in results, excluded from filtering"
                    .format(i=":".join(interaction)))
                continue
            # find pose ids for ligands with desired interactions
            sql_string += " AND Pose_ID IN ({interaction_str})".format(
                interaction_str=self._generate_interaction_filtering_query(
                    interaction_filter_indices))

        # add ligand filters
        if ligand_filters_list != []:
            sql_string += " AND LigName IN ({ligand_str})".format(
                ligand_str=self._generate_ligand_filtering_query(
                    ligand_filters_list))

        # adding if we only want to keep
        # one pose per ligand (will keep first entry)
        if self.log_distinct_ligands:
            sql_string += " GROUP BY LigName"

        # add how to order results
        if self.order_results is not None:
            try:
                sql_string += " ORDER BY " + self.field_to_column_name[
                    self.order_results]
            except KeyError:
                print(
                    "Please ensure you are only requesting one option for --order_results and have written it correctly"
                )
                raise KeyError

        return sql_string

    def _generate_interaction_index_filtering_query(self, interaction_list):
        """takes list of interaction info for a given ligand,
            looks up corresponding interaction index

        Args:
            interaction_list (List): List containing interaction info
                in format [<interaction_type>, <rec_chain>, <rec_resname>,
                <rec_resid>, <rec_atom>]

        Returns:
            String: SQLite-formated query on Interaction_indices table
        """
        interaction_info = [
            "interaction_type", "rec_chain", "rec_resname", "rec_resid",
            "rec_atom"
        ]
        len_interaction_info = len(interaction_info)
        sql_string = "SELECT interaction_id FROM Interaction_indices WHERE "

        sql_string += " AND ".join([
            "{column} LIKE '%{value}%'".format(column=interaction_info[i],
                                               value=interaction_list[i])
            for i in range(len_interaction_info) if interaction_list[i] != ""
        ])

        return sql_string

    def _generate_interaction_filtering_query(self, interaction_index_list):
        """takes list of interaction indices and searches for ligand ids
            which have those interactions

        Args:
            interaction_index_list (list): List of interaction indices

        Returns:
            String: SQLite-formatted query
        """

        return "SELECT Pose_id FROM Interaction_bitvectors WHERE " + " OR ".join(
            [
                "Interaction_{index_n} = 1".format(index_n=index)
                for index in interaction_index_list
            ])

    def _generate_ligand_filtering_query(self, ligand_filters):
        """write string to select from ligand table

        Args:
            ligand_filters (list): List of filters on ligand table

        Returns:
            String: SQLite-formatted query
        """

        sql_ligand_string = "SELECT LigName FROM Ligands WHERE "

        substruct_flag = ligand_filters['F'][0].upper()
        for kw in ligand_filters.keys():
            fils = ligand_filters[kw]
            if kw == 'N':
                for name in fils:
                    name_sql_str = "LigName LIKE '%{value}%' OR ".format(
                        value=name)
                    sql_ligand_string += name_sql_str
            if kw == "S":
                for substruct in fils:
                    substruct_sql_str = "ligand_smile LIKE '%{value}%' {flag} ".format(
                        value=substruct, flag=substruct_flag)
                    sql_ligand_string += substruct_sql_str

        if sql_ligand_string.endswith("AND "):
            sql_ligand_string = sql_ligand_string.rstrip("AND ")
        if sql_ligand_string.endswith("OR "):
            sql_ligand_string = sql_ligand_string.rstrip("OR ")

        return sql_ligand_string

    def _generate_interaction_bitvectors(self, interactions_list):
        """takes string of interactions and makes bitvector

        Args:
            interactions_list (list): list of list of tuples.
                Inner lists contain interaction tuples for
                saved poses for a single ligand

        Returns:
            List: List of bitvectors for saved poses
        """
        bitvectors_list = []
        for pose_interactions in interactions_list:
            pose_bitvector = [None] * len(self.unique_interactions)
            for interaction_tuple in pose_interactions:
                interaction_idx = self.unique_interactions[interaction_tuple]
                # index corrected for interaction_indices starting at 1
                pose_bitvector[interaction_idx - 1] = 1

            bitvectors_list.append(pose_bitvector)

        return bitvectors_list

    def _insert_interaction_bitvectors(self, bitvectors):
        """Insert bitvectors of interaction data into database

        Args:
            bitvectors (List): List of lists With inner list representing
                interaction bitvector for a pose
        """

        interaction_columns = range(len(self.unique_interactions))
        column_str = ""
        filler_str = ""
        for i in interaction_columns:
            column_str += "Interaction_" + str(i + 1) + ", "
            filler_str += "?,"
        column_str = column_str.rstrip(", ")
        filler_str = filler_str.rstrip(",")

        sql_insert = '''INSERT INTO Interaction_bitvectors ({columns}) VALUES ({fillers})'''.format(
            columns=column_str, fillers=filler_str)

        try:
            cur = self.conn.cursor()
            cur.executemany(sql_insert, bitvectors)
            self.conn.commit()
            cur.close()

        except Exception as e:
            print("Error while inserting bitvectors")
            raise e

    def _generate_percentile_rank_window(self):
        """makes window with percentile ranks for percentile filtering

        Returns:
            String: SQLite-formatted string for creating
                percent ranks on energies_binding and leff
        """
        column_names = ",".join(self._fetch_results_column_names())
        return "SELECT {columns}, PERCENT_RANK() OVER (ORDER BY energies_binding) energy_percentile_rank, PERCENT_RANK() OVER (ORDER BY leff) leff_percentile_rank FROM Results".format(
            columns=column_names)

    def _fetch_results_column_names(self):
        """Fetches list of string for column names in results table

        Returns:
            List: List of strings of results table column names
        """
        return [
            column_tuple[1]
            for column_tuple in self.conn.execute("PRAGMA table_info(Results)")
        ]

    def _delete_from_results(self):
        """Remove rows from results table if they did not pass filtering
        """
        cur = self.conn.cursor()
        cur.execute("DELETE FROM Results WHERE Pose_ID NOT IN {view}".format(
            view=self.passing_results_view_name))
        self.conn.commit()
        cur.close()

    def _delete_from_ligands(self):
        """Remove rows from ligands table if they did not pass filtering
        """
        cur = self.conn.cursor()
        cur.execute("DELETE FROM Ligands WHERE LigName NOT IN {view}".format(
            view=self.passing_results_view_name))
        self.conn.commit()
        cur.close()

    def _delete_from_interactions(self):
        """Remove rows from interactions bitvector table
            if they did not pass filtering
        """
        cur = self.conn.cursor()
        cur.execute(
            "DELETE FROM Interaction_bitvectors WHERE Pose_ID NOT IN {view}".
            format(view=self.passing_results_view_name))
        self.conn.commit()
        cur.close()

    def _generate_view_names_query(self):
        """Generate string to return names of views in database
        """
        return "SELECT name FROM sqlite_schema WHERE type = 'view'"

