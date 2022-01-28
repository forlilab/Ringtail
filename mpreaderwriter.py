import multiprocessing
from time import sleep
import parsers
import numpy as np

class DockingFileReader(multiprocessing.Process):
    """ this class is the individual worker for processing dlgs"""
    def __init__(self, queueIn, queueOut, mode, numclusters):
        #set mode for which file parser to use
        self.mode = mode
        #set number of clusters to write
        self.num_clusters = numclusters
        # initialize the parent class to inherit all multiprocessing methods
        multiprocessing.Process.__init__(self)
        # each worker knows about the queue in (where data to process comes from)...
        self.queueIn = queueIn
        # ...and a queue out (where to send the results)
        self.queueOut = queueOut

    def get_best_cluster_poses(self, ligand_dict):
        """takes input ligand dictionary, reads run pose clusters, adds "cluster_best_run" entry with the top scoring run for each cluster"""
        top_poses = []
        cluster_dict = ligand_dict["clusters"]
        for cluster in cluster_dict:
            top_poses.append(cluster_dict[cluster][0])
        ligand_dict["cluster_top_poses"] = top_poses
        return ligand_dict

    def write_results_interaction_rows(self, ligand_dict):
        """writes list of lists of ligand values to be inserted into sqlite database"""

        ligand_rows = []
        interaction_rows = []

        #initialize list for sql row with name and smile string
        ligand_name = ligand_dict["ligname"]
        ligand_smile = ligand_dict["ligand_smile_string"]

        ligand_data_keys = ["cluster_rmsds",
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

        ligand_interaction_keys = ["type",
        "chain",
        "residue",
        "resid",
        "recname",
        "recid"]

        stateVar_keys = ["pose_about",
        "pose_translations",
        "pose_quarternions"]
        
        ######get pose-specific data
        for i in range(len(ligand_dict["sorted_runs"])):
            #check if run is best for a cluster. We are only saving the top pose for each cluster
            pose_rank = i
            run_number = ligand_dict["sorted_runs"][pose_rank]
            try:
                cluster_top_pose_runs = ligand_dict["cluster_top_poses"][:self.num_clusters] #will only select top n clusters. Default 3
            except IndexError:
                cluster_top_pose_runs = ligand_dict["cluster_top_poses"] #catch indexerror if not enough clusters for given ligand
            if run_number in cluster_top_pose_runs:
                ligand_data_list = [ligand_name, ligand_smile, pose_rank+1, run_number]
                #get energy data
                for key in ligand_data_keys:
                    ligand_data_list.append(ligand_dict[key][pose_rank])

                pose_interactions_dict = ligand_dict["interactions"][pose_rank]
                num_interactions = pose_interactions_dict["count"][0]
                ligand_data_list.append(num_interactions) 
                interaction_strings_list = []
                hb_flag = False
                for key in ligand_interaction_keys:
                    interaction_data = pose_interactions_dict[key]
                    if type(interaction_data) == list:
                        interaction_data_string = ""
                        for interaction in interaction_data:
                            interaction_data_string = interaction_data_string + interaction + ", "
                    interaction_strings_list.append(interaction_data_string)
                    if key == 'type' and not hb_flag: #add hb_count only once
                        hb_count = interaction_data_string.count("H")
                        ligand_data_list.append(hb_count)
                        hb_flag = True
                hb_flag = False
                #reformat interaction data to have all information for each interaction together
                interactions_string = ""
                for i in range(int(num_interactions)):
                    single_interaction_string = ""
                    for line in interaction_strings_list:
                        line_list = line.split(",")
                        single_interaction_string += line_list[i] + ":"
                    interactions_string += single_interaction_string.replace(" ", "") + ", "
                interaction_rows.append(interactions_string)


                for key in stateVar_keys:
                    stateVar_data = ligand_dict[key][pose_rank]
                    for dim in stateVar_data:
                        ligand_data_list.append(dim)
                pose_dihedrals = ligand_dict["pose_dihedrals"][pose_rank]
                dihedral_string = ""
                for dihedral in pose_dihedrals:
                    dihedral_string = dihedral_string + str(dihedral) + ", "
                ligand_data_list.append(dihedral_string)

                ligand_rows.append(ligand_data_list)

        return ligand_rows, interaction_rows

    def write_ligand_row(self, ligand_dict):
        """writes row to be inserted into ligand table"""
        ligand_name = ligand_dict["ligname"]
        ligand_smile = ligand_dict["ligand_smile_string"]
        input_pdbqt = "\n".join(ligand_dict["ligand_input_pdbqt"])
        best_binding = ligand_dict["scores"][0]
        best_run = ligand_dict["sorted_runs"][0]

        return [ligand_name, ligand_smile, input_pdbqt, best_binding, best_run]

    def run(self):
        # method overload from parent class
        #
        # this is where the task of this class is performed
        #
        # each multiprocessing.Process class must have a "run" method which
        # is called by the initialization (see below) with start()
        #
        proc_name = multiprocessing.current_process().name
        while True:
            # retrieve from the queue in the next task to be done
            next_task = self.queueIn.get()
            # if a poison pill is received, this worker's job is done, quit
            if next_task is None:
                # the poison pill can be anything
                print('%s: Exiting' % proc_name)
                # before leaving, pass the poison pill back in the queue (for the writer, see below)
                self.queueOut.put(None)
                break
            print('%s: %s' % (proc_name, next_task))
            # generate CPU LOAD
            if self.mode == "dlg":
                parsed_file_dict = parsers.parse_single_dlg(next_task)
                parsed_file_dict = self.get_best_cluster_poses(parsed_file_dict)
                resultsAndInteractions = self.write_results_interaction_rows(parsed_file_dict)
                results_rows = resultsAndInteractions[0]
                interaction_rows = resultsAndInteractions[1]
                ligand_row = self.write_ligand_row(parsed_file_dict)
                file_packet = (results_rows, ligand_row, interaction_rows)
            # put the result in the out queue
            self.queueOut.put(file_packet)
        return

class Writer(multiprocessing.Process):
    # this class is a listener that retrieves data from the queue and writes it
    # into datbase
    def __init__(self, queue, maxProcesses, chunksize, db_obj):
        multiprocessing.Process.__init__(self)
        self.queue = queue
        # this class knows about how many multi-processing workers there are
        self.maxProcesses=maxProcesses
        # assign pointer to db object, set chunksize
        self.db = db_obj
        self.chunksize = chunksize
        # initialize data arrays
        self.results_array = []
        self.ligands_array = []
        self.interaction_rows_list = []
        self.bitvectors_list = []
        self.counter = 0
        self.pose_id_list = []

        self.current_pose_id = 1
        self.unique_interactions = set()
        self.unique_interactions_list = []
        self.unique_interactions_split = []

    def run(self):
        # method overload from parent class
        #
        # this is where the task of this class is performed
        #
        # each multiprocessing.Process class must have a "run" method which
        # is called by the initialization (see below) with start()
        #
        while True:
            # retrieve the next task from the queue
            next_task = self.queue.get()
            if next_task == None:
                # if a poison pill is found, it means one of the workers quit
                self.maxProcesses -= 1
            else:
                # if not a poison pill, process the task item

                # after every n (chunksize) files, write to db and clear data-holders
                if self.counter >= self.chunksize:
                    self.write2db()

                #process next file    
                self.process_file(next_task)

            if self.maxProcesses == 0:
                # received as many poison pills as workers
                #perform final db write
                self.write2db()
                # no workers left, no job to do
                break
        return

    def get_unique_interactions(self):
        self.new_interactions = []
        for pose in self.interaction_rows_list:
            interaction_list = pose.replace(":,",",").split(",")
            for interaction in interaction_list:
                if interaction.endswith(":"):
                        interaction = interaction.rstrip(":") #remove trailing colon
                if interaction not in self.unique_interactions and interaction != "":
                    interaction_attributes= interaction.split(":")
                    self.unique_interactions.add(interaction)
                    self.unique_interactions_list.append(interaction)
                    self.unique_interactions_split.append(interaction_attributes) #type,ligand atom, ligand id, chain, residue, resid, residue atom, residue atom id
                    self.new_interactions.append(interaction_attributes)

    def write_interaction_bitvectors(self):
        """takes string of interactions and all possible interactions and makes bitvector"""
        for i in range(len(self.pose_id_list)):
            pose_id = self.pose_id_list[i]
            interactions_string = self.interaction_rows_list[i]
            pose_bitvector = [pose_id]
            for interaction in self.unique_interactions_list:
                if interaction in interactions_string:
                    pose_bitvector.append(1) #true
                else:
                    pose_bitvector.append(None) #false

            self.bitvectors_list.append(pose_bitvector)

    def write2db(self):
        #insert result and ligand data
        self.db.insert_results(np.array(self.results_array))
        self.db.insert_ligands(np.array(self.ligands_array))

        #perform operations for inserting iteraction data
        self.new_interaction_idx_counter = len(self.unique_interactions) + 1 #save index of next unique interaction before adding new ones
        self.get_unique_interactions()
        try:
            self.unique_interactions_split.remove([' ']) #make sure we are not passing any empty lists to dbman
        except ValueError:
            pass

        if not self.db.bv_table_flag:
            self.db.insert_interactions(unique_interactions = self.unique_interactions, unique_interactions_split = self.unique_interactions_split)

        else:
            self.db.insert_interactions(new_interactions = self.new_interactions, interaction_idx_counter = self.new_interaction_idx_counter)

        self.write_interaction_bitvectors()
        self.db.insert_interaction_BVs(self.bitvectors_list, len(self.unique_interactions))

        #reset all data-holders for next chunk
        self.results_array = []
        self.ligands_array = []
        self.interaction_rows_list = []
        self.pose_id_list = []
        self.bitvectors_list = []
        self.counter = 0

    def process_file(self, next_task):
        file_packet = next_task
        results_rows = file_packet[0]
        ligand_row = file_packet[1]
        interaction_rows = file_packet[2]
        for pose in results_rows:
            self.results_array.append(pose)
        for pose in interaction_rows:
            self.interaction_rows_list.append(pose)
            self.pose_id_list.append(self.current_pose_id)
            self.current_pose_id += 1
        self.ligands_array.append(ligand_row)









