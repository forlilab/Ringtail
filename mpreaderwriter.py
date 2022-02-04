import multiprocessing
from time import sleep
import parsers
import numpy as np

class DockingFileReader(multiprocessing.Process):
    """ this class is the individual worker for processing dlgs"""
    def __init__(self, queueIn, queueOut, dbman, mode, numclusters, no_print):
        #set mode for which file parser to use
        self.mode = mode
        #set number of clusters to write
        self.num_clusters = numclusters
        #set dbmanager
        self.dbman = dbman
        #set flag for printing
        self.no_print = no_print
        # initialize the parent class to inherit all multiprocessing methods
        multiprocessing.Process.__init__(self)
        # each worker knows about the queue in (where data to process comes from)...
        self.queueIn = queueIn
        # ...and a queue out (where to send the results)
        self.queueOut = queueOut

    def find_best_cluster_poses(self, ligand_dict):
        """takes input ligand dictionary, reads run pose clusters, adds "cluster_best_run" entry with the top scoring run for each cluster"""
        top_poses = []
        cluster_dict = ligand_dict["clusters"]
        for cluster in cluster_dict:
            top_poses.append(cluster_dict[cluster][0])
        ligand_dict["cluster_top_poses"] = top_poses
        return ligand_dict

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
                if not self.no_print:
                    print('%s: Exiting' % proc_name)
                # before leaving, pass the poison pill back in the queue (for the writer, see below)
                self.queueOut.put(None)
                break
            if not self.no_print:
                print('%s: %s' % (proc_name, next_task))
            # generate CPU LOAD
            if self.mode == "dlg":
                parsed_file_dict = parsers.parse_single_dlg(next_task)
                parsed_file_dict = self.find_best_cluster_poses(parsed_file_dict)
                file_packet = self.dbman.format_rows_from_dict(parsed_file_dict)
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
        self.interactions_list = []
        self.counter = 0
        self.pose_id_list = []

        self.current_pose_id = 1

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
                    self.write_to_db()

                #process next file    
                self.process_file(next_task) #stack rows onto corresponding arrays
            if self.maxProcesses == 0:
                # received as many poison pills as workers
                #perform final db write
                self.write_to_db()
                # no workers left, no job to do
                print("Break")
                self.close()
                break

        return

    def write_to_db(self):
        #insert result and ligand data
        self.db.insert_results(filter(None, self.results_array))#filter out stray Nones
        self.db.insert_ligands(filter(None,self.ligands_array))

        #perform operations for inserting iteraction data
        self.db.insert_interactions(self.pose_id_list, self.interactions_list)

        #reset all data-holders for next chunk
        self.results_array = []
        self.ligands_array = []
        self.interactions_list = []
        self.pose_id_list = []
        self.counter = 0

    def process_file(self, file_packet):
        results_rows, ligand_row, interaction_rows = file_packet
        for pose in results_rows:
            self.results_array.append(pose)
        for pose in interaction_rows:
            self.interactions_list.append(pose)
            self.pose_id_list.append(self.current_pose_id)
            self.current_pose_id += 1
        self.ligands_array.append(ligand_row)









