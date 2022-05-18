#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail multiprocessing workers
#

import multiprocessing
import time
import sys
from .parsers import parse_single_dlg


class DockingFileReader(multiprocessing.Process):
    """ this class is the individual worker for processing dlgs"""

    def __init__(self, queueIn, queueOut, dbman, mode, numclusters, no_print):
        # set mode for which file parser to use
        self.mode = mode
        # set number of clusters to write
        self.num_clusters = numclusters
        # set dbmanager
        self.dbman = dbman
        # set flag for printing
        self.no_print = no_print
        # initialize the parent class to inherit all multiprocessing methods
        multiprocessing.Process.__init__(self)
        # each worker knows the queue in (where data to process comes from)
        self.queueIn = queueIn
        # ...and a queue out (where to send the results)
        self.queueOut = queueOut

    def find_best_cluster_poses(self, ligand_dict):
        """takes input ligand dictionary, reads run pose clusters,
            adds "cluster_best_run" entry with the top scoring
            run for each cluster"""
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
        while True:
            # retrieve from the queue in the next task to be done
            next_task = self.queueIn.get()
            # if a poison pill is received, this worker's job is done, quit
            if next_task is None:
                # before leaving, pass the poison pill back in the queue
                self.queueOut.put(None)
                break

            # generate CPU LOAD
            # parser depends on requested mode
            if self.mode == "dlg":
                parsed_file_dict = parse_single_dlg(next_task)
            # future: Vina parser, etc
            # elif self.mode == "vina":
            #    parsed_file_dict = parser.parse_single_vina_log(next_task)
            parsed_file_dict = self.find_best_cluster_poses(parsed_file_dict)
            file_packet = self.dbman.format_rows_from_dict(parsed_file_dict)
            # put the result in the out queue
            self.queueOut.put(file_packet)
        return


class Writer(multiprocessing.Process):
    # this class is a listener that retrieves data from the queue and writes it
    # into datbase
    def __init__(self, queue, maxProcesses, chunksize, db_obj, num_files, single_receptor):
        multiprocessing.Process.__init__(self)
        self.queue = queue
        # this class knows about how many multi-processing workers there are
        self.maxProcesses = maxProcesses
        # assign pointer to db object, set chunksize, single_receptor flag
        self.db = db_obj
        self.chunksize = chunksize
        self.single_receptor = single_receptor
        # initialize data arrays
        self.results_array = []
        self.ligands_array = []
        self.interactions_list = []
        self.receptor_array = []
        # progress tracking instance variables
        self.first_insert = True
        self.counter = 0
        self.num_files_remaining = num_files
        self.total_num_files = num_files
        self.num_files_written = 0
        self.time0 = time.perf_counter()
        # based on estimated 3 seconds per 1000 files, converted to minutes
        self.est_time_remaining = 3 * num_files / (1000 * 60)
        self.last_write_time = 0

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
            if next_task is None:
                # if a poison pill is found, it means one of the workers quit
                self.maxProcesses -= 1
                print("Closing process. Remaining open processes:",
                      self.maxProcesses)
            else:
                # if not a poison pill, process the task item

                # after every n (chunksize) files, write to db
                if self.counter >= self.chunksize:
                    self.write_to_db()

                # process next file
                self.process_file(next_task)
                # print info about files and time remaining
                sys.stdout.write('\r')
                sys.stdout.write(
                    "{n} files remaining to process. Estimated time remaining: {est_time:.2f} minutes. Last chunk took {write_time:.2f} seconds to write."
                    .format(n=self.num_files_remaining,
                            est_time=self.est_time_remaining,
                            write_time=self.last_write_time))
                sys.stdout.flush()
                self.num_files_remaining -= 1

            if self.maxProcesses == 0:
                # received as many poison pills as workers
                print("Performing final database write")
                # perform final db write
                self.write_to_db()
                # no workers left, no job to do
                print("File processing completed")
                self.close()
                break

        return

    def write_to_db(self):
        time1 = time.perf_counter()
        # insert result, ligand, and receptor data
        self.db.insert_results(filter(
            None, self.results_array))  # filter out stray Nones
        self.db.insert_ligands(filter(None, self.ligands_array))
        # if this is the first insert or we have multiple receptors, insert the receptor array
        if self.first_insert or not self.single_receptor:
            self.db.insert_receptors(filter(None, self.receptor_array))
            self.first_insert = False

        # perform operations for inserting iteraction data
        self.db.insert_interactions(self.interactions_list)

        # calulate time for processing/writing previous chunk
        self.num_files_written += self.chunksize
        total_runtime = time.perf_counter() - self.time0
        self.est_time_remaining = (total_runtime / self.num_files_written) * (
            self.total_num_files - self.num_files_written) / 60  # converted to minutes
        self.last_write_time = time.perf_counter() - time1

        # reset all data-holders for next chunk
        self.results_array = []
        self.ligands_array = []
        self.interactions_list = []
        self.counter = 0

        return self.est_time_remaining

    def process_file(self, file_packet):
        results_rows, ligand_row, interaction_rows, receptor_row = file_packet
        for pose in results_rows:
            self.results_array.append(pose)
        for pose in interaction_rows:
            self.interactions_list.append(pose)
        self.ligands_array.append(ligand_row)
        self.receptor_array.append(receptor_row)
        self.counter += 1
