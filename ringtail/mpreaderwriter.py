#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail multiprocessing workers
#

import multiprocessing
import time
import sys
import logging
import traceback
from .parsers import parse_single_dlg, parse_vina_pdbqt
from .exceptions import FileParsingError, WriteToDatabaseError
from .interactions import InteractionFinder


class DockingFileReader(multiprocessing.Process):
    """this class is the individual worker for processing dlgs"""

    def __init__(self, queueIn, queueOut, pipe_conn, dbman, mode, max_poses, interaction_tolerance, store_all_poses, target, add_interactions, interaction_cutoffs, receptor_file):
        # set mode for which file parser to use
        self.mode = mode
        # set number of clusters to write
        self.max_poses = max_poses
        self.store_all_poses_flag = store_all_poses
        # set interaction_tolerance cutoff
        self.interaction_tolerance = interaction_tolerance
        # set options for finding interactions
        self.add_interactions = add_interactions
        self.interaction_cutoffs = interaction_cutoffs
        self.receptor_file = receptor_file
        # set dbmanager
        self.dbman = dbman
        # set target name to check against
        self.target = target
        # initialize the parent class to inherit all multiprocessing methods
        multiprocessing.Process.__init__(self)
        # each worker knows the queue in (where data to process comes from)
        self.queueIn = queueIn
        # ...and a queue out (where to send the results)
        self.queueOut = queueOut
        # ...and a pipe to the parent
        self.pipe = pipe_conn

        self.interaction_finder = None

        self.exception = None

    def _find_best_cluster_poses(self, ligand_dict):
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
        try:
            while True:
                # retrieve from the queue in the next task to be done
                next_task = self.queueIn.get()
                # logging.debug("Next Task: " + str(next_task))
                # if a poison pill is received, this worker's job is done, quit
                if next_task is None:
                    # before leaving, pass the poison pill back in the queue
                    self.queueOut.put(None)
                    break

                # generate CPU LOAD
                # parser depends on requested mode
                if self.mode == "dlg":
                    parsed_file_dict = parse_single_dlg(next_task)
                elif self.mode == "vina":
                    parsed_file_dict = parse_vina_pdbqt(next_task)
                # future: NG parser, etc

                # check receptor name from file against that which we expect
                if (
                    parsed_file_dict["receptor"] != self.target
                    and self.target is not None
                    and self.mode != "vina"
                ):
                    raise FileParsingError(
                        "Receptor name {0} in {1} does not match given target name {2}. Please ensure that this file belongs to the current virtual screening.".format(
                            parsed_file_dict["receptor"], next_task, self.target
                        )
                    )
                # find the run number for the best pose in each cluster for adgpu
                if self.mode == "dlg":
                    parsed_file_dict = self._find_best_cluster_poses(parsed_file_dict)
                # find run numbers for poses we want to save
                parsed_file_dict["poses_to_save"] = self._find_poses_to_save(parsed_file_dict)
                # Calculate interactions if requested
                if self.add_interactions:
                    if self.interaction_finder is None:
                        self.interaction_finder = InteractionFinder(self.receptor_file)
                    if parsed_file_dict["interactions"] == []:
                        for pose in parsed_file_dict["pose_coordinates"]:
                            parsed_file_dict["interactions"].append(self.interaction_finder.find_pose_interactions(parsed_file_dict["ligand_atomtypes"], pose))
                            parsed_file_dict["num_interactions"].append(int(parsed_file_dict["interactions"][-1]["count"][0]))
                            parsed_file_dict["num_hb"].append(len([1 for i in parsed_file_dict["interactions"][-1]["type"] if i == "H"]))
                # find poses we want to save tolerated interactions for
                if self.interaction_tolerance is not None:
                    parsed_file_dict["tolerated_interaction_runs"] = self._find_tolerated_interactions(parsed_file_dict)
                else:
                    parsed_file_dict["tolerated_interaction_runs"] = []
                file_packet = self.dbman.format_rows_from_dict(parsed_file_dict)
                # put the result in the out queue
                self.queueOut.put(file_packet)
        except Exception:
            tb = traceback.format_exc()
            self.pipe.send((FileParsingError("Error while parsing file"), tb))
        finally:
            return

    def _find_poses_to_save(self, ligand_dict: dict) -> list:
        """returns list of the run numbers for the top run in the
        top self.max_pose clusters (ADGPU) or just the top poses overall

        Args:
            ligand_dict (Dictionary): Dictionary of ligand data from parser

        Returns:
            List: List of run numbers to save
        """
        if self.store_all_poses_flag:
            poses_to_save = ligand_dict["sorted_runs"]
        elif self.mode == "dlg":
            # will only select top n clusters. Default 3
            poses_to_save = ligand_dict["cluster_top_poses"][
                : self.max_poses
            ]
        # if not adgpu, save top n poses
        else:
            poses_to_save = ligand_dict["sorted_runs"][:self.max_poses]

        return poses_to_save

    def _find_tolerated_interactions(self, ligand_dict):
        """take ligand dict and finds which poses we should save the
        interactions for as tolerated interactions for the top pose
        of the cluster. These runs are within the
        <self.interaction_tolerance> angstroms RMSD of the top pose
        for a given cluster. All data for the cluster's top pose is saved.

        Args:
            ligand_dict (Dictionary): Dictionary of ligand data from parser

        Returns:
            List: List of run numbers of tolerated runs
        """
        tolerated_runs = []
        for idx, run in enumerate(ligand_dict["sorted_runs"]):
            if (
                float(ligand_dict["cluster_rmsds"][idx])
                <= self.interaction_tolerance
            ):
                tolerated_runs.append(run)
        return tolerated_runs


class Writer(multiprocessing.Process):
    # this class is a listener that retrieves data from the queue and writes it
    # into datbase
    def __init__(
        self, queue, maxProcesses, pipe_conn, chunksize, db_obj, num_files, mode="dlg"
    ):
        multiprocessing.Process.__init__(self)
        self.queue = queue
        # this class knows about how many multi-processing workers there are
        self.maxProcesses = maxProcesses
        # assign pointer to db object, set chunksize
        self.mode = mode
        self.db = db_obj
        self.chunksize = chunksize
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

        self.pipe = pipe_conn

    def run(self):
        # method overload from parent class
        #
        # this is where the task of this class is performed
        #
        # each multiprocessing.Process class must have a "run" method which
        # is called by the initialization (see below) with start()
        #
        try:
            while True:
                # retrieve the next task from the queue
                next_task = self.queue.get()
                if next_task is None:
                    # if a poison pill is found, it means one of the workers quit
                    self.maxProcesses -= 1
                    logging.info(
                        "Closing process. Remaining open processes: " + str(self.maxProcesses)
                    )
                else:
                    # if not a poison pill, process the task item

                    # after every n (chunksize) files, write to db
                    if self.counter >= self.chunksize:
                        self.write_to_db()

                    # process next file
                    self.process_file(next_task)
                    # print info about files and time remaining
                    sys.stdout.write("\r")
                    sys.stdout.write(
                        "{n} files remaining to process. Estimated time remaining: {est_time:.2f} minutes. Last chunk took {write_time:.2f} seconds to write.".format(
                            n=self.num_files_remaining,
                            est_time=self.est_time_remaining,
                            write_time=self.last_write_time,
                        )
                    )
                    sys.stdout.flush()
                    self.num_files_remaining -= 1

                if self.maxProcesses == 0:
                    # received as many poison pills as workers
                    logging.info("Performing final database write")
                    # perform final db write
                    self.write_to_db()
                    # no workers left, no job to do
                    logging.info("File processing completed")
                    self.close()
                    break
        except Exception:
            tb = traceback.format_exc()
            self.pipe.send(
                (WriteToDatabaseError("Error occured while writing database"), tb)
            )
        finally:
            return

    def write_to_db(self):
        time1 = time.perf_counter()
        # insert result, ligand, and receptor data
        self.db.insert_results(
            filter(None, self.results_array)
        )  # filter out stray Nones
        self.db.insert_ligands(filter(None, self.ligands_array))
        # if this is the first insert or we have multiple receptors, insert the receptor array
        if self.first_insert and filter(None, self.receptor_array) != []:
            self.db.insert_receptors(filter(None, self.receptor_array))
            self.first_insert = False

        # perform operations for inserting iteraction data
        if self.interactions_list != []:
            self.db.insert_interactions(self.interactions_list)

        # calulate time for processing/writing previous chunk
        self.num_files_written += self.chunksize
        total_runtime = time.perf_counter() - self.time0
        self.est_time_remaining = (
            (total_runtime / self.num_files_written)
            * (self.total_num_files - self.num_files_written)
            / 60
        )  # converted to minutes
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
