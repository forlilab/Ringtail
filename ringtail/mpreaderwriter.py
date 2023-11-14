#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail multiprocessing workers
#

import platform
import time
import sys
import logging
import traceback
import queue
from .parsers import parse_single_dlg, parse_vina_pdbqt
from .exceptions import FileParsingError, WriteToStorageError, MultiprocessingError
from .interactions import InteractionFinder

os_string = platform.system()
if os_string == "Darwin":  # mac
    import multiprocess as multiprocessing
else:
    import multiprocessing


class DockingFileReader(multiprocessing.Process):
    """this class is the individual worker for processing dlgs"""

    def __init__(
        self,
        queueIn,
        queueOut,
        pipe_conn,
        storageman,
        storageman_class,
        mode,
        max_poses,
        interaction_tolerance,
        store_all_poses,
        target,
        add_interactions,
        interaction_cutoffs,
        receptor_file,
    ):
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
        # set storagemanager and class
        self.storageman = storageman
        self.storageman_class = storageman_class
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
        while True:
            try:
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
                    # find the run number for the best pose in each cluster for adgpu
                    parsed_file_dict = self._find_best_cluster_poses(parsed_file_dict)
                elif self.mode == "vina":
                    parsed_file_dict = parse_vina_pdbqt(next_task)
                # Example code for calling user-implemented mode
                # elif self.mode == "my_mode":
                #     parsed_file_dict = myparser(next_task)
                else:
                    raise NotImplementedError(f"Parser for input file mode {self.mode} not implemented!")

                # check receptor name from file against that which we expect
                if (
                    parsed_file_dict["receptor"] != self.target
                    and self.target is not None
                    and self.mode == "dlg"
                ):
                    raise FileParsingError(
                        "Receptor name {0} in {1} does not match given target name {2}. Please ensure that this file belongs to the current virtual screening.".format(
                            parsed_file_dict["receptor"], next_task, self.target
                        )
                    )
                
                # find run numbers for poses we want to save
                parsed_file_dict["poses_to_save"] = self._find_poses_to_save(
                    parsed_file_dict
                )
                # Calculate interactions if requested
                if self.add_interactions:
                    if self.interaction_finder is None:
                        self.interaction_finder = InteractionFinder(self.receptor_file)
                    if parsed_file_dict["interactions"] == []:
                        for pose in parsed_file_dict["pose_coordinates"]:
                            parsed_file_dict["interactions"].append(
                                self.interaction_finder.find_pose_interactions(
                                    parsed_file_dict["ligand_atomtypes"], pose
                                )
                            )
                            parsed_file_dict["num_interactions"].append(
                                int(parsed_file_dict["interactions"][-1]["count"][0])
                            )
                            parsed_file_dict["num_hb"].append(
                                len(
                                    [
                                        1
                                        for i in parsed_file_dict["interactions"][-1][
                                            "type"
                                        ]
                                        if i == "H"
                                    ]
                                )
                            )
                # find poses we want to save tolerated interactions for
                if self.interaction_tolerance is not None:
                    parsed_file_dict[
                        "tolerated_interaction_runs"
                    ] = self._find_tolerated_interactions(parsed_file_dict)
                else:
                    parsed_file_dict["tolerated_interaction_runs"] = []
                # put the result in the out queue
                data_packet = self.storageman_class.format_for_storage(parsed_file_dict)
                self._add_to_queueout(data_packet)
            except Exception:
                tb = traceback.format_exc()
                self.pipe.send(
                    (FileParsingError(f"Error while parsing {next_task}"), tb, next_task)
                )

    def _add_to_queueout(self, obj):
        max_attempts = 750
        timeout = 0.5  # seconds
        attempts = 0
        while True:
            if attempts >= max_attempts:
                raise MultiprocessingError(
                    "Something is blocking the progressing of file writing. Exiting program."
                ) from queue.Full
            try:
                self.queueOut.put(obj, block=True, timeout=timeout)
                break
            except queue.Full:
                # logging.debug(f"Queue full: queueOut.put attempt {attempts} timed out. {max_attempts - attempts} put attempts remaining.")
                attempts += 1

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
            poses_to_save = ligand_dict["cluster_top_poses"][: self.max_poses]
        # if not adgpu, save top n poses
        else:
            poses_to_save = ligand_dict["sorted_runs"][: self.max_poses]

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
            if float(ligand_dict["cluster_rmsds"][idx]) <= self.interaction_tolerance:
                tolerated_runs.append(run)
        return tolerated_runs


class Writer(multiprocessing.Process):
    # this class is a listener that retrieves data from the queue and writes it
    # into datbase
    def __init__(
        self, queue, num_readers, pipe_conn, chunksize, storageman, mode="dlg"
    ):
        multiprocessing.Process.__init__(self)
        self.queue = queue
        # this class knows about how many multi-processing workers there are and where the pipe to the parent is
        self.num_readers = num_readers
        self.pipe = pipe_conn
        # assign pointer to storage object, set chunksize
        self.mode = mode
        self.storageman = storageman
        self.chunksize = chunksize
        # initialize data array (stack of dictionaries)
        self.results_array = []
        self.ligands_array = []
        self.interactions_list = []
        self.receptor_array = []
        # progress tracking instance variables
        self.first_insert = True
        self.counter = 0
        self.num_files_written = 0
        self.time0 = time.perf_counter()
        self.last_write_time = 0

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
                    self.num_readers -= 1
                    logging.info(
                        "Closing process. Remaining open processes: "
                        + str(self.num_readers)
                    )
                else:
                    # if not a poison pill, process the task item

                    # after every n (chunksize) files, write to storage
                    if self.counter >= self.chunksize:
                        self.write_to_storage()
                        # print info about files and time remaining
                        sys.stdout.write("\r")
                        sys.stdout.write(
                            "{0} files written to database. Writing {1:.0f} files/minute. Elapsed time {2:.0f} seconds.".format(
                                self.num_files_written,
                                self.num_files_written * 60 / self.total_runtime,
                                self.total_runtime,
                            )
                        )
                        sys.stdout.flush()

                    # process next file
                    self.process_file(next_task)

                if self.num_readers == 0:
                    # received as many poison pills as workers
                    logging.info("Performing final database write")
                    # perform final storage write
                    self.write_to_storage(final=True)
                    # no workers left, no job to do
                    logging.info("File processing completed")
                    self.close()
                    break
        except Exception:
            tb = traceback.format_exc()
            self.pipe.send(
                (
                    WriteToStorageError("Error occured while writing database"),
                    tb,
                    "Database",
                )
            )

    def write_to_storage(self, final=False):
        # insert result, ligand, and receptor data
        self.storageman.insert_data(
            self.results_array,
            self.ligands_array,
            self.interactions_list,
            self.receptor_array,
            self.first_insert,
        )
        if self.first_insert:  # will only insert receptor for first insertion
            self.first_insert = False

        # calulate time for processing/writing speed
        self.num_files_written += self.chunksize
        self.total_runtime = time.perf_counter() - self.time0

        # reset data holder for next chunk
        self.results_array = []
        self.ligands_array = []
        self.interactions_list = []
        self.receptor_array = []
        self.counter = 0

        if final:
            # if final write, tell storageman to index
            self.storageman.create_indices()
            self.storageman.set_ringtaildb_version()

    def process_file(self, file_packet):
        results_rows, ligand_row, interaction_rows, receptor_row = file_packet
        for pose in results_rows:
            self.results_array.append(pose)
        for pose in interaction_rows:
            self.interactions_list.append(pose)
        self.ligands_array.append(ligand_row)
        self.receptor_array.append(receptor_row)

        self.counter += 1
