#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail multiprocess workers
#

import platform
import time
import sys
from .logutils import LOGGER as logger
import traceback
import queue
from .parsers import parse_single_dlg, parse_vina_result
from .exceptions import (
    FileParsingError,
    WriteToStorageError,
    MultiprocessingError,
    ResultsProcessingError,
)
from .interactions import InteractionFinder

import multiprocess


class DockingFileReader(multiprocess.Process):
    """This class is the individual worker for processing docking results.
    One instance of this class is instantiated for each available processor.

    Attributes:
        queueIn (multiprocess.Queue): current queue for the processor/file reader
        queueOut (multiprocess.Queue): queue for the processor/file reader after adding or removing an item
        pipe_conn (multiprocess.Pipe): pipe connection to the reader
        storageman (StorageManager): storageman object
        storageman_class (StorageManager): storagemanager child class/database type
        docking_mode (str): describes what docking engine was used to produce the results
        max_poses (int): max number of poses to store for each ligand
        interaction_tolerance (float): Will add the interactions for poses within some tolerance RMSD range of the top pose in a cluster to that top pose."
        store_all_poses (bool): Store all poses from docking results
        add_interactions (bool): find and save interactions between ligand poses and receptor
        interaction_cutoffs (list(float)): cutoff for interactions of hydrogen bonds and VDW interactions, in ångströms
        target (str): receptor name
    """

    def __init__(
        self,
        queueIn,
        queueOut,
        pipe_conn,
        storageman,
        storageman_class,
        docking_mode,
        max_poses,
        interaction_tolerance,
        store_all_poses,
        target,
        add_interactions,
        interaction_cutoffs,
        receptor_file,
    ):
        # set docking_mode for which file parser to use (and for vina, '_string' if parsing string output directly)
        self.docking_mode = docking_mode
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
        # initialize the parent class to inherit all multiprocess methods
        multiprocess.Process.__init__(self)
        # each worker knows the queue in (where data to process comes from)
        self.queueIn = queueIn
        # ...and a queue out (where to send the results)
        self.queueOut = queueOut
        # ...and a pipe to the parent
        self.pipe = pipe_conn
        self.interaction_finder = None
        self.exception = None

    def _find_best_cluster_poses(self, ligand_dict):
        """Takes input ligand dictionary, reads run pose clusters, adds "cluster_best_run"
        entry with the top scoring run for each cluster

        Args:
            ligand_dict (dict): dictionary of ligands

        Returns:
            dict: top poses in cluster for each ligand
        """

        top_poses = []
        cluster_dict = ligand_dict["clusters"]
        for cluster in cluster_dict:
            top_poses.append(cluster_dict[cluster][0])
        ligand_dict["cluster_top_poses"] = top_poses
        return ligand_dict

    def run(self):
        """Method overload from parent class .This is where the task of this class is performed.
        Each multiprocess.Process class must have a "run" method which is called by the
        initialization (see below) with start()

        Raises:
            NotImplementedError: if parser for specific docking result type is not implemented
            FileParsingError
        """

        while True:
            try:
                # retrieve from the queue in the next task to be done
                next_task = self.queueIn.get()
                if type(next_task) == dict:
                    text = list(next_task.keys())[0]
                else:
                    text = next_task
                logger.debug("Next Task: " + str(text))
                # if a poison pill is received, this worker's job is done, quit
                if next_task is None:
                    # before leaving, pass the poison pill back in the queue
                    self.queueOut.put(None)
                    break

                # generate CPU LOAD
                # parser depends on requested docking_mode
                if self.docking_mode == "dlg":
                    parsed_file_dict = parse_single_dlg(next_task)
                    # find the run number for the best pose in each cluster for adgpu
                    parsed_file_dict = self._find_best_cluster_poses(parsed_file_dict)
                elif self.docking_mode == "vina":
                    parsed_file_dict = parse_vina_result(next_task)

                # Example code for calling user-implemented docking_mode
                # elif self.docking_mode == "my_docking_mode":
                #     parsed_file_dict = myparser(next_task)
                else:
                    raise NotImplementedError(
                        f"Parser for input file docking_mode {self.docking_mode} not implemented!"
                    )
                # check receptor name from file against that which we expect
                if (
                    parsed_file_dict["receptor"] != self.target
                    and self.target is not None
                    and self.docking_mode == "dlg"
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
                    try:
                        from .receptormanager import ReceptorManager as rm

                        with self.storageman:
                            # grab receptor info from database, this assumes there is only one receptor in the database
                            receptor_blob = self.storageman.fetch_receptor_objects()[0][
                                1
                            ]  # method returns and iter of tuples, blob is the second tuple element in the first list element
                            # convert receptor blob to string
                            receptor_string = rm.blob2str(receptor_blob)
                    except:
                        raise ResultsProcessingError(
                            "add_interactions was requested, but cannot find the receptor in the database. Please ensure to include the receptor_file and save_receptor if the receptor has not already been added to the database."
                        )
                    if self.interaction_finder is None:
                        self.interaction_finder = InteractionFinder(
                            receptor_string, self.interaction_cutoffs
                        )
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
                    parsed_file_dict["tolerated_interaction_runs"] = (
                        self._find_tolerated_interactions(parsed_file_dict)
                    )
                else:
                    parsed_file_dict["tolerated_interaction_runs"] = []
                # put the result in the out queue
                data_packet = self.storageman_class.format_for_storage(parsed_file_dict)
                self._add_to_queueout(data_packet)
            except Exception:
                tb = traceback.format_exc()
                self.pipe.send(
                    (
                        FileParsingError(f"Error while parsing {next_task}"),
                        tb,
                        next_task,
                    )
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
                logger.debug(
                    f"Queue full: queueOut.put attempt {attempts} timed out. {max_attempts - attempts} put attempts remaining."
                )
                attempts += 1

    def _find_poses_to_save(self, ligand_dict: dict) -> list:
        """Returns list of the run numbers for the top run in the
        top self.max_pose clusters (ADGPU) or just the top poses overall

        Args:
            ligand_dict (dict): Dictionary of ligand data from parser

        Returns:
            list: List of run numbers to save
        """
        if self.store_all_poses_flag:
            poses_to_save = ligand_dict["sorted_runs"]
        elif self.docking_mode == "dlg":
            # will only select top n clusters. Default 3
            poses_to_save = ligand_dict["cluster_top_poses"][: self.max_poses]
        # if not adgpu, save top n poses
        else:
            poses_to_save = ligand_dict["sorted_runs"][: self.max_poses]

        return poses_to_save

    def _find_tolerated_interactions(self, ligand_dict):
        """Take ligand dict and finds which poses we should save the
        interactions for as tolerated interactions for the top pose
        of the cluster. These runs are within the
        <self.interaction_tolerance> angstroms RMSD of the top pose
        for a given cluster. All data for the cluster's top pose is saved.

        Args:
            ligand_dict (dict): Dictionary of ligand data from parser

        Returns:
            list: run numbers of tolerated runs
        """
        tolerated_runs = []
        for idx, run in enumerate(ligand_dict["sorted_runs"]):
            if float(ligand_dict["cluster_rmsds"][idx]) <= self.interaction_tolerance:
                tolerated_runs.append(run)
        return tolerated_runs


class Writer(multiprocess.Process):
    """This class is a listener that retrieves data from the queue and writes it
    into datbase"""

    def __init__(
        self, queue, num_readers, pipe_conn, chunksize, storageman, docking_mode
    ):
        multiprocess.Process.__init__(self)
        self.queue = queue
        # this class knows about how many multi-processing workers there are and where the pipe to the parent is
        self.num_readers = num_readers
        self.pipe = pipe_conn
        # assign pointer to storage object, set chunksize
        self.docking_mode = docking_mode
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
        """Method overload from parent class. This is where the task of this class
        is performed. Each multiprocess.Process class must have a "run" method which
        is called by the initialization (see below) with start()

        Raises:
            WriteToStorageError
        """

        try:
            while True:
                # retrieve the next task from the queue
                next_task = self.queue.get()
                if next_task is None:
                    # if a poison pill is found, it means one of the workers quit
                    self.num_readers -= 1
                    logger.debug(
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
                                self.num_files_written + 1,
                                self.num_files_written * 60 / self.total_runtime,
                                self.total_runtime,
                            )
                        )
                        sys.stdout.flush()

                    # process next file
                    self.process_data(next_task)
                if self.num_readers == 0:
                    # received as many poison pills as workers
                    logger.info("Performing final database write")
                    # perform final storage write
                    self.write_to_storage()
                    # no workers left, no job to do
                    logger.info("File processing completed")
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

    def write_to_storage(self):
        """Inserting data to the database through the designated storagemanager."""
        # insert result, ligand, and receptor data
        self.storageman.insert_data(
            self.results_array,
            self.ligands_array,
            self.interactions_list,
            self.receptor_array,
            self.first_insert,
        )
        # So at this point the ligand array is empty
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

    def process_data(self, data_packet):
        """Breaks up the data in the data_packet to distribute between
        the different arrays to be inserted in the database.

        Args:
            data_packet (any): File packet to be processed
        """
        results_rows, ligand_row, interaction_rows, receptor_row = data_packet
        for pose in results_rows:
            self.results_array.append(pose)
        for pose in interaction_rows:
            self.interactions_list.append(pose)
        self.ligands_array.append(ligand_row)
        self.receptor_array.append(receptor_row)

        self.counter += 1
