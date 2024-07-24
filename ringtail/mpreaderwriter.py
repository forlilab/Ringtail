#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail multiprocessing workers
#

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
from .storagemanager import *
import multiprocessing


class DockingFileReader(
    multiprocessing.Process,
):
    """This class is the individual worker for processing docking results.
    One instance of this class is instantiated for each available processor.

    Attributes:
        queueIn (multiprocessing.Manager.Queue): current queue for the processor/file reader
        queueOut (multiprocessing.Manager.Queue): queue for the processor/file reader after adding or removing an item
        pipe_conn (multiprocessing.Pipe): pipe connection to the reader
        docking_mode (str): describes what docking engine was used to produce the results
        storageman_class (StorageManager): storagemanager child class/database type
        max_poses (int): max number of poses to store for each ligand
        store_all_poses (bool): Store all poses from docking results
        add_interactions (bool): find and save interactions between ligand poses and receptor
        interaction_tolerance (float): Will add the interactions for poses within some tolerance RMSD range of the top pose in a cluster to that top pose."
        interaction_cutoffs (list(float)): cutoff for interactions of hydrogen bonds and VDW interactions, in ångströms
        interaction_finder (InteractionFinder): object that processes ligand and receptor data to find interactions
        receptor_blob (str): the complete receptor description
        target (str): receptor name
    """

    def __init__(
        self,
        queueIn,
        queueOut,
        pipe_conn,
        storageman_class,
        docking_mode,
        max_poses,
        interaction_tolerance,
        store_all_poses,
        target,
        add_interactions,
        interaction_cutoffs,
        string_processing=False,
        receptor_blob=None,
    ):
        # initialize parent Process class
        multiprocessing.Process.__init__(self)
        # attributes related to multiprocessing
        self.queueIn = queueIn
        self.queueOut = queueOut
        self.pipe = pipe_conn
        # attributes related to the database
        self.docking_mode = docking_mode
        self.storageman_class = storageman_class
        self.string_processing = string_processing
        # attributes related to data processing
        self.max_poses = max_poses
        self.store_all_poses_flag = store_all_poses
        self.interaction_tolerance = interaction_tolerance
        self.interaction_cutoffs = interaction_cutoffs
        self.add_interactions = add_interactions
        # receptor information, needed if adding interactions
        self.receptor_blob = receptor_blob
        self.target = target
        self.interaction_finder = None
        # self.exception = None

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
        Each multiprocessing.Process class must have a "run" method which is called by the
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
                    parsed_file_dict = parse_vina_result(
                        next_task, self.string_processing
                    )

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
                    from .receptormanager import ReceptorManager as rm

                    receptor_string = rm.blob2str(self.receptor_blob)

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
        """
        Method that adds processed data packet to queue out, which will then be written to database by the Writer class

        Args:
            obj (any): Data packet that is ready to be added to out queue and written to db

        Raises:
            MultiprocessingError
        """
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


class Writer(multiprocessing.Process):
    """This class is a listener that retrieves data from the queue and writes it
    into datbase

    Attributes:
        queue (multiprocessing.Manager.Queue): incoming queue of objects to be written to db
        num_readers (int): number of CPUs on which where Writer class is active
        pipe_conn (multiprocessing.Pipe connection):
        chunksize (int): decides how many docking results are processed before writing to database
        docking_mode (str): docking mode used for the results
        storageman_class (StorageManager): type of database used so data is correctly formatted before writing
        duplicate_handling (bool): how to handle duplicate entries in the database
        db_file (str): database file for which to connect and write results to
        results_array (list): holds db-specific formatted data for the results table
        ligands_array (list): holds db-specific formatted data for the ligands table
        interactions_list (list): holds db-specific formatted data for the interactions table
        receptor_array (list): holds db-specific formatted data for the receptor table
        first_insert (bool): first db writing is slightly different than the remainders
        self.counter (int): keeps track of how many docking results are in the processed data lists
        num_files_written (int): how many files/docking results have been written to the database
        time0 (time): start time for the first file written to db
        last_write_time (time): end time for the last file written to db
    """

    def __init__(
        self,
        queue,
        num_readers,
        pipe_conn,
        chunksize,
        docking_mode,
        duplicate_handling,
        db_file,
        storageman_class,
    ):
        # initialize parent Process class
        multiprocessing.Process.__init__(self)
        # attributes related to multiprocessing
        self.queue = queue
        self.num_readers = num_readers
        self.pipe = pipe_conn
        self.chunksize = chunksize
        # attributes related to the database
        self.docking_mode = docking_mode
        self.storageman_class = storageman_class
        self.duplicate_handling = duplicate_handling
        self.db_file = db_file
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
        is performed. Each multiprocessing.Process class must have a "run" method which
        is called by the initialization (see below) with start()

        Raises:
            WriteToStorageError
        """
        self.storageman = self.storageman_class(self.db_file)
        self.storageman.duplicate_handling = self.duplicate_handling
        with self.storageman:
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
                        logger.info("Performing final database write")
                        # perform final storage write
                        self.write_to_storage(final=True)
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

    def write_to_storage(self, final=False):
        """Inserting data to the database through the designated storagemanager.

        Args:
            final (bool): if last data entry, finalize database
        """
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

        if final:
            # if final write, tell storageman to index
            self.storageman.create_indices()
            self.storageman.set_ringtail_db_schema_version()

    def process_file(self, file_packet):
        """Breaks up the data in the file_packet to distribute between
        the different arrays to be inserted in the database.

        Args:
            file_packet (any): File packet to be processed
        """
        results_rows, ligand_row, interaction_rows, receptor_row = file_packet
        for pose in results_rows:
            self.results_array.append(pose)
        for pose in interaction_rows:
            self.interactions_list.append(pose)
        self.ligands_array.append(ligand_row)
        self.receptor_array.append(receptor_row)

        self.counter += 1
