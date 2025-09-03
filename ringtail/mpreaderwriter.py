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
)
from .interactions import InteractionFinder
import multiprocessing as mp


class DockingFileReader(mp.Process):
    """This class is the individual worker for processing docking results.
    One instance of this class is instantiated for each available processor.

    Attributes:
        queueIn (multiprocess.Queue): current queue for the processor/file reader
        queueOut (multiprocess.Queue): queue for the processor/file reader after adding or removing an item
        pipe_conn (multiprocess.Pipe): pipe connection to the reader
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
        shared_dict,
    ):
        self.shared = shared_dict
        # initialize the parent class to inherit all multiprocess methods
        mp.Process.__init__(self)
        # each worker knows the queue in (where data to process comes from)
        self.queueIn = queueIn
        # ...and a queue out (where to send the results)
        self.queueOut = queueOut
        # ...and a pipe to the parent
        self.pipe = pipe_conn
        self.interaction_finder = None
        self.exception = None
        self.docking_mode = self.shared.get("docking_mode")

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
                if self.docking_mode == "adgpu":
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
                    parsed_file_dict["receptor"] != self.shared["target"]
                    and self.shared.get("target") is not None
                    and self.docking_mode == "adgpu"
                ):
                    raise FileParsingError(
                        "Receptor name {0} in {1} does not match given target name {2}. Please ensure that this file belongs to the current virtual screening.".format(
                            parsed_file_dict["receptor"],
                            next_task,
                            self.shared.get("target"),
                        )
                    )

                # find run numbers for poses we want to save
                # NOTE where num poses plays in
                parsed_file_dict["poses_to_save"] = self._find_poses_to_save(
                    parsed_file_dict
                )
                # Calculate interactions if requested
                if self.shared.get("add_interactions"):
                    if self.interaction_finder is None:
                        self.interaction_finder = InteractionFinder(
                            self.shared.get("receptor_string"),
                            self.shared.get("interaction_cutoffs"),
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
                if self.shared.get("interaction_tolerance") is not None:
                    parsed_file_dict["tolerated_interaction_runs"] = (
                        self._find_tolerated_interactions(parsed_file_dict)
                    )
                else:
                    parsed_file_dict["tolerated_interaction_runs"] = []
                # put the result in the out queue
                # the result is a dict where the singular key is the ligand name
                # TODO this is not working right now, the data_packet should be expected format for writer
                data_packet = self.shared.get("format_method")(parsed_file_dict)
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
        store_all_poses = self.shared.get("store_all_poses")
        max_poses = self.shared.get("max_poses")
        if max_poses:
            if max_poses > len(ligand_dict["sorted_runs"]):
                store_all_poses = True
        if store_all_poses:
            poses_to_save = ligand_dict["sorted_runs"]
        elif self.docking_mode == "adgpu":
            # will only select top n clusters.
            poses_to_save = ligand_dict["cluster_top_poses"][:max_poses]
        # if not adgpu, save top n poses
        else:
            poses_to_save = ligand_dict["sorted_runs"][:max_poses]

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
            if float(ligand_dict["cluster_rmsds"][idx]) <= self.shared.get(
                "interaction_tolerance"
            ):
                tolerated_runs.append(run)
        return tolerated_runs


class Writer(mp.Process):
    """This class is a listener that retrieves data from the queue and writes it
    into datbase"""

    def __init__(self, queue, num_readers: int, options: dict):
        mp.Process.__init__(self)
        self.queue = queue
        # this class knows about how many multi-processing workers there are and where the pipe to the parent is
        self.num_readers = num_readers
        # assign pointer to storage object, set chunksize
        db_file = options.pop("db_file")
        storage_class = options.pop("storageman_class")
        self.storageman = storage_class(db_file)
        self.chunk_size = options.pop("chunk_size")
        self.options = options
        # initialize data array (stack of dictionaries)
        self.docked_ligands = {"ligands": [], "poses": [], "interactions": []}
        self.receptor_written_to_db = False
        self.receptor_row = None
        # progress tracking instance variables
        self.counter = 0
        self.num_files_written = 0
        self.time0 = time.perf_counter()
        self.total_runtime = 0
        self.last_write_time = 0

    def run(self):
        self.time0 = time.perf_counter()

        try:
            while True:
                next_task = self.queue.get()
                if next_task is None:
                    self.num_readers -= 1
                    logger.debug(
                        f"Closing process. Remaining open processes: {self.num_readers}"
                    )
                    if self.num_readers == 0:
                        logger.info("Performing final database write")
                        self.write_to_storage()
                        logger.info("File processing completed")
                        break
                    continue

                if self.receptor_row is None and not self.receptor_written_to_db:
                    self.receptor_row = list(next_task.get("receptor_row"))

                self.docked_ligands["ligands"].append(next_task["ligand"])
                self.docked_ligands["poses"].extend(next_task["poses"])
                self.docked_ligands["interactions"].extend(next_task["interactions"])
                self.counter += 1

                # After every n (chunk size) files, write to storage
                if self.counter >= self.chunk_size:
                    self._log_progress()
                    self.write_to_storage()
            if self.counter > 0:
                logger.info("Performing final database write")
                self.write_to_storage()
                logger.info("File processing completed")

        except Exception:
            tb = traceback.format_exc()
            logger.error("Exception during writing:\n" + tb)
            raise WriteToStorageError("Error occurred while writing to the database.")

    def write_to_storage(self):
        """Inserting data to the database through the designated storagemanager."""
        # insert result, ligand, and receptor data
        with self.storageman as sm:
            if not self.receptor_written_to_db and self.receptor_row:
                sm.insert_receptor(self.receptor_row)
                self.receptor_written_to_db = True
                self.receptor_row = None
            if self.docked_ligands["ligands"]:
                sm.insert_data(
                    self.docked_ligands,
                    self.options,
                )

        # calulate time for processing/writing speed
        self.num_files_written += self.counter
        self.total_runtime = time.perf_counter() - self.time0

        # reset data holder for next chunk
        self.docked_ligands = {"ligands": [], "poses": [], "interactions": []}
        self.counter = 0

    def _log_progress(self):
        rate = self.num_files_written * 60 / (self.total_runtime or 1)
        sys.stdout.write(
            f"\r{self.num_files_written} files written. {rate:.0f} files/min. "
            f"Elapsed: {self.total_runtime:.0f}s."
        )
        sys.stdout.flush()
