#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail multiprocess workers
#

import time
import sys
from .logutils import get_logger

logger = get_logger(__name__)
import traceback
import queue
from .parsers import docking_file_parsers
from .exceptions import (
    FileParsingErrorAdgpu,
    FileParsingErrorPdbqt,
    FileParsingErrorSdf,
    FileParsingError,
    WriteToStorageError,
    MultiprocessingError,
)
import multiprocessing as mp
from .storagemanager import StorageManager
from .interactions import InteractionFinder


class DockingFileReader(mp.Process):
    """This class is the individual worker for processing docking results.
    One instance of this class is instantiated for each available processor.

    Attributes:
        queueIn (multiprocess.Queue): current queue for the processor/file reader
        queueOut (multiprocess.Queue): queue for the processor/file reader after adding or removing an item
        pipe_conn (multiprocess.Pipe): pipe connection to the reader
        interaction_finder (InteractionFinder): class that calculates interactions
        exception:  for surfacing up an excpetion without crashing the program
        docking_mode (str): docking mode which determines which parser to use
        target_name(str): receptor name to check against certain docking mode files to ensure they belong to the correct target
        num_poses (int): number of poses to store, -1 means all
        interaction_tolerance (float): for docking clusters (eg adgpu) if wanting to store interactions for poses that are part of a cluster
        interaction_cutoffs (list[float,float]): interaction cutoff distances for hb and vdw if calculating interactions
        calculate_interactions (bool):
        receptor_string (str): string representation of receptor to use for calculating interactions
    """

    def __init__(
        self,
        queueIn: mp.Queue,
        queueOut: mp.Queue,
        pipe_conn,
        docking_mode: str,
        target_name: str,
        num_poses: int,
        interaction_tolerance,
        calculate_interactions: bool,
        interaction_cutoffs: list[float, float],
        receptor_string: str,
    ):

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
        self.docking_mode = docking_mode
        self.target_name = target_name
        self.num_poses = num_poses
        self.interaction_tolerance = interaction_tolerance
        self.interaction_cutoffs = interaction_cutoffs
        self.calculate_interactions = calculate_interactions
        self.receptor_string = receptor_string

    def run(self):
        """Method overload from parent class .This is where the task of this class is performed.
        Each multiprocess.Process class must have a "run" method which is called by the
        initialization (see below) with start()

        Raises:
            NotImplementedError: if parser for specific docking result type is not implemented
            FileParsingError
        """
        common_processing_vars = {}
        if self.docking_mode == "adgpu":
            common_processing_vars.update(
                {
                    "target": self.target_name,
                    "interaction_tolerance": self.interaction_tolerance,
                }
            )
        if self.calculate_interactions:
            try:
                interaction_finder = InteractionFinder(
                    self.receptor_string,
                    *self.interaction_cutoffs,
                )
                common_processing_vars.update(
                    {
                        "calculate_interactions": True,
                        "interaction_finder": interaction_finder,
                    }
                )
            except Exception as e:
                logger.warning(
                    f"InteractionFinder setup failed; interactions will not be calculated. Reason: {e}"
                )
                common_processing_vars.update(
                    {
                        "calculate_interactions": False,
                    }
                )

        while True:
            try:
                # retrieve from the queue in the next task to be done
                next_task = self.queueIn.get()
                if isinstance(next_task, dict):
                    text = list(next_task.keys())[0]
                else:
                    text = next_task
                logger.debug("Next Task: " + str(text))
                # if a poison pill is received, this worker's job is done, quit
                if next_task is None:
                    # before leaving, pass the poison pill back in the queue
                    self.queueOut.put(None)
                    break

                # initialize a parser for each process with kw-args
                parser_class = docking_file_parsers.get(self.docking_mode)
                if parser_class is None:
                    raise NotImplementedError(
                        f"Parser for docking_mode {self.docking_mode} not implemented!"
                    )
                parser = parser_class(self.num_poses, **common_processing_vars)

                try:
                    # generate CPU LOAD
                    for data_packet in parser(next_task):
                        self._add_to_queueout(data_packet)
                except FileParsingErrorAdgpu as e:
                    raise FileParsingError(
                        f"Problems when parsing the ADGPU docking log file: {str(e)}"
                    )
                except FileParsingErrorPdbqt as e:
                    raise FileParsingError(
                        f"Problems when parsing the vina docking file: {str(e)}"
                    )
                except FileParsingErrorSdf as e:
                    raise FileParsingError(
                        f"Problems when parsing the SDF docking file: {str(e)}"
                    )

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
        Adds processed object (string, dict, etc) to the output queue which will shuttle it to the Writer

        Args:
            obj (_type_): _description_

        Raises:
            MultiprocessingError: _description_
        """
        max_attempts = 750
        timeout = 0.5  # seconds
        attempts = 0
        while True:
            if attempts >= max_attempts:
                raise MultiprocessingError(
                    "Something is blocking the progressing of file writing. Exiting program."
                ) from queue.Full()
            try:
                self.queueOut.put(obj, block=True, timeout=timeout)
                break
            except queue.Full:
                logger.debug(
                    f"Queue full: queueOut.put attempt {attempts} timed out. {max_attempts - attempts} put attempts remaining."
                )
                attempts += 1


class Writer(mp.Process):
    """This class is a listener that retrieves data from the queue and writes it
    into datbase"""

    def __init__(
        self,
        queue,
        num_readers: int,
        db_file: str,
        storageman_class: StorageManager,
        chunk_size: int,
        duplicate_handling: str,
    ):
        mp.Process.__init__(self)
        self.queue = queue
        # this class knows about how many multi-processing workers there are and where the pipe to the parent is
        self.num_readers = num_readers
        # assign pointer to storage object, set chunksize
        self.storageman: StorageManager = storageman_class(db_file)
        self.chunk_size = chunk_size
        self.duplicate_handling = duplicate_handling
        # initialize data array (stack of dictionaries)
        self.docked_ligands = {"ligands": [], "poses": [], "interactions": []}
        self.receptor_written_to_db = False
        self.receptor_row = None
        # progress tracking instance variables
        self.counter = 0
        self.num_files_written = 0
        self.time0 = time.perf_counter()
        self.total_runtime = 0
        self.last_print_time = 0

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
                        sys.stdout.write(
                            f"\nWrote {self.num_files_written} docking results to the database.\n"
                        )
                        sys.stdout.flush()
                        break
                    continue

                if self.receptor_row is None and not self.receptor_written_to_db:
                    self.receptor_row = list(next_task.get("receptor"))
                self.docked_ligands["ligands"].extend(next_task["ligands"])
                self.docked_ligands["poses"].extend(next_task["poses"])
                self.docked_ligands["interactions"].extend(next_task["interactions"])
                self.counter += 1
                now = time.perf_counter()
                if now - self.last_print_time >= 2.0:
                    self._log_progress()
                    self.last_print_time = now
                if self.counter >= self.chunk_size:
                    self.write_to_storage()

        except Exception:
            tb = traceback.format_exc()
            logger.error("Exception during writing:\n" + tb)
            raise WriteToStorageError("Error occurred while writing to the database.")

    def write_to_storage(self):
        """Inserting data to the database through the designated storagemanager."""
        # insert result, ligand, and receptor data
        with self.storageman as sm:
            if not self.receptor_written_to_db and self.receptor_row:
                sm.insert_receptor_basic_info(self.receptor_row)
                self.receptor_written_to_db = True
                self.receptor_row = None
            sm.insert_data(
                self.docked_ligands,
                self.duplicate_handling,
            )

        # calulate time for processing/writing speed
        self.num_files_written += self.counter
        self.total_runtime = time.perf_counter() - self.time0

        # reset data holder for next chunk
        self.docked_ligands = {"ligands": [], "poses": [], "interactions": []}
        self.counter = 0

    def _log_progress(self):
        current = self.num_files_written + self.counter
        elapsed = time.perf_counter() - self.time0
        rate = current * 60 / (elapsed or 1)
        sys.stdout.write(
            f"\r{current} files processed. {rate:.0f} files/min. "
            f"Elapsed: {elapsed:.0f}s."
        )
        sys.stdout.flush()
