#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail multiprocess workers
#

import time
import sys
from .logutils import LOGGER
import traceback
import queue
from .parsers import docking_file_parsers
from .exceptions import (
    FileParsingError,
    WriteToStorageError,
    MultiprocessingError,
)
import multiprocessing as mp


class DockingFileReader(mp.Process):
    """This class is the individual worker for processing docking results.
    One instance of this class is instantiated for each available processor.

    Attributes:
        queueIn (multiprocess.Queue): current queue for the processor/file reader
        queueOut (multiprocess.Queue): queue for the processor/file reader after adding or removing an item
        pipe_conn (multiprocess.Pipe): pipe connection to the reader
        shared_dict (dict): docking opitions including number of poses, docking mode, duplicate handling, and interaction calculation information
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
                LOGGER.debug("Next Task: " + str(text))
                # if a poison pill is received, this worker's job is done, quit
                if next_task is None:
                    # before leaving, pass the poison pill back in the queue
                    self.queueOut.put(None)
                    break

                # generate CPU LOAD
                # parser depends on requested docking_mode
                data_packet = docking_file_parsers.get(
                    self.docking_mode, "missing_parser"
                )(
                    next_task,
                    self.shared.get("num_poses"),
                    **{
                        "interaction_tolerance": self.shared.get(
                            "interaction_tolerance", None
                        ),
                        "calculate_interactions": not self.shared.get(
                            "no_interactions", True
                        ),
                    },
                )

                if data_packet == "missing_parser":
                    raise NotImplementedError(
                        f"Parser for input file docking_mode {self.docking_mode} not implemented!"
                    )
                # for adgpu only: check receptor name from file against that which we expect
                if (
                    rec_name := data_packet.get("receptor_row")[0]
                    != (target := self.shared.get("target", None))
                    and target is not None
                    and self.docking_mode == "adgpu"
                ):
                    raise FileParsingError(
                        "Receptor name {0} in {1} does not match given target name {2}. Please ensure that this file belongs to the current virtual screening.".format(
                            rec_name,
                            next_task,
                            target,
                        )
                    )

                # put the result in the out queue (will be sent to Writer)
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
                ) from queue.Full
            try:
                self.queueOut.put(obj, block=True, timeout=timeout)
                break
            except queue.Full:
                LOGGER.debug(
                    f"Queue full: queueOut.put attempt {attempts} timed out. {max_attempts - attempts} put attempts remaining."
                )
                attempts += 1


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
                    LOGGER.debug(
                        f"Closing process. Remaining open processes: {self.num_readers}"
                    )
                    if self.num_readers == 0:
                        LOGGER.info("Performing final database write")
                        self.write_to_storage()
                        LOGGER.info("File processing completed")
                        break
                    continue

                if self.receptor_row is None and not self.receptor_written_to_db:
                    self.receptor_row = list(next_task.get("receptor_row"))

                self.docked_ligands["ligands"].append(next_task["ligand_row"])
                self.docked_ligands["poses"].extend(next_task["results_rows"])
                self.docked_ligands["interactions"].extend(
                    next_task["interaction_rows"]
                )
                self.counter += 1

                # After every n (chunk size) files, write to storage
                if self.counter >= self.chunk_size:
                    self._log_progress()
                    self.write_to_storage()
            if self.counter > 0:
                LOGGER.info("Performing final database write")
                self.write_to_storage()
                LOGGER.info("File processing completed")

        except Exception:
            tb = traceback.format_exc()
            LOGGER.error("Exception during writing:\n" + tb)
            raise WriteToStorageError("Error occurred while writing to the database.")

    def write_to_storage(self):
        """Inserting data to the database through the designated storagemanager."""
        # insert result, ligand, and receptor data
        with self.storageman as sm:
            if not self.receptor_written_to_db and self.receptor_row:
                sm.insert_receptor_basic_info(self.receptor_row)
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
