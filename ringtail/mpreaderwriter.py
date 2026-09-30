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

# tag the writer puts on pipe messages, so the parent treats them as fatal
WRITER_ERROR_SOURCE = "Database"


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
        super().__init__()
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
        # ad6 has multiple ligands per file, so need to keep processing file if one bad ligand/record
        if self.docking_mode == "ad6":
            common_processing_vars["report_error"] = self._report_bad_record
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
        Adds a parsed result to the output queue for the Writer, waiting while it is full.

        Args:
            obj (dict): parsed docking result

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
                ) from queue.Full()
            try:
                self.queueOut.put(obj, block=True, timeout=timeout)
                break
            except queue.Full:
                logger.debug(
                    f"Queue full: queueOut.put attempt {attempts} timed out. {max_attempts - attempts} put attempts remaining."
                )
                attempts += 1

    def _report_bad_record(self, name: str, tb: str):
        """Sends a skipped record to the parent for the failed-files log."""
        self.pipe.send((FileParsingError(f"Error while parsing {name}"), tb, name))


class Writer(mp.Process):
    """Listener that writes parsed docking results from the queue into the database.

    Args:
        queue (mp.Queue): parsed results, ending with one None per reader
        pipe_conn: connection to the parent, used to report a failed write
        num_readers (int): number of readers to wait for
        db_file (str): database file
        storageman_class (type[StorageManager]): storage manager class for the database
        chunk_size (int): number of results to buffer before each write
        duplicate_handling (str): how to handle duplicate Results rows
    """

    def __init__(
        self,
        queue,
        pipe_conn,
        num_readers: int,
        db_file: str,
        storageman_class: type[StorageManager],
        chunk_size: int,
        duplicate_handling: str,
    ):
        super().__init__()
        self.queue = queue
        self.pipe = pipe_conn
        # this class knows about how many multi-processing workers there are and where the pipe to the parent is
        self.num_readers = num_readers
        # assign pointer to storage object, set chunksize
        self.storageman: StorageManager = storageman_class(db_file)
        self.chunk_size = chunk_size
        self.duplicate_handling = duplicate_handling
        # parsed results of the current chunk, one per docking file or SDF record
        self.packets = []
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
                        if self.num_files_written:
                            sys.stdout.write(
                                f"\nWrote {self.num_files_written} docking results to the database.\n"
                            )
                        else:
                            # Parse failures are logged and skipped rather than fatal, so
                            # when every file fails this was the only thing the user saw,
                            # phrased as though it had succeeded.
                            sys.stdout.write(
                                "\nNo docking results were written to the database: none "
                                "of the provided files could be parsed. See "
                                "ringtail_failed_files.log and the log above for why.\n"
                            )
                        sys.stdout.flush()
                        break
                    continue

                if self.receptor_row is None and not self.receptor_written_to_db:
                    self.receptor_row = list(next_task.get("receptor"))
                self.packets.append(next_task)
                self.counter += 1
                now = time.perf_counter()
                if now - self.last_print_time >= 2.0:
                    self._log_progress()
                    self.last_print_time = now
                if self.counter >= self.chunk_size:
                    self.write_to_storage()

        except Exception as e:
            tb = traceback.format_exc()
            error = WriteToStorageError(f"Error occurred while writing to the database: {e}")
            # child logging is not visible under spawn, so report to the parent
            self.pipe.send((error, tb, WRITER_ERROR_SOURCE))
            raise error

    def write_to_storage(self):
        """Inserting data to the database through the designated storagemanager."""
        # insert receptor data
        with self.storageman as sm:
            if not self.receptor_written_to_db and self.receptor_row:
                sm.insert_receptor_basic_info(self.receptor_row)
                self.receptor_written_to_db = True
                self.receptor_row = None
            last_ids = sm.last_row_ids()

        # insert ligand, result and interaction data, one ligand at a time if the chunk fails
        try:
            self._insert(self.packets)
            self.num_files_written += len(self.packets)
        except Exception:
            self._undo_write(last_ids)
            self.num_files_written += self._insert_one_by_one()

        # calulate time for processing/writing speed
        self.total_runtime = time.perf_counter() - self.time0

        # reset data holder for next chunk
        self.packets = []
        self.counter = 0

    def _insert(self, packets: list):
        """Writes the ligands, poses and interactions of the given parsed results in one go."""
        data = {
            key: [row for packet in packets for row in packet[key]]
            for key in ("ligands", "poses", "interactions")
        }
        with self.storageman as sm:
            sm.insert_data(data, self.duplicate_handling)

    def _undo_write(self, last_ids: dict):
        """Deletes what a failed write left in the database."""
        with self.storageman as sm:
            sm.delete_rows_after(last_ids)

    def _insert_one_by_one(self) -> int:
        """Writes the chunk's parsed results one at a time, reporting the ones that fail.

        Raises:
            WriteToStorageError: if every result of a chunk of several fails

        Returns:
            int: number of results written
        """
        written = 0
        for packet in self.packets:
            with self.storageman as sm:
                last_ids = sm.last_row_ids()
            try:
                self._insert([packet])
                written += 1
            except Exception as e:
                self._undo_write(last_ids)
                if packet["ligands"]:
                    name = packet["ligands"][0][0]
                elif packet["poses"]:
                    name = packet["poses"][0].ligname
                else:
                    name = "unknown"
                self.pipe.send(
                    (
                        WriteToStorageError(f"Error while writing {name} to the database: {e}"),
                        traceback.format_exc(),
                        f"ligand {name}",
                    )
                )
        if not written and len(self.packets) > 1:
            raise WriteToStorageError(
                "No ligand of the chunk could be written, the database may not be writable."
            )
        return written

    def _log_progress(self):
        current = self.num_files_written + self.counter
        elapsed = time.perf_counter() - self.time0
        rate = current * 60 / (elapsed or 1)
        sys.stdout.write(
            f"\r{current} files processed. {rate:.0f} files/min. "
            f"Elapsed: {elapsed:.0f}s."
        )
        sys.stdout.flush()
