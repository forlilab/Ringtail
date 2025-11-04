#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail multiprocess manager
#

from time import sleep
import sys
import queue
import fnmatch
import os
import glob
from .mpreaderwriter import DockingFileReader
from .mpreaderwriter import Writer
from .logutils import LOGGER as logger
from .exceptions import MultiprocessingError, RTCoreError
import traceback
from datetime import datetime
import multiprocessing as mp

from .util import iterate_nested
from .ringtailoptions import ResultsObject, validate_file_pattern


def _set_start_method():
    preferred = "fork" if sys.platform != "win32" else "spawn"
    try:
        mp.set_start_method(preferred, force=True)
    except RuntimeError:
        # Already set elsewhere; ignore
        pass


_set_start_method()


class MPManager:
    """Manager that orchestrates paralell processing of docking results data, using one of the supported
    multiprocessors.

    Attributes:
        results (ResultsObject): has all docking results data such as file paths and directories to be written
        num_readers (int): numer of processors
        chunk_size (int): how many tasks ot send to the write processor at the time, essentially how many ligands files to read before writing them per processor
        num_files (int): number of files processed at any given time
        queueIn (mp.Queue): queue for files about to be processed by file reader
        queueOut (mp.Queue): queue for files that have been processed and should be written
        writer_options (dict): incoming results writing options
        shared (dict): limited write_options for the docking file readers
        file_pattern (str): file extension for selected docking_mode
        p_conn (one side of mp.Pipe): parent connection for piping information between processes
        c_conn (one side of mp.Pipe): child connection for piping information between processes
        receptor_file: path to receptor file
        workers (list): reader and writer processes

    """

    def __init__(
        self,
        results: ResultsObject,  # also has receptor information
        writer_options: dict,  # has mostly settings for Writer
        max_proc: int,
    ):
        self.results = results
        # initialize multiprocessing variables
        # number of processors
        if max_proc is None:
            max_proc = mp.cpu_count()
        self.num_readers = max_proc - 1
        # shared queue
        self.queueIn = mp.Queue(maxsize=2 * max_proc)
        self.queueOut = mp.Queue(maxsize=2 * max_proc)
        # Start reader processes
        self.workers = []
        self.p_conn, self.c_conn = mp.Pipe(duplex=True)
        self.num_files = 0

        self.writer_options = writer_options
        self.shared = {}
        # populate read only dict
        self.shared["docking_mode"] = self.writer_options.pop("docking_mode")
        self.shared["receptor_string"] = results.receptor_string
        self.shared["target"] = results.target_name

        # these two variables are needed for queueing files
        self.receptor_file = results.receptor_file_path
        self.file_pattern = validate_file_pattern(
            self.shared["docking_mode"], results.file_pattern
        )

    def process_results(self, reader_options: dict):
        logger.info(f"Starting {self.num_readers} docking results readers")
        self.shared.update(reader_options)

        for _ in range(self.num_readers):
            reader = DockingFileReader(
                self.queueIn,
                self.queueOut,
                self.c_conn,
                self.shared,
            )
            reader.start()
            self.workers.append(reader)
        # Start writer in a thread (pulls from queueOut)
        writer = Writer(self.queueOut, self.num_readers, self.writer_options)
        writer.start()
        self.workers.append(writer)
        try:
            self._process_data_sources()
        except Exception as e:
            tb = traceback.format_exc()
            self._kill_all_workers(e, "results sources processing", tb)

        # Now safe to shut down readers
        for _ in range(self.num_readers):
            self.queueIn.put(None)

        # check for exceptions
        while writer.is_alive():
            sleep(0.5)
            self._check_for_worker_exceptions()

        writer.join()

        logger.info(f"Wrote {self.num_files} docking results to the database")

    def _process_data_sources(self):
        """Adds each docking result item to the queue, including files and data provided as string/dict.
        For files, processes lists of files, recursively traveresed filepaths, and individually listed file paths.
        """

        if self.results.has_file_results:
            # add individual file(s)
            files = list(iterate_nested(self.results.file))
            if files:
                for file in files:
                    if (
                        fnmatch.fnmatch(file, self.file_pattern)
                        and file != self.receptor_file
                    ):
                        self._add_to_queue(file)

            # add files from file path(s)
            path_list = list(iterate_nested(self.results.file_path))
            if path_list:
                for path in path_list:
                    # scan for ligand dlgs
                    for files in self._scan_dir(
                        path, self.file_pattern, self.results.recursive_path_traverse
                    ):
                        for file in files:
                            self._add_to_queue(file)

            # add files from file list(s)
            file_lists = list(iterate_nested(self.results.file_list))
            if file_lists:
                for file_list in file_lists:
                    self._scan_file_list(file_list, self.file_pattern.replace("*", ""))

        # add docking data from input strings
        if self.results.strings:
            try:
                for (
                    ligand_name,
                    docking_result,
                ) in self.results.strings.items():
                    string_data = {ligand_name: docking_result}
                    self._add_to_queue(string_data)
            except:
                raise RTCoreError(
                    "There was an error while reading the results string input."
                )

    def _add_to_queue(self, results_data):
        """_summary_

        Args:
            results_data (string or dict): results data provided as a file path or a dictionary kw pair

        Raises:
            MultiprocessingError
        """
        # adds result file to the multiprocess queue
        max_attempts = 750
        timeout = 0.5  # seconds
        if not isinstance(results_data, dict) and self.receptor_file is not None:
            if (
                os.path.split(results_data)[-1] == os.path.split(self.receptor_file)[-1]
            ):  # check that we don't try to add the receptor
                return
        attempts = 0
        while True:
            if attempts >= max_attempts:
                raise MultiprocessingError(
                    "Something is blocking the progressing of results data reading. Exiting program."
                ) from queue.Full
            try:
                self.queueIn.put(results_data, block=True, timeout=timeout)
                self.num_files += 1
                self._check_for_worker_exceptions()
                break
            except queue.Full:
                attempts += 1
                self._check_for_worker_exceptions()

    def _check_for_worker_exceptions(self):
        if self.p_conn.poll():
            error, tb, filename = self.p_conn.recv()
            logger.error(f"Caught error in multiprocess from {filename}:")
            logger.error(f"{tb}")
            # don't kill parser errors, only database error
            if filename == "Database":
                self._kill_all_workers(error, filename, tb)
            else:
                with open("ringtail_failed_files.log", "a") as f:
                    f.write(
                        str(datetime.now()) + f"\tRingtail failed to parse {filename}\n"
                    )
                    self.num_files -= 1
                    logger.debug(tb)

    def _kill_all_workers(self, error, filename, tb):
        for s in self.workers:
            s.kill()
        logger.debug(f"Error encountered while handling {filename}")
        logger.debug(tb)
        raise error

    def _scan_dir(self, path, pattern, recursive=False):
        """scan for valid output files in a directory the pattern is used
        to glob files optionally, a recursive search is performed

        Args:
            path (str): folder path
            pattern (str): file extension
            recursive (bool, optional): look for files and folders recursively

        Yields:
            list: of file paths found in the search
        """
        logger.info(
            "Scanning directory [%s] for files (pattern:|%s|)" % (path, pattern)
        )
        if recursive:
            path = os.path.normpath(path)
            path = os.path.expanduser(path)
            for dirpath, _, filenames in os.walk(path):
                yield (  # <----
                    os.path.join(dirpath, f)
                    for f in fnmatch.filter(filenames, "*" + pattern)
                )
        else:
            yield glob.glob(os.path.join(path, pattern))  # <----

    def _scan_file_list(self, filename, pattern):
        """read file names from file list and ensures they exist,
        then adding them to the list of files to be processed

        Args:
            filename (str): filename provided in list
            pattern (str): file extension

        Raises:
            MultiprocessingError
        """

        lig_accepted = []
        c = 0
        with open(filename, "r") as fp:
            for line in fp.readlines():
                line = line.strip()
                c += 1
                if os.path.isfile(line):
                    if line.endswith(pattern) or line.endswith(
                        pattern + ".gz"
                    ):  # NOTE if adding zip option change here
                        lig_accepted.append(line)
                else:
                    logger.warning("Warning! file |%s| does not exist" % line)
        if len(lig_accepted) == 0:
            raise MultiprocessingError(
                "*ERROR* No valid files were found when reading from |%s|" % filename
            )
        logger.info(
            "# [ %5.3f%% files in list accepted (%d) ]"
            % ((len(lig_accepted) / c * 100, c))
        )

        for file in lig_accepted:
            if file != self.receptor_file:
                self._add_to_queue(file)
