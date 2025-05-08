#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail multiprocess manager
#

from time import sleep
import queue
import fnmatch
import os
import glob
from .mpreaderwriter import DockingFileReader
from .mpreaderwriter import Writer
from .logutils import LOGGER as logger
from .exceptions import MultiprocessingError, RTCoreError, ResultsProcessingError
import traceback
from datetime import datetime
import pickle
import multiprocessing as mp
import threading

mp.set_start_method("spawn", force=True)


class MPManager:
    """Manager that orchestrates paralell processing of docking results data, using one of the supported
    multiprocessors.

    Attributes:
        docking_mode (str): describes what docking engine was used to produce the results
        max_poses (int): max number of poses to store for each ligand
        interaction_tolerance (float): Will add the interactions for poses within some tolerance RMSD range of the top pose in a cluster to that top pose."
        store_all_poses (bool): Store all poses from docking results
        add_interactions (bool): find and save interactions between ligand poses and receptor
        interaction_cutoffs (list(float)): cutoff for interactions of hydrogen bonds and VDW interactions, in ångströms
        max_proc (int): Maximum number of processes to create during parallel file parsing.
        storageman (StorageManager): storageman object
        chunk_size (int): how many tasks ot send to a processor at the time
        target (str): name of receptor
        receptor_file (str): file path to receptor
        file_pattern (str, optional): file pattern to look for if recursively finding results files to process
        file_sources (InputFiles, optional): RingtailOption object that holds all attributes related to results files
        string_sources (InputStrings, optional): RingtailOption object that holds all attributes related to results strings
        num_files (int): number of files processed at any given time

    """

    def __init__(
        self,
        docking_mode,
        max_poses,
        interaction_tolerance,
        store_all_poses,
        add_interactions,
        interaction_cutoffs,
        max_proc,
        storageman,
        chunk_size,
        target,
        receptor_file,
        file_pattern=None,
        file_sources=None,
        string_sources=None,
    ):
        self.storageman = storageman
        self.chunk_size = chunk_size
        self.docking_mode = docking_mode

        # shared memory
        manager = mp.Manager()
        self.shared = manager.dict()
        # general docking info
        self.shared["docking_mode"] = docking_mode
        self.shared["max_poses"] = max_poses
        self.shared["store_all_poses"] = store_all_poses
        self.shared["format_method"] = self.storageman.format_for_storage

        # interaction info
        self.shared["add_interactions"] = add_interactions
        if add_interactions:
            self.shared["receptor_string"] = self._get_receptor_string()
        self.shared["interaction_tolerance"] = interaction_tolerance
        self.shared["interaction_cutoffs"] = interaction_cutoffs
        self.shared["target"] = target
        self.receptor_file = receptor_file
        self.file_sources = file_sources
        self.file_pattern = file_pattern
        self.string_sources = string_sources
        self.num_files = 0
        self.max_proc = max_proc

    def process_results(self):
        if self.max_proc is None:
            self.max_proc = mp.cpu_count()
        self.num_readers = self.max_proc - 1
        self.queueIn = mp.Queue(maxsize=2 * self.max_proc)
        self.queueOut = mp.Queue(maxsize=2 * self.max_proc)

        logger.info(f"Starting {self.num_readers} docking results readers")

        # Start reader processes
        self.workers = []
        self.p_conn, self.c_conn = mp.Pipe(duplex=True)
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
        dbfile = self.storageman.db_file
        writer = Writer(
            self.queueOut,
            self.num_readers,
            self.chunk_size,
            dbfile,
            self.docking_mode,
        )
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

    def _get_receptor_string(self) -> str:
        try:
            from .receptormanager import ReceptorManager as rm

            with self.storageman:
                # grab receptor info from database, this assumes there is only one receptor in the database
                receptor_blob = self.storageman.fetch_receptor_objects()[0][
                    1
                ]  # method returns an iter of tuples, blob is the second tuple element in the first list element
                # convert receptor blob to string
                return rm.blob2str(receptor_blob)
        except:
            raise ResultsProcessingError(
                "add_interactions was requested, but cannot find the receptor in the database. Please ensure to include the receptor_file and save_receptor if the receptor has not already been added to the database."
            )

    def _process_data_sources(self):
        """Adds each docking result item to the queue, including files and data provided as string/dict.
        For files, processes lists of files, recursively traveresed filepaths, and individually listed file paths.
        """

        def _iterate_nested(obj):
            """
            File inputs can come in multiple levels of nested lists, this method unpacks them

            Args:
                obj (list[list[list[etc]]]): None or nested lists

            Returns:
                None: if input is None

            Yields:
                str: should be unpacked paths to docking results
            """
            if obj is None:
                return None
            elif isinstance(obj, list):
                for item in obj:
                    yield from _iterate_nested(item)
            else:
                yield obj

        if self.file_sources:
            # add individual file(s)
            files = list(_iterate_nested(self.file_sources.file))
            if files:
                for file in files:
                    if (
                        fnmatch.fnmatch(file, self.file_pattern)
                        and file != self.receptor_file
                    ):
                        self._add_to_queue(file)

            # add files from file path(s)
            path_list = list(_iterate_nested(self.file_sources.file_path))
            if path_list:
                for path in path_list:
                    # scan for ligand dlgs
                    for files in self._scan_dir(
                        path, self.file_pattern, recursive=True
                    ):
                        for file in files:
                            self._add_to_queue(file)

            # add files from file list(s)
            file_lists = list(_iterate_nested(self.file_sources.file_list))
            if file_lists:
                for file_list in file_lists:
                    self._scan_file_list(file_list, self.file_pattern.replace("*", ""))

        # add docking data from input strings
        if self.string_sources:
            try:
                for (
                    ligand_name,
                    docking_result,
                ) in self.string_sources.results_strings.items():
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
            string (bool, optional): switch if results provided as a string

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
            logger.error(f"Caught error in multiprocess from {filename}")
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
            yield glob(os.path.join(path, pattern))  # <----

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
