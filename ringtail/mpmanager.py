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
from .logutils import LOGGER
from .exceptions import MultiprocessingError, RTCoreError
import traceback
from datetime import datetime

import multiprocess


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
        storageman_class (StorageManager): storagemanager child class/database type
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
        storageman_class,
        chunk_size,
        target,
        receptor_file,
        file_pattern=None,
        file_sources=None,
        string_sources=None,
    ):
        self.docking_mode = docking_mode
        self.chunk_size = chunk_size
        self.max_poses = max_poses
        self.store_all_poses = store_all_poses
        self.interaction_tolerance = interaction_tolerance
        self.target = target
        self.add_interactions = add_interactions
        self.interaction_cutoffs = interaction_cutoffs
        self.receptor_file = receptor_file
        self.file_sources = file_sources
        self.file_pattern = file_pattern
        self.string_sources = string_sources
        self.storageman = storageman
        self.storageman_class = storageman_class
        self.num_files = 0
        self.max_proc = max_proc
        self.logger = LOGGER

    def process_results(self):
        """Processes results data (files or string sources) by adding them to the queue
        and starting their processing in multiprocess.
        """
        if self.max_proc is None:
            self.max_proc = multiprocess.cpu_count()
        self.num_readers = self.max_proc - 1
        self.queueIn = multiprocess.Queue(maxsize=2 * self.max_proc)
        self.queueOut = multiprocess.Queue(maxsize=2 * self.max_proc)
        # start the workers in background
        self.workers = []
        self.p_conn, self.c_conn = multiprocess.Pipe(True)
        self.logger.info(
            "Starting {0} docking results readers".format(self.num_readers)
        )
        for i in range(self.num_readers):
            # one worker is started for each processor to be used
            s = DockingFileReader(
                self.queueIn,
                self.queueOut,
                self.c_conn,
                self.storageman,
                self.storageman_class,
                self.docking_mode,
                self.max_poses,
                self.interaction_tolerance,
                self.store_all_poses,
                self.target,
                self.add_interactions,
                self.interaction_cutoffs,
                self.receptor_file,
            )
            # this method calls .run() internally
            s.start()
            self.workers.append(s)

        # start the writer to process the data from the workers
        w = Writer(
            self.queueOut,
            self.num_readers,
            self.c_conn,
            self.chunk_size,
            self.storageman,
            self.docking_mode,
        )

        w.start()
        self.workers.append(w)

        # process items in the queue
        try:
            self._process_data_sources()
        except Exception as e:
            tb = traceback.format_exc()
            self._kill_all_workers(e, "results sources processing", tb)
        # put as many poison pills in the queue as there are workers
        for i in range(self.num_readers):
            self.queueIn.put(None)

        # check for exceptions
        while w.is_alive():
            sleep(0.5)
            self._check_for_worker_exceptions()

        w.join()

        self.logger.info(
            "Wrote {0} docking results to the database".format(self.num_files)
        )

    def _process_data_sources(self):
        """Adds each docking result item to the queue, including files and data provided as string/dict.
        For files, processes lists of files, recursively traveresed filepaths, and individually listed file paths.
        """

        if self.file_sources:
            # add individual file(s)
            if self.file_sources.file != (None and [[]]):
                for file_list in self.file_sources.file:
                    for file in file_list:
                        if (
                            fnmatch.fnmatch(file, self.file_pattern)
                            and file != self.receptor_file
                        ):
                            self._add_to_queue(file)

            # add files from file path(s)
            if self.file_sources.file_path != (None and [[]]):
                for path_list in self.file_sources.file_path:
                    for path in path_list:
                        # scan for ligand dlgs
                        for files in self._scan_dir(
                            path, self.file_pattern, recursive=True
                        ):
                            for file in files:
                                self._add_to_queue(file)

            # add files from file list(s)
            if self.file_sources.file_list != (None and [[]]):
                for filelist_list in self.file_sources.file_list:
                    for filelist in filelist_list:
                        self._scan_file_list(
                            filelist, self.file_pattern.replace("*", "")
                        )

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
            self.logger.error(f"Caught error in multiprocess from {filename}")
            # don't kill parser errors, only database error
            if filename == "Database":
                self._kill_all_workers(error, filename, tb)
            else:
                with open("ringtail_failed_files.log", "a") as f:
                    f.write(
                        str(datetime.now()) + f"\tRingtail failed to parse {filename}\n"
                    )
                    self.logger.debug(tb)

    def _kill_all_workers(self, error, filename, tb):
        for s in self.workers:
            s.kill()
        self.logger.debug(f"Error encountered while handling {filename}")
        self.logger.debug(tb)
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
        self.logger.info(
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
                    self.logger.warning("Warning! file |%s| does not exist" % line)
        if len(lig_accepted) == 0:
            raise MultiprocessingError(
                "*ERROR* No valid files were found when reading from |%s|" % filename
            )
        self.logger.info(
            "# [ %5.3f%% files in list accepted (%d) ]"
            % ((len(lig_accepted) / c * 100, c))
        )

        for file in lig_accepted:
            if file != self.receptor_file:
                self._add_to_queue(file)
