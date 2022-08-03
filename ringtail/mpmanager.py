#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail multiprocessing manager
#

import platform
from time import sleep
import logging
import queue
import fnmatch
import os
import glob
import warnings
from .mpreaderwriter import DockingFileReader
from .mpreaderwriter import Writer
from .exceptions import MultiprocessingError

os_string = platform.system()
if os_string == "Darwin":  # mac
    import multiprocess as multiprocessing
else:
    import multiprocessing


class MPManager:
    def __init__(
        self,
        db_obj,
        opts={
            "mode": "dlg",
            "chunk_size": 1,
            "max_poses": 3,
            "interaction_tolerance": None,
            "store_all_poses": False,
            "target": None,
            "add_interactions": False,
            "interaction_cutoffs": [3.7, 4.0],
            "receptor_file": None,
            "file_sources": None,
            "file_pattern": None,
        },
    ):
        # confirm that requested parser mode is implemented
        self.implemented_modes = ["dlg", "vina"]
        if opts["mode"] not in self.implemented_modes:
            raise NotImplementedError(
                "Requested file parsing mode {0} not yet implemented".format(
                    opts["mode"]
                )
            )
        self.mode = opts["mode"]
        self.db = db_obj
        self.chunksize = opts["chunk_size"]
        self.max_poses = opts["max_poses"]
        self.store_all_poses = opts["store_all_poses"]
        self.interaction_tolerance = opts["interaction_tolerance"]
        self.target = opts["target"]
        self.add_interactions = opts["add_interactions"]
        self.interaction_cutoffs = opts["interaction_cutoffs"]
        self.receptor_file = opts["receptor_file"]
        self.file_sources = opts["file_sources"]
        self.file_pattern = opts["file_pattern"]
        self.num_files = 0

        self.max_proc = multiprocessing.cpu_count()
        self.queueIn = multiprocessing.Queue(maxsize=2 * self.max_proc)
        self.queueOut = multiprocessing.Queue()

    def process_files(self):
        # start the workers in background
        self.workers = []
        self.p_conn, self.c_conn = multiprocessing.Pipe(True)
        logging.info("Starting {0} file readers".format(self.max_proc))
        for i in range(self.max_proc):
            # one worker is started for each processor to be used
            s = DockingFileReader(
                self.queueIn,
                self.queueOut,
                self.c_conn,
                self.db,
                self.mode,
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
            self.max_proc,
            self.c_conn,
            self.chunksize,
            self.db,
            self.mode,
        )

        w.start()
        self.workers.append(w)

        # process items in the queue
        self._process_sources()
        # put as many poison pills in the queue as there are workers
        for i in range(self.max_proc):
            self.queueIn.put(None)

        # check for exceptions
        while w.is_alive():
            sleep(0.5)
            self._check_for_worker_exceptions()

        w.join()

        logging.info("Wrote {0} files to database".format(self.num_files))

    def _process_sources(self):
        # add individual files
        if self.file_sources["file"] is not None:
            for file_list in self.file_sources["file"]:
                for file in file_list:
                    if fnmatch.fnmatch(file, self.file_pattern):
                        self._add_to_queue(file)

        if self.file_sources["file_path"] is not None:
            for path_list in self.file_sources["file_path"]["path"]:
                for path in path_list:
                    # scan for ligand dlgs
                    for files in self._scan_dir(
                        path, self.file_pattern, recursive=True
                    ):
                        for f in files:
                            self._add_to_queue(f)

        if self.file_sources["file_list"] is not None:
            for filelist_list in self.file_sources["file_list"]:
                for filelist in filelist_list:
                    self._scan_file_list(filelist, self.file_pattern.replace("*", ""))

    def _add_to_queue(self, file):
        max_attempts = 750
        timeout = 0.5  # seconds
        if file == self.receptor_file:
            return
        attempts = 0
        while True:
            if attempts >= max_attempts:
                raise MultiprocessingError(
                    "Something is blocking the progressing of file reading. Exiting program."
                ) from queue.Full
            try:
                self.queueIn.put(file, block=True, timeout=timeout)
                self.num_files += 1
                self._check_for_worker_exceptions()
                break
            except queue.Full:
                attempts += 1
                self._check_for_worker_exceptions()

    def _check_for_worker_exceptions(self):
        if self.p_conn.poll():
            logging.debug("Caught error in multiprocessing")
            error, tb, file_name = self.p_conn.recv()
            for s in self.workers:
                s.kill()
            logging.debug(f"Error encountered while parsing {file_name}")
            logging.debug(tb)
            raise MultiprocessingError("Error occurred during file parsing!")

    def _scan_dir(self, path, pattern, recursive=False):
        """scan for valid output files in a directory
        the pattern is used to glob files
        optionally, a recursive search is performed
        """
        logging.info(
            "-Scanning directory [%s] for files (pattern:|%s|)" % (path, pattern)
        )
        if recursive:
            path = os.path.normpath(path)
            path = os.path.expanduser(path)
            for dirpath, dirnames, filenames in os.walk(path):
                yield (  # <----
                    os.path.join(dirpath, f)
                    for f in fnmatch.filter(filenames, "*" + pattern)
                )
        else:
            yield glob(os.path.join(path, pattern))  # <----

    def _scan_file_list(self, filename, pattern=".dlg"):
        """read file names from file list"""
        lig_accepted = []
        c = 0
        with open(filename, "r") as fp:
            for line in fp.readlines():
                line = line.strip()
                c += 1
                if os.path.isfile(line):
                    if line.endswith(pattern) or line.endswith(pattern + ".gz"):
                        lig_accepted.append(line)
                else:
                    warnings.warn("Warning! file |%s| does not exist" % line)
        if len(lig_accepted) == 0:
            raise MultiprocessingError(
                "*ERROR* No valid files were found when reading from |%s|" % filename
            )
        logging.info(
            "# [ %5.3f%% files in list accepted (%d) ]"
            % ((len(lig_accepted) / c * 100, c))
        )

        for file in lig_accepted:
            self._add_to_queue(file)
