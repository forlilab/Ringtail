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
from .storagemanager import StorageManager, StorageManagerSQLite
from .mpreaderwriter import DockingFileReader
from .mpreaderwriter import Writer
from .exceptions import MultiprocessingError
import traceback
from datetime import datetime

os_string = platform.system()
if os_string == "Darwin":  # mac
    import multiprocess as multiprocessing
else:
    import multiprocessing


class MPManager:
    def __init__(
        self,
        storageman,
        storageman_class=StorageManagerSQLite,
        mode="dlg",
        chunk_size=1,
        max_poses=3,
        interaction_tolerance=None,
        store_all_poses=False,
        target=None,
        add_interactions=False,
        interaction_cutoffs=[3.7, 4.0],
        receptor_file=None,
        file_sources={
            "file": [[]],
            "file_path": {"path": [[]], "pattern": "*.dlg*", "recursive": None},
            "file_list": [[]],
        },
        file_pattern="*.dlg*",
        max_proc=None,
    ):

        # confirm that requested parser mode is implemented
        self.implemented_modes = ["dlg", "vina"]
        if mode not in self.implemented_modes:
            raise NotImplementedError(
                f"Requested file parsing mode {mode} not yet implemented"
            )
        self.mode = mode
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

        self.storageman = storageman
        self.storageman_class = storageman_class
        self.num_files = 0
        self.max_proc = max_proc

    def process_files(self):
        if self.max_proc is None:
            self.max_proc = multiprocessing.cpu_count()
        self.num_readers = self.max_proc - 1
        self.queueIn = multiprocessing.Queue(maxsize=2 * self.max_proc)
        self.queueOut = multiprocessing.Queue(maxsize=2 * self.max_proc)
        # start the workers in background
        self.workers = []
        self.p_conn, self.c_conn = multiprocessing.Pipe(True)
        logging.info("Starting {0} file readers".format(self.num_readers))
        for i in range(self.num_readers):
            # one worker is started for each processor to be used
            s = DockingFileReader(
                self.queueIn,
                self.queueOut,
                self.c_conn,
                self.storageman,
                self.storageman_class,
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
            self.num_readers,
            self.c_conn,
            self.chunk_size,
            self.storageman,
            self.mode,
        )

        w.start()
        self.workers.append(w)

        # process items in the queue
        try:
            self._process_sources()
        except Exception as e:
            tb = traceback.format_exc()
            self._kill_all_workers(e, "file sources", tb)
        # put as many poison pills in the queue as there are workers
        for i in range(self.num_readers):
            self.queueIn.put(None)

        # check for exceptions
        while w.is_alive():
            sleep(0.5)
            self._check_for_worker_exceptions()

        w.join()

        logging.info("Wrote {0} files to database".format(self.num_files))

    def _process_sources(self):
        # add individual file(s)
        if self.file_sources["file"] != [[]]:
            for file_list in self.file_sources["file"]:
                for file in file_list:
                    if (
                        fnmatch.fnmatch(file, self.file_pattern)
                        and file != self.receptor_file
                    ):
                        self._add_to_queue(file)
        # add files from file path(s)
        if self.file_sources["file_path"]["path"] != [[]]:
            for path_list in self.file_sources["file_path"]["path"]:
                for path in path_list:
                    # scan for ligand dlgs
                    for files in self._scan_dir(
                        path, self.file_pattern, recursive=True
                    ):
                        for f in files:
                            self._add_to_queue(f)
        # add files from file list(s)
        if self.file_sources["file_list"] != [[]]:
            for filelist_list in self.file_sources["file_list"]:
                for filelist in filelist_list:
                    self._scan_file_list(filelist, self.file_pattern.replace("*", ""))

    def _add_to_queue(self, file):
        max_attempts = 750
        timeout = 0.5  # seconds
        if self.receptor_file is not None:
            if (
                os.path.split(file)[-1] == os.path.split(self.receptor_file)[-1]
            ):  # check that we don't try to add the receptor
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
                # logging.debug(f"Queue full: queueIn.put attempt {attempts} timed out. {max_attempts - attempts} put attempts remaining.")
                attempts += 1
                self._check_for_worker_exceptions()

    def _check_for_worker_exceptions(self):
        if self.p_conn.poll():
            error, tb, filename = self.p_conn.recv()
            logging.error(f"Caught error in multiprocessing from {filename}")
            # don't kill parser errors, only database error
            if filename == "Database":
                self._kill_all_workers(error, filename, tb)
            else:
                with open("ringtail_failed_files.log", 'a') as f:
                    f.write(str(datetime.now()) + f"\tRingtail failed to parse {filename}\n")
                    logging.debug(tb)

    def _kill_all_workers(self, error, filename, tb):
        for s in self.workers:
            s.kill()
        logging.debug(f"Error encountered while handling {filename}")
        logging.debug(tb)
        raise error

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
                    logging.warning("Warning! file |%s| does not exist" % line)
        if len(lig_accepted) == 0:
            raise MultiprocessingError(
                "*ERROR* No valid files were found when reading from |%s|" % filename
            )
        logging.info(
            "# [ %5.3f%% files in list accepted (%d) ]"
            % ((len(lig_accepted) / c * 100, c))
        )

        for file in lig_accepted:
            if file != self.receptor_file:
                self._add_to_queue(file)
