#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail multiprocessing manager
#

import multiprocessing
from time import sleep
import logging
from .mpreaderwriter import DockingFileReader
from .mpreaderwriter import Writer
from .exceptions import MultiprocessingError


class MPManager:
    def __init__(
        self,
        filelist,
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
        self.filelist = filelist
        self.db = db_obj
        self.chunksize = opts["chunk_size"]
        self.max_poses = opts["max_poses"]
        self.store_all_poses = opts["store_all_poses"]
        self.interaction_tolerance = opts["interaction_tolerance"]
        self.target = opts["target"]
        self.add_interactions = opts["add_interactions"]
        self.interaction_cutoffs = opts["interaction_cutoffs"]
        self.receptor_file = opts["receptor_file"]

        self.num_files = len(self.filelist)

        self.max_proc = multiprocessing.cpu_count()
        self.queueIn = multiprocessing.Queue(maxsize=self.max_proc)
        self.queueOut = multiprocessing.Queue()

    def process_files(self):
        # start the workers in background
        workers = []
        p_conn, c_conn = multiprocessing.Pipe(False)
        for i in range(self.max_proc):
            # one worker is started for each processor to be used
            s = DockingFileReader(
                self.queueIn,
                self.queueOut,
                c_conn,
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
            workers.append(s)

        # start the writer to process the data from the workers
        w = Writer(
            self.queueOut,
            self.max_proc,
            c_conn,
            self.chunksize,
            self.db,
            self.num_files,
            self.mode,
        )

        w.start()
        workers.append(w)

        # process items in the queue
        for file in self.filelist:
            self.queueIn.put(file, block=True)
        # put as many poison pills in the queue as there are workers
        for i in range(self.max_proc):
            self.queueIn.put(None)

        # check for exceptions
        while w.is_alive():
            sleep(0.5)
            if p_conn.poll():
                logging.debug("Caught error in multiprocessing")
                error, tb = p_conn.recv()
                for s in workers:
                    s.kill()
                logging.debug(tb)
                raise MultiprocessingError("Error occurred during file parsing!")

        w.join()
