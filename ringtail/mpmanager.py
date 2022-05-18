#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail multiprocessing manager
#

import multiprocessing


class MPManager():

    def __init__(self, filelist, mode, db_obj, chunksize, numclusters,
                 no_print_flag, single_receptor):
        # confirm that requested parser mode is implemented
        self.implemented_modes = ["dlg"]
        if mode not in self.implemented_modes:
            raise NotImplementedError(
                "Requested file parsing mode {0} not yet implemented".format(
                    mode))
        self.mode = mode
        self.filelist = filelist
        self.db = db_obj
        self.chunksize = chunksize
        self.numclusters = numclusters
        self.no_print = no_print_flag
        self.single_receptor = single_receptor
        self.num_files = len(self.filelist)

        self.max_proc = multiprocessing.cpu_count()
        self.queueIn = multiprocessing.Queue(maxsize=self.max_proc)
        self.queueOut = multiprocessing.Queue()

    def process_files(self):
        # start the workers in background
        for i in range(self.max_proc):
            # one worker is started for each processor to be used
            s = DockingFileReader(self.queueIn, self.queueOut, self.db,
                                  self.mode, self.numclusters, self.no_print)
            # this method calls .run() internally
            s.start()

        # start the writer to process the data from the workers
        w = Writer(self.queueOut, self.max_proc, self.chunksize, self.db,
                   self.num_files, self.single_receptor)
        w.start()

        # process items in the queue
        for file in self.filelist:

            self.queueIn.put(file, block=True)
        # put as many poison pills in the queue as there are workers
        for i in range(self.max_proc):
            self.queueIn.put(None)

        w.join()
