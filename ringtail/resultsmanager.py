#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail results manager
#

from ringtail import MPManager
from .exceptions import MultiprocessingError, ResultsProcessingError


class ResultsManager():

    def __init__(self,
                 mode='dlg',
                 dbman=None,
                 filelist=None,
                 chunk_size=1000,
                 numclusters=3,
                 no_print_flag=False,
                 target=None):
        self.dbman = dbman
        self.filelist = filelist
        self.num_result_files = len(filelist)
        self.no_print_flag = no_print_flag
        self.parser = MPManager(filelist=self.filelist,
                                db_obj=self.dbman,
                                chunksize=chunk_size,
                                mode=mode,
                                numclusters=numclusters,
                                no_print_flag=self.no_print_flag,
                                target=target)

    def process_results(self):
        try:
            # start MP process
            self.parser.process_files()
        except MultiprocessingError as e:
            raise ResultsProcessingError("Error occurred while processing results") from e
