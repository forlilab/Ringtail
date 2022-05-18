#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail results manager
#

from ringtail import MPManager


class ResultsManager():

    def __init__(self,
                 mode='dlg',
                 dbman=None,
                 filelist=None,
                 chunk_size=1000,
                 numclusters=3,
                 no_print_flag=False,
                 single_receptor=False):
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
                                single_receptor=single_receptor)

    def process_results(self):
        # start MP process
        self.parser.process_files()
