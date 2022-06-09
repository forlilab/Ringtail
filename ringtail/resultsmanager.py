#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail results manager
#

from ringtail import MPManager
from .exceptions import MultiprocessingError, ResultsProcessingError


class ResultsManager:
    def __init__(
        self,
        mode="dlg",
        dbman=None,
        filelist=None,
        chunk_size=1000,
        numclusters=3,
        interaction_tolerance_cutoff=None,
        store_all_poses=False,
        target=None,
    ):
        self.dbman = dbman
        self.filelist = filelist
        self.num_result_files = len(filelist)
        self.parser = MPManager(
            filelist=self.filelist,
            db_obj=self.dbman,
            chunksize=chunk_size,
            mode=mode,
            numclusters=numclusters,
            interaction_tolerance_cutoff=interaction_tolerance_cutoff,
            store_all_poses=store_all_poses,
            target=target,
        )

    def process_results(self):
        try:
            # start MP process
            self.parser.process_files()
        except MultiprocessingError as e:
            raise ResultsProcessingError(
                "Error occurred while processing results"
            ) from e
