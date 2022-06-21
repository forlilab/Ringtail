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
        opts={
            "mode": "dlg",
            "filelist": None,
            "chunk_size": 1,
            "max_poses": 3,
            "interaction_tolerance": None,
            "store_all_poses": False,
            "target": None,
            "add_interactions": False,
            "interaction_cutoffs": [3.7, 4.0],
            "receptor_file": None,
        },
        dbman=None,
    ):
        self.dbman = dbman
        self.filelist = opts["filelist"]
        self.parser = MPManager(
            filelist=self.filelist,
            db_obj=self.dbman,
            opts=opts,
        )

    def process_results(self):
        try:
            # start MP process
            self.parser.process_files()
        except MultiprocessingError as e:
            raise ResultsProcessingError(
                "Error occurred while processing results"
            ) from e
