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
        dbman=None,
    ):
        self.dbman = dbman
        self.parser = MPManager(
            db_obj=self.dbman,
            opts=opts,
        )

    def process_results(self):
        # start MP process
        self.parser.process_files()
