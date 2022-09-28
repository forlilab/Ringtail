#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail results manager
#

from .mpmanager import MPManager
from .exceptions import ResultsProcessingError


class ResultsManager:
    def __init__(
        self,
        dbman,
        parser_manager = "multiprocessing",
        mode="dlg",
        chunk_size=1,
        max_poses=3,
        interaction_tolerance=None,
        store_all_poses=False,
        target=None,
        add_interactions=False,
        interaction_cutoffs=[3.7, 4.0],
        receptor_file=None,
        file_sources={'file': [[]],
                      'file_path': {
                          'path': [[]],
                          'pattern': '*.dlg*',
                          'recursive': None},
                      'file_list': [[]]},
        file_pattern="*.dlg*",
        _stop_at_defaults=False
    ):
        parser_managers = {'multiprocessing': MPManager,}
        
        self.parser_manager = parser_manager
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
        if _stop_at_defaults:
            return

        self.dbman = dbman
        parser_opts = {}
        for k,v in self.__dict__.items():
            if k == "parser_manager":
                continue
            parser_opts[k] = v
        self.parser = parser_managers[self.parser_manager](**parser_opts)

    def process_results(self):
        # check that we have file source(s)
        if self.file_sources["file"] == [[]] and self.file_sources["file_path"]["path"] == [[]] and self.file_sources["file_list"] == [[]]:
            raise ResultsProcessingError("No file sources given. File sources must be given for writing results to database.")
        # start MP process
        self.parser.process_files()

    @classmethod
    def get_defaults(cls):
        return cls(None, _stop_at_defaults=True).__dict__
