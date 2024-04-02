#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail results manager
#

from .mpmanager import MPManager
from .exceptions import ResultsProcessingError
from .storagemanager import StorageManager
from .logmanager import logger


class ResultsManager:
    def __init__(
        self,
        mode: str = None,
        max_poses: int = None,
        interaction_tolerance: float = None,
        store_all_poses: bool = None,
        add_interactions: bool = None,
        interaction_cutoffs: list = None,
        file_sources = None, #TODO needs to accommodate results strings
        string_sources = None,
        max_proc: int = None,
        storageman: StorageManager = None,
        storageman_class = None,
        chunk_size: int = 1, 
        parser_manager: str = "multiprocessing",
    ):
        self.parser_manager = parser_manager
        self.mode = mode
        self.chunk_size = chunk_size
        self.max_poses = max_poses
        self.store_all_poses = store_all_poses
        self.interaction_tolerance = interaction_tolerance
        self.add_interactions = add_interactions
        self.interaction_cutoffs = interaction_cutoffs
        self.target = None
        self.receptor_file = None
        self.file_pattern = None
        self.max_proc = max_proc
        self.storageman_class = storageman_class
        self.storageman = storageman

        self.file_sources = file_sources
        if file_sources is not None:
            self.file_pattern = file_sources.file_pattern
            self.target = file_sources.target
            self.receptor_file = file_sources.receptor_file

        self.string_sources = string_sources

    def process_files(self):
        # check that we have file source(s)
        if (
            self.file_sources.file == (None and [[]])
            and self.file_sources.file_path == (None and [[]])
            and self.file_sources.file_list == (None and [[]])
        ):
            raise ResultsProcessingError(
                "No file sources given. File sources must be given for writing results to database."
            )
        if self.mode == "vina" and self.add_interactions and self.receptor_file is None:
                raise ResultsProcessingError(
                    "Gave --add_interactions with Vina mode but did not specify receptor name. Please give receptor pdbqt name with --receptor_file.")
        
        # start MP process
        logger.debug(f'These are the file options being procesed: {str(self.file_sources.todict())}.')

        # NOTE: if implementing a new parser manager (i.e. serial) must add it to this dict
        implemented_parser_managers = {
            "multiprocessing": MPManager,
        }
        parser_opts = {}
        for k, v in self.__dict__.items():
            if k == "parser_manager":
                continue
            parser_opts[k] = v
        self.parser = implemented_parser_managers[self.parser_manager](**parser_opts)
        self.parser.process_files()

    def process_strings(self):
        # check that we have file source(s)
        if not self.string_sources:
            raise ResultsProcessingError(
                "No string sources given. String sources must be given for writing results to database."
            )
        if self.mode == "vina_string" and self.add_interactions and self.receptor_file is None:
                raise ResultsProcessingError(
                    "Gave 'add_interactions' with Vina mode but did not specify receptor name. Please give receptor pdbqt name with 'receptor_file'.")
        self.mode = "vina_string" # quick and dirty way to reroute mpreaderwriter, change it in here so no negative effects elsewhere
        # start MP process
        logger.debug(f'These are the strings options being procesed: {self.string_sources}.') 
        print(f'\n\n past log msg\n\n')
        #NOTE get this far
        # NOTE: if implementing a new parser manager (i.e. serial) must add it to this dict
        implemented_parser_managers = {
            "multiprocessing": MPManager,
        }
        print(f'\n\n setting up parsed options\n\n')
        parser_opts = {}
        for k, v in self.__dict__.items():
            if k == "parser_manager":
                continue
            parser_opts[k] = v
        self.parser = implemented_parser_managers[self.parser_manager](**parser_opts)
        print(f'\n\n next up is processing strings, and I have a parser: {self.parser}\n\n')
        self.parser.process_strings()

