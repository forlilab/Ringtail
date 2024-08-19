#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail results manager
#

from .mpmanager import MPManager
from .exceptions import ResultsProcessingError
from .storagemanager import StorageManager
from .logutils import LOGGER as logger


class ResultsManager:
    """Class that handles the processing of the results, including passing on the docking results to the appropriate paralell/multi-processing unit

    Args:
        max_poses (int): max number of poses to store for each ligand
        interaction_tolerance (float): Will add the interactions for poses within some tolerance RMSD range of the top pose in a cluster to that top pose."
        store_all_poses (bool): Store all poses from docking results
        add_interactions (bool): find and save interactions between ligand poses and receptor
        interaction_cutoffs (list(float)): cutoff for interactions of hydrogen bonds and VDW interactions, in ångströms
        max_proc (int): Maximum number of processes to create during parallel file parsing.
        storageman (StorageManager): storageman object
        storageman_class (StorageManager): storagemanager child class/database type
        chunk_size (int): how many tasks ot send to a processor at the time
        parser_manager (str, optional): what paralellization or multiprocessing package to use
        file_sources (InputFiles, optional): given file sources including the receptor file
        string_sources (InputStrings, optional): given string sources including the path to the receptor

    Raises:
        ResultsProcessingError
    """

    def __init__(
        self,
        docking_mode: str = None,
        max_poses: int = None,
        interaction_tolerance: float = None,
        store_all_poses: bool = None,
        add_interactions: bool = None,
        interaction_cutoffs: list = None,
        max_proc: int = None,
        storageman: StorageManager = None,
        storageman_class: StorageManager = None,
        chunk_size: int = 1,
        parser_manager: str = "multiprocess",
        file_sources=None,
        string_sources=None,
    ):
        self.parser_manager = parser_manager
        self.docking_mode = docking_mode
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
        # if results are provided as files
        self.file_sources = file_sources
        if file_sources is not None:
            self.file_pattern = file_sources.file_pattern
            self.target = file_sources.target
            self.receptor_file = file_sources.receptor_file
        # if results are provided as strings
        self.string_sources = string_sources
        if self.string_sources is not None:
            self.target = self.string_sources.target
            self.receptor_file = self.string_sources.receptor_file

    def process_docking_data(self):
        """Processes docking data in the form of files or strings

        Raises:
            ResultsProcessingError: if no file or string sources are provided, or if both are provided
        """
        # check that we have results source(s)
        files_sources = bool(self.file_sources)
        if files_sources:
            files_present = not bool(
                self.file_sources.file == (None and [[]])
                and self.file_sources.file_path == (None and [[]])
                and self.file_sources.file_list == (None and [[]])
            )

        strings_present = bool(self.string_sources)

        # if no results are given
        if not files_sources and not strings_present:
            raise ResultsProcessingError(
                "No results sources given. Docking results sources must be given for writing results to database."
            )
        if files_sources is True and files_present is False:
            raise ResultsProcessingError(
                "Indicated docking results files would be processed, but the results file object is empty. Please add results files."
            )
        # if results are given as both types of results
        if files_sources and strings_present:
            raise ResultsProcessingError(
                "Docking results were provided as both file sources and string sources. Currently only one results type is accepeted at the time."
            )

        # start MP process
        if files_sources:
            logmsg = f"These are the file sources being processed: {str(self.file_sources.todict())}"
        else:
            logmsg = f'This is the list of ligands whos strings ware being procssed: {str(self.string_sources.todict()["results_strings"].keys())}'
        logger.debug(logmsg)

        # NOTE: if implementing a new parser manager (i.e. serial) must add it to this dict
        implemented_parser_managers = {
            "multiprocess": MPManager,
        }
        parser_opts = {}
        for k, v in self.__dict__.items():
            if k == "parser_manager":
                continue
            parser_opts[k] = v
        self.parser = implemented_parser_managers[self.parser_manager](**parser_opts)
        self.parser.process_results()
