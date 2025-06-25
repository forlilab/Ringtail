#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail results manager
#

from .mpmanager import MPManager
from .exceptions import ResultsProcessingError
from .storagemanager import StorageManager
from .logutils import LOGGER as logger
from .ringtailoptions import ResultsObject


class ResultsManager:
    """Class that handles the processing of the results, including passing on the docking results to the appropriate paralell/multi-processing unit

    Attributes:
        write_options (dict): incoming results writing options
        parser (MPManager): manager to handle reading and writing docking results

    Raises:
        ResultsProcessingError
    """

    def __init__(
        self,
        db_file: str,
        docking_mode: str,
        storageman_class: StorageManager,
        duplicate_handling: str,
        chunk_size: int = 1,
    ):

        self.writer_options = {
            "db_file": db_file,
            "docking_mode": docking_mode,
            "storageman_class": storageman_class,
            "chunk_size": chunk_size,
            "duplicate_handling": duplicate_handling,
        }

    def process_docking_data(
        self,
        results: ResultsObject,
        processing_options: dict = {},
        parser_manager: str = "multiprocessing",
    ):
        """Processes docking data in the form of files or strings

        Raises:
            ResultsProcessingError: if no file or string sources are provided, or if both are provided
        """

        # if no results are given
        if not results.has_results:
            raise ResultsProcessingError(
                "No results sources given. Docking results sources must be given for writing results to database."
            )
        # if results are given as both types of results
        if results.has_file_results and results.strings:
            raise ResultsProcessingError(
                "Docking results were provided as both file sources and string sources. Currently only one results type is accepeted at the time."
            )

        if results.has_file_results:
            logmsg = f"These are the provided files: {results.file}, directories: {results.file_path}, and file lists: {results.file_list} provided for database storage."
        else:
            logmsg = f"This is the list of ligands whos strings are being procssed: {str(results.strings.keys())}"
        logger.debug(logmsg)

        if processing_options.get("store_all_poses") and processing_options.get(
            "interaction_tolerance"
        ):
            logger.warning(
                "Cannot use 'interaction_tolerance' with 'store_all_poses'. Removing 'interaction_tolerance'."
            )
            processing_options["interaction_tolerance"] = None

        # NOTE: if implementing a new parser manager (i.e. serial) must add it to this dict
        implemented_parser_managers = {
            "multiprocessing": MPManager,
        }
        self.parser: MPManager = implemented_parser_managers[parser_manager](
            results, self.writer_options, processing_options.pop("max_proc")
        )
        # start MP process
        self.parser.process_results(processing_options)
