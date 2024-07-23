#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail multiprocessing manager
#

import platform
from time import sleep
import queue
import fnmatch
import os
import glob

# from .mpreaderwriter import DockingFileReader
# from .mpreaderwriter import Writer
from .logutils import LOGGER
from .exceptions import MultiprocessingError, RTCoreError
import traceback
from datetime import datetime
import multiprocessing


class MPManager:
    """Manager that orchestrates paralell processing of docking results data, using one of the supported
    multiprocessors.

    Attributes:
        docking_mode (str): describes what docking engine was used to produce the results
        max_poses (int): max number of poses to store for each ligand
        interaction_tolerance (float): Will add the interactions for poses within some tolerance RMSD range of the top pose in a cluster to that top pose."
        store_all_poses (bool): Store all poses from docking results
        add_interactions (bool): find and save interactions between ligand poses and receptor
        interaction_cutoffs (list(float)): cutoff for interactions of hydrogen bonds and VDW interactions, in ångströms
        max_proc (int): Maximum number of processes to create during parallel file parsing.
        storageman (StorageManager): storageman object
        storageman_class (StorageManager): storagemanager child class/database type
        chunk_size (int): how many tasks ot send to a processor at the time
        target (str): name of receptor
        receptor_file (str): file path to receptor
        file_pattern (str, optional): file pattern to look for if recursively finding results files to process
        file_sources (InputFiles, optional): RingtailOption object that holds all attributes related to results files
        string_sources (InputStrings, optional): RingtailOption object that holds all attributes related to results strings
        num_files (int): number of files processed at any given time

    """

    def __init__(
        self,
        docking_mode,
        max_poses,
        interaction_tolerance,
        store_all_poses,
        add_interactions,
        interaction_cutoffs,
        max_proc,
        db_file,
        storageman_class,
        chunk_size,
        target,
        receptor_file,
        file_pattern=None,
        file_sources=None,
        string_sources=None,
    ):
        self.docking_mode = docking_mode
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
        self.string_sources = string_sources
        self.db_file = db_file
        self.storageman_class = storageman_class
        self.num_files = 0
        self.max_proc = max_proc
        self.logger = LOGGER

        # set up queue
        if self.max_proc is None:
            self.max_proc = multiprocessing.cpu_count()
        self.num_readers = self.max_proc - 1
        # start the workers in background
        self.workers = []
        self.p_conn, self.c_conn = multiprocessing.Pipe(True)
        self.logger.info(
            "Starting {0} docking results readers".format(self.num_readers)
        )
        self.managed_queue_in = multiprocessing.Manager().Queue(
            maxsize=2 * self.max_proc
        )
        self.managed_queue_out = multiprocessing.Manager().Queue(
            maxsize=2 * self.max_proc
        )

    def process_results(self, string_processing=False):
        """Processes results data (files or string sources) by adding them to the queue
        and starting their processing in multiprocessing.

        Args:
            string_sources (bool, optional): Switch for processing results that are provided as strings instead of files.
        """

        # exception comment: "Queue objects should only be shared between processes through inheritance"
        for _ in range(self.num_readers):
            # one worker is started for each processor to be used
            s = DockingFileReader(
                self.managed_queue_in,
                self.managed_queue_out,
                self.c_conn,
                self.db_file,
                self.storageman_class,
                self.docking_mode,
                self.max_poses,
                self.interaction_tolerance,
                self.store_all_poses,
                self.target,
                self.add_interactions,
                self.interaction_cutoffs,
                self.receptor_file,
                string_processing,
            )
            print(s.__dict__)
            # this method calls .run() internally
            s.start()
            self.workers.append(s)

        # start the writer to process the data from the workers
        w = Writer(
            self.managed_queue_out,
            self.num_readers,
            self.c_conn,
            self.chunk_size,
            self.docking_mode,
            self.db_file,
            self.storageman_class,
        )

        w.start()
        self.workers.append(w)

        # process items in the queue
        try:
            if string_processing == True:
                self._process_string_sources()
            else:
                self._process_file_sources()
        except Exception as e:
            tb = traceback.format_exc()
            self._kill_all_workers(e, "results sources processing", tb)
        # put as many poison pills in the queue as there are workers
        for i in range(self.num_readers):
            self.managed_queue_in.put(None)

        # check for exceptions
        while w.is_alive():
            sleep(0.5)
            self._check_for_worker_exceptions()

        w.join()

        self.logger.info(
            "Wrote {0} docking results to the database".format(self.num_files)
        )

    def _process_file_sources(self):
        """Adds each results file item to the queue, and processes lists of files,
        recursively traveresed filepaths, and individually listed file paths.
        """
        # add individual file(s)
        if self.file_sources.file != (None and [[]]):
            for file_list in self.file_sources.file:
                for file in file_list:
                    if (
                        fnmatch.fnmatch(file, self.file_pattern)
                        and file != self.receptor_file
                    ):
                        self._add_to_queue(file)

        # add files from file path(s)
        if self.file_sources.file_path != (None and [[]]):
            for path_list in self.file_sources.file_path:
                for path in path_list:
                    # scan for ligand dlgs
                    for files in self._scan_dir(
                        path, self.file_pattern, recursive=True
                    ):
                        for f in files:
                            self._add_to_queue(f)
        # add files from file list(s)
        if self.file_sources.file_list != (None and [[]]):
            for filelist_list in self.file_sources.file_list:
                for filelist in filelist_list:
                    self._scan_file_list(filelist, self.file_pattern.replace("*", ""))

    def _process_string_sources(self):
        """Adds each results string dictionary item to the queue

        Raises:
            RTCoreError
        """
        if self.string_sources != None:
            for (
                ligand_name,
                docking_result,
            ) in self.string_sources.results_strings.items():
                string_data = {ligand_name: docking_result}
                self._add_to_queue(string_data, string=True)
        else:
            raise RTCoreError(
                "There was an error while reading the results string input."
            )

    def _add_to_queue(self, results_data, string=False):
        """_summary_

        Args:
            results_data (string or dict): results data provided as a file path or a dictionary kw pair
            string (bool, optional): switch if results provided as a string

        Raises:
            MultiprocessingError
        """
        # adds result file to the multiprocessing queue
        max_attempts = 750
        timeout = 0.5  # seconds
        if string == False and self.receptor_file is not None:
            if (
                os.path.split(results_data)[-1] == os.path.split(self.receptor_file)[-1]
            ):  # check that we don't try to add the receptor
                return
        attempts = 0
        while True:
            if attempts >= max_attempts:
                raise MultiprocessingError(
                    "Something is blocking the progressing of results data reading. Exiting program."
                ) from queue.Full
            try:
                self.managed_queue_in.put(results_data, block=True, timeout=timeout)
                self.num_files += 1
                self._check_for_worker_exceptions()
                break
            except queue.Full:
                attempts += 1
                self._check_for_worker_exceptions()

    def _check_for_worker_exceptions(self):
        if self.p_conn.poll():
            error, tb, filename = self.p_conn.recv()
            self.logger.error(f"Caught error in multiprocessing from {filename}")
            # don't kill parser errors, only database error
            if filename == "Database":
                self._kill_all_workers(error, filename, tb)
            else:
                with open("ringtail_failed_files.log", "a") as f:
                    f.write(
                        str(datetime.now()) + f"\tRingtail failed to parse {filename}\n"
                    )
                    self.logger.debug(tb)

    def _kill_all_workers(self, error, filename, tb):
        for s in self.workers:
            s.kill()
        self.logger.debug(f"Error encountered while handling {filename}")
        self.logger.debug(tb)
        raise error

    def _scan_dir(self, path, pattern, recursive=False):
        """scan for valid output files in a directory the pattern is used
        to glob files optionally, a recursive search is performed

        Args:
            path (str): folder path
            pattern (str): file extension
            recursive (bool, optional): look for files and folders recursively

        Yields:
            list: of file paths found in the search
        """
        self.logger.info(
            "Scanning directory [%s] for files (pattern:|%s|)" % (path, pattern)
        )
        if recursive:
            path = os.path.normpath(path)
            path = os.path.expanduser(path)
            for dirpath, dirnames, filenames in os.walk(path):
                yield (  # <----
                    os.path.join(dirpath, f)
                    for f in fnmatch.filter(filenames, "*" + pattern)
                )
        else:
            yield glob(os.path.join(path, pattern))  # <----

    def _scan_file_list(self, filename, pattern):
        """read file names from file list and ensures they exist,
        then adding them to the list of files to be processed

        Args:
            filename (str): filename provided in list
            pattern (str): file extension

        Raises:
            MultiprocessingError
        """

        lig_accepted = []
        c = 0
        with open(filename, "r") as fp:
            for line in fp.readlines():
                line = line.strip()
                c += 1
                if os.path.isfile(line):
                    if line.endswith(pattern) or line.endswith(
                        pattern + ".gz"
                    ):  # NOTE if adding zip option change here
                        lig_accepted.append(line)
                else:
                    self.logger.warning("Warning! file |%s| does not exist" % line)
        if len(lig_accepted) == 0:
            raise MultiprocessingError(
                "*ERROR* No valid files were found when reading from |%s|" % filename
            )
        self.logger.info(
            "# [ %5.3f%% files in list accepted (%d) ]"
            % ((len(lig_accepted) / c * 100, c))
        )

        for file in lig_accepted:
            if file != self.receptor_file:
                self._add_to_queue(file)


import platform
import time
import sys
from .logutils import LOGGER as logger
import traceback
import queue
from .parsers import parse_single_dlg, parse_vina_result
from .exceptions import (
    FileParsingError,
    WriteToStorageError,
    MultiprocessingError,
    ResultsProcessingError,
)
from .interactions import InteractionFinder
from .storagemanager import *
import multiprocessing


class DockingFileReader(
    multiprocessing.Process,
):
    """This class is the individual worker for processing docking results.
    One instance of this class is instantiated for each available processor.

    Attributes:
        queueIn (multiprocessing.Queue): current queue for the processor/file reader
        queueOut (multiprocessing.Queue): queue for the processor/file reader after adding or removing an item
        pipe_conn (multiprocessing.Pipe): pipe connection to the reader
        storageman (StorageManager): storageman object
        storageman_class (StorageManager): storagemanager child class/database type
        docking_mode (str): describes what docking engine was used to produce the results
        max_poses (int): max number of poses to store for each ligand
        interaction_tolerance (float): Will add the interactions for poses within some tolerance RMSD range of the top pose in a cluster to that top pose."
        store_all_poses (bool): Store all poses from docking results
        add_interactions (bool): find and save interactions between ligand poses and receptor
        interaction_cutoffs (list(float)): cutoff for interactions of hydrogen bonds and VDW interactions, in ångströms
        target (str): receptor name
        string_processing (bool, optional): switch for processing result strings
    """

    def __init__(
        self,
        queueIn,
        queueOut,
        pipe_conn,
        db_file,
        storageman_class,
        docking_mode,
        max_poses,
        interaction_tolerance,
        store_all_poses,
        target,
        add_interactions,
        interaction_cutoffs,
        receptor_file,
        string_processing=False,
    ):
        print("         I doubt I make it this far")
        # set docking_mode for which file parser to use (and for vina, '_string' if parsing string output directly)
        self.docking_mode = docking_mode
        # set number of clusters to write
        self.max_poses = max_poses
        self.store_all_poses_flag = store_all_poses
        # set interaction_tolerance cutoff
        self.interaction_tolerance = interaction_tolerance
        # set options for finding interactions
        self.add_interactions = add_interactions
        self.interaction_cutoffs = interaction_cutoffs
        self.receptor_file = receptor_file
        multiprocessing.Process.__init__(self)
        self.storageman_class = storageman_class
        self.db_file = db_file
        # set target name to check against
        self.target = target
        # initialize the parent class to inherit all multiprocessing methods

        # # each worker knows the queue in (where data to process comes from)
        # # ...and a queue out (where to send the results)
        self.queueIn = queueIn
        self.queueOut = queueOut
        # # ...and a pipe to the parent
        self.pipe = pipe_conn
        self.interaction_finder = None
        self.exception = None
        # if the results being processed comes as a string instead of a file (currently only implemented for vina)
        self.string_processing = string_processing

    def _find_best_cluster_poses(self, ligand_dict):
        """Takes input ligand dictionary, reads run pose clusters, adds "cluster_best_run"
        entry with the top scoring run for each cluster

        Args:
            ligand_dict (dict): dictionary of ligands

        Returns:
            dict: top poses in cluster for each ligand
        """

        top_poses = []
        cluster_dict = ligand_dict["clusters"]
        for cluster in cluster_dict:
            top_poses.append(cluster_dict[cluster][0])
        ligand_dict["cluster_top_poses"] = top_poses
        return ligand_dict

    def run(self):
        """Method overload from parent class .This is where the task of this class is performed.
        Each multiprocessing.Process class must have a "run" method which is called by the
        initialization (see below) with start()

        Raises:
            NotImplementedError: if parser for specific docking result type is not implemented
            FileParsingError
        """

        while True:
            try:
                # retrieve from the queue in the next task to be done
                next_task = self.queueIn.get()
                if type(next_task) == dict:
                    text = list(next_task.keys())[0]
                else:
                    text = next_task
                logger.debug("Next Task: " + str(text))
                # if a poison pill is received, this worker's job is done, quit
                if next_task is None:
                    # before leaving, pass the poison pill back in the queue
                    self.queueOut.put(None)
                    break

                # generate CPU LOAD
                # parser depends on requested docking_mode
                if self.docking_mode == "dlg":
                    parsed_file_dict = parse_single_dlg(next_task)
                    # find the run number for the best pose in each cluster for adgpu
                    parsed_file_dict = self._find_best_cluster_poses(parsed_file_dict)
                elif self.docking_mode == "vina":
                    parsed_file_dict = parse_vina_result(
                        next_task, self.string_processing
                    )

                # Example code for calling user-implemented docking_mode
                # elif self.docking_mode == "my_docking_mode":
                #     parsed_file_dict = myparser(next_task)
                else:
                    raise NotImplementedError(
                        f"Parser for input file docking_mode {self.docking_mode} not implemented!"
                    )
                # check receptor name from file against that which we expect
                if (
                    parsed_file_dict["receptor"] != self.target
                    and self.target is not None
                    and self.docking_mode == "dlg"
                ):
                    raise FileParsingError(
                        "Receptor name {0} in {1} does not match given target name {2}. Please ensure that this file belongs to the current virtual screening.".format(
                            parsed_file_dict["receptor"], next_task, self.target
                        )
                    )

                # find run numbers for poses we want to save
                parsed_file_dict["poses_to_save"] = self._find_poses_to_save(
                    parsed_file_dict
                )
                # Calculate interactions if requested
                if self.add_interactions:
                    try:
                        from .receptormanager import ReceptorManager as rm

                        with self.storageman_class(self.db_file) as self.storageman:
                            # grab receptor info from database, this assumes there is only one receptor in the database
                            receptor_blob = self.storageman.fetch_receptor_objects()[0][
                                1
                            ]  # method returns and iter of tuples, blob is the second tuple element in the first list element
                            # convert receptor blob to string
                            receptor_string = rm.blob2str(receptor_blob)
                    except:
                        raise ResultsProcessingError(
                            "add_interactions was requested, but cannot find the receptor in the database. Please ensure to include the receptor_file and save_receptor if the receptor has not already been added to the database."
                        )
                    if self.interaction_finder is None:
                        self.interaction_finder = InteractionFinder(
                            receptor_string, self.interaction_cutoffs
                        )
                    if parsed_file_dict["interactions"] == []:
                        for pose in parsed_file_dict["pose_coordinates"]:
                            parsed_file_dict["interactions"].append(
                                self.interaction_finder.find_pose_interactions(
                                    parsed_file_dict["ligand_atomtypes"], pose
                                )
                            )
                            parsed_file_dict["num_interactions"].append(
                                int(parsed_file_dict["interactions"][-1]["count"][0])
                            )
                            parsed_file_dict["num_hb"].append(
                                len(
                                    [
                                        1
                                        for i in parsed_file_dict["interactions"][-1][
                                            "type"
                                        ]
                                        if i == "H"
                                    ]
                                )
                            )
                # find poses we want to save tolerated interactions for
                if self.interaction_tolerance is not None:
                    parsed_file_dict["tolerated_interaction_runs"] = (
                        self._find_tolerated_interactions(parsed_file_dict)
                    )
                else:
                    parsed_file_dict["tolerated_interaction_runs"] = []
                # put the result in the out queue
                data_packet = self.storageman_class.format_for_storage(parsed_file_dict)
                self._add_to_queueout(data_packet)
            except Exception:
                tb = traceback.format_exc()
                self.pipe.send(
                    (
                        FileParsingError(f"Error while parsing {next_task}"),
                        tb,
                        next_task,
                    )
                )

    def _add_to_queueout(self, obj):
        max_attempts = 750
        timeout = 0.5  # seconds
        attempts = 0
        while True:
            if attempts >= max_attempts:
                raise MultiprocessingError(
                    "Something is blocking the progressing of file writing. Exiting program."
                ) from queue.Full
            try:
                self.queueOut.put(obj, block=True, timeout=timeout)
                break
            except queue.Full:
                logger.debug(
                    f"Queue full: queueOut.put attempt {attempts} timed out. {max_attempts - attempts} put attempts remaining."
                )
                attempts += 1

    def _find_poses_to_save(self, ligand_dict: dict) -> list:
        """Returns list of the run numbers for the top run in the
        top self.max_pose clusters (ADGPU) or just the top poses overall

        Args:
            ligand_dict (dict): Dictionary of ligand data from parser

        Returns:
            list: List of run numbers to save
        """
        if self.store_all_poses_flag:
            poses_to_save = ligand_dict["sorted_runs"]
        elif self.docking_mode == "dlg":
            # will only select top n clusters. Default 3
            poses_to_save = ligand_dict["cluster_top_poses"][: self.max_poses]
        # if not adgpu, save top n poses
        else:
            poses_to_save = ligand_dict["sorted_runs"][: self.max_poses]

        return poses_to_save

    def _find_tolerated_interactions(self, ligand_dict):
        """Take ligand dict and finds which poses we should save the
        interactions for as tolerated interactions for the top pose
        of the cluster. These runs are within the
        <self.interaction_tolerance> angstroms RMSD of the top pose
        for a given cluster. All data for the cluster's top pose is saved.

        Args:
            ligand_dict (dict): Dictionary of ligand data from parser

        Returns:
            list: run numbers of tolerated runs
        """
        tolerated_runs = []
        for idx, run in enumerate(ligand_dict["sorted_runs"]):
            if float(ligand_dict["cluster_rmsds"][idx]) <= self.interaction_tolerance:
                tolerated_runs.append(run)
        return tolerated_runs


class Writer(multiprocessing.Process):
    """This class is a listener that retrieves data from the queue and writes it
    into datbase"""

    def __init__(
        self,
        queue,
        num_readers,
        pipe_conn,
        chunksize,
        docking_mode,
        db_file,
        storageman_class,
    ):
        multiprocessing.Process.__init__(self)
        self.queue = queue
        # this class knows about how many multi-processing workers there are and where the pipe to the parent is
        self.num_readers = num_readers
        self.pipe = pipe_conn
        # assign pointer to storage object, set chunksize
        self.docking_mode = docking_mode
        self.chunksize = chunksize
        self.storageman_class = storageman_class
        self.db_file = db_file
        # initialize data array (stack of dictionaries)
        self.results_array = []
        self.ligands_array = []
        self.interactions_list = []
        self.receptor_array = []
        # progress tracking instance variables
        self.first_insert = True
        self.counter = 0
        self.num_files_written = 0
        self.time0 = time.perf_counter()
        self.last_write_time = 0

    def run(self):
        """Method overload from parent class. This is where the task of this class
        is performed. Each multiprocessing.Process class must have a "run" method which
        is called by the initialization (see below) with start()

        Raises:
            WriteToStorageError
        """

        try:
            while True:
                # retrieve the next task from the queue
                next_task = self.queue.get()
                if next_task is None:
                    # if a poison pill is found, it means one of the workers quit
                    self.num_readers -= 1
                    logger.debug(
                        "Closing process. Remaining open processes: "
                        + str(self.num_readers)
                    )
                else:
                    # if not a poison pill, process the task item

                    # after every n (chunksize) files, write to storage
                    if self.counter >= self.chunksize:
                        self.write_to_storage()
                        # print info about files and time remaining
                        sys.stdout.write("\r")
                        sys.stdout.write(
                            "{0} files written to database. Writing {1:.0f} files/minute. Elapsed time {2:.0f} seconds.".format(
                                self.num_files_written,
                                self.num_files_written * 60 / self.total_runtime,
                                self.total_runtime,
                            )
                        )
                        sys.stdout.flush()

                    # process next file
                    self.process_file(next_task)
                if self.num_readers == 0:
                    # received as many poison pills as workers
                    logger.info("Performing final database write")
                    # perform final storage write
                    self.write_to_storage(final=True)
                    # no workers left, no job to do
                    logger.info("File processing completed")
                    self.close()
                    break
        except Exception:
            tb = traceback.format_exc()
            self.pipe.send(
                (
                    WriteToStorageError("Error occured while writing database"),
                    tb,
                    "Database",
                )
            )

    def write_to_storage(self, final=False):
        """Inserting data to the database through the designated storagemanager.

        Args:
            final (bool): if last data entry, finalize database
        """
        # insert result, ligand, and receptor data
        with self.storageman_class(self.db_file) as self.storageman:
            self.storageman.insert_data(
                self.results_array,
                self.ligands_array,
                self.interactions_list,
                self.receptor_array,
                self.first_insert,
            )
        # So at this point the ligand array is empty
        if self.first_insert:  # will only insert receptor for first insertion
            self.first_insert = False

        # calulate time for processing/writing speed
        self.num_files_written += self.chunksize
        self.total_runtime = time.perf_counter() - self.time0

        # reset data holder for next chunk
        self.results_array = []
        self.ligands_array = []
        self.interactions_list = []
        self.receptor_array = []
        self.counter = 0

        if final:
            # if final write, tell storageman to index
            with self.storageman_class(self.db_file) as self.storageman:
                self.storageman.create_indices()
                self.storageman.set_ringtail_db_schema_version()

    def process_file(self, file_packet):
        """Breaks up the data in the file_packet to distribute between
        the different arrays to be inserted in the database.

        Args:
            file_packet (any): File packet to be processed
        """
        results_rows, ligand_row, interaction_rows, receptor_row = file_packet
        for pose in results_rows:
            self.results_array.append(pose)
        for pose in interaction_rows:
            self.interactions_list.append(pose)
        self.ligands_array.append(ligand_row)
        self.receptor_array.append(receptor_row)

        self.counter += 1
