#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail virtual screening manager
#

import matplotlib.pyplot as plt
import json
from meeko import RDKitMolCreate
from .storagemanager import StorageManager, StorageManagerSQLite
from .resultsmanager import ResultsManager
from .receptormanager import ReceptorManager
from .outputmanager import OutputManager
from .ringtailoptions import *
from .exceptions import RTCoreError, OutputError
from rdkit import Chem
import itertools
import os
from os import path
from .logmanager import logger

class RingtailCore:
    """Core class for coordinating different actions on virtual screening
    i.e. adding results to storage, filtering, output options

    #TODO
    Attributes:
        generaloptions (GeneralOptions object): sets logging level and general processing options, including read/write, and outputting summary to console
        db_file (str): name of database file being operated on
        storageman (StorageManager object): Interface module with database
        storageopts (StorageOptions object): options for storageman
        file_sources (InputFiles object): object containing all specified ligand result filtes and receptor file
        resultsman (ResultsManager object): Module to deal with result file processing before adding to databae
        resultsmanopts (ResultsProcessingOptions object): settings for processing result files
        outputman (OutputManager object): Manager for output tasks of log-writting, plotting, ligand SDF writing
        readopts (ReadOptions object): options for outputman and methods to process data with
        filterobj (Filters object): object holding all optional filters
    """
#-#-#- Base methods -#-#-#
    
    def __init__(self, db_file = "output.db", storage_type = "sqlite", process_mode="read"):
        """
        Initialize RingtailCore object and create a storageman object with the db file.
        Does not open access to the storage. Future option will include opening database as readonly.
        """
        self.db_file = db_file
        self.process_mode = process_mode
        storageman = StorageManager.check_storage_compatibility(storage_type) 
        self.storageman = storageman(db_file)
        self.set_general_options()
             
    def update_database_version(self, consent=False):
        # Method to update database version
        return self.storageman.update_database_version(consent)
    
#-#-#- Private methods -#-#-#
    
    def _add_poses(
        self,
        atom_indices,
        poses,
        mol,
        flexres_mols,
        flexres_info,
        ligand_saved_coords,
        flexres_saved_coords,
        properties,
    ):
        """Add poses from given cursor to rdkit mols for ligand and flexible residues

        Args:
            atom_indices (List): List of ints indicating mapping of coordinate indices to smiles indices
            poses (iterable): iterable containing ligand_pose, flexres_pose, flexres_names
            mol (RDKit Mol): RDKit molecule for ligand
            flexres_mols (list): list of rdkit molecules for flexible residues
            flexres_info (list): list of tuples containing info for each flexible residue (res_smiles, res_index_map, res_h_parents)
            ligand_saved_coords (list): list of coordinates to save for adding hydrogens later
            flexres_saved_coords (list): list of lists of flexres coords to save for adding hydrogens later
            properties (dict): Dictionary of lists of properties, with each element corresponding to that conformer in the rdkit mol

        """
        for (
            Pose_ID,
            docking_score,
            leff,
            ligand_pose,
            flexres_pose,
        ) in poses:
            # fetch info about pose interactions and format into string with format <type>-<chain>:<resname>:<resnum>:<atomname>:<atomnumber>, joined by commas
            pose_bitvector = self.storageman.fetch_interaction_bitvector(Pose_ID)
            if pose_bitvector is not None:
                interaction_indices = []
                interactions_list = []
                for idx, bit in enumerate(pose_bitvector):
                    if bit == 1:
                        interaction_indices.append(
                            idx + 1
                        )  # adjust for indexing starting at 1
                for int_idx in interaction_indices:
                    interaction_info = self.storageman.fetch_interaction_info_by_index(
                        int_idx
                    )
                    interaction = (
                        interaction_info[0] + "-" + ":".join(interaction_info[1:])
                    )
                    interactions_list.append(interaction)
                interactions_str = ", ".join(interactions_list)
                properties["Interactions"].append(interactions_str)
            # add properties to dictionary lists
            properties["Binding energies"].append(docking_score)
            properties["Ligand effiencies"].append(leff)
            # get pose coordinate info
            ligand_pose = json.loads(ligand_pose)
            flexres_pose = json.loads(flexres_pose)
            mol = RDKitMolCreate.add_pose_to_mol(mol, ligand_pose, atom_indices)
            for fr_idx, fr_mol in enumerate(flexres_mols):
                flexres_mols[fr_idx] = RDKitMolCreate.add_pose_to_mol(
                    fr_mol, flexres_pose[fr_idx], flexres_info[fr_idx][1]
                )
                flexres_saved_coords[fr_idx].append(flexres_pose[fr_idx])
            ligand_saved_coords.append(ligand_pose)
        return mol, flexres_mols, ligand_saved_coords, flexres_saved_coords, properties

    def _produce_summary(self, columns=["docking_score", "leff"], percentiles=[1, 10]) -> None:
        """Print summary of data in storage
        """
        with self.storageman: summary_data = self.storageman.fetch_summary_data(columns, percentiles)
        print("Total Stored Ligands          :", summary_data.pop("num_ligands"))
        print("Total Stored Poses            :", summary_data.pop("num_poses"))
        print("Total Unique Interactions     :", summary_data.pop("num_unique_interactions"))
        print("Number Interacting Residues   :", summary_data.pop("num_interacting_residues"))

        colon_col = 18
        print("\nEnergy statistics:")
        print("=======================================")
        for col in columns:
            if col == "docking_score":
                min_e = summary_data["min_docking_score"]
                max_e = summary_data["max_docking_score"]
                print(f"Energy (min)      : {min_e:.2f} kcal/mol")
                print(f"Energy (max)      : {max_e:.2f} kcal/mol")
            elif col == "leff":
                min_le = summary_data["min_leff"]
                max_le = summary_data["max_leff"]
                print(f"LE     (min)      : {min_le:.2f} kcal/mol/heavyatom")
                print(f"LE     (max)      : {max_le:.2f} kcal/mol/heavyatom")
            else:
                min_col = summary_data[f"min_{col}"]
                max_col = summary_data[f"max_{col}"]
                print(f"{col} (min) : {min_col}")
                print(f"{col} (max) : {max_col}")
        if percentiles != [] and percentiles is not None:
            print("----------------------------------------")

            for col in columns:
                for p in percentiles:
                    if col == "docking_score":
                        p_string = f"Energy (top {p}% )"
                        p_string += ' ' * (colon_col - len(p_string))
                        print(f"{p_string}: {summary_data[f'{p}%_docking_score']:.2f} kcal/mol")
                    elif col == "leff":
                        p_string = f"LE     (top {p}% )"
                        p_string += ' ' * (colon_col - len(p_string))
                        print(f"{p_string}: {summary_data[f'{p}%_leff']:.2f} kcal/mol/heavyatom")
                    else:
                        p_string = f"{col} (top {p}%)"
                        p_string += ' ' * (colon_col - len(p_string))
                        print(f"{p_string}: {summary_data[f'{p}%_{col}']:.2f}")

    def _generate_interaction_combinations(self, max_miss=0):
        """Recursive function to list of tuples of possible interaction filter combinations, excluding up to max_miss interactions per filtering round

        Args:
            max_miss (int): Maximum number of interactions to be excluded
        """

        all_interactions = []
        for _type in Filters.get_filter_keys("interaction"):
            interactions = getattr(self.filterobj, _type)
            for interact in interactions:
                all_interactions.append(_type + "-" + interact[0])
        # warn if max_miss greater than number of interactions
        if max_miss > len(all_interactions):
            logger.warning(
                "Requested max_miss options greater than number of interaction filters given. Defaulting to max_miss = number interaction filters"
            )
            max_miss = len(all_interactions)

        # BASE CASE:
        if max_miss == 0:
            return [tuple(all_interactions)]
        else:
            combinations = list(
                itertools.combinations(
                    all_interactions, len(all_interactions) - max_miss
                )
            )
            return combinations + self._generate_interaction_combinations(
                max_miss=max_miss - 1
            )

    def _prepare_filters_for_storageman(self, interaction_combination):
        """Takes desired interaction combination, formats Filter object to dict, removes interactions not in given interaction_combination

        Args:
            interaction_combination (list): list of interactions to be included in this round of filtering
        Returns:
            dict: dictionary of filters for storageman
        """

        filters_dict = self.filterobj.todict()
        for itype in Filters.get_filter_keys("interaction"):
            itype_interactions = filters_dict[itype]
            for interaction in itype_interactions:
                if itype + "-" + interaction[0] not in interaction_combination:
                    filters_dict[itype].remove(interaction)

        return filters_dict

    def _docking_mode_from_file_extension(self, file_pattern: str) -> str:
        if file_pattern.lower() == "*.dlg*":
            mode = "dlg"
        elif file_pattern.lower() == "*.pdbqt*":
            mode = "vina"
        else:
            raise OptionError(f'{file_pattern} is not a valid docking file format for ringtail.')
        return mode
    
#-#-#- API -#-#-#
    #-#- Processing methods -#-# 

    def add_options_from_file(self, options_file: str ="options.json"):
        """
        Provide ringtail options from file, *not currently in use
        Args:
            config_file: json formatted file containing ringtail and filter options
        """
        
        if options_file is None: #or not json compatible
            raise OptionError("No option file was found in the Ringtail/util_files directory.")
        
        filepath = self._options_file_path(options_file)
        with open(filepath, "r") as f:
            logger.info("Reading Ringtail options from options file")
            options: dict = json.load(f)

        # Set each given object option to dict to respective class or manager
        optmap = {
            "generalopts": self.set_general_options,
            "writeopts": self.set_results_processing_options,
            "storageopts": self.set_storage_options,
            "readopts": self.set_read_options,
            "filterobj": self.set_filters,
            "fileobj": self.set_file_sources
        }
        for k, v in options.items():
            optmap[k](dict=v)
            logger.debug(f'{optmap[k]} was ran with these options: {v}')
    
    def add_results_from_files(self,
                               file = None, 
                               file_path = None, 
                               file_list = None, 
                               file_pattern = None, 
                               recursive = None, 
                               receptor_file=None,
                               save_receptor = None,
                               file_source_object = None,
                               store_all_poses: bool = None,
                               max_poses: int = None,
                               add_interactions: bool = None,
                               interaction_tolerance: float = None,
                               interaction_cutoffs: list = None,
                               max_proc: int = None,
                               optionsdict=None
                               ):
        #TODO Can files be added as lists and as single paths? 
        """
        Call storage manager to process result files and add to database.
        Options can be provided as a dict or as individual options.
        Creates a database, or adds to an existing one if using "append_results" in storageman_opts

        Args:
            file (str): ligand result file
            file_path (list(str)): list of folders containing one or more result files
            file_list (list(str)): list of ligand result file(s)
            file_pattern (str): file pattern to use with recursive search in a file_path, "*.dlg*" for AutoDock-GDP and "*.pdbqt*" for vina
            recursive (bool): used to recursively search file_path for folders inside folders
            receptor_file (str): string containing the receptor .pdbqt
            save_receptor (bool): whether or not to store the full receptor details in the database (needed for some things)
                optional: file_source_object (InputFiles): file sources already as an object #TODO how do I handle file objects in all of this
            store_all_poses (bool): store all ligand poses, does it take precedence over max poses? 
            max_poses (int): how many poses to save (ordered by soem score?)
            add_interactions (bool): add ligand-receptor interaction data, only in vina mode
            interaction_tolerance (float): longest ångström distance that is considered interaction?
            interaction_cutoffs (list): ångström distance cutoffs for x and y interaction
            max_proc (int): max number of computer processors to use for file reading
            optionsdict (dict): write options as a dict

        """

        self.process_mode = "write"

        if file_source_object is not None:
            files = file_source_object
        else: 
            self.set_file_sources(file, file_path, file_list, file_pattern, recursive, receptor_file, save_receptor)
            files = self.files
    

        results_files_given = (files.file is not None or files.file_path is not None or files.file_list is not None)

        if not results_files_given and not files.save_receptor:
            raise OptionError("At least one input option needs to be used: --file, --file_path, --file_list, or --input_db and --save_receptor")

        if not hasattr(self, "storageopts"):
            self.set_storage_options()

        # If there are any ligand files, process ligand data
        if results_files_given: 
            with self.storageman:
                if self.storageopts.append_results and not RTOptions.is_valid_path(self.db_file):
                    raise OptionError("The provided --input_db is not a valid path, please check the provided path.")
                
                self.set_results_processing_options(store_all_poses, max_poses, add_interactions, interaction_tolerance, interaction_cutoffs, max_proc, optionsdict)
                
                # Docking mode compatibility check
                if self.generalopts.docking_mode == "vina" and self.resultsmanopts.interaction_tolerance is not None:
                    logger.warning("Cannot use interaction_tolerance with Vina mode. Removing interaction_tolerance.")
                    self.interaction_tolerance = None

                # Prepare the results manager object
                if not hasattr(self, "resultsman"):self.resultsman = ResultsManager(file_sources=files)
                self.resultsman.storageman = self.storageman
                self.resultsman.storageman_class = self.storageman.__class__
                for k,v in self.resultsmanopts.todict().items():  
                    setattr(self.resultsman, k, v)
                self.resultsman.mode = self.generalopts.docking_mode
                logger.debug("Results manager object has been initialized.")

                # Process results files and handle database versioning 
                self.storageman.check_storage_ready()
                logger.info("Adding results...")
                self.resultsman.process_results() 
                self.storageman.set_ringtaildb_version()
                #-#-#- I get this far
                if self.generalopts.summary: self._produce_summary()

        
        if files.save_receptor: 
            self.save_receptor(files.receptor_file)
    
    def save_receptor(self, receptor_file):
            """
            Add receptor to database. Context managed by self.storageman

            Args:
                receptors (list): list of receptor blobs to add to database
                * currently only one receptor allowed per database
            """
            receptor_list = ReceptorManager.make_receptor_blobs([receptor_file])
            with self.storageman:
                for rec, rec_name in receptor_list:
                    # NOTE: in current implementation, only one receptor allowed per database
                    # Check that any receptor row is incomplete (needs receptor blob) before inserting
                    filled_receptor_rows = self.storageman.count_receptors_in_db()
                    if filled_receptor_rows != 0:
                        raise RTCoreError(
                            "Expected Receptors table to have no receptor objects present, already has {0} receptor present. Cannot add more than 1 receptor to a database.".format(
                                filled_receptor_rows
                            )
                        )
                    self.storageman.save_receptor(rec)
                    logger.info("Receptor data was added to the database.")

    def filter(self, enumerate_interaction_combs=False, return_iter=False):
        """
        Prepare list of filters, then hand it off to storageManager to
            perform filtering. Create log of passing results.
        Args:
            enumerate_interaction_combs (bool): inherently handled atm, might need to be depreceated?
        """
        #NOTE This is clunky for now, but sets options to default if they have not yet been set.
        self.process_mode = "read"

        if not hasattr(self, "filterobj"):
            logger.debug("No filters have been set, using default values found in Filters class") 
            self.set_filters()
        if not hasattr(self, "storageopts"):
            logger.debug("No storage options have been set, using default values found in 'set_storage_options'") 
            self.set_storage_options()

        if not hasattr(self, "readopts"):
            logger.debug("No read options have been set, using default values found in 'set_read_options'") 
            self.set_read_options()
        # make sure enumerate_interaction_combs always true if max_miss = 0, since we don't ever worry about the union in this case
        if self.filterobj.max_miss == 0:
            self.readopts.enumerate_interaction_combs = True

        # guard against unsing percentile filter with all_poses
        if self.storageopts.output_all_poses and not (self.filterobj.score_percentile is None or self.filterobj.le_percentile is None):
            logger.warning(
                "Cannot return all passing poses with percentile filter. Will only log best pose."
            )
            self.storageopts.output_all_poses = False

        logger.info("Filtering results...")

        # get possible permutations of interaction with max_miss excluded
        interaction_combs = self._generate_interaction_combinations(
                            self.filterobj.max_miss)
        ligands_passed = 0
        '''This for comprehension takes all combinations represented in one union of one or multiple, and filters, and goes around until all combinations have been used to filter'''
        with self.storageman:
            for ic_idx, combination in enumerate(interaction_combs):
                # prepare Filter object with only desired interaction combination for storageManager
                filters_dict = self._prepare_filters_for_storageman(combination)
                # set storageMan's internal ic_counter to reflect current ic_idx
                if len(interaction_combs) > 1:
                    self.storageman.set_view_suffix(ic_idx)
                # ask storageManager to fetch results
                filtered_results = self.storageman.filter_results(
                    filters_dict, not self.readopts.enumerate_interaction_combs
                )
                if filtered_results is not None:
                    if return_iter:
                        return filtered_results
                    result_bookmark_name = self.storageman.get_current_view_name()
                    with self.outputman: 
                        self.outputman.write_filters_to_log(self.filterobj.todict(), combination, f"Morgan Fingerprints butina clustering cutoff: {self.storageman.mfpt_cluster}\nInteraction Fingerprints clustering cutoff: {self.storageman.interaction_cluster}")
                        self.outputman.write_results_bookmark_to_log(result_bookmark_name)
                        number_passing = self.outputman.write_log(filtered_results)
                        self.outputman.log_num_passing_ligands(number_passing)
                        print("\nNumber of ligands passing filters:", number_passing)
                        ligands_passed = number_passing
                else:
                    logger.warning("WARNING: No ligands found passing filters")
            if len(interaction_combs) > 1:
                maxmiss_union_results = self.storageman.get_maxmiss_union(len(interaction_combs))
                with self.outputman:
                    self.outputman.write_maxmiss_union_header()
                    self.outputman.write_results_bookmark_to_log(self.storageman.results_view_name + "_union")
                    number_passing_union = self.outputman.write_log(maxmiss_union_results)
                    self.outputman.log_num_passing_ligands(number_passing_union)
                    print("\nNumber passing ligands in max_miss union:", number_passing_union)
                    ligands_passed = number_passing_union
            
        return ligands_passed

    #-#-#- Methods to set ringtail options explicitly -#-#-#
    def set_general_options(self, docking_mode=None,
                                summary=None,
                                verbose=None,
                                debug=None,
                                rtopts=None,
                                dict: dict=None):
        """
        Settings for ringtail

        Args:
            docking_mode (str): dlg (autodock-gpu) or vina 
            summary (bool): print database summary to terminal
            verbose (bool): set logging level to "info" (second lowest level)
            debug (bool): set logging level to "debug" (lowest level)
            dict (dict): dictionary of one or more of the above args, is overwritten by individual args
        """
        # Dict of individual arguments
        indiv_options: dict = vars(); 
        del indiv_options["self"]; del indiv_options["dict"]

        # Create option object with default values if needed
        if not hasattr(self, "generalopts"): self.generalopts = GeneralOptions()
            
        # Set options from dict if provided
        if dict is not None:
            for k,v in dict.items():
                setattr(self.generalopts, k, v) 

        # Set additional options from individual arguments
        #NOTE Will overwrite config file
        for k,v in indiv_options.items():
            if v is not None: setattr(self.generalopts, k, v)

        # Assign attributes to storage manager
        for k,v in self.generalopts.todict().items():
            setattr(self, k, v)

        # Set logger level
        if self.generalopts.debug:
            logger.setLevel("DEBUG")
        elif self.generalopts.verbose:
            logger.setLevel("INFO")

    def set_storage_options(self, 
                            filter_bookmark = None,
                            append_results = None,
                            duplicate_handling = None,
                            overwrite_log_file = None,
                            order_results_by = None,
                            outfields = None,
                            output_all_poses = None,
                            mfpt_cluster = None,
                            interaction_cluster = None,
                            results_view_name =None,
                            dict=None):
        """
        Options or storagemanager, used both in read and write.

        Args:
            filter_bookmark (str): Perform filtering over specified bookmark. (in output group in CLI)
            append_results (bool): Add new results to an existing database, specified by database choice in ringtail initialization or --input_db in cli
            duplicate_handling (str, options): specify how duplicate Results rows should be handled when inserting into database. Options are "ignore" or "replace". Default behavior will allow duplicate entries.
            overwrite_logfile (bool): by default, if a log file exists, it doesn't get overwritten and an error is returned; this option enable overwriting existing log files. Will also overwrite existing database
            order_results_by (str): Stipulates how to order the results when written to the log file. By default will be ordered by order results were added to the database. ONLY TAKES ONE OPTION."
                    "available fields are:  "
                    '"e" (docking_score), '
                    '"le" (ligand efficiency), '
                    '"delta" (delta energy from best pose), '
                    '"ref_rmsd" (RMSD to reference pose), '
                    '"e_inter" (intermolecular energy), '
                    '"e_vdw" (van der waals energy), '
                    '"e_elec" (electrostatic energy), '
                    '"e_intra" (intermolecular energy), '
                    '"n_interact" (number of interactions), '
                    '"rank" (rank of ligand pose), '
                    '"run" (run number for ligand pose), '
                    '"hb" (hydrogen bonds); '
            outfields (str): defines which fields are used when reporting the results (to stdout and to the log file); fields are specified as comma-separated values, e.g. "--outfields=e,le,hb"; by default, docking_score (energy) and ligand name are reported; ligand always reported in first column available fields are:  '
                    '"Ligand_name" (Ligand name), '
                    '"e" (docking_score), '
                    '"le" (ligand efficiency), '
                    '"delta" (delta energy from best pose), '
                    '"ref_rmsd" (RMSD to reference pose), '
                    '"e_inter" (intermolecular energy), '
                    '"e_vdw" (van der waals energy), '
                    '"e_elec" (electrostatic energy), '
                    '"e_intra" (intermolecular energy), '
                    '"n_interact" (number of interactions), '
                    '"ligand_smile" , '
                    '"rank" (rank of ligand pose), '
                    '"run" (run number for ligand pose), '
                    '"hb" (hydrogen bonds), '
                    '"receptor" (receptor name); '
                    "Fields are printed in the order in which they are provided. Ligand name will always be returned and will be added in first position if not specified.
            output_all_poses (bool): By default, will output only top-scoring pose passing filters per ligand. This flag will cause each pose passing the filters to be logged.
            mfpt_cluster (float): Cluster filered ligands by Tanimoto distance of Morgan fingerprints with Butina clustering and output ligand with lowest ligand efficiency from each cluster. Default clustering cutoff is 0.5. Useful for selecting chemically dissimilar ligands.
            interaction_cluster (float): Cluster filered ligands by Tanimoto distance of interaction fingerprints with Butina clustering and output ligand with lowest ligand efficiency from each cluster. Default clustering cutoff is 0.5. Useful for enhancing selection of ligands with diverse interactions.
            results_view_name (str): name for resulting book mark file. Default value is "passing_results"
            dict (dict): dictionary of one or more of the above args, is overwritten by individual args
        """
        # Dict of individual arguments
        indiv_options: dict = vars(); 
        del indiv_options["self"]; del indiv_options["dict"]

        # Create option object with default values if needed
        if not hasattr(self, "storageopts"): self.storageopts = StorageOptions()
            
        # Set options from dict if provided
        if dict is not None:
            for k,v in dict.items():
                setattr(self.storageopts, k, v) 

        # Set additional options from individual arguments
        #NOTE Will overwrite config file
        for k,v in indiv_options.items():
            if v is not None: setattr(self.storageopts, k, v)

        # Assign attributes to storage manager
        for k,v in self.storageopts.todict().items():
            setattr(self.storageman, k, v)

    def set_results_processing_options(self, store_all_poses = None,
                                            max_poses = None,
                                            add_interactions = None,
                                            interaction_tolerance = None,
                                            interaction_cutoffs = None,
                                            max_proc = None,
                                            dict: dict =None):
            
        #TODO problem with how options are assigned to manager
        """
        Create resultsmanager object if needed and set options.

        Args:
            store_all_poses (bool): store all ligand poses, does it take precedence over max poses? 
            max_poses (int): how many poses to save (ordered by soem score?)
            add_interactions (bool): add ligand-receptor interaction data, only in vina mode
            interaction_tolerance (float): longest ångström distance that is considered interaction?
            interaction_cutoffs (list): ångström distance cutoffs for x and y interaction
            max_proc (int): max number of computer processors to use for file reading
            dict (dict): dictionary of one or more of the above args, is overwritten by individual args
        """
        # Dict of individual arguments
        indiv_options: dict = vars(); 
        del indiv_options["self"]; del indiv_options["dict"]

        # Create option object with default values if needed
        if not hasattr(self, "resultsmanopts"): self.resultsmanopts = ResultsProcessingOptions()
            
        # Set options from dict if provided
        if dict is not None:
            for k,v in dict.items():
                setattr(self.resultsmanopts, k, v)

        # Set additional options from individual arguments
        #NOTE Will overwrite config file
        for k,v in indiv_options.items():
            if v is not None: setattr(self.resultsmanopts, k, v)
        
    def set_read_options(self, 
                         filtering = None,
                         plot = None,
                         find_similar_ligands = None,
                         export_bookmark_csv =  None,
                         export_bookmark_db =  None,
                         export_query_csv =  None,
                         export_receptor =  None,
                         data_from_bookmark =  None,
                         pymol =  None,
                         enumerate_interaction_combs =  None,
                         log_file = None,
                         export_sdf_path = None,
                         dict: dict =None):
        #TODO this requires some work because some of these options are really used to call functions
        """ Class that holds options related to reading from the database, including format for
        result export and alternate ways of displaying the data (plotting),

        Args:
            filtering (bool): implicit argument set if there are any optional filters present
            plot (bool): Makes scatterplot of LE vs Best Energy, saves as scatter.png.
            find_similar_ligands (str): Allows user to find similar ligands to given ligand name based on previously performed morgan fingerprint or interaction clustering.
            export_bookmark_csv (str): Create csv of the bookmark given with bookmark_name. Output as <bookmark_name>.csv. Can also export full database tables
            export_bookmark_db (bool): Export a database containing only the results found in the bookmark specified by --bookmark_name. Will save as <input_db>_<bookmark_name>.db
            export_query_csv (str): Create csv of the requested SQL query. Output as query.csv. MUST BE PRE-FORMATTED IN SQL SYNTAX e.g. SELECT [columns] FROM [table] WHERE [conditions]
            export_receptor (bool): Export stored receptor pdbqt. Will write to current directory.
            data_from_bookmark (bool): Write log of --outfields data for bookmark specified by --bookmark_name. Must use without any filters.
            pymol (bool): Lauch PyMOL session and plot of ligand efficiency vs docking score for molecules in bookmark specified with --bookmark_name. Will display molecule in PyMOL when clicked on plot. Will also open receptor if given.
            enumerate_interaction_combs (bool): #TODO
            log_file (str): by default, results are saved in "output_log.txt"; if this option is used, ligands and requested info passing the filters will be written to specified file
            export_sdf_path (str): specify the path where to save poses of ligands passing the filters (SDF format); if the directory does not exist, it will be created; if it already exist, it will throw an error, unless the --overwrite is used  NOTE: the log file will be automatically saved in this path. Ligands will be stored as SDF files in the order specified.
                optional: readopts (ReadOptions): provide options as object, will overwrite other options
            dict (dict): dictionary of one or more of the above args, is overwritten by individual args
        """
        # Dict of individual arguments
        indiv_options: dict = vars(); 
        del indiv_options["self"]; del indiv_options["dict"]

        # Create option object with default values if needed
        if not hasattr(self, "readopts"): self.readopts = ReadOptions()
            
        # Set options from dict if provided
        if dict is not None:
            for k,v in dict.items():
                setattr(self.readopts, k, v)

        # Set additional options from individual arguments
        #NOTE Will overwrite config file
        for k,v in indiv_options.items():
            if v is not None: setattr(self.readopts, k, v)

        # Creates output man with attributes if needed
        if hasattr(self, "outputman"):
            self.outputman.export_sdf_path = self.readopts.export_sdf_path
            self.outputman.log_file = self.readopts.log_file
        else: self.outputman = OutputManager(self.readopts.log_file, self.readopts.export_sdf_path)

    def set_filters(self,
                    eworst=None, 
                    ebest=None, 
                    leworst=None, 
                    lebest=None, 
                    score_percentile=None,
                    le_percentile=None,
                    vdw_interactions=None,
                    hb_interactions=None,
                    reactive_interactions=None,
                    interactions_count=None,
                    react_any=None,
                    max_miss=None,
                    ligand_name=None,
                    ligand_substruct=None,
                    ligand_substruct_pos=None,
                    ligand_max_atoms=None,
                    ligand_operator=None,     
                    dict: dict = None):
        """
        Object holding all specific filters, options, types.
        Args:
            #TODO
            dict (dict): dictionary of one or more of the above args, is overwritten by individual args
        """
        # Dict of individual arguments
        indiv_options: dict = vars(); 
        del indiv_options["self"]; del indiv_options["dict"]

        # Create option object with default values if needed
        if not hasattr(self, "filterobj"):
            self.filterobj = Filters()
            
        # Set options from dict if provided
        if dict is not None:
            for k,v in dict.items():
                if v is not None: setattr(self.filterobj, k, v) 

        # Set additional options from individual arguments
        #NOTE Will overwrite config file
        for k,v in indiv_options.items():
            if v is not None: setattr(self.filterobj, k, v)

        # Compatibility check with docking mode
        if self.generalopts.docking_mode == "vina" and self.filterobj.react_any:
            logger.warning("Cannot use reaction filters with Vina mode. Removing react_any filter.")
            self.filterobj.react_any = False
    
    def set_file_sources(self,
                    file=None, 
                    file_path=None, 
                    file_list=None, 
                    file_pattern=None, 
                    recursive=None, 
                    receptor_file=None,
                    save_receptor=None,
                    dict: dict = None):
        """
        Object holding all ligand docking results files.
        Args:
            dict (dict): dictionary of one or more of the above args, is overwritten by individual args
        """
        # Dict of individual arguments
        indiv_options: dict = vars(); 
        del indiv_options["self"]; del indiv_options["dict"]

        # Create option object with default values if needed
        if not hasattr(self, "file"): files = InputFiles()
            
        # Set options from dict if provided
        if dict is not None:
            for k,v in dict.items():
                setattr(files, k, v) 

        # Set additional options from individual arguments
        #NOTE Will overwrite config file
        for k,v in indiv_options.items():
            if v is not None: setattr(files, k, v)
        
        self.files = files

    #-#-#- Output API -#-#-#
    def create_ligand_rdkit_mol(self, ligname, smiles, atom_indices, h_parent_line, flexible_residues, flexres_atomnames, pose_ID=None, write_nonpassing=False):
        """creates rdkit molecule for given ligand, either for a specific pose_ID or for all passing (and nonpassing?) poses

        Args:
            ligname (string): ligand name
            smiles (string): ligand smiles string
            atom_indices (list): list of atom indices converting pdbqt to rdkit mol
            h_parent_line (list): list of atom indices for heteroatoms with attached hydrogens
            flexible_residues (list): list of flexible residue names
            flexres_atomnames (list): list of atomtypes in flexible residue
            pose_ID (int, optional): pose_ID for single pose to return. Defaults to None.
            write_nonpassing (bool, optional): _description_. Defaults to False.

        Raises:
            OutputError: raises error if there is an issue with determining a flexible residue identity

        Returns:
            tuple: (ligand rdkit mol, [flexres rdkit mols], {properties for ligand})
        """        
        mol = Chem.MolFromSmiles(smiles)
        flexres_mols = []
        flexres_info = []
        atom_indices = json.loads(atom_indices)
        ligand_saved_coords = []
        flexres_saved_coords = []
        # make flexible residue molecules
        for res, res_ats in zip(flexible_residues, flexres_atomnames):
            flexres_saved_coords.append([])
            resname = res[:3]
            res_ats = [
                at.strip() for at in res_ats
            ]  # strip out whitespace around atom names
            (
                res_smiles,
                res_index_map,
                res_h_parents,
            ) = RDKitMolCreate.guess_flexres_smiles(resname, res_ats)
            if res_smiles is None:  # catch error in guessing smiles
                raise OutputError(
                    f"Error while creating Mol for flexible residue {res}: unrecognized residue or incorrect atomtypes"
                )
            frm = Chem.MolFromSmiles(res_smiles)
            frm.SetProp("resinfo", res)
            flexres_mols.append(frm)
            flexres_info.append((res_smiles, res_index_map, res_h_parents))

        # fetch coordinates for passing poses and add to
        # rdkit ligand mol, add flexible residues
        properties = {
            "Binding energies": [],
            "Ligand effiencies": [],
            "Interactions": [],
        }
        with self.storageman:
            if pose_ID is None:  # get all passing and nonpassing poses (if requested)
                passing_properties = self.storageman.fetch_passing_pose_properties(
                    ligname
                )
                (
                    mol,
                    flexres_mols,
                    ligand_saved_coords,
                    flexres_saved_coords,
                    properties,
                ) = self._add_poses(
                    atom_indices,
                    passing_properties,
                    mol,
                    flexres_mols,
                    flexres_info,
                    ligand_saved_coords,
                    flexres_saved_coords,
                    properties,
                )

                # fetch coordinates for non-passing poses
                # and add to ligand mol, flexible residue mols
                if write_nonpassing:
                    nonpassing_properties = (
                        self.storageman.fetch_nonpassing_pose_properties(ligname)
                    )
                    (
                        mol,
                        flexres_mols,
                        ligand_saved_coords,
                        flexres_saved_coords,
                        properties,
                    ) = self._add_poses(
                        atom_indices,
                        nonpassing_properties,
                        mol,
                        flexres_mols,
                        flexres_info,
                        ligand_saved_coords,
                        flexres_saved_coords,
                        properties,
                    )
            else:
                pose_properties = (self.storageman.fetch_single_pose_properties(pose_ID))
                (
                    mol,
                    flexres_mols,
                    ligand_saved_coords,
                    flexres_saved_coords,
                    properties,
                ) = self._add_poses(
                    atom_indices,
                    pose_properties,
                    mol,
                    flexres_mols,
                    flexres_info,
                    ligand_saved_coords,
                    flexres_saved_coords,
                    properties,
                )

        # add hydrogens to mols
        lig_h_parents = [int(idx) for idx in json.loads(h_parent_line)]
        mol = RDKitMolCreate.add_hydrogens(
            mol, ligand_saved_coords, lig_h_parents
        )
        flexres_hparents = []
        for idx, res in enumerate(flexres_mols):
            flexres_hparents = flexres_info[idx][2]
            flexres_mols[idx] = RDKitMolCreate.add_hydrogens(
                res, flexres_saved_coords[idx], flexres_hparents
            )

        return mol, flexres_mols, properties

    def find_similar_ligands(self, query_ligname: str):
        """Find ligands in cluster with query_ligname
        """
        
        with self.storageman: similar_ligands, bookmark_name, cluster_name = self.storageman.fetch_clustered_similars(query_ligname)

        if similar_ligands is not None:
            with self.outputman:
                self.outputman.write_find_similar_header(query_ligname, cluster_name)
                self.outputman.write_results_bookmark_to_log(bookmark_name)
                number_similar = self.outputman.write_log(similar_ligands)
                self.outputman.log_num_passing_ligands(number_similar)
                print("Number similar ligands:", number_similar)

    def get_previous_filter_data(self):
        """Get data requested in self.out_opts['outfields'] from the
        results view of a previous filtering
        """
        with self.storageman: new_data = self.storageman.fetch_data_for_passing_results()
        with self.outputman: self.outputman.write_log(new_data)

    def plot(self, save=True):
        """
        Get data needed for creating Ligand Efficiency vs
        Energy scatter plot from storageManager. Call OutputManager to create plot.
        """
        if self.filter.max_miss > 0:
            raise OptionError("Cannot use --plot with --max_miss > 0. Can plot for desired bookmark with --bookmark_name.")
        
        logger.info("Creating plot of results")
        # get data from storageMan
        with self.storageman: all_data, passing_data = self.storageman.get_plot_data()
        all_plot_data_binned = dict()
        # bin the all_ligands data by 1000ths to make plotting faster
        for line in all_data:
            # add to dictionary as bin of energy and le
            if None in line:
                continue
                #raise OutputError("Detected empty data line when plotting. Please check that database and bookmarks are not empty.")
            data_bin = (round(line[0], 3), round(line[1], 3))
            if data_bin not in all_plot_data_binned:
                all_plot_data_binned[data_bin] = 1
            else:
                all_plot_data_binned[data_bin] += 1
        # plot the data
        self.outputman.plot_all_data(all_plot_data_binned)
        if passing_data != []:  # handle if no passing ligands
            for line in passing_data:
                self.outputman.plot_single_point(
                    line[0], line[1], "red"
                )  # energy (line[0]) on x axis, le (line[1]) on y axis
        if save:
            self.outputman.save_scatterplot()
        else:
            plt.show()

    def write_molecule_sdfs(self, write_nonpassing=False, return_rdmol_dict=False):
        """
        Have output manager write sdf molecules for passing results in given results bookmark

        Args:
            write_nonpassing (bool, optional): Option to include non-passing poses for passing ligands
            return_rdmol_dict (bool, optional): Suppresses SDF file writing, returns dictionary of rdkit mols
        """

        with self.storageman:
            if self.filterobj.max_miss > 0:
                logger.warning("WARNING: Requested --export_sdf_path with --max_miss. Exported SDFs will be for union of interaction combinations.")
                self.storageman.results_view_name = self.storageman.results_view_name + "_union"
            if not self.storageman.check_passing_view_exists():
                logger.warning(
                    "Given results bookmark does not exist in database. Cannot write passing molecule SDFs"
                )
                return None
            # make temp table for SDF writing
            self.storageman.create_temp_passing_table()
            passing_molecule_info = self.storageman.fetch_passing_ligand_output_info()
            flexible_residues, flexres_atomnames = self.storageman.fetch_flexres_info()

        if flexible_residues != []:
            flexible_residues = json.loads(flexible_residues)
            flexres_atomnames = json.loads(flexres_atomnames)
        all_mols = {}
        for (ligname, smiles, atom_indices, h_parent_line) in passing_molecule_info:
            logger.info("Writing " + ligname.split(".")[0] + ".sdf")
            # create rdkit ligand molecule and flexible residue container
            if smiles == "":
                logger.warning(
                    f"No SMILES found for {ligname}. Cannot create SDF."
                )
                continue
            mol, flexres_mols, properties = self.create_ligand_rdkit_mol(ligname, smiles, atom_indices, h_parent_line, flexible_residues, flexres_atomnames, write_nonpassing=write_nonpassing)

            # write out mol
            if not return_rdmol_dict:
                self.outputman.write_out_mol(
                    ligname, mol, flexres_mols, properties
                )
            else:
                all_mols[ligname] = {"ligand": mol, "flex_residues": flexres_mols}

        if return_rdmol_dict:
            return all_mols

    def display_pymol(self):
        """
        Launch pymol session and plot of LE vs docking score. Displays molecules when clicked
        """

        import subprocess
        from rdkit.Chem import PyMol

        # launch pymol session
        p = subprocess.Popen(
        ["pymol", "-R"],
        stdout=subprocess.PIPE,
        )

        poseIDs = {}
        with self.storageman:
            # fetch data for passing ligands
            _, passing_data = self.storageman.get_plot_data(only_passing=True)
            for line in passing_data:
                plt.plot(line[0], line[1], '.r', mfc='None', picker=5)
                poseIDs[(line[0], line[1])] = (line[2], line[3])  # line[0] is LE, line[1] is docking score, line[2] is pose_ID, line[3] is LigName
            plt.ylabel("Ligand Efficiency (kcal/mol/heavy atom)")
            plt.xlabel("Docking Score (kcal/mol)")
            plt.title("Passing Docking Poses")

            try:
                pymol = PyMol.MolViewer()
            except ConnectionRefusedError as e:
                raise RTCoreError("Error establishing connection with PyMol. Try manually launching PyMol with `pymol -R` in another terminal window.") from e
            
            def onpick(event):
                line = event.artist
                coords = tuple([c[0] for c in line.get_data()])
                chosen_pose = poseIDs[coords]
                logger.info(f"LigName: {chosen_pose[1]}; Pose_ID: {chosen_pose[0]}")

                # make rdkit mol for poseid
                ligname, ligand_smile, atom_index_map, hydrogen_parents = self.storageman.fetch_single_ligand_output_info(chosen_pose[1])
                flexible_residues, flexres_atomnames = self.storageman.fetch_flexres_info()
                if flexible_residues != []:  # converts string to list
                    flexible_residues = json.loads(flexible_residues)
                    flexres_atomnames = json.loads(flexres_atomnames)

                mol, flexres_mols, _ = self.create_ligand_rdkit_mol(ligname, ligand_smile, atom_index_map, hydrogen_parents, flexible_residues, flexres_atomnames, pose_ID=chosen_pose[0])
                logger.debug(Chem.MolToSmiles(mol))
                pymol.ShowMol(mol, name=ligname, showOnly=False)
                for idx, resmol in enumerate(flexres_mols):
                    pymol.ShowMol(resmol, name=ligname + "_" + flexible_residues[idx], showOnly=False)

            fig = plt.gcf()
            cid = fig.canvas.mpl_connect('pick_event', onpick)
            plt.show()

    def export_csv(self, requested_data: str, csv_name: str, table=False):
        """Get requested data from database, export as CSV

        Args:
            requested_data (string): Table name or SQL-formatted query
            csv_name (string): Name for exported CSV file
            table (bool): flag indicating is requested data is a table name
        """
        with self.storageman: 
            df = self.storageman.to_dataframe(requested_data, table=table)
            df.to_csv(csv_name)
        
    def export_bookmark_db(self, bookmark_db_name: str):
        """Export database containing data from bookmark

        Args:
            bookmark_db_name (str): name for bookmark_db
        """
        
        logger.info("Exporting bookmark database")
        if os.path.exists(bookmark_db_name):
            logger.warning(
                "Requested export DB name already exists. Please rename or remove existing database. New database not exported."
            )
            return
        with self.storageman: self.storageman.clone(bookmark_db_name)
        # connect to cloned database
        self.db_file = bookmark_db_name
        with StorageManagerSQLite(**self.storageopts) as db_clone:
            db_clone.prune()
            db_clone.close_storage(vacuum=True)

    def export_receptors(self):
        """
        Export receptor to pdbqt
        """
        with self.storageman: receptor_tuples = self.storageman.fetch_receptor_objects()
        for recname, recblob in receptor_tuples:
            if recblob is None:
                logger.warning(f"No receptor pdbqt stored for {recname}. Export failed.")
                continue
            self.outputman.write_receptor_pdbqt(recname, recblob)

    def get_previous_filter_data(self):
        """Get data requested in self.out_opts['outfields'] from the
        results view of a previous filtering
        """
        with self.storageman: new_data = self.storageman.fetch_data_for_passing_results()
        with self.outputman: self.outputman.write_log(new_data)
    
    @staticmethod
    def generate_options_json_template(to_file=True):
        """
        Creates a dict of all Ringtail option classes, and their 
        key-default value pairs. Outputs to options.json in 
        "util_files" of to_file = true, else it returns the dict 
        of default option values.
        Args:
            to_file (bool): if true writes file to standard options json path, if false returns a json string with values
        """
        readopts = ReadOptions().todict()
        generalopts = GeneralOptions().todict()
        fileobj = InputFiles().todict()
        writeopts = ResultsProcessingOptions().todict()
        storageopts = StorageOptions().todict()
        filterobj = Filters().todict()

        json_string = {"generalopts": generalopts,
                       "writeopts": writeopts,
                       "storageopts": storageopts,
                       "readopts": readopts,
                       "filterobj": filterobj,
                       "fileobj": fileobj}
        if to_file:
            filepath= RingtailCore._options_file_path()
            with open(filepath, 'w') as f: 
                f.write(json.dumps(json_string, indent=4))
            logger.debug(f"Default ringtail option values written to file {filepath}")
            return filepath
        else:
            logger.debug("Default ringtail option values prepared as a string.")
            return json_string

    @staticmethod       
    def get_defaults(object="all") -> dict:
        """
        Gets default values from RingtailOptions and returns dict of all,
        or specific object. 
        
        Args:
            object (str): ["all", "generalopts", "writeopts", "storageopts", "readopts", "filterobj", "fileobj"]
        """
        all_defaults = RingtailCore.generate_options_json_template(to_file=False)

        if object.lower() not in ["all", "generalopts", "writeopts", "storageopts", "readopts", "filterobj", "fileobj"]:
            raise OptionError(f'The options object {object.lower()} does not exist. Please choose amongst \n ["all", "generalopts", "writeopts", "storageopts", "readopts", "filterobj", "fileobj"]')
        
        if object.lower() == "all":
            logger.debug("All ringtail default values have been fetched.")
            return all_defaults
        else:
            logger.debug(f"Ringtail default values for {object} have been fetched.")
            return all_defaults[object.lower()]

    @staticmethod
    def _options_file_path(filename="options.json"):
        utilfolder = path.abspath(__file__ + "/../../util_files/")
        if not os.path.exists(utilfolder):
            os.makedirs(utilfolder) 
        return utilfolder + "/" + filename

