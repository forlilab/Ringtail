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
    including adding results to storage, filtering and clusteirng, and outputting data as
    rdkit molecules, plotting docking results, and visualizing select ligands in pymol.

    Attributes:
        db_file (str): name of database file being operated on
        docking_mode (str): 
        storageman (StorageManager object): Interface module with database
        resultsman (ResultsManager object): Module to deal with results processing before adding to database
        outputman (OutputManager object): Manager for output tasks of log-writting, plotting, ligand SDF writing, starting pymol sessions
        filters (Filters object): object holding all optional filters
    """
#-#-#- Base methods -#-#-#
    
    def __init__(self, 
                 db_file: str = "output.db", 
                 storage_type: str = "sqlite", 
                 logging_level: str = None):
        """
        Initialize RingtailCore object and create a storageman object with the db file.
        Does not open access to the storage. Future option will include opening database as readonly.

        _run_mode refers to whether ringtail is ran from the command line or through direct API use, 
        where the former is more restrictive. 
        """
        if logging_level is not None: logger.setLevel(logging_level)
        self.db_file = db_file
        storageman = StorageManager.check_storage_compatibility(storage_type) 
        self.storageman = storageman(db_file)
        self.set_storageman_attributes()
        self._run_mode = "api"
        self._docking_mode = "dlg"

    def update_database_version(self, consent=False):
        # Method to update database version from 1.0.0 to 1.1.0
        return self.storageman.update_database_version(consent)
    

#-#-#- Private methods -#-#-#
    
    def _validate_docking_mode(self, docking_mode: str):
        """ Method that validates specified AutoDock program used to generate results.
        Args:
            docking_mode (str): string that describes docking mode
        Raises:
            RTCoreError if docking_mode is not supported
        """
        if type(docking_mode) is not str:
            logger.warning('The given docking mode was not given as a string, it will be set to default value "dlg".')
            self._docking_mode = "dlg"
        elif docking_mode.lower() not in ["dlg", "vina"]:
            raise NotImplementedError(f'Docking mode {docking_mode} is not supported. Please choose between "dlg" and "vina".')
        else:
            self._docking_mode = docking_mode.lower()
            logger.debug(f"Docking mode set to {self.docking_mode}.")
    
    def _get_docking_mode(self):
        return self._docking_mode

    docking_mode = property(fget = _get_docking_mode, fset = _validate_docking_mode)

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

    def _generate_interaction_combinations(self, max_miss=0):
        """Recursive function to list of tuples of possible interaction filter combinations, excluding up to max_miss interactions per filtering round

        Args:
            max_miss (int): Maximum number of interactions to be excluded
        """

        all_interactions = []
        for _type in Filters.get_filter_keys("interaction"):
            interactions = getattr(self.filters, _type)
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

        filters_dict = self.filters.todict()
        for itype in Filters.get_filter_keys("interaction"):
            itype_interactions = filters_dict[itype]
            for interaction in itype_interactions:
                if itype + "-" + interaction[0] not in interaction_combination:
                    filters_dict[itype].remove(interaction)

        return filters_dict

    def _create_rdkit_mol(self, ligname, smiles, atom_indices, h_parent_line, flexible_residues, flexres_atomnames, pose_ID=None, write_nonpassing=False):
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
        
        Note: needs to be ran inside a storageman context manager, will not be able to access the temporary table otherwise.
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

    def _set_file_sources(
            self,
            file=None, 
            file_path=None, 
            file_list=None, 
            file_pattern=None, 
            recursive=None, 
            receptor_file=None,
            save_receptor=None,
            dict: dict = None) -> InputFiles:
        """
        Object holding all ligand docking results files.
        Args:
            file (str, optional: list(str)): ligand result file
            file_path (str, optional: list(str)): list of folders containing one or more result files
            file_list (str, optional: list(str)): list of ligand result file(s)
            file_pattern (str): file pattern to use with recursive search in a file_path, "*.dlg*" for AutoDock-GDP and "*.pdbqt*" for vina
            recursive (bool): used to recursively search file_path for folders inside folders
            receptor_file (str): string containing the receptor .pdbqt
            save_receptor (bool): whether or not to store the full receptor details in the database (needed for some things)
            dict (dict): dictionary of one or more of the above args, is overwritten by individual args
        Returns:
            InputFiles object
        """

        def ensure_double_list(object) -> list:
            """Most of ringtail is set up to handle files, file paths, and file lists as double lists [[items]].
            Instead of changing that for now, the input is checked to ensure they are presented as, or can be 
            converted to a double list. Really only matters for using API."""
            if type(object)==list:
                if type(object[0])==list: 
                    if type(object[0][0])==str: pass
                    else: raise OptionError( f"error, object is more than two encapsulated lists: '{object[0][0]}' should be a string.")
                elif type(object[0])==str:
                    object=[object]
                else: logger.error("Unable to parse file input.")
            elif type(object)==str:
                object=[[object]]
            else: logger.error("Unable to parse file input.")
            
            return object
        
        # Ensure files are in current writeable format
        if file is not None: file = ensure_double_list(file)
        if file_path is not None: file_path = ensure_double_list(file_path)
        if file_list is not None: file_list = ensure_double_list(file_list)

        # Set file format 
        if file_pattern is None:
            if file is not None and "pdbqt" in file[0][0]:
                file_pattern = "*.pdbqt*"
            else:
                file_pattern = "*.dlg*"
                logger.warning("File pattern was not specified, set to default '*.dlg*'.")

        # Dict of individual arguments
        indiv_options: dict = vars(); 
        del indiv_options["self"]; del indiv_options["dict"]; del indiv_options["ensure_double_list"]

        # Create option object with default values if needed
        files = InputFiles()
            
        # Set options from dict if provided
        if dict is not None:
            for k,v in dict.items():
                setattr(files, k, v) 
                logger.debug(f'File attribute {k} was set to {v}.')

        # Set additional options from individual arguments
        #NOTE Will overwrite config file
        for k,v in indiv_options.items():
            if v is not None: 
                setattr(files, k, v)
                logger.debug(f'File attribute {k} was set to {v}.')
        
        # set docking mode based on file pattern
        if files.file_pattern != None: 
            if "pdbqt" in files.file_pattern.lower() and self.docking_mode != "vina": self.docking_mode = "vina"
            elif "dlg" in files.file_pattern.lower() and self.docking_mode != "dlg": self.docking_mode = "dlg"

        return files
    
    def _set_results_sources(
            self,
            results_strings: dict = None,
            receptor_file: str =None,
            save_receptor: bool =None,
            dict: dict = None) -> InputStrings:
        """
        Object holding all ligand vina docking results string and corresponding receptor.
        Args:
            results_string (dict): string containing the ligand identified and docking results as a dictionary
            receptor_file (str): string containing the receptor .pdbqt
            save_receptor (bool): whether or not to store the full receptor details in the database (needed for some things)
            dict (dict): dictionary of one or more of the above args, is overwritten by individual args
        Return:
            InputStrings object
        """
        # Dict of individual arguments
        indiv_options: dict = vars(); 
        del indiv_options["self"]; del indiv_options["dict"]

        # Create option object with default values if needed
        strings = InputStrings()
            
        # Set options from dict if provided
        if dict is not None:
            for k,v in dict.items():
                setattr(strings, k, v) 
                logger.debug(f'Docking string results attribute {k} was set to {v}.')

        # Set additional options from individual arguments
        #NOTE Will overwrite config file
        for k,v in indiv_options.items():
            if v is not None: setattr(strings, k, v)
            logger.debug(f'Docking string results attribute {k} was set to {v}.')
        
        return strings

    def _create_resultsmanager(self, file_sources: InputFiles = None, string_sources: InputStrings = None) -> ResultsManager:
        """Creates a results manager object based on results provided either as files or as strings (currently only for vina).
        Will create a new object each time results are added.
        In its current state it assumes only one source of results will be provided. 
        If both are provided, it will only process the strings. 
        Args:
            file_sources (InputFiles): used if docking results are provided through files
            string_sources (InputStrings): used if docking results are provided through strings
        """
        if file_sources is not None:
            self.resultsman = ResultsManager(file_sources=file_sources)
            logger.debug("Results manager object has been created with results files.")
        elif string_sources is not None:
            self.resultsman = ResultsManager(string_sources=string_sources)
            logger.debug("Results manager object has been created with results strings.")
        else:
            raise RTCoreError("No results sources were provided, a results manager object could not be created.")


#-#-#- Core attribute setting methods -#-#-#
    """ These methods are used internally to assing values to all ringtail options. 
    This ensures:   - that all options are set to specific types through RingtailOptions
                    - that internal consistency checks are performed on a group of options
                    - these methods ensure options are assigned to the appropriate ringtail manager classes""" 
    
    def set_storageman_attributes(self, 
                            filter_bookmark: str = None,
                            append_results: bool = None,
                            duplicate_handling: str = None,
                            overwrite: bool = None,
                            order_results: str = None,
                            outfields: str = None,
                            output_all_poses: str = None,
                            mfpt_cluster: float = None,
                            interaction_cluster: float = None,
                            bookmark_name: str = None,
                            dict: dict = None):
        """
        Create storage_manager_options object if needed, sets options, and assigns them to the storage manager object.

        Args:
            filter_bookmark (str): Perform filtering over specified bookmark. (in output group in CLI)
            append_results (bool): Add new results to an existing database, specified by database choice in ringtail initialization or --input_db in cli
            duplicate_handling (str, options): specify how duplicate Results rows should be handled when inserting into database. Options are "ignore" or "replace". Default behavior will allow duplicate entries.
            overwrite (bool): by default, if a log file exists, it doesn't get overwritten and an error is returned; this option enable overwriting existing log files. Will also overwrite existing database
            order_results (str): Stipulates how to order the results when written to the log file. By default will be ordered by order results were added to the database. ONLY TAKES ONE OPTION."
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
            bookmark_name (str): name for resulting book mark file. Default value is "passing_results"
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
                if v is not None: setattr(self.storageopts, k, v) 
                logger.debug(f'Storage manager attribute {k} was set to {v}.')

        # Set additional options from individual arguments
        #NOTE Will overwrite config file
        for k,v in indiv_options.items():
            if v is not None: setattr(self.storageopts, k, v)
            logger.debug(f'Storage manager attribute {k} was set to {v}.')

        # Assign attributes to storage manager
        for k,v in self.storageopts.todict().items():
            setattr(self.storageman, k, v)
        logger.debug("Options for storage manager have been changed.")      

    def set_resultsman_attributes(self,
                                    store_all_poses: bool = None,
                                    max_poses: int = None,
                                    add_interactions: bool = None,
                                    interaction_tolerance: float = None,
                                    interaction_cutoffs = None,
                                    max_proc: int = None,
                                    dict: dict = None):
            
        """
        Create results_manager_options object if needed, sets options, and assigns them to the results manager object.

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
                logger.debug(f'Results manager attribute {k} was set to {v}.')

        # Set additional options from individual arguments
        #NOTE Will overwrite config file
        for k,v in indiv_options.items():
            if v is not None: 
                setattr(self.resultsmanopts, k, v)
            logger.debug(f'Results manager attribute {k} was set to {v}.')

        # Assigns options to the results manager object
        for k,v in self.resultsmanopts.todict().items():  
            if v is not None: setattr(self.resultsman, k, v)    
        logger.debug("Options for results manager have been changed.")       
    
    def set_output_options(self, 
                         log_file: str = None,
                         export_sdf_path: str = None,
                         enumerate_interaction_combs: bool =  None,
                         find_similar_ligands: str = None,
                         export_bookmark_csv: str =  None,
                         export_query_csv: str =  None,
                         dict: dict =None):
        """ Creates output options object that holds attributes related to reading and outputting results.
        Will assign log_file name and export_sdf_path to the output_manager object.

        Args:
            log_file (str): by default, results are saved in "output_log.txt"; if this option is used, ligands and requested info passing the filters will be written to specified file
            export_sdf_path (str): specify the path where to save poses of ligands passing the filters (SDF format); if the directory does not exist, it will be created; if it already exist, it will throw an error, unless the --overwrite is used  NOTE: the log file will be automatically saved in this path. Ligands will be stored as SDF files in the order specified.
            enumerate_interaction_combs (bool): When used with `max_miss` > 0, will log ligands/poses passing each separate interaction filter combination as well as union of combinations. Can significantly increase runtime.
            find_similar_ligands (str): Allows user to find similar ligands to given ligand name based on previously performed morgan fingerprint or interaction clustering.
            export_bookmark_csv (str): Create csv of the bookmark given with bookmark_name. Output as <bookmark_name>.csv. Can also export full database tables
            export_query_csv (str): Create csv of the requested SQL query. Output as query.csv. MUST BE PRE-FORMATTED IN SQL SYNTAX e.g. SELECT [columns] FROM [table] WHERE [conditions]
            dict (dict): dictionary of one or more of the above args, is overwritten by individual args
        """
        # Dict of individual arguments
        indiv_options: dict = vars(); 
        del indiv_options["self"]; del indiv_options["dict"]

        # Create option object with default values if needed
        if not hasattr(self, "outputopts"): self.outputopts = OutputOptions()
            
        # Set options from dict if provided
        if dict is not None:
            for k,v in dict.items():
                if v is not None: setattr(self.outputopts, k, v)
                logger.debug(f'Output options {k} was set to {v}.')

        # Set additional options from individual arguments
        #NOTE Will overwrite config file
        for k,v in indiv_options.items():
            if v is not None: setattr(self.outputopts, k, v)
            logger.debug(f'Output options {k} was set to {v}.')

        # Creates output man with attributes if needed
        self.outputman = OutputManager(self.outputopts.log_file, self.outputopts.export_sdf_path)
        logger.debug("Options for output manager have been changed.")   

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
        Create a filter object containing all numerical and string filters. 
        """

        # Dict of individual arguments
        indiv_options: dict = vars(); 
        del indiv_options["self"]; del indiv_options["dict"]

        # Create a filter object
        self.filters = Filters()
            
        # Set options from dict if provided
        if dict is not None:
            for k,v in dict.items():
                if v is not None: setattr(self.filters, k, v) 
                logger.debug(f'Filter {k} was set to {v}.')

        # Set additional options from individual arguments
        #NOTE Will overwrite config file
        for k,v in indiv_options.items():
            if v is not None: setattr(self.filters, k, v)
            logger.debug(f'Filter {k} was set to {v}.')
    

    #-#-#- API -#-#-#
    def add_results_from_files(self,
                               file: str = None, 
                               file_path: str = None, 
                               file_list: str = None, 
                               file_pattern: str = None, 
                               recursive: bool = None, 
                               receptor_file: str = None,
                               save_receptor: bool = None,
                               filesources_dict: dict = None,
                               append_results: bool = None,
                               # results processing options
                               duplicate_handling: str = None,
                               overwrite: bool = None,
                               store_all_poses: bool = None,
                               max_poses: int = None,
                               add_interactions: bool = None,
                               interaction_tolerance: float = None,
                               interaction_cutoffs: list = None,
                               max_proc: int = None,
                               summary: bool = None,
                               options_dict: dict =None
                               ):
        """
        Call storage manager to process result files and add to database. Creates a database, or adds to an existing one if using "append_results".
        Options can be provided as a dict or as individual options. If both are provided, individual options will overwrite those from the dictionary. 

        Args:
            file (str, optional: list(str)): ligand result file
            file_path (str, optional: list(str)): list of folders containing one or more result files
            file_list (str, optional: list(str)): list of ligand result file(s)
            file_pattern (str): file pattern to use with recursive search in a file_path, "*.dlg*" for AutoDock-GDP and "*.pdbqt*" for vina
            recursive (bool): used to recursively search file_path for folders inside folders
            receptor_file (str): string containing the receptor .pdbqt
            save_receptor (bool): whether or not to store the full receptor details in the database (needed for some things)
            filesources_dict (dict): file sources already as an object 
            append_results (bool): Add new results to an existing database, specified by database choice in ringtail initialization or --input_db in cli
            duplicate_handling (str, options): specify how duplicate Results rows should be handled when inserting into database. Options are "ignore" or "replace". Default behavior will allow duplicate entries.
            store_all_poses (bool): store all ligand poses, does it take precedence over max poses? 
            max_poses (int): how many poses to save (ordered by soem score?)
            add_interactions (bool): add ligand-receptor interaction data, only in vina mode
            interaction_tolerance (float): longest ångström distance that is considered interaction?
            interaction_cutoffs (list): ångström distance cutoffs for x and y interaction
            max_proc (int): max number of computer processors to use for file reading
            options_dict (dict): write options as a dict
        """

        files = self._set_file_sources(file, file_path, file_list, file_pattern, recursive, receptor_file, save_receptor, filesources_dict)
        results_files_given = (files.file is not None or files.file_path is not None or files.file_list is not None)
        if not results_files_given and not files.save_receptor:
            raise OptionError("At least one input option needs to be used: --file, --file_path, --file_list, or --input_db and --save_receptor")

        # if dictionary of options provided, attribute to appropriate managers
        if options_dict is not None:
            results_dict, storage_dict = RingtailCore.split_dict(options_dict, ["append_results", "duplicate_handling", "overwrite"])
        else:
            storage_dict = None
            results_dict = None
        self.set_storageman_attributes(append_results=append_results, duplicate_handling=duplicate_handling, overwrite=overwrite, dict=storage_dict)

        # If there are ligand files present, process ligand data
        if results_files_given: 
            with self.storageman:

                # check storage exist and can be appended to if specified
                if self.storageopts.append_results and not RTOptions.is_valid_path(self.db_file):
                    raise OptionError("The provided --input_db is not a valid path, please check the provided path.")
                
                # Prepare the results manager 
                self._create_resultsmanager(file_sources=files)
                self.resultsman.storageman = self.storageman               
                self.resultsman.storageman_class = self.storageman.__class__
                self.set_resultsman_attributes(store_all_poses, max_poses, add_interactions, interaction_tolerance, interaction_cutoffs, max_proc, results_dict)

                # Docking mode compatibility check
                if self.docking_mode == "vina" and self.resultsman.interaction_tolerance is not None:
                    logger.warning("Cannot use interaction_tolerance with Vina mode. Removing interaction_tolerance.")
                    self.resultsman.interaction_tolerance = None
                self.resultsman.mode = self.docking_mode

                # Process results files and handle database versioning 
                self.storageman.check_storage_ready(self._run_mode, self.docking_mode, self.resultsman.store_all_poses, self.resultsman.max_poses)
                logger.info("Adding results...")
                self.resultsman.process_files()
                self.storageman.set_ringtail_db_schema_version()
                if summary: self.produce_summary()

        if files.save_receptor: 
            self.save_receptor(files.receptor_file)
    
    def add_results_from_vina_string(self,
                                    results_strings: dict = None,
                                    receptor_file: str = None,
                                    save_receptor: bool = None,
                                    resultsources_dict: dict = None,
                                    append_results: bool = None,
                                    # result processing options
                                    duplicate_handling: str = None,
                                    overwrite: bool = None,
                                    store_all_poses: bool = None,
                                    max_poses: int = None,
                                    add_interactions: bool = None,
                                    interaction_cutoffs: list = None,
                                    max_proc: int = None,
                                    summary: bool = None,
                                    options_dict=None
                                    ):
        """
        Call storage manager to process the given vina output string and add to database.
        Options can be provided as a dict or as individual options.
        Creates a database, or adds to an existing one if using "append_results".
        #TODO a lot of overlap with add_results_from_files
        Args:
            results_string (dict): string containing the ligand identified and docking results as a dictionary
            receptor_file (str): string containing the receptor .pdbqt
            save_receptor (bool): whether or not to store the full receptor details in the database (needed for some things)
            resultsources_dict (dict): file sources already as an object 

            append_results (bool): Add new results to an existing database, specified by database choice in ringtail initialization or --input_db in cli
            duplicate_handling (str, options): specify how duplicate Results rows should be handled when inserting into database. Options are "ignore" or "replace". Default behavior will allow duplicate entries.
            store_all_poses (bool): store all ligand poses, does it take precedence over max poses? 
            max_poses (int): how many poses to save (ordered by soem score?)
            add_interactions (bool): add ligand-receptor interaction data, only in vina mode
            interaction_cutoffs (list): ångström distance cutoffs for x and y interaction
            max_proc (int): max number of computer processors to use for file reading
            options_dict (dict): write options as a dict

        """
        # Method currently only works with vina output, set automatically
        self.docking_mode = "vina"

        # create results string object
        results = self._set_results_sources(results_strings, receptor_file, save_receptor, resultsources_dict)
        results_strings_given = bool(results.results_strings)
        if not results_strings_given and not results.save_receptor:
            raise OptionError("At least one input option needs to be used: 'results_strings', or 'save_receptor'")

        # if dictionary of options provided, attribute to appropriate managers
        if options_dict is not None:
            results_dict, storage_dict = RingtailCore.split_dict(options_dict, ["append_results", "duplicate_handling", "overwrite"])
        else:
            storage_dict = None
            results_dict = None
        self.set_storageman_attributes(append_results=append_results, duplicate_handling=duplicate_handling, overwrite=overwrite, dict=storage_dict)
        # If there are any ligand files, process ligand data
        if results_strings_given: 
            with self.storageman:
                if self.storageopts.append_results and not RTOptions.is_valid_path(self.db_file):
                    raise OptionError("The provided 'input_db' is not a valid path, please check the provided path.")
                # Prepare the results manager 
                
                self._create_resultsmanager(string_sources=results)
                self.resultsman.storageman = self.storageman               
                self.resultsman.storageman_class = self.storageman.__class__
                self.set_resultsman_attributes(store_all_poses, max_poses, add_interactions, None, interaction_cutoffs, max_proc, results_dict)

                self.resultsman.mode = self.docking_mode

                # Process results files and handle database versioning 
                self.storageman.check_storage_ready(self._run_mode, self.docking_mode, self.resultsman.store_all_poses, self.resultsman.max_poses)
                logger.info("Adding results...")
                self.resultsman.process_strings()
                self.storageman.set_ringtail_db_schema_version()
                if summary: self.produce_summary()

        if results.save_receptor: 
            self.save_receptor(results.receptor_file)

    def save_receptor(self, receptor_file):
            """
            Add receptor to database. Context managed by self.storageman

            Args:
                receptors (list): list of receptor blobs to add to database
                * currently only one receptor allowed per database
            """
            receptor_list = ReceptorManager.make_receptor_blobs([receptor_file])
            with self.storageman:
                for rec, _ in receptor_list:
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
    
    def produce_summary(self, columns=["docking_score", "leff"], percentiles=[1, 10]) -> None:
        """Print summary of data in storage
        Args:
            columns (list(str)): data columns used to prepare summary
            percentiles (list(int)): cutoff percentiles for the summary
        """
        with self.storageman: summary_data = self.storageman.fetch_summary_data(columns, percentiles)
        if 'summary_data' not in locals():
            raise RTCoreError(f'Summary data is empty, please check the database.')
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

    def filter(self, 
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
               filters_dict: dict = None,
               # other processing options: 
               enumerate_interaction_combs=False, 
               output_all_poses: bool = None, 
               mfpt_cluster = None, 
               interaction_cluster = None, 
               log_file: str = None, 
               overwrite: bool = None, 
               order_results: str = None, 
               outfields: str = None, 
               bookmark_name: str =None, 
               options_dict: dict = None,
               return_iter=False):
        """Prepare list of filters, then hand it off to storageman to
            perform filtering. Create log of passing results.
        Args:
            Filters:
                eworst (float): specify the worst energy value accepted
                ebest (float): specify the best energy value accepted
                leworst (float): specify the worst ligand efficiency value accepted
                lebest (float): specify the best ligand efficiency value accepted
                score_percentile (float): specify the worst energy percentile accepted. Express as percentage e.g. 1 for top 1 percent.
                le_percentile (float): specify the worst ligand efficiency percentile accepted. Express as percentage e.g. 1 for top 1 percent.
                vdw_interactions (list[tuple]): define van der Waals interactions with residue as [-][CHAIN]:[RES]:[NUM]:[ATOM_NAME]. E.g., [('A:VAL:279:', True), ('A:LYS:162:', True)] -> [('chain:resname:resid:atomname', <wanted (bool)>), ('chain:resname:resid:atomname', <wanted (bool)>)]
                hb_interactions (list[tuple]): define HB (ligand acceptor or donor) interaction as [-][CHAIN]:[RES]:[NUM]:[ATOM_NAME]. E.g., [('A:VAL:279:', True), ('A:LYS:162:', True)] -> [('chain:resname:resid:atomname', <wanted (bool)>), ('chain:resname:resid:atomname', <wanted (bool)>)]
                reactive_interactions (list[tuple]): check if ligand reacted with specified residue as [-][CHAIN]:[RES]:[NUM]:[ATOM_NAME]. E.g., [('A:VAL:279:', True), ('A:LYS:162:', True)] -> [('chain:resname:resid:atomname', <wanted (bool)>), ('chain:resname:resid:atomname', <wanted (bool)>)]
                interactions_count (list[tuple]): accept ligands with at least the requested number of HB interactions. If a negative number is provided, then accept ligands with no more than the requested number of interactions. E.g., [('hb_count', 5)]
                react_any (bool): check if ligand reacted with any residue
                max_miss (int): Will compute all possible combinations of interaction filters excluding up to max_miss numer of interactions from given set. Default will only return union of poses interaction filter combinations. Use with 'enumerate_interaction_combs' for enumeration of poses passing each individual combination of interaction filters.
                ligand_name (list[str]): specify ligand name(s). Will combine name filters with OR, e.g., ["lig1", "lig2"]
                ligand_substruct (list[str]): MARTS, index of atom in SMARTS, cutoff dist, and target XYZ coords, e.g., ["ccc", "CN"]
                ligand_substruct_pos (list[str]): SMARTS pattern(s) for substructure matching, e.g., ['"[Oh]C" 0 1.2 -5.5 10.0 15.5'] -> ["smart_string index_of_positioned_atom cutoff_distance x y z"]
                ligand_max_atoms (int): Maximum number of heavy atoms a ligand may have
                ligand_operator (str): logical join operator for multiple SMARTS (default: OR), either AND or OR
                filters_dict (dict): provide filters as a dictionary
            Ligand results options:
                enumerate_interaction_combs (bool): When used with `max_miss` > 0, will log ligands/poses passing each separate interaction filter combination as well as union of combinations. Can significantly increase runtime.
                output_all_poses (bool): By default, will output only top-scoring pose passing filters per ligand. This flag will cause each pose passing the filters to be logged.
                mfpt_cluster (float): Cluster filtered ligands by Tanimoto distance of Morgan fingerprints with Butina clustering and output ligand with lowest ligand efficiency from each cluster. Default clustering cutoff is 0.5. Useful for selecting chemically dissimilar ligands.
                interaction_cluster (float): Cluster filtered ligands by Tanimoto distance of interaction fingerprints with Butina clustering and output ligand with lowest ligand efficiency from each cluster. Default clustering cutoff is 0.5. Useful for enhancing selection of ligands with diverse interactions.
                log_file (str): by default, results are saved in "output_log.txt"; if this option is used, ligands and requested info passing the filters will be written to specified file
                overwrite (bool): by default, if a log file exists, it doesn't get overwritten and an error is returned; this option enable overwriting existing log files. Will also overwrite existing database
                order_results (str): Stipulates how to order the results when written to the log file. By default will be ordered by order results were added to the database. ONLY TAKES ONE OPTION.
                                        available fields are:  
                                        "e" (docking_score), 
                                        "le" (ligand efficiency), 
                                        "delta" (delta energy from best pose), 
                                        "ref_rmsd" (RMSD to reference pose), 
                                        "e_inter" (intermolecular energy), 
                                        "e_vdw" (van der waals energy), 
                                        "e_elec" (electrostatic energy), 
                                        "e_intra" (intermolecular energy), 
                                        "n_interact" (number of interactions), 
                                        "rank" (rank of ligand pose), 
                                        "run" (run number for ligand pose), 
                                        "hb" (hydrogen bonds); 
                outfields (str): defines which fields are used when reporting the results (to stdout and to the log file); fields are specified as comma-separated values, e.g. "--outfields=e,le,hb"; by default, docking_score (energy) and ligand name are reported; ligand always reported in first column available fields are: \n
                                    "Ligand_name" (Ligand name), 
                                    "e" (docking_score), 
                                    "le" (ligand efficiency), 
                                    "delta" (delta energy from best pose), 
                                    "ref_rmsd" (RMSD to reference pose), 
                                    "e_inter" (intermolecular energy), 
                                    "e_vdw" (van der waals energy), 
                                    "e_elec" (electrostatic energy), 
                                    "e_intra" (intermolecular energy), 
                                    "n_interact" (number of iteractions), 
                                    "ligand_smile" , 
                                    "rank" (rank of ligand pose), 
                                    "run" (run number for ligand pose), 
                                    "hb" (hydrogen bonds), 
                                    "receptor" (receptor name)
                bookmark_name (str): name for resulting book mark file. Default value is 'passing_results'
                options_dict (dict): write options as a dict
        Returns:
            int: number of ligands passing filter
        """
        
        self.set_filters(eworst=eworst, 
                        ebest=ebest, 
                        leworst=leworst, 
                        lebest=lebest, 
                        score_percentile=score_percentile,
                        le_percentile=le_percentile,
                        vdw_interactions=vdw_interactions,
                        hb_interactions=hb_interactions,
                        reactive_interactions=reactive_interactions,
                        interactions_count=interactions_count,
                        react_any=react_any,
                        max_miss=max_miss,
                        ligand_name=ligand_name,
                        ligand_substruct=ligand_substruct,
                        ligand_substruct_pos=ligand_substruct_pos,
                        ligand_max_atoms=ligand_max_atoms,
                        ligand_operator=ligand_operator,     
                        dict = filters_dict)

        if options_dict is not None:
            storage_dict, output_dict = RingtailCore.split_dict(options_dict, ["log_file", "enumerate_interaction_combs"])
        else:
            storage_dict = None
            output_dict = None
        self.set_storageman_attributes(output_all_poses = output_all_poses, 
                                        mfpt_cluster = mfpt_cluster, 
                                        interaction_cluster = interaction_cluster, 
                                        overwrite = overwrite, 
                                        order_results = order_results, 
                                        outfields = outfields, 
                                        bookmark_name = bookmark_name,
                                        dict=storage_dict)
        self.set_output_options(log_file=log_file,enumerate_interaction_combs=enumerate_interaction_combs, dict = output_dict)

        # Compatibility check with docking mode
        if self.docking_mode == "vina" and self.filters.react_any:
            logger.warning("Cannot use reaction filters with Vina mode. Removing react_any filter.")
            self.filters.react_any = False

        # make sure enumerate_interaction_combs always true if max_miss = 0, since we don't ever worry about the union in this case
        if self.filters.max_miss == 0:
            self.outputopts.enumerate_interaction_combs = True

        # guard against unsing percentile filter with all_poses
        if self.storageopts.output_all_poses and not (self.filters.score_percentile is None or self.filters.le_percentile is None):
            logger.warning(
                "Cannot return all passing poses with percentile filter. Will only log best pose."
            )
            self.storageopts.output_all_poses = False
            
        logger.info("Filtering results...")

        # get possible permutations of interaction with max_miss excluded
        interaction_combs = self._generate_interaction_combinations(
                            self.filters.max_miss)
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
                    filters_dict, not self.outputopts.enumerate_interaction_combs
                )
                if filtered_results is not None:
                    if return_iter:
                        return filtered_results
                    result_bookmark_name = self.storageman.get_current_view_name()
                    with self.outputman: 
                        self.outputman.write_filters_to_log(self.filters.todict(), combination, f"Morgan Fingerprints butina clustering cutoff: {self.storageman.mfpt_cluster}\nInteraction Fingerprints clustering cutoff: {self.storageman.interaction_cluster}")
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
                    self.outputman.write_results_bookmark_to_log(self.storageman.bookmark_name + "_union")
                    number_passing_union = self.outputman.write_log(maxmiss_union_results)
                    self.outputman.log_num_passing_ligands(number_passing_union)
                    print("\nNumber passing ligands in max_miss union:", number_passing_union)
                    ligands_passed = number_passing_union
            
        return ligands_passed
   
    def write_molecule_sdfs(self, sdf_path = None, bookmark_name = None,  write_nonpassing = None):
        """
        Have output manager write molecule sdf files for passing results in given results bookmark

        Args:
            sdf_path (str, optional): Optional path existing or to be created in cd where SDF files will be saved
            bookmark_name (str, optional): Option to run over specified bookmark other than that just used for filtering
            write_nonpassing (bool, optional): Option to include non-passing poses for passing ligands
        """

        if sdf_path is not None:
            self.set_output_options(export_sdf_path=sdf_path)
        
        all_mols = self.ligands_rdkit_mol(bookmark_name=bookmark_name, write_nonpassing=write_nonpassing)
        
        for ligname, info in all_mols.items():
            logger.info("Writing " + ligname + ".sdf")
            self.outputman.write_out_mol(
                ligname, info["ligand"], info["flex_residues"], info["properties"]
            )
        
    def ligands_rdkit_mol(self, bookmark_name = None, write_nonpassing=False) -> dict:
        """
        Creates a dictionary of RDKit mols of all ligands specified from a bookmark, either excluding (default) or including 
        those ligands that did not pass the filter(s).
        
        Args:
            bookmark_name (str, optional): Option to run over specified bookmark other than that just used for filtering
            write_nonpassing (bool, optional): Option to include non-passing poses for passing ligands
        Returns:
            all_mols (dict): containing ligand names, RDKit mols, flexible residue bols, and other ligand properties
        """

        if bookmark_name is not None:
            self.set_storageman_attributes(bookmark_name=bookmark_name)

        with self.storageman:        
            if self.filters.max_miss > 0:
                logger.warning("WARNING: Requested 'export_sdf_path' with 'max_miss'. Exported SDFs will be for union of interaction combinations.")
                self.storageman.bookmark_name = self.storageman.bookmark_name + "_union"
            if not self.storageman.check_passing_view_exists():
                logger.warning(
                    "Given results bookmark does not exist in database. Cannot write passing molecule SDFs"
                )
                return None
            
            # make temp table 
            self.storageman.create_temp_passing_table()
            passing_molecule_info = self.storageman.fetch_passing_ligand_output_info()
            flexible_residues, flexres_atomnames = self.storageman.fetch_flexres_info()

            if flexible_residues != []:
                flexible_residues = json.loads(flexible_residues)
                flexres_atomnames = json.loads(flexres_atomnames)

            all_mols = {}
            for (ligname, smiles, atom_indices, h_parent_line) in passing_molecule_info:
                logger.info("Creating an RDKIT mol for ligand: " + ligname + ".")
                # create rdkit ligand molecule and flexible residue container
                if smiles == "":
                    logger.warning(
                        f"No SMILES found for {ligname}. Cannot create SDF."
                    )
                    continue
                # some work needed bc of info needed for creating rd kit, can I remove more than flex stuff? 
                mol, flexres_mols, properties = self._create_rdkit_mol(ligname, smiles, atom_indices, h_parent_line, flexible_residues, flexres_atomnames, write_nonpassing=write_nonpassing)
                
                all_mols[ligname] = {"ligand": mol, "flex_residues": flexres_mols, "properties": properties}
        
        return all_mols

    def find_similar_ligands(self, query_ligname: str):
        """
        Find ligands in cluster with query_ligname
        Args:
            query_ligname (str): name of the ligand in the ligand table to look for similars to
        Returns:
            int: number of ligands that are similar
        """

        number_similar = 0
        with self.storageman: similar_ligands, bookmark_name, cluster_name = self.storageman.fetch_clustered_similars(query_ligname)

        if similar_ligands is not None:
            with self.outputman:
                self.outputman.write_find_similar_header(query_ligname, cluster_name)
                self.outputman.write_results_bookmark_to_log(bookmark_name)
                number_similar = self.outputman.write_log(similar_ligands)
                self.outputman.log_num_passing_ligands(number_similar)
                print("Number similar ligands:", number_similar)
        return number_similar
        
    def plot(self, save=True):
        """
        Get data needed for creating Ligand Efficiency vs
        Energy scatter plot from storageManager. Call OutputManager to create plot.
        """
        if self.filters.max_miss > 0:
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

    def display_pymol(self, bookmark_name = None):
        """
        Launch pymol session and plot of LE vs docking score. Displays molecules when clicked.
        Args:
            bookmark_name (str): bookmark name to use in pymol. 'None' uses the whole db? 
        """

        import subprocess
        from rdkit.Chem import PyMol

        # launch pymol session
        p = subprocess.Popen(
        ["pymol", "-R"],
        stdout=subprocess.PIPE,
        )
        # ensure pymol was opened
        import time
        time.sleep(2)
        
        if bookmark_name is not None:
            self.set_storageman_attributes(bookmark_name=bookmark_name)
            
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

                mol, flexres_mols, _ = self._create_rdkit_mol(ligname, ligand_smile, atom_index_map, hydrogen_parents, flexible_residues, flexres_atomnames, pose_ID=chosen_pose[0])
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
        
    def export_bookmark_db(self, bookmark_name: str = None):
        """Export database containing data from bookmark

        Args:
            bookmark_db_name (str): name for bookmark_db
        """
        if bookmark_name is not None:
            self.set_storageman_attributes(bookmark_name=bookmark_name)
        bookmark_db_name = (
                        self.db_file.rstrip(".db")
                        + "_"
                        + self.storageman.bookmark_name
                        + ".db"
                    )
        logger.info("Exporting bookmark database")
        if os.path.exists(bookmark_db_name):
            logger.warning(
                "Requested export DB name already exists. Please rename or remove existing database. New database not exported."
            )
            return
        with self.storageman: self.storageman.clone(bookmark_db_name)
        # connect to cloned database
        self.db_file = bookmark_db_name
        #TODO needs rejiggering so agnostic to db engine. It also rewires core to new db, is that ok? 
        dictionary = self.storageopts.todict()
        dictionary["db_file"] = self.db_file
        with StorageManagerSQLite(**dictionary) as db_clone:
            db_clone.prune()
            db_clone.close_storage(vacuum=True)

    def export_receptors(self):
        """
        Export receptor in database to pdbqt
        """
        with self.storageman: receptor_tuples = self.storageman.fetch_receptor_objects()
        for recname, recblob in receptor_tuples:
            if recblob is None:
                logger.warning(f"No receptor pdbqt stored for {recname}. Export failed.")
                continue
            if not hasattr(self, "outputman"):
                self.set_output_options()
            self.outputman.write_receptor_pdbqt(recname, recblob)

    def get_previous_filter_data(self, outfields = None, bookmark_name = None):
        """Get data requested in self.out_opts['outfields'] from the
        results view of a previous filtering
        Args:
            outfields (str): use outfields as described in RingtailOptions > StorageOptions
            bookmark_name (str): bookmark for which the filters were used
        """
        if bookmark_name is not None:
            self.set_storageman_attributes(bookmark_name=bookmark_name)
        if outfields is not None:
            self.set_storageman_attributes(outfields=outfields)

        with self.storageman: new_data = self.storageman.fetch_data_for_passing_results()
        with self.outputman: self.outputman.write_log(new_data)

    def drop_bookmark(self, bookmark_name: str):
        """Drops specified bookmark from the database
        Args:
            bookmark_name (str): name of bookmark to be dropped."""
        
        with self.storageman: self.storageman._drop_bookmark(bookmark_name=bookmark_name)
        logger.info("Bookmark {0} was dropped from the database {1}".format(bookmark_name, self.storageman.db_file))

    def add_config_from_file(self, config_file: str ="config.json"):
        """
        Provide ringtail config from file, will directly set storage manager settings, 
        and return dictionaries for the remaining options. 
        Args:
            config_file: json formatted file containing ringtail and filter options
        Returns:
            file_dict: dictionary of files to use in add results
            write_dict: dictionary of options to use in add results
            output_dict: dictionary of options to use for filtering
            filters_dict: dictionary of files containing filters
        """
        (file_dict, write_dict, output_dict, filters_dict, storageman_dict) = RingtailCore.read_config_file(config_file=config_file)

        self.set_storageman_attributes(dict= storageman_dict)
        logger.info("A dictionary containing storage options was extracted from config file and assigned storagemanager.")

        return (file_dict, write_dict, output_dict, filters_dict)

    @staticmethod
    def generate_config_json_template(to_file: bool = True) -> str:
        """
        Creates a dict of all Ringtail option classes, and their 
        key-default value pairs. Outputs to options.json in 
        "util_files" of to_file = true, else it returns the dict 
        of default option values.
        Args:
            to_file (bool): if true writes file to standard options json path, if false returns a json string with values
        Return:
            str: filename or json string with options
        """
        
        json_string = {"outputopts":OutputOptions().todict(),
                       "resultsmanopts":ResultsProcessingOptions().todict(),
                       "storageopts":StorageOptions().todict(),
                       "filters":Filters().todict(),
                       "fileobj": InputFiles().todict()}
        if to_file:
            filename="config.json"
            with open(filename, 'w') as f: 
                f.write(json.dumps(json_string, indent=4))
            logger.debug(f"Default ringtail option values written to file {filename}.")
            return filename
        else:
            logger.debug("Default ringtail option values prepared as a string.")
            return json_string

    @staticmethod
    def read_config_file(config_file: str = "config.json", return_as_string = False):
        """
        Will read and parse a file containing ringtail options following the format 
        given from RingtailCore.generate_config_json_template
        Args:
            config_file (str, file name in cd): json formatted file containing ringtail and filter options
            return_as_string (bool): will return all dictionaries as dict without sections if true
        Returns:
            file_dict: dictionary of files to use in add results
            write_dict: dictionary of options to use in add results
            output_dict: dictionary of options to use for filtering
            filters_dict: dictionary of files containing filters
            storageman_dict: dictionary of storage manager options
        """

        try: 
            if config_file is None: 
                raise OptionError("No config file was found in the Ringtail/util_files directory.")
            
            filepath = "config.json"
            with open(filepath, "r") as f:
                logger.info("Reading Ringtail options from config file")
                options: dict = json.load(f)

            file_dict = options["fileobj"]
            logger.info("A dictionary containing results files was extracted from config file.")
            write_dict = {**options["resultsmanopts"], 
                        "append_results": options["storageopts"]["append_results"], 
                        "duplicate_handling":options["storageopts"]["duplicate_handling"], 
                        "overwrite":options["storageopts"]["overwrite"]}
            logger.info("A dictionary containing write options was extracted from config file.")

            storageman_dict = options["storageopts"]
            logger.info("A dictionary containing storage options was extracted from config file.")

            output_dict = options["outputopts"]
            logger.info("A dictionary containing database read options was extracted from config file.")

            filters_dict = options["filters"]
            logger.info("A dictionary containing filters was extracted from config file.")

            if return_as_string:
                return file_dict | write_dict | output_dict | filters_dict | storageman_dict
            else:
                return (file_dict, write_dict, output_dict, filters_dict, storageman_dict)
            
        except FileNotFoundError:
            logger.error("Please ensure config file is in the working directory.")
        except Exception as e:
            OptionError(f"There were issues with the configuration file: {e}")

    @staticmethod       
    def get_defaults(object: str = "all") -> dict:
        """
        Gets default values from RingtailOptions and returns dict of all options,
        or options belonging to a specific group. 
        
        Args:
            object (str): ["all", "resultsmanopts", "storageopts", "outputopts", "filters", "fileobj"]
        """

        all_defaults = RingtailCore.generate_config_json_template(to_file=False)

        if object.lower() not in ["all", "resultsmanopts", "storageopts", "outputopts", "filterobj", "fileobj"]:
            raise OptionError(f'The options object {object.lower()} does not exist. Please choose amongst \n ["all", "writeopts", "storageopts", "outputopts", "filterobj", "fileobj"]')
        
        if object.lower() == "all":
            logger.debug("All ringtail default values have been fetched.")
            return all_defaults
        else:
            logger.debug(f"Ringtail default values for {object} have been fetched.")
            return all_defaults[object.lower()]
        
    @staticmethod       
    def get_options_info(object: str = "all") -> dict:
        """
        Gets default values from RingtailOptions and returns dict of all,
        or specific object. 
        
        Args:
            object (str): ["all", "writeopts", "storageopts", "outputopts", "filters", "fileobj"]
        """
        all_info = {"outputopts": OutputOptions.options,
                    "fileobj": InputFiles.options,
                    "writeopts": ResultsProcessingOptions.options,
                    "storageopts": StorageOptions.options,
                    "filters":Filters.options,}
        
        if object.lower() not in ["all", "writeopts", "storageopts", "outputopts", "filters", "fileobj"]:
            raise OptionError(f'The options object {object.lower()} does not exist. Please choose amongst \n ["all", "writeopts", "storageopts", "outputopts", "filters", "fileobj"]')
        if object.lower() == "all":
            all_info_one_dict = {}
            for _,v in all_info.items():
                all_info_one_dict.update(v)
            logger.debug("All ringtail default values have been fetched.")
            return all_info_one_dict
        else:
            logger.debug(f"Ringtail default values for {object} have been fetched.")
            return all_info[object.lower()]
    
#-#-#- Util method -#-#-# 
    @staticmethod
    def split_dict(dict: dict, items: list) -> tuple:
        """ Utility method that takes one dictionary and splits it into two based on the listed keys 
            Args: 
                dict (dict): original dictionary
                items (list): list of keys to use for separation
            Returns:
                tuple of:
                    dict (dict): original dict minus the removed items
                    new_dict (dict): dict containing the items removed from the original dict
        """
        new_dict = {}

        for key in items:
            new_dict[key] = dict.pop(key)
        
        return dict, new_dict
