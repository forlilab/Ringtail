#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail virtual screening manager
#

import matplotlib.pyplot as plt
import json
from meeko import RDKitMolCreate
from meeko.export_flexres import pdb_updated_flexres_from_rdkit
from .storagemanager import StorageManager
from .resultsmanager import ResultsManager
from .receptormanager import ReceptorManager
from .outputmanager import OutputManager
from .ringtailoptions import *
from .util import *
from .exceptions import RTCoreError, OutputError, StorageError
from rdkit import Chem
import itertools
import os
from .logutils import LOGGER


class RingtailCore:
    """Core class for coordinating different actions on virtual screening
    including adding results to storage, filtering and clusteirng, and outputting data as
    rdkit molecules, plotting docking results, and visualizing select ligands in pymol.

    Attributes:
        db_file (str): name of database file being operated on
        docking_mode (str): specifies what docking mode has been used for the results in the database
        storageman (StorageManager): Interface module with database
        resultsman (ResultsManager): Module to deal with results processing before adding to database
        outputman (OutputManager): Manager for output tasks of log-writting, plotting, ligand SDF writing, starting pymol sessions
        filters (Filters): object holding all optional filters
        _run_mode (str): refers to whether ringtail is ran from the command line or through direct API use, where the former is more restrictive
    """

    # region #-#-#- Base methods -#-#-#

    def __init__(
        self,
        db_file: str = "output.db",
        storage_type: str = "sqlite",
        docking_mode: str = "dlg",
        logging_level: str = "WARNING",
    ):
        """Initialize ringtail core, and create a storageman object with the db file.
        Can set logger level here, otherwise change it by logger.setLevel("level")

        Args:
            db_file (str): Database file to initialize core with. Defaults to "output.db".
            storage_type (str, optional): Database setup to use for interacting with the database. Defaults to "sqlite".
            docking_mode (str, optional): Docking mode for which the database will be used
            logging_level (str, optional): Global logger level. Defaults to "DEBUG".
        """

        # Initiate logging
        self.logger = LOGGER
        self.logger.set_level(logging_level)
        # create log file handler if log level is debug
        if self.logger.level() == "DEBUG":
            self.logger.add_filehandler()

        self.logger.info(
            f"[     New RingtailCore object initialized with database file {db_file}    ]"
        )
        # Check if storage type is implemented
        try:
            storageman = StorageManager.check_storage_compatibility(storage_type)
        except NotImplementedError as e:
            self.logger.error(e)
            raise e

        # initialize remaining variables and storagemanager
        self.storagetype = storage_type
        self.db_file = db_file
        self.storageman = storageman(db_file)
        self._run_mode = "api"
        self._docking_mode = docking_mode
        self.set_storageman_attributes()

    def update_database_version(self, consent=False, new_version="2.0.0"):
        """Method to update database version from earlier versions to either 1.1.0 or 2.0.0"""
        with self.storageman:
            self.storageman.update_database_version(new_version, consent)
        # return self.storageman.update_database_version(new_version, consent)

    # -#-#- Private methods -#-#-#

    def _validate_docking_mode(self, docking_mode: str):
        """Method that validates specified AutoDock program used to generate results.

        Args:
            docking_mode (str): string that describes docking mode

        Raises:
            RTCoreError: if docking_mode is not supported
        """
        if type(docking_mode) is not str:
            self.logger.warning(
                'The given docking mode was not given as a string, it will be set to default value "dlg".'
            )
            self._docking_mode = "dlg"
        elif docking_mode.lower() not in ["dlg", "vina"]:
            raise NotImplementedError(
                f'Docking mode {docking_mode} is not supported. Please choose between "dlg" and "vina".'
            )
        else:
            self._docking_mode = docking_mode.lower()
            self.logger.debug(f"Docking mode set to {self.docking_mode}.")

    def _get_docking_mode(self):
        """
        Private method to retrieve docking mode

        Returns:
            str: docking mode
        """
        return self._docking_mode

    # making docking mode a property of the ringtail core object
    docking_mode = property(fget=_get_docking_mode, fset=_validate_docking_mode)

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
            atom_indices (list): List of ints indicating mapping of coordinate indices to smiles indices
            poses (iterable): iterable containing ligand_pose, flexres_pose, flexres_names
            mol (RDKit.Chem.Mol): RDKit molecule for ligand
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
            interactions = self.storageman.fetch_pose_interactions(Pose_ID)
            # if that pose id has interactions
            if interactions is not None:
                # make a list of all of them
                interactions_list = []
                # for each interaction row, make into a string according to format above
                for interaction_info in interactions:
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
            self.logger.warning(
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
            dict: filters for storageman
        """

        filters_dict = self.filters.todict()
        for itype in Filters.get_filter_keys("interaction"):
            itype_interactions = filters_dict[itype]
            for interaction in itype_interactions:
                if itype + "-" + interaction[0] not in interaction_combination:
                    filters_dict[itype].remove(interaction)
        # we want a match for all interactions when enumerating combinations
        filters_dict["max_miss"] = 0

        return filters_dict

    def _create_rdkit_mol(
        self,
        ligname,
        smiles,
        atom_indices,
        h_parent_line,
        flexible_residues,
        flexres_atomnames,
        pose_ID=None,
        write_nonpassing=False,
    ):
        """Creates rdkit molecule for given ligand, either for a specific pose_ID or for all passing (and nonpassing?) poses

        Args:
            ligname (string): ligand name
            smiles (string): ligand smiles string
            atom_indices (list): list of atom indices converting pdbqt to rdkit mol
            h_parent_line (list): list of atom indices for heteroatoms with attached hydrogens
            flexible_residues (list): list of flexible residue names
            flexres_atomnames (list): list of atomtypes in flexible residue
            pose_ID (int, optional): pose_ID for single pose to return. Defaults to None.
            write_nonpassing (bool, optional): whether or not to write ligands not passing filter

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
            passing_properties = self.storageman.fetch_passing_pose_properties(ligname)
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
            pose_properties = self.storageman.fetch_single_pose_properties(pose_ID)
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
        # add ligand name to properties
        properties["_Name"] = ligname
        # add hydrogens to mols
        lig_h_parents = [int(idx) for idx in json.loads(h_parent_line)]
        mol = RDKitMolCreate.add_hydrogens(mol, ligand_saved_coords, lig_h_parents)
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
        dict: dict = None,
    ) -> InputFiles:
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
            InputFiles
        """

        # Method to deal with the file read manager currently expecting file data to be presented as double [[]] lists
        def ensure_double_list(object) -> list:
            """Most of ringtail is set up to handle files, file paths, and file lists as double lists [[items]].
            Instead of changing that for now, the input is checked to ensure they are presented as, or can be
            converted to a double list. Really only matters for using API."""
            if type(object) == list:
                if type(object[0]) == list:
                    if type(object[0][0]) == str:
                        pass
                    else:
                        raise OptionError(
                            f"error, object is more than two encapsulated lists: '{object[0][0]}' should be a string."
                        )
                elif type(object[0]) == str:
                    object = [object]
                else:
                    self.logger.error("Unable to parse file input.")
            elif type(object) == str:
                object = [[object]]
            else:
                self.logger.error("Unable to parse file input.")

            return object

        # keywords this pertains to
        need_to_be_double_list = ["file", "file_path", "file_list"]

        # Set file format automatically if not specified
        if file_pattern is None:
            if self.docking_mode == "dlg":
                file_pattern = "*.dlg*"

            elif self.docking_mode == "vina":
                file_pattern = "*.pdbqt*"
            else:
                self.logger.error(
                    "Docking mode and file pattern was not specified, can not continue."
                )

        # Dict of individual arguments
        individual_options = {
            "file": file,
            "file_path": file_path,
            "file_list": file_list,
            "file_pattern": file_pattern,
            "recursive": recursive,
            "receptor_file": receptor_file,
            "save_receptor": save_receptor,
        }

        # Create option object with default values if needed
        files = InputFiles()

        # Set options from dict if provided
        if dict is not None:
            for k, v in dict.items():
                if k in need_to_be_double_list and v is not None:
                    v = ensure_double_list(v)
                setattr(files, k, v)
            self.logger.debug(
                f"File attributes {dict} were assigned from provided option dictionary."
            )

        # Set additional options from individual arguments
        # NOTE Will overwrite provided dictionary
        for k, v in individual_options.items():
            if v is not None:
                if k in need_to_be_double_list:
                    v = ensure_double_list(v)
                setattr(files, k, v)
        self.logger.debug(
            f"File attributes {individual_options} were assigned from provided individual options."
        )

        return files

    def _set_results_sources(
        self,
        results_strings: dict = None,
        receptor_file: str = None,
        save_receptor: bool = None,
        dict: dict = None,
    ) -> InputStrings:
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
        individual_options = {
            "results_strings": results_strings,
            "receptor_file": receptor_file,
            "save_receptor": save_receptor,
        }
        # Create option object with default values if needed
        strings = InputStrings()

        # Set options from dict if provided
        if dict is not None:
            for k, v in dict.items():
                setattr(strings, k, v)
            self.logger.debug(f"Docking string results attributes {dict} were set.")

        # Set additional options from individual arguments
        # NOTE Will overwrite config file
        for k, v in individual_options.items():
            if v is not None:
                setattr(strings, k, v)
        self.logger.debug(
            f"Docking string results attributes {individual_options} were set."
        )

        return strings

    def _create_resultsmanager(
        self, file_sources: InputFiles = None, string_sources: InputStrings = None
    ) -> ResultsManager:
        """Creates a results manager object based on results provided either as files or as strings (currently only for vina).
        Will create a new object each time results are added, and use the docking mode specified as an attribute of the core.
        In its current state it assumes only one source of results will be provided.
        If both are provided, it will only process the strings.

        Args:
            file_sources (InputFiles): used if docking results are provided through files
            string_sources (InputStrings): used if docking results are provided through strings

        Raises:
            RTCoreError
        """
        if file_sources is not None:
            self.resultsman = ResultsManager(
                file_sources=file_sources, docking_mode=self.docking_mode
            )
            self.logger.debug(
                "Results manager object has been created with results files."
            )
        elif string_sources is not None:
            self.resultsman = ResultsManager(
                string_sources=string_sources, docking_mode=self.docking_mode
            )
            self.logger.debug(
                "Results manager object has been created with results strings."
            )
        else:
            raise RTCoreError(
                "No results sources were provided, a results manager object could not be created."
            )

    def _add_results(
        self,
        results_sources,
        strings=False,
        duplicate_handling: str = None,
        overwrite: bool = None,
        store_all_poses: bool = None,
        max_poses: int = None,
        add_interactions: bool = None,
        interaction_tolerance: float = None,
        interaction_cutoffs: list = None,
        max_proc: int = None,
        options_dict: dict | None = None,
        finalize: bool = True,
    ):
        """Method that is agnostic of results type, and will do the actual call to storage manager to process result files and add to database.

        Args:
            results_sources(InputFiles or InputStrings): type checked and validated results object
            strings (bool): whether or not results are provided as strings or files
            duplicate_handling (str): specify how duplicate Results rows should be handled when inserting into database. Options are "ignore" or "replace". Default behavior will allow duplicate entries.
            store_all_poses (bool): store all ligand poses, does it take precedence over max poses?
            max_poses (int): how many poses to save (ordered by soem score?)
            add_interactions (bool): add ligand-receptor interaction data, only in vina mode
            interaction_tolerance (float): longest ångström distance that is considered interaction?
            interaction_cutoffs (list): ångström distance cutoffs for x and y interaction
            max_proc (int): max number of computer processors to use for file reading
            options_dict (dict): write options as a dict

        Raises:
            OptionError
        """
        # if dictionary of options provided, attribute to appropriate managers
        if options_dict is not None:
            write_dict, storage_dict = split_dict(
                options_dict, ["duplicate_handling", "overwrite"]
            )
        else:
            storage_dict = None
            write_dict = None
        self.set_storageman_attributes(
            duplicate_handling=duplicate_handling,
            overwrite=overwrite,
            dict=storage_dict,
        )

        if results_sources.save_receptor:
            self.save_receptor(results_sources.receptor_file)

        with self.storageman:
            # Prepare the results manager with the provided docking results sources
            if strings == False:
                self._create_resultsmanager(file_sources=results_sources)
            elif strings == True:
                self._create_resultsmanager(string_sources=results_sources)
            self.resultsman.storageman = self.storageman
            self.resultsman.storageman_class = self.storageman.__class__
            self.set_resultsman_attributes(
                store_all_poses,
                max_poses,
                add_interactions,
                interaction_tolerance,
                interaction_cutoffs,
                max_proc,
                write_dict,
            )

            # Docking mode compatibility check
            if (
                self.docking_mode == "vina"
                and self.resultsman.interaction_tolerance is not None
            ):
                self.logger.warning(
                    "Cannot use interaction_tolerance with Vina mode. Removing interaction_tolerance."
                )
                self.resultsman.interaction_tolerance = None

            # Process results files and handle database versioning
            self.storageman.check_storage_ready(
                self._run_mode,
                self.docking_mode,
                self.resultsman.store_all_poses,
                self.resultsman.max_poses,
            )
            self.logger.info("Adding results...")
            self.resultsman.process_docking_data()
            if finalize:
                self.storageman.finalize_database_write()

    # endregion

    # region #-#-#- Core attribute setting methods -#-#-#
    """ These methods are used internally to assing values to all ringtail options. 
    This ensures:   - that all options are set to specific types through RingtailOptions
                    - that internal consistency checks are performed on a group of options
                    - these methods ensure options are assigned to the appropriate ringtail manager classes"""

    def set_storageman_attributes(
        self,
        filter_bookmark: str = None,
        duplicate_handling: str = None,
        overwrite: bool = None,
        order_results: str = None,
        outfields: str = None,
        output_all_poses: str = None,
        mfpt_cluster: float = None,
        interaction_cluster: float = None,
        bookmark_name: str = None,
        dict: dict = None,
    ):
        """
        Create storage_manager_options object if needed, sets options, and assigns them to the storage manager object.

        Args:
            filter_bookmark (str): Perform filtering over specified bookmark. (in output group in CLI)
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
        individual_options = {
            "filter_bookmark": filter_bookmark,
            "duplicate_handling": duplicate_handling,
            "overwrite": overwrite,
            "order_results": order_results,
            "outfields": outfields,
            "output_all_poses": output_all_poses,
            "mfpt_cluster": mfpt_cluster,
            "interaction_cluster": interaction_cluster,
            "bookmark_name": bookmark_name,
        }

        # Create option object with default values if needed
        if not hasattr(self, "storageopts"):
            self.storageopts = StorageOptions()

        # Set options from dict if provided
        if dict is not None:
            for k, v in dict.items():
                if v is not None:
                    setattr(self.storageopts, k, v)
            self.logger.debug(f"Storage manager attributes {dict} were set.")

        # Set additional options from individual arguments
        # NOTE Will overwrite config file
        for k, v in individual_options.items():
            if v is not None:
                setattr(self.storageopts, k, v)
        self.logger.debug(f"Storage manager attributes {individual_options} were.")

        # Assign attributes to storage manager
        for k, v in self.storageopts.todict().items():
            setattr(self.storageman, k, v)
        self.logger.info("Options for storage manager have been changed.")

    def set_resultsman_attributes(
        self,
        store_all_poses: bool = None,
        max_poses: int = None,
        add_interactions: bool = None,
        interaction_tolerance: float = None,
        interaction_cutoffs: list = None,
        max_proc: int = None,
        dict: dict = None,
    ):
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
        individual_options = {
            "store_all_poses": store_all_poses,
            "max_poses": max_poses,
            "add_interactions": add_interactions,
            "interaction_tolerance": interaction_tolerance,
            "interaction_cutoffs": interaction_cutoffs,
            "max_proc": max_proc,
        }

        # Create option object with default values if needed
        if not hasattr(self, "resultsmanopts"):
            self.resultsmanopts = ResultsProcessingOptions()

        # Set options from dict if provided
        if dict is not None:
            for k, v in dict.items():
                if v is not None:
                    setattr(self.resultsmanopts, k, v)
            self.logger.debug(f"Results manager attributes {dict} were set.")

        # Set additional options from individual arguments
        # NOTE Will overwrite config file
        for k, v in individual_options.items():
            if v is not None:
                setattr(self.resultsmanopts, k, v)
        self.logger.debug(f"Results manager attributes {individual_options} were set.")

        # Assigns options to the results manager object
        for k, v in self.resultsmanopts.todict().items():
            if v is not None:
                setattr(self.resultsman, k, v)
        self.logger.info("Options for results manager have been changed.")

    def set_output_options(
        self,
        log_file: str = None,
        export_sdf_path: str = None,
        enumerate_interaction_combs: bool = None,
        dict: dict = None,
    ):
        """Creates output options object that holds attributes related to reading and outputting results.
        Will assign log_file name and export_sdf_path to the output_manager object.

        Args:
            log_file (str): by default, results are saved in "output_log.txt"; if this option is used, ligands and requested info passing the filters will be written to specified file
            export_sdf_path (str): specify the path where to save poses of ligands passing the filters (SDF format); if the directory does not exist, it will be created; if it already exist, it will throw an error, unless the --overwrite is used  NOTE: the log file will be automatically saved in this path. Ligands will be stored as SDF files in the order specified.
            enumerate_interaction_combs (bool): When used with `max_miss` > 0, will log ligands/poses passing each separate interaction filter combination as well as union of combinations. Can significantly increase runtime.
            dict (dict): dictionary of one or more of the above args, is overwritten by individual args

        """
        # Dict of individual arguments
        individual_options = {
            "log_file": log_file,
            "export_sdf_path": export_sdf_path,
            "enumerate_interaction_combs": enumerate_interaction_combs,
        }
        # Create option object with default values if needed
        if not hasattr(self, "outputopts"):
            self.outputopts = OutputOptions()

        # Set options from dict if provided
        if dict is not None:
            for k, v in dict.items():
                if v is not None:
                    setattr(self.outputopts, k, v)
            self.logger.debug(f"Output options {dict} were set.")

        # Set additional options from individual arguments
        # NOTE Will overwrite config file
        for k, v in individual_options.items():
            if v is not None:
                setattr(self.outputopts, k, v)
        self.logger.debug(f"Output options {individual_options} were set.")

        # Creates output man with attributes if needed
        self.outputman = OutputManager(
            self.outputopts.log_file, self.outputopts.export_sdf_path
        )
        self.logger.info("Options for output manager have been changed.")

    def set_filters(
        self,
        eworst=None,
        ebest=None,
        leworst=None,
        lebest=None,
        score_percentile=None,
        le_percentile=None,
        vdw_interactions=None,
        hb_interactions=None,
        reactive_interactions=None,
        hb_count=None,
        react_any=None,
        max_miss=None,
        ligand_name=None,
        ligand_operator=None,
        ligand_substruct=None,
        ligand_substruct_pos=None,
        ligand_max_atoms=None,
        dict: dict = None,
    ):
        """
        Create a filter object containing all numerical and string filters.

        Args:
            eworst (float): specify the worst energy value accepted
            ebest (float): specify the best energy value accepted
            leworst (float): specify the worst ligand efficiency value accepted
            lebest (float): specify the best ligand efficiency value accepted
            score_percentile (float): specify the worst energy percentile accepted. Express as percentage e.g. 1 for top 1 percent.
            le_percentile (float): specify the worst ligand efficiency percentile accepted. Express as percentage e.g. 1 for top 1 percent.
            vdw_interactions (list[tuple]): define van der Waals interactions with residue as [-][CHAIN]:[RES]:[NUM]:[ATOM_NAME]. E.g., [('A:VAL:279:', True), ('A:LYS:162:', True)] -> [('chain:resname:resid:atomname', <wanted (bool)>), ('chain:resname:resid:atomname', <wanted (bool)>)]
            hb_interactions (list[tuple]): define HB (ligand acceptor or donor) interaction as [-][CHAIN]:[RES]:[NUM]:[ATOM_NAME]. E.g., [('A:VAL:279:', True), ('A:LYS:162:', True)] -> [('chain:resname:resid:atomname', <wanted (bool)>), ('chain:resname:resid:atomname', <wanted (bool)>)]
            reactive_interactions (list[tuple]): check if ligand reacted with specified residue as [-][CHAIN]:[RES]:[NUM]:[ATOM_NAME]. E.g., [('A:VAL:279:', True), ('A:LYS:162:', True)] -> [('chain:resname:resid:atomname', <wanted (bool)>), ('chain:resname:resid:atomname', <wanted (bool)>)]
            hb_count (list[tuple]): accept ligands with at least the requested number of HB interactions. If a negative number is provided, then accept ligands with no more than the requested number of interactions. E.g., [('hb_count', 5)]
            react_any (bool): check if ligand reacted with any residue
            max_miss (int): Will compute all possible combinations of interaction filters excluding up to max_miss numer of interactions from given set. Default will only return union of poses interaction filter combinations. Use with 'enumerate_interaction_combs' for enumeration of poses passing each individual combination of interaction filters.
            ligand_name (list[str]): specify ligand name(s). Will combine name filters with OR, e.g., ["lig1", "lig2"]
            ligand_substruct (list[str]): SMARTS, index of atom in SMARTS, cutoff dist, and target XYZ coords, e.g., ["ccc", "CN"]
            ligand_substruct_pos (list[str]): SMARTS pattern(s) for substructure matching, e.g., ['"[Oh]C" 0 1.2 -5.5 10.0 15.5'] -> ["smart_string index_of_positioned_atom cutoff_distance x y z"]
            ligand_max_atoms (int): Maximum number of heavy atoms a ligand may have
            ligand_operator (str): logical join operator for multiple SMARTS (default: OR), either AND or OR
            dict (dict): dictionary of one or more of the above args, is overwritten by individual args
        """

        # Dict of individual arguments
        individual_options = {
            "eworst": eworst,
            "ebest": ebest,
            "leworst": leworst,
            "lebest": lebest,
            "score_percentile": score_percentile,
            "le_percentile": le_percentile,
            "vdw_interactions": vdw_interactions,
            "hb_interactions": hb_interactions,
            "reactive_interactions": reactive_interactions,
            "hb_count": hb_count,
            "react_any": react_any,
            "max_miss": max_miss,
            "ligand_operator": ligand_operator,
            "ligand_name": ligand_name,
            "ligand_substruct": ligand_substruct,
            "ligand_substruct_pos": ligand_substruct_pos,
            "ligand_max_atoms": ligand_max_atoms,
        }

        # Create a filter object
        self.filters = Filters()

        # Set options from dict if provided
        if dict is not None:
            for k, v in dict.items():
                if v is not None:
                    setattr(self.filters, k, v)
            self.logger.debug(f"Filter {dict} were set.")

        # Set additional options from individual arguments
        # NOTE Will overwrite config file
        for k, v in individual_options.items():
            if v is not None:
                setattr(self.filters, k, v)
        self.logger.debug(f"Filter {individual_options} were set.")

        self.logger.info("A filter object has been prepared.")

    # endregion

    # region # -#-#- API -#-#-#
    def add_results_from_files(
        self,
        file: str = None,
        file_path: str = None,
        file_list: str = None,
        file_pattern: str = None,
        recursive: bool = None,
        receptor_file: str = None,
        save_receptor: bool = None,
        filesources_dict: dict = None,
        duplicate_handling: str = None,
        overwrite: bool = None,
        store_all_poses: bool = None,
        max_poses: int = None,
        add_interactions: bool = None,
        interaction_tolerance: float = None,
        interaction_cutoffs: list = None,
        max_proc: int = None,
        options_dict: dict = None,
        finalize: bool = True,
    ):
        """
        Call storage manager to process result files and add to database. Creates or adds to an existing a database.
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
            duplicate_handling (str, options): specify how duplicate Results rows should be handled when inserting into database. Options are "ignore" or "replace". Default behavior will allow duplicate entries.
            store_all_poses (bool): store all ligand poses, does it take precedence over max poses?
            max_poses (int): how many poses to save (ordered by soem score?)
            add_interactions (bool): add ligand-receptor interaction data, only in vina mode
            interaction_tolerance (float): longest ångström distance that is considered interaction?
            interaction_cutoffs (list): ångström distance cutoffs for x and y interaction
            max_proc (int): max number of computer processors to use for file reading
            options_dict (dict): write options as a dict

        Raises:
            OptionError

        """

        files = self._set_file_sources(
            file,
            file_path,
            file_list,
            file_pattern,
            recursive,
            receptor_file,
            save_receptor,
            filesources_dict,
        )
        results_files_given = (
            files.file is not None
            or files.file_path is not None
            or files.file_list is not None
        )

        if not results_files_given and not files.save_receptor:
            raise OptionError(
                "At least one input option needs to be used: file, file_path, file_list, or save_receptor"
            )

        # If there are ligand files present, process ligand data
        if results_files_given or files.save_receptor:
            self._add_results(
                files,
                False,
                duplicate_handling,
                overwrite,
                store_all_poses,
                max_poses,
                add_interactions,
                interaction_tolerance,
                interaction_cutoffs,
                max_proc,
                options_dict,
                finalize,
            )

    def add_results_from_vina_string(
        self,
        results_strings: dict = None,
        receptor_file: str = None,
        save_receptor: bool = None,
        resultsources_dict: dict = None,
        duplicate_handling: str = None,
        overwrite: bool = None,
        store_all_poses: bool = None,
        max_poses: int = None,
        add_interactions: bool = None,
        interaction_cutoffs: list = None,
        max_proc: int = None,
        options_dict: dict = None,
        finalize: bool = True,
    ):
        """
        Call storage manager to process the given vina output string and add to database.
        Options can be provided as a dict or as individual options.
        Creates or adds to an existing a database.

        Args:
            results_string (dict): string containing the ligand identified and docking results as a dictionary
            receptor_file (str): string containing the receptor .pdbqt
            save_receptor (bool): whether or not to store the full receptor details in the database (needed for some things)
            resultsources_dict (dict): file sources already as an object
            duplicate_handling (str, options): specify how duplicate Results rows should be handled when inserting into database. Options are "ignore" or "replace". Default behavior will allow duplicate entries.
            store_all_poses (bool): store all ligand poses, does it take precedence over max poses?
            max_poses (int): how many poses to save (ordered by soem score?)
            add_interactions (bool): add ligand-receptor interaction data, only in vina mode
            interaction_cutoffs (list): ångström distance cutoffs for x and y interaction
            max_proc (int): max number of computer processors to use for file reading
            options_dict (dict): write options as a dict

        Raises:
            OptionError

        """
        # Method currently only works with vina output, set automatically
        if self.docking_mode != "vina":
            self.docking_mode = "vina"

        # create results string object
        results = self._set_results_sources(
            results_strings, receptor_file, save_receptor, resultsources_dict
        )
        results_strings_given = bool(results.results_strings)

        if not results_strings_given and not results.save_receptor:
            raise OptionError(
                "At least one input option needs to be used: 'results_strings', or 'save_receptor'"
            )

        # If there are any docking data strings, process docking results
        if results_strings_given or results.save_receptor:
            self._add_results(
                results,
                True,
                duplicate_handling,
                overwrite,
                store_all_poses,
                max_poses,
                add_interactions,
                None,
                interaction_cutoffs,
                max_proc,
                options_dict,
                finalize,
            )

    def finalize_write(self):
        """
        Finalize database write by creating interaction tables and setting database version
        """
        with self.storageman:
            self.storageman.finalize_database_write()

    def save_receptor(self, receptor_file):
        """
        Add receptor to database.

        Args:
            receptor_file (str): path to receptor file

        """
        receptor_list = ReceptorManager.make_receptor_blobs([receptor_file])
        with self.storageman:
            for rec, rec_name in receptor_list:
                # NOTE: in current implementation, only one receptor allowed per database
                # Check that any receptor row is incomplete (needs receptor blob) before inserting
                filled_receptor_rows = self.storageman.count_receptors_in_db()
                if (
                    filled_receptor_rows != 0
                ):  # throw error if receptor is already present
                    raise RTCoreError(
                        "Expected Receptors table to have no receptor objects present, already has {0} receptor present. Cannot add more than 1 receptor to a database.".format(
                            filled_receptor_rows
                        )
                    )
                self.storageman.insert_receptor_blob(rec, rec_name)
                self.logger.info("Receptor data was added to the database.")

    def produce_summary(
        self, columns=["docking_score", "leff"], percentiles=[1, 10]
    ) -> None:
        """Print summary of data in storage to sdout

        Args:
            columns (list(str)): data columns used to prepare summary
            percentiles (list(int)): cutoff percentiles for the summary

        """
        with self.storageman:
            summary_data = self.storageman.fetch_summary_data(columns, percentiles)
        if "summary_data" not in locals():
            raise RTCoreError(f"Summary data is empty, please check the database.")
        print("Total Stored Ligands          :", summary_data.pop("num_ligands"))
        print("Total Stored Poses            :", summary_data.pop("num_poses"))
        print(
            "Total Unique Interactions     :",
            summary_data.pop("num_unique_interactions"),
        )
        print(
            "Number Interacting Residues   :",
            summary_data.pop("num_interacting_residues"),
        )

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
                        p_string += " " * (colon_col - len(p_string))
                        print(
                            f"{p_string}: {summary_data[f'{p}%_docking_score']:.2f} kcal/mol"
                        )
                    elif col == "leff":
                        p_string = f"LE     (top {p}% )"
                        p_string += " " * (colon_col - len(p_string))
                        print(
                            f"{p_string}: {summary_data[f'{p}%_leff']:.2f} kcal/mol/heavyatom"
                        )
                    else:
                        p_string = f"{col} (top {p}%)"
                        p_string += " " * (colon_col - len(p_string))
                        print(f"{p_string}: {summary_data[f'{p}%_{col}']:.2f}")

    def filter(
        self,
        eworst=None,
        ebest=None,
        leworst=None,
        lebest=None,
        score_percentile=None,
        le_percentile=None,
        vdw_interactions=None,
        hb_interactions=None,
        reactive_interactions=None,
        hb_count=None,
        react_any=None,
        max_miss=None,
        ligand_name=None,
        ligand_operator=None,
        ligand_substruct=None,
        ligand_substruct_pos=None,
        ligand_max_atoms=None,
        filters_dict: dict | None = None,
        # other processing options:
        enumerate_interaction_combs: bool = False,
        output_all_poses: bool = None,
        mfpt_cluster: float = None,
        interaction_cluster: float = None,
        log_file: str = None,
        overwrite: bool = None,
        order_results: str = None,
        outfields: str = None,
        bookmark_name: str = None,
        filter_bookmark: str = None,
        options_dict: dict | None = None,
        return_iter=False,
    ):
        """Prepare list of filters, then hand it off to storageman to perform filtering. Creates log of all ligand docking results that passes.

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
                hb_count (list[tuple]): accept ligands with at least the requested number of HB interactions. If a negative number is provided, then accept ligands with no more than the requested number of interactions. E.g., [('hb_count', 5)]
                react_any (bool): check if ligand reacted with any residue
                max_miss (int): Will compute all possible combinations of interaction filters excluding up to max_miss numer of interactions from given set. Default will only return union of poses interaction filter combinations. Use with 'enumerate_interaction_combs' for enumeration of poses passing each individual combination of interaction filters.
                ligand_name (list[str]): specify ligand name(s). Will combine name filters with OR, e.g., [["lig1", "lig2"]]
                ligand_substruct (list[str]): SMARTS, index of atom in SMARTS, cutoff dist, and target XYZ coords, e.g., [["ccc", "CN"]]
                ligand_substruct_pos (list[list[type]]): SMARTS pattern(s) for substructure matching, e.g., [["[Oh]C", 0, 1.2, -5.5, 10.0, 15.5]] -> [["smart_string", index_of_positioned_atom, cutoff_distance, x, y, z]]
                ligand_max_atoms (int): Maximum number of heavy atoms a ligand may have
                ligand_operator (str): logical join operator for multiple SMARTS (default: OR), either AND or OR
                filters_dict (dict): provide filters as a dictionary
            Ligand results options:
                enumerate_interaction_combs (bool): When used with `max_miss` > 0, will log ligands/poses passing each separate interaction filter combination as well as union of combinations. Can significantly increase runtime.
                output_all_poses (bool): By default, will output only top-scoring pose passing filters per ligand. This flag will cause each pose passing the filters to be logged.
                mfpt_cluster (float): Cluster filtered ligands by Tanimoto distance of Morgan fingerprints with Butina clustering and output ligand with lowest ligand efficiency from each cluster. Default clustering cutoff is 0.5. Useful for selecting chemically dissimilar ligands.
                interaction_cluster (float): Cluster filtered ligands by Tanimoto distance of interaction fingerprints with Butina clustering and output ligand with lowest ligand efficiency from each cluster. Default clustering cutoff is 0.5. Useful for enhancing selection of ligands with diverse interactions.
                log_file (str): by default, results are saved in `output_log.txt`; if this option is used, ligands and requested info passing the filters will be written to specified file
                overwrite (bool): by default, if a log file exists, it doesn't get overwritten and an error is returned; this option enable overwriting existing log files. Will also overwrite existing database
                order_results (str): Stipulates how to order the results when written to the log file. By default will be ordered by order results were added to the database. ONLY TAKES ONE OPTION. Available fields are:\n
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
                outfields (str): defines which fields are used when reporting the results (to stdout and to the log file); fields are specified as comma-separated values, e.g. `--outfields=e,le,hb`; by default, docking_score (energy) and ligand name are reported; ligand always reported in first column available fields are: \n
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
                filter_bookmark (str): name of bookmark to perform filtering over
                options_dict (dict): write options as a dict
                return_inter (bool): return an iterable of all of the filtering results

        Returns:
            int: number of ligands passing filter
            iter (optional): an iterable of all of the filtering results

        """

        self.set_filters(
            eworst=eworst,
            ebest=ebest,
            leworst=leworst,
            lebest=lebest,
            score_percentile=score_percentile,
            le_percentile=le_percentile,
            vdw_interactions=vdw_interactions,
            hb_interactions=hb_interactions,
            reactive_interactions=reactive_interactions,
            hb_count=hb_count,
            react_any=react_any,
            max_miss=max_miss,
            ligand_name=ligand_name,
            ligand_operator=ligand_operator,
            ligand_substruct=ligand_substruct,
            ligand_substruct_pos=ligand_substruct_pos,
            ligand_max_atoms=ligand_max_atoms,
            dict=filters_dict,
        )

        if options_dict is not None:
            storage_dict, output_dict = split_dict(
                options_dict, ["log_file", "enumerate_interaction_combs"]
            )
        else:
            storage_dict = None
            output_dict = None
        self.set_storageman_attributes(
            output_all_poses=output_all_poses,
            mfpt_cluster=mfpt_cluster,
            interaction_cluster=interaction_cluster,
            overwrite=overwrite,
            order_results=order_results,
            outfields=outfields,
            bookmark_name=bookmark_name,
            filter_bookmark=filter_bookmark,
            dict=storage_dict,
        )
        self.set_output_options(
            log_file=log_file,
            enumerate_interaction_combs=enumerate_interaction_combs,
            dict=output_dict,
        )

        # Compatibility check with docking mode
        if self.docking_mode == "vina" and self.filters.react_any:
            self.logger.warning(
                "Cannot use reaction filters with Vina mode. Removing react_any filter."
            )
            self.filters.react_any = False

        # make sure enumerate_interaction_combs always true if max_miss = 0, since we don't ever worry about the union in this case
        if self.filters.max_miss == 0:
            self.outputopts.enumerate_interaction_combs = True

        # guard against unsing percentile filter with all_poses
        if self.storageopts.output_all_poses and not (
            self.filters.score_percentile is None or self.filters.le_percentile is None
        ):
            self.logger.warning(
                "Cannot return all passing poses with percentile filter. Will only log best pose."
            )
            self.storageopts.output_all_poses = False

        self.logger.info("Filtering results...")
        ligands_passed = 0
        # get possible permutations of interaction with max_miss excluded
        if self.filters.max_miss > 0 and self.outputopts.enumerate_interaction_combs:
            write_one_bookmark = False
        else:
            write_one_bookmark = True

        with self.storageman:
            # pre-process if filtering to multiple bookmark combinations
            if write_one_bookmark:
                filtered_results = self.storageman.filter_results(
                    self.filters.todict(),
                )
                # if there were results of the filtering
                if filtered_results:
                    # if retuning an iterable with the resulting pose ids
                    if return_iter:
                        return filtered_results
                    result_bookmark_name = self.storageman.get_current_bookmark_name()
                    # write output log file
                    with self.outputman:
                        self.outputman.write_filters_to_log(
                            self.filters.todict(),
                            [],
                            f"Morgan Fingerprints butina clustering cutoff: {self.storageman.mfpt_cluster}\nInteraction Fingerprints clustering cutoff: {self.storageman.interaction_cluster}",
                        )
                        self.outputman.write_results_bookmark_to_log(
                            result_bookmark_name
                        )
                        number_passing = self.outputman.write_filter_log(
                            filtered_results
                        )
                        self.outputman.log_num_passing_ligands(number_passing)
                        print("\nNumber of ligands passing filters:", number_passing)
                        ligands_passed = number_passing
                else:
                    self.logger.warning(f"WARNING: No ligands found passing filter.")
                    self.storageman.drop_bookmark(self.storageman.bookmark_name)
            # else produce a bookmark for each interaction combination
            elif not write_one_bookmark:
                interaction_combs = self._generate_interaction_combinations(
                    self.filters.max_miss
                )
                for ic_idx, combination in enumerate(interaction_combs):
                    # prepare Filter object with only desired interaction combination for storageManager
                    filters_dict = self._prepare_filters_for_storageman(combination)
                    # set storageMan's internal ic_counter to reflect current ic_idx
                    if len(interaction_combs) > 1:
                        self.storageman.set_bookmark_suffix(ic_idx)
                    # ask storageManager to fetch results
                    filtered_results = self.storageman.filter_results(
                        filters_dict,
                        not self.outputopts.enumerate_interaction_combs,
                    )
                    if filtered_results:
                        if return_iter:
                            return filtered_results
                        result_bookmark_name = (
                            self.storageman.get_current_bookmark_name()
                        )
                        with self.outputman:
                            self.outputman.write_filters_to_log(
                                self.filters.todict(),
                                combination,
                                f"Morgan Fingerprints butina clustering cutoff: {self.storageman.mfpt_cluster}\nInteraction Fingerprints clustering cutoff: {self.storageman.interaction_cluster}",
                            )
                            self.outputman.write_results_bookmark_to_log(
                                result_bookmark_name
                            )
                            number_passing = self.outputman.write_filter_log(
                                filtered_results
                            )
                            self.outputman.log_num_passing_ligands(number_passing)
                            print(
                                "\nNumber of ligands passing filters:", number_passing
                            )
                            ligands_passed = number_passing
                    elif len(interaction_combs) > 1:
                        self.logger.warning(
                            f"WARNING: No ligands found passing given interaction combination {combination}"
                        )
                        self.storageman.drop_bookmark(self.storageman.bookmark_name)
                if len(interaction_combs) > 1:
                    maxmiss_union_results = self.storageman.get_maxmiss_union(
                        len(interaction_combs)
                    )
                with self.outputman:
                    self.outputman.write_maxmiss_union_header()
                    self.outputman.write_results_bookmark_to_log(
                        self.storageman.bookmark_name + "_union"
                    )
                    number_passing_union = self.outputman.write_filter_log(
                        maxmiss_union_results
                    )
                    self.outputman.log_num_passing_ligands(number_passing_union)
                    print(
                        "\nNumber passing ligands in max_miss union:",
                        number_passing_union,
                    )
                    ligands_passed = number_passing_union

        return ligands_passed

    def write_flexres_pdb(
        self, receptor_polymer, ligname: str, filename: str, bookmark_name: str = None
    ):
        """
        Writes a receptor pdb with flexible residues based on the ligand provided

        Args:
            receptor_polymer (Polymer): version of receptor produced by meeko
            ligname (str): ligand name for which the receptor flexible residue info should be collected
            filename (str): name of the output pdb, extension is optional, will default to '.pdb'
            bookmark_name (str, optional): will use last used bookmark if not specified, will not work in a db without any filtering performed
        """
        # make flexres rdkit mols for ligand-receptor docking
        if bookmark_name is not None:
            self.set_storageman_attributes(bookmark_name=bookmark_name)

        with self.storageman:
            self.storageman.create_temp_table_from_bookmark()
            ligname, ligand_smile, atom_index_map, hydrogen_parents = (
                self.storageman.fetch_single_ligand_output_info(ligname)
            )
            flexible_residues, flexres_atomnames = self.storageman.fetch_flexres_info()
            if flexible_residues != []:  # converts string to list
                flexible_residues = json.loads(flexible_residues)
                flexres_atomnames = json.loads(flexres_atomnames)

            ligand_mol, flexres_mols, _ = self._create_rdkit_mol(
                ligname,
                ligand_smile,
                atom_index_map,
                hydrogen_parents,
                flexible_residues,
                flexres_atomnames,
            )
            if filename:
                # if providing filename, make sure it has .pdb extension
                root, ext = os.path.splitext(filename)
                if not ext:
                    ext = ".pdb"
                path = root + ext
            else:
                # name if after the receptor
                receptor_name, _ = self.storageman.fetch_receptor_objects()[0]
                path = receptor_name + ".pdb"

        flexmoldict = {}
        # string in list of strings
        for index, flexres in enumerate(flexible_residues):
            # res id is chain:resnum
            flexmoldict[f"{flexres[4]}:{flexres[-3:]}"] = flexres_mols[index]

        pdb_str = pdb_updated_flexres_from_rdkit(receptor_polymer, flexmoldict)
        # write pdb string to file
        with open(path, "w") as file:
            file.write(pdb_str)

        return ligand_mol, flexmoldict 

    def write_molecule_sdfs(
        self,
        sdf_path: str | None = None,
        all_in_one: bool = True,
        bookmark_name: str = None,
        write_nonpassing: bool = None,
    ):
        """
        Have output manager write molecule sdf files for passing results in given results bookmark

        Args:
            sdf_path (str, optional): Optional path existing or to be created in cd where SDF files will be saved
            all_in_one (bool, optional): If True will write all molecules to one SDF (separated by $$$$), if False will write one molecule pre SDF
            bookmark_name (str, optional): Option to run over specified bookmark other than that just used for filtering
            write_nonpassing (bool, optional): Option to include non-passing poses for passing ligands

        Raises:
            StorageError: if bookmark or data not found
        """

        if sdf_path is not None:
            self.set_output_options(export_sdf_path=sdf_path)
        else:
            self.set_output_options(export_sdf_path=".")
        try:
            all_mols = self.ligands_rdkit_mol(
                bookmark_name=bookmark_name, write_nonpassing=write_nonpassing
            )
        except StorageError as e:
            self.logger.error(str(e))
            return

        if all_mols is None:
            self.logger.error(
                "Selected bookmark {0} does not exist or does not have any data, cannot write molecule SDFS.".format(
                    bookmark_name
                )
            )
            return

        for ligname, info in all_mols.items():
            # determine filename
            if all_in_one:
                # will write one SDF file for all molecules in bookmark (_None if no bookmark present)
                db_file_name = os.path.splitext(self.db_file)[0]
                sdf_file_name = ("{0}_{1}.sdf").format(
                    db_file_name, str(self.storageman.bookmark_name)
                )
                self.logger.info("Writing " + ligname + " to {0}".format(sdf_file_name))
            else:
                # filename is name of ligand
                sdf_file_name = ligname + ".sdf"
                self.logger.info("Writing " + ligname + ".sdf")

            self.outputman.write_out_mol(
                sdf_file_name,
                info["ligand"],
                info["flex_residues"],
                info["properties"],
            )

    def ligands_rdkit_mol(self, bookmark_name=None, write_nonpassing=False) -> dict:
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
            # Ensure bookmarks exist and have data
            all_bookmarks = self.storageman.get_all_bookmark_names()

            # is bookmark name actually in database
            if self.storageman.bookmark_name in all_bookmarks:
                # check if has max_miss filter
                bookmark_filters = self.storageman.fetch_filters_from_bookmark(
                    self.storageman.bookmark_name
                )
                try:
                    max_miss_present = bool(
                        bookmark_filters["max_miss"] > 0
                        and not bookmark_filters["enumerate_interaction_combs"]
                        and not "_union" in self.storageman.bookmark_name
                    )
                except:
                    #  in case bookmark query string does not contain the phrase 'max_miss', carry on
                    pass
                else:
                    if max_miss_present:
                        self.logger.warning(
                            "'max_miss' used in filtering, but the bookmark used for sdfs writing is not the union of the search"
                        )
            # if bookmark name is not in the database
            elif not self.storageman.bookmark_name in all_bookmarks:
                # does bookmark name + _union resolve the issue
                if self.storageman.check_passing_bookmark_exists(
                    self.storageman.bookmark_name + "_union"
                ):
                    self.storageman.bookmark_name = (
                        self.storageman.bookmark_name + "_union"
                    )
                    self.logger.warning(
                        "Requested 'export_sdf_path' with 'max_miss' and 'enumerate_interaction_combs' used in the filtering process. Exported SDFs will be for union of interaction combinations."
                    )
                # if not, raise error
                else:
                    raise StorageError(
                        "Filtering bookmark {0} does not exist in database. Cannot write passing molecule SDFs".format(
                            self.storageman.bookmark_name
                        )
                    )

            if not self.storageman.bookmark_has_rows(self.storageman.bookmark_name):
                raise StorageError(
                    "Given results bookmark exists but does not have any data. Cannot write passing molecule SDFs"
                )

            # make temp table
            self.storageman.create_temp_table_from_bookmark()
            passing_molecule_info = self.storageman.fetch_passing_ligand_output_info()
            flexible_residues, flexres_atomnames = self.storageman.fetch_flexres_info()

            if flexible_residues != []:
                flexible_residues = json.loads(flexible_residues)
                flexres_atomnames = json.loads(flexres_atomnames)

            all_mols = {}
            for ligname, smiles, atom_indices, h_parent_line in passing_molecule_info:
                self.logger.info("Creating an RDKIT mol for ligand: " + ligname + ".")
                # create rdkit ligand molecule and flexible residue container
                if smiles == "":
                    self.logger.warning(
                        f"No SMILES found for {ligname}. Cannot create SDF."
                    )
                    continue
                # some work needed bc of info needed for creating rd kit, can I remove more than flex stuff?
                mol, flexres_mols, properties = self._create_rdkit_mol(
                    ligname,
                    smiles,
                    atom_indices,
                    h_parent_line,
                    flexible_residues,
                    flexres_atomnames,
                    write_nonpassing=write_nonpassing,
                )

                all_mols[ligname] = {
                    "ligand": mol,
                    "flex_residues": flexres_mols,
                    "properties": properties,
                }

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
        with self.storageman:
            similar_ligands, bookmark_name, cluster_name = (
                self.storageman.fetch_clustered_similars(query_ligname)
            )
            if similar_ligands is not None:
                if not hasattr(self, "outputman"):
                    self.set_output_options()
                with self.outputman:
                    self.outputman.write_find_similar_header(
                        query_ligname, cluster_name
                    )
                    self.outputman.write_results_bookmark_to_log(bookmark_name)
                    number_similar = self.outputman.write_filter_log(similar_ligands)
                    self.outputman.log_num_passing_ligands(number_similar)
                    print("Number similar ligands:", number_similar)
        return number_similar

    def plot(self, save=True, bookmark_name: str = None):
        """
        Get data needed for creating Ligand Efficiency vs
        Energy scatter plot from storageManager. Call OutputManager to create plot.

        Args:
            save (bool): whether to save plot to cd
        """
        if bookmark_name is not None:
            self.set_storageman_attributes(bookmark_name=bookmark_name)
        with self.storageman:
            bookmark_filters = (
                self.storageman.fetch_filters_from_bookmark()
            )  # fetches the filters used to produce the bookmark
        if bookmark_filters:
            max_miss = bookmark_filters["max_miss"]
            if max_miss > 0:
                raise OptionError(
                    "Cannot use --plot with --max_miss > 0. Can plot for desired bookmark with --bookmark_name."
                )

        logger.info("Creating plot of results")
        # get data from storageMan
        with self.storageman:
            all_data, passing_data = self.storageman.get_plot_data()
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

    def display_pymol(self, bookmark_name=None):
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
                plt.plot(line[0], line[1], ".r", mfc="None", picker=5)
                poseIDs[(line[0], line[1])] = (
                    line[2],
                    line[3],
                )  # line[0] is LE, line[1] is docking score, line[2] is pose_ID, line[3] is LigName
            plt.ylabel("Ligand Efficiency (kcal/mol/heavy atom)")
            plt.xlabel("Docking Score (kcal/mol)")
            plt.title("Passing Docking Poses")

            try:
                pymol = PyMol.MolViewer()
            except ConnectionRefusedError as e:
                raise RTCoreError(
                    "Error establishing connection with PyMol. Try manually launching PyMol with `pymol -R` in another terminal window."
                ) from e

            def onpick(event):
                line = event.artist
                coords = tuple([c[0] for c in line.get_data()])
                chosen_pose = poseIDs[coords]
                self.logger.info(
                    f"LigName: {chosen_pose[1]}; Pose_ID: {chosen_pose[0]}"
                )

                # make rdkit mol for poseid
                ligname, ligand_smile, atom_index_map, hydrogen_parents = (
                    self.storageman.fetch_single_ligand_output_info(chosen_pose[1])
                )
                flexible_residues, flexres_atomnames = (
                    self.storageman.fetch_flexres_info()
                )
                if flexible_residues != []:  # converts string to list
                    flexible_residues = json.loads(flexible_residues)
                    flexres_atomnames = json.loads(flexres_atomnames)

                mol, flexres_mols, _ = self._create_rdkit_mol(
                    ligname,
                    ligand_smile,
                    atom_index_map,
                    hydrogen_parents,
                    flexible_residues,
                    flexres_atomnames,
                    pose_ID=chosen_pose[0],
                )
                self.logger.debug(Chem.MolToSmiles(mol))
                pymol.ShowMol(mol, name=ligname, showOnly=False)
                for idx, resmol in enumerate(flexres_mols):
                    pymol.ShowMol(
                        resmol,
                        name=ligname + "_" + flexible_residues[idx],
                        showOnly=False,
                    )

            fig = plt.gcf()
            cid = fig.canvas.mpl_connect("pick_event", onpick)
            plt.show()

    def export_csv(self, requested_data: str, csv_name: str, table=False):
        """Get requested data from database, export as CSV

        Args:
            requested_data (str): Table name or SQL-formatted query
            csv_name (str): Name for exported CSV file
            table (bool): flag indicating is requested data is a table name
        """
        with self.storageman:
            df = self.storageman.to_dataframe(requested_data, table=table)
            df.to_csv(csv_name)

    def export_bookmark_db(self, bookmark_name: str = None) -> str:
        """Export database containing data from bookmark

        Args:
            bookmark_name (str): name for bookmark_db

        Returns:
            str: name of the new, exported database
        """
        if bookmark_name is not None:
            self.set_storageman_attributes(bookmark_name=bookmark_name)
        bookmark_db_name = (
            self.db_file.rstrip(".db") + "_" + self.storageman.bookmark_name + ".db"
        )
        self.logger.info("Exporting bookmark database")
        if os.path.exists(bookmark_db_name):
            self.logger.warning(
                "Requested export DB name already exists. Please rename or remove existing database. New database not exported."
            )
            return
        with self.storageman:
            self.storageman.clone(bookmark_db_name)
        # connect to cloned database
        dictionary = self.storageopts.todict()
        dictionary["db_file"] = bookmark_db_name
        temp_storageman = StorageManager.check_storage_compatibility(self.storagetype)
        with temp_storageman(**dictionary) as db_clone:
            db_clone.prune()
            db_clone.close_storage(vacuum=True)

        return bookmark_db_name

    def export_receptors(self):
        """
        Export receptor in database to pdbqt
        """
        with self.storageman:
            receptor_tuples = self.storageman.fetch_receptor_objects()
        for recname, recblob in receptor_tuples:
            if recblob is None:
                self.logger.warning(
                    f"No receptor pdbqt stored for {recname}. Export failed."
                )
                continue
            if not hasattr(self, "outputman"):
                self.set_output_options()
            self.outputman.write_receptor_pdbqt(recname, recblob)

    def get_previous_filter_data(
        self, outfields=None, bookmark_name=None, log_file=None
    ):
        """Get data requested in self.out_opts['outfields'] from the
        results bookmark of a previous filtering

        Args:
            outfields (str): use outfields as described in RingtailOptions > StorageOptions
            bookmark_name (str): bookmark for which the filters were used
        """
        if bookmark_name is not None:
            self.set_storageman_attributes(bookmark_name=bookmark_name)
        if outfields is not None:
            self.set_storageman_attributes(outfields=outfields)
        if log_file is not None:
            self.set_output_options(log_file=log_file)

        with self.storageman:
            new_data = self.storageman.fetch_data_for_passing_results()
        with self.outputman:
            self.outputman.write_filter_log(new_data)

    def drop_bookmark(self, bookmark_name: str):
        """Drops specified bookmark from the database

        Args:
            bookmark_name (str): name of bookmark to be dropped.
        """

        with self.storageman:
            self.storageman.drop_bookmark(bookmark_name)
        self.logger.info(
            "Bookmark {0} was dropped from the database {1}".format(
                bookmark_name, self.storageman.db_file
            )
        )

    def get_bookmark_names(self):
        """
        Method to retrieve all bookmark names in a database

        Returns:
            list: of all bookmarks in a database
        """
        with self.storageman:
            return self.storageman.get_all_bookmark_names()

    @staticmethod
    def default_dict() -> dict:
        """
        Creates a dict of all Ringtail options.

        Return:
            str: json string with options
        """
        defaults = {}
        defaults.update(OutputOptions().todict())
        defaults.update(ResultsProcessingOptions().todict())
        defaults.update(StorageOptions().todict())
        defaults.update(Filters().todict())
        defaults.update(InputFiles().todict())
        defaults.update(ReadOptions().todict())
        defaults.update(GeneralOptions().todict())

        return defaults

    @staticmethod
    def generate_config_file_template():
        """Outputs to "config.json in current working directory if to_file = true,
        else it returns the dict of default option values used for API (for command
        line a few more options are included that are always used explicitly when using API)

        Args:
            to_file (bool): whether to produce the template as a json string or as a file "config.json"

        Returns:
            str: file name of config file or json string with template including default values
        """

        json_string = RingtailCore.default_dict()

        filename = "config.json"
        with open(filename, "w") as f:
            f.write(json.dumps(json_string, indent=4))
        return filename

    @staticmethod
    def get_options_info() -> dict:
        """
        Gets names, default values, and meta data for all Ringtail options.
        """
        options = {}
        options.update(OutputOptions.options)
        options.update(ResultsProcessingOptions.options)
        options.update(StorageOptions.options)
        options.update(Filters.options)
        options.update(InputFiles.options)
        options.update(ReadOptions.options)
        options.update(GeneralOptions.options)

        return options

    # endergion
