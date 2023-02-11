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
from .filters import Filters
from .exceptions import RTCoreError, OutputError
from rdkit import Chem
import itertools
import logging
import os


class RingtailCore:
    """Core class for coordinating different actions on virtual screening
    i.e. adding results to storage, filtering, output options

    Attributes:
        storageman (storageManager): Interface module with database
        eworst (float): The worst scoring energy filter value requested by user
        filter_file (string): Name of file containing filters provided by user
        filtered_results (storage cursor object): Cursor object
            containing results passing requested filters (iterable)
        filters (dictionary): Dictionary containing user-specified filters
        no_print_flag (boolean): Flag specifying whether passing results
            should be printed to terminal
        out_opts (dictionary): Specified output options including data fields
            to output, export_sdf_path, log file name
        output_manager (OutputManager object): Manager for output tasks of
            log-writting, plotting, ligand SDF writing
        results_filters_list (List): List of tuples of filter option and value
        results_man (ResultsManager object): Manager for processing result
            files for insertion into database
    """

    def __init__(self, storage_type = "sqlite", opts_dict: dict = None):
        """Initialize RingtailCore object. Will create storageManager object to serve
        as interface with database (currently implemented in SQLite).
        Will create ResultsManager to process result files.
        Will create OutputManager object to assist in creating output files.

        Args:
            storage_opts (dictionary): dictionary of options required by storageManager
            rman_opts (dictionary): dictionary of options required by
                results manager
            filters (dictionary): Dictionary containing user-specified filters
            out_opts (dictionary): Specified output options including data
                fields to output, export_sdf_path, log file name
        """
        # check if given opt_dict and set attrbutes if so
        if opts_dict is None:
            self.storage_opts = {}
            self.rman_opts = {}
            self.filters_dict = {}
            self.out_opts = {}
        else:
            self.storage_opts = opts_dict["storage_opts"]["values"]
            self.rman_opts = opts_dict["rman_opts"]["values"]
            self.filters_dict = opts_dict["filters"]["values"]
            self.out_opts = opts_dict["out_opts"]["values"]

        # set storage type
        self.storage_type = storage_type
        storage_types = {
            "sqlite": StorageManagerSQLite,
        }
        # confirm given storage type is implemented
        if self.storage_type not in storage_types:
            raise NotImplementedError(f"Given storage type {self.storage_type} is not implemented.")
        
        # initialize "worker" classes
        self.storageman = storage_types[storage_type](**self.storage_opts)
        self.results_man = ResultsManager(storageman=self.storageman, storageman_class=storage_types[storage_type], **self.rman_opts)
        self.output_manager = OutputManager(**self.out_opts)
        self.filters = Filters(**self.filters_dict)

    @classmethod
    def get_defaults(cls, storage_type="sqlite", terse=False) -> dict:

        cls.__default_options = {
            "storage_opts": {
                "values": StorageManager.get_defaults(storage_type),
                "ignore": [],
                "types": StorageManager.get_default_types(storage_type),
            },
            "rman_opts": {
                "values": ResultsManager.get_defaults(),
                "ignore": [],
                "types": ResultsManager.get_default_types(),
            },
            "out_opts": {
                "values": OutputManager.get_defaults(),
                "ignore": [],
                "types": OutputManager.get_default_types(),
            },
            "filters_dict": {
                "values": Filters.get_defaults(),
                "ignore": [],
                "types": Filters.get_default_types(),
            },
        }

        defaults = cls.__default_options.copy()
        if terse:
            for group, opt in defaults.items():
                defaults[group] = {k: v for k, v in opt.items() if not k == "ignore"}

        return defaults

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close_database()

    def add_results(self):
        """
        Call results manager to process result files and add to database
        """
        # check that we want to overwrite or add results if storage already exists
        if (
            not self.storage_opts["overwrite"]
            and not self.storage_opts["append_results"]
        ):
            self.storageman.check_storage_empty()
        logging.info("Adding results...")
        self.results_man.process_results()

    def save_receptors(self, receptor_file):
        """Add receptor to database

        Args:
            receptors (list): list of receptor blobs to add to database
        """
        receptor_list = ReceptorManager.make_receptor_blobs([receptor_file])
        for rec, rec_name in receptor_list:
            # NOTE: in current implementation, only one receptor allowed per database
            # Check that any receptor row is incomplete (needs receptor blob) before inserting
            filled_receptor_rows = self.storageman.check_receptors_saved()
            if filled_receptor_rows != 0:
                raise RTCoreError(
                    "Expected Receptors table to have no receptor objects present, already has {0} receptor present. Cannot add more than 1 receptor to a database.".format(
                        filled_receptor_rows
                    )
                )
            self.storageman.save_receptor(rec)

    def produce_summary(self, columns=["docking_score", "leff"], percentiles=[1, 10]) -> None:
        """Print summary of data in storage
        """
        summary_data = self.storageman.fetch_summary_data(columns, percentiles)
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

    def display_pymol(self):
        """launch pymol session and plot of LE vs docking score. Displays molecules when clicked
        """

        import subprocess
        from rdkit.Chem import PyMol

        # launch pymol session
        p = subprocess.Popen(
        ["pymol", "-R"],
        stdout=subprocess.PIPE,
        )

        poseIDs = {}

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
            logging.info(f"LigName: {chosen_pose[1]}; Pose_ID: {chosen_pose[0]}")

            # make rdkit mol for poseid
            ligname, ligand_smile, atom_index_map, hydrogen_parents = self.storageman.fetch_single_ligand_output_info(chosen_pose[1])
            flexible_residues, flexres_atomnames = self.storageman.fetch_flexres_info()
            if flexible_residues != []:  # converts string to list
                flexible_residues = json.loads(flexible_residues)
                flexres_atomnames = json.loads(flexres_atomnames)

            mol, flexres_mols, _ = self.create_ligand_rdkit_mol(ligname, ligand_smile, atom_index_map, hydrogen_parents, flexible_residues, flexres_atomnames, pose_ID=chosen_pose[0])
            logging.debug(Chem.MolToSmiles(mol))
            pymol.ShowMol(mol, name=ligname, showOnly=False)
            for idx, resmol in enumerate(flexres_mols):
                pymol.ShowMol(resmol, name=ligname + "_" + flexible_residues[idx], showOnly=False)

        fig = plt.gcf()
        cid = fig.canvas.mpl_connect('pick_event', onpick)
        plt.show()
    
    def filter(self):
        """
        Prepare list of filters, then hand it off to storageManager to
            perform filtering. Create log of passing results.
        """

        logging.info("Filtering results...")
        self.output_manager.create_log_file()
        # get possible permutations of interaction with max_miss excluded
        interaction_combs = self._generate_interaction_combinations(
            self.filters.max_miss
        )

        for ic_idx, combination in enumerate(interaction_combs):
            # prepare list of filter values and keys for storageManager
            results_filters_list = self.prepare_results_filter_list(combination)

            # set storageMan's internal ic_counter to reflect current ic_idx
            if len(interaction_combs) > 1:
                self.storageman.set_view_suffix(str(ic_idx))
            # ask storageManager to fetch results
            self.filtered_results = self.storageman.filter_results(
                self.filters,
                results_filters_list,
                self.filters["ligand_filters"],
                self.out_opts["outfields"],
                self.out_opts["filter_bookmark"],
            )
            number_passing_ligands = self.storageman.get_number_passing_ligands()
            print("Number passing Ligands:", number_passing_ligands)
            result_bookmark_name = self.storageman.get_current_view_name()
            self.output_manager.write_filters_to_log(self.filters, combination)
            self.output_manager.write_results_bookmark_to_log(result_bookmark_name)
            self.output_manager.log_num_passing_ligands(number_passing_ligands)
            self.output_manager.write_log(self.filtered_results)

    def get_previous_filter_data(self):
        """Get data requested in self.out_opts['outfields'] from the
        results view of a previous filtering
        """
        self.output_manager.create_log_file()
        new_data = self.storageman.fetch_data_for_passing_results(
            self.out_opts["outfields"]
        )
        self.output_manager.write_log(new_data)

    def plot(self, save=True):
        """
        Get data needed for creating Ligand Efficiency vs
        Energy scatter plot from storageManager. Call OutputManager to create plot.
        """

        logging.info("Creating plot of results")
        # get data from storageMan
        all_data, passing_data = self.storageman.get_plot_data()
        all_plot_data_binned = dict()
        # bin the all_ligands data by 1000ths to make plotting faster
        for line in all_data:
            # add to dictionary as bin of energy and le
            data_bin = (round(line[0], 3), round(line[1], 3))
            if data_bin not in all_plot_data_binned:
                all_plot_data_binned[data_bin] = 1
            else:
                all_plot_data_binned[data_bin] += 1
        # plot the data
        self.output_manager.plot_all_data(all_plot_data_binned)
        if passing_data != []:  # handle if no passing ligands
            for line in passing_data:
                self.output_manager.plot_single_point(
                    line[0], line[1], "red"
                )  # energy (line[0]) on x axis, le (line[1]) on y axis
        if save:
            self.output_manager.save_scatterplot()
        else:
            plt.show()

    def prepare_results_filter_list(self, included_interactions):
        """takes filters dictionary from option parser.
        Output list of tuples to be inserted into sql call string

        Args:
            included_interactions (tuple): Tuple of interactions to include in filter
        """

        filters_list = []

        # get property filters
        properties_keys = [
            "eworst",
            "ebest",
            "leworst",
            "lebest",
            "score_percentile",
            "le_percentile",
        ]

        for key in properties_keys:
            if getattr(self.filters, key) is not None:
                filters_list.append(getattr(self.filters, key))

        # interaction filters
        interaction_filters = {"V": self.filters.vdw_interactions, "H": self.filters.hb_interactions, "R": self.filters.reactive_interactions}
        for key in interaction_filters:
            if interaction_filters[key] is not None:
                kept_interactions = []
                for interaction in interaction_filters[key]:
                    # only keep interactions specified by included_interactions
                    if key + "-" + interaction[0] in included_interactions:
                        kept_interactions.append(interaction)
                filters_list.append((key, kept_interactions))

        # get interaction count filters
        for count in self.filters.interactions_count:
            filters_list.append(count)  # already a tuple, don't need to reformat

        # add react_any flag
        filters_list.append(("react_any", self.filters.react_any))

        return filters_list

    def write_molecule_sdfs(self, write_nonpassing=False):
        """have output manager write sdf molecules for passing results in given results bookmark

        Args:
            write_nonpassing (bool, optional): Option to include non-passing poses for passing ligands
        """

        if not self.storageman.check_passing_view_exists():
            logging.warning(
                "Given results bookmark does not exist in database. Cannot write passing molecule SDFs"
            )
            return
        passing_molecule_info = self.storageman.fetch_passing_ligand_output_info()
        flexible_residues, flexres_atomnames = self.storageman.fetch_flexres_info()
        if flexible_residues != []:
            flexible_residues = json.loads(flexible_residues)
            flexres_atomnames = json.loads(flexres_atomnames)
        for (ligname, smiles, atom_indices, h_parent_line) in passing_molecule_info:
            logging.info("Writing " + ligname.split(".")[0] + ".sdf")
            # create rdkit ligand molecule and flexible residue container
            if smiles == "":
                logging.warning(
                    f"No SMILES found for {ligname}. Cannot create SDF."
                )
                continue

            mol, flexres_mols, properties = self.create_ligand_rdkit_mol(ligname, smiles, atom_indices, h_parent_line, flexible_residues, flexres_atomnames, write_nonpassing=write_nonpassing)

            # write out mol
            self.output_manager.write_out_mol(
                ligname, mol, flexres_mols, properties
            )

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
            flexres_mols.append(Chem.MolFromSmiles(res_smiles))
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
    
    def export_csv(self, requested_data: str, csv_name: str, table=False):
        """Get requested data from database, export as CSV

        Args:
            requested_data (string): Table name or SQL-formatted query
            csv_name (string): Name for exported CSV file
            table (bool): flag indicating is requested data is a table name
        """
        df = self.storageman.to_dataframe(requested_data, table=table)
        df.to_csv(csv_name)

    def export_bookmark_db(self, bookmark_db_name: str):
        """Export database containing data from bookmark

        Args:
            bookmark_db_name (str): name for bookmark_db
        """
        logging.info("Exporting bookmark database")
        if os.path.exists(bookmark_db_name):
            logging.warning(
                "Requested export DB name already exists. Please rename or remove existing database. New database not exported."
            )
            return
        self.storageman.clone(bookmark_db_name)
        # connect to cloned database
        db_clone = StorageManagerSQLite(bookmark_db_name, self.storage_opts)
        db_clone.prune()
        db_clone.close_connection(vacuum=True)

    def close_database(self):
        """Tell database we are done and it can close the connection"""
        self.storageman.close_connection()

    def _clean_storage_string(self, input_str):
        """take a storage string representing a list,
        strips unwanted characters

        Args:
            input_str (str)

        Returns:
            String: cleaned string
        """
        return (
            input_str.replace("[", "")
            .replace("]", "")
            .replace("'", "")
            .replace('"', "")
            .replace("\n", "")
            .replace("\\n", "")
            .replace(" \n", "")
        )

    def _generate_pdbqt_block(self, pdbqt_lines):
        """Generate pdbqt block from given lines from a pdbqt

        Args:
            flexres_lines (TYPE): Description
        """

        return "\n".join(
            [
                line.lstrip(" ")
                for line in list(
                    filter(
                        None, [self._clean_storage_string(line) for line in pdbqt_lines]
                    )
                )
            ]
        )

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
            print(ligand_pose)
            mol = RDKitMolCreate.add_pose_to_mol(mol, ligand_pose, atom_indices)
            print(Chem.MolToMolBlock(mol))
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
        for _type, interactions in {"V": self.filters.vdw_interactions, "H": self.filters.hb_interactions, "R": self.filters.reactive_interactions}.items():
            for interact in interactions:
                all_interactions.append(_type + "-" + interact[0])

        # warn if max_miss greater than number of interactions
        if max_miss > len(all_interactions):
            logging.warning(
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
