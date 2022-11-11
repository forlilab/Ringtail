#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail virtual screening manager
#

from cProfile import label
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import colors
import numpy as np
import json
from meeko import RDKitMolCreate
from .storagemanager import StorageManager, StorageManagerSQLite
from .resultsmanager import ResultsManager
from .receptormanager import ReceptorManager
from .exceptions import StorageError
from .exceptions import RTCoreError, OutputError
from rdkit import Chem
from rdkit.Chem import SDWriter
import itertools
import logging
import os
import typing


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

    def __init__(self, **opts):
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
        defaults = self.get_defaults()

        # unpack given opts dictionaries
        for k, v in opts.items():
            if k in defaults.keys():
                defaults[k] = v

        self.storage_opts = defaults["storage_opts"]["values"]
        self.rman_opts = defaults["rman_opts"]["values"]
        self.filters = defaults["filters"]["values"]
        self.out_opts = defaults["out_opts"]["values"]

        storage_type = self.storage_opts.pop("storage_type")
        # confirm storage_type is not None
        if storage_type is None:
            raise RTCoreError(
                "storage_type not defined and must be defined at runtime."
            )

        storage_types = {
            "sqlite": StorageManagerSQLite,
        }
        self.storageman = storage_types[storage_type](**self.storage_opts)
        self.rman_opts["storageman"] = self.storageman
        self.rman_opts["storageman_class"] = storage_types[storage_type]
        self.output_manager = None

    @classmethod
    def get_defaults(cls, storage_type="sqlite", terse=False):

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
        }

        # add storage type default
        cls.__default_options["storage_opts"]["values"]["storage_type"] = storage_type
        cls.__default_options["storage_opts"]["types"]["storage_type"] = str

        out_opts = {
            "log": "output_log.txt",
            "export_sdf_path": None,
            "plot": None,
            "outfields": "e",
            "no_print": True,
            "export_bookmark_csv": None,
            "export_query_csv": None,
            "export_bookmark_db": None,
            "data_from_bookmark": None,
            "filter_bookmark": None,
        }

        out_types = {
            "log": str,
            "export_sdf_path": str,
            "plot": bool,
            "outfields": str,
            "no_print": bool,
            "export_bookmark_csv": str,
            "export_query_csv": str,
            "export_bookmark_db": str,
            "data_from_bookmark": str,
            "filter_bookmark": str,
        }

        filters = {
            "properties": {
                "eworst": None,
                "ebest": None,
                "leworst": None,
                "lebest": None,
                "energy_percentile": None,
                "le_percentile": None,
            },
            "interactions": {"V": [], "H": [], "R": []},
            "interactions_count": [],
            "ligand_filters": {"N": [], "S": [], "F": "OR"},
            "filter_ligands_flag": False,
            "max_miss": 0,
            "react_any": None,
        }

        filter_types = {
            "properties": {
                "eworst": float,
                "ebest": float,
                "leworst": float,
                "lebest": float,
                "energy_percentile": float,
                "le_percentile": float,
            },
            "interactions": dict,
            "interactions_count": list,
            "ligand_filters": dict,
            "filter_ligands_flag": bool,
            "max_miss": int,
            "react_any": bool,
        }

        defaults = cls.__default_options.copy()
        if terse:
            for group, opt in defaults.items():
                defaults[group] = {k: v for k, v in opt.items() if not k == "ignore"}
        defaults["out_opts"] = {"values": out_opts, "ignore": [], "types": out_types}
        defaults["filters"] = {"values": filters, "ignore": [], "types": filter_types}

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
        self.results_man = ResultsManager(**self.rman_opts)
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
            try:
                filled_receptor_rows = self.storageman.check_receptors_saved()
                if filled_receptor_rows != 0:
                    raise RTCoreError(
                        "Expected Receptors table to have no receptor objects present, already has {0} receptor present. Cannot add more than 1 receptor to a database.".format(
                            filled_receptor_rows
                        )
                    )
                self.storageman.save_receptor(rec)
            except Exception as e:
                raise e

    def filter(self):
        """
        Prepare list of filters, then hand it off to storageManager to
            perform filtering. Create log of passing results.
        """

        logging.info("Filtering results...")
        self._prepare_output_manager()
        # get possible permutations of interaction with max_miss excluded
        interaction_combs = self._generate_interaction_combinations(
            self.filters["max_miss"]
        )

        for ic_idx, combination in enumerate(interaction_combs):
            # prepare list of filter values and keys for storageManager
            results_filters_list = self.prepare_results_filter_list(combination)

            # make sure we have ligand filter list
            if not self.filters["filter_ligands_flag"]:
                self.filters["ligand_filters"] = []
            # set storageMan's internal ic_counter to reflect current ic_idx
            if len(interaction_combs) > 1:
                self.storageman.set_view_suffix(str(ic_idx))
            # ask storageManager to fetch results
            try:
                self.filtered_results = self.storageman.filter_results(
                    results_filters_list,
                    self.filters["ligand_filters"],
                    self.out_opts["outfields"],
                    self.out_opts["filter_bookmark"],
                )
                number_passing_ligands = self.storageman.get_number_passing_ligands()
                result_bookmark_name = self.storageman.get_current_view_name()
                self.output_manager.write_filters_to_log(self.filters, combination)
                self.output_manager.write_results_bookmark_to_log(result_bookmark_name)
                self.output_manager.log_num_passing_ligands(number_passing_ligands)
                self.output_manager.write_log(self.filtered_results)
            except StorageError as e:
                logging.exception("Database error occurred while filtering")
                raise e
            except OutputError as e:
                logging.exception("Logging error occurred after filtering")
                raise e
            except Exception as e:
                raise e

    def get_previous_filter_data(self):
        """Get data requested in self.out_opts['outfields'] from the
        results view of a previous filtering
        """
        try:
            self._prepare_output_manager()
            new_data = self.storageman.fetch_data_for_passing_results(
                self.out_opts["outfields"]
            )
            self.output_manager.write_log(new_data)
        except StorageError as e:
            logging.exception(
                "Database error occurred while fetching data for bookmark"
            )
            raise e
        except OutputError as e:
            logging.exception("Error occurred while writing log")
            raise e
        except Exception as e:
            raise e

    def plot(self):
        """
        Get data needed for creating Ligand Efficiency vs
        Energy scatter plot from storageManager. Call OutputManager to create plot.
        """
        try:
            logging.info("Creating plot of results")
            self._prepare_output_manager()
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
            self.output_manager.save_scatterplot()
        except StorageError as e:
            logging.exception("Database error occurred while fetching data for plot")
            raise e
        except Exception as e:
            logging.exception("Error occurred during plotting")
            raise e

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
            "energy_percentile",
            "le_percentile",
        ]

        property_filters = self.filters["properties"]
        for key in properties_keys:
            if property_filters[key] is not None:
                filters_list.append((key, property_filters[key]))

        # interaction filters
        interaction_filters = self.filters["interactions"]
        for key in interaction_filters:
            if interaction_filters[key] is not None:
                kept_interactions = []
                for interaction in interaction_filters[key]:
                    # only keep interactions specified by included_interactions
                    if key + "-" + interaction[0] in included_interactions:
                        kept_interactions.append(interaction)
                filters_list.append((key, kept_interactions))

        # get interaction count filters
        interact_count_filters = self.filters["interactions_count"]
        for count in interact_count_filters:
            filters_list.append(count)  # already a tuple, don't need to reformat

        # add react_any flag
        filters_list.append(("react_any", self.filters["react_any"]))

        return filters_list

    def write_molecule_sdfs(self, write_nonpassing=False):
        """have output manager write sdf molecules for passing results in given results bookmark

        Args:
            write_nonpassing (bool, optional): Option to include non-passing poses for passing ligands
        """
        try:
            self._prepare_output_manager()
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
                        passing_properties,
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

                # write out mol
                self.output_manager.write_out_mol(
                    ligname, mol, flexres_mols, properties
                )

        except StorageError as e:
            logging.exception(
                "Error occurred while fetching database information for SDF output"
            )
            raise e
        except OutputError as e:
            logging.exception(f"Error occured while writing {ligname} to an SDF")
            raise e
        except Exception as e:
            logging.exception("Error occurred during SDF output")
            raise e

    def export_csv(self, requested_data: str, csv_name: str, table=False):
        """Get requested data from database, export as CSV

        Args:
            requested_data (string): Table name or SQL-formatted query
            csv_name (string): Name for exported CSV file
            table (bool): flag indicating is requested data is a table name
        """
        try:
            df = self.storageman.to_dataframe(requested_data, table=table)
            df.to_csv(csv_name)
        except StorageError as e:
            logging.exception(
                f"Error occured while getting data for exporting CSV of {requested_data}"
            )
            raise e
        except Exception as e:
            logging.exception(f"Error occured while exporting CSV of {requested_data}")
            raise e

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

    def _prepare_output_manager(self):
        if self.output_manager is None:
            self.output_manager = OutputManager(
                self.out_opts["log"], self.out_opts["export_sdf_path"]
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
            energies_binding,
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
            properties["Binding energies"].append(energies_binding)
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
        for _type, interactions in self.filters["interactions"].items():
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


# # # # # # # # # # # # # #
# # # Output class # # #
# # # # # # # # # # # # # #


class OutputManager:
    """Class for creating outputs

    Attributes:
        ad_to_std_atomtypes (dictionary): Mapping of autodock atomtypes (key)
            to standard PDB atomtypes (value). Initialized as None, will be
            loaded from AD_to_STD_ATOMTYPES.json if required
        ax (pyplot axis): Axis for scatterplot of LE vs Energy
        conformer_indices (list): Description
        fig_base_name (str): Name for figures excluding file extension
        flex_residue_smiles (Dictionary): Contains smiles used to create
            RDKit objects for flexible residues
        log (string): name for log file
        vsman (VSManager): VSManager object that created output_manager

    """

    def __init__(self, log_file="output_log.txt", export_sdf_path=""):
        """Initialize OutputManager object and create log file

        Args:
            log_file (string): name for log file
            export_sdf_path (string): path for exporting sdf files
        """

        self.log = log_file
        self.export_sdf_path = export_sdf_path

        self._create_log_file()

    def plot_all_data(self, binned_data):
        """takes dictionary of binned data where key is the
        coordinates of the bin and value is the number of points in that bin.
        Adds to scatter plot colored by value

        Args:
            binned_data (dict): Keys are tuples of key and y value for bin.
                Value is the count of points falling into that bin.
        """
        try:
            # gather data
            energies = []
            leffs = []
            bin_counts = []
            for data_bin in binned_data.keys():
                energies.append(data_bin[0])
                leffs.append(data_bin[1])
                bin_counts.append(binned_data[data_bin])

            # start with a square Figure
            fig = plt.figure()

            gs = fig.add_gridspec(
                2,
                2,
                width_ratios=(7, 2),
                height_ratios=(2, 7),
                left=0.1,
                right=0.9,
                bottom=0.1,
                top=0.9,
                wspace=0.05,
                hspace=0.05,
            )

            self.ax = fig.add_subplot(gs[1, 0])
            ax_histx = fig.add_subplot(gs[0, 0], sharex=self.ax)
            ax_histy = fig.add_subplot(gs[1, 1], sharey=self.ax)
            fig.colorbar(
                mappable=cm.ScalarMappable(
                    colors.Normalize(vmin=min(bin_counts), vmax=max(bin_counts))
                ),
                label="Scatterplot bin count",
            )
            self.ax.set_xlabel("Best Binding Energy / kcal/mol")
            self.ax.set_ylabel("Best Ligand Efficiency")
        except Exception as e:
            raise OutputError("Error occurred while initializing plot") from e

        self.scatter_hist(energies, leffs, bin_counts, self.ax, ax_histx, ax_histy)

    def plot_single_point(self, x, y, color="black"):
        """Add point to scatter plot with given x and y coordinates and color.

        Args:
            x (float): x coordinate
            y (float): y coordinate
            color (str, optional): Color for point. Default black.
        """
        try:
            self.ax.scatter([x], [y], c=color)
        except Exception as e:
            raise OutputError("Error occurred while plotting") from e

    def save_scatterplot(self):
        """
        Saves current figure as scatter.png
        """
        try:
            plt.savefig("scatter.png", bbox_inches="tight")
            plt.close()
        except Exception as e:
            raise OutputError("Error while saving figure") from e

    def write_log(self, lines):
        """Writes lines from results iterable into log file

        Args:
            lines (iterable): Iterable with tuples of data for
                writing into log
        """
        try:
            for line in lines:
                logging.info(line)
                self._write_log_line(
                    str(line).replace("(", "").replace(")", "")
                )  # strip parens from line, which is natively a tuple
            self._write_log_line("***************\n")
        except Exception as e:
            raise OutputError("Error occurred during log writing") from e

    def _write_log_line(self, line):
        """write a single row to the log file

        Args:
            line (string): Line to write to log
        """
        try:
            with open(self.log, "a") as f:
                f.write(line)
                f.write("\n")
        except Exception as e:
            raise OutputError(f"Error writing line {line} to log") from e

    def log_num_passing_ligands(self, number_passing_ligands):
        """
        Write the number of ligands which pass given filter to log file

        Args:
            number_passing_ligands (int): number of ligands that passed filter
        """
        try:
            with open(self.log, "a") as f:
                f.write("\n")
                f.write(
                    "Number passing ligands: {num} \n".format(
                        num=str(number_passing_ligands)
                    )
                )
                f.write("---------------\n")
        except Exception as e:
            raise OutputError("Error writing number of passing ligands in log") from e

    def write_results_bookmark_to_log(self, bookmark_name):
        """Write the name of the result bookmark into log

        Args:
            bookmark_name (string): name of current results' bookmark in db
        """
        try:
            with open(self.log, "a") as f:
                f.write("\n")
                f.write(f"Result bookmark name: {bookmark_name}\n")
        except Exception as e:
            raise OutputError("Error writing bookmark name to log") from e

    def write_out_mol(self, ligname, mol, flexres_mols, properties):
        """writes out given mol as sdf

        Args:
            ligname (string): name of ligand that will be used to
                name output SDF file
            mol (RDKit mol object): RDKit molobject to be written to SDF
            flexres_mols (list): dictionary of rdkit molecules for flexible residues
            properties (dict): dictionary of list of properties to add to mol before writing
        """
        try:
            filename = self.export_sdf_path + ligname + ".sdf"
            mol_flexres_list = [mol]
            mol_flexres_list += flexres_mols
            mol = RDKitMolCreate.combine_rdkit_mols(mol_flexres_list)
            # convert properties to strings as needed
            for k, v in properties.items():
                if isinstance(v, list):
                    v = json.dumps(v)
                elif not isinstance(v, str):
                    v = str(v)
                mol.SetProp(k, v)

            with SDWriter(filename) as w:
                for conf in mol.GetConformers():
                    w.write(mol, conf.GetId())

        except Exception as e:
            raise OutputError("Error occurred while writing SDF from RDKit Mol") from e

    def _create_log_file(self):
        """
        Initializes log file
        """
        try:
            with open(self.log, "w") as f:
                f.write("Filtered poses:\n")
                f.write("***************\n")
        except Exception as e:
            raise OutputError("Error while creating log file") from e

    def scatter_hist(self, x, y, z, ax, ax_histx, ax_histy):
        """
        Makes scatterplot with a histogram on each axis

        Args:
            x (list): x coordinates for data
            y (list): y coordinates for data
            z (list): z coordinates for data
            ax (matplotlib axis): scatterplot axis
            ax_histx (matplotlib axis): x histogram axis
            ax_histy (matplotlib axis): y histogram axis
        """
        try:
            # no labels
            ax_histx.tick_params(axis="x", labelbottom=False)
            ax_histy.tick_params(axis="y", labelleft=False)

            # the scatter plot:
            ax.scatter(x, y, c=z, cmap="viridis")

            # now determine nice limits by hand:
            xbinwidth = 0.25
            ybinwidth = 0.01
            xminlim = (int(min(x) / xbinwidth) + 3) * xbinwidth
            xmaxlim = (int(max(x) / xbinwidth) + 3) * xbinwidth
            yminlim = (int(min(y) / ybinwidth) + 3) * ybinwidth
            ymaxlim = (int(max(y) / ybinwidth) + 3) * ybinwidth

            xbins = np.arange(xminlim, xmaxlim + xbinwidth, xbinwidth)
            ybins = np.arange(yminlim, ymaxlim + ybinwidth, ybinwidth)

            ax_histx.hist(x, bins=xbins)
            ax_histy.hist(y, bins=ybins, orientation="horizontal")
        except Exception as e:
            raise OutputError("Error occurred while adding all data to plot") from e

    def write_filters_to_log(self, filters_dict, included_interactions):
        """Takes dictionary of filters, formats as string and writes to log file

        Args:
            filters_dict (dict): dictionary of filtering options
        """
        try:
            buff = ["##### PROPERTIES"]
            for k, v in filters_dict["properties"].items():
                if v is not None:
                    v = "%2.3f" % v
                else:
                    v = " [ none ]"
                buff.append("#  % 7s : %s" % (k, v))
            if filters_dict["filter_ligands_flag"]:
                buff.append("#### LIGAND FILTERS")
                for k, v in filters_dict["ligand_filters"].items():
                    if v is not None:
                        if isinstance(v, list):
                            v = ", ".join([f for f in v if f != ""])
                    else:
                        v = " [ none ]"
                    buff.append("#  % 7s : %s" % (k, v))
            buff.append("#### INTERACTIONS")
            labels = ["~", ""]
            for _type, info in filters_dict["interactions"].items():
                kept_interactions = []
                if len(info) == 0:
                    buff.append("#  % 7s :  [ none ]" % (_type))
                    continue
                for interact in info:
                    if _type + "-" + interact[0] not in included_interactions:
                        continue
                    else:
                        kept_interactions.append(interact)
                res_str = ", ".join(
                    ["(%s)%s" % (labels[int(x[1])], x[0]) for x in kept_interactions]
                )
                l_str = "#  % 7s : %s" % (_type, res_str)
                buff.append(l_str)

            for line in buff:
                self._write_log_line(line)

        except Exception as e:
            raise OutputError("Error occurred while writing filters to log") from e
