#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail virtual screening manager
#

import matplotlib.pyplot as plt
import numpy as np
import json
import warnings
from meeko import RDKitMolCreate
from ringtail import DBManagerSQLite
from ringtail import ResultsManager
from .exceptions import (
    DatabaseConnectionError,
    DatabaseTableCreationError,
    DatabaseError,
)
from .exceptions import VirtualScreeningError, ResultsProcessingError, OutputError
from rdkit import Chem
from rdkit.Chem import SDWriter
import itertools
import logging


class VSManager:
    """Manager for coordinating different actions on virtual screening
    i.e. adding results to db, filtering, output options

    Attributes:
        dbman (DBManager): Interface module with database
        eworst (float): The worst scoring energy filter value requested by user
        filter_file (string): Name of file containing filters provided by user
        filtered_results (DB cursor object): Cursor object
            containing results passing requested filters (iterable)
        filters (dictionary): Dictionary containing user-specified filters
        no_print_flag (boolean): Flag specifying whether passing results
            should be printed to terminal
        out_opts (dictionary): Specified output options including data fields
            to output, export_poses_path, log file name
        output_manager (Outputter object): Manager for output tasks of
            log-writting, plotting, ligand SDF writing
        results_filters_list (List): List of tuples of filter option and value
        results_man (ResultsManager object): Manager for processing result
            files for insertion into database
    """

    def __init__(self, db_opts, rman_opts, filters, out_opts):
        """Initialize VSManager object. Will create DBManager object to serve
        as interface with database (currently implemented in SQLite).
        Will create ResultsManager to process result files.
        Will create Outputter object to assist in creating output files.

        Args:
            db_opts (dictionary): dictionary of options required by DBManager
            rman_opts (dictionary): dictionary of options required by
                results manager
            filters (dictionary): Dictionary containing user-specified filters
            out_opts (dictionary): Specified output options including data
                fields to output, export_poses_path, log file name
        """
        self.filters = filters
        self.out_opts = out_opts
        self.rman_opts = rman_opts
        self.db_opts = db_opts

    def __enter__(self):
        try:
            self.dbman = DBManagerSQLite(self.db_opts["dbFile"], self.db_opts)
        except DatabaseConnectionError as e:
            raise VirtualScreeningError(
                "Error encountered while connecting to database. Please ensure that given database file name is correct."
            ) from e
        except DatabaseTableCreationError as e:
            raise VirtualScreeningError(
                "Error encountered while creating database tables. If database already exists, use --add_results or --overwrite."
            ) from e
        except DatabaseError as e:
            raise VirtualScreeningError(
                "Error occurred while initializing database."
            ) from e

        # if requested, write database or add results to an existing one
        if self.dbman.write_db_flag or self.db_opts["add_results"]:
            logging.info("Adding results...")
            try:
                self.results_man = ResultsManager(opts=self.rman_opts, dbman=self.dbman)
                self.add_results()
            except ResultsProcessingError as e:
                raise VirtualScreeningError("Error occured while adding results") from e

        else:
            try:
                self.output_manager = Outputter(
                    self.out_opts["log"], self.out_opts["export_poses_path"]
                )
            except OutputError as e:
                raise VirtualScreeningError(
                    "Error occured while creating output manager"
                ) from e

        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close_database()

    def add_results(self):
        """
        Call results manager to process result files and add to database
        """
        self.results_man.process_results()

    def add_receptors_to_db(self, receptors):
        """Add receptor to database

        Args:
            receptors (list): list of receptor blobs to add to database
        """
        for rec, rec_name in receptors:
            # NOTE: in current implementation, only one receptor allowed per database
            # Check that any receptor row is incomplete (needs receptor blob) before inserting
            try:
                filled_receptor_rows = self.dbman.get_number_filled_receptor_rows()
                if filled_receptor_rows != 0:
                    raise VirtualScreeningError(
                        "Expected Receptors table to have no receptor objects present, already has {0} receptor present. Cannot add more than 1 receptor to a database.".format(
                            filled_receptor_rows
                        )
                    )
                self.dbman.add_receptor_object_to_row(rec)
            except DatabaseError as e:
                raise VirtualScreeningError(
                    "Error occurred while adding receptor to database"
                ) from e

    def filter(self):
        """
        Prepare list of filters, then hand it off to DBManager to
            perform filtering. Create log of passing results.
        """

        logging.info("Filtering results")
        # get possible permutations of interaction with max_miss excluded
        interaction_combs = self._generate_interaction_combinations(
            self.filters["max_miss"]
        )

        for ic_idx, combination in enumerate(interaction_combs):
            # prepare list of filter values and keys for DBManager
            results_filters_list = self.prepare_results_filter_list(combination)

            # make sure we have ligand filter list
            if not self.filters["filter_ligands_flag"]:
                self.filters["ligand_filters"] = []
            # set DBMan's internal ic_counter to reflect current ic_idx
            if len(interaction_combs) > 1:
                self.dbman.set_view_suffix(str(ic_idx))
            # ask DBManager to fetch results
            try:
                self.filtered_results = self.dbman.filter_results(
                    results_filters_list,
                    self.filters["ligand_filters"],
                    self.out_opts["outfields"],
                )
                number_passing_ligands = self.dbman.get_number_passing_ligands()
                result_bookmark_name = self.dbman.get_current_view_name()
                self.output_manager.write_filters_to_log(self.filters, combination)
                self.output_manager.write_results_bookmark_to_log(result_bookmark_name)
                self.output_manager.log_num_passing_ligands(number_passing_ligands)
                self.output_manager.write_log(self.filtered_results)
            except DatabaseError as e:
                raise VirtualScreeningError(
                    "Database error occurred while filtering"
                ) from e
            except OutputError as e:
                raise VirtualScreeningError(
                    "Logging error occurred after filtering"
                ) from e
            except Exception as e:
                raise VirtualScreeningError("Error occurred while filtering") from e

    def get_previous_filter_data(self):
        """Get data requested in self.out_opts['outfields'] from the
        results view of a previous filtering
        """
        try:
            new_data = self.dbman.fetch_data_for_passing_results(
                self.out_opts["outfields"]
            )
            self.output_manager.write_log(new_data)
        except DatabaseError as e:
            raise VirtualScreeningError(
                "Database error occurred while fetching data for bookmark"
            ) from e
        except OutputError as e:
            raise VirtualScreeningError("Error occurred while writing log") from e
        except Exception as e:
            raise VirtualScreeningError(
                "Error occurred while fetching data for bookmark"
            ) from e

    def plot(self):
        """
        Get data needed for creating Ligand Efficiency vs
        Energy scatter plot from DBManager. Call Outputter to create plot.
        """
        try:
            logging.info("Creating plot of results")
            # get data from DBMan
            all_data, passing_data = self.dbman.get_plot_data()
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
        except DatabaseError as e:
            raise VirtualScreeningError(
                "Database error occurred while fetching data for plot"
            ) from e
        except Exception as e:
            raise VirtualScreeningError("Error occurred during plotting") from e

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
        """have output manager write sdf molecules for passing results

        Args:
            write_nonpassing (bool, optional): Option to include non-passing poses for passing ligands
        """
        try:
            if not self.dbman.check_passing_view_exists():
                warnings.warn(
                    "Passing results view does not exist in database. Cannot write passing molecule SDFs"
                )
                return
            passing_molecule_info = self.dbman.fetch_passing_ligand_output_info()
            for (ligname, smiles, atom_indices, h_parent_line) in passing_molecule_info:
                logging.info("Writing " + ligname.split(".")[0] + ".sdf")
                # create rdkit ligand molecule and flexible residue container
                if smiles == "":
                    warnings.warn(f"No SMILES found for {ligname}. Cannot create SDF.")
                    continue
                mol = Chem.MolFromSmiles(smiles)
                atom_indices = self._db_string_to_list(atom_indices)
                flexres_mols = {}
                saved_coords = []

                # fetch coordinates for passing poses and add to
                # rdkit ligand mol, add flexible residues
                properties = {
                    "Binding energies": [],
                    "Ligand effiencies": [],
                    "Interactions": [],
                }
                passing_properties = self.dbman.fetch_passing_pose_properties(ligname)
                mol, flexres_mols, saved_coords, properties = self._add_poses(
                    atom_indices,
                    passing_properties,
                    mol,
                    flexres_mols,
                    saved_coords,
                    properties,
                )

                # fetch coordinates for non-passing poses
                # and add to ligand mol, flexible residue mols
                if write_nonpassing:
                    nonpassing_properties = self.dbman.fetch_nonpassing_pose_properties(
                        ligname
                    )
                    mol, flexres_mols, saved_coords, properties = self._add_poses(
                        atom_indices,
                        nonpassing_properties,
                        mol,
                        flexres_mols,
                        saved_coords,
                        properties,
                    )

                # write out molecule
                # h_parents = self._format_h_parents(h_parent_line)
                h_parents = [int(idx) for idx in self._db_string_to_list(h_parent_line)]
                self.output_manager.write_out_mol(
                    ligname, mol, flexres_mols, saved_coords, h_parents, properties
                )

        except DatabaseError as e:
            raise VirtualScreeningError(
                "Error occurred while fetching database information for SDF output"
            ) from e
        except OutputError as e:
            raise VirtualScreeningError(
                f"Error occured while writing {ligname} to an SDF"
            ) from e
        except Exception as e:
            raise VirtualScreeningError("Error occurred during SDF output") from e

    def export_csv(self, requested_data, csv_name, table=False):
        """Get requested data from database, export as CSV

        Args:
            requested_data (string): Table name or SQL-formatted query
            csv_name (string): Name for exported CSV file
            table (bool): flag indicating is requested data is a table name
        """
        try:
            df = self.dbman.to_dataframe(requested_data, table=table)
            df.to_csv(csv_name)
        except DatabaseError as e:
            raise VirtualScreeningError(
                f"Error occured while getting data for exporting CSV of {requested_data}"
            ) from e
        except Exception as e:
            raise VirtualScreeningError(
                f"Error occured while exporting CSV of {requested_data}"
            ) from e

    def close_database(self):
        """Tell database we are done and it can close the connection"""
        self.dbman.close_connection()

    def _clean_db_string(self, input_str):
        """take a db string representing a list,
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

    def _db_string_to_list(self, input_str):
        """Convert string form of list from database to list

        Args:
            input_str (TYPE): Description
        """

        return json.loads(input_str)

    def _generate_pdbqt_block(self, pdbqt_lines):
        """Generate pdbqt block from given lines from a pdbqt

        Args:
            flexres_lines (TYPE): Description
        """

        return "\n".join(
            [
                line.lstrip(" ")
                for line in list(
                    filter(None, [self._clean_db_string(line) for line in pdbqt_lines])
                )
            ]
        )

    def _add_poses(
        self, atom_indices, poses, mol, flexres_mols, saved_coords, properties
    ):
        """Add poses from given cursor to rdkit mols for ligand and flexible residues

        Args:
            atom_indices (List): List of ints indicating mapping of coordinate indices to smiles indices
            poses (iterable): iterable containing ligand_pose, flexres_pose, flexres_names
            mol (RDKit Mol): RDKit molecule for ligand
            flexres_mols (Dict): Dictionary of rdkit molecules for flexible residues
            saved_coords (list): list of coordinates to save for adding hydrogens later
            properties (dict): Dictionary of lists of properties, with each element corresponding to that conformer in the rdkit mol

        """
        for (
            Pose_ID,
            energies_binding,
            leff,
            ligand_pose,
            flexres_pose,
            flexres_names,
        ) in poses:
            # fetch info about pose interactions and format into string with format <type>-<chain>:<resname>:<resnum>:<atomname>:<atomnumber>, joined by commas
            pose_bitvector = self.dbman.fetch_interaction_bitvector(Pose_ID)
            if pose_bitvector is not None:
                interaction_indices = []
                interactions_list = []
                for idx, bit in enumerate(pose_bitvector):
                    if bit == 1:
                        interaction_indices.append(
                            idx + 1
                        )  # adjust for indexing starting at 1
                for int_idx in interaction_indices:
                    interaction_info = self.dbman.fetch_interaction_info_by_index(int_idx)
                    interaction = interaction_info[0] + "-" + ":".join(interaction_info[1:])
                    interactions_list.append(interaction)
                interactions_str = ", ".join(interactions_list)
                properties["Interactions"].append(interactions_str)
            # add properties to dictionary lists
            properties["Binding energies"].append(energies_binding)
            properties["Ligand effiencies"].append(leff)
            # get pose coordinate info
            ligand_pose = self._db_string_to_list(ligand_pose)
            flexres_pose = self._db_string_to_list(flexres_pose)
            flexres_names = [
                name for idx, name in enumerate(self._db_string_to_list(flexres_names))
            ]
            flexres_pdbqts = [self._generate_pdbqt_block(res) for res in flexres_pose]
            mol, flexres_mols = RDKitMolCreate.add_pose_to_mol(
                mol,
                ligand_pose,
                atom_indices,
                flexres_mols=flexres_mols,
                flexres_poses=flexres_pdbqts,
                flexres_names=flexres_names,
            )
            saved_coords.append(ligand_pose)
        return mol, flexres_mols, saved_coords, properties

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
            warnings.warn(
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


class Outputter:
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
        vsman (VSManager): VSManager object that created outputter

    """

    def __init__(self, log_file, export_sdf_path=""):
        """Initialize Outputter object and create log file

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
        Saves current figure as [self.fig_base_name]_scatter.png
        """
        try:
            plt.savefig(self.fig_base_name + "_scatter.png", bbox_inches="tight")
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

    def write_out_mol(
        self, ligname, mol, flexres_mols, saved_coords, h_parents, properties
    ):
        """writes out given mol as sdf

        Args:
            ligname (string): name of ligand that will be used to
                name output SDF file
            mol (RDKit mol object): RDKit molobject to be written to SDF
            flexres_mols (dict): dictionary of rdkit molecules for flexible residues
            saved_coords (list): list of coordinates that have been used already
            h_parents (list): list of atom indices of hydrogens and their parents
            properties (dict): dictionary of list of properties to add to mol before writing
        """
        try:
            filename = self.export_sdf_path + ligname + ".sdf"
            mol = RDKitMolCreate.export_combined_rdkit_mol(
                mol, flexres_mols, saved_coords, h_parents
            )
            # convert properties to strings as needed
            for k, v in properties.items():
                if isinstance(v, list):
                    properties[k] = json.dumps(v)
                elif not isinstance(v, str):
                    properties[k] = str(v)
            RDKitMolCreate.add_properties_to_mol(mol, properties)

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
            ax.scatter(x, y, c=z, cmap="Blues")

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
