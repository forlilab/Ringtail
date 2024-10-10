#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail output manager
#

from .receptormanager import ReceptorManager
from .exceptions import OutputError
from .ringtailoptions import Filters
import os
import json
import numpy as np
import time
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import colors
from rdkit.Chem import SDWriter
from meeko import RDKitMolCreate
from .logutils import LOGGER


class OutputManager:
    """Class for creating outputs, can be a context manager to handle log files

    Attributes:
        log_file (str): name for log file
        export_sdf_path (str): path for exporting SDF molecule files
        _log_open (bool): if log file is open or not
    """

    def __init__(self, log_file=None, export_sdf_path=None):
        self.log_file = log_file
        self.export_sdf_path = export_sdf_path
        self._log_open = False
        self.logger = LOGGER

    def __enter__(self):
        """Opening outputmanager as a context manager"""
        self.open_logfile()

    def __exit__(self, exc_type, exc_value, traceback):
        """Closing outputmanager as a context manager"""
        self.close_logfile()

    # -#-#- Log file methods -#-#-#
    def open_logfile(self, write_filters_header=True):
        """
        Opens log file and creates it if needed

        Args:
            write_filters_header (bool): only used because one method does not take the same headers

        Raises:
            OutputError
        """
        self.log_file = open(
            self.log_file, "w"
        )  # makes log_file attribute a file pointer from the string path name
        self._log_open = True
        try:
            if write_filters_header:
                self.log_file.write("Filters:\n")
                self.log_file.write("***************\n")
        except Exception as e:
            raise OutputError("Error while creating log file") from e

    def close_logfile(self):
        """Closes the log file properly and reset file pointer to filename"""
        if self._log_open:
            self.log_file.close()
            self.log_file = os.path.basename(
                self.log_file.name
            )  # returns log_file attribute from file pointer to string path name
            self._log_open = False

    def write_filter_log(self, lines):
        """Writes lines from results iterable into log file

        Args:
            lines (iterable): Iterable with tuples of data for
                writing into log

        Raises:
            OutputError

        Returns:
            int: number of ligands passing that are written to log file
        """
        try:
            time0 = time.perf_counter()
            num_passing = 0
            for line in lines:
                self._write_log_line(
                    str(line).replace("(", "").replace(")", "")
                )  # strip parens from line, which is natively a tuple
                num_passing += 1
            self._write_log_line("***************\n")
            self.logger.debug(
                f"Time to write log: {time.perf_counter() - time0:.2f} seconds"
            )
            return num_passing
        except Exception as e:
            raise OutputError("Error occurred during log writing") from e

    def _write_log_line(self, line: str):
        """write a single row to the log file

        Args:
            line (str): Line to write to log

        Raises:
            OutputError
        """
        try:
            self.log_file.write(line)
            self.log_file.write("\n")
        except Exception as e:
            raise OutputError(f"Error writing line {line} to log") from e

    def log_num_passing_ligands(self, number_passing_ligands: int):
        """
        Write the number of ligands which pass given filter to log file

        Args:
            number_passing_ligands (int): number of ligands that passed filter

        Raises:
            OutputError
        """
        try:
            self.log_file.write("\n")
            self.log_file.write(
                f"Number passing ligands: {str(number_passing_ligands)} \n"
            )
            self.log_file.write("---------------\n")
        except Exception as e:
            raise OutputError("Error writing number of passing ligands in log") from e

    def write_results_bookmark_to_log(self, bookmark_name):
        """Write the name of the result bookmark into log

        Args:
            bookmark_name (str): name of current results' bookmark in db

        Raises:
            OutputError
        """
        try:
            self.log_file.write("\n")
            self.log_file.write(f"Result bookmark name: {bookmark_name}\n")
            self.log_file.write("***************\n")
        except Exception as e:
            raise OutputError("Error writing bookmark name to log") from e

    def write_filters_to_log(
        self, filters_dict, included_interactions, additional_info=""
    ):
        """Takes dictionary of filters, formats as string and writes to log file

        Args:
            filters_dict (dict): dictionary with filtering options
            included_interactions (list): types of interactions to include in the filtering
            additional_info (str): any additional information to write to top of log file

        Raises:
            OutputError
        """
        try:
            buff = [additional_info, "##### PROPERTIES"]
            for k in Filters.get_filter_keys("property"):
                v = filters_dict.pop(k)
                if v is not None:
                    v = "%2.3f" % v
                else:
                    v = " [ none ]"
                buff.append("#  % 7s : %s" % (k, v))
            buff.append("#### LIGAND FILTERS")
            for k in Filters.get_filter_keys("ligand"):
                v = filters_dict.pop(k)
                if v is not None:
                    if isinstance(v, list):
                        v = ", ".join([str(f) for f in v if f != ""])
                else:
                    v = " [ none ]"
                buff.append("#  % 7s : %s" % (k, v))
            buff.append("#### INTERACTIONS")
            labels = ["~", ""]
            for _type in Filters.get_filter_keys("interaction"):
                info = filters_dict.pop(_type)
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

            buff.append("#### OTHER FILTERS")
            for k, v in filters_dict.items():
                if v is None:
                    v = " [ none ]"
                buff.append("#  % 7s : %s" % (k, v))

            for line in buff:
                self._write_log_line(line)

        except Exception as e:
            raise OutputError("Error occurred while writing filters to log") from e

    def write_maxmiss_union_header(self):
        """
        Properly formats header for the log file if using max_miss and enumerate_interaction_combs
        """
        self.log_file.write("\n---------------\n")
        self.log_file.write("Max Miss Union:\n")

    def write_find_similar_header(self, query_ligname, cluster_name):
        """
        Properly formats header for the log file find_similar_ligands
        """
        if not self._log_open:
            self.open_logfile(write_filters_header=False)
        self.log_file.write("\n---------------\n")
        self.log_file.write(
            f"Found ligands similar to {query_ligname} in clustering {cluster_name}:\n"
        )

    # -#-#- Non-logfile methods -#-#-#

    def write_out_mol(self, filename, mol, flexres_mols, properties):
        """Writes out given mol as sdf. Will create the specified sdf folder in
        current working directory if needed.

        Args:
            filename (str): name of SDF file that will be written to
            mol (RDKit.Chem.Mol): RDKit molobject to be written to SDF
            flexres_mols (list): dictionary of rdkit molecules for flexible residues
            properties (dict): dictionary of list of properties to add to mol before writing

        Raises:
            OutputError
        """
        if (
            self.export_sdf_path is not None
            and not self.export_sdf_path == ""
            and not os.path.isdir(self.export_sdf_path)
        ):
            os.makedirs(self.export_sdf_path)
            self.logger.info(
                "Specified directory for SDF files was created in current working directory."
            )
        filename = self.export_sdf_path + "/" + filename
        try:
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

            # open/create file so it can be appended to if requested
            with open(filename, "a") as sdf_file:
                w = SDWriter(sdf_file)
                for conf in mol.GetConformers():
                    w.write(mol, conf.GetId())
                w.close()

        except Exception as e:
            raise OutputError("Error occurred while writing SDF from RDKit Mol") from e

    def write_receptor_pdbqt(self, recname: str, receptor_compbytes):
        """
        Writes a pdbqt file from receptor "blob"

        Args:
            recname (str): name of receptor to use in output filename
            receptor_compbytes (blob): receptor blob
        """

        if not recname.endswith(".pdbqt"):
            recname = recname + ".pdbqt"
        receptor_str = ReceptorManager.blob2str(receptor_compbytes)
        with open(recname, "w") as f:
            f.write(receptor_str)

    def scatter_hist(self, x, y, z, ax, ax_histx, ax_histy):
        """
        Makes scatterplot with a histogram on each axis

        Args:
            x (list): x coordinates for data
            y (list): y coordinates for data
            z (list): z coordinates for data
            ax (matplotlib.axis): scatterplot axis
            ax_histx (matplotlib.axis): x histogram axis
            ax_histy (matplotlib.axis): y histogram axis

        Raises:
            OutputError
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

    def plot_all_data(self, binned_data):
        """Takes dictionary of binned data where key is the
        coordinates of the bin and value is the number of points in that bin.
        Adds to scatter plot colored by value

        Args:
            binned_data (dict): Keys are tuples of key and y value for bin.
                Value is the count of points falling into that bin.

        Raises:
            OutputError
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
                    colors.Normalize(vmin=min(bin_counts), vmax=max(bin_counts)),
                ),
                cax=self.ax.inset_axes([0.85, 0.1, 0.05, 0.8]),
                label="Scatterplot bin count",
            )
            self.ax.set_xlabel("Best docking score / kcal/mol")
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

        Raises:
            OutputError
        """
        try:
            self.ax.scatter([x], [y], c=color)
        except Exception as e:
            raise OutputError("Error occurred while plotting") from e

    def save_scatterplot(self):
        """
        Saves current figure as scatter.png

        Raises:
            OutputError
        """
        try:
            plt.savefig("scatter.png", bbox_inches="tight")
            plt.close()
        except Exception as e:
            raise OutputError("Error while saving figure") from e
