from dbmanager import DBManagerSQLite
from resultsmanager import ResultsManager
import matplotlib.pyplot as plt
import numpy as np
import json
import warnings
from meeko import RDKitMolCreate
from rdkit import Chem
from rdkit.Chem import SDWriter
import itertools


class VSManager():
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

    def __init__(self,
                 db_opts,
                 rman_opts,
                 filters,
                 out_opts,
                 filter_fname=None):
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
            filter_fname (None/string, optional): Name of file
                containing filters provided by user. None by default.
        """
        self.filters = filters
        self.out_opts = out_opts
        self.eworst = self.filters['properties'][
            'eworst']  # has default -3 kcal/mol
        self.filter_file = filter_fname
        self.no_print_flag = self.out_opts["no_print"]

        self.dbman = DBManagerSQLite(db_opts)
        self.results_man = ResultsManager(
            mode=rman_opts['mode'],
            dbman=self.dbman,
            chunk_size=rman_opts['chunk_size'],
            filelist=rman_opts['filelist'],
            numclusters=rman_opts['num_clusters'],
            no_print_flag=self.no_print_flag,
            single_receptor=out_opts['single_receptor'])
        self.output_manager = Outputter(self, self.out_opts['log'])

        # if requested, write database or add results to an existing one
        if self.dbman.write_db_flag or db_opts["add_results"]:
            print("Adding results...")
            self.add_results()

    def add_results(self):
        """
        Call results manager to process result files and add to database
        """
        self.results_man.process_results()

    def filter(self):
        """
        Prepare list of filters, then hand it off to DBManager to
            perform filtering. Create log of passing results.
        """

        print("Filtering results")
        # get possible permutations of interaction with max_miss excluded
        interaction_combs = self._generate_interaction_combinations(self.filters["max_miss"])

        for ic_idx, combination in enumerate(interaction_combs):
            print(combination)
            # prepare list of filter values and keys for DBManager
            results_filters_list = self.prepare_results_filter_list(combination)

            # make sure we have ligand filter list
            if not self.filters['filter_ligands_flag']:
                self.filters["ligand_filters"] = []
            # set DBMan's internal ic_counter to reflect current ic_idx
            if len(interaction_combs) > 1:
                self.dbman.set_view_suffix(str(ic_idx))
            # ask DBManager to fetch results
            self.filtered_results = self.dbman.filter_results(
                results_filters_list, self.filters["ligand_filters"],
                self.out_opts['outfields'])
            number_passing_ligands = self.dbman.get_number_passing_ligands()
            self.output_manager.write_filters_to_log(self.filters, combination)
            self.output_manager.log_num_passing_ligands(number_passing_ligands)
            self.write_log(self.filtered_results)

    def get_previous_filter_data(self):
        """Get data requested in self.out_opts['outfields'] from the
        results view of a previous filtering
        """
        new_data = self.dbman.fetch_data_for_passing_results(
            self.out_opts['outfields'])
        self.write_log(new_data)

    def plot(self):
        """
        Get data needed for creating Ligand Efficiency vs
        Energy scatter plot from DBManager. Call Outputter to create plot.
        """
        print("Creating plot of results")
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

    def write_log(self, lines):
        """Writes lines from results cursor into log file

        Args:
            lines (DB cursor): Iterable cursor with tuples of data for
                writing into log
        """

        for line in lines:
            if not self.no_print_flag:
                print(line)
            self.output_manager.write_log_line(
                str(line).replace("(", "").replace(
                    ")",
                    ""))  # strip parens from line, which is natively a tuple
        self.output_manager.write_log_line("***************\n")

    def prepare_results_filter_list(self, included_interactions):
        """takes filters dictionary from option parser.
        Output list of tuples to be inserted into sql call string

        Args:
            included_interactions (tuple): Tuple of interactions to include in filter
        """

        filters_list = []

        # get property filters
        properties_keys = [
            'eworst', 'ebest', 'leworst', 'lebest', 'epercentile',
            'leffpercentile'
        ]

        property_filters = self.filters['properties']
        for key in properties_keys:
            if property_filters[key] is not None:
                filters_list.append((key, property_filters[key]))

        interaction_filters = self.filters['interactions']
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
            filters_list.append(
                count)  # already a tuple, don't need to reformat

        # add react_any flag
        filters_list.append(("react_any", self.filters["react_any"]))

        return filters_list

    def write_molecule_sdfs(self):
        """have output manager write sdf molecules for passing results
        """

        if not self.dbman.check_passing_view_exists():
            warnings.warn(
                "Passing results view does not exist in database. Cannot write passing molecule SDFs"
            )
            return
        passing_molecule_info = self.dbman.fetch_passing_ligand_output_info()
        for (ligname, smiles, atom_indices,
             h_parent_line) in passing_molecule_info:
            print("Writing " + ligname.split(".")[0] + ".sdf")
            # create rdkit ligand molecule and flexible residue container
            if smiles == '':
                warnings.warn(f"No SMILES found for {ligname}. Cannot create SDF.")
                continue
            mol = Chem.MolFromSmiles(smiles)
            atom_indices = self._db_string_to_list(atom_indices)
            flexres_mols = {}
            saved_coords = []

            # fetch coordinates for passing poses and add to
            # rdkit ligand mol, add flexible residues
            passing_coordinates = self.dbman.fetch_passing_pose_coordinates(
                ligname)
            mol, flexres_mols, saved_coords = self._add_poses(atom_indices, passing_coordinates,
                                                              mol, flexres_mols, saved_coords)

            # fetch coordinates for non-passing poses
            # and add to ligand mol, flexible residue mols
            nonpassing_coordinates = self.dbman.fetch_nonpassing_pose_coordinates(
                ligname)
            mol, flexres_mols, saved_coords = self._add_poses(atom_indices, nonpassing_coordinates,
                                                              mol, flexres_mols, saved_coords)

            # write out molecule
            # h_parents = self._format_h_parents(h_parent_line)
            h_parents = [int(idx) for idx in self._db_string_to_list(h_parent_line)]
            self.output_manager.write_out_mol(ligname, mol, flexres_mols, saved_coords, h_parents)

    def export_csv(self, requested_data, csv_name):
        """Get requested data from database, export as CSV

        Args:
            requested_data (string): Table name or SQL-formatted query
            csv_name (string): Name for exported CSV file
        """
        df = self.dbman.fetch_dataframe_from_db(requested_data)
        df.to_csv(csv_name)

    def close_database(self):
        """Tell database we are done and it can close the connection
        """
        self.dbman.close_connection()

    def _clean_db_string(self, input_str):
        """take a db string representing a list,
        strips unwanted characters

        Args:
            input_str (str)

        Returns:
            String: cleaned string
        """
        return input_str.replace("[", "").replace("]", "").replace(
            "'", "").replace('"',
                             '').replace("\n",
                                         "").replace("\\n",
                                                     "").replace(" \n", "")

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

        return "\n".join([
            line.lstrip(" ") for line in list(
                filter(None, [self._clean_db_string(line) for line in pdbqt_lines]))
        ])

    def _add_poses(self, atom_indices, poses, mol, flexres_mols, saved_coords):
        """Add poses from given cursor to rdkit mols for ligand and flexible residues

        Args:
            atom_indices (List): List of ints indicating mapping of coordinate indices to smiles indices
            poses (iterable): iterable containing ligand_pose, flexres_pose, flexres_names
            mol (RDKit Mol): RDKit molecule for ligand
            flexres_mols (Dict): Dictionary of rdkit molecules for flexible residues
            saved_coords (list): list of coordinates to save for adding hydrogens later
        """
        for ligand_pose, flexres_pose, flexres_names in poses:
            ligand_pose = self._db_string_to_list(ligand_pose)
            flexres_pose = self._db_string_to_list(flexres_pose)
            flexres_names = [name for idx, name in enumerate(self._db_string_to_list(flexres_names))]
            flexres_pdbqts = [self._generate_pdbqt_block(res) for res in flexres_pose]
            mol, flexres_mols = RDKitMolCreate.add_pose_to_mol(mol, ligand_pose, atom_indices,
                                                               flexres_mols=flexres_mols,
                                                               flexres_poses=flexres_pdbqts,
                                                               flexres_names=flexres_names)
            saved_coords.append(ligand_pose)
        return mol, flexres_mols, saved_coords

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
            warnings.warn("Requested max_miss options greater than number of interaction filters given. Defaulting to max_miss = number interaction filters")
            max_miss = len(all_interactions)

        # BASE CASE:
        if max_miss == 0:
            return [tuple(all_interactions)]
        else:
            combinations = list(itertools.combinations(all_interactions, len(all_interactions) - max_miss))
            return combinations + self._generate_interaction_combinations(max_miss=max_miss - 1)

# # # # # # # # # # # # # #
# # # Output class # # #
# # # # # # # # # # # # # #

class Outputter():
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

    def __init__(self, vsman, log_file):
        """Initialize Outputter object and create log file

        Args:
            vsman (VSManager): VSManager object that created outputter
            log_file (string): name for log file
        """

        self.log = log_file
        self.vsman = vsman

        if self.vsman.filter_file is not None:
            self.fig_base_name = self.vsman.filter_file.split(".")[0]
        else:
            self.fig_base_name = "all_ligands"

        self._create_log_file()

    def plot_all_data(self, binned_data):
        """takes dictionary of binned data where key is the
        coordinates of the bin and value is the number of points in that bin.
        Adds to scatter plot colored by value

        Args:
            binned_data (dict): Keys are tuples of key and y value for bin.
                Value is the count of points falling into that bin.
        """
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

        gs = fig.add_gridspec(2,
                              2,
                              width_ratios=(7, 2),
                              height_ratios=(2, 7),
                              left=0.1,
                              right=0.9,
                              bottom=0.1,
                              top=0.9,
                              wspace=0.05,
                              hspace=0.05)

        self.ax = fig.add_subplot(gs[1, 0])
        ax_histx = fig.add_subplot(gs[0, 0], sharex=self.ax)
        ax_histy = fig.add_subplot(gs[1, 1], sharey=self.ax)
        self.ax.set_xlabel("Best Binding Energy / kcal/mol")
        self.ax.set_ylabel("Best Ligand Efficiency")

        self.scatter_hist(energies, leffs, bin_counts, self.ax, ax_histx,
                          ax_histy)

    def plot_single_point(self, x, y, color="black"):
        """Add point to scatter plot with given x and y coordinates and color.

        Args:
            x (float): x coordinate
            y (float): y coordinate
            color (str, optional): Color for point. Default black.
        """
        self.ax.scatter([x], [y], c=color)

    def save_scatterplot(self):
        """
        Saves current figure as [self.fig_base_name]_scatter.png
        """
        plt.savefig(self.fig_base_name + "_scatter.png", bbox_inches="tight")
        plt.close()

    def write_log_line(self, line):
        """write a single row to the log file

        Args:
            line (string): Line to write to log
        """
        with open(self.log, "a") as f:
            f.write(line)
            f.write("\n")

    def log_num_passing_ligands(self, number_passing_ligands):
        """
        Write the number of ligands which pass given filter to log file

        Args:
            number_passing_ligands (int): number of ligands that passed filter
        """
        with open(self.log, "a") as f:
            f.write("\n")
            f.write("Number passing ligands: {num} \n".format(
                num=str(number_passing_ligands)))
            f.write("---------------\n")

    def write_out_mol(self, ligname, mol, flexres_mols, saved_coords, h_parents):
        """writes out given mol as sdf

        Args:
            ligname (string): name of ligand that will be used to
                name output SDF file
            mol (meeko PDBQTMolecule object): Meeko PDBQTMolecule object to be written to SDF
        """
        filename = self.vsman.out_opts["export_poses_path"] + ligname + ".sdf"
        mol = RDKitMolCreate.export_combined_rdkit_mol(mol, flexres_mols, saved_coords, h_parents)
        with SDWriter(filename) as w:
            for conf in mol.GetConformers():
                w.write(mol, conf.GetId())

    def _create_log_file(self):
        """
        Initializes log file
        """
        with open(self.log, 'w') as f:
            f.write("Filtered poses:\n")
            f.write("***************\n")

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
        ax_histy.hist(y, bins=ybins, orientation='horizontal')

    def write_filters_to_log(self, filters_dict, included_interactions):
        """Takes dictionary of filters, formats as string and writes to log file

        Args:
            filters_dict (dict): dictionary of filtering options
        """

        buff = ['##### PROPERTIES']
        for k, v in filters_dict["properties"].items():
            if v is not None:
                v = "%2.3f" % v
            else:
                v = " [ none ]"
            buff.append("#  % 7s : %s" % (k, v))
        if filters_dict['filter_ligands_flag']:
            buff.append("#### LIGAND FILTERS")
            for k, v in filters_dict["ligand_filters"].items():
                if v is not None:
                    v = "%2.3f" % v
                else:
                    v = " [ none ]"
                buff.append("#  % 7s : %s" % (k, v))
        buff.append("#### INTERACTIONS")
        labels = ['-', '+']
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
            res_str = ", ".join(['(%s)%s' % (labels[int(x[1])], x[0]) for x in kept_interactions])
            l_str = "#  % 7s : %s" % (_type, res_str)
            buff.append(l_str)

        for line in buff:
            self.write_log_line(line)
