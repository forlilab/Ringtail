from dbmanager import DBManager, DBManagerSQLite
from resultsmanager import ResultsManager
import matplotlib.pyplot as plt
import numpy as np
from rdkit import Chem
from rdkit.Geometry import Point3D
import re

class VSManager():
    """ DOCUMENTATION GOES HERE """
    def __init__(self, db_opts, rman_opts, filters, out_opts, filter_fname=None):
        self.filters = filters
        self.out_opts = out_opts
        self.eworst = self.filters['properties']['eworst'] # has default -3 kcal/mol
        self.filter_file = filter_fname
        self.no_print_flag = self.out_opts["no_print"]

        self.dbman = DBManagerSQLite(db_opts)
        self.results_man = ResultsManager(mode=rman_opts['mode'], dbman = self.dbman, chunk_size=rman_opts['chunk_size'], filelist=rman_opts['filelist'], numclusters=rman_opts['num_clusters'], no_print_flag = self.out_opts["no_print"])
        self.output_manager = Outputter(self, self.out_opts['log'])

        #if requested, write database or add results to an existing one
        if self.dbman.write_db_flag or db_opts["add_results"]:
            print("Adding results...")
            self.add_results()

    def add_results(self):
        """"""
        self.results_man.process_results()
        return 

    def filter(self):
        """"""
        #check if we need to calculate percentiles. If we need to, do so.
        if self.filters['properties']['epercentile'] is not None:
            self.calculate_energy_percentile(self.filters['properties']['epercentile'])
            if self.energy_percentile < self.eworst:
                self.eworst = self.energy_percentile #only keep the most negative (most stringent) filter
        if self.filters['properties']['leffpercentile'] is not None:
            self.calculate_energy_percentile(self.filters['properties']['leffpercentile'])
            if self.leff_percentile < self.eworst:
                self.eworst = self.leff_percentile #only keep the most negative (most stringent) filter
        self.filters['properties']['eworst'] = self.eworst #reset to new value

        #prepare list of filter values and keys for DBManager
        self.prepare_results_filter_list()

        
        print("Filtering results")
        #make sure we have ligand filter list
        if not self.filters['filter_ligands_flag']:
            self.filters["ligand_filters"] = []
        #ask DBManager to fetch results
        self.filtered_results = self.dbman.filter_results(self.results_filters_list, self.filters["ligand_filters"], self.out_opts['outfields'])
        number_passing_ligands = self.dbman.get_number_passing_ligands()[0]
        self.output_manager.log_num_passing_ligands(number_passing_ligands)
        for line in self.filtered_results:
            if not self.no_print_flag:
                print(line)
            self.output_manager.write_log_line(str(line).replace("()",""))#strip parens from line, which is natively a tuple

    def plot(self):
        print("Creating plot of results")
        #plot as requested
        all_data, passing_data = self.dbman.get_plot_data()
        all_plot_data_binned = {}
        #bin the all_ligands data by 1000ths to make plotting faster
        for line in all_data:
            #add to dictionary as bin of energy and le
            data_bin = (round(line[0], 3), round(line[1],3))
            if data_bin not in all_plot_data_binned:
                all_plot_data_binned[data_bin] = 1
            else:
                all_plot_data_binned[data_bin] += 1
        #plot the data
        self.output_manager.plot_all_data(all_plot_data_binned)
        for line in passing_data:
            self.output_manager.plot_single_point(line[0],line[1],"red") #energy (line[0]) on x axis, le (line[1]) on y axis
        self.output_manager.save_scatterplot()

    def prepare_results_filter_list(self):
        """takes filters dictionary from option parser. Output list of tuples to be inserted into sql call string"""

        filters_list = []

        #get property filters
        properties_keys = ['eworst',
        'ebest',
        'leworst',
        'lebest']

        property_filters = self.filters['properties']
        for key in properties_keys:
            if property_filters[key] != None:
                filters_list.append((key, property_filters[key]))

        #get interaction filters
        interaction_keys = ['V',
        'H',
        'R']

        interaction_filters = self.filters['interactions']
        for key in interaction_filters:
            if interaction_filters[key] != None:
                filters_list.append((key, interaction_filters[key]))

        #get interaction count filters
        interact_count_filters = self.filters["interactions_count"]
        for count in interact_count_filters:
            filters_list.append(count) #already a tuple, don't need to reformat

        #add react_any flag
        filters_list.append(("react_any",self.filters["react_any"]))

        self.results_filters_list = filters_list

    def write_molecule_sdfs(self):
        """have output manager write sdf molecules for passing results"""
        passing_molecule_info = self.dbman.fetch_passing_ligand_output_info()
        for (ligname, smiles, atom_indices, hparents) in passing_molecule_info:
            #create rdkit molecule
            mol = self.output_manager.create_molecule(ligname, smiles)
            #fetch coordinates for passing poses and add to rdkit mol
            passing_coordinates = self.dbman.fetch_passing_pose_coordinates(ligname)
            for pose in passing_coordinates:
                mol = self.output_manager.add_pose_to_mol(mol, pose[0], atom_indices)
            #fetch coordinates for non-passing poses and add to mol
            nonpassing_coordinates = self.dbman.fetch_nonpassing_pose_coordinate(ligname)
            for pose in nonpassing_coordinates:
                mol = self.output_manager.add_pose_to_mol(mol, pose[0], atom_indices)

            #write out molecule
            self.output_manager.write_out_mol(ligname, mol)

    def close_database(self):
        """Tell database we are done and it can close the connection"""
        self.dbman.close_connection()

###################################################
############ plotting class #######################
###################################################

class Outputter():

    def __init__(self, vsman, log_file):
        self.log = log_file
        self.vsman = vsman
        self.filter_ligands_flag = self.vsman.filters["filter_ligands_flag"]
        self.num_ligands = self.vsman.results_man.num_result_files
        if self.filter_ligands_flag:
            self.passing_ligand = vsman.filtered_ligands

        if self.vsman.filter_file != None:
            self.fig_base_name = self.vsman.filter_file.split(".")[0]
        else:
            self.fig_base_name = "all_ligands"

        self._create_log_file()

    def plot_all_data(self, binned_data):
        """takes dictionary of binned data where key is the coordinates of the bin and value is the number of points in that bin. Adds to scatter plot colored by value"""
        #gather data
        energies = []
        leffs = []
        bin_counts = []
        for data_bin in binned_data.keys():
            energies.append(data_bin[0])
            leffs.append(data_bin[1])
            bin_counts.append(binned_data[data_bin])

        # start with a square Figure
        fig = plt.figure()

        gs = fig.add_gridspec(2, 2,  width_ratios=(7, 2), height_ratios=(2, 7),
                      left=0.1, right=0.9, bottom=0.1, top=0.9,
                      wspace=0.05, hspace=0.05)

        self.ax = fig.add_subplot(gs[1, 0])
        ax_histx = fig.add_subplot(gs[0, 0], sharex=self.ax)
        ax_histy = fig.add_subplot(gs[1, 1], sharey=self.ax)
        self.ax.set_xlabel("Best Binding Energy / kcal/mol")
        self.ax.set_ylabel("Best Ligand Efficiency")

        self.scatter_hist(energies, leffs, bin_counts, self.ax, ax_histx, ax_histy)


    def plot_single_point(self,x,y,color="black"):
        #self.scatter_plot.scatter(x,y,color)
        self.ax.scatter([x],[y],c=color)

    def save_scatterplot(self):
        plt.savefig(self.fig_base_name + "_scatter.png", bbox_inches="tight")
        plt.close()

    def write_log_line(self, line):
        """write a single row to the log file"""
        with open(self.log, "a") as f:
          f.write(line)
          f.write("\n")

    def log_num_passing_ligands(self, number_passing_ligands):
        with open(self.log, "a") as f:
            f.write("\n")
            f.write("Number passing ligands: {num} \n".format(num=number_passing_ligands))
            f.write("-----------------\n")

    def create_molecule(self, ligname, smiles):
        """creates rdkit molecule from given ligand information"""

        if smiles == "":
            raise RuntimeError("Need SMILES for {molname}".format(molname=ligname))
        return Chem.MolFromSmiles(smiles)

    def add_pose_to_mol(self, mol, coordinates, index_map):
        """add given coordinates to given molecule as new conformer. Index_map maps order of coordinates to order in smile string used to generate rdkit mol"""
        n_atoms = mol.GetNumAtoms()
        conf = Chem.Conformer(n_atoms)
        #split string from database into list, with each element containing the x,y,and z coordinates for one atom
        atom_coordinates = coordinates.split("],")
        if len(atom_coordinates) != n_atoms: #confirm we have the right number of coordinates
            raise RuntimeError("ERROR! Incorrect number of coordinates! Given {n_coords} coordinates for {n_at} atoms!".format(n_coords = len(atom_coordinates), n_at = n_atoms))
        for i in range(n_atoms):
            pdbqt_index = index_map[i+1] - 1
            x, y, z = [float(coord.replace("[","").replace("]","").replace("'","").replace('"','')) for coord in atom_coordinates.split(",")]
            conf.SetAtomPosition(i, Point3D(x, y, z))
        conf_id = mol.AddConformer(conf)
        mol = Chem.AddHs(mol, addCoords=True)
        return mol

    def write_out_mol(self, ligname, mol):
        """writes out given mol as sdf"""
        filename = self.vsman.out_opts["export_poses_path"] + ligname.replace(".pdbqt", ".sdf")
        with SDWriter(filename) as w:
            w.write(mol)


    def _create_log_file(self):
        with open(self.log, 'w') as f:
            f.write("Filtered poses:\n")
            f.write("***************\n")

    def scatter_hist(self, x, y, z, ax, ax_histx, ax_histy):
        # no labels
        ax_histx.tick_params(axis="x", labelbottom=False)
        ax_histy.tick_params(axis="y", labelleft=False)

        # the scatter plot:
        ax.scatter(x, y, c=z, cmap = "Blues")

        # now determine nice limits by hand:
        xbinwidth = 0.25
        ybinwidth = 0.01
        xminlim = (int(min(x)/xbinwidth) + 3) * xbinwidth
        xmaxlim = (int(max(x)/xbinwidth) + 3) * xbinwidth
        yminlim = (int(min(y)/ybinwidth) + 3) * ybinwidth
        ymaxlim = (int(max(y)/ybinwidth) + 3) * ybinwidth

        xbins = np.arange(xminlim, xmaxlim + xbinwidth, xbinwidth)
        ybins = np.arange(yminlim, ymaxlim + ybinwidth, ybinwidth)

        ax_histx.hist(x, bins=xbins)
        ax_histy.hist(y, bins=ybins, orientation='horizontal')
