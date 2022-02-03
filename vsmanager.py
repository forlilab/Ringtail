from dbmanager import DBManager, DBManagerSQLite
from resultsmanager import ResultsManager

class VSManager():
    """ DOCUMENTATION GOES HERE """
    def __init__(self, db_opts, rman_opts, filters, out_opts, filter_fname=None):
        self.dbman = DBManagerSQLite(db_opts)
        self.results_man = ResultsManager(mode=rman_opts['mode'], dbman = self.dbman, chunk_size=rman_opts['chunk_size'], filelist=rman_opts['filelist'], numclusters=rman_opts['num_clusters'])
        self.filters = filters
        self.out_opts = out_opts
        self.eworst = self.filters['properties']['eworst'] # has default -3 kcal/mol
        self.filter_file = filter_fname
        self.output_manager = Outputter(self, self.out_opts['log'], self.out_opts['plot'])

        #if requested, write database
        if self.dbman.write_db_flag:
            print("adding results")
            self.add_results()

    def add_results(self):
        """"""
        self.results_man.process_results()

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

        #ask DBManager to fetch results
        print("filtering results")
        self.filtered_results = self.dbman.filter_results(self.results_filters_list, self.out_opts['outfields'])
        for line in self.filtered_results:
            self.output_manager.write_log_line(line)
            self.output_manager.plot_single_point(line[0],line[1],"blue") #energy (line[0]) on x axis, le (line[1]) on y axis

        #perform ligand filtering if requested
        if not self.filters['filter_ligands_flag']:
            return
        ligand_filters = self.filters["ligand_filters"]
        print("filtering ligands")
        self.filtered_ligands = self.dbman.filter_ligands(ligand_filters)
        self.output_manager.write_filtered_ligand_header()
        for line in self.filtered_ligands:
            self.output_manager.write_log_line(line)
            self.output_manager.plot_single_point(line[0],line[1],"red") #energy (line[0]) on x axis, le (line[1]) on y axis

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
        'HA',
        'HD',
        'R']

        interaction_filters = self.filters['interactions']
        for key in interaction_filters:
            if interaction_filters[key] != None:
                filters_list.append((key, interaction_filters[key]))

        #get interaction count filters
        interact_count_filters = self.filters["interactions_count"]
        for count in interact_count_filters:
            filters_list.append(count)

        self.results_filters_list = filters_list

    def calculate_energy_percentile(self, percent):
        self.energy_percentile = np.percentile(self.energies, percent)

    def calculate_leff_percentile(self, percent):
        self.leff_percentile = np.percentile(self.leffs, percent)

    def close_database(self):
        """Tell database we are done and it can close the connection"""
        self.dbman.close_connection()

###################################################
############ plotting class #######################
###################################################

class Outputter():

    def __init__(self, vsman, log_file, plot_flag):
        self.log = log_file
        self.vsman = vsman
        self.filter_ligands_flag = self.vsman.filters["filter_ligands_flag"]
        self.plot_flag = plot_flag
        self.num_ligands = self.vsman.results_man.num_result_files
        if self.filter_ligands_flag:
            self.passing_ligand = vsman.filtered_ligands

        if self.vsman.filter_file != None:
            self.fig_base_name = self.vsman.filter_file.split(".")[0]
        else:
            self.fig_base_name = "all_ligands"

        self._create_log_file()
        if self.plot_flag:
            self._initialize_scatter_plot()

    def _initialize_scatter_plot(self):
        self.scatter_plot = plt.subplot(1, 1, 1)
        self.scatter_plot.scatter([],[])
        plot_data = self.vsman.dbman.get_plot()

        for line in plot_data:
            self.plot_single_point(line[0], line[1]) #energy (line[0]) on x axis, le (line[1]) on y axis

    def plot_single_point(self,x,y,color="black"):
        if not self.plot_flag: #make sure user wants plot
            return
        self.scatter_plot.plot(x,y,color)

    def save_scatterplot(self):
        plt.savefig(self.fig_base_name + "_scatter.png")
        plt.close()

    def make_histograms(self):
        self.histFig, (en, le) = plt.subplots(2)
        en.set_title("Histogram of energies of best ligand poses")
        le.set_title("Histogram of ligand efficiencies of best ligand poses")
        en.hist(self.energies, bins = len(self.energies)/150, histtype = "stepfilled")
        le.hist(self.leffs, bins = len(self.leffs)/150, histtype = "stepfilled")

        #check if there were energy or leff percentile filters, add to plot if so
        if self.filtering.energy_percentile_flag:
            en.axvline(c="red", x = self.filtering.energy_percentile)
        if self.filtering.leff_percentile_flag:
            le.axvline(c="red", x=self.filtering.leff_percentile)

        plt.savefig(self.fig_base_name + "_hist.png")
        plt.close(self.histFig)

    def make_scatterplot(self):
        plot = plt.scatter(self.energies, self.leffs, c=self.colors)

        #check if there were energy or leff percentile filters, add to plot if so
        if self.filtering.energy_percentile_flag:
            plot.axvline(c="orange", x = self.filtering.energy_percentile)
        if self.filtering.leff_percentile_flag:
            plot.axhline(c="orange", y=self.filtering.leff_percentile)

        plt.savefig(self.fig_base_name + "_scatter.png")

    def write_log_line(self, line):
        """write a single row to the log file"""
        with open(self.log, "a") as f:
          f.write(new_line)
          f.write("\n")

    def write_filtered_ligand_header(self):
        with open(self.log, "a") as f:
            f.write("\n")
            f.write("***************\n")
            f.write("\n")
            f.write("Filtered Ligands:\n")
            f.write("-----------------\n")

    def write_log(self):

        passing_results_ligands = {}
        passing_ligand_count = 0

        with open(self.log, 'w') as f:
            f.write("Filtered poses:\n")
            f.write("---------------\n")
            for line in self.passing_results:
                first_el_str = str(line[0])
                if first_el_str.endswith(".pdbqt"): #make sure we are only counting pdbqts
                    passing_results_ligands.add(first_el_str) #will only add unique pdbqts
                f.write(" ".join(map(str,line)))
                f.write("\n")
            f.write("\n")
            f.write("***************\n")
            f.write("\n")
            f.write("Filtered Ligands:\n")
            f.write("-----------------\n")
            if self.filtered_ligands != None:
                for line in self.passing_ligands:
                    passing_ligand_count += 1
                    f.write(line + "\n")
            f.write("***************\n")
            f.write("\n")
            f.write("Unique ligands passing result-based filters:\n")
            f.write("-----------------\n")
            for unique_ligand in passing_results_ligands:
                f.write(unique_ligand + "\n")
            f.write("***************\n")
            f.write("\n")
            f.write("Total Ligands:\n")
            f.write("-----------------\n")
            f.write("Ligands passing results-based filters: ")
            f.write(str(len(passing_results_ligands)) + "\n")
            f.write("Ligands passing ligand-based filters: ")
            f.write(str(passing_ligand_count))

    def _create_log_file(self):
        with open(self.log, 'w') as f:
            f.write("Filtered poses:\n")
            f.write("---------------\n")
            
            






