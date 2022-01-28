from dbmanager import DBManager
from resultsmanager import ResultsManager

class VSManager():
    """ DOCUMENTATION GOES HERE """
    def __init__(self, db_opts, rman_opts, filters, out_opts):
        self.dbman = DBManager(opts=db_opts)
        self.results_man = ResultsManager(mode=rman_opts['mode'], dbman = self.dbman, chunk_size=rman_opts['chunk_size'], file_list=rman_opts['filelist'], numclusters=rman_opts['num_clusters'])
        self.filters = filters
        self.out_opts = out_opts
        self.eworst = self.filters['properties']['eworst'] # has default -3 kcal/mol

        #if requested, write database
        if dbman.write_flag:
            self.add_results()

        if self.filters['properties']['epercentile'] not None or self.filters['properties']['leffpercentile'] not None or self.out_opts['plot']:
            self.top_energies, self.top_leffs, self.top_data = self.dbman.get_top_energies_leffs()

    def add_results(self):
        """"""
        self.results_man.process_results()

    def filter(self):
        """"""
        #check if we need to calculate percentiles. If we need to, do so.
        if self.filters['properties']['epercentile'] not None:
            self.calculate_energy_percentile(self.filters['properties']['epercentile'])
            if self.energy_percentile < self.eworst:
                self.eworst = self.energy_percentile #only keep the most negative (most stringent) filter
        if self.filters['properties']['leffpercentile'] not None:
            self.calculate_energy_percentile(self.filters['properties']['leffpercentile'])
            if self.leff_percentile < self.eworst:
                self.eworst = self.leff_percentile #only keep the most negative (most stringent) filter
        self.filters['properties']['eworst'] = self.eworst #reset to new value

        #prepare list of filter values and keys for DBManager
        self.prepare_results_filter_list()

        #ask DBManager to fetch results
        self.filtered_results = dbman.filter_results(self.results_filters_list, self.out_opts['outfields'])

        #perform ligand filtering if requested
        if self.filters['filter_ligands_flag']:
            ligand_filters = self.filters["ligand_filters"]
            self.filtered_ligands = dbman.filter_ligands(ligand_filters)

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

###################################################
############ plotting class #######################
###################################################

class Outputter():

    def __init__(self, vsman, log_file):
        self.log = log_file
        self.filter_ligands_flag = self.vsman.filters["filter_ligands_flag"]
        self.energies = vsman.top_energies
        self.leffs = vsman.top_leffs
        self.plot_data = vsman.top_data
        self.passing_results = vsman.filtered_results
        self.passing_ligand = vsman.filtered_ligands

        if self.vsman.filter_file != None:
            self.fig_base_name = self.vsman.filter_file.split(".")[0]
        else:
            self.fig_base_name = "all_ligands"

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

    def assign_scatter_colors(self):
        colors = []
        for row in self.plot_data:
            if int(row[3]) == 1:
                ligand_name = row[0]
                if not self.filter_ligands_flag:
                    filtered_ligands = []

                if ligand_name in [row[0] for row in self.filtered_results] and ligand_name in [row[0] for row in self.filtered_ligands]:
                    colors.append("purple")
                elif ligand_name in [row[0] for row in self.filtered_ligands]:
                    colors.append("red")
                elif ligand_name in [row[0] for row in self.filtered_results]:
                    colors.append("blue")
                else:
                    colors.append("black")

        self.colors = colors

    def make_scatterplot(self):
        plot = plt.scatter(self.energies, self.leffs, c=self.colors)

        #check if there were energy or leff percentile filters, add to plot if so
        if self.filtering.energy_percentile_flag:
            plot.axvline(c="orange", x = self.filtering.energy_percentile)
        if self.filtering.leff_percentile_flag:
            plot.axhline(c="orange", y=self.filtering.leff_percentile)

        plt.savefig(self.fig_base_name + "_scatter.png")

    def write_log(self):

        with open(self.log, 'w') as f:
            f.write("Filtered poses:\n")
            f.write("---------------\n")
            for line in self.passing_results:
                f.write(" ".join(map(str,line)))
                f.write("\n")
            f.write("\n")
            f.write("***************\n")
            f.write("\n")
            f.write("Filtered Ligands:\n")
            f.write("-----------------\n")
            if self.filtered_ligands != None:
                for line in self.passing_ligands:
                    f.write(line + "\n")






