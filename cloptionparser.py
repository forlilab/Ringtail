import sys
import argparse
from glob import glob
import os

class CLOptionParser():

    def __init__(self):
        self.write_db_flag = False
        self.description = """description something of something"""
        self.usage = "READ THE MANUAL!"
        self.name= "vs_results_display"
        self.epilog="""

        REQUIRED PACKAGES
                Requires numpy, multiprocessing, matplotlib, sqlite3.\n

        AUTHOR
                Written by Althea Hansel-Harris. Based on code by Stefano Forli, PhD and Andreas Tillack, PhD.\n

        REPORTING BUGS
                Please report bugs to:
                AutoDock mailing list   http://autodock.scripps.edu/mailing_list\n

        COPYRIGHT
                Copyright (C) 2021 Stefano Forli Laboratory, Center for Computational Structural Biology, 
                             The Scripps Research Institute.
                GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>
        """ 

        self.option_groups = [
                { "INPUT" : {
                    'desc': 'Specify input data',
                    'args': [ 
                        ('--file', {
                            'help':'ligand DLG or DLG.gz (compressed) file to filter', 
                            'action':'append', 'type':str, 'metavar': "FILENAME.DLG[.gz]", 'required':False},
                            ),
                        ('--file_path', {
                            'help': 'directory containing DLG or DLG.gz (compressed) files to filter', 
                            'action':'store', 'type':str, 'metavar':"DIRNAME"},
                            ),
                        ('--file_list', {'help': 'file containing the list of DLG or DLG.gz (compressed) files to filter; relative or absolute paths are allowed', 
                            'action':'store', 'type':str, 'metavar':"FILENAME"},
                            ),
                        ('--recursive', {'help': 'enable recursive directory scan when --dir is used', 
                            'action':'store_true', 'default': False },
                            ),
                        ('--pattern', {'help': 'specify which pattern to use when searching for files to process in directories [only with "--dir", default "*.dlg*"]', 
                            'action':'store', 'type': str, 'metavar': 'PATTERN','default': "*.dlg*" },
                            ),
                        ('--filters_file', {
                            'help':'specify a file containing the filter definitions; each line in the file must contain the definition of a single filter, using the keywords of each filter as variable names, e.g.: "eworst=-3.0", "vdw=A:THR:276". NOTE: properties defined here have the precedence on the other CLI options!', 
                            'action':'store', 'type':str, 'metavar': "FILTERS_FILE", 'default':None},
                         ),
                        ('--sql_db', {
                            'help':'specify a file containing a sqlite database', 
                            'action':'store', 'type':str, 'metavar': "SQL_DB", 'default':None},
                         ),
                            ],
                },},

                {'OUTPUT': {
                    'desc': 'Manage the type of data reported and where it is written.',
                    'args': [ 
                        ('--output_sql', {
                            'help':('Name for output SQLite file'), 
                            'action':'store', 'type':str, 'metavar': "[FILE_NAME].DB", 'default':"output.db"},
                            ),
                        ('--num_clusters', {
                            'help':('n: Store top pose for top n clusters'), 
                            'action':'store', 'type':int, 'default':3, 'metavar':'INT'},
                            ),
                        ('--log', {
                            'help':('by default, results are printed in the terminal (STDOUT); '
                                'if this option is used, ligands passing the filters will be written '
                                'to this file'), 
                            'action':'store', 'type':str, 'metavar': "[FILE_NAME].TXT", 'default':"output_log.txt"},
                            ),
                        ('--out_fields', {
                            'help':'defines which fields are used when reporting '
                            'the results (either to stdout or to the log file); in the output, values '
                            'are separated by a "|" symbol; fields are specified ' 
                            'as comma-separated values, e.g. "--out_fields=e,le,fname"; by '
                            'default, energies_binding (energy) and file name (fname) are reported; fname always reported in first column'
                            'available fields are:  '
                            '"e" (energies_binding), '
                            '"le" (ligand efficiency), '
                            '"delta" (delta energy from best pose), '
                            '"ref_rmsd" (RMSD to reference pose), '
                            '"e_inter" (intermolecular energy), '
                            '"e_vdw" (van der waals energy), '
                            '"e_elec" (electrostatic energy), '
                            '"e_intra" (intermolecular energy), '
                            '"n_interact" (number of interactions), '
                            '"interactions" (all interactions), '
                            '"fname" (full path filename), '
                            '"ligand_smile" , '
                            '"rank" (rank of ligand pose), '
                            '"run" (run number for ligand pose), '
                            #'"mname" (molecule name from the \'move\' kw in the DLG); '
                            #'"name" (basename of the dlg[.gz] file); '
                            '"hb" (hydrogen bonds); '
                            #'"hba" (hydrogen bond, ligand acceptor); '
                            #'"hbd" (hydrogen bond, ligand donor); '
                            #'"r" (reacted: 1 (yes), 0 (no));'
                            #'"ri" (reacted residue info); '
                            '"all" (print all poses passing filters). Fields are '
                            'printed in the order in which they are provided.',
                            'action':'store', 'type':str, 'metavar': "FIELD1,FIELD2,...", 
                            'default':'fname,e'},
                            ),
                        ('--export_poses_path', {
                            'help':('specify the path where to save poses of ligands passing the filters (PDBQT format); '
                                'if the directory does not exist, it will be created; if it already exist, it will throw '
                                'an error, unless the --overwrite is used  [QUESTION????: NOTE: the log file will be automatically saved in this path.]'),
                            'action':'store', 'type':str, 'metavar': "DIRECTORY_NAME", 'default':None},
                            ),
                        ('--no_header', {
                            'help':('by default, a commented header ("#") with a summary of the filters used '
                                'is written to both STDOUT and the log file; this option suppresses the header'),
                            'action':'store_true', 'default':False },
                            ),
                        ('--no_print', {
                            'help':('suppress printing the results to STDOUT. NOTE: runtime may be faster if no_print option used.'),
                            'action':'store_true', 'default':False },
                            ),
                        ('--plot', {
                            'help':('Makes best energy and le histograms and displays in command line. Makes scatterplot of LE vs Best Energy, saves as [filters_file].png or out.png if no filters_file given.'),
                            'action':'store_true', 'default':False},
                            ),
                        ('--overwrite', {
                            'help':('by default, if a log file exists, it doesn\'t get '
                                'overwritten and an error is returned; this option enable overwriting existing log files. Will also overwrite existing database'),
                            'action':'store_true', 'default':False },
                            ),
                        ('--log_distinct_ligands', {
                            'help':('by default, will output all poses passing filters, including multiple poses for the same ligand. '
                                'This flag will cause each ligand passing the filters to only be logged once, with the best pose.'),
                            'action':'store_true', 'default':False },
                            ),
                        ('--order_results', {
                            'help':'Stipulates how to order the results when written to the log file. By default will be ordered by order results were added to the database. ONLY TAKES ONE OPTION.'
                            'available fields are:  '
                            '"e" (energies_binding), '
                            '"le" (ligand efficiency), '
                            '"delta" (delta energy from best pose), '
                            '"ref_rmsd" (RMSD to reference pose), '
                            '"e_inter" (intermolecular energy), '
                            '"e_vdw" (van der waals energy), '
                            '"e_elec" (electrostatic energy), '
                            '"e_intra" (intermolecular energy), '
                            '"n_interact" (number of interactions), '
                            '"interactions" (all interactions), '
                            '"fname" (full path filename), '
                            '"rank" (rank of ligand pose), '
                            '"run" (run number for ligand pose), '
                            #'"mname" (molecule name from the \'move\' kw in the DLG); '
                            #'"name" (basename of the dlg[.gz] file); '
                            '"hb" (hydrogen bonds); ',
                            #'"hba" (hydrogen bond, ligand acceptor); '
                            #'"hbd" (hydrogen bond, ligand donor); '
                            'action':'store', 'type':str, 'metavar': "STRING", 
                            'default':None},
                            ),


                        ],
                    },},

                {'PROPERTY FILTERS': {
                    'desc': ('Specify energy and ligand efficiency filters'),
                    'args': [ 
                        ('--eworst', {
                            'help':'specify the worst energy value accepted', 
                            'action':'store', 'type':float, 'metavar': "FLOAT", 'default':-3.0},
                            ),
                        ('--ebest', {
                            'help':'specify the best energy value accepted', 
                            'action':'store', 'type':float, 'metavar': "FLOAT",},
                            ),
                        ('--leworst', {
                            'help':'specify the worst ligand efficiency value accepted', 
                            'action':'store', 'type':float, 'metavar': "FLOAT"},
                            ),
                        ('--lebest', {
                            'help':'specify the best ligand efficiency value accepted', 
                            'action':'store', 'type':float, 'metavar': "FLOAT",},
                            ),
                        ('--epercentile', {
                            'help':'specify the worst energy percentile accepted', 
                            'action':'store', 'type':float, 'metavar': "FLOAT"},
                            ),
                        ('--leffpercentile', {
                            'help':'specify the worst ligand efficiency percentile accepted', 
                            'action':'store', 'type':float, 'metavar': "FLOAT"},
                            )

                            ],}},

                {'LIGAND FILTERS': {
                    'desc': ('Specify filters on ligands, including substructures or names'),
                    'args': [ 
                        ('--name', {
                            'help':'specify ligand name(s). Will seach OR', 
                            'action':'store', 'type':str, 'default':None, 'metavar': "STRING"},
                            ),
                        ('--substructure', {
                            'help':'specify SMILES substring(s) to search for', 
                            'action':'store', 'type':str, 'default':None, 'metavar': "STRING"},
                            ),
                        ('--substruct_flag', {
                            'help':'specify whether to search AND or OR for substructures. Default OR', 
                            'action':'store', 'type':str, 'default':"OR", 'metavar': "STRING"},
                            )
                            ],}},
                {'INTERACTION FILTERS': {
                    'desc': ('Specify interaction filters, either by count or by specific residue interaction. '
                             'Residue specifications are described using CHAIN:RES:NUM:ATOM_NAME, '
                             'and any combination is allowed, e.g.: CHAIN:::, :RES::, ::NUM:, :::ATOM_NAME, :RES:NUM:, etc... '
                             'Unwanted interactions can be defined by prepending "-" to the residue specification, e.g. "-B:THR:276:". '
                             'Multiple residues can be specified in a single option by separating them with a comma (e.g.: --vdw=B:THR:276:,B:HIS:226:), '
                             'or by repeating the interaction options (e.g. --vdw=B:THR:276: --vdw=B:HIS:226: ).'),
                    'args': [ 
                        ('--vdw', {
                            'help':'define van der Waals interactions with residue', 
                            'action':'append', 'type':str, 'metavar': "[-][CHAIN]:[RES]:[NUM]:[ATOM_NAME]"},
                            ),
                        ('--hb', {
                            'help':'define HB (ligand acceptor or donor) interaction',
                            'action':'append', 'type':str, 'metavar': "[-][CHAIN]:[RES]:[NUM]:[ATOM_NAME]"},
                            ),
                        ('--hba', {
                            'help':'define HB acceptor (ligand as acceptor) interaction',
                            'action':'append', 'type':str, 'metavar': "[-][CHAIN]:[RES]:[NUM]:[ATOM_NAME]"},
                            ),
                        ('--hbd', {
                            'help':'define HB donor (ligand as donor) interaction',
                            'action':'append', 'type':str, 'metavar': "[-][CHAIN]:[RES]:[NUM]:[ATOM_NAME]"},
                            ),
                        ('--hb_count', {
                            'help':('accept ligands with at least the requested number of HB interactions. '
                                'If a negative number is provided, then accept ligands with no more than the '
                                ' requested number of interactions'),
                            'action':'store', 'type':int, 'metavar': "NUMBER"},
                            ),
                        # ('--react_any', {
                        #     'help':'check if ligand reacted with any residue',
                        #     'action':'store_true'},
                        #     ),
                        ('--react_res', {
                            'help':'check if ligand reacted with specified residue',
                            'action':'append', 'type':str, 'metavar': "[-][CHAIN]:[RES]:[NUM]:[ATOM_NAME]"},
                            ),
                        ('--react_count', {
                            'help':('accept ligands that reacted with at least the requested number of residues; '
                                'if a negative number is provided, then accept ligands that do not react '
                                'with more than the specified number of residues'),
                            'action':'store', 'type':int, 'metavar': "NUMBER"},
                            ),
                        ],
                    },},
            ]

        self._initialize_parser()
        self._process_sources()

    def _initialize_parser(self):
        # create parser
        parser = argparse.ArgumentParser(description=self.description, usage=self.usage, epilog=self.epilog,
                formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                )
        # populate options menu
        for group_class in self.option_groups:
            for group_name, group_info in group_class.items():
                group_desc = group_info['desc']
                group_args = group_info['args']
                # create new group 
                group = parser.add_argument_group(group_name,group_desc)
                # add options to the group
                for name, args in group_args:
                    group.add_argument(name, **args)
        # parse options
        cmdline_opts = sys.argv[1:]
        # hack to allow defining options from a file (argparse does not allow to do it in a clean way)
        if "--filters_file" in cmdline_opts:
            idx = cmdline_opts.index('--filters_file')
            cmdline_opts.pop(idx)
            ffile = cmdline_opts.pop(idx)
            cmdline_opts += self.read_filter_file(ffile)
            self.filter_file = ffile
        # add a function here to validate the cmdline (repeated options?)
        # validate policy of file > cmdline? (or vice versa?)
        parsed_opts = parser.parse_args(cmdline_opts)
        #print("PARSED OPTIONS", parsed_opts)
        self.process_options(parsed_opts)

    def read_filter_file(self, fname):
            """ parse the filter file to define filters """
            opts = []
            with open(fname, 'r') as fp:
                for l in fp.readlines():
                    l = l.strip()
                    # remove empty lines
                    if not len(l):
                        continue
                    # remove comments
                    if l[0] == "#":
                        continue
                    opts.append("--%s" % l.strip())
            self.cmdline_opts += opts

    def process_options(self, parsed_opts):
        """ convert command line options to the dict of filters """
        # check that required input options are provided
        file_sources = {} #'file':None, 'files_path':None, 'file_list':None}
        file_sources['file'] = parsed_opts.file
        if not parsed_opts.file_path is None:
            file_sources['file_path'] = {'path':parsed_opts.file_path, 'pattern':parsed_opts.pattern, 'recursive': parsed_opts.recursive}
        else:
            file_sources['file_path'] = None
        file_sources['file_list'] = parsed_opts.file_list
        print("")
        if (file_sources['file'] is None) and (file_sources['file_path'] is None) and (file_sources['file_list'] is None) and (parsed_opts.sql_db is None):
            print("*ERROR* at least one input option needs to be used:  --file, --file_path, --file_list, --sql_db")
            sys.exit(1)
        ##### parse output options
        output = {'fields': parsed_opts.out_fields,
                  'log': parsed_opts.log,
                  'header': not parsed_opts.no_header,
                  'stdout': not parsed_opts.no_print,
                  'overwrite': parsed_opts.overwrite,
                  'export_poses_path': parsed_opts.export_poses_path,
                  'plot': parsed_opts.plot,
                  'outfields': parsed_opts.out_fields
                }
        db_opts = {'num_clusters': parsed_opts.num_clusters, "order_results":parsed_opts.order_results, "log_distinct_ligands":parsed_opts.log_distinct_ligands}

        # if a path for saving poses is specified, then the log will be written there
        if not parsed_opts.export_poses_path is None:
            output['log'] = "%s%s%s" % (output['export_poses_path'], os.sep, output['log'])
        #print(">>>>>>>>>>>>>chekc that fields are recognized")
        if parsed_opts.log is None:
            print("*ERROR* print to STDOUT is disabled and no log file has been specified; at least one output source needs to be used.")
            sys.exit(1)
        ##### filters
        # property filters
        properties = {'eworst':None, 'ebest':None, 'leworst':None, 'lebest':None, 'epercentile':None, 'leffpercentile':None}
        for kw, _ in properties.items():
            properties[kw] = getattr(parsed_opts, kw)
        # interaction filters (residues)
        interactions = {}
        res_interactions_kw = [('vdw', 'V'), ('hb','H'), ('hba','HA'), ('hbd','HD'), ('react_res','R')]
        for opt, _type in res_interactions_kw:
            interactions[_type] = []
            res_list = getattr(parsed_opts, opt)
            if res_list is None:
                continue
            found_res = []
            for res in res_list:
                if "," in res:
                    for r in res.split(','):
                        found_res.append(r)
                else:
                    found_res.append(res)
            for res in found_res:
                wanted = True
                if not res.count(":") == 3:
                    print( ('*ERROR* [%s]: to specify a residue use '
                        'the format CHAIN:RES:NUM:ATOM_NAME. Any item can be omitted, '
                        'as far as the number of semicolons is always 3 '
                        '(e.g.: CHAIN:::, :RES::, CHAIN::NUM:, etc.)') % res )
                    sys.exit()
                if res[0] == '-':
                    res=res[1:]
                    wanted = False
                interactions[_type].append( (res, wanted) )
        # count interactions
        interactions_count = []
        count_kw = [('hb_count', ("hb_count")), ('react_count',('R'))]
        for kw, pool in count_kw:
            c = getattr(parsed_opts, kw, None)
            if c is None:
                continue
            interactions_count.append((pool, c))
        #make dictionary for ligand filters
        ligand_filters_kw = [('name', 'N'), ('substructure', 'S'), ('substruct_flag', 'F')]
        ligand_filters = {}
        filter_ligands_flag = True
        for kw, _type in ligand_filters_kw:
            ligand_filters[_type] = []
            ligand_filter_list = getattr(parsed_opts, kw)
            if ligand_filter_list is None:
                filter_ligands_flag = False
            else:
                ligand_filter_list = ligand_filter_list.split(",")
                for fil in ligand_filter_list:
                    ligand_filters[_type].append(fil)

        # TODO combine everything into a single parameters dictionary to be passed to DockingResultManager (with all the kws in its init)
        filters = {'properties':properties, 'interactions':interactions, 'interactions_count':interactions_count, "ligand_filters":ligand_filters, 'filter_ligands_flag':filter_ligands_flag, 'max_miss':0} 
        self.file_sources = file_sources 
        self.filters = filters
        self.out_opts = output
        self.num_clusters = parsed_opts.num_clusters
        if parsed_opts.sql_db != None:
            sqlFile = parsed_opts.sql_db
        else:
            sqlFile = parsed_opts.output_sql
            db_opts['write_db_flag'] = True
            if parsed_opts.overwrite: #confirm user wants to overwrite
                if os.path.exists(sqlFile): #check if database file already exists
                    os.remove(sqlFile)
        db_opts['sqlFile'] = sqlFile
        self.db_opts = db_opts

    def _process_sources(self):
        """ process the options for input files (parse dictionary) """
        sources = self.file_sources
        self.files_pool = []
        if sources['file'] is not None:
            # print("DRM> initialized with %d individual files" % len(sources['file']))
            self.files_pool = sources['file']
        # update the files pool with the all the files found in the path
        if sources['file_path'] is not None:
            # print("DRM> scanning path [%s] (recursive: %s, pattern '%s')" % (sources['file_path']['path'],
            #         str(sources['file_path']['recursive']), sources['file_path']['pattern']))
            self.scan_dir(sources['file_path']['path'], 
                          sources['file_path']['pattern'], 
                          sources['file_path']['recursive'])
        # update the files pool with the files specified in the files list
        if sources['file_list'] is not None:
            # print("DRM> searching for files listed in [%s]" % sources['file_list'])
            self.scan_file_list(sources['file_list'])

    def scan_dir(self, path, pattern, recursive=False):
        """ scan for valid output files in a directory 
            the pattern is used to glob files
            optionally, a recursive search is performed
        """
        print("-Scanning directory [%s] for DLG files (pattern:|%s|)" % (path, pattern))
        files = []
        if recursive:
            path = os.path.normpath(path)
            path = os.path.expanduser(path)
            for dirpath, dirnames, filenames in os.walk(path):
                files.extend(os.path.join(dirpath,f) for f in fnmatch.filter(filenames,'*'+pattern))
        else:
            files = glob(os.path.join(path, pattern))
        print("-Found %d files." % len(files))
        self.files_pool.extend(files)

    def scan_file_list(self, filename):
        """ read file names from file list """
        accepted = []
        c = 0
        with open(filename, 'r') as fp:
            for l in fp.readlines():
                l = l.strip()
                c+=1
                if os.path.isfile(l):
                    accepted.append(l)
                else:
                    print("Warning! file |%s| does not exist" % l)
        if len(accepted) == 0:
            print("*ERROR* No valid files were found when reading from |%s|" % filename )
            #sys.exit(1)
        print("# [ %5.3f%% files in list accepted (%d) ]" % ( len(accepted)/c *100, c ))
        self.files_pool.extend(accepted)

