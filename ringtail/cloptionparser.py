#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail command line option parser
#

import sys
import argparse
from glob import glob
import os
import fnmatch
import warnings


class CLOptionParser():

    def __init__(self):
        self.write_db_flag = False
        self.description = """Package for creating SQLite database from virtual screening DLGs and performing filtering on results."""
        self.usage = "Please see GitHub for full usage details."
        self.name = "Ringtail"
        self.epilog = """

        REQUIRED PACKAGES
                Requires RDkit, SciPy, Meeko.\n

        AUTHOR
                Written by Althea Hansel-Harris. Based on code by Stefano Forli, PhD, Andreas Tillack, PhD, and Diogo Santos-Martins, PhD.\n

        REPORTING BUGS
                Please report bugs to:
                AutoDock mailing list   http://autodock.scripps.edu/mailing_list\n

        COPYRIGHT
                Copyright (C) 2022 Stefano Forli Laboratory, Center for Computational Structural Biology, 
                             The Scripps Research Institute.
                GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>
        """

        self.option_groups = [
            {
                "INPUT": {
                    'desc':'Specify input data',
                    'args': [
                        (
                            '--file',
                            {
                                'help': 'ligand DLG(s) and receptor PDBQT file to save and filter. Compressed (.gz) files allowed. Only 1 receptor allowed.',
                                'action': 'append',
                                'type': str,
                                'metavar': "FILENAME.[DLG/PDBQT][.gz]",
                                'required': False
                            },
                        ),
                        (
                            '--file_path',
                            {
                                'help': 'directory(s) containing DLG and PDBQT files to save and filter. Compressed (.gz) files allowed',
                                'action': 'append',
                                'type': str,
                                'metavar': "DIRNAME"
                            },
                        ),
                        (
                            '--file_list',
                            {
                                'help': 'file(s) containing the list of DLG and PDBQT files to filter; relative or absolute paths are allowed. Compressed (.gz) files allowed',
                                'action': 'append',
                                'type': str,
                                'metavar': "FILENAME"
                            },
                        ),
                        (
                            '--save_receptors',
                            {
                                'help': 'Saves receptor PDBQTs to database. Receptor location(s) must be specied with in --file, --file_path directory or --file_list file',
                                'action': 'store_true',
                                'default': False
                            },
                        ),
                        (
                            '--recursive',
                            {
                                'help': 'enable recursive directory scan when --file_path is used',
                                'action': 'store_true',
                                'default': False
                            },
                        ),
                        (
                            '--pattern',
                            {
                                'help': 'specify which pattern to use when searching for DLG files to process in directories [only with "--file_path", default "*.dlg*"]',
                                'action': 'store',
                                'type': str,
                                'metavar': 'PATTERN',
                                'default': "*.dlg*"
                            },
                        ),
                        (
                            '--filters_file',
                            {
                                'help': 'specify a file containing the filter definitions; each line in the file must contain the definition of a single filter, using the keywords of each filter as variable names, e.g.: "eworst=-3.0", "vdw=A:THR:276". NOTE: properties defined here have the precedence on the other CLI options!',
                                'action': 'store',
                                'type': str,
                                'metavar': "FILTERS_FILE",
                                'default': None
                            },
                        ),
                        (
                            '--input_db',
                            {
                                'help': 'specify a database file to perform actions with',
                                'action': 'store',
                                'type': str,
                                'metavar': "DATABASE",
                                'default': None
                            },
                        ),
                        (
                            '--add_results',
                            {
                                'help': 'Add new results to an existing database, specified by --input_db',
                                'action': 'store_true',
                                'default': False
                            },
                        ),
                        (
                            '--conflict_handling',
                            {
                                'help': 'specify how conflicting Results rows should be handled when inserting into database. Options are "ignore" or "replace". Default behavior will duplicate entries.',
                                'action': 'store',
                                'type': str,
                                'metavar': "'ignore' or 'replace'",
                                'default': None
                            },
                        ),
                    ],
                },
            },
            {
                'OUTPUT': {
                    'desc':
                    'Manage the type of data reported and where it is written.',
                    'args': [
                        (
                            '--output_db',
                            {
                                'help': ('Name for output database file'),
                                'action': 'store',
                                'type': str,
                                'metavar': "[FILE_NAME].DB",
                                'default': "output.db"
                            },
                        ),
                        (
                            '--export_table_csv',
                            {
                                'help': ('Create csv of the requested database table. Output as <table_name>.csv'),
                                'action': 'store',
                                'type': str,
                                'metavar': "TABLE_NAME",
                                'default': None
                            },
                        ),
                        (
                            '--export_query_csv',
                            {
                                'help': ('Create csv of the requested SQL query. Output as query.csv. MUST BE PRE-FORMATTED IN SQL SYNTAX e.g. SELECT [columns] FROM [table] WHERE [conditions]'),
                                'action': 'store',
                                'type': str,
                                'metavar': "[VALID SQL QUERY]",
                                'default': None
                            },
                        ),
                        (
                            '--max_poses',
                            {
                                'help': ('n: Store top pose for top n clusters'),
                                'action': 'store',
                                'type': int,
                                'default': 3,
                                'metavar': 'INT'
                            },
                        ),
                        (
                            '--store_all_poses',
                            {
                                'help': ('Store all poses from input files. Overrides --max_poses'),
                                'action': 'store_true',
                                'default': False
                            },
                        ),
                        (
                            '--log',
                            {
                                'help': ('by default, results are saved in "output_log.txt"; '
                                 'if this option is used, ligands passing the filters will be written '
                                 'to specified file'),
                                'action': 'store',
                                'type': str,
                                'metavar': "[FILE_NAME].TXT",
                                'default': "output_log.txt"
                            },
                        ),
                        (
                            '--subset_name',
                            {
                                'help': ('Specify name for db view of passing results'),
                                'action': 'store',
                                'type': str,
                                'default': "passing_results",
                                'metavar': 'STRING'
                            },
                        ),
                        (
                            '--out_fields',
                            {
                                'help':
                                'defines which fields are used when reporting '
                                'the results (to stdout and to the log file); '
                                'fields are specified '
                                'as comma-separated values, e.g. "--out_fields=e,le,hb"; by '
                                'default, energies_binding (energy) and ligand name are reported; ligand always reported in first column'
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
                                '"ligand_smile" , '
                                '"rank" (rank of ligand pose), '
                                '"run" (run number for ligand pose), '
                                '"hb" (hydrogen bonds); '
                                'Fields are '
                                'printed in the order in which they are provided. Ligand name will always be returned and should not be specified',
                                'action':'store',
                                'type': str,
                                'metavar': "FIELD1,FIELD2,...",
                                'default':'e'
                            },
                        ),
                        (
                            '--data_from_subset',
                            {
                                'help': ('Write log of --out_fields data for subset specified by --subset_name. Must use with --no_filter.'),
                                'action': 'store_true',
                                'default': False
                            },
                        ),
                        (
                            '--export_poses_path',
                            {
                                'help':
                                ('specify the path where to save poses of ligands passing the filters (SDF format); '
                                 'if the directory does not exist, it will be created; if it already exist, it will throw '
                                 'an error, unless the --overwrite is used  NOTE: the log file will be automatically saved in this path.'
                                 'Ligands will be stored as SDF files, with the poses passing the filtering criteria first, followed the non-passing poses, in the order specified.'
                                 ),
                                'action': 'store',
                                'type': str,
                                'metavar': "DIRECTORY_NAME",
                                'default': None
                            },
                        ),
                        (
                            '--no_print',
                            {
                                'help':
                                ('suppress printing of results passing filtering criteria. NOTE: runtime may be faster if no_print option used.'
                                 ),
                                'action':
                                'store_true',
                                'default':
                                False
                            },
                        ),
                        (
                            '--no_filter',
                            {
                                'help': ('Do not perform any filtering.'),
                                'action': 'store_true',
                                'default': False
                            },
                        ),
                        (
                            '--plot',
                            {
                                'help':
                                ('Makes scatterplot of LE vs Best Energy, saves as [filters_file].png or out.png if no filters_file given.'
                                 ),
                                'action': 'store_true',
                                'default': False
                            },
                        ),
                        (
                            '--overwrite',
                            {
                                'help':
                                ('by default, if a log file exists, it doesn\'t get '
                                 'overwritten and an error is returned; this option enable overwriting existing log files. Will also overwrite existing database'
                                 ),
                                'action': 'store_true',
                                'default': False
                            },
                        ),
                        (
                            '--all_poses',
                            {
                                'help':
                                ('by default, will output only top-scoring pose passing filters per ligand. '
                                 'This flag will cause each pose passing the filters to be logged.'
                                 ),
                                'action': 'store_false',
                                'default': True
                            },
                        ),
                        (
                            '--order_results',
                            {
                                'help':
                                'Stipulates how to order the results when written to the log file. By default will be ordered by order results were added to the database. ONLY TAKES ONE OPTION.'
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
                                '"rank" (rank of ligand pose), '
                                '"run" (run number for ligand pose), '
                                '"hb" (hydrogen bonds); ',
                                'action': 'store',
                                'type': str,
                                'metavar': "STRING",
                                'default': None
                            },
                        ),
                    ],
                },
            },
            {
                'PROPERTY FILTERS': {
                    'desc': ('Specify energy and ligand efficiency filters'),
                    'args':
                    [(
                        '--eworst',
                        {
                            'help': 'specify the worst energy value accepted',
                            'action': 'store',
                            'type': float,
                            'metavar': "FLOAT"
                        },
                    ),
                     (
                         '--ebest',
                         {
                             'help': 'specify the best energy value accepted',
                             'action': 'store',
                             'type': float,
                             'metavar': "FLOAT",
                         },
                     ),
                     (
                         '--leworst',
                         {
                             'help':
                             'specify the worst ligand efficiency value accepted',
                             'action': 'store',
                             'type': float,
                             'metavar': "FLOAT"
                         },
                     ),
                     (
                         '--lebest',
                         {
                             'help':
                             'specify the best ligand efficiency value accepted',
                             'action': 'store',
                             'type': float,
                             'metavar': "FLOAT",
                         },
                     ),
                     (
                         '--epercentile',
                         {
                             'help':
                             'specify the worst energy percentile accepted. Express as percentage e.g. 1 for top 1 percent.',
                             'action': 'store',
                             'type': float,
                             'metavar': "FLOAT",
                             'default': 1.0
                         },
                     ),
                     (
                         '--leffpercentile',
                         {
                             'help':
                             'specify the worst ligand efficiency percentile accepted. Express as percentage e.g. 1 for top 1 percent.',
                             'action': 'store',
                             'type': float,
                             'metavar': "FLOAT"
                         },
                     )],
                }
            },
            {
                'LIGAND FILTERS': {
                    'desc':
                    ('Specify filters on ligands, including substructures or names'
                     ),
                    'args': [
                        (
                            '--name',
                            {
                                'help':
                                'specify ligand name(s). Will seach OR',
                                'action': 'store',
                                'type': str,
                                'default': None,
                                'metavar': "STRING"
                            },
                        ),
                        (
                            '--substructure',
                            {
                                'help':
                                'specify SMILES substring(s) to search for. Join multiple strings with "," ',
                                'action': 'store',
                                'type': str,
                                'default': None,
                                'metavar': "STRING"
                            },
                        ),
                        (
                            '--substruct_join',
                            {
                                'help':
                                'specify whether to search AND or OR for substructures. Default OR',
                                'action': 'store',
                                'type': str,
                                'default': "OR",
                                'metavar': "STRING"
                            },
                        )
                    ],
                }
            },
            {
                'INTERACTION FILTERS': {
                    'desc':
                    ('Specify interaction filters, either by count or by specific residue interaction. '
                     'Residue specifications are described using CHAIN:RES:NUM:ATOM_NAME, '
                     'and any combination is allowed, e.g.: CHAIN:::, :RES::, ::NUM:, :::ATOM_NAME, :RES:NUM:, etc... '
                     'Unwanted interactions can be defined by prepending "~" to the residue specification, e.g. "~B:THR:276:". '
                     'Multiple residues can be specified in a single option by separating them with a comma (e.g.: --vdw=B:THR:276:,B:HIS:226:), '
                     'or by repeating the interaction options (e.g. --vdw=B:THR:276: --vdw=B:HIS:226: ).'
                     ),
                    'args': [
                        (
                            '--vdw',
                            {
                                'help':
                                'define van der Waals interactions with residue',
                                'action': 'append',
                                'type': str,
                                'metavar': "[-][CHAIN]:[RES]:[NUM]:[ATOM_NAME]"
                            },
                        ),
                        (
                            '--hb',
                            {
                                'help':
                                'define HB (ligand acceptor or donor) interaction',
                                'action': 'append',
                                'type': str,
                                'metavar': "[-][CHAIN]:[RES]:[NUM]:[ATOM_NAME]"
                            },
                        ),
                        (
                            '--max_miss',
                            {
                                'help':
                                'Will separately log all possible combinations of interaction filters in log file excluding up to max_miss numer of interactions from given set. Cannot be used with --plot or --export_poses_path.',
                                'action': 'store',
                                'type': int,
                                'metavar': "INTEGER",
                                'default': 0
                            },
                        ),
                        (
                            '--hb_count',
                            {
                                'help':
                                ('accept ligands with at least the requested number of HB interactions. '
                                 'If a negative number is provided, then accept ligands with no more than the '
                                 ' requested number of interactions'),
                                'action':
                                'store',
                                'type':
                                int,
                                'metavar':
                                "NUMBER"
                            },
                        ),
                        (
                            '--react_any',
                            {
                                'help':
                                'check if ligand reacted with any residue',
                                'action': 'store_true',
                                'default': False
                            },
                        ),
                        (
                            '--react_res',
                            {
                                'help':
                                'check if ligand reacted with specified residue',
                                'action': 'append',
                                'type': str,
                                'metavar': "[-][CHAIN]:[RES]:[NUM]:[ATOM_NAME]"
                            },
                        ),
                        (
                            '--interaction_tolerance',
                            {
                                'help': ('Will add the interactions for poses within some tolerance RMSD range of the top pose in a cluster to that top pose. Can use as flag with default tolerance of 0.8, or give other value as desired'),
                                'action': 'store',
                                'type': float,
                                'metavar': "FLOAT",
                                'default': None,
                                'const': 0.8,
                                'nargs': '?'
                            },
                        ),
                    ],
                },
            },
        ]

        self.outfield_options = [
            "e", "le", "delta", "ref_rmsd", "e_inter", "e_vdw", "e_elec",
            "e_intra", "n_interact", "interactions", "fname", "ligand_smile",
            "rank", "run", "hb", "source_file"
        ]

        self._initialize_parser()
        self._process_sources()

        # confirm that files were found, else throw error
        # if only receptor files found and --save_receptor, assume we just want to add receptor and not modify the rest of the db, so turn off write_db_flag
        if len(self.lig_files_pool) == 0 and len(self.rec_files_pool) != 0 and self.save_receptors:
            self.db_opts["write_db_flag"] = False
            # raise error if not input db not given
            if self.input_db is None:
                raise RuntimeError("No input database given for saving receptor(s)")
        if len(self.lig_files_pool) == 0 and (self.db_opts["write_db_flag"] or self.db_opts["add_results"]):
            raise RuntimeError(
                "No ligand files found. Please check file source.")
        if len(self.rec_files_pool) > 1:
            raise RuntimeError("Found more than 1 receptor PDBQTs. Please check input files and only include receptor associated with DLGs")

    def _initialize_parser(self):
        # create parser
        parser = argparse.ArgumentParser(
            description=self.description,
            usage=self.usage,
            epilog=self.epilog,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        )
        # populate options menu
        for group_class in self.option_groups:
            for group_name, group_info in group_class.items():
                group_desc = group_info['desc']
                group_args = group_info['args']
                # create new group
                group = parser.add_argument_group(group_name, group_desc)
                # add options to the group
                for name, args in group_args:
                    group.add_argument(name, **args)
        # parse options
        cmdline_opts = sys.argv[1:]
        # allows defining options from file (argparse does not allow to do it in a clean way)
        if "--filters_file" in cmdline_opts:
            idx = cmdline_opts.index('--filters_file')
            cmdline_opts.pop(idx)
            ffile = cmdline_opts.pop(idx)
            cmdline_opts += self.read_filter_file(ffile)
            self.filter_file = ffile
        # add a function here to validate the cmdline (repeated options?)
        # validate policy of file > cmdline? (or vice versa?)
        parsed_opts = parser.parse_args(cmdline_opts)
        self.process_options(parsed_opts)

    def read_filter_file(self, fname):
        """ parse the filter file to define filters """
        opts = []
        with open(fname, 'r') as fp:
            for line in fp.readlines():
                line = line.strip()
                # remove empty lines
                if not len(line):
                    continue
                # remove comments
                if line[0] == "#":
                    continue
                opts.append("--%s" % line.strip())
        self.cmdline_opts += opts

    def process_options(self, parsed_opts):
        """ convert command line options to the dict of filters """
        # check that required input options are provided
        file_sources = {}  # 'file':None, 'files_path':None, 'file_list':None}
        file_sources['file'] = parsed_opts.file
        self.save_receptors = parsed_opts.save_receptors
        if parsed_opts.file_path is not None:
            file_sources['file_path'] = {
                'path': parsed_opts.file_path,
                'pattern': parsed_opts.pattern,
                'recursive': parsed_opts.recursive
            }
        else:
            file_sources['file_path'] = None
        file_sources['file_list'] = parsed_opts.file_list
        print("")
        if (file_sources['file'] is
                None) and (file_sources['file_path'] is
                           None) and (file_sources['file_list'] is
                                      None) and (parsed_opts.input_db is None):
            raise FileNotFoundError(
                "*ERROR* at least one input option needs to be used:  --file, --file_path, --file_list, --input_db"
            )
        if parsed_opts.add_results and parsed_opts.input_db is None:
            raise RuntimeError(
                "ERRROR! Must specify --input_db if adding results to an existing database"
            )
        if parsed_opts.max_miss < 0:
            raise RuntimeError("--max_miss must be greater than or equal to 0")
        if parsed_opts.max_miss > 0:
            if parsed_opts.plot:
                raise RuntimeError("Cannot use --plot with --max_miss > 0. Can plot for desired subset with --no_filter,--data_from_subset and, --subset_name.")
            if parsed_opts.export_poses_path is not None:
                raise RuntimeError("Cannot use --export_poses_path with --max_miss > 0. Can export poses for desired subset with --no_filter, --data_from_subset, and --subset_name")
        parsed_opts.out_fields = parsed_opts.out_fields.split(",")
        for outfield in parsed_opts.out_fields:
            if outfield not in self.outfield_options:
                raise RuntimeError(
                    "WARNING: {out_f} is not a valid output option. Please see --help or documentation"
                    .format(out_f=outfield))
        # parse output options
        # Make sure that export_poses_path has trailing /, is directory
        if parsed_opts.export_poses_path is not None:
            if not parsed_opts.export_poses_path.endswith("/"):
                parsed_opts.export_poses_path += "/"
            if not os.path.isdir(parsed_opts.export_poses_path):
                raise FileNotFoundError(
                    "--export_poses_path directory does not exist. Please create directory first"
                )
        # confirm that conflict_handling is an allowed option
        conflict_options = {"IGNORE", "REPLACE"}
        conflict_handling = parsed_opts.conflict_handling
        if conflict_handling is not None:
            conflict_handling = parsed_opts.conflict_handling.upper()
            if conflict_handling not in conflict_options:
                warnings.warn(f"--conflict_handing option {parsed_opts.conflict_handling} not allowed. Reverting to default behavior.")
                conflict_handling = None

        output = {
            'log': parsed_opts.log,
            'overwrite': parsed_opts.overwrite,
            'export_poses_path': parsed_opts.export_poses_path,
            'plot': parsed_opts.plot,
            'outfields': parsed_opts.out_fields,
            'no_print': parsed_opts.no_print,
            'no_filter': parsed_opts.no_filter,
            'data_from_subset': parsed_opts.data_from_subset,
            'export_table': parsed_opts.export_table_csv,
            'export_query': parsed_opts.export_query_csv,
        }
        db_opts = {
            'num_clusters': parsed_opts.max_poses,
            "order_results": parsed_opts.order_results,
            "log_distinct_ligands": parsed_opts.all_poses,
            "interaction_tolerance": parsed_opts.interaction_tolerance,
            "results_view_name": parsed_opts.subset_name,
            "store_all_poses": parsed_opts.store_all_poses,
            "overwrite": parsed_opts.overwrite,
            "add_results": parsed_opts.add_results,
            "conflict_opt": conflict_handling
        }

        # if a path for saving poses is specified, then the log will be written there
        if not parsed_opts.export_poses_path is None:
            output['log'] = os.path.join(output["export_poses_path"],
                                         output['log'])
        if parsed_opts.log is None:
            print(
                "*ERROR* print to STDOUT is disabled and no log file has been specified; at least one output source needs to be used."
            )
            sys.exit(1)
        # # # filters
        # property filters
        properties = {
            'eworst': None,
            'ebest': None,
            'leworst': None,
            'lebest': None,
            'epercentile': None,
            'leffpercentile': None
        }
        for kw, _ in properties.items():
            properties[kw] = getattr(parsed_opts, kw)
        # interaction filters (residues)
        interactions = {}
        res_interactions_kw = [('vdw', 'V'), ('hb', 'H'), ('react_res', 'R')]
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
                    print((
                        '*ERROR* [%s]: to specify a residue use '
                        'the format CHAIN:RES:NUM:ATOM_NAME. Any item can be omitted, '
                        'as long as the number of semicolons is always 3 '
                        '(e.g.: CHAIN:::, :RES::, CHAIN::NUM:, etc.)') % res)
                    sys.exit()
                if res[0] == '~':
                    res = res[1:]
                    wanted = False
                interactions[_type].append((res, wanted))
        # count interactions
        interactions_count = []
        count_kw = [('hb_count', ("hb_count")), ('react_count', ('R'))]
        for kw, pool in count_kw:
            c = getattr(parsed_opts, kw, None)
            if c is None:
                continue
            interactions_count.append((pool, c))
        # make dictionary for ligand filters
        ligand_filters_kw = [('name', 'N'), ('substructure', 'S'),
                             ('substruct_join', 'F')]
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

        filters = {
            'properties': properties,
            'interactions': interactions,
            'interactions_count': interactions_count,
            "ligand_filters": ligand_filters,
            'filter_ligands_flag': filter_ligands_flag,
            'max_miss': parsed_opts.max_miss,
            "react_any": parsed_opts.react_any
        }
        self.file_sources = file_sources
        self.input_db = parsed_opts.input_db
        if parsed_opts.input_db is not None:
            sqlFile = parsed_opts.input_db
            if not os.path.exists(sqlFile):
                print("WARNING: input database does not exist!")
                sys.exit(1)
            db_opts['write_db_flag'] = False
        else:
            sqlFile = parsed_opts.output_db
            db_opts['write_db_flag'] = True
        db_opts['sqlFile'] = sqlFile

        # make attributes for parsed opts
        self.db_opts = db_opts
        self.filters = filters
        self.out_opts = output

    def _process_sources(self):
        """ process the options for input files (parse dictionary) """
        sources = self.file_sources
        self.lig_files_pool = []
        self.rec_files_pool = []
        if sources['file'] is not None:
            self.lig_files_pool = list(filter(lambda file: file.endswith(".dlg") or file.endswith(".dlg.gz"), sources['file']))
            self.rec_files_pool = list(filter(lambda file: file.endswith(".pdbqt") or file.endswith(".pdbqt.gz"), sources['file']))
        # update the files pool with the all the files found in the path
        if sources['file_path'] is not None:
            for path in sources['file_path']['path']:
                # scan for ligand dlgs
                self.scan_dir(path,
                              sources['file_path']['pattern'],
                              sources['file_path']['recursive'])
                # scan for receptor pdbqts
                if self.save_receptors:
                    self.scan_dir(path,
                                  "*.pdbqt*",
                                  sources['file_path']['recursive'])
        # update the files pool with the files specified in the files list
        if sources['file_list'] is not None:
            for file in sources['file_list']:
                self.scan_file_list(file)

        if len(self.lig_files_pool) > 0 or len(self.rec_files_pool) > 0:
            print("-Found %d ligand files." % len(self.lig_files_pool))
            print("-Found %d receptor files." % len(self.rec_files_pool))

        # raise error if --save_receptors and none found
        if self.save_receptors and len(self.rec_files_pool) == 0:
            raise FileNotFoundError("--save_receptors flag specified but no receptor PDBQT found. Please check location of receptor file and file source options")

    def scan_dir(self, path, pattern, recursive=False):
        """ scan for valid output files in a directory
            the pattern is used to glob files
            optionally, a recursive search is performed
        """
        print("-Scanning directory [%s] for files (pattern:|%s|)" %
              (path, pattern))
        files = []
        if recursive:
            path = os.path.normpath(path)
            path = os.path.expanduser(path)
            for dirpath, dirnames, filenames in os.walk(path):
                files.extend(
                    os.path.join(dirpath, f)
                    for f in fnmatch.filter(filenames, '*' + pattern))
        else:
            files = glob(os.path.join(path, pattern))
        if "dlg" in pattern:
            self.lig_files_pool.extend(files)
        elif "pdbqt" in pattern:
            self.rec_files_pool.extend(files)

    def scan_file_list(self, filename):
        """ read file names from file list """
        lig_accepted = []
        rec_accepted = []
        c = 0
        with open(filename, 'r') as fp:
            for line in fp.readlines():
                line = line.strip()
                c += 1
                if os.path.isfile(line):
                    if line.endswith(".dlg") or line.endswith(".dlg.gz"):
                        lig_accepted.append(line)
                    if line.endswith(".pdbqt"):
                        rec_accepted.append(line)
                else:
                    print("Warning! file |%s| does not exist" % line)
        if len(lig_accepted) + len(rec_accepted) == 0:
            raise FileNotFoundError("*ERROR* No valid files were found when reading from |%s|" % filename)
        print("# [ %5.3f%% files in list accepted (%d) ]" %
              ((len(lig_accepted) + len(rec_accepted)) / c * 100, c))
        self.lig_files_pool.extend(lig_accepted)
        self.rec_files_pool.extend(rec_accepted)
