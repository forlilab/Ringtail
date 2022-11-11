#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail command line option parser
#

import sys
import argparse
import json
from glob import glob
import os
import fnmatch
import logging
from .exceptions import OptionError
import __main__
from .ringtailcore import RingtailCore


def cmdline_parser(defaults={}):

    conf_parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        add_help=False,
    )
    conf_parser.add_argument("-c", "--config")
    confargs, remaining_argv = conf_parser.parse_known_args()

    defaults = RingtailCore.get_defaults()
    """defaults = {
            'append_results': False,
            'order_results': None,
            'output_all_poses': None,
            'results_view_name': 'passing_results',
            'overwrite': None,
            'conflict_opt': None,
            'mode': 'dlg',
            'parser_manager': 'multiprocessing',
            'chunk_size': 1,
            'max_poses': 3,
            'store_all_poses': False,
            'interaction_tolerance': None,
            'target': None,
            'add_interactions': False,
            'interaction_cutoffs': [3.7, 4.0],
            'receptor_file': None,
            'file_sources': {
                'file': [[]],
                'file_path': {'path': [[]], 'pattern': '*.dlg*', 'recursive': None},
                'file_list': [[]]},
            'file_pattern': '*.dlg*',
            'max_proc': None,
            'log': 'output_log.txt',
            'export_sdf_path': None,
            'plot': None,
            'outfields': 'e',
            'no_print': True,
            'export_bookmark_csv': None,
            'export_query_csv': None,
            'export_bookmark_db': None,
            'data_from_bookmark': None,
            'filter_bookmark': None,
            'eworst': None,
            'ebest': None,
            'leworst': None,
            'lebest': None,
            'energy_percentile': None,
            'le_percentile': None, 
            'van_der_waals': [],
            'hydrogen_bond': [],
            'reactive_res': [],
            'interactions_count': [],
            'name': [],
            'substructure': [],
            'substructure_join': 'OR',
            'filter_ligands_flag': False,
            'max_miss': 0,
            'react_any': None,
            'save_receptor': False,
            'filter': False,
    }"""
    config = {}
    for section in defaults:
        for k, v in defaults[section]["values"].items():
            config[k] = v
    # separate filters into separate arguments
    filter_keys = ["properties", "interactions", "ligand_filters"]
    replace_filter_keys = {
        "V": "van_der_waals",
        "H": "hydrogen_bond",
        "R": "reactive_res",
        "N": "name",
        "S": "substructure",
        "F": "substructure_join",
    }

    for fk in filter_keys:
        filterdict = config.pop(fk)
        for k, v in filterdict.items():
            if k in replace_filter_keys:
                config[replace_filter_keys[k]] = v
            else:
                config[k] = v

    if confargs.config is not None:
        logging.info("Reading options from config file")
        with open(confargs.config) as f:
            c = json.load(f)
            config.update(c)

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        usage="""rt_process_vs.py [-c CONFIG] write|read OPTS  \n
                 Script to process virtual screening data by writing/filtering/exporting from database. In 'write' mode, provide docking files for writing to database with --file, --file_list, and/or --file_path. In 'read' mode, provide input database with --input_db. Please see GitHub for full usage details.""",
        description="Package for creating, managing, and filtering databases of virtual screening results.",
        epilog="""
        AUTHORS
                Written by Althea Hansel-Harris. Based on code by Stefano Forli, PhD, Andreas Tillack, PhD, Diogo Santos-Martins, PhD, and Matthew Holcomb, PhD.\n

        REPORTING BUGS
                Please report bugs to:
                AutoDock mailing list   http://autodock.scripps.edu/mailing_list\n

        COPYRIGHT
                Copyright (C) 2022 Forli Lab, Center for Computational Structural Biology,
                             The Scripps Research Institute.
                GNU L-GPL version 2.1 or later <https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html> \n
        """,
        exit_on_error=False,
    )

    subparsers = parser.add_subparsers(
        help="Specify if should write to or read from database. To show options of each mode use the in-line help, e.g.: %s read -h"
        % os.path.basename(__main__.__file__),
        dest="rr_mode",
    )

    write_parser = subparsers.add_parser("write")
    write_parser.add_argument(
        "-i",
        "--input_db",
        help="specify a database file to perform actions with",
        action="store",
        type=str,
        metavar="DATABASE",
    )
    write_parser.add_argument(
        "-s",
        "--bookmark_name",
        help="Specify name for db view of passing results to create (write mode) or export from (read mode)",
        action="store",
        type=str,
        metavar="STRING",
    )
    write_parser.add_argument(
        "-m",
        "--mode",
        help='specify AutoDock program used to generate results. Available options are "DLG" and "Vina". Vina mode will automatically change --pattern to *.pdbqt',
        action="store",
        type=str,
        metavar="[dlg] or [vina]",
    )
    write_parser.add_argument(
        "-v",
        "--verbose",
        help="Print results passing filtering criteria to STDOUT. NOTE: runtime may be slower option used.",
        action="store_true",
    )
    write_parser.add_argument(
        "-d",
        "--debug",
        help="Print additional error information to STDOUT.",
        action="store_true",
    )
    write_parser.add_argument(
        "-f",
        "--file",
        help="ligand docking output file to save. Compressed (.gz) files allowed. Only 1 receptor allowed.",
        action="append",
        type=str,
        metavar="FILENAME.[DLG/PDBQT][.gz]",
        nargs="+",
    )
    write_parser.add_argument(
        "-fp",
        "--file_path",
        help="directory(s) containing docking output files to save. Compressed (.gz) files allowed",
        action="append",
        type=str,
        metavar="DIRNAME",
        nargs="+",
    )
    write_parser.add_argument(
        "-fl",
        "--file_list",
        help="file(s) containing the list of docking output files to save; relative or absolute paths are allowed. Compressed (.gz) files allowed",
        action="append",
        type=str,
        metavar="FILENAME",
        nargs="+",
    )
    write_parser.add_argument(
        "-p",
        "--pattern",
        help='specify which pattern to use when searching for result files to process [only with "--file_path"]',
        action="store",
        type=str,
        metavar="PATTERN",
        dest="file_pattern",
    )
    write_parser.add_argument(
        "-r",
        "--recursive",
        help="enable recursive directory scan when --file_path is used",
        action="store_true",
    )
    write_parser.add_argument(
        "-a",
        "--append_results",
        help="Add new results to an existing database, specified by --input_db",
        action="store_true",
    )
    write_parser.add_argument(
        "-dh",
        "--duplicate_handling",
        help='specify how duplicate Results rows should be handled when inserting into database. Options are "ignore" or "replace". Default behavior will allow duplicate entries.',
        action="store",
        type=str,
        metavar="'ignore' or 'replace'",
    )
    write_parser.add_argument(
        "-sr",
        "--save_receptor",
        help="Saves receptor PDBQT to database. Receptor location must be specied with in --file, --file_path directory or --file_list file",
        action="store_true",
    )
    write_parser.add_argument(
        "-o",
        "--output_db",
        help="Name for output database file",
        action="store",
        type=str,
        metavar="[FILE_NAME].DB",
        default="output.db",
    )
    write_parser.add_argument(
        "-ov",
        "--overwrite",
        help="by default, if a log file exists, it doesn't get overwritten and an error is returned; this option enable overwriting existing log files. Will also overwrite existing database",
        action="store_true",
    )
    write_parser.add_argument(
        "-mp",
        "--max_poses",
        help="n: Store top pose for top n clusters",
        action="store",
        type=int,
        metavar="INT",
    )
    write_parser.add_argument(
        "-sa",
        "--store_all_poses",
        help="Store all poses from input files. Overrides --max_poses",
        action="store_true",
    )
    write_parser.add_argument(
        "-it",
        "--interaction_tolerance",
        help="Will add the interactions for poses within some tolerance RMSD range of the top pose in a cluster to that top pose. Can use as flag with default tolerance of 0.8, or give other value as desired. Only compatible with ADGPU mode",
        action="store",
        type=float,
        metavar="FLOAT",
        const=0.8,
        nargs="?",
    )
    write_parser.add_argument(
        "-ai",
        "--add_interactions",
        help="Find interactions between ligand poses and receptor and save to database. Requires receptor PDBQT to be given with input files (all modes) and --receptor_file to be specified with Vina mode. SIGNIFICANTLY INCREASES DATBASE WRITE TIME.",
        action="store_true",
    )
    write_parser.add_argument(
        "-ic",
        "--interaction_cutoffs",
        help="Use with --add_interactions, specify distance cutoffs for measuring interactions between ligand and receptor in angstroms. Give as string, separating cutoffs for hydrogen bonds and VDW with comma (in that order). E.g. '-ic 3.7,4.0' will set the cutoff for hydrogen bonds to 3.7 angstroms and for VDW to 4.0. These are the default cutoffs.",
        action="store",
        type=str,
        metavar="[HB CUTOFF],[VDW CUTOFF]",
    )
    write_parser.add_argument(
        "-rf",
        "--receptor_file",
        help="Use with Vina mode. Give file for receptor PDBQT.",
        action="store",
        type=str,
        metavar="STRING",
    )
    write_parser.add_argument(
        "-mpr",
        "--max_proc",
        help="Maximum number of processes to create during parallel file parsing. Defaults to number of CPU processors.",
        action="store",
        type=int,
    )

    read_parser = subparsers.add_parser("read")
    read_parser.add_argument(
        "-i",
        "--input_db",
        help="specify a database file to perform actions with",
        action="store",
        type=str,
        metavar="DATABASE",
    )
    read_parser.add_argument(
        "-s",
        "--bookmark_name",
        help="Specify name for db view of passing results to create or export from",
        action="store",
        type=str,
        metavar="STRING",
        dest="results_view_name",
    )
    read_parser.add_argument(
        "-m",
        "--mode",
        help='specify AutoDock program used to generate results. Available options are "DLG" and "Vina". Vina mode will automatically change --pattern to *.pdbqt',
        action="store",
        type=str,
        metavar="[dlg] or [vina]",
    )
    read_parser.add_argument(
        "-v",
        "--verbose",
        help="Print results passing filtering criteria to STDOUT. NOTE: runtime may be slower option used.",
        action="store_true",
    )
    read_parser.add_argument(
        "-d",
        "--debug",
        help="Print additional error information to STDOUT.",
        action="store_true",
    )

    output_group = read_parser.add_argument_group("Output options")
    output_group.add_argument(
        "-l",
        "--log",
        help='by default, results are saved in "output_log.txt"; if this option is used, ligands and requested info passing the filters will be written to specified file',
        action="store",
        type=str,
        metavar="[FILE_NAME].TXT",
    )
    output_group.add_argument(
        "-of",
        "--outfields",
        help=(
            'defines which fields are used when reporting the results (to stdout and to the log file); fields are specified as comma-separated values, e.g. "--outfields=e,le,hb"; by default, energies_binding (energy) and ligand name are reported; ligand always reported in first column available fields are:  '
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
            '"hb" (hydrogen bonds), '
            '"receptor" (receptor name); '
            "Fields are "
            "printed in the order in which they are provided. Ligand name will always be returned and should not be specified"
        ),
        action="store",
        type=str,
        metavar="FIELD1,FIELD2,...",
    )
    output_group.add_argument(
        "-ord",
        "--order_results",
        help="Stipulates how to order the results when written to the log file. By default will be ordered by order results were added to the database. ONLY TAKES ONE OPTION."
        "available fields are:  "
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
        action="store",
        type=str,
        metavar="STRING",
    )
    output_group.add_argument(
        "-oap",
        "--output_all_poses",
        help="By default, will output only top-scoring pose passing filters per ligand. This flag will cause each pose passing the filters to be logged.",
        action="store_true",
    )
    output_group.add_argument(
        "-xs",
        "--export_bookmark_csv",
        help="Create csv of the bookmark given with bookmark_name. Output as <bookmark_name>.csv. Can also export full database tables",
        action="store",
        type=str,
        metavar="BOOKMARK_NAME",
    )
    output_group.add_argument(
        "-xq",
        "--export_query_csv",
        help="Create csv of the requested SQL query. Output as query.csv. MUST BE PRE-FORMATTED IN SQL SYNTAX e.g. SELECT [columns] FROM [table] WHERE [conditions]",
        action="store",
        type=str,
        metavar="[VALID SQL QUERY]",
    )
    output_group.add_argument(
        "-sdf",
        "--export_sdf_path",
        help="specify the path where to save poses of ligands passing the filters (SDF format); if the directory does not exist, it will be created; if it already exist, it will throw an error, unless the --overwrite is used  NOTE: the log file will be automatically saved in this path. Ligands will be stored as SDF files in the order specified.",
        action="store",
        type=str,
        metavar="DIRECTORY_NAME",
    )
    output_group.add_argument(
        "-xdb",
        "--export_bookmark_db",
        help="Export a database containing only the results found in the bookmark specified by --bookmark_name. Will save as <input_db>_<bookmark_name>.db",
        action="store_true",
    )
    output_group.add_argument(
        "-nd",
        "--data_from_bookmark",
        help="Write log of --outfields data for bookmark specified by --bookmark_name. Must use without any filters.",
        action="store_true",
    )
    output_group.add_argument(
        "-fb",
        "--filter_bookmark",
        help="Perform filtering over specified bookmark.",
        action="store",
        type=str,
    )
    output_group.add_argument(
        "-p",
        "--plot",
        help="Makes scatterplot of LE vs Best Energy, saves as scatter.png.",
        action="store_true",
    )

    properties_group = read_parser.add_argument_group(
        "Property Filters", "Specify energy and ligand efficiency filters"
    )
    properties_group.add_argument(
        "-e",
        "--eworst",
        help="specify the worst energy value accepted",
        action="store",
        type=float,
        metavar="FLOAT",
    )
    properties_group.add_argument(
        "-eb",
        "--ebest",
        help="specify the best energy value accepted",
        action="store",
        type=float,
        metavar="FLOAT",
    )
    properties_group.add_argument(
        "-le",
        "--leworst",
        help="specify the worst ligand efficiency value accepted",
        action="store",
        type=float,
        metavar="FLOAT",
    )
    properties_group.add_argument(
        "-leb",
        "--lebest",
        help="specify the best ligand efficiency value accepted",
        action="store",
        type=float,
        metavar="FLOAT",
    )
    properties_group.add_argument(
        "-pe",
        "--energy_percentile",
        help="specify the worst energy percentile accepted. Express as percentage e.g. 1 for top 1 percent.",
        action="store",
        type=float,
        metavar="FLOAT",
    )
    properties_group.add_argument(
        "-ple",
        "--le_percentile",
        help="specify the worst ligand efficiency percentile accepted. Express as percentage e.g. 1 for top 1 percent.",
        action="store",
        type=float,
        metavar="FLOAT",
    )

    ligand_group = read_parser.add_argument_group(
        "Ligand Filters", "Specify ligand name filter(s)"
    )
    ligand_group.add_argument(
        "-n",
        "--name",
        help="specify ligand name(s). Will combine name filters with OR",
        action="store",
        type=str,
        metavar="STRING",
        nargs="+",
    )
    # SUBSTRUCTURE SEARCH: To be implemented in future version via RDKit
    # ligand_group.add_argument(
    #    "-st",
    #    "--substructure",
    #    help="specify SMILES substring(s) to search for.",
    #    action="store",
    #    type=str,
    #    metavar="STRING",
    #    nargs="+",
    # )
    # ligand_group.add_argument(
    #    "-sj",
    #    "--substructure_join",
    #    help="specify whether to join substructures filters with AND or OR.",
    #    action="store",
    #    type=str,
    #    metavar="STRING",
    # )

    interaction_group = read_parser.add_argument_group(
        "Interaction Filters",
        'Specify interaction filters, either by count or by specific residue interaction. Residue specifications are described using CHAIN:RES:NUM:ATOM_NAME, and any combination is allowed, e.g.: CHAIN:::, :RES::, ::NUM:, :::ATOM_NAME, :RES:NUM:, etc... Unwanted interactions can be defined by prepending "~" to the residue specification, e.g. "~B:THR:276:". Multiple residues can be specified in a single option by separating them with a space (e.g.: --vdw=B:THR:276: B:HIS:226:), or by repeating the interaction options (e.g. --vdw=B:THR:276: --vdw=B:HIS:226: ).',
    )
    interaction_group.add_argument(
        "-vdw",
        "--van_der_waals",
        help="define van der Waals interactions with residue",
        action="append",
        type=str,
        metavar="[-][CHAIN]:[RES]:[NUM]:[ATOM_NAME]",
    )
    interaction_group.add_argument(
        "-hb",
        "--hydrogen_bond",
        help="define HB (ligand acceptor or donor) interaction",
        action="append",
        type=str,
        metavar="[-][CHAIN]:[RES]:[NUM]:[ATOM_NAME]",
    )
    interaction_group.add_argument(
        "-r",
        "--reactive_res",
        help="check if ligand reacted with specified residue",
        action="append",
        type=str,
        metavar="[-][CHAIN]:[RES]:[NUM]:[ATOM_NAME]",
    )
    interaction_group.add_argument(
        "-hc",
        "--hb_count",
        help="accept ligands with at least the requested number of HB interactions. If a negative number is provided, then accept ligands with no more than the requested number of interactions",
        action="store",
        type=int,
        metavar="NUMBER",
    )
    interaction_group.add_argument(
        "-ra",
        "--react_any",
        help="check if ligand reacted with any residue",
        action="store_true",
    )
    interaction_group.add_argument(
        "-mm",
        "--max_miss",
        help="Will separately log all possible combinations of interaction filters in log file excluding up to max_miss numer of interactions from given set. Cannot be used with --plot or --export_sdf_path.",
        action="store",
        type=int,
        metavar="INTEGER",
    )

    # catch if running with no options
    if len(sys.argv) == 1:
        parser.print_help()
        raise OptionError("Script called with no commandline options. Aborting.")

    parser.set_defaults(**config)
    write_parser.set_defaults(**config)
    read_parser.set_defaults(**config)
    args = parser.parse_args(remaining_argv)

    return args, parser, conf_parser, write_parser, read_parser


class CLOptionParser:
    def __init__(self):

        self.outfield_options = [
            "e",
            "le",
            "delta",
            "ref_rmsd",
            "e_inter",
            "e_vdw",
            "e_elec",
            "e_intra",
            "n_interact",
            "interactions",
            "fname",
            "ligand_smile",
            "rank",
            "run",
            "hb",
            "source_file",
        ]
        self.order_options = {
            "e",
            "le",
            "delta",
            "ref_rmsd",
            "e_inter",
            "e_vdw",
            "e_elec",
            "e_intra",
            "n_interact",
            "rank",
            "run",
            "hb",
        }

        self._initialize_parser()

        # confirm that receptor file was found if needed, else throw error
        # if only receptor files found and --save_receptor, assume we just want to
        # add receptor and not modify the rest of the db
        if self.save_receptor:
            if self.storage_opts["append_results"]:
                raise OptionError(
                    "Cannot use --append_results with --save_receptor. Please remove the --save_receptor flag"
                )
            if self.receptor_file is None:
                raise OptionError(
                    "Must provide path for receptor PDBQT file is using --save_receptor"
                )

        # raise error if --save_receptor or --add_interactions and receptor not found
        if (self.save_receptor or self.add_interactions) and not os.path.exists(
            self.receptor_file
        ):
            raise OptionError(
                "--save_receptor or --add_interaction flag specified but no receptor PDBQT found. Please check location of receptor file and --receptor_file option"
            )

    def _initialize_parser(self):
        # create parser
        try:
            (
                parsed_opts,
                self.parser,
                self.conf_parser,
                self.write_parser,
                self.read_parser,
            ) = cmdline_parser()
            self.process_options(parsed_opts)
        except argparse.ArgumentError as e:
            self.parser.print_help()
            logging.error("\n")
            raise OptionError(
                "Invalid option or option ordering. Be sure to put read/write mode before any other arguments"
            ) from e
        except Exception as e:
            try:
                if self.rr_mode == "write":
                    self.write_parser.print_help()
                elif self.rr_mode == "read":
                    self.read_parser.print_help()
            finally:
                logging.error("\n")
                raise e

    def read_filter_file(self, fname):
        """parse the filter file to define filters"""
        opts = []
        with open(fname, "r") as fp:
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
        """convert command line options to the dict of filters"""
        self.debug = parsed_opts.debug
        if self.debug:
            logging.getLogger().setLevel(logging.DEBUG)
        self.filter = False  # set flag indicating if any filters given
        conflict_handling = None
        file_sources = None
        self.rec_files_pool = []
        filters = {}
        storage_opts = {}
        out_opts = {}
        # make sure mode is allowed
        allowed_modes = {"dlg", "vina"}
        self.mode = parsed_opts.mode.lower()
        # set receptor name, add_interactions if given
        self.rec_name = parsed_opts.receptor_file
        self.add_interactions = parsed_opts.add_interactions
        if self.mode not in allowed_modes:
            raise OptionError(
                "Given mode {0} not allowed. Please be sure that requested mode is 'vina' or 'dlg'".format(
                    self.mode
                )
            )
        if self.mode == "vina":
            # Guard against non-compatible options being called in Vina mode
            if parsed_opts.react_any:
                logging.warning(
                    "Cannot use reaction filters with Vina mode. Removing react_any filter."
                )
                parsed_opts.react_any = False
            if parsed_opts.interaction_tolerance is not None:
                logging.warning(
                    "Cannot use interaction_tolerance with Vina mode. Removing interaction_tolerance."
                )
                parsed_opts.interaction_tolerance = None
            if parsed_opts.add_interactions and parsed_opts.receptor_file is None:
                raise OptionError(
                    "Gave --add_interactions with Vina mode but did not specify receptor name. Please give receptor pdbqt name with --receptor_file."
                )
            # set pattern to .pdbqt
            parsed_opts.file_pattern = "*.pdbqt*"

        # guard against unsing interaction tolerance with storing all poses
        if (
            parsed_opts.interaction_tolerance is not None
            and parsed_opts.store_all_poses
        ):
            logging.warning(
                "Cannot use interaction_tolerance with store_all_poses. Removing interaction_tolerance."
            )
            parsed_opts.interaction_tolerance = False

        # guard against unsing percentile filter with all_poses
        if parsed_opts.output_all_poses and (
            parsed_opts.energy_percentile is not None
            or parsed_opts.le_percentile is not None
        ):
            logging.warning(
                "Cannot return all passing poses with percentile filter. Will only log best pose."
            )
            parsed_opts.output_all_poses = False

        # check options for write mode
        self.rr_mode = parsed_opts.rr_mode
        self.pattern = parsed_opts.file_pattern

        if self.rr_mode is None:
            raise OptionError(
                "No mode specified for rt_process_vs.py. Please specify mode (write/read)."
            )
        if self.rr_mode == "write":
            self.save_receptor = parsed_opts.save_receptor
            self.receptor_file = parsed_opts.receptor_file
            # check that required input options are provided
            file_sources = {}  # 'file':None, 'files_path':None, 'file_list':None}
            if parsed_opts.file is None:
                file_sources["file"] = [[]]
            else:
                file_sources["file"] = parsed_opts.file
            if parsed_opts.file_path is not None:
                file_sources["file_path"] = {
                    "path": parsed_opts.file_path,
                    "pattern": parsed_opts.file_pattern,
                    "recursive": parsed_opts.recursive,
                }
            else:
                file_sources["file_path"] = {
                    "path": [[]],
                    "pattern": "*.dlg*",
                    "recursive": None,
                }
            if parsed_opts.file_list is None:
                file_sources["file_list"] = [[]]
            else:
                file_sources["file_list"] = parsed_opts.file_list
            if (
                (file_sources["file"] == [[]])
                and (file_sources["file_path"] == [[]])
                and (file_sources["file_list"] == [[]])
                and (parsed_opts.input_db is None)
            ):
                raise OptionError(
                    "At least one input option needs to be used:  --file, --file_path, --file_list, --input_db"
                )
            if parsed_opts.append_results and parsed_opts.input_db is None:
                raise OptionError(
                    "Must specify --input_db if adding results to an existing database"
                )
            # confirm that conflict_handling is an allowed option
            conflict_options = {"IGNORE", "REPLACE"}
            if parsed_opts.duplicate_handling is not None:
                conflict_handling = parsed_opts.duplicate_handling.upper()
                if conflict_handling not in conflict_options:
                    logging.warning(
                        f"--conflict_handing option {parsed_opts.duplicate_handling} not allowed. Reverting to default behavior."
                    )
                    conflict_handling = None
        else:
            self.save_receptor = False
            self.receptor_file = None
            if parsed_opts.input_db is None:
                raise OptionError(
                    "No input database specified in read mode. Please specify database with --input_db"
                )
            if parsed_opts.max_miss < 0:
                raise OptionError("--max_miss must be greater than or equal to 0")
            if parsed_opts.max_miss > 0:
                if parsed_opts.plot:
                    raise OptionError(
                        "Cannot use --plot with --max_miss > 0. Can plot for desired bookmark with --bookmark_name."
                    )
                if parsed_opts.export_sdf_path is not None:
                    raise OptionError(
                        "Cannot use --export_sdf_path with --max_miss > 0. Can export poses for desired bookmark --bookmark_name"
                    )
            outfields_list = parsed_opts.outfields.split(",")
            for outfield in outfields_list:
                if outfield not in self.outfield_options:
                    raise OptionError(
                        "WARNING: {out_f} is not a valid output option. Please see --help".format(
                            out_f=outfield
                        )
                    )
            if parsed_opts.order_results is not None:
                if parsed_opts.order_results not in self.order_options:
                    raise OptionError(
                        "Requested ording option that is not available. Please see --help for available options."
                    )
            # parse output options
            # Make sure that export_sdf_path has trailing /, is directory
            if parsed_opts.export_sdf_path is not None:
                if not parsed_opts.export_sdf_path.endswith("/"):
                    parsed_opts.export_sdf_path += "/"
                if not os.path.isdir(parsed_opts.export_sdf_path):
                    raise OptionError(
                        "--export_sdf_path directory does not exist. Please create directory first"
                    )

            # # # filters
            if parsed_opts.filter_bookmark == parsed_opts.results_view_name:
                logging.error(
                    f"Specified filter_bookmark and bookmark_name are the same: {parsed_opts.results_view_name}"
                )
                raise OptionError(
                    "--filter_bookmark and --bookmark_name cannot be the same! Please rename --bookmark_name"
                )
            # property filters
            properties = {
                "eworst": None,
                "ebest": None,
                "leworst": None,
                "lebest": None,
                "energy_percentile": None,
                "le_percentile": None,
            }
            for kw, _ in properties.items():
                properties[kw] = getattr(parsed_opts, kw)
                if properties[kw] is not None:
                    self.filter = True
            # Cannot use energy/le cuttoffs with percentiles. Override percentile with given cutoff
            if (
                properties["eworst"] is not None
                and properties["energy_percentile"] is not None
            ):
                logging.warning(
                    "Cannot use --eworst cutoff with --energy_percentile. Overiding energy_percentile with eworst."
                )
                properties["energy_percentile"] = None
            if (
                properties["leworst"] is not None
                and properties["le_percentile"] is not None
            ):
                logging.warning(
                    "Cannot use --leworst cutoff with --le_percentile. Overiding le_percentile with leworst."
                )
                properties["le_percentile"] = None
            if (
                properties["energy_percentile"] is not None
                or properties["le_percentile"] is not None
            ) and parsed_opts.filter_bookmark is not None:
                raise OptionError(
                    "Cannot use --energy_percentile or --le_percentile with --filter_bookmark."
                )
            # interaction filters (residues)
            interactions = {}
            res_interactions_kw = [
                ("van_der_waals", "V"),
                ("hydrogen_bond", "H"),
                ("reactive_res", "R"),
            ]
            for opt, _type in res_interactions_kw:
                interactions[_type] = []
                res_list = getattr(parsed_opts, opt)
                if res_list is None:
                    continue
                found_res = []
                for res in res_list:
                    if "," in res:
                        for r in res.split(","):
                            found_res.append(r)
                    else:
                        found_res.append(res)
                for res in found_res:
                    wanted = True
                    if not res.count(":") == 3:
                        raise OptionError(
                            (
                                "[%s]: to specify a residue use "
                                "the format CHAIN:RES:NUM:ATOM_NAME. Any item can be omitted, "
                                "as long as the number of semicolons is always 3 "
                                "(e.g.: CHAIN:::, :RES::, CHAIN::NUM:, etc.)"
                            )
                            % res
                        )
                    if res[0] == "~":
                        res = res[1:]
                        wanted = False
                    interactions[_type].append((res, wanted))
                    self.filter = True
            # count interactions
            interactions_count = []
            count_kw = [("hb_count", ("hb_count")), ("react_count", ("R"))]
            for kw, pool in count_kw:
                c = getattr(parsed_opts, kw, None)
                if c is None:
                    continue
                interactions_count.append((pool, c))
                self.filter = True
            if parsed_opts.react_any is not None:
                self.filter = True
            # make dictionary for ligand filters
            ligand_filters_kw = [("name", "N"), ("substructure", "S")]
            ligand_filters = {}
            filter_ligands_flag = True
            ligand_filter_list = []
            for kw, _type in ligand_filters_kw:
                ligand_filters[_type] = []
                ligand_filter_list = getattr(parsed_opts, kw)
                if ligand_filter_list is None:
                    continue
                for fil in ligand_filter_list:
                    ligand_filters[_type].append(fil)
            ligand_filters["F"] = getattr(parsed_opts, "substructure_join")
            if ligand_filters["N"] == [] and ligand_filters["S"] == []:
                filter_ligands_flag = False
            if filter_ligands_flag:
                self.filter = True

            filters = {
                "properties": properties,
                "interactions": interactions,
                "interactions_count": interactions_count,
                "ligand_filters": ligand_filters,
                "filter_ligands_flag": filter_ligands_flag,
                "max_miss": parsed_opts.max_miss,
                "react_any": parsed_opts.react_any,
            }

        out_opts = {
            "log": parsed_opts.log,
            "export_sdf_path": parsed_opts.export_sdf_path,
            "plot": parsed_opts.plot,
            "outfields": parsed_opts.outfields,
            "no_print": not parsed_opts.verbose,
            "export_bookmark_csv": parsed_opts.export_bookmark_csv,
            "export_query_csv": parsed_opts.export_query_csv,
            "export_bookmark_db": parsed_opts.export_bookmark_db,
            "data_from_bookmark": parsed_opts.data_from_bookmark,
            "filter_bookmark": parsed_opts.filter_bookmark,
        }

        storage_opts = {
            "order_results": parsed_opts.order_results,
            "output_all_poses": parsed_opts.output_all_poses,
            "results_view_name": parsed_opts.results_view_name,
            "overwrite": parsed_opts.overwrite,
            "append_results": parsed_opts.append_results,
            "conflict_opt": conflict_handling,
            "mode": self.mode,
            "storage_type": "sqlite",
        }

        if isinstance(parsed_opts.interaction_tolerance, str):
            parsed_opts.interaction_tolerance = [
                float(val) for val in parsed_opts.interaction_cutoffs.split(",")
            ]
        rman_opts = {
            "chunk_size": 1,
            "mode": self.mode,
            "store_all_poses": parsed_opts.store_all_poses,
            "max_poses": parsed_opts.max_poses,
            "interaction_tolerance": parsed_opts.interaction_tolerance,
            "add_interactions": parsed_opts.add_interactions,
            "interaction_cutoffs": parsed_opts.interaction_cutoffs,
            "receptor_file": parsed_opts.receptor_file,
            "file_sources": file_sources,
            "file_pattern": self.pattern,
            "parser_manager": "multiprocessing",
            "max_proc": parsed_opts.max_proc,
        }

        # save target name
        if rman_opts["receptor_file"] is not None:
            receptor = (
                rman_opts["receptor_file"].split(".")[0].split("/")[-1]
            )  # remove file extension and path
        else:
            receptor = None
        rman_opts["target"] = receptor

        self.file_sources = file_sources
        self.input_db = parsed_opts.input_db
        if parsed_opts.input_db is not None:
            dbFile = parsed_opts.input_db
            if not os.path.exists(dbFile):
                raise OptionError("WARNING: input database does not exist!")
        else:
            dbFile = parsed_opts.output_db
        storage_opts["db_file"] = dbFile

        # make attributes for parsed opts
        self.storage_opts = storage_opts
        self.filters = filters
        self.out_opts = out_opts
        self.rman_opts = rman_opts
