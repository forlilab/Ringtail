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
    config = {}
    for section in defaults:
        for k, v in defaults[section]["values"].items():
            config[k] = v
    # separate filters into separate arguments
    replace_filter_keys = {
        "vdw_interactions": "van_der_waals",
        "hb_interactions": "hydrogen_bond",
        "reactive_interactions": "reactive_res",
        "ligand_name": "name",
        "ligand_substruct": "smarts",
        "ligand_operator": "smarts_join",
        "ligand_substruct_pos": "smarts_idxyz",
        "ligand_max_atoms": "max_nr_atoms",
    }


    for k, v in replace_filter_keys.items():
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

        CITATION
                If using Ringtail in your work, please cite the following publication:\n

                Ringtail: A Python Tool for Efficient Management and Storage of Virtual Screening Results
                Althea T. Hansel-Harris, Diogo Santos-Martins, Niccolò Bruciaferri, Andreas F. Tillack, Matthew Holcomb, and Stefano Forli
                Journal of Chemical Information and Modeling 2023 63 (7), 1858-1864
                DOI: 10.1021/acs.jcim.3c00166

        REPORTING BUGS
                Please report bugs to:
                AutoDock mailing list   http://autodock.scripps.edu/mailing_list\n

        COPYRIGHT
                Copyright (C) 2023 Forli Lab, Center for Computational Structural Biology,
                             The Scripps Research Institute.
                GNU L-GPL version 2.1 or later <https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html> \n
        """,
        exit_on_error=False,
    )

    subparsers = parser.add_subparsers(
        help=f"Specify if should write to or read from database. To show options of each mode use the in-line help, e.g.: {os.path.basename(__main__.__file__)} read -h",
        dest="process_mode",
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
        "-su",
        "--summary",
        help='Prints summary information about stored data to STDOUT. Includes number of stored ligands and poses, min and max docking score and ligand efficiency, and 1%% (percentile) and 10%% (percentile) energy and ligand efficiency.',
        action="store_true",
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
        "-su",
        "--summary",
        help='Prints summary information about stored data to STDOUT. Includes number of stored ligands and poses, min and max docking score and ligand efficiency, and 1%% (percentile) and 10%% (percentile) energy and ligand efficiency.',
        action="store_true",
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
        "--log_file",
        help='by default, results are saved in "output_log.txt"; if this option is used, ligands and requested info passing the filters will be written to specified file',
        action="store",
        type=str,
        metavar="[FILE_NAME].TXT",
    )
    output_group.add_argument(
        "-of",
        "--outfields",
        help=(
            'defines which fields are used when reporting the results (to stdout and to the log file); fields are specified as comma-separated values, e.g. "--outfields=e,le,hb"; by default, docking_score (energy) and ligand name are reported; ligand always reported in first column available fields are:  '
            '"Ligand_name" (Ligand name), '
            '"e" (docking_score), '
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
            "printed in the order in which they are provided. Ligand name will always be returned and will be added in first position if not specified."
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
        '"e" (docking_score), '
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
        "-mfpc",
        "--mfpt_cluster",
        help="Cluster filered ligands by Tanimoto distance of Morgan fingerprints with Butina clustering and output ligand with lowest ligand efficiency from each cluster. Default clustering cutoff is 0.5. Useful for selecting chemically dissimilar ligands.",
        action="store",
        type=float,
        metavar="FLOAT",
        const=0.5,
        nargs="?",
    )
    output_group.add_argument(
        "-ifpc",
        "--interaction_cluster",
        help="Cluster filered ligands by Tanimoto distance of interaction fingerprints with Butina clustering and output ligand with lowest ligand efficiency from each cluster. Default clustering cutoff is 0.5. Useful for enhancing selection of ligands with diverse interactions.",
        action="store",
        type=float,
        metavar="FLOAT",
        const=0.5,
        nargs="?",
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
        "-xr",
        "--export_receptor",
        help="Export stored receptor pdbqt. Will write to current directory.",
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
        "-fsl",
        "--find_similar_ligands",
        help="Allows user to find similar ligands to given ligand name based on previously performed morgan fingerprint or interaction clustering.",
        action="store",
        type=str,
    )
    output_group.add_argument(
        "-p",
        "--plot",
        help="Makes scatterplot of LE vs Best Energy, saves as scatter.png.",
        action="store_true",
    )
    output_group.add_argument(
        '-py',
        '--pymol',
        help="Lauch PyMOL session and plot of ligand efficiency vs docking score for molecules in bookmark specified with --bookmark_name. Will display molecule in PyMOL when clicked on plot. Will also open receptor if given.",
        action="store_true"
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
        "--score_percentile",
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
    ligand_group.add_argument(
        "-mna",
        "--max_nr_atoms",
        help="Maximum number of heavy atoms a ligand may have",
        action="store",
        type=int,
        metavar="INT",
    )
    ligand_group.add_argument(
        "--smarts",
        help="SMARTS pattern(s) for substructure matching",
        action="store",
        type=str,
        metavar="STRING",
        nargs="+",
    )
    ligand_group.add_argument(
        "--smarts_idxyz",
        help="SMARTS, index of atom in SMARTS, cutoff dist, and target XYZ coords",
        action="store",
        type=str,
        metavar="STRING",
        nargs="+",
    )
    ligand_group.add_argument(
        "-sj",
        "--smarts_join",
        choices=["AND", "OR"],
        help="logical operator for multiple SMARTS (default: OR)",
        action="store",
        type=str,
        metavar="STRING",
    )

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
        help="Will compute all possible combinations of interaction filters excluding up to max_miss numer of interactions from given set. Default will only return union of poses interaction filter combinations. Use with --enumerate_interaction_combs for enumeration of poses passing each individual combination of interaction filters.",
        action="store",
        type=int,
        metavar="INTEGER",
    )
    interaction_group.add_argument(
        "-eic",
        "--enumerate_interaction_combs",
        help="Use with max_miss. If used, will output ligands passing each individual combination of interaction filters with max_miss.",
        action="store_true",
    )

    # catch if running with no options
    if len(sys.argv) == 1:
        parser.print_help()
        raise OptionError("Script called with no commandline options. Please call with either 'read' or 'write'. See --help for details.")

    parser.set_defaults(**config)
    write_parser.set_defaults(**config)
    read_parser.set_defaults(**config)
    args = parser.parse_args(remaining_argv)

    return args, parser, conf_parser, write_parser, read_parser


class CLOptionParser:
    def __init__(self, ringtail_core):

        self.rtcore = ringtail_core
        self._initialize_parser()

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
            logging.error("\n")
            raise OptionError(
                "Invalid option or option ordering. Be sure to put read/write mode before any other arguments"
            ) from e
        except Exception as e:
            try:
                if self.process_mode == "write":
                    self.write_parser.print_help()
                elif self.process_mode == "read":
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

        self.process_mode = parsed_opts.process_mode

        filter_flag = False  # set flag indicating if any filters given
        conflict_handling = None
        file_sources = None
        self.rec_files_pool = []
        all_filters = {}
        storage_opts = {}
        parsed_opts.mode = parsed_opts.mode.lower()

        if parsed_opts.mode == "vina":
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
            parsed_opts.score_percentile is not None
            or parsed_opts.le_percentile is not None
        ):
            logging.warning(
                "Cannot return all passing poses with percentile filter. Will only log best pose."
            )
            parsed_opts.output_all_poses = False

        if self.process_mode is None:
            raise OptionError(
                "No mode specified for rt_process_vs.py. Please specify mode (write/read)."
            )
        if self.process_mode == "write":
            # initialize rt_process read-only options to prevent errors
            parsed_opts.plot = None
            parsed_opts.find_similar_ligands=None
            parsed_opts.export_bookmark_csv = None
            parsed_opts.export_query_csv = None
            parsed_opts.export_bookmark_db = None
            parsed_opts.export_receptor = None
            parsed_opts.data_from_bookmark = None
            parsed_opts.pymol = None
            parsed_opts.log_file = None
            parsed_opts.enumerate_interaction_combs = None
            # confirm that receptor file was found if needed, else throw error
            # if only receptor files found and --save_receptor, assume we just want to
            # add receptor and not modify the rest of the db
            # set receptor name, add_interactions if given
            if parsed_opts.save_receptor:
                if parsed_opts.append_results:
                    raise OptionError(
                        "Cannot use --append_results with --save_receptor. Please remove the --save_receptor flag"
                    )
                if parsed_opts.receptor_file is None:
                    raise OptionError(
                        "Must provide path for receptor PDBQT file is using --save_receptor"
                    )

            # raise error if --save_receptor or --add_interactions and receptor not found
            if (parsed_opts.save_receptor or parsed_opts.add_interactions) and not os.path.exists(
                parsed_opts.receptor_file
            ):
                raise OptionError(
                    "--save_receptor or --add_interaction flag specified but no receptor PDBQT found. Please check location of receptor file and --receptor_file option"
                )
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
                and (file_sources["file_path"] == {'path': [[]], 'pattern': '*.dlg*', 'recursive': None})
                and (file_sources["file_list"] == [[]])
                and (parsed_opts.input_db is None)
            ):
                raise OptionError(
                    "At least one input option needs to be used:  --file, --file_path, --file_list, --input_db"
                )
            if (
                (file_sources["file"] == [[]])
                and (file_sources["file_path"] == {'path': [[]], 'pattern': '*.dlg*', 'recursive': None})
                and (file_sources["file_list"] == [[]])
                and parsed_opts.save_receptor
            ):
                self.process_mode = "add_receptor"

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
        elif self.process_mode == "read":
            # initialize write-only rt_process options to prevent errors
            parsed_opts.receptor_file = None
            parsed_opts.save_receptor = None

            # option compatibility checking
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
                if parsed_opts.export_sdf_path:
                    logging.warning("WARNING: Requested --export_sdf_path with --max_miss. Exported SDFs will be for union of interaction combinations.")
            if parsed_opts.order_results:
                if parsed_opts.order_results not in self.order_options:
                    raise OptionError(
                        "Requested ording option that is not available. Please see --help for available options."
                    )
            # parse output options
            # Make sure that export_sdf_path has trailing /, is directory
            if parsed_opts.export_sdf_path:
                if not parsed_opts.export_sdf_path.endswith("/"):
                    parsed_opts.export_sdf_path += "/"
                if not os.path.isdir(parsed_opts.export_sdf_path):
                    raise OptionError(
                        "--export_sdf_path directory does not exist. Please create directory first"
                    )

            # # # filters
            optional_filters = ["eworst", "ebest", "leworst", "lebest", "score_percentile", "le_percentile", "van_der_waals", "hydrogen_bond", "reactive_res", "name", "smarts", "smarts_idxyz", "max_nr_atoms"]
            for f in optional_filters:
                if getattr(parsed_opts, f) is not None:
                    filter_flag = True
            # Cannot use energy/le cuttoffs with percentiles. Override percentile with given cutoff
            if (
                parsed_opts.eworst is not None
                and parsed_opts.score_percentile is not None
            ):
                logging.warning(
                    "Cannot use --eworst cutoff with --score_percentile. Overiding score_percentile with eworst."
                )
                parsed_opts.score_percentile = None
            if (
                parsed_opts.leworst is not None
                and parsed_opts.le_percentile is not None
            ):
                logging.warning(
                    "Cannot use --leworst cutoff with --le_percentile. Overiding le_percentile with leworst."
                )
                parsed_opts.le_percentile = None
            if (
                parsed_opts.score_percentile is not None
                or parsed_opts.le_percentile is not None
            ) and parsed_opts.filter_bookmark is not None:
                raise OptionError(
                    "Cannot use --score_percentile or --le_percentile with --filter_bookmark."
                )
            # property filters
            property_filters = {
                "eworst": None,
                "ebest": None,
                "leworst": None,
                "lebest": None,
                "score_percentile": None,
                "le_percentile": None,
            }
            for kw, _ in property_filters.items():
                property_filters[kw] = getattr(parsed_opts, kw)
            # interaction filters (residues)
            interactions = {}
            res_interactions_kw = [
                ("van_der_waals", "vdw_interactions"),
                ("hydrogen_bond", "hb_interactions"),
                ("reactive_res", "reactive_interactions"),
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
            # count interactions
            interactions_count = []
            count_kw = [("hb_count", ("hb_count"))]
            for kw, pool in count_kw:
                c = getattr(parsed_opts, kw, None)
                if c is None:
                    continue
                interactions_count.append((pool, c))
            # make dictionary for ligand filters
            ligand_filters_kw = [("name", "ligand_name"), ("smarts", "ligand_substruct"), ("smarts_idxyz", "ligand_substruct_pos"), ("max_nr_atoms", "ligand_max_atoms")]
            ligand_filters = {}
            ligand_filter_list = []
            for kw, _type in ligand_filters_kw:
                ligand_filter_list = getattr(parsed_opts, kw)
                if kw == "max_nr_atoms":
                    ligand_filters[_type] = ligand_filter_list
                    continue
                ligand_filters[_type] = []
                if ligand_filter_list is None:
                    continue
                for fil in ligand_filter_list:
                    ligand_filters[_type].append(fil)
            if ligand_filters["ligand_max_atoms"] is not None and len(ligand_filters["ligand_max_atoms"]) % 6 != 0:
                msg = "--smarts_idxyz needs groups of 6 values:\n"
                msg += "  1. SMARTS\n"
                msg += "  2. index of atom in SMARTS (0 based)\n"
                msg += "  3. distance cutoff\n"
                msg += "  4. X\n"
                msg += "  5. Y\n"
                msg += "  6. Z\n"
                msg += "For example --smarts_idxyz \"[C][Oh]\" 1 1.5 -20. 42. -7.1"
                raise OptionError(msg)
            ligand_filters["ligand_operator"] = getattr(parsed_opts, "smarts_join")

            all_filters = property_filters | ligand_filters | interactions
            all_filters["max_miss"] = parsed_opts.max_miss

        # set all object options for both read and write mode

        self.rt_process_options = {
            "filter": filter_flag,
            "verbose": parsed_opts.verbose,  # both modes
            "debug": parsed_opts.debug,  # both modes
            "summary": parsed_opts.summary,  # both modes
            "plot": parsed_opts.plot,
            "find_similar_ligands": parsed_opts.find_similar_ligands,
            "export_bookmark_csv": parsed_opts.export_bookmark_csv,
            "export_query_csv": parsed_opts.export_query_csv,
            "export_bookmark_db": parsed_opts.export_bookmark_db,
            "export_receptor": parsed_opts.export_receptor,
            "data_from_bookmark": parsed_opts.data_from_bookmark,
            "pymol": parsed_opts.pymol,
            "enumerate_interaction_combs": parsed_opts.enumerate_interaction_combs,
            "export_sdf_path": parsed_opts.export_sdf_path,
            "receptor_file": parsed_opts.receptor_file,  # write only
            "save_receptor": parsed_opts.save_receptor,  # write only
        }

        # set core filter attributes
        all_filters["react_any"] = parsed_opts.react_any
        for k,v in all_filters.items():
            setattr(self.rtcore.filters, k, v)

        # set storageman opts
        if parsed_opts.input_db is not None:
            dbFile = parsed_opts.input_db
            if not os.path.exists(dbFile):
                raise OptionError("WARNING: input database does not exist!")
        else:
            dbFile = parsed_opts.output_db
        storage_opts = {
            "db_file": dbFile,
            "order_results": parsed_opts.order_results,
            "outfields": parsed_opts.outfields,
            "filter_bookmark": parsed_opts.filter_bookmark,
            "output_all_poses": parsed_opts.output_all_poses,
            "mfpt_cluster": parsed_opts.mfpt_cluster,
            "interaction_cluster": parsed_opts.interaction_cluster,
            "results_view_name": parsed_opts.results_view_name,
            "overwrite": parsed_opts.overwrite,
            "append_results": parsed_opts.append_results,
            "conflict_opt": conflict_handling,
        }
        for k,v in storage_opts.items():
            setattr(self.rtcore.storageman, k, v)

        # prepare and set results_man opts
        if isinstance(parsed_opts.interaction_tolerance, str):
            parsed_opts.interaction_tolerance = [
                float(val) for val in parsed_opts.interaction_cutoffs.split(",")
            ]
         # save target name
        if parsed_opts.receptor_file is not None:
            receptor = (
                os.path.basename(parsed_opts.receptor_file).split(".")[0]
            )  # remove file extension and path
        else:
            receptor = None
        rman_opts = {
            "chunk_size": 1,
            "mode": parsed_opts.mode,
            "store_all_poses": parsed_opts.store_all_poses,
            "max_poses": parsed_opts.max_poses,
            "interaction_tolerance": parsed_opts.interaction_tolerance,
            "add_interactions": parsed_opts.add_interactions,
            "interaction_cutoffs": parsed_opts.interaction_cutoffs,
            "receptor_file": parsed_opts.receptor_file,
            "target": receptor,
            "file_sources": file_sources,
            "file_pattern": parsed_opts.file_pattern,
            "parser_manager": "multiprocessing",
            "max_proc": parsed_opts.max_proc,
        }
        for k,v in rman_opts.items():
            setattr(self.rtcore.results_man, k, v)

        # set outputman options
        outman_opts = {
            "log_file": parsed_opts.log_file,
            "export_sdf_path": parsed_opts.export_sdf_path,
        }
        for k,v in outman_opts.items():
            setattr(self.rtcore.output_manager, k, v)
