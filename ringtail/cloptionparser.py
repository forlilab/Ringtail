#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail command line option parser
#

import sys
import argparse
import os
from .exceptions import OptionError
import __main__
from .logutils import get_logger, setup_logging

logger = get_logger(__name__)
from .ringtailoptions import (
    Filters,
    validate_docking_mode,
    ringtail_defaults,
    RingtailDefaults,
)
from .schema import OUTFIELD_SCHEMA, ORDER_RESULT_SCHEMA

_outfields_opts = "; ".join(
    f'"{name}" ({col.description})' for name, col in OUTFIELD_SCHEMA.items()
)
_outfields_help = (
    "Defines which fields are used when reporting results (to eg output csv or output_log). "
    "Comma-separated, e.g. --outfields=docking_score,leff,num_hb. "
    f"Available fields: {_outfields_opts}. "
    "Ligand name is always returned first."
)

_order_opts = "; ".join(
    f'"{name}" ({col.description})' for name, col in ORDER_RESULT_SCHEMA.items()
)
_order_help = (
    "Order results by this single field when e.g., writing to the output log file."
    "By default ordered by database insertion order. "
    f"Available fields: {_order_opts}."
)
from types import SimpleNamespace


def cmdline_parser(defaults: dict = {}):
    """Parses options provided using the command line.
    All arguments are first populated with default values.
    If a config file is provided, these will overwrite default values.
    Any single arguments provided using the argument parser will overwrite
    default and config file values.

    Args:
        defaults (dict): default argument values

    """

    conf_parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        add_help=False,
    )
    conf_parser.add_argument("-c", "--config")
    confargs, remaining_argv = conf_parser.parse_known_args()

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
        "-st",
        "--storage_type",
        help="specify a database file to perform actions with",
        action="store",
        type=str,
        metavar="sqlite",
    )
    write_parser.add_argument(
        "-m",
        "--docking_mode",
        help='specify AutoDock program used to generate results. Available options are "adgpu"/"dlg" and "vina"/"pdbqt". Vina mode will automatically change --pattern to *.pdbqt',
        action="store",
        type=str,
        metavar="'adgpu' or 'vina'",
    )
    write_parser.add_argument(
        "-su",
        "--print_summary",
        help="Prints summary information about stored data to STDOUT. Includes number of stored ligands and poses, min and max docking score and ligand efficiency, and 1%% (percentile) and 10%% (percentile) energy and ligand efficiency.",
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
        "--logfile",
        help="Write log output to this file (useful with --debug).",
        type=str,
        metavar="FILE.log",
        default=None,
    )
    write_parser.add_argument(
        "-dr",
        "--docking_results",
        help="docking result file(s), director(ies) to scan, or text file(s) listing result paths. "
             "Type is auto-detected per item. Replaces --file / --file_path / --file_list.",
        action="append",
        type=str,
        metavar="FILE_OR_DIR_OR_LIST",
        nargs="+",
    )
    write_parser.add_argument(
        "-f",
        "--file",
        help="(deprecated, use --docking_results) ligand docking output file to save.",
        action="append",
        type=str,
        metavar="FILENAME.[DLG/PDBQT][.gz]",
        nargs="+",
    )
    write_parser.add_argument(
        "-fp",
        "--file_path",
        help="(deprecated, use --docking_results) directory(s) containing docking output files to save.",
        action="append",
        type=str,
        metavar="DIRNAME",
        nargs="+",
    )
    write_parser.add_argument(
        "-fl",
        "--file_list",
        help="(deprecated, use --docking_results) file(s) containing the list of docking output files to save.",
        action="append",
        type=str,
        metavar="FILENAME",
        nargs="+",
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
        help="by default, if a database or log file exists, it doesn't get overwritten and an error is returned; this option enable overwriting existing database and/or log file.",
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
        "-ni",
        "--no_interactions",
        help="If interactions for vina results should not be calculated and stored",
        action="store_true",
    )
    write_parser.add_argument(
        "-ic",
        "--interaction_cutoffs",
        help="Specify distance cutoffs for measuring interactions between ligand and receptor in angstroms if other than default. Give as string, separating cutoffs for hydrogen bonds and VDW with comma (in that order). E.g. '-ic 3.7,4.0' will set the cutoff for hydrogen bonds to 3.7 angstroms and for VDW to 4.0. These are the default cutoffs.",
        action="store",
        type=str,
        metavar="[HB CUTOFF],[VDW CUTOFF]",
    )
    write_parser.add_argument(
        "-rf",
        "--receptor_file",
        help="Give file for receptor .pdbqt OR .json file for a receptor meeko.Polymer.",
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
        help="Specify name for db view of passing results to create or export from. Can only contain lower case letters, numbers, and underscore (_)",
        action="store",
        type=str,
        metavar="STRING",
    )
    read_parser.add_argument(
        "-su",
        "--print_summary",
        help="Prints summary information about stored data to STDOUT. Includes number of stored ligands and poses, min and max docking score and ligand efficiency, and 1%% (percentile) and 10%% (percentile) energy and ligand efficiency.",
        action="store_true",
    )
    read_parser.add_argument(
        "-pb",
        "--print_bookmarks",
        help="Prints all bookmark names stored in the database to STDOUT.",
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
    read_parser.add_argument(
        "--logfile",
        help="Write log output to this file (useful with --debug).",
        type=str,
        metavar="FILE.log",
        default=None,
    )

    output_group = read_parser.add_argument_group("Output options")
    output_group.add_argument(
        "-l",
        "--output_log",
        help="if this option is used, ligands and requested info passing the filters will be written to specified file. Can be used as a flag, which will write to output_log.txt.",
        action="store",
        type=str,
        const="output_log.txt",
        nargs="?",
        metavar="[FILE_NAME].TXT",
    )
    output_group.add_argument(
        "-of",
        "--outfields",
        help=_outfields_help,
        action="store",
        type=str,
        metavar="FIELD1,FIELD2,...",
    )
    output_group.add_argument(
        "-ord",
        "--order_results",
        help=_order_help,
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
        help="Create csv of the bookmark given with bookmark_name. Output as <filename>.csv. Can also export full database tables",
        action="store",
        type=str,
        metavar="<filename>.csv",
        const=True,
        default=None,
        nargs="?",
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
        help="specify the path where to save poses of ligands passing the filters (SDF format); if the directory does not exist, it will be created; if it already exist, it will throw an error, unless the --overwrite is used  NOTE: the log file will be automatically saved in this path. Ligands will be stored as one large SDF file unless using '--individual_sdf_files'.",
        action="store",
        type=str,
        metavar="DIRECTORY_NAME",
    )
    output_group.add_argument(
        "--individual_sdf_files",
        help="Use if you like to print chosen molecules to individual SDF files, as opposed to one big SDF.",
        action="store_true",
    )
    output_group.add_argument(
        "-xdb",
        "--export_bookmark_db",
        help="Export a database containing only the results found in the bookmark specified by --bookmark_name. Will save as <input_db>_<bookmark_name>.db",
        action="store_true",
    )
    output_group.add_argument(
        "-xr",
        "--export_receptor_pdbqt",
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
        "--input_bookmark",
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
        "--ligand_name",
        help="specify ligand name(s). Will combine name filters with OR",
        action="append",
        type=str,
        metavar="STRING",
        nargs="+",
    )
    ligand_group.add_argument(
        "--ligand_name_file",
        help=".csv file containing list of ligand names to be used as a filter",
        action="store",
        type=str,
        metavar="STRING",
    )
    ligand_group.add_argument(
        "-mna",
        "--ligand_max_atoms",
        help="Maximum number of heavy atoms (non-hydrogens) a ligand may have",
        action="store",
        type=int,
        metavar="INT",
    )
    ligand_group.add_argument(
        "-lmw",
        "--ligand_min_molweight",
        help="min molweight (inclusive, g/mol) for the ligand",
        action="store",
        type=float,
        metavar="float",
    )
    ligand_group.add_argument(
        "-hmw",
        "--ligand_max_molweight",
        help="max molweight (inclusive, g/mol) for the ligand",
        action="store",
        type=float,
        metavar="float",
    )
    ligand_group.add_argument(
        "--ligand_substruct",
        help="SMARTS pattern(s) for substructure matching, if error delimit each substructure with ''.",
        action="append",
        type=str,
        metavar="STRING",
        nargs="+",
    )
    ligand_group.add_argument(
        "--ligand_substruct_pos",
        help='"SMARTS, index of atom in SMARTS, cutoff dist, and target XYZ coords". Group each set of six values with "".',
        action="append",
        type=str,
        metavar="STRING",
        nargs="+",
    )
    ligand_group.add_argument(
        "-sj",
        "--ligand_operator",
        choices=["AND", "OR"],
        help="logical join operator for multiple SMARTS (default: OR)",
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
        "--vdw_interactions",
        help="define van der Waals interactions with residue",
        action="append",
        type=str,
        metavar="[-][CHAIN]:[RES]:[NUM]:[ATOM_NAME]",
    )
    interaction_group.add_argument(
        "-hb",
        "--hb_interactions",
        help="define HB (ligand acceptor or donor) interaction",
        action="append",
        type=str,
        metavar="[-][CHAIN]:[RES]:[NUM]:[ATOM_NAME]",
    )
    interaction_group.add_argument(
        "-r",
        "--reactive_interactions",
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
        raise OptionError(
            "Script called with no commandline options. Please call with either 'read' or 'write'. See --help for details."
        )

    if confargs.config is not None:
        # update default values with those given in the config file
        try:
            with open(confargs.config) as json_file:
                import json

                config_dict = json.load(json_file)
                defaults.update(config_dict)
        except FileNotFoundError as e:
            raise OptionError("Config file not found in current directory.") from e

    cmd_line_prompt = str(sys.argv)
    parser.set_defaults(**defaults)
    write_parser.set_defaults(**defaults)
    read_parser.set_defaults(**defaults)
    args = parser.parse_args(remaining_argv)
    return args, parser, confargs, write_parser, read_parser, cmd_line_prompt


class CLOptionParser:
    """Command line option/argument parser. Options and switches are utilized in the script
    'rt_process_vs.py'.

    Attributes:
        process_mode (str): operating in 'write' or 'read' mode
        db_file
        logging_level
        storage_type
        file_sources (dict): paths to docking results and receptor files
        print_summary (bool): switch to print database summary
        filtering (bool): switch to run filtering method
        in write mode:
            write_options (SimpleNamespace)
        in read mode:
            filters (dict): parsed and organized filters
            filter_options (SimpleNamespace)
            output_options (SimpleNamespace)

    Raises:
        OptionError: Error when an option cannot be parsed correctly
    """

    def __init__(self):
        # create parser
        parsed_opts = None
        try:
            # add default values from ringtailoptions
            defaults_dict = ringtail_defaults()

            (
                parsed_opts,
                self.parser,
                self.confargs,
                self.write_parser,
                self.read_parser,
                self.cmd_line_prompt,
            ) = cmdline_parser(defaults_dict)
            self.process_options(parsed_opts)

        except argparse.ArgumentError as e:
            raise OptionError(
                "Invalid option or option ordering. Be sure to put read/write mode before any other arguments"
            ) from e
        except Exception as e:
            try:
                if parsed_opts is not None and parsed_opts.process_mode == "write":
                    self.write_parser.print_help()
                elif parsed_opts is not None and parsed_opts.process_mode == "read":
                    self.read_parser.print_help()
            finally:
                raise e

    def process_options(self, parsed_opts):
        """
        Process and organize command line options to into ringtail options and filter dictionaries and ringtail core attributes

        Args:
            parsed_opts (argparse.Namespace): arguments provided through the cmdline_parser method.

        """
        if parsed_opts.debug:
            log_level = "DEBUG"
        elif parsed_opts.verbose:
            log_level = "INFO"
        else:
            log_level = "WARNING"
        setup_logging(level=log_level, logfile=getattr(parsed_opts, "logfile", None))

        # make sure we log the command line prompt
        logger.info("Command line prompt: " + self.cmd_line_prompt)

        if parsed_opts.process_mode is None:
            raise OptionError(
                "No mode specified for rt_process_vs.py. Please specify mode (write/read)."
            )
        self.process_mode = parsed_opts.process_mode.lower()
        # initialize flag indicating if any filters given
        self.filtering = False
        self.clustering_only = False
        # Check database options
        if parsed_opts.input_db is not None:
            if not os.path.exists(parsed_opts.input_db):
                raise OptionError(
                    f"WARNING: input database -{parsed_opts.input_db}- does not exist or cannot be found!"
                )
            db_file = parsed_opts.input_db
        elif self.process_mode == "read" and parsed_opts.input_db is None:
            raise OptionError(
                "No input database specified in read mode. Please specify database with --input_db"
            )
        else:
            db_file = parsed_opts.output_db

        docking_mode = validate_docking_mode(parsed_opts.docking_mode)

        # Ringtail Core attributes
        self.db_file = db_file
        self.storage_type = (
            parsed_opts.storage_type.lower()
            if parsed_opts.storage_type
            else RingtailDefaults.storage_type
        )

        self.logging_level = log_level

        self.print_summary = parsed_opts.print_summary
        self.print_bookmarks = (
            parsed_opts.print_bookmarks if self.process_mode == "read" else False
        )

        if self.process_mode == "write":
            # Check if writing to an existing database
            only_adding_receptor = (
                parsed_opts.receptor_file
                and not parsed_opts.docking_results
                and not parsed_opts.file
                and not parsed_opts.file_path
                and not parsed_opts.file_list
            )
            if (
                os.path.exists(db_file)
                and not parsed_opts.append_results
                and not parsed_opts.overwrite
                and not only_adding_receptor
            ):
                raise OptionError(
                    f"The database {db_file} exists but the user has not specified to --append_results or --overwrite. Please include one of these options if writing to an existing database."
                )

            import itertools

            def _flatten(arg):
                return list(itertools.chain.from_iterable(arg)) if arg else []

            all_inputs = (
                _flatten(parsed_opts.docking_results)
                + _flatten(parsed_opts.file)
                + _flatten(parsed_opts.file_path)
                + _flatten(parsed_opts.file_list)
            )
            self.file_sources = {
                "docking_results": all_inputs or None,
                "recursive": parsed_opts.recursive,
                "receptor_file": parsed_opts.receptor_file,
                "save_receptor": parsed_opts.save_receptor,
            }

            if isinstance(parsed_opts.interaction_tolerance, str):
                parsed_opts.interaction_tolerance = [
                    float(val) for val in parsed_opts.interaction_cutoffs.split(",")
                ]

            self.write_options = SimpleNamespace(
                docking_mode=docking_mode,
                duplicate_handling=parsed_opts.duplicate_handling,
                overwrite=parsed_opts.overwrite,
                store_all_poses=parsed_opts.store_all_poses,
                max_poses=parsed_opts.max_poses,
                calculate_interactions=not parsed_opts.no_interactions,
                interaction_tolerance=parsed_opts.interaction_tolerance,
                interaction_cutoffs=parsed_opts.interaction_cutoffs,
                max_proc=parsed_opts.max_proc,
            )

        elif self.process_mode == "read":
            # set up filters
            filters = {}
            # empty lists are counted as filters, and the default keyword "OR" too
            optional_filters = Filters.get_filter_keys("all")
            optional_filters.append(
                "hb_count"
            )  # hb/interaction count is not part of the three formal list but works as a filter on its own
            optional_filters.append(
                "react_any"
            )  # react_any is not part of the three formal list but works as a filter on its own
            for f in optional_filters:
                if (
                    getattr(parsed_opts, f) is not None
                    and getattr(parsed_opts, f) != []
                ):
                    if f == "ligand_operator":
                        pass
                    else:
                        self.filtering = True

            if self.filtering:
                # property filters
                property_list = Filters.get_filter_keys("property")
                for kw in property_list:
                    filters[kw] = getattr(parsed_opts, kw)

                # interaction filters (residues)
                interactions = {}
                interactions_kw = Filters.get_filter_keys("interaction")
                for _type in interactions_kw:
                    interactions[_type] = []
                    res_list = getattr(parsed_opts, _type)
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
                        if type(res) == str:
                            logger.debug(
                                "cloptionparser: interaction filters provided as string"
                            )
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
                        else:
                            logger.debug(
                                "cloptionparser: interaction filters provided as list"
                            )
                            if (
                                not res[0].count(":") == 3
                            ):  # first element of list is the interaction
                                raise OptionError(
                                    (
                                        "[%s]: to specify a residue use "
                                        "the format CHAIN:RES:NUM:ATOM_NAME. Any item can be omitted, "
                                        "as long as the number of semicolons is always 3 "
                                        "(e.g.: CHAIN:::, :RES::, CHAIN::NUM:, etc.)"
                                    )
                                    % res
                                )
                            interactions[_type].append(res)

                for k, v in interactions.items():
                    filters[k] = v

                # count interactions
                hb_count = []
                count_kw = [("hb_count", ("hb_count"))]
                for kw, pool in count_kw:
                    c = getattr(parsed_opts, kw, None)
                    if c is None:
                        continue
                    hb_count.append((pool, c))
                filters["hb_count"] = hb_count

                # make dictionary for ligand filters
                ligand_kw = Filters.get_filter_keys("ligand")
                ligand_filters = {}
                # can only use one type of ligand name specification
                if (
                    parsed_opts.ligand_name is not None
                    and parsed_opts.ligand_name_file is not None
                ):
                    raise OptionError(
                        "Cannot use --ligand_name and --ligand_name_file together, please choose one."
                    )
                # parse the ligand filters, depending on how the keywords are used they will be a list of list or list of lists
                for _type in ligand_kw:
                    ligand_filter_value = getattr(parsed_opts, _type)
                    # just a simple string
                    if _type in [
                        "ligand_max_atoms",
                        "ligand_min_molweight",
                        "ligand_max_molweight",
                    ]:
                        ligand_filters[_type] = ligand_filter_value
                        continue
                    # don't include None values
                    if ligand_filter_value is (None):
                        continue
                    ligand_filters[_type] = []
                    # the other ligand filters can come as [[filter1,filter2,filter3]] or [[filter1],[filter2, filter3]]
                    for filter_list in ligand_filter_value:
                        # if more than one filter in list, go through each
                        if len(filter_list) > 1:
                            for filter in filter_list:
                                if _type == "ligand_substruct_pos":
                                    # make a lits of the six values
                                    ligand_filters[_type].append(
                                        [i for i in filter.split(" ")]
                                    )
                                else:
                                    ligand_filters[_type].append(filter)
                        else:
                            if _type == "ligand_substruct_pos":
                                # if only one item in list, append to ligand list
                                ligand_filters[_type].append(
                                    [i for i in filter_list[0].split(" ")]
                                )
                            else:
                                ligand_filters[_type].append(filter_list[0])

                ligand_filters["ligand_operator"] = parsed_opts.ligand_operator
                # ligand substruct pos needs six items
                if (
                    "ligand_substruct_pos" in ligand_filters
                    and len(ligand_filters["ligand_substruct_pos"][0]) % 6
                ):
                    msg = "--ligand_substruct_pos needs groups of 6 values:\n"
                    msg += "  1. Ligand SMARTS\n"
                    msg += "  2. index of atom in ligand SMARTS (0 based)\n"
                    msg += "  3. distance cutoff\n"
                    msg += "  4. X\n"
                    msg += "  5. Y\n"
                    msg += "  6. Z\n"
                    msg += (
                        'For example --ligand_substruct_pos "[C][Oh] 1 1.5 -20 42 -7.1"'
                    )
                    raise OptionError(msg)

                for k, v in ligand_filters.items():
                    filters[k] = v
                filters["max_miss"] = parsed_opts.max_miss
                filters["react_any"] = parsed_opts.react_any

            self.filters = filters

            # another form of filtering is clustering, which can in theory be performed without filtering (such as clustering over a filtered bookmark)
            if (
                parsed_opts.mfpt_cluster or parsed_opts.interaction_cluster
            ) and not self.filtering:
                self.clustering_only = True

            self.filter_options = SimpleNamespace(
                output_bookmark=parsed_opts.bookmark_name,
                input_bookmark=parsed_opts.input_bookmark,
                enumerate_interaction_combs=parsed_opts.enumerate_interaction_combs,
                mfpt_cluster=parsed_opts.mfpt_cluster,
                interaction_cluster=parsed_opts.interaction_cluster,
                ligand_name_file=parsed_opts.ligand_name_file,
                output_log=parsed_opts.output_log,
                outfields=parsed_opts.outfields,
                order_results=parsed_opts.order_results,
                output_all_poses=parsed_opts.output_all_poses,
            )

            self.output_options = SimpleNamespace(
                export_bookmark_db=parsed_opts.export_bookmark_db,
                export_receptor_pdbqt=parsed_opts.export_receptor_pdbqt,
                data_from_bookmark=parsed_opts.data_from_bookmark,
                individual_sdf_files=parsed_opts.individual_sdf_files,
                output_log=parsed_opts.output_log,
                export_sdf_path=parsed_opts.export_sdf_path,
                export_query_csv=parsed_opts.export_query_csv,
                find_similar_ligands=parsed_opts.find_similar_ligands,
                export_bookmark_csv=parsed_opts.export_bookmark_csv,
                outfields=parsed_opts.outfields,
                order_results=parsed_opts.order_results,
                output_all_poses=parsed_opts.output_all_poses,
            )
