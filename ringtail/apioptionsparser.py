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

class APIOptionParser:
    def __init__(self, ringtail_core, opts):

        self.rtcore = ringtail_core
        self.process_options(opts)
        self._initialize_options()

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

    def _initialize_options(self):
        RingtailCore.get_defaults() #I changed this, could be bad

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
            #TODO I won't need this if we use the objects and restructure the code
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

            #TODO This might be OK to keep here, or put in whatever code runs the main argument logic
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
                if getattr(parsed_opts, f) is not None: #TODO this just says if either of these values is not None, then there are filters available
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

        #TODO ok so here I think RT took cmd line inputs, then "erased" the read options if write was 
        # used, and then now the dict draws values from the opt object, where cmd line opts had been
        # erased at the appropriate places. So some of these values (read only) are sat twice basically
        self.rt_process_options = {
            "filter": filter_flag, #TODO what are you
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
        #TODO this for example is not needed if you are writing, right? 
        all_filters["react_any"] = parsed_opts.react_any
        for k,v in all_filters.items():
            setattr(self.rtcore.filters, k, v)

        # set storageman opts
        #TODO only one database treated at the time, so it writes only to bookmarks after write, it 
        # does not write back to the database? 
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
            receptor = None #TODO seems redundant

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
