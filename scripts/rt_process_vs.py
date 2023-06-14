#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#

import sys
import time
from ringtail import CLOptionParser
from ringtail import RingtailCore
import logging
import traceback

if __name__ == "__main__":

    time0 = time.perf_counter()
    level = logging.INFO
    logging.basicConfig(
        level=level, stream=sys.stdout, filemode="w", format="%(message)s"
    )

    try:
        rt_core = RingtailCore()
        cl_opts = CLOptionParser(rt_core)
    except Exception as e:
        tb = traceback.format_exc()
        logging.debug(tb)
        logging.critical("ERROR: " + str(e))
        sys.exit(1)

    # set logging level
    levels = {10: "debug", 20: "info", 30: "warning"}
    debug = cl_opts.rt_process_options["debug"]
    verbose = cl_opts.rt_process_options["verbose"]
    if debug:
        level = logging.DEBUG
    elif verbose:
        level = logging.INFO
    else:
        level = logging.WARNING
    logging.getLogger().setLevel(level)
    logging.info(f"Logging level set to {levels[level]}")

    # create manager object for virtual screening. Will make database if needed
    try:
        defaults = RingtailCore.get_defaults()
        with rt_core:

                # parse command line options and filters file (if given)

            if cl_opts.process_mode == "write":
                rt_core.add_results()

            # Add receptors to database if requested
            if cl_opts.rt_process_options["save_receptor"]:
                rt_core.save_receptors(cl_opts.rt_process_options["receptor_file"])

            # print summary
            if cl_opts.rt_process_options["summary"]:
                rt_core.produce_summary()

            time1 = time.perf_counter()

            if cl_opts.process_mode == "read":

                # perform filtering
                if cl_opts.rt_process_options["filter"]:
                    rt_core.filter(cl_opts.rt_process_options["enumerate_interaction_combs"])

                # Write log with new data for previous filtering results
                if cl_opts.rt_process_options["data_from_bookmark"] and not cl_opts.rt_process_options["filter"]:
                    rt_core.get_previous_filter_data()
                
                if cl_opts.rt_process_options["find_similar_ligands"]:
                    rt_core.find_similar_ligands(cl_opts.rt_process_options["find_similar_ligands"])

                # plot if requested
                if cl_opts.rt_process_options["plot"]:
                    rt_core.plot()

                # def open pymol viewer
                if cl_opts.rt_process_options["pymol"]:
                    rt_core.display_pymol()

                # write out molecules if requested
                if cl_opts.rt_process_options["export_sdf_path"]:
                    rt_core.write_molecule_sdfs()

                # write out requested CSVs
                if cl_opts.rt_process_options["export_bookmark_csv"]:
                    rt_core.export_csv(
                        cl_opts.rt_process_options["export_bookmark_csv"],
                        cl_opts.rt_process_options["export_bookmark_csv"] + ".csv",
                        table=True,
                    )
                if cl_opts.rt_process_options["export_query_csv"]:
                    rt_core.export_csv(cl_opts.rt_process_options["export_query_csv"], "query.csv")
                if cl_opts.rt_process_options["export_bookmark_db"]:
                    bookmark_name = (
                        rt_core.storageman.db_file.rstrip(".db")
                        + "_"
                        + rt_core.storageman.results_view_name
                        + ".db"
                    )
                    rt_core.export_bookmark_db(bookmark_name)

                if cl_opts.rt_process_options["export_receptor"]:
                    rt_core.export_receptors()

    except Exception as e:
        tb = traceback.format_exc()
        logging.debug(tb)
        logging.critical("ERROR: " + str(e))
        sys.exit(1)

    # print performance times
    time2 = time.perf_counter()
    logging.info(
        "Time to initialize/write database: "
        + str(round(time1 - time0, 2))
        + " seconds"
    )
    logging.info(
        "Time to perform filtering: " + str(round(time2 - time1, 2)) + " seconds "
    )
    logging.info(cl_opts.parser.epilog)
