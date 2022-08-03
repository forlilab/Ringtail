#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#

import sys
import time
from ringtail import CLOptionParser
from ringtail import VSManager
from ringtail import ReceptorManager
from ringtail import VirtualScreeningError, OptionError
import logging
import traceback

if __name__ == "__main__":

    time0 = time.perf_counter()
    debug = True
    if debug:
        level = logging.DEBUG
    else:
        level = logging.INFO
    logging.basicConfig(
        level=level, stream=sys.stdout, filemode="w", format="%(message)s"
    )
    # parse command line options and filters file (if given)
    try:
        cl_opts = CLOptionParser()
    except OptionError as e:
        tb = traceback.format_exc()
        logging.debug(tb)
        logging.critical(e)
        sys.exit(1)

    # prepare option dictionaries for VSManager
    dbman_opts = cl_opts.db_opts
    rman_opts = cl_opts.rman_opts
    # save target name
    if rman_opts["receptor_file"] is not None:
        receptor = (
            rman_opts["receptor_file"].split(".")[0].split("/")[-1]
        )  # remove file extension and path
    else:
        receptor = None
    rman_opts = {
        "chunk_size": 1,
        "mode": rman_opts["mode"],
        "max_poses": rman_opts["max_poses"],
        "interaction_tolerance": rman_opts["interaction_tolerance"],
        "store_all_poses": rman_opts["store_all_poses"],
        "target": receptor,
        "receptor_file": rman_opts["receptor_file"],
        "add_interactions": rman_opts["add_interactions"],
        "interaction_cutoffs": rman_opts["interaction_cutoffs"],
        "file_sources": rman_opts["file_sources"],
        "file_pattern": rman_opts["file_pattern"],
    }
    filters = cl_opts.filters
    out_opts = cl_opts.out_opts

    # set logging level
    levels = {10: "debug", 20: "info", 30: "error"}
    no_print = out_opts["no_print"]
    if no_print is False and debug is False:
        level = logging.INFO
    elif debug is False:
        level = logging.WARNING
    logging.getLogger().setLevel(level)
    logging.info(f"Logging level set to {levels[level]}")

    # create manager object for virtual screening. Will make database if needed
    try:
        with VSManager(
            db_opts=dbman_opts, rman_opts=rman_opts, filters=filters, out_opts=out_opts
        ) as vsman:

            # Add receptors to database if requested
            if cl_opts.save_receptor:
                receptor_list = ReceptorManager.make_receptor_blobs(
                    [rman_opts["receptor_file"]]
                )
                vsman.add_receptors_to_db(receptor_list)

            time1 = time.perf_counter()

            if cl_opts.rr_mode == "read":
                # perform filtering
                if cl_opts.filter:
                    vsman.filter()

                # Write log with new data for previous filtering results
                if out_opts["data_from_bookmark"] and not cl_opts.filter:
                    vsman.get_previous_filter_data()

                # plot if requested
                if out_opts["plot"]:
                    vsman.plot()

                # write out molecules if requested
                if out_opts["export_poses_path"] is not None:
                    vsman.write_molecule_sdfs()

                # write out requested CSVs
                if out_opts["export_table"] is not None:
                    vsman.export_csv(
                        out_opts["export_table"],
                        out_opts["export_table"] + ".csv",
                        table=True,
                    )
                if out_opts["export_query"] is not None:
                    vsman.export_csv(out_opts["export_query"], "query.csv")
                if out_opts["export_bookmark_db"] is not None:
                    bookmark_db_name = (
                        dbman_opts["dbFile"].rstrip(".db")
                        + "_"
                        + dbman_opts["results_view_name"]
                        + ".db"
                    )
                    vsman.export_bookmark_db(bookmark_db_name)

    except VirtualScreeningError as e:
        tb = traceback.format_exc()
        logging.debug(tb)
        logging.critical(e)
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
