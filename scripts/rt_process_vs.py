#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#

import sys
import time
from ringtail import CLOptionParser
from ringtail import RingtailCore
from ringtail import ReceptorManager
import logging
import traceback

if __name__ == "__main__":

    time0 = time.perf_counter()
    level = logging.DEBUG
    logging.basicConfig(
        level=level, stream=sys.stdout, filemode="w", format="%(message)s"
    )
    # parse command line options and filters file (if given)
    try:
        cl_opts = CLOptionParser()
        debug = cl_opts.debug
    except Exception as e:
        tb = traceback.format_exc()
        logging.debug(tb)
        logging.critical("ERROR: " + str(e))
        sys.exit(1)

    # prepare option dictionaries for RingtailCore
    dbman_opts = cl_opts.db_opts
    rman_opts = cl_opts.rman_opts
    filters = cl_opts.filters
    out_opts = cl_opts.out_opts

    # set logging level
    levels = {10: "debug", 20: "info", 30: "warning"}
    no_print = out_opts.pop("no_print")
    if debug:
        level = logging.DEBUG
    elif no_print is False:
        level = logging.INFO
    else:
        level = logging.WARNING
    logging.getLogger().setLevel(level)
    logging.info(f"Logging level set to {levels[level]}")

    # create manager object for virtual screening. Will make database if needed
    try:
        with RingtailCore(
            db_opts=dbman_opts, rman_opts=rman_opts, filters=filters, out_opts=out_opts
        ) as rt_core:

            # Add receptors to database if requested
            # TODO change this
            if cl_opts.save_receptor:
                receptor_list = ReceptorManager.make_receptor_blobs(
                    [rman_opts["receptor_file"]]
                )
                rt_core.add_receptors_to_db(receptor_list)

            time1 = time.perf_counter()

            if cl_opts.rr_mode == "read":
                # perform filtering
                if cl_opts.filter:
                    rt_core.filter()

                # Write log with new data for previous filtering results
                if out_opts["data_from_bookmark"] and not cl_opts.filter:
                    rt_core.get_previous_filter_data()

                # plot if requested
                if out_opts["plot"]:
                    rt_core.plot()

                # write out molecules if requested
                if out_opts["export_sdf_path"] is not None:
                    rt_core.write_molecule_sdfs()

                # write out requested CSVs
                if out_opts["export_bookmark_csv"] is not None:
                    rt_core.export_csv(
                        out_opts["export_bookmark_csv"],
                        out_opts["export_bookmark_csv"] + ".csv",
                        table=True,
                    )
                if out_opts["export_query_csv"] is not None:
                    rt_core.export_csv(out_opts["export_query_csv"], "query.csv")
                if out_opts["export_bookmark_db"] is not None:
                    bookmark_db_name = (
                        dbman_opts["dbFile"].rstrip(".db")
                        + "_"
                        + dbman_opts["results_view_name"]
                        + ".db"
                    )
                    rt_core.export_bookmark_db(bookmark_db_name)

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
