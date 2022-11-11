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
    storage_opts = cl_opts.storage_opts
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
        defaults = RingtailCore.get_defaults()
        with RingtailCore(
            storage_opts={
                "values": storage_opts,
                "types": defaults["storage_opts"]["types"],
            },
            rman_opts={"values": rman_opts, "types": defaults["rman_opts"]["types"]},
            filters={"values": filters, "types": defaults["filters"]["types"]},
            out_opts={"values": out_opts, "types": defaults["out_opts"]["types"]},
        ) as rt_core:

            file_sources_empty = rt_core.get_defaults()["rman_opts"]["values"][
                "file_sources"
            ]

            if cl_opts.rr_mode == "write" and (
                rman_opts["file_sources"] != file_sources_empty
            ):
                rt_core.add_results()

            # Add receptors to database if requested
            if cl_opts.save_receptor:
                rt_core.save_receptors(rman_opts["receptor_file"])

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
                    bookmark_name = (
                        storage_opts["dbFile"].rstrip(".db")
                        + "_"
                        + storage_opts["results_view_name"]
                        + ".db"
                    )
                    rt_core.export_bookmark_db(bookmark_name)

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
