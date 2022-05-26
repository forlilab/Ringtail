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

if __name__ == '__main__':

    time0 = time.perf_counter()
    # parse command line options and filters file (if given)
    try:
        cl_opts = CLOptionParser()
    except OptionError as e:
        logging.critical(e)
        sys.exit(1)

    # save target name
    if len(cl_opts.rec_files_pool) != 0:
        receptor = cl_opts.rec_files_pool[0].split(".")[0].split("/")[-1]  # remove file extension and path
    else:
        receptor = None

    # prepare option dictionaries for VSManager
    dbman_opts = cl_opts.db_opts
    rman_opts = {'chunk_size': 1,
                 'filelist': cl_opts.lig_files_pool,
                 'mode': cl_opts.mode,
                 'num_clusters': dbman_opts["num_clusters"],
                 'target': receptor}
    filters = cl_opts.filters
    out_opts = cl_opts.out_opts

    # set logging level
    no_print = out_opts["no_print"]
    debug = False
    if debug:
        level = logging.DEBUG
    elif not no_print:
        level = logging.INFO
    else:
        level = logging.WARNING
    logging.basicConfig(level=level)

    # create manager object for virtual screening. Will make database if needed
    try:
        with VSManager(db_opts=dbman_opts,
                       rman_opts=rman_opts,
                       filters=filters,
                       out_opts=out_opts) as vsman:

            # Add receptors to database if requested
            if cl_opts.save_receptor:
                recman = ReceptorManager(cl_opts.rec_files_pool, vsman.dbman)
                recman.add_receptors_to_db()

            time1 = time.perf_counter()

            # perform filtering
            if cl_opts.filter:
                vsman.filter()

            # Write log with new data for previous filtering results
            if out_opts["data_from_subset"] and not cl_opts.filter:
                vsman.get_previous_filter_data()

            # plot if requested
            if out_opts["plot"]:
                vsman.plot()

            # write out molecules if requested
            if out_opts["export_poses_path"] is not None:
                vsman.write_molecule_sdfs()

            # write out requested CSVs
            if out_opts["export_table"] is not None:
                vsman.export_csv(out_opts["export_table"], out_opts["export_table"] + ".csv", table=True)
            if out_opts["export_query"] is not None:
                vsman.export_csv(out_opts["export_query"], "query.csv")

    except VirtualScreeningError as e:
        logging.critical(e)
        sys.exit(1)

    # print performance times
    time2 = time.perf_counter()
    print("Time to initialize/write database: " + str(round(time1 - time0, 2)) + " seconds")
    print("Time to perform filtering: " + str(round(time2 - time1, 2)) + " seconds ")
