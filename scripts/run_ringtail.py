#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#

import time
from ringtail import CLOptionParser
from ringtail import VSManager
from ringtail import ReceptorManager

if __name__ == '__main__':

    time0 = time.perf_counter()
    # parse command line options and filters file (if given)
    cl_opts = CLOptionParser()

    # save target name
    if len(cl_opts.rec_files_pool) != 0:
        receptor = cl_opts.rec_files_pool[0].split(".")[0].split("/")[-1]  # remove file extension and path
    else:
        receptor = None

    # prepare option dictionaries for VSManager
    dbman_opts = cl_opts.db_opts
    rman_opts = {'chunk_size': 1,
                 'filelist': cl_opts.lig_files_pool,
                 'mode': 'dlg',
                 'num_clusters': dbman_opts["num_clusters"],
                 'target': receptor}
    filters = cl_opts.filters
    out_opts = cl_opts.out_opts

    # create manager object for virtual screening. Will make database if needed
    vsman = VSManager(db_opts=dbman_opts,
                      rman_opts=rman_opts,
                      filters=filters,
                      out_opts=out_opts)

    # Add receptors to database if requested
    if cl_opts.save_receptor:
        recman = ReceptorManager(cl_opts.rec_files_pool, vsman.dbman)
        recman.add_receptors_to_db()

    time1 = time.perf_counter()

    # perform filtering
    if not out_opts["no_filter"]:
        vsman.filter()

    # Write log with new data for previous filtering results
    if out_opts["no_filter"] and out_opts["data_from_subset"]:
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

    # close database
    vsman.close_database()
    time2 = time.perf_counter()

    # print performance times
    print("Time to initialize/write database: " + str(round(time1 - time0, 2)) + " seconds")
    print("Time to perform filtering: " + str(round(time2 - time1, 2)) + " seconds ")
