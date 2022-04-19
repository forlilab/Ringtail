import time
from cloptionparser import CLOptionParser
from vsmanager import VSManager

if __name__ == '__main__':

    time0 = time.perf_counter()
    #parse command line options and filters file (if given)
    cl_opts = CLOptionParser()

    #prepare option dictionaries for VSManager
    dbman_opts = cl_opts.db_opts
    rman_opts = {'chunk_size': 10000,
            'filelist': cl_opts.files_pool,
            'mode' : 'dlg',
            'num_clusters':dbman_opts["num_clusters"]
        }
    filters = cl_opts.filters
    out_opts = cl_opts.out_opts

    #create manager object for virtual screening. Will write database if needed
    vsman = VSManager(db_opts = dbman_opts, rman_opts = rman_opts, filters=filters, out_opts = out_opts)
    time1 = time.perf_counter()

    #perform filtering
    if not out_opts["no_filter"]:
        vsman.filter()

    # Write log with new data for previous filtering results if we did not filter
    if out_opts["no_filter"] and out_opts["data_from_subset"]:
        vsman.get_previous_filter_data()

    #plot if requested
    if out_opts["plot"]:
        vsman.plot()

    #write out molecules if requested
    if out_opts["export_poses_path"] is not None:
        vsman.write_molecule_sdfs()

    #close database
    vsman.close_database()
    time2 = time.perf_counter()

    #print performance times
    print("Time to initialize/write database: " + str(round(time1-time0,2)) + " seconds")
    print("Time to perform filtering: " + str(round(time2-time1,2)) + " seconds ")