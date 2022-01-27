from cloptionparser import CLOptionParser
from vsmanager import VSManager

if __name__ == '__main__':

    #parse command line options and filters file (if given)
    cl_opts = CLOptionParser()

    #prepare option dictionaries for VSManager
    dbman_opts = cl_opts.output
    dbman_opts["write_db_flag"] = cl_opts.write_db_flag
    rman_opts = {'chunk_size': 1000,
            'filelist': cl_opts.files_pool,
            'mode' : 'dlg',
            'num_clusters':cl_opts.num_clusters
        }

    vsman = VSManager(db_fname = cl_opts.sqlFile, db_opts = dbman_opts, rman_opts = rman_opts)
    pass