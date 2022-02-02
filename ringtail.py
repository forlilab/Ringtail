

from cloptionparser import CLOptionParser
from vsmanager import VSManager, Outputter

if __name__ == '__main__':

    #parse command line options and filters file (if given)
    cl_opts = CLOptionParser()

    #prepare option dictionaries for VSManager
    dbman_opts = cl_opts.db_opts
    rman_opts = {'chunk_size': 1000,
            'filelist': cl_opts.files_pool,
            'mode' : 'dlg',
            'num_clusters':cl_opts.num_clusters
        }
    filters = cl_opts.filters
    out_opts = cl_opts.out_opts

    #create manager object for virtual screening. Will write database if needed
    vsman = VSManager(db_opts = dbman_opts, rman_opts = rman_opts, filters=filters, out_opts = out_opts)

    #perform filtering
    vsman.filter()
    
    #close database
    vsman.close_database()