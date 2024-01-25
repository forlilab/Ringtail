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
    ''' The altered rt_process_vs now makes use of the new ringtail api. The CLOptionParser just makes the option objects, 
    and this script assigns the objects to the relevant managers then runs the ringtail methods with the new API.
    I am not using ringtail core as a context anymore'''
    time0 = time.perf_counter()
    level = logging.INFO
    logging.basicConfig(
        level=level, stream=sys.stdout, filemode="w", format="%(message)s"
    )


    try:
        # parse command line options and filters file (if given)
        cmdinput = CLOptionParser()
        rtopts = cmdinput.rtopts
    except Exception as e:
        tb = traceback.format_exc()
        logging.debug(tb)
        logging.critical("ERROR: " + str(e))
        sys.exit(1)

    # set logging level
    levels = {10: "debug", 20: "info", 30: "warning"}
    debug = rtopts.debug
    verbose = rtopts.verbose
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
        rtcore = RingtailCore(rtopts.db_file)
        read_opts = cmdinput.read_opts

        #-#-#- Universal/shared options
        rtcore.set_storage_options(storageopts = cmdinput.storageman_opts)
        rtcore.set_general_options(process_mode=rtopts.process_mode, rtopts = rtopts)
        rtcore.open()
        if rtopts.process_mode == "write":
            #-#-#- Set write options to the write managers, and processes results
            #-#-#- Will add receptor if "save_receptor" is true
            rtcore.add_results_from_files(file_source_object=cmdinput.file_sources, 
                                          writeopts=cmdinput.write_opts)

        time1 = time.perf_counter()
        if rtopts.process_mode == "read":
                #-#-#- Set read options to the read managers, inc result and output man right? 
                # perform filtering
                if rtcore.summary:
                    rtcore._produce_summary()
                    
                rtcore.set_read_options(readopts=read_opts)
                rtcore.filterobj = cmdinput.filters
                if read_opts.filtering:
                    rtcore.filter(read_opts.enumerate_interaction_combs)
                # Write log with new data for previous filtering results
                if read_opts.data_from_bookmark and not read_opts.filtering:
                    rtcore.get_previous_filter_data()
                
                if read_opts.find_similar_ligands:
                    rtcore.find_similar_ligands(read_opts.find_similar_ligands)

                # plot if requested
                if read_opts.plot:
                    rtcore.plot()

                # def open pymol viewer
                if read_opts.pymol:
                    rtcore.display_pymol()

                # write out molecules if requested
                if read_opts.export_sdf_path:
                    rtcore.write_molecule_sdfs()

                # write out requested CSVs
                if read_opts.export_bookmark_csv:
                    rtcore.export_csv(
                        read_opts.export_bookmark_csv,
                        read_opts.export_bookmark_csv + ".csv",
                        table=True,)

                if read_opts.export_query_csv:
                    rtcore.export_csv(read_opts.export_query_csv, "query.csv")

                if read_opts.export_bookmark_db:
                    bookmark_name = (
                        rtcore.storageman.db_file.rstrip(".db")
                        + "_"
                        + rtcore.storageman.results_view_name
                        + ".db"
                    )
                    rtcore.export_bookmark_db(bookmark_name)

                if read_opts.export_receptor:
                    rtcore.export_receptors()
        rtcore.close_storage()

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
    logging.info(cmdinput.parser.epilog)
