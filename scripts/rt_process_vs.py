#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#

import sys
import time
from ringtail import CLOptionParser
from ringtail import RingtailCore, logger
import traceback

if __name__ == "__main__":
    ''' The altered rt_process_vs now makes use of the new ringtail api. The CLOptionParser just makes the option objects, 
    and this script assigns the objects to the relevant managers then runs the ringtail methods with the new API.
    I am not using ringtail core as a context anymore'''
    time0 = time.perf_counter()

    try:
        # parse command line options and filters file (if given)
        cmdinput = CLOptionParser()
        rtcore: RingtailCore = cmdinput.rtcore
        rtopts = rtcore.generalopts
        # Set logging level
        if rtopts.debug:
            logger.setLevel("DEBUG")
        elif rtopts.verbose:
            logger.setLevel("INFO")
        else:
            logger.setLevel("WARNING")
        logger.warning(f"Logger set to level {logger.level()}")
    except Exception as e:
        tb = traceback.format_exc()
        logger.debug(tb)
        logger.critical("ERROR: " + str(e))
        sys.exit(1)

    # create manager object for virtual screening. Will make database if needed
    try:
        readopts = rtcore.readopts
        if rtcore.process_mode == "write":
            logger.debug("Starting write process")
            #-#-#- Processes results, will add receptor if "save_receptor" is true
            rtcore.add_results_from_files(file_source_object=rtcore.files)

        time1 = time.perf_counter()
        if  rtcore.process_mode == "read":
                logger.debug("Starting read process")
                #-#-#- Perform filtering
                if rtcore.generalopts.summary:
                    rtcore._produce_summary()
                    
                if readopts.filtering:
                    rtcore.filter(readopts.enumerate_interaction_combs)
                # Write log with new data for previous filtering results
                if readopts.data_from_bookmark and not readopts.filtering:
                    rtcore.get_previous_filter_data()
                
                if readopts.find_similar_ligands:
                    rtcore.find_similar_ligands(readopts.find_similar_ligands)

                # plot if requested
                if readopts.plot:
                    rtcore.plot()

                # def open pymol viewer
                if readopts.pymol:
                    rtcore.display_pymol()

                # write out molecules if requested
                if readopts.export_sdf_path:
                    rtcore.write_molecule_sdfs()

                # write out requested CSVs
                if readopts.export_bookmark_csv:
                    rtcore.export_csv(
                        readopts.export_bookmark_csv,
                        readopts.export_bookmark_csv + ".csv",
                        table=True,)

                if readopts.export_query_csv:
                    rtcore.export_csv(readopts.export_query_csv, "query.csv")

                if readopts.export_bookmark_db:
                    bookmark_name = (
                        rtcore.storageman.db_file.rstrip(".db")
                        + "_"
                        + rtcore.storageman.results_view_name
                        + ".db"
                    )
                    rtcore.export_bookmark_db(bookmark_name)

                if readopts.export_receptor:
                    rtcore.export_receptors()

    #TODO can depreciate use of traceback in this file
    except Exception as e:
        tb = traceback.format_exc()
        logger.debug(tb)
        logger.critical("ERROR: " + str(e))
        sys.exit(1)

    # print performance times
    time2 = time.perf_counter()
    logger.info(
        "Time to initialize/write database: "
        + str(round(time1 - time0, 2))
        + " seconds"
    )
    logger.info(
        "Time to perform filtering: " + str(round(time2 - time1, 2)) + " seconds "
    )
    if rtopts.debug or rtopts.verbose: print(cmdinput.parser.epilog)
