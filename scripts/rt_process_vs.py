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
        # parse command line options and config file (if given)
        cmdinput = CLOptionParser()
        rtcore: RingtailCore = cmdinput.rtcore
        rtcore._set_general_options(dict=cmdinput.generalopts)
    except Exception as e:
        tb = traceback.format_exc()
        logger.debug(tb)
        logger.critical("ERROR: " + str(e))
        sys.exit(1)

    # create manager object for virtual screening. Will make database if needed
    try:
        rtcore.set_read_options(dict=cmdinput.readopts)
        rtcore.set_storageman_attributes(dict=cmdinput.storageopts)
        readopts = rtcore.readopts 
        if rtcore.process_mode == "write":
            logger.debug("Starting write process")
            #-#-#- Processes results, will add receptor if "save_receptor" is true
            rtcore.add_results_from_files(filesources_dict=cmdinput.file_sources, 
                                          options_dict=cmdinput.writeopts, 
                                          summary=rtcore.print_summary)
        time1 = time.perf_counter()
        if  rtcore.process_mode == "read":
                logger.debug("Starting read process")
                
                #-#-#- Print database summary
                if rtcore.print_summary:
                    rtcore.produce_summary()
                
                #-#-#- Perform filtering
                if readopts.filtering:
                    rtcore.filter(filters_dict=cmdinput.filters,
                                  enumerate_interaction_combs=readopts.enumerate_interaction_combs)

                # Write log with new data for previous filtering results
                if readopts.data_from_bookmark and not readopts.filtering:
                    rtcore.get_previous_filter_data()
                
                # find similar ligands to that specified, if specified (i.e., not None)
                if readopts.find_similar_ligands:
                    rtcore.find_similar_ligands(readopts.find_similar_ligands)

                # plot if requested
                if readopts.plot:
                    rtcore.plot()

                # open pymol viewer
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

                # export query as csv
                if readopts.export_query_csv:
                    rtcore.export_csv(readopts.export_query_csv, "query.csv")

                # export bookmark as database
                if readopts.export_bookmark_db:
                    rtcore.export_bookmark_db(rtcore.storageman.bookmark_name)

                # export receptor as .pdbqt
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
    if logger.level() in ["DEBUG", "INFO"]: print(cmdinput.parser.epilog)
