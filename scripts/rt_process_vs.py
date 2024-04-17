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
    except Exception as e:
        tb = traceback.format_exc()
        logger.debug(tb)
        logger.critical("ERROR: " + str(e))
        sys.exit(1)

    # create manager object for virtual screening. Will make database if needed
    try:
        rtcore.set_output_options(dict=cmdinput.outputopts)
        rtcore.set_storageman_attributes(dict=cmdinput.storageopts)
        outopts = rtcore.outputopts 
        if cmdinput.process_mode == "write":
            logger.debug("Starting write process")
            #-#-#- Processes results, will add receptor if "save_receptor" is true
            rtcore.add_results_from_files(filesources_dict=cmdinput.file_sources, 
                                          options_dict=cmdinput.writeopts, 
                                          summary=cmdinput.print_summary)
        time1 = time.perf_counter()
        if  cmdinput.process_mode == "read":
                logger.debug("Starting read process")
                
                #-#-#- Print database summary
                if cmdinput.print_summary:
                    rtcore.produce_summary()
                
                #-#-#- Perform filtering
                if cmdinput.filtering:
                    rtcore.filter(filters_dict=cmdinput.filters,
                                  enumerate_interaction_combs=outopts.enumerate_interaction_combs)

                # Write log with new data for previous filtering results
                if cmdinput.data_from_bookmark and not cmdinput.filtering:
                    rtcore.get_previous_filter_data()
                
                # find similar ligands to that specified, if specified (i.e., not None)
                if outopts.find_similar_ligands:
                    rtcore.find_similar_ligands(outopts.find_similar_ligands)

                # plot if requested
                if cmdinput.plot:
                    rtcore.plot()

                # open pymol viewer
                if cmdinput.pymol:
                    rtcore.display_pymol()

                # write out molecules if requested
                if outopts.export_sdf_path:
                    rtcore.write_molecule_sdfs()

                # write out requested CSVs
                if outopts.export_bookmark_csv:
                    rtcore.export_csv(
                        outopts.export_bookmark_csv,
                        outopts.export_bookmark_csv + ".csv",
                        table=True,)

                # export query as csv
                if outopts.export_query_csv:
                    rtcore.export_csv(outopts.export_query_csv, "query.csv")

                # export bookmark as database
                if cmdinput.export_bookmark_db:
                    rtcore.export_bookmark_db(rtcore.storageman.bookmark_name)

                # export receptor as .pdbqt
                if cmdinput.export_receptor:
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
