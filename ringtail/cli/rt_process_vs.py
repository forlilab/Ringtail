#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#

import sys
import time
from ringtail import CLOptionParser
from ringtail import RingtailCore
from ringtail import logutils
import traceback


def main():
    time0 = time.perf_counter()

    try:
        # set up the logger
        logger = logutils.LOGGER
        # parse command line options and config file (if given)
        cmdinput = CLOptionParser()
        rtcore: RingtailCore = cmdinput.rtcore
    except Exception as e:
        logger.critical("ERROR: " + str(e))
        sys.exit(1)

    # create manager object for virtual screening. Will make database if needed
    try:
        rtcore.set_output_options(dict=cmdinput.outputopts)
        rtcore.set_storageman_attributes(dict=cmdinput.storageopts)
        readopts = cmdinput.readopts
        outopts = rtcore.outputopts
        if cmdinput.process_mode == "write":
            logger.debug("Starting write process")
            # -#-#- Processes results, will add receptor if "save_receptor" is true
            rtcore.add_results_from_files(
                filesources_dict=cmdinput.file_sources, options_dict=cmdinput.writeopts
            )
        time1 = time.perf_counter()

        # -#-#- Print database summary
        if cmdinput.print_summary:
            rtcore.produce_summary()

        if cmdinput.process_mode == "read":
            logger.debug("Starting read process")

            # -#-#- Perform filtering
            if cmdinput.filtering:
                rtcore.filter(
                    filters_dict=cmdinput.filters,
                    enumerate_interaction_combs=outopts.enumerate_interaction_combs,
                )

            # Write log with new data for previous filtering results
            if cmdinput.data_from_bookmark and not cmdinput.filtering:
                rtcore.get_previous_filter_data()

            # find similar ligands to that specified, if specified (i.e., not None)
            if readopts["find_similar_ligands"]:
                rtcore.find_similar_ligands(readopts["find_similar_ligands"])

            # write out molecules if requested
            if outopts.export_sdf_path:
                rtcore.write_molecule_sdfs(sdf_path = outopts.export_sdf_path)

            # write out requested CSVs
            if readopts["export_bookmark_csv"]:
                rtcore.export_csv(
                    readopts["export_bookmark_csv"],
                    readopts["export_bookmark_csv"] + ".csv",
                    table=True,
                )

            # export query as csv
            if readopts["export_query_csv"]:
                rtcore.export_csv(readopts["export_query_csv"], "query.csv")

            # export bookmark as database
            if cmdinput.export_bookmark_db:
                rtcore.export_bookmark_db(rtcore.storageman.bookmark_name)

            # export receptor as .pdbqt
            if cmdinput.export_receptor:
                rtcore.export_receptors()

            # plot if requested
            if cmdinput.plot:
                rtcore.plot()

            # open pymol viewer
            if cmdinput.pymol:
                rtcore.display_pymol()

    except Exception as e:
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
    if logger.level() in ["DEBUG", "INFO"]:
        print(cmdinput.parser.epilog)
    return


if __name__ == "__main__":
    """Script that sets up a command line option parser (cloptionparser) and processes all arguments into dictionaries
    and options that are then used with the ringtail core api.
    This script will allow either a write or a read session at the time.
    Available database operations are described in the readme.md document of this codebase.
    """
    sys.exit(main())
