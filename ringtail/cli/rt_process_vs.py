#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#

import sys
import time
from ringtail import CLOptionParser
from ringtail import RingtailCore
from ringtail import logutils
from ringtail.exceptions import NoInputError
import traceback


def main():
    time0 = time.perf_counter()

    try:
        # set up the logger
        logger = logutils.LOGGER
        # parse command line options and config file (if given)
        cmdinput = CLOptionParser()
        rtcore: RingtailCore = cmdinput.rtcore
    except NoInputError as e:
        sys.tracebacklimit = 0
        logger.critical("ERROR: " + str(e))
        sys.exit(1)
    except BaseException as e:
        traceback.print_exc()
        logger.critical("ERROR: " + str(e))
        sys.exit(1)

    # create manager object for virtual screening. Will make database if needed
    try:
        if cmdinput.process_mode == "write":
            logger.debug("Starting write process")
            # -#-#- Processes results, will add receptor if "save_receptor" is true
            rtcore.add_results_from_files(**cmdinput.file_sources, **cmdinput.writeopts)
        time1 = time.perf_counter()

        # -#-#- Print database summary
        if cmdinput.print_summary:
            rtcore.produce_summary()

        if cmdinput.process_mode == "read":
            logger.debug("Starting read process")
            bookmark_name = cmdinput.filter_options["bookmark_name"]
            # -#-#- Perform filtering
            if cmdinput.filtering:
                _, bookmark_name = rtcore.filter(
                    **cmdinput.filters,
                    **cmdinput.filter_options,
                )

            # Write log with new data for previous filtering results
            if cmdinput.output_options["data_from_bookmark"] and not cmdinput.filtering:
                rtcore.get_previous_filter_data(
                    cmdinput.filter_options["bookmark_name"]
                )

            # find similar ligands to that specified, if specified (i.e., not None)
            if ligname := cmdinput.output_options["find_similar_ligands"]:
                rtcore.find_similar_ligands(
                    ligname, cmdinput.output_options["output_log"]
                )

            # write out molecules if requested
            if cmdinput.output_options["export_sdf_path"]:
                rtcore.write_molecule_sdfs(
                    sdf_path=cmdinput.output_options["export_sdf_path"],
                    all_in_one=not cmdinput.output_options["individual_sdf_files"],
                    bookmark_name=bookmark_name,
                )

            # write out requested CSVs
            if input := cmdinput.output_options["export_bookmark_csv"]:
                if input is True:
                    table_name = bookmark_name
                else:
                    table_name = input
                rtcore.export_csv(
                    table_name,
                    table_name + ".csv",
                    table=True,
                )

            # export query as csv
            if cmdinput.output_options["export_query_csv"]:
                rtcore.export_csv(
                    cmdinput.output_options["export_query_csv"], "query.csv"
                )

            # export bookmark as database
            if cmdinput.output_options["export_bookmark_db"]:
                rtcore.export_bookmark_db(bookmark_name)

            # export receptor as .pdbqt
            if cmdinput.output_options["export_receptor_pdbqt"]:
                rtcore.export_receptors()

            # plot if requested
            if cmdinput.output_options["plot"]:
                rtcore.plot(bookmark_name=bookmark_name)

            # open pymol viewer
            if cmdinput.output_options["pymol"]:
                rtcore.display_pymol(bookmark_name=bookmark_name)

    except BaseException as e:
        traceback.print_exc()
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
