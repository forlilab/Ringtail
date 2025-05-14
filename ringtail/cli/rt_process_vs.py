#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#

import sys
import time
from ringtail import CLOptionParser
from ringtail import RingtailCore
from ringtail import logutils


def main():
    time0 = time.perf_counter()

    try:
        # set up the logger
        logger = logutils.LOGGER
        # parse command line options and config file (if given)
        cli = CLOptionParser()
        rtcore: RingtailCore = cli.rtcore
    except Exception as e:
        logger.critical("ERROR: " + str(e))
        sys.exit(1)

    # create manager object for virtual screening. Will make database if needed
    try:
        if cli.process_mode == "write":
            logger.debug("Starting write process")
            # -#-#- Processes results, will add receptor if "save_receptor" is true
            rtcore.add_results_from_files(**cli.file_sources, **vars(cli.write_options))
        time1 = time.perf_counter()

        # -#-#- Print database summary
        if cli.print_summary:
            rtcore.produce_summary()

        if cli.process_mode == "read":
            bookmark_name = cli.filter_options.bookmark_name
            logger.debug("Starting read process")

            # -#-#- Perform filtering
            if cli.filtering:
                filtered_data = rtcore.filter(**cli.filters, **vars(cli.filter_options))
                if cli.output_options.log_file:
                    rtcore.write_filter_log(
                        filtered_data,
                        cli.output_options.outfields,
                        cli.output_options.order_results,
                        not cli.output_options.output_all_poses,
                        cli.output_options.log_file,
                    )

            # Write log with new data for previous filtering results
            if cli.output_options.data_from_bookmark and not cli.filtering:
                rtcore.get_previous_filter_data(
                    bookmark_name,
                    cli.output_options.outfields,
                    cli.output_options.order_results,
                    cli.output_options.log_file,
                )

            # find similar ligands to that specified, if specified (i.e., not None)
            if cli.output_options.find_similar_ligands:
                rtcore.find_similar_ligands(cli.output_options.find_similar_ligands)

            # write out molecules if requested
            if cli.output_options.export_sdf_path:
                rtcore.write_molecule_sdfs(
                    bookmark_name,
                    sdf_path=cli.output_options.export_sdf_path,
                    all_in_one=not cli.output_options.individual_sdf_files,
                )

            # write out requested CSVs
            if cli.output_options.export_bookmark_csv:
                rtcore.export_csv(
                    cli.output_options.export_bookmark_csv,
                    cli.output_options.export_bookmark_csv + ".csv",
                    table=True,
                )

            # export query as csv
            if cli.output_options.export_query_csv:
                rtcore.export_csv(cli.output_options.export_query_csv, "query.csv")

            # export bookmark as database
            if cli.output_options.export_bookmark_db:
                rtcore.export_bookmark_db(bookmark_name)

            # export receptor as .pdbqt
            if cli.output_options.export_receptor:
                rtcore.export_receptors()

            # plot if requested
            if cli.output_options.plot:
                rtcore.plot(bookmark_name)

            # open pymol viewer
            if cli.output_options.pymol:
                rtcore.display_pymol(bookmark_name)

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
        print(cli.parser.epilog)
    return


if __name__ == "__main__":
    """Script that sets up a command line option parser (cloptionparser) and processes all arguments into dictionaries
    and options that are then used with the ringtail core api.
    This script will allow either a write or a read session at the time.
    Available database operations are described in the readme.md document of this codebase.
    """
    sys.exit(main())
