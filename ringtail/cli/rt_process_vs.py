#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""Script that sets up a command line option parser (cloptionparser) and processes all arguments into dictionaries
and options that are then used with the ringtail core api.
This script will allow either a write or a read session at the time.
Available database operations are described in the readme.md document of this codebase.
"""


import sys
import time
from ringtail import CLOptionParser, RingtailCore, setup_logging, get_logger

logger = get_logger(__name__)


def print_db_summary(summary_data: dict, requested_data: dict):
    """
    Formats and prints to stdout a database summary

    Args:
        summary_data (dict): data from database
        requested_data (dict): columns and percentiles requested
    """
    columns = requested_data.get("columns")
    percentiles = requested_data.get("percentiles")
    print("Total Stored Ligands          :", summary_data.pop("num_ligands"))
    print("Total Stored Poses            :", summary_data.pop("num_poses"))
    print(
        "Total Unique Interactions     :",
        summary_data.pop("num_unique_interactions"),
    )
    print(
        "Number Interacting Residues   :",
        summary_data.pop("num_interacting_residues"),
    )

    colon_col = 18
    print("\nEnergy statistics:")
    print("=======================================")
    for col in columns:
        if col == "docking_score":
            min_e = summary_data["min_docking_score"]
            max_e = summary_data["max_docking_score"]
            print(f"Energy (min)      : {min_e:.2f} kcal/mol")
            print(f"Energy (max)      : {max_e:.2f} kcal/mol")
        elif col == "leff":
            min_le = summary_data["min_leff"]
            max_le = summary_data["max_leff"]
            print(f"LE     (min)      : {min_le:.2f} kcal/mol/heavyatom")
            print(f"LE     (max)      : {max_le:.2f} kcal/mol/heavyatom")
        else:
            min_col = summary_data[f"min_{col}"]
            max_col = summary_data[f"max_{col}"]
            print(f"{col} (min) : {min_col}")
            print(f"{col} (max) : {max_col}")
    if percentiles != [] and percentiles is not None:
        print("----------------------------------------")

        for col in columns:
            for p in percentiles:
                if col == "docking_score":
                    p_string = f"Energy (top {p}% )"
                    p_string += " " * (colon_col - len(p_string))
                    print(
                        f"{p_string}: {summary_data[f'{p}%_docking_score']:.2f} kcal/mol"
                    )
                elif col == "leff":
                    p_string = f"LE     (top {p}% )"
                    p_string += " " * (colon_col - len(p_string))
                    print(
                        f"{p_string}: {summary_data[f'{p}%_leff']:.2f} kcal/mol/heavyatom"
                    )
                else:
                    p_string = f"{col} (top {p}%)"
                    p_string += " " * (colon_col - len(p_string))
                    print(f"{p_string}: {summary_data[f'{p}%_{col}']:.2f}")


def main():
    time0 = time.perf_counter()

    setup_logging(level="WARNING")
    try:
        cli = CLOptionParser()
        rtcore = RingtailCore(
            db_file=cli.db_file, storage_type=cli.storage_type, access_mode="cli"
        )
    except Exception as e:
        logger.critical("ERROR: " + str(e))
        return 1

    # create manager object for virtual screening. Will make database if needed
    try:
        if cli.process_mode == "write":
            logger.debug("Starting write process")
            # -#-#- Processes results, will add receptor if "save_receptor" is true
            rtcore.add_results_from_files(**cli.file_sources, **vars(cli.write_options))
        time1 = time.perf_counter()

        # -#-#- Print database summary
        if cli.print_summary:
            summary_data, requested_data = rtcore.db_summary_data()
            print_db_summary(summary_data, requested_data)

        if cli.process_mode == "read":
            bookmark_name = cli.filter_options.bookmark_name
            logger.debug("Starting read process")

            # -#-#- Perform filtering
            if cli.filtering:
                num_passing, bookmark_name = rtcore.filter(
                    **cli.filters, **vars(cli.filter_options)
                )
                print(f"\nNumber of passing ligands: {num_passing}")

            if cli.clustering_only:
                cluster_data = {
                    "ifp": cli.filter_options.interaction_cluster,
                    "mfp": cli.filter_options.mfpt_cluster,
                }
                bookmark_name, _ = rtcore._parse_clustering(
                    cluster_data, cli.filter_options.bookmark_name
                )

            # Write log with new data for previous filtering results
            if (
                cli.output_options.data_from_bookmark and not cli.filtering
            ) or cli.clustering_only:
                rtcore.get_previous_filter_data(
                    bookmark_name,
                    cli.output_options.outfields,
                    cli.output_options.order_results,
                    cli.output_options.output_log,
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
                rtcore.export_table_as_csv(
                    cli.output_options.export_bookmark_csv,
                    cli.output_options.export_bookmark_csv + ".csv",
                )

            # export query as csv
            if cli.output_options.export_query_csv:
                rtcore.export_sql_as_csv(
                    cli.output_options.export_query_csv, "query.csv"
                )

            # export bookmark as database
            if cli.output_options.export_bookmark_db:
                rtcore.export_bookmark_db(bookmark_name)

            # export receptor as .pdbqt
            if cli.output_options.export_receptor_pdbqt:
                rtcore.export_receptor_pdbqt()

            # plot if requested
            if cli.output_options.plot:
                rtcore.plot(bookmark_name)

            # open pymol viewer
            if cli.output_options.pymol:
                rtcore.display_pymol(bookmark_name)

    except Exception as e:
        logger.critical("ERROR: " + str(e))
        return 1

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
    logger.info(cli.parser.epilog)
    return 0


if __name__ == "__main__":
    sys.exit(main())
