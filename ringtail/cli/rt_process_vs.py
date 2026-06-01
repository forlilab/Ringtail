#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""Script that sets up a command line option parser (cloptionparser) and processes all arguments that are
 used with the ringtail core api.
This script will allow either a write or a read session at the time.
Available database operations are described in ringtail.readthedocs.io.
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

        # -#-#- Print bookmark names
        if cli.print_bookmarks:
            bookmark_names = rtcore.get_bookmark_names()
            print("\nBookmarks in database:")
            for name in bookmark_names:
                print(f"  {name}")

        if cli.process_mode == "read":
            bookmark_name = cli.filter_options.output_bookmark
            logger.debug("Starting read process")

            # -#-#- Perform filtering
            if cli.filtering:
                num_passing, bookmark_name = rtcore.filter(
                    **cli.filters, **vars(cli.filter_options)
                )
                if (
                    cli.filter_options.mfpt_cluster
                    or cli.filter_options.interaction_cluster
                ):
                    print(f"\nNumber of cluster representatives: {num_passing}")
                else:
                    print(f"\nNumber of passing ligands: {num_passing}")

            if cli.clustering_only:
                cluster_data = {}
                if cli.filter_options.mfpt_cluster:
                    cluster_data["mfp"] = cli.filter_options.mfpt_cluster
                if cli.filter_options.interaction_cluster:
                    cluster_data["ifp"] = cli.filter_options.interaction_cluster
                filter_bm = cli.filter_options.input_bookmark
                user_bm = cli.filter_options.output_bookmark
                _, count_reps = rtcore._parse_clustering(
                    cluster_data, user_bm, filter_bm
                )
                bookmark_name = user_bm
                print(f"\nNumber of cluster representative ligands: {count_reps}")

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
                ligname = cli.output_options.find_similar_ligands
                options = rtcore.fetch_cluster_options(ligname)
                print(
                    "Here are the existing clustering groups. Please ensure that your query ligand(s) is a part of the group you select."
                )
                print(
                    "   Choice number   |   Underlying filter bookmark   |   Morgan or interaction fingerprint_cutoff   "
                )
                print(
                    "----------------------------------------------------------------------------------------------------------"
                )
                ids = []
                for cluster_id, filter_window, name in options:
                    ids.append(cluster_id)
                    print(f"{cluster_id}             |    {filter_window}             |    {name}")
                cluster_choice = int(
                    input("Please specify choice number for the cluster you would like to return similar ligands from: ")
                )
                if cluster_choice not in ids:
                    raise ValueError(
                        f"Given cluster number {cluster_choice} does not exist in the database. Please be sure you are specifying an integer in the given cluster options."
                    )
                rtcore.fetch_clustered_similars(ligname, cluster_choice, output_log="cluster_log.txt")

            # write out molecules if requested
            if cli.output_options.export_sdf_path:
                rtcore.write_molecule_sdfs(
                    bookmark_name,
                    sdf_path=cli.output_options.export_sdf_path,
                    all_in_one=not cli.output_options.individual_sdf_files,
                )

            # write out requested CSVs
            if val := cli.output_options.export_bookmark_csv:
                # export whole table/all results it not specifying outfields
                filename = (
                    (val if val.endswith(".csv") else val + ".csv")
                    if isinstance(val, str)
                    else bookmark_name + ".csv"
                )
                if not cli.output_options.outfields:
                    rtcore.export_table_as_csv(
                        bookmark_name,
                        filename,
                    )
                # If export csv and outfields are both used
                if cli.output_options.outfields:
                    outs = cli.output_options.outfields
                    # ensure outfields is list of strings:
                    if isinstance(outs, str):
                        outs = [item.strip() for item in outs.split(",")]
                    rtcore.export_columns_as_csv(
                        outs,
                        bookmark_name,
                        filename,
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
