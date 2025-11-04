#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail results parsers
#

import os
import gzip
import bz2
import json
import numpy as np
from .exceptions import FileParsingError
from .logutils import LOGGER
from rdkit import Chem
from meeko import PDBQTMolecule, RDKitMolCreate


def make_ringtail_data_dict(
    ligand_row: list, receptor_row: list, data_rows: list, interaction_rows: list
) -> dict:
    """
    Simple wrapper method to unify dict keywords

    Args:
        ligand_row (list): _description_
        receptor_row (list): _description_
        data_rows (list): _description_
        interaction_rows (list): _description_

    Returns:
        _type_: _description_
    """
    return {
        "ligand_row": ligand_row,
        "receptor_row": receptor_row,
        "results_rows": data_rows,
        "interaction_rows": interaction_rows,
    }


def parse_adgpu_results(filename: str, num_poses: int, **kwargs) -> dict:
    """
    Parses an ADGPU dlg (docking log file) and returns dict of lists of rows ready to be inserted in database

    Args:
        filename (str): dlg file name to parse
        num_poses (int): number of poses to store
        interaction_tolerance (float, optional): if wanting to include interactions for poses clustered with passing poses. Defaults to None.

    Returns:
        dict: containing a
            ligand row (list of ligand name, smiles, and binary rdkit mol)
            receptor row (list of receptor data)
            resulst rows* (list of docking results such as energies, coordinates, and other relevant per pose data)
            interaction rows (list of interaction tuples that includes ligname, pose rank, run number to uniqueliy identify
                to which pose an interaction belongs, which makes the interaction list independent of other results
                when it comes to e.g., inserting data in the database)

            * results_row looks like this
                In same order as expected in insert_results:
                receptor, [2]
                pose_rank, [3]
                run_number, [4]
                cluster_rmsd, [5]
                reference_rmsd, [6]
                docking_score, [7]
                leff, [8]
                deltas, [9]
                energies_inter, [10]
                energies_vdw, [11]
                energies_electro, [12]
                energies_flexLig, [13]
                energies_flexLR, [14]
                energies_intra, [15]
                energies_torsional, [16]
                unbound_energy, [17]
                nr_interactions, [18]
                num_hb, [19]
                cluster_size, [20]
                ligand_coordinates, [21]
                flexible_res_coordinates [22]
                ligname [business key]
    """

    def _find_best_cluster_poses(clusters):
        """Takes input ligand dictionary, reads run pose clusters, adds "cluster_best_run"
        entry with the top scoring run for each cluster

        Args:
            ligand_dict (dict): dictionary of ligands

        Returns:
            dict: top poses in cluster for each ligand
        """
        # this can take just clusters and return top poses
        top_poses = []

        for cluster in clusters:
            top_poses.append(clusters[cluster][0])
        return top_poses

    def _parse_docking_file_dlg(fname: str, num_poses: int):
        """Parse an ADGPU DLG file uncompressed or gzipped

        Args:
            fname (str): ligand docking result file name

        Raises:
            ValueError
            FileParsingError

        Returns:
            dict: minimally parsed results dictionaries
        """

        STD_END = "DOCKED: ENDMDL"
        STD_KW = "DOCKED: "

        INPUT_KW = "INPUT LIGAND PDBQT FILE"
        INPUT_END = "FINAL DOCKED STATE"

        # split the first name/extension
        fname_clean = os.path.basename(fname)
        name, ext = os.path.splitext(fname_clean)
        ext = ext[1:].lower()
        if ext == "gz":
            open_fn = gzip.open
            # split the second name/extension
            name, ext = os.path.splitext(name)
            ext = ext[1:].lower()
        else:
            open_fn = open

        # intialize containers for pose data
        interactions = []
        scores = []
        intermolecular_energy = []
        vdw_hb_desolv = []
        electrostatic = []
        flex_ligand = []
        flexLigand_flexReceptor = []
        internal_energy = []
        torsion = []
        unbound_energy = []
        clusters = {}
        cluster_sizes = {}
        cluster_list = []  # list indicating cluster each run belongs to
        pose_interact_count = []
        pose_hb_counts = []
        pose_coordinates = []
        flexible_residues = []
        flexres_startlines = set()
        flexible_res_coords = []

        # Define empty center list for backwards compatibility with DLGs without grid centers
        center = [None, None, None]

        # read poses
        heavy_at_count = 0
        heavy_at_count_complete = False
        file_as_string = ""
        with open_fn(fname, "rb") as fp:
            inside_header = True
            inside_pose = False
            inside_input = False
            inside_res = False
            flexres_atomnames = []
            for line in fp.readlines():
                line = line.decode("utf-8")
                file_as_string += line
                if inside_header:
                    # store ligand file name
                    if line[0:11] == "Ligand file":
                        ligname = (
                            line.split(":", 1)[1].split("/")[-1].split(".")[0].strip()
                        )  # remove path and file extension
                    # store receptor name and grid parameters
                    elif line[:13] == "Receptor name":
                        receptor = line.split()[2]
                    elif line[:21] == "Number of grid points":
                        npts = [
                            pts.rstrip("\n").replace(" ", "")
                            for pts in line.split(":")[1].split(",")
                        ]
                    elif line[:12] == "Grid spacing":
                        spacing = line.split()[2].rstrip(
                            "A"
                        )  # remove A unit from string
                    elif line[:11] == "Grid center":
                        center = [
                            coord.rstrip("\n").replace(" ", "")
                            for coord in line.split(":")[1].split(",")
                        ]
                    # store flexible residue identities and atomtyps
                    if "INPUT-FLEXRES-PDBQT:" in line:
                        if "ATOM" in line or "HETATM" in line:
                            flexres_atomnames[-1].append(line[33:37])
                            if (
                                line[38:41] + ":" + line[42] + line[44:47]
                                in flexible_residues
                            ):
                                continue
                            flexible_residues.append(
                                line[38:41] + ":" + line[42] + line[44:47]
                            )  # RES:<chain><resnum>
                            flexres_startlines.add(line[21:53])  # save startline
                        # add new list for new flexres
                        elif "INPUT-FLEXRES-PDBQT: ROOT" in line:
                            flexres_atomnames.append([])
                    # store number of runs
                    elif "Number of runs:" in line:
                        nruns = int(line.split()[3])
                        cluster_list = list(range(nruns))
                        cluster_rmsds = list(range(nruns))
                        ref_rmsds = list(range(nruns))
                    # store input pdbqt lines
                    elif INPUT_KW in line:
                        inside_input = True
                    elif INPUT_END in line:
                        inside_input = False
                    if inside_input is True:
                        if line.startswith("INPUT-LIGAND-PDBQT"):
                            if (
                                " UNK " in line
                            ):  # replace ligand atoms ATOM flag with HETATM
                                line = line.replace("ATOM", "HETATM")
                    if "FINAL DOCKED STATE" in line:
                        inside_header = False
                    if inside_header:
                        continue

                if "FINAL DOCKED STATE" in line:
                    # first time inside a pose block
                    inside_pose = True
                    interactions.append({})
                    pose_coordinates.append([])
                    flexible_res_coords.append([])
                # store pose anaylsis
                elif line[0:9] == "ANALYSIS:":
                    # storing interactions
                    line = line.split("ANALYSIS:")[1]
                    kw, info = line.split(None, 1)
                    info = info.replace("{", "")
                    info = info.replace("}", "")
                    info = info.replace('"', "")
                    interactions[-1][kw.lower()] = [x.strip() for x in info.split(",")]
                    if "COUNT" in line:
                        interact_count = int(line.split()[1])
                        pose_interact_count.append(str(interact_count))
                        if interact_count == 0:
                            pose_hb_counts.append("0")
                    else:
                        if "TYPE" in line:
                            hb_count = line.count("H")
                            pose_hb_counts.append(hb_count)

                # make new flexible residue list if in the coordinates for a flexible residue
                if line[8:40] in flexres_startlines:
                    flexible_res_coords[-1].append([])
                    inside_res = True
                elif STD_END in line:
                    inside_pose = False
                    inside_res = False
                    heavy_at_count_complete = True
                elif "DOCKED: ROOT" in line:
                    inside_res = False
                if (line[: len(STD_KW)] == STD_KW) and inside_pose:
                    # store the pose raw data
                    line = line.split(STD_KW)[1]
                    # store pose coordinates
                    if "ATOM" in line or "HETATM" in line:
                        if inside_res:
                            flexible_res_coords[-1][-1].append(
                                [line[30:38], line[38:46], line[46:54]]
                            )
                        else:
                            pose_coordinates[-1].append(
                                [line[30:38], line[38:46], line[46:54]]
                            )
                    # store pose data
                    elif "Estimated Free Energy of Binding" in line:
                        try:
                            e = float(line.split()[7])
                        except (
                            ValueError
                        ):  # catch off-by-one error if number is next to =
                            try:
                                e = float(line.split()[6].lstrip("="))
                            except ValueError:
                                raise ValueError(
                                    "ERROR! Cannot parse {0} in {1}".format(line, fname)
                                )
                        finally:
                            if np.isnan(e):
                                raise ValueError(
                                    "Error! File contains NaN value for energy."
                                )
                        scores.append(e)
                    elif "Final Intermolecular Energy" in line:
                        try:
                            e = float(line.split()[6])
                        except (
                            ValueError
                        ):  # catch off-by-one error if number is next to =
                            try:
                                e = float(line.split()[5].lstrip("="))
                            except ValueError:
                                raise ValueError(
                                    "ERROR! Cannot parse {0} in {1}".format(line, fname)
                                )
                        intermolecular_energy.append(e)
                    elif "vdW + Hbond + desolv Energy" in line:
                        try:
                            e = float(line.split()[8])
                        except (
                            ValueError
                        ):  # catch off-by-one error if number is next to =
                            try:
                                e = float(line.split()[7].lstrip("="))
                            except ValueError:
                                raise ValueError(
                                    "ERROR! Cannot parse {0} in {1}".format(line, fname)
                                )
                        vdw_hb_desolv.append(e)
                    elif "Electrostatic Energy" in line:
                        try:
                            e = float(line.split()[4])
                        except (
                            ValueError
                        ):  # catch off-by-one error if number is next to =
                            try:
                                e = float(line.split()[3].lstrip("="))
                            except ValueError:
                                raise ValueError(
                                    "ERROR! Cannot parse {0} in {1}".format(line, fname)
                                )
                        electrostatic.append(e)
                    elif "Moving Ligand-Fixed Receptor" in line:
                        try:
                            e = float(line.split()[5])
                        except (
                            ValueError
                        ):  # catch off-by-one error if number is next to =
                            try:
                                e = float(line.split()[4].lstrip("="))
                            except ValueError:
                                raise ValueError(
                                    "ERROR! Cannot parse {0} in {1}".format(line, fname)
                                )
                        flex_ligand.append(e)
                    elif "Moving Ligand-Moving Receptor" in line:
                        try:
                            e = float(line.split()[5])
                        except (
                            ValueError
                        ):  # catch off-by-one error if number is next to =
                            try:
                                e = float(line.split()[4].lstrip("="))
                            except ValueError:
                                raise ValueError(
                                    "ERROR! Cannot parse {0} in {1}".format(line, fname)
                                )
                        flexLigand_flexReceptor.append(e)
                    elif "Final Total Internal Energy" in line:
                        try:
                            e = float(line.split()[7])
                        except (
                            ValueError
                        ):  # catch off-by-one error if number is next to =
                            try:
                                e = float(line.split()[6].lstrip("="))
                            except ValueError:
                                raise ValueError(
                                    "ERROR! Cannot parse {0} in {1}".format(line, fname)
                                )
                        internal_energy.append(e)
                    elif "Torsional Free Energy" in line:
                        try:
                            e = float(line.split()[6])
                        except (
                            ValueError
                        ):  # catch off-by-one error if number is next to =
                            try:
                                e = float(line.split()[5].lstrip("="))
                            except ValueError:
                                raise ValueError(
                                    "ERROR! Cannot parse {0} in {1}".format(line, fname)
                                )
                        torsion.append(e)
                    elif "Unbound System's Energy" in line:
                        try:
                            e = float(line.split()[6])
                        except (
                            ValueError
                        ):  # catch off-by-one error if number is next to =
                            try:
                                e = float(line.split()[5].lstrip("="))
                            except ValueError:
                                raise ValueError(
                                    "ERROR! Cannot parse {0} in {1}".format(line, fname)
                                )
                        unbound_energy.append(e)

                    # update heavy atom count
                    if heavy_at_count_complete:
                        continue
                    elif line[0:4] == "ATOM" or line[0:6] == "HETATM":
                        # count heavy atoms
                        if not line[-2] == "HD":
                            heavy_at_count += 1

                    continue

                # store poses in each cluster in dictionary as list of ordered runs
                elif "RANKING" in line:
                    cluster_num = line.split()[0]
                    run = line.split()[2]
                    cluster_list[int(run) - 1] = cluster_num
                    if cluster_num in clusters:
                        clusters[cluster_num].append(int(run))
                        cluster_sizes[cluster_num] += 1
                    else:
                        clusters[cluster_num] = [int(run)]
                        cluster_sizes[cluster_num] = 1

                    cluster_rmsds[int(run) - 1] = float(
                        line.split()[4]
                    )  # will be stored in order of runs
                    ref_rmsds[int(run) - 1] = float(line.split()[5])

        # ensure adta is complete
        if (
            len(scores) == 0
            or len(intermolecular_energy) == 0
            or len(vdw_hb_desolv) == 0
            or len(electrostatic) == 0
            or len(internal_energy) == 0
            or len(torsion) == 0
            or len(unbound_energy) == 0
        ):
            raise FileParsingError("Incomplete data in " + fname)

        # sort runs based on docking score
        sorted_idx = np.argsort(scores)

        results = {
            # results columns
            "cluster_rmsds": cluster_rmsds,
            "ref_rmsds": ref_rmsds,
            "scores": scores,
            "leff": [round(x / heavy_at_count, 2) for x in scores],
            "delta": [round(x - scores[0], 2) for x in scores],
            "intermolecular_energy": intermolecular_energy,
            "vdw_hb_desolv": vdw_hb_desolv,
            "electrostatics": electrostatic,
            "flex_ligand": flex_ligand,
            "flexLigand_flexReceptor": flexLigand_flexReceptor,
            "internal_energy": internal_energy,
            "torsional_energy": torsion,
            "unbound_energy": unbound_energy,
            "num_interactions": pose_interact_count,
            "num_hb": pose_hb_counts,
            "cluster_sizes": cluster_sizes,
            "pose_coordinates": pose_coordinates,
            "flexible_res_coordinates": flexible_res_coords,
            # interactions
            "interactions": interactions,
            # for parsing only
            "sorted_runs": [int(x + 1) for x in sorted_idx],
            "clusters": clusters,
            "cluster_list": cluster_list,
        }

        # sort all lists by ranked index
        sorted_results = {
            key: (
                sort_list_by_sorted_idx(value, sorted_idx)
                if isinstance(value, list)
                else value
            )
            for key, value in results.items()
        }
        top_cluster_poses = _find_best_cluster_poses(sorted_results["clusters"])
        sorted_results["poses_to_save"] = pick_best_poses(top_cluster_poses, num_poses)

        receptor_dict = {
            "receptor": receptor,
            "grid_center": center,
            "grid_dim": npts,
            "grid_spacing": spacing,
            "flexible_residues": flexible_residues,
            "flexres_atomnames": flexres_atomnames,
        }

        ligand_dict = {"name": ligname, "file_str": file_as_string}

        return ligand_dict, receptor_dict, sorted_results

    def _find_tolerated_interactions(
        sorted_runs: list, clusters: dict, cluster_rmsds: list, int_tol
    ) -> list[int]:
        """Take ligand dict and finds which poses we should save the
        interactions for as tolerated interactions for the top pose
        of the cluster. These runs are within the
        <self.interaction_tolerance> angstroms RMSD of the top pose
        for a given cluster. All data for the cluster's top pose is saved.

        Args:
            ligand_dict (dict): Dictionary of ligand data from parser

        Returns:
            list: run numbers of tolerated runs
        """
        tolerated_interaction_cluster_runs = {}
        tolerated_runs = []
        for idx, run in enumerate(sorted_runs):
            if float(cluster_rmsds[idx]) <= int_tol:
                tolerated_runs.append(run)
                for key, value in clusters.items():
                    if run in value:
                        tolerated_interaction_cluster_runs.setdefault(key, []).append(
                            run
                        )

        return tolerated_interaction_cluster_runs

    ligand_dict, receptor_dict, data_dict = _parse_docking_file_dlg(
        fname=filename, num_poses=num_poses
    )

    ligand_row = generate_ligand_data_list_from_pdbqt_dlg(
        ligand_dict["name"], ligand_dict["file_str"], is_dlg=True
    )
    ligname = ligand_row[0]

    receptor_row = generate_receptor_row(receptor_dict)
    receptor = receptor_row[0]

    if isinstance(kwargs["interaction_tolerance"], float):
        tolerated_interaction_clusters = _find_tolerated_interactions(
            data_dict["sorted_runs"],
            data_dict["clusters"],
            data_dict["cluster_rmsds"],
            kwargs["interaction_tolerance"],
        )

    passing_interaction_dicts = []
    results_rows = []
    interactions = data_dict.get("interactions", [])
    for index, run_number in enumerate(data_dict["poses_to_save"]):
        # TODO check these indices, I might have mixede something up
        pose_rank = index + 1
        data_idx = run_number - 1

        if interactions[data_idx] not in [[], {}, None]:
            # make shallow copy
            current_interactions: dict = interactions[data_idx].copy()
            # add identifiers
            unique_pose_id = {
                "ligand_name": ligname,
                "run_number": run_number,
                "pose_rank": pose_rank,
            }
            current_interactions.update(unique_pose_id)
            passing_interaction_dicts.append(current_interactions)

        if isinstance(kwargs["interaction_tolerance"], float):
            # check if there are runs tolerated for only interactions in this cluster
            current_cluster_id = data_dict["cluster_list"][data_idx]
            if tol_int_runs := tolerated_interaction_clusters[current_cluster_id] != []:
                # these interactions are technically a different run/pose, but have been clustered with
                # ligname, run number, and pose rank used above (best ranked pose in the cluster)
                for tol_int_run_number in tol_int_runs:
                    tolerated_interactions: dict = interactions[
                        tol_int_run_number - 1
                    ].copy()
                    tolerated_interactions.update(unique_pose_id)
                    passing_interaction_dicts.append(current_interactions)
        cluster_size = data_dict["cluster_sizes"][data_dict["cluster_list"][pose_rank]]
        results_rows.append(
            [
                receptor,
                pose_rank,
                run_number,
                data_dict["cluster_rmsds"][data_idx],
                data_dict["ref_rmsds"][data_idx],
                data_dict["scores"][data_idx],
                data_dict["leff"][data_idx],
                data_dict["delta"][data_idx],
                data_dict["intermolecular_energy"][data_idx],
                data_dict["vdw_hb_desolv"][data_idx],
                data_dict["electrostatics"][data_idx],
                data_dict["flex_ligand"][data_idx],
                data_dict["flexLigand_flexReceptor"][data_idx],
                data_dict["internal_energy"][data_idx],
                data_dict["torsional_energy"][data_idx],
                data_dict["unbound_energy"][data_idx],
                data_dict["num_interactions"][data_idx],
                data_dict["num_hb"][data_idx],
                cluster_size,
                json.dumps(data_dict["pose_coordinates"][data_idx]),
                json.dumps(data_dict["flexible_res_coordinates"][data_idx]),
                ligname,
            ]
        )

    # prepare data in ringtail recognizable format
    return make_ringtail_data_dict(
        ligand_row,
        receptor_row,
        results_rows,
        generate_interaction_tuples(passing_interaction_dicts),
    )


def parse_vina_results(data_pointer, num_poses: int, **kwargs) -> dict:
    """Parser for vina docking results, supporting either pdbqt or gzipped (.gz) files, or with the
    docking results provided as a string.

    Args:
        data_pointer (any): either filename or dictionary of string docking results

    Returns:
        dict: parsed results ready to be inserted in database
    """

    def _read_vina_results_lines(data_object, name) -> dict:
        """Reads vina docking results line by line

        Args:
            data_object (any): filename or dict containing docking results
            name (str): given ligand name/file name

        Raises:
            ValueError: if a line cannot be parsed

        Returns:
            dict: parsed results ready to be inserted in database
        """
        pose_coordinates = []
        scores = []
        sorted_runs = []
        leff = []
        delta = []
        intermolecular_energy = []
        internal_energy = []
        unbound_energy = []
        num_heavy_atoms = 0
        flexible_res_coords = []
        inside_res = False
        flexible_residues = []
        ligand_atomtypes = []
        flexres_atomnames = []
        first_model = True
        cluster = 1  # treat every pose in vina like new cluster
        cluster_list = []
        clusters = {}
        cluster_sizes = {}
        cluster_rmsds = []

        for line in data_object:
            try:
                if line.startswith("MODEL"):
                    cluster_list.append(cluster)
                    clusters[cluster] = [cluster]
                    cluster_sizes[cluster] = 1
                    cluster_rmsds.append(0.0)
                    cluster += 1
                    pose_coordinates.append([])
                    flexible_res_coords.append([])
                    sorted_runs.append(line.split()[1])
                if line.startswith("REMARK VINA RESULT:"):
                    scores.append(float(line.split()[3]))
                if line.startswith("REMARK INTER:"):
                    intermolecular_energy.append(float(line.split()[2]))
                if line.startswith("REMARK INTRA:"):
                    internal_energy.append(float(line.split()[2]))
                if line.startswith("REMARK UNBOUND:"):
                    unbound_energy.append(float(line.split()[2]))
                if line.startswith("HETATM") or line.startswith("ATOM"):
                    if inside_res:
                        flexible_res_coords[-1][-1].append(
                            [line[30:38], line[38:46], line[46:54]]
                        )
                        if first_model:
                            flexres_atomnames[-1].append(line[12:16].strip())
                    else:
                        pose_coordinates[-1].append(
                            [line[30:38], line[38:46], line[46:54]]
                        )
                        if first_model:
                            ligand_atomtypes.append(line[77:].strip())
                            if line[13] != "H":
                                num_heavy_atoms += 1
                if line.startswith("ENDMDL") and first_model:
                    first_model = False
                # make new flexible residue list if in the coordinates for a flexible residue
                if line.startswith("BEGIN_RES"):
                    flexible_res_coords[-1].append([])
                    if first_model:
                        flexres_atomnames.append([])
                    inside_res = True
                if line.startswith("END_RES"):
                    inside_res = False
                # store flexible residue identities
                if line.startswith("BEGIN_RES") and first_model:
                    res = line[10:13].strip()
                    chain = line[14].strip()
                    resnum = line[15:19].strip()
                    res_string = "%s:%s%s" % (res, chain, resnum)
                    flexible_residues.append(res_string)
            except ValueError:
                raise ValueError("ERROR! Cannot parse {0} in {1}".format(line, name))

        # calculate ligand efficiency and deltas from the best pose
        leff = [round(x / num_heavy_atoms, 2) for x in scores]
        delta = [round(x - scores[0], 2) for x in scores]
        ligand_dict = {"name": ligname, "file_str": data_object}

        results = {
            "receptor": None,
            "pose_rank": [1] * len(scores),
            "run_number": sorted_runs,
            "pose_coordinates": pose_coordinates,  # list
            "flexible_res_coordinates": flexible_res_coords,
            "flexible_residues": flexible_residues,
            "flexres_atomnames": flexres_atomnames,
            "clusters": clusters,
            "cluster_rmsds": cluster_rmsds,
            "cluster_sizes": cluster_sizes,
            "cluster_list": cluster_list,
            "scores": scores,  # list
            "leff": leff,  # list
            "delta": delta,  # list
            "intermolecular_energy": intermolecular_energy,  # list
            "internal_energy": internal_energy,  # list
            "unbound_energy": unbound_energy,  # list
            "poses_to_save": pick_best_poses(sorted_runs, num_poses),
            "ligand_atomtypes": ligand_atomtypes,
        }
        return ligand_dict, results

    # check if input is file or data string
    try:
        os.path.splitext(data_pointer)
        from_file = True
    except TypeError:
        from_file = False

    if from_file:
        # read and decode contents
        fname_clean = os.path.basename(data_pointer)
        # split the first name/extension
        name, ext = os.path.splitext(fname_clean)
        ext = ext[1:].lower()
        if ext == "gz":
            open_fn = gzip.open
            # split the second name/extension
            name, ext = os.path.splitext(name)
            ext = ext[1:].lower()
        else:
            open_fn = open
        ligname = data_pointer.split(".pdbqt")[0].split("/")[-1]
        LOGGER.debug("Parsing vina docking file")
        with open_fn(data_pointer, "rb") as fp:
            # this should decode lines into iterable object
            data_object = [line.decode("utf-8") for line in fp]
    else:
        # input provided as string and convert from '\n' separated to to iterable lines
        data_object = list(data_pointer.values())[0].splitlines()
        # get the first (and only, probably not optimal) key which should be the ligand name
        ligname = list(data_pointer.keys())[0]
        LOGGER.debug("Parsing vina docking string")

    ligand_dict, data_dict = _read_vina_results_lines(data_object, ligname)

    ligand_row = generate_ligand_data_list_from_pdbqt_dlg(
        ligand_dict["name"], ligand_dict["file_str"], is_dlg=False
    )
    return make_ringtail_data_dict(ligand_row, [], data_rows, interaction_rows)


def parse_docking_file_sdf(fname: str, num_poses: int) -> dict:
    mols = []
    suppl = Chem.SDMolSupplier(fname)
    for mol in suppl:
        mols.append(mol)
    scores = [mol.GetProp("docking_score") for mol in mols]

    # name them sorted_runs and sort based on docking_score
    sorted_mols = sort_list_by_sorted_idx(mols, np.argsort(scores))

    return {"poses_to_save": pick_best_poses(sorted_mols, num_poses)}


def sort_list_by_sorted_idx(input_list: list, sorted_idx: list[int]) -> list:
    """
    Sort given input list to match order of sorted indices list

    Args:
        input_list (list): list to be sorted
        sorted_idx (list[int]): sorted indices

    Returns:
        list: sorted list
    """
    if input_list == []:
        return input_list
    return [input_list[i] for i in sorted_idx]


def pick_best_poses(pose_based_list: list, num_poses: int = -1) -> list:
    """
    Given max number of poses to save (-1 = all), return given list
    with only the items to be stored

    Args:
        pose_based_list (list): _description_
        num_poses (int, optional): _description_. Defaults to -1.

    Returns:
        list: _description_
    """
    if num_poses == -1:
        return pose_based_list
    else:
        return pose_based_list[:num_poses]


def calculate_interactions():
    # TODO I think it might be better to calculate interactions in storage manager and not in file parser
    # Calculate interactions if requested
    if self.shared.get("no_interactions") != True:
        if self.interaction_finder is None:
            if receptor_string := self.shared.get("receptor_string"):
                self.interaction_finder = InteractionFinder(
                    receptor_string,
                    self.shared.get(
                        "interaction_cutoffs",
                        RingtailDefaults.interaction_cutoffs,
                    ),
                )
            else:
                LOGGER.error(
                    "Cannot calculate interactions, missing receptor representation"
                )

        if parsed_file_dict["interactions"] == []:
            for pose in parsed_file_dict["pose_coordinates"]:
                parsed_file_dict["interactions"].append(
                    self.interaction_finder.find_pose_interactions(
                        parsed_file_dict["ligand_atomtypes"], pose
                    )
                )
                parsed_file_dict["num_interactions"].append(
                    int(parsed_file_dict["interactions"][-1]["count"][0])
                )
                parsed_file_dict["num_hb"].append(
                    len(
                        [
                            1
                            for i in parsed_file_dict["interactions"][-1]["type"]
                            if i == "H"
                        ]
                    )
                )


def generate_ligand_data_list_from_pdbqt_dlg(
    ligname: str, file_str: str, is_dlg=True
) -> list:
    """writes row to be inserted into ligand table

    Args:
        ligand_dict (dict): Dictionary of ligand data from parser

    Returns:
        list: List of data to be written as row in ligand table. Format:
            [ligand_name, ligand_smile, ligand_rdbin]
    """

    # use meeko prepare molecule
    pdbqt_mol = PDBQTMolecule(file_str, name=ligname, is_dlg=is_dlg, skip_typing=True)
    rdkit_mol = RDKitMolCreate.from_pdbqt_mol(pdbqt_mol)[0]
    # remove conformer data
    rdkit_mol.RemoveAllConformers()
    ligand_rdbin = rdkit_mol.ToBinary()

    return [
        ligname,
        Chem.MolToSmiles(rdkit_mol),
        ligand_rdbin,
    ]


def generate_receptor_row(receptor_data: dict) -> list:
    """Writes row to be inserted into receptor table

    Args:
        ligand_dict (dict): Dictionary of ligand data from parser

    Returns:
        list: receptor row columns
    """
    rec_name = receptor_data["receptor"]
    box_dim = json.dumps(receptor_data["grid_dim"])
    box_center = json.dumps(receptor_data["grid_center"])
    grid_spacing = receptor_data["grid_spacing"]
    if grid_spacing != "":
        grid_spacing = float(grid_spacing)
    flexible_residues = json.dumps(receptor_data["flexible_residues"])
    flexres_atomnames = json.dumps(receptor_data["flexres_atomnames"])

    return [
        rec_name,
        box_dim,
        box_center,
        grid_spacing,
        flexible_residues,
        flexres_atomnames,
    ]


def generate_interaction_tuples(interaction_dictionaries: list, unique_id=True) -> list:
    """Takes dictionary of file results, formats as list of tuples for interactions.
    To each interaction description is added business keys/columns that identifies
    which results row/pose each interaction belongs to

    Args:
        interaction_dictionaries (list): List of pose interaction
            dictionaries from parser

    Returns:
        list: of tuples of interaction data
    """
    interaction_keywords = ["type", "chain", "residue", "resid", "recname", "recid"]
    interactions = {
        (
            ((pose["ligand_name"], pose["run_number"], pose["pose_rank"]))
            if unique_id
            else ()
        )
        + tuple(pose[kw][i] for kw in interaction_keywords)
        for pose in interaction_dictionaries
        for i in range(int(pose["count"][0]))
    }

    return list(interactions)


docking_file_parsers: dict[str, callable] = {
    "adgpu": parse_adgpu_results,
    "vina": parse_vina_results,
    "adng": parse_docking_file_sdf,
}
