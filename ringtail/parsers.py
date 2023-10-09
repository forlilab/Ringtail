#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail file parsers
#

import os
import gzip
import numpy as np
from .exceptions import FileParsingError


def parse_single_dlg(fname):
    """parse an ADGPU DLG file uncompressed or gzipped"""
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
    pose_about = []
    pose_trans = []
    pose_quarternions = []
    pose_dihedrals = []
    clusters = {}
    cluster_sizes = {}
    cluster_list = []  # list indicating cluster each run belongs to
    pose_interact_count = []
    pose_hb_counts = []
    pose_coordinates = []
    flexible_residues = []
    flexres_startlines = set()
    flexible_res_coords = []
    smile_string = ""
    index_map = ""
    h_parents = ""

    # Define empty center list for backwards compatibility with DLGs without grid centers
    center = [None, None, None]

    # read poses
    heavy_at_count = 0
    heavy_at_count_complete = False
    with open_fn(fname, "rb") as fp:
        inside_header = True
        inside_pose = False
        inside_input = False
        inside_res = False
        smile_string = ""
        input_model = []
        index_map = []
        h_parents = []
        ligand_atomtypes = []
        flexres_atomnames = []
        for line in fp.readlines():
            line = line.decode("utf-8")
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
                    spacing = line.split()[2].rstrip("A")  # remove A unit from string
                elif line[:11] == "Grid center":
                    center = [
                        coord.rstrip("\n").replace(" ", "")
                        for coord in line.split(":")[1].split(",")
                    ]
                # store smile string
                elif (
                    "REMARK SMILES" in line and "IDX" not in line and smile_string == ""
                ):
                    smile_string = line.split("REMARK SMILES")[-1]
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
                        input_model.append(line[20:])
                        # save ligand atomtypes
                        if line.startswith("INPUT-LIGAND-PDBQT") and (
                            "ATOM" in line or "HETATM" in line
                        ):
                            ligand_atomtypes.append(line.strip()[97:])
                    if line.startswith("INPUT-LIGAND-PDBQT: REMARK SMILES IDX"):
                        index_map += (
                            line.lstrip("INPUT-LIGAND-PDBQT: REMARK SMILES IDX")
                            .rstrip("\n")
                            .split()
                        )
                    if line.startswith("INPUT-LIGAND-PDBQT: REMARK H PARENT"):
                        h_parents += (
                            line.lstrip("INPUT-LIGAND-PDBQT: REMARK H PARENT")
                            .rstrip("\n")
                            .split()
                        )
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
                    except ValueError:  # catch off-by-one error if number is next to =
                        try:
                            e = float(line.split()[6].lstrip("="))
                        except ValueError:
                            raise ValueError(
                                "ERROR! Cannot parse {0} in {1}".format(line, fname)
                            )
                    finally:
                        if np.isnan(e):
                            raise ValueError("Error! File contains NaN value for energy.")
                    scores.append(e)
                elif "Final Intermolecular Energy" in line:
                    try:
                        e = float(line.split()[6])
                    except ValueError:  # catch off-by-one error if number is next to =
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
                    except ValueError:  # catch off-by-one error if number is next to =
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
                    except ValueError:  # catch off-by-one error if number is next to =
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
                    except ValueError:  # catch off-by-one error if number is next to =
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
                    except ValueError:  # catch off-by-one error if number is next to =
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
                    except ValueError:  # catch off-by-one error if number is next to =
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
                    except ValueError:  # catch off-by-one error if number is next to =
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
                    except ValueError:  # catch off-by-one error if number is next to =
                        try:
                            e = float(line.split()[5].lstrip("="))
                        except ValueError:
                            raise ValueError(
                                "ERROR! Cannot parse {0} in {1}".format(line, fname)
                            )
                    unbound_energy.append(e)
                # store state variables
                elif "NEWDPF about" in line:
                    ind_pose_abt = [float(i) for i in line.split()[3:]]
                    pose_about.append(ind_pose_abt)
                elif "NEWDPF tran0" in line:
                    translations = [float(i) for i in line.split()[3:]]
                    pose_trans.append(translations)
                elif "NEWDPF axisangle0" in line:
                    quaternion = [float(i) for i in line.split()[3:]]
                    pose_quarternions.append(quaternion)
                elif "NEWDPF dihe0" in line:
                    dihedrals = [float(i) for i in line.split()[3:]]
                    pose_dihedrals.append(dihedrals)

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

    sorted_idx = np.argsort(scores)
    # sort poses, scores, and interactions

    def sort_list_by_sorted_idx(input_list):
        """Sort given input list to match order of sorted indices list"""
        if input_list == []:
            return input_list
        return [input_list[i] for i in sorted_idx]

    pose_coordinates = sort_list_by_sorted_idx(pose_coordinates)
    flexible_res_coords = sort_list_by_sorted_idx(flexible_res_coords)
    scores = sort_list_by_sorted_idx(scores)
    interactions = sort_list_by_sorted_idx(interactions)
    intermolecular_energy = sort_list_by_sorted_idx(intermolecular_energy)
    vdw_hb_desolv = sort_list_by_sorted_idx(vdw_hb_desolv)
    electrostatic = sort_list_by_sorted_idx(electrostatic)
    flex_ligand = sort_list_by_sorted_idx(flex_ligand)
    flexLigand_flexReceptor = sort_list_by_sorted_idx(flexLigand_flexReceptor)
    internal_energy = sort_list_by_sorted_idx(internal_energy)
    torsion = sort_list_by_sorted_idx(torsion)
    unbound_energy = sort_list_by_sorted_idx(unbound_energy)
    pose_about = sort_list_by_sorted_idx(pose_about)
    pose_trans = sort_list_by_sorted_idx(pose_trans)
    pose_quarternions = sort_list_by_sorted_idx(pose_quarternions)
    pose_dihedrals = sort_list_by_sorted_idx(pose_dihedrals)
    cluster_rmsds = sort_list_by_sorted_idx(cluster_rmsds)
    ref_rmsds = sort_list_by_sorted_idx(ref_rmsds)
    pose_interact_count = sort_list_by_sorted_idx(pose_interact_count)
    pose_hb_counts = sort_list_by_sorted_idx(pose_hb_counts)
    cluster_list = sort_list_by_sorted_idx(cluster_list)

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
    # calculate ligand efficiency and deltas from the best pose
    leff = [x / heavy_at_count for x in scores]
    delta = [x - scores[0] for x in scores]

    return {
        "ligname": ligname,
        "receptor": receptor,
        "grid_center": center,
        "grid_dim": npts,
        "grid_spacing": spacing,
        "ligand_input_model": input_model,
        "ligand_index_map": index_map,
        "ligand_h_parents": h_parents,
        "pose_coordinates": pose_coordinates,
        "flexible_res_coordinates": flexible_res_coords,
        "flexible_residues": flexible_residues,
        "ligand_smile_string": smile_string,
        "clusters": clusters,
        "cluster_rmsds": cluster_rmsds,
        "cluster_sizes": cluster_sizes,
        "cluster_list": cluster_list,
        "ref_rmsds": ref_rmsds,
        "scores": scores,
        "leff": leff,
        "delta": delta,
        "intermolecular_energy": intermolecular_energy,
        "vdw_hb_desolv": vdw_hb_desolv,
        "electrostatics": electrostatic,
        "flex_ligand": flex_ligand,
        "flexLigand_flexReceptor": flexLigand_flexReceptor,
        "internal_energy": internal_energy,
        "torsional_energy": torsion,
        "unbound_energy": unbound_energy,
        "interactions": interactions,
        "num_interactions": pose_interact_count,
        "num_hb": pose_hb_counts,
        "sorted_runs": [int(x + 1) for x in sorted_idx],
        "pose_about": pose_about,
        "pose_translations": pose_trans,
        "pose_quarternions": pose_quarternions,
        "pose_dihedrals": pose_dihedrals,
        "fname": fname,
        "ligand_atomtypes": ligand_atomtypes,
        "flexres_atomnames": flexres_atomnames,
    }


def parse_vina_pdbqt(fname):

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

    ligname = fname.split(".pdbqt")[0].split("/")[-1]
    pose_coordinates = []
    scores = []
    sorted_runs = []
    leff = []
    delta = []
    intermolecular_energy = []
    internal_energy = []
    unbound_energy = []
    num_heavy_atoms = 0
    smile_string = ""
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
    smile_idx_map = []
    ligand_h_parents = []

    with open_fn(fname, "rb") as fp:
        for line in fp.readlines():
            line = line.decode("utf-8")
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
                        flexible_res_coords[-1][-1].append([line[30:38], line[38:46], line[46:54]])
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
                if line.startswith("REMARK SMILES IDX") and first_model:
                    smile_idx_map += (
                        line.lstrip("REMARK SMILES IDX").rstrip("\n").split()
                    )
                elif line.startswith("REMARK SMILES") and first_model:
                    smile_string = line.split()[2]
                if line.startswith("REMARK H PARENT") and first_model:
                    ligand_h_parents += (
                        line.lstrip("REMARK H PARENT").rstrip("\n").split()
                    )
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
                raise ValueError("ERROR! Cannot parse {0} in {1}".format(line, fname))

    # calculate ligand efficiency and deltas from the best pose
    leff = [x / num_heavy_atoms for x in scores]
    delta = [x - scores[0] for x in scores]

    return {
        "ligname": ligname,  # string
        "receptor": "",
        "grid_center": "",
        "grid_dim": "",
        "grid_spacing": "",
        "ligand_input_model": "",
        "ligand_index_map": smile_idx_map,
        "ligand_h_parents": ligand_h_parents,
        "pose_coordinates": pose_coordinates,  # list
        "flexible_res_coordinates": flexible_res_coords,
        "flexible_residues": flexible_residues,
        "flexres_atomnames": flexres_atomnames,
        "ligand_smile_string": smile_string,
        "clusters": clusters,
        "cluster_rmsds": cluster_rmsds,
        "cluster_sizes": cluster_sizes,
        "cluster_list": cluster_list,
        "ref_rmsds": [],
        "scores": scores,  # list
        "leff": leff,  # list
        "delta": delta,  # list
        "intermolecular_energy": intermolecular_energy,  # list
        "vdw_hb_desolv": [],
        "electrostatics": [],
        "flex_ligand": [],
        "flexLigand_flexReceptor": [],
        "internal_energy": internal_energy,  # list
        "torsional_energy": [],
        "unbound_energy": unbound_energy,  # list
        "interactions": [],  # list of dictionaries
        "num_interactions": [],
        "num_hb": [],
        "sorted_runs": sorted_runs,
        "pose_about": [],
        "pose_translations": [],
        "pose_quarternions": [],
        "pose_dihedrals": [],
        "fname": fname,
        "ligand_atomtypes": ligand_atomtypes,
    }


def receptor_pdbqt_parser(fname):
    """Parse receptor PDBQT file into list of dictionary with dictionary containing data for a single atom line

    Args:
        fname (string): name of receptor pdbqt file to parse
    """

    with open(fname, "rb") as f:
        raw_lines = f.readlines()

    lines = []
    for line in raw_lines:
        line = line.decode("utf-8")
        if not line.startswith("ATOM") and not line.startswith("HETATM"):
            continue
        line_dict = {}
        # column magic numbers from PDBQT format
        line_dict["atm_flag"] = line[:6]
        line_dict["atm_num"] = int(line[6:11])
        line_dict["atm_name"] = line[12:16]
        line_dict["res_name"] = line[17:20]
        line_dict["chain"] = line[21]
        line_dict["res_num"] = int(line[22:26])
        line_dict["x"] = float(line[30:38])
        line_dict["y"] = float(line[38:46])
        line_dict["z"] = float(line[46:54])
        line_dict["occupancy"] = float(line[54:60])
        line_dict["b_iso"] = float(line[60:66])
        line_dict["q"] = float(line[70:76])
        line_dict["atomtype"] = line[77:]

        lines.append(line_dict)

    return lines
