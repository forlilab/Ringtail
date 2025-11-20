#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail results parsers
#

import os
import gzip
import bz2
import json
from typing import Union
from collections import defaultdict
import numpy as np
from .exceptions import (
    FileParsingErrorAdgpu,
    FileParsingErrorPdbqt,
    FileParsingErrorSdf,
)
from .logutils import LOGGER
from rdkit import Chem, Geometry
from meeko import PDBQTMolecule, RDKitMolCreate, MoleculePreparation
from .interactions import calculate_interactions


def results_row() -> dict:
    """
    Returns a dict with expected keys and None values for inserting into ringtail

    Returns:
        dict:
    """
    return {
        "receptor": None,
        "pose_rank": None,
        "run_number": None,
        "cluster_rmsd": 0,
        "reference_rmsd": 0,
        "docking_score": 0,
        "leff": 0,
        "deltas": 0,
        "energies_inter": 0,
        "energies_vdw": 0,
        "energies_electro": 0,
        "energies_flexLig": 0,
        "energies_flexLR": 0,
        "energies_intra": 0,
        "energies_torsional": 0,
        "unbound_energy": 0,
        "num_interactions": 0,
        "num_hb": 0,
        "cluster_size": 1,
        "pose_coordinates": None,
        "flexible_res_coordinates": None,
        "ligname": None,
    }


RESULTS_TEMPLATE = results_row()


def make_ringtail_data_dict(
    ligand_row: list,
    receptor_row: list,
    data_rows: list,
    interaction_rows: list,
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


class VinaMoleculeSupplier:
    def __init__(self, num_poses: int, **kwargs):
        self.kwargs = kwargs
        self.num_poses = num_poses

    def __call__(self, data_pointer):
        yield self.parse_vina_results(data_pointer)

    def parse_vina_results(self, data_pointer) -> dict:
        """Parser for vina docking results, supporting either pdbqt or gzipped (.gz) files, or with the
        docking results provided as a string.

        Args:
            data_pointer (any): either filename or dictionary of string docking results

        Returns:
            dict: parsed results ready to be inserted in database
        """

        def _read_vina_results_lines(data_object, name) -> tuple[dict, dict, dict]:
            """Reads vina docking results line by line

            Args:
                data_object (any): filename or dict containing docking results
                name (str): given ligand name/file name

            Raises:
                ValueError: if a line cannot be parsed

            Returns:
                dict: parsed results ready to be inserted in database
            """
            scores = []
            run_number = []
            num_heavy_atoms = 0
            intermolecular_energy = []
            internal_energy = []
            unbound_energy = []
            num_heavy_atoms = 0
            flexible_res_coords = []
            inside_res = False
            flexible_residues = []
            flexres_atomnames = []
            first_model = True

            for line in data_object:
                try:
                    if line.startswith("MODEL"):
                        flexible_res_coords.append([])
                        run_number.append(line.split()[1])
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
                            if first_model:
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
                    raise FileParsingErrorPdbqt(
                        "ERROR! Cannot parse {0} in {1}".format(line, name)
                    )

            # calculate ligand efficiency and deltas from the best pose
            leff = [round(x / num_heavy_atoms, 2) for x in scores]
            delta = [round(x - scores[0], 2) for x in scores]

            receptor_dict = {
                "flexible_residues": flexible_residues,
                "flexres_atomnames": flexres_atomnames,
            }

            results = {
                "pose_rank": [1] * len(scores),  # list
                "run_number": run_number,  # list
                "flexible_res_coordinates": [
                    json.dumps(coordinate) for coordinate in flexible_res_coords
                ],  # list of lists
                "docking_score": scores,  # list
                "leff": leff,  # list
                "delta": delta,  # list
                "energies_inter": intermolecular_energy,  # list
                "energies_intra": internal_energy,  # list
                "unbound_energy": unbound_energy,  # list
            }

            # pick desired number of poses
            selected_results = {
                key: (
                    pick_best_poses(value, self.num_poses)
                    if isinstance(value, list)
                    else value
                )
                for key, value in results.items()
            }
            return selected_results, receptor_dict

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
                pdbqt_str = "".join(data_object)
        else:
            # input provided as string and convert from '\n' separated to to iterable lines
            data_object = list(data_pointer.values())[0].splitlines()
            pdbqt_str = data_pointer
            # get the first (and only, probably not optimal) key which should be the ligand name
            ligname = list(data_pointer.keys())[0]
            LOGGER.debug("Parsing vina docking string")

        results_dict, receptor_dict = _read_vina_results_lines(data_object, ligname)

        data_length = len(results_dict["docking_score"])
        # add fields used by other docking engine for db data insert compatibility
        empty_columns = {
            "receptor": ["receptor"] * data_length,
            "cluster_rmsd": [None] * data_length,
            "deltas": [None] * data_length,
            "reference_rmsd": [None] * data_length,
            "energies_inter": [None] * data_length,
            "energies_vdw": [None] * data_length,
            "energies_electro": [None] * data_length,
            "energies_flexLig": [None] * data_length,
            "energies_flexLR": [None] * data_length,
            "energies_torsional": [None] * data_length,
            # these two are calculated with interactions
            "cluster_size": [1] * data_length,
            "ligname": [ligname] * data_length,
        }

        # add empty columns and interaction data
        results_dict.update(empty_columns)
        ligand_row, coordinates, mol = generate_ligand_data_list_from_pdbqt_dlg(
            ligname,
            pdbqt_str,
            is_dlg=False,
        )
        pose_coordinates = pick_best_poses(coordinates, self.num_poses)
        if (
            "calculate_interactions" in self.kwargs
            and self.kwargs["calculate_interactions"] == True
        ):
            try:
                receptor_string = self.kwargs["receptor_string"]
            except:
                LOGGER.error(
                    "Cannot calculate interactions, missing receptor representation. Interactions can be calculated at a later time if receptor is provided."
                )
                interaction_info = {
                    "nr_interactions": [0] * data_length,
                    "num_hb": [0] * data_length,
                }
                interactions = []
            interaction_fields = {
                k: v
                for k, v in results_dict.items()
                if k in ["ligname", "run_number", "pose_rank"]
            }
            interaction_fields.update({"pose_coordinates": pose_coordinates})

            interactionless_poses = [
                dict(zip(interaction_fields.keys(), values))
                for values in zip(*interaction_fields.values())
            ]
            """
            Each value is a list, and for each value I need to zip item n with other item ns into
            a dict with key 
            """
            interactions, num_hb, num_interactions = calculate_interactions(
                interactionless_poses,
                mol,
                receptor_string,
                *self.kwargs["interaction_cutoffs"],
            )
            interaction_rows = generate_interaction_tuples(interactions)
            # use interaction stuff here
            interaction_info = {
                "num_interactions": num_interactions,
                "num_hb": num_hb,
            }
        else:
            interaction_info = {
                "num_interactions": [0] * data_length,
                "num_hb": [0] * data_length,
            }
            interaction_rows = []

        results_dict.update(interaction_info)
        results_dict.update(
            {
                "pose_coordinates": [
                    json.dumps(coordinate.tolist()) for coordinate in pose_coordinates
                ]
            }
        )
        # zip all the rows together so there is one dict per row/pose
        results_rows = [
            dict(zip(results_dict.keys(), values))
            for values in zip(*results_dict.values())
        ]

        return {
            "ligands": [ligand_row],
            "receptor": generate_receptor_row(receptor_dict),
            "poses": results_rows,
            "interactions": interaction_rows,
        }


class ADGPUMoleculeSupplier:
    def __init__(self, num_poses: int, **kwargs):
        self.kwargs = kwargs
        self.num_poses = num_poses

    def __call__(self, filename: str):
        yield self.parse_adgpu_results(filename)

    def parse_adgpu_results(self, filename: str) -> dict:
        """
        Parses an ADGPU dlg (docking log file) and returns dict of lists of rows ready to be inserted in database

        Args:
            filename (str): dlg file name to parse

        Returns:
            dict: containing a
                ligand row (list of ligand name, smiles, and binary rdkit mol)
                receptor row (list of receptor data)
                resulst rows* (list of docking results such as energies, coordinates, and other relevant per pose data)
                interaction rows (list of interaction tuples that includes ligname, pose rank, run number to uniqueliy identify
                    to which pose an interaction belongs, which makes the interaction list independent of other results
                    when it comes to e.g., inserting data in the database)
        """

        ligand_dict, receptor_dict, results_dict = self._parse_docking_file_dlg(
            fname=filename
        )

        # make sure receptor is correct
        receptor_row = generate_receptor_row(receptor_dict)
        receptor_name = receptor_row[0]
        if target := self.kwargs.get("target"):
            if receptor_name != target:
                raise FileParsingErrorAdgpu(
                    f"Receptor name {receptor_name} in {filename} does not match given target name {target}. Please ensure that this file belongs to the current virtual screening."
                )

        # create rdkit mol from file, get coordinates for each pose
        ligand_row, poses_coordinates, _ = generate_ligand_data_list_from_pdbqt_dlg(
            ligand_dict["name"], ligand_dict["file_str"], is_dlg=True
        )
        ligname = ligand_row[0]
        results_dict.update(
            {
                "pose_coordinates": [
                    json.dumps(coordinate.tolist()) for coordinate in poses_coordinates
                ]
            }
        )

        # sort clusters and select based on number of poses to store
        clusters: dict = results_dict.get("clusters")
        top_cluster_runs = [cluster[0] for cluster in clusters.values()]
        sorted_selected_pose_clusters = pick_best_poses(
            np.argsort(
                [results_dict.get("scores")[run - 1] for run in top_cluster_runs]
            ),
            self.num_poses,
        )

        # delete clusters not considered after filtering for num_poses
        clusters = {
            k: v for k, v in clusters.items() if k in sorted_selected_pose_clusters
        }

        cluster_sizes = {k: len(v) for k, v in clusters.items()}
        # remove runs/clustered poses that we are not storing
        cluster_representatives = {k: v[0] for k, v in clusters.items()}

        # deal with request to include more interactions for poses that are clustered with "main" poses
        if isinstance(self.kwargs["interaction_tolerance"], float):
            int_tol = self.kwargs["interaction_tolerance"]
            rmsds = results_dict["cluster_rmsds"]
            cluster_rmsd_dict = {
                k: [rmsds[run - 1] for run in v] for k, v in clusters.items()
            }
            tolerated_cluster_run_numbers = defaultdict(list)
            # go through the cluster rmsd dict and determine which indices pass muster
            for cluster, rmsds in cluster_rmsd_dict.items():
                for index, rmsd in enumerate(rmsds):
                    if rmsd <= int_tol:
                        tolerated_cluster_run_numbers[cluster].append(
                            clusters[cluster][index]
                        )
        else:
            tolerated_cluster_run_numbers = {
                k: [v] for k, v in cluster_representatives.items()
            }

        interaction_dicts = []
        results_rows = []
        interactions = results_dict.get("interactions", [])

        for pose_rank, cluster_number in enumerate(sorted_selected_pose_clusters):
            run_number = cluster_representatives[cluster_number]
            data_index = run_number - 1
            # parse Results data
            results_rows.append(
                {
                    "receptor": receptor_name,
                    "pose_rank": pose_rank,
                    "run_number": run_number,
                    "cluster_rmsd": results_dict["cluster_rmsds"][data_index],
                    "reference_rmsd": results_dict["ref_rmsds"][data_index],
                    "docking_score": results_dict["scores"][data_index],
                    "leff": results_dict["leff"][data_index],
                    "deltas": results_dict["delta"][data_index],
                    "energies_inter": results_dict["intermolecular_energy"][data_index],
                    "energies_vdw": results_dict["vdw_hb_desolv"][data_index],
                    "energies_electro": results_dict["electrostatics"][data_index],
                    "energies_flexLig": results_dict["flex_ligand"][data_index],
                    "energies_flexLR": results_dict["flexLigand_flexReceptor"][
                        data_index
                    ],
                    "energies_intra": results_dict["internal_energy"][data_index],
                    "energies_torsional": results_dict["torsional_energy"][data_index],
                    "unbound_energy": results_dict["unbound_energy"][data_index],
                    "num_interactions": results_dict["num_interactions"][data_index],
                    "num_hb": results_dict["num_hb"][data_index],
                    "cluster_size": cluster_sizes[cluster_number],
                    "pose_coordinates": results_dict["pose_coordinates"][data_index],
                    "flexible_res_coordinates": json.dumps(
                        results_dict["flexible_res_coordinates"][data_index]
                    ),
                    "ligname": ligname,
                }
            )

            # parse interactions
            # define unique id for interaction rows
            unique_pose_id = {
                "ligand_name": ligname,
                "run_number": run_number,
                "pose_rank": pose_rank,
            }
            for interaction_run_number in tolerated_cluster_run_numbers[cluster_number]:
                current_interactions: dict = interactions[
                    interaction_run_number - 1
                ].copy()
                current_interactions.update(unique_pose_id)
                interaction_dicts.append(current_interactions)

        # prepare data in ringtail recognizable format
        return {
            "ligands": [ligand_row],
            "poses": results_rows,
            "receptor": receptor_row,
            "interactions": generate_interaction_tuples(interaction_dicts),
        }

    def _parse_docking_file_dlg(self, fname: str) -> tuple[dict, dict, dict]:
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
        clusters = defaultdict(list)
        pose_interact_count = []
        pose_hb_counts = []
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
                    flexible_res_coords.append([])
                # store pose anaylsis
                elif line[0:9] == "ANALYSIS:":
                    # storing interactions
                    line = line.split("ANALYSIS:")[1]
                    kw, info = line.split(None, 1)
                    info = info.replace("{", "")
                    info = info.replace("}", "")
                    info = info.replace('"', "")
                    if kw.lower() == "count":
                        interactions[-1][kw.lower()] = int(
                            [x.strip() for x in info.split(",")][0]
                        )
                    else:
                        interactions[-1][kw.lower()] = [
                            x.strip() for x in info.split(",")
                        ]

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
                    cluster_num = int(line.split()[0]) - 1  # make it zero based integer
                    run = line.split()[2]
                    clusters[cluster_num].append(int(run))

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
            raise FileParsingErrorAdgpu("Incomplete data in " + fname)

        results = {
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
            "flexible_res_coordinates": flexible_res_coords,
            "cluster_rmsds": cluster_rmsds,
            "ref_rmsds": ref_rmsds,
            "interactions": interactions,
            "clusters": clusters,
        }

        receptor_dict = {
            "receptor": receptor,
            "grid_center": center,
            "grid_dim": npts,
            "grid_spacing": spacing,
            "flexible_residues": flexible_residues,
            "flexres_atomnames": flexres_atomnames,
        }

        ligand_dict = {"name": ligname, "file_str": file_as_string}

        return ligand_dict, receptor_dict, results


class SDFMoleculeSupplier:
    def __init__(self, num_poses: int, **kwargs):
        self.kwargs = kwargs
        if num_poses == -1:
            self.num_poses = float("inf")
        else:
            self.num_poses = num_poses

        if (
            "calculate_interactions" in self.kwargs
            and self.kwargs["calculate_interactions"] == True
        ):
            try:
                self.receptor_string = self.kwargs["receptor_string"]
            except:
                LOGGER.error(
                    "Cannot calculate interactions, missing receptor representation. "
                )
                self.calc_interactions = False
            else:
                self.calc_interactions = True

    def __call__(self, fname: str):
        # keep track of parsed ligands
        processed_ligands = set()
        last_ligand = None

        file_stream = open(fname, "rb")
        suppl = Chem.ForwardSDMolSupplier(file_stream, removeHs=False)
        for mol in suppl:
            if mol is None:
                continue
            ligname = mol.GetProp("_Name")
            # only grab num_poses
            if ligname != last_ligand:
                last_ligand = ligname
                pose_count = 1
            else:
                pose_count += 1
            if pose_count > self.num_poses:
                continue

            results_dict = RESULTS_TEMPLATE.copy()
            # remove coordinates and docking properties (retains other custom properties)
            mol, mol_properties = prepare_mol_for_database(
                mol, store_properties=["ligname", "docking_score", "pose_rank"]
            )
            # ligname = mol_properties.get("ligname", None)

            if ligname not in processed_ligands:
                smiles = Chem.MolToSmiles(mol)
                processed_ligands.add(ligname)
                ligand_row = [ligname, smiles, mol.ToBinary()]

            results_dict.update(mol_properties)
            results_dict.update(
                {
                    "ligname": ligname,
                    "run_number": 1,
                    "leff": round(results_dict["docking_score"] / mol.GetNumAtoms(), 2),
                }
            )
            single_pose_coordinate = results_dict["pose_coordinates"][0]
            results_dict.update({"pose_coordinates": single_pose_coordinate})
            # calculate interactions
            if self.calc_interactions:
                interaction_fields = {
                    key: results_dict[key]
                    for key in [
                        "ligname",
                        "run_number",
                        "pose_rank",
                        "pose_coordinates",
                    ]
                }

                interactions, num_hb, num_interactions = calculate_interactions(
                    [interaction_fields],
                    mol,
                    self.receptor_string,
                    *self.kwargs["interaction_cutoffs"],
                )
                interaction_rows = generate_interaction_tuples(interactions)
                # use interaction stuff here
                results_dict.update(
                    {
                        "num_interactions": num_interactions[0],
                        "num_hb": num_hb[0],
                        "pose_coordinates": json.dumps(
                            results_dict["pose_coordinates"].tolist()
                        ),
                    }
                )

            # add non-specified fields as list of None
            yield {
                "ligands": [ligand_row],
                "poses": [results_dict],
                "interactions": interaction_rows,
                "receptor": [],
            }


def process_docked_mol(mol, **kwargs):
    """
    #TODO part of the poin t of this method is that it takes one mol with perhaps multiple
    conformers, it is not one mol one pose as in the file reading
    #TODO maybe make this some sort of streaming service?
    """

    ligname = mol.GetProp("_Name")
    smiles = Chem.MolToSmiles(mol)
    ligand_row = [ligname, smiles, mol.ToBinary()]
    interaction_rows = []

    results_dict = RESULTS_TEMPLATE.copy()
    # remove coordinates and docking properties (retains other custom properties)
    # I think conformers are properly parsed but not sure about scores
    mol, mol_properties = prepare_mol_for_database(
        mol, store_properties=["docking_score", "pose_rank"]
    )

    results_dict.update(mol_properties)
    results_dict.update(
        {
            "ligname": ligname,
            "run_number": 1,
            "leff": round(results_dict["docking_score"] / mol.GetNumAtoms(), 2),
        }
    )
    single_pose_coordinate = results_dict["pose_coordinates"][0]
    results_dict.update({"pose_coordinates": single_pose_coordinate})
    # calculate interactions
    if "calculate_interactions" in kwargs and kwargs["calculate_interactions"] == True:
        interaction_fields = {
            key: results_dict[key]
            for key in [
                "ligname",
                "run_number",
                "pose_rank",
                "pose_coordinates",
            ]
        }

        interactions, num_hb, num_interactions = calculate_interactions(
            [interaction_fields],
            mol,
            kwargs["receptor_string"],
            *kwargs["interaction_cutoffs"],
        )
        interaction_rows = generate_interaction_tuples(interactions)
        # use interaction stuff here
        results_dict.update(
            {
                "num_interactions": num_interactions[0],
                "num_hb": num_hb[0],
                "pose_coordinates": json.dumps(
                    results_dict["pose_coordinates"].tolist()
                ),
            }
        )

    # add non-specified fields as list of None
    return {
        "ligands": [ligand_row],
        "poses": [results_dict],
        "interactions": interaction_rows,
        "receptor": [],
    }


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


def generate_ligand_data_list_from_pdbqt_dlg(
    ligname: str, file_str: str, is_dlg=True
) -> tuple[list, list, Union[Chem.Mol, None]]:
    """writes row to be inserted into ligand table

    Args:
        ligand_dict (dict): Dictionary of ligand data from parser

    Returns:
        list: List of data to be written as row in ligand table. Format:
            [ligand_name, ligand_rdbin]
    """

    # use meeko prepare molecule
    pdbqt_mol = PDBQTMolecule(file_str, name=ligname, is_dlg=is_dlg, skip_typing=True)
    # return the whole list with conformers if requested
    rdkit_mol = RDKitMolCreate.from_pdbqt_mol(pdbqt_mol)[0]

    rdkit_mol, properties = prepare_mol_for_database(rdkit_mol)
    pose_coordinates = properties.get("pose_coordinates")
    smiles = Chem.MolToSmiles(rdkit_mol)
    ligand_rdbin = rdkit_mol.ToBinary()

    ligand_row = [ligname, smiles, ligand_rdbin]

    return (ligand_row, pose_coordinates, rdkit_mol)


def prepare_mol_for_database(
    mol: Chem.Mol, store_properties: list = None
) -> tuple[Chem.Mol, dict]:
    """
    In this method right now the pose coordinates assumes one or more conformers,
    while the property assumes only one mol with non-list property
    #TODO probably make two different methods

    Args:
        mol (Chem.Mol): _description_
        store_properties (dict, optional): _description_. Defaults to None.

    Returns:
        tuple[Chem.Mol, dict]: _description_
    """
    pose_coordinates = []
    for conf in mol.GetConformers():
        pose_coordinates.append(conf.GetPositions())
    # remove conformer data and store one binary version of the rdkit
    mol.RemoveAllConformers()
    # Get a list of property names
    properties = {"pose_coordinates": pose_coordinates}
    if store_properties:
        for db_column in store_properties:
            prop_name = db_to_adng(db_column)
            if mol.HasProp(prop_name):
                properties.update(
                    {db_column: type_casting[db_column](mol.GetProp(prop_name))}
                )
                mol.ClearProp(prop_name)
            else:
                # TODO raise error if missing ligand name?
                properties.update({db_column: None})
    return mol, properties


def adng_to_db(adng_property: str) -> str:
    return adng_aliases.get(adng_property, adng_property)


def db_to_adng(db_column: str) -> str:
    return db_alias_to_adng.get(db_column, db_column)


def generate_receptor_row(receptor_data: dict) -> list:
    """Writes row to be inserted into receptor table

    Args:
        ligand_dict (dict): Dictionary of ligand data from parser

    Returns:
        list: receptor row columns
    """
    rec_name = receptor_data.get("receptor", None)
    box_dim = json.dumps(receptor_data.get("grid_dim", None))
    box_center = json.dumps(receptor_data.get("grid_center", None))
    grid_spacing = receptor_data.get("grid_spacing", "")
    if grid_spacing != "":
        grid_spacing = float(grid_spacing)
    flexible_residues = json.dumps(receptor_data.get("flexible_residues", None))
    flexres_atomnames = json.dumps(receptor_data.get("flexres_atomnames", None))

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
        for i in range(pose["count"])
    }

    return list(interactions)


docking_file_parsers: dict[str, callable] = {
    "adgpu": ADGPUMoleculeSupplier,
    "vina": VinaMoleculeSupplier,
    "adng": SDFMoleculeSupplier,
}


adng_aliases = {
    "adng_free_energy": "docking_score",
    "pose_id": "pose_rank",
    "_Name": "ligname",
}
type_casting = {"docking_score": float, "pose_rank": int, "ligname": str}
db_alias_to_adng = {v: k for k, v in adng_aliases.items()}
