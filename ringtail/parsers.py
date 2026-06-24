#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail results parsers
#

import os
import gzip
import json
from typing import Union, NamedTuple
from collections import defaultdict
import numpy as np
from .exceptions import (
    FileParsingErrorAdgpu,
    FileParsingErrorPdbqt,
    FileParsingErrorSdf,
)
from .logutils import get_logger

logger = get_logger(__name__)
from rdkit import Chem
from meeko import PDBQTMolecule, RDKitMolCreate
from .interactions import find_interactions, make_interaction_finder


class PoseRecord(NamedTuple):
    receptor: object = None
    pose_rank: object = None
    run_number: object = None
    cluster_rmsd: object = 0
    reference_rmsd: object = 0
    docking_score: object = 0
    leff: object = 0
    delta: object = 0
    energies_inter: object = 0
    energies_vdw: object = 0
    energies_electro: object = 0
    energies_flexLig: object = 0
    energies_flexLR: object = 0
    energies_intra: object = 0
    energies_torsional: object = 0
    unbound_energy: object = 0
    num_interactions: object = 0
    num_hb: object = 0
    cluster_size: object = 1
    pose_coordinates: object = None
    flexible_res_coordinates: object = "[]"
    ligname: object = None


class InteractionRecord(NamedTuple):
    ligname: object
    run_number: object
    pose_rank: object
    type: object
    chain: object
    residue: object
    resid: object
    recname: object
    recid: object


def _dicts_to_interaction_records(interactions: list) -> list:
    """Convert deduplicated interaction dicts to InteractionRecords."""
    return [InteractionRecord(**d) for d in interactions]


def _dict_to_pose_record(d: dict) -> PoseRecord:
    """Convert a pose dict to a PoseRecord, ignoring any extra keys."""
    return PoseRecord(*(d[f] for f in PoseRecord._fields))


def results_row() -> dict:
    """Return a dict of default pose field values, derived from PoseRecord."""
    return PoseRecord()._asdict()


def _open_fn_and_name(fname: str) -> tuple:
    """Return (open_fn, stem_name) for a plain or .gz file."""
    clean = os.path.basename(fname)
    name, ext = os.path.splitext(clean)
    if ext[1:].lower() == "gz":
        name, _ = os.path.splitext(name)
        return gzip.open, name
    return open, name


class VinaMoleculeSupplier:
    def __init__(self, num_poses: int, **kwargs):
        self.kwargs = kwargs
        self.num_poses = num_poses

    def __call__(self, data_pointer):
        yield self.parse_vina_results(data_pointer, **self.kwargs)

    def parse_vina_results(
        self,
        data_pointer,
        is_file=True,
        calculate_interactions: bool = None,
        interaction_finder=None,
        interaction_cutoffs: list[float] = None,
        receptor_string: str = None,
    ) -> dict:
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
                    elif line.startswith("REMARK VINA RESULT:"):
                        scores.append(float(line.split()[3]))
                    elif line.startswith("REMARK INTER:"):
                        intermolecular_energy.append(float(line.split()[2]))
                    elif line.startswith("REMARK INTRA:"):
                        internal_energy.append(float(line.split()[2]))
                    elif line.startswith("REMARK UNBOUND:"):
                        unbound_energy.append(float(line.split()[2]))
                    elif line.startswith("HETATM") or line.startswith("ATOM"):
                        if inside_res:
                            flexible_res_coords[-1][-1].append(
                                [line[30:38], line[38:46], line[46:54]]
                            )
                            if first_model:
                                flexres_atomnames[-1].append(line[12:16].strip())
                        elif first_model and line[13] != "H":
                            num_heavy_atoms += 1
                    elif line.startswith("ENDMDL") and first_model:
                        first_model = False
                    # make new flexible residue list if in the coordinates for a flexible residue
                    elif line.startswith("BEGIN_RES"):
                        flexible_res_coords[-1].append([])
                        inside_res = True
                        if first_model:
                            flexres_atomnames.append([])
                            res = line[10:13].strip()
                            chain = line[14].strip()
                            resnum = line[15:19].strip()
                            flexible_residues.append("%s:%s%s" % (res, chain, resnum))
                    elif line.startswith("END_RES"):
                        inside_res = False
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

        if is_file:
            open_fn, _ = _open_fn_and_name(data_pointer)
            ligname = os.path.basename(data_pointer.split(".pdbqt")[0])
            logger.debug("Parsing vina docking file")
            with open_fn(data_pointer, "rb") as fp:
                pdbqt_str_list = fp.read().decode("utf-8")
            data_object = pdbqt_str_list.splitlines()
        else:
            # input provided as string and convert from '\n' separated to to iterable lines
            pdbqt_str_list = list(data_pointer.values())[0]
            data_object = pdbqt_str_list.splitlines()
            # get the first (and only, probably not optimal) key which should be the ligand name
            ligname = list(data_pointer.keys())[0]
            logger.debug("Parsing vina docking string")

        results_dict, receptor_dict = _read_vina_results_lines(data_object, ligname)

        data_length = len(results_dict["docking_score"])
        # add fields used by other docking engine for db data insert compatibility
        empty_columns = {
            "receptor": ["receptor"] * data_length,
            "cluster_rmsd": [None] * data_length,
            "reference_rmsd": [None] * data_length,
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
            pdbqt_str_list,
            is_dlg=False,
        )
        pose_coordinates = pick_best_poses(coordinates, self.num_poses)

        if calculate_interactions and not interaction_finder:
            interaction_finder = make_interaction_finder(
                receptor_string, interaction_cutoffs
            )

        if interaction_finder:
            pose_identifiers = [
                {"ligname": lig, "run_number": run, "pose_rank": rank}
                for lig, run, rank in zip(
                    results_dict["ligname"],
                    results_dict["run_number"],
                    results_dict["pose_rank"],
                )
            ]
            poseids_coordinates = list(zip(pose_identifiers, pose_coordinates))
            interactions, num_hb, num_interactions = find_interactions(
                poseids_coordinates, mol, interaction_finder=interaction_finder
            )
            interaction_rows = generate_interaction_tuples(interactions)
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
                # plain float lists; the storage manager serializes per dialect
                # (DuckDB FLOAT[][], SQLite packed float32 BLOB)
                "pose_coordinates": [
                    coordinate.tolist() for coordinate in pose_coordinates
                ]
            }
        )
        results_rows = [
            _dict_to_pose_record(dict(zip(results_dict.keys(), values)))
            for values in zip(*results_dict.values())
        ]

        return {
            "ligands": [ligand_row],
            "receptor": generate_receptor_row(receptor_dict),
            "poses": results_rows,
            "interactions": _dicts_to_interaction_records(interaction_rows),
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
                # plain float lists; serialized per dialect by the storage manager
                "pose_coordinates": [
                    coordinate.tolist() for coordinate in poses_coordinates
                ]
            }
        )

        # sort clusters and select based on number of poses to store
        if self.num_poses > -1:
            clusters: dict = results_dict.get("clusters")
            top_cluster_runs = [cluster[0] for cluster in clusters.values()]
            sorted_selected_pose_clusters = pick_best_poses(
                np.argsort(
                    [
                        results_dict.get("docking_score")[run - 1]
                        for run in top_cluster_runs
                    ]
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
                rmsds = results_dict["cluster_rmsd"]
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

            for sorted_idx, cluster_number in enumerate(sorted_selected_pose_clusters):
                pose_rank = sorted_idx + 1
                run_number = cluster_representatives[cluster_number]
                data_index = run_number - 1
                # parse Results data
                results_rows.append(
                    self._build_pose_record(
                        results_dict,
                        data_index,
                        run_number,
                        pose_rank,
                        receptor_name,
                        ligname,
                        cluster_sizes[cluster_number],
                    )
                )

                # parse interactions
                # define unique id for interaction rows
                unique_pose_id = {
                    "ligname": ligname,
                    "run_number": run_number,
                    "pose_rank": pose_rank,
                }
                for interaction_run_number in tolerated_cluster_run_numbers[
                    cluster_number
                ]:
                    current_interactions: dict = interactions[
                        interaction_run_number - 1
                    ].copy()
                    current_interactions.update({"id": unique_pose_id})
                    interaction_dicts.append(current_interactions)
        elif self.num_poses == -1:
            interaction_dicts = []
            results_rows = []
            interactions = results_dict.get("interactions", [])

            sorted_indices_all = np.argsort(results_dict.get("docking_score"))
            # results_dict needs to be made into list of dicts for each index
            for sorted_idx, run_number in enumerate(sorted_indices_all):
                data_index = run_number - 1
                pose_rank = sorted_idx + 1
                # parse Results data
                results_rows.append(
                    self._build_pose_record(
                        results_dict,
                        data_index,
                        run_number,
                        pose_rank,
                        receptor_name,
                        ligname,
                        1,
                    )
                )

                # parse interactions
                # define unique id for interaction rows
                unique_pose_id = {
                    "ligname": ligname,
                    "run_number": run_number,
                    "pose_rank": pose_rank,
                }

                current_interactions = interactions[data_index].copy()
                current_interactions.update({"id": unique_pose_id})
                interaction_dicts.append(current_interactions)

        # prepare data in ringtail recognizable format
        return {
            "ligands": [ligand_row],
            "poses": results_rows,
            "receptor": receptor_row,
            "interactions": _dicts_to_interaction_records(
                generate_interaction_tuples(interaction_dicts)
            ),
        }

    @staticmethod
    def _parse_energy(line: str, idx: int, fname: str) -> float:
        """
        Helper method to parse energies from .dlg strings

        Args:
            line (str): _description_
            idx (int): _description_
            fname (str): _description_

        Raises:
            ValueError: _description_

        Returns:
            float: _description_
        """
        parts = line.split()
        try:
            return float(parts[idx])
        except ValueError:
            try:
                return float(parts[idx - 1].lstrip("="))
            except ValueError:
                raise ValueError(f"ERROR! Cannot parse {line!r} in {fname}")

    def _build_pose_record(
        self,
        results_dict: dict,
        data_index: int,
        run_number: int,
        pose_rank: int,
        receptor_name: str,
        ligname: str,
        cluster_size: int,
    ) -> PoseRecord:
        """
        Convert parsed results to PoseRecord

        Args:
            results_dict (dict): _description_
            data_index (int): _description_
            run_number (int): _description_
            pose_rank (int): _description_
            receptor_name (str): _description_
            ligname (str): _description_
            cluster_size (int): _description_

        Returns:
            PoseRecord: _description_
        """
        return PoseRecord(
            receptor=receptor_name,
            pose_rank=pose_rank,
            run_number=run_number,
            cluster_rmsd=results_dict["cluster_rmsd"][data_index],
            reference_rmsd=results_dict["ref_rmsd"][data_index],
            docking_score=results_dict["docking_score"][data_index],
            leff=results_dict["leff"][data_index],
            delta=results_dict["delta"][data_index],
            energies_inter=results_dict["energies_inter"][data_index],
            energies_vdw=results_dict["energies_vdw"][data_index],
            energies_electro=results_dict["energies_electro"][data_index],
            energies_flexLig=results_dict["energies_flexLig"][data_index],
            energies_flexLR=results_dict["energies_flexLR"][data_index],
            energies_intra=results_dict["energies_intra"][data_index],
            energies_torsional=results_dict["energies_torsional"][data_index],
            unbound_energy=results_dict["unbound_energy"][data_index],
            num_interactions=results_dict["num_interactions"][data_index],
            num_hb=results_dict["num_hb"][data_index],
            cluster_size=cluster_size,
            pose_coordinates=results_dict["pose_coordinates"][data_index],
            flexible_res_coordinates=json.dumps(
                results_dict["flexible_res_coordinates"][data_index]
            ),
            ligname=ligname,
        )

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

        open_fn, name = _open_fn_and_name(fname)

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
        flexres_atomnames = []

        # Define empty center list for backwards compatibility with DLGs without grid centers
        center = [None, None, None]

        heavy_at_count = 0
        heavy_at_count_complete = False
        # read file contents
        with open_fn(fname, "rb") as fp:
            file_as_string = fp.read().decode("utf-8")

        # divide
        header_str, sep, body_str = file_as_string.partition("FINAL DOCKED STATE")

        # parse header
        for line in header_str.splitlines():
            if line[0:11] == "Ligand file":
                ligname = (
                    line.split(":", 1)[1].strip().replace("\\", "/").split("/")[-1].split(".")[0]
                )  # remove path and file extension (separator-agnostic)
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
            # store flexible residue identities and atomtypes
            if "INPUT-FLEXRES-PDBQT:" in line:
                if "ATOM" in line or "HETATM" in line:
                    flexres_atomnames[-1].append(line[33:37])
                    if (
                        line[38:41] + ":" + line[42] + line[44:47]
                        not in flexible_residues
                    ):
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

        # energy to variable table with nan ok/not ok
        _ENERGY_SPECS = [
            ("Estimated Free Energy of Binding", 7, scores, True),
            ("Final Intermolecular Energy", 6, intermolecular_energy, False),
            ("vdW + Hbond + desolv Energy", 8, vdw_hb_desolv, False),
            ("Electrostatic Energy", 4, electrostatic, False),
            ("Moving Ligand-Fixed Receptor", 5, flex_ligand, False),
            ("Moving Ligand-Moving Receptor", 5, flexLigand_flexReceptor, False),
            ("Final Total Internal Energy", 7, internal_energy, False),
            ("Torsional Free Energy", 6, torsion, False),
            ("Unbound System's Energy", 6, unbound_energy, False),
        ]

        # parse poses
        inside_pose = False
        inside_res = False
        for line in (sep + body_str).splitlines():
            if "FINAL DOCKED STATE" in line:
                inside_pose = True
                interactions.append({})
                flexible_res_coords.append([])
            # store pose analysis
            elif line[0:9] == "ANALYSIS:":
                line = line.split("ANALYSIS:")[1]
                kw, info = line.split(None, 1)
                info = info.replace("{", "").replace("}", "").replace('"', "")
                if kw.lower() == "count":
                    interactions[-1][kw.lower()] = int(
                        [x.strip() for x in info.split(",")][0]
                    )
                else:
                    interactions[-1][kw.lower()] = [x.strip() for x in info.split(",")]
                if "COUNT" in line:
                    interact_count = int(line.split()[1])
                    pose_interact_count.append(str(interact_count))
                    if interact_count == 0:
                        pose_hb_counts.append("0")
                elif "TYPE" in line:
                    pose_hb_counts.append(line.count("H"))

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
                    for marker, idx, target, check_nan in _ENERGY_SPECS:
                        if marker in line:
                            e = self._parse_energy(line, idx, fname)
                            if check_nan and np.isnan(e):
                                raise ValueError(
                                    "Error! File contains NaN value for energy."
                                )
                            target.append(e)
                            break
                # update heavy atom count
                if not heavy_at_count_complete:
                    if line[0:4] == "ATOM" or line[0:6] == "HETATM":
                        if not line[-2] == "HD":
                            heavy_at_count += 1

            # store poses in each cluster in dictionary as list of ordered runs
            elif "RANKING" in line:
                cluster_num = int(line.split()[0]) - 1  # make it zero based integer
                run = line.split()[2]
                clusters[cluster_num].append(int(run))
                cluster_rmsds[int(run) - 1] = float(
                    line.split()[4]
                )  # will be stored in order of runs
                ref_rmsds[int(run) - 1] = float(line.split()[5])

        # ensure data is complete
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
            "docking_score": scores,
            "leff": [round(x / heavy_at_count, 2) for x in scores],
            "delta": [round(x - scores[0], 2) for x in scores],
            "energies_inter": intermolecular_energy,
            "energies_vdw": vdw_hb_desolv,
            "energies_electro": electrostatic,
            "energies_flexLig": flex_ligand,
            "energies_flexLR": flexLigand_flexReceptor,
            "energies_intra": internal_energy,
            "energies_torsional": torsion,
            "unbound_energy": unbound_energy,
            "num_interactions": pose_interact_count,
            "num_hb": pose_hb_counts,
            "flexible_res_coordinates": flexible_res_coords,
            "cluster_rmsd": cluster_rmsds,
            "ref_rmsd": ref_rmsds,
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
        self.num_poses = num_poses

    def __call__(self, fname: str):
        return self.parse_adng_results(fname, **self.kwargs)

    def parse_adng_results(
        self,
        fname: str,
        calculate_interactions: bool = None,
        interaction_finder=None,
        interaction_cutoffs: list[float] = None,
        receptor_string: str = None,
    ):
        """
        Method to parse docking results from SDF files, uses Chem.ForwardSDMolSupplier

        Args:
            fname (str): sdf file name
            calculate_interactions (bool, optional): _description_. Defaults to None.
            interaction_finder (_type_, optional): _description_. Defaults to None.
            interaction_cutoffs (list[float], optional): _description_. Defaults to None.
            receptor_string (str, optional): _description_. Defaults to None.

        Yields:
            _type_: _description_
        """
        if calculate_interactions and not interaction_finder:
            interaction_finder = make_interaction_finder(
                receptor_string, interaction_cutoffs
            )

        processed_ligands = set()
        last_ligand = None

        with open(fname, "rb") as file_stream:
            suppl = Chem.ForwardSDMolSupplier(file_stream, removeHs=False)
            for mol in suppl:
                if mol is None:
                    continue
                ligname = mol.GetProp("_Name")
                if ligname != last_ligand:
                    last_ligand = ligname
                    pose_count = 1
                else:
                    pose_count += 1
                if self.num_poses != -1 and pose_count > self.num_poses:
                    continue

                try:
                    result = process_docked_mol(
                        mol,
                        calculate_interactions=bool(interaction_finder),
                        interaction_finder=interaction_finder,
                    )
                except Exception as e:
                    raise FileParsingErrorSdf(
                        f"Failed to process mol '{ligname}' in {fname}: {e}"
                    ) from e
                # Only emit the ligand row on first occurrence — subsequent poses
                # of the same ligand rely on ON CONFLICT IGNORE in the DB, but
                # avoid the smiles/binary cost by clearing the list here.
                if ligname in processed_ligands:
                    result["ligands"] = []
                else:
                    processed_ligands.add(ligname)

                yield result


def process_docked_mol(
    mol,
    calculate_interactions: bool,
    interaction_finder=None,
    add_default_columns: bool = True,
):
    """
    Processes a single docked Mol, which includes removing certain properties and preserving others,
    as well as parsing some like coordinates, as needed for the Ringtail database

    Args:
        mol (_type_): _description_
        calculate_interactions (bool): _description_
        interaction_finder (_type_, optional): _description_. Defaults to None.
        add_default_columns (bool, optional): _description_. Defaults to True.

    Returns:
        _type_: _description_
    """
    ligname = mol.GetProp("_Name")
    run_number = 1
    smiles = Chem.MolToSmiles(mol)
    interaction_rows = []
    results_dict = results_row() if add_default_columns else {}
    n_heavy = mol.GetNumHeavyAtoms()
    mol, mol_properties = prepare_mol_for_database(
        mol, store_properties=["docking_score", "pose_rank"]
    )
    ligand_row = [
        ligname,
        smiles,
        mol.ToBinary(
            Chem.PropertyPickleOptions.MolProps
            | Chem.PropertyPickleOptions.AtomProps
            | Chem.PropertyPickleOptions.BondProps
        ),
    ]

    results_dict.update(mol_properties)
    results_dict.update(
        {
            "ligname": ligname,
            "run_number": run_number,
            "leff": round(results_dict["docking_score"] / n_heavy, 2),
        }
    )
    single_pose_coordinate = results_dict["pose_coordinates"][0]
    results_dict["pose_coordinates"] = (
        single_pose_coordinate.tolist()
        if isinstance(single_pose_coordinate, np.ndarray)
        else single_pose_coordinate
    )

    if calculate_interactions and interaction_finder:
        poseids_coordinates = [
            (
                {"ligname": ligname, "run_number": run_number, "pose_rank": 1},
                single_pose_coordinate,
            )
        ]
        interactions, num_hb, num_interactions = find_interactions(
            poseids_coordinates, mol, interaction_finder=interaction_finder
        )
        interaction_rows = generate_interaction_tuples(interactions)
        results_dict.update(
            {"num_interactions": num_interactions[0], "num_hb": num_hb[0]}
        )

    return {
        "ligands": [ligand_row],
        "poses": [_dict_to_pose_record(results_dict)],
        "interactions": _dicts_to_interaction_records(interaction_rows),
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
    rdkit_mol = RDKitMolCreate.from_pdbqt_mol(pdbqt_mol, only_hs_with_coords=True)[0]

    rdkit_mol, properties = prepare_mol_for_database(rdkit_mol)
    pose_coordinates = properties.get("pose_coordinates")
    smiles = Chem.MolToSmiles(rdkit_mol)
    ligand_rdbin = rdkit_mol.ToBinary(
        Chem.PropertyPickleOptions.MolProps
        | Chem.PropertyPickleOptions.AtomProps
        | Chem.PropertyPickleOptions.BondProps
    )

    ligand_row = [ligname, smiles, ligand_rdbin]

    return (ligand_row, pose_coordinates, rdkit_mol)


def prepare_mol_for_database(
    mol: Chem.Mol, store_properties: list = None
) -> tuple[Chem.Mol, dict]:
    """
    In this method right now the pose coordinates assumes one or more conformers,
    while the property assumes only one mol with non-list property

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
                properties.update({db_column: None})
    return mol, properties


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
    else:
        grid_spacing = None
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


def generate_interaction_tuples(interaction_dictionaries: list) -> list:
    """Format pose interaction dicts into deduplicated InteractionRecord-ready dicts.

    Each entry includes the pose identity keys (ligname, run_number, pose_rank)
    plus the interaction descriptor fields.

    Args:
        interaction_dictionaries (list): List of pose interaction dicts from parser.

    Returns:
        list: Deduplicated dicts ready to be converted to InteractionRecords.
    """
    if not interaction_dictionaries:
        return []

    interaction_keywords = ["type", "chain", "residue", "resid", "recname", "recid"]
    id_keys = list(interaction_dictionaries[0]["id"].keys())
    all_keys = id_keys + interaction_keywords

    tuples = dict.fromkeys(
        tuple(pose["id"][k] if k in id_keys else pose[k][i] for k in all_keys)
        for pose in interaction_dictionaries
        for i in range(pose["count"])
    )

    return [dict(zip(all_keys, t)) for t in tuples]


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
