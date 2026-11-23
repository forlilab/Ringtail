#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail types
#


import os
from typing import NamedTuple
from dataclasses import dataclass
from rdkit import Chem


class PoseMol(NamedTuple):
    """One RDKit mol per pose, returned by RingtailCore.create_rdkit_mols.

    flexres_mols is [] unless flexres_data was passed to the call.
    """

    pose_id: int
    mol: Chem.Mol
    ligname: str
    flexres_mols: list
    properties: dict


class LigandMol(NamedTuple):
    """One multi-conformer RDKit mol per ligand, from create_rdkit_mols_by_ligand."""

    ligname: str
    mol: object
    flexres_per_pose: list
    properties: dict
    flexres_residues: list


class PoseData(NamedTuple):
    """Row returned by StorageManager.fetch_rdkit_pose_properties.

    Optional fields are populated only when the corresponding flag is passed to
    the fetch: flexres_coordinates is None unless include_flexres=True, and
    comment is '' unless include_comment=True (and the Pose_comments table exists).
    """

    pose_id: int
    rdmol: object
    ligname: str
    docking_score: object
    leff: object
    pose_coordinates: object
    pose_rank: object
    flexres_coordinates: object = None
    comment: str = ""


@dataclass
class ResultsObject:
    """Class that handles sources of data to be written including ligand data paths and how
    to traverse them, and options to store receptor.
    """

    docking_results: object = None
    recursive_path_traverse: bool | None = None
    receptor_file_path: str | None = None
    save_receptor: bool | None = None
    receptor_string: str | None = None

    @property
    def target_name(self):
        if self.receptor_file_path and os.path.exists(self.receptor_file_path):
            return os.path.basename(self.receptor_file_path).split(".")[0]
        return None

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
    energies_flexlig: object = 0
    energies_flexlr: object = 0
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


def dicts_to_interaction_records(interactions: list) -> list:
    """Convert deduplicated interaction dicts to InteractionRecords."""
    return [InteractionRecord(**d) for d in interactions]


def dict_to_pose_record(d: dict) -> PoseRecord:
    """Convert a pose dict to a PoseRecord, ignoring any extra keys."""
    return PoseRecord(*(d[f] for f in PoseRecord._fields))


def results_row() -> dict:
    """Return a dict of default pose field values, derived from PoseRecord."""
    return PoseRecord()._asdict()
