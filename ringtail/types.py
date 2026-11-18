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
