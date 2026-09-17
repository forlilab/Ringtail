#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail interaction finder
#

import json

import numpy as np
from typing import Union
from meeko import MoleculePreparation, Polymer
from .receptormanager import (
    receptor_atoms_from_pdbqt_string,
    receptor_atoms_from_polymer
)
from rdkit import Chem
from .logutils import get_logger

logger = get_logger(__name__)



# Ligand H-bond donor/acceptor SMARTS, from meeko's vina_params.json.
# Ordered cascade, last match wins, so file order must be preserved.
# Flags land on match[0] -- for the donor rules ("[#8][#1]") that is the heavy atom.
_MEEKO_TYPING_FILE = "vina_params"
_MEEKO_TYPING_GROUP = "vina_typing"
_ROLE_FLAGS = {"vina_donor": "donor", "vina_acceptor": "acceptor"}


def _load_interaction_roles() -> list[tuple[str, dict]]:
    """[(smarts, {"donor": bool, "acceptor": bool})] from meeko's parameter file."""
    path = MoleculePreparation.packaged_params[_MEEKO_TYPING_FILE]
    with open(path) as handle:
        rules = json.load(handle)[_MEEKO_TYPING_GROUP]

    roles = []
    for rule in rules:
        flags = {
            role: rule[key] for key, role in _ROLE_FLAGS.items() if key in rule
        }
        if flags:
            roles.append((rule["smarts"], flags))
    if not roles:
        raise RuntimeError(
            f"no donor/acceptor rules found in meeko's {_MEEKO_TYPING_FILE}.json"
        )
    return roles


_INTERACTION_ROLES = _load_interaction_roles()


class InteractionFinder:
    """Class for handling and calculating ligand-receptor interactions.

    The receptor is read into the three objects used by the batched distance
    query: a _RECEPTOR_DTYPE atom array, a scipy KDTree over atom coordinates,
    and per-property annotation index sets, all keyed by row position in the
    array. The receptor may be supplied as a meeko Polymer (or polymer-JSON
    string) — read natively, no pdbqt — or as a legacy pdbqt string. Both
    sources produce the same triple.

    Attributes:
        hb_cutoff (float): cutoff for hydrogen-bond interactions in ångströms
        vdw_cutoff (float): cutoff for van der Waals interactions in ångströms
    """

    def __init__(self, receptor, hb_cutoff: float, vdw_cutoff: float):
        self.hb_cutoff = hb_cutoff
        self.vdw_cutoff = vdw_cutoff
        self.pdbqt_rec = None
        self._compiled_roles = [
            (Chem.MolFromSmarts(s), flags) for s, flags in _INTERACTION_ROLES
        ]
        if isinstance(receptor, Polymer):
            self._atoms_arr, self._annotations, self._kdtree = receptor_atoms_from_polymer(receptor)
        elif self._looks_like_polymer_json(receptor):
            self._atoms_arr, self._annotations, self._kdtree = receptor_atoms_from_polymer(
                Polymer.from_json(receptor)
            )
        else:
            self._atoms_arr, self._annotations, self._kdtree = receptor_atoms_from_pdbqt_string(receptor)

    @staticmethod
    def _looks_like_polymer_json(rec_string: str) -> bool:
        return isinstance(rec_string, str) and rec_string.lstrip().startswith("{")


    @staticmethod
    def _empty_result() -> dict:
        return {
            "type": [],
            "recid": [],
            "recname": [],
            "residue": [],
            "resid": [],
            "chain": [],
            "count": 0,
            "hb_count": 0,
        }

    def _atoms_for_neighbors(self, neighbor_indices, atom_property):
        """Keep the neighbour indices (from a batched query_ball_point call) that
        carry `atom_property`, then return their atom records. 
        """
        if len(neighbor_indices) == 0:
            return []
        selected = set(neighbor_indices)
        selected.intersection_update(self._annotations[atom_property])
        if not selected:
            return []
        return self._atoms_arr[list(selected)].copy()

    def ligand_interaction_atoms(self, mol: Chem.Mol) -> tuple[set[int],set[int],set[int]]:
        """
        Identify H bond donors, acceptors, and vdw atoms in ligand

        Heavy atoms only, indexed to match the pose coordinates. Topology only, so one
        call serves every pose. vdw is the complement: heavy, not H-bond capable.

        Args:
            mol (Chem.Mol)

        Returns:
            tuple[set[int],set[int],set[int]]: donor, acceptor and vdw atom indices
        """
        n = mol.GetNumAtoms()
        donor = [False]*n
        acceptor = [False]*n
        for smarts, flags in self._compiled_roles:
            for match in mol.GetSubstructMatches(smarts):
                if "donor" in flags:
                    donor[match[0]] = flags["donor"]
                if "acceptor" in flags:
                    acceptor[match[0]] = flags["acceptor"]
        donors = {i for i, v in enumerate(donor) if v}
        acceptors = {i for i, v in enumerate(acceptor) if v}
        heavy = {a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1}

        return donors, acceptors, heavy - acceptors - donors


    def find_pose_interactions(
        self, donor_idxs: set[int], acceptor_idxs: set[int], vdw_idxs: set[int], lig_coordinates: list
    ) -> dict:
        """Identify interactions for a pose within the cutoff distances.

        Args:
            donor_idxs (set[int]): ligand atom indices for hydrogen donors
            acceptor_idxs (set[int]): ligand atom indices for hydrogen acceptors
            vdw_idxs (set[int]): ligand atom indices for atoms that may participate in vdw interactions
            lig_coordinates (list): coordinates for the atoms in the ligand

        Returns:
            dict: all interaction details for a given ligand pose
        """
        type_list = []
        recid_list = []
        recname_list = []
        residue_list = []
        resid_list = []
        chain_list = []

        def append_rec_atom_info(rec_at):
            recid_list.append(str(rec_at["idx"]))
            recname_list.append(str(rec_at["name"]))
            residue_list.append(rec_at["resname"])
            resid_list.append(str(rec_at["resid"]))
            chain_list.append(rec_at["chain"])

        hb_idxs = sorted(donor_idxs | acceptor_idxs)
        h_coord_arr = np.array(
            [[float(c) for c in lig_coordinates[i]] for i in hb_idxs]
        ).reshape(-1, 3)
        vdw_coord_arr = np.array(
                    [[float(c) for c in lig_coordinates[i]] for i in sorted(vdw_idxs)]
                ).reshape(-1, 3)

        # one query per radius, over only the atoms eligible for it
        # _atoms_for_neighbors in the for loop keeps the receptor atoms of the right kind
        hb_neighbors = self._kdtree.query_ball_point(
            h_coord_arr, self.hb_cutoff, p=2, return_sorted=True
        )
        vdw_neighbors = self._kdtree.query_ball_point(
            vdw_coord_arr, self.vdw_cutoff, p=2, return_sorted=True
        )
        # enumerate possible hydrogen bonders
        for k, idx in enumerate(hb_idxs):
            if idx in acceptor_idxs:
                for rec_at in self._atoms_for_neighbors(hb_neighbors[k], "hb_don"):
                    append_rec_atom_info(rec_at)
                    type_list.append("H")
            if idx in donor_idxs:
                for rec_at in self._atoms_for_neighbors(hb_neighbors[k], "hb_acc"):
                    append_rec_atom_info(rec_at)
                    type_list.append("H")
        # enumerate across possible van der Waal bonders
        for neighbors in vdw_neighbors:
            for rec_at in self._atoms_for_neighbors(neighbors, "vdw"):
                append_rec_atom_info(rec_at)
                type_list.append("V")

        # get full interaction counts (count here will not match rows in interaction table
        # because ligand atoms are not counted, and results deduplicated on receptor atom)
        # (move the counts after unique/dedup to have counts match)
        int_count = len(type_list)
        h_count = type_list.count("H")

        unique = dict.fromkeys(
            zip(
                type_list,
                recid_list,
                recname_list,
                residue_list,
                resid_list,
                chain_list,
            )
        )
        types, recids, recnames, residues, resids, chains = (
            [list(col) for col in zip(*unique)] if unique else [[]] * 6
        )

        return {
            "type": types,
            "recid": recids,
            "recname": recnames,
            "residue": residues,
            "resid": resids,
            "chain": chains,
            "count": int_count,
            "hb_count": h_count,
        }


def make_interaction_finder(
    receptor_string: str,
    interaction_cutoffs: list,
) -> Union[InteractionFinder, None]:
    """Construct an InteractionFinder; return None on failure."""
    try:
        return InteractionFinder(receptor_string, *interaction_cutoffs)
    except Exception as e:
        logger.warning("Could not create InteractionFinder: %s", e)
        return None


def find_interactions(
    poses_coordinates: list[tuple[dict, list]],
    mol: Chem.Mol,
    receptor_string: str = None,
    hb_cutoff: float = None,
    vdw_cutoff: float = None,
    interaction_finder: InteractionFinder = None,
):
    if not interaction_finder:
        interaction_finder = make_interaction_finder(
            receptor_string, [hb_cutoff, vdw_cutoff]
        )
    if not interaction_finder:
        n = len(poses_coordinates)
        return [], [0] * n, [0] * n

    interactions = []
    num_hb = []
    num_interactions = []

    # calculate interactions for each pose
    h_don_idx,h_acc_idx,vdw_idx = interaction_finder.ligand_interaction_atoms(mol)

    for pose_meta, coords in poses_coordinates:
        pose_interactions = interaction_finder.find_pose_interactions(
             h_don_idx,h_acc_idx,vdw_idx, coords
        )
        logger.debug(
            f"Ligand {pose_meta.get('ligname', '?')} pose {pose_meta.get('pose_rank', '?')}: {pose_interactions.get('count', 0)} interactions found"
        )
        pose_interactions.update({"id": pose_meta})
        num_hb.append(pose_interactions.pop("hb_count"))
        num_interactions.append(pose_interactions["count"])
        interactions.append(pose_interactions)
    return interactions, num_hb, num_interactions
