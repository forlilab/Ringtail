#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail interaction finder
#

import numpy as np
from typing import Union
from meeko import MoleculePreparation, Polymer
from .receptormanager import (
    receptor_atoms_from_pdbqt_string,
    receptor_atoms_from_polymer
)
from rdkit import Chem, Geometry
from .logutils import get_logger

logger = get_logger(__name__)



def _looks_like_polymer_json(rec_string: str) -> bool:
    return isinstance(rec_string, str) and rec_string.lstrip().startswith("{")


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
        if isinstance(receptor, Polymer):
            self._atoms_arr, self._annotations, self._kdtree = receptor_atoms_from_polymer(receptor)
        elif _looks_like_polymer_json(receptor):
            self._atoms_arr, self._annotations, self._kdtree = receptor_atoms_from_polymer(
                Polymer.from_json(receptor)
            )
        else:
            self._atoms_arr, self._annotations, self._kdtree = receptor_atoms_from_pdbqt_string(receptor)

    def __call__(self, lig_atomtype_list: list, lig_coordinates: list) -> dict:
        """Identify interactions for a pose; delegates to find_pose_interactions."""
        return self.find_pose_interactions(lig_atomtype_list, lig_coordinates)

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

    def find_pose_interactions(
        self, lig_atomtype_list: list, lig_coordinates: list
    ) -> dict:
        """Identify interactions for a pose within the cutoff distances.

        Args:
            lig_atomtype_list (list): list of atom types in the ligand
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

        # The original code ran the hb_don/hb_acc queries for every atom but
        # discarded every result unless the ligand atom was an acceptor ("A") /
        # donor ("D"). Checking the type before querying is equivalent and skips
        # the KDTree query entirely for non-acceptor/non-donor atoms.
        valid = [i for i, at in enumerate(lig_atomtype_list) if at is not None]
        if not valid:
            return self._empty_result()

        # One KDTree query per radius for the whole pose, then per-point
        # post-processing (set ∩ annotation -> atom records) in
        # _atoms_for_neighbors. Equivalent to per-atom queries; only vectorized.
        coords_arr = np.asarray(
            [[float(c) for c in lig_coordinates[i]] for i in valid], dtype=float
        )
        hb_neighbors = self._kdtree.query_ball_point(
            coords_arr, self.hb_cutoff, p=2, return_sorted=True
        )
        vdw_neighbors = self._kdtree.query_ball_point(
            coords_arr, self.vdw_cutoff, p=2, return_sorted=True
        )
        for k, idx in enumerate(valid):
            atomtype = lig_atomtype_list[idx]
            if atomtype.endswith("A"):
                for rec_at in self._atoms_for_neighbors(hb_neighbors[k], "hb_don"):
                    append_rec_atom_info(rec_at)
                    type_list.append("H")
            if atomtype.endswith("D"):
                for rec_at in self._atoms_for_neighbors(hb_neighbors[k], "hb_acc"):
                    append_rec_atom_info(rec_at)
                    type_list.append("H")
            for rec_at in self._atoms_for_neighbors(vdw_neighbors[k], "vdw"):
                append_rec_atom_info(rec_at)
                type_list.append("V")

        # Deduplicate on the receptor atom. Several ligand atoms can fall inside one
        # receptor atom's cutoff sphere, which produced one entry each — but the
        # ligand atom is not part of what gets stored, so those collapsed to a single
        # Interactions row while count/hb_count still reported the pairs. Deduplicating
        # here, where the counts are derived, is what keeps Results.num_hb and
        # Results.num_interactions equal to the rows the database actually holds.
        # (Remove this block to go back to counting ligand-atom/receptor-atom pairs.)
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
            "count": len(types),
            "hb_count": types.count("H"),
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

    num_atoms = mol.GetNumAtoms()
    conf = Chem.Conformer(num_atoms)
    # add some conformer so meeko is happy
    for i in range(num_atoms):
        conf.SetAtomPosition(i, Geometry.Point3D(0, 0, 0))
    # Add conformer to molecule
    mol.AddConformer(conf, assignId=True)
    mol = Chem.AddHs(mol, addCoords=True)

    # The molsetup (and thus atom_types) depends only on mol's topology, which is
    # identical for every pose here — pose coordinates are applied separately in
    # find_pose_interactions, never to `mol`. So prepare once instead of per pose
    # (MoleculePreparation was ~3/4 of this function's runtime when looped).
    mk_prep = MoleculePreparation(rigid_macrocycles=True)
    molsetup_list = mk_prep(mol)
    if not molsetup_list:
        ligname = (
            poses_coordinates[0][0].get("ligname", "?") if poses_coordinates else "?"
        )
        logger.warning(
            f"MoleculePreparation returned no setups for ligand {ligname} — skipping interaction calculation for all of its poses"
        )

        def empty(pose_meta):
            return {**InteractionFinder._empty_result(), "id": pose_meta}

        n = len(poses_coordinates)
        return [empty(pm) for pm, _ in poses_coordinates], [0] * n, [0] * n

    molsetup = molsetup_list[0]
    # Full length, with None where meeko ignores an atom (merged non-polar hydrogens).
    # find_pose_interactions indexes pose coordinates by position in this list, so a
    # compacted list silently shifts every atom after the first ignored one onto another
    # atom's coordinates. Heavy atoms survive that (they lead the stored molecule) but
    # polar hydrogens do not, and they are the ligand's only H-bond donors.
    atom_types = [None] * num_atoms
    for i, atom in enumerate(molsetup.atoms):
        if not atom.is_ignore and i < num_atoms:
            atom_types[i] = atom.atom_type
    if mol.GetNumAtoms() != num_atoms:
        ligname = (
            poses_coordinates[0][0].get("ligname", "?") if poses_coordinates else "?"
        )
        logger.warning(
            f"Ligand {ligname}: AddHs added {mol.GetNumAtoms() - num_atoms} hydrogen(s), "
            "which have no stored coordinates and are excluded from interactions. "
            "Hydrogen bonds donated by the ligand may be incomplete."
        )

    # calculate interactions for each pose
    for pose_meta, coords in poses_coordinates:
        pose_interactions = interaction_finder.find_pose_interactions(
            atom_types, coords
        )
        logger.debug(
            f"Ligand {pose_meta.get('ligname', '?')} pose {pose_meta.get('pose_rank', '?')}: {pose_interactions.get('count', 0)} interactions found"
        )
        pose_interactions.update({"id": pose_meta})
        num_hb.append(pose_interactions.pop("hb_count"))
        num_interactions.append(pose_interactions["count"])
        interactions.append(pose_interactions)
    return interactions, num_hb, num_interactions
