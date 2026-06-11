#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail interaction finder
#

import numpy as np
from typing import Union
from scipy import spatial
from meeko import PDBQTReceptor, MoleculePreparation, Polymer
from .exceptions import InteractionError
from rdkit import Chem, Geometry
from .logutils import get_logger

logger = get_logger(__name__)

# AutoDock atom type -> interaction property. Ported verbatim from meeko's
# receptor_pdbqt.atom_property_definitions so Ringtail can classify receptor
# atoms itself (from either a Polymer or a pdbqt blob) without reaching into
# meeko's private _atom_annotations. Only hb_don/hb_acc/vdw are used by the
# interaction finder; the rest are kept for parity with meeko's mapping.
ATOM_PROPERTY_DEFINITIONS = {
    "H": "vdw", "C": "vdw", "A": "vdw", "N": "vdw", "P": "vdw", "S": "vdw",
    "Br": "vdw", "I": "vdw", "F": "vdw", "Cl": "vdw",
    "NA": "hb_acc", "OA": "hb_acc", "SA": "hb_acc", "OS": "hb_acc", "NS": "hb_acc",
    "HD": "hb_don", "HS": "hb_don",
    "Mg": "metal", "Ca": "metal", "Fe": "metal", "Zn": "metal", "Mn": "metal",
    "MG": "metal", "CA": "metal", "FE": "metal", "ZN": "metal", "MN": "metal",
    "W": "water",
    "G0": "glue", "G1": "glue", "G2": "glue", "G3": "glue",
    "CG0": "glue", "CG1": "glue", "CG2": "glue", "CG3": "glue",
}

# Structured-array layout for receptor atoms built from a Polymer. Field names
# match the subset of meeko's PDBQTReceptor._atoms dtype that the interaction
# finder reads, so the pdbqt path (which reuses meeko's array via the public
# PDBQTReceptor.atoms()) and the Polymer path are interchangeable downstream.
_RECEPTOR_DTYPE = [
    ("idx", "i4"),
    ("serial", "i4"),
    ("name", "U4"),
    ("resid", "i4"),
    ("resname", "U4"),
    ("chain", "U2"),
    ("xyz", "f4", (3,)),
    ("atom_type", "U2"),
]


def _build_kdtree_and_annotations(atoms: np.ndarray):
    """Build a KDTree over receptor atom coordinates and the per-property
    annotation index sets, from any structured array carrying 'xyz' and
    'atom_type' fields. Replaces meeko's private _KDTree / _atom_annotations
    with Ringtail-owned equivalents (identical construction).

    Returns:
        tuple[scipy.spatial.cKDTree, dict[str, set[int]]]
    """
    kdtree = spatial.cKDTree(atoms["xyz"])
    annotations = {"hb_acc": set(), "hb_don": set(), "vdw": set(), "metal": set()}
    for i, atom_type in enumerate(atoms["atom_type"]):
        prop = ATOM_PROPERTY_DEFINITIONS.get(str(atom_type))
        if prop in annotations:
            annotations[prop].add(i)
    return kdtree, annotations


def _receptor_atoms_from_polymer(polymer: Polymer) -> np.ndarray:
    """Read rigid receptor atoms directly from a meeko Polymer, returning a
    structured array of _RECEPTOR_DTYPE. Mirrors PDBQTWriterLegacy.write_from_polymer's
    atom selection (skip movable monomers, ignored/padding atoms, and flexible
    sidechain atoms) so the resulting atom set matches the legacy pdbqt path.
    """
    mk_prep = MoleculePreparation(load_atom_params=["ad4_types"])
    polymer.parameterize(mk_prep)

    records = []
    serial = 0
    for res_id, monomer in polymer.get_valid_monomers().items():
        # Mirror write_from_polymer exactly: every valid monomer contributes its
        # non-ignored, non-flexres atoms to the rigid receptor — including the
        # backbone of flexible residues. Only per-atom is_flexres_atom is skipped,
        # NOT the whole movable monomer.
        chain, resnum = res_id.split(":")
        if resnum and resnum[-1].isalpha():
            resnum = resnum[:-1]
        resid = int(resnum)
        molsetup = monomer.molsetup
        resname = monomer.input_resname
        is_flexres_atom = monomer.is_flexres_atom or []
        for atom_idx, atom in enumerate(molsetup.atoms):
            if atom.is_ignore:
                continue
            if atom_idx < len(is_flexres_atom) and is_flexres_atom[atom_idx]:
                continue
            xyz = atom.coord
            records.append(
                (
                    serial,
                    serial,
                    atom.pdbinfo.name,
                    resid,
                    resname,
                    chain,
                    (xyz[0], xyz[1], xyz[2]),
                    atom.atom_type,
                )
            )
            serial += 1
    if not records:
        raise InteractionError("No rigid receptor atoms found in polymer")
    return np.array(records, dtype=_RECEPTOR_DTYPE)


def _receptor_atoms_from_pdbqt(pdbqt_string: str) -> np.ndarray:
    """Parse a legacy pdbqt receptor via meeko's PUBLIC PDBQTReceptor.atoms().
    Returns meeko's structured array (which carries 'xyz'/'atom_type' plus the
    name/resid/resname/chain fields the finder reads). No private attributes."""
    return PDBQTReceptor(pdbqt_string).atoms()


def _looks_like_polymer_json(rec_string: str) -> bool:
    return isinstance(rec_string, str) and rec_string.lstrip().startswith("{")


class InteractionFinder:
    """Class for handling and calculating ligand-receptor interactions.

    The receptor is read into a Ringtail-owned set of three objects used by the
    batched distance query: an atom array, a scipy KDTree over atom coordinates,
    and per-property annotation index sets. The receptor may be supplied as a
    meeko Polymer (or polymer-JSON string) — read natively, no pdbqt — or as a
    legacy pdbqt string, parsed via meeko's public PDBQTReceptor.atoms().

    Attributes:
        hb_cutoff (float): cutoff for hydrogen-bond interactions in ångströms
        vdw_cutoff (float): cutoff for van der Waals interactions in ångströms
    """

    def __init__(self, receptor, hb_cutoff: float, vdw_cutoff: float):
        self.hb_cutoff = hb_cutoff
        self.vdw_cutoff = vdw_cutoff
        try:
            if isinstance(receptor, Polymer):
                self._atoms_arr = _receptor_atoms_from_polymer(receptor)
            elif _looks_like_polymer_json(receptor):
                self._atoms_arr = _receptor_atoms_from_polymer(
                    Polymer.from_json(receptor)
                )
            else:
                self._atoms_arr = _receptor_atoms_from_pdbqt(receptor)
        except InteractionError:
            raise
        except Exception as e:
            raise InteractionError("No valid receptor option given") from e

        self._kdtree, self._annotations = _build_kdtree_and_annotations(
            self._atoms_arr
        )

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
        carry `atom_property`, then return their atom records. This is the
        post-processing meeko's closest_atoms_from_positions does per point:
        set(neighbours) ∩ annotation -> atom records.
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

        return {
            "type": type_list,
            "recid": recid_list,
            "recname": recname_list,
            "residue": residue_list,
            "resid": resid_list,
            "chain": chain_list,
            "count": len(type_list),
            "hb_count": type_list.count("H"),
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
        interaction_finder = make_interaction_finder(receptor_string, [hb_cutoff, vdw_cutoff])
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
        ligname = poses_coordinates[0][0].get("ligname", "?") if poses_coordinates else "?"
        logger.warning(
            f"MoleculePreparation returned no setups for ligand {ligname} — skipping interaction calculation for all of its poses"
        )
        def empty(pose_meta):
            return {**InteractionFinder._empty_result(), "id": pose_meta}

        n = len(poses_coordinates)
        return [empty(pm) for pm, _ in poses_coordinates], [0] * n, [0] * n

    molsetup = molsetup_list[0]
    atom_types = [atom.atom_type for atom in molsetup.atoms if not atom.is_ignore]

    # calculate interactions for each pose
    for pose_meta, coords in poses_coordinates:
        pose_interactions = interaction_finder.find_pose_interactions(atom_types, coords)
        logger.debug(
            f"Ligand {pose_meta.get('ligname', '?')} pose {pose_meta.get('pose_rank', '?')}: {pose_interactions.get('count', 0)} interactions found"
        )
        pose_interactions.update({"id": pose_meta})
        num_hb.append(pose_interactions.pop("hb_count"))
        num_interactions.append(pose_interactions["count"])
        interactions.append(pose_interactions)
    return interactions, num_hb, num_interactions
