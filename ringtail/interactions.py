#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail interaction finder
#

import numpy as np
import tempfile
from typing import Union
from meeko import PDBQTReceptor, MoleculePreparation
from .receptormanager import ReceptorManager
from .exceptions import InteractionError
from rdkit import Chem, Geometry
from .logutils import get_logger

logger = get_logger(__name__)


class InteractionFinder:
    """Class for handling and calculating ligand-receptor interactions.

    Attributes:
        rec_string (str): string describing the receptor
        hb_cutoff (float) cutoff for interactions of hydrogen bonds ain ångströms
        vdw_cutoff (float) cutoff for interactions of van der Waals interactions in ångströms
    """

    def __init__(self, rec_string: str, hb_cutoff: float, vdw_cutoff: float):
        try:
            self.pdb = PDBQTReceptor(rec_string)
        except OSError:
            with tempfile.NamedTemporaryFile(mode="wt", delete=False, suffix=".pdbqt") as f:
                f.write(rec_string)
                fname = f.name
            self.pdb = PDBQTReceptor(fname)
        except Exception:
            try:
                # assume it is a polymer json string
                pdbqt_str = ReceptorManager.polymer_json2pdbqt_str(rec_string)
                self.pdb = PDBQTReceptor(pdbqt_str)
            except Exception as e:
                raise InteractionError("No valid receptor option given") from e
        self.hb_cutoff = hb_cutoff
        self.vdw_cutoff = vdw_cutoff

        # other interaction pre-calculations

    def __call__(self, lig_atomtype_list: list, lig_coordinates: list) -> dict:
        """Identify interactions for a pose; delegates to find_pose_interactions."""
        return self.find_pose_interactions(lig_atomtype_list, lig_coordinates)

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
        needs_conversion = not isinstance(lig_coordinates, np.ndarray)

        type_list = []
        recid_list = []
        recname_list = []
        residue_list = []
        resid_list = []
        chain_list = []

        def append_rec_atom_info(rec_at):
            # rec_at: (atom_id, ?, atom_name, resid, resname, chain, xyz, q, t)
            recid_list.append(str(rec_at[0]))
            recname_list.append(str(rec_at[2]))
            residue_list.append(rec_at[4])
            resid_list.append(str(rec_at[3]))
            chain_list.append(rec_at[5])

        for idx, atomtype in enumerate(lig_atomtype_list):
            if atomtype is None:
                continue
            coords = (
                np.array([float(c) for c in lig_coordinates[idx]])
                if needs_conversion
                else lig_coordinates[idx]
            )

            for rec_at in self.pdb.closest_atoms_from_positions(
                coords, self.hb_cutoff, atom_properties="hb_don"
            ):
                if not atomtype.endswith("A"):
                    continue
                append_rec_atom_info(rec_at)
                type_list.append("H")

            for rec_at in self.pdb.closest_atoms_from_positions(
                coords, self.hb_cutoff, atom_properties="hb_acc"
            ):
                if not atomtype.endswith("D"):
                    continue
                append_rec_atom_info(rec_at)
                type_list.append("H")

            for rec_at in self.pdb.closest_atoms_from_positions(
                coords, self.vdw_cutoff, atom_properties="vdw"
            ):
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
    # calculate interactions for each pose
    mk_prep = MoleculePreparation(rigid_macrocycles=True)
    for pose_meta, coords in poses_coordinates:
        molsetup_list = mk_prep(mol)
        if not molsetup_list:
            logger.warning(
                f"MoleculePreparation returned no setups for ligand {pose_meta.get('ligname', '?')} — skipping interaction calculation for this pose"
            )
            interactions.append(
                {
                    "type": [],
                    "recid": [],
                    "recname": [],
                    "residue": [],
                    "resid": [],
                    "chain": [],
                    "count": 0,
                    "hb_count": 0,
                    "id": pose_meta,
                }
            )
            num_hb.append(0)
            num_interactions.append(0)
            continue
        molsetup = molsetup_list[0]
        atom_types = [atom.atom_type for atom in molsetup.atoms if not atom.is_ignore]
        pose_interactions = interaction_finder.find_pose_interactions(atom_types, coords)
        logger.debug(
            f"Ligand {pose_meta.get('ligname', '?')} pose {pose_meta.get('pose_rank', '?')}: {pose_interactions.get('count', 0)} interactions found"
        )
        pose_interactions.update({"id": pose_meta})
        num_hb.append(pose_interactions.pop("hb_count"))
        num_interactions.append(pose_interactions["count"])
        interactions.append(pose_interactions)
    return interactions, num_hb, num_interactions
