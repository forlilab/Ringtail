#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail interaction finder
#

import numpy as np
import tempfile
from meeko import PDBQTReceptor, MoleculePreparation
from .receptormanager import ReceptorManager
from rdkit import Chem, Geometry


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
        except OSError as e:
            with tempfile.NamedTemporaryFile(mode="wt") as f:
                f.write(rec_string)
                self.pdb = PDBQTReceptor(f.name)
        except:
            try:
                # assume it is a polymer json string
                pdbqt_str = ReceptorManager.polymer_json2pdbqt_str(rec_string)
                self.pdb = PDBQTReceptor(pdbqt_str)
            except:
                raise Exception("No valid receptor option given")
        self.hb_cutoff = hb_cutoff
        self.vdw_cutoff = vdw_cutoff

    def find_pose_interactions(
        self, lig_atomtype_list: list, lig_coordinates: list
    ) -> dict:
        """Method that identifies interactions for a pose within th given cutoff distances in the main class.

        Args:
            lig_atomtype_list (list): list of atoms in the ligand
            lig_coordinates (list): coordinates for the atoms in the ligand

        Returns:
            dict: all interaction details for a given ligand pose
        """

        def append_rec_atom_info(rec_at):
            # rec_at is array of format (atom_id, ?, atom_name, resid, resname, chain, xyz, q, t)
            recid_list.append(str(rec_at[0]))
            recname_list.append(str(rec_at[2]))
            residue_list.append(rec_at[4])
            resid_list.append(str(rec_at[3]))
            chain_list.append(rec_at[5])

        type_list = []
        recid_list = []
        recname_list = []
        residue_list = []
        resid_list = []
        chain_list = []
        # check coordinate type:
        if not type(lig_coordinates) == np.ndarray:
            needs_conversion = True
        else:
            needs_conversion = False

        for idx, atomtype in enumerate(lig_atomtype_list):
            if needs_conversion:
                coords = np.array([float(coord) for coord in lig_coordinates[idx]])
            else:
                coords = lig_coordinates[idx]

            hbd_neighbors = self.pdb.closest_atoms_from_positions(
                coords, self.hb_cutoff, atom_properties="hb_don"
            )
            for rec_at in hbd_neighbors:
                # rec_at is array of format (atom_id, atom_name, resname, resid, chainid, xyz, q, t)
                if not atomtype.endswith("A"):
                    continue
                append_rec_atom_info(rec_at)
                type_list.append("H")

            hba_neighbors = self.pdb.closest_atoms_from_positions(
                coords, self.hb_cutoff, atom_properties="hb_acc"
            )
            for rec_at in hba_neighbors:
                # rec_at is array of format (atom_id, atom_name, resname, resid, chainid, xyz, q, t)
                if not atomtype.endswith("D"):
                    continue
                append_rec_atom_info(rec_at)
                type_list.append("H")

            vdw_neighbors = self.pdb.closest_atoms_from_positions(
                coords, self.vdw_cutoff, atom_properties="vdw"
            )
            for rec_at in vdw_neighbors:
                # rec_at is array of format (atom_id, atom_name, resname, resid, chainid, xyz, q, t)
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


def calculate_interactions(
    poses_coordinates: list[tuple[dict, list]],
    mol: Chem.Mol,
    receptor_string: str,
    hb_cutoff: float,
    vdw_cutoff: float,
):
    interaction_finder = InteractionFinder(receptor_string, hb_cutoff, vdw_cutoff)
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
    # calculate interactions for each pose
    for id, coords in poses_coordinates:
        # make a molsetup for the Mol which includes atom types needed for interaction calculations
        mk_prep = MoleculePreparation(rigid_macrocycles=True)
        molsetup_list = mk_prep(mol)
        molsetup = molsetup_list[0]
        atom_types = []
        for _, atom in enumerate(molsetup.atoms):
            if atom.is_ignore:
                continue
            atom_types.append(atom.atom_type)
        # engage interaction finder
        pose_interactions = interaction_finder.find_pose_interactions(
            atom_types, coords
        )
        # add unique
        pose_interactions.update({"id": id})
        num_hb.append(pose_interactions.pop("hb_count"))
        num_interactions.append(pose_interactions["count"])
        interactions.append(pose_interactions)
    return interactions, num_hb, num_interactions
