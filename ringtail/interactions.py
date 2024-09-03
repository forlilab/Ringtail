#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail interaction finder
#

import numpy as np
import tempfile
import meeko
from meeko import PDBQTReceptor


class InteractionFinder:
    """Class for handling and calculating ligand-receptor interactions.

    Attributes:
        rec_string (str): string describing the receptor
        interaction_cutoff_radii (list(float)): cutoff for interactions of hydrogen bonds and VDW interactions, in ångströms
    """

    def __init__(self, rec_string, interaction_cutoff_radii):
        self.rec_string = rec_string
        try:
            self.pdb = PDBQTReceptor(rec_string)
        except OSError as e:
            with tempfile.NamedTemporaryFile(dir="/dev/shm", mode="wt") as f:
                f.write(rec_string)
                self.pdb = PDBQTReceptor(f.name)
        self.interaction_cutoff_radii = interaction_cutoff_radii

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

        for idx, atomtype in enumerate(lig_atomtype_list):
            coords = np.array([float(coord) for coord in lig_coordinates[idx]])

            hbd_neighbors = self.pdb.closest_atoms_from_positions(
                coords, self.interaction_cutoff_radii[0], atom_properties="hb_don"
            )
            for rec_at in hbd_neighbors:
                # rec_at is array of format (atom_id, atom_name, resname, resid, chainid, xyz, q, t)
                if not atomtype.endswith("A"):
                    continue
                append_rec_atom_info(rec_at)
                type_list.append("H")

            hba_neighbors = self.pdb.closest_atoms_from_positions(
                coords, self.interaction_cutoff_radii[0], atom_properties="hb_acc"
            )
            for rec_at in hba_neighbors:
                # rec_at is array of format (atom_id, atom_name, resname, resid, chainid, xyz, q, t)
                if not atomtype.endswith("D"):
                    continue
                append_rec_atom_info(rec_at)
                type_list.append("H")

            vdw_neighbors = self.pdb.closest_atoms_from_positions(
                coords, self.interaction_cutoff_radii[1], atom_properties="vdw"
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
            "count": [str(len(type_list))],
            "ligid": [],
            "ligname": [],
        }
