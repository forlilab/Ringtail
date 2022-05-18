#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail receptor dataclass
#

import dataclasses
from dataclasses import dataclass, field
import os
import json


@dataclass(frozen=True)
class Receptor:
    name: str
    sequence: str = field(default_factory=str, hash=False)
    residues: dict = field(default_factory=dict, compare=False, hash=False, repr=False)

    def add_AA_to_sequence(self, new_AA: str) -> "Receptor":
        return dataclasses.replace(self, sequence=self.sequence + new_AA)


@dataclass(frozen=True)
class Residue:
    name: str
    num: int
    chain: str
    atoms: dict = field(default_factory=dict, compare=False, hash=False, repr=False)


@dataclass(frozen=True)
class Atom:
    name: str
    num: int
    atomtype: str
    atm_flag: str
    x: float
    y: float
    z: float
    occupancy: float
    b_iso: float
    q: float


def make_receptor_object(fname, name):
    """Make receptor object from list of pdbqt line dictionaries

    Args:
        pdbqt_lines (list): list of dicitonaries, with each dictionary containing the data for one atom
    """

    pdbqt_lines = receptor_pdbqt_parser(fname)
    receptor = Receptor(name)

    with open(os.path.join(os.path.dirname(__file__),
                           'AA_CODES_3TO1.json'), 'r') as f:
        aminoacid_codes_dict = json.load(f)

    for atom in pdbqt_lines:
        atom_name = atom["atm_name"]
        atom_num = atom["atm_num"]
        atomtype = atom["atomtype"]
        atm_flag = atom["atm_flag"]
        x = atom["x"]
        y = atom["y"]
        z = atom["z"]
        occupancy = atom["occupancy"]
        b_iso = atom["b_iso"]
        q = atom["q"]
        new_atom = Atom(atom_name, atom_num, atomtype, atm_flag, x, y, z, occupancy, b_iso, q)

        resname = atom["res_name"]
        resnum = atom["res_num"]
        chain = atom["chain"]
        current_res = Residue(resname, resnum, chain)
        if current_res in receptor.residues:
            receptor.residues[current_res].atoms[new_atom] = new_atom
        else:
            current_res.atoms[new_atom] = new_atom
            receptor.residues[current_res] = current_res
            receptor = receptor.add_AA_to_sequence(aminoacid_codes_dict[current_res.name.replace(" ", "")])

    return receptor

def write_receptor_pdbqt(receptor_object, fname=""):
    """Writes pdbqt file from receptor object 
    
    Args:
        receptor_object (Receptor): Reptor object to be written to pdbqt
        fname (string): name of output pdbqt file. By default will be <receptor_object.name>.pdbqt
    """
    raise NotImplementedError
