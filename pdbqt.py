import dataclasses
from dataclasses import dataclass, field
from parsers import pdbqt_parser
import os
import json


@dataclass(frozen=True)
class PDBQT:
    name: str
    sequence: str = field(default_factory=str, hash=False)
    residues: dict = field(default_factory=dict, compare=False, hash=False, repr=False)

    def add_AA_to_sequence(self, new_AA: str) -> "PDBQT":
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


class PDBQTReader():

    def __init__(self, fname):
        self.fname = fname
        self.name = fname.split(".")[0]

    def make_pdbqt_object(self):
        """Make PDBQT object from list of pdbqt line dictionaries

        Args:
            pdbqt_lines (list): list of dicitonaries, with each dictionary containing the data for one atom
        """

        pdbqt_lines = pdbqt_parser(self.fname)
        pdbqt = PDBQT(self.name)

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
            if current_res in pdbqt.residues:
                pdbqt.residues[current_res].atoms[new_atom] = new_atom
            else:
                current_res.atoms[new_atom] = new_atom
                pdbqt.residues[current_res] = current_res
                pdbqt = pdbqt.add_AA_to_sequence(aminoacid_codes_dict[current_res.name.replace(" ", "")])

        return pdbqt
