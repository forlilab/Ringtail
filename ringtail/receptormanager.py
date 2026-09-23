#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail receptor manager
#

from pathlib import Path
import gzip
from dataclasses import dataclass
from typing import Optional
from .logutils import get_logger
from .exceptions import ReceptorError

logger = get_logger(__name__)
from meeko import PDBQTWriterLegacy, Polymer, MoleculePreparation

import numpy as np
from scipy import spatial



@dataclass(frozen=True)
class ReceptorData:
    """Receptor data retrieved from the database.

    A receptor may be stored as a pdbqt blob (blob_str), a meeko polymer JSON
    string (polymer_json), or both. Use the accessor methods to get the right
    representation for each use case rather than reading fields directly.
    """

    name: Optional[str]
    blob_str: Optional[str]
    polymer_json: Optional[str]

    def receptor_string(self) -> Optional[str]:
        """String for interaction calculation: meeko Polymer JSON preferred, else
        pdbqt blob str. The interaction finder reads receptor atoms natively from a
        Polymer (no pdbqt round-trip) and only parses the pdbqt blob for legacy
        receptors that have no Polymer JSON."""
        return self.polymer_json or self.blob_str

    def pdbqt_str(self) -> Optional[str]:
        """Pdbqt-format string, converting polymer JSON if needed. Use for .pdbqt file export."""
        if self.blob_str:
            return self.blob_str
        if self.polymer_json:
            return polymer_json2pdbqt_str(self.polymer_json)
        return None




def make_receptor_blob(receptor_file: str) -> tuple[str, bytes]:
    """Creates compressed receptor info (blob)

    Args:
        receptor_file (str): path to receptor file

    Returns:
        tuple[str, bytes]: rec_name and blob (compressed receptor)
    """
    # lstrip(".") so a dot-leading filename doesn't yield an empty name, and
    # split rather than Path().stem so "4j8m.pdbqt.gz" still gives "4j8m"
    rec_name = Path(receptor_file).name.lstrip(".").split(".")[0]
    if receptor_file.endswith(".gz"):
        with open(receptor_file, "rb") as r:
            receptor = r.read()
    else:
        with open(receptor_file, "r") as r:
            receptor = gzip.compress(r.read().encode())
    logger.debug(f"Receptor blob for receptor {rec_name} prepared successfully.")
    return rec_name, receptor


def blob2str(receptor_blob):
    """Decompresses a receptor blob to a string.

    Args:
        receptor_blob (bytes): gzip-compressed receptor blob

    Returns:
        str: receptor string, or None if receptor_blob is None
    """
    if receptor_blob is None:
        return None
    return gzip.decompress(receptor_blob).decode()


def receptor_str_from_file(receptor_file: str) -> str:
    if receptor_file.endswith(".gz"):
        with gzip.open(receptor_file, "rt") as r:
            return r.read()
    else:
        with open(receptor_file, "r") as r:
            return r.read()


def _parse_polymer_json(polymer_json: str) -> tuple[str, dict]:
    """
    Makes a polymer object from a receptor polymer json, and uses
    meeko method PDBQTWriterLegacy to create a string and dict representation
    of the receptor in the pdbqt format

    Args:
        polymer_json (str): json string (not dict) representation of receptor

    Returns:
        tuple[str, dict]: _description_
    """
    polymer = Polymer.from_json(polymer_json)
    mk_prep = MoleculePreparation(load_atom_params=["ad4_types"])
    polymer.parameterize(mk_prep)
    return PDBQTWriterLegacy.write_from_polymer(polymer)


def polymer_json2pdbqt_str(polymer_json: str) -> str:
    """
    Returns pdbqt string representation of a polymer json

    Args:
        polymer_json (str): _description_

    Returns:
        str: _description_
    """
    return _parse_polymer_json(polymer_json)[0]


def polymer_json2pdbqt_dict(polymer_json: str) -> dict:
    """
    Returns pdbqt dict representation of a polymer json

    Args:
        polymer_json (str): _description_

    Returns:
        dict: dict repr of the pdbqt string
    """
    return _parse_polymer_json(polymer_json)[1]




_ATOM_PROPERTY_DEFINITIONS = {
    "H": "vdw",
    "C": "vdw",
    "A": "vdw",
    "N": "vdw",
    "P": "vdw",
    "S": "vdw",
    "Br": "vdw",
    "I": "vdw",
    "F": "vdw",
    "Cl": "vdw",
    "NA": "hb_acc",
    "OA": "hb_acc",
    "SA": "hb_acc",
    "OS": "hb_acc",
    "NS": "hb_acc",
    "HD": "hb_don",
    "HS": "hb_don",
    "Mg": "metal",
    "Ca": "metal",
    "Fe": "metal",
    "Zn": "metal",
    "Mn": "metal",
    "MG": "metal",
    "CA": "metal",
    "FE": "metal",
    "ZN": "metal",
    "MN": "metal",
}

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

def _annotate_receptor_atoms(atoms: np.ndarray) -> dict[str, set[int]]:
    """Bucket receptor atoms by AD4 type, keyed by row position in `atoms`.

    One bucket per atom type, straight out of _ATOM_PROPERTY_DEFINITIONS, so the sets
    are disjoint -> an HD is a donor and nothing else, an OA an acceptor and nothing else.
    From meeko's classification 

    Unknown atom types land in vdw 
    """
    annotations = {"hb_acc": set(), "hb_don": set(), "vdw": set(), "metal": set()}
    for row, atom_type in enumerate(atoms["atom_type"]):
        annotations[_ATOM_PROPERTY_DEFINITIONS.get(str(atom_type), "vdw")].add(row)
    return annotations


def receptor_atoms_from_pdbqt_string(
    pdbqt_string: str,
) -> tuple[np.ndarray, dict[str, set[int]], spatial.cKDTree]:
    """Parse a legacy pdbqt receptor into the same three objects the Polymer
    path produces: a _RECEPTOR_DTYPE atom array, per-property annotation index
    lists, and a KDTree over the coordinates.

    Only the fields the interaction finder reads are kept. The pdbqt columns for
    partial charge, altloc, insertion code, occupancy, b-factor and record type
    are parsed by nothing downstream, so they are not read.

    Returns:
        tuple[np.ndarray, dict[str, set[int]], spatial.cKDTree]: atoms,
        annotations (keys 'hb_acc', 'hb_don', 'vdw', 'metal'), and a KDTree.
        All indices are row positions in atoms, not its 'idx' field.
    """

    atoms = []
    atom_annotations = {"hb_acc": set(), "hb_don": set(), "vdw": set(), "metal": set()}
    # TZ is a pseudo atom for AutoDock4Zn FF
    pseudo_atom_types = ["TZ"]

    idx = 0
    for line in pdbqt_string.split("\n"):
        if line.startswith("ATOM") or line.startswith("HETATM"):
            serial = int(line[6:11].strip())
            name = line[12:16].strip()
            resname = line[17:20].strip()
            chainid = line[21].strip()
            resid = int(line[22:26].strip())
            xyz = np.array(
                [line[30:38].strip(), line[38:46].strip(), line[46:54].strip()],
                dtype=np.float32,
            )
            atom_type = line[77:79].strip()

            if not atom_type in pseudo_atom_types:
                # key annotations by ROW POSITION in `atoms`, not by `idx`
                atoms.append(
                    (
                        idx,
                        serial,
                        name,
                        resid,
                        resname,
                        chainid,
                        xyz,
                        atom_type,
                    )
                )

            idx += 1
    if not atoms:
        raise ReceptorError("no receptor atoms found in pdbqt string")
    atoms = np.array(atoms, dtype=_RECEPTOR_DTYPE)

    KDTree = spatial.cKDTree(atoms["xyz"])

    return atoms, _annotate_receptor_atoms(atoms), KDTree


def receptor_atoms_from_polymer(
    polymer: Polymer,
) -> tuple[np.ndarray, dict[str, set[int]], spatial.cKDTree]:
    """Read rigid receptor atoms directly from a meeko Polymer. Mirrors
    PDBQTWriterLegacy.write_from_polymer's atom selection (skip movable monomers,
    ignored/padding atoms, and flexible sidechain atoms) so the resulting atom
    set matches the legacy pdbqt path.

    Returns:
        tuple[np.ndarray, dict[str, set[int]], spatial.cKDTree]: atoms,
        annotations (keys 'hb_acc', 'hb_don', 'vdw', 'metal'), and a KDTree.
        All indices are row positions in atoms.
    """
    mk_prep = MoleculePreparation(load_atom_params=["ad4_types"])
    polymer.parameterize(mk_prep)

    records = []
    serial = 0
    for res_id, monomer in polymer.get_valid_monomers().items():
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
        raise ReceptorError("No rigid receptor atoms found in polymer")
    atoms = np.array(records, dtype=_RECEPTOR_DTYPE)
    kdtree = spatial.cKDTree(atoms["xyz"])

    return atoms, _annotate_receptor_atoms(atoms), kdtree
