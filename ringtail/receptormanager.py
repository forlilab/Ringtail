#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail receptor manager
#

import gzip
from .logutils import get_logger

logger = get_logger(__name__)
from meeko import PDBQTWriterLegacy, Polymer, MoleculePreparation


class ReceptorManager:
    """Class with methods dealing with formatting of receptor information"""

    @staticmethod
    def make_receptor_blob(receptor_file: str) -> tuple[str, bytes]:
        """Creates compressed receptor info (blob)

        Args:
            receptor_file (str): path to receptor file

        Returns:
            tuple[str, bytes]: rec_name and blob (compressed receptor)
        """
        rec_name = receptor_file.split(".")[0].split("/")[-1]
        if receptor_file.endswith(".gz"):
            with open(receptor_file, "rb") as r:
                receptor = r.read()
        else:
            with open(receptor_file, "r") as r:
                receptor = gzip.compress(r.read().encode())
        logger.debug(f"Receptor blob for receptor {rec_name} prepared successfully.")
        return rec_name, receptor

    @staticmethod
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

    @staticmethod
    def receptor_str_from_file(receptor_file: str) -> str:
        if receptor_file.endswith(".gz"):
            with gzip.open(receptor_file, "rt") as r:
                return r.read()
        else:
            with open(receptor_file, "r") as r:
                return r.read()

    @staticmethod
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

    @staticmethod
    def polymer_json2pdbqt_str(polymer_json: str) -> str:
        """
        Returns pdbqt string representation of a polymer json

        Args:
            polymer_json (str): _description_

        Returns:
            str: _description_
        """
        return ReceptorManager._parse_polymer_json(polymer_json)[0]

    @staticmethod
    def polymer_json2pdbqt_dict(polymer_json: str) -> dict:
        """
        Returns pdbqt dict representation of a polymer json

        Args:
            polymer_json (str): _description_

        Returns:
            dict: dict repr of the pdbqt string
        """
        return ReceptorManager._parse_polymer_json(polymer_json)[1]
