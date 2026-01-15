#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail receptor manager
#

import gzip
from .logutils import LOGGER
from meeko import PDBQTWriterLegacy, Polymer, MoleculePreparation


class ReceptorManager:
    """Class with methods dealing with formatting of receptor information"""

    @staticmethod
    def make_receptor_blob(receptor_file: str) -> tuple[str, bytes]:
        """Creates compressed receptor info (blob)

        Args:
            file_list (str): path to receptor file

        Returns:
            dict: rec_name and blob (compressed receptor)
        """
        # check file extension, compress to bytes if needed
        rec_name = receptor_file.split(".")[0].split("/")[
            -1
        ]  # remove file extension and path
        if receptor_file.endswith(".gz"):
            with open(receptor_file, "rb") as r:
                receptor = r.read()
        else:
            with open(receptor_file, "r") as r:
                receptor = gzip.compress(r.read().encode())
        LOGGER.debug(f"Receptor blob for receptor {rec_name} parepared successfully.")
        return rec_name, receptor

    @staticmethod
    def blob2str(receptor_blob):
        """Creates blob of compresser receptor file info

        Args:
            receptor_blob (blob): zipped receptor blob

        Returns:
            str: receptor string
        """
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
    def _parse_polymer_json(polymer_json) -> tuple[str, dict]:
        """
        Makes a polymer object from a receptor polymer json, and uses
        meeko method PDBQTWriterLegacy to create a string and dict representation
        of the receptor in the pdbqt format

        Args:
            polymer_json (_type_): _description_

        Returns:
            tuple[str, dict]: _description_
        """
        polymer = Polymer.from_json(polymer_json)
        mk_prep = MoleculePreparation(load_atom_params=["ad4_types"])
        polymer.parameterize(mk_prep)
        return PDBQTWriterLegacy.write_from_polymer(polymer)

    @staticmethod
    def polymer_json2pdbqt_str(polymer_json) -> str:
        """
        Returns pdbqt string representation of a polymer json

        Args:
            polymer_json (_type_): _description_

        Returns:
            str: _description_
        """
        return ReceptorManager._parse_polymer_json(polymer_json)[0]

    @staticmethod
    def polymer_json2pdbqt_dict(polymer_json) -> dict:
        """
        Returns pdbqt dict representation of a polymer json

        Args:
            polymer_json (_type_): _description_

        Returns:
            dict: _description_
        """
        return ReceptorManager._parse_polymer_json(polymer_json)[1]
