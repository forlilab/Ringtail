#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail receptor manager
#

import gzip
from .logutils import LOGGER


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
