#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail receptor manager
#

import gzip
from .logutils import LOGGER as logger


class ReceptorManager:
    # TODO add b2z method too?
    """Class with methods dealing with formatting of receptor information"""

    @staticmethod
    def make_receptor_blobs(file_list):
        """Creates compressed receptor info

        Args:
            file_list (str): path to receptor file

        Returns:
            blob: compressed receptor
        """
        receptors = []
        for rec_file in file_list:
            # check file extension, compress to bytes if needed
            rec_name = rec_file.split(".")[0].split("/")[
                -1
            ]  # remove file extension and path
            if rec_file.endswith(".gz"):
                with open(rec_file, "rb") as r:
                    receptors.append((r.read(), rec_name))
            else:
                with open(rec_file, "r") as r:
                    receptors.append((gzip.compress(r.read().encode()), rec_name))
        logger.debug("Receptor blob parepared successfully.")
        return receptors

    @staticmethod
    def blob2str(receptor_blob):
        """Creates blob of compresser receptor file info

        Args:
            receptor_blob (blob): zipped receptor blob

        Returns:
            str: receptor string
        """
        return gzip.decompress(receptor_blob).decode()
