#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail receptor manager
#

import gzip


class ReceptorManager:
    @staticmethod
    def make_receptor_blobs(file_list):
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

        return receptors

    @staticmethod
    def blob2str(receptor_blob):
        return gzip.decompress(receptor_blob).decode()
