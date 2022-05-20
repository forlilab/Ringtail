#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail receptor manager
#

import gzip

class ReceptorManager():

    def __init__(self, rec_file_list, dbmanager):

        self.file_list = rec_file_list
        self.receptors = []
        self.dbman = dbmanager

        self._make_receptor_blobs()

    def _make_receptor_blobs(self):
        for rec_file in self.file_list:
            # check file extension, compress to bytes if needed
            rec_name = rec_file.split(".")[0].split("/")[-1]  # remove file extension and path
            if rec_file.endswith(".gz"):
                with open(rec_file, 'rb') as r:
                    self.receptors.append((r.read(), rec_name))
            else:
                with open(rec_file, 'r') as r:
                    self.receptors.append((gzip.compress(r.read().encode()), rec_name))

    def add_receptors_to_db(self):
        for rec, rec_name in self.receptors:
            # NOTE: in current implementation, only one receptor allowed per database
            # Check that any receptor row is incomplete (needs receptor blob) before inserting
            filled_receptor_rows = self.dbman.get_number_filled_receptor_rows()
            if filled_receptor_rows != 0:
                raise RuntimeError("Expected Receptors table to have no receptor objects present, already has {0} receptor present. Cannot add more than 1 receptor to a database.".format(filled_receptor_rows))
            self.dbman.add_receptor_object_to_row(rec, rec_name)
