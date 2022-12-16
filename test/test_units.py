#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail unit testing
#
from ringtail import StorageManagerSQLite

class Test_StorageManSQLite:

    def test_fetch_summary_data(self):
        with StorageManagerSQLite("output.db") as dbman:
            summ_dict = dbman.fetch_summary_data()
            assert summ_dict == {'num_ligands': 287, 'num_poses': 645, 'num_unique_interactions': 183, 'min_energies_binding': -7.93, 'max_energies_binding': -2.03, '1%_energies_binding': -7.43, '10%_energies_binding': -6.46, 'min_leff': -0.6183333333333333, 'max_leff': -0.13277777777777777, '1%_leff': -0.581, '10%_leff': -0.4653846153846154}