#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail unit testing
#
from ringtail import StorageManagerSQLite, RingtailCore, Filters
import sqlite3
import os
import json

class Test_StorageManSQLite:

    def test_fetch_summary_data(self):
        with StorageManagerSQLite("output.db") as dbman:
            summ_dict = dbman.fetch_summary_data()
            assert summ_dict == {'num_ligands': 287, 'num_poses': 645, 'num_unique_interactions': 183, 'min_docking_score': -7.93, 'max_docking_score': -2.03, '1%_docking_score': -7.43, '10%_docking_score': -6.46, 'min_leff': -0.6183333333333333, 'max_leff': -0.13277777777777777, '1%_leff': -0.581, '10%_leff': -0.4653846153846154, 'num_interacting_residues': 82}

    def test_bookmark_info(self):
        os.system("rm output.db")

        opts = RingtailCore.get_defaults()
        opts["storage_opts"]["values"]["storage_type"] = "sqlite"
        opts["rman_opts"]["values"]["file_sources"]["file_path"]["path"] = [["test_data/"]]
        opts["rman_opts"]["values"]["file_sources"]["file_path"]["recursive"] = True

        filters = {'properties': {'eworst': -15.0, 'ebest': -16.0, 'leworst': -0.4, 'lebest': -0.5, 'score_percentile': None, 'le_percentile': None}, 'interactions': {'V': [('A:VAL:279:', True)], 'H': [('A:LYS:162:', True)], 'R': [('A:TYR:169:', True)]}, 'interactions_count': [('hb_count', 5)], 'ligand_filters': {'N': ['127458'], 'S': [], 'F': 'OR', 'X': []}, 'filter_ligands_flag': True, 'max_miss': 0, 'react_any': True}
        opts["filters"]["values"] = filters

        with RingtailCore(**opts) as rt_core:
            rt_core.add_results()
            rt_core.filter()

        conn = sqlite3.connect("output.db")
        cur = conn.cursor()
        bookmark = cur.execute("SELECT filters FROM Bookmarks WHERE Bookmark_name LIKE 'passing_results'")
        bookmark_filters_db_str = bookmark.fetchone()[0]
        assert bookmark_filters_db_str == json.dumps(filters)
        cur.close()
        conn.close()

class Test_RingtailCore:

    def test_prepare_filters_for_storageman(self):
        f = Filters(hb_interactions=[("A:ARG:123:", True), ("A:VAL:124:", True)], vdw_interactions=[("A:ARG:123:", True), ("A:VAL:124:", True)])
        test_filters = []
        with RingtailCore() as rtc:
            rtc.filters = f
            interaction_combs = rtc._generate_interaction_combinations(1)
            for ic in interaction_combs:
                test_filters.append(rtc._prepare_filters_for_storageman(ic))

        assert {'eworst': None, 'ebest': None, 'leworst': None, 'lebest': None, 'score_percentile': None, 'le_percentile': None, 'vdw_interactions': [('A:ARG:123:', True), ('A:VAL:124:', True)], 'hb_interactions': [('A:ARG:123:', True)], 'reactive_interactions': [], 'interactions_count': [], 'react_any': None, 'max_miss': 0, 'ligand_name': [], 'ligand_substruct': [], 'ligand_substruct_pos': [], 'ligand_max_atoms': None, 'ligand_operator': 'OR'} in test_filters

        assert {'eworst': None, 'ebest': None, 'leworst': None, 'lebest': None, 'score_percentile': None, 'le_percentile': None, 'vdw_interactions': [('A:ARG:123:', True), ('A:VAL:124:', True)], 'hb_interactions': [('A:VAL:124:', True)], 'reactive_interactions': [], 'interactions_count': [], 'react_any': None, 'max_miss': 0, 'ligand_name': [], 'ligand_substruct': [], 'ligand_substruct_pos': [], 'ligand_max_atoms': None, 'ligand_operator': 'OR'} in test_filters

        assert {'eworst': None, 'ebest': None, 'leworst': None, 'lebest': None, 'score_percentile': None, 'le_percentile': None, 'vdw_interactions': [('A:ARG:123:', True)], 'hb_interactions': [('A:ARG:123:', True), ('A:VAL:124:', True)], 'reactive_interactions': [], 'interactions_count': [], 'react_any': None, 'max_miss': 0, 'ligand_name': [], 'ligand_substruct': [], 'ligand_substruct_pos': [], 'ligand_max_atoms': None, 'ligand_operator': 'OR'} in test_filters

        assert {'eworst': None, 'ebest': None, 'leworst': None, 'lebest': None, 'score_percentile': None, 'le_percentile': None, 'vdw_interactions': [('A:VAL:124:', True)], 'hb_interactions': [('A:ARG:123:', True), ('A:VAL:124:', True)], 'reactive_interactions': [], 'interactions_count': [], 'react_any': None, 'max_miss': 0, 'ligand_name': [], 'ligand_substruct': [], 'ligand_substruct_pos': [], 'ligand_max_atoms': None, 'ligand_operator': 'OR'} in test_filters

        assert {'eworst': None, 'ebest': None, 'leworst': None, 'lebest': None, 'score_percentile': None, 'le_percentile': None, 'vdw_interactions': [('A:ARG:123:', True), ('A:VAL:124:', True)], 'hb_interactions': [('A:ARG:123:', True), ('A:VAL:124:', True)], 'reactive_interactions': [], 'interactions_count': [], 'react_any': None, 'max_miss': 0, 'ligand_name': [], 'ligand_substruct': [], 'ligand_substruct_pos': [], 'ligand_max_atoms': None, 'ligand_operator': 'OR'} in test_filters

        assert len(test_filters) == 5