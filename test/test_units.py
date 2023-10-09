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
        os.system("rm output.db")
        status1 = os.system(
            "python ../scripts/rt_process_vs.py write -d --file_list filelist1.txt"
        )
        with StorageManagerSQLite("output.db") as dbman:
            summ_dict = dbman.fetch_summary_data()
            print(summ_dict)
            assert summ_dict == {'num_ligands': 3, 'num_poses': 7, 'num_unique_interactions': 57, 'num_interacting_residues': 30, 'min_docking_score': -6.66, 'max_docking_score': -4.98, '1%_docking_score': -6.66, '10%_docking_score': -6.66, 'min_leff': -0.444, 'max_leff': -0.35000000000000003, '1%_leff': -0.444, '10%_leff': -0.444}

    def test_bookmark_info(self):
        os.system("rm output.db")

        opts = RingtailCore.get_defaults()
        opts["rman_opts"]["values"]["file_sources"]["file_path"]["path"] = [["test_data/"]]
        opts["rman_opts"]["values"]["file_sources"]["file_path"]["recursive"] = True

        opts["filters"]["values"] = {'eworst': -3, 'ebest': None, 'leworst': None, 'lebest': None, 'score_percentile': None, 'le_percentile': None, 'vdw_interactions': [('A:ARG:123:', True), ('A:VAL:124:', True)], 'hb_interactions': [('A:ARG:123:', True)], 'reactive_interactions': [], 'interactions_count': [], 'react_any': None, 'max_miss': 0, 'ligand_name': [], 'ligand_substruct': [], 'ligand_substruct_pos': [], 'ligand_max_atoms': None, 'ligand_operator': 'OR'}

        with RingtailCore(opts_dict=opts) as rt_core:
            rt_core.add_results()
            rt_core.filter()

        conn = sqlite3.connect("output.db")
        cur = conn.cursor()
        bookmark = cur.execute("SELECT filters FROM Bookmarks WHERE Bookmark_name LIKE 'passing_results'")
        bookmark_filters_db_str = bookmark.fetchone()[0]
        print(bookmark_filters_db_str)
        filters = {'eworst': -3, 'ebest': None, 'leworst': None, 'lebest': None, 'score_percentile': None, 'le_percentile': None, 'vdw_interactions': [('A:ARG:123:', True), ('A:VAL:124:', True)], 'hb_interactions': [('A:ARG:123:', True)], 'reactive_interactions': [], 'interactions_count': [], 'react_any': None, 'max_miss': 0, 'ligand_name': [], 'ligand_substruct': [], 'ligand_substruct_pos': [], 'ligand_max_atoms': None, 'ligand_operator': 'OR'}
        assert bookmark_filters_db_str == json.dumps(filters)
        cur.close()
        conn.close()

    def test_version_info(self):
        os.system("rm output.db")

        opts = RingtailCore.get_defaults()
        opts["rman_opts"]["values"]["file_sources"]["file_path"]["path"] = [["test_data/"]]
        opts["rman_opts"]["values"]["file_sources"]["file_path"]["recursive"] = True

        with RingtailCore(opts_dict=opts) as rt_core:
            rt_core.add_results()
            versionmatch, version = rt_core.storageman.check_ringtaildb_version()
        
        assert versionmatch
        assert int(version) == 110  # TODO: update for new versions

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