#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail unit testing
#
from ringtail import StorageManagerSQLite, RingtailCore
import sqlite3
import os

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
        bookmark = cur.execute("SELECT * FROM Bookmarks WHERE Bookmark_name LIKE 'passing_results'")
        bookmark_tuple = list(bookmark.fetchone())
        assert tuple([bookmark_tuple[0]] + bookmark_tuple[3:]) == ('passing_results', -15.0, -16.0, -0.4, -0.5, None, None, '["127458"]', None, None, None, '[["V", "A", "VAL", "279", "", true], ["H", "A", "LYS", "162", "", true], ["R", "A", "TYR", "169", "", true], ["R", "", "", "", "", true]]', 5)
        cur.close()
        conn.close()

    def test_set_opts(self):
        opts = RingtailCore.get_defaults()
        RingtailCore.set_opts(opts, "storage_type", "sqlite")
        RingtailCore.set_opts(opts, "path", [["test_data/"]])
        RingtailCore.set_opts(opts, "recursive", True)

        assert opts["storage_opts"]["values"]["storage_type"] == "sqlite"
        assert opts["rman_opts"]["values"]["file_sources"]["file_path"]["path"] == [["test_data/"]]
        assert opts["rman_opts"]["values"]["file_sources"]["file_path"]["recursive"] == True

        opts2 = RingtailCore.get_defaults()
        RingtailCore.set_opts(opts2, ["storage_type","path","recursive"], ["sqlite", [["test_data/"]], True])

        assert opts2["storage_opts"]["values"]["storage_type"] == "sqlite"
        assert opts2["rman_opts"]["values"]["file_sources"]["file_path"]["path"] == [["test_data/"]]
        assert opts2["rman_opts"]["values"]["file_sources"]["file_path"]["recursive"] == True

    def test_get_opts(self):
        opts = RingtailCore.get_defaults()
        opts["storage_opts"]["values"]["storage_type"] = "sqlite"
        opts["rman_opts"]["values"]["file_sources"]["file_path"]["path"] = [["test_data/"]]
        opts["rman_opts"]["values"]["file_sources"]["file_path"]["recursive"] = True

        assert RingtailCore.get_opts(opts, "storage_type")[0] == "sqlite"
        assert RingtailCore.get_opts(opts, "path")[0] == [["test_data/"]]
        assert RingtailCore.get_opts(opts, "recursive")[0] == True

        opts2 = RingtailCore.get_defaults()
        opts2["storage_opts"]["values"]["storage_type"] = "sqlite"
        opts2["rman_opts"]["values"]["file_sources"]["file_path"]["path"] = [["test_data/"]]
        opts2["rman_opts"]["values"]["file_sources"]["file_path"]["recursive"] = True

        assert RingtailCore.get_opts(opts2, ["storage_type", "path", "recursive"]) == ["sqlite", [["test_data/"]], True]

    def test_set_filters(self):
        opts = RingtailCore.get_defaults()
        RingtailCore.set_filters(opts, "eworst", -3.0)
        assert opts["filters"]["values"]["properties"]["eworst"] == -3.0
        RingtailCore.set_filters(opts, "V", ["cucumber"])
        assert opts["filters"]["values"]["interactions"]["V"] == ["cucumber"]
        RingtailCore.set_filters(opts, "N", ["ligand1"])
        assert opts["filters"]["values"]["ligand_filters"]["N"] == ["ligand1"]
        assert opts["filters"]["values"]["filter_ligands_flag"] == True
        RingtailCore.set_filters(opts, "max_miss", 3)
        assert opts["filters"]["values"]["max_miss"] == 3

        opts2 = RingtailCore.get_defaults()
        RingtailCore.set_filters(opts2, ["eworst", "V", "N", "max_miss"], [-3.0, ["cucumber"], ["ligand1"], 3])
        assert opts2["filters"]["values"]["properties"]["eworst"] == -3.0
        assert opts2["filters"]["values"]["interactions"]["V"] == ["cucumber"]
        assert opts2["filters"]["values"]["ligand_filters"]["N"] == ["ligand1"]
        assert opts2["filters"]["values"]["filter_ligands_flag"] == True
        assert opts2["filters"]["values"]["max_miss"] == 3

    def test_get_filters(self):
        opts = RingtailCore.get_defaults()
        
        opts["filters"]["values"]["properties"]["eworst"] = -3.0
        opts["filters"]["values"]["interactions"]["V"] = ["cucumber"]
        opts["filters"]["values"]["ligand_filters"]["N"] = ["ligand1"]
        opts["filters"]["values"]["max_miss"] = 3

        assert RingtailCore.get_filters(opts, "eworst")[0] == -3.0
        assert RingtailCore.get_filters(opts, "V")[0] == ["cucumber"]
        assert RingtailCore.get_filters(opts, "N")[0] == ["ligand1"]
        assert RingtailCore.get_filters(opts, "max_miss")[0] == 3

        assert RingtailCore.get_filters(opts, ["eworst", "V", "N", "max_miss"]) == [-3.0, ["cucumber"], ["ligand1"], 3]
