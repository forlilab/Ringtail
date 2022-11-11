import sqlite3
import os
from ringtail import RingtailCore
from ringtail import RTCoreError
from ringtail import StorageManagerSQLite
import pytest

class TestWriteOpts:
    def test_set_opts(self):
        os.system("rm output.db")

        opts = RingtailCore.get_defaults()
        opts["storage_opts"]["values"]["storage_type"] = "sqlite"
        opts["rman_opts"]["values"]["file_sources"]["file_path"]["path"] = [["test_data/"]]
        opts["rman_opts"]["values"]["file_sources"]["file_path"]["recursive"] = True
        
        expected_opts = {'storage_opts': {'values': {'append_results': False, 'order_results': None, 'output_all_poses': None, 'results_view_name': 'passing_results', 'overwrite': None, 'conflict_opt': None, 'mode': 'ADGPU', 'db_file': 'output.db',}, 'ignore': []}, 'rman_opts': {'values': {'parser_manager': 'multiprocessing', "storageman_class": StorageManagerSQLite, 'mode': 'dlg', 'chunk_size': 1, 'max_poses': 3, 'store_all_poses': False, 'interaction_tolerance': None, 'target': None, 'add_interactions': False, 'interaction_cutoffs': [3.7, 4.0], 'receptor_file': None, 'file_sources': {'file': [[]], 'file_path': {'path': [['test_data/']], 'pattern': '*.dlg*', 'recursive': True}, 'file_list': [[]]}, 'file_pattern': '*.dlg*', 'max_proc': None}, 'ignore': []}, 'out_opts': {'values': {'log': 'output_log.txt', 'export_sdf_path': None, 'plot': None, 'outfields': 'e', 'no_print': True, 'export_bookmark_csv': None, 'export_query_csv': None, 'export_bookmark_db': None, 'data_from_bookmark': None, 'filter_bookmark': None}, 'ignore': []}, 'filters': {'values': {'properties': {'eworst': None, 'ebest': None, 'leworst': None, 'lebest': None, 'energy_percentile': None, 'le_percentile': None}, 'interactions': {'V': [], 'H': [], 'R': []}, 'interactions_count': [], 'ligand_filters': {'N': [], 'S': [], 'F': 'OR'}, 'filter_ligands_flag': False, 'max_miss': 0, 'react_any': None}, 'ignore': []}}

        with RingtailCore(**opts) as rt_core:
            assert rt_core.storage_opts == expected_opts["storage_opts"]["values"]
            dbman = rt_core.rman_opts["storageman"]
            expected_opts["rman_opts"]["values"]["storageman"] = dbman
            assert rt_core.rman_opts == expected_opts["rman_opts"]["values"]
            assert rt_core.filters == expected_opts["filters"]["values"]
            assert rt_core.out_opts == expected_opts["out_opts"]["values"]

            rt_core.add_results()

        conn = sqlite3.connect("output.db")
        cur = conn.cursor()
        cur.execute("SELECT COUNT(*) FROM Ligands")
        ligcount = cur.fetchone()[0]

        cur.execute("SELECT COUNT(*) FROM Results")
        posecount = cur.fetchone()[0]

        cur.close()
        conn.close()

        assert ligcount == 287
        assert posecount == 645
   
class TestReadOpts:
    def test_all_filters_python(self):
        os.system("rm output.db")

        opts = RingtailCore.get_defaults()
        opts["storage_opts"]["values"]["storage_type"] = "sqlite"
        opts["rman_opts"]["values"]["file_sources"]["file_path"]["path"] = [["test_data/"]]
        opts["rman_opts"]["values"]["file_sources"]["file_path"]["recursive"] = True

        """filters = {'properties': {'eworst': -15.0, 'ebest': -16.0, 'leworst': -0.4, 'lebest': -0.5, 'energy_percentile': None, 'le_percentile': None}, 'interactions': {'V': [('A:VAL:279:', True)], 'H': [('A:LYS:162:', True)], 'R': [('A:TYR:169:', True)]}, 'interactions_count': [('hb_count', 5)], 'ligand_filters': {'N': ['127458'], 'S': [], 'F': 'OR'}, 'filter_ligands_flag': True, 'max_miss': 0, 'react_any': True}"""
        filters = {'properties': {'eworst': -6.0, 'ebest': None, 'leworst': None, 'lebest': None, 'energy_percentile': None, 'le_percentile': None}, 'interactions': {'V': [], 'H': [], 'R': []}, 'interactions_count': [], 'ligand_filters': {'N': [], 'S': [], 'F': 'OR'}, 'filter_ligands_flag': False, 'max_miss': 0, 'react_any': False}
        opts["filters"]["values"] = filters

        with RingtailCore(**opts) as rt_core:
            rt_core.add_results()
            rt_core.filter()

        conn = sqlite3.connect("output.db")
        cur = conn.cursor()
        cur.execute("SELECT COUNT(*) FROM passing_results")
        bookmarkcount = cur.fetchone()[0]

        cur.close()
        conn.close()
        
        assert bookmarkcount == 65

