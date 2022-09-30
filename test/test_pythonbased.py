import sqlite3
import os
from ringtail import RingtailCore
from ringtail import RTCoreError
import pytest

class TestWriteOpts:
    def test_set_opts(self):
        os.system("rm output.db")

        opts = RingtailCore.get_defaults()
        opts["storage_opts"]["values"]["storage_type"] = "sqlite"
        opts["rman_opts"]["values"]["file_sources"]["file_path"]["path"] = [["test_data/"]]
        opts["rman_opts"]["values"]["file_sources"]["file_path"]["recursive"] = True
        
        expected_opts = {'storage_opts': {'values': {'append_results': False, 'order_results': None, 'output_all_poses': None, 'results_view_name': 'passing_results', 'overwrite': None, 'conflict_opt': None, 'mode': 'ADGPU', 'db_file': 'output.db',}, 'ignore': []}, 'rman_opts': {'values': {'parser_manager': 'multiprocessing', 'mode': 'dlg', 'chunk_size': 1, 'max_poses': 3, 'store_all_poses': False, 'interaction_tolerance': None, 'target': None, 'add_interactions': False, 'interaction_cutoffs': [3.7, 4.0], 'receptor_file': None, 'file_sources': {'file': [[]], 'file_path': {'path': [['test_data/']], 'pattern': '*.dlg*', 'recursive': True}, 'file_list': [[]]}, 'file_pattern': '*.dlg*', 'max_proc': None}, 'ignore': []}, 'out_opts': {'values': {'log': 'output_log.txt', 'overwrite': None, 'export_sdf_path': None, 'plot': None, 'outfields': 'e', 'no_print': True, 'export_bookmark_csv': None, 'export_query_csv': None, 'export_bookmark_db': None, 'data_from_bookmark': None, 'filter_bookmark': None}, 'ignore': []}, 'filters': {'values': {'properties': {'eworst': None, 'ebest': None, 'leworst': None, 'lebest': None, 'energy_percentile': None, 'le_percentile': None}, 'interactions': {'V': [], 'H': [], 'R': []}, 'interactions_count': [], 'ligand_filters': {'N': [], 'S': [], 'F': 'OR'}, 'filter_ligands_flag': False, 'max_miss': 0, 'react_any': None}, 'ignore': []}}

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

    def test_set_max_proc(self):
        with pytest.raises(RTCoreError):
            os.system("rm output.db")

            opts = RingtailCore.get_defaults()
            opts["storage_opts"]["values"]["storage_type"] = "sqlite"
            opts["rman_opts"]["values"]["file_sources"]["file_path"]["path"] = [["test_data/"]]
            opts["rman_opts"]["values"]["file_sources"]["file_path"]["recursive"] = True
            opts["rman_opts"]["values"]["max_proc"] = "1"
            with RingtailCore(**opts) as rt_core:
                rt_core.add_results()

    def test_hb_count_python(self):
        pass

