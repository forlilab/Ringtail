#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail command line tool end-to-end testing
#

import sqlite3
import os
import pytest

@pytest.fixture
def countrows():
    def __dbconnect(query):
        conn = sqlite3.connect("output.db")
        curs = conn.cursor()
        curs.execute(query)
        count = curs.fetchone()[0]
        curs.close()
        conn.close()
        os.system("rm output.db")
        return count
    return __dbconnect

class TestInputs:
    os.system("rm output.db")

    #TODO
    def test_cmdline_default(self):
        # run command line with minimal stuff, and read parsed opts
        assert True
    
    #TODO
    def test_cmdline_config_file(self):
        # run command line with config file only and it should pass (may need read or write), and read parsed_opts
        assert True
    
    #TODO
    def test_organization_of_parsed_opts(self):
        # print parsed dictionaries and assert composition
        assert True

    #TODO
    def test_db_file(self):
        # add file that exists, check for append
        # if append, check poses and docking mode consistent (read error message)
        assert True
       
    def test_files(self, countrows):
        os.system(
                "python ../scripts/rt_process_vs.py write -d --file test_data/group1/127458.dlg.gz --file test_data/group1/173101.dlg.gz --file test_data/group1/100729.dlg.gz"
            )
        count1 = countrows("SELECT COUNT(*) FROM Ligands")
        
        os.system(
            "python ../scripts/rt_process_vs.py write -d --file test_data/group1/127458.dlg.gz test_data/group1/173101.dlg.gz --file test_data/group1/100729.dlg.gz"
        )
        count2 = countrows("SELECT COUNT(*) FROM Ligands")

        assert count1 == count2 == 3

    def test_file_paths(self, countrows):
        os.system(
            "python ../scripts/rt_process_vs.py write -d --file_path test_data/group1 --file_path test_data/group2"
        )
        count1 = countrows("SELECT COUNT(*) FROM Ligands")

        os.system(
            "python ../scripts/rt_process_vs.py write -d --file_path test_data/group1 test_data/group2"
        )
        count2 = countrows("SELECT COUNT(*) FROM Ligands")

        assert count1 == count2 == 217

    def test_file_list(self, countrows):
        os.system(
            "python ../scripts/rt_process_vs.py write -d --file_list filelist1.txt --file_list filelist2.txt"
        )
        count1 = countrows("SELECT COUNT(*) FROM Ligands")

        os.system(
            "python ../scripts/rt_process_vs.py write -d --file_list filelist1.txt filelist2.txt"
        )
        count2 = countrows("SELECT COUNT(*) FROM Ligands")

        assert count1 == count2 == 5

    def test_all_file_inputs(self, countrows):
        os.system(
            "python ../scripts/rt_process_vs.py write -d --file_list filelist1.txt --file test_data/group2/361056.dlg.gz test_data/group2/53506.dlg.gz --file_path test_data/group3"
        )
        count = countrows("SELECT COUNT(*) FROM Ligands")

        assert count == 75

    #TODO this is a biggie
    def test_append_results(self, countrows):
        os.system(
            "python ../scripts/rt_process_vs.py write -d --file_path test_data/group1"
        )
        os.system(
            "python ../scripts/rt_process_vs.py write -d --input_db output.db --file_path test_data/group2 --append_results"
        )
        count = countrows("SELECT COUNT(*) FROM Ligands")

        assert count == 217

    def test_duplicate_handling(self, countrows):
        os.system(
            "python ../scripts/rt_process_vs.py write -d --file_path test_data/group1"
        )
        os.system(
            "python ../scripts/rt_process_vs.py write -d --input_db output.db --file_path test_data/group1 --append_results --duplicate_handling ignore"
        )
        count = countrows("SELECT COUNT(*) FROM Ligands")
        assert count == 138

    def test_save_rec_file(self, countrows):
        
        os.system(
            "python ../scripts/rt_process_vs.py write -d --file_list filelist1.txt --receptor_file test_data/4j8m.pdbqt --save_receptor"
        )
        count = countrows("SELECT COUNT(*) FROM Receptors WHERE receptor_object NOT NULL")
        assert count == 1
    
    def test_save_rec_file_gz(self, countrows):
        os.system(
            "python ../scripts/rt_process_vs.py write -d --file_list filelist1.txt --receptor_file test_data/4j8m.pdbqt.gz --save_receptor"
        )
        count = countrows("SELECT COUNT(*) FROM Receptors WHERE receptor_object NOT NULL")
        assert count == 1

class TestOutputs:
    def test_export_bookmark_csv(self):
        status1 = os.system(
            "python ../scripts/rt_process_vs.py write -d --file_list filelist1.txt"
        )
        status2 = os.system(
            "python ../scripts/rt_process_vs.py read -d --input_db output.db --export_bookmark_csv Ligands"
        )

        assert status1 == status2 == 0
        assert os.path.exists("Ligands.csv")
        os.system("rm output.db")
        os.system("rm Ligands.csv")

    def test_export_query_csv(self):
        status1 = os.system(
            "python ../scripts/rt_process_vs.py write -d --file_list filelist1.txt"
        )

        status2 = os.system(
            "python ../scripts/rt_process_vs.py read -d --input_db output.db --export_query_csv 'SELECT * FROM Results'"
        )

        assert status1 == status2 == 0
        assert os.path.exists("query.csv")
        os.system("rm output.db")
        os.system("rm query.csv")

    def test_interaction_tolerance(self):
        status_notol = os.system(
            "python ../scripts/rt_process_vs.py write -d --file test_data/group1/127458.dlg.gz"
        )

        conn = sqlite3.connect("output.db")
        cur = conn.cursor()

        cur.execute(
            "SELECT * FROM Interaction_bitvectors WHERE Pose_ID in (SELECT Pose_ID FROM Results WHERE LigName LIKE '127458' AND run_number = 13)"
        )
        count_notol = sum([1 for interaction in cur.fetchone() if interaction == 1])

        cur.close()
        conn.close()

        os.system("rm output.db")

        status_tol = os.system(
            "python ../scripts/rt_process_vs.py write -d --file test_data/group1/127458.dlg.gz --interaction_tolerance"
        )

        conn = sqlite3.connect("output.db")
        cur = conn.cursor()
        cur.execute(
            "SELECT * FROM Interaction_bitvectors WHERE Pose_ID in (SELECT Pose_ID FROM Results WHERE LigName LIKE '127458' AND run_number = 13)"
        )
        count_tol = sum([1 for interaction in cur.fetchone() if interaction == 1])

        cur.close()
        conn.close()

        os.system("rm output.db")

        status_tol2 = os.system(
            "python ../scripts/rt_process_vs.py write -d --file test_data/group1/127458.dlg.gz --interaction_tolerance 2.0"
        )

        conn = sqlite3.connect("output.db")
        cur = conn.cursor()
        cur.execute(
            "SELECT * FROM Interaction_bitvectors WHERE Pose_ID in (SELECT Pose_ID FROM Results WHERE LigName LIKE '127458' AND run_number = 13)"
        )
        count_tol2 = sum([1 for interaction in cur.fetchone() if interaction == 1])
        
        cur.close()
        conn.close()

        os.system("rm output.db")

        assert status_notol == 0
        assert status_tol == 0
        assert status_tol2 == 0
        assert (
            count_notol != count_tol
            or count_tol2 != count_tol
            or count_tol2 != count_notol
        )

    def test_max_poses(self):
        status3 = os.system(
            "python ../scripts/rt_process_vs.py write -d --file_list filelist1.txt"
        )
        conn = sqlite3.connect("output.db")
        cur = conn.cursor()
        cur.execute("SELECT COUNT(*) FROM Results")
        count3 = cur.fetchone()[0]

        cur.close()
        conn.close()

        os.system("rm output.db")

        status1 = os.system(
            "python ../scripts/rt_process_vs.py write -d --file_list filelist1.txt --max_poses 1"
        )
        conn = sqlite3.connect("output.db")
        cur = conn.cursor()
        cur.execute("SELECT COUNT(*) FROM Results")
        count1 = cur.fetchone()[0]

        cur.execute("SELECT COUNT(*) FROM Ligands")
        ligcount1 = cur.fetchone()[0]

        cur.close()
        conn.close()

        os.system("rm output.db")

        status5 = os.system(
            "python ../scripts/rt_process_vs.py write -d --file_list filelist1.txt --max_poses 5"
        )
        conn = sqlite3.connect("output.db")
        cur = conn.cursor()
        cur.execute("SELECT COUNT(*) FROM Results")
        count5 = cur.fetchone()[0]

        cur.close()
        conn.close()

        os.system("rm output.db")

        assert status1 == 0
        assert status3 == 0
        assert status5 == 0
        assert (count1 < count3) or (count1 < count5)
        assert count1 == ligcount1

    def test_store_all(self):
        status = os.system(
            "python ../scripts/rt_process_vs.py write -d --file_list filelist1.txt --store_all_poses"
        )
        conn = sqlite3.connect("output.db")
        cur = conn.cursor()
        cur.execute("SELECT COUNT(*) FROM Results")
        count = cur.fetchone()[0]

        cur.execute("SELECT COUNT(*) FROM Ligands")
        ligcount = cur.fetchone()[0]

        cur.close()
        conn.close()

        os.system("rm output.db")

        assert status == 0
        assert count == ligcount * 20

class TestFilters:
    #TODO
    def test_complex_filters_error(self):
        # check that filters are parsed right
        # check that wrong filter formats causes the correct error
        assert True
    
    #TODO
    def test_interaction_count(self):
        # this was not tested before, and was not working
        assert True
    
    def test_eworst(self):
        status1 = os.system(
            "python ../scripts/rt_process_vs.py write -d --file_list filelist1.txt"
        )
        status2 = os.system(
            "python ../scripts/rt_process_vs.py read -d --input_db output.db --eworst -15"
        )

        os.system("rm output.db")
        assert status1 == status2 == 0

    def test_ebest(self):
        status1 = os.system(
            "python ../scripts/rt_process_vs.py write -d --file_list filelist1.txt"
        )
        status2 = os.system(
            "python ../scripts/rt_process_vs.py read -d --input_db output.db --ebest -15"
        )

        os.system("rm output.db")

        assert status1 == status2 == 0

    def test_leworst(self):
        status1 = os.system(
            "python ../scripts/rt_process_vs.py write -d --file_list filelist1.txt"
        )
        status2 = os.system(
            "python ../scripts/rt_process_vs.py read -d --input_db output.db --leworst -0.4"
        )

        os.system("rm output.db")

        assert status1 == status2 == 0

    def test_lebest(self):
        status1 = os.system(
            "python ../scripts/rt_process_vs.py write -d --file_list filelist1.txt"
        )
        status2 = os.system(
            "python ../scripts/rt_process_vs.py read -d --input_db output.db --leworst -0.4"
        )

        os.system("rm output.db")

        assert status1 == status2 == 0

    def test_epercentile(self):
        status1 = os.system(
            "python ../scripts/rt_process_vs.py write -d --file_list filelist1.txt"
        )
        status2 = os.system(
            "python ../scripts/rt_process_vs.py read -d --input_db output.db --score_percentile 0.1"
        )

        os.system("rm output.db")

        assert status1 == status2 == 0

    def test_lepercentile(self):
        status1 = os.system(
            "python ../scripts/rt_process_vs.py write -d --file_list filelist1.txt"
        )
        status2 = os.system(
            "python ../scripts/rt_process_vs.py read -d --input_db output.db --le_percentile 0.1"
        )

        os.system("rm output.db")

        assert status1 == status2 == 0

    def test_epercentile_eworst(self):
        status1 = os.system(
            "python ../scripts/rt_process_vs.py write -d --file_list filelist1.txt"
        )
        status2 = os.system(
            "python ../scripts/rt_process_vs.py read -d --input_db output.db --score_percentile 0.1 --eworst -14"
        )

        os.system("rm output.db")

        assert status1 == status2 == 0

    def test_lepercentile_leworst(self):
        status1 = os.system(
            "python ../scripts/rt_process_vs.py write -d --file_list filelist1.txt"
        )
        status2 = os.system(
            "python ../scripts/rt_process_vs.py read -d --input_db output.db --le_percentile 0.1 --leworst -0.4"
        )

        os.system("rm output.db")

        assert status1 == status2 == 0

    def test_name(self):
        status1 = os.system(
            "python ../scripts/rt_process_vs.py write -d --file_list filelist1.txt"
        )
        status2 = os.system(
            "python ../scripts/rt_process_vs.py read -d --input output.db --ligand_name 127458"
        )

        os.system("rm output.db")

        assert status1 == status2 == 0

    def test_hbcount(self):
        status1 = os.system(
            "python ../scripts/rt_process_vs.py write -d --file_list filelist1.txt"
        )
        status2 = os.system(
            "python ../scripts/rt_process_vs.py read -d --input_db output.db --hb_count 5"
        )

        assert status1 == status2 == 0

    """def test_react_any(self):
        status = os.system(
            "python ../scripts/rt_process_vs.py read -d --input_db output.db --react_any"
        )

        assert status == 0

    def test_react1(self):
        status = os.system(
            "python ../scripts/rt_process_vs.py read -d --input_db output.db  --reactive_interactions A:TYR:169:"
        )

        assert status == 0

    def test_react2(self):
        status = os.system(
            "python ../scripts/rt_process_vs.py read -d --input_db output.db  --reactive_interactions :TYR:169:"
        )

        assert status == 0

    def test_react3(self):
        status = os.system(
            "python ../scripts/rt_process_vs.py read -d --input_db output.db --reactive_interactions :TYR::"
        )

        assert status == 0

    def test_react4(self):
        status = os.system(
            "python ../scripts/rt_process_vs.py read -d --input_db output.db --reactive_interactions A:TYR::"
        )

        assert status == 0

    def test_react5(self):
        status = os.system(
            "python ../scripts/rt_process_vs.py read -d --input_db output.db --reactive_interactions A::169:"
        )

        assert status == 0

    def test_react6(self):
        status = os.system(
            "python ../scripts/rt_process_vs.py read -d --input_db output.db --reactive_interactions A:::"
        )

        assert status == 0

    def test_react7(self):
        status = os.system(
            "python ../scripts/rt_process_vs.py read -d --input_db output.db --reactive_interactions ::169:"
        )

        assert status == 0"""

    def test_hb1(self):
        status = os.system(
            "python ../scripts/rt_process_vs.py read -d --input_db output.db -hb A:LYS:162:"
        )

        assert status == 0

    def test_hb2(self):
        status = os.system(
            "python ../scripts/rt_process_vs.py read -d --input_db output.db -hb :LYS:162:"
        )

        assert status == 0

    def test_hb3(self):
        status = os.system(
            "python ../scripts/rt_process_vs.py read -d --input_db output.db -hb :LYS::"
        )

        assert status == 0

    def test_hb4(self):
        status = os.system(
            "python ../scripts/rt_process_vs.py read -d --input_db output.db -hb A:LYS::"
        )

        assert status == 0

    def test_hb5(self):
        status = os.system(
            "python ../scripts/rt_process_vs.py read -d --input_db output.db -hb A::162:"
        )

        assert status == 0

    def test_hb6(self):
        status = os.system(
            "python ../scripts/rt_process_vs.py read -d --input_db output.db -hb A:::"
        )

        assert status == 0

    def test_hb7(self):
        status = os.system(
            "python ../scripts/rt_process_vs.py read -d --input_db output.db -hb ::162:"
        )

        assert status == 0

    def test_vdw1(self):
        status = os.system(
            "python ../scripts/rt_process_vs.py read -d --input_db output.db -vdw A:VAL:279:"
        )

        assert status == 0

    def test_vdw2(self):
        status = os.system(
            "python ../scripts/rt_process_vs.py read -d --input_db output.db -vdw :VAL:279:"
        )

        assert status == 0

    def test_vdw3(self):
        status = os.system(
            "python ../scripts/rt_process_vs.py read -d --input_db output.db -vdw :VAL::"
        )

        assert status == 0

    def test_vdw4(self):
        status = os.system(
            "python ../scripts/rt_process_vs.py read -d --input_db output.db -vdw A:VAL::"
        )

        assert status == 0

    def test_vdw5(self):
        status = os.system(
            "python ../scripts/rt_process_vs.py read -d --input_db output.db -vdw A::279:"
        )

        assert status == 0

    def test_vdw6(self):
        status = os.system(
            "python ../scripts/rt_process_vs.py read -d --input_db output.db -vdw A:::"
        )

        assert status == 0

    def test_vdw7(self):
        status = os.system(
            "python ../scripts/rt_process_vs.py read -d --input_db output.db -vdw ::279:"
        )

        assert status == 0
    
    def test_all_filters(self):
        status = os.system(
            "python ../scripts/rt_process_vs.py read -d --input_db output.db --eworst -15 --ebest -16 --leworst -0.4 --lebest -0.5 --score_percentile 99 --le_percentile 99 --ligand_name 127458 --hb_count 5 --react_any -hb A:LYS:162: -vdw A:VAL:279: --reactive_interactions A:TYR:169:"
        )

        assert status == 0
    
    def test_export_sdf(self):
        status1 = os.system(
            "python ../scripts/rt_process_vs.py write -d --file_list filelist1.txt"
        )
        status2 = os.system(
            "python ../scripts/rt_process_vs.py read --input_db output.db -e -4 -sdf . -d "
        )

        import glob
        sdf_files = glob.glob("*.sdf")
        expected = ['127458.sdf', '100729.sdf', '173101.sdf']
        assert len(sdf_files) == len(expected)
        for f in sdf_files:
            assert f in expected
            os.remove(f)

        assert status2 == 0
    
        os.system("rm output_log.txt output.db")
