#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail end-to-end testing
#

import sqlite3
import os


class TestInputs:

    def test_multiple_files1(self):
        os.system(
            "python ../scripts/run_ringtail.py --file test_data/TEST_0000003-results/OB3Z3759450716_RX1--4xfx_mon_prep--tyr169.dlg.gz --file test_data/TEST_0000003-results/OB3Z3759440327_RX1--4xfx_mon_prep--tyr169.dlg.gz --file test_data/TEST_0000003-results/OB3Z3759457587_RX1--4xfx_mon_prep--tyr169.dlg.gz --file test_data/TEST_0000003-results/OB3Z3759444505_1_RX1--4xfx_mon_prep--tyr169.dlg.gz --file test_data/TEST_0000001-results/OB3Z3759305928_RX1--4xfx_mon_prep--tyr169.dlg.gz --no_print"
        )

        conn = sqlite3.connect("output.db")
        cur = conn.cursor()
        cur.execute("SELECT COUNT(*) FROM Ligands")
        count = cur.fetchone()[0]

        cur.close()
        conn.close()

        os.system("rm output.db")

        assert count == 5

    def test_multiple_files2(self):
        os.system(
            "python ../scripts/run_ringtail.py --file test_data/TEST_0000003-results/OB3Z3759450716_RX1--4xfx_mon_prep--tyr169.dlg.gz test_data/TEST_0000003-results/OB3Z3759440327_RX1--4xfx_mon_prep--tyr169.dlg.gz test_data/TEST_0000003-results/OB3Z3759457587_RX1--4xfx_mon_prep--tyr169.dlg.gz test_data/TEST_0000003-results/OB3Z3759444505_1_RX1--4xfx_mon_prep--tyr169.dlg.gz --file test_data/TEST_0000001-results/OB3Z3759305928_RX1--4xfx_mon_prep--tyr169.dlg.gz --no_print"
        )

        conn = sqlite3.connect("output.db")
        cur = conn.cursor()
        cur.execute("SELECT COUNT(*) FROM Ligands")
        count = cur.fetchone()[0]

        cur.close()
        conn.close()

        os.system("rm output.db")

        assert count == 5

    def test_multiple_paths1(self):
        os.system(
            "python ../scripts/run_ringtail.py --file_path test_data/TEST_0000003-results --file_path test_data/TEST_0000005-results --no_print"
        )

        conn = sqlite3.connect("output.db")
        cur = conn.cursor()
        cur.execute("SELECT COUNT(*) FROM Ligands")
        count = cur.fetchone()[0]

        cur.close()
        conn.close()

        os.system("rm output.db")

        assert count == 19999

    def test_multiple_paths2(self):
        os.system(
            "python ../scripts/run_ringtail.py --file_path test_data/TEST_0000003-results test_data/TEST_0000005-results --no_print"
        )

        conn = sqlite3.connect("output.db")
        cur = conn.cursor()
        cur.execute("SELECT COUNT(*) FROM Ligands")
        count = cur.fetchone()[0]

        cur.close()
        conn.close()

        os.system("rm output.db")

        assert count == 19999

    def test_filelist1(self):
        os.system(
            "python ../scripts/run_ringtail.py --file_list filelist1.txt --file_list filelist2.txt --no_print"
        )

        conn = sqlite3.connect("output.db")
        cur = conn.cursor()
        cur.execute("SELECT COUNT(*) FROM Ligands")
        count = cur.fetchone()[0]

        cur.close()
        conn.close()

        os.system("rm output.db")

        assert count == 5

    def test_filelist2(self):
        os.system(
            "python ../scripts/run_ringtail.py --file_list filelist1.txt filelist2.txt --no_print"
        )

        conn = sqlite3.connect("output.db")
        cur = conn.cursor()
        cur.execute("SELECT COUNT(*) FROM Ligands")
        count = cur.fetchone()[0]

        cur.close()
        conn.close()

        os.system("rm output.db")

        assert count == 5

    def test_all_input_opts(self):
        os.system(
            "mv test_data/TEST_0000003-results/OB3Z3759450716_RX1--4xfx_mon_prep--tyr169.dlg.gz test_data/OB3Z3759450716_RX1--4xfx_mon_prep--tyr169.dlg.gz"
        )
        os.system(
            "mv test_data/TEST_0000003-results/OB3Z3759440327_RX1--4xfx_mon_prep--tyr169.dlg.gz test_data/OB3Z3759440327_RX1--4xfx_mon_prep--tyr169.dlg.gz"
        )
        os.system(
            "mv test_data/TEST_0000003-results/OB3Z3759457587_RX1--4xfx_mon_prep--tyr169.dlg.gz test_data/OB3Z3759457587_RX1--4xfx_mon_prep--tyr169.dlg.gz"
        )
        os.system(
            "mv test_data/TEST_0000003-results/OB3Z3759444505_1_RX1--4xfx_mon_prep--tyr169.dlg.gz test_data/OB3Z3759444505_1_RX1--4xfx_mon_prep--tyr169.dlg.gz"
        )
        os.system(
            "mv test_data/TEST_0000001-results/OB3Z3759305928_RX1--4xfx_mon_prep--tyr169.dlg.gz test_data/OB3Z3759305928_RX1--4xfx_mon_prep--tyr169.dlg.gz"
        )

        os.system(
            "python ../scripts/run_ringtail.py --file_list filelist3.txt --file test_data/OB3Z3759444505_1_RX1--4xfx_mon_prep--tyr169.dlg.gz test_data/OB3Z3759305928_RX1--4xfx_mon_prep--tyr169.dlg.gz --file_path test_data/TEST_0000003-results --no_print"
        )

        conn = sqlite3.connect("output.db")
        cur = conn.cursor()
        cur.execute("SELECT COUNT(*) FROM Ligands")
        count = cur.fetchone()[0]

        cur.close()
        conn.close()

        os.system("rm output.db")
        os.system(
            "mv test_data/OB3Z3759450716_RX1--4xfx_mon_prep--tyr169.dlg.gz test_data/TEST_0000003-results/OB3Z3759450716_RX1--4xfx_mon_prep--tyr169.dlg.gz"
        )
        os.system(
            "mv test_data/OB3Z3759440327_RX1--4xfx_mon_prep--tyr169.dlg.gz test_data/TEST_0000003-results/OB3Z3759440327_RX1--4xfx_mon_prep--tyr169.dlg.gz"
        )
        os.system(
            "mv test_data/OB3Z3759457587_RX1--4xfx_mon_prep--tyr169.dlg.gz test_data/TEST_0000003-results/OB3Z3759457587_RX1--4xfx_mon_prep--tyr169.dlg.gz"
        )
        os.system(
            "mv test_data/OB3Z3759444505_1_RX1--4xfx_mon_prep--tyr169.dlg.gz test_data/TEST_0000003-results/OB3Z3759444505_1_RX1--4xfx_mon_prep--tyr169.dlg.gz"
        )
        os.system(
            "mv test_data/OB3Z3759305928_RX1--4xfx_mon_prep--tyr169.dlg.gz test_data/TEST_0000001-results/OB3Z3759305928_RX1--4xfx_mon_prep--tyr169.dlg.gz"
        )

        assert count == 10001

    def test_add_results(self):
        os.system(
            "python ../scripts/run_ringtail.py --file_path test_data/TEST_0000005-results --no_print"
        )
        os.system(
            "python ../scripts/run_ringtail.py --input_db output.db --file_path test_data/TEST_0000003-results --add_results --no_print"
        )

        conn = sqlite3.connect("output.db")
        cur = conn.cursor()
        cur.execute("SELECT COUNT(*) FROM Ligands")
        count = cur.fetchone()[0]

        cur.close()
        conn.close()

        os.system("rm output.db")

        assert count == 19999

    def test_conflict_handling_ign(self):
        os.system(
            "python ../scripts/run_ringtail.py --file_path test_data/TEST_0000005-results --no_print"
        )
        os.system(
            "python ../scripts/run_ringtail.py --input_db output.db --file_path test_data/TEST_0000005-results --add_results --no_print --conflict_handling ignore"
        )

        conn = sqlite3.connect("output.db")
        cur = conn.cursor()
        cur.execute("SELECT COUNT(*) FROM Ligands")
        count = cur.fetchone()[0]

        cur.close()
        conn.close()

        os.system("rm output.db")

        assert count == 9999

    def test_conflict_handling_rpl(self):
        os.system(
            "python ../scripts/run_ringtail.py --file_path test_data/TEST_0000005-results --no_print"
        )
        os.system(
            "python ../scripts/run_ringtail.py --input_db output.db --file_path test_data/TEST_0000005-results --add_results --no_print --conflict_handling replace"
        )

        conn = sqlite3.connect("output.db")
        cur = conn.cursor()
        cur.execute("SELECT COUNT(*) FROM Ligands")
        count = cur.fetchone()[0]

        cur.close()
        conn.close()

        os.system("rm output.db")

        assert count == 9999

    def test_save_rec_file(self):
        os.system(
            "python ../scripts/run_ringtail.py --file_list filelist1.txt --file test_data/4xfx_mon_prep--tyr169_rigid.pdbqt --save_receptor --no_print"
        )

        conn = sqlite3.connect("output.db")
        cur = conn.cursor()
        cur.execute(
            "SELECT COUNT(*) FROM Receptors WHERE receptor_object NOT NULL")
        count = cur.fetchone()[0]

        cur.close()
        conn.close()

        os.system("rm output.db")

        assert count == 1

    def test_save_rec_filepath(self):
        os.system(
            "python ../scripts/run_ringtail.py --file_list filelist1.txt --file_path test_data/rec --save_receptor --no_print"
        )

        conn = sqlite3.connect("output.db")
        cur = conn.cursor()
        cur.execute(
            "SELECT COUNT(*) FROM Receptors WHERE receptor_object NOT NULL")
        count = cur.fetchone()[0]

        cur.close()
        conn.close()

        os.system("rm output.db")

        assert count == 1

    def test_save_rec_filelist(self):
        os.system(
            "python ../scripts/run_ringtail.py --file_list filelist4.txt --save_receptor --no_print"
        )

        conn = sqlite3.connect("output.db")
        cur = conn.cursor()
        cur.execute(
            "SELECT COUNT(*) FROM Receptors WHERE receptor_object NOT NULL")
        count = cur.fetchone()[0]

        cur.close()
        conn.close()

        os.system("rm output.db")

        assert count == 1

    def test_save_rec_file_gz(self):
        os.system(
            "python ../scripts/run_ringtail.py --file_list filelist1.txt --file test_data/4xfx_mon_prep--tyr169_rigid.pdbqt.gz --save_receptor --no_print"
        )

        conn = sqlite3.connect("output.db")
        cur = conn.cursor()
        cur.execute(
            "SELECT COUNT(*) FROM Receptors WHERE receptor_object NOT NULL")
        count = cur.fetchone()[0]

        cur.close()
        conn.close()

        os.system("rm output.db")

        assert count == 1

    def test_save_rec_filepath_gz(self):
        os.system(
            "python ../scripts/run_ringtail.py --file_list filelist1.txt --file_path test_data/rec_gz --save_receptor --no_print"
        )

        conn = sqlite3.connect("output.db")
        cur = conn.cursor()
        cur.execute(
            "SELECT COUNT(*) FROM Receptors WHERE receptor_object NOT NULL")
        count = cur.fetchone()[0]

        cur.close()
        conn.close()

        os.system("rm output.db")

        assert count == 1

    def test_save_rec_filelist_gz(self):
        os.system(
            "python ../scripts/run_ringtail.py --file_list filelist5.txt --save_receptor --no_print"
        )

        conn = sqlite3.connect("output.db")
        cur = conn.cursor()
        cur.execute(
            "SELECT COUNT(*) FROM Receptors WHERE receptor_object NOT NULL")
        count = cur.fetchone()[0]

        cur.close()
        conn.close()

        os.system("rm output.db")

        assert count == 1


class TestOutputs:

    def test_export_table_csv(self):
        status = os.system(
            "python ../scripts/run_ringtail.py --file_list filelist1.txt --export_table_csv Ligands --no_print"
        )

        os.system("rm output.db")

        assert status == 0

    def test_export_query_csv(self):
        status = os.system(
            "python ../scripts/run_ringtail.py --file_list filelist1.txt --export_query_csv 'SELECT * FROM Results' --no_print"
        )

        os.system("rm output.db")

        assert status == 0

    def test_interaction_tolerance(self):
        status_notol = os.system(
            "python ../scripts/run_ringtail.py --file_list filelist1.txt")

        conn = sqlite3.connect("output.db")
        cur = conn.cursor()
        cur.execute(
            "SELECT * FROM Interaction_bitvectors WHERE Pose_ID in (SELECT Pose_ID FROM Results WHERE LigName LIKE 'OB3Z3759440327_RX1--4xfx_mon_prep--tyr169' AND run_number = 27)"
        )
        count_notol = sum(
            [1 for interaction in cur.fetchone() if interaction == 1])

        cur.close()
        conn.close()

        os.system("rm output.db")

        status_tol = os.system(
            "python ../scripts/run_ringtail.py --file_list filelist1.txt --interaction_tolerance"
        )

        conn = sqlite3.connect("output.db")
        cur = conn.cursor()
        cur.execute(
            "SELECT * FROM Interaction_bitvectors WHERE Pose_ID in (SELECT Pose_ID FROM Results WHERE LigName LIKE 'OB3Z3759440327_RX1--4xfx_mon_prep--tyr169' AND run_number = 27)"
        )
        count_tol = sum(
            [1 for interaction in cur.fetchone() if interaction == 1])

        cur.close()
        conn.close()

        os.system("rm output.db")

        status_tol2 = os.system(
            "python ../scripts/run_ringtail.py --file_list filelist1.txt --interaction_tolerance 2.0"
        )

        conn = sqlite3.connect("output.db")
        cur = conn.cursor()
        cur.execute(
            "SELECT * FROM Interaction_bitvectors WHERE Pose_ID in (SELECT Pose_ID FROM Results WHERE LigName LIKE 'OB3Z3759440327_RX1--4xfx_mon_prep--tyr169' AND run_number = 27)"
        )
        count_tol2 = sum(
            [1 for interaction in cur.fetchone() if interaction == 1])

        cur.close()
        conn.close()

        os.system("rm output.db")

        assert status_notol == 0
        assert status_tol == 0
        assert status_tol2 == 0
        assert count_notol != count_tol
        assert count_tol2 != count_tol
        assert count_tol2 != count_notol

    def test_max_poses(self):
        status3 = os.system(
            "python ../scripts/run_ringtail.py --file_list filelist1.txt")
        conn = sqlite3.connect("output.db")
        cur = conn.cursor()
        cur.execute("SELECT COUNT(*) FROM Results")
        count3 = cur.fetchone()[0]

        cur.close()
        conn.close()

        os.system("rm output.db")

        status1 = os.system(
            "python ../scripts/run_ringtail.py --file_list filelist1.txt --max_poses 1"
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
            "python ../scripts/run_ringtail.py --file_list filelist1.txt --max_poses 5"
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
        assert count1 < count3 < count5
        assert count1 == ligcount1

    def test_store_all(self):
        status = os.system(
            "python ../scripts/run_ringtail.py --file_list filelist1.txt --store_all_poses"
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
        assert count == ligcount * 50


class TestFilters():

    def test_eworst(self):
        status = os.system(
            "python ../scripts/run_ringtail.py --file_list filelist1.txt --eworst -15"
        )

        os.system("rm output.db")

        assert status == 0

    def test_ebest(self):
        status = os.system(
            "python ../scripts/run_ringtail.py --file_list filelist1.txt --ebest -15"
        )

        os.system("rm output.db")

        assert status == 0

    def test_leworst(self):
        status = os.system(
            "python ../scripts/run_ringtail.py --file_list filelist1.txt --leworst -0.4"
        )

        os.system("rm output.db")

        assert status == 0

    def test_lebest(self):
        status = os.system(
            "python ../scripts/run_ringtail.py --file_list filelist1.txt --leworst -0.4"
        )

        os.system("rm output.db")

        assert status == 0

    def test_epercentile(self):
        status = os.system(
            "python ../scripts/run_ringtail.py --file_list filelist1.txt --epercentile 0.1"
        )

        os.system("rm output.db")

        assert status == 0

    def test_lepercentile(self):
        status = os.system(
            "python ../scripts/run_ringtail.py --file_list filelist1.txt --leffpercentile 0.1"
        )

        os.system("rm output.db")

        assert status == 0

    def test_epercentile_eworst(self):
        status = os.system(
            "python ../scripts/run_ringtail.py --file_list filelist1.txt --epercentile 0.1 --eworst -14"
        )

        os.system("rm output.db")

        assert status == 0

    def test_lepercentile_leworst(self):
        status = os.system(
            "python ../scripts/run_ringtail.py --file_list filelist1.txt --leffpercentile 0.1 --leworst -0.4"
        )

        os.system("rm output.db")

        assert status == 0

    def test_name(self):
        status = os.system(
            "python ../scripts/run_ringtail.py --file_list filelist1.txt --name OB3Z3759440327_RX1--4xfx_mon_prep--tyr169"
        )

        os.system("rm output.db")

        assert status == 0

    def test_substruct(self):
        status = os.system(
            "python ../scripts/run_ringtail.py --file_list filelist1.txt --substructure CO"
        )

        os.system("rm output.db")

        assert status == 0

    def test_substruct_join_or(self):
        status = os.system(
            "python ../scripts/run_ringtail.py --file_list filelist1.txt --substructure CO,OC"
        )

        os.system("rm output.db")

        assert status == 0

    def test_substruct_join_and(self):
        status = os.system(
            "python ../scripts/run_ringtail.py --file_list filelist1.txt --substructure CO,C=C --substruct_join AND"
        )

        os.system("rm output.db")

        assert status == 0

    def test_hbcount(self):
        status = os.system(
            "python ../scripts/run_ringtail.py --file_list filelist1.txt --hb_count 5"
        )

        os.system("rm output.db")

        assert status == 0

    def test_react_any(self):
        status = os.system(
            "python ../scripts/run_ringtail.py --file_list filelist1.txt --react_any"
        )

        os.system("rm output.db")

        assert status == 0

    def test_react1(self):
        status = os.system(
            "python ../scripts/run_ringtail.py --file_list filelist1.txt --react_res A:TYR:169:"
        )

        os.system("rm output.db")

        assert status == 0

    def test_react2(self):
        status = os.system(
            "python ../scripts/run_ringtail.py --file_list filelist1.txt --react_res :TYR:169:"
        )

        os.system("rm output.db")

        assert status == 0

    def test_react3(self):
        status = os.system(
            "python ../scripts/run_ringtail.py --file_list filelist1.txt --react_res :TYR::"
        )

        os.system("rm output.db")

        assert status == 0

    def test_react4(self):
        status = os.system(
            "python ../scripts/run_ringtail.py --file_list filelist1.txt --react_res A:TYR::"
        )

        os.system("rm output.db")

        assert status == 0

    def test_react5(self):
        status = os.system(
            "python ../scripts/run_ringtail.py --file_list filelist1.txt --react_res A::169:"
        )

        os.system("rm output.db")

        assert status == 0

    def test_react6(self):
        status = os.system(
            "python ../scripts/run_ringtail.py --file_list filelist1.txt --react_res A:::"
        )

        os.system("rm output.db")

        assert status == 0

    def test_react7(self):
        status = os.system(
            "python ../scripts/run_ringtail.py --file_list filelist1.txt --react_res ::169:"
        )

        os.system("rm output.db")

        assert status == 0

    def test_hb1(self):
        status = os.system(
            "python ../scripts/run_ringtail.py --file_list filelist1.txt --hb A:MET:214:"
        )

        os.system("rm output.db")

        assert status == 0

    def test_hb2(self):
        status = os.system(
            "python ../scripts/run_ringtail.py --file_list filelist1.txt --hb :MET:214:"
        )

        os.system("rm output.db")

        assert status == 0

    def test_hb3(self):
        status = os.system(
            "python ../scripts/run_ringtail.py --file_list filelist1.txt --hb :MET::"
        )

        os.system("rm output.db")

        assert status == 0

    def test_hb4(self):
        status = os.system(
            "python ../scripts/run_ringtail.py --file_list filelist1.txt --hb A:MET::"
        )

        os.system("rm output.db")

        assert status == 0

    def test_hb5(self):
        status = os.system(
            "python ../scripts/run_ringtail.py --file_list filelist1.txt --hb A::214:"
        )

        os.system("rm output.db")

        assert status == 0

    def test_hb6(self):
        status = os.system(
            "python ../scripts/run_ringtail.py --file_list filelist1.txt --hb A:::"
        )

        os.system("rm output.db")

        assert status == 0

    def test_hb7(self):
        status = os.system(
            "python ../scripts/run_ringtail.py --file_list filelist1.txt --hb ::214:"
        )

        os.system("rm output.db")

        assert status == 0

    def test_vdw1(self):
        status = os.system(
            "python ../scripts/run_ringtail.py --file_list filelist1.txt --vdw A:PRO:207:"
        )

        os.system("rm output.db")

        assert status == 0

    def test_vdw2(self):
        status = os.system(
            "python ../scripts/run_ringtail.py --file_list filelist1.txt --vdw :PRO:207:"
        )

        os.system("rm output.db")

        assert status == 0

    def test_vdw3(self):
        status = os.system(
            "python ../scripts/run_ringtail.py --file_list filelist1.txt --vdw :PRO::"
        )

        os.system("rm output.db")

        assert status == 0

    def test_vdw4(self):
        status = os.system(
            "python ../scripts/run_ringtail.py --file_list filelist1.txt --vdw A:PRO::"
        )

        os.system("rm output.db")

        assert status == 0

    def test_vdw5(self):
        status = os.system(
            "python ../scripts/run_ringtail.py --file_list filelist1.txt --vdw A::207:"
        )

        os.system("rm output.db")

        assert status == 0

    def test_vdw6(self):
        status = os.system(
            "python ../scripts/run_ringtail.py --file_list filelist1.txt --vdw A:::"
        )

        os.system("rm output.db")

        assert status == 0

    def test_vdw7(self):
        status = os.system(
            "python ../scripts/run_ringtail.py --file_list filelist1.txt --vdw ::207:"
        )

        os.system("rm output.db")

        assert status == 0

    def test_all_filters(self):
        status = os.system(
            "python ../scripts/run_ringtail.py --file_list filelist1.txt --eworst -15 --ebest -16 --leworst -0.4 --lebest -0.5 --epercentile 99 --leffpercentile 99 --name OB3Z3759440327_RX1--4xfx_mon_prep--tyr169,OB3Z3759440327_RX1--4xfx_mon_prep--tyr169 --substructure CO,C=C --substruct_join AND, --hb_count 5 --react_any --hb A:MET:214: --vdw A:PRO:207: --react_res A:TYR:169:"
        )

        os.system("rm output.db")

        assert status == 0
