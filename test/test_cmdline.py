#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail command line tool end-to-end testing
#

import sqlite3
import os
import pytest
from ringtail import RingtailCore


@pytest.fixture
def passingcount():
    def __dbconnect(bookmark):
        rtc = RingtailCore("output.db")
        with rtc.storageman as sm:
            count = sm.get_passing_poses_count(bookmark, True)
        return count

    return __dbconnect


@pytest.fixture
def tablecount():
    def __dbconnect(table):
        rtc = RingtailCore("output.db")
        with rtc.storageman as sm:
            count = sm.db_query(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
        return count

    return __dbconnect


class TestInputs:
    os.system("rm output.db")

    def test_files(self, tablecount):
        os.system(
            "python ../ringtail/cli/rt_process_vs.py write --file test_data/adgpu/group1/127458.dlg.gz --file test_data/adgpu/group1/173101.dlg.gz --file test_data/adgpu/group1/100729.dlg.gz"
        )
        count1 = tablecount("Ligands")

        os.system(
            "python ../ringtail/cli/rt_process_vs.py write --file test_data/adgpu/group1/127458.dlg.gz test_data/adgpu/group1/173101.dlg.gz --file test_data/adgpu/group1/100729.dlg.gz --append_results"
        )
        count2 = tablecount("Ligands")

        os.system("rm output.db")

        assert count1 == count2 == 3

    def test_file_paths(self, tablecount):
        os.system(
            "python ../ringtail/cli/rt_process_vs.py write --file_path test_data/adgpu/group1 --file_path test_data/adgpu/group2"
        )
        count1 = tablecount("Ligands")

        os.system("rm output.db")

        os.system(
            "python ../ringtail/cli/rt_process_vs.py write --file_path test_data/adgpu/group1 test_data/adgpu/group2"
        )
        count2 = tablecount("Ligands")

        os.system("rm output.db")

        assert count1 == count2 == 217

    def test_file_list(self, tablecount):
        os.system(
            "python ../ringtail/cli/rt_process_vs.py write --file_list test_data/filelist1.txt --file_list test_data/filelist2.txt"
        )
        count1 = tablecount("Ligands")

        os.system("rm output.db")

        os.system(
            "python ../ringtail/cli/rt_process_vs.py write --file_list test_data/filelist1.txt test_data/filelist2.txt"
        )
        count2 = tablecount("Ligands")

        os.system("rm output.db")

        assert count1 == count2 == 5

    def test_all_file_inputs(self, tablecount):
        os.system(
            "python ../ringtail/cli/rt_process_vs.py write --file_list test_data/filelist1.txt --file test_data/adgpu/group2/361056.dlg.gz test_data/adgpu/group2/53506.dlg.gz --file_path test_data/adgpu/group3"
        )
        count = tablecount("Ligands")

        os.system("rm output.db")

        assert count == 75

    def test_vina_input(self, tablecount):
        os.system(
            "python ../ringtail/cli/rt_process_vs.py write -m vina --file_path test_data/vina -rf test_data/vina/receptor.pdbqt -sr"
        )

        count = tablecount("Results")
        assert count == 6

    def test_overwrite(self, tablecount):
        # count result rows in database to be overwritten

        count_old_db = tablecount("Ligands")
        os.system(
            "python ../ringtail/cli/rt_process_vs.py write --file_list test_data/filelist1.txt --overwrite"
        )
        count_new_db = tablecount("Ligands")
        assert count_old_db == 2
        assert count_new_db == 3

    def test_overwrite_false(self, tablecount):
        # count result rows in database to be overwritten

        count_old_db = tablecount("Ligands")
        assert count_old_db == 3

        code = os.system(
            "python ../ringtail/cli/rt_process_vs.py write --file_list test_data/filelist1.txt"
        )
        assert (
            code == 256
        )  # indicates failure of rt_process_vs.py, log file will have error w traceback

    def test_cmdline_config_file(self, tablecount):
        from ringtail import RingtailCore
        import json

        filepath = RingtailCore.generate_config_file_template()

        with open(filepath, "r") as f:
            data = json.load(f)
        # all fields to be changed
        data["file_list"] = [["test_data/filelist1.txt"]]

        with open(filepath, "w") as f:
            f.write(json.dumps(data, indent=4))

        os.system("python ../ringtail/cli/rt_process_vs.py write --config config.json")

        count = tablecount("Ligands")

        os.system("rm output.db config.json")

        assert count == 3

    def test_duplicate_handling(self, tablecount):
        os.system(
            "python ../ringtail/cli/rt_process_vs.py write --file_path test_data/adgpu/group1"
        )
        os.system(
            "python ../ringtail/cli/rt_process_vs.py write --input_db output.db --file_path test_data/adgpu/group1 --append_results --duplicate_handling ignore"
        )
        count = tablecount("Ligands")
        assert count == 138

    def test_append_results(self, tablecount):
        os.system(
            "python ../ringtail/cli/rt_process_vs.py write --input_db output.db --file_path test_data/adgpu/group2 --append_results"
        )
        count = tablecount("Ligands")

        assert count == 217

    def test_save_rec_file(self, tablecount):

        os.system(
            "python ../ringtail/cli/rt_process_vs.py write --input_db output.db --receptor_file test_data/adgpu/4j8m.pdbqt --save_receptor --append_results"
        )
        count = tablecount("Receptors")

        os.system("rm output.db")

        assert count == 1

    def test_save_rec_file_gz(self, tablecount):
        os.system(
            "python ../ringtail/cli/rt_process_vs.py write --file_list test_data/filelist1.txt --receptor_file test_data/adgpu/4j8m.pdbqt.gz --save_receptor"
        )
        count = tablecount("Receptors")

        os.system("rm output.db")

        assert count == 1


class TestOutputs:
    def test_export_bookmark_csv(self):
        status1 = os.system(
            "python ../ringtail/cli/rt_process_vs.py write --file_list test_data/filelist1.txt"
        )
        status2 = os.system(
            "python ../ringtail/cli/rt_process_vs.py read --input_db output.db --export_bookmark_csv Ligands"
        )

        assert status1 == status2 == 0
        assert os.path.exists("Ligands.csv")

        os.system("rm Ligands.csv")

    def test_export_query_csv(self):

        status = os.system(
            "python ../ringtail/cli/rt_process_vs.py read --input_db output.db --export_query_csv 'SELECT * FROM Results'"
        )

        assert status == 0
        assert os.path.exists("query.csv")

        os.system("rm output.db")
        os.system("rm query.csv")

    def test_interaction_tolerance(self):
        status_notol = os.system(
            "python ../ringtail/cli/rt_process_vs.py write --file test_data/adgpu/group1/127458.dlg.gz"
        )
        query_to_check = """
        SELECT * FROM Interactions 
        WHERE Pose_ID in 
        (SELECT Pose_ID FROM Results 
            WHERE ligand_id = 
                (SELECT ligand_id from Ligands WHERE LigName = '127458')
            AND run_number = 13)"""

        conn = sqlite3.connect("output.db")
        cur = conn.cursor()

        cur.execute(query_to_check)
        count_notol = len(cur.fetchall())

        cur.close()
        conn.close()

        os.system("rm output.db")

        status_tol = os.system(
            "python ../ringtail/cli/rt_process_vs.py write --file test_data/adgpu/group1/127458.dlg.gz --interaction_tolerance"
        )

        conn = sqlite3.connect("output.db")
        cur = conn.cursor()
        cur.execute(query_to_check)
        count_tol = len(cur.fetchall())

        cur.close()
        conn.close()

        os.system("rm output.db")

        status_tol2 = os.system(
            "python ../ringtail/cli/rt_process_vs.py write --file test_data/adgpu/group1/127458.dlg.gz --interaction_tolerance 2.0"
        )

        conn = sqlite3.connect("output.db")
        cur = conn.cursor()
        cur.execute(query_to_check)
        count_tol2 = len(cur.fetchall())

        cur.close()
        conn.close()

        assert status_notol == 0
        assert status_tol == 0
        assert status_tol2 == 0
        assert (
            count_notol != count_tol
            or count_tol2 != count_tol
            or count_tol2 != count_notol
        )

    def test_max_poses(self):
        os.system("rm output.db")
        status3 = os.system(
            "python ../ringtail/cli/rt_process_vs.py write --file_list test_data/filelist1.txt"
        )
        conn = sqlite3.connect("output.db")
        cur = conn.cursor()
        cur.execute("SELECT COUNT(*) FROM Results")
        count3 = cur.fetchone()[0]

        cur.close()
        conn.close()

        os.system("rm output.db")

        status1 = os.system(
            "python ../ringtail/cli/rt_process_vs.py write --file_list test_data/filelist1.txt --max_poses 1"
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
            "python ../ringtail/cli/rt_process_vs.py write --file_list test_data/filelist1.txt --max_poses 5"
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
            "python ../ringtail/cli/rt_process_vs.py write --file_list test_data/filelist1.txt --store_all_poses"
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

    def test_eworst(self):
        status1 = os.system(
            "python ../ringtail/cli/rt_process_vs.py write --file_list test_data/filelist1.txt"
        )
        status2 = os.system(
            "python ../ringtail/cli/rt_process_vs.py read --input_db output.db --eworst -15"
        )

        assert status1 == status2 == 0

    def test_ebest(self):
        status = os.system(
            "python ../ringtail/cli/rt_process_vs.py read --input_db output.db --ebest -15"
        )

        assert status == 0

    def test_leworst(self):
        status = os.system(
            "python ../ringtail/cli/rt_process_vs.py read --input_db output.db --leworst -0.4"
        )

        assert status == 0

    def test_lebest(self):
        status = os.system(
            "python ../ringtail/cli/rt_process_vs.py read --input_db output.db --leworst -0.4"
        )

        assert status == 0

    def test_epercentile(self):
        status = os.system(
            "python ../ringtail/cli/rt_process_vs.py read --input_db output.db --score_percentile 0.1"
        )

        assert status == 0

    def test_lepercentile(self):
        status = os.system(
            "python ../ringtail/cli/rt_process_vs.py read --input_db output.db --le_percentile 0.1"
        )

        assert status == 0

    def test_epercentile_eworst(self):
        status = os.system(
            "python ../ringtail/cli/rt_process_vs.py read --input_db output.db --score_percentile 0.1 --eworst -14"
        )

        assert status == 0

    def test_lepercentile_leworst(self):
        status = os.system(
            "python ../ringtail/cli/rt_process_vs.py read --input_db output.db --le_percentile 0.1 --leworst -0.4"
        )

        assert status == 0

    def test_name(self):
        status = os.system(
            "python ../ringtail/cli/rt_process_vs.py read --input output.db --ligand_name 127458"
        )

        assert status == 0

    def test_ligand_filters(self):
        status = os.system(
            """python ../ringtail/cli/rt_process_vs.py read --input output.db --ligand_substruct "NC" --ligand_operator AND"""
        )

        assert status == 0

    def test_hbcount(self, passingcount):
        status = os.system(
            "python ../ringtail/cli/rt_process_vs.py read -s hb_count --input_db output.db --hb_count 5"
        )
        count = passingcount("hb_count")

        assert status == 0
        assert count == 1

    def test_hb1(self):
        status = os.system(
            "python ../ringtail/cli/rt_process_vs.py read --input_db output.db -hb A:LYS:162:"
        )

        assert status == 0

    def test_hb2(self):
        status = os.system(
            "python ../ringtail/cli/rt_process_vs.py read --input_db output.db -hb :LYS:162:"
        )

        assert status == 0

    def test_hb3(self):
        status = os.system(
            "python ../ringtail/cli/rt_process_vs.py read --input_db output.db -hb :LYS::"
        )

        assert status == 0

    def test_hb4(self):
        status = os.system(
            "python ../ringtail/cli/rt_process_vs.py read --input_db output.db -hb A:LYS::"
        )

        assert status == 0

    def test_hb5(self):
        status = os.system(
            "python ../ringtail/cli/rt_process_vs.py read --input_db output.db -hb A::162:"
        )

        assert status == 0

    def test_hb6(self):
        status = os.system(
            "python ../ringtail/cli/rt_process_vs.py read --input_db output.db -hb A:::"
        )

        assert status == 0

    def test_hb7(self):
        status = os.system(
            "python ../ringtail/cli/rt_process_vs.py read --input_db output.db -hb ::162:"
        )

        assert status == 0

    def test_vdw1(self):
        status = os.system(
            "python ../ringtail/cli/rt_process_vs.py read --input_db output.db -vdw A:VAL:279:"
        )

        assert status == 0

    def test_vdw2(self):
        status = os.system(
            "python ../ringtail/cli/rt_process_vs.py read --input_db output.db -vdw :VAL:279:"
        )

        assert status == 0

    def test_vdw3(self):
        status = os.system(
            "python ../ringtail/cli/rt_process_vs.py read --input_db output.db -vdw :VAL::"
        )

        assert status == 0

    def test_vdw4(self):
        status = os.system(
            "python ../ringtail/cli/rt_process_vs.py read --input_db output.db -vdw A:VAL::"
        )

        assert status == 0

    def test_vdw5(self):
        status = os.system(
            "python ../ringtail/cli/rt_process_vs.py read --input_db output.db -vdw A::279:"
        )

        assert status == 0

    def test_vdw6(self):
        status = os.system(
            "python ../ringtail/cli/rt_process_vs.py read --input_db output.db -vdw A:::"
        )

        assert status == 0

    def test_vdw7(self):
        status = os.system(
            "python ../ringtail/cli/rt_process_vs.py read --input_db output.db -vdw ::279:"
        )

        assert status == 0

    def test_all_filters(self):
        status = os.system(
            "python ../ringtail/cli/rt_process_vs.py read --input_db output.db --eworst -15 --ebest -16 --leworst -0.4 --lebest -0.5 --score_percentile 99 --le_percentile 99 --ligand_name 127458 --hb_count 5 --react_any -hb A:LYS:162: -vdw A:VAL:279:"
        )

        assert status == 0

    def test_export_sdf(self):
        status = os.system(
            "python ../ringtail/cli/rt_process_vs.py read --input_db output.db -s sdf_test -e -4 -sdf . "
        )

        import glob

        sdf_files = glob.glob("*.sdf")

        assert len(sdf_files) == 1
        os.remove(sdf_files[0])

        assert status == 0

    def test_clustering(self, tablecount):
        os.system("rm output.db")
        os.system(
            "python ../ringtail/cli/rt_process_vs.py write --file_path test_data/adgpu/group1 --file_path test_data/adgpu/group2"
        )
        status = os.system(
            "python ../ringtail/cli/rt_process_vs.py read --input_db output.db --eworst -6"
        )
        count1 = tablecount("filtered_poses")
        status = os.system(
            "python ../ringtail/cli/rt_process_vs.py read --input_db output.db --bookmark_name passing_results -mfpc 0.9"
        )
        count2 = tablecount("filtered_poses")

        assert status == 0
        assert count1 == 68
        assert count2 == 68 + 5

    def test_filters_value_error(self):

        status = os.system(
            "python ../ringtail/cli/rt_process_vs.py read --input_db output.db --score_percentile 109"
        )
        # checking that code exited with error since a percentile cannot be above 100
        assert status != 0
        os.system("rm output_log.txt output.db")

    def test_react_any(self):
        # write new db with reactive data
        os.system(
            "python ../ringtail/cli/rt_process_vs.py write --output_db output.db --file_path test_data/reactive --receptor_file test_data/reactive/4j8m_m_rigid.pdbqt"
        )

        status = os.system(
            "python ../ringtail/cli/rt_process_vs.py read --input_db output.db --react_any"
        )

        assert status == 0

    def test_react1(self):
        status = os.system(
            "python ../ringtail/cli/rt_process_vs.py read --input_db output.db  --reactive_interactions A:TYR:212:"
        )

        assert status == 0

    def test_react2(self):
        status = os.system(
            "python ../ringtail/cli/rt_process_vs.py read --input_db output.db  --reactive_interactions :TYR:212:"
        )

        assert status == 0

    def test_react3(self):
        status = os.system(
            "python ../ringtail/cli/rt_process_vs.py read --input_db output.db --reactive_interactions :TYR::"
        )

        assert status == 0

    def test_react4(self):
        status = os.system(
            "python ../ringtail/cli/rt_process_vs.py read --input_db output.db --reactive_interactions A:TYR::"
        )

        assert status == 0

    def test_react5(self):
        status = os.system(
            "python ../ringtail/cli/rt_process_vs.py read --input_db output.db --reactive_interactions A::212:"
        )

        assert status == 0

    def test_react6(self):
        status = os.system(
            "python ../ringtail/cli/rt_process_vs.py read --input_db output.db --reactive_interactions A:::"
        )

        assert status == 0

    def test_react7(self):
        status = os.system(
            "python ../ringtail/cli/rt_process_vs.py read --input_db output.db --reactive_interactions ::212:"
        )

        assert status == 0
        os.system("rm output_log.txt output.db")


class TestOtherScripts:
    # TODO failing
    def test_rt_compare(self):
        # first database
        os.system(
            "python ../ringtail/cli/rt_process_vs.py write --file_path test_data/adgpu/group1"
        )
        # second database
        os.system(
            "python ../ringtail/cli/rt_process_vs.py write --output_db output2.db --file_path test_data/adgpu/group1"
        )
        # filter producing 30 ligands
        os.system(
            "python ../ringtail/cli/rt_process_vs.py read --input_db output.db --eworst -6"
        )
        # filter producing 5 ligands
        os.system(
            "python ../ringtail/cli/rt_process_vs.py read --input_db output2.db --eworst -7"
        )
        # should produce 25 ligands
        os.system(
            "python ../ringtail/cli/rt_compare.py --wanted output.db --unwanted output2.db --log compared_ligands.txt"
        )
        with open("compared_ligands.txt") as f:
            for pos, line in enumerate(f):
                if pos + 1 == 4:  # zero based line indexing
                    assert line == "Number passing ligands: 25 \n"
                    break

        os.system("rm output.db output2.db compared_ligands.txt output_log.txt")
