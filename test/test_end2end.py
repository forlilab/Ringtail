#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail end-to-end testing
#

import sqlite3
import os
import pytest

"""class TestInputs:
	def test_multiple_files1(self):
		os.system("python ../scripts/run_ringtail.py --file test_data/TEST_0000003-results/OB3Z3759450716_RX1--4xfx_mon_prep--tyr169.dlg.gz --file test_data/TEST_0000003-results/OB3Z3759440327_RX1--4xfx_mon_prep--tyr169.dlg.gz --file test_data/TEST_0000003-results/OB3Z3759457587_RX1--4xfx_mon_prep--tyr169.dlg.gz --file test_data/TEST_0000003-results/OB3Z3759444505_1_RX1--4xfx_mon_prep--tyr169.dlg.gz --file test_data/TEST_0000001-results/OB3Z3759305928_RX1--4xfx_mon_prep--tyr169.dlg.gz --no_print")

		conn = sqlite3.connect("output.db")
		cur = conn.cursor()
		cur.execute("SELECT COUNT(*) FROM Ligands")
		count = cur.fetchone()[0]

		cur.close()
		conn.close()

		os.system("rm output.db")

		assert count == 5

	def test_multiple_files2(self):
		os.system("python ../scripts/run_ringtail.py --file test_data/TEST_0000003-results/OB3Z3759450716_RX1--4xfx_mon_prep--tyr169.dlg.gz test_data/TEST_0000003-results/OB3Z3759440327_RX1--4xfx_mon_prep--tyr169.dlg.gz test_data/TEST_0000003-results/OB3Z3759457587_RX1--4xfx_mon_prep--tyr169.dlg.gz test_data/TEST_0000003-results/OB3Z3759444505_1_RX1--4xfx_mon_prep--tyr169.dlg.gz --file test_data/TEST_0000001-results/OB3Z3759305928_RX1--4xfx_mon_prep--tyr169.dlg.gz --no_print")

		conn = sqlite3.connect("output.db")
		cur = conn.cursor()
		cur.execute("SELECT COUNT(*) FROM Ligands")
		count = cur.fetchone()[0]

		cur.close()
		conn.close()

		os.system("rm output.db")

		assert count == 5

	def test_multiple_paths1(self):
		os.system("python ../scripts/run_ringtail.py --file_path test_data/TEST_0000003-results --file_path test_data/TEST_0000005-results --no_print")

		conn = sqlite3.connect("output.db")
		cur = conn.cursor()
		cur.execute("SELECT COUNT(*) FROM Ligands")
		count = cur.fetchone()[0]

		cur.close()
		conn.close()

		os.system("rm output.db")

		assert count == 19999

	def test_multiple_paths2(self):
		os.system("python ../scripts/run_ringtail.py --file_path test_data/TEST_0000003-results test_data/TEST_0000005-results --no_print")

		conn = sqlite3.connect("output.db")
		cur = conn.cursor()
		cur.execute("SELECT COUNT(*) FROM Ligands")
		count = cur.fetchone()[0]

		cur.close()
		conn.close()

		os.system("rm output.db")

		assert count == 19999

	def test_filelist1(self):
		os.system("python ../scripts/run_ringtail.py --file_list filelist1.txt --file_list filelist2.txt --no_print")

		conn = sqlite3.connect("output.db")
		cur = conn.cursor()
		cur.execute("SELECT COUNT(*) FROM Ligands")
		count = cur.fetchone()[0]

		cur.close()
		conn.close()

		os.system("rm output.db")

		assert count == 5

	def test_filelist2(self):
		os.system("python ../scripts/run_ringtail.py --file_list filelist1.txt filelist2.txt --no_print")

		conn = sqlite3.connect("output.db")
		cur = conn.cursor()
		cur.execute("SELECT COUNT(*) FROM Ligands")
		count = cur.fetchone()[0]

		cur.close()
		conn.close()

		os.system("rm output.db")

		assert count == 5

	def test_all_input_opts(self):
		os.system("mv test_data/TEST_0000003-results/OB3Z3759450716_RX1--4xfx_mon_prep--tyr169.dlg.gz test_data/OB3Z3759450716_RX1--4xfx_mon_prep--tyr169.dlg.gz")
		os.system("mv test_data/TEST_0000003-results/OB3Z3759440327_RX1--4xfx_mon_prep--tyr169.dlg.gz test_data/OB3Z3759440327_RX1--4xfx_mon_prep--tyr169.dlg.gz")
		os.system("mv test_data/TEST_0000003-results/OB3Z3759457587_RX1--4xfx_mon_prep--tyr169.dlg.gz test_data/OB3Z3759457587_RX1--4xfx_mon_prep--tyr169.dlg.gz")
		os.system("mv test_data/TEST_0000003-results/OB3Z3759444505_1_RX1--4xfx_mon_prep--tyr169.dlg.gz test_data/OB3Z3759444505_1_RX1--4xfx_mon_prep--tyr169.dlg.gz")
		os.system("mv test_data/TEST_0000001-results/OB3Z3759305928_RX1--4xfx_mon_prep--tyr169.dlg.gz test_data/OB3Z3759305928_RX1--4xfx_mon_prep--tyr169.dlg.gz")

		os.system("python ../scripts/run_ringtail.py --file_list filelist3.txt --file test_data/OB3Z3759444505_1_RX1--4xfx_mon_prep--tyr169.dlg.gz test_data/OB3Z3759305928_RX1--4xfx_mon_prep--tyr169.dlg.gz --file_path test_data/TEST_0000003-results --no_print")

		conn = sqlite3.connect("output.db")
		cur = conn.cursor()
		cur.execute("SELECT COUNT(*) FROM Ligands")
		count = cur.fetchone()[0]

		cur.close()
		conn.close()

		os.system("rm output.db")
		os.system("mv test_data/OB3Z3759450716_RX1--4xfx_mon_prep--tyr169.dlg.gz test_data/TEST_0000003-results/OB3Z3759450716_RX1--4xfx_mon_prep--tyr169.dlg.gz")
		os.system("mv test_data/OB3Z3759440327_RX1--4xfx_mon_prep--tyr169.dlg.gz test_data/TEST_0000003-results/OB3Z3759440327_RX1--4xfx_mon_prep--tyr169.dlg.gz")
		os.system("mv test_data/OB3Z3759457587_RX1--4xfx_mon_prep--tyr169.dlg.gz test_data/TEST_0000003-results/OB3Z3759457587_RX1--4xfx_mon_prep--tyr169.dlg.gz")
		os.system("mv test_data/OB3Z3759444505_1_RX1--4xfx_mon_prep--tyr169.dlg.gz test_data/TEST_0000003-results/OB3Z3759444505_1_RX1--4xfx_mon_prep--tyr169.dlg.gz")
		os.system("mv test_data/OB3Z3759305928_RX1--4xfx_mon_prep--tyr169.dlg.gz test_data/TEST_0000001-results/OB3Z3759305928_RX1--4xfx_mon_prep--tyr169.dlg.gz")

		assert count == 10001

	def test_add_results(self):
		os.system("python ../scripts/run_ringtail.py --file_path test_data/TEST_0000005-results --no_print")
		os.system("python ../scripts/run_ringtail.py --input_db output.db --file_path test_data/TEST_0000003-results --add_results --no_print")

		conn = sqlite3.connect("output.db")
		cur = conn.cursor()
		cur.execute("SELECT COUNT(*) FROM Ligands")
		count = cur.fetchone()[0]

		cur.close()
		conn.close()

		os.system("rm output.db")

		assert count == 19999

	def test_conflict_handling_ign(self):
		os.system("python ../scripts/run_ringtail.py --file_path test_data/TEST_0000005-results --no_print")
		os.system("python ../scripts/run_ringtail.py --input_db output.db --file_path test_data/TEST_0000005-results --add_results --no_print --conflict_handling ignore")

		conn = sqlite3.connect("output.db")
		cur = conn.cursor()
		cur.execute("SELECT COUNT(*) FROM Ligands")
		count = cur.fetchone()[0]

		cur.close()
		conn.close()

		os.system("rm output.db")

		assert count == 9999

	def test_conflict_handling_rpl(self):
		os.system("python ../scripts/run_ringtail.py --file_path test_data/TEST_0000005-results --no_print")
		os.system("python ../scripts/run_ringtail.py --input_db output.db --file_path test_data/TEST_0000005-results --add_results --no_print --conflict_handling replace")

		conn = sqlite3.connect("output.db")
		cur = conn.cursor()
		cur.execute("SELECT COUNT(*) FROM Ligands")
		count = cur.fetchone()[0]

		cur.close()
		conn.close()

		os.system("rm output.db")

		assert count == 9999

	def test_save_rec_file(self):
		os.system("python ../scripts/run_ringtail.py --file_list filelist1.txt --file test_data/4xfx_mon_prep--tyr169_rigid.pdbqt --save_receptor --no_print")

		conn = sqlite3.connect("output.db")
		cur = conn.cursor()
		cur.execute("SELECT COUNT(*) FROM Receptors WHERE receptor_object NOT NULL")
		count = cur.fetchone()[0]

		cur.close()
		conn.close()

		os.system("rm output.db")

		assert count == 1

	def test_save_rec_filepath(self):
		os.system("python ../scripts/run_ringtail.py --file_list filelist1.txt --file_path test_data/rec --save_receptor --no_print")

		conn = sqlite3.connect("output.db")
		cur = conn.cursor()
		cur.execute("SELECT COUNT(*) FROM Receptors WHERE receptor_object NOT NULL")
		count = cur.fetchone()[0]

		cur.close()
		conn.close()

		os.system("rm output.db")

		assert count == 1

	def test_save_rec_filelist(self):
		os.system("python ../scripts/run_ringtail.py --file_list filelist4.txt --save_receptor --no_print")

		conn = sqlite3.connect("output.db")
		cur = conn.cursor()
		cur.execute("SELECT COUNT(*) FROM Receptors WHERE receptor_object NOT NULL")
		count = cur.fetchone()[0]

		cur.close()
		conn.close()

		os.system("rm output.db")

		assert count == 1

	def test_save_rec_file_gz(self):
		os.system("python ../scripts/run_ringtail.py --file_list filelist1.txt --file test_data/4xfx_mon_prep--tyr169_rigid.pdbqt.gz --save_receptor --no_print")

		conn = sqlite3.connect("output.db")
		cur = conn.cursor()
		cur.execute("SELECT COUNT(*) FROM Receptors WHERE receptor_object NOT NULL")
		count = cur.fetchone()[0]

		cur.close()
		conn.close()

		os.system("rm output.db")

		assert count == 1

	def test_save_rec_filepath_gz(self):
		os.system("python ../scripts/run_ringtail.py --file_list filelist1.txt --file_path test_data/rec_gz --save_receptor --no_print")

		conn = sqlite3.connect("output.db")
		cur = conn.cursor()
		cur.execute("SELECT COUNT(*) FROM Receptors WHERE receptor_object NOT NULL")
		count = cur.fetchone()[0]

		cur.close()
		conn.close()

		os.system("rm output.db")

		assert count == 1

	def test_save_rec_filelist_gz(self):
		os.system("python ../scripts/run_ringtail.py --file_list filelist5.txt --save_receptor --no_print")

		conn = sqlite3.connect("output.db")
		cur = conn.cursor()
		cur.execute("SELECT COUNT(*) FROM Receptors WHERE receptor_object NOT NULL")
		count = cur.fetchone()[0]

		cur.close()
		conn.close()

		os.system("rm output.db")

		assert count == 1"""


class TestOutputs:
	def test_export_table_csv(self):
		status = os.system("python ../scripts/run_ringtail.py --file_list filelist1.txt --export_table_csv Ligands --no_print")

		os.system("rm output.db")

		assert status == 0

	def test_export_query_csv(self):
		status = os.system("python ../scripts/run_ringtail.py --file_list filelist1.txt --export_query_csv 'SELECT * FROM Results' --no_print")

		os.system("rm output.db")

		assert status == 0







