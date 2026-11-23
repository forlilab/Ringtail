#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail CLI end-to-end tests
#
import json
import shutil
import subprocess
from pathlib import Path

import pytest
from ringtail import RingtailCore
from ringtail.schema import RESULTS_SCHEMA

# all tests here are slow since each one runs a subprocess
pytestmark = pytest.mark.slow

TEST_DIR = Path(__file__).parent
TEST_DATA = TEST_DIR / "test_data"


class TestInputs:
    """Write-mode CLI input variations: files, paths, lists, docking modes."""

    def test_file_paths(self, cli):
        rc = cli.write(
            "-m",
            "adgpu",
            "--docking_results",
            str(TEST_DATA / "adgpu/group1"),
            "--docking_results",
            str(TEST_DATA / "adgpu/group2"),
        )
        assert rc == 0
        assert cli.count("Ligands") == 217

    def test_file_list(self, cli):
        rc = cli.write(
            "-m",
            "adgpu",
            "--docking_results",
            str(TEST_DATA / "adgpu/filelist1.txt"),
            "--docking_results",
            str(TEST_DATA / "adgpu/filelist2.txt"),
        )
        assert rc == 0
        assert cli.count("Ligands") == 5

    def test_all_file_inputs(self, cli):
        rc = cli.write(
            "-m",
            "adgpu",
            "--file_list",
            str(TEST_DATA / "adgpu/filelist1.txt"),
            "--file",
            str(TEST_DATA / "adgpu/group2/361056.dlg.gz"),
            str(TEST_DATA / "adgpu/group2/53506.dlg.gz"),
            "--file_path",
            str(TEST_DATA / "adgpu/group3"),
        )
        assert rc == 0
        assert cli.count("Ligands") == 75

    def test_ad6_input(self, cli):
        rc = cli.write(
            "-m",
            "ad6",
            "--docking_results",
            str(TEST_DATA / "ad6"),
            "-rf",
            str(TEST_DATA / "ad6/helix--scofu01.json"),
            "-sr",
        )
        assert rc == 0
        assert cli.count("Results") == 9
        assert cli.count("Interactions") == 53

    def test_vina_input(self, cli):
        rc = cli.write(
            "-m",
            "vina",
            "--docking_results",
            str(TEST_DATA / "vina"),
            "-rf",
            str(TEST_DATA / "vina/receptor.pdbqt"),
        )
        assert rc == 0
        assert cli.count("Results") == 6

    def test_overwrite(self, cli):
        cli.write(
            "-m", "adgpu", "--docking_results", str(TEST_DATA / "adgpu/filelist2.txt")
        )
        count_before = cli.count("Ligands")
        rc = cli.write(
            "-m",
            "adgpu",
            "--docking_results",
            str(TEST_DATA / "adgpu/filelist1.txt"),
            "--overwrite",
        )
        count_after = cli.count("Ligands")
        assert rc == 0
        assert count_before == 2
        assert count_after == 3

    def test_cmdline_config_file(self, cli, tmp_path):
        config_data = RingtailCore.defaults()
        config_data["docking_results"] = [[str(TEST_DATA / "adgpu/filelist1.txt")]]
        config_data["docking_mode"] = "adgpu"
        config_path = tmp_path / "config.json"
        config_path.write_text(json.dumps(config_data, indent=4))

        rc = cli.write("--config", str(config_path))
        assert rc == 0
        assert cli.count("Ligands") == 3

    def test_append_results(self, cli):
        cli.write("-m", "adgpu", "--docking_results", str(TEST_DATA / "adgpu/group1"))
        rc = cli.write(
            "-m",
            "adgpu",
            "--input_db",
            str(cli.db),
            "--docking_results",
            str(TEST_DATA / "adgpu/group2"),
            "--append_results",
        )
        assert rc == 0
        assert cli.count("Ligands") == 217

    def test_append_results_with_bookmarks(self, cli):
        """Appending to a db with bookmarks asks for consent, and refusing writes nothing"""
        cli.write("-m", "adgpu", "--docking_results", str(TEST_DATA / "adgpu/group1"))
        cli.read("--eworst", "-6", "-s", "hits")
        append = (
            "-m",
            "adgpu",
            "--input_db",
            str(cli.db),
            "--docking_results",
            str(TEST_DATA / "adgpu/group2"),
            "--append_results",
        )
        ligands_before = cli.count("Ligands")

        # no answer to the prompt counts as refusing
        assert cli.write(*append) == 1
        assert cli.count("Ligands") == ligands_before
        assert cli.count("hits") > 0

        assert cli.write(*append, "--yes") == 0
        assert cli.count("Ligands") == 217
        assert cli.count("Filters") == 0

    
    def test_max_poses(self, cli):
        rc = cli.write(
            "-m",
            "adgpu",
            "--docking_results",
            str(TEST_DATA / "adgpu/filelist1.txt"),
            "--max_poses",
            "1",
        )
        assert rc == 0
        assert cli.count("Results") == cli.count("Ligands")

    def test_store_all(self, cli):
        rc = cli.write(
            "-m",
            "adgpu",
            "--docking_results",
            str(TEST_DATA / "adgpu/filelist1.txt"),
            "--store_all_poses",
        )
        assert rc == 0
        assert cli.count("Results") == cli.count("Ligands") * 20


    def test_save_rec_file(self, cli):
        rc = cli.write(
            "-m",
            "adgpu",
            "--docking_results",
            str(TEST_DATA / "adgpu/group1"),
            "--receptor_file",
            str(TEST_DATA / "adgpu/4j8m.pdbqt"),
            "--save_receptor",
        )
        assert rc == 0
        assert cli.count("Receptors") == 1


class TestOutputs:
    def test_output_log_and_csv(self, cli, tmp_path):
        import csv

        cli.write(
            "-m", "adgpu", "--docking_results", str(TEST_DATA / "adgpu/filelist1.txt")
        )
        log_file = tmp_path / "filter_log.txt"
        csv_file = tmp_path / "hits.csv"
        rc = cli.read(
            "-e",
            "-6",
            "-s",
            "hits",
            "--output_log",
            str(log_file),
            "--export_bookmark_csv",
            str(csv_file),
        )
        assert rc == 0

        # log: line 32 is first data row
        with open(log_file) as f:
            lines = f.readlines()
        assert lines[31].strip() == "127458, -6.66"

        # without --outfields the csv has every Results column
        with open(csv_file) as f:
            reader = csv.reader(f)
            headers = next(reader)
            rows = list(reader)
        assert set(headers) == set(RESULTS_SCHEMA.columns)
        score_idx = headers.index("docking_score")
        assert len(rows) == 1
        assert float(rows[0][score_idx]) == pytest.approx(-6.66, abs=0.01)

        # with --outfields the csv only has those columns
        csv_file = tmp_path / "hits_outfields.csv"
        rc = cli.read(
            "-e",
            "-6",
            "-s",
            "hits_outfields",
            "--outfields",
            "ligname,docking_score",
            "--export_bookmark_csv",
            str(csv_file),
        )
        assert rc == 0
        with open(csv_file) as f:
            reader = csv.reader(f)
            headers = next(reader)
            rows = list(reader)
        assert set(headers) - {""} == {"ligname", "docking_score"}
        ligname_idx = headers.index("ligname")
        score_idx = headers.index("docking_score")
        assert len(rows) == 1
        assert rows[0][ligname_idx] == "127458"
        assert float(rows[0][score_idx]) == pytest.approx(-6.66, abs=0.01)

    def test_export_table_csv(self, cli, tmp_path):
        cli.write(
            "-m", "adgpu", "--docking_results", str(TEST_DATA / "adgpu/filelist1.txt")
        )
        rc = cli.read("-s", "Ligands", "-xs")
        assert rc == 0
        assert Path(tmp_path / "Ligands.csv").exists()

@pytest.mark.parametrize("storage_type", ["duckdb"], indirect=True)
class TestFilters:
    """Filter options on a copy of populated_db, same data as the filter tests in test_units.py"""

    @pytest.fixture(autouse=True)
    def setup_db(self, cli, populated_db):
        shutil.copy2(populated_db.db_file, cli.db)
        self.cli = cli

    def test_score_percentile(self, ad6_db):
        shutil.copy2(ad6_db.db_file, self.cli.db)
        assert self.cli.read("--score_percentile", "50", "-s", "percentile") == 0
        assert self.cli.count("percentile", "ligand_id") == 2

    def test_ligand_filters(self):
        rc = self.cli.read(
            "--ligand_substruct",
            "C=O",
            "CC(C)(C)",
            "--ligand_operator",
            "AND",
            "-s",
            "substruct_and",
        )
        assert rc == 0
        assert self.cli.count("substruct_and", "ligand_id") == 18

    @pytest.mark.parametrize("hb_spec", ["A:LYS:162:", ":LYS::", "A:::"])
    def test_hb_interaction_formats(self, hb_spec):
        assert self.cli.read("-hb", hb_spec) == 0

    def test_all_filters(self):
        """CLI accepts every filter type at once"""
        rc = self.cli.read(
            "--eworst",
            "-15",
            "--ebest",
            "-16",
            "--leworst",
            "-0.4",
            "--lebest",
            "-0.5",
            "--ligand_name",
            "127458",
            "--hb_count",
            "5",
            "--react_any",
            "-hb",
            "A:LYS:162:",
            "-vdw",
            "A:VAL:279:",
        )
        assert rc == 0

    def test_export_sdf(self, tmp_path):
        rc = self.cli.read("-s", "sdf_test", "-e", "-4", "-sdf", str(tmp_path))
        sdf_files = list(tmp_path.glob("*.sdf"))
        assert len(sdf_files) == 1
        assert rc == 0

    def test_clustering(self):
        cli = self.cli
        cli.read("--eworst", "-6")
        count1 = cli.count("filtered_poses")
        rc = cli.read(
            "--input_bookmark",
            "passing_results",
            "--bookmark_name",
            "clustered_results",
            "-mfpc",
            "0.9",
        )
        count2 = cli.count("filtered_poses")
        assert rc == 0
        assert count1 == 68
        assert count2 == 68 + 5


class TestOtherScripts:
    def test_rt_compare_same_file_names(self, tmp_path):
        """rt_compare works on dbs with the same file name in different folders"""
        for target in ["target_a", "target_b", "target_c"]:
            (tmp_path / target).mkdir()
            db = str(tmp_path / target / "output.db")
            subprocess.run(
                [
                    "rt_process_vs",
                    "write",
                    "-o",
                    db,
                    "-st",
                    "duckdb",
                    "-m",
                    "adgpu",
                    "--docking_results",
                    str(TEST_DATA / "adgpu/filelist1.txt"),
                ],
                cwd=str(TEST_DIR),
                capture_output=True,
            )
            subprocess.run(
                ["rt_process_vs", "read", "--input_db", db, "--eworst", "-6"],
                cwd=str(tmp_path),
                capture_output=True,
            )
        rc = subprocess.run(
            [
                "rt_compare",
                "-w",
                "target_a/output.db",
                "passing_results",
                "-w",
                "target_b/output.db",
                "passing_results",
                "-uw",
                "target_c/output.db",
                "passing_results",
            ],
            cwd=str(tmp_path),
            capture_output=True,
        ).returncode
        assert rc == 0
