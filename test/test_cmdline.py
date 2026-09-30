#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail CLI end-to-end tests
#
import json
import subprocess
from pathlib import Path

import pytest
from ringtail import RingtailCore
from ringtail.schema import RESULTS_SCHEMA

# This whole module is skipped for a quick run: `pytest -m "not slow"`, because each test
# costs a subprocess
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
        # same dataset and count as test_units.py::TestAD6Handling::test_file_add
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

        rc = subprocess.run(
            ["rt_process_vs", "write", "-o", str(cli.db), "--config", str(config_path)],
            cwd=str(TEST_DIR),
            capture_output=True,
        ).returncode
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
        """Appending to a filtered database asks first; refusing writes nothing."""
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

        # the prompt gets no answer (stdin is closed), which counts as refusing
        assert cli.write(*append) == 1
        assert cli.count("Ligands") == ligands_before
        assert cli.passing("hits") > 0

        assert cli.write(*append, "--yes") == 0
        assert cli.count("Ligands") == 217
        assert cli.count("Filters") == 0

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

        # Without --outfields, CSV export keeps every Results column.
        with open(csv_file) as f:
            reader = csv.reader(f)
            headers = next(reader)
            rows = list(reader)
        assert set(headers) == set(RESULTS_SCHEMA.columns)
        score_idx = headers.index("docking_score")
        assert len(rows) == 1
        assert float(rows[0][score_idx]) == pytest.approx(-6.66, abs=0.01)

        invalid_csv = tmp_path / "invalid.csv"
        rc = cli.read(
            "--bookmark_name",
            "does_not_exist",
            "--outfields",
            "ligname,docking_score",
            "--export_bookmark_csv",
            str(invalid_csv),
        )
        assert rc != 0
        assert not invalid_csv.exists()

    def test_output_log_and_csv_with_outfields(self, cli, tmp_path):
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
            "--outfields",
            "ligname,docking_score",
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

        # With --outfields, CSV export has only the requested columns
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

    def test_max_poses(self, cli, tmp_path):
        cli.write(
            "-m", "adgpu", "--docking_results", str(TEST_DATA / "adgpu/filelist1.txt")
        )
        count_default = cli.count("Results")

        db_max1 = str(tmp_path / "max1.db")
        subprocess.run(
            [
                "rt_process_vs",
                "write",
                "-o",
                db_max1,
                "-st",
                "duckdb",
                "-m",
                "adgpu",
                "--docking_results",
                str(TEST_DATA / "adgpu/filelist1.txt"),
                "--max_poses",
                "1",
            ],
            cwd=str(TEST_DIR),
            capture_output=True,
        )
        assert RingtailCore(db_max1).table_length("Results") == RingtailCore(
            db_max1
        ).table_length("Ligands")
        assert RingtailCore(db_max1).table_length("Results") < count_default

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


class TestFilters:
    """
    Filter tests: use adgpu filelist1.txt setup per test via autouse fixture.
    Storage type is parametrized from conftest to cover both backends.
    """

    @pytest.fixture(autouse=True)
    def setup_db(self, cli):
        cli.write(
            "-m", "adgpu", "--docking_results", str(TEST_DATA / "adgpu/filelist1.txt")
        )
        self.cli = cli

    def test_epercentile(self):
        assert self.cli.read("--score_percentile", "0.1") == 0

    def test_ligand_filters(self):
        assert (
            self.cli.read("--ligand_substruct", "NC", "--ligand_operator", "AND") == 0
        )

    @pytest.mark.parametrize("hb_spec", ["A:LYS:162:", ":LYS::", "A:::"])
    def test_hb_interaction_formats(self, hb_spec):
        assert self.cli.read("-hb", hb_spec) == 0

    def test_all_filters(self):
        """Smoke test that the CLI accepts a filter of every kind at once."""
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

    def test_export_sdf(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        rc = self.cli.read("-s", "sdf_test", "-e", "-4", "-sdf", str(tmp_path))
        sdf_files = list(tmp_path.glob("*.sdf"))
        assert len(sdf_files) == 1
        assert rc == 0

    def test_clustering(self, cli):
        cli.write(
            "-m",
            "adgpu",
            "--overwrite",
            "--docking_results",
            str(TEST_DATA / "adgpu/group1"),
            "--docking_results",
            str(TEST_DATA / "adgpu/group2"),
        )
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

    def test_react_any(self, cli):
        cli.write(
            "-m",
            "adgpu",
            "--overwrite",
            "--docking_results",
            str(TEST_DATA / "reactive"),
            "--receptor_file",
            str(TEST_DATA / "reactive/4j8m_m_rigid.pdbqt"),
        )
        assert cli.read("--react_any") == 0


class TestOtherScripts:
    def test_rt_compare(self, tmp_path):
        db1 = str(tmp_path / "output.db")
        db2 = str(tmp_path / "output2.db")
        log = "compared_ligands.txt"

        subprocess.run(
            [
                "rt_process_vs",
                "write",
                "-o",
                db1,
                "-st",
                "duckdb",
                "-m",
                "adgpu",
                "--docking_results",
                str(TEST_DATA / "adgpu/group1"),
            ],
            cwd=str(TEST_DIR),
            capture_output=True,
        )
        subprocess.run(
            [
                "rt_process_vs",
                "write",
                "-o",
                db2,
                "-st",
                "duckdb",
                "-m",
                "adgpu",
                "--docking_results",
                str(TEST_DATA / "adgpu/group1"),
            ],
            cwd=str(TEST_DIR),
            capture_output=True,
        )
        subprocess.run(
            ["rt_process_vs", "read", "--input_db", db1, "--eworst", "-6"],
            cwd=str(TEST_DIR),
            capture_output=True,
        )
        subprocess.run(
            ["rt_process_vs", "read", "--input_db", db2, "--eworst", "-7"],
            cwd=str(TEST_DIR),
            capture_output=True,
        )
        subprocess.run(
            [
                "rt_compare",
                "--wanted",
                db1,
                "passing_results",
                "--unwanted",
                db2,
                "passing_results",
                "--output_log",
                log,
            ],
            cwd=str(TEST_DIR),
            capture_output=True,
        )

        output_log = tmp_path / f"output_{Path(log).name}"
        with open(output_log) as f:
            for pos, line in enumerate(f):
                if pos + 1 == 94:
                    assert line == "Number passing ligands: 24 \n"
                    break

    @pytest.mark.parametrize("storage_type", ["duckdb"])
    def test_rt_compare_same_file_names(self, tmp_path, storage_type):
        """Screens are often each named output.db, in separate folders."""
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
                    storage_type,
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
