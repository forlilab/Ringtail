#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail unit tests
#
import json
import logging
from pathlib import Path

import pytest
from ringtail import Filters, QueryBuilder, RingtailCore

TEST_DATA = Path(__file__).parent / "test_data"


class TestCoreOperations:
    """Basic write and read operations using adgpu (default mode)."""

    def test_add_file(self, tmp_db):
        tmp_db.add_results_from_files(
            docking_results=str(TEST_DATA / "adgpu/group1/1451.dlg.gz"),
            max_poses=3,
            docking_mode="adgpu",
        )
        assert tmp_db.table_length("Results") == 3

    def test_store_all_poses(self, tmp_db):
        tmp_db.add_results_from_files(
            docking_results=str(TEST_DATA / "adgpu/group1/1451.dlg.gz"),
            store_all_poses=True,
            docking_mode="adgpu",
        )
        assert tmp_db.table_length("Results") == 20

    def test_add_folder(self, tmp_db):
        tmp_db.add_results_from_files(
            docking_results=str(TEST_DATA / "adgpu/group1"), docking_mode="adgpu"
        )
        assert tmp_db.table_length("Ligands") == 138

    def test_append_to_database(self, tmp_db):
        tmp_db.add_results_from_files(
            docking_results=str(TEST_DATA / "adgpu/group1"), docking_mode="adgpu"
        )
        tmp_db.add_results_from_files(
            docking_results=str(TEST_DATA / "adgpu/group2"), docking_mode="adgpu"
        )
        assert tmp_db.table_length("Ligands") == 217

    def test_save_receptor(self, tmp_db):
        tmp_db.add_results_from_files(
            docking_results=str(TEST_DATA / "adgpu/group1"), docking_mode="adgpu"
        )
        count_before = tmp_db.db_query(
            "SELECT COUNT(*) FROM Receptors WHERE receptor_object NOT NULL"
        )[0][0]
        assert count_before == 0

        tmp_db.save_receptor(receptor=str(TEST_DATA / "adgpu/4j8m.pdbqt"))
        count_after = tmp_db.db_query(
            "SELECT COUNT(*) FROM Receptors WHERE receptor_object NOT NULL"
        )[0][0]
        assert count_after == 1

    def test_db_summary_data(self, tmp_db):
        from ringtail import exceptions as e

        with pytest.raises(e.StorageError):
            tmp_db.db_summary_data()

        tmp_db.add_results_from_files(
            docking_results=str(TEST_DATA / "adgpu/group1"), docking_mode="adgpu"
        )
        tmp_db.add_results_from_files(
            docking_results=str(TEST_DATA / "adgpu/group2"), docking_mode="adgpu"
        )
        data, _ = tmp_db.db_summary_data()
        assert isinstance(data, dict)
        assert len(data.keys()) == 12
        assert all(isinstance(v, (float, int)) for v in data.values())

    def test_duplicate_handling(self, tmp_db):
        f = str(TEST_DATA / "adgpu/group1/1451.dlg.gz")
        tmp_db.add_results_from_files(docking_results=f, docking_mode="adgpu")
        result_count = tmp_db.table_length("Results")
        inter_count = tmp_db.table_length("Interactions")

        tmp_db.add_results_from_files(
            docking_results=f, docking_mode="adgpu", duplicate_handling="replace"
        )
        assert tmp_db.table_length("Results") == result_count
        assert tmp_db.table_length("Interactions") == inter_count

        tmp_db.add_results_from_files(
            docking_results=f, docking_mode="adgpu", duplicate_handling="ignore"
        )
        assert tmp_db.table_length("Results") == result_count
        assert tmp_db.table_length("Interactions") == inter_count

        tmp_db.add_results_from_files(docking_results=f, docking_mode="adgpu")
        assert tmp_db.table_length("Results") == result_count * 2
        assert tmp_db.table_length("Interactions") == inter_count * 2

    def test_db_num_poses_warning(self, tmp_db, tmp_path):
        from ringtail import setup_logging

        logfile = str(tmp_path / "poses_warning.log")
        setup_logging(level="DEBUG", logfile=logfile)

        tmp_db.add_results_from_files(
            docking_results=str(TEST_DATA / "adgpu/group1/1451.dlg.gz"),
            max_poses=1,
            docking_mode="adgpu",
        )
        tmp_db.add_results_from_files(
            docking_results=str(TEST_DATA / "adgpu/group1/1620.dlg.gz"),
            max_poses=4,
            docking_mode="adgpu",
        )
        warning_string = "The following database properties do not agree with the properties last used for this database: \nCurrent number of poses saved is 4 but database was previously set to 1."
        with open(logfile) as f:
            assert warning_string in f.read()


class TestFiltering:
    """Filter operations on the full 217-ligand adgpu dataset."""

    def test_filter(self, populated_db):
        count, _ = populated_db.filter(
            eworst=-6,
            hb_interactions=[("A:VAL:279:", True), ("A:LYS:162:", True)],
            vdw_interactions=[("A:VAL:279:", True), ("A:LYS:162:", True)],
            max_miss=1,
            output_bookmark="union_bookmark",
        )
        assert count == 33
        bookmarks = populated_db.get_bookmark_names()
        assert len(bookmarks) == 1
        assert bookmarks[0] == "union_bookmark"

    def test_return_iter(self, populated_db):
        iterable = populated_db.filter(
            eworst=-7, output_bookmark="iterable", return_iter=True
        )
        assert len(iterable) == 8

    def test_enumerate_interaction_combinations(self, populated_db):
        bookmarks_before = populated_db.get_bookmark_names()
        count, _ = populated_db.filter(
            eworst=-6,
            hb_interactions=[("A:VAL:279:", True), ("A:LYS:162:", True)],
            vdw_interactions=[("A:VAL:279:", True), ("A:LYS:162:", True)],
            max_miss=1,
            enumerate_interaction_combs=True,
            output_bookmark="enumerated_bookmark",
        )
        assert count == 33
        new_bookmarks = [
            b for b in populated_db.get_bookmark_names() if b not in bookmarks_before
        ]
        assert len(new_bookmarks) == 6
        assert "enumerated_bookmark_0" in new_bookmarks
        assert "enumerated_bookmark_union" in new_bookmarks

    def test_filter_from_bookmark(self, populated_db):
        count1, _ = populated_db.filter(eworst=-6, output_bookmark="filter_window")
        count2, _ = populated_db.filter(
            eworst=-7, output_bookmark="bookmark", input_bookmark="filter_window"
        )
        assert count1 > count2

    def test_ligand_filters(self, populated_db):
        count_name, _ = populated_db.filter(
            ligand_name=["88"], output_bookmark="ligname"
        )
        assert count_name == 7

        count_or, _ = populated_db.filter(
            ligand_substruct=["C=O", "CC(C)(C)"], output_bookmark="substruct_or"
        )
        assert count_or == 90

        count_and, _ = populated_db.filter(
            ligand_substruct=["C=O", "CC(C)(C)"],
            ligand_operator="AND",
            output_bookmark="substruct_and",
        )
        assert count_and == 18

        count_pos, _ = populated_db.filter(
            ligand_substruct_pos=[
                ["[C][Oh]", 1, 10, 102, 106, 154],
                ["C=O", 1, 10, 102, 106, 154],
            ],
            output_bookmark="substruct_pos",
        )
        assert count_pos == 12

        count_file, _ = populated_db.filter(
            ligand_name_file=str(TEST_DATA / "adgpu/ligand_names.csv"),
        )
        assert count_file == 16

    def test_hb_count_boundary(self, populated_db):
        """hb_count is inclusive: N admits poses with exactly N hydrogen bonds.
        The previous exclusive comparison yielded the "at least N+1" counts here
        (106 for hb_count=4, 81 for hb_count=5)."""
        count_4, _ = populated_db.filter(hb_count=4, output_bookmark="hb_least_4")
        assert count_4 == 141

        count_5, _ = populated_db.filter(hb_count=5, output_bookmark="hb_least_5")
        assert count_5 == 106

    def test_hb_count_at_most(self, populated_db):
        """A negative value means "no more than", also inclusive."""
        count, _ = populated_db.filter(hb_count=-1, output_bookmark="hb_most_1")
        assert count == 36

    def test_hb_count_zero(self, populated_db):
        """0 means "no more than 0": only poses with no hydrogen bonds. The sign
        carries the direction, so there is no -0 to express this with otherwise."""
        count, _ = populated_db.filter(hb_count=0, output_bookmark="hb_zero")
        assert count == 7

    def test_all_filters(self, populated_db):
        count, _ = populated_db.filter(
            eworst=-6,
            hb_interactions=[("A:VAL:279:", True), ("A:LYS:162:", True)],
            vdw_interactions=[("A:VAL:279:", True), ("A:LYS:162:", True)],
            max_miss=1,
            output_bookmark="big_query",
            ligand_name=["88"],
        )
        assert count == 1

    def test_generate_interaction_combinations(self, tmp_db):
        filters = Filters(
            hb_interactions=[("A:ARG:123:", True), ("A:VAL:124:", True)],
            vdw_interactions=[("A:ARG:123:", True), ("A:VAL:124:", True)],
        )
        combos = tmp_db._generate_interaction_combinations(filters.asdict(), 1)
        result_filters = [
            tmp_db._prepare_interaction_combo_filters(filters.asdict(), c)
            for c in combos
        ]
        assert len(result_filters) == 5
        assert (
            Filters(
                hb_interactions=[("A:ARG:123:", True), ("A:VAL:124:", True)],
                vdw_interactions=[("A:ARG:123:", True)],
            ).asdict()
            in result_filters
        )
        assert (
            Filters(
                hb_interactions=[("A:ARG:123:", True)],
                vdw_interactions=[("A:ARG:123:", True), ("A:VAL:124:", True)],
            ).asdict()
            in result_filters
        )


class TestOutput:
    """Output operations: SDFs, CSVs, logs, bookmark exports."""

    def test_get_filterdata(self, populated_db, tmp_path):
        populated_db.filter(eworst=-7, output_bookmark="has_filterdata")
        log_file = str(tmp_path / "filterdata.txt")
        populated_db.get_previous_filter_data(
            "has_filterdata", "delta, reference_rmsd", output_log=log_file
        )
        with open(log_file) as f:
            contents = f.read()
        assert "11991, 0.0, 226.06" in contents
        assert "3961, -0.02, 215.96" in contents

    def test_export_csv_and_log(self, populated_db, tmp_path):
        log_file = str(tmp_path / "filter_log.txt")
        populated_db.filter(
            eworst=-7,
            output_log=log_file,
            output_bookmark="export_csv",
            outfields=["ligname", "docking_score"],
        )

        # verify log content
        target_line = None
        with open(log_file) as f:
            lines = f.readlines()
        for i, line in enumerate(lines):
            if "bookmark" in line and i + 2 < len(lines):
                target_line = lines[i + 2].strip()
        assert target_line == "11128, -7.25"

        csv_ligands = str(tmp_path / "Ligands.csv")
        populated_db.export_table_as_csv("Ligands", csv_ligands)
        assert Path(csv_ligands).exists()

        csv_bookmark = str(tmp_path / "export_csv.csv")
        populated_db.export_table_as_csv("export_csv", csv_bookmark)
        assert Path(csv_bookmark).exists()

    def test_create_rdkitmol(self, populated_db):
        populated_db.filter(ebest=-3, output_bookmark="rdkit_test")
        ligands_poses = populated_db.fetch_select_ligands_poses(
            ligand_names=["14303"], bookmark_name="rdkit_test"
        )
        _, mol, _, _, _ = populated_db.create_rdkit_mols(ligands_poses["14303"])[0]
        assert mol.GetNumAtoms() == 10

    def test_write_sdfs(self, populated_db, tmp_path):
        sdf_dir = str(tmp_path / "sdf_files")
        populated_db.filter(eworst=-7, output_bookmark="sdf_bookmark")
        populated_db.write_molecule_sdfs("sdf_bookmark", sdf_dir, all_in_one=False)

        sdf_files = {f.stem for f in Path(sdf_dir).iterdir()}
        expected_lignames = {
            "3961",
            "5995",
            "11128",
            "11991",
            "13974",
            "15776",
            "136065",
            "127947",
        }
        assert sdf_files == expected_lignames

        with open(Path(sdf_dir) / "136065.sdf") as f:
            lines = f.readlines()
        assert lines[3] == " 27 28  0  0  0  0  0  0  0  0999 V2000\n"

    def test_export_bookmark_db(self, populated_db, tmp_path):
        populated_db.filter(eworst=-7, output_bookmark="export_db")
        bookmark_db_name = populated_db.export_bookmark_db("export_db")
        assert Path(bookmark_db_name).exists()

        rtc_bm = RingtailCore(db_file=bookmark_db_name)
        assert rtc_bm.table_length("Results") == 8

    def test_compress_decompress_db(self, populated_db, tmp_path):
        from ringtail.util import compress_file, decompress_file, detect_db_type
        import shutil as _sh

        src = populated_db.db_file
        full = populated_db.table_length("Results")

        # filtered subset -> compress (gzip) -> decompress -> open and verify count
        populated_db.filter(eworst=-7, output_bookmark="exp")
        subset = populated_db.export_bookmark_db("exp", str(tmp_path / "subset.db"))
        populated_db.delete_bookmark("exp")
        sub_count = RingtailCore(db_file=subset).table_length("Results")
        assert sub_count == 8

        art = compress_file(
            subset, str(tmp_path / "subset.db.gz"), method="gzip", level=6
        )
        assert art.endswith(".gz") and Path(art).exists()
        back = decompress_file(art, str(tmp_path / "back.db"))
        assert detect_db_type(back) in ("sqlite", "duckdb")
        assert RingtailCore(db_file=back).table_length("Results") == sub_count
        assert Path(src).exists()  # source database is never destroyed

        # no-filter: compress the whole db -> decompress -> full count preserved
        art2 = compress_file(src, str(tmp_path / "whole.db.gz"), method="gzip")
        back2 = decompress_file(art2, str(tmp_path / "whole_back.db"))
        assert RingtailCore(db_file=back2).table_length("Results") == full

        # zstd -> gzip fallback when the zstd binary is unavailable
        orig_which = _sh.which
        _sh.which = lambda x: None if x == "zstd" else orig_which(x)
        try:
            fb = compress_file(src, str(tmp_path / "fb.db.zst"), method="zstd")
        finally:
            _sh.which = orig_which
        assert fb.endswith(".gz")

    def test_similar_ligands_interaction(self, populated_db, tmp_path):
        populated_db.filter(ebest=-6, interaction_cluster=0.5)
        options = populated_db.fetch_cluster_options("28837")
        assert len(options) > 0
        cluster_id = options[0][0]
        ligands, bookmark_name, cluster_name = populated_db.fetch_clustered_similars(
            "28837", cluster_id, output_log=str(tmp_path / "cluster_log.txt")
        )
        assert len(ligands) == 13

    def test_similar_ligands_mfpt(self, populated_db, tmp_path):
        populated_db.filter(ebest=-6, mfpt_cluster=0.5)
        options = populated_db.fetch_cluster_options("287065")
        assert len(options) > 0
        cluster_id = options[0][0]
        ligands, bookmark_name, cluster_name = populated_db.fetch_clustered_similars(
            "287065", cluster_id, output_log=str(tmp_path / "cluster_log.txt")
        )
        assert len(ligands) == 8


class TestStorageMan:
    def test_bookmark_info(self, tmp_db: RingtailCore):
        tmp_db.add_results_from_files(
            docking_results=str(TEST_DATA / "adgpu/group2"), docking_mode="adgpu"
        )
        tmp_db.filter(
            eworst=-3,
            hb_interactions=[("A:VAL:279:", True), ("A:LYS:162:", True)],
            vdw_interactions=[("A:VAL:279:", True)],
            output_bookmark="bookmark_info",
        )
        qb = QueryBuilder()
        query = (
            qb.SELECT("filters")
            .FROM("Filters")
            .WHERE("name='bookmark_info'")
            .build()[0]
        )
        bookmark_filters_db_str = tmp_db.db_query(query)[0][0]
        assert (
            json.loads(bookmark_filters_db_str)
            == Filters(
                eworst=-3.0,
                vdw_interactions=[["A:VAL:279:", True]],
                hb_interactions=[["A:VAL:279:", True], ["A:LYS:162:", True]],
            ).asdict()
        )

    def test_version_info(self, tmp_db):
        from importlib.metadata import version

        tmp_db.add_results_from_files(
            docking_results=str(TEST_DATA / "adgpu/group1/1451.dlg.gz"),
            max_poses=1,
            docking_mode="adgpu",
        )
        with tmp_db.storageman:
            versionmatch, db_version = tmp_db.storageman.check_ringtaildb_version()
        assert versionmatch
        assert db_version == version("ringtail")


class TestMergeDB:
    def test_merge_workflow(self, tmp_path, storage_type):
        db1 = str(tmp_path / "primary.db")
        db2 = str(tmp_path / "secondary.db")
        db3 = str(tmp_path / "tertiary.db")

        rtc1 = RingtailCore(db1, storage_type=storage_type)
        rtc1.add_results_from_files(
            str(TEST_DATA / "adgpu/group1/1451.dlg.gz"), docking_mode="adgpu"
        )
        rtc2 = RingtailCore(db2, storage_type=storage_type)
        rtc2.add_results_from_files(
            str(TEST_DATA / "adgpu/group1/1620.dlg.gz"), docking_mode="adgpu"
        )
        rtc3 = RingtailCore(db3, storage_type=storage_type)
        rtc3.add_results_from_files(
            str(TEST_DATA / "adgpu/group1/1751.dlg.gz"), docking_mode="adgpu"
        )

        assert (
            rtc1.table_length("Ligands")
            == rtc2.table_length("Ligands")
            == rtc3.table_length("Ligands")
            == 1
        )

        # before merge: no poses in the tight interval, one in the loose interval
        assert rtc1.filter(eworst=-2, ebest=-5)[0] == 0
        assert rtc1.filter(eworst=-5)[0] == 1

        # merge secondary and tertiary into primary
        rtc1 = RingtailCore(db1)
        rtc1.merge_databases(db2, False)
        rtc1.merge_databases(db3, False)
        assert rtc1.table_length("Ligands") == 3
        assert rtc1.filter(eworst=-2, ebest=-5)[0] == 2

        # PKs in secondary should be reassigned in the merged db
        secondary_pose_in_own_db = RingtailCore(db2).db_query(
            "SELECT pose_id FROM Results WHERE pose_rank = 1 AND ligand_id = "
            "(SELECT ligand_id FROM Ligands WHERE ligname = '1620')"
        )[0][0]
        assert secondary_pose_in_own_db == 1

        secondary_pose_in_merged = rtc1.db_query(
            "SELECT pose_id FROM Results WHERE pose_rank = 1 AND ligand_id = "
            "(SELECT ligand_id FROM Ligands WHERE ligname = '1620')"
        )[0][0]
        assert secondary_pose_in_merged != secondary_pose_in_own_db


class TestCrossref:
    """Cross-referencing ligands across databases by bookmark and by status table."""

    @staticmethod
    def _build_db(path, storage_type, ligands):
        rtc = RingtailCore(path, storage_type=storage_type)
        rtc.add_results_from_files(
            docking_results=[
                str(TEST_DATA / f"adgpu/group1/{lig}.dlg.gz") for lig in ligands
            ],
            docking_mode="adgpu",
        )
        return rtc

    @staticmethod
    def _accept_ligand(rtc, ligname):
        """Mark the best pose of a ligand as Accepted (status 1)."""
        pose_id = rtc.db_query(
            "SELECT pose_id FROM Results WHERE pose_rank = 1 AND ligand_id = "
            f"(SELECT ligand_id FROM Ligands WHERE ligname = '{ligname}')"
        )[0][0]
        rtc.update_pose_status(pose_id, 1)

    @staticmethod
    def _bookmark_exists(path, name):
        return name in RingtailCore(path).get_bookmark_names()

    def test_crossref_status_tables(self, tmp_path, storage_type):
        """Crossref two dbs scoped on the Accepted status table; only the ligand
        accepted in BOTH databases should pass."""
        dbA = str(tmp_path / "targetA.db")
        dbB = str(tmp_path / "targetB.db")
        # shared ligand "1451"; "1620"/"1751" are unique to one db each
        rtcA = self._build_db(dbA, storage_type, ["1451", "1620"])
        rtcB = self._build_db(dbB, storage_type, ["1451", "1751"])

        self._accept_ligand(rtcA, "1451")
        self._accept_ligand(rtcA, "1620")
        self._accept_ligand(rtcB, "1451")
        self._accept_ligand(rtcB, "1751")

        rtcA = RingtailCore(dbA, storage_type=storage_type)
        count, new_bookmarks, _ = rtcA.cross_reference_databases(
            wanted_dbs=[(dbA, "accepted"), (dbB, "accepted")],
        )

        assert count == 1  # only "1451" is accepted in both
        assert new_bookmarks[dbA] == "crossref_accepted"
        assert new_bookmarks[dbB] == "crossref_accepted"
        assert self._bookmark_exists(dbA, "crossref_accepted")
        assert self._bookmark_exists(dbB, "crossref_accepted")

    def test_crossref_mixed_scope(self, tmp_path, storage_type):
        """A status table in wanted_dbs interoperates with a bookmark in
        unwanted_dbs: the shared accepted ligand is excluded by the bookmark."""
        dbA = str(tmp_path / "targetA.db")
        dbB = str(tmp_path / "targetB.db")
        dbC = str(tmp_path / "offtarget.db")
        rtcA = self._build_db(dbA, storage_type, ["1451", "1620"])
        rtcB = self._build_db(dbB, storage_type, ["1451", "1751"])
        rtcC = self._build_db(dbC, storage_type, ["1451"])

        self._accept_ligand(rtcA, "1451")
        self._accept_ligand(rtcB, "1451")
        # bookmark in the off-target db that captures the shared ligand
        rtcC = RingtailCore(dbC, storage_type=storage_type)
        assert rtcC.filter(eworst=0, output_bookmark="offtarget_hits")[0] >= 1

        rtcA = RingtailCore(dbA, storage_type=storage_type)
        count, _, _ = rtcA.cross_reference_databases(
            wanted_dbs=[(dbA, "accepted"), (dbB, "accepted")],
            unwanted_dbs=[(dbC, "offtarget_hits")],
        )

        assert count == 0  # "1451" passes wanted intersect but is excluded


class TestADGPUHandling:
    def test_reactive_filtering(self, tmp_db):
        tmp_db.add_results_from_files(
            docking_results=str(TEST_DATA / "reactive"),
            store_all_poses=True,
            receptor_file=str(TEST_DATA / "reactive/4j8m_m_rigid.pdbqt"),
            docking_mode="adgpu",
        )
        count, _ = tmp_db.filter(reactive_interactions=[("A:TYR:212:", True)])
        assert count == 10


class TestVinaHandling:
    def test_file_add(self, tmp_db):
        vina_path = TEST_DATA / "vina"
        tmp_db.add_results_from_files(
            docking_results=str(vina_path),
            receptor_file=str(vina_path / "receptor.pdbqt"),
            save_receptor=True,
            docking_mode="vina",
        )
        assert tmp_db.table_length("Results") == 6

    def test_string_add(self, tmp_db):
        vina_path = TEST_DATA / "vina"
        sample1 = (vina_path / "sample-result.pdbqt").read_text()
        sample2 = (vina_path / "sample-result-2.pdbqt").read_text()
        tmp_db.save_receptor(str(vina_path / "receptor.pdbqt"))
        tmp_db.add_results_from_vina_string(
            results={"sample1": sample1, "sample2": sample2}
        )
        assert tmp_db.table_length("Results") == 6

    def test_add_interactions(self, tmp_db):
        vina_path = TEST_DATA / "vina"
        tmp_db.add_results_from_files(
            docking_results=str(vina_path),
            receptor_file=str(vina_path / "receptor.pdbqt"),
            save_receptor=True,
            docking_mode="vina",
        )
        assert tmp_db.table_length("Interaction_indices") == 32
        assert tmp_db.table_length("Interactions") == 77

    def test_add_interactions_from_polymer(self, tmp_path, storage_type):
        data_path = TEST_DATA / "flexres"
        rtc_json = RingtailCore(
            str(tmp_path / "flexres_json.db"), storage_type=storage_type
        )
        rtc_json.add_results_from_files(
            docking_results=str(data_path / "ligand.pdbqt"),
            docking_mode="vina",
            receptor_file=str(data_path / "receptor.json"),
            save_receptor=True,
        )
        rtc_pdbqt = RingtailCore(
            str(tmp_path / "flexres_pdbqt.db"), storage_type=storage_type
        )
        rtc_pdbqt.add_results_from_files(
            docking_results=str(data_path / "ligand.pdbqt"),
            docking_mode="vina",
            receptor_file=str(data_path / "receptor.pdbqt"),
            save_receptor=True,
        )
        assert (
            rtc_json.table_length("Ligands") == rtc_pdbqt.table_length("Ligands") == 1
        )
        assert (
            rtc_json.table_length("Interactions")
            == rtc_pdbqt.table_length("Interactions")
            == 38
        )

    def test_polymer_receptor(self, tmp_db):
        data_path = TEST_DATA / "flexres"
        tmp_db.add_results_from_files(
            docking_results=str(data_path / "ligand.pdbqt"),
            docking_mode="vina",
            receptor_file=str(data_path / "receptor.json"),
            save_receptor=True,
        )
        receptor_items = tmp_db.get_receptor_object()
        assert receptor_items.name == "receptor"
        assert not receptor_items.blob_str
        assert receptor_items.polymer_json is not None

    def test_write_flexres_pdb(self, tmp_db, tmp_path):
        pytest.importorskip("meeko")
        import meeko

        data_path = TEST_DATA / "flexres"
        tmp_db.add_results_from_files(
            docking_results=str(data_path / "ligand.pdbqt"),
            receptor_file=str(data_path / "receptor.pdbqt"),
            recursive=True,
            docking_mode="vina",
        )
        tmp_db.filter(eworst=-1, output_bookmark="flexres")
        polymer = meeko.Polymer.from_json((data_path / "receptor.json").read_text())
        export_base = str(tmp_path / "exported_flex_rec")
        tmp_db.write_flexres_pdb(polymer, "ligand", "flexres", export_base)

        expected = Path(f"{export_base}_ligand.pdb")
        assert expected.exists()
        content = expected.read_text()
        assert (
            "ATOM     11  C   HIS A   2       2.368   0.239  -0.349                       C"
            in content
        )
        assert (
            "ATOM     44 HD2  HIS A   3      -0.400  -5.308  -1.507                       H"
            in content
        )

    def test_various_filters(self, tmp_db):
        vina_path = TEST_DATA / "vina"
        tmp_db.add_results_from_files(
            docking_results=str(vina_path),
            receptor_file=str(vina_path / "receptor.pdbqt"),
            save_receptor=True,
            docking_mode="vina",
        )
        count, _ = tmp_db.filter(eworst=-6, ligand_substruct=["[N]"])
        assert count == 1

    def test_db_dockingmode_warning(self, tmp_db, tmp_path):
        from ringtail import setup_logging

        logfile = str(tmp_path / "dockingmode_warning.log")
        setup_logging(level="DEBUG", logfile=logfile)

        tmp_db.add_results_from_files(
            docking_results=str(TEST_DATA / "adgpu/group1/1451.dlg.gz"),
            docking_mode="adgpu",
        )
        rtc2 = RingtailCore(tmp_db.db_file, storage_type=tmp_db.storagetype)
        rtc2.add_results_from_files(
            docking_results=str(TEST_DATA / "vina/sample-result.pdbqt"),
            docking_mode="vina",
        )
        warning = (
            "The following database properties do not agree with the properties last used for this database: \n"
            "Current docking mode is vina but last used docking mode of database is adgpu."
        )
        with open(logfile) as f:
            assert warning in f.read()


class TestAD6Handling:
    def test_stream(self, tmp_db):
        rdkit = pytest.importorskip("rdkit")
        from rdkit import Chem

        tmp_db.save_receptor(str(TEST_DATA / "ad6/helix--scofu01.json"))
        suppl = Chem.SDMolSupplier(
            str(TEST_DATA / "ad6/docked_ligands.sdf"), removeHs=False
        )
        tmp_db.add_mol(suppl)
        assert tmp_db.table_length("Results") == 9
        assert tmp_db.table_length("Interactions") == 60

    def test_file_add(self, tmp_db):
        ad6_path = TEST_DATA / "ad6"
        tmp_db.add_results_from_files(
            docking_results=str(ad6_path),
            receptor_file=str(ad6_path / "helix--scofu01.json"),
            save_receptor=True,
            docking_mode="ad6",
        )
        assert tmp_db.table_length("Results") == 9
        assert tmp_db.table_length("Interactions") == 60

    def test_file_add_no_interactions(self, tmp_db):
        ad6_path = TEST_DATA / "ad6"
        tmp_db.add_results_from_files(
            docking_results=str(ad6_path),
            calculate_interactions=False,
            docking_mode="ad6",
        )
        assert tmp_db.table_length("Results") == 9
        assert tmp_db.table_length("Interactions") == 0

    def test_calc_interactions_deferred(self, tmp_db):
        ad6_path = TEST_DATA / "ad6"
        tmp_db.add_results_from_files(
            docking_results=str(ad6_path),
            calculate_interactions=False,
            docking_mode="ad6",
        )
        assert tmp_db.table_length("Interactions") == 0

        tmp_db.save_receptor(str(ad6_path / "helix--scofu01.json"))
        tmp_db.add_interactions()
        assert tmp_db.table_length("Results") == 9
        assert tmp_db.table_length("Interactions") == 60

    def test_add_interactions_recalc_larger_cutoffs(self, tmp_db):
        # populate with interactions at the default cutoffs (3.7 HB, 4.0 VDW)
        ad6_path = TEST_DATA / "ad6"
        tmp_db.add_results_from_files(
            docking_results=str(ad6_path),
            receptor_file=str(ad6_path / "helix--scofu01.json"),
            save_receptor=True,
            docking_mode="ad6",
        )
        interactions_before = tmp_db.table_length("Interactions")
        hb_before = tmp_db.db_query("SELECT SUM(num_hb) FROM Results")[0][0]
        assert interactions_before > 0

        # recalculating over existing interactions requires consent; without it the
        # call is a no-op and the existing interactions are kept
        tmp_db.add_interactions(hb_cutoff=6.0, vdw_cutoff=7.0)
        assert tmp_db.table_length("Interactions") == interactions_before

        # with consent, existing interactions are deleted and recomputed with the
        # larger cutoffs, which captures more contacts -> more interactions and more
        # hydrogen bonds (num_hb is recomputed in the Results table)
        tmp_db.add_interactions(hb_cutoff=6.0, vdw_cutoff=7.0, consent=True)
        interactions_after = tmp_db.table_length("Interactions")
        hb_after = tmp_db.db_query("SELECT SUM(num_hb) FROM Results")[0][0]

        assert tmp_db.table_length("Results") == 9  # poses themselves unchanged
        assert interactions_after > interactions_before
        assert hb_after > hb_before

    def test_filtering(self, tmp_db):
        ad6_path = TEST_DATA / "ad6"
        tmp_db.add_results_from_files(
            docking_results=str(ad6_path),
            receptor_file=str(ad6_path / "helix--scofu01.json"),
            save_receptor=True,
            docking_mode="ad6",
        )
        count, _ = tmp_db.filter(
            eworst=-13, ligand_substruct=["C=O"], vdw_interactions=[(":VAL::", True)]
        )
        assert count == 1


class TestLogger:
    def test_set_log_level(self):
        from ringtail import LOGGER, setup_logging

        setup_logging(level="INFO")
        assert LOGGER.level == logging.INFO


class TestOptions:
    def test_filter_option_checks(self, tmp_db, tmp_path):
        from ringtail.ringtailoptions import Filters

        filelist = tmp_path / "filelist.txt"
        filelist.write_text(
            "\n".join(
                str(TEST_DATA / "adgpu/group1" / f)
                for f in ["127458.dlg.gz", "173101.dlg.gz", "100729.dlg.gz"]
            )
        )
        tmp_db.add_results_from_files(
            docking_results=str(filelist), docking_mode="adgpu"
        )
        tmp_db.filters = Filters(score_percentile=20)
        assert tmp_db.filters.eworst is None
        assert tmp_db.filters.score_percentile == 20

        tmp_db.filters = Filters(score_percentile=20, eworst=-6)
        tmp_db.filters.checks()
        assert tmp_db.filters.eworst == -6
        assert tmp_db.filters.score_percentile is None

    def test_overwrite_db(self, tmp_db, tmp_path):
        list1 = tmp_path / "list1.txt"
        list1.write_text(
            "\n".join(
                str(TEST_DATA / "adgpu/group1" / f)
                for f in ["127458.dlg.gz", "173101.dlg.gz", "100729.dlg.gz"]
            )
        )
        list2 = tmp_path / "list2.txt"
        list2.write_text(
            str(TEST_DATA / "adgpu/group1/272275.dlg.gz")
            + "\n"
            + str(TEST_DATA / "adgpu/group3/60239.dlg.gz")
            + "\n"
        )
        tmp_db.add_results_from_files(docking_results=str(list1), docking_mode="adgpu")
        count_old = tmp_db.table_length("Ligands")
        tmp_db.add_results_from_files(
            docking_results=str(list2), docking_mode="adgpu", overwrite=True
        )
        count_new = tmp_db.table_length("Ligands")
        assert count_old == 3
        assert count_new == 2
