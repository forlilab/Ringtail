#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail unit tests
#
import json
from pathlib import Path

import pytest
from ringtail import (
    Filters,
    QueryBuilder,
    RECALC_TRACKING_TABLE,
    RingtailCore,
)
from ringtail.exceptions import OptionError

TEST_DATA = Path(__file__).parent / "test_data"


def duplicate_pairs(rtc) -> int:
    """How many (pose_id, interaction_id) pairs are stored more than once."""
    return rtc.db_query(
        """SELECT COUNT(*) FROM (SELECT pose_id, interaction_id, COUNT(*) c
           FROM Interactions GROUP BY pose_id, interaction_id HAVING c > 1)"""
    )[0][0]


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
        rows = tmp_db.db_query(
            "SELECT pose_rank, run_number, docking_score FROM Results ORDER BY pose_rank"
        )
        assert sorted(run for _, run, _ in rows) == list(range(1, 21))
        scores = [score for _, _, score in rows]
        assert scores == sorted(scores)

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

    def test_save_receptor(self, populated_db):
        count_before = populated_db.db_query(
            "SELECT COUNT(*) FROM Receptors WHERE receptor_object NOT NULL"
        )[0][0]
        assert count_before == 0

        populated_db.save_receptor(receptor=str(TEST_DATA / "adgpu/4j8m.pdbqt"))
        count_after = populated_db.db_query(
            "SELECT COUNT(*) FROM Receptors WHERE receptor_object NOT NULL"
        )[0][0]
        assert count_after == 1

    def test_db_summary_data(self, tmp_db, populated_db):
        from ringtail import exceptions as e

        with pytest.raises(e.StorageError):
            tmp_db.db_summary_data()

        data, _ = populated_db.db_summary_data()
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

    def test_adding_results_with_filters_needs_consent(self, ad6_db):
        """Adding results to a db with bookmarks needs consent, and deletes all filter data"""
        sdf = str(TEST_DATA / "ad6/docked_ligands.sdf")
        ad6_db.filter(eworst=0, output_bookmark="before")
        ad6_db.filter(eworst=0, mfpt_cluster=0.6, output_bookmark="clustered")
        result_count = ad6_db.table_length("Results")

        with pytest.raises(OptionError, match="consent=True"):
            ad6_db.add_results_from_files(docking_results=sdf)
        assert ad6_db.table_length("Results") == result_count
        assert "before" in ad6_db.get_bookmark_names()

        ad6_db.add_results_from_files(docking_results=sdf, consent=True)
        assert ad6_db.table_length("Results") == result_count * 2
        assert ad6_db.get_bookmark_names() == []
        assert not ad6_db.has_filter_data()
        for table in ("Filters", "Filtered_poses", "Clusters", "Pose_clusters"):
            assert ad6_db.table_length(table) == 0


class TestFiltering:
    """Filter operations on the full 217-ligand adgpu dataset."""

    @pytest.mark.parametrize(
        ("percentile", "expected_count"),
        [(0, 0), (5, 0), (25, 1), (50, 2), (100, 4)],
    )
    def test_score_percentile_floor_rounding(self, ad6_db, percentile, expected_count):
        """Score percentile rounds down and includes poses at the cutoff"""
        count, _ = ad6_db.filter(
            score_percentile=percentile,
            output_bookmark=f"score_percentile_{percentile}",
        )
        assert count == expected_count

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
        """hb_count includes poses with exactly that many hydrogen bonds"""
        count_4, _ = populated_db.filter(hb_count=4, output_bookmark="hb_least_4")
        assert count_4 == 141

        count_5, _ = populated_db.filter(hb_count=5, output_bookmark="hb_least_5")
        assert count_5 == 106

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

    def test_tiered_filter_or_of_groups(self, populated_db):
        """OR of two filter groups gives the union of each group filtered on its own"""
        expr = {
            "op": "or",
            "children": [
                {
                    "eworst": -6,
                    "hb_interactions": [("A:VAL:279:", True), ("A:LYS:162:", True)],
                    "vdw_interactions": [("A:VAL:279:", True), ("A:LYS:162:", True)],
                    "max_miss": 1,
                },
                {"eworst": -7},
            ],
        }
        count_or, _ = populated_db.filter(filters=expr, output_bookmark="tier_or")
        # each tier on its own
        populated_db.filter(
            eworst=-6,
            hb_interactions=[("A:VAL:279:", True), ("A:LYS:162:", True)],
            vdw_interactions=[("A:VAL:279:", True), ("A:LYS:162:", True)],
            max_miss=1,
            output_bookmark="tier_a",
        )
        populated_db.filter(eworst=-7, output_bookmark="tier_b")

        ligs_a = set(
            populated_db.fetch_select_ligands_poses(bookmark_name="tier_a").keys()
        )
        ligs_b = set(
            populated_db.fetch_select_ligands_poses(bookmark_name="tier_b").keys()
        )
        ligs_or = set(
            populated_db.fetch_select_ligands_poses(bookmark_name="tier_or").keys()
        )
        assert ligs_or == ligs_a | ligs_b
        assert count_or == len(ligs_a | ligs_b)

    def test_smarts_inside_group(self, populated_db):
        """A SMARTS filter in a filter group is combined with the other criteria in that group"""
        score_only = {"op": "and", "children": [{"eworst": -6}]}
        with_smarts = {
            "op": "and",
            "children": [{"eworst": -6}, {"ligand_substruct": ["C=O"]}],
        }
        c_score, _ = populated_db.filter(
            filters=score_only, output_bookmark="sig_score"
        )
        c_both, _ = populated_db.filter(filters=with_smarts, output_bookmark="sig_both")
        assert c_both > 0
        assert c_both <= c_score  # adding SMARTS only narrows


class TestOutput:
    """Output operations: SDFs, CSVs, logs, bookmark exports."""

    def test_best_passing_pose_and_order_results(self, ad6_db, tmp_path):
        rows = ad6_db.filter(
            eworst=0,
            output_bookmark="ordered_best_poses",
            output_log=str(tmp_path / "ordered_best_poses.txt"),
            outfields=["ligname", "docking_score", "pose_rank"],
            order_results="docking_score",
            return_iter=True,
        )

        assert [row[0] for row in rows] == [
            "second_mol",
            "third_mol",
            "first_mol",
            "fourth_mol",
        ]
        assert [row[1] for row in rows] == pytest.approx(
            [-16.74, -14.77, -13.04, -12.44]
        )
        assert [row[2] for row in rows] == [1, 1, 1, 1]

    def test_get_filterdata(self, populated_db, tmp_path):
        populated_db.filter(eworst=-7, output_bookmark="has_filterdata")
        log_file = str(tmp_path / "filterdata.txt")
        populated_db.get_previous_filter_data(
            "has_filterdata", "delta, reference_rmsd", output_log=log_file
        )
        with open(log_file) as f:
            contents = f.read()
        assert "11991, 0.0, 226.06" in contents
        assert "3961, 0.0, 215.96" in contents

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
        populated_db.export_table_to_csv("Ligands", csv_ligands)
        assert Path(csv_ligands).exists()

        csv_bookmark = str(tmp_path / "export_csv.csv")
        populated_db.export_table_to_csv("export_csv", csv_bookmark)
        assert Path(csv_bookmark).exists()

    def test_export_sql_to_csv(self, populated_db, tmp_path):
        csv_file = tmp_path / "query.csv"
        populated_db.export_sql_to_csv(
            "SELECT pose_id, docking_score FROM Results", str(csv_file)
        )
        lines = csv_file.read_text().splitlines()
        assert lines[0] == "pose_id,docking_score"
        assert len(lines) - 1 == populated_db.table_length("Results")

    def test_pose_comments(self, ad6_db):
        pose_id = ad6_db.db_query("SELECT pose_id FROM Results LIMIT 1")[0][0]

        assert ad6_db.get_pose_comment(pose_id) is None
        ad6_db.set_pose_comment(pose_id, "first comment")
        assert ad6_db.get_pose_comment(pose_id) == "first comment"

        ad6_db.set_pose_comment(pose_id, "updated comment")
        assert ad6_db.get_pose_comment(pose_id) == "updated comment"

        ad6_db.set_pose_comment(pose_id, "")
        assert ad6_db.get_pose_comment(pose_id) is None

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

    def test_histogram_matches_numpy(self, populated_db):
        import numpy as np

        scores = [
            row[0] for row in populated_db.db_query("SELECT docking_score FROM Results")
        ]
        edges, counts = populated_db.get_histogram_data("docking_score", 5)
        expected, _ = np.histogram(scores, bins=edges)
        assert counts == expected.tolist()

    def test_compress_decompress_db(self, populated_db, tmp_path, monkeypatch):
        import shutil
        from ringtail.util import compress_file, decompress_file, detect_db_type

        src = populated_db.db_file
        art = compress_file(src, str(tmp_path / "whole.db.gz"), method="gzip")
        assert art.endswith(".gz")
        back = decompress_file(art, str(tmp_path / "back.db"))
        assert detect_db_type(back) in ("sqlite", "duckdb")
        assert RingtailCore(db_file=back).table_length(
            "Results"
        ) == populated_db.table_length("Results")
        assert Path(src).exists()

        # zstd falls back to gzip if the zstd binary is missing
        real_which = shutil.which
        monkeypatch.setattr(
            shutil, "which", lambda name: None if name == "zstd" else real_which(name)
        )
        fallback = compress_file(src, str(tmp_path / "fallback.db.zst"), method="zstd")
        assert fallback.endswith(".gz")

    @pytest.mark.parametrize(
        ("cluster", "ligname", "expected"),
        [("interaction_cluster", "28837", 13), ("mfpt_cluster", "287065", 8)],
    )
    def test_similar_ligands(self, populated_db, cluster, ligname, expected):
        populated_db.filter(ebest=-6, **{cluster: 0.5})
        options = populated_db.fetch_cluster_options(ligname)
        assert len(options) > 0
        ligands, _, _ = populated_db.fetch_clustered_similars(ligname, options[0][0])
        assert len(ligands) == expected


class TestInteractionAnalysis:
    """Polymer JSON and pdbqt receptors give the same interactions"""

    @staticmethod
    def _from_db(rtc: RingtailCore) -> dict:
        """Interactions in the db as {(ligname, run_number): set of interactions}"""
        rows = rtc.db_query("""SELECT L.ligname, R.run_number, II.interaction_type,
                      II.rec_resname, II.rec_resid, II.rec_chain, II.rec_atom
               FROM Interactions I
               JOIN Interaction_indices II ON II.interaction_id = I.interaction_id
               JOIN Results R ON R.pose_id = I.pose_id
               JOIN Ligands L ON L.ligand_id = R.ligand_id""")
        poses = {}
        for ligname, run, itype, resname, resid, chain, recname in rows:
            poses.setdefault((ligname, run), set()).add(
                (itype, resname, str(resid), chain, recname)
            )
        return poses

    def test_receptor_formats_agree(self, flexres_json_db, flexres_pdbqt_db):
        assert (
            flexres_json_db.table_length("Interactions")
            == flexres_pdbqt_db.table_length("Interactions")
            == 25
        )
        assert self._from_db(flexres_json_db) == self._from_db(flexres_pdbqt_db)


class TestStorageMan:
    def test_bookmark_info(self, populated_db: RingtailCore):
        populated_db.filter(
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
        bookmark_filters_db_str = populated_db.db_query(query)[0][0]
        assert (
            json.loads(bookmark_filters_db_str)
            == Filters(
                eworst=-3.0,
                vdw_interactions=[["A:VAL:279:", True]],
                hb_interactions=[["A:VAL:279:", True], ["A:LYS:162:", True]],
            ).to_dict()
        )

    def test_bookmark_dependents(self, populated_db: RingtailCore):
        populated_db.filter(eworst=-6, output_bookmark="parent_bm")
        populated_db.filter(
            eworst=-7, input_bookmark="parent_bm", output_bookmark="child_bm"
        )
        populated_db.filter(eworst=-7, output_bookmark="unrelated_bm")

        assert populated_db.get_bookmark_dependents("parent_bm") == ["child_bm"]
        assert populated_db.get_bookmark_dependents("child_bm") == []
        assert populated_db.get_bookmark_dependents("unrelated_bm") == []

        # deleting the parent keeps the child bookmark and its poses
        child_poses = populated_db.table_length("child_bm")
        assert child_poses > 0
        populated_db.delete_bookmark("parent_bm")

        assert "parent_bm" not in populated_db.get_bookmark_names()
        assert "child_bm" in populated_db.get_bookmark_names()
        assert populated_db.table_length("child_bm") == child_poses

    def test_version_info(self, ad6_db):
        from ringtail.schema import SCHEMA_VERSION

        with ad6_db.storageman:
            versionmatch, db_version = ad6_db.storageman.check_ringtaildb_version()

        assert versionmatch
        assert db_version == SCHEMA_VERSION

    def test_bookmark_paging(self, populated_db):
        """Page through a bookmark the way the GUI table does"""
        populated_db.build_gui_tables()
        populated_db.filter(eworst=-7, output_bookmark="page")
        start = populated_db.get_starting_rowid("page")
        page = populated_db.get_scrolling_table_data("page", 3, start)
        pose_col = page["headers"].index("pose_id")
        poses = [row[pose_col] for row in page["data"]]
        assert len(poses) == 3

        row = populated_db.get_row_from_pose("page", poses[1])
        following = populated_db.get_scrolling_table_data("page", 1, row + 1)["data"]
        assert following[0][pose_col] == poses[2]


@pytest.mark.slow  # clones a database and merges another into it
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
        assert rtc1.merge_databases([db2, db3], False) == []
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

    def test_merge_wipes_statuses_and_comments(self, tmp_path, storage_type):
        dbs = []
        for name, ligand in [("primary", "1451"), ("secondary", "1620")]:
            db = str(tmp_path / f"{name}.db")
            RingtailCore(db, storage_type=storage_type).add_results_from_files(
                str(TEST_DATA / f"adgpu/group1/{ligand}.dlg.gz"), docking_mode="adgpu"
            )
            dbs.append(db)
        rtc1 = RingtailCore(dbs[0])
        rtc1.build_gui_tables()
        rtc1.update_pose_status([1, 2], 1)
        rtc1.set_pose_comment(1, "wiped by the merge")
        rtc1.filter(eworst=-5, output_bookmark="stale")

        assert rtc1.merge_databases(dbs[1:], False) == []
        assert rtc1.table_length("Accepted") == 0
        assert rtc1.get_pose_comment(1) is None
        assert rtc1.get_bookmark_names() == []
        # status tables come back empty, so statuses can be assigned again
        rtc1.update_pose_status(3, 2)
        assert rtc1.table_length("Maybe") == 1


@pytest.mark.slow  # attaches extra databases and cross-references them
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
        """Crossref on the accepted status table only passes ligands accepted in both dbs"""
        dbA = str(tmp_path / "targetA.db")
        dbB = str(tmp_path / "targetB.db")
        # 1451 is in both dbs, 1620 and 1751 are only in one each
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
        """A bookmark in unwanted_dbs excludes ligands from a status table in wanted_dbs"""
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

    def test_crossref_gui_style_call(self, tmp_path, storage_type):
        """The GUI passes '' for this database and an alias for each attached one."""
        dbA = str(tmp_path / "targetA.db")
        dbB = str(tmp_path / "targetB.db")
        rtcA = self._build_db(dbA, storage_type, ["1451", "1620"])
        rtcB = self._build_db(dbB, storage_type, ["1451", "1751"])
        rtcA.filter(eworst=0, output_bookmark="all_a")
        rtcB.filter(eworst=0, output_bookmark="all_b")

        count, bookmarks, _ = rtcA.cross_reference_databases(
            [("", "all_a"), (dbB, "all_b")], [], "crossref_gui", {dbB: "target_b"}
        )
        assert count == 1
        assert bookmarks == {"": "crossref_gui_all_a", dbB: "crossref_gui_all_b"}
        # every pose of the shared ligand, not only its best one
        shared_poses = rtcA.db_query(
            "SELECT COUNT(*) FROM Results R JOIN Ligands L USING(ligand_id) "
            "WHERE L.ligname = '1451'"
        )[0][0]
        assert rtcA.table_length("crossref_gui_all_a") == shared_poses > 1


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

    def test_interaction_tolerance(self, tmp_db):
        tmp_db.add_results_from_files(
            docking_results=str(TEST_DATA / "adgpu/group1/127458.dlg.gz"),
            docking_mode="adgpu",
            interaction_tolerance=2.0,
        )
        # the top poses have 53 interactions, the tolerated poses add 4 more
        assert tmp_db.table_length("Interactions") == 57


class TestVinaHandling:
    def test_file_add(self, vina_db):
        assert vina_db.table_length("Results") == 6
        assert vina_db.table_length("Interaction_indices") == 25
        assert vina_db.table_length("Interactions") == 60

    def test_string_add(self, tmp_db):
        vina_path = TEST_DATA / "vina"
        sample1 = (vina_path / "sample-result.pdbqt").read_text()
        sample2 = (vina_path / "sample-result-2.pdbqt").read_text()
        tmp_db.save_receptor(str(vina_path / "receptor.pdbqt"))
        tmp_db.add_results_from_vina_string(
            results={"sample1": sample1, "sample2": sample2}
        )
        assert tmp_db.table_length("Results") == 6

    def test_string_add_without_receptor_or_interactions(self, tmp_db):
        vina_path = TEST_DATA / "vina/sample-result.pdbqt"
        tmp_db.add_results_from_vina_string(
            results={"sample": vina_path.read_text()},
            calculate_interactions=False,
        )
        assert tmp_db.table_length("Results") == 3
        assert tmp_db.table_length("Interactions") == 0

    def test_string_add_with_supplied_receptor(self, tmp_db):
        vina_path = TEST_DATA / "vina"
        tmp_db.add_results_from_vina_string(
            results={"sample": (vina_path / "sample-result.pdbqt").read_text()},
            receptor_string=(vina_path / "receptor.pdbqt").read_text(),
        )
        assert tmp_db.table_length("Results") == 3
        assert tmp_db.table_length("Interactions") > 0

    def test_polymer_receptor(self, flexres_json_db):
        receptor_items = flexres_json_db.get_receptor_object()
        assert receptor_items.name == "receptor"
        assert not receptor_items.blob_str
        assert receptor_items.polymer_json is not None

    def test_write_flexres_pdb(self, tmp_db, tmp_path):
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
            calculate_interactions=False,
        )
        warning = (
            "The following database properties do not agree with the properties last used for this database: \n"
            "Current docking mode is vina but last used docking mode of database is adgpu."
        )
        with open(logfile) as f:
            assert warning in f.read()


class TestAD6Handling:
    def test_stream(self, tmp_db):
        from rdkit import Chem

        tmp_db.save_receptor(str(TEST_DATA / "ad6/helix--scofu01.json"))
        suppl = Chem.SDMolSupplier(
            str(TEST_DATA / "ad6/docked_ligands.sdf"), removeHs=False
        )
        tmp_db.add_mol(suppl)
        assert tmp_db.table_length("Results") == 9
        assert tmp_db.table_length("Interactions") == 53

    def test_stream_without_receptor_or_interactions(self, tmp_db):
        from rdkit import Chem

        suppl = Chem.SDMolSupplier(
            str(TEST_DATA / "ad6/docked_ligands.sdf"), removeHs=False
        )
        tmp_db.add_mol(suppl, calculate_interactions=False)
        assert tmp_db.table_length("Results") == 9
        assert tmp_db.table_length("Interactions") == 0

    def test_file_add(self, ad6_db):
        assert ad6_db.table_length("Results") == 9
        assert ad6_db.table_length("Interactions") == 53

    def test_file_add_no_interactions(self, ad6_db_no_interactions):
        assert ad6_db_no_interactions.table_length("Results") == 9
        assert ad6_db_no_interactions.table_length("Interactions") == 0

    def test_gzipped_sdf(self, tmp_db, tmp_path):
        import gzip

        sdf_gz = tmp_path / "docked_ligands.sdf.gz"
        sdf_gz.write_bytes(
            gzip.compress((TEST_DATA / "ad6/docked_ligands.sdf").read_bytes())
        )
        tmp_db.add_results_from_files(
            docking_results=str(sdf_gz),
            docking_mode="ad6",
            calculate_interactions=False,
        )
        assert tmp_db.table_length("Results") == 9

    def test_bad_record_is_logged_and_skipped(self, tmp_db, tmp_path, monkeypatch):
        from rdkit import Chem

        # the failed files log is written to the working directory
        monkeypatch.chdir(tmp_path)
        mols = list(
            Chem.SDMolSupplier(
                str(TEST_DATA / "ad6/docked_ligands.sdf"), removeHs=False
            )
        )
        mols[1].ClearProp("autodock_free_energy")
        sdf = tmp_path / "one_bad.sdf"
        with Chem.SDWriter(str(sdf)) as writer:
            for mol in mols:
                writer.write(mol)
        tmp_db.add_results_from_files(
            docking_results=str(sdf), docking_mode="ad6", calculate_interactions=False
        )
        assert tmp_db.table_length("Results") == 8
        log = (tmp_path / "ringtail_failed_files.log").read_text()
        assert f"{mols[1].GetProp('_Name')} in {sdf}" in log

    def test_calc_interactions_deferred(self, ad6_db_no_interactions, ad6_db):
        """Adding interactions later gives the same interactions as adding them at ingest"""
        db = ad6_db_no_interactions
        db.save_receptor(str(TEST_DATA / "ad6" / "helix--scofu01.json"))
        progress = []
        db.add_interactions(
            chunk_size=2, progress_callback=lambda done, total: progress.append((done, total))
        )
        assert progress[-1] == (9, 9)
        assert db.table_length("Interactions") == 53
        assert TestInteractionAnalysis._from_db(db) == TestInteractionAnalysis._from_db(
            ad6_db
        )

    def test_add_interactions_recalc_larger_cutoffs(self, ad6_db):
        """Recalculating needs consent, backs up the db, and larger cutoffs give more interactions"""
        before = ad6_db.table_length("Interactions")
        hb_before = ad6_db.db_query("SELECT SUM(num_hb) FROM Results")[0][0]

        # without consent nothing changes
        ad6_db.add_interactions(hb_cutoff=6.0, vdw_cutoff=7.0)
        assert ad6_db.table_length("Interactions") == before

        ad6_db.add_interactions(hb_cutoff=6.0, vdw_cutoff=7.0, consent=True, backup=True)
        assert ad6_db.table_length("Results") == 9
        assert ad6_db.table_length("Interactions") > before
        assert ad6_db.db_query("SELECT SUM(num_hb) FROM Results")[0][0] > hb_before

        backup = RingtailCore(ad6_db.db_file + ".bk", storage_type=ad6_db.storagetype)
        assert backup.table_length("Interactions") == before

    def test_recalc_cancel_stays_resumable(self, ad6_db_no_interactions):
        """A cancelled recalculation stops after a committed batch and can be resumed"""
        db = ad6_db_no_interactions
        db.save_receptor(str(TEST_DATA / "ad6" / "helix--scofu01.json"))

        calls = {"n": 0}

        def cancel_after_one_batch():
            calls["n"] += 1
            return calls["n"] > 1

        cancelled = db.add_interactions(
            consent=True, chunk_size=2, should_cancel=cancel_after_one_batch
        )

        assert cancelled["completed"] is False
        assert cancelled["poses_done"] == 2
        # tracking table is kept so the run can be resumed
        assert RECALC_TRACKING_TABLE in db.all_database_tables()
        status = db.interaction_recalc_status()
        assert status["pending"] is True
        assert status["poses_done"] == 2
        assert status["poses_total"] == 9
        assert status["cutoffs"] == (3.7, 4.0)

        # resuming finishes the run with the same interactions
        finished = db.add_interactions(consent=True, chunk_size=2)
        assert finished["completed"] is True
        assert db.table_length("Interactions") == 53
        assert RECALC_TRACKING_TABLE not in db.all_database_tables()
        assert db.interaction_recalc_status()["pending"] is False
        assert duplicate_pairs(db) == 0

    def test_bookmarks_with_interaction_filters(self, ad6_db):
        """Only bookmarks with interaction filters are returned"""
        assert ad6_db.bookmarks_with_interaction_filters() == []

        ad6_db.filter(eworst=-13, output_bookmark="score_only")
        ad6_db.filter(vdw_interactions=[(":VAL::", True)], output_bookmark="uses_vdw")
        ad6_db.filter(hb_count=0, output_bookmark="uses_hb_count")

        found = set(ad6_db.bookmarks_with_interaction_filters())
        assert "uses_vdw" in found
        # hb_count=0 is still an interaction filter
        assert "uses_hb_count" in found
        assert "score_only" not in found


class TestOptions:
    def test_filter_option_checks(self):
        leaf = Filters(score_percentile=20).leaf
        assert leaf.eworst is None
        assert leaf.score_percentile == 20

        # cutoff and percentile on the same column is not allowed
        with pytest.raises(OptionError):
            Filters(score_percentile=20, eworst=-6)
        with pytest.raises(OptionError):
            Filters(le_percentile=20, leworst=-0.4)

        # each alone is fine
        assert Filters(eworst=-6).leaf.eworst == -6
        assert Filters(le_percentile=20).leaf.le_percentile == 20

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


class TestCustomFields:

    def test_custom_fields(self, tmp_db: RingtailCore, make_ad6_sdf):
        """Make sure
        1. SDFs with same columns populates correctly
        2. SDFs with some shared some different columns populate none or value as appropriate
        """
        basic = make_ad6_sdf()
        filler = make_ad6_sdf(
            smiles=" O(C(=O)CC)C(=O)CC",
            ligname="filler_ligand",
            docking_score=-7,
            pose_rank=1,
        )

        custom = make_ad6_sdf(
            smiles="O=C(O)c1ccccc1C(=O)O",
            ligname="custom_ligand",
            docking_score=-13,
            pose_rank=1,
            custom={"supplier": "fischer", "special_score": 2.1},
        )
        custom_same = make_ad6_sdf(
            smiles="O=C(O)c1ccccc1C(=O)O",
            ligname="custom_ligand",
            docking_score=-6.3,
            pose_rank=2,
            custom={"supplier": "sigma-aldrich", "special_score": 8},
        )
        custom_other = make_ad6_sdf(
            smiles="CCO",
            ligname="custom_ligand_other",
            docking_score=-8.9,
            pose_rank=1,
            custom={"supplier": "emolecule", "category": "inactive"},
        )

        tmp_db.add_results_from_files(
            docking_results=[basic, filler, custom, custom_same, custom_other],
            docking_mode="ad6",
            calculate_interactions=False,
        )
        # test that the fields are as specified for each ligand and pose
        rows = tmp_db.db_query("""
            SELECT L.ligname, R.pose_rank,
                C.supplier, C.special_score, C.category
            FROM Results AS R
            JOIN Ligands AS L ON L.ligand_id = R.ligand_id
            LEFT JOIN custom_pose_data AS C ON C.pose_id = R.pose_id
        """)

        actual = {
            (ligname, rank): (supplier, score, category)
            for ligname, rank, supplier, score, category in rows
        }

        assert actual[("test_ligand", 1)] == (None, None, None)
        assert actual[("custom_ligand", 1)] == ("fischer", 2.1, None)
        assert actual[("custom_ligand", 2)] == ("sigma-aldrich", 8.0, None)
        assert actual[("custom_ligand_other", 1)] == ("emolecule", None, "inactive")
        assert actual[("filler_ligand", 1)] == (None, None, None)
        assert len(rows) == 5

    @pytest.mark.parametrize("duplicate_handling", ["ignore", "replace"])
    def test_custom_fields_deduplication(
        self, tmp_db: RingtailCore, make_ad6_sdf, duplicate_handling
    ):
        """
        Make sure custom fields deduplicate as requested
        """
        original = make_ad6_sdf(custom={"supplier": "fischer", "special_score": 3.5})
        replacement = make_ad6_sdf(
            custom={"supplier": "sigma-aldrich", "category": "active"}
        )
        filler = make_ad6_sdf(
            smiles="O=C(O)c1ccccc1C(=O)O",
            ligname="filler_ligand",
            docking_score=-13,
            pose_rank=1,
            custom={"supplier": "fischer", "special_score": 2.1},
        )

        tmp_db.add_results_from_files(
            docking_results=[original, filler],
            docking_mode="ad6",
            calculate_interactions=False,
        )
        test_query = """
            SELECT L.ligname, R.pose_rank, R.run_number,
                R.pose_id, C.supplier, C.special_score
            FROM Results AS R
            JOIN Ligands AS L ON L.ligand_id = R.ligand_id
            JOIN custom_pose_data AS C ON C.pose_id = R.pose_id
        """
        before = {
            (name, rank, run): (pose_id, supplier, score)
            for name, rank, run, pose_id, supplier, score in tmp_db.db_query(test_query)
        }
        target = ("test_ligand", 1, 1)
        filler_key = ("filler_ligand", 1, 1)
        assert len(before) == 2
        assert before[target][1:] == ("fischer", 3.5)

        tmp_db.add_results_from_files(
            docking_results=replacement,
            docking_mode="ad6",
            duplicate_handling=duplicate_handling,
            calculate_interactions=False,
        )
        after = {
            (name, rank, run): (pose_id, supplier, score)
            for name, rank, run, pose_id, supplier, score in tmp_db.db_query(test_query)
        }
        assert tmp_db.table_length("Results") == 2
        assert tmp_db.table_length("custom_pose_data") == 2
        assert after.keys() == before.keys()
        assert after[filler_key] == before[filler_key]

        if duplicate_handling == "ignore":
            assert after == before
        else:
            old_id = before[target][0]
            new_id = after[target][0]

            assert new_id != old_id

            assert after[target][1:] == ("sigma-aldrich", None)

            assert tmp_db.db_query(
                "SELECT category FROM custom_pose_data WHERE pose_id = ?",
                (new_id,),
            ) == [("active",)]

            assert (
                tmp_db.db_query(
                    "SELECT pose_id FROM custom_pose_data WHERE pose_id = ?",
                    (old_id,),
                )
                == []
            )

    def test_custom_fields_export(self, tmp_db, make_ad6_sdf, tmp_path):
        """Test that custom fields are exported"""
        import csv

        custom = make_ad6_sdf(
            ligname="with_custom",
            custom={
                "supplier": "Fischer, Inc.",
                "special score": 3.5,
            },
        )
        basic = make_ad6_sdf(ligname="without_custom")

        tmp_db.add_results_from_files(
            docking_results=[custom, basic],
            docking_mode="ad6",
            calculate_interactions=False,
        )

        csv_path = tmp_path / "custom_fields.csv"
        tmp_db.export_columns_to_csv(
            columns={
                "pose_id": None,
                "ligname": None,
                "supplier": None,
                "special score": "custom_score",
            },
            csv_name=str(csv_path),
        )

        with csv_path.open(newline="") as handle:
            reader = csv.DictReader(handle)
            assert reader.fieldnames == [
                "pose_id",
                "ligname",
                "supplier",
                "custom_score",
            ]
            rows = list(reader)

        assert len(rows) == 2
        by_name = {row["ligname"]: row for row in rows}
        assert set(by_name) == {"with_custom", "without_custom"}

        assert by_name["with_custom"]["supplier"] == "Fischer, Inc."
        assert float(by_name["with_custom"]["custom_score"]) == 3.5

        assert by_name["without_custom"]["supplier"] == ""
        assert by_name["without_custom"]["custom_score"] == ""

        expected_ids = {name: pose_id for name, pose_id in tmp_db.db_query("""
                SELECT L.ligname, R.pose_id
                FROM Results AS R
                JOIN Ligands AS L ON L.ligand_id = R.ligand_id
            """)}
        assert {
            name: int(row["pose_id"]) for name, row in by_name.items()
        } == expected_ids

    def test_custom_fields_subset(self, tmp_db: RingtailCore, make_ad6_sdf, tmp_path):
        """Make sure a subset database keeps the custom fields of the poses it copies"""
        custom = make_ad6_sdf(
            ligname="custom_ligand",
            docking_score=-13,
            custom={"supplier": "fischer", "special_score": 2.1},
        )
        plain = make_ad6_sdf(
            smiles="O=C(O)c1ccccc1C(=O)O", ligname="plain_ligand", docking_score=-12
        )
        left_out = make_ad6_sdf(
            smiles="CCO",
            ligname="left_out_ligand",
            docking_score=-6,
            custom={"supplier": "sigma-aldrich", "special_score": 8.5},
        )
        tmp_db.add_results_from_files(
            docking_results=[custom, plain, left_out],
            docking_mode="ad6",
            calculate_interactions=False,
        )
        tmp_db.filter(eworst=-10, output_bookmark="best")

        subset = tmp_db.export_bookmark_db("best", str(tmp_path / "subset.db"))
        rows = RingtailCore(subset).db_query("""
            SELECT L.ligname, C.supplier, C.special_score
            FROM Results AS R
            JOIN Ligands AS L ON L.ligand_id = R.ligand_id
            LEFT JOIN custom_pose_data AS C ON C.pose_id = R.pose_id
        """)
        assert {ligname: (supplier, score) for ligname, supplier, score in rows} == {
            "custom_ligand": ("fischer", 2.1),
            "plain_ligand": (None, None),
        }

    @pytest.mark.slow  # merges one database into another
    def test_custom_fields_merge(self, tmp_path, storage_type, make_ad6_sdf):
        """Make sure merged poses keep their custom fields under their new pose_ids,
        and custom columns only the merging database has are added"""
        primary = RingtailCore(str(tmp_path / "primary.db"), storage_type=storage_type)
        primary.add_results_from_files(
            docking_results=make_ad6_sdf(
                ligname="primary_ligand",
                custom={"supplier": "fischer", "special_score": 2.1},
            ),
            docking_mode="ad6",
            calculate_interactions=False,
        )
        secondary_db = str(tmp_path / "secondary.db")
        RingtailCore(secondary_db, storage_type=storage_type).add_results_from_files(
            docking_results=[
                make_ad6_sdf(
                    smiles="CCO",
                    ligname="secondary_ligand",
                    custom={"supplier": "sigma-aldrich", "category": "active"},
                ),
                make_ad6_sdf(smiles="O=C(O)c1ccccc1C(=O)O", ligname="secondary_plain"),
            ],
            docking_mode="ad6",
            calculate_interactions=False,
        )

        assert primary.merge_databases([secondary_db], False) == []
        rows = primary.db_query("""
            SELECT L.ligname, C.supplier, C.special_score, C.category
            FROM Results AS R
            JOIN Ligands AS L ON L.ligand_id = R.ligand_id
            LEFT JOIN custom_pose_data AS C ON C.pose_id = R.pose_id
        """)
        assert {
            ligname: (supplier, score, category)
            for ligname, supplier, score, category in rows
        } == {
            "primary_ligand": ("fischer", 2.1, None),
            "secondary_ligand": ("sigma-aldrich", None, "active"),
            "secondary_plain": (None, None, None),
        }
