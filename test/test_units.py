#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail unit tests
#
import json
import logging
from pathlib import Path

import pytest
from ringtail import (
    Filters,
    QueryBuilder,
    RECALC_TRACKING_TABLE,
    RingtailCore,
)
from ringtail.exceptions import OptionError, RTCoreError

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

    def test_molweight_and_max_atoms_together(self, populated_db):
        """Both ligand size filters may be combined through the API, and intersect.
        Contradictory bounds simply return nothing, as with crossing energy bounds.
        The CLI still rejects the combination as a likely mistake."""
        mw, _ = populated_db.filter(ligand_min_molweight=190, output_bookmark="mw_only")
        atoms, _ = populated_db.filter(
            ligand_max_atoms=13, output_bookmark="atoms_only"
        )
        both, _ = populated_db.filter(
            ligand_min_molweight=190,
            ligand_max_atoms=13,
            output_bookmark="mw_and_atoms",
        )
        assert (mw, atoms, both) == (74, 140, 15)

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
        combos = tmp_db._generate_interaction_combinations(filters.leaf.asdict(), 1)
        result_filters = [
            tmp_db._prepare_interaction_combo_filters(filters.leaf.asdict(), c)
            for c in combos
        ]
        assert len(result_filters) == 5
        assert (
            Filters(
                hb_interactions=[("A:ARG:123:", True), ("A:VAL:124:", True)],
                vdw_interactions=[("A:ARG:123:", True)],
            ).leaf.asdict()
            in result_filters
        )
        assert (
            Filters(
                hb_interactions=[("A:ARG:123:", True)],
                vdw_interactions=[("A:ARG:123:", True), ("A:VAL:124:", True)],
            ).leaf.asdict()
            in result_filters
        )

    def test_tiered_filter_or_of_groups(self, populated_db):
        """Nested specification: OR of two AND-groups. The passing ligand set must
        equal the union of each tier filtered on its own (true OR-of-AND semantics)."""
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
        count_or, _ = populated_db.filter(
            filters=expr, output_bookmark="tier_or"
        )
        # each tier on its own
        count_a, _ = populated_db.filter(
            eworst=-6,
            hb_interactions=[("A:VAL:279:", True), ("A:LYS:162:", True)],
            vdw_interactions=[("A:VAL:279:", True), ("A:LYS:162:", True)],
            max_miss=1,
            output_bookmark="tier_a",
        )
        count_b, _ = populated_db.filter(eworst=-7, output_bookmark="tier_b")
        assert count_a == 33  # matches TestFiltering.test_filter

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

    def test_tiered_query_sql_shape(self, tmp_db):
        """Property-only tiers render one parenthesized OR-of-AND WHERE predicate."""
        expr = {
            "op": "or",
            "children": [
                {"eworst": -9, "ebest": -12},
                {"eworst": -7},
            ],
        }
        with tmp_db.storageman as sm:
            sql = sm._generate_filtering_query(Filters.from_dict(expr), "out")
        assert (
            "WHERE ((R.docking_score <= -9 AND R.docking_score >= -12) "
            "OR (R.docking_score <= -7))" in sql
        )

    def test_tiered_query_nested_depth(self, tmp_db):
        """The renderer is fully recursive: arbitrary AND/OR nesting depth works."""
        expr = {
            "op": "and",
            "children": [
                {"eworst": -8},
                {
                    "op": "or",
                    "children": [
                        {"ebest": -12},
                        {"op": "and", "children": [{"eworst": -9}, {"lebest": -0.5}]},
                    ],
                },
            ],
        }
        with tmp_db.storageman as sm:
            sql = sm._generate_filtering_query(Filters.from_dict(expr), "out")
        assert (
            "((R.docking_score <= -8) AND ((R.docking_score >= -12) "
            "OR ((R.docking_score <= -9) AND (R.leff >= -0.5))))" in sql
        )

    def test_smarts_inside_group(self, populated_db):
        """A SMARTS criterion works inside a filter group (RDKit leaf ->
        pose_id IN (...)), ANDing with SQL criteria in the same group."""
        score_only = {"op": "and", "children": [{"eworst": -6}]}
        with_smarts = {
            "op": "and",
            "children": [{"eworst": -6}, {"ligand_substruct": ["C=O"]}],
        }
        c_score, _ = populated_db.filter(
            filters=score_only, output_bookmark="sig_score"
        )
        c_both, _ = populated_db.filter(
            filters=with_smarts, output_bookmark="sig_both"
        )
        assert c_both > 0
        assert c_both <= c_score  # adding SMARTS only narrows

    def test_expr_rdkit_group_count(self):
        """Cross-group SMARTS detection (drives the API warning / GUI consent dialog)."""
        one = {"op": "and", "children": [{"eworst": -6}, {"ligand_substruct": ["C=O"]}]}
        two = {
            "op": "or",
            "children": [
                {"op": "and", "children": [{"ligand_substruct": ["C=O"]}]},
                {
                    "op": "and",
                    "children": [{"ligand_substruct": ["CN"]}, {"eworst": -7}],
                },
            ],
        }
        assert Filters.from_dict(one).rdkit_group_count() == 1
        assert Filters.from_dict(two).rdkit_group_count() == 2

    def test_tiered_with_global_smarts(self, populated_db):
        """A SMARTS criterion applying to a whole expression is its own group ANDed on."""
        expr = {"op": "or", "children": [{"eworst": -6}]}
        count_expr, _ = populated_db.filter(
            filters=expr, output_bookmark="tier_smarts_a"
        )
        count_both, _ = populated_db.filter(
            filters={
                "op": "and",
                "children": [expr, {"ligand_substruct": ["C=O"]}],
            },
            output_bookmark="tier_smarts_b",
        )
        assert count_both > 0
        assert count_both <= count_expr  # SMARTS only narrows the tree result

    def test_filters_and_flat_args_are_exclusive(self, tmp_db):
        """Passing both a specification and flat criteria is a mistake, not a merge."""
        with pytest.raises(OptionError):
            tmp_db.filter(filters={"eworst": -6}, ligand_substruct=["C=O"])


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

    def test_export_bookmark_db_carries_statuses_and_comments(
        self, populated_db, tmp_path
    ):
        """Statuses and comments are hand curation from visual inspection, so losing
        them on export silently throws away the most expensive data in the database.
        Pose_comments is the sneaky half: it is created by ensure_gui_tables rather than
        _create_tables, so the subset had no comments table at all to copy into."""
        populated_db.filter(eworst=-7, output_bookmark="curated")
        with populated_db.storageman as sm:
            sm.ensure_gui_tables()
            kept = [
                row[0]
                for row in sm.db_query(
                    "SELECT pose_id FROM Filtered_poses FP "
                    "JOIN Filters F ON F.filter_id = FP.filter_id "
                    "WHERE F.name = 'curated' ORDER BY pose_id"
                ).fetchall()
            ]
        assert len(kept) >= 3, "need a few poses to curate"
        accepted, maybe, rejected = kept[0], kept[1], kept[2]

        populated_db.update_pose_status(accepted, 1)
        populated_db.update_pose_status(maybe, 2)
        populated_db.update_pose_status(rejected, 3)
        populated_db.set_pose_comment(accepted, "clear H-bond, keep")

        subset_path = populated_db.export_bookmark_db(
            "curated", str(tmp_path / "curated.db")
        )
        subset = RingtailCore(subset_path)

        with subset.storageman as sm:
            for table, expected in (
                ("Accepted", accepted),
                ("Maybe", maybe),
                ("Rejected", rejected),
            ):
                rows = sm.db_query(f"SELECT pose_id FROM {table}").fetchall()
                assert [r[0] for r in rows] == [expected], f"{table} did not transfer"
        assert subset.get_pose_comment(accepted) == "clear H-bond, keep"

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


class TestInteractionAnalysis:
    """Interactions must pair the right atoms, not merely the right number of them.

    Comparing the two receptor formats does that without a fixed expectation: the
    Polymer JSON and pdbqt paths read different files through different code and must
    arrive at the same interactions, atom for atom. Counts alone cannot show that.

    An earlier version of this class compared against the ANALYSIS: blocks AutoDock-GPU
    writes into its own .dlg files. That was backed out. Reaching agreement required
    adopting AutoDock-GPU's receptor conventions wholesale -- H-bonds reported on the
    donor's heavy atom, every H-bond-capable atom excluded from van der Waals -- and
    Ringtail is not AutoDock-GPU. The implementations are allowed to differ.
    """

    @staticmethod
    def _from_db(rtc: RingtailCore) -> dict:
        """{(ligname, run_number): {(type, residue, resid, chain, recname), ...}}

        Interactions are held as sets: a receptor atom reached by several ligand atoms
        is stored once, so the set is what the database actually holds.
        """
        rows = rtc.db_query(
            """SELECT L.ligname, R.run_number, II.interaction_type,
                      II.rec_resname, II.rec_resid, II.rec_chain, II.rec_atom
               FROM Interactions I
               JOIN Interaction_indices II ON II.interaction_id = I.interaction_id
               JOIN Results R ON R.pose_id = I.pose_id
               JOIN Ligands L ON L.ligand_id = R.ligand_id"""
        )
        poses = {}
        for ligname, run, itype, resname, resid, chain, recname in rows:
            poses.setdefault((ligname, run), set()).add(
                (itype, resname, str(resid), chain, recname)
            )
        return poses

    def test_receptor_formats_agree(self, flexres_json_db, flexres_pdbqt_db):
        """A Polymer JSON receptor and a pdbqt one must give identical interactions."""
        json_poses = self._from_db(flexres_json_db)
        pdbqt_poses = self._from_db(flexres_pdbqt_db)
        assert json_poses, "the flexres fixture stored no interactions to compare"
        assert set(json_poses) == set(pdbqt_poses), (
            "the two receptor formats resolved different poses: "
            f"json={sorted(json_poses)} pdbqt={sorted(pdbqt_poses)}"
        )
        for pose in sorted(json_poses):
            assert json_poses[pose] == pdbqt_poses[pose], (
                f"{pose[0]} run {pose[1]} differs between receptor formats:\n"
                f"  only from json : {sorted(json_poses[pose] - pdbqt_poses[pose])}\n"
                f"  only from pdbqt: {sorted(pdbqt_poses[pose] - json_poses[pose])}"
            )


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

        # a dependent's poses are its own, so deleting the parent leaves its data
        # intact and only orphans the recorded lineage
        child_poses = populated_db.table_length("child_bm")
        assert child_poses > 0
        populated_db.delete_bookmark("parent_bm")

        assert "parent_bm" not in populated_db.get_bookmark_names()
        assert "child_bm" in populated_db.get_bookmark_names()
        assert populated_db.table_length("child_bm") == child_poses

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

    def test_merge_several_databases_in_one_call(self, tmp_path, storage_type):
        """Each secondary must be detached before the next one is attached."""
        dbs = []
        for name, ligand in [
            ("primary", "1451"),
            ("secondary", "1620"),
            ("tertiary", "1751"),
            ("quaternary", "10091"),
        ]:
            db = str(tmp_path / f"{name}.db")
            rtc = RingtailCore(db, storage_type=storage_type)
            rtc.add_results_from_files(
                str(TEST_DATA / f"adgpu/group1/{ligand}.dlg.gz"), docking_mode="adgpu"
            )
            assert rtc.table_length("Ligands") == 1
            dbs.append(db)

        rtc1 = RingtailCore(dbs[0])
        assert rtc1.merge_databases(dbs[1:], False) == []
        assert rtc1.table_length("Ligands") == 4


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
    def test_file_add(self, vina_db):
        assert vina_db.table_length("Results") == 6

    def test_string_add(self, tmp_db):
        vina_path = TEST_DATA / "vina"
        sample1 = (vina_path / "sample-result.pdbqt").read_text()
        sample2 = (vina_path / "sample-result-2.pdbqt").read_text()
        tmp_db.save_receptor(str(vina_path / "receptor.pdbqt"))
        tmp_db.add_results_from_vina_string(
            results={"sample1": sample1, "sample2": sample2}
        )
        assert tmp_db.table_length("Results") == 6

    def test_add_interactions(self, vina_db):
        assert vina_db.table_length("Interaction_indices") == 25
        assert vina_db.table_length("Interactions") == 60

    def test_add_interactions_from_polymer(self, flexres_json_db, flexres_pdbqt_db):
        """The two receptor formats must produce identical interactions.

        The count is pinned here; that the two agree atom for atom is checked in
        TestInteractionAnalysis::test_receptor_formats_agree.
        """
        assert (
            flexres_json_db.table_length("Ligands")
            == flexres_pdbqt_db.table_length("Ligands")
            == 1
        )
        assert (
            flexres_json_db.table_length("Interactions")
            == flexres_pdbqt_db.table_length("Interactions")
            == 25
        )

    def test_polymer_receptor(self, flexres_json_db):
        receptor_items = flexres_json_db.get_receptor_object()
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

    def test_various_filters(self, vina_db):
        count, _ = vina_db.filter(eworst=-6, ligand_substruct=["[N]"])
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
        assert tmp_db.table_length("Interactions") == 53

    def test_file_add(self, ad6_db):
        assert ad6_db.table_length("Results") == 9
        assert ad6_db.table_length("Interactions") == 53

    def test_file_add_no_interactions(self, ad6_db_no_interactions):
        assert ad6_db_no_interactions.table_length("Results") == 9
        assert ad6_db_no_interactions.table_length("Interactions") == 0

    def test_calc_interactions_deferred(self, ad6_db_no_interactions, ad6_db):
        """Calculating later must land on the same interactions as calculating at ingest.

        Compared against the ingest-time fixture rather than a hardcoded count, so the
        test keeps checking the two paths agree even when perception changes.
        """
        db = ad6_db_no_interactions
        assert db.table_length("Interactions") == 0

        db.save_receptor(str(TEST_DATA / "ad6" / "helix--scofu01.json"))
        db.add_interactions()
        assert db.table_length("Results") == 9
        assert db.table_length("Interactions") == 53
        # and the same interactions, not merely the same number of them
        assert TestInteractionAnalysis._from_db(db) == TestInteractionAnalysis._from_db(
            ad6_db
        )

    def test_recalc_in_batches_matches_one_pass(self, ad6_db, ad6_db_no_interactions):
        """chunk_size smaller than the pose count must give the same answer.

        The default chunk_size (500) exceeds any test dataset, so the batching path
        used to go unexercised — and it was broken: the pending-poses query was a live
        cursor over the tracking table that the batch commit wrote to, which on duckdb
        replaced the rows still being read.
        """
        one_pass = ad6_db
        one_pass.add_interactions(consent=True, chunk_size=500)  # no mid-run flush
        batched = ad6_db_no_interactions
        batched.save_receptor(str(TEST_DATA / "ad6" / "helix--scofu01.json"))
        batched.add_interactions(consent=True, chunk_size=2)  # 9 poses -> 5 batches

        assert batched.table_length("Interactions") == one_pass.table_length(
            "Interactions"
        )
        query = """SELECT I.pose_id, II.interaction_type, II.rec_chain, II.rec_resid,
                          II.rec_atom
                   FROM Interactions I
                   JOIN Interaction_indices II ON II.interaction_id=I.interaction_id"""
        assert set(batched.db_query(query)) == set(one_pass.db_query(query))

    def test_recalc_resumes_after_interruption(self, ad6_db_no_interactions):
        """An interrupted run picks up from the tracking table, without duplicating."""
        import ringtail.ringtailcore as rtc_module

        db = ad6_db_no_interactions
        db.save_receptor(str(TEST_DATA / "ad6" / "helix--scofu01.json"))

        real = rtc_module.find_interactions
        calls = {"n": 0}

        def fail_partway(*args, **kwargs):
            calls["n"] += 1
            if calls["n"] > 4:
                raise KeyboardInterrupt("simulated interruption")
            return real(*args, **kwargs)

        rtc_module.find_interactions = fail_partway
        try:
            with pytest.raises(BaseException):
                db.add_interactions(consent=True, chunk_size=2)
        finally:
            rtc_module.find_interactions = real

        partial = db.table_length("Interactions")
        assert partial > 0, "some work should have been committed before the failure"

        db.add_interactions(consent=True, chunk_size=2)
        assert db.table_length("Results") == 9
        # every pose accounted for, and none of the committed work redone
        assert (
            db.db_query("SELECT COUNT(DISTINCT pose_id) FROM Interactions")[0][0] == 9
        )
        assert duplicate_pairs(db) == 0

    def test_recalc_refuses_to_resume_with_different_cutoffs(
        self, ad6_db_no_interactions
    ):
        """Resuming at new cutoffs would mix two calculations with no record of which."""
        import ringtail.ringtailcore as rtc_module

        db = ad6_db_no_interactions
        db.save_receptor(str(TEST_DATA / "ad6" / "helix--scofu01.json"))

        real = rtc_module.find_interactions
        calls = {"n": 0}

        def fail_partway(*args, **kwargs):
            calls["n"] += 1
            if calls["n"] > 4:
                raise KeyboardInterrupt("simulated interruption")
            return real(*args, **kwargs)

        rtc_module.find_interactions = fail_partway
        try:
            with pytest.raises(BaseException):
                db.add_interactions(consent=True, chunk_size=2)
        finally:
            rtc_module.find_interactions = real

        with pytest.raises(OptionError, match="cutoffs"):
            db.add_interactions(hb_cutoff=6.0, vdw_cutoff=7.0, consent=True)

        # the original cutoffs still resume, so the guard is not a dead end
        db.add_interactions(consent=True, chunk_size=2)
        assert (
            db.db_query("SELECT COUNT(DISTINCT pose_id) FROM Interactions")[0][0] == 9
        )

    def test_interactions_land_on_their_own_pose(self, ad6_db):
        """Each pose keeps its own interactions.

        Interactions resolve to a pose through (ligname, run_number, pose_rank). The sdf
        path reported a fixed pose_rank of 1, so every pose of a ligand resolved to that
        ligand's first pose: rank 1 accumulated everyone's interactions, as duplicates,
        and the other poses were stored with none.
        """
        rows = ad6_db.db_query(
            """SELECT R.pose_id, R.pose_rank,
                      (SELECT COUNT(*) FROM Interactions I
                       WHERE I.pose_id = R.pose_id) AS n
               FROM Results R"""
        )
        assert rows
        assert all(n > 0 for _, _, n in rows), (
            f"poses stored with no interactions: {[p for p, _, n in rows if n == 0]}"
        )
        assert len({rank for _, rank, _ in rows}) > 1, "need >1 rank to be meaningful"

        assert duplicate_pairs(ad6_db) == 0

    def test_add_interactions_recalc_larger_cutoffs(self, ad6_db):
        # the fixture is populated at the default cutoffs (3.7 HB, 4.0 VDW)
        tmp_db = ad6_db
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

    def test_recalc_reports_progress(self, ad6_db_no_interactions):
        """Progress is reported per committed batch and ends at the pose count."""
        db = ad6_db_no_interactions
        db.save_receptor(str(TEST_DATA / "ad6" / "helix--scofu01.json"))

        seen = []
        result = db.add_interactions(
            consent=True, chunk_size=2, progress_callback=lambda d, t: seen.append((d, t))
        )

        assert seen, "a batched run should report at least once"
        assert all(total == 9 for _, total in seen)
        assert [done for done, _ in seen] == sorted(done for done, _ in seen)
        assert seen[-1][0] == 9
        assert result == {"completed": True, "poses_done": 9, "poses_total": 9}

    def test_recalc_progress_resumes_from_where_it_stopped(self, ad6_db_no_interactions):
        """A resumed run counts the whole database, not just what is left.

        Reporting only the remaining poses would send a progress bar back to zero on
        every resume, which reads as the work being thrown away.
        """
        db = ad6_db_no_interactions
        db.save_receptor(str(TEST_DATA / "ad6" / "helix--scofu01.json"))

        stop_after = {"n": 0}

        def cancel_after_two_batches():
            stop_after["n"] += 1
            return stop_after["n"] > 2

        first = db.add_interactions(
            consent=True, chunk_size=2, should_cancel=cancel_after_two_batches
        )
        assert not first["completed"]
        assert 0 < first["poses_done"] < 9

        seen = []
        second = db.add_interactions(
            consent=True, chunk_size=2, progress_callback=lambda d, t: seen.append(d)
        )
        assert second["completed"]
        assert seen[0] > first["poses_done"], "resumed progress must not restart at zero"
        assert seen[-1] == 9

    def test_recalc_cancel_stays_resumable(self, ad6_db_no_interactions, ad6_db):
        """Cancelling stops on a committed boundary and keeps the tracking table."""
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
        # the tracking table is what makes the run resumable, so it must survive
        assert RECALC_TRACKING_TABLE in db.all_database_tables()
        status = db.interaction_recalc_status()
        assert status["pending"] is True
        assert status["poses_done"] == 2
        assert status["poses_total"] == 9
        assert status["cutoffs"] == (3.7, 4.0)

        # finishing gives the same answer as never having been interrupted
        finished = db.add_interactions(consent=True, chunk_size=2)
        assert finished["completed"] is True
        assert db.table_length("Interactions") == 53
        assert RECALC_TRACKING_TABLE not in db.all_database_tables()
        assert db.interaction_recalc_status()["pending"] is False
        assert duplicate_pairs(db) == 0

    def test_recalc_resumes_without_consent(self, ad6_db_no_interactions):
        """Finishing an interrupted run needs no consent, because it deletes nothing.

        Consent exists to guard the delete-and-recompute. A resume only computes the
        poses that were never reached, so requiring it there made the plain
        ``add_interactions()`` call the docs show return a silent no-op -- reporting
        poses_done: 0 and logging "Consent not given for deleting and re-calculating
        interactions" about a call that would not have deleted anything.
        """
        db = ad6_db_no_interactions
        db.save_receptor(str(TEST_DATA / "ad6" / "helix--scofu01.json"))

        calls = {"n": 0}

        def cancel_after_one_batch():
            calls["n"] += 1
            return calls["n"] > 1

        first = db.add_interactions(chunk_size=2, should_cancel=cancel_after_one_batch)
        assert first["completed"] is False
        assert db.table_length("Interactions") > 0, "partial work must be present"

        second = db.add_interactions(chunk_size=2)
        assert second == {"completed": True, "poses_done": 9, "poses_total": 9}
        assert db.table_length("Interactions") == 53

        # with the tracking table gone, consent is required again
        assert db.add_interactions() == {
            "completed": False,
            "poses_done": 0,
            "poses_total": 0,
        }
        assert db.table_length("Interactions") == 53

    def test_recalc_without_a_receptor_raises_before_deleting(self, ad6_db):
        """No receptor must stop the run, not silently empty the database.

        make_interaction_finder reports failure by returning None, and find_interactions
        then reports zero interactions for every pose. Checked after the tables were
        cleared, that deleted every interaction, zeroed num_hb/num_interactions and
        reported completed: True, with only a logged warning to say otherwise.
        """
        before = TestInteractionAnalysis._from_db(ad6_db)
        counts_before = ad6_db.db_query(
            "SELECT pose_id, num_interactions, num_hb FROM Results ORDER BY pose_id"
        )
        assert before

        ad6_db.db_query("DELETE FROM Receptors", commit=True)
        with pytest.raises(RTCoreError, match="no receptor"):
            ad6_db.add_interactions(consent=True)

        assert TestInteractionAnalysis._from_db(ad6_db) == before
        assert (
            ad6_db.db_query(
                "SELECT pose_id, num_interactions, num_hb FROM Results ORDER BY pose_id"
            )
            == counts_before
        )
        assert RECALC_TRACKING_TABLE not in ad6_db.all_database_tables()

    def test_recalc_with_an_unreadable_receptor_raises_before_deleting(self, ad6_db):
        """Same guard, for a receptor that is stored but cannot be parsed."""
        before = TestInteractionAnalysis._from_db(ad6_db)
        assert before

        # valid JSON, so the column accepts it on duckdb, but not a Polymer
        ad6_db.db_query(
            "UPDATE Receptors SET receptor_object=NULL, polymer=?",
            ['{"not_a_polymer": 1}'],
            commit=True,
        )
        with pytest.raises(RTCoreError, match="could not be read"):
            ad6_db.add_interactions(consent=True)

        assert TestInteractionAnalysis._from_db(ad6_db) == before
        assert RECALC_TRACKING_TABLE not in ad6_db.all_database_tables()

    def test_recalc_groups_poses_by_ligand(self, ad6_db, ad6_db_no_interactions):
        """Grouping a ligand's poses into one call must not change the answer.

        find_interactions is handed a whole ligand at a time so meeko prepares the
        molsetup once instead of once per pose. That preparation produces the atom-type
        list every pose is indexed against, so a grouping mistake reintroduces exactly
        the index-shift bug this method exists to repair — silently, since it would
        still produce plausible-looking interactions.
        """
        import ringtail.ringtailcore as rtc_module

        one_pass = ad6_db
        one_pass.add_interactions(consent=True, chunk_size=500)

        grouped = ad6_db_no_interactions
        grouped.save_receptor(str(TEST_DATA / "ad6" / "helix--scofu01.json"))

        real = rtc_module.find_interactions
        poses_per_call = []

        def counting_find_interactions(poses_coordinates, *args, **kwargs):
            poses_per_call.append(len(poses_coordinates))
            return real(poses_coordinates, *args, **kwargs)

        rtc_module.find_interactions = counting_find_interactions
        try:
            grouped.add_interactions(consent=True, chunk_size=500)
        finally:
            rtc_module.find_interactions = real

        assert max(poses_per_call) > 1, (
            "no ligand's poses were batched together, so the grouping is not in effect"
        )
        assert sum(poses_per_call) == 9

        query = """SELECT I.pose_id, II.interaction_type, II.rec_chain, II.rec_resid,
                          II.rec_atom
                   FROM Interactions I
                   JOIN Interaction_indices II ON II.interaction_id=I.interaction_id"""
        assert set(grouped.db_query(query)) == set(one_pass.db_query(query))
        counts = "SELECT pose_id, num_interactions, num_hb FROM Results ORDER BY pose_id"
        assert grouped.db_query(counts) == one_pass.db_query(counts)

    def test_recalc_backup_leaves_the_original_untouched(self, ad6_db):
        """The backup is taken before anything is deleted, not after."""
        before = ad6_db.table_length("Interactions")
        assert before > 0

        ad6_db.add_interactions(
            hb_cutoff=6.0, vdw_cutoff=7.0, consent=True, backup=True
        )

        backup_file = Path(ad6_db.db_file + ".bk")
        assert backup_file.is_file()
        # a backup taken after clear_interaction_tables would be a copy of the damage
        backup_db = RingtailCore(str(backup_file), storage_type=ad6_db.storagetype)
        assert backup_db.table_length("Interactions") == before
        assert ad6_db.table_length("Interactions") > before

    def test_filtering(self, ad6_db):
        count, _ = ad6_db.filter(
            eworst=-13, ligand_substruct=["C=O"], vdw_interactions=[(":VAL::", True)]
        )
        assert count == 1

    def test_bookmarks_with_interaction_filters(self, ad6_db):
        """Only bookmarks whose filters touch interactions are reported.

        These are the ones a recalculation invalidates the meaning of, so the caller
        can warn about them by name before deleting anything.
        """
        assert ad6_db.bookmarks_with_interaction_filters() == []

        ad6_db.filter(eworst=-13, output_bookmark="score_only")
        ad6_db.filter(
            vdw_interactions=[(":VAL::", True)], output_bookmark="uses_vdw"
        )
        ad6_db.filter(hb_count=0, output_bookmark="uses_hb_count")

        found = set(ad6_db.bookmarks_with_interaction_filters())
        assert "uses_vdw" in found
        # 0 is a real criterion for hb_count ("no hydrogen bonds"), so a truthiness
        # test here would drop this bookmark
        assert "uses_hb_count" in found
        assert "score_only" not in found

    def test_bookmarks_with_interaction_filters_sees_nested_filters(self, populated_db):
        """A tiered filter stores its criteria nested, so the scan must recurse."""
        populated_db.filter(
            filters={
                "op": "or",
                "children": [
                    {"eworst": -7},
                    {"hb_interactions": [("A:VAL:279:", True)]},
                ],
            },
            output_bookmark="nested_tier",
        )
        assert "nested_tier" in populated_db.bookmarks_with_interaction_filters()


class TestLogger:
    def test_set_log_level(self):
        from ringtail import LOGGER, setup_logging

        setup_logging(level="INFO")
        assert LOGGER.level == logging.INFO


class TestOptions:
    def test_filter_option_checks(self, tmp_db, tmp_path):
        # No results are ingested here on purpose: Filter.checks() runs at construction and
        # touches no data, so building a database would only make this test slow.

        # criteria live on the leaf; a flat specification has exactly one
        leaf = Filters(score_percentile=20).leaf
        assert leaf.eworst is None
        assert leaf.score_percentile == 20

        # An absolute cutoff and a percentile on the same column are two different
        # requests, so they are refused rather than one silently overriding the other.
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



