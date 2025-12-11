#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail unit testing
#
from ringtail import RingtailCore, Filters, QueryBuilder
import os
import json
import pytest


def _create_test_db(db_name: str = "output.db"):
    rtc = RingtailCore(db_file=db_name)
    rtc.add_results_from_files(file_path="test_data/adgpu/group1/")
    rtc.add_results_from_files(file_path="test_data/adgpu/group2/")


def _db_exists(db_name: str = "output.db"):
    return os.path.exists(db_name)


class TestRingtailCore:
    def test_add_file(self):
        os.system("rm output.db*")
        rtc = RingtailCore(db_file="output.db")
        rtc.add_results_from_files(
            file="test_data/adgpu/group1/1451.dlg.gz", max_poses=3
        )
        count = rtc.table_length("Results")
        assert count == 3

    def test_storeallposes(self):
        rtc = RingtailCore(db_file="output.db")
        rtc.add_results_from_files(
            file="test_data/adgpu/group1/1451.dlg.gz",
            store_all_poses=True,
            overwrite=True,
        )
        count = rtc.table_length("Results")
        assert count == 20

    def test_add_folder(self):
        rtc = RingtailCore(db_file="output.db")
        rtc.add_results_from_files(file_path="test_data/adgpu/group1", overwrite=True)
        count = rtc.table_length("Ligands")
        assert count == 138

    def test_save_receptor(self):
        rtc = RingtailCore(db_file="output.db", logging_level="DEBUG")
        count0 = rtc.db_query(
            "SELECT COUNT(*) FROM Receptors WHERE receptor_object NOT NULL"
        )[0][0]

        assert count0 == 0

        rtc.save_receptor(receptor="test_data/adgpu/4j8m.pdbqt")
        count = rtc.db_query(
            "SELECT COUNT(*) FROM Receptors WHERE receptor_object NOT NULL"
        )[0][0]

        assert count == 1

    def test_produce_summary(self):
        # Ensure storage error thrown if no data in database
        from ringtail import exceptions as e

        with pytest.raises(e.StorageError):
            fake_rtc = RingtailCore("nodata.db")
            fake_rtc.produce_summary()
        os.system("rm nodata.db")

        import sys

        class ListStream:
            def __init__(self):
                self.data = []

            def write(self, s):
                self.data.append(s)

        sys.stdout = summary_items = ListStream()
        rtc = RingtailCore(db_file="output.db")
        rtc.produce_summary()
        sys.stdout = sys.__stdout__

        assert len(summary_items.data) == 38

    def test_append_to_database(self):
        rtc = RingtailCore(db_file="output.db")
        rtc.add_results_from_files(file_path="test_data/adgpu/group2/")
        count = rtc.table_length("Ligands")

        assert count == 217

    def test_filter(self):
        if not _db_exists():
            _create_test_db()
        rtc = RingtailCore(db_file="output.db")
        count_ligands_passing = rtc.filter(
            eworst=-6,
            hb_interactions=[("A:VAL:279:", True), ("A:LYS:162:", True)],
            vdw_interactions=[("A:VAL:279:", True), ("A:LYS:162:", True)],
            max_miss=1,
            bookmark_name="union_bookmark",
        )
        # make sure correct number of ligands passing
        assert count_ligands_passing == 33
        # make sure only one bookmark was created
        bookmarks = rtc.get_bookmark_names()
        assert len(bookmarks) == 1
        assert bookmarks[0] == "union_bookmark"
        rtc.delete_bookmark("union_bookmark")

    def test_return_iter(self):
        if not _db_exists():
            _create_test_db()
        rtc = RingtailCore(db_file="output.db")
        iterable = rtc.filter(eworst=-7, bookmark_name="iterable", return_iter=True)

        assert len(iterable) == 8

    def test_enumerate_interaction_combinations(self):
        if not _db_exists():
            _create_test_db()
        # first, test without enumerate, check number of passing union as well as number of bookmarks
        rtc = RingtailCore(db_file="output.db")
        # get current bookmark count
        bookmarks_old = rtc.get_bookmark_names()
        count_ligands_passing = rtc.filter(
            eworst=-6,
            hb_interactions=[("A:VAL:279:", True), ("A:LYS:162:", True)],
            vdw_interactions=[("A:VAL:279:", True), ("A:LYS:162:", True)],
            max_miss=1,
            enumerate_interaction_combs=True,
            bookmark_name="enumerated_bookmark",
        )
        # make sure correct number of ligands passing
        assert count_ligands_passing == 33

        # make sure additional bookmarks were created for the enumerated combinations
        bookmarks = rtc.get_bookmark_names()
        bookmarks = [element for element in bookmarks if element not in bookmarks_old]
        # This filtering session should produce 6 bookmarks
        assert len(bookmarks) == 6

        # check that naming works properly
        assert "enumerated_bookmark_0" in bookmarks
        assert "enumerated_bookmark_union" in bookmarks

    def test_filter_from_bookmark(self):
        if not _db_exists():
            _create_test_db()
        rtc = RingtailCore(db_file="output.db")
        count_passing_ligands1 = rtc.filter(eworst=-6, bookmark_name="filter_window")
        count_passing_ligands2 = rtc.filter(
            eworst=-7, bookmark_name="bookmark", filter_bookmark="filter_window"
        )
        assert count_passing_ligands1 > count_passing_ligands2

    def test_ligand_filters(self):
        if not _db_exists():
            _create_test_db()
        rtc = RingtailCore(db_file="output.db")

        # tests for partial names
        count_ligname = rtc.filter(ligand_name=["88"], bookmark_name="ligname")
        assert count_ligname == 7

        # test substructure search (default 'OR' ligand_operator)
        count_substruct_or = rtc.filter(
            ligand_substruct=["C=O", "CC(C)(C)"], bookmark_name="substruct_or"
        )
        assert count_substruct_or == 90

        # test substructure search ('AND' ligand_operator)
        count_substruct_and = rtc.filter(
            ligand_substruct=["C=O", "CC(C)(C)"],
            ligand_operator="AND",
            bookmark_name="substruct_and",
        )
        assert count_substruct_and == 18

        count_substruct_pos = rtc.filter(
            ligand_substruct_pos=[
                ["[C][Oh]", 1, 10, 102, 106, 154],
                ["C=O", 1, 10, 102, 106, 154],
            ],
            bookmark_name="substruct_pos",
        )
        assert count_substruct_pos == 12

    def test_all_filters(self):
        if not _db_exists():
            _create_test_db()
        rtc = RingtailCore(db_file="output.db")
        count_ligands_passing = rtc.filter(
            eworst=-6,
            hb_interactions=[("A:VAL:279:", True), ("A:LYS:162:", True)],
            vdw_interactions=[("A:VAL:279:", True), ("A:LYS:162:", True)],
            max_miss=1,
            bookmark_name="big_query",
            ligand_name=["88"],
        )

        assert count_ligands_passing == 1

    def test_get_filterdata(self):
        if not _db_exists():
            _create_test_db()
        rtc = RingtailCore(db_file="output.db")
        rtc.filter(eworst=-7, bookmark_name="has_filterdata")
        log_file_name = "output_log_test.txt"
        rtc.get_previous_filter_data(
            "has_filterdata", "delta, reference_rmsd", log_file=log_file_name
        )

        with open(log_file_name) as f:
            file_contents = f.read()
        import linecache

        final_line = linecache.getline(log_file_name, 11)

        assert "11991, 0.0, 226.06" in file_contents
        assert "3961, -0.02, 215.96" in file_contents
        assert final_line == "***************\n"

        os.system(("rm " + log_file_name))

    def test_similar_ligands_interaction(self, monkeypatch):

        rtc = RingtailCore(db_file="ifp_cluster.db")
        rtc.add_results_from_files(
            file_path=["test_data/adgpu/group1/", "test_data/adgpu/group2/"],
        )
        ligand_name = "28837"
        rtc.filter(ebest=-6, interaction_cluster=0.5)
        monkeypatch.setattr("builtins.input", lambda _: 1)  # provides terminal input
        number_similar = rtc.find_similar_ligands(ligand_name)
        os.system("rm ifp_cluster.db")

        assert number_similar == 13

    def test_similar_ligands_mfpt(self, monkeypatch):
        rtc = RingtailCore(db_file="mfpt_cluster.db")
        rtc.add_results_from_files(
            file_path=["test_data/adgpu/group1/", "test_data/adgpu/group2/"],
        )
        ligand_name = "287065"
        rtc.filter(ebest=-6, mfpt_cluster=0.5)
        monkeypatch.setattr("builtins.input", lambda _: 1)  # provides terminal input
        number_similar = rtc.find_similar_ligands(ligand_name)
        os.system("rm mfpt_cluster.db")
        assert number_similar == 8

    def test_create_rdkitmol(self):
        if not _db_exists():
            _create_test_db()
        bookmark_name = "rdkit_test"
        rtc = RingtailCore(db_file="output.db")
        ligname = "14303"
        rtc.filter(ebest=-3, bookmark_name=bookmark_name)
        ligands_poses = rtc._fetch_select_ligands_poses(
            ligand_names=[ligname], bookmark_name=bookmark_name
        )

        mol = rtc.create_rdkit_mol(ligname, ligands_poses[ligname])[0]
        # grab one molecule from bookmark and check number of atoms
        num_of_atoms = mol.GetNumAtoms()
        assert num_of_atoms == 10

    def test_write_sdfs(self):
        if not _db_exists():
            _create_test_db()
        sdf_path = "sdf_files"
        rtc = RingtailCore(db_file="output.db")
        rtc.filter(eworst=-7, bookmark_name="sdf_bookmark")
        rtc.write_molecule_sdfs("sdf_bookmark", sdf_path, all_in_one=False)

        # ensure correct number of files written
        sdf_files = os.listdir(sdf_path)
        expected = [
            "3961.sdf",
            "5995.sdf",
            "11128.sdf",
            "11991.sdf",
            "13974.sdf",
            "15776.sdf",
            "136065.sdf",
            "127947.sdf",
        ]
        assert len(sdf_files) == len(expected)

        # ensure contents is correct
        with open("sdf_files/136065.sdf") as sdf:
            sdf.readline()
            sdf.readline()
            sdf.readline()
            fourth_line = sdf.readline()
        assert fourth_line == " 27 28  0  0  0  0  0  0  0  0999 V2000\n"

        # ensure the correct files were written
        for f in sdf_files:
            assert f in expected
            os.remove(sdf_path + "/" + f)
        os.rmdir(sdf_path)

    def test_pymol(self):
        # will not add a test for now, as I cannot figure out an unambiguous, lightweight way to test
        pass

    def test_export_csv(self):
        if not _db_exists():
            _create_test_db()
        rtc = RingtailCore(db_file="output.db")
        rtc.filter(eworst=-7, log_file="different_log.txt", bookmark_name="export_csv")
        rtc.export_table_as_csv("Ligands", "Ligands.csv")

        assert os.path.exists("Ligands.csv")
        os.system("rm Ligands.csv")
        rtc.export_table_as_csv("export_csv", "export_csv.csv")
        assert os.path.exists("export_csv.csv")
        os.system("rm export_csv.csv")

    def export_receptor_pdbqt(self):
        if not _db_exists():
            _create_test_db()
        rtc = RingtailCore(db_file="output.db")
        rtc.export_receptor_pdbqt()
        receptor_name = rtc.db_query("SELECT RecName FROM Receptors;")[0][0]
        receptor_file = receptor_name + ".pdbqt"

        assert os.path.exists(receptor_file)

        os.system("rm " + receptor_file)

    def test_generate_interactions_prepare_filters(self):
        test_filters = []
        rtc = RingtailCore()
        filters = Filters(
            {
                "hb_interactions": [("A:ARG:123:", True), ("A:VAL:124:", True)],
                "vdw_interactions": [("A:ARG:123:", True), ("A:VAL:124:", True)],
            }
        )
        interaction_combs = rtc._generate_interaction_combinations(filters.asdict(), 1)
        for ic in interaction_combs:
            nufilter = rtc._prepare_interaction_combo_filters(filters.asdict(), ic)
            test_filters.append(nufilter)

        assert (
            Filters(
                {
                    "hb_interactions": [("A:ARG:123:", True), ("A:VAL:124:", True)],
                    "vdw_interactions": [("A:ARG:123:", True)],
                }
            ).asdict()
            in test_filters
        )

        assert (
            Filters(
                {
                    "hb_interactions": [("A:ARG:123:", True), ("A:VAL:124:", True)],
                    "vdw_interactions": [("A:VAL:124:", True)],
                }
            ).asdict()
            in test_filters
        )

        assert (
            Filters(
                {
                    "hb_interactions": [("A:ARG:123:", True)],
                    "vdw_interactions": [("A:ARG:123:", True), ("A:VAL:124:", True)],
                }
            ).asdict()
            in test_filters
        )

        assert (
            Filters(
                {
                    "hb_interactions": [("A:VAL:124:", True)],
                    "vdw_interactions": [("A:ARG:123:", True), ("A:VAL:124:", True)],
                }
            ).asdict()
            in test_filters
        )

        assert (
            Filters(
                {
                    "hb_interactions": [("A:ARG:123:", True), ("A:VAL:124:", True)],
                    "vdw_interactions": [("A:ARG:123:", True), ("A:VAL:124:", True)],
                }
            ).asdict()
            in test_filters
        )

        assert len(test_filters) == 5

    def test_logfile_write(self):
        assert os.path.exists("different_log.txt")

        with open("different_log.txt") as f:
            target_line_no = None
            for line_no, line_content in enumerate(f):
                if "bookmark" in line_content:
                    target_line_no = line_no + 2
                if line_no == target_line_no:
                    break

        assert line_content == "11128, -7.25\n"

    def test_plot(self):
        if not _db_exists():
            _create_test_db()
        rtcore = RingtailCore(db_file="output.db")
        rtcore.filter(eworst=-7, bookmark_name="plot_data")
        rtcore.plot("plot_data")
        assert os.path.isfile("scatter.png") == True
        os.system("rm scatter.png")

    def test_export_bookmark_db(self):
        if not _db_exists():
            _create_test_db()
        rtc = RingtailCore(db_file="output.db")
        rtc.filter(eworst=-7, bookmark_name="export_db")
        bookmark_db_name = rtc.export_bookmark_db("export_db")

        assert os.path.exists(bookmark_db_name)
        rtc_bm = RingtailCore(db_file=bookmark_db_name)
        count = rtc_bm.table_length("Results")

        assert count == 8

        os.system("rm " + bookmark_db_name)

    def test_duplicate_handling(self):
        os.system("rm output.db*")

        rtc = RingtailCore(db_file="output.db")
        file = "test_data/adgpu/group1/1451.dlg.gz"
        rtc.add_results_from_files(file=file)
        # ensure three results rows were added
        result_count = rtc.table_length("Results")
        inter_count = rtc.table_length("Interactions")
        # add same file but replace the duplicate
        rtc.add_results_from_files(file=file, duplicate_handling="replace")
        result_count_replace = rtc.table_length("Results")
        inter_count_replace = rtc.table_length("Interactions")
        # add same file but ignore the duplicate
        rtc.add_results_from_files(file=file, duplicate_handling="ignore")
        result_count_ignore = rtc.table_length("Results")
        inter_count_ignore = rtc.table_length("Interactions")

        os.system("rm output.db*")
        # add same file but allow the duplicate
        rtc = RingtailCore(db_file="output.db")
        rtc.add_results_from_files(file=file)
        rtc.add_results_from_files(file=file)
        result_count_dupl = rtc.table_length("Results")
        inter_count_dupl = rtc.table_length("Interactions")

        assert (
            result_count
            == result_count_replace
            == result_count_ignore
            == result_count_dupl / 2
        )
        assert (
            inter_count
            == inter_count_replace
            == inter_count_ignore
            == inter_count_dupl / 2
        )

        os.system("rm output.db*")

    def test_db_num_poses_warning(self):
        from ringtail import LOGGER

        # make sure we make ringtail core object with log file
        rtc = RingtailCore(db_file="output.db", logging_level="DEBUG")

        # add results with max poses = 1
        rtc.add_results_from_files(
            file="test_data/adgpu/group1/1451.dlg.gz", max_poses=1
        )
        # add results with different max poses
        rtc.add_results_from_files(
            file="test_data/adgpu/group1/1620.dlg.gz", max_poses=4
        )
        warning_string = "The following database properties do not agree with the properties last used for this database: \nCurrent number of poses saved is 4 but database was previously set to 1."

        log_file = LOGGER._log_fp.baseFilename
        with open(log_file) as f:
            if warning_string in f.read():
                warning_worked = True
            else:
                warning_worked = False

        os.system("rm output.db*")

        assert warning_worked

    def test_reactive_filtering(self):
        rtc = RingtailCore(db_file="output.db")
        rtc.add_results_from_files(
            file_path="test_data/reactive/",
            store_all_poses=True,
            receptor_file="test_data/reactive/4j8m_m_rigid.pdbqt",
        )
        count_ligands_passing = rtc.filter(reactive_interactions=[("A:TYR:212:", True)])

        os.system("rm output.db*")

        assert count_ligands_passing == 10

    def test_polymer_receptor(self):
        rtc = RingtailCore(db_file="flexres.db")
        data_path = "test_data/flexres"
        rtc.add_results_from_files(
            file=data_path + "/ligand.pdbqt",
            docking_mode="vina",
            receptor_file=data_path + "/receptor.json",
            save_receptor=True,
        )
        receptor_items = rtc.get_receptor_object()

        os.system("rm flexres.db")

        assert receptor_items[0] == "receptor"
        assert not receptor_items[1]
        assert receptor_items[2] != None


class TestADNGHandling:
    def test_adng_stream(self):
        try:
            from rdkit import Chem
        except:
            return

        rtc = RingtailCore()
        rtc.save_receptor("test_data/adng/helix--scofu01.json")
        suppl = Chem.SDMolSupplier(
            "test_data/adng/docked_ligands.sdf",
            removeHs=False,
        )
        rtc.add_mol(suppl)

        results_count = rtc.table_length("Results")
        interaction_count = rtc.table_length("Interactions")
        os.system("rm output.db")

        assert results_count == 9
        assert interaction_count == 156

    def test_adng_file_add(self):
        adng_path = "test_data/adng"
        rtc = RingtailCore("output.db")
        rtc.add_results_from_files(
            file_path=adng_path,
            receptor_file=adng_path + "/helix--scofu01.json",
            save_receptor=True,
            docking_mode="adng",
        ),
        results_count = rtc.table_length("Results")
        interaction_count = rtc.table_length("Interactions")
        os.system("rm output.db")

        assert results_count == 9
        assert interaction_count == 156

    def test_adng_filtering(self):
        adng_path = "test_data/adng"
        rtc = RingtailCore("output.db")
        rtc.add_results_from_files(
            file_path=adng_path,
            receptor_file=adng_path + "/helix--scofu01.json",
            save_receptor=True,
            docking_mode="adng",
        ),
        count = rtc.filter(
            eworst=-13, ligand_substruct=["C=O"], vdw_interactions=[(":VAL::", True)]
        )
        os.system("rm output.db")
        assert count == 1


class TestVinaHandling:

    def test_vina_file_add(self):
        vina_path = "test_data/vina"
        rtc = RingtailCore("output.db")
        rtc.add_results_from_files(
            file_path=vina_path,
            receptor_file=vina_path + "/receptor.pdbqt",
            save_receptor=True,
            docking_mode="vina",
        ),
        count = rtc.table_length("Results")
        os.system("rm output.db*")

        assert count == 6

    def test_vina_string_add(self):
        vina_path = "test_data/vina"
        with open("test_data/vina/sample-result.pdbqt") as f:
            sample1 = f.read()
        with open("test_data/vina/sample-result-2.pdbqt") as f:
            sample2 = f.read()
        rtc = RingtailCore("output.db")
        rtc.save_receptor(
            vina_path + "/receptor.pdbqt",
        )
        rtc.add_results_from_vina_string(
            results={"sample1": sample1, "sample2": sample2},
        )
        count = rtc.table_length("Results")
        os.system("rm output.db*")

        assert count == 6

    def test_add_interactions(self):
        vina_path = "test_data/vina"
        rtc = RingtailCore("output.db", logging_level="DEBUG")
        rtc.add_results_from_files(
            file_path=vina_path,
            receptor_file=vina_path + "/receptor.pdbqt",
            save_receptor=True,
            docking_mode="vina",
        )
        unique_definition_count = rtc.table_length("Interaction_indices")
        interaction_count = rtc.table_length("Interactions")
        os.system("rm output.db*")

        assert unique_definition_count == 32
        assert interaction_count == 162

    def test_add_interactions_from_polymer(self):
        rtc = RingtailCore(db_file="flexres.db")
        data_path = "test_data/flexres"
        rtc.add_results_from_files(
            file=data_path + "/ligand.pdbqt",
            docking_mode="vina",
            receptor_file=data_path + "/receptor.json",
            save_receptor=True,
        )

        ligands_1 = rtc.table_length("Ligands")
        interactions_1 = rtc.table_length("Interactions")

        rtc = RingtailCore(db_file="flexres2.db")
        data_path = "test_data/flexres"
        rtc.add_results_from_files(
            file=data_path + "/ligand.pdbqt",
            docking_mode="vina",
            receptor_file=data_path + "/receptor.pdbqt",
            save_receptor=True,
        )
        ligands_2 = rtc.table_length("Ligands")
        interactions_2 = rtc.table_length("Interactions")

        os.system("rm flexres*.db")

        assert ligands_1 == ligands_2 == 1
        assert interactions_1 == interactions_2 == 85

    def test_db_dockingmode_warning(self):
        from ringtail import RingtailDefaults
        from ringtail import LOGGER

        rtc = RingtailCore(db_file="output.db", logging_level="DEBUG")
        rtc.add_results_from_files(file="test_data/adgpu/group1/1451.dlg.gz")
        rtc = RingtailCore(db_file="output.db", logging_level="DEBUG")
        rtc.add_results_from_files(
            file="test_data/vina/sample-result.pdbqt", docking_mode="vina"
        )

        warning_string = f"The following database properties do not agree with the properties last used for this database: \nCurrent docking mode is vina but last used docking mode of database is {RingtailDefaults.docking_mode}."
        log_file = LOGGER._log_fp.baseFilename
        with open(log_file, "r") as f:
            if warning_string in f.read():
                warning_worked = True
            else:
                warning_worked = False

        os.system("rm output.db*")

        assert warning_worked

    def test_various_filters_vina(self):
        vina_path = "test_data/vina"
        rtc = RingtailCore("output.db")
        rtc.add_results_from_files(
            file_path=vina_path,
            receptor_file=vina_path + "/receptor.pdbqt",
            save_receptor=True,
            docking_mode="vina",
        ),
        count = rtc.filter(eworst=-6, ligand_substruct=["[N]"])
        os.system("rm output.db")
        assert count == 1


class TestStorageMan:
    def test_fetch_summary_data(self):
        rtc = RingtailCore("output.db")
        rtc.add_results_from_files(
            file_list="test_data/filelist1.txt",
            receptor_file="test_data/adgpu/4j8m.pdbqt",
            save_receptor=True,
        )
        with rtc.storageman:
            summ_dict = rtc.storageman.fetch_summary_data()
        rounded_dict = {
            k: round(v, 2) if isinstance(v, float) else v for k, v in summ_dict.items()
        }
        assert rounded_dict == {
            "num_ligands": 3,
            "num_poses": 7,
            "num_unique_interactions": 57,
            "num_interacting_residues": 30,
            "min_docking_score": -6.66,
            "max_docking_score": -4.98,
            "1%_docking_score": -6.66,
            "10%_docking_score": -6.66,
            "min_leff": -0.44,
            "max_leff": -0.35,
            "1%_leff": -0.44,
            "10%_leff": -0.44,
        }

    def test_bookmark_info(self):
        rtc = RingtailCore("output.db")
        rtc.add_results_from_files(
            file_path="test_data/adgpu/group2",
        )
        rtc.filter(
            eworst=-3,
            hb_interactions=[("A:VAL:279:", True), ("A:LYS:162:", True)],
            vdw_interactions=[("A:VAL:279:", True)],
            bookmark_name="bookmark_info",
        )
        qb = QueryBuilder()
        query_string = (
            qb.SELECT("filters")
            .FROM("Filters")
            .WHERE("name='bookmark_info'")
            .build()[0]
        )
        bookmark_filters_db_str = rtc.db_query(query_string)[0][0]

        assert (
            json.loads(bookmark_filters_db_str)
            == Filters(
                {
                    "eworst": -3.0,
                    "vdw_interactions": [["A:VAL:279:", True]],
                    "hb_interactions": [["A:VAL:279:", True], ["A:LYS:162:", True]],
                }
            ).asdict()
        )

    def test_version_info(self):

        from importlib.metadata import version

        rtc = RingtailCore("output.db")
        with rtc.storageman:
            versionmatch, db_version = rtc.storageman.check_ringtaildb_version()
        os.system("rm output.db*")
        assert versionmatch
        assert db_version == version("ringtail")


class TestMergeDB:
    def test_db_write(self):
        rtc1 = RingtailCore("primary.db")
        rtc1.add_results_from_files("test_data/adgpu/group1/1451.dlg.gz")

        rtc2 = RingtailCore("secondary.db")
        rtc2.add_results_from_files("test_data/adgpu/group1/1620.dlg.gz")

        rtc3 = RingtailCore("tertiary.db")
        rtc3.add_results_from_files("test_data/adgpu/group1/1751.dlg.gz")

        # they should all have one ligand each
        assert (
            rtc1.table_length("Ligands")
            == rtc2.table_length("Ligands")
            == rtc3.table_length("Ligands")
            == 1
        )

    def test_before_merge(self):
        rtc1 = RingtailCore("primary.db")
        # should not be any poses in this interval
        assert rtc1.filter(eworst=-2, ebest=-5) == 0
        assert rtc1.filter(eworst=-5) == 1
        assert rtc1.table_length("filtered_poses") == 3

    def test_after_merge(self):
        rtc1 = RingtailCore("primary.db")
        rtc1.merge_databases("secondary.db", False)
        rtc1.merge_databases("tertiary.db", False)
        # this should add two more ligands
        assert rtc1.table_length("Ligands") == 3
        # should now be data in this interval
        assert rtc1.filter(eworst=-2, ebest=-5) == 2

    def test_check_PKs(self):
        rtc2 = RingtailCore("secondary.db")

        # get best ranked pose id in secondary database for ligand 1620
        secondary_db_pose_as_main = rtc2.db_query(
            "SELECT Pose_ID FROM Results WHERE pose_rank = 1 AND ligand_id = (SELECT ligand_id FROM Ligands WHERE LigName = '1620')"
        )[0][0]
        assert secondary_db_pose_as_main == 1

        # compare to pose id for best ranked pose for same ligand in the merged database
        rtc1 = RingtailCore("primary.db")
        secondary_db_pose_as_merged = rtc1.db_query(
            "SELECT Pose_ID FROM Results WHERE pose_rank = 1 AND ligand_id = (SELECT ligand_id FROM Ligands WHERE LigName = '1620')"
        )[0][0]
        assert secondary_db_pose_as_merged != secondary_db_pose_as_main

        os.system("rm primary.db secondary.db tertiary.db")


class TestLogger:

    def test_set_log_level(self):
        from ringtail.logutils import RaccoonLogger

        logger = RaccoonLogger()
        logger.set_level("info")
        log_level = logger.level()
        assert log_level == "INFO"


class TestOptions:
    def test_object_checks(self):
        # checking that incompatible options are handled
        from ringtail.ringtailoptions import Filters

        rtc = RingtailCore()
        rtc.add_results_from_files(file_list="test_data/filelist1.txt")
        rtc.filters = Filters({"score_percentile": 20})
        assert rtc.filters.eworst == None
        assert rtc.filters.score_percentile == 20

        # conflicting options, score percentile should be set to none
        rtc.filters = Filters({"score_percentile": 20, "eworst": -6})
        rtc.filters.checks()

        assert rtc.filters.eworst == -6
        assert rtc.filters.score_percentile == None

    def test_overwrite_db(self):
        rtc = RingtailCore()
        rtc.add_results_from_files(file_list="test_data/filelist1.txt")
        count_old_db = rtc.table_length("Ligands")

        rtc.add_results_from_files(file_list="test_data/filelist2.txt", overwrite=True)
        count_new_db = rtc.table_length("Ligands")

        assert count_old_db == 3
        assert count_new_db == 2

    def test_cleanup(self):
        # Alter this method if you wish to not delete all log files after testing automatically
        os.system("rm *_ringtail.log")
        os.system("rm output.db* output2.db")
        os.system("rm different_log.txt")
        os.system("rm cluster_log.txt output_log.txt")
