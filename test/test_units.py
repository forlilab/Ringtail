#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail unit testing
#
from ringtail import RingtailCore
import sqlite3
import os
import json
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
        return count
    return __dbconnect

@pytest.fixture(scope='class')
def dbquery():
    conn = sqlite3.connect("output.db")
    curs = conn.cursor()
    def __dbconnect(query):
        curs.execute(query)
        return curs
    yield __dbconnect
    curs.close()
    conn.close()

class Test_StorageManSQLite:
    # Common setup of ringtail core
    rtcore = RingtailCore("output.db")
    rtcore.add_results_from_files(file_list = [['filelist1.txt']],
                                     recursive = True,
                                     receptor_file="test_data/4j8m.pdbqt",
                                     save_receptor=True)

    def test_fetch_summary_data(self):
        with self.rtcore.storageman: summ_dict = self.rtcore.storageman.fetch_summary_data()
        assert summ_dict == {'num_ligands': 3, 'num_poses': 7, 'num_unique_interactions': 57, 'num_interacting_residues': 30, 'min_docking_score': -6.66, 'max_docking_score': -4.98, '1%_docking_score': -6.66, '10%_docking_score': -6.66, 'min_leff': -0.444, 'max_leff': -0.35000000000000003, '1%_leff': -0.444, '10%_leff': -0.444}

    def test_bookmark_info(self):
        self.rtcore.filter(eworst = -3, 
                    vdw_interactions=[('A:ARG:123:', True), ('A:VAL:124:', True)],
                    hb_interactions=[('A:ARG:123:', True)],
                    ligand_operator='OR',)
        conn = sqlite3.connect("output.db")
        cur = conn.cursor()
        bookmark = cur.execute("SELECT filters FROM Bookmarks WHERE Bookmark_name LIKE 'passing_results'")
        bookmark_filters_db_str = bookmark.fetchone()[0]
        filters = {"eworst": -3.0, "ebest": None, "leworst": None, "lebest": None, "score_percentile": None, "le_percentile": None, "vdw_interactions": [["A:ARG:123:", True], ["A:VAL:124:", True]], "hb_interactions": [["A:ARG:123:", True]], "reactive_interactions": [], "interactions_count": [], "react_any": None, "max_miss": 0, "ligand_name": [], "ligand_substruct": [], "ligand_substruct_pos": [], "ligand_max_atoms": None, "ligand_operator": "OR"}
        assert bookmark_filters_db_str == json.dumps(filters)
        cur.close()
        conn.close()

    def test_version_info(self):
        with self.rtcore.storageman: versionmatch, version = self.rtcore.storageman.check_ringtaildb_version()
        
        assert versionmatch
        assert int(version) == 110  # NOTE: update for new versions

class Test_RingtailCore:

    def test_get_defaults(self):
        from ringtail import ringtailoptions
        defaults = RingtailCore.get_defaults("resultsmanopts")
        object_dict = ringtailoptions.ResultsProcessingOptions().todict()
        assert defaults == object_dict
    
    def test_db_setup(self):
        RingtailCore(db_file="output.db", logging_level="debug")
        assert os.path.isfile("output.db") == True
    
    def test_save_receptor(self, countrows):
        rtc = RingtailCore(db_file="output.db")
        rtc.add_results_from_files(receptor_file="test_data/4j8m.pdbqt.gz",
                                        save_receptor=True)
        count = countrows("SELECT COUNT(*) FROM Receptors WHERE receptor_object NOT NULL")
        assert count == 1

    def test_folder1(self, countrows):
        rtc = RingtailCore(db_file="output.db")
        rtc.add_results_from_files(file_path = [['test_data/group1']],
                                            append_results=True)
        count = countrows("SELECT COUNT(*) FROM Ligands")
        assert count == 138

    def test_produce_summary(self):        
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

        assert len(summary_items.data) == 40

    def test_create_rdkitmol(self):
        pass

    def test_append_to_database(self, countrows):
        rtc = RingtailCore(db_file="output.db")
        rtc.add_results_from_files(file_path=[['test_data/group2/']], append_results=True)
        count = countrows("SELECT COUNT(*) FROM Ligands")

        assert count == 217

    def test_filter(self, countrows):
        rtc = RingtailCore(db_file="output.db")
        count_ligands_passing = rtc.filter(eworst = -6, hb_interactions=[('A:VAL:279:', True), ('A:LYS:162:', True)], vdw_interactions=[('A:VAL:279:', True), ('A:LYS:162:', True)], max_miss = 1)

        assert count_ligands_passing == 36

    def test_similar_ligands(self, monkeypatch):
        rtc = RingtailCore(db_file="output.db")
        ligand_name = "287065"
        rtc.filter(ebest = -6, mfpt_cluster=0.5)
        monkeypatch.setattr('builtins.input', lambda _: 0) # provide terminal input
        number_similar = rtc.find_similar_ligands(ligand_name)

        assert number_similar == 8

    #TODO not working
    def test_plot(self):
        return
        
        rtc = RingtailCore(db_file="output.db")
        rtc.add_results_from_files(file_path = [['test_data/group3']],
                                            append_results=True)
        rtc.filter(ebest = -6)
        rtc.plot()

        assert os.path.isfile("scatter.png") == True

    def test_write_sdfs(self):
        rtc = RingtailCore(db_file="output.db")
        rtc.filter(eworst = -7)
        rtc.write_molecule_sdfs(".")

        import glob
        sdf_files = glob.glob("*.sdf")
        expected = ['3961.sdf', '5995.sdf', '11128.sdf', '11991.sdf', '13974.sdf', '15776.sdf', '136065.sdf']
        assert len(sdf_files) == len(expected)
        for f in sdf_files:
            assert f in expected
            os.remove(f)

    def test_pymol(self):
        pass

    def test_export_csv(self):
        rtc = RingtailCore(db_file="output.db")
        rtc.filter(eworst = -7)
        rtc.export_csv("Ligands", "Ligands.csv", True)
        
        assert os.path.exists("Ligands.csv")
        os.system("rm Ligands.csv")

    def test_export_bookmark_db(self):
        pass

    def export_receptor(self):
        pass

    def test_get_filterdata(self):
        pass

    def test_generate_interactions_prepare_filters(self):
        test_filters = []
        rtc = RingtailCore()
        rtc.process_mode = "read"
        rtc._set_general_options(docking_mode="dlg", logging_level="DEBUG")
        rtc.set_filters(hb_interactions=[("A:ARG:123:", True), ("A:VAL:124:", True)], vdw_interactions=[("A:ARG:123:", True), ("A:VAL:124:", True)])
        interaction_combs = rtc._generate_interaction_combinations(1)
        for ic in interaction_combs:
            nufilter=rtc._prepare_filters_for_storageman(ic)
            test_filters.append(nufilter)

        assert {'eworst': None, 'ebest': None, 'leworst': None, 'lebest': None, 'score_percentile': None, 'le_percentile': None, 'vdw_interactions': [('A:ARG:123:', True), ('A:VAL:124:', True)], 'hb_interactions': [('A:ARG:123:', True)], 'reactive_interactions': [], 'interactions_count': [], 'react_any': None, 'max_miss': 0, 'ligand_name': [], 'ligand_substruct': [], 'ligand_substruct_pos': [], 'ligand_max_atoms': None, 'ligand_operator': 'OR'} in test_filters

        assert {'eworst': None, 'ebest': None, 'leworst': None, 'lebest': None, 'score_percentile': None, 'le_percentile': None, 'vdw_interactions': [('A:ARG:123:', True), ('A:VAL:124:', True)], 'hb_interactions': [('A:VAL:124:', True)], 'reactive_interactions': [], 'interactions_count': [], 'react_any': None, 'max_miss': 0, 'ligand_name': [], 'ligand_substruct': [], 'ligand_substruct_pos': [], 'ligand_max_atoms': None, 'ligand_operator': 'OR'} in test_filters

        assert {'eworst': None, 'ebest': None, 'leworst': None, 'lebest': None, 'score_percentile': None, 'le_percentile': None, 'vdw_interactions': [('A:ARG:123:', True)], 'hb_interactions': [('A:ARG:123:', True), ('A:VAL:124:', True)], 'reactive_interactions': [], 'interactions_count': [], 'react_any': None, 'max_miss': 0, 'ligand_name': [], 'ligand_substruct': [], 'ligand_substruct_pos': [], 'ligand_max_atoms': None, 'ligand_operator': 'OR'} in test_filters

        assert {'eworst': None, 'ebest': None, 'leworst': None, 'lebest': None, 'score_percentile': None, 'le_percentile': None, 'vdw_interactions': [('A:VAL:124:', True)], 'hb_interactions': [('A:ARG:123:', True), ('A:VAL:124:', True)], 'reactive_interactions': [], 'interactions_count': [], 'react_any': None, 'max_miss': 0, 'ligand_name': [], 'ligand_substruct': [], 'ligand_substruct_pos': [], 'ligand_max_atoms': None, 'ligand_operator': 'OR'} in test_filters

        assert {'eworst': None, 'ebest': None, 'leworst': None, 'lebest': None, 'score_percentile': None, 'le_percentile': None, 'vdw_interactions': [('A:ARG:123:', True), ('A:VAL:124:', True)], 'hb_interactions': [('A:ARG:123:', True), ('A:VAL:124:', True)], 'reactive_interactions': [], 'interactions_count': [], 'react_any': None, 'max_miss': 0, 'ligand_name': [], 'ligand_substruct': [], 'ligand_substruct_pos': [], 'ligand_max_atoms': None, 'ligand_operator': 'OR'} in test_filters

        os.system("rm output_log.txt output.db")
        assert len(test_filters) == 5

class Test_logger:
    def test_set_log_level(self):
        from ringtail import logger
        logger.setLevel("info")
        log1 = logger.level()
        logger.setLevel(50)
        log2 = logger.level()
        assert log1 == log2 == "INFO" 
    
    def test_loggerfile_format(self):
        from ringtail import logger, exceptions
        try:
            raise(exceptions.OptionError("This is a test error."))
        except Exception as e:
            pass

        with open(logger.filename()) as f:
            for line in f:
                pass
        last_line = line
        keywords = ["ERROR", "test_units.py[", "ringtail.exceptions:This is a test error."]
        assert all(x in last_line for x in keywords)
        

#TODO
class Test_exceptions:
    def test_option_error(self):
        # make sure errors are caught and surfaced properly
            # in command line
            # in API
            # write to stout when necessary
        assert True

#TODO 
class Test_outputmanager:
    def test_logfile_write(self):
        # if making new
        # if overwriting
        # is there an option to add?
        assert True

#TODO
class Test_options:
    def test_options_type_check(self):
        # if corrrect type
        # if wrong type
        # if re-assigning
        # does it initialize right? 
        # does it re-initialize and default right?
        assert True
    
    def test_file_format(self):
        # that files are double lists ready for results manager
        # file extension is correct
        # folders are proper directories
        assert True
    
    def test_object_checks(self):
        # take one object or two and make sure checks are performed when appropriate
        assert True