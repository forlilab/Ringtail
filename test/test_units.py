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

class TestRingtailCore:

    def test_get_defaults(self):
        os.system("rm output.db output_log.txt")
        from ringtail import ringtailoptions
        defaults = RingtailCore.get_defaults("resultsmanopts")
        object_dict = ringtailoptions.ResultsProcessingOptions().todict()
        assert defaults == object_dict

    def test_add_folder(self, countrows):
        rtc = RingtailCore(db_file="output.db")
        rtc.add_results_from_files(file_path = [['test_data/group1']])
        count = countrows("SELECT COUNT(*) FROM Ligands")
        assert count == 138
    
    def test_save_receptor(self, countrows):
        rtc = RingtailCore(db_file="output.db", logging_level="debug")
        count0 = countrows("SELECT COUNT(*) FROM Receptors WHERE receptor_object NOT NULL")
        
        assert count0 == 0
        
        rtc.save_receptor(receptor_file="test_data/4j8m.pdbqt")
        count = countrows("SELECT COUNT(*) FROM Receptors WHERE receptor_object NOT NULL")
        
        assert count == 1
    
    def test_produce_summary(self):        
        
        # Ensure storage error thrown if no data in database
        from ringtail import exceptions as e
        with pytest.raises(e.RTCoreError):
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

    def test_append_to_database(self, countrows):
        rtc = RingtailCore(db_file="output.db")
        rtc.add_results_from_files(file_path=[['test_data/group2/']], append_results=True)
        count = countrows("SELECT COUNT(*) FROM Ligands")

        assert count == 217

    def test_filter(self):
        rtc = RingtailCore(db_file="output.db")
        count_ligands_passing = rtc.filter(eworst = -6, hb_interactions=[('A:VAL:279:', True), ('A:LYS:162:', True)], vdw_interactions=[('A:VAL:279:', True), ('A:LYS:162:', True)], max_miss = 1)

        assert count_ligands_passing == 36

    def test_get_filterdata(self):
        rtc = RingtailCore(db_file="output.db")
        rtc.filter(eworst = -7) 
        log_file_name = "output_log_test.txt"
        rtc.set_output_options(log_file=log_file_name)
        rtc.get_previous_filter_data("delta, ref_rmsd", bookmark_name="passing_results")


        with open(log_file_name) as f:
            file_contents = f.read()
        import linecache
        final_line = linecache.getline(log_file_name, 10)

        assert "'11991', '11991', 0.0, 226.06" in file_contents
        assert "'3961', '3961', 0.0, 215.96" in file_contents
        assert final_line == "***************\n"

        os.system(("rm " + log_file_name))

    def test_similar_ligands(self, monkeypatch):
        rtc = RingtailCore(db_file="output.db")
        ligand_name = "287065"
        rtc.filter(ebest = -6, mfpt_cluster=0.5)
        monkeypatch.setattr('builtins.input', lambda _: 0) # provides terminal input
        number_similar = rtc.find_similar_ligands(ligand_name)

        assert number_similar == 8

    def test_create_rdkitmol(self):
        bookmark_name = "rdkit_test"
        rtc = RingtailCore(db_file="output.db")
        rtc.filter(ebest = -3, bookmark_name=bookmark_name)
        rdkit_dict = rtc.ligands_rdkit_mol(bookmark_name=bookmark_name)
        assert len(rdkit_dict) == 8
        # grab one molecule from bookmark and check number of atoms
        num_of_atoms = rdkit_dict["14303"]["ligand"].GetNumAtoms()
        assert num_of_atoms == 10

    def test_write_sdfs(self):
        sdf_path = "sdf_files"
        rtc = RingtailCore(db_file="output.db")
        rtc.filter(eworst = -7)
        rtc.write_molecule_sdfs(sdf_path)

        # ensure correct number of files written
        sdf_files = os.listdir(sdf_path)
        expected = ['3961.sdf', '5995.sdf', '11128.sdf', '11991.sdf', '13974.sdf', '15776.sdf', '136065.sdf']
        assert len(sdf_files) == len(expected)

        # ensure contents is correct
        with open('sdf_files/136065.sdf') as sdf:
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
        rtc = RingtailCore(db_file="output.db")
        rtc.filter(eworst = -7)
        rtc.export_csv("Ligands", "Ligands.csv", True)
        
        assert os.path.exists("Ligands.csv")
        os.system("rm Ligands.csv")

    def test_export_receptor(self, dbquery):
        rtc = RingtailCore(db_file="output.db")
        rtc.export_receptors()
        curs = dbquery("SELECT RecName FROM Receptors;")
        receptor_name = curs.fetchone()[0]
        receptor_file = receptor_name + ".pdbqt"
        
        assert os.path.exists(receptor_file)

        os.system("rm " + receptor_file)
    
    def test_generate_interactions_prepare_filters(self):
        test_filters = []
        rtc = RingtailCore()
        rtc.docking_mode="dlg"
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

        assert len(test_filters) == 5

    def test_logfile_write(self):
        rtc = RingtailCore("output.db")
        assert os.path.exists("output_log.txt")

        with open("output_log.txt") as f:
            for line_no, line_content in enumerate(f):
                if line_no == 28:
                    break

        assert line_content == "'11128', -7.25\n"

    def test_plot(self):     
        rtcore = RingtailCore(db_file="output.db")
        rtcore.filter(eworst = -7)
        rtcore.plot()
        assert os.path.isfile("scatter.png") == True
        os.system("rm scatter.png")

    def test_export_bookmark_db(self):
        rtc = RingtailCore(db_file="output.db")
        rtc.filter(eworst = -7)
        rtc.export_bookmark_db()
        
        assert os.path.exists(rtc.db_file)

        conn = sqlite3.connect(rtc.db_file)
        curs = conn.cursor()
        curs.execute("SELECT COUNT(*) FROM Results")
        count = curs.fetchone()[0]
        curs.close()
        conn.close()

        assert count == 7

        os.system("rm " + rtc.db_file)

    def test_duplicate_handling(self, countrows):
        os.system("rm output.db output_log.txt")

        rtc = RingtailCore(db_file="output.db")
        file = 'test_data/group1/1451.dlg.gz'
        rtc.add_results_from_files(file = file, duplicate_handling="replace")
        # ensure three results rows were added
        count = countrows("SELECT COUNT(*) FROM Results")
        # add same file but replace the duplicate
        rtc.add_results_from_files(file = file, append_results=True)
        count_replace = countrows("SELECT COUNT(*) FROM Results")

        os.system("rm output.db")
        rtc = RingtailCore(db_file="output.db")
        rtc.add_results_from_files(file = file, duplicate_handling="ignore")
        # add same file but ignore the duplicate
        rtc.add_results_from_files(file = file, append_results=True)
        count_ignore = countrows("SELECT COUNT(*) FROM Results")

        os.system("rm output.db")
        rtc = RingtailCore(db_file="output.db")
        rtc.add_results_from_files(file = file)
        # add same file but allow the duplicate
        rtc.add_results_from_files(file = file, append_results=True)
        count_dupl = countrows("SELECT COUNT(*) FROM Results")
        assert count == count_replace == count_ignore == count_dupl/2

        os.system("rm output.db")
        
    def test_db_num_poses_warning(self):
        rtc = RingtailCore(db_file="output.db")
        rtc.add_results_from_files(file = [['test_data/group1/1451.dlg.gz']], max_poses=1)
        rtc.add_results_from_files(file = [['test_data/group1/1620.dlg.gz']], max_poses=4, append_results=True)
        from ringtail import logger
        warning_string = 'The following database properties do not agree with the properties last used for this database: \nCurrent number of poses saved is 4 but database was previously set to 1.'
        with open(logger.filename(), "r") as f:
            if warning_string in f.read():
                warning_worked = True
            else:
                warning_worked = False
                
        os.system("rm output.db")

        assert warning_worked

class TestVinaHandling:

    def test_vina_file_add(self, countrows):
        vina_path = 'test_data/vina'
        rtc = RingtailCore("output.db")
        rtc.docking_mode="vina"
        rtc.add_results_from_files(file_path=[[vina_path]], file_pattern = "*.pdbqt*", receptor_file=vina_path+"/receptor.pdbqt", save_receptor=True),
        count = countrows("SELECT COUNT(*) FROM Results")
        os.system("rm output.db")

        assert count == 6

    def test_vina_string_add(self, countrows):
        vina_path = 'test_data/vina'
        with open('test_data/vina/sample-result.pdbqt') as f:
            sample1 = f.read()
        with open('test_data/vina/sample-result-2.pdbqt') as f:
            sample2 = f.read()
        rtc = RingtailCore("output.db")
        rtc.add_results_from_vina_string(results_strings={"sample1":sample1, "sample2": sample2}, receptor_file=vina_path+"/receptor.pdbqt", save_receptor=True)
        count = countrows("SELECT COUNT(*) FROM Results")
        os.system("rm output.db")

        assert count == 6

    def test_add_interactions(self, countrows):
        vina_path = 'test_data/vina'
        rtc = RingtailCore("output.db")
        rtc.docking_mode="vina"
        rtc.add_results_from_files(file_path=vina_path, file_pattern = "*.pdbqt*", receptor_file=vina_path+"/receptor.pdbqt", save_receptor=True, add_interactions=True)
        count = countrows("SELECT COUNT(*) FROM Interaction_indices")
        os.system("rm output.db")

        assert count == 45

    def test_db_dockingmode_warning(self):
        rtc = RingtailCore(db_file="output.db")
        rtc.add_results_from_files(file = 'test_data/group1/1451.dlg.gz')
        rtc.add_results_from_files(file = 'test_data/vina/sample-result.pdbqt', append_results=True)
        
        from ringtail import logger
        warning_string = 'The following database properties do not agree with the properties last used for this database: \nCurrent docking mode is vina but last used docking mode of database is dlg.'
        with open(logger.filename(), "r") as f:
            if warning_string in f.read():
                warning_worked = True
            else:
                warning_worked = False
                
        os.system("rm output.db")

        assert warning_worked

class TestStorageMan:
   
    def test_storageman_setup(self):
        rtc = RingtailCore("output.db")
        rtc.add_results_from_files(file_list = 'filelist1.txt',
                                    recursive = True,
                                    receptor_file="test_data/4j8m.pdbqt",
                                    save_receptor=True)
        
        storageman_attributes = {'append_results': rtc.storageman.append_results, 
                                 'duplicate_handling': rtc.storageman.duplicate_handling, 
                                 'filter_bookmark': rtc.storageman.filter_bookmark, 
                                 'overwrite': rtc.storageman.overwrite, 
                                 'order_results': rtc.storageman.order_results, 
                                 'outfields': rtc.storageman.outfields, 
                                 'output_all_poses': rtc.storageman.output_all_poses, 
                                 'mfpt_cluster': rtc.storageman.mfpt_cluster, 
                                 'interaction_cluster': rtc.storageman.interaction_cluster, 
                                 'bookmark_name': rtc.storageman.bookmark_name}
        storageman_attributes_defaults = RingtailCore.get_defaults("storageopts")
        # ensure defaults values are set correctly and do not change during processing
        assert storageman_attributes == storageman_attributes_defaults

    def test_fetch_summary_data(self):
        rtc = RingtailCore("output.db")
        with rtc.storageman: summ_dict = rtc.storageman.fetch_summary_data()
        assert summ_dict == {'num_ligands': 3, 'num_poses': 7, 'num_unique_interactions': 57, 'num_interacting_residues': 30, 'min_docking_score': -6.66, 'max_docking_score': -4.98, '1%_docking_score': -6.66, '10%_docking_score': -6.66, 'min_leff': -0.444, 'max_leff': -0.35000000000000003, '1%_leff': -0.444, '10%_leff': -0.444}

    def test_bookmark_info(self, dbquery):
        rtc = RingtailCore("output.db")
        rtc.filter(eworst = -3, 
                    vdw_interactions=[('A:ARG:123:', True), ('A:VAL:124:', True)],
                    hb_interactions=[('A:ARG:123:', True)],
                    ligand_operator='OR',)
        curs = dbquery("SELECT filters FROM Bookmarks WHERE Bookmark_name LIKE 'passing_results'")
        bookmark_filters_db_str = curs.fetchone()[0]
        filters = {"eworst": -3.0, "ebest": None, "leworst": None, "lebest": None, "score_percentile": None, "le_percentile": None, "vdw_interactions": [["A:ARG:123:", True], ["A:VAL:124:", True]], "hb_interactions": [["A:ARG:123:", True]], "reactive_interactions": [], "interactions_count": [], "react_any": None, "max_miss": 0, "ligand_name": [], "ligand_substruct": [], "ligand_substruct_pos": [], "ligand_max_atoms": None, "ligand_operator": "OR"}
        assert bookmark_filters_db_str == json.dumps(filters)

    def test_version_info(self):
        rtc = RingtailCore("output.db")
        with rtc.storageman: versionmatch, version = rtc.storageman.check_ringtaildb_version()
        os.system("rm output.db output_log.txt")
        assert versionmatch
        assert int(version) == 110  # NOTE: update for new database schema versions

class TestConfigFile:

    def test_generate_config_file(self):
        RingtailCore.generate_config_json_template()
        filepath = "config.json"

        assert os.path.exists(filepath)

        # Assure file is created, and populate it with non-default values
        with open(filepath, "r") as f:
            data = json.load(f)
        # all fields I want to change
        data["fileobj"]["file_path"] = [['test_data/']]
        data["fileobj"]["file_path"]
        data["fileobj"]["recursive"] = True
        data["fileobj"]["receptor_file"] = "test_data/4j8m.pdbqt"
        data["fileobj"]["save_receptor"] = True
        data["filters"]["eworst"] = -6
        data["filters"]["vdw_interactions"] = [('A:VAL:279:', True), ('A:LYS:162:', True)]
        data["filters"]["hb_interactions"] = [('A:VAL:279:', True), ('A:LYS:162:', True)]
        data["filters"]["max_miss"] = 1
        with open(filepath, "w") as f:
            f.write(json.dumps(data, indent=4))

    def test_adding_results(self, dbquery):
        rtcore = RingtailCore(db_file="output.db", logging_level = "DEBUG")
        (file_dict, write_dict, _, _) = rtcore.add_config_from_file(config_file="config.json")
        rtcore.add_results_from_files(filesources_dict= file_dict, options_dict=write_dict)
        curs = dbquery("""SELECT COUNT(*) FROM Results;""")
        count = curs.fetchone()[0]

        assert count == 645
        
    def test_filter(self):
        rtcore = RingtailCore(db_file="output.db", logging_level = "DEBUG")
        (_, _, _, filters_dict) = rtcore.add_config_from_file(config_file="config.json")
        count_ligands_passing = rtcore.filter(filters_dict=filters_dict)
        os.system("rm output_log.txt output.db config.json")
        assert count_ligands_passing == 51 

class TestLogger:

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
        
class TestOptions:
    def test_option_error(self):
        from ringtail import exceptions as e
        with pytest.raises(e.OptionError):
            RingtailCore.get_defaults("not_a_ringtail_object")

    def test_options_type_check(self):
        rtc = RingtailCore()
        # set to wrong type, ensure reset to default value
        rtc.set_storageman_attributes(outfields=5)
        assert rtc.storageman.outfields == "Ligand_name,e"
    
    def test_object_checks(self):
        # checking that incompatible options are handled
        rtc = RingtailCore()
        rtc.set_filters(score_percentile = 20)
        assert rtc.filters.eworst == None
        assert rtc.filters.score_percentile == 20

        # conflicting options, score percentile should be set to none
        rtc.set_filters(eworst =-6)
        assert rtc.filters.eworst == -6
        assert rtc.filters.score_percentile == None

    def test_set_order(self):
        rtc = RingtailCore()
        rtc.set_filters(dict = {"eworst": -5})
        assert rtc.filters.eworst == -5
        # ensure single options overwrite dict options
        rtc.set_filters(eworst =-6, dict = {"eworst": -5})
        assert rtc.filters.eworst == -6
