#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail unit testing
#
from ringtail import StorageManagerSQLite, RingtailCore, Filters
import sqlite3
import os
import json
import pytest

class Test_StorageManSQLite:
    #TODO rewrite tests to go better with API update? 
    def test_fetch_summary_data(self):
        
        status1 = os.system(
            "python ../scripts/rt_process_vs.py write -d --file_list filelist1.txt"
        )
        rtcore = RingtailCore("output.db")
        with rtcore.storageman: summ_dict = rtcore.storageman.fetch_summary_data()
        assert summ_dict == {'num_ligands': 3, 'num_poses': 7, 'num_unique_interactions': 57, 'num_interacting_residues': 30, 'min_docking_score': -6.66, 'max_docking_score': -4.98, '1%_docking_score': -6.66, '10%_docking_score': -6.66, 'min_leff': -0.444, 'max_leff': -0.35000000000000003, '1%_leff': -0.444, '10%_leff': -0.444}
        
        os.system("rm output.db")

    def test_bookmark_info(self):
        rtcore = RingtailCore("output.db")
        rtcore.add_results_from_files(file_path=[["test_data/"]], recursive=True)

        rtcore.filter(eworst = -3, 
                    vdw_interactions=[('A:ARG:123:', True), ('A:VAL:124:', True)],
                    hb_interactions=[('A:ARG:123:', True)],
                    ligand_operator='OR')

        conn = sqlite3.connect("output.db")
        cur = conn.cursor()
        bookmark = cur.execute("SELECT filters FROM Bookmarks WHERE Bookmark_name LIKE 'passing_results'")
        bookmark_filters_db_str = bookmark.fetchone()[0]
        filters = {'eworst': -3.0, 'ebest': None, 'leworst': None, 'lebest': None, 'score_percentile': None, 'le_percentile': None, 'vdw_interactions': [('A:ARG:123:', True), ('A:VAL:124:', True)], 'hb_interactions': [('A:ARG:123:', True)], 'reactive_interactions': [], 'interactions_count': [], 'react_any': None, 'max_miss': 0, 'ligand_name': [], 'ligand_substruct': [], 'ligand_substruct_pos': [], 'ligand_max_atoms': None, 'ligand_operator': 'OR'}
        assert bookmark_filters_db_str == json.dumps(filters)
        cur.close()
        conn.close()

        os.system("rm output.db")   

    def test_version_info(self):
        rtcore = RingtailCore("output.db")
        rtcore.add_results_from_files(file_path=[["test_data/"]], recursive=True)
        with rtcore.storageman: versionmatch, version = rtcore.storageman.check_ringtaildb_version()
        
        assert versionmatch
        assert int(version) == 110  # NOTE: update for new versions
        os.system("rm output.db")

    #TODO
    def test_context_manager(self):
        # ensure database open correctly
        # ensure database closes correctly
        assert True

class Test_RingtailCore:

    def test_prepare_filters_for_storageman(self):
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



        assert len(test_filters) == 5
    
        os.system("rm output_log.txt output.db")

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