#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail API testing
#

import sqlite3
import os
import pytest
from ringtail import RingtailCore 
import json

@pytest.fixture
def cur():
    conn = sqlite3.connect("outputapi.db")
    curs = conn.cursor()
    yield curs
    curs.close()
    conn.close()

@pytest.fixture
def countrows():
    def __dbconnect(query):
        conn = sqlite3.connect("outputapi.db")
        curs = conn.cursor()
        curs.execute(query)
        count = curs.fetchone()[0]
        curs.close()
        conn.close()
        return count
    return __dbconnect

@pytest.fixture(scope='class')
def dbquery():
    conn = sqlite3.connect("outputapi.db")
    curs = conn.cursor()
    def __dbconnect(query):
        curs.execute(query)
        return curs
    yield __dbconnect
    curs.close()
    conn.close()

class TestAPI:
    """
    Notes:
    Tests to ensure the API operates properly, although this is really more the RingtailCore 
    tests I am writing in unit tests, since the API is used by everything the command line tool.
    Maybe depreceate this test file after I complete the unit tests
    """
    # Setup
    rtstorage = RingtailCore(db_file="outputapi.db", logging_level="debug")

    rtstorage.add_results_from_files(file_path = [['test_data/']],
                                     recursive = True,
                                     receptor_file="test_data/4j8m.pdbqt",
                                     save_receptor=True)

    def test_receptor_save(self, dbquery):
        curs = dbquery("""SELECT COUNT(*) FROM Receptors WHERE receptor_object NOT NULL;""" )
        count = curs.fetchone()[0]
        assert count == 1

    def test_ligand_write(self, dbquery):
        curs = dbquery("""SELECT COUNT(*) FROM Results;""" )
        count = curs.fetchone()[0]
        assert count == 645

    def test_one_filter(self):
        count_ligands_passing = self.rtstorage.filter(eworst = -6)
        assert count_ligands_passing == 65
        os.system("rm output_log.txt ")

    def test_two_filters(self):
        count_ligands_passing = self.rtstorage.filter(eworst = -6, hb_interactions= [('A:VAL:279:', True), ('A:LYS:162:', True)])
        assert count_ligands_passing == 18
        os.system("rm output_log.txt")

    def test_three_filters(self, dbquery):
        count_ligands_passing = self.rtstorage.filter(eworst = -6, hb_interactions=[('A:VAL:279:', True), ('A:LYS:162:', True)], vdw_interactions=[('A:VAL:279:', True), ('A:LYS:162:', True)], max_miss = 1)
        os.system("rm outputapi.db output_log.txt")
        assert count_ligands_passing == 51 


        