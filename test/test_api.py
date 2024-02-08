import sqlite3
import os
import pytest
from ringtail import RingtailCore 

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

def write_standard_setup():
    rtstorage = RingtailCore(db_file="outputapi.db")
    rtstorage.set_general_options(process_mode="write", summary=False, debug=True)
    rtstorage.add_results_from_files(file_path = [['test_data/']],
                                     recursive = True,
                                     receptor_file="test_data/4j8m.pdbqt",
                                     save_receptor=True)
    return rtstorage

class TestRingtailWrite:
    rtstorage = write_standard_setup()

    def test_receptor_save(self, dbquery):
        curs = dbquery("""SELECT COUNT(*) FROM Receptors WHERE receptor_object NOT NULL;""" )
        count = curs.fetchone()[0]
        assert count == 1

    def test_ligand_write(self, dbquery):
        curs = dbquery("""SELECT COUNT(*) FROM Results;""" )
        count = curs.fetchone()[0]
        assert count == 645

    def test_one_filter(self):
        self.rtstorage.set_read_options(log_file="outputapi1_log.txt")
        self.rtstorage.set_filters(eworst = -6)
        count_ligands_passing = self.rtstorage.filter()
        assert count_ligands_passing == 65
        os.system("rm outputapi1_log.txt ")

    def test_two_filters(self):
        self.rtstorage.set_read_options(log_file="outputapi2_log.txt")
        self.rtstorage.set_filters(eworst = -6, hb_interactions= [('A:VAL:279:', True), ('A:LYS:162:', True)])
        count_ligands_passing = self.rtstorage.filter()
        assert count_ligands_passing == 18
        os.system("rm outputapi2_log.txt")

    def test_three_filters(self, dbquery):
        self.rtstorage.set_read_options(log_file="outputapi3_log.txt")
        #self.rtstorage.set_storage_options(results_view_name="multiplechoice")
        self.rtstorage.set_filters(eworst = -6, hb_interactions=[('A:VAL:279:', True), ('A:LYS:162:', True)], vdw_interactions=[('A:VAL:279:', True), ('A:LYS:162:', True)], max_miss = 1)
        count_ligands_passing = self.rtstorage.filter()

        assert count_ligands_passing == 51 
        os.system("rm outputapi.db outputapi3_log.txt")

        
