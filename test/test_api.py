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
    # Setup
    rtstorage = RingtailCore(db_file="outputapi.db")

    rtstorage.set_general_options(summary=False, debug=True)
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
        self.rtstorage.set_filters(eworst = -6)
        count_ligands_passing = self.rtstorage.filter()
        assert count_ligands_passing == 65
        os.system("rm output_log.txt ")

    def test_two_filters(self):
        self.rtstorage.set_filters(eworst = -6, hb_interactions= [('A:VAL:279:', True), ('A:LYS:162:', True)])
        count_ligands_passing = self.rtstorage.filter()
        assert count_ligands_passing == 18
        os.system("rm output_log.txt")

    def test_three_filters(self, dbquery):
        self.rtstorage.set_filters(eworst = -6, hb_interactions=[('A:VAL:279:', True), ('A:LYS:162:', True)], vdw_interactions=[('A:VAL:279:', True), ('A:LYS:162:', True)], max_miss = 1)
        count_ligands_passing = self.rtstorage.filter()
        assert count_ligands_passing == 51 
        os.system("rm outputapi.db output_log.txt")

class TestConfigFile:
    # Setup
    rtcore = RingtailCore(db_file="outputapi.db")
    rtcore.generate_config_json_template()
    filepath = rtcore._config_file_path()

    with open(filepath, "r") as f:
        data = json.load(f)
    # all fields I want to change
    data["fileobj"]["file_path"] = [['test_data/']]
    data["fileobj"]["file_path"]
    data["fileobj"]["recursive"] = True
    data["fileobj"]["receptor_file"] = "test_data/4j8m.pdbqt"
    data["fileobj"]["save_receptor"] = True
    data["filterobj"]["eworst"] = -6
    data["filterobj"]["vdw_interactions"] = [('A:VAL:279:', True), ('A:LYS:162:', True)]
    data["filterobj"]["hb_interactions"] = [('A:VAL:279:', True), ('A:LYS:162:', True)]
    data["filterobj"]["max_miss"] = 1

    with open(filepath, "w") as f:
        f.write(json.dumps(data, indent=4))
    rtcore.add_config_from_file()

    def test_adding_results(self, dbquery):
        self.rtcore.add_results_from_files(file_source_object= self.rtcore.files)
        curs = dbquery("""SELECT COUNT(*) FROM Results;""" )
        count = curs.fetchone()[0]

        assert count == 645
        
    def test_filter(self):
        count_ligands_passing = self.rtcore.filter()

        assert count_ligands_passing == 51 
        os.system("rm outputapi.db")

class TestOptionsHandling:
    pass
    # Test the order of operations, single options overwrite passing through a dict

    # Test that type is checked even if you set a value later on

    # Test that you can set an individual value on an object after it has been set
        # Using the method of course, but type must be checked. Can test this by seeing what is written to the output log or something? 

    # Test that any values in an options object can be accessed without error (and issues with _ referencing etc)

    # Test that if you set an invalid option (including None) it reverts to the default value
        # So I really need to have type and default and name in one space
