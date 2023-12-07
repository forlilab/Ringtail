import sqlite3
import os
import pytest
from ringtail import RingtailCore, RingtailArguments as RTArgs, APIOptionParser as api

@pytest.fixture
def cur():
    conn = sqlite3.connect("output.db")
    curs = conn.cursor()
    yield curs
    curs.close()
    conn.close()
    os.system("rm output.db")

@pytest.fixture
def countrows():
    def __dbconnect(query):
        conn = sqlite3.connect("output.db")
        curs = conn.cursor()
        curs.execute(query)
        count = curs.fetchone()[0]
        curs.close()
        conn.close()
        os.system("rm output.db")
        return count
    return __dbconnect

def test_new_api(countrows):
    
    rtopt = RTArgs()
    rtopt.process_mode = "write"
    rtopt.file_path = [['test_data/group1/']]
    rtopt.recursive = True
    rtopt.save_receptor = True
    rtopt.receptor_file = "test_data/4j8m.pdbqt"

    rtcore = RingtailCore()
    api(ringtail_core=rtcore, opts =rtopt)
    with rtcore: rtcore.add_results()
    with rtcore: rtcore.save_receptors(rtopt.receptor_file)

    count = countrows("""SELECT COUNT(*) FROM Receptors WHERE receptor_object NOT NULL;""" )

    assert count == 1

def test_new_api2(countrows):
    
    rtopt = RTArgs()
    rtopt.process_mode = "write"
    rtopt.file_path = [['test_data/group1/']]
    rtopt.recursive = True
    rtopt.save_receptor = False
    rtopt.receptor_file = "test_data/4j8m.pdbqt"

    rtcore = RingtailCore()
    api(ringtail_core=rtcore, opts =rtopt)
    with rtcore: rtcore.add_results()

    count = countrows("""SELECT COUNT(*) FROM Receptors WHERE receptor_object NOT NULL;""" )
    assert count == 0

def test_new_api3(countrows):
    
    rtopt = RTArgs()
    rtopt.process_mode = "write"
    rtopt.file_path = [['test_data/group1/']]
    rtopt.recursive = True
    rtopt.save_receptor = False
    rtopt.receptor_file = "test_data/4j8m.pdbqt"

    rtcore = RingtailCore()
    api(ringtail_core=rtcore, opts =rtopt)
    with rtcore: rtcore.add_results()

    count = countrows("""SELECT COUNT(*) FROM Results;""" )
    assert count == 370
