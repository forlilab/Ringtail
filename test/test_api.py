import sqlite3
import os
import pytest
from ringtail import RingtailCore, RingtailArguments as RTArgs, APIOptionParser as api

@pytest.fixture(scope='function')
def cur():
    conn = sqlite3.connect("output.db")
    curs = conn.cursor()
    yield curs
    curs.close()
    conn.close()
    os.system("rm output.db")

def test_new_api(cur):
    
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

    query = """SELECT COUNT(*) FROM Receptors WHERE receptor_object NOT NULL;""" 
    cur.execute(query)
    count = cur.fetchone()[0]

    assert count == 1

def test_new_api2(cur):
    
    rtopt = RTArgs()
    rtopt.process_mode = "write"
    rtopt.file_path = [['test_data/group1/']]
    rtopt.recursive = True
    rtopt.save_receptor = False
    rtopt.receptor_file = "test_data/4j8m.pdbqt"

    rtcore = RingtailCore()
    api(ringtail_core=rtcore, opts =rtopt)
    with rtcore: rtcore.add_results()

    query = """SELECT COUNT(*) FROM Receptors WHERE receptor_object NOT NULL;""" 
    cur.execute(query)
    count = cur.fetchone()[0]

    assert count == 0

def test_new_api3(cur):
    
    rtopt = RTArgs()
    rtopt.process_mode = "write"
    rtopt.file_path = [['test_data/group1/']]
    rtopt.recursive = True
    rtopt.save_receptor = False
    rtopt.receptor_file = "test_data/4j8m.pdbqt"

    rtcore = RingtailCore()
    api(ringtail_core=rtcore, opts =rtopt)
    with rtcore: rtcore.add_results()

    query = """SELECT COUNT(*) FROM Results;""" 
    cur.execute(query)
    count = cur.fetchone()[0]

    assert count == 370
