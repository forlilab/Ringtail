import sqlite3
import os
import pytest
from ringtail import RingtailCore, RTArgs, RTWrite, RTRead, APIOptionParser as api

if os.path.isfile("output.db"): os.system("rm output.db")

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
    os.system("rm output.db")

def write_standard_setup():
    print("\nSetting up Ringtail Core\n")
    rtopt = RTWrite()
    rtopt.process_mode = "write"
    rtopt.file_path = [['test_data/group1/']]
    rtopt.recursive = True
    rtopt.save_receptor = True
    rtopt.receptor_file = "test_data/4j8m.pdbqt"
    rtopt.summary = True
    rtopt.max_poses = 2
    rtcore = RingtailCore()
    api(ringtail_core=rtcore, opts =rtopt)
    with rtcore: rtcore.add_results()
    with rtcore: rtcore.save_receptors(rtopt.receptor_file)
    print("\n\nRingtail Core set up\n")

write_standard_setup()
class TestRingtailWrite:

    def test_receptor_save(self, dbquery):
        curs = dbquery("""SELECT COUNT(*) FROM Receptors WHERE receptor_object NOT NULL;""" )
        count = curs.fetchone()[0]
        assert count == 1

    def test_ligand_write(self, dbquery):
        curs = dbquery("""SELECT COUNT(*) FROM Results;""" )
        count = curs.fetchone()[0]
        assert count == 242

# write_standard_setup()
# class TestRingtailRead:
#     print("\nSetting up Ringtail Core\n")
#     rtread = RTWrite()
#     rtread.process_mode = "read"
    
#     rtread.summary = True
    
#     rtcore = RingtailCore()
#     api(ringtail_core=rtcore, opts =rtread)
#     print("\n\nRingtail Core set up\n")
