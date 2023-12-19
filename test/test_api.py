import sqlite3
import os
import pytest
from ringtail import RingtailCore as ringtail

if os.path.isfile("output.db"): os.system("rm output.db")




@pytest.fixture
def cur():
    conn = sqlite3.connect("output.db")
    curs = conn.cursor()
    yield curs
    curs.close()
    conn.close()
    # os.system("rm output.db")

@pytest.fixture
def countrows():
    def __dbconnect(query):
        conn = sqlite3.connect("output.db")
        curs = conn.cursor()
        curs.execute(query)
        count = curs.fetchone()[0]
        curs.close()
        conn.close()
        # os.system("rm output.db")
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
    # os.system("rm output.db")


def write_standard_setup():
    print("\nSetting up Ringtail Core\n")
    rtstorage = ringtail(db_file="output.db")
    rtstorage.open()
    rtstorage.set_general_options(process_mode="write", summary=True)
    # rtstorage.file_writer_options(max_poses=2)
    rtstorage.add_results_from_files(file_path = [['test_data/group2/']],
                                     recursive = True,
                                     receptor_file="test_data/4j8m.pdbqt")

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
        assert count == 176 #242
