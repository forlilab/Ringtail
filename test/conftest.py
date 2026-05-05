#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Shared pytest fixtures for Ringtail test suite
#
import subprocess
from pathlib import Path
from types import SimpleNamespace

import pytest
from ringtail import RingtailCore

TEST_DIR = Path(__file__).parent
TEST_DATA = TEST_DIR / "test_data"


@pytest.fixture(params=["duckdb", "sqlite"])
def storage_type(request):
    return request.param


@pytest.fixture
def tmp_db(tmp_path, storage_type) -> RingtailCore:
    """Fresh empty RingtailCore db for each test, parametrized over storage backends."""
    return RingtailCore(db_file=str(tmp_path / "test.db"), storage_type=storage_type)


@pytest.fixture
def populated_db(tmp_path, storage_type):
    """Fresh 217-ligand adgpu db (group1 + group2) for each test."""
    rtc = RingtailCore(str(tmp_path / "populated.db"), storage_type=storage_type)
    rtc.add_results_from_files(
        file_path=[str(TEST_DATA / "adgpu/group1"), str(TEST_DATA / "adgpu/group2")],
        docking_mode="adgpu",
    )
    return rtc


@pytest.fixture
def cli(tmp_path):
    """
    Fixture providing write/read helpers that route CLI output to tmp_path.
    Always uses duckdb (unit tests cover both backends via storage_type).
    Returns a SimpleNamespace with:
      .write(*args)  -> returncode
      .read(*args)   -> returncode
      .count(table)  -> int (row count in table)
      .passing(bookmark) -> int (passing pose count)
      .db            -> Path to the output db
    """
    db = tmp_path / "output.db"

    def write(*args):
        cmd = ["rt_process_vs", "write", "-o", str(db), "-st", "duckdb", *args]
        return subprocess.run(cmd, cwd=str(TEST_DIR), capture_output=True).returncode

    def read(*args):
        cmd = ["rt_process_vs", "read", "--input_db", str(db), *args]
        return subprocess.run(cmd, cwd=str(tmp_path), capture_output=True).returncode

    def count(table):
        return RingtailCore(str(db)).table_length(table)

    def passing(bookmark):
        rtc = RingtailCore(str(db))
        with rtc.storageman as sm:
            return sm.get_passing_poses_count(bookmark, True)

    return SimpleNamespace(write=write, read=read, count=count, passing=passing, db=db)
