#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Shared pytest fixtures for Ringtail test suite
#
# The expensive thing in this suite is ingesting docking results: one 217-ligand adgpu
# import takes ~7 s, and 40-odd tests want the same database. So it is built **once per
# backend per session** and each test gets a private copy of the file — measured at 3 ms
# against 6.7 s, i.e. ~2000x cheaper, for exactly the same isolation.
#
import shutil
import subprocess
from pathlib import Path
from types import SimpleNamespace

import pytest
from ringtail import RingtailCore

TEST_DIR = Path(__file__).parent
TEST_DATA = TEST_DIR / "test_data"

TEST_MAX_PROC = 3


@pytest.fixture(autouse=True, scope="session")
def _limit_ingest_workers():
    """Cap ingest parallelism for the whole session."""
    from _pytest.monkeypatch import MonkeyPatch

    import ringtail.mpmanager as mpmanager

    patch = MonkeyPatch()
    patch.setattr(mpmanager.mp, "cpu_count", lambda: TEST_MAX_PROC)
    yield
    patch.undo()


#: The adgpu directories the populated fixtures ingest, and the ligand count they produce.
POPULATED_SOURCES = ["adgpu/group1", "adgpu/group2"]
POPULATED_LIGANDS = 217


# Session-scoped so the template below can depend on it: a fixture may depend on one of
# equal or wider scope, never narrower. Pytest also groups tests by parameter, so each
# backend's template is built once rather than once per switch.
@pytest.fixture(scope="session", params=["duckdb", "sqlite"])
def storage_type(request):
    return request.param


@pytest.fixture
def tmp_db(tmp_path, storage_type) -> RingtailCore:
    """Fresh empty RingtailCore db for each test, parametrized over storage backends."""
    return RingtailCore(db_file=str(tmp_path / "test.db"), storage_type=storage_type)


@pytest.fixture(scope="session")
def _populated_template(tmp_path_factory, storage_type) -> Path:
    """The ingested database, built once per backend per session. Never handed to a test.

    ``add_results_from_files`` closes its connection before returning, which matters here:
    DuckDB keeps a single-writer lock and a write-ahead log beside the file, so copying
    while a connection is open could hand a test a locked or half-written database.
    """
    path = tmp_path_factory.mktemp(f"template_{storage_type}") / "populated.db"
    rtc = RingtailCore(str(path), storage_type=storage_type)
    rtc.add_results_from_files(
        docking_results=[str(TEST_DATA / src) for src in POPULATED_SOURCES],
        docking_mode="adgpu",
    )
    assert rtc.table_length("Ligands") == POPULATED_LIGANDS, (
        "the shared template is not what tests expect; every populated_db test will be "
        "wrong in the same way, so fail loudly here instead"
    )
    return path


@pytest.fixture
def populated_db(tmp_path, storage_type, _populated_template) -> RingtailCore:
    """A private 217-ligand adgpu db (group1 + group2) for each test.

    Copied from the session template rather than re-ingested. Each test still gets its own
    file, so tests may freely write bookmarks, statuses and comments.
    """
    dest = tmp_path / "populated.db"
    shutil.copy2(_populated_template, dest)
    return RingtailCore(str(dest), storage_type=storage_type)


def _template(tmp_path_factory, storage_type, label, **ingest) -> Path:
    path = tmp_path_factory.mktemp(f"{label}_{storage_type}") / f"{label}.db"
    rtc = RingtailCore(str(path), storage_type=storage_type)
    rtc.add_results_from_files(**ingest)
    return path


def _copy_of(template: Path, tmp_path, storage_type, name) -> RingtailCore:
    dest = tmp_path / f"{name}.db"
    shutil.copy2(template, dest)
    return RingtailCore(str(dest), storage_type=storage_type)


@pytest.fixture(scope="session")
def _vina_template(tmp_path_factory, storage_type) -> Path:
    vina_path = TEST_DATA / "vina"
    return _template(
        tmp_path_factory,
        storage_type,
        "vina",
        docking_results=str(vina_path),
        receptor_file=str(vina_path / "receptor.pdbqt"),
        save_receptor=True,
        docking_mode="vina",
    )


@pytest.fixture
def vina_db(tmp_path, storage_type, _vina_template) -> RingtailCore:
    """The 6-pose vina dataset with its receptor and interactions, as a private copy."""
    return _copy_of(_vina_template, tmp_path, storage_type, "vina")


@pytest.fixture(scope="session")
def _ad6_template(tmp_path_factory, storage_type) -> Path:
    ad6_path = TEST_DATA / "ad6"
    return _template(
        tmp_path_factory,
        storage_type,
        "ad6",
        docking_results=str(ad6_path),
        receptor_file=str(ad6_path / "helix--scofu01.json"),
        save_receptor=True,
        docking_mode="ad6",
    )


@pytest.fixture
def ad6_db(tmp_path, storage_type, _ad6_template) -> RingtailCore:
    """The 9-pose ad6 dataset with its receptor and interactions, as a private copy."""
    return _copy_of(_ad6_template, tmp_path, storage_type, "ad6")


@pytest.fixture(scope="session")
def _ad6_no_interactions_template(tmp_path_factory, storage_type) -> Path:
    return _template(
        tmp_path_factory,
        storage_type,
        "ad6_plain",
        docking_results=str(TEST_DATA / "ad6"),
        calculate_interactions=False,
        docking_mode="ad6",
    )


@pytest.fixture
def ad6_db_no_interactions(
    tmp_path, storage_type, _ad6_no_interactions_template
) -> RingtailCore:
    """The ad6 dataset ingested with interactions switched off, as a private copy."""
    return _copy_of(_ad6_no_interactions_template, tmp_path, storage_type, "ad6_plain")


def _flexres_template(tmp_path_factory, storage_type, receptor: str) -> Path:
    data_path = TEST_DATA / "flexres"
    return _template(
        tmp_path_factory,
        storage_type,
        f"flexres_{Path(receptor).suffix.lstrip('.')}",
        docking_results=str(data_path / "ligand.pdbqt"),
        docking_mode="vina",
        receptor_file=str(data_path / receptor),
        save_receptor=True,
    )


@pytest.fixture(scope="session")
def _flexres_json_template(tmp_path_factory, storage_type) -> Path:
    return _flexres_template(tmp_path_factory, storage_type, "receptor.json")


@pytest.fixture(scope="session")
def _flexres_pdbqt_template(tmp_path_factory, storage_type) -> Path:
    return _flexres_template(tmp_path_factory, storage_type, "receptor.pdbqt")


@pytest.fixture
def flexres_json_db(tmp_path, storage_type, _flexres_json_template) -> RingtailCore:
    """One flexible-residue pose, receptor supplied as meeko Polymer JSON."""
    return _copy_of(_flexres_json_template, tmp_path, storage_type, "flexres_json")


@pytest.fixture
def flexres_pdbqt_db(tmp_path, storage_type, _flexres_pdbqt_template) -> RingtailCore:
    """The same pose, receptor supplied as pdbqt. Interactions must come out identical to
    the JSON route, which is what makes the pair worth having."""
    return _copy_of(_flexres_pdbqt_template, tmp_path, storage_type, "flexres_pdbqt")


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

    Note these shell out to the installed console script, so the tests using them need
    Ringtail's entry points on PATH. They are marked ``slow`` for that reason: a subprocess
    per call, each doing its own import and ingest.
    """
    db = tmp_path / "output.db"

    def write(*args):
        # --max_proc for the same reason as TEST_MAX_PROC above
        cmd = [
            "rt_process_vs",
            "write",
            "-o",
            str(db),
            "-st",
            "duckdb",
            "--max_proc",
            str(TEST_MAX_PROC),
            *args,
        ]
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
