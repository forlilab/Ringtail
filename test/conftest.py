#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Shared pytest fixtures for Ringtail test suite
#
# Populated dbs are built once per backend per session, and each test gets its own copy
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


# Session scope so the session templates can use it
@pytest.fixture(
    scope="session",
    # quick run (-m "not slow") only tests duckdb
    params=["duckdb", pytest.param("sqlite", marks=pytest.mark.slow)],
)
def storage_type(request):
    return request.param


@pytest.fixture
def tmp_db(tmp_path, storage_type) -> RingtailCore:
    """Fresh empty RingtailCore db for each test, parametrized over storage backends."""
    return RingtailCore(db_file=str(tmp_path / "test.db"), storage_type=storage_type)


@pytest.fixture(scope="session")
def _populated_template(tmp_path_factory, storage_type) -> Path:
    """Populated db built once per backend, only used as a template for populated_db"""
    return _template(
        tmp_path_factory,
        storage_type,
        "populated",
        docking_results=[str(TEST_DATA / "adgpu/group1"), str(TEST_DATA / "adgpu/group2")],
        docking_mode="adgpu",
    )


@pytest.fixture
def populated_db(tmp_path, storage_type, _populated_template) -> RingtailCore:
    """Copy of the 217-ligand adgpu db (group1 + group2) for each test"""
    return _copy_of(_populated_template, tmp_path, storage_type, "populated")


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
    """Copy of the vina db (6 poses) with receptor and interactions"""
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
    """Copy of the ad6 db (9 poses) with receptor and interactions"""
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
    """Copy of the ad6 db without interactions"""
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
    """Copy of the flexres db (1 pose) with receptor from meeko Polymer JSON"""
    return _copy_of(_flexres_json_template, tmp_path, storage_type, "flexres_json")


@pytest.fixture
def flexres_pdbqt_db(tmp_path, storage_type, _flexres_pdbqt_template) -> RingtailCore:
    """Copy of the flexres db with receptor from pdbqt, same interactions as the JSON one"""
    return _copy_of(_flexres_pdbqt_template, tmp_path, storage_type, "flexres_pdbqt")


@pytest.fixture
def make_ad6_sdf(tmp_path):
    from rdkit import Chem

    counter = 0

    def make(
        smiles="OCC(N)(CO)CO",
        ligname="test_ligand",
        docking_score=-10,
        pose_rank=1,
        custom=None,
    ):
        nonlocal counter
        mol = Chem.MolFromSmiles(smiles)

        mol.AddConformer(Chem.Conformer(mol.GetNumAtoms()))

        for name, value in (custom or {}).items():
            if value is None:
                continue
            if isinstance(value, int):
                mol.SetIntProp(name, value)
            elif isinstance(value, float):
                mol.SetDoubleProp(name, value)
            else:
                mol.SetProp(name, str(value))

        mol.SetProp("_Name", ligname)
        mol.SetDoubleProp("autodock_free_energy", docking_score)
        mol.SetIntProp("pose_id", pose_rank)

        path = tmp_path / f"ligand_{counter}.sdf"
        counter += 1
        with Chem.SDWriter(str(path)) as writer:
            writer.write(mol)

        return str(path)

    return make


@pytest.fixture
def cli(tmp_path):
    """Helpers that run the rt_process_vs CLI with a duckdb output db in tmp_path

    Returns:
        SimpleNamespace: write(*args), read(*args), count(table, distinct_column), db
    """
    db = tmp_path / "output.db"

    def write(*args):
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
        # no stdin so a consent prompt is refused instead of hanging
        return subprocess.run(
            cmd, cwd=str(TEST_DIR), capture_output=True, stdin=subprocess.DEVNULL
        ).returncode

    def read(*args):
        cmd = ["rt_process_vs", "read", "--input_db", str(db), *args]
        return subprocess.run(cmd, cwd=str(tmp_path), capture_output=True).returncode

    def count(table, distinct_column=None):
        return RingtailCore(str(db)).table_length(table, distinct_column)

    return SimpleNamespace(write=write, read=read, count=count, db=db)
