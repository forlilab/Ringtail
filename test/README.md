# Ringtail test suite

## Running the tests

```bash
pytest -m "not slow"     # quick run, duckdb only, no CLI tests
pytest                   # full run, adds sqlite and the CLI tests
```

The CLI tests call `rt_process_vs`, so run pytest with your Ringtail environment activated.

```bash
conda activate <your_ringtail_env>
pytest
```

The 10 slowest tests are printed at the end of each run.

## Running specific tests

From the Ringtail folder, add a file, class or test to the path with `::`.

```bash
pytest test/test_units.py
pytest test/test_units.py::TestFiltering
pytest test/test_units.py::TestCoreOperations::test_append_to_database
```

Most tests run once per storage backend, and the backend is part of the test id. Pick one with `-k`, or add the id in brackets.

```bash
pytest test/test_units.py::TestFiltering -k duckdb
pytest test/test_units.py -k "TestAD6Handling and sqlite"
pytest "test/test_units.py::TestCoreOperations::test_append_to_database[sqlite]"
```

`-k` also matches test and class names, so you can run every test with a given word in the name.

```bash
pytest -k bookmark
pytest -m "not slow" -k TestCustomFields
```

The sqlite cases are marked slow, so `-m "not slow"` will skip them even if you ask for sqlite with `-k`.

Other useful options:

| option | what it does |
|---|---|
| `-x` | stop at the first failure |
| `-v` | print each test id and its result |
| `-s` | show print output from the tests |
| `--lf` | rerun only the tests that failed last time |
| `--collect-only -q` | list the test ids without running them |

## Layout

| file | what it covers |
|---|---|
| `conftest.py` | shared fixtures, see below |
| `test_units.py` | the `RingtailCore` API, for both storage backends |
| `test_cmdline.py` | the command line scripts, all marked slow |
| `test_data/` | adgpu, vina and ad6 docking results, receptors, file lists |

## Fixtures

Every fixture that uses a database runs once for duckdb and once for sqlite, through `storage_type`.

| fixture | what you get |
|---|---|
| `tmp_db` | empty db |
| `populated_db` | 217 adgpu ligands (group1 + group2) |
| `vina_db` | 6 vina poses with receptor and interactions |
| `ad6_db` | 9 ad6 poses with receptor and interactions |
| `ad6_db_no_interactions` | the ad6 poses without interactions |
| `flexres_json_db` / `flexres_pdbqt_db` | 1 flexible-residue pose, receptor as meeko Polymer JSON or pdbqt |
| `make_ad6_sdf` | writes an ad6 sdf for one ligand, with optional custom properties |
| `cli` | runs `rt_process_vs` write and read on a duckdb output db |

The populated databases are built once per backend per session, and each test gets its own copy. To add a new one, build it with `_template` and hand out copies with `_copy_of` in `conftest.py`. Tests about adding results should use `tmp_db` and add the results in the test.

Ingest is capped at `TEST_MAX_PROC` processes for the whole session.

## Slow tests

Tests marked slow are left out of the quick run:

- all of `test_cmdline.py`
- `TestMergeDB` and `TestCrossref`
- `test_custom_fields_merge` in `TestCustomFields`
- every sqlite case

New tests are in the quick run unless you mark them with `@pytest.mark.slow`. A mistyped marker is an error, since `--strict-markers` is set in `pyproject.toml`.

The ingest tests in `TestCoreOperations`, `TestVinaHandling` and `TestAD6Handling` are not marked slow, so the quick run still covers adding results.
