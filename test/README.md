# Ringtail test suite

## Running

```bash
pytest -m "not slow"     # quick: 106 cases, ~2m20s — API level only
pytest                   # deep:  175 cases, ~9m20s — adds the CLI integration tests
```

`pytest` with no arguments is the full run, on purpose — see "Why `slow` and not `quick`"
below.

The deep run shells out to Ringtail's console scripts, so **the entry points have to be on
`PATH`**. Inside the environment that is automatic; if you invoke pytest by absolute path
without activating the env, every CLI test errors with
`FileNotFoundError: 'rt_process_vs'`.

```bash
micromamba activate saltyringtail && pytest      # entry points on PATH
```

## Layout

| File | What it covers |
|---|---|
| `conftest.py` | shared fixtures: `tmp_db`, `populated_db`, `cli`, `storage_type` |
| `test_units.py` | the `RingtailCore` API, both storage backends |
| `test_cmdline.py` | the same ground through `rt_process_vs` and friends — all `slow` |
| `test_data/` | adgpu, vina and ad6 docking results, receptors, file lists |

Every fixture that touches a database is parametrized over **duckdb and sqlite**, so a
dialect-specific bug in one backend fails the suite rather than waiting for a user to find
it.

## Why `slow` and not `quick`

Marking the *expensive* tests, and selecting with `-m "not slow"`, means a newly added test
is in the quick run **by default**. The inverse — marking the essential ones `quick` — looks
equivalent but fails in the worse direction: a new test is silently excluded until someone
remembers to mark it, so the quick run slowly stops meaning anything, invisibly. Here the
worst case is a quick run that gets slower, which you notice immediately.

`--strict-markers` is set in `pyproject.toml` so a typo'd marker is an error rather than a
decorator that does nothing.

## What is `slow`

- **all of `test_cmdline.py`** — a subprocess per call, each with its own import and ingest,
  re-covering through the CLI what `test_units.py` covers through the API.
- **`TestMergeDB`, `TestCrossref`** — clone databases and attach extra ones.

Deliberately *not* marked slow: the ingest tests in `TestCoreOperations`,
`TestVinaHandling` and `TestAD6Handling`. They do cost seconds each, but the write path is
where the bugs are, and a quick run with no ingest coverage would be false comfort.

## The expensive thing is ingest

Two changes account for the 17 minutes becoming 9, and they are worth understanding before
adding tests.

**1. Shared templates.** Many tests ingest byte-identical data and differ only in what they
assert. Those datasets are built **once per backend per session** and each test gets a copy
of the file — measured at 3 ms against 6.7 s, roughly 2000x cheaper, with identical
isolation, since every test still gets a private file it can write bookmarks into.

| fixture | dataset |
|---|---|
| `populated_db` | 217-ligand adgpu (group1 + group2) |
| `vina_db` | 6-pose vina with receptor and interactions |
| `ad6_db` | 9-pose ad6 with receptor and interactions |
| `ad6_db_no_interactions` | the same, ingested with `calculate_interactions=False` |
| `flexres_json_db` / `flexres_pdbqt_db` | one flexible-residue pose, receptor as Polymer JSON vs pdbqt |

If you add one, note the ordering constraint: a session-scoped fixture may only depend on
fixtures of session scope or wider, which is why `storage_type` is session-scoped. And let
the connection close before copying — DuckDB holds a single-writer lock and keeps a
write-ahead log beside the file.

Tests that are *about* ingest still ingest for real. `test_add_folder` and
`test_append_to_database` would test nothing if handed a prebuilt database.

**2. Capped ingest parallelism** (`TEST_MAX_PROC`, and `--max_proc` for the CLI helper). See
the gotcha below — this one turned out to be the larger of the two, taking the quick run
from 4m06 to 2m22 on its own.

## Gotchas

**Ingest cannot run from a piped script.** Ringtail forces the `spawn` start method, and
spawned children re-import `__main__`; from `python - <<EOF` that is `<stdin>`, which does
not exist, and you get a wall of `FileNotFoundError` from the child processes. Put the code
in a real file under `if __name__ == "__main__":`.

**`max_proc` defaults to `cpu_count()` regardless of input size.** A 3-file ingest measured
4.67 s with the default 11 readers against 1.93 s with one reader — the same work, the
difference being 11 interpreter startups under the forced `spawn` method. `conftest.py` caps
this session-wide via `TEST_MAX_PROC`, and passes `--max_proc` to the CLI helper because a
subprocess cannot see the in-process patch.

This is worth knowing beyond the tests: a user importing a handful of files pays the same
overhead, so scaling the default by input size would be a real improvement to Ringtail
itself.

**`--durations=10` is on by default** so the slowest tests are always printed. If that list
starts growing, something has regressed.
