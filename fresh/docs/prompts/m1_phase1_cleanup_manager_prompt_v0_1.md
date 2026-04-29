# Claude M1 Phase 1 Prompt — Cleanup Manager v0.1

You are working in the `fresh/` workflow area of the EGFR-MYO1D repository on branch `claude/task10` (or successor). This is **M1 Phase 1** of the M1 completion plan stored at `C:\Users\admin\.claude\plans\1-enchanted-pumpkin.md`.

## 1. Project context

The fresh EGFR-MYO1D workflow is being completed strictly to the `milestone1_foundation_codex_handoff_v0_5.md` spec. Tasks 1-9 already exist; this phase fills in a missing M1 module (cleanup manager).

This phase must **not** perform docking, pocket discovery, receptor normalization, MYO1D slicing, ligand prep, scoring, or candidate nomination. It must not run qsub.

## 2. Absolute rules

Do not modify the old workflow:

```text
run_production.py
main.py
egfr_pipeline/**
config/**
docs/runbook.md
output/**
results_export/**
```

All new code must be under `fresh/`. All run outputs must stay under `fresh/runs/<run_id>/`. Cleanup must NEVER delete outside the active run directory.

Maintain Python 2.7.11 / 3.9 syntax compatibility per existing project pattern (commit 40336dd). No f-strings, use `from __future__ import` headers, `.format()` for strings, `# type: (...) -> ...` comments instead of inline annotations where existing modules do.

## 3. Scope

In scope:
- Implement `fresh/src/egfr_myo1d/core/cleanup.py` per handoff §18
- Add `cleanup` CLI subcommand in `fresh/src/egfr_myo1d/cli.py` following `_cmd_init_run` pattern
- Replace `fresh/scripts/cleanup_run.py` placeholder with a real entrypoint that delegates to `core/cleanup.py`
- Add tests under `fresh/tests/test_m1_phase1_cleanup_manager.py`
- Add fixtures under `fresh/tests/fixtures/m1_phase1_cleanup/`
- Add docs `fresh/docs/m1_phase1_cleanup_manager.md` and `fresh/docs/m1_phase1_changes.md`

Out of scope:
- PBS generation (Phase 6)
- Receptor/MYO1D/ligand prep (Phases 3-5, 7)
- M2 docking, pocket discovery
- Reading `gates.yaml` cleanup_policy beyond the test_policy / production_cleanup_default keys (already in config)

## 4. Required CLI behavior

Add subcommand `cleanup` to `fresh/src/egfr_myo1d/cli.py`:

```bash
python -m egfr_myo1d.cli cleanup --run-id RUN --mode test [--dry-run true|false]
python -m egfr_myo1d.cli cleanup --run-id RUN --mode production [--dry-run true|false]
```

Behavior:

- `--run-id`: required; resolved via `RunContext.for_existing(run_id)`. Path-traversal IDs (containing `..`, `/`, `\`) are rejected with nonzero exit.
- `--mode`:
  - `test`: deletes intermediate/scratch files (see §6 below)
  - `production`: defaults to dry-run unless `--dry-run false` is explicitly passed
- `--dry-run`: accepts `true` or `false`. For `--mode test`, default is `false`. For `--mode production`, default is `true`. Explicit `--dry-run` overrides defaults.
- Process exit: `0` for PASS or PASS_WITH_WARNINGS, `1` for FAIL.
- Stdout: short summary printing run_id, mode, dry_run, deleted_count, preserved_count, cleanup_report path.

## 5. Files to create / modify

Create:

```text
fresh/src/egfr_myo1d/core/cleanup.py
fresh/tests/test_m1_phase1_cleanup_manager.py
fresh/tests/fixtures/m1_phase1_cleanup/sample_run_layout.txt   # describes synthetic run dir layout used by tests
fresh/docs/m1_phase1_cleanup_manager.md
fresh/docs/m1_phase1_changes.md
```

Modify:

```text
fresh/src/egfr_myo1d/cli.py            # add cleanup subparser + _cmd_cleanup handler
fresh/scripts/cleanup_run.py           # replace placeholder; delegate to cli `cleanup` or directly to core/cleanup.py
```

Do not modify any old workflow files.

## 6. Cleanup deletion / preservation policy (handoff §18)

Deletion candidates (only inside the active run_dir):

```text
scratch/
tmp/
*.tmp
*.pdbqt.tmp
*.vina.tmp
*.silent
per-pose scratch files
test docking intermediates
```

Preservation list (must NEVER be deleted):

```text
manifest/                          # all manifest files
logs/                              # master.log, phase_status.jsonl, job_status.jsonl, jobs/, errors/
qc/                                # all QC CSVs
reports/                           # all reports
manifest/cleanup_report.json       # written by this command
```

Safety rules:

```text
- Cleanup must refuse to delete outside ctx.run_dir (use ctx.require_within_run_dir for every candidate).
- Cleanup must support dry-run (no deletions, but still write cleanup_report.json with what would be deleted).
- Cleanup must always write manifest/cleanup_report.json.
- Production cleanup default is dry_run=true.
- Test cleanup default is dry_run=false.
```

## 7. Required output files

After `cleanup` runs:

```text
fresh/runs/<run_id>/manifest/cleanup_report.json
fresh/runs/<run_id>/logs/phase_status.jsonl                # appended (PASS or WARN)
fresh/runs/<run_id>/logs/master.log                        # appended
```

`cleanup_report.json` schema:

```json
{
  "run_id": "...",
  "run_dir": "fresh/runs/.../",
  "mode": "test|production",
  "dry_run": true,
  "timestamp": "ISO-8601",
  "deleted_files": [
    {"path": "fresh/runs/.../scratch/foo.tmp", "size_bytes": 1234, "deleted": false}
  ],
  "preserved_files": [
    {"path": "fresh/runs/.../manifest/run_manifest.json", "reason": "manifest"}
  ],
  "deleted_count": 0,
  "preserved_count": 6,
  "errors": [],
  "status": "PASS|WARN|FAIL"
}
```

Notes:
- `deleted` boolean per file: `true` if actually deleted, `false` if dry-run or preservation
- `errors[]` entries record any individual file deletion failure (permission, in-use, etc.) without aborting the whole command — these elevate status to WARN

## 8. Severity rules

```text
PASS:  successful cleanup with zero errors
WARN:  cleanup completed but some files could not be deleted (recorded in errors[])
FAIL:  attempted to delete outside run_dir, or run_id rejected, or required preservation directory was missing
```

## 9. Tests required

Create `fresh/tests/test_m1_phase1_cleanup_manager.py` with at least these tests:

```text
test_cleanup_refuses_outside_run_dir
test_cleanup_preserves_manifest_logs_qc_reports
test_cleanup_dry_run_makes_no_changes
test_cleanup_test_mode_removes_scratch_and_tmp
test_cleanup_production_default_is_dry_run
test_cleanup_writes_cleanup_report_json
test_cleanup_appends_phase_status
test_cli_help_includes_cleanup
```

Test conventions (match `test_task9_pocket_candidate_prioritization.py` style):
- Local helpers: `unique_run_id()`, `make_tmp_run_context(tmp_path)`, `read_csv()`, `read_json()`, `assert_under_run_dir()`
- Direct API call for most tests
- One `subprocess.run([sys.executable, "-m", "egfr_myo1d.cli", "--help"])` for the help test
- pytest tmp_path fixture for isolation
- No requirement for real PDB / ligand files

## 10. Acceptance commands

Run all of these from repo root before reporting completion:

```bash
export PYTHONPATH="$PWD/fresh/src:${PYTHONPATH:-}"

# Full pytest including new tests
pytest -q fresh/tests

# Targeted phase test
pytest -q fresh/tests/test_m1_phase1_cleanup_manager.py

# CLI smoke
python -m egfr_myo1d.cli --help | grep cleanup
python -m egfr_myo1d.cli init-run --mode smoke_env --run-id m1_phase1_local
python -m egfr_myo1d.cli cleanup --run-id m1_phase1_local --mode test --dry-run true
python -m egfr_myo1d.cli cleanup --run-id m1_phase1_local --mode test --dry-run false
python -m egfr_myo1d.cli cleanup --run-id m1_phase1_local --mode production
python -m egfr_myo1d.cli status --run-id m1_phase1_local

# Path traversal negative test
python -m egfr_myo1d.cli cleanup --run-id ../bad_run --mode test

# Old workflow protection
git diff --name-only -- run_production.py main.py egfr_pipeline/ config/ docs/runbook.md output/ results_export/
```

Expected:
- Path traversal exits nonzero with no outside directory created
- Old workflow diff prints nothing
- All cleanup invocations write `manifest/cleanup_report.json`
- Production default invocation has `dry_run: true` in the report

## 11. Mandatory end-of-task self-test block

Codex / Claude must run the §10 acceptance commands and record actual results. No claims of "tests pass" without evidence.

Record per command:
- Command line
- Exit code
- Short result summary

For tests: number of tests collected, passed, failed, warnings.

## 12. Final response format

```text
M1 Phase 1 status: PASS / PASS WITH WARNINGS / FAIL

Files created:
- ...

Files modified:
- ...

Files deleted: (none in this phase)

Commands run and results:
- ...

Test summary:
- prior tests passing: 98 (or current count)
- new phase 1 tests: 8 passing
- total: 106

Cleanup report verification:
- manifest/cleanup_report.json schema validated
- preservation list intact
- dry-run vs real run outputs differ only in deleted boolean

Old workflow protection:
- git diff prints nothing for protected paths

Acceptance closure:
- M1 §23 #8 closed (cleanup safe + writes cleanup_report.json)
- HPC-pending items: none for this phase

Known limitations / not implemented by design:
- No PBS generation (Phase 6)
- No receptor / MYO1D / ligand handling (later phases)
- No M2 work
```

Do not end with vague claims. Either show command output or state precisely why a command could not be run.
