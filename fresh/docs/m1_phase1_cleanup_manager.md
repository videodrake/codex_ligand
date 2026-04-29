# M1 Phase 1 — Cleanup Manager

Implements the cleanup policy from `milestone1_foundation_codex_handoff_v0_5.md` §18, closing M1 §23 #8.

## What it does

`fresh/src/egfr_myo1d/core/cleanup.py` walks the active run directory and:

1. Identifies deletion candidates per handoff §18.2:
   - Every file under `scratch/` and `tmp/` (recursively)
   - Files anywhere under `run_dir` whose name matches `*.tmp`, `*.pdbqt.tmp`, `*.vina.tmp`, or `*.silent` (excluding the preserved subdirectories below)
2. Identifies preserved files per handoff §18.3:
   - Every file under `manifest/`, `logs/`, `qc/`, and `reports/`
3. Either reports (dry-run) or deletes the candidate files (real mode)
4. Always writes `manifest/cleanup_report.json` with the deleted/preserved file lists
5. Appends a phase status entry to `logs/phase_status.jsonl` and `logs/master.log`
6. Refuses to touch anything outside `ctx.run_dir` — every candidate is double-checked via `_is_under` before deletion

## Modes

| `--mode` | Default `--dry-run` | Notes |
| --- | --- | --- |
| `test` | `false` | Used between development runs to clean scratch artifacts |
| `production` | `true` | Defensive default; explicit `--dry-run false` required to actually delete |

## CLI

```bash
python -m egfr_myo1d.cli cleanup --run-id RUN --mode test
python -m egfr_myo1d.cli cleanup --run-id RUN --mode test --dry-run true
python -m egfr_myo1d.cli cleanup --run-id RUN --mode production
python -m egfr_myo1d.cli cleanup --run-id RUN --mode production --dry-run false
```

`fresh/scripts/cleanup_run.py` is a thin wrapper that sets `PYTHONPATH` and invokes the CLI.

## Output

```text
fresh/runs/<run_id>/manifest/cleanup_report.json
fresh/runs/<run_id>/logs/phase_status.jsonl     (appended)
fresh/runs/<run_id>/logs/master.log             (appended)
```

`cleanup_report.json` schema:

```json
{
  "run_id": "...",
  "run_dir": "...",
  "mode": "test|production",
  "dry_run": true,
  "profile": "codex_dev|hpc_strict",
  "timestamp": "ISO-8601",
  "deleted_files": [
    {"path": "scratch/foo.tmp", "size_bytes": 123, "deleted": false}
  ],
  "preserved_files": [
    {"path": "manifest/run_manifest.json", "reason": "preservation_list"}
  ],
  "deleted_count": 7,
  "preserved_count": 6,
  "errors": [],
  "status": "PASS|WARN|FAIL"
}
```

## Severity

| Status | Conditions |
| --- | --- |
| `PASS` | All eligible candidates handled with zero errors |
| `WARN` | Cleanup completed but some files could not be deleted (recorded in `errors[]`) |
| `FAIL` | Path-traversal `run_id`, missing `run_dir`, or attempt to delete outside `run_dir` |

## Safety guarantees

- Cleanup operates only inside `ctx.run_dir`. `RunContext.for_existing` validates the run id; the candidate enumerator only walks paths under `run_dir`; and an additional `_is_under` check before each `unlink()` raises `RunContextError` if anything escapes.
- `manifest/cleanup_report.json` is always preservation-listed; the cleanup never deletes its own report.
- Per-file deletion errors (permission, in-use, etc.) elevate status to `WARN`; they do not abort the run.

## What is intentionally not in this phase

- PBS/qsub generation (Phase 6)
- Receptor / MYO1D / membrane frame / ligand work (Phases 3-5, 7)
- M2 docking, scoring, candidate work
