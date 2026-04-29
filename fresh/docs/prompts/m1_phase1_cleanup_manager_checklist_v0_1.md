# M1 Phase 1 Acceptance Checklist v0.1 — Cleanup Manager

Use this after the implementer applies M1 Phase 1.

## 1. Confirm pre-Phase state preserved

Old workflow files/directories must not be modified:

```text
run_production.py
main.py
egfr_pipeline/**
config/**
docs/runbook.md
output/**
results_export/**
```

Tasks 1-9 modules under `fresh/src/egfr_myo1d/` must not be modified except `cli.py` (subparser addition only).

## 2. Import and CLI smoke

Bash:

```bash
export PYTHONPATH="$PWD/fresh/src:${PYTHONPATH:-}"
python -m egfr_myo1d.cli --help
python -m egfr_myo1d.cli version
```

PowerShell:

```powershell
$env:PYTHONPATH = "$PWD\fresh\src;$env:PYTHONPATH"
python -m egfr_myo1d.cli --help
python -m egfr_myo1d.cli version
```

Expected:

```text
--help output includes the new `cleanup` subcommand line
version exits 0
no docking starts, no old workflow starts
```

## 3. Cleanup test mode (dry-run)

```bash
python -m egfr_myo1d.cli init-run --mode smoke_env --run-id m1_phase1_local
# Seed scratch/tmp files for the cleanup target
mkdir -p fresh/runs/m1_phase1_local/scratch
echo dummy > fresh/runs/m1_phase1_local/scratch/foo.tmp
echo dummy > fresh/runs/m1_phase1_local/tmp/bar.pdbqt.tmp

python -m egfr_myo1d.cli cleanup --run-id m1_phase1_local --mode test --dry-run true
```

Expected:

```text
Exit 0
fresh/runs/m1_phase1_local/manifest/cleanup_report.json exists
report.dry_run == true
report.deleted_count >= 2
foo.tmp and bar.pdbqt.tmp still exist on disk (dry-run did not delete)
```

## 4. Cleanup test mode (real delete)

```bash
python -m egfr_myo1d.cli cleanup --run-id m1_phase1_local --mode test --dry-run false
```

Expected:

```text
Exit 0
report.dry_run == false
report.deleted[].deleted == true for scratch/foo.tmp and tmp/bar.pdbqt.tmp
foo.tmp and bar.pdbqt.tmp no longer exist on disk
manifest/, logs/, qc/, reports/ directories still exist and are populated
```

## 5. Cleanup production mode (default dry-run)

```bash
python -m egfr_myo1d.cli cleanup --run-id m1_phase1_local --mode production
```

Expected:

```text
Exit 0
report.dry_run == true (default for production mode)
no files deleted regardless of scratch/tmp content
```

## 6. Path traversal safety

```bash
python -m egfr_myo1d.cli cleanup --run-id ../bad_run --mode test
python -m egfr_myo1d.cli cleanup --run-id /absolute/bad --mode test
```

Expected:

```text
Both exit nonzero
No directory created or modified outside fresh/runs/
No cleanup_report.json written outside fresh/runs/
```

## 7. Status command after cleanup

```bash
python -m egfr_myo1d.cli status --run-id m1_phase1_local
```

Expected stdout:

```text
run_id
run_dir
last phase status (one entry should be cleanup PASS or WARN)
WARN/FAIL counts
master.log path
cleanup_report.json path
```

## 8. Tests

```bash
pytest -q fresh/tests/test_m1_phase1_cleanup_manager.py
pytest -q fresh/tests
```

Required tests must pass (see prompt §9).

## 9. Old workflow protection

```bash
git diff --name-only -- run_production.py main.py egfr_pipeline/ config/ docs/runbook.md output/ results_export/
```

Expected:

```text
No output (empty diff)
```

## 10. Run-output containment

```bash
find . -path './fresh/runs/m1_phase1_local' -prune -o -name 'cleanup_report.json' -print
```

Expected:

```text
No matches outside fresh/runs/m1_phase1_local/
```

## 11. What must not be implemented yet (in this phase)

```text
PBS generation
receptor normalization
MYO1D construct slicing
membrane frame generation
ligand manifest
prepare-inputs orchestrator
M2 work (PyRosetta / Vina / fpocket / pose QC / consensus / pocket / candidate)
```

## 12. Phase 1 accepted if

```text
- core/cleanup.py exists with the public API in §6 of the prompt.
- cli.py exposes `cleanup` subcommand.
- scripts/cleanup_run.py is no longer a placeholder.
- Test mode dry-run leaves files; real delete removes scratch/*.tmp.
- Production default is dry_run=true.
- cleanup_report.json schema matches prompt §7.
- preservation list (manifest/, logs/, qc/, reports/) is never deleted.
- Path-traversal run_id rejected with nonzero exit and no outside writes.
- All 8+ phase 1 tests pass.
- Existing 98 tests pass.
- Old workflow files unmodified.
- M1 §23 #8 acceptance item is closed.
```

## 13. Implementer final response must include

```text
M1 Phase 1 status: PASS / PASS WITH WARNINGS / FAIL
Files created
Files modified
Commands run and results (with exit codes and short summaries)
Test counts (prior + new + total, all passing)
Cleanup report verification
Old workflow protection
Acceptance closure note: M1 §23 #8 closed
Known limitations / not implemented by design
```
