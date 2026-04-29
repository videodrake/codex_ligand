# M1 Phase 1 — Changes

Closes M1 §23 #8 (cleanup safe + writes cleanup_report.json) per `milestone1_foundation_codex_handoff_v0_5.md` §18 and v1.0 plan §16 M1 Task 4.

## Files created

```text
fresh/src/egfr_myo1d/core/cleanup.py
fresh/tests/test_m1_phase1_cleanup_manager.py
fresh/docs/m1_phase1_cleanup_manager.md
fresh/docs/m1_phase1_changes.md
```

## Files modified

```text
fresh/src/egfr_myo1d/cli.py        # added `cleanup` subparser + _cmd_cleanup handler
fresh/scripts/cleanup_run.py       # replaced placeholder with real CLI wrapper
```

## Files deleted

None.

## Public API additions

```python
from egfr_myo1d.core.cleanup import (
    run_cleanup,            # primary entry point
    CleanupReport,          # dataclass
    resolve_dry_run_default,
    DELETION_FILE_PATTERNS,
    DELETION_DIR_BASENAMES,
    PRESERVED_DIR_BASENAMES,
)
```

## CLI surface additions

```bash
python -m egfr_myo1d.cli cleanup --run-id RUN --mode test [--dry-run true|false] [--profile codex_dev|hpc_strict]
python -m egfr_myo1d.cli cleanup --run-id RUN --mode production [--dry-run true|false] [--profile codex_dev|hpc_strict]
```

Total CLI subparser count after Phase 1: 12 (was 11).

## Acceptance closure

- M1 §23 #8 closed: cleanup is safe, refuses to delete outside `run_dir`, supports dry-run, preserves `manifest/`, `logs/`, `qc/`, `reports/`, writes `manifest/cleanup_report.json`, defaults to dry-run for `--mode production`.

## Out of scope (next phases)

- Phase 2: relocate `preparation/{constructs,pdb_writer}.py` → `myo1d/`
- Phase 3: MYO1D construct + QC (`prepare-myo1d`)
- Phase 4: receptor normalization (`prepare-receptor`)
- Phase 5: membrane frame generation (`compute-membrane-frame`)
- Phase 6: PBS generator (`prepare-pbs`)
- Phase 7: ligand manifest shell (`manifest-ligands`)
- Phase 8: `prepare-inputs` orchestrator + M1 integration test
- Phase 9: Tasks 4-9 schema realignment
