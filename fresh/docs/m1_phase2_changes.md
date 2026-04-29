# M1 Phase 2 — Changes

Pure structural move; no logic change. Aligns module tree with handoff §4 spec; unblocks Phase 3.

## Files created

```text
fresh/src/egfr_myo1d/myo1d/__init__.py
fresh/docs/m1_phase2_myo1d_relocation.md
fresh/docs/m1_phase2_changes.md
```

## Files moved (git mv)

```text
fresh/src/egfr_myo1d/preparation/constructs.py  →  fresh/src/egfr_myo1d/myo1d/construct.py
fresh/src/egfr_myo1d/preparation/pdb_writer.py  →  fresh/src/egfr_myo1d/myo1d/pdb_writer.py
```

## Files modified (imports only)

```text
fresh/src/egfr_myo1d/myo1d/construct.py            # 1 self-import: preparation.pdb_writer → myo1d.pdb_writer
fresh/src/egfr_myo1d/validation/prepared_inputs.py # 2 import statements consolidated and rewritten to myo1d.*
fresh/src/egfr_myo1d/validation/real_inputs.py     # 1 import statement rewritten to myo1d.construct
fresh/tests/test_task4_ppi_input_preparation.py    # 1 test import rewritten to myo1d.construct
```

## Files deleted

None directly — the two relocated files are tracked by git as renames.

## Logic changes

None.

## Verification

- All 98 prior tests + 16 Phase 1 tests = 114 tests pass.
- `python -c "from egfr_myo1d.myo1d import construct, pdb_writer"` works.
- `grep -rn "from egfr_myo1d.preparation.constructs\|from egfr_myo1d.preparation.pdb_writer" fresh/src fresh/tests --include="*.py"` returns no matches.
- Old workflow files unmodified: `git diff --name-only -- run_production.py main.py egfr_pipeline/ config/ docs/runbook.md output/ results_export/` returns empty.

## Acceptance closure

Phase 2 closes nothing new in M1 §23 by itself — it is a prerequisite restructure for Phase 3, which closes M1 §23 #13 (MYO1D 955-1006 construct QC).

## Out of scope (next phase)

Phase 3:
- Extend `myo1d/construct.py` with `slice_myo1d_construct`, `emit_myo1d_construct_pdb`
- Create `myo1d/qc.py` with `run_myo1d_qc`
- Add `prepare-myo1d` CLI subcommand
- New tests under `fresh/tests/test_m1_phase3_myo1d_construct_qc.py`
