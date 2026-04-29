# Task 7 Changes

Created:

- `fresh/src/egfr_myo1d/analysis/__init__.py`
- `fresh/src/egfr_myo1d/analysis/ppi_consensus.py`
- `fresh/src/egfr_myo1d/validation/ppi_consensus.py`
- `fresh/tests/fixtures/task7_ppi_consensus/accepted_ppi_contacts.csv`
- `fresh/tests/test_task7_ppi_consensus_patch.py`
- `fresh/docs/task7_ppi_consensus_patch.md`
- `fresh/docs/task7_changes.md`

Modified:

- `fresh/src/egfr_myo1d/cli.py`

Summary:

- Added `summarize-ppi-consensus`.
- Added schema validation for supplied PPI contact tables.
- Added EGFR and MYO1D residue-list parsing audits.
- Added active-face, sheet-12, tail/noise, ATP-overlap, membrane-proximal, and
  convergence audits.
- Added deterministic EGFR-side consensus patch CSV generation.
- Preserved Task 6 spec-only posture and no-score-bonus policy.

Validation:

- `python -m py_compile fresh/src/egfr_myo1d/analysis/ppi_consensus.py fresh/src/egfr_myo1d/validation/ppi_consensus.py fresh/tests/test_task7_ppi_consensus_patch.py`
- `pytest -q fresh/tests/test_task7_ppi_consensus_patch.py`
- `pytest -q fresh/tests`


---

## M1 dependency (added by Phase 9 alignment)

This task module was originally designed to run before the M1 foundation
modules (cleanup, receptor normalize, membrane frame, MYO1D construct, ligand
manifest) existed. After M1 Phases 1-8 the canonical input artifacts now live
under `fresh/runs/<run_id>/normalized/`, `fresh/runs/<run_id>/manifest/membrane_frame.json`,
and `fresh/runs/<run_id>/qc/<state>_receptor_mapping.csv`.

Phase 9 takes an **additive** approach: this module's logic and output paths
are unchanged. The alignment between M1 outputs and this task's consumption
points is recorded by `validation/m1_alignment.py`'s `record_m1_alignment(ctx)`
helper, which writes `manifest/m1_alignment.json` listing which M1 artifacts
are present and which Task 4-9 modules would naturally consume them in a future
M2 actual-execution phase.

See `fresh/docs/m1_phase9_tasks4to9_realignment.md` for the full rationale.
