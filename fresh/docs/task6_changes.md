# Task 6 Changes

Created:

- `fresh/src/egfr_myo1d/planning/__init__.py`
- `fresh/src/egfr_myo1d/planning/ppi_sampling.py`
- `fresh/src/egfr_myo1d/planning/pose_qc_policy.py`
- `fresh/src/egfr_myo1d/validation/ppi_sampling_plan.py`
- `fresh/tests/test_task6_ppi_sampling_plan.py`
- `fresh/docs/task6_ppi_sampling_plan.md`
- `fresh/docs/task6_changes.md`

Modified:

- `fresh/src/egfr_myo1d/cli.py`

Summary:

- Added the `plan-ppi-sampling` CLI command.
- Added spec-only PPI job planning for future PyRosetta and LightDock engines.
- Added future pose-acceptance policy serialization.
- Added strict gating through Task 5 real-input readiness.
- Added quarantine/blocker reporting for fixture-only artifacts and strict-mode
  blockers.
- Added tests for job specs, active-face propagation, score-bonus prohibition,
  strict blockers, path safety, old-workflow protection, and heavy-tool
  non-import behavior.

Tests:

- Run `pytest -q fresh/tests`.


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
