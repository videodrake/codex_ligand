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

