# Task 8 Changes

Created:

- `fresh/src/egfr_myo1d/analysis/pocket_selection.py`
- `fresh/src/egfr_myo1d/validation/pocket_discovery.py`
- `fresh/tests/fixtures/task8_pocket_planning/ppi_consensus_patch.csv`
- `fresh/tests/test_task8_pocket_discovery_plan.py`
- `fresh/docs/task8_pocket_discovery_plan.md`
- `fresh/docs/task8_changes.md`

Modified:

- `fresh/src/egfr_myo1d/cli.py`

Summary:

- Added `plan-pocket-discovery`.
- Added Task 7 consensus patch schema intake.
- Added deterministic pocket-selection planning records.
- Added ATP, membrane, dimer/protomer, and PPI-evidence audits.
- Added `NO_GO` behavior when no accepted PPI patch is available.
- Preserved zero docking/scoring/candidate-nomination execution.

Validation:

- `python -m py_compile fresh/src/egfr_myo1d/analysis/pocket_selection.py fresh/src/egfr_myo1d/validation/pocket_discovery.py fresh/tests/test_task8_pocket_discovery_plan.py`
- `pytest -q fresh/tests/test_task8_pocket_discovery_plan.py`
- `pytest -q fresh/tests`

