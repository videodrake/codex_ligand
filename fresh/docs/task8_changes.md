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
- Task 8.1 added `qc/task7_consensus_schema_audit.csv`.
- Task 8.1 added `pockets/egfr_myo1d_ppi_guided_pocket_plan_records.csv` as a clearer alias for the legacy planning CSV.
- Task 8.1 added explicit plan-vs-detected semantic flags to planning rows, plan JSON, and manifests.
- Added ATP, membrane, dimer/protomer, and PPI-evidence audits.
- Added `NO_GO` behavior when no accepted PPI patch is available.
- Preserved zero docking/scoring/candidate-nomination execution.

Validation:

- `python -m py_compile fresh/src/egfr_myo1d/analysis/pocket_selection.py fresh/src/egfr_myo1d/validation/pocket_discovery.py fresh/tests/test_task8_pocket_discovery_plan.py`
- `pytest -q fresh/tests/test_task8_pocket_discovery_plan.py`
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
