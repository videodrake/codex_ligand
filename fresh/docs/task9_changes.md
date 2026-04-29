# Task 9 Changes

Created:

- `fresh/src/egfr_myo1d/analysis/pocket_candidate_prioritization.py`
- `fresh/src/egfr_myo1d/validation/pocket_candidate_prioritization.py`
- `fresh/tests/fixtures/task9_pocket_candidates/ppi_consensus_patch.csv`
- `fresh/tests/fixtures/task9_pocket_candidates/pocket_discovery_plan.json`
- `fresh/tests/fixtures/task9_pocket_candidates/pocket_detector_candidates.csv`
- `fresh/tests/test_task9_pocket_candidate_prioritization.py`
- `fresh/docs/task9_pocket_candidate_prioritization.md`
- `fresh/docs/task9_changes.md`

Modified:

- `fresh/src/egfr_myo1d/cli.py`
- `fresh/src/egfr_myo1d/analysis/pocket_selection.py`
- `fresh/src/egfr_myo1d/validation/pocket_discovery.py`
- `fresh/tests/test_task8_pocket_discovery_plan.py`
- `fresh/docs/task8_pocket_discovery_plan.md`
- `fresh/docs/task8_changes.md`

Summary:

- Added `prioritize-pocket-candidates`.
- Added detector-agnostic pocket candidate schema validation.
- Added deterministic PPI proximity, ATP exclusion, membrane compatibility, and dimer accessibility classification.
- Added prioritized pocket candidate and family outputs.
- Added explicit no-runtime/no-docking/no-scoring/no-nomination manifest flags.
- Preserved protomer identity and EGFR residue numbering from supplied records.

Validation:

- `python -m py_compile fresh/src/egfr_myo1d/analysis/pocket_candidate_prioritization.py fresh/src/egfr_myo1d/validation/pocket_candidate_prioritization.py fresh/tests/test_task9_pocket_candidate_prioritization.py`
- `pytest -q fresh/tests/test_task8_pocket_discovery_plan.py fresh/tests/test_task9_pocket_candidate_prioritization.py`
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
