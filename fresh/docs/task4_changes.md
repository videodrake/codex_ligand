# Task 4 Changes

Task 4 added pure-Python PPI input preparation and QC reports for future EGFR-MYO1D docking tasks.

## Created

- `fresh/src/egfr_myo1d/preparation/__init__.py`
- `fresh/src/egfr_myo1d/preparation/pdb_writer.py`
- `fresh/src/egfr_myo1d/preparation/constructs.py`
- `fresh/src/egfr_myo1d/preparation/restraints.py`
- `fresh/src/egfr_myo1d/preparation/masks.py`
- `fresh/src/egfr_myo1d/validation/prepared_inputs.py`
- `fresh/tests/test_task4_ppi_input_preparation.py`
- `fresh/tests/fixtures/task4_inputs/ppi_input_contract.json`
- `fresh/tests/fixtures/task4_inputs/egfr_dimer_AB_valid.pdb`
- `fresh/tests/fixtures/task4_inputs/egfr_single_chain_X_ambiguous.pdb`
- `fresh/tests/fixtures/task4_inputs/egfr_v924r_warn.pdb`
- `fresh/tests/fixtures/task4_inputs/myo1d_955_1001_source.pdb`
- `fresh/tests/fixtures/task4_inputs/myo1d_955_1006_tail_masked_source.pdb`
- `fresh/tests/fixtures/task4_inputs/myo1d_962_1006_terminal_artifact_bad.pdb`
- `fresh/tests/fixtures/task4_inputs/myo1d_capped_hetatm_955_1001.pdb`
- `fresh/tests/fixtures/task4_inputs/membrane_frame_valid.json`
- `fresh/docs/task4_ppi_input_preparation.md`
- `fresh/docs/task4_changes.md`

## Modified

- `fresh/src/egfr_myo1d/cli.py`

## Guardrails

The implementation preserves chain IDs, residue numbering, insertion codes, record types, and source mappings. It does not mutate EGFR residues, normalize ambiguous receptors into production inputs, drop caps/HETATM biological residues, run docking, generate PBS/qsub jobs, delete cleanup targets, score candidates, or nominate compounds.


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
