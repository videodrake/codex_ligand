# Task 3 Changes

Task 3 added structural input contracts and QC reporting for EGFR receptor inputs, MYO1D beta-meander partner inputs, and membrane-frame metadata.

## Created

- `fresh/src/egfr_myo1d/structure/__init__.py`
- `fresh/src/egfr_myo1d/structure/pdb_parser.py`
- `fresh/src/egfr_myo1d/structure/contracts.py`
- `fresh/src/egfr_myo1d/structure/geometry.py`
- `fresh/src/egfr_myo1d/structure/myo1d_annotation.py`
- `fresh/src/egfr_myo1d/validation/structure_inputs.py`
- `fresh/tests/test_task3_structure_input_contracts.py`
- `fresh/tests/fixtures/task3_inputs/contract.json`
- `fresh/tests/fixtures/task3_inputs/egfr_valid_dimer_AB.pdb`
- `fresh/tests/fixtures/task3_inputs/egfr_single_chain_dimer_ambiguous_X.pdb`
- `fresh/tests/fixtures/task3_inputs/egfr_residue_reset_bad.pdb`
- `fresh/tests/fixtures/task3_inputs/egfr_v924r_warn.pdb`
- `fresh/tests/fixtures/task3_inputs/myo1d_955_1001_valid.pdb`
- `fresh/tests/fixtures/task3_inputs/myo1d_962_1006_bad_terminal.pdb`
- `fresh/tests/fixtures/task3_inputs/myo1d_955_1006_tail_masked_warn.pdb`
- `fresh/tests/fixtures/task3_inputs/membrane_frame_valid.json`
- `fresh/tests/fixtures/task3_inputs/membrane_frame_zero_normal_bad.json`
- `fresh/docs/task3_structure_input_contracts.md`
- `fresh/docs/task3_changes.md`

## Modified

- `fresh/src/egfr_myo1d/cli.py`

## Guardrails

Task 3 does not implement docking, pocket discovery, receptor normalization, MYO1D production slicing or relaxation, qsub/PBS generation, cleanup deletion, scoring, or candidate nomination.

All generated QC outputs are restricted to `fresh/runs/<run_id>/`.

