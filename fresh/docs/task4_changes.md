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

