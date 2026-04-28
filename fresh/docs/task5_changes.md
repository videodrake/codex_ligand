# Task 5 Changes

Task 5 added the real-input readiness bridge for the fresh EGFR-MYO1D workflow.

## Created

- `fresh/src/egfr_myo1d/validation/real_inputs.py`
- `fresh/tests/test_task5_real_input_readiness.py`
- `fresh/docs/task5_real_input_readiness.md`
- `fresh/docs/task5_changes.md`

## Modified

- `fresh/src/egfr_myo1d/cli.py`

## New CLI

```bash
python -m egfr_myo1d.cli validate-real-inputs \
  --run-id <run_id> \
  --mode smoke_input \
  --profile hpc_strict \
  --input-root fresh/data/raw \
  --contract fresh/data/raw/ppi_input_contract.json
```

## Guardrails

The command writes only readiness manifests and QC CSVs under `fresh/runs/<run_id>/`. It does not run docking, mutate structures, repair EGFR mutations, renumber residues, normalize receptors, generate PBS/qsub files, delete cleanup targets, score compounds, or nominate candidates.

