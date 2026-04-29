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
