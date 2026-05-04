# M2.2 PyRosetta Adapter Harness

M2.2 prepares a PyRosetta adapter harness from M2.1 input specs. In this implementation, the harness is dry-run only: it validates wiring, writes job manifests, and records commands, but it does not import PyRosetta or run docking/relaxation.

## Inputs

```text
fresh/runs/<run_id>/manifest/m2_1_ppi_input_manifest.json
fresh/runs/<run_id>/prepared/m2_1_ppi_inputs/<state>/specs/*pyrosetta-global-ppi_input_spec.json
```

The M2.1 spec carries the receptor pack, MYO1D primary partner pack, receptor mapping CSV, and membrane context.

## Outputs

```text
fresh/runs/<run_id>/phase1_ppi/pyrosetta_adapter/pyrosetta_job_manifest.csv
fresh/runs/<run_id>/phase1_ppi/pyrosetta_adapter/pyrosetta_job_manifest.jsonl
fresh/runs/<run_id>/phase1_ppi/pyrosetta_adapter/launch_pyrosetta_dry_runs.sh
fresh/runs/<run_id>/manifest/m2_2_pyrosetta_adapter_manifest.json
fresh/runs/<run_id>/qc/m2_2_pyrosetta_adapter_qc.csv
fresh/runs/<run_id>/reports/m2_2_pyrosetta_adapter.md
```

Dry-run job execution writes:

```text
fresh/runs/<run_id>/phase1_ppi/pyrosetta_adapter/outputs/<job_name>/dry_run_status.json
```

## CLI

```bash
python -m egfr_myo1d.cli prepare-m2-pyrosetta-harness \
  --run-id <run_id> \
  --mode smoke_input \
  --profile codex_dev

python -m egfr_myo1d.ppi.run_ppi_job \
  --run-id <run_id> \
  --job-name <job_name> \
  --dry-run true
```

Smoke mode emits one state/seed dry-run job. Mini emits primary-state jobs for seeds 0 and 1. Production emits manifests only and still does not submit or execute anything.

## Guardrails

- PyRosetta is not imported in M2.2.
- Docking and relaxation are not executed.
- PyRosetta output numbering is not trusted or consumed in this phase.
- The M1 receptor mapping CSV remains the source of truth.
- stdout/stderr paths are under `fresh/runs/<run_id>/logs/jobs/`.
- All raw or dry-run outputs remain under `fresh/runs/<run_id>/`.
- No qsub/PBS submission is performed.
- No scoring, pocket discovery, compound docking, or candidate nomination is performed.

`fresh/scripts/pbs/run_ppi_state_seed.pbs.template` is a dry-run template for later HPC wrapping. It does not call `qsub`.

