# M2 Phase 2 Changes

Implements the M2.2 PyRosetta adapter harness as a dry-run-only stage.

## Files Created

```text
fresh/src/egfr_myo1d/ppi/__init__.py
fresh/src/egfr_myo1d/ppi/pyrosetta_adapter.py
fresh/src/egfr_myo1d/ppi/run_ppi_job.py
fresh/scripts/pbs/run_ppi_state_seed.pbs.template
fresh/tests/test_m2_phase2_pyrosetta_adapter.py
fresh/docs/m2_phase2_pyrosetta_adapter.md
fresh/docs/m2_phase2_changes.md
```

## Files Modified

```text
fresh/src/egfr_myo1d/cli.py
```

## CLI Additions

```bash
python -m egfr_myo1d.cli prepare-m2-pyrosetta-harness --run-id RUN
python -m egfr_myo1d.ppi.run_ppi_job --run-id RUN --job-name JOB --dry-run true
```

## Behavior

- Consumes M2.1 PyRosetta input specs.
- Writes dry-run job manifests under `fresh/runs/<run_id>/phase1_ppi/pyrosetta_adapter/`.
- Writes centralized job stdout/stderr paths under `fresh/runs/<run_id>/logs/jobs/`.
- Dry-run jobs create `dry_run_status.json` under the run directory and update job/phase logs.
- Records `execution_allowed=false`, `pyrosetta_imported=false`, `docking_executed=false`, and `relaxation_executed=false`.

## Not Implemented By Design

- PyRosetta import, docking, or relaxation
- LightDock execution
- qsub/PBS/sbatch submission
- PPI pose collection/contact extraction
- MYO1D pose QC or artifact filtering
- Pocket discovery, compound docking, scoring, or candidate nomination

