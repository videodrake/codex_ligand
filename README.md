# EGFR-MYO1D Pipeline (`codex_ligand`)

This workspace contains the active EGFR-MYO1D state-comparison pipeline. The routine baseline is Vina-centered ligand analysis with Phase 1 PPI evidence from PyRosetta and secondary validation support from LightDock.

## What This README Covers

- what this workspace is
- where to start reading
- how to run main commands
- where outputs are written

For deep scientific or architecture detail, use `docs/`.

## Quick Onboarding

Read in this order:

1. [docs/onboarding/README.md](docs/onboarding/README.md)
2. [docs/AI_START_HERE.md](docs/AI_START_HERE.md)
3. [docs/current_pipeline_status.md](docs/current_pipeline_status.md)
4. [docs/first_time_environment_setup.md](docs/first_time_environment_setup.md)
5. [docs/runbook.md](docs/runbook.md)
6. [docs/manual_execution.md](docs/manual_execution.md)
7. [config/README.md](config/README.md)
8. [output/README.md](output/README.md)

## Current Output Reading Order

For production runs driven by `run_production.py` or `qsub config/run_production.pbs`, start interpretation at `output/{project}/step_index.md`.

Recommended reading order:

1. `output/{project}/step_index.md`
2. `output/{project}/step6_report/project_report.txt`
3. `output/{project}/step5_verdict/valid_sites.csv`
4. `output/{project}/step4_vina_postprocess/vina_pocket_table.csv`
5. `output/{project}/step3_ppi_postprocess/ppi_pyrosetta_residues.csv`

Canonical runtime outputs remain under the existing project root. The `step1_vina_raw/` through `step7_validate/` folders are derived interpretation views that can be regenerated from canonical outputs; they do not replace the root artifacts.

## Repository Layout

- `main.py`: Unified CLI entry point.
- `egfr_pipeline/`: Core implementation package.
- `config/`: YAML, INI, and PBS wrappers.
- `docs/`: Onboarding, runbooks, architecture, and phase plans.
- `input/`: Receptor and ligand inputs.
- `output/`: Baseline and phase-separated outputs.
- `tests/`: Validation and test suite.
- `scripts/`: Utility scripts used by workflows.

## Command Quickstart

Run from `codex_ligand/` after activating the expected environment.

Prerequisites:

- `conda activate pyrosetta`
- main config: `config/example-project.yaml`
- run precheck before production or heavy submissions

```bash
python main.py --help
python main.py -c config/example-project.yaml validate --help
qsub config/run_pre_qsub_checks.pbs
```

Routine baseline lane (execution order):

```bash
python main.py -c config/example-project.yaml vina
python main.py -c config/example-project.yaml postprocess
python main.py -c config/example-project.yaml verdict
python main.py -c config/example-project.yaml report
python main.py -c config/example-project.yaml validate
```

Additional commands:

```bash
python main.py -c config/example-project.yaml pyrosetta
python main.py -c config/example-project.yaml md
python main.py -c config/example-project.yaml ppi-postprocess
python main.py -c config/example-project.yaml full
```

`md` opens the MD analysis submenu; the downstream analysis tools still take their own CLI arguments after that entry point.

Scheduler wrappers:

```bash
qsub config/run_pre_qsub_checks.pbs
qsub config/run_production.pbs
```

Expected output checkpoints after the routine baseline lane:

- `output/egfr_myo1d_vina/step_index.md`
- `output/egfr_myo1d_vina/vina_pose_table.csv`
- `output/egfr_myo1d_vina/vina_pocket_table.csv`
- `output/egfr_myo1d_vina/valid_sites.csv`
- `output/egfr_myo1d_vina/project_report.txt`

## Output Entry Points

- `output/{project}/step_index.md`: First human-readable entry point for completed production runs.
- [output/README.md](output/README.md): Output root index.
- [output/phase1_ppi/README.md](output/phase1_ppi/README.md)
- [output/phase2_pockets/README.md](output/phase2_pockets/README.md)
- [output/phase3_docking/README.md](output/phase3_docking/README.md)
- [output/phase4_perturbation/README.md](output/phase4_perturbation/README.md)

Routine baseline project output root:

- `output/egfr_myo1d_vina/`

## Documentation Indexes

- [docs/README.md](docs/README.md): Full docs index.
- [docs/onboarding/README.md](docs/onboarding/README.md): New-contributor package.
- [config/README.md](config/README.md): Config semantics and wrapper roles.

## Current Baseline Guardrails

- Treat AFM as legacy optional unless explicitly re-enabled.
- Keep the three receptor states separated in interpretation and reporting.
- Treat `max_workers = 16` as the safe routine operating bound.
- Prefer active code/config over older prose when conflicts appear.

