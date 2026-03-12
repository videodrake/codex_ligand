# Config Guide

## Status Note

Current documentation baseline:

- Phase 1 primary engine: PyRosetta
- Phase 1 active secondary validation: LightDock
- AFM parser code still exists, but AFM is not part of the active routine workflow

Read [current_pipeline_status.md](/Users/admin/Desktop/hwang/codex/codex_ligand/docs/current_pipeline_status.md) before using older planning notes.

This directory contains runtime configuration and PBS submission files for the
current EGFR-MYO1D workflow.

## What Lives Here

| File | Type | Purpose |
|------|------|---------|
| `example-project.yaml` | YAML | Project-level Vina config example |
| `full_test.yaml` | YAML | Full-pipeline test config |
| `ppi_test_beta_meander.ini` | INI | PyRosetta PPI test config for beta-meander |
| `ppi_test_TH1.ini` | INI | PyRosetta PPI test config for TH1 |
| `ppi_prod_beta_meander.ini` | INI | PyRosetta PPI production config for beta-meander |
| `ppi_prod_TH1.ini` | INI | PyRosetta PPI production config for TH1 |
| `run_pre_qsub_checks.pbs` | PBS | Scheduler-side validation lane before production |
| `run_production.pbs` | PBS | Main production pipeline submission |
| `run_ppi_test.pbs` | PBS | PyRosetta PPI test submission |
| `run_ppi_prod.pbs` | PBS | PyRosetta PPI production submission |
| `run_full_test.pbs` | PBS | Full test submission helper |

## Current Separation of Config Styles

Two config styles are in use right now.

- Vina and the unified CLI use YAML project configs.
- Legacy PyRosetta PPI entrypoints still use INI configs.

This split is intentional for the current repository state. Do not assume that
all runtime paths share one config schema yet.

## Vina / Unified CLI Usage

Use YAML configs for the current project-level flow:

```bash
python main.py vina --config config/example-project.yaml
python main.py postprocess --config config/example-project.yaml
python main.py validate
```

The YAML model is where current project-wide settings belong:
- receptor definitions
- ligand definitions
- Vina parameters
- postprocess parameters
- worker count
- output root

## PyRosetta PPI Usage

PyRosetta PPI still uses INI configs and PBS wrappers.

Interactive or direct run:

```bash
python main.py pyrosetta
python -m egfr_pipeline.pyrosetta_docking.pipeline_manager config/ppi_test_TH1.ini
```

PBS test submission:

```bash
qsub config/run_ppi_test.pbs
qsub -v CONFIG_FILE=config/ppi_test_TH1.ini config/run_ppi_test.pbs
qsub -v RUN_MODE=both config/run_ppi_test.pbs
```

PBS production submission:

```bash
qsub config/run_ppi_prod.pbs
qsub -v CONFIG_FILE=config/ppi_prod_TH1.ini config/run_ppi_prod.pbs
qsub -v RUN_MODE=both config/run_ppi_prod.pbs
```

## Pre-Qsub Validation Lane

Before heavy scheduler runs, use the precheck lane:

```bash
qsub config/run_pre_qsub_checks.pbs
```

On success, it writes:

```bash
output/pre_qsub_status/last_pass.json
```

`run_production.pbs` now checks for that success marker by default. The safest
submission pattern is:

```bash
PRECHECK_JOB=$(qsub config/run_pre_qsub_checks.pbs)
qsub -W depend=afterok:${PRECHECK_JOB} config/run_production.pbs
```

If you intentionally need to bypass the guard:

```bash
qsub -v SKIP_PRECHECK_GUARD=1 config/run_production.pbs
```

The current server-side baseline is the shared `pyrosetta` conda environment.
Repository scripts do not create a separate test environment anymore.

- environment name: `pyrosetta`
- package reference list: [requirements-test.txt](/Users/admin/Desktop/hwang/codex/codex_ligand/requirements-test.txt)
- manual setup note: [server_environment_setup.md](/Users/admin/Desktop/hwang/codex/codex_ligand/docs/server_environment_setup.md)

## Production Modes

`run_production.pbs` forwards mode selection into `run_production.py`.

Examples:

```bash
qsub config/run_production.pbs
qsub -v MODE=force config/run_production.pbs
qsub -v MODE=from,FROM=4 config/run_production.pbs
qsub -v MODE=status config/run_production.pbs
qsub -v MODE=vina-only config/run_production.pbs
qsub -v MODE=ppi-only config/run_production.pbs
qsub -v MODE=post-only config/run_production.pbs
```

## Notes on the INI Files

The current `ppi_*.ini` files now carry explicit metadata fields used by the
Phase 1 traceability path, including:
- `receptor_id`
- `partner_id`
- `construct_type`
- `receptor_construct`
- `partner_construct`
- chain IDs
- numbering system

They also use metadata-tagged output directory naming so test and production
runs do not overwrite each other as easily.

## Recommended Reading

For operational details, also see:
- [pre_qsub_test_line.md](/Users/admin/Desktop/hwang/codex/codex_ligand/docs/pre_qsub_test_line.md)
- [phase1_output_chain_note.md](/Users/admin/Desktop/hwang/codex/codex_ligand/docs/phase1_output_chain_note.md)
- [phase1_pyrosetta_execution_note.md](/Users/admin/Desktop/hwang/codex/codex_ligand/docs/phase1_pyrosetta_execution_note.md)
