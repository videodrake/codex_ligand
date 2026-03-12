# Runbook

## EGFR-MYO1D Pipeline

This runbook is the operator-facing execution guide for the current repository.
It describes how the project is expected to be run and interpreted today.

## 1. Current Operating Baseline

Use this runbook together with:

1. `docs/current_pipeline_status.md`
2. `README.md`
3. `docs/architecture.md`
4. `config/README.md`

Current active baseline:

- Vina-centered ligand workflow
- PyRosetta-centered Phase 1 PPI workflow
- LightDock as the active secondary Phase 1 validation axis
- MD as a downstream stability gate

AlphaFold-Multimer is not part of the active routine workflow.

## 2. Fixed Inputs

Current receptor states:

- `3GT8_raw`
- `3GT8_cl38_48`
- `3GT8_cl85_100`

Current ligand set is defined in `config/example-project.yaml`.

PPI input structures live under `input/PPI/`.

## 3. Worker Policy

The server may expose 32 CPU cores, but routine safe operation should assume 16
workers unless there is explicit reason to go higher.

Treat 16 as the normal operating upper bound.

## 4. Normal Execution Paths

### Local or interactive orchestration

- `python main.py`
- `python main.py vina -c config/example-project.yaml`
- `python main.py postprocess -c config/example-project.yaml`
- `python main.py verdict -c config/example-project.yaml`
- `python main.py report -c config/example-project.yaml`
- `python main.py validate -c config/example-project.yaml`

### PyRosetta / server-side

- `qsub config/run_ppi_test.pbs`
- `qsub config/run_ppi_prod.pbs`

### Production orchestration

- `qsub config/run_pre_qsub_checks.pbs`
- `qsub config/run_production.pbs`

Safest submission pattern:

```bash
PRECHECK_JOB=$(qsub config/run_pre_qsub_checks.pbs)
qsub -W depend=afterok:${PRECHECK_JOB} config/run_production.pbs
```

## 5. Execution Order

### A. Pre-qsub validation

Run the lightweight validation lane before heavy server work.

Outputs:

- `output/pre_qsub_status/last_pass.json`

### B. Vina docking

Run docking for the three receptor states against the configured ligand set.

Main output family:

- raw docking pose files under the project output root

### C. Vina postprocess

Run:

- pose parsing
- contact extraction
- pocket clustering
- pocket summarization
- cross-receptor comparison
- optional bootstrap

Main outputs:

- `vina_pose_table.csv`
- `vina_pocket_table.csv`
- `vina_drug_pocket_map.csv`
- `vina_pocket_comparison.csv`
- `vina_pocket_bootstrap.csv` (optional)

### D. Phase 1 PyRosetta branch

Run PyRosetta PPI workflow and preserve traceable outputs by receptor state and
partner construct.

Main outputs include:

- `pyrosetta_run_metadata.json`
- `pyrosetta_decoy_scores.csv`
- `phase1_input_validation_report.json`
- `phase1_input_validation_summary.md`
- `ppi_pyrosetta_residues.csv`
- `ppi_pyrosetta_summary.csv`
- `ppi_cluster_summary.csv`
- `ppi_hotspot_residues.csv`
- `ppi_interface_patch_table.csv`

### E. LightDock secondary validation

Use LightDock as the active independent secondary validation axis for Phase 1.

Current LightDock outputs include:

- `lightdock_run_metadata.json`
- `lightdock_interface_support_table.csv`
- `lightdock_model_summary.csv`
- `cross_method_convergence.csv`

LightDock remains secondary evidence only. LightDock-only residues should not
be treated as primary patch truth without PyRosetta support.

### F. MD stability gate

Use MD outputs to classify stability before advancing to stronger downstream
interpretation.

### G. Integration and reporting

Run:

- `python main.py verdict -c config/example-project.yaml`
- `python main.py report -c config/example-project.yaml`
- `python main.py validate -c config/example-project.yaml`

Main outputs:

- `cross_method_agreement.csv`
- `valid_sites.csv`
- `vina_consensus_sites.csv`
- `project_report.txt`
- `combined_residue_evidence.csv`

## 6. Interpretation Rules

- Preserve receptor-state separation
- do not hard-code historical site names
- trust newly generated structured outputs over old residue labels
- treat LightDock as supporting method-independence evidence
- treat AFM as inactive unless explicitly re-enabled

## 7. Common Failure Mode To Avoid

The main documentation failure to avoid is using older AFM-heavy planning notes
as if they were the active baseline.

Current rule:

- use LightDock for active Phase 1 secondary validation
- do not plan normal execution around AFM
