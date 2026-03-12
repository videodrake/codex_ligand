# Execution Manual
## EGFR-MYO1D Pipeline

Last updated: 2026-03-12

## Status Note

This manual reflects the current active baseline:

- Vina-centered ligand workflow
- PyRosetta-centered Phase 1 PPI workflow
- LightDock as the active secondary Phase 1 validation axis
- MD as a downstream stability gate

AlphaFold-Multimer is not part of the active routine workflow.
AFM parser code may still exist in the repository, but it should be treated as
legacy optional support only.

## 0. Environment

Recommended server-side environment:

```bash
cd ~/codex_ligand
conda activate pyrosetta
```

Minimal packages for local validation:

```bash
pip install pyyaml numpy pandas scipy matplotlib pytest
```

Optional packages:

```bash
pip install vina rdkit
pip install MDAnalysis
```

## 1. Project Config

Use `config/example-project.yaml` as the baseline project config.

Important points:

- routine safe worker baseline is 16
- receptor states are fixed to the three active states
- ligand list is config-driven
- PyRosetta result directories are tracked under `ppi.pyrosetta_result_dirs`
- AFM config fields may still appear, but they are not part of the active default workflow

## 2. Recommended Reading Before Running

1. `docs/current_pipeline_status.md`
2. `README.md`
3. `docs/architecture.md`
4. `docs/runbook.md`
5. `config/README.md`

For Phase 1:

6. `docs/phase1_pyrosetta_execution_note.md`
7. `docs/phase1_ppi_handoff_note.md`
8. `docs/phase1_lightdock_validation_note.md`
9. `docs/phase1_output_chain_note.md`

## 3. Pre-qsub Validation

Before heavy server submission:

```bash
qsub config/run_pre_qsub_checks.pbs
```

Success marker:

```text
output/pre_qsub_status/last_pass.json
```

Safest chained submission:

```bash
PRECHECK_JOB=$(qsub config/run_pre_qsub_checks.pbs)
qsub -W depend=afterok:${PRECHECK_JOB} config/run_production.pbs
```

## 4. Main Command Paths

### Interactive

```bash
python main.py
```

### Non-interactive

```bash
python main.py vina -c config/example-project.yaml
python main.py postprocess -c config/example-project.yaml
python main.py verdict -c config/example-project.yaml
python main.py report -c config/example-project.yaml
python main.py validate -c config/example-project.yaml
python main.py full -c config/example-project.yaml
```

### PyRosetta / PBS

```bash
qsub config/run_ppi_test.pbs
qsub config/run_ppi_prod.pbs
```

## 5. Current Execution Order

### Step 1. Precheck

Run the lightweight validation lane.

### Step 2. Vina docking

Generate raw docking poses for the three receptor states and configured ligands.

### Step 3. Vina postprocess

Run:

- parse
- contacts
- cluster
- summarize
- compare
- bootstrap (optional)

Main outputs:

- `vina_pose_table.csv`
- `vina_pocket_table.csv`
- `vina_drug_pocket_map.csv`
- `vina_pocket_comparison.csv`
- `vina_pocket_bootstrap.csv` (optional)

### Step 4. PyRosetta Phase 1

Run the PyRosetta PPI workflow and preserve receptor/partner metadata.

Main outputs:

- `pyrosetta_run_metadata.json`
- `pyrosetta_decoy_scores.csv`
- `phase1_input_validation_report.json`
- `phase1_input_validation_summary.md`
- `ppi_pyrosetta_residues.csv`
- `ppi_pyrosetta_summary.csv`
- `ppi_cluster_summary.csv`
- `ppi_hotspot_residues.csv`
- `ppi_interface_patch_table.csv`

### Step 5. LightDock secondary validation

Use LightDock for independent secondary validation of Phase 1 patch evidence.

Main outputs:

- `lightdock_run_metadata.json`
- `lightdock_interface_support_table.csv`
- `lightdock_model_summary.csv`
- `cross_method_convergence.csv`

Important rule:

- LightDock remains secondary evidence only
- LightDock-only residues should not be promoted as primary Phase 2 patch truth by themselves

### Step 6. MD stability gate

Use MD outputs to classify state stability before stronger downstream claims.

### Step 7. Integration

Run:

```bash
python main.py verdict -c config/example-project.yaml
python main.py report -c config/example-project.yaml
python main.py validate -c config/example-project.yaml
```

Main outputs:

- `cross_method_agreement.csv`
- `valid_sites.csv`
- `vina_consensus_sites.csv`
- `project_report.txt`
- `combined_residue_evidence.csv`

## 6. Interpretation Rules

- preserve receptor-state separation
- do not hard-code old site labels
- trust current structured outputs over historical residue labels
- treat LightDock as the active secondary validation path
- treat AFM as inactive unless explicitly re-enabled

## 7. Common Mistake To Avoid

Do not use older AFM-heavy documents as the current execution baseline.

If there is any conflict:

1. trust `docs/current_pipeline_status.md`
2. trust the current code
3. trust the current LightDock notes over older AFM-oriented prose
