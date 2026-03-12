# EGFR-MYO1D Pipeline

This repository is a research pipeline for EGFR-MYO1D state-comparison work.
It combines a Vina-centered ligand workflow with a PyRosetta-centered Phase 1
PPI workflow across three fixed EGFR receptor states:

- `3GT8_raw`
- `3GT8_cl38_48`
- `3GT8_cl85_100`

## Current Baseline

The current active baseline is:

- Vina-centered ligand evidence for pocket and pose analysis
- PyRosetta as the primary Phase 1 PPI engine
- LightDock as the active secondary Phase 1 validation axis
- MD as a downstream stability gate

Important clarification:

- AlphaFold-Multimer is not part of the active routine workflow
- `egfr_pipeline/ppi/afm_extract.py` still exists as a legacy optional parser
- do not plan new work around AFM unless the user explicitly asks to re-enable it

## Read First

Before using or extending the project, read these files in order:

1. [docs/current_pipeline_status.md](/Users/admin/Desktop/hwang/codex/codex_ligand/docs/current_pipeline_status.md)
2. [docs/project_context.md](/Users/admin/Desktop/hwang/codex/codex_ligand/docs/project_context.md)
3. [docs/architecture.md](/Users/admin/Desktop/hwang/codex/codex_ligand/docs/architecture.md)
4. [docs/runbook.md](/Users/admin/Desktop/hwang/codex/codex_ligand/docs/runbook.md)
5. [config/README.md](/Users/admin/Desktop/hwang/codex/codex_ligand/config/README.md)

For Phase 1 PPI details, then read:

6. [docs/phase1_pyrosetta_execution_note.md](/Users/admin/Desktop/hwang/codex/codex_ligand/docs/phase1_pyrosetta_execution_note.md)
7. [docs/phase1_ppi_handoff_note.md](/Users/admin/Desktop/hwang/codex/codex_ligand/docs/phase1_ppi_handoff_note.md)
8. [docs/phase1_lightdock_validation_note.md](/Users/admin/Desktop/hwang/codex/codex_ligand/docs/phase1_lightdock_validation_note.md)
9. [docs/phase1_output_chain_note.md](/Users/admin/Desktop/hwang/codex/codex_ligand/docs/phase1_output_chain_note.md)

## Main Entry Points

- `python main.py`
- `python main.py vina -c config/example-project.yaml`
- `python main.py postprocess -c config/example-project.yaml`
- `python main.py verdict -c config/example-project.yaml`
- `python main.py report -c config/example-project.yaml`
- `python main.py validate -c config/example-project.yaml`
- `python main.py pyrosetta`
- `qsub config/run_pre_qsub_checks.pbs`
- `qsub config/run_production.pbs`

## Repository Structure

Key folders:

- `egfr_pipeline/vina/`: docking, pose parsing, contacts, pocket clustering, summaries, comparison
- `egfr_pipeline/pyrosetta_docking/`: PyRosetta Phase 1 docking and scoring
- `egfr_pipeline/phase1/`: Phase 1 consensus, cross-state comparison, LightDock validation, review report
- `egfr_pipeline/ppi/`: PPI preparation and residue extraction
- `egfr_pipeline/md/`: MD analysis helpers
- `config/`: YAML configs, INI configs, PBS submission files
- `input/`: receptor, ligand, and PPI inputs
- `output/`: generated outputs
- `tests/`: test suite including the pre-qsub lane

## Inputs

Current active inputs include:

- receptors: `input/receptors/*.pdb`
- ligands: `input/ligands/*.sdf` and matching PDBQT files
- PPI inputs: `input/PPI/`
- project config: `config/example-project.yaml`

## Core Outputs

Vina/postprocess outputs:

- `vina_pose_table.csv`
- `vina_pocket_table.csv`
- `vina_drug_pocket_map.csv`
- `vina_pocket_comparison.csv`
- `vina_pocket_bootstrap.csv` (optional)

Phase 1 / PPI outputs:

- `ppi_pyrosetta_residues.csv`
- `ppi_pyrosetta_summary.csv`
- `ppi_cluster_summary.csv`
- `ppi_hotspot_residues.csv`
- `ppi_interface_patch_table.csv`
- `ppi_patch_cross_state_comparison.csv`
- `ppi_patch_state_robustness.csv`
- `cross_method_convergence.csv`
- `phase1_downstream_patch_reference.csv`

Final integration outputs:

- `cross_method_agreement.csv`
- `valid_sites.csv`
- `vina_consensus_sites.csv`
- `project_report.txt`
- `combined_residue_evidence.csv`

## Operating Notes

- The main server has 32 CPU cores, but the routine operating baseline is 16 workers.
- Keep receptor states separated in all outputs.
- Do not hard-code legacy site names or old residue labels into logic.
- Treat old AFM-centric docs as historical unless they explicitly say otherwise.

## Historical Note

Some older documents still contain AFM-centered or AFM-available language.
Those documents are retained for historical context, but the current active
Phase 1 secondary validation path is LightDock.
