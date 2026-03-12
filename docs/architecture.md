# Architecture Overview

Last updated: 2026-03-12

This document describes the current active architecture.
If older documents disagree, use this file and the current code as the baseline.

## Active Entry Points

Primary entry points:

- `python main.py`
- `python main.py vina`
- `python main.py postprocess`
- `python main.py verdict`
- `python main.py report`
- `python main.py validate`
- `python main.py full`
- `python main.py pyrosetta`
- `qsub config/run_pre_qsub_checks.pbs`
- `qsub config/run_production.pbs`

## High-Level Data Flow

```text
config/example-project.yaml
        |
        +--> Vina branch
        |      dock.py
        |      -> parse_poses.py
        |      -> contacts.py
        |      -> cluster.py
        |      -> summarize.py
        |      -> compare.py
        |      -> bootstrap.py (optional)
        |      -> vina_pose_table.csv
        |      -> vina_pocket_table.csv
        |      -> vina_drug_pocket_map.csv
        |      -> vina_pocket_comparison.csv
        |
        +--> Phase 1 PPI branch
        |      pyrosetta_docking/pipeline_manager.py
        |      -> ppi/postprocess_ppi.py
        |      -> ppi/pyrosetta_extract.py
        |      -> phase1/cluster_consensus.py
        |      -> phase1/compare_states.py
        |      -> phase1/lightdock_validation.py
        |      -> phase1/review_report.py
        |      -> ppi_interface_patch_table.csv
        |      -> ppi_cluster_summary.csv
        |      -> ppi_hotspot_residues.csv
        |      -> cross_method_convergence.csv
        |
        +--> MD branch
        |      md/gromacs_analysis.py
        |      md/ligand_contacts.py
        |
        +--> Integration branch
               verdict.py
               -> cross_method_agreement.csv
               -> valid_sites.csv
               -> vina_consensus_sites.csv
               report.py
               -> project_report.txt
               -> combined_residue_evidence.csv
               validate.py
```

## Current Evidence Roles

- Vina branch: primary ligand and pocket evidence
- PyRosetta branch: primary Phase 1 receptor-side PPI evidence
- LightDock: active secondary Phase 1 validation axis
- MD: downstream stability gate
- verdict/report/validate: integration and interpretation support

## AFM Status In Architecture

`egfr_pipeline/ppi/afm_extract.py` still exists, but AFM is not part of the
active routine architecture.

That means:

- AFM is not an active branch in the default data-flow diagram
- AFM should not be treated as the current secondary-validation baseline
- LightDock is the active secondary-validation path for Phase 1

## Package Map

### `egfr_pipeline/vina/`

- `dock.py`: Vina execution
- `parse_poses.py`: pose-level parsing into `vina_pose_table.csv`
- `contacts.py`: receptor contact extraction
- `cluster.py`: pocket assignment
- `summarize.py`: pocket and ligand-pocket summaries
- `compare.py`: cross-receptor pocket comparison
- `bootstrap.py`: optional pocket stability resampling

### `egfr_pipeline/pyrosetta_docking/`

- `pipeline_manager.py`: main Phase 1 PyRosetta execution manager
- `docking.py`: docking stages
- `analysis.py`: score and interface analysis
- `metadata.py`: run metadata and output naming helpers

### `egfr_pipeline/ppi/`

- `prepare_dimer_pdb.py`: receptor/partner preparation
- `postprocess_ppi.py`: restored-output registration and postprocessing
- `pyrosetta_extract.py`: PyRosetta residue extraction
- `afm_extract.py`: legacy optional AFM parser, not active baseline
- `submit.py`: PBS helper logic

### `egfr_pipeline/phase1/`

- `cluster_consensus.py`: Phase 1 cluster and hotspot consensus
- `compare_states.py`: cross-state comparison
- `lightdock_validation.py`: LightDock setup, extraction, and convergence
- `review_report.py`: Phase 1 review and downstream handoff

### `egfr_pipeline/md/`

- `gromacs_analysis.py`: MD analysis helpers
- `ligand_contacts.py`: ligand-contact analysis

### Root-level integration

- `main.py`: unified CLI
- `run_production.py`: production orchestration
- `egfr_pipeline/verdict.py`: site evidence integration
- `egfr_pipeline/report.py`: human-readable report generation
- `egfr_pipeline/validate.py`: schema and handoff validation
