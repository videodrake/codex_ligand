# Architecture Overview

Last updated: 2026-03-13

This document describes the current data-flow architecture of the repository. It is not the run procedure guide and it is not the plan-vs-gap document. Use [runbook.md](runbook.md) for operator steps, [current_vs_plan_matrix.md](current_vs_plan_matrix.md) for plan differences, and [output_artifact_map.md](output_artifact_map.md) for artifact-by-artifact meaning.

## Architecture Snapshot

The repository currently has two architectural layers that coexist:

1. A routine Vina-centered baseline flow used by the default CLI and production path
2. A newer phase-separated scientific flow for Phase 1 through Phase 4

These layers share concepts and some outputs, but they are not yet one unified default execution path.

## Operational Phase Vs Step Interpretation Layer

`run_production.py` now exposes an additive step output layer on top of the canonical project root.

| Operational phase | Canonical root role | Derived step folder | Read first |
|------|------|------|------|
| Phase 1 | raw Vina pose `.pdbqt` generation | `step1_vina_raw/` | `raw_pose_index.csv` |
| Phase 2 | PyRosetta ranking outputs | `step2_ppi_raw/` | `TH1_final_ranking.csv`, `beta_meander_final_ranking.csv` |
| Phase 3 | PPI residue extraction | `step3_ppi_postprocess/` | `ppi_pyrosetta_residues.csv` |
| Phase 4 | Vina pocket summarization | `step4_vina_postprocess/` | `vina_pocket_table.csv` |
| Phase 5 | verdict integration | `step5_verdict/` | `valid_sites.csv` |
| Phase 6 | report generation | `step6_report/` | `project_report.txt` |
| Phase 7 | validation | `step7_validate/` | `validation_status.json` |

Canonical runtime outputs remain under `output/{project}/` and stay the runtime source of truth. `egfr_pipeline/output_steps.py` builds the step folders, `step_index.md`, and `current_run_manifest.json` as a derived interpretation view. Large raw PyRosetta directories are referenced by manifest rather than duplicated.

## Layer 1: Routine Baseline Flow

This is the current default operational architecture.

```text
config/example-project.yaml
        |
        +--> input/receptors/*.pdb
        +--> input/ligands/*.sdf + runtime-prepared *.pdbqt
        +--> ppi.pyrosetta_result_dirs / legacy PPI registrations
        |
        +--> Vina execution branch
        |      egfr_pipeline/vina/dock.py
        |      -> output/{project}/raw pose files
        |
        +--> Vina postprocess branch
        |      parse_poses.py
        |      -> contacts.py
        |      -> cluster.py
        |      -> summarize.py
        |      -> compare.py
        |      -> bootstrap.py (optional)
        |      -> output/{project}/
        |         vina_pose_table.csv
        |         vina_pocket_table.csv
        |         vina_drug_pocket_map.csv
        |         vina_pocket_comparison.csv
        |
        +--> PPI support branch
        |      pyrosetta_docking/pipeline_manager.py
        |      -> ppi/postprocess_ppi.py
        |      -> ppi/pyrosetta_extract.py
        |      -> output/{project}/
        |
        +--> Integration branch
               verdict.py
               -> output/{project}/cross_method_agreement.csv
               -> output/{project}/valid_sites.csv
               -> output/{project}/vina_consensus_sites.csv
               report.py
               -> output/{project}/project_report.txt
               -> output/{project}/combined_residue_evidence.csv
               validate.py
               -> output_steps.py
               -> output/{project}/step1_vina_raw ... step7_validate
               -> output/{project}/step_index.md
               -> output/{project}/current_run_manifest.json
```

Key point:

- This layer is why the repo is still best described as Vina-centered in routine operation.

## Layer 2: Phase-Separated Scientific Flow

This is the newer architecture reflected in `egfr_pipeline/phase1` through `egfr_pipeline/phase4`.

```text
input/PPI/phase1/
        |
        +--> Phase 1: receptor-side patch definition
        |      prepare_inputs.py
        |      -> launch_docking.py / pyrosetta_docking/*
        |      -> extract_interface.py
        |      -> orientation_filter.py
        |      -> cluster_consensus.py
        |      -> compare_states.py
        |      -> lightdock_validation.py
        |      -> review_report.py
        |      -> output/phase1_ppi/
        |         phase1_downstream_patch_reference.csv
        |
        +--> Phase 2: candidate pocket proposal
        |      patch_ingestion.py
        |      -> pocket_proposal.py
        |      -> pocket_merge.py
        |      -> patch_relationship.py
        |      -> druggability_confidence.py
        |      -> cross_state_alignment.py
        |      -> phase3_export.py
        |      -> output/phase2_pockets/
        |         phase3_candidate_pocket_reference.csv
        |
        +--> Phase 3: diversity-aware docking
        |      pocket_reference_ingestion.py
        |      -> job_construction.py
        |      -> budget_policy.py
        |      -> run_diverse_docking.py
        |      -> pose_attribution.py
        |      -> diversity_validation.py
        |      -> phase4_export.py
        |      -> output/phase3_docking/
        |         phase4_docking_evidence_reference.csv
        |
        +--> Phase 4: perturbation ranking
               evidence_ingestion.py
               -> score_framework.py
               -> mechanistic_classification.py
               -> perturbation_scoring.py
               -> state_interpretation.py
               -> review_output.py
               -> final_report.py
               -> presentation_summary.py
               -> output/phase4_perturbation/
```

Key point:

- This layer reflects the scientific target structure more directly than the default operational path does.

## Current Cross-Branch Handoffs

These are the most important machine-readable handoff points in the current architecture.

| From | Handoff artifact | To |
|------|------|------|
| Phase 1 | `output/phase1_ppi/phase1_downstream_patch_reference.csv` | Phase 2 patch ingestion |
| Phase 2 | `output/phase2_pockets/phase3_candidate_pocket_reference.csv` | Phase 3 pocket reference ingestion |
| Phase 3 | `output/phase3_docking/phase4_docking_evidence_reference.csv` | Phase 4 evidence ingestion |
| Routine Vina + PPI baseline | `vina_pocket_table.csv` plus PPI residue exports | `verdict.py`, `report.py`, `validate.py` |

## Current Data Domains

| Domain | Main location | Architectural role |
|------|------|------|
| Receptor and ligand inputs | `input/receptors/`, `input/ligands/` | Routine Vina baseline inputs |
| Structured Phase 1 inputs | `input/PPI/phase1/` | Current full-kinase-domain Phase 1 preparation layer |
| Legacy prepared PPI inputs | `input/PPI/prepared/` | Older dimer-centered preparation layer kept for provenance |
| Routine baseline outputs | `output/{project}/` | Default operational canonical output root plus additive derived step view |
| Structured Phase outputs | `output/phase1_ppi/` through `output/phase4_perturbation/` | Phase-separated scientific output roots |

## Package Responsibilities

| Package or module area | Responsibility |
|------|------|
| `egfr_pipeline/vina/` | Vina execution and postprocess chain for the routine ligand workflow |
| `egfr_pipeline/pyrosetta_docking/` | Core PyRosetta docking execution, scoring, and metadata |
| `egfr_pipeline/ppi/` | PPI preparation, chain restoration, residue extraction, and legacy AFM parsing |
| `egfr_pipeline/phase1/` | Structured Phase 1 interface mapping, filtering, convergence, and handoff |
| `egfr_pipeline/phase2/` | Pocket proposal, merge logic, patch relationship, druggability, and Phase 3 export |
| `egfr_pipeline/phase3/` | Pocket-guided docking setup, budget control, diversity tracking, and Phase 4 export |
| `egfr_pipeline/phase4/` | Multi-phase evidence ingestion, mechanistic classification, scoring, and final review outputs |
| `egfr_pipeline/verdict.py` | Current routine baseline site-judgment layer |
| `egfr_pipeline/report.py` | Current routine baseline text-report layer |
| `egfr_pipeline/validate.py` | Current routine baseline schema and handoff validation |
| `main.py` | Unified interactive and CLI entry surface |
| `run_production.py` | Production orchestration for the routine baseline flow |

## Architectural Rules

- Treat Vina as the center of gravity for the current routine baseline.
- Treat PyRosetta as the primary Phase 1 receptor-side evidence layer.
- Treat LightDock as the active Phase 1 secondary validation branch.
- Treat AFM as a legacy optional parser path, not as an active architecture branch.
- Treat canonical root outputs as the runtime source of truth and the step folders as a derived interpretation view.
- Treat phase-separated outputs as real architectural components, but do not assume they drive the default CLI path unless the user explicitly works on them.

## Current Architectural Cautions

- `run_production.py` phase numbers are operational stages, not the same as the scientific Phase 1-4 structure.
- `output/egfr_myo1d_vina/` mixes payload areas and pointer stub files, so file location alone is not enough to infer artifact meaning.
- The routine baseline final-decision layer is still `verdict/report/validate`, not the advanced Phase 4 perturbation stack.

## Read Next

- [current_pipeline_status.md](current_pipeline_status.md): short current baseline summary
- [runbook.md](runbook.md): operator-facing run sequence
- [data_inventory.md](data_inventory.md): where current inputs and outputs physically live
- [output_artifact_map.md](output_artifact_map.md): which output files are handoff files, review files, or trace files
- [current_vs_plan_matrix.md](current_vs_plan_matrix.md): where this architecture still differs from the intended 4-phase plan
