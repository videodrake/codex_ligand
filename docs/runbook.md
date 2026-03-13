# Runbook

Last updated: 2026-03-13

This document is the operator-facing run procedure for the current repository. Use it to understand what order to run things in, which checkpoints matter, and when to stop or escalate. For exact commands, use [manual_execution.md](manual_execution.md). For the current baseline summary, use [current_pipeline_status.md](current_pipeline_status.md). For artifact meaning, use [output_artifact_map.md](output_artifact_map.md).

## Operating Frame

Run this repository with the current baseline in mind:

| Topic | Operational rule |
|------|------|
| Receptor states | Treat `3GT8_raw`, `3GT8_cl38_48`, and `3GT8_cl85_100` as the fixed comparison set |
| Ligand workflow center | Routine ligand evidence is Vina-centered |
| Phase 1 primary evidence | PyRosetta |
| Phase 1 secondary validation | LightDock |
| AFM | Legacy optional only; do not include it in routine execution unless explicitly re-enabled |
| Runtime environment | Production and pre-qsub lanes reuse the shared `pyrosetta` conda environment |
| Worker policy | Treat `max_workers = 16` as the routine safe operating bound |

## Interpretation Start Point

When a production run has completed, start with `output/{project}/step_index.md`.

Recommended reading order:

1. `output/{project}/step_index.md`
2. `output/{project}/step6_report/project_report.txt`
3. `output/{project}/step5_verdict/valid_sites.csv`
4. `output/{project}/step4_vina_postprocess/vina_pocket_table.csv`
5. `output/{project}/step3_ppi_postprocess/ppi_pyrosetta_residues.csv`

Canonical runtime outputs remain under `output/{project}/` and remain the source of truth. The step folders are a derived interpretation view that can be regenerated from canonical outputs.

## Reference Output Layout

```text
output/egfr_myo1d_vina/
  vina_pose_table.csv
  vina_pocket_table.csv
  vina_drug_pocket_map.csv
  ppi_pyrosetta_residues.csv
  ppi_pyrosetta_summary.csv
  valid_sites.csv
  cross_method_agreement.csv
  combined_residue_evidence.csv
  project_report.txt
  step_index.md
  current_run_manifest.json
  step1_vina_raw/
  step2_ppi_raw/
  step3_ppi_postprocess/
  step4_vina_postprocess/
  step5_verdict/
  step6_report/
  step7_validate/
```

## Before You Run

Confirm these conditions before submitting heavier work:

1. The active config is the intended project config, normally `config/example-project.yaml` or a direct derivative.
2. The three receptor states are present and still mapped explicitly.
3. Required ligands and prepared inputs are available for the intended lane.
4. AFM-dependent fields are not being treated as active requirements.
5. Worker count does not exceed the routine safe bound without an explicit reason.
6. You know whether you are running the routine Vina-centered baseline, a Phase 1-focused branch, or both.

If you need the exact shell or PBS commands for these checks, use [manual_execution.md](manual_execution.md).

## Standard Operator Sequence

### 1. Run Pre-qsub Validation First

Always run the lightweight precheck lane before the heavier production submission path.

Required checkpoint:

- `output/pre_qsub_status/last_pass.json`

If the precheck does not pass, stop and fix the configuration, input registration, or environment issue before moving on.

### 2. Choose The Execution Lane

Pick the run lane that matches the task instead of assuming every run should execute the full repository.

| Lane | Use when | Primary outputs to review first |
|------|------|------|
| Routine baseline lane | You need the current default ligand-facing evidence flow | `output/egfr_myo1d_vina/results/valid_sites.csv`, `cross_method_agreement.csv`, `project_report.txt` |
| Phase 1 lane | You need receptor-side interface evidence or Phase 2 handoff material | `output/phase1_ppi/phase1_downstream_patch_reference.csv`, `phase1_interface_report.md` |
| Combined review lane | You need routine baseline outputs plus current Phase 1 evidence for interpretation | Routine baseline result files plus the Phase 1 handoff and review artifacts |

For the exact command surface of each lane, use [manual_execution.md](manual_execution.md).

### 3. Submit Heavy Work In A Safe Order

The default safe pattern is:

1. Run pre-qsub validation.
2. Submit the production lane only after precheck success.
3. Run or review the Phase 1 PyRosetta branch when receptor-side evidence is required.
4. Run LightDock only as the secondary validation path for Phase 1, not as a replacement for PyRosetta.
5. Run verdict, report, and validate after the routine baseline outputs are available.

Important rule:

- Do not treat the scientific Phase 1-4 documents as proof that the whole repository should always run in a single unified end-to-end chain.

## Required Checkpoints

Review these checkpoints before moving to interpretation.

| Stage | Required checkpoint | Stop condition |
|------|------|------|
| Precheck | `output/pre_qsub_status/last_pass.json` exists and reflects a pass | Missing or failed precheck |
| Routine Vina postprocess | Core Vina tables are populated | Missing pose or pocket summary tables |
| Routine integration | `valid_sites.csv`, `cross_method_agreement.csv`, `project_report.txt` exist | Final decision files missing or obviously stale |
| Phase 1 PyRosetta | Phase 1 residue, patch, and review outputs exist | No structured Phase 1 exports for downstream use |
| Phase 1 LightDock | Cross-method convergence outputs exist when LightDock was requested | LightDock run incomplete but being cited as supporting evidence |

For exact file names and which ones are handoff artifacts, use [output_artifact_map.md](output_artifact_map.md).

## Step Folder Mode Semantics

`run_production.py` keeps step-folder behavior aligned with canonical phase semantics:

- `--status`: read-only for the step layer. It reports canonical phase status and derived step status without regenerating step folders, `step_index.md`, or `current_run_manifest.json`.
- `--from N`: reruns canonical phases `N` and above. Earlier step folders remain untouched, while steps `N` through `7` are treated as stale until rebuilt.
- `--only N[,M]`: rebuilds only the explicitly selected phases and matching step folders. Unselected steps are left unchanged.
- `--force`: reruns the selected scope even when canonical outputs already exist and refreshes the matching derived step views.
- Fresh run: there is no dedicated fresh-run flag in `run_production.py`. Use `python scripts/reset_production_outputs.py --execute` to clear the old project output root before a clean production rerun. This removes the step folders and root step files together with the canonical project outputs.

Operational cautions:

- Step folders are additive derived views, not replacements for canonical files.
- `step2_ppi_raw/` records PyRosetta raw run paths and metadata, but it does not duplicate large raw directories.

## Interpretation Rules During Operation

- Preserve receptor-state separation in every review step.
- Treat PyRosetta as the primary Phase 1 structural evidence layer.
- Treat LightDock as independent secondary validation, not as standalone primary truth.
- Treat AFM as inactive unless the user explicitly asks to re-enable it.
- Treat `verdict`, `report`, and `validate` as the routine final interpretation layer for the default baseline.
- Do not promote advanced Phase 4 perturbation outputs as the default final layer unless the task is explicitly Phase 4-oriented.

## Common Operator Mistakes

- Running with historical AFM expectations even though AFM is not in the routine baseline.
- Assuming production stage numbers and scientific Phase 1-4 numbers mean the same thing.
- Treating pointer stub files in `output/egfr_myo1d_vina/` as the actual payload files.
- Running more than 16 workers by default because the machine exposes more cores.
- Collapsing the three receptor states too early in summaries or handoff interpretation.

## When To Stop And Escalate

Stop and resolve the issue before continuing if any of these are true:

- pre-qsub validation fails
- receptor state registration is incomplete or ambiguous
- the run depends on AFM inputs that are currently unset
- LightDock is being cited without PyRosetta support in a context that expects primary evidence
- the expected result files are missing but downstream interpretation is already starting
- output files appear to be stale pointers rather than current payloads

## Use These Docs Next

- [manual_execution.md](manual_execution.md): exact commands and execution surfaces
- [architecture.md](architecture.md): data-flow map and handoff structure
- [data_inventory.md](data_inventory.md): physical input and output locations
- [output_artifact_map.md](output_artifact_map.md): artifact meaning, priority, and downstream consumption
- [current_vs_plan_matrix.md](current_vs_plan_matrix.md): current implementation gaps relative to the 4-phase plan
