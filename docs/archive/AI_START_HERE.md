# AI Start Here

This repository is an EGFR-MYO1D state-comparison research pipeline that combines a Vina-centered ligand evidence layer with a Phase 1 receptor-side PPI evidence layer. The current routine baseline is built around exactly three receptor states (`3GT8_raw`, `EGFR_160-185`, `EGFR_170-200`), uses PyRosetta as the primary Phase 1 structural engine, uses LightDock as the active secondary Phase 1 validation path, treats AlphaFold-Multimer only as legacy optional support, and feeds integrated interpretation through verdict, report, and validation outputs.

## Current Operational Baseline

| Topic | Current baseline |
|------|------|
| Project focus | EGFR-MYO1D state-comparison pipeline |
| Receptor states | `3GT8_raw`, `EGFR_160-185`, `EGFR_170-200` |
| Central evidence layer | Vina-centered docking, pose parsing, pocket summarization, and cross-state comparison |
| Phase 1 primary evidence | PyRosetta |
| Phase 1 active secondary validation | LightDock |
| AFM status | Legacy optional parser only, inactive unless explicitly re-enabled |
| MD status | Downstream stability gate, not the first-trust onboarding surface |
| Runtime environment | Shared `pyrosetta` conda environment for production and pre-qsub lanes |
| Routine worker policy | Treat `max_workers = 16` as the safe routine upper bound |

## Source-of-Truth Hierarchy

Use the repository in this order of trust.

1. Executable reality for current behavior
   - `main.py`
   - `run_production.py`
   - `config/example-project.yaml`
   - `config/run_pre_qsub_checks.pbs`
   - `config/run_production.pbs`
   - active modules under `egfr_pipeline/`
2. Current-state documents for interpretation of current behavior
   - `docs/current_pipeline_status.md`
   - `docs/architecture.md`
   - `docs/runbook.md`
   - `docs/project_context.md`
   - `config/README.md`
3. Planning documents for intended phase design
   - `docs/brief_egfr_myo_1_d_pipeline_v_2.md`
   - `docs/prd_phase_1_ppi_first_interface_mapping_v2.md`
   - `docs/tasks_phase_1_ppi_first_interface_mapping_v2.md`
   - `docs/prd_phase_2_pocket_proposal_and_druggability_mapping.md`
   - `docs/tasks_phase_2_pocket_proposal_and_druggability_mapping.md`
   - `docs/prd_phase_3_diversity_aware_pocket_guided_docking.md`
   - `docs/tasks_phase_3_diversity_aware_pocket_guided_docking.md`
   - `docs/prd_phase_4_perturbation_relevance_scoring.md`
   - `docs/tasks_phase_4_perturbation_relevance_scoring.md`
4. Supporting notes and manual references
   - `README.md`
   - Phase 1 execution notes under `docs/phase1_*`
   - `docs/manual_execution.md`
   - `docs/manual_vina.md`
   - `docs/manual_pyrosetta.md`
5. Historical material
   - older handoff docs
   - older AFM-first planning docs
   - older technical summary documents

Important rule: planning documents explain where the project is trying to go; they do not prove that the default CLI and production flow already behave that way today.

## Recommended Read Order

Read in this order when onboarding a new GPT.

1. `docs/AI_START_HERE.md`
2. `docs/current_pipeline_status.md`
3. `docs/project_context.md`
4. `docs/architecture.md`
5. `docs/runbook.md`
6. `config/README.md`
7. `config/example-project.yaml`
8. `README.md`
9. `docs/brief_egfr_myo_1_d_pipeline_v_2.md`
10. Phase PRD and task files for Phases 1-4
11. Current entry-point and validation code:
    - `main.py`
    - `run_production.py`
    - `egfr_pipeline/config.py`
    - `egfr_pipeline/validate.py`
    - `egfr_pipeline/report.py`
    - `egfr_pipeline/verdict.py`
12. Phase-specific implementation files only after the above

## What Not To Trust First

Do not use these as your first mental model of the repository.

- Any text that implies AFM is part of the default active Phase 1 workflow
- Any document that treats LightDock as absent or secondary only in theory
- Any planning doc that implies the default operational pipeline already runs as a clean end-to-end Phase 1 -> Phase 2 -> Phase 3 -> Phase 4 stack
- Older handoff and technical documents that predate the current PyRosetta + LightDock baseline
- Any environment note that assumes a separate test-only conda environment is the baseline
- Any performance assumption that normalizes worker counts above the routine safe bound of 16

## Current Position of Key Methods

| Method | Current position |
|------|------|
| Vina | Current center of gravity for ligand evidence and for the default output backbone |
| PyRosetta | Primary Phase 1 structural evidence for receptor-side PPI mapping |
| LightDock | Active secondary validation path for Phase 1 and source of cross-method convergence tables |
| AFM | Legacy optional parser in `egfr_pipeline/ppi/afm_extract.py`; not part of the routine baseline |
| MD | Downstream stability gate mentioned in current docs, but not the first place to anchor onboarding |

## Current Code Reality Check

These repository facts are directly visible in current code and config.

- `config/example-project.yaml` defines the three active receptor states, three ligands, `max_workers: 16`, and keeps `ppi.afm_models: null`
- `main.py` exposes the active default CLI around `vina`, `postprocess`, `pyrosetta`, `verdict`, `report`, `validate`, and `full`
- `run_production.py` still orchestrates the routine production flow around Vina, PPI postprocessing, verdict, report, and validation
- `egfr_pipeline/phase1/lightdock_validation.py` implements LightDock setup, extraction, and cross-method convergence outputs
- `egfr_pipeline/ppi/afm_extract.py` still exists, but it is an auxiliary parser path rather than the default active Phase 1 route
- `config/run_production.pbs`, `config/run_pre_qsub_checks.pbs`, and `scripts/setup_test_env.sh` all align to the shared `pyrosetta` environment rather than a separate baseline test-only env

## How To Resolve Document Conflicts

When two sources disagree, use this rule.

1. For current operational truth, trust executable code and active runtime config first.
2. Use `docs/current_pipeline_status.md`, `docs/architecture.md`, and `docs/runbook.md` to interpret the current baseline.
3. Use the brief and phase PRD/task files as target design documents, not as automatic proof of implemented behavior.
4. Treat AFM references as legacy unless the user explicitly asks to revive AFM support.
5. If code exposes an optional path but the active config leaves it unset or null, treat that path as inactive in the routine baseline.

## Practical Starting Assumptions For A New GPT

Start with these assumptions unless the user says otherwise.

- The project is currently operated as a Vina-centered evidence pipeline with a PyRosetta-based Phase 1 layer.
- The three receptor states are fixed and must remain explicitly separated in outputs.
- LightDock is the active Phase 1 secondary validation path.
- AFM is not part of the default routine plan.
- Production and pre-qsub lanes reuse the `pyrosetta` environment.
- `max_workers = 16` is the routine safe operating assumption.
