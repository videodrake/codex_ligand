# Current Pipeline Status

Last updated: 2026-03-13

This is a short derived summary of the repository's current baseline. Use it for quick orientation, not as the only source of truth. For onboarding order and conflict resolution, start with [AI_START_HERE.md](AI_START_HERE.md). For plan-vs-implementation gaps, use [current_vs_plan_matrix.md](current_vs_plan_matrix.md).

## Current Summary

The repository currently behaves as an EGFR-MYO1D state-comparison pipeline with a Vina-centered routine ligand workflow and a PyRosetta-centered Phase 1 PPI branch. The active receptor ensemble is fixed to three states, LightDock is the active secondary validation path for Phase 1, AFM is legacy optional only, and routine final integration still flows through `verdict`, `report`, and `validate` rather than the advanced Phase 4 perturbation stack.

## Current Interpretation Path

For completed production runs, start at `output/{project}/step_index.md`.

The current interpretation order is:

1. `step_index.md`
2. `step6_report/project_report.txt`
3. `step5_verdict/valid_sites.csv`
4. `step4_vina_postprocess/vina_pocket_table.csv`
5. `step3_ppi_postprocess/ppi_pyrosetta_residues.csv`

Canonical runtime outputs remain under `output/{project}/` and stay the source of truth. The step folders are a derived view layered on top of those canonical outputs; they are additive, not a migration.

## Active Baseline

| Topic | Current state |
|------|------|
| Receptor states | `3GT8_raw`, `EGFR_160-185`, `EGFR_170-200` |
| Routine ligand evidence | Vina docking and Vina postprocess outputs |
| Phase 1 primary evidence | PyRosetta |
| Phase 1 active secondary validation | LightDock |
| AFM status | Legacy optional parser only, inactive unless explicitly re-enabled |
| Runtime environment baseline | Shared `pyrosetta` conda environment |
| Routine worker policy | Treat `max_workers = 16` as the safe routine operating bound |

## Current Default Execution Surface

The default code-facing workflow currently centers on:

- `python main.py vina`
- `python main.py postprocess`
- `python main.py pyrosetta`
- `python main.py verdict`
- `python main.py report`
- `python main.py validate`
- `python main.py full`
- `qsub config/run_pre_qsub_checks.pbs`
- `qsub config/run_production.pbs`

Important interpretation:

- The default `full` flow is still Vina-centered.
- The refactored scientific Phase 1-4 stack exists in code and outputs, but it is not yet the single default orchestration path for routine runs.

## Current Evidence Hierarchy

- Primary ligand evidence: Vina pose and pocket outputs
- Primary Phase 1 PPI evidence: PyRosetta
- Secondary independent Phase 1 validation: LightDock
- Downstream stability gate: MD
- Current routine final integration: `verdict`, `report`, and `validate`

## What Is Not The Routine Baseline

- AFM is not part of the active routine workflow.
- The advanced `output/phase2_pockets/`, `output/phase3_docking/`, and `output/phase4_perturbation/` trees should not be assumed to drive the default CLI path unless the user explicitly works on those phases.
- `valid_sites.csv` is the routine baseline judgment table; it is not the same thing as the advanced Phase 4 perturbation-ranking outputs.

## Current Output Checkpoints

If you need the current baseline outputs first, open:

1. `output/egfr_myo1d_vina/step_index.md`
2. `output/egfr_myo1d_vina/step6_report/project_report.txt`
3. `output/egfr_myo1d_vina/step5_verdict/valid_sites.csv`

If you need the canonical files behind that derived view, open:

1. `output/egfr_myo1d_vina/vina_pocket_table.csv`
2. `output/egfr_myo1d_vina/valid_sites.csv`
3. `output/egfr_myo1d_vina/project_report.txt`

If you need the current structured Phase 1 handoff first, open:

1. `output/phase1_ppi/phase1_downstream_patch_reference.csv`
2. `output/phase1_ppi/phase1_interface_report.md`

For a fuller artifact map, use [output_artifact_map.md](output_artifact_map.md).

## Use These Docs Next

- [AI_START_HERE.md](AI_START_HERE.md): onboarding order, source-of-truth hierarchy, and conflict rules
- [data_inventory.md](data_inventory.md): current input and output inventory
- [current_vs_plan_matrix.md](current_vs_plan_matrix.md): where current implementation still diverges from the 4-phase plan
- [output_artifact_map.md](output_artifact_map.md): what each major artifact means and which files are handoff files
- [glossary_and_assumptions.md](glossary_and_assumptions.md): project-specific vocabulary, numbering, chain, and field semantics
- [architecture.md](architecture.md): current data flow
- [runbook.md](runbook.md): operator-facing execution procedure

## Interpretation Rule

When you need a quick rule for current work:

1. Trust active code and active config over stale prose.
2. Treat this file as a quick summary, not as the document that overrides all others.
3. Treat AFM as inactive unless explicitly requested.
4. Treat LightDock as the active Phase 1 secondary validation path.
