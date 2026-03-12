# Current Pipeline Status

Last updated: 2026-03-12

This is the canonical current-state document for the repository.
If any older document disagrees with this file, use this file plus the current
code as the source of truth.

## Active Scientific Baseline

The repository currently operates as a state-comparison pipeline for the
EGFR-MYO1D project.

The active receptor ensemble is fixed to exactly three receptor states:

- `3GT8_raw`
- `3GT8_cl38_48`
- `3GT8_cl85_100`

The active workflow is:

1. Vina-centered ligand docking and postprocessing
2. PyRosetta-centered Phase 1 PPI mapping
3. LightDock secondary validation for Phase 1
4. MD-based stability gate
5. verdict / report / validate integration

## Evidence Hierarchy

Current evidence roles are:

- Primary ligand evidence: Vina pose and pocket outputs
- Primary Phase 1 PPI evidence: PyRosetta
- Secondary independent Phase 1 validation: LightDock
- Downstream stability gate: MD
- Final integration layer: verdict, report, validate

## AFM Status

AlphaFold-Multimer is not part of the active routine workflow.

Important implications:

- `egfr_pipeline/ppi/afm_extract.py` still exists in the repository
- some older docs still mention AFM
- AFM must be treated as a legacy optional parser, not as an active planning baseline
- new implementation or interpretation work must not assume AFM is required

Unless the user explicitly asks to revive AFM support, future work should use:

- PyRosetta as the primary Phase 1 engine
- LightDock as the active secondary validation axis

## Current Document Priority

Read documents in this order:

1. `docs/current_pipeline_status.md`
2. `README.md`
3. `docs/project_context.md`
4. `docs/architecture.md`
5. `docs/runbook.md`
6. `config/README.md`
7. `docs/phase1_pyrosetta_execution_note.md`
8. `docs/phase1_ppi_handoff_note.md`
9. `docs/phase1_lightdock_validation_note.md`
10. `docs/phase1_output_chain_note.md`

## Historical Documents

The following classes of documents may contain older AFM-oriented assumptions:

- `docs/prd_phase_1_ppi_first_interface_mapping.md`
- `docs/tasks_phase_1_ppi_first_interface_mapping.md`
- `docs/manual_execution.md`
- `docs/EGFR_MYO1D_Pipeline_Technical_Document.md`
- other earlier planning notes that mention AFM as a normal part of Phase 1

These files are kept for history and context, not as the default planning
baseline.

## Current Code Reality

What is currently true in code:

- Vina branch is active and structured under `egfr_pipeline/vina/`
- PyRosetta Phase 1 execution is active under `egfr_pipeline/pyrosetta_docking/`
- Phase 1 downstream integration and review are active under `egfr_pipeline/phase1/`
- LightDock integration exists under `egfr_pipeline/phase1/lightdock_validation.py`
- AFM extraction code exists but is not the active secondary-validation path
- pre-qsub validation and production guard flow are implemented

## Operating Constraints

- Routine worker baseline: 16
- Do not infer production performance from the local Codex workspace
- Keep outputs traceable by receptor state, ligand, and construct metadata
- Do not let legacy labels outrank newly generated structured outputs

## Rule For Future AI Agents

When in doubt:

1. trust current code over stale prose
2. trust this file over older planning docs
3. treat AFM as inactive unless explicitly requested
4. treat LightDock as the active Phase 1 secondary validation path
