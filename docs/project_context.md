# Project Context

## EGFR-MYO1D Pipeline

## 1. Why this repository exists

This repository exists to standardize the computational workflow used in the
EGFR-MYO1D study.

The problem is not the absence of computation. The problem is the lack of
stable, comparable, and reviewable outputs across receptor states and across
evidence types.

This repository is therefore a structured research pipeline, not a collection
of one-off scripts.

## 2. Scientific framing

The current scientific focus is the EGFR kinase-domain surface in the context
of MYO1D interaction and perturbation.

The pipeline is designed to support:

- receptor-state-specific ligand docking behavior
- pose-level and pocket-level evidence
- receptor-side PPI patch definition
- cross-state comparison
- downstream perturbation-oriented ranking

This is not a generic docking project. It is a state-comparison research
pipeline.

## 3. Fixed receptor states

The current receptor ensemble is limited to exactly these three structures:

1. `3GT8_raw`
2. `3GT8_cl38_48`
3. `3GT8_cl85_100`

All major outputs and comparisons should preserve these state IDs explicitly.
Residue numbering consistency across them remains a high-priority concern.

## 4. Current computational baseline

The active baseline is:

1. Vina-centered ligand workflow
2. PyRosetta-centered Phase 1 PPI workflow
3. LightDock secondary validation for Phase 1
4. MD as a downstream stability gate

This means the repository is no longer AFM-centered for Phase 1.

## 5. AFM status

AlphaFold-Multimer is not part of the active routine workflow.

The repository still contains AFM parsing code, but that code should be treated
as legacy optional support, not as the default planning baseline.

In practice:

- do not treat AFM as part of the current core workflow
- do not create new requirements that depend on AFM unless explicitly requested
- do not let older AFM-oriented docs outrank the current LightDock-based Phase 1 baseline

## 6. Current evidence philosophy

The pipeline is evidence-driven, not conclusion-driven.

That means the system should preserve:

- raw pose-level evidence
- pocket-level summaries
- receptor-side interface evidence
- cross-state overlap metrics
- method-agreement information

It should not over-interpret every overlap or every patch relationship as fixed
truth.

## 7. Interpretation rule

Legacy residue labels, site names, and older report labels are historical
reference only.

The current repository should prioritize:

- newly generated pose data
- newly generated pocket assignments
- current receptor-state-specific overlap evidence
- current PyRosetta and LightDock outputs for Phase 1

Do not hard-code old site names into logic.

## 8. Operating constraint

The main execution environment has 32 CPU cores, but the routine safe operating
assumption is 16 workers.

That means:

- parallel execution must stay configurable
- 16 workers should be treated as the practical routine upper bound
- local Codex workspace behavior must not be treated as server-performance truth

## 9. Development style

This repository should be improved by refactoring and standardizing the current
codebase, not by discarding it and starting over.

Preferred style:

- inspect the current repo first
- preserve what works
- tighten traceability and schema consistency
- expand in small safe steps

## 10. Current output goals

The repository should provide structured outputs for:

- receptor metadata
- ligand metadata
- Vina pose-level parsed output
- Vina pocket-level summaries
- ligand-to-pocket mapping
- cross-receptor pocket comparison
- PyRosetta receptor-side residue summaries
- LightDock cross-method convergence outputs
- markdown summary reports
- validation outputs

## 11. What this repository is not

This is not:

- a web application
- a SaaS platform
- a cloud deployment project
- a generic docking toolkit

It is a focused computational research pipeline for EGFR-MYO1D analysis.

## 12. Read-first rule for future contributors

Anyone entering the repository should read, in order:

1. `docs/current_pipeline_status.md`
2. `README.md`
3. `docs/project_context.md`
4. `docs/architecture.md`
5. `docs/runbook.md`

Then, if they are working on Phase 1:

6. `docs/phase1_pyrosetta_execution_note.md`
7. `docs/phase1_ppi_handoff_note.md`
8. `docs/phase1_lightdock_validation_note.md`
9. `docs/phase1_output_chain_note.md`
