# Project Context

Last updated: 2026-03-12

This document explains why this repository exists, what scientific question it is trying to support, and which project assumptions should shape future documentation and implementation work. It is not the read-order guide, the operator runbook, or the detailed data-flow map. Start with [AI_START_HERE.md](AI_START_HERE.md) for onboarding order, [current_pipeline_status.md](current_pipeline_status.md) for the present baseline, and [architecture.md](architecture.md) for flow-level structure.

## Project Identity

This repository is an EGFR-MYO1D state-comparison research pipeline. Its purpose is to turn docking, receptor-side interface mapping, and downstream interpretation into a traceable workflow whose outputs stay comparable across receptor states and across evidence types.

The repository is therefore not a generic docking toolkit and not a collection of one-off analysis scripts. It is a structured research system for asking whether receptor state changes alter pocket behavior, interface evidence, and perturbation relevance in ways that remain reviewable after the run is finished.

## Scientific Question

The central question is how the EGFR kinase-domain surface behaves across receptor states in the context of MYO1D interaction and ligand perturbation. In practice, the project tries to connect four kinds of evidence:

- ligand pose and pocket behavior
- receptor-side interface patch evidence
- cross-state agreement or divergence
- downstream perturbation-oriented interpretation

This framing matters because the project is not only looking for strong docking scores. It is trying to preserve enough structure in the outputs to support state-specific comparison, mechanistic review, and later decision-making.

## Fixed Scope Assumptions

The current research scope is intentionally narrow.

| Topic | Current project assumption |
|------|------|
| Receptor ensemble | Exactly three states: `3GT8_raw`, `EGFR_160-185`, `EGFR_170-200` |
| Comparison unit | State-specific evidence should remain explicit through the pipeline |
| Ligand evidence center | Vina-centered outputs remain the current routine backbone |
| Receptor-side Phase 1 evidence | PyRosetta is the primary structural source |
| Secondary Phase 1 validation | LightDock is the active independent validation path |
| AFM status | Legacy optional support only, not part of the routine scientific baseline |
| MD role | Downstream stability gate, not the first-trust onboarding surface |

These assumptions should be treated as the current project frame unless a user explicitly asks to reopen the scope.

## Why The Repository Exists

The problem this repository addresses is not the lack of computational methods. The problem is the lack of stable, comparable, and reviewable outputs when multiple receptor states, multiple evidence layers, and multiple analysis passes are involved.

The repository exists to standardize:

- how receptor states are named and compared
- how ligand-side and receptor-side evidence are kept traceable
- how handoff artifacts feed later phases
- how final summaries remain tied back to machine-readable outputs

Without this structure, it becomes too easy for historical labels, ad hoc notebooks, or one-time interpretation choices to outrank the actual run evidence.

## Evidence Philosophy

This project is evidence-driven rather than conclusion-driven. The pipeline should preserve enough intermediate structure that later reviewers can see where a claim came from.

That means the repository should keep:

- pose-level evidence before aggressive summarization
- pocket-level summaries that still preserve state identity
- receptor-side residue and patch evidence
- cross-method agreement where it exists
- cross-state comparison tables rather than single-state claims only

It should avoid treating every overlap, patch match, or classification label as settled truth without the supporting artifacts.

## Interpretation Boundaries

Several constraints keep the project grounded:

- Historical residue labels and older site names are reference material, not authoritative truth.
- Legacy AFM-oriented material should not define the current scientific baseline.
- Planning documents can describe a desirable future end-to-end phase system without proving that the default operational path already behaves that way.
- State comparison is part of the core question, so outputs that collapse receptor identity too early are lower-value artifacts.

## What Counts As A Good Output

A useful output in this project usually has most of these properties:

- it is tied to one or more explicit receptor states
- it preserves machine-readable provenance
- it can be handed to the next stage without manual reinterpretation
- it can be reviewed by a human without re-running the whole analysis
- it helps distinguish orthosteric, rim, allosteric, or low-relevance behavior more clearly than raw scores alone

This is why the project values structured CSV handoffs and review reports together instead of choosing only one style of artifact.

## What This Repository Is Not

This repository is not:

- a web application
- a SaaS product
- a generic molecular modeling framework
- an AFM-first workflow
- a repo that should be restarted from scratch whenever the design evolves

The preferred development style is to inspect the current codebase, preserve the working baseline, and improve traceability in small safe steps.

## Use These Docs For The Rest

- [AI_START_HERE.md](AI_START_HERE.md): onboarding order, trust hierarchy, and conflict rules
- [current_pipeline_status.md](current_pipeline_status.md): short summary of the current baseline
- [architecture.md](architecture.md): current data flow and package responsibilities
- [data_inventory.md](data_inventory.md): current input and output locations
- [output_artifact_map.md](output_artifact_map.md): meaning and priority of major artifacts
- [current_vs_plan_matrix.md](current_vs_plan_matrix.md): differences between planned behavior and current implementation
- [runbook.md](runbook.md): operator-facing execution procedure
