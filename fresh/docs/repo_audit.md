# Repository Audit Template

This template is for later donor-module review. It is not a completed audit.

## Goals

- Identify old workflow modules that may be useful as references.
- Classify each module as direct reuse, wrapper-only reuse, logic-only reuse, or do-not-reuse.
- Record dependencies, risks, and tests required before any future reuse.

## Guardrails

- Do not modify old workflow files during Milestone 1 Task 1.
- Do not treat old Workflow A/B outputs as fresh workflow source-of-truth.
- Do not copy confidential ligand identifiers into public tracked files.

## Pending Audit Areas

- PDB parsing and residue mapping utilities.
- PBS/qsub helper patterns.
- Logging and report patterns.
- Vina, fpocket, and PyRosetta code paths for later donor review only.
