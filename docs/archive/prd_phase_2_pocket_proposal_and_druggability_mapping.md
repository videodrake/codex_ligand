## Context Summary
- Project: EGFR–MYO1D interface-centric perturbation discovery pipeline
- Current Subphase: Phase 2 PRD
- Purpose: Enumerate ligandable pockets and classify them relative to the Phase 1 MYO1D interface patch
- Primary Tools: fpocket, P2Rank, optionally FTMap-like hotspot mapping

---

# PRD — Phase 2: Pocket Proposal and Druggability Mapping

## Goal
Generate a receptor-state-specific catalog of candidate ligandable pockets and classify each candidate pocket by its structural relationship to the Phase 1 receptor-side MYO1D patch.

## Why this phase exists
A receptor-side PPI patch does not automatically imply a ligandable pocket. This phase separates:
- PPI importance
nfrom
- small-molecule tractability.

The goal is to enumerate candidate pockets first, then determine whether they are:
- directly overlapping the PPI patch,
- sitting at the interface rim,
- or acting as possible allosteric modulators.

## Inputs
- Phase 1 receptor-side interface patch outputs rooted in PyRosetta primary mapping, with LightDock secondary validation support and AFM only as optional legacy auxiliary context if present
- Three receptor-state structures
- Pocket proposal tools and/or externally supplied pocket predictions

## Phase 1 Handoff Contract
- The structured Phase 1 handoff into Phase 2 should preserve `construct_type` and `orientation_validation_status` in the machine-readable patch reference.
- Current operational default: Phase 2 remains in compatibility mode while legacy and skipped-filter runs still exist.
- `construct_type` should be validated and carried through normalization because downstream pocket interpretation depends on the receptor construct context.
- `orientation_validation_status` should be preserved as-is. Calibrated classes (`pass`, `fail`, `ambiguous`) are preferred, but `not_available` is still allowed for legacy or skipped-filter runs and must surface as a compatibility warning rather than a silent normalization.

## Core Requirements

### Feature 1. Multi-tool pocket proposal
Support candidate pocket generation from more than one method, such as:
- fpocket
- P2Rank
- future hotspot mapping tools

### Feature 2. Candidate pocket normalization
Merge or normalize closely overlapping candidate pockets within each receptor state.

### Feature 3. Patch relationship classification
For every candidate pocket, classify its relation to the Phase 1 PPI patch as:
- orthosteric candidate
- rim candidate
- allosteric candidate
- structurally distant / low-relevance candidate

### Feature 4. Pocket metadata persistence
Store candidate pocket metadata in a standardized format, including:
- receptor_id
- candidate_pocket_id
- centroid
- optional box sizes
- proposal source
- proposal score
- patch relationship class

### Feature 5. Druggability-oriented summary
Generate a summary that separates:
- geometric pocket presence
- predicted ligandability
- and PPI relevance

## User Story
As the researcher, I want to enumerate candidate ligandable pockets before running diversity-aware docking, so that docking budget is spent on scientifically relevant pockets rather than repeatedly rediscovering dominant pockets from a giant blind box.

## Acceptance Criteria
- [ ] Candidate pockets can be generated or ingested for each receptor state.
- [ ] Pocket proposals from multiple methods can be merged or normalized.
- [ ] Each candidate pocket is assigned a structured receptor-local ID.
- [ ] Each candidate pocket is classified relative to the Phase 1 PPI patch.
- [ ] Candidate pockets are stored in a machine-readable file.
- [ ] Pocket catalogs remain receptor-specific and are not mixed globally.
- [ ] Phase 1 patch handoff validation preserves `construct_type` and `orientation_validation_status`.
- [ ] Compatibility-mode ingestion reports when `orientation_validation_status` is `not_available` instead of treating it as an automatic failure.

## Primary Outputs
- `candidate_pockets.csv`
- `candidate_pocket_merge_table.csv`
- `pocket_patch_relationship.csv`
- `druggability_proposal_summary.csv`

## Non-Goals for Phase 2
- Running the final ligand docking campaign
- Ranking perturbation candidates
- Declaring final MYO1D-disruptive sites

## Open Questions
- Which pocket proposal methods should be mandatory in the MVP?
- What merge threshold should define overlapping candidate pockets?
- How conservative should orthosteric vs rim vs allosteric classification be?

---

## Korean Summary

이 Phase의 목표는 **candidate pocket catalog**를 만드는 것이다. fpocket/P2Rank 같은 포켓 탐지 결과를 합쳐서 receptor별 포켓 목록을 만들고, 각 포켓이 Phase 1의 MYO1D patch와 직접 겹치는지, rim인지, allosteric인지 분류한다.

