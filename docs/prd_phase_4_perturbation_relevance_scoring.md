## Context Summary
- Project: EGFR–MYO1D interface-centric perturbation discovery pipeline
- Current Subphase: Phase 4 PRD
- Purpose: Rank pockets and ligands by their likelihood of disrupting MYO1D attachment rather than by affinity alone
- Core Concept: Perturbation relevance > raw docking strength

---

# PRD — Phase 4: Perturbation Relevance Scoring

## Goal
Rank candidate pockets and ligand-supported sites by their **likelihood of disrupting MYO1D attachment** to EGFR, using a structured multi-axis evidence framework.

## Why this phase exists
The project’s real endpoint is not “best pocket” or “best affinity.”
The real endpoint is:
- which pocket is most relevant to MYO1D attachment disruption,
- and whether that disruption is likely orthosteric, interface-rim-mediated, or allosteric.

This phase converts prior outputs into biologically meaningful prioritization.

## Inputs
- Phase 1 receptor-side interface patch outputs
- Phase 2 candidate pocket and patch relationship outputs
- Phase 3 diversity-aware ligand docking outputs
- Cross-state receptor evidence

## Core Requirements

### Feature 1. Multi-axis perturbation scoring
Score candidate sites using at least four evidence axes:
1. PPI Interface Confidence
2. Druggability Confidence
3. Perturbation Relevance
4. State Robustness / Accessibility

### Feature 2. Mechanistic classification
Every high-priority site should be classified as one of:
- orthosteric disruptor candidate
- interface-rim modulator candidate
- allosteric modulator candidate
- ligandable but PPI-irrelevant candidate

### Feature 3. Separation of affinity from relevance
Affinity must contribute to ranking, but it must not dominate the final interpretation if a pocket is structurally irrelevant to the MYO1D patch.

### Feature 4. Structured final verdict outputs
The system must generate a final structured table that preserves:
- raw axis scores
- mechanistic class
- receptor-state support
- supporting ligand evidence
- key residue notes

### Feature 5. Review-first reporting logic
The final ranking must remain auditable.
It should support human review instead of hiding conclusions inside one opaque score.

## User Story
As the researcher, I want final candidate sites ranked by how likely they are to perturb MYO1D attachment, so that downstream validation effort is spent on biologically meaningful candidates rather than on generic strong binders.

## Acceptance Criteria
- [ ] Final ranking uses more than affinity alone.
- [ ] PPI patch information is directly incorporated into final ranking.
- [ ] Orthosteric, rim, and allosteric interpretations can be separated.
- [ ] Cross-state support affects ranking rather than being treated as an afterthought.
- [ ] Final outputs preserve raw score components and not only a collapsed verdict.
- [ ] A final review table can be generated for manual inspection.

## Primary Outputs
- `perturbation_candidate_table.csv`
- `perturbation_axis_scores.csv`
- `final_candidate_classes.csv`
- `integrated_phase4_report.md`

## Non-Goals for Phase 4
- Declaring absolute biological truth
- Replacing experimental validation
- Automatically rejecting all low-affinity allosteric candidates

## Open Questions
- What score weights should be used in the first release?
- How conservative should orthosteric classification be?
- Should allosteric candidates require stronger state-robustness evidence than rim candidates?

---

## Korean Summary

이 Phase의 목표는 포켓이나 리간드를 단순 affinity 순으로 정렬하는 것이 아니라, **MYO1D 부착 방해 가능성** 기준으로 순위화하는 것이다. 최종 산출물은 orthosteric / rim / allosteric / irrelevant 분류와 함께 각 축의 점수를 모두 보존해야 한다.

