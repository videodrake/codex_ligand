## Context Summary
- Project: EGFR–MYO1D interface-centric perturbation discovery pipeline
- Current Subphase: Phase 4 Task Breakdown
- Purpose: Convert the Phase 4 PRD into implementation-facing tasks for perturbation relevance scoring
- Upstream Dependencies: Phase 1 receptor-side patch reference, Phase 2 candidate pocket classification, and Phase 3 docking evidence package
- Key Principle: Final ranking must reflect MYO1D attachment disruption relevance, not affinity alone

---

# Task Breakdown
## Phase 4: Perturbation Relevance Scoring

This document breaks Phase 4 into implementation-facing task groups. The purpose of this phase is to transform structured outputs from Phases 1–3 into a final evidence-based ranking of candidate pockets and ligand-supported sites by their **likelihood of disrupting MYO1D attachment**.

Phase 4 is **not** a generic final scoring layer and **not** a simple affinity-ranking stage.
Its function is to answer the real biological question of the project:

> Which candidate sites are most plausible for disrupting MYO1D attachment to EGFR, and by what mechanism?

This phase must therefore preserve mechanistic meaning, traceability, and human interpretability.

Each task group below includes:
- Objective
- Priority
- Zone
- Main tasks
- Subtasks
- Test tasks
- Dependencies
- Deliverables

---

# Task Group 4.0: Multi-Phase Evidence Ingestion and Validation
**Priority:** Must-Have  
**Zone:** 🟢 Green

## Objective
Load and validate the structured evidence outputs from Phases 1–3 so that final scoring is based on complete, traceable, and internally consistent inputs.

## Main Tasks

### 4.0.1 Phase 1 evidence ingestion
- Load receptor-side patch reference
- Load patch robustness information
- Load patch confidence labels
- Load hotspot residue summaries

### 4.0.2 Phase 2 evidence ingestion
- Load candidate pocket reference
- Load patch relationship classes
- Load druggability support levels
- Load cross-state pocket class labels

### 4.0.3 Phase 3 evidence ingestion
- Load docking evidence reference
- Load ligand support information
- Load search-policy and budget provenance
- Load diversity/occupancy summaries if needed

### 4.0.4 Consistency validation
- Confirm receptor_id consistency across all phases
- Confirm candidate_pocket_id consistency across Phase 2 and Phase 3
- Confirm biological class fields are preserved and interpretable
- Confirm missing required fields are surfaced clearly

## Subtasks
- Define a unified Phase 4 evidence schema.
- Define validation rules for each upstream file.
- Add field-level compatibility checks.
- Record warnings for incomplete or partially missing evidence.

## Test Tasks
- Confirm all upstream evidence files can be loaded.
- Confirm receptor IDs, patch IDs, and candidate pocket IDs align correctly.
- Confirm missing evidence is surfaced as warnings or errors.
- Confirm final scoring does not silently proceed with broken mappings.

## Dependencies
Depends on successful completion of Phases 1–3.

## Deliverables
- `phase4_evidence_validation.md`
- `phase4_evidence_normalized.csv`

---

# Task Group 4.1: Multi-Axis Score Framework Definition
**Priority:** Must-Have  
**Zone:** 🟢 Green

## Objective
Define a final scoring framework that reflects MYO1D perturbation relevance rather than simple ligand affinity.

## Main Tasks

### 4.1.1 Define the core evidence axes
At minimum, define the following score axes:
1. **PPI Interface Confidence**
2. **Druggability Confidence**
3. **Perturbation Relevance**
4. **State Robustness / Accessibility**

### 4.1.2 Define score semantics
Each axis must have a clear interpretation.
For example:
- PPI Interface Confidence = how strong and reliable the receptor-side patch evidence is
- Druggability Confidence = how plausible the site is for small-molecule engagement
- Perturbation Relevance = how directly the site can plausibly alter MYO1D attachment
- State Robustness / Accessibility = whether the site is persistent, shifted, or state-dependent

### 4.1.3 Preserve raw subcomponents
Do not collapse everything into one opaque score too early.
Each axis should preserve enough raw submetrics to remain reviewable.

## Subtasks
- Define axis names and meanings.
- Define raw submetric mapping for each axis.
- Define optional normalization or scaling strategy.
- Define whether any axis is conditional or optional in early versions.

## Test Tasks
- Confirm each axis can be computed from existing upstream evidence.
- Confirm axis meanings are not redundant with one another.
- Confirm final scoring remains interpretable rather than black-box.
- Confirm low-affinity but biologically relevant sites are not automatically eliminated.

## Dependencies
Depends on Task Group 4.0.

## Deliverables
- `phase4_score_framework.md`
- `phase4_axis_definition_table.csv`

---

# Task Group 4.2: Mechanistic Classification Logic
**Priority:** Must-Have  
**Zone:** 🟢 Green

## Objective
Assign mechanistic interpretation labels to candidate sites so the final output distinguishes different modes of MYO1D perturbation.

## Main Tasks

### 4.2.1 Define mechanistic classes
Every sufficiently supported site should be classified as one of:
- `orthosteric_disruptor_candidate`
- `interface_rim_modulator_candidate`
- `allosteric_modulator_candidate`
- `ligandable_but_ppi_irrelevant_candidate`
- `uncertain_mechanism_candidate`

### 4.2.2 Map Phase 2 classes into Phase 4 interpretation
Use:
- patch relationship class
- hotspot overlap
- pocket support
- docking evidence
- receptor-state support
as inputs to final mechanistic labeling.

### 4.2.3 Preserve uncertainty
Do not force all candidates into strong mechanistic labels when evidence is weak or contradictory.

## Subtasks
- Define mechanistic classification rules.
- Define the minimum evidence required for each class.
- Define uncertainty fallback behavior.
- Define how contradictory evidence affects class assignment.

## Test Tasks
- Confirm a candidate can be classified as orthosteric, rim, allosteric, irrelevant, or uncertain.
- Confirm contradictory evidence does not silently default to a strong label.
- Confirm low-support candidates are not overcalled.
- Confirm mechanistic class can be traced back to upstream evidence.

## Dependencies
Depends on Task Group 4.1.

## Deliverables
- `final_candidate_classes.csv`
- `phase4_mechanistic_classification_note.md`

---

# Task Group 4.3: Perturbation-Relevance Scoring Logic
**Priority:** Must-Have  
**Zone:** 🟢 Green

## Objective
Implement the actual candidate ranking logic so final prioritization reflects MYO1D perturbation relevance, not just geometric or docking strength.

## Main Tasks

### 4.3.1 Combine axis scores into candidate-level ranking
- Compute per-candidate axis scores
- Preserve raw axis values
- Generate a final perturbation score or ranked evidence summary

### 4.3.2 Prevent affinity domination
- Ensure affinity is only one component of the ranking
- Do not allow high affinity at a biologically irrelevant site to dominate the final ranking

### 4.3.3 Preserve evidence provenance
For each final candidate, retain references to:
- receptor-side patch evidence
- patch relationship class
- pocket support level
- docking support count
- budget/search-policy provenance if relevant

## Subtasks
- Define scoring combination strategy.
- Define ranking output schema.
- Add hooks for later tuning of score weights.
- Preserve raw evidence columns in output.

## Test Tasks
- Confirm final ranking can be produced from upstream evidence.
- Confirm high-affinity but irrelevant pockets do not outrank biologically relevant moderate-affinity pockets by default.
- Confirm final ranking remains explainable from axis scores.
- Confirm provenance fields survive into final outputs.

## Dependencies
Depends on Task Groups 4.1 and 4.2.

## Deliverables
- `perturbation_axis_scores.csv`
- `perturbation_candidate_table.csv`
- `phase4_ranking_method_note.md`

---

# Task Group 4.4: State-Robustness and Accessibility Interpretation
**Priority:** Must-Have  
**Zone:** 🟡 Yellow

## Objective
Ensure final ranking properly reflects whether a candidate site is consistently accessible, conditionally accessible, or strongly state-dependent.

## Main Tasks

### 4.4.1 Cross-state support interpretation
- Interpret Phase 2 and Phase 3 state evidence jointly
- Distinguish robust, shifted, and state-specific pockets

### 4.4.2 Accessibility-aware scoring influence
- Allow state robustness to strengthen confidence
- Allow state-specificity to be preserved as potentially meaningful rather than automatically penalized

### 4.4.3 Cryptic/allosteric caution handling
- Preserve candidates that may only appear in certain receptor states
- Flag them appropriately rather than discarding them as weak by default

## Subtasks
- Define state-interpretation logic.
- Define how state robustness modifies overall confidence.
- Define how state-specific pockets should be represented in final ranking.

## Test Tasks
- Confirm state-robust pockets can be recognized.
- Confirm state-specific but biologically interesting pockets are retained.
- Confirm cryptic-like candidates remain visible rather than silently filtered out.

## Dependencies
Depends on Task Group 4.3.

## Deliverables
- `phase4_state_interpretation.csv`
- `phase4_accessibility_note.md`

---

# Task Group 4.5: Final Review-First Output Design
**Priority:** Must-Have  
**Zone:** 🟡 Yellow

## Objective
Design the final output tables so they remain useful for human review, presentation, and downstream scientific discussion.

## Main Tasks

### 4.5.1 Candidate review table design
Create a final review-oriented table that includes:
- candidate identifier
- mechanistic class
- axis scores
- ligand support summary
- receptor-state support summary
- key residues or patch notes

### 4.5.2 Raw-evidence preservation
Ensure that a final condensed table exists, but also preserve a deeper table where the full evidence trail is still accessible.

### 4.5.3 Traceability support
- Make it possible to trace every top-ranked site back to the upstream files and metrics that generated it.

## Subtasks
- Define review table schema.
- Define minimal and expanded output views.
- Add traceability fields for upstream evidence IDs.

## Test Tasks
- Confirm the final table can be read without opening all upstream files.
- Confirm top-ranked candidates remain traceable to raw evidence.
- Confirm the condensed view does not hide critical uncertainties.

## Dependencies
Depends on Task Groups 4.2, 4.3, and 4.4.

## Deliverables
- `phase4_final_review_table.csv`
- `phase4_expanded_evidence_table.csv`

---

# Task Group 4.6: Final Report and Interpretation Guide
**Priority:** Should-Have  
**Zone:** 🟡 Yellow

## Objective
Generate a final Phase 4 report that summarizes how the project reached its final candidate ranking and how the results should be interpreted.

## Main Tasks

### 4.6.1 Final candidate ranking summary
- Present the final ranked candidates
- Summarize their mechanistic classes
- Highlight top orthosteric, rim, and allosteric candidates separately

### 4.6.2 Interpretation guidance
- Explain why affinity alone was not used as the final criterion
- Explain how PPI patch relevance changed the ranking
- Explain how state robustness and uncertainty should be read

### 4.6.3 Validation and caution section
- Summarize what the final ranking can and cannot claim
- Preserve the distinction between computational prioritization and biological proof

## Subtasks
- Define final report sections.
- Add candidate summary tables.
- Add interpretation notes for mechanistic classes.
- Add caution statements about experimental validation needs.

## Test Tasks
- Confirm the report can be read without deep upstream file inspection.
- Confirm the report explains the rationale of the ranking clearly.
- Confirm uncertainty and non-finality are preserved appropriately.

## Dependencies
Depends on Task Groups 4.3, 4.4, and 4.5.

## Deliverables
- `integrated_phase4_report.md`

---

# Task Group 4.7: Presentation-Ready Summary Layer
**Priority:** Should-Have  
**Zone:** 🟡 Yellow

## Objective
Prepare a compact presentation-ready view of the final Phase 4 results for internal discussion, presentation, or publication planning.

## Main Tasks

### 4.7.1 Short summary table
- Generate a compact top-candidate table for presentation use

### 4.7.2 Mechanism-oriented shortlist
- Separate shortlists by mechanism class:
  - orthosteric
  - rim
  - allosteric

### 4.7.3 Confidence-oriented notes
- Mark which candidates are high-confidence, moderate-confidence, or exploratory
- Keep notes on why each candidate is prioritized

## Subtasks
- Define presentation summary schema.
- Define shortlist cutoffs.
- Define compact rationale text fields.

## Test Tasks
- Confirm presentation tables are compact and readable.
- Confirm each shortlisted candidate still retains a rationale.
- Confirm the presentation layer does not erase uncertainty.

## Dependencies
Depends on Task Groups 4.5 and 4.6.

## Deliverables
- `phase4_presentation_shortlist.csv`
- `phase4_top_candidates_brief.md`

---

# Recommended Initial Execution Order for Phase 4

The recommended order is:

1. **Task Group 4.0** — Multi-phase evidence ingestion and validation  
2. **Task Group 4.1** — Multi-axis score framework definition  
3. **Task Group 4.2** — Mechanistic classification logic  
4. **Task Group 4.3** — Perturbation-relevance scoring logic  
5. **Task Group 4.4** — State-robustness and accessibility interpretation  
6. **Task Group 4.5** — Final review-first output design  
7. **Task Group 4.6** — Final report and interpretation guide  
8. **Task Group 4.7** — Presentation-ready summary layer

### Why this order
- Phase 4 must start by validating that all upstream evidence is internally consistent.
- A score framework must be defined before candidates can be ranked.
- Mechanistic class assignment must be available before final interpretation is meaningful.
- Ranking must be complete before final reports and presentation layers are generated.
- Final outputs should prioritize interpretability and scientific defensibility.

---

## Korean Summary

이 문서는 Phase 4를 구현 단위로 쪼갠 task 문서다. 핵심은 **좋은 포켓**을 고르는 것이 아니라, **MYO1D 부착 방해 가능성이 가장 높은 site를 mechanistic class와 함께 최종 우선순위화하는 것**이다. 주요 흐름은 다중 phase 증거 검증 → 4축 점수 체계 정의 → orthosteric/rim/allosteric/irrelevant 분류 → perturbation relevance ranking → state robustness 해석 → 최종 review table → 최종 보고서 → 발표용 요약 정리 순서다.

