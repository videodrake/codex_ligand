## Context Summary
- Project: EGFR–MYO1D interface-centric perturbation discovery pipeline
- Current Subphase: Phase 2 Task Breakdown
- Purpose: Convert the Phase 2 PRD into implementation-facing tasks for candidate pocket proposal and druggability mapping
- Upstream Dependency: Phase 1 must provide a machine-readable receptor-side patch reference
- Key Principle: Pocket discovery is not the goal by itself; candidate pockets must be classified by relevance to the MYO1D receptor-side attachment patch

---

# Task Breakdown
## Phase 2: Pocket Proposal and Druggability Mapping

This document breaks Phase 2 into implementation-facing task groups. The purpose of this phase is to build a receptor-state-specific catalog of candidate ligandable pockets and classify each candidate pocket by its relationship to the Phase 1 MYO1D receptor-side patch.

This phase is intentionally **not** the final ligand docking phase.
Its final responsibility is to produce a pocket catalog that is not only biologically annotated, but also ready for budget allocation in Phase 3.
It is a **candidate pocket enumeration and annotation phase**.

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

# Task Group 2.0: Phase 1 Patch Reference Ingestion and Validation
**Priority:** Must-Have  
**Zone:** 🟢 Green

## Objective
Load and validate the machine-readable receptor-side patch reference produced by Phase 1 so that pocket relevance can be anchored to a biologically defined interface patch.

## Main Tasks

### 2.0.1 Phase 1 patch reference ingestion
- Read the Phase 1 downstream patch reference file.
- Preserve receptor_id, patch_id, hotspot residues, patch centroid, robustness label, and confidence class.

### 2.0.2 Patch-reference validation
- Confirm that each receptor state has a corresponding patch reference.
- Confirm residue numbering and chain labels remain interpretable.
- Confirm that Phase 1 confidence labels are present.

### 2.0.3 Patch reference standardization
- Convert Phase 1 patch reference into a stable internal format that Phase 2 tools can consume without ambiguity.

## Subtasks
- Define patch reference schema for Phase 2 ingestion.
- Add validation rules for missing or malformed patch fields.
- Add compatibility checks between receptor metadata and patch metadata.

## Test Tasks
- Confirm the patch reference file can be loaded for all three receptor states.
- Confirm receptor IDs and patch IDs remain traceable.
- Confirm malformed or incomplete patch references trigger warnings or errors.

## Dependencies
Depends on Phase 1 completion.

## Deliverables
- `phase2_patch_reference_validation.md`
- `phase2_patch_reference_normalized.csv`

---

# Task Group 2.1: Multi-Tool Candidate Pocket Proposal
**Priority:** Must-Have  
**Zone:** 🟢 Green

## Objective
Enumerate candidate ligandable pockets for each receptor state using one or more pocket proposal methods, while keeping the workflow open to multiple tools rather than hard-coding a single source of truth.

## Main Tasks

### 2.1.1 Pocket proposal source integration
Support candidate pocket generation from multiple possible sources, such as:
- fpocket
- P2Rank
- optionally other future hotspot/druggability tools

### 2.1.2 Receptor-local pocket proposal
- Generate candidate pockets independently for each receptor state.
- Preserve receptor identity in every candidate pocket record.

### 2.1.3 Pocket metadata extraction
For each candidate pocket, extract at least:
- receptor_id
- candidate_pocket_id
- centroid_x / centroid_y / centroid_z
- proposal source
- proposal score if available
- optional box size estimates

## Subtasks
- Define a tool-agnostic candidate pocket schema.
- Implement ingestion for pocket proposal outputs.
- Add source labels and raw proposal score fields.
- Preserve proposal provenance by tool.

## Test Tasks
- Confirm at least one pocket proposal source can generate pockets for each receptor state.
- Confirm candidate pocket records remain receptor-local.
- Confirm output records preserve source labels and coordinates.
- Confirm proposal outputs are machine-readable and reviewable.

## Dependencies
Depends on Task Group 2.0 and receptor structure availability.

## Deliverables
- `candidate_pockets_raw.csv`
- `candidate_pocket_source_summary.csv`
- `phase2_pocket_proposal_note.md`

---

# Task Group 2.2: Candidate Pocket Normalization and Merge Logic
**Priority:** Must-Have  
**Zone:** 🟢 Green

## Objective
Normalize and merge closely overlapping candidate pockets within each receptor state so that downstream docking budget is not wasted on duplicate pockets proposed by different tools.

## Main Tasks

### 2.2.1 Merge-threshold definition
- Define how close two candidate pockets must be to count as overlapping proposals.
- Preserve raw metrics rather than silently merging without trace.

### 2.2.2 Receptor-local merge operation
- Merge or group candidate pockets only within the same receptor state.
- Never merge across receptor states.

### 2.2.3 Pocket provenance preservation
- Keep track of which original tool-specific proposals contributed to a merged candidate pocket.

### 2.2.4 Stable candidate pocket IDs
- Generate stable, receptor-local candidate pocket IDs that downstream phases can reuse.

## Subtasks
- Define centroid-based merge metrics.
- Define optional residue-overlap merge logic if available later.
- Build source-contribution tracking.
- Create merged pocket ID naming convention.

## Test Tasks
- Confirm nearby proposals from different tools can be merged.
- Confirm distant proposals are not merged incorrectly.
- Confirm merged pockets retain provenance to original proposal sources.
- Confirm candidate pocket IDs remain stable and reproducible.

## Dependencies
Depends on Task Group 2.1.

## Deliverables
- `candidate_pocket_merge_table.csv`
- `candidate_pockets.csv`
- `candidate_pocket_provenance.csv`

---

# Task Group 2.3: Patch Relationship Classification
**Priority:** Must-Have  
**Zone:** 🟢 Green

## Objective
Classify each candidate pocket by its structural relationship to the Phase 1 MYO1D receptor-side patch.

## Main Tasks

### 2.3.1 Define relationship classes
Each candidate pocket must be assigned one of the following working classes:
- `orthosteric_candidate`
- `rim_candidate`
- `allosteric_candidate`
- `low_relevance_candidate`

### 2.3.2 Relationship metrics
Use one or more of the following to support classification:
- patch-to-pocket centroid distance
- hotspot residue overlap
- contact-residue overlap if available later
- patch coverage / edge proximity

### 2.3.3 Classification transparency
- Preserve raw metrics used for classification.
- Do not reduce classification to a hidden black-box label.

## Subtasks
- Define relationship thresholds.
- Define a raw metrics table.
- Define relationship label assignment logic.
- Add room for future refinement without breaking current outputs.

## Test Tasks
- Confirm every candidate pocket receives a relationship class.
- Confirm raw classification metrics are preserved.
- Confirm orthosteric and rim candidates can be distinguished where appropriate.
- Confirm the system does not merge all near-patch candidates into one label without evidence.

## Dependencies
Depends on Task Groups 2.0 and 2.2.

## Deliverables
- `pocket_patch_relationship.csv`
- `pocket_patch_relationship_metrics.csv`
- `phase2_relationship_classification_note.md`

---

# Task Group 2.4: Druggability Confidence Layer
**Priority:** Should-Have  
**Zone:** 🟡 Yellow

## Objective
Separate simple geometric pocket presence from stronger evidence of ligandability/druggability.

## Main Tasks

### 2.4.1 Proposal-score normalization
- Preserve raw proposal scores from each tool.
- Add a simple normalized confidence field or source confidence annotation if possible.

### 2.4.2 Multi-source support scoring
- Indicate whether a pocket is supported by more than one proposal source.
- Distinguish single-tool pockets from consensus pockets.

### 2.4.3 Optional hotspot/druggability support
- Allow future integration of hotspot-based or fragment-based support layers without restructuring the core file schema.
-Preferred hotspot support tool: FTMap should be treated as the primary hotspot-support method when available, because Phase 2 needs stronger evidence for PPI-site ligandability than geometric pocket proposal alone.

## Subtasks
- Define a minimal druggability support schema.
- Define a consensus-support flag.
- Keep tool-specific raw fields intact.

## Test Tasks
- Confirm single-tool and multi-tool candidate pockets can be distinguished.
- Confirm raw proposal confidence is preserved.
- Confirm future tool additions would not break the schema.

## Dependencies
Depends on Task Group 2.2.

## Deliverables
- `druggability_proposal_summary.csv`
- `candidate_pocket_support_flags.csv`

---

# Task Group 2.5: Cross-State Pocket Proposal Alignment
**Priority:** Should-Have  
**Zone:** 🟡 Yellow

## Objective
Compare candidate pocket proposals across receptor states before ligand docking, in order to identify pockets that are conserved, shifted, or state-specific.
This comparison is critical not only for robustness assessment, but also for distinguishing persistent pockets from state-dependent or potentially cryptic pockets.

## Main Tasks

### 2.5.1 Receptor-state pocket alignment
- Compare candidate pocket locations across receptor states.
- Preserve receptor-state identity explicitly.

### 2.5.2 Cross-state proposal categories
Classify candidate pockets as:
- `state_robust_pocket`
- `state_shifted_pocket`
- `state_specific_pocket`
- `uncertain_alignment`

### 2.5.3 Preserve raw metrics
Store cross-state comparison metrics rather than only final labels.

## Subtasks
- Define cross-state pocket comparison metrics.
- Define alignment labels.
- Preserve uncertainty where alignment is weak or numbering is problematic.

## Test Tasks
- Confirm at least one pocket proposal can be compared across receptor states.
- Confirm robust vs state-specific candidate pockets can be distinguished.
- Confirm raw comparison metrics remain visible.

## Dependencies
Depends on Task Groups 2.2 and 2.3.

## Deliverables
- `candidate_pocket_cross_state_comparison.csv`
- `candidate_pocket_state_classes.csv`

---

# Task Group 2.6: Phase 3 Docking Preparation Export
**Priority:** Must-Have  
**Zone:** 🟡 Yellow

## Objective
Export a clean candidate pocket catalog and metadata package that Phase 3 can use directly for diversity-aware pocket-guided docking.

## Main Tasks

### 2.6.1 Candidate pocket docking export
For each final candidate pocket, export at least:
- receptor_id
- candidate_pocket_id
- centroid
- optional docking box estimate
- patch relationship class
- pocket support annotation
- state class if available

### 2.6.2 Phase 3 compatibility checks
- Confirm that the exported pocket catalog can be consumed without additional manual cleanup.
- Confirm that receptor IDs and candidate pocket IDs are stable.

### 2.6.3 Downstream readiness note
- Summarize what Phase 3 should treat as the primary pocket list.
- Explicitly mark low-confidence candidates where relevant.

## Subtasks
- Define export schema for Phase 3.
- Add final field validation.
- Add machine-readable and human-readable export notes.

## Test Tasks
- Confirm the Phase 3 export file is complete and loadable.
- Confirm candidate pocket IDs remain stable.
- Confirm every exported pocket retains receptor identity and patch relationship class.

## Dependencies
Depends on Task Groups 2.3 and 2.4, and optionally 2.5.

## Deliverables
- `phase3_candidate_pocket_reference.csv`
- `phase2_to_phase3_handoff_note.md`

---

# Task Group 2.7: Phase 2 Review Report
**Priority:** Should-Have  
**Zone:** 🟡 Yellow

## Objective
Generate a readable Phase 2 review package that summarizes candidate pockets, patch relationship classes, and the resulting Phase 3 docking reference set.

## Main Tasks

### 2.7.1 Candidate pocket summary report
- Summarize candidate pockets per receptor state.
- Summarize proposal sources.
- Summarize merged pocket counts and classes.

### 2.7.2 Patch relationship summary
- Show which pockets are orthosteric, rim, allosteric, or low relevance.
- Preserve uncertainty notes where classification is weak.

### 2.7.3 Downstream handoff summary
- Explain which pockets should be used in Phase 3 docking.
- Explain which pockets are experimental, optional, or low confidence.

## Subtasks
- Define report sections.
- Add summary tables for receptors and pocket classes.
- Add a final candidate-pockets-for-docking section.

## Test Tasks
- Confirm the report can be read without opening raw proposal files.
- Confirm the report clearly distinguishes geometric pockets from biologically relevant pockets.
- Confirm the report is sufficient to hand off into Phase 3.

## Dependencies
Depends on Task Groups 2.3, 2.4, and 2.6.

## Deliverables
- `phase2_candidate_pocket_report.md`

---

# Recommended Initial Execution Order for Phase 2

The recommended order is:

1. **Task Group 2.0** — Phase 1 patch reference ingestion and validation  
2. **Task Group 2.1** — Multi-tool candidate pocket proposal  
3. **Task Group 2.2** — Candidate pocket normalization and merge logic  
4. **Task Group 2.3** — Patch relationship classification  
5. **Task Group 2.4** — Druggability confidence layer  
6. **Task Group 2.5** — Cross-state pocket proposal alignment  
7. **Task Group 2.6** — Phase 3 docking preparation export  
8. **Task Group 2.7** — Phase 2 review report

### Why this order
- Phase 2 must begin with the validated Phase 1 patch reference.
- Candidate pockets must be proposed before they can be merged or classified.
- Candidate pockets must be normalized before patch relationship classes become stable.
- Druggability confidence is useful but secondary to basic pocket identification and classification.
- Phase 3 should only begin after a clean candidate pocket reference exists.

---

## Korean Summary

이 문서는 Phase 2를 구현 단위로 쪼갠 task 문서다. 핵심은 **candidate pocket을 많이 찾는 것 자체가 목표가 아니라, Phase 1에서 정의한 MYO1D receptor-side patch와의 관계를 기준으로 포켓을 분류하는 것**이다. 주요 흐름은 Phase 1 patch reference 검증 → pocket proposal → pocket merge → orthosteric/rim/allosteric 분류 → druggability confidence → cross-state pocket 비교 → Phase 3용 pocket reference export 순서다.

