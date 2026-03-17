## Context Summary
- Project: EGFR–MYO1D interface-centric perturbation discovery pipeline
- Current Subphase: Phase 3 Task Breakdown
- Purpose: Convert the Phase 3 PRD into implementation-facing tasks for diversity-aware, pocket-guided ligand docking
- Upstream Dependency: Phase 2 must provide a budget-ready candidate pocket reference with machine-readable patch-relationship classification, druggability annotations, and receptor-local priority fields
- Key Principle: Docking should maximize biologically useful pocket diversity rather than repeatedly over-sampling the same dominant pocket

---

# Task Breakdown
## Phase 3: Diversity-Aware Pocket-Guided Docking

This document breaks Phase 3 into implementation-facing task groups. The purpose of this phase is **not** to run naive repeated giant-box blind docking. The purpose is to run ligand docking in a way that deliberately spreads search effort across biologically relevant candidate pockets, while preserving traceability and structured outputs.

This phase is governed by the following project rules:

- Phase 3 must consume the **Phase 2 budget-ready candidate pocket reference** rather than rediscovering pockets from scratch.
- The goal is not simply better affinity ranking.
- The goal is to allocate docking effort efficiently across candidate pockets that may matter for MYO1D perturbation.
- Dominant pockets must not consume unlimited search budget once they are sufficiently sampled.
- Search-budget decisions must remain transparent and auditable.

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

# Task Group 3.0: Phase 2 Candidate Pocket Reference Ingestion
**Priority:** Must-Have  
**Zone:** 🟢 Green

## Objective
Load and validate the structured Phase 2 candidate pocket reference, including patch-relationship classification, druggability annotations, and receptor-local priority fields, so that all docking jobs are anchored to explicitly prioritized receptor-local candidate pockets.

## Main Tasks

### 3.0.1 Candidate pocket reference ingestion
- Read the Phase 2 export file containing the structured candidate pocket catalog with patch-relationship classification, druggability annotations, and receptor-local priority fields.
- Preserve receptor_id, candidate_pocket_id, centroid, patch relationship class, druggability support, state class, phase3_priority_tier, and recommended budget fields.

### 3.0.2 Reference validation
- Confirm that candidate pockets are available for each receptor state.
- Confirm that required fields for docking budget allocation are present.
- Confirm that receptor IDs remain aligned with receptor metadata.

### 3.0.3 Internal normalization
- Convert the Phase 2 pocket reference into a stable internal representation that can drive receptor-local docking loops.

## Subtasks
- Define Phase 3 pocket-reference schema.
- Add validation logic for missing candidate pocket fields.
- Add warnings for malformed box definitions or incomplete priority annotations.

## Test Tasks
- Confirm the Phase 2 pocket reference can be loaded for all receptor states.
- Confirm candidate pocket IDs and receptor IDs remain traceable.
- Confirm malformed or incomplete pocket references trigger warnings or errors.

## Dependencies
Depends on Phase 2 completion.

## Deliverables
- `phase3_candidate_reference_validation.md`
- `phase3_candidate_reference_normalized.csv`

---

# Task Group 3.1: Pocket-Guided Docking Job Construction
**Priority:** Must-Have  
**Zone:** 🟢 Green

## Objective
Build receptor- and pocket-specific docking jobs from the Phase 2 candidate pocket reference.

## Main Tasks

### 3.1.1 Receptor-local job generation
- For each receptor state, generate docking jobs per candidate pocket.
- Preserve receptor-local separation.
- Never mix pockets across receptor states.

### 3.1.2 Pocket-local box generation
- Define or ingest docking box parameters for each candidate pocket.
- Preserve box metadata per candidate pocket.
- Allow later refinement without changing pocket identity.

### 3.1.3 Ligand dispatch matrix construction
- Combine ligand list with receptor-local candidate pockets.
- Support pocket-level job creation rather than only receptor-level blind jobs.

## Subtasks
- Define a docking job schema.
- Define receptor/pocket/ligand job naming rules.
- Add support for optional per-pocket box size overrides.
- Preserve output path conventions for later parsing.

## Test Tasks
- Confirm jobs can be generated for all receptor-pocket-ligand combinations.
- Confirm receptor-local and pocket-local identities are preserved.
- Confirm job naming is stable and deterministic.
- Confirm no cross-receptor pocket mixing occurs.

## Dependencies
Depends on Task Group 3.0.

## Deliverables
- `phase3_docking_job_table.csv`
- `phase3_job_box_table.csv`

---

# Task Group 3.2: Search Budget Policy and Saturation Rules
**Priority:** Must-Have  
**Zone:** 🟢 Green

## Objective
Define and implement an explicit search-budget policy so that already well-sampled pockets stop consuming disproportionate docking effort.

## Main Tasks

### 3.2.1 Budget parameter definition
Support fields such as:
- `recommended_seed_budget`
- `recommended_max_poses_per_pocket`
- `max_rounds`
- `pocket_cutoff`
- `saturation_affinity_window`

### 3.2.2 Saturation rule definition
A pocket should be marked as saturated when it has accumulated enough acceptable poses under the current round’s policy.
This must be configurable and recorded explicitly.

### 3.2.3 Budget reallocation policy
Once a pocket is saturated, remaining search effort should be redirected to unsaturated pockets according to priority tier and remaining budget.

### 3.2.4 Saturation transparency
The pipeline must preserve:
- why a pocket was marked saturated,
- how much budget it consumed,
- and what budget remained for other pockets.

## Subtasks
- Define minimum viable saturation rule.
- Define pocket status labels:
  - `open`
  - `saturated`
  - `skipped`
  - `exhausted`
- Define budget accounting fields.
- Define budget reallocation logic by priority tier.

## Test Tasks
- Confirm pockets can transition from open to saturated.
- Confirm saturated pockets stop receiving new search budget in later rounds.
- Confirm unsaturated pockets can receive reallocated budget.
- Confirm budget usage is logged in structured outputs.

## Dependencies
Depends on Task Group 3.1.

## Deliverables
- `phase3_budget_policy.md`
- `pocket_search_status.csv`
- `phase3_budget_tracking.csv`

---

# Task Group 3.3: Diversity-Aware Vina Execution Layer
**Priority:** Must-Have  
**Zone:** 🟢 Green

## Objective
Execute ligand docking in a way that follows the pocket-guided, budget-aware policy rather than naive repeated blind docking.

## Main Tasks

### 3.3.1 New execution entry point
Add a dedicated execution path or script for diversity-aware docking, such as:
- `run_diverse_docking.py`

This should coexist with legacy docking entry points rather than immediately replacing them.

### 3.3.2 Pocket-guided local docking execution
- Run docking per receptor and per candidate pocket.
- Support multiple seeds per pocket if configured.
- Preserve deterministic metadata for each job.

### 3.3.3 Controlled round-based execution
- Allow search to proceed in rounds.
- Reevaluate pocket status after each round.
- Stop allocating budget to saturated pockets.

### 3.3.4 Workspace vs server rule preservation
- Keep `max_workers=16` as the intended server-side default.
- Do not overinterpret current workspace performance.
- Treat this environment as functional-validation-only.

## Subtasks
- Define the diversity-aware docking runner interface.
- Add job-level metadata persistence.
- Add round-aware job dispatch.
- Preserve compatibility with current Vina execution utilities where possible.

## Test Tasks
- Confirm diversity-aware jobs can be launched per receptor and per pocket.
- Confirm round metadata is preserved.
- Confirm saturated pockets are not repeatedly resubmitted.
- Confirm functional behavior does not depend on current workspace performance assumptions.

## Dependencies
Depends on Task Group 3.2.

## Deliverables
- `run_diverse_docking.py`
- `phase3_run_metadata.json`
- `phase3_round_log.csv`

---

# Task Group 3.4: Pose Parsing Compatibility and Pocket Attribution
**Priority:** Must-Have  
**Zone:** 🟢 Green

## Objective
Ensure Phase 3 docking outputs remain compatible with the existing pose parsing flow while preserving pocket-guided provenance.

## Main Tasks

### 3.4.1 Pose provenance extension
Extend or preserve pose-level outputs so each pose can be traced to:
- receptor_id
- ligand_id
- candidate_pocket_id
- round_id
- seed
- docking_mode

### 3.4.2 Compatibility with existing parsing logic
- Reuse existing parse_vina_results.py and extract_contacts.py where possible.
- Avoid breaking current pose-table generation.

### 3.4.3 Pocket-attribution-aware pose table
Ensure pose-level outputs can later support:
- receptor-local clustering validation
- dominant pocket assignment
- perturbation scoring in Phase 4

## Subtasks
- Define added pose provenance fields.
- Update or wrap the current parsing flow if needed.
- Preserve backward-compatible outputs where reasonable.

## Test Tasks
- Confirm diversity-aware docking outputs can still be parsed.
- Confirm every pose retains candidate_pocket_id provenance.
- Confirm pose table remains machine-readable and consistent across runs.
- Confirm no information required for Phase 4 is lost.

## Dependencies
Depends on Task Group 3.3.

## Deliverables
- `vina_pose_table.csv` (extended or regenerated)
- `phase3_pose_provenance_note.md`

---

# Task Group 3.5: Pocket Occupancy and Diversity Validation
**Priority:** Must-Have  
**Zone:** 🟡 Yellow

## Objective
Validate that the new Phase 3 workflow actually improves site diversity and reduces pathological over-concentration into a small number of dominant pockets.

## Main Tasks

### 3.5.1 Pocket occupancy summary
Summarize, per receptor and ligand:
- how many poses landed in each candidate pocket
- how many pockets remained unsampled
- how many pockets became saturated

### 3.5.2 Diversity metrics
Define simple diversity-validation metrics such as:
- number of distinct candidate pockets sampled
- pose concentration ratio in the most dominant pockets
- fraction of total budget consumed by top N pockets

### 3.5.3 Comparison against naive blind baseline
Where possible, compare diversity-aware results against the older naive blind approach to show whether over-concentration was reduced.

## Subtasks
- Define occupancy summary schema.
- Define diversity-validation metrics.
- Add optional naive-vs-diverse comparison logic.
- Preserve summary outputs for later report use.

## Test Tasks
- Confirm pocket occupancy can be summarized for a completed Phase 3 run.
- Confirm diversity metrics can be computed from structured outputs.
- Confirm the system can show whether search became more distributed across candidate pockets.
- Confirm dominant pockets no longer consume unrestricted budget once saturated.

## Dependencies
Depends on Task Groups 3.3 and 3.4.

## Deliverables
- `phase3_pocket_occupancy_summary.csv`
- `phase3_diversity_metrics.csv`
- `phase3_blind_vs_diverse_comparison.csv` (optional if baseline available)

---

# Task Group 3.6: Phase 4-Ready Export
**Priority:** Must-Have  
**Zone:** 🟡 Yellow

## Objective
Export a Phase 4-ready docking evidence package that preserves ligand support, pocket provenance, and diversity-aware search history.

## Main Tasks

### 3.6.1 Final docking evidence export
For each candidate pocket and ligand-supported site, preserve:
- receptor_id
- candidate_pocket_id
- ligand_id
- pose support count
- best affinity
- mean affinity
- contact residues if available
- round and budget history summary

### 3.6.2 Perturbation-relevance handoff preparation
Ensure the export preserves enough information for Phase 4 to evaluate:
- orthosteric relevance
- rim relevance
- allosteric relevance
- ligand support strength
- state support context

### 3.6.3 Handoff quality checks
- Confirm the export is complete enough for Phase 4 scoring.
- Confirm no essential provenance fields are missing.

## Subtasks
- Define Phase 4 handoff schema.
- Add field completeness checks.
- Add a short downstream note explaining how to use the export.

## Test Tasks
- Confirm the export file is loadable by downstream tools.
- Confirm ligand support information is preserved.
- Confirm budget history is not lost before Phase 4.
- Confirm pocket provenance remains explicit.

## Dependencies
Depends on Task Groups 3.4 and 3.5.

## Deliverables
- `phase4_docking_evidence_reference.csv`
- `phase3_to_phase4_handoff_note.md`

---

# Task Group 3.7: Phase 3 Review Report
**Priority:** Should-Have  
**Zone:** 🟡 Yellow

## Objective
Generate a readable Phase 3 review package that shows how docking budget was distributed and whether candidate pockets were explored in a more balanced and biologically useful way.

## Main Tasks

### 3.7.1 Search-budget summary report
- Summarize how much budget each receptor-local pocket received.
- Show which pockets saturated and when.
- Show which pockets remained open or underexplored.

### 3.7.2 Ligand support summary
- Summarize ligand-supported pockets.
- Show multimodal vs dominant-pocket behavior.
- Preserve uncertainty where support is sparse.

### 3.7.3 Diversity outcome summary
- Show whether the new workflow reduced over-concentration into a few pockets.
- Preserve both raw metrics and readable interpretation.

## Subtasks
- Define report sections.
- Add budget summary tables.
- Add pocket occupancy and diversity summary sections.
- Add a final Phase 4 handoff summary section.

## Test Tasks
- Confirm the report can be read without opening raw docking files.
- Confirm budget behavior is understandable.
- Confirm ligand support summaries remain interpretable.
- Confirm the report is sufficient to hand off into Phase 4.

## Dependencies
Depends on Task Groups 3.5 and 3.6.

## Deliverables
- `phase3_diverse_docking_report.md`

---

# Recommended Initial Execution Order for Phase 3

The recommended order is:

1. **Task Group 3.0** — Phase 2 candidate pocket reference ingestion  
2. **Task Group 3.1** — Pocket-guided docking job construction  
3. **Task Group 3.2** — Search budget policy and saturation rules  
4. **Task Group 3.3** — Diversity-aware Vina execution layer  
5. **Task Group 3.4** — Pose parsing compatibility and pocket attribution  
6. **Task Group 3.5** — Pocket occupancy and diversity validation  
7. **Task Group 3.6** — Phase 4-ready export  
8. **Task Group 3.7** — Phase 3 review report

### Why this order
- Phase 3 must start from a validated Phase 2 candidate pocket reference.
- Jobs must be defined before budget policy can be applied.
- Budget policy must exist before diversity-aware execution becomes meaningful.
- Parsed output compatibility must be restored immediately after the new execution path is added.
- Diversity claims should only be made after structured occupancy validation exists.
- Phase 4 should only begin after a complete docking evidence reference is exported.

---

## Korean Summary

이 문서는 Phase 3를 구현 단위로 쪼갠 task 문서다. 핵심은 **naive giant-box blind docking 반복을 버리고, Phase 2에서 준비한 candidate pocket reference를 기준으로 search budget을 통제하면서 receptor-local pocket 다양성을 확보하는 것**이다. 주요 흐름은 pocket reference 검증 → pocket-guided job 생성 → saturation/budget policy → diversity-aware 실행 → pose provenance 유지 → pocket occupancy/diversity 검증 → Phase 4용 docking evidence export 순서다.

