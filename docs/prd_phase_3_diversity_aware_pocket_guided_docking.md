## Context Summary
- Project: EGFR–MYO1D interface-centric perturbation discovery pipeline
- Current Subphase: Phase 3 PRD
- Purpose: Replace naive repeated giant-box blind docking with diversity-aware pocket-guided docking
- Primary Engine: AutoDock Vina
- Key Control Concept: Search budget management with pocket saturation

---

# PRD — Phase 3: Diversity-Aware Pocket-Guided Docking

## Goal
Run ligand docking in a way that maximizes **scientifically useful pocket diversity** rather than over-allocating search effort to already dominant pockets.

## Why this phase exists
Repeated giant-box blind docking tends to over-concentrate results into strong pockets. That is useful for affinity discovery, but not sufficient for this project, which needs broad exploration of candidate pockets relevant to MYO1D perturbation.

Therefore, this phase must shift docking from:
- naive repeated blind docking

to
- candidate-pocket-guided docking with explicit search-budget control.

## Inputs
- Candidate pocket catalog from Phase 2
- Receptor metadata
- Ligand metadata
- Existing Vina execution machinery

## Core Requirements

### Feature 1. Candidate-pocket-driven docking
Dock ligands per receptor and per candidate pocket rather than relying only on one giant blind box.

### Feature 2. Configurable search budget
Support explicit budget parameters such as:
- seeds_per_pocket
- max_poses_per_pocket
- max_rounds
- pocket_cutoff
- saturation_affinity_window

### Feature 3. Saturation rule
A pocket that already accumulated enough acceptable poses should stop consuming the same search budget round.
This rule must be configurable and visible in outputs.

### Feature 4. Budget reallocation
Once pockets are saturated, remaining search effort should be redirected to unsaturated pockets.

### Feature 5. Structured diversity tracking
The pipeline must record:
- which pockets were searched
- which pockets saturated
- how much budget each pocket consumed
- which pocket each pose came from

## User Story
As the researcher, I want docking effort to spread across meaningful candidate pockets instead of repeatedly collapsing into 1–2 dominant pockets, so that I can evaluate a broader and more relevant set of perturbation candidates.

## Acceptance Criteria
- [ ] Docking can be run per receptor and per candidate pocket.
- [ ] Search budget parameters are configurable.
- [ ] Pocket saturation is tracked in a structured way.
- [ ] Remaining search effort can be redirected to unsaturated pockets.
- [ ] Pose outputs remain compatible with the existing pose table structure.
- [ ] Receptor-local diversity is preserved without mixing pockets across receptors.
- [ ] Behavior is deterministic enough to support repeated analysis with the same settings.

## Primary Outputs
- `vina_pose_table.csv` (extended if needed)
- `vina_pocket_table.csv`
- `vina_drug_pocket_map.csv`
- `pocket_search_status.csv`
- `candidate_pocket_run_log.csv`

## Non-Goals for Phase 3
- Final perturbation ranking
- Cross-method final verdict
- Full report generation

## Open Questions
- Should the MVP still allow one fallback giant-box mode for exploratory use?
- What exact saturation rule should be considered “sufficient” in the first release?
- How should multimodal ligands be handled in budget allocation?

---

## Korean Summary

이 Phase의 목표는 Vina를 단순 blind docking 반복이 아니라 **pocket-guided + saturation-controlled docking**으로 바꾸는 것이다. 즉, 한 포켓에 포즈가 충분히 쌓이면 그 포켓에는 더 이상 예산을 쓰지 않고, 다른 unsaturated pocket으로 계산 예산을 돌리게 만드는 단계다.

