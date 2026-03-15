> Status note (2026-03-12): This brief remains part of the active source-of-truth set. `docs/current_pipeline_status.md` is a derived summary only and must not override it.
> Current active Phase 1 baseline is PyRosetta + LightDock.
> AlphaFold-Multimer references in older planning language are historical or optional only.

## Context Summary
- Project: EGFR-MYO1D interface-centric perturbation discovery pipeline
- Current Phase: Phase 0 (Refactored Project Brief)
- Status: Refactored brief drafted to replace the earlier Vina-first framing
- Key Shift: Pipeline logic is reorganized from Vina-first to PPI-first -> pocket proposal -> diversity-aware docking -> perturbation scoring
- Intended Next Step: Create phase-specific PRDs and then phase-specific task files after review

---

# Project Brief (Refactored)

## One-Line Description
A research pipeline that first maps the EGFR receptor-side MYO1D attachment interface, then identifies ligandable pockets that can directly or indirectly disrupt that attachment across three EGFR receptor states.

## Core Scientific Goal
The goal is **not** to simply find good Vina pockets on EGFR C-lobe.
The goal is to answer two linked questions in the correct order:

1. **Where does MYO1D attach on the EGFR receptor-side surface?**
2. **Which ligandable sites can disrupt that attachment, either directly (orthosteric), at the interface rim, or allosterically?**

Therefore, the pipeline must prioritize:
- defining the receptor-side PPI interface patch first,
- mapping ligandable pockets second,
- and ranking pockets by **MYO1D-perturbation relevance**, not by affinity alone.

## Fixed Receptor States
The current receptor ensemble is explicitly fixed to:
1. **3GT8_raw**
2. **3GT8_cl38_48**
3. **3GT8_cl85_100**

These receptor states must remain comparable by receptor identity, residue numbering, and metadata.

## Refactored 4-Phase Scientific Architecture

### Phase 1 ??PPI-first Interface Mapping
Define the receptor-side MYO1D attachment patch using PyRosetta global docking as the primary engine, LightDock as the active secondary validation path, and AlphaFold-Multimer only as optional legacy auxiliary evidence.

### Phase 2 ??Pocket Proposal and Druggability Mapping
Enumerate ligandable pockets using structure-based pocket proposal tools and classify their spatial relationship to the Phase 1 PPI patch.

### Phase 3 ??Diversity-Aware Pocket-Guided Ligand Docking
Replace naive repeated giant-box blind docking with candidate-pocket-guided docking plus search-budget control and pocket saturation rules.

### Phase 4 ??Perturbation Relevance Scoring
Rank candidate pockets and ligands by their likelihood of disrupting MYO1D attachment, using orthosteric/rim/allosteric interpretation rather than affinity-only logic.

## Why This Refactor Is Necessary
The older Vina-first framing was useful for building a docking infrastructure, but it does not align perfectly with the real research question.

The actual question is not:
- ??Which pocket is strongest???

It is:
- ??Which ligandable site is most relevant to MYO1D attachment disruption???

That requires a PPI-first architecture.

## Explicitly Out of Scope
- Public web application features
- Generic docking platform behavior unrelated to EGFR-MYO1D
- Treating old residue/site labels as fixed truth
- Fully automated scientific conclusions without raw evidence
- Performance tuning based only on the current non-server Codex workspace

## Success Criteria
The refactored architecture is successful when:
- [ ] A receptor-side MYO1D interface patch can be defined and summarized across receptor states.
- [ ] Candidate ligandable pockets can be proposed independently of naive repeated giant-box docking.
- [ ] Ligand docking becomes diversity-aware rather than dominant-pocket-driven.
- [ ] Final outputs classify pockets by perturbation relevance to MYO1D attachment.
- [ ] All major outputs remain structured, reviewable, and compatible with future Codex-assisted development.

## Definition of Done for This Brief
This brief is complete when the project is clearly reframed as a 4-phase, PPI-first perturbation-discovery workflow and approved as the new basis for phase-specific PRDs.

---

## Korean Summary (간단 ?-약)

????로??트??????심?? ??음??다.
- 먼?? MYO1D가 EGFR ??디??붙는지 ?*의??다.
- ?????음 ???부착을 방해??????는 ligandable pocket??찾는??
- ??라??구조??**PPI-first ??pocket proposal ??diverse docking ??perturbation scoring**??로 바뀐다.
- 최종 목표????좋?? ??켓??이 ??니????MYO1D 부???방해 가????켓??을 찾는 것이??


