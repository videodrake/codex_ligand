# CODEX HANDOFF
## EGFR-MYO1D Pipeline

---

## 1. What this document is

This is the **master handoff document** for Codex and future technical contributors.

Its role is **not** to contain the full scientific or implementation detail of the project.
Instead, its purpose is to:

- explain what the project is,
- identify which documents are currently active,
- define the current scientific and technical direction,
- state the working rules for repository modification,
- and prevent Codex from following outdated Vina-first assumptions.

This file should remain **short, current, and authoritative**.
Detailed planning belongs in the phase-specific documents listed below.

---

## 2. Project in one sentence

This repository supports an **EGFR-MYO1D interface-centric perturbation discovery pipeline** whose goal is to first define the receptor-side MYO1D attachment patch on EGFR, and then identify ligandable pockets that may disrupt that attachment directly or indirectly.

---

## 3. Current active scientific framing

The project is **no longer organized around a Vina-first pocket search**.

The current active scientific framing is:

1. **Phase 1 ??PPI-first interface mapping**  
   Define the receptor-side MYO1D attachment patch.

2. **Phase 2 ??Pocket proposal and druggability mapping**  
   Enumerate candidate ligandable pockets and classify them relative to the Phase 1 patch.

3. **Phase 3 ??Diversity-aware pocket-guided docking**  
   Replace naive repeated giant-box blind docking with diversity-aware docking guided by candidate pockets.

4. **Phase 4 ??Perturbation relevance scoring**  
   Rank pockets and ligand-supported sites by their likelihood of disrupting MYO1D attachment.

This means the project should always be interpreted as:

> **PPI-first ??pocket proposal ??diversity-aware docking ??perturbation relevance**

not as:

> **blind docking first, then try to explain biology later**

---

## 4. Core scientific conclusions already adopted as project rules

These points should be treated as active project assumptions unless explicitly changed later.

### 4.1 AlphaFold-Multimer is not part of the Phase 1 core workflow
Project-specific experience suggests that AlphaFold-Multimer was not sufficiently reliable in this use case.
It may be kept as historical or optional reference, but it is **not part of the active core Phase 1 plan**.

### 4.2 Extended beta-meander is the primary MYO1D docking input
The full TH1 domain generated too much noise when used as the main blind-search docking input.
The **extended beta-meander (~residue 955??006)** is currently the preferred primary input for interface mapping.
The earlier truncated construct (962??006) produced a VAL962 N-terminal artifact and must be replaced.
The extension adds at least 7 upstream residues to eliminate terminal charge and backbone freedom artifacts.

### 4.3 TH1 is a downstream plausibility envelope, not the main search input
TH1 should be used later to assess whether top beta-meander-derived poses remain structurally plausible in the larger MYO1D domain context.
It should **not** be treated as the primary blind-search input in the current active workflow.

### 4.4 Sheet 8 and sheet 9 are treated as the primary active face
Current project experience and prior literature support the interpretation that sheets **8 and 9** are the main functional interface face for docking purposes.

### 4.5 Sheet 12 is treated as structural support, not the main direct-contact face
Sheet 12 may matter for PPI disruption when mutated, but the current working interpretation is that it is more likely a **structural support element** than the primary direct-contact face.
This interpretation is supported by project-specific MD simulation data, in which sheet 12 did not form direct contacts with the receptor surface. The functional essentiality of sheet 12 (Ko et al. alanine substitution) is therefore attributed to its role in stabilizing the beta-meander fold, which indirectly enables sheets 8/9 to bind.

### 4.6 Orientation-aware filtering is mandatory, not optional
A pose is not acceptable simply because sheet 8/9 residues contact the receptor.
The beta-meander has a thin geometry and often produces **face-flipped poses**.
Therefore, the active workflow must be **orientation-aware**, not merely contact-aware.
Orientation filtering must be implemented as a mandatory pass/fail gate before any model enters consensus building.
This means computing the face direction of sheet 8/9 relative to the receptor surface and rejecting poses where the active face points away from the receptor.

### 4.7 New computational outputs outrank legacy site labels
Older residue/site labels from previous reports are reference material only.
They must never be hard-coded as truth if new structured outputs disagree with them.

### 4.8 Full kinase domain replaces C-lobe fragment for Phase 1 docking
The earlier C-lobe fragment docking (45 residues) produced useful pilot data but has known limitations:
- N-lobe absence distorts the C-lobe electrostatic landscape and steric environment
- N-lobe steric occlusion was only checked post-hoc, not enforced during docking
- The fragment context may create artificial surface pockets or expose buried residues
Therefore, Phase 1 core docking must use the full kinase domain (~280 residues).
C-lobe fragment results are preserved as pilot/reference data for comparison.

### 4.9 Pilot data is historical reference only, not a validation target
The C-lobe fragment docking results (C02, C04, C07 sites) were produced with a system that has known structural deficiencies.
New full-kinase-domain docking results must be interpreted on their own merit.
New results that differ from pilot data are expected improvements, not failures.
Pilot site names (C02, C04, C07) must not be used as expectations or validation targets for the new system.
Phase 2 should not begin until the new system has produced its own defensible patch definition ??independent of any pilot data.

---

## 5. Current receptor ensemble

The active receptor ensemble is fixed to these three states:

- `3GT8_raw`
- `3GT8_cl38_48`
- `3GT8_cl85_100`

These states must remain explicitly labeled and separated in all outputs.
Direct cross-state comparison requires careful residue-numbering and chain validation.

### 5.1 Receptor construct requirement (v2)
All three receptor states must be prepared as **full kinase domain**, including both N-lobe and C-lobe.
The exact residue range depends on the numbering system (see PRD v2 for detailed numbering caution). The active kinase core spans approximately UniProt 696??79 (excluding C-terminal tail); 3GT8 PDB numbering differs by an offset of +24.
The earlier C-lobe fragment (45 residues) is retained as historical reference only and must not be used for new Phase 1 core docking runs.
Full kinase domain preparation preserves N-lobe steric occlusion, hinge region context, and realistic electrostatic landscape during docking.

---

## 6. Current compute rule

The main server has **32 CPU cores**, but for routine use this project assumes only **16 cores are safely available**.

Important additional constraint:

- The current Codex workspace is **not the same as the real server environment**.
- Code may be developed and structurally validated here.
- Actual high-load parallel performance validation must be treated as a **server-side task**, not as something proven in the Codex workspace.

Therefore:
- `max_workers=16` should remain the intended server-side default,
- but no environment-specific performance claims should be made based only on the Codex workspace.

---

## 7. Repository modification rules for Codex

Codex should follow these rules unless the user explicitly overrides them.

### 7.1 Do not rewrite the whole repository
Refactor and extend the existing codebase whenever possible.
Prefer narrow, inspectable changes over broad rewrites.

### 7.2 Do not use outdated Vina-first documents as active design truth
Older broad-scope Vina-first documents should be treated as legacy reference only.
The active plan is the new 4-phase structure.

### 7.3 Keep raw evidence visible
Do not hide logic behind opaque conclusions.
Outputs should remain reviewable and traceable.

### 7.4 Separate primary and secondary evidence
- Primary evidence in Phase 1: PyRosetta + orientation-aware filtering
- Secondary independent support in Phase 1: LightDock
- TH1: plausibility envelope only

### 7.5 Avoid hard-coding biological assumptions into code unless explicitly approved
For example:
- old site names,
- old residue labels,
- legacy patch IDs,
- or unvalidated receptor mapping shortcuts
should not be silently embedded in code.

### 7.6 Preserve machine-readable intermediate outputs
Every phase should produce reusable tables, not just terminal logs or final narrative summaries.

---

## 8. Current active documents (read these first)

Codex should read these active documents in the following order.

### 8.1 Entry documents
1. `README.md`
2. `docs/project_context.md`
3. `docs/runbook.md`

### 8.2 Active planning documents
4. `docs/brief_egfr_myo_1_d_pipeline_v_2.md`
5. `docs/prd_phase_1_ppi_first_interface_mapping_v2.md`
6. `docs/tasks_phase_1_ppi_first_interface_mapping_v2.md`
7. `docs/prd_phase_2_pocket_proposal_and_druggability_mapping.md`
8. `docs/prd_phase_3_diversity_aware_pocket_guided_docking.md`
9. `docs/prd_phase_4_perturbation_relevance_scoring.md`

If phase-specific task files for Phases 2?? exist later, those should be read after the corresponding phase PRDs.

---

## 9. Current active implementation priority

The current implementation priority is:

### Priority 0 (prerequisite for all else)
**Structural input upgrade**
- Full kinase domain receptor preparation (3 states)
- Extended beta-meander (~955??006) preparation
- Input validation and metadata

### Priority 1
**Phase 1 refinement and stabilization (with upgraded inputs)**
- extended beta-meander input definition
- active face / forbidden face definition
- **orientation-aware filtering implementation (mandatory)**
- face-flip filtering
- PyRosetta global docking with full kinase domain
- receptor-side interface consensus (orientation-validated only)
- LightDock-based secondary validation
- TH1 plausibility evaluation
- pilot data comparison (C-lobe fragment vs full kinase domain)
- VAL962 artifact assessment
- Phase 2 patch reference export

### Priority 2
**Phase 2 pocket proposal**
- candidate pocket catalog
- patch relationship classification
- druggability mapping

### Priority 3
**Phase 3 diversity-aware docking**
- candidate-pocket-guided docking
- saturation rule
- budget reallocation

### Priority 4
**Phase 4 perturbation relevance scoring**
- orthosteric / rim / allosteric / irrelevant classification
- multi-axis ranking

---

## 10. What Phase 1 must achieve before Phase 2 can begin

Phase 2 should not begin until Phase 1 can export a defensible receptor-side patch reference.

At minimum, Phase 1 must produce:

- a receptor-side attachment patch candidate set **derived from full-kinase-domain docking (not C-lobe fragment)**,
- cross-state patch robustness information,
- a confidence classification for each patch,
- **orientation-filtered evidence only** (no face-flipped poses in consensus),
- and a machine-readable Phase 2 patch reference file.

All of the above must be derived entirely from the new full-kinase-domain + extended-beta-meander system. Pilot data (C-lobe fragment results) is not part of the Phase 2 gate criteria.

Additionally, MD validation (100??00 ns) of the top 1?? cluster representatives is strongly recommended as a Phase 1?? gate, though it is not included in the Phase 1 core task list. If MD validation reveals that the top patch is dynamically unstable, Phase 2 pocket proposal would be premature.

If those are not available, pocket relevance classification in later phases will be biologically underconstrained.

---

## 11. What Codex should never assume

Codex must **not** assume any of the following without explicit confirmation:

- that the old Vina-first project documents are still active,
- that TH1 should be the primary search input,
- that AlphaFold-Multimer is still part of the core workflow,
- that sheet 12 is the primary direct-contact face,
- that sheet 8/9 contact alone is sufficient to accept a pose (orientation filtering is mandatory),
- **that the C-lobe fragment (45 res) is still the active receptor input (full kinase domain is now required),**
- **that the truncated beta-meander (962??006) is still the active partner input (extended ~955??006 is now required),**
- **that VAL962 is a confirmed anchor residue (it is under artifact assessment),**
- that current workspace performance reflects real server performance,
- or that old residue/site labels are automatically correct.

---

## 12. Current expected document split

This repository should now be understood as having two document layers:

### A. Master/index documents
These explain the repository and how to work in it.
- `README.md`
- `docs/project_context.md`
- `docs/runbook.md`
- `docs/codex_handoff_egfr_myo_1_d_pipeline_v2.md`

### B. Active scientific planning documents
These contain the real project logic.
- `docs/brief_egfr_myo_1_d_pipeline_v_2.md`
- phase-specific PRDs
- phase-specific task files

This handoff file belongs to layer A.
It should remain short and should not become another giant technical specification.

---

## 13. If Codex is uncertain what to do next

Codex should follow this rule:

1. First determine which active phase is being worked on.
2. Then read the corresponding phase PRD and task document.
3. Then make the smallest safe change needed for that phase.
4. Then report changes, assumptions, and next blockers clearly.

If there is ambiguity, Codex should prefer:
- preserving traceability,
- keeping outputs structured,
- and avoiding premature architectural expansion.

---

## 14. Korean summary (간단 ?-약)

??문서??Codex??**최종 마스??handoff 문서**??

??심 ??용?? ??음???같다.

- ????로??트????제 **Vina-first**가 ??니??**PPI-first / 4-phase 구조**??
- 최신 active ??름??:
  1. MYO1D receptor-side patch 규명
  2. candidate pocket ??안
  3. diversity-aware docking
  4. perturbation relevance ranking
- AlphaFold-Multimer??core??서 ??외??다.
- beta-meander???primary input??로 ??용??다.
- TH1?? main blind input????니??plausibility check??이??
- sheet 8/9??primary active face, sheet 12??support face???본다.
- orientation-aware docking???face-flip filtering????수??
- LightDock????...립 보조 검증축??로 ??용??다.
- ??계산 결과가 기존 보고????벨보다 ??선??다.
- ??재 Codex workspace????제 ??버?? ??르므??? ??능 검증?? ??버??서 ??야 ??다.

??문서??**짧?? ??덱??문서**???? ????고, ??제 ??용?? phase???PRD/task 문서?????는 구조?????용??다.



