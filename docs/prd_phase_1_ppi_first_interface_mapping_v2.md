## Context Summary
- Project: EGFR–MYO1D interface-centric perturbation discovery pipeline
- Current Subphase: Phase 1 PRD (v2 — Structural Input Upgrade)
- Purpose: Define the receptor-side MYO1D attachment patch before any pocket-level ligand prioritization
- Primary Engine: PyRosetta global docking
- Supporting Engine: LightDock (independent secondary validation)
- Auxiliary Engine: AlphaFold-Multimer (auxiliary structural evidence only, not core)
- Key Change in v2: Receptor and partner structural inputs are upgraded to address known fragment-docking limitations

---

# PRD — Phase 1: PPI-First Interface Mapping (v2)

## Goal
Identify and prioritize the **EGFR receptor-side MYO1D attachment patch** across the three receptor states, using structurally adequate inputs that minimize known fragment-docking artifacts.

## Why this phase comes first
The project's primary question is where MYO1D attaches to EGFR. Ligandable pockets only become scientifically meaningful after this receptor-side interface patch is defined.

## What changed from v1 to v2

### Problem statement
The v1 Phase 1 plan did not specify structural input requirements with enough precision. As a result, the existing pilot docking campaign was conducted with:
- **Receptor:** C-lobe fragment only (45 residues), missing the entire N-lobe (~235 residues)
- **Partner:** beta-meander truncated at residue 962, making VAL962 an N-terminal first residue

This produced usable pilot data but introduced two systematic artifacts:

1. **N-lobe absence artifact:** The C-lobe surface electrostatic landscape, solvent accessibility, and steric environment differ from the full kinase domain context. N-lobe steric occlusion was checked post-hoc via PyMOL superposition, but N-lobe presence during docking would change the sampled energy landscape itself, not merely filter results afterward.

2. **N-terminal truncation artifact:** VAL962 appeared as a 100% contact residue across all 5 valid structures. However, as the first residue in the chain, it has artificial N-terminal charge and excess backbone freedom. The project already recognized this issue and flagged it as requiring validation by extending the construct to ~955.

### v2 resolution
Phase 1 now explicitly requires:
- **Full kinase domain** (~280 residues, covering both N-lobe and C-lobe) as the receptor input
- **Extended beta-meander** (starting from ~residue 955 instead of 962) as the partner input
- **Orientation-aware filtering** as a mandatory filtering step, not an optional future improvement
- **Pilot data preservation** — existing C-lobe fragment results are retained as reference/pilot data for comparison, not discarded

---

## Inputs

### Receptor inputs (upgraded)
- Three receptor states, each as **full kinase domain**:
  - 3GT8_raw (full kinase domain)
  - 3GT8_cl38_48 (full kinase domain)
  - 3GT8_cl85_100 (full kinase domain)
- **Residue range caution:** The exact residue range depends on the numbering system used:
  - 3GT8 PDB numbering: ~681–990 (as deposited in PDB)
  - UniProt numbering (P00533): PDB + 24, yielding ~705–1014 (including C-terminal tail)
  - PPI fragment numbering (legacy): offset −16 from PDB
  - The active kinase core (excluding C-terminal tail) spans approximately UniProt 696–979
  - **All documents must explicitly declare which numbering system is in use.** Cross-referencing between systems must use the verified conversion formulas (see PyRosetta_PPI_Handoff numbering section).
- The full kinase domain includes both N-lobe and C-lobe, preserving:
  - N-lobe steric occlusion during docking (not just post-hoc filtering)
  - hinge region conformational context
  - activation loop influence on C-lobe surface accessibility
  - realistic electrostatic landscape across the entire receptor surface
- **Activation loop note:** 3GT8 has a disordered activation loop (residues ~831–852, UniProt numbering) that was previously restored using SwissModel for the Vina docking receptor. The same loop modeling approach or an equivalent must be applied for the full kinase domain preparation. An unresolved activation loop creates an artificial surface cavity that can attract docking poses to a non-physiological site.

### Partner inputs (upgraded)
- Primary partner: **extended beta-meander** (~residue 955–1006)
  - Extends at least 7 residues N-terminal to the current construct (962–1006)
  - Eliminates VAL962 N-terminal artifact
  - Preserves all 5 beta-sheets (8th–12th) with additional upstream context
- TH1 domain: plausibility envelope check only, not primary search input
- AlphaFold-Multimer: auxiliary structural evidence only, not core workflow

### Pilot data (historical reference only)
- Existing C-lobe fragment × beta-meander(962–1006) results exist as historical records.
- These are retained for optional methodological comparison but carry known systematic artifacts.
- Legacy site labels and residue assignments are historical reference material only and must not constrain interpretation of new results.

---

## Core Requirements

### Feature 1. Full kinase domain receptor preparation
For each of the three receptor states, prepare a full kinase domain structure suitable for PyRosetta global docking. This includes:
- N-lobe + C-lobe as a single chain
- Activation loop modeled or confirmed as resolved
- Chain ID and residue numbering explicitly documented
- Numbering system (PDB vs UniProt) declared and conversion recorded

### Feature 2. Extended beta-meander partner preparation
Prepare the extended beta-meander construct (~955–1006) with:
- At least 7 additional residues upstream of the current 962 start
- Structural integrity confirmed (no broken backbone, reasonable geometry)
- Sheet boundaries (8th–12th) explicitly annotated
- Active face (sheets 8, 9) and structural support face (sheet 12) documented

### Feature 3. Orientation-aware filtering (mandatory)
Every surviving PPI model must pass orientation-aware filtering before it is accepted as valid interface evidence. This means:
- Sheet 8/9 active face must be oriented toward the receptor surface, not flipped away
- A geometric or vector-based face orientation check must be implemented
- Poses where sheet 8/9 residues contact the receptor only via backbone or edge-on contacts must be flagged
- Face-flip detection must be explicit and logged, not implicit in other filters

### Feature 4. Receptor-side interface extraction
For every surviving, orientation-validated PPI model, extract receptor-side and partner-side interface residues separately.

### Feature 5. Cluster-level interface consensus
Cluster PPI models and compute receptor-side residue occupancy/frequency per cluster. Only orientation-validated models contribute to consensus.

### Feature 6. Multi-state interface comparison
Compare receptor-side interface patches across the three receptor states.

### Feature 7. Pilot data comparison layer (optional reference only)
If useful for methodological comparison, the new full-kinase-domain results may be compared against existing C-lobe fragment pilot data. However:
- The pilot data (C-lobe fragment × truncated beta-meander) was generated with a fundamentally different system and carries known systematic artifacts.
- New results should be interpreted on their own merit, without any expectation that they reproduce pilot site locations or rankings.
- The absence of a pilot site in the new data is not a failure; it may indicate the pilot site was an artifact of the fragment system.
- This comparison is optional and informational. It must not constrain the interpretation of new results.

### Feature 8. LightDock secondary validation
Use LightDock as an independent secondary validation axis. LightDock results should be compared with PyRosetta results to identify convergent receptor-side patches. LightDock is not a replacement for PyRosetta but provides method-independence evidence.

### Feature 9. TH1 plausibility envelope
After receptor-side patches are defined using the beta-meander, assess whether the top patches remain structurally plausible when the full TH1 domain context is considered. TH1 is not used as the primary docking input.

### Feature 10. Interface patch summary table
Create a standardized table that summarizes:
- receptor_id
- receptor_construct (full_kinase_domain)
- partner_construct (extended_beta_meander)
- cluster_id
- receptor_interface_residues
- receptor_residue_occupancy
- orientation_validation_status
- dG / dSASA / sc / packstat if available
- partner-side residue summary
- confidence notes

---

## User Story
As the researcher, I want to define the receptor-side MYO1D attachment patch using structurally adequate inputs (full kinase domain + extended beta-meander) with orientation-aware quality control, so that later ligand ranking is biologically anchored to a defensible PPI interface definition, free from known fragment-docking artifacts.

## Acceptance Criteria
- [ ] Full kinase domain structures are prepared and validated for all three receptor states.
- [ ] Extended beta-meander (~955–1006) is prepared and validated.
- [ ] Orientation-aware filtering is implemented and applied to all surviving models.
- [ ] Receptor-side interface residues are extracted for PyRosetta models that pass orientation filtering.
- [ ] Cluster-level receptor-side residue occupancy can be computed from orientation-validated models.
- [ ] The three receptor states can be compared at the receptor-side patch level.
- [ ] LightDock secondary validation is available for convergence comparison.
- [ ] Pilot data comparison shows which prior conclusions are stable vs artifact-dependent. *(optional — new results stand on their own)*
- [ ] A standardized Phase 1 interface summary file is generated with all required fields.
- [ ] Old site labels are not treated as fixed truth.
- [ ] The extended beta-meander construct eliminates N-terminal truncation artifacts. Any residue's significance is assessed fresh from the new data, not by comparison to pilot results.

## Primary Outputs
- `ppi_interface_patch_table.csv`
- `ppi_cluster_summary.csv`
- `ppi_hotspot_residues.csv`
- `orientation_filter_log.csv`
- `lightdock_interface_support_table.csv`
- `afm_interface_support_table.csv` (if available)
- `pilot_comparison_table.csv` (optional, only if methodological comparison performed)

## Non-Goals for Phase 1
- Ligand docking prioritization
- Final perturbation ranking
- Full integrated report generation
- MD simulation (deferred to post-Phase-1 validation, but recognized as needed before Phase 2 commitment)

## Open Questions (updated from v1)
- ~~Which MYO1D construct should be treated as the primary docking partner?~~ → **Resolved: extended beta-meander (~955–1006)**
- What minimum cluster occupancy threshold should define a receptor-side hotspot?
- What geometric criterion defines an acceptable sheet 8/9 face orientation?
- How should the orientation filter handle ambiguous edge cases?
- Should a minimum dSASA threshold be applied to full-kinase-domain models differently than fragment models?
- What C-terminal tail handling policy should be adopted for the full kinase domain receptor? (include tail, exclude tail, or test both)

## Deferred but recognized needs
- **MD validation (100–200 ns)** of top cluster representatives from the full-kinase-domain docking is recognized as scientifically necessary before committing to Phase 2. It is not part of Phase 1's core implementation but should be planned as a Phase 1→2 gate.

---

## Korean Summary (간단 요약)

이 Phase의 목표는 ligand보다 먼저 **MYO1D가 EGFR 어디에 붙는지 receptor-side patch를 정의하는 것**이다.

v2의 핵심 변경:
- **Receptor를 C-lobe fragment(45 res)에서 전장 kinase domain(~280 res)으로 업그레이드**한다. N-lobe가 도킹 과정에서 steric occlusion과 electrostatic landscape에 직접 영향을 미치기 때문이다.
- **Beta-meander를 ~955부터 시작하도록 N-terminal 확장**한다. VAL962 말단 artifact를 제거하기 위함이다.
- **Orientation-aware filtering을 필수화**한다. Beta-meander의 얇은 β-sheet 구조는 face-flip이 쉽게 발생하므로, sheet 8/9의 active face가 receptor를 향하고 있는지 기하학적으로 검증해야 한다.
- 기존 C-lobe fragment 결과는 **역사적 참고용으로만 보존**한다. 새 결과의 해석을 제약하는 데 사용하지 않는다.
- LightDock를 독립 보조 검증축으로 사용한다.
- MD validation은 Phase 1 core에는 포함하지 않지만, Phase 2 진입 전 gate로 인식한다.

