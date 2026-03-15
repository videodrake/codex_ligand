> Status note (2026-03-12): This v2 PRD is the active Phase 1 PRD for the current baseline, but AFM should still
> be treated as legacy optional support only. The active Phase 1 secondary-validation path is LightDock.
> Use this PRD for Phase 1 target design details, and use `docs/current_pipeline_status.md` for current operational status.

## Context Summary
- Project: EGFR-MYO1D interface-centric perturbation discovery pipeline
- Current Subphase: Phase 1 PRD (v2 - Structural Input Upgrade)
- Purpose: Define the receptor-side MYO1D attachment patch before any pocket-level ligand prioritization
- Primary Engine: PyRosetta global docking
- Supporting Engine: LightDock (independent secondary validation)
- Auxiliary Engine: AlphaFold-Multimer (auxiliary structural evidence only, not core)
- Key Change in v2: Receptor and partner structural inputs are upgraded to address known fragment-docking limitations

---

# PRD ??Phase 1: PPI-first Interface Mapping (v2)

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
- **Pilot data preservation** ??existing C-lobe fragment results are retained as reference/pilot data for comparison, not discarded

---

## Inputs

### Receptor inputs (upgraded)
- Three receptor states, each as **full kinase domain**:
  - 3GT8_raw (full kinase domain)
  - 3GT8_cl38_48 (full kinase domain)
  - 3GT8_cl85_100 (full kinase domain)
- **Residue range caution:** The exact residue range depends on the numbering system used:
  - 3GT8 PDB numbering: ~681??90 (as deposited in PDB)
  - UniProt numbering (P00533): PDB + 24, yielding ~705??014 (including C-terminal tail)
  - PPI fragment numbering (legacy): offset ??6 from PDB
  - The active kinase core (excluding C-terminal tail) spans approximately UniProt 696??79
  - **All documents must explicitly declare which numbering system is in use.** Cross-referencing between systems must use the verified conversion formulas (see PyRosetta_PPI_Handoff numbering section).
- The full kinase domain includes both N-lobe and C-lobe, preserving:
  - N-lobe steric occlusion during docking (not just post-hoc filtering)
  - hinge region conformational context
  - activation loop influence on C-lobe surface accessibility
  - realistic electrostatic landscape across the entire receptor surface
- **Activation loop note:** 3GT8 has a disordered activation loop (residues ~831??52, UniProt numbering) that was previously restored using SwissModel for the Vina docking receptor. The same loop modeling approach or an equivalent must be applied for the full kinase domain preparation. An unresolved activation loop creates an artificial surface cavity that can attract docking poses to a non-physiological site.

### Partner inputs (upgraded)
- Primary partner: **extended beta-meander** (~residue 955??006)
  - Extends at least 7 residues N-terminal to the current construct (962??006)
  - Eliminates VAL962 N-terminal artifact
  - Preserves all 5 beta-sheets (8th??2th) with additional upstream context
- TH1 domain: plausibility envelope check only, not primary search input
- AlphaFold-Multimer: auxiliary structural evidence only, not core workflow

### Pilot data (historical reference only)
- Existing C-lobe fragment × beta-meander(962??006) results exist as historical records.
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
Prepare the extended beta-meander construct (~955??006) with:
- At least 7 additional residues upstream of the current 962 start
- Structural integrity confirmed (no broken backbone, reasonable geometry)
- Sheet boundaries (8th??2th) explicitly annotated
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
- [ ] Extended beta-meander (~955??006) is prepared and validated.
- [ ] Orientation-aware filtering is implemented and applied to all surviving models.
- [ ] Receptor-side interface residues are extracted for PyRosetta models that pass orientation filtering.
- [ ] Cluster-level receptor-side residue occupancy can be computed from orientation-validated models.
- [ ] The three receptor states can be compared at the receptor-side patch level.
- [ ] LightDock secondary validation is available for convergence comparison.
- [ ] Pilot data comparison shows which prior conclusions are stable vs artifact-dependent. *(optional ??new results stand on their own)*
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
- ~~Which MYO1D construct should be treated as the primary docking partner?~~ ??**Resolved: extended beta-meander (~955??006)**
- What minimum cluster occupancy threshold should define a receptor-side hotspot?
- What geometric criterion defines an acceptable sheet 8/9 face orientation?
- How should the orientation filter handle ambiguous edge cases?
- Should a minimum dSASA threshold be applied to full-kinase-domain models differently than fragment models?
- What C-terminal tail handling policy should be adopted for the full kinase domain receptor? (include tail, exclude tail, or test both)

## Deferred but recognized needs
- **MD validation (100??00 ns)** of top cluster representatives from the full-kinase-domain docking is recognized as scientifically necessary before committing to Phase 2. It is not part of Phase 1's core implementation but should be planned as a Phase 1?? gate.

---

## Korean Summary (간단 ?-약)

??Phase??목표??ligand보다 먼?? **MYO1D가 EGFR ??디??붙는지 receptor-side patch????*의??는 ???*??다.

v2????심 변???
- **Receptor???C-lobe fragment(45 res)??서 ??장 kinase domain(~280 res)??로 ?...그??이??*??다. N-lobe가 ??킹 과정??서 steric occlusion???electrostatic landscape??직접 ??향??미치?????문??다.
- **Beta-meander???~955부????작??도???N-terminal ?*장**??다. VAL962 말단 artifact?????거??기 ??함??다.
- **Orientation-aware filtering????수??*??다. Beta-meander?????? β-sheet 구조??face-flip????게 발생??????? sheet 8/9??active face가 receptor?????하?????는지 기하??적??로 검증해????다.
- 기존 C-lobe fragment 결과??**??????참고??으로만 보존**??다. ??결과????석????약??는 ????용???? ??는??
- LightDock????...립 보조 검증축??로 ??용??다.
- MD validation?? Phase 1 core??는 ??함???? ??????? Phase 2 진입 ??gate?????식??다.



