> Historical document: superseded by the current LightDock-based Phase 1 baseline.
> Read `docs/current_pipeline_status.md` and `docs/prd_phase_1_ppi_first_interface_mapping_v2.md` first.
> Do not use this file as the default planning baseline for new work.

## Context Summary
- Project: EGFR-MYO1D interface-centric perturbation discovery pipeline
- Current Subphase: Phase 1 PRD
- Purpose: Define the receptor-side MYO1D attachment patch before any pocket-level ligand prioritization
- Primary Engine: PyRosetta global docking
- Supporting Engine: AlphaFold-Multimer (auxiliary evidence only)

---

# PRD ??Phase 1: PPI-first Interface Mapping

## Goal
Identify and prioritize the **EGFR receptor-side MYO1D attachment patch** across the three receptor states.

## Why this phase comes first
The project??s primary question is where MYO1D attaches to EGFR. Ligandable pockets only become scientifically meaningful after this receptor-side interface patch is defined.

## Inputs
- Three receptor states:
  - 3GT8_raw
  - EGFR_160-185
  - EGFR_170-200
- MYO1D partner structures or fragments:
  - TH1 domain and/or beta-meander constructs
- Existing PyRosetta global docking pipeline
- Existing or future AlphaFold-Multimer runs

## Core Requirements

### Feature 1. Receptor-side interface extraction
For every surviving PPI model, extract receptor-side and partner-side interface residues separately.

### Feature 2. Cluster-level interface consensus
Cluster PPI models and compute receptor-side residue occupancy/frequency per cluster.

### Feature 3. Multi-state interface comparison
Compare receptor-side interface patches across the three receptor states.

### Feature 4. Auxiliary AFM support
Use AlphaFold-Multimer only as supporting evidence.
AFM outputs must be stored, scored, and compared, but must not replace the primary PyRosetta evidence layer.

### Feature 5. Interface patch summary table
Create a standardized table that summarizes:
- receptor_id
- cluster_id
- receptor_interface_residues
- receptor_residue_occupancy
- dG / dSASA / sc / packstat if available
- partner-side residue summary
- confidence notes

## User Story
As the researcher, I want to define the receptor-side MYO1D attachment patch before looking at ligand pockets, so that later ligand ranking is biologically anchored to the real PPI question.

## Acceptance Criteria
- [ ] Receptor-side interface residues are extracted for PyRosetta models.
- [ ] Cluster-level receptor-side residue occupancy can be computed.
- [ ] The three receptor states can be compared at the receptor-side patch level.
- [ ] AFM outputs can be included as auxiliary patch evidence.
- [ ] A standardized Phase 1 interface summary file is generated.
- [ ] Old site labels are not treated as fixed truth.

## Primary Outputs
- `ppi_interface_patch_table.csv`
- `ppi_cluster_summary.csv`
- `ppi_hotspot_residues.csv`
- `afm_interface_support_table.csv`

## Non-Goals for Phase 1
- Ligand docking prioritization
- Final perturbation ranking
- Full integrated report generation

## Open Questions
- Which MYO1D construct should be treated as the primary docking partner?
- What minimum cluster occupancy threshold should define a receptor-side hotspot?
- How should AFM low-confidence but spatially convergent results be weighted?

---

## Korean Summary

??Phase??목표??ligand보다 먼?? **MYO1D가 EGFR ??디??붙는지 receptor-side patch????*의??는 ???*??다.
????-진?? PyRosetta??고, AFM?? 보조 증거?? ??출물?? receptor-side interface residue?? cluster-level hotspot ?-약??다.


