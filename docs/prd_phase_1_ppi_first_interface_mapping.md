## Context Summary
- Project: EGFR–MYO1D interface-centric perturbation discovery pipeline
- Current Subphase: Phase 1 PRD
- Purpose: Define the receptor-side MYO1D attachment patch before any pocket-level ligand prioritization
- Primary Engine: PyRosetta global docking
- Supporting Engine: AlphaFold-Multimer (auxiliary evidence only)

---

# PRD — Phase 1: PPI-First Interface Mapping

## Goal
Identify and prioritize the **EGFR receptor-side MYO1D attachment patch** across the three receptor states.

## Why this phase comes first
The project’s primary question is where MYO1D attaches to EGFR. Ligandable pockets only become scientifically meaningful after this receptor-side interface patch is defined.

## Inputs
- Three receptor states:
  - 3GT8_raw
  - 3GT8_cl38_48
  - 3GT8_cl85_100
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

이 Phase의 목표는 ligand보다 먼저 **MYO1D가 EGFR 어디에 붙는지 receptor-side patch를 정의하는 것**이다.
주 엔진은 PyRosetta이고, AFM은 보조 증거다. 산출물은 receptor-side interface residue와 cluster-level hotspot 요약이다.

