> Historical document: superseded by the current LightDock-based Phase 1 baseline.
> Read `docs/current_pipeline_status.md` and `docs/tasks_phase_1_ppi_first_interface_mapping_v2.md` first.
> Do not use this file as the default task baseline for new work.

## Context Summary
- Project: EGFR-MYO1D interface-centric perturbation discovery pipeline
- Current Subphase: Phase 1 Task Breakdown
- Purpose: Convert the Phase 1 PRD into implementation-facing tasks for receptor-side MYO1D interface mapping
- Primary Engine: PyRosetta global docking
- Supporting Engine: AlphaFold-Multimer (auxiliary evidence only)
- Key Principle: Define the receptor-side MYO1D attachment patch before ligand-site prioritization

---

# Task Breakdown
## Phase 1: PPI-first Interface Mapping

This document breaks Phase 1 into implementation-facing task groups. The purpose of this phase is to define the receptor-side MYO1D attachment patch in a reproducible, structured, and reviewable format before any ligand-pocket prioritization begins.

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

# Task Group 1.0: Receptor and Partner Input Validation
**Priority:** Must-Have  
**Zone:** ??? Green

## Objective
Validate that receptor structures and MYO1D partner inputs are complete, comparable, and suitable for PPI-first interface mapping.

## Main Tasks

### 1.0.1 Receptor ensemble validation
- Confirm the three receptor states are available and correctly labeled:
  - 3GT8_raw
  - EGFR_160-185
  - EGFR_170-200
- Confirm chain IDs are known
- Confirm residue numbering ranges are recorded
- Detect obvious numbering mismatches across receptor states

### 1.0.2 MYO1D partner input validation
- Identify the current MYO1D constructs used for docking
- Record whether TH1, beta-meander, or another fragment is treated as the primary docking partner
- Confirm input file availability and identity

### 1.0.3 Metadata persistence
- Build or update receptor metadata output
- Build or update partner metadata output
- Preserve source paths and structural notes

## Subtasks
- Define receptor metadata schema for Phase 1
- Define MYO1D partner metadata schema
- Add numbering/chain validation checks
- Add warnings for non-comparable residue numbering

## Test Tasks
- Confirm all three receptor states can be loaded and identified
- Confirm chain and residue metadata are emitted in a stable format
- Confirm partner input files are discovered and labeled consistently
- Confirm numbering mismatches surface as warnings or explicit validation output

## Dependencies
Depends on the repository/document baseline and existing receptor input availability.

## Deliverables
- `receptor_metadata.csv`
- `partner_metadata.csv`
- `phase1_input_validation_report.md`

---

# Task Group 1.1: PyRosetta Global Docking Standardization
**Priority:** Must-Have  
**Zone:** ??? Green

## Objective
Standardize PyRosetta global docking execution and output handling so receptor-side interface evidence can be extracted reproducibly.

## Main Tasks

### 1.1.1 Audit current PyRosetta execution path
- Identify current entry points
- Identify how decoys are generated
- Identify current score outputs
- Identify current cluster outputs, if any

### 1.1.2 Normalize execution inputs
- Ensure receptor ID and partner ID are propagated into run metadata
- Ensure run parameters are recorded
- Ensure random seed or reproducibility settings are stored when possible

### 1.1.3 Standardize raw output placement
- Separate outputs by receptor and partner input
- Preserve decoy, score, and clustering artifacts in stable locations

### 1.1.4 Standardize score extraction
- Extract or normalize score tables containing at least:
  - decoy_id
  - total score
  - interface-related metrics if available
  - receptor_id
  - partner_id

## Subtasks
- Inventory current PyRosetta score file conventions
- Define normalized output path convention
- Add or update run metadata file generation
- Implement a standardized PyRosetta score export if missing

## Test Tasks
- Confirm a PyRosetta run can be executed or replayed with stable metadata
- Confirm receptor_id and partner_id appear in downstream outputs
- Confirm score tables can be parsed without manual file reconstruction
- Confirm output placement does not mix receptor states

## Dependencies
Depends on Task Group 1.0.

## Deliverables
- `pyrosetta_run_metadata.json`
- `pyrosetta_decoy_scores.csv`
- `phase1_pyrosetta_execution_note.md`

---

# Task Group 1.2: Receptor-Side Interface Residue Extraction
**Priority:** Must-Have  
**Zone:** ??? Green

## Objective
Extract receptor-side and partner-side interface residues from PyRosetta models in a structured, reusable format.

## Main Tasks

### 1.2.1 Define interface extraction rule
- Choose and document the receptor?-partner contact rule
- Keep receptor-side and partner-side residues separate
- Keep residue string formatting consistent across outputs

### 1.2.2 Extract per-model interface residues
- For each decoy or selected model, extract:
  - receptor_interface_residues
  - partner_interface_residues
  - residue counts
  - optional per-residue interface metrics if available

### 1.2.3 Build receptor-side residue table
- Emit a long-form table that can be aggregated later by residue frequency and occupancy

## Subtasks
- Define residue string standard
- Implement or refactor interface extraction utility
- Store receptor-side and partner-side outputs in separate columns or files
- Preserve decoy/model IDs for traceability

## Test Tasks
- Confirm interface residues can be extracted from at least one known model
- Confirm receptor-side and partner-side residues are not mixed
- Confirm residue string format is stable and reusable
- Confirm outputs can be joined back to score tables by model/decoy ID

## Dependencies
Depends on Task Group 1.1.

## Deliverables
- `pyrosetta_interface_residue_table.csv`
- `pyrosetta_interface_models.csv`

---

# Task Group 1.3: Cluster-Level Interface Consensus
**Priority:** Must-Have  
**Zone:** ??? Green

## Objective
Aggregate receptor-side interface evidence across PyRosetta clusters to identify stable receptor-side interface patches and hotspot residues.

## Main Tasks

### 1.3.1 Standardize cluster identity
- Ensure cluster IDs are stable and traceable
- Link cluster summaries back to member models

### 1.3.2 Compute receptor-side residue occupancy
- For each receptor-side residue, compute:
  - occupancy across models within a cluster
  - occupancy across selected top clusters if needed

### 1.3.3 Compute cluster summary fields
- cluster_id
- n_members
- representative_model
- receptor hotspot residues
- optional interface metrics such as mean dG, mean dSASA, mean sc if available

### 1.3.4 Define hotspot candidate rule
- Add a configurable rule for hotspot-like residues
- Start with occupancy-based logic rather than overcomplicated energetic ranking

## Subtasks
- Link cluster membership to interface residue table
- Build cluster-wise receptor residue counts
- Convert counts to occupancy/frequency
- Define top residue selection logic for cluster summaries

## Test Tasks
- Confirm cluster summaries can be generated for at least one receptor state
- Confirm receptor-side occupancy values are computed reproducibly
- Confirm hotspot residue lists remain traceable back to raw models
- Confirm cluster summary files are readable without raw-file inspection

## Dependencies
Depends on Task Group 1.2.

## Deliverables
- `ppi_cluster_summary.csv`
- `ppi_hotspot_residues.csv`
- `ppi_interface_patch_table.csv`

---

# Task Group 1.4: Auxiliary AlphaFold-Multimer Support
**Priority:** Should-Have  
**Zone:** ??? Yellow

## Objective
Standardize AFM outputs as auxiliary receptor-side interface evidence that can be compared with PyRosetta but does not replace it.

## Main Tasks

### 1.4.1 Audit current AFM outputs
- Identify available result files
- Identify whether ranking, pLDDT, pTM, ipTM, and PAE are accessible

### 1.4.2 Standardize AFM model summary
- Model ID
- rank
- ipTM
- pTM
- mean pLDDT
- optional confidence notes

### 1.4.3 Extract receptor-side AFM contacts
- Extract receptor-side contact/interface residues from AFM structures
- Store them in a separate auxiliary evidence table

### 1.4.4 Add confidence-aware labeling
- Keep low-confidence models if useful
- Mark them clearly rather than silently discarding them

## Subtasks
- Inventory AFM output files
- Define AFM summary schema
- Define receptor-side residue extraction for AFM models
- Define auxiliary evidence label conventions

## Test Tasks
- Confirm at least one AFM run can produce a standardized summary
- Confirm receptor-side contact residues can be extracted from AFM models
- Confirm low-confidence outputs remain flagged rather than disappearing
- Confirm AFM outputs are clearly separated from PyRosetta primary evidence

## Dependencies
Depends on availability of AFM outputs and can proceed in parallel after Task Group 1.0.

## Deliverables
- `afm_model_summary.csv`
- `afm_interface_support_table.csv`
- `phase1_afm_support_note.md`

---

# Task Group 1.5: Multi-State Interface Patch Comparison
**Priority:** Must-Have  
**Zone:** ??? Yellow

## Objective
Compare receptor-side interface patches across the three receptor states to determine which receptor-side residues or patches are state-robust and which are state-specific.

## Main Tasks

### 1.5.1 Normalize cross-state residue comparison inputs
- Ensure receptor-side residue tables from each state can be compared
- Surface numbering mismatch warnings if comparison is unsafe

### 1.5.2 Compare receptor-side patch overlap
- Compare hotspot residue overlap across receptor states
- Compare patch centroids or patch summaries if available
- Preserve raw overlap metrics instead of only yes/no conclusions

### 1.5.3 Build state-robustness summary
- Identify residues or local patches that recur across states
- Distinguish robust patch candidates from state-specific candidates

## Subtasks
- Define cross-state comparison metrics for interface patches
- Build residue overlap table
- Add cross-state patch summary logic
- Preserve warnings where residue mapping is uncertain

## Test Tasks
- Confirm at least one receptor-side patch comparison can be produced across states
- Confirm raw overlap metrics are preserved
- Confirm state-robust vs state-specific patch candidates can be distinguished
- Confirm numbering mismatch warnings are visible where needed

## Dependencies
Depends on Task Group 1.3 and, if used, Task Group 1.4.

## Deliverables
- `ppi_patch_cross_state_comparison.csv`
- `ppi_patch_state_robustness.csv`
- `phase1_interface_comparison_report.md`

---

# Task Group 1.6: Phase 1 Review Report
**Priority:** Should-Have  
**Zone:** ??? Yellow

## Objective
Generate a readable Phase 1 review package that summarizes receptor-side MYO1D attachment evidence before moving into pocket proposal.

## Main Tasks

### 1.6.1 Build a Phase 1 summary report
- Summarize receptor inputs
- Summarize PyRosetta evidence
- Summarize receptor-side hotspot residues
- Summarize cross-state receptor-side patch behavior

### 1.6.2 Keep primary vs auxiliary evidence explicit
- PyRosetta receptor-side patch = primary evidence
- AFM receptor-side support = auxiliary evidence

### 1.6.3 Prepare downstream handoff to Phase 2
- Provide the receptor-side patch definition in a form that Phase 2 can use as reference input

## Subtasks
- Define minimum report sections
- Add hotspot summary tables
- Add state-robustness summary section
- Add downstream-ready patch export or reference note

## Test Tasks
- Confirm the report can be read without opening raw decoy files
- Confirm primary vs auxiliary evidence are labeled correctly
- Confirm the report is sufficient to hand off into Phase 2 candidate pocket mapping

## Dependencies
Depends on Task Groups 1.3, 1.4, and 1.5.

## Deliverables
- `phase1_interface_report.md`
- `phase1_downstream_patch_reference.csv`

---

# Recommended Initial Execution Order for Phase 1

If Phase 1 is implemented incrementally, the recommended order is:

1. **Task Group 1.0** ??receptor and partner input validation  
2. **Task Group 1.1** ??PyRosetta execution standardization  
3. **Task Group 1.2** ??receptor-side interface residue extraction  
4. **Task Group 1.3** ??cluster-level interface consensus  
5. **Task Group 1.5** ??multi-state interface patch comparison  
6. **Task Group 1.4** ??auxiliary AFM support  
7. **Task Group 1.6** ??Phase 1 review report

### Why this order
- You must validate receptor and partner inputs before interface mapping.
- You must standardize PyRosetta outputs before residue extraction.
- You must extract receptor-side residues before consensus building.
- You must build receptor-side consensus before comparing receptor states.
- AFM is useful but not the primary engine, so it can lag behind the PyRosetta core.
- The review report should come last.

---

## Korean Summary

??문서??Phase 1??구현 ??위???쪼갠 task 문서?? ??심?? ligand보다 먼?? **MYO1D가 EGFR ??디??붙는지 receptor-side patch????*의??는 ???*??다. 주요 ??름?? ?...력 검?????PyRosetta ??행 ????????interface residue 추출 ??cluster-level hotspot ??receptor ??태 ???비교 ??AFM 보조 증거 ??Phase 1 ?-약 보고????서??


