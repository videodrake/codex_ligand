> Status note (2026-03-12): This v2 task file is the active Phase 1 task breakdown for the current baseline,
> but current active secondary validation is LightDock. AFM should be treated as legacy optional support only.
> Use this task file for Phase 1 implementation intent, and use `docs/current_pipeline_status.md` for current operational status.

## Context Summary
- Project: EGFR-MYO1D interface-centric perturbation discovery pipeline
- Current Subphase: Phase 1 Task Breakdown (v2 - Structural Input Upgrade)
- Purpose: Convert the Phase 1 PRD v2 into implementation-facing tasks for receptor-side MYO1D interface mapping
- Primary Engine: PyRosetta global docking
- Supporting Engine: LightDock (independent secondary validation)
- Auxiliary Engine: AlphaFold-Multimer (auxiliary evidence only, not core)
- Key Principle: Define the receptor-side MYO1D attachment patch before ligand-site prioritization
- Key Change in v2: Receptor upgraded to full kinase domain, partner extended to ~955, orientation-aware filtering mandatory

---

# Task Breakdown
## Phase 1: PPI-first Interface Mapping (v2)

This document breaks Phase 1 into implementation-facing task groups. The purpose of this phase is to define the receptor-side MYO1D attachment patch in a reproducible, structured, and reviewable format before any ligand-pocket prioritization begins.

**v2 changes from v1:**
- Task Group 1.0 expanded to include full kinase domain preparation and extended beta-meander preparation
- New Task Group 1.2A added: Orientation-Aware Filtering (mandatory, inserted between residue extraction and consensus)
- New Task Group 1.7 added: Pilot Data Comparison Layer
- Task Group 1.4 refocused from AFM to LightDock as primary secondary validation; AFM moved to optional sub-task
- All downstream task groups updated to require orientation-validated models only

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

# Task Group 1.0: Receptor and Partner Input Preparation and Validation
**Priority:** Must-Have  
**Zone:** ??? Green

## Objective
Prepare and validate structurally adequate receptor and partner inputs that eliminate known fragment-docking artifacts before PPI-first interface mapping begins.

## Main Tasks

### 1.0.1 Full kinase domain receptor preparation
- For each of the three receptor states (3GT8_raw, EGFR_160-185, EGFR_170-200):
  - Extract or prepare the full kinase domain (~residues 696??79, UniProt numbering)
  - Include both N-lobe and C-lobe as a single chain
  - Confirm activation loop is modeled or resolved
  - Confirm chain ID assignment
  - Record residue numbering system and conversion formula
  - Validate that no missing residues or chain breaks exist in the docking-ready structure

### 1.0.2 Full kinase domain receptor validation
- Confirm all three receptor states are comparable:
  - Same chain ID convention
  - Same residue numbering range
  - Residue sequence identity verified
- Detect and record any numbering mismatches across receptor states
- Confirm that the N-lobe region is structurally complete (no truncated helices or missing loops that would create artificial cavities)

### 1.0.3 Extended beta-meander partner preparation
- Prepare the extended beta-meander construct starting from ~residue 955 (instead of 962)
- Confirm that the extension includes at least 7 residues upstream of the current construct
- Verify structural integrity: no backbone breaks, reasonable phi/psi angles, no steric clashes in the extended region
- Explicitly annotate:
  - Sheet 8 boundaries and residues
  - Sheet 9 boundaries and residues
  - Sheet 10, 11, 12 boundaries and residues
  - Inter-sheet loops
  - Active face definition (sheets 8, 9 ??primary contact face)
  - Structural support face definition (sheet 12 ??not primary direct-contact face per current working assumption)

### 1.0.4 Extended beta-meander partner validation
- Confirm VAL962 is no longer the first residue in the construct
- Confirm that the new N-terminal region does not introduce artifactual contacts (e.g., the new terminal residue should not become a spurious anchor)
- Record partner metadata: source structure, residue range, sheet annotations

### 1.0.5 Metadata persistence
- Build or update receptor metadata output (now including N-lobe/C-lobe boundary info)
- Build or update partner metadata output (now including extended range and sheet annotations)
- Preserve source paths, structural notes, and preparation method

### 1.0.6 Pilot data input registration
- Register the existing C-lobe fragment × beta-meander(962??006) results as pilot/reference data
- Record pilot data construct specifications for later comparison
- Do not modify or overwrite pilot data files

## Subtasks
- Define receptor metadata schema for Phase 1 v2 (including construct type field)
- Define MYO1D partner metadata schema (including sheet annotation fields)
- Add numbering/chain validation checks for full kinase domain
- Add structural integrity checks for extended beta-meander
- Define pilot data reference schema

## Test Tasks
- Confirm all three full-kinase-domain receptor states can be loaded and identified
- Confirm N-lobe and C-lobe are both present in each receptor structure
- Confirm extended beta-meander can be loaded and VAL962 is not the first residue
- Confirm chain and residue metadata are emitted in a stable format
- Confirm partner input files are discovered and labeled consistently
- Confirm numbering mismatches surface as warnings or explicit validation output
- Confirm pilot data files are registered but not modified

## Dependencies
Depends on the repository/document baseline and existing receptor input availability.
Full kinase domain preparation may require structure extraction from 3GT8 PDB + activation loop modeling.
Extended beta-meander preparation may require upstream sequence modeling from the MYO1D TH1 domain structure.

## Deliverables
- `receptor_metadata.csv` (v2: includes construct_type = full_kinase_domain)
- `partner_metadata.csv` (v2: includes extended range, sheet annotations)
- `pilot_data_reference.csv`
- `phase1_input_validation_report.md`

---

# Task Group 1.1: PyRosetta Global Docking Standardization
**Priority:** Must-Have  
**Zone:** ??? Green

## Objective
Standardize PyRosetta global docking execution and output handling so receptor-side interface evidence can be extracted reproducibly from full-kinase-domain docking runs.

## Main Tasks

### 1.1.1 Audit current PyRosetta execution path
- Identify current entry points
- Identify how decoys are generated
- Identify current score outputs
- Identify current cluster outputs, if any
- Assess whether current execution path can handle full-kinase-domain inputs without modification

### 1.1.2 Normalize execution inputs
- Ensure receptor ID and partner ID are propagated into run metadata
- Ensure **construct type** (full_kinase_domain vs legacy_clobe_fragment) is recorded
- Ensure run parameters are recorded
- Ensure random seed or reproducibility settings are stored when possible

### 1.1.3 Assess compute scaling for full kinase domain
- Full kinase domain (~280 res) vs C-lobe fragment (45 res) will significantly increase compute time per decoy
- **Scaling estimate:** RosettaDock scoring scales approximately as O(N × M) where N and M are receptor and partner residue counts. For rigid-body perturbation:
  - C-lobe fragment system: 45 × 45 = 2,025 residue-pair evaluations
  - Full kinase domain system: 280 × 52 = ~14,560 residue-pair evaluations (~7× increase)
  - Including solvation terms (which scale worse): expect **8??5× slower per decoy**
- **Target decoy count adjustment:**
  - If per-decoy time increases ~10×, then 1M decoys ??~100K decoys is more realistic for initial run
  - Alternatively, 200K??00K with multi-seed approach (5 seeds × 50K??00K each)
  - Literature suggests 100K decoys is sufficient for global docking with proper filtering (Comprehensive Filtering Strategies guide: 10,000??00,000)
- **Multi-seed strategy preferred:** 5??0 independent seeds × 50K??00K decoys each provides better sampling diversity than single-seed mega-run
- Document the compute trade-off and chosen target
- **All production runs are server-side only** ??current workspace cannot validate performance

### 1.1.4 Standardize raw output placement
- Separate outputs by receptor state and partner construct
- Preserve decoy, score, and clustering artifacts in stable locations
- Add construct-type labels to output directories

### 1.1.5 Standardize score extraction
- Extract or normalize score tables containing at least:
  - decoy_id
  - total score
  - I_sc (interface score) ??preferred primary ranking metric
  - dG_separated
  - dSASA
  - sc (shape complementarity)
  - packstat
  - nres_int
  - delta_unsatHbonds
  - receptor_id
  - partner_id
  - construct_type

## Subtasks
- Inventory current PyRosetta score file conventions
- Define normalized output path convention for v2
- Add or update run metadata file generation
- Implement a standardized PyRosetta score export if missing
- Add compute time estimation for full-kinase-domain runs

## Test Tasks
- Confirm a PyRosetta run can be initiated with full-kinase-domain input
- Confirm receptor_id, partner_id, and construct_type appear in downstream outputs
- Confirm score tables can be parsed without manual file reconstruction
- Confirm output placement does not mix receptor states or construct types

## Dependencies
Depends on Task Group 1.0.

## Deliverables
- `pyrosetta_run_metadata.json` (v2: includes construct_type)
- `pyrosetta_decoy_scores.csv` (v2: includes I_sc, expanded metrics)
- `phase1_pyrosetta_execution_note.md`
- `compute_scaling_estimate.md`

---

# Task Group 1.2: Receptor-Side Interface Residue Extraction
**Priority:** Must-Have  
**Zone:** ??? Green

## Objective
Extract receptor-side and partner-side interface residues from PyRosetta models in a structured, reusable format.

## Main Tasks

### 1.2.1 Define interface extraction rule
- Choose and document the receptor?-partner contact rule (recommended: 8 Å Cα distance or InterfaceAnalyzerMover-based)
- Keep receptor-side and partner-side residues separate
- Keep residue string formatting consistent across outputs
- For full kinase domain models, receptor-side residues may include N-lobe residues ??these must be preserved, not filtered out

### 1.2.2 Extract per-model interface residues
- For each decoy or selected model, extract:
  - receptor_interface_residues (may include N-lobe and C-lobe residues)
  - partner_interface_residues
  - residue counts
  - per-residue interface metrics if available (per_residue_dG, per_residue_dSASA)
  - receptor_lobe_label per residue (N-lobe vs C-lobe, based on residue number boundary)

### 1.2.3 Build receptor-side residue table
- Emit a long-form table that can be aggregated later by residue frequency and occupancy
- Include lobe annotation for each receptor residue

## Subtasks
- Define residue string standard
- Define N-lobe/C-lobe boundary residue number for annotation
- Implement or refactor interface extraction utility for full kinase domain
- Store receptor-side and partner-side outputs in separate columns or files
- Preserve decoy/model IDs for traceability

## Test Tasks
- Confirm interface residues can be extracted from at least one full-kinase-domain model
- Confirm receptor-side and partner-side residues are not mixed
- Confirm N-lobe residues are captured if present in interface (not silently discarded)
- Confirm residue string format is stable and reusable
- Confirm outputs can be joined back to score tables by model/decoy ID

## Dependencies
Depends on Task Group 1.1.

## Deliverables
- `pyrosetta_interface_residue_table.csv` (v2: includes lobe_label column)
- `pyrosetta_interface_models.csv`

---

# Task Group 1.2A: Orientation-Aware Filtering
**Priority:** Must-Have  
**Zone:** ??? Green

## Objective
Implement and apply mandatory orientation-aware filtering to all PyRosetta PPI models before they are accepted as valid interface evidence. This eliminates face-flipped poses where sheet 8/9 residues contact the receptor but the active face points away from the binding surface.

## Rationale
The beta-meander has a thin, flat β-sheet geometry. In global docking, this means the structure can easily land on the receptor in two orientations: active-face-down (correct) or active-face-up (flipped). A simple residue contact count cannot distinguish these because the same residues may make edge-on or backbone contacts in either orientation. Prior pilot data did not include this filter, which may have allowed face-flipped poses to contaminate cluster-level consensus.

## Main Tasks

### 1.2A.1 Define active face and forbidden face
- **Active face:** The surface of sheets 8 and 9 that should face the receptor. This is the face where the functional side chains of the experimentally validated residues (8th, 9th sheet alanine mutants that abolished decoy function) are exposed.
- **Forbidden face (back face):** The opposite surface of the beta-meander, primarily the back of sheets 10, 11 and some of sheet 12's structural scaffold.
- Define these faces using:
  - Cα?-Cβ vector direction for representative residues on each face
  - Or: surface normal vector computed from sheet 8/9 plane
  - Or: centroid-to-residue vectors for key active-face residues

### 1.2A.2 Implement face orientation metric
- For each docked pose, compute a face orientation score that indicates whether the active face is directed toward the receptor surface or away from it
- Recommended approach:
  1. Compute the centroid of sheet 8/9 Cα atoms
  2. Compute the mean Cα?-Cβ vector for active-face residues (pointing outward from the sheet plane toward the functional side)
  3. Compute the vector from the sheet 8/9 centroid toward the receptor centroid (or nearest receptor surface)
  4. The dot product of (mean Cα?-Cβ vector) and (sheet-to-receptor vector) indicates orientation:
     - Positive: active face toward receptor (correct)
     - Negative: active face away from receptor (flipped)
- Store the raw orientation score for each model

### 1.2A.3 Define pass/fail threshold
- Define a configurable threshold for the orientation score
- Conservative default: reject models where the orientation dot product is negative (face clearly flipped)
- Edge cases (near-zero dot product) should be flagged as `orientation_ambiguous` rather than silently accepted or rejected
- The threshold must be reviewable and adjustable

### 1.2A.4 Apply orientation filter and log results
- Apply the orientation filter to all models that passed basic energy/interface filters
- Log:
  - model_id
  - orientation_score (raw)
  - orientation_class (pass / fail / ambiguous)
  - active_face_residues_in_contact (count)
  - back_face_residues_in_contact (count)
- Models that fail orientation filtering are excluded from consensus building but preserved in raw data

### 1.2A.5 Validate filter against pilot data
- Apply the orientation filter retroactively to the existing 5 valid pilot structures (C02_M01, C02_M03, C04_M01, C04_M02, C07_M03)
- Report how many pilot structures would have passed/failed/ambiguous
- This provides a calibration baseline for the new filter

## Subtasks
- Define representative active-face residue set for Cα?-Cβ vector computation
- Implement orientation metric computation function
- Define orientation class labels and threshold
- Add orientation columns to interface residue table
- Build retroactive pilot validation report

## Test Tasks
- Confirm orientation score can be computed for at least one known model
- Confirm a clearly face-flipped pose can be detected and rejected
- Confirm a clearly correct-orientation pose passes the filter
- Confirm ambiguous cases are flagged, not silently binned
- Confirm the filter does not reject all models (sanity check)
- Confirm pilot structures produce orientation scores consistent with their known geometry

## Dependencies
Depends on Task Group 1.2 (interface residues must be extracted first).

## Deliverables
- `orientation_filter_log.csv`
- `orientation_filter_pilot_validation.csv`
- `phase1_orientation_filter_note.md`

---

# Task Group 1.3: Cluster-Level Interface Consensus
**Priority:** Must-Have  
**Zone:** ??? Green

## Objective
Aggregate receptor-side interface evidence across PyRosetta clusters to identify stable receptor-side interface patches and hotspot residues. **Only orientation-validated models contribute to consensus.**

## Main Tasks

### 1.3.1 Standardize cluster identity
- Ensure cluster IDs are stable and traceable
- Link cluster summaries back to member models
- Record how many members per cluster passed orientation filtering vs total

### 1.3.2 Compute receptor-side residue occupancy
- For each receptor-side residue, compute:
  - occupancy across **orientation-validated** models within a cluster
  - occupancy across selected top clusters if needed
- Separately record N-lobe vs C-lobe residue occupancy

### 1.3.3 Compute cluster summary fields
- cluster_id
- n_members_total
- n_members_orientation_valid
- representative_model (must be orientation-validated)
- receptor hotspot residues
- receptor_lobe_distribution (fraction of interface on N-lobe vs C-lobe)
- optional interface metrics: mean dG, mean dSASA, mean sc, mean I_sc

### 1.3.4 Define hotspot candidate rule
- Add a configurable rule for hotspot-like residues
- Start with occupancy-based logic across orientation-validated models
- Hotspot residues must meet minimum occupancy threshold (configurable, default suggestion: ??50% within cluster)
- Separate hotspot lists for receptor-side and partner-side

## Subtasks
- Link cluster membership to interface residue table (orientation-filtered)
- Build cluster-wise receptor residue counts (orientation-validated only)
- Convert counts to occupancy/frequency
- Define top residue selection logic for cluster summaries
- Add N-lobe/C-lobe distribution tracking

## Test Tasks
- Confirm cluster summaries can be generated for at least one receptor state
- Confirm receptor-side occupancy values are computed only from orientation-validated models
- Confirm hotspot residue lists remain traceable back to raw models
- Confirm cluster summary files are readable without raw-file inspection
- Confirm N-lobe interface residues are captured if they appear

## Dependencies
Depends on Task Group 1.2A (orientation filtering must be complete).

## Deliverables
- `ppi_cluster_summary.csv` (v2: includes orientation-valid counts, lobe distribution)
- `ppi_hotspot_residues.csv`
- `ppi_interface_patch_table.csv`

---

# Task Group 1.4: LightDock Secondary Validation
**Priority:** Must-Have  
**Zone:** ??? Yellow

## Objective
Use LightDock as an independent secondary validation method to provide method-independence evidence for receptor-side patch identification. LightDock results are compared with PyRosetta results but do not replace them.

## Main Tasks

### 1.4.1 LightDock execution
- Run LightDock with the same full-kinase-domain receptor and extended beta-meander partner inputs
- Use appropriate LightDock scoring function and swarm configuration
- Record run parameters and metadata

### 1.4.2 LightDock interface extraction
- Extract receptor-side contact/interface residues from LightDock top models
- Store them in a separate secondary evidence table
- Apply the same orientation-aware filtering logic (or equivalent) to LightDock poses

### 1.4.3 Cross-method convergence analysis
- Compare receptor-side interface residues from PyRosetta and LightDock
- Identify convergent residues (appearing in top clusters of both methods)
- Identify method-specific residues (appearing only in one method)
- Quantify convergence using residue overlap metrics

### 1.4.4 (Optional) Auxiliary AlphaFold-Multimer support
- If AFM outputs are available, standardize and compare as a third auxiliary evidence layer
- AFM is not part of the core workflow but may provide additional structural context
- AFM outputs must be clearly labeled as auxiliary and not mixed with primary evidence

## Subtasks
- Define LightDock run configuration for this system
- Define interface extraction procedure for LightDock outputs
- Define convergence metric (e.g., Jaccard index of top-N receptor residues)
- (Optional) Inventory and standardize any available AFM outputs

## Test Tasks
- Confirm LightDock can be executed with full-kinase-domain inputs
- Confirm receptor-side residues can be extracted from LightDock models
- Confirm convergence comparison can be generated between PyRosetta and LightDock
- Confirm LightDock evidence is clearly separated from PyRosetta primary evidence

## Dependencies
Depends on Task Group 1.0 (same inputs as PyRosetta). Can proceed in parallel after 1.0.

## Deliverables
- `lightdock_run_metadata.json`
- `lightdock_interface_support_table.csv`
- `cross_method_convergence.csv`
- `phase1_lightdock_validation_note.md`
- (Optional) `afm_interface_support_table.csv`, `afm_model_summary.csv`

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
- Confirm that all states used the same receptor construct type (full_kinase_domain)

### 1.5.2 Compare receptor-side patch overlap
- Compare hotspot residue overlap across receptor states
- Compare patch centroids or patch summaries if available
- Preserve raw overlap metrics instead of only yes/no conclusions
- Track whether patches span only C-lobe, or include N-lobe residues

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
- Summarize receptor inputs (full kinase domain, three states)
- Summarize partner inputs (extended beta-meander)
- Summarize PyRosetta evidence (orientation-filtered)
- Summarize receptor-side hotspot residues
- Summarize cross-state receptor-side patch behavior
- Summarize LightDock convergence findings

### 1.6.2 Keep primary vs secondary vs auxiliary evidence explicit
- PyRosetta receptor-side patch = primary evidence
- LightDock receptor-side support = secondary independent validation
- AFM receptor-side support = auxiliary evidence (if available)

### 1.6.3 Prepare downstream handoff to Phase 2
- Provide the receptor-side patch definition in a form that Phase 2 can use as reference input
- Include patch confidence classification
- Include cross-state robustness information
- All patch definitions derived entirely from new full-kinase-domain data

## Subtasks
- Define minimum report sections
- Add hotspot summary tables
- Add state-robustness summary section
- Add downstream-ready patch export or reference note

## Test Tasks
- Confirm the report can be read without opening raw decoy files
- Confirm primary vs secondary vs auxiliary evidence are labeled correctly
- Confirm the report is sufficient to hand off into Phase 2 candidate pocket mapping

## Dependencies
Depends on Task Groups 1.3, 1.4, and 1.5. Optionally includes 1.7 if pilot comparison was performed.

## Deliverables
- `phase1_interface_report.md`
- `phase1_downstream_patch_reference.csv`
  - **v2 schema requirements for Phase 2 compatibility:** This file must include `construct_type` and `orientation_validation_status` fields. Phase 2 Task Group 2.0 now consumes these fields; legacy or skipped-filter runs may still carry `orientation_validation_status = not_available`, which should surface as a compatibility warning at ingestion time.

---

# Task Group 1.7: Pilot Data Comparison Layer (Optional Reference)
**Priority:** Should-Have  
**Zone:** ??? Yellow

## Objective
If methodologically useful, compare new full-kinase-domain docking results against existing C-lobe fragment pilot data. This comparison is **informational only** and must not constrain the interpretation of new results.

**Critical framing rule:** The new full-kinase-domain + extended-beta-meander docking is a fresh start with a fundamentally different system. Pilot data (C-lobe fragment 45 res × truncated beta-meander 962??006) carries known systematic artifacts (N-lobe absence, VAL962 terminal artifact, no orientation filtering). New results that differ from pilot data are not failures ??they are expected improvements.

## Main Tasks

### 1.7.1 Register pilot data as historical reference
- Record pilot data construct specifications
- Label all pilot outputs as `legacy_pilot_fragment_system`
- Do not integrate pilot data into new consensus or new patch definitions

### 1.7.2 Optional methodological comparison
- If useful, note whether new top clusters occupy similar or different receptor surface regions compared to pilot clusters
- This comparison is for understanding how system improvements changed results, not for validating new results against old ones

### 1.7.3 N-terminal artifact resolution
- Note whether any residue that was the first residue in the pilot construct (VAL962) still appears as a dominant contact in the extended construct
- If it does not, record this as expected resolution of the truncation artifact
- If it does, note this as independent confirmation (but do not over-interpret based on pilot precedent)

## Subtasks
- Register pilot data files with legacy labels
- Define optional comparison metrics (if comparison is performed)

## Test Tasks
- Confirm pilot data files are registered but not mixed into new consensus
- Confirm new results are interpreted independently of pilot expectations

## Dependencies
Depends on Task Group 1.3 (new results must exist first).

## Deliverables
- `pilot_data_reference.csv` (registration only)
- `phase1_pilot_comparison_note.md` (optional, if comparison performed)

---

# Recommended Initial Execution Order for Phase 1 v2

If Phase 1 is implemented incrementally, the recommended order is:

1. **Task Group 1.0** ??receptor and partner input preparation and validation
2. **Task Group 1.1** ??PyRosetta execution standardization (including compute scaling assessment)
3. **Task Group 1.2** ??receptor-side interface residue extraction
4. **Task Group 1.2A** ??orientation-aware filtering (mandatory)
5. **Task Group 1.3** ??cluster-level interface consensus (orientation-filtered)
6. **Task Group 1.5** ??multi-state interface patch comparison
7. **Task Group 1.4** ??LightDock secondary validation (can overlap with 1.2??.3)
8. **Task Group 1.7** ??pilot data comparison layer
9. **Task Group 1.6** ??Phase 1 review report

### Why this order
- You must prepare and validate full-kinase-domain + extended-beta-meander inputs before any docking.
- You must standardize PyRosetta execution before residue extraction.
- You must extract interface residues before orientation filtering.
- **Orientation filtering must happen before consensus building** ??this is the key v2 insertion point.
- You must build orientation-validated consensus before comparing receptor states.
- LightDock can run in parallel after inputs are prepared, but convergence analysis needs PyRosetta results.
- Pilot comparison requires both new results and registered pilot data.
- The review report comes last, incorporating all evidence layers.

### Compute budget note
Full-kinase-domain docking will be substantially more expensive per decoy than C-lobe fragment docking. The project should plan for:
- Reduced decoy count if needed (100K??00K may still be sufficient with improved filtering)
- Multi-seed approach (5??0 seeds) rather than single-seed mega-run
- Server-side execution only (current Codex workspace is not suitable for production runs)

---

## Korean Summary (간단 ?-약)

??문서??Phase 1??구현 ??위???쪼갠 task 문서 v2?? v1 ???????심 변경사??

- **Task 1.0:** Receptor?????장 kinase domain??로, partner???~955부????작??는 ?*장 beta-meander????...그??이??
- **Task 1.2A (? 규):** Orientation-aware filtering????수 ??계?????입. Sheet 8/9 active face가 receptor?????하?-?? 기하??적??로 검???
- **Task 1.4:** AFM ????LightDock????...립 secondary validation??로 ??배???
- **Task 1.7 (? 규):** 기존 C-lobe fragment pilot data?? ??결과??체계??비교. VAL962 artifact ??정 ??함
- 모든 consensus ???downstream ?'업?? orientation-validated model?????용
- ??장 kinase domain ??킹??compute cost 증???????????고, decoy count 조정???multi-seed ??략??계획

??심 ??름: ?...력 준?????PyRosetta ????????interface residue 추출 ??**orientation filtering** ??cluster consensus ??state 비교 ??LightDock 검?????pilot 비교 ??최종 보고??



