# Data Flow Guide

This document explains the current repository from a scientist's point of
view.

The goal is not to describe the code module-by-module, but to answer the
questions a structural biology or computational chemistry user usually asks:

- what are the real inputs?
- which tool is used at each stage?
- what files are produced?
- what does each output mean scientifically?
- what can we conclude, and what should we still treat cautiously?

## 0. Step Interpretation Layer

Production runs now expose two output layers under the same project root:

- canonical runtime outputs under `output/{project}/`
- additive step folders that act as a derived interpretation view

Start reading a completed production run in this order:

1. `output/{project}/step_index.md`
2. `output/{project}/step6_report/project_report.txt`
3. `output/{project}/step5_verdict/valid_sites.csv`
4. `output/{project}/step4_vina_postprocess/vina_pocket_table.csv`
5. `output/{project}/step3_ppi_postprocess/ppi_pyrosetta_residues.csv`

Canonical runtime outputs remain the source of truth. The step layer is additive and can be regenerated from canonical outputs. Large raw PyRosetta directories are indexed by manifest rather than duplicated into `step2_ppi_raw/`.

### 0.1 Phase-To-Step Mapping

| Production phase | Derived step folder | Interpretation role |
|------|------|------|
| Phase 1 | `step1_vina_raw/` | Raw Vina pose inventory |
| Phase 2 | `step2_ppi_raw/` | PyRosetta ranking summaries and raw path index |
| Phase 3 | `step3_ppi_postprocess/` | Receptor-side residue evidence |
| Phase 4 | `step4_vina_postprocess/` | Pocket-level Vina interpretation |
| Phase 5 | `step5_verdict/` | Final site prioritization |
| Phase 6 | `step6_report/` | Narrative report-first view |
| Phase 7 | `step7_validate/` | Persisted validation status |

### 0.2 Reference Output Layout

```text
output/egfr_myo1d_vina/
  vina_pose_table.csv
  vina_pocket_table.csv
  vina_drug_pocket_map.csv
  ppi_pyrosetta_residues.csv
  ppi_pyrosetta_summary.csv
  valid_sites.csv
  cross_method_agreement.csv
  combined_residue_evidence.csv
  project_report.txt
  step_index.md
  current_run_manifest.json
  step1_vina_raw/
  step2_ppi_raw/
  step3_ppi_postprocess/
  step4_vina_postprocess/
  step5_verdict/
  step6_report/
  step7_validate/
```

## 1. Scientific Question

This project is an EGFR-MYO1D state-comparison pipeline.

The central question is not just "where can a ligand bind on EGFR?".
The actual goal is broader:

- compare three EGFR receptor states
- define receptor-side MYO1D interaction evidence
- identify pocket proposals relative to that receptor-side patch
- test ligand-pose diversity across receptor states
- integrate these layers into a final mechanistic interpretation

The current routine baseline combines:

- a Vina-centered ligand evidence layer
- a PyRosetta-centered Phase 1 PPI evidence layer
- LightDock as active secondary validation for Phase 1
- downstream integration through verdict, report, and validation outputs

AlphaFold-Multimer is not part of the current routine baseline.

## 2. Real Inputs

### 2.1 Receptor states

The current active receptor set contains exactly three states:

- `input/receptors/3GT8_raw.pdb`
- `input/receptors/EGFR_160-185.pdb`
- `input/receptors/EGFR_170-200.pdb`

These correspond to:

- `3GT8_raw`
- `EGFR_160-185`
- `EGFR_170-200`

For Vina runs, matching receptor `.pdbqt` files are also expected or prepared:

- `input/receptors/3GT8_raw_receptor.pdbqt`
- `input/receptors/EGFR_160-185_receptor.pdbqt`
- `input/receptors/EGFR_170-200_receptor.pdbqt`

### 2.2 Ligands

The current ligand set is defined in:

- `config/example-project.yaml`

Current baseline ligand inputs are:

- `input/ligands/173940_ligand.sdf`
- `input/ligands/97806_ligand.sdf`
- `input/ligands/VAX-C12_0_ligand.sdf`

Matching ligand `.pdbqt` files are expected or prepared for docking.

### 2.3 PPI inputs

The receptor-side PPI branch uses MYO1D partner inputs under:

- `input/PPI/`

Important prepared structures include:

- `input/PPI/prepared/EGFR_dimer_TH1.pdb`
- `input/PPI/prepared/EGFR_dimer_beta_meander.pdb`

These are the docking-ready complexes or paired inputs used for Phase 1
receptor-side patch definition.

### 2.4 Runtime config

The current project-level routine config is:

- `config/example-project.yaml`

Key current baseline values include:

- 3 receptor states
- ligand list
- `max_workers: 16`
- Vina mode and scoring parameters
- postprocess thresholds
- references to PyRosetta PPI result directories

## 3. The Pipeline Has Two Major Evidence Branches

Scientifically, the repository is easiest to understand as two major evidence
branches that later meet:

1. Vina-centered ligand and pocket evidence
2. Phase 1 receptor-side PPI evidence

Those branches are later integrated into:

3. verdict / report / validation outputs

MD exists as a downstream stability-gate concept, but it is not currently the
first automatic production path.

## 4. Branch A: Vina-Centered Ligand Evidence

### 4.1 Docking input

Tool:

- AutoDock Vina Python API

Code path:

- `egfr_pipeline/vina/dock.py`
- `main.py vina`
- `run_production.py` Phase 1

What goes in:

- receptor PDBQT for each receptor state
- ligand PDBQT for each ligand
- Vina parameters from `config/example-project.yaml`

What is produced first:

- raw docking pose files (`.pdbqt`)

Scientific meaning:

- these are the direct ligand-pose hypotheses
- they tell us where each ligand can be placed in each receptor state
- by themselves they are still raw pose collections, not pocket-level claims

### 4.2 Pose parsing

Tool:

- `egfr_pipeline/vina/parse_poses.py`

Input:

- raw Vina pose `.pdbqt` files

Output:

- `vina_pose_table.csv`

What is extracted:

- receptor ID
- ligand ID
- pose rank
- affinity
- RMSD bounds
- ligand pose centroid

Scientific meaning:

- this converts raw docking text into a pose-level table
- each row is now a machine-readable ligand pose hypothesis
- we can compare pose counts and affinities across receptor states

### 4.3 Contact extraction

Tool:

- `egfr_pipeline/vina/contacts.py`

Input:

- `vina_pose_table.csv`
- receptor PDB coordinates

Output:

- updated `vina_pose_table.csv` with contact residue fields

Key fields:

- `contact_residues`
- `n_contact_residues`
- `contact_distances` (per-residue minimum distance in semicolon-delimited format)

Additional output:

- `vina_contact_distances.csv` (long-form: receptor_id, ligand_id, pose_rank, residue_id, min_distance_A)

Scientific meaning:

- this step tells us which receptor residues physically surround a pose
- the distance information enables distance-weighted analysis and threshold tuning
- this is the first place where a pose becomes biologically interpretable
- we can move from "a pose exists here" to "this pose touches these residues at known distances"

### 4.4 Pocket clustering

Tool:

- `egfr_pipeline/vina/cluster.py`

Input:

- pose-level rows with centroids and contact residues

Output:

- updated `vina_pose_table.csv` with `pocket_id`
- `vina_clustering_merge_log.csv` (pocket merge provenance with merge reasons)
- `vina_clustering_parameters.json` (clustering parameter snapshot)

Scientific meaning:

- multiple poses are grouped into recurring spatial pockets
- merge provenance is preserved so clustering decisions can be audited
- this reduces pose-level clutter into pocket-level patterns
- the clustering is receptor-local, so each receptor state keeps its own pocket
  organization

### 4.5 Pocket summarization

Tool:

- `egfr_pipeline/vina/summarize.py`

Outputs:

- `vina_pocket_table.csv`
- `vina_drug_pocket_map.csv`
- `vina_pocket_residue_occupancy.csv`

What these mean:

- `vina_pocket_table.csv`
  - pocket-level summary per receptor
  - centroid
  - number of poses
  - number of ligands
  - best and mean affinity
  - union of contact residues
- `vina_drug_pocket_map.csv`
  - ligand-to-dominant-pocket mapping
  - tells us which pocket each ligand prefers in each receptor state
- `vina_pocket_residue_occupancy.csv`
  - per-pocket residue occupancy counts and fractions
  - hotspot flag for residues appearing in a high fraction of poses

Scientific meaning:

- the pose cloud is now summarized into interpretable pockets
- residue occupancy reveals which residues are most consistently contacted
- we can ask whether a pocket is broad or narrow
- we can ask whether multiple ligands converge on the same pocket
- we can ask whether a ligand changes dominant pocket across receptor states

### 4.6 Cross-receptor comparison

Tool:

- `egfr_pipeline/vina/compare.py`

Output:

- `vina_pocket_comparison.csv`

Scientific meaning:

- this compares pockets across `3GT8_raw`, `EGFR_160-185`, and
  `EGFR_170-200`
- it tells us whether a pocket is conserved, shifted, or state-specific
- this is central for state-comparison interpretation

### 4.7 Optional bootstrap stability check

Tool:

- `egfr_pipeline/vina/bootstrap.py`

Output:

- `vina_pocket_bootstrap.csv`

Scientific meaning:

- this asks whether the pocket pattern is stable under resampling
- it does not create a new biological result by itself
- it adds confidence or caution around pocket reproducibility

## 5. Branch B: Phase 1 Receptor-Side PPI Evidence

This branch addresses a different scientific question:

- where does MYO1D likely engage the receptor side of EGFR?

That information is later used to judge whether ligand pockets are relevant to
the receptor-side attachment patch.

### 5.1 PyRosetta global docking

Tool:

- PyRosetta / RosettaDock

Code path:

- `egfr_pipeline/pyrosetta_docking/pipeline_manager.py`
- `egfr_pipeline/pyrosetta_docking/docking.py`
- `config/phase1/*.ini` (Phase 1 PyRosetta configs)

Primary outputs include:

- `pyrosetta_run_metadata.json`
- `pyrosetta_decoy_scores.csv`
- `scored_all_models.csv` (all models with Pass 1 metrics + filter_status for post-hoc re-analysis)
- `scored_stage2_models.csv` (Stage 2 candidates with expensive metrics)
- `filter_thresholds.csv` (filter parameters and pass/reject counts)
- `cluster_membership.csv` (full model-to-cluster mapping)
- `*_ContactPairs.csv` (per-residue-pair minimum distances for final models)
- decoy structural outputs inside Phase 1 PPI output directories

Scientific meaning:

- this generates many receptor-partner docking hypotheses
- it is the primary structural evidence layer for receptor-side patch mapping
- the comprehensive metric capture (dSASA polar/hydrophobic decomposition, I_RMSD,
  per-residue energy with 9 terms, contact pair distances) enables post-hoc
  re-analysis without re-docking
- it is not a ligand docking step; it is a protein-protein interface search

### 5.2 Interface extraction and filtering

Tools:

- `egfr_pipeline/pyrosetta_docking/analysis.py`
- `egfr_pipeline/ppi/pyrosetta_extract.py`
- Phase 1 consensus modules

Outputs include:

- `ppi_pyrosetta_residues.csv`
- `ppi_pyrosetta_summary.csv`
- `ppi_cluster_summary.csv` (with centroid_x/y/z, centroid_spread_A, energy distribution stats)
- `ppi_hotspot_residues.csv`
- `ppi_interface_patch_table.csv`

Scientific meaning:

- this turns many raw PPI docking decoys into residue-level receptor-side patch
  evidence
- instead of asking "which docking model is best?", we ask:
  - which receptor residues repeatedly participate in the interface?
  - which residues are robust across states?
  - which residues look like stable candidate patch residues?

### 5.3 Cross-state Phase 1 interpretation

Tools:

- `egfr_pipeline/phase1/cluster_consensus.py`
- `egfr_pipeline/phase1/compare_states.py`
- `egfr_pipeline/phase1/review_report.py`

Outputs include:

- `ppi_patch_cross_state_comparison.csv`
- `ppi_patch_state_robustness.csv`
- `phase1_downstream_patch_reference.csv`

Scientific meaning:

- this is where receptor-side patch evidence becomes state-comparison evidence
- the patch can now be described as:
  - robust across states
  - shifted across states
  - state-specific
- this is the bridge from protein-protein docking to later pocket relevance

## 6. LightDock as Phase 1 Secondary Validation

Tool:

- LightDock

Code path:

- `egfr_pipeline/phase1/lightdock_validation.py`

Outputs include:

- `lightdock_run_metadata.json`
- `lightdock_interface_support_table.csv`
- `lightdock_model_summary.csv`
- `cross_method_convergence.csv`

Scientific meaning:

- LightDock is not the primary truth source
- it is an independent secondary validation axis
- it asks whether a separate docking method supports similar receptor-side
  residues or patch regions

Interpretation rule:

- PyRosetta remains primary evidence
- LightDock adds method-independence support
- LightDock-only residues should be treated cautiously and not promoted as
  primary truth on their own

## 7. How Branch A and Branch B Meet

Once we have:

- Vina pockets and ligand-pocket maps
- Phase 1 receptor-side patch evidence

the project can ask the biologically meaningful question:

- which ligand pockets are relevant to the MYO1D receptor-side patch?

This is where simple docking becomes perturbation relevance.

## 8. Later Phases: Patch-Relevant Pocket Interpretation

### 8.1 Phase 2 pocket proposal and relationship mapping

Tools:

- `egfr_pipeline/phase2/pocket_proposal.py`
- `egfr_pipeline/phase2/pocket_merge.py`
- `egfr_pipeline/phase2/patch_ingestion.py`
- `egfr_pipeline/phase2/patch_relationship.py`
- `egfr_pipeline/phase2/druggability_confidence.py`
- `egfr_pipeline/phase2/cross_state_alignment.py`
- `egfr_pipeline/phase2/phase3_export.py`

Outputs include:

- `candidate_pockets.csv`
- `pocket_patch_relationship.csv`
- `druggability_proposal_summary.csv`
- `candidate_pocket_state_classes.csv`
- `phase3_candidate_pocket_reference.csv`

Scientific meaning:

- this layer asks whether a candidate pocket is:
  - orthosteric-like relative to the patch
  - rim-like
  - allosteric-like
  - low relevance
- it also turns pocket geometry into docking-ready targets for the next phase

### 8.2 Phase 3 diversity-aware docking

Tools:

- `egfr_pipeline/phase3/run_diverse_docking.py`
- `egfr_pipeline/phase3/pose_attribution.py`
- `egfr_pipeline/phase3/review_report.py`

Scientific meaning:

- this phase moves from broad pocket detection to diversity-aware docking over
  curated pocket targets
- the question becomes:
  - which ligands repeatedly support the same mechanistically relevant pocket?
  - which pockets recruit diverse ligands?

### 8.3 Phase 4 perturbation relevance scoring

Tools:

- `egfr_pipeline/phase4/perturbation_scoring.py`
- `egfr_pipeline/phase4/score_framework.py`
- `egfr_pipeline/phase4/mechanistic_classification.py`
- `egfr_pipeline/phase4/final_report.py`

Scientific meaning:

- this is where all evidence layers are integrated
- the output is no longer just "best affinity"
- it becomes a final ranking shaped by:
  - ligand evidence
  - receptor-side patch relevance
  - cross-state behavior
  - mechanistic class

## 9. Final Integration Outputs

Current integrated outputs include:

- `cross_method_agreement.csv`
- `valid_sites.csv` (now includes per-component score breakdown: vina_affinity_pts,
  vina_convergence_pts, vina_consensus_pts, ppi_spatial_pts, ppi_overlap_pts,
  cross_receptor_pts, score_denominator)
- `vina_consensus_sites.csv`
- `combined_residue_evidence.csv`
- `project_report.txt`

Scientific meaning:

- these files summarize the project-level conclusion
- the component-level score breakdown in `valid_sites.csv` enables transparent
  auditing of which evidence axes drove each site's final verdict
- they are the first place to look when asking:
  - which sites look biologically meaningful?
  - which pockets are supported by multiple evidence types?
  - which residues appear repeatedly across methods and states?

## 10. MD as a Downstream Stability Gate

MD is part of the broader scientific design, but it is important to describe it
honestly.

Current reality:

- MD is documented as a downstream stability gate
- MD is not the first automatic production path in the current routine baseline

Scientific meaning:

- MD is used to test whether a proposed complex or interface remains stable
  over time
- this adds dynamic confidence to a structure-derived hypothesis
- it should be interpreted as a later-stage strengthening or rejection step,
  not as the first entry point for repository understanding

## 11. What We Can Learn From This Pipeline

When the current baseline is working properly, the project can tell us:

- whether different EGFR receptor states expose different ligand-accessible
  pockets
- whether those pockets are conserved or state-specific
- whether those pockets lie near, overlap with, or remain distant from the
  receptor-side MYO1D patch
- whether a ligand repeatedly prefers one mechanistically relevant pocket
- whether multiple methods support the same receptor-side residues
- whether the final evidence favors a stronger, weaker, or more uncertain site
  interpretation

## 12. What Still Requires Caution

Even with the current structured outputs, some things still require careful
scientific judgment.

- LightDock is secondary evidence, not replacement primary truth
- AFM code may still exist, but it is not part of the active routine baseline
- residue numbering and chain mapping must remain explicit and traceable
- Vina pockets are useful hypotheses, not direct proof of biological binding
- MD is not automatically part of the current production path
- planning documents may describe a cleaner end-to-end architecture than the
  current code already executes by default

## 13. If You Want The Shortest Possible Summary

This project starts from:

- three receptor states
- a small ligand panel
- PPI partner structures

It then builds two evidence layers:

1. ligand-pocket evidence from Vina
2. receptor-side patch evidence from PyRosetta, checked secondarily by LightDock

Those layers are combined into:

- pocket relevance
- state comparison
- final verdict and report outputs

So the real scientific product is not just docking scores.
It is a structured interpretation of which EGFR surface regions remain
mechanistically interesting across receptor states and how well those claims
are supported by multiple layers of evidence.
