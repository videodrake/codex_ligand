# Data Inventory

This document inventories the current repository data surfaces that a new GPT needs to recognize before touching any analysis or documentation. It focuses on what is actually present in the workspace, what the active config expects at runtime, and where derived outputs currently land.

## Inventory Rules

- Treat paths under `input/` as repository-managed inputs or prepared handoff assets.
- Treat paths under `output/` as derived artifacts.
- Distinguish between files that are physically present in the workspace and files that are only referenced by config or by downstream expectations.
- Treat AFM-related assets as legacy optional context, not as routine baseline input.

## Active Dataset Summary

| Category | Current baseline |
|------|------|
| Receptor states | `3GT8_raw`, `3GT8_cl38_48`, `3GT8_cl85_100` |
| Ligand set in active config | `173940`, `97806`, `VAX-C12_0` |
| Phase 1 primary partner baseline | extended beta-meander derived from `input/PPI/TH1 domain.pdb` |
| Phase 1 secondary validation | LightDock outputs under `output/phase1_ppi/<state>/lightdock/` |
| Legacy optional partner context | TH1 full-domain and AFM-related references |

## 1. Raw Inputs

### 1.1 Receptor structures

These are the current checked-in receptor PDB sources used across the repo.

| Path | Meaning | Current status |
|------|------|------|
| `input/receptors/3GT8_raw.pdb` | Raw 3GT8-derived receptor state | active |
| `input/receptors/3GT8_cl38_48.pdb` | MD cluster representative for 38-48 ns state | active |
| `input/receptors/3GT8_cl85_100.pdb` | MD cluster representative for 85-100 ns state | active |

Notes:

- These three receptor IDs are the fixed state ensemble for current onboarding and interpretation.
- `input/PPI/phase1/receptor_metadata.csv` records the current full-kinase-domain interpretation for all three states as range `699-1007`, with boundary residue `838` used for N-lobe/C-lobe split.

### 1.2 Ligand structures

Current ligand source files present in the repository:

| Path | Meaning | Current status |
|------|------|------|
| `input/ligands/173940_ligand.sdf` | Ligand `173940` structure source | active |
| `input/ligands/97806_ligand.sdf` | Ligand `97806` structure source | active |
| `input/ligands/VAX-C12_0_ligand.sdf` | Ligand `VAX-C12_0` structure source | active |

### 1.3 Raw PPI partner inputs

| Path | Meaning | Current status |
|------|------|------|
| `input/PPI/beta_meander.pdb` | beta-meander partner source | legacy/raw partner source |
| `input/PPI/TH1 domain.pdb` | TH1-domain source file; current extended beta-meander is derived from this | active source for current Phase 1 preparation |

Important interpretation:

- `TH1 domain.pdb` still exists as a full-domain source, but the current Phase 1 baseline centers on the extended beta-meander construct, not on TH1 as the primary search target.
- `beta_meander.pdb` remains relevant as historical raw input and legacy comparison context.

## 2. Prepared Inputs

Prepared inputs are repository-side transformations that sit between raw structures and derived pipeline outputs.

### 2.1 Phase 1 current preparation assets

These files represent the current prepared Phase 1 input layer for the full-kinase-domain plus extended-beta-meander baseline.

| Path | Meaning |
|------|------|
| `input/PPI/phase1/receptor_metadata.csv` | per-state receptor metadata for the current full-kinase-domain Phase 1 baseline |
| `input/PPI/phase1/partner_metadata.csv` | metadata for the current extended beta-meander partner, including sheet annotations |
| `input/PPI/phase1/docking_pair_metadata.csv` | state-specific docking pair assembly inventory |
| `input/PPI/phase1/phase1_input_validation_report.md` | human-readable validation of the current prepared Phase 1 inputs |
| `input/PPI/phase1/pilot_data_reference.csv` | references to older fragment-based pilot results retained only for comparison |

What these files mean:

- `receptor_metadata.csv` is the quickest structured source for receptor numbering, construct type, and lobe counts.
- `partner_metadata.csv` defines the current partner as `extended_beta_meander`, spanning residues `955-1006`, with `VAL962` no longer the first residue.
- `docking_pair_metadata.csv` records the actual receptor-partner assembly units used for current Phase 1 docking preparation.
- `pilot_data_reference.csv` points to legacy fragment-based runs and should not be read as the current scientific baseline.

### 2.2 Legacy prepared dimer assets

These files are prepared PPI assets that still matter for reference and historical traceability.

| Path | Meaning | Current status |
|------|------|------|
| `input/PPI/prepared/EGFR_dimer_beta_meander.pdb` | prepared EGFR dimer plus beta-meander assembly | legacy prepared asset |
| `input/PPI/prepared/EGFR_dimer_beta_meander_info.json` | metadata for the legacy dimer-plus-beta-meander assembly | legacy prepared asset |
| `input/PPI/prepared/EGFR_dimer_beta_meander_mapping.csv` | chain/residue remapping for that legacy prepared assembly | legacy prepared asset |
| `input/PPI/prepared/EGFR_dimer_TH1.pdb` | prepared EGFR dimer plus TH1 assembly | legacy prepared asset |
| `input/PPI/prepared/EGFR_dimer_TH1_info.json` | metadata for the legacy dimer-plus-TH1 assembly | legacy prepared asset |
| `input/PPI/prepared/EGFR_dimer_TH1_mapping.csv` | chain/residue remapping for that legacy prepared assembly | legacy prepared asset |

Important interpretation:

- The `prepared/` directory documents older dimer-centered preparation paths and is still useful for provenance.
- These assets should not outrank the current `input/PPI/phase1/` preparation layer when describing the active baseline.

### 2.3 Config-declared prepared inputs not present in the current workspace snapshot

`config/example-project.yaml` expects these runtime-prepared files:

| Config field pattern | Expected file type | Present in current `input/` tree |
|------|------|------|
| `receptors[*].pdbqt` | receptor PDBQT for Vina | not present |
| `ligands[*].pdbqt` | ligand PDBQT for Vina | not present |

Current interpretation:

- The active config references `.pdbqt` inputs for Vina runs.
- Those `.pdbqt` files are not present in the current checked-in `input/` tree of this workspace snapshot.
- For onboarding purposes, treat them as runtime-required prepared inputs that may be generated or staged outside the current git-tracked snapshot.

## 3. Derived Outputs

## 3.1 Current project root

The active config defines:

- `project_name: egfr_myo1d_vina`
- `output_root: ./output`

So the main project output root is:

- `output/egfr_myo1d_vina/`

Observed contents:

- root-level summary file names such as `vina_pose_table.csv`, `cross_method_agreement.csv`, and `project_report.txt`
- subdirectories `vina/`, `ppi/`, and `results/`

Important caveat:

- Several root-level files in `output/egfr_myo1d_vina/` are thin path-pointer stubs rather than full payload files.
- Example: `output/egfr_myo1d_vina/vina_pose_table.csv` contains the string `vina/vina_pose_table.csv`, which points to the actual file under `output/egfr_myo1d_vina/vina/`.
- New GPTs should not assume every root-level `.csv` there is the full dataset body.

### 3.2 Vina-derived outputs

Primary Vina output location:

- `output/egfr_myo1d_vina/vina/`

Key files:

| Path | Meaning |
|------|------|
| `output/egfr_myo1d_vina/vina/vina_pose_table.csv` | pose-level docking table; each row is a receptor-ligand pose with affinity, centroid, raw pose file path, and pocket assignment |
| `output/egfr_myo1d_vina/vina/vina_pocket_table.csv` | pocket-level aggregation of Vina poses within each receptor |
| `output/egfr_myo1d_vina/vina/vina_drug_pocket_map.csv` | mapping from ligands to inferred pockets |
| `output/egfr_myo1d_vina/vina/vina_pocket_comparison.csv` | cross-receptor pocket comparison output |
| `output/egfr_myo1d_vina/vina/vina_pocket_bootstrap.csv` | optional bootstrap stability summary for pocket-level patterns |
| `output/egfr_myo1d_vina/vina/vina_contact_distances.csv` | long-form per-pose contact distances (receptor_id, ligand_id, pose_rank, residue_id, min_distance_A) |
| `output/egfr_myo1d_vina/vina/vina_pocket_residue_occupancy.csv` | per-pocket residue occupancy and hotspot flags |
| `output/egfr_myo1d_vina/vina/vina_clustering_merge_log.csv` | pocket merge provenance log with merge reasons |
| `output/egfr_myo1d_vina/vina/vina_clustering_parameters.json` | clustering parameter snapshot for reproducibility |

### 3.3 Integrated verdict/report outputs

Primary integrated summary location:

- `output/egfr_myo1d_vina/results/`

Key files:

| Path | Meaning |
|------|------|
| `output/egfr_myo1d_vina/results/valid_sites.csv` | final rule-based site verdicts by pocket |
| `output/egfr_myo1d_vina/results/cross_method_agreement.csv` | Vina-PPI overlap and agreement table |
| `output/egfr_myo1d_vina/results/vina_consensus_sites.csv` | consensus pocket/site grouping across receptors |
| `output/egfr_myo1d_vina/results/combined_residue_evidence.csv` | combined residue-level evidence view |
| `output/egfr_myo1d_vina/results/project_report.txt` | human-readable top-level project report |

### 3.4 Legacy or earlier PPI-derived outputs under the project root

Observed PPI location under the project root:

- `output/egfr_myo1d_vina/ppi/`

This contains older or parallel PPI run products such as:

- `ppi_pyrosetta_residues.csv`
- `ppi_pyrosetta_summary.csv`
- legacy run folders like `beta_meander/`

Interpretation:

- This area preserves previous PPI-derived products and historical run outputs.
- For the current structured Phase 1 baseline, prefer the dedicated `output/phase1_ppi/` tree described below.

### 3.5 PyRosetta PPI docking outputs (per-PDB run)

Each PPI docking run produces the following enhanced output files under `<PDB_NAME>/`:

| File | Meaning |
|------|------|
| `scored_all_models.csv` | All models with Pass 1 metrics (dG, dSASA, dSASA_polar, dSASA_hydrophobic, sc, total_score) + filter_status |
| `scored_stage2_models.csv` | Stage 2 candidate models with expensive metrics (packstat, unsatHb, nres_int, hbonds_int) |
| `filter_thresholds.csv` | Filter thresholds applied per stage with input/output counts |
| `cluster_results/cluster_membership.csv` | Full model-to-cluster mapping including non-representative members |
| `final_result/*_ContactPairs.csv` | Per-residue-pair minimum distances for final ranked models |
| `final_ranking.csv` | Comprehensive ranking with dSASA_polar, dSASA_hydrophobic, I_RMSD, center_x/y/z |
| `cluster_results/cluster_summary.csv` | Cluster summaries with centroid coordinates, dSASA decomposition, energy distribution stats |

## 4. Phase-Specific Derived Trees

These output roots reflect the newer phase-separated document and code organization.

### 4.1 Phase 1 structured PPI outputs

Root:

- `output/phase1_ppi/`

Key cross-state files:

| Path | Meaning |
|------|------|
| `output/phase1_ppi/phase1_interface_report.md` | human-readable Phase 1 review report |
| `output/phase1_ppi/phase1_downstream_patch_reference.csv` | Phase 2 handoff patch reference |
| `output/phase1_ppi/ppi_patch_cross_state_comparison.csv` | residue-level cross-state comparison |
| `output/phase1_ppi/ppi_patch_state_robustness.csv` | state robustness labeling for Phase 1 residues |

Per-state subtrees:

- `output/phase1_ppi/3GT8_raw/`
- `output/phase1_ppi/3GT8_cl38_48/`
- `output/phase1_ppi/3GT8_cl85_100/`

Typical per-state files:

| File | Meaning |
|------|------|
| `ppi_cluster_summary.csv` | cluster-level interface summary |
| `ppi_hotspot_residues.csv` | hotspot residue summary |
| `ppi_interface_patch_table.csv` | structured residue occupancy table for the state |
| `pyrosetta_interface_models.csv` | per-model interface metadata |
| `pyrosetta_interface_residue_table.csv` | per-model residue extraction table |
| `orientation_filter_log.csv` | orientation-aware filter log when available |
| `cross_method_convergence.csv` | PyRosetta vs LightDock residue convergence table |
| `lightdock/lightdock_interface_support_table.csv` | LightDock residue support table |
| `lightdock/lightdock_run_metadata.json` | LightDock run metadata |

### 4.2 Phase 2 pocket proposal outputs

Root:

- `output/phase2_pockets/`

Representative files:

| Path | Meaning |
|------|------|
| `candidate_pockets_raw.csv` | unmerged pocket proposals by source tool |
| `candidate_pockets.csv` | normalized candidate pocket catalog |
| `candidate_pocket_merge_table.csv` | merge relationships among raw proposals |
| `candidate_pocket_provenance.csv` | provenance from source proposal tools |
| `pocket_patch_relationship.csv` | relationship of each candidate pocket to the Phase 1 patch |
| `druggability_proposal_summary.csv` | higher-level druggability confidence layer |
| `phase3_candidate_pocket_reference.csv` | handoff file intended for Phase 3 |

### 4.3 Phase 3 diverse docking outputs

Root:

- `output/phase3_docking/`

Representative files:

| Path | Meaning |
|------|------|
| `phase3_docking_job_table.csv` | generated receptor-pocket-ligand job inventory |
| `phase3_job_box_table.csv` | per-pocket docking box definitions |
| `pocket_search_status.csv` | pocket saturation/open status |
| `phase3_budget_tracking.csv` | docking budget accounting |
| `phase3_round_log.csv` | round-level execution log |
| `vina_pose_table.csv` | Phase 3 pose output in a Vina-compatible shape |
| `phase4_docking_evidence_reference.csv` | handoff file intended for Phase 4 |

### 4.4 Phase 4 perturbation outputs

Root:

- `output/phase4_perturbation/`

Representative files:

| Path | Meaning |
|------|------|
| `perturbation_candidate_table.csv` | final ranked candidate table |
| `perturbation_axis_scores.csv` | axis-level scoring breakdown |
| `final_candidate_classes.csv` | mechanistic classification labels |
| `integrated_phase4_report.md` | Phase 4 integrated narrative report |
| `phase4_final_review_table.csv` | review-first presentation table |
| `phase4_presentation_shortlist.csv` | shortlist for presentation or follow-up review |

## 5. Fast Interpretation Guide

Use this shortcut when onboarding.

| If you need to know... | Look here first |
|------|------|
| Which receptor and ligand inputs are active | `config/example-project.yaml` plus this document |
| What the current prepared Phase 1 input set is | `input/PPI/phase1/` |
| What older PPI prepared assets still exist | `input/PPI/prepared/` |
| Where the main Vina outputs live | `output/egfr_myo1d_vina/vina/` |
| Where the final verdict/report outputs live | `output/egfr_myo1d_vina/results/` |
| Where current structured Phase 1 outputs live | `output/phase1_ppi/` |
| Where future phase-separated pocket/docking/scoring outputs live | `output/phase2_pockets/`, `output/phase3_docking/`, `output/phase4_perturbation/` |

## 6. Current Cautions

- Do not confuse root-level pointer files inside `output/egfr_myo1d_vina/` with the actual payload files under `vina/`, `ppi/`, or `results/`.
- Do not treat `input/PPI/prepared/` as the current primary Phase 1 preparation layer when `input/PPI/phase1/` already provides the newer full-kinase-domain baseline metadata.
- Do not assume AFM inputs exist just because AFM parser code exists.
- Do not assume `.pdbqt` inputs are checked into this workspace just because the active YAML expects them at runtime.
