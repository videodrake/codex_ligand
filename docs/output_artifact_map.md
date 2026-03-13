# Output Artifact Map

This document maps the main output artifacts in the repository by phase and by operational role. It is designed for onboarding: a new GPT should be able to tell which files are presentation-facing, which files are machine handoff files, and which files are mostly trace or provenance records.

## Reading Rules

- Start with `output/{project}/step_index.md` when it exists. That is the first human-readable entry point for the step output layer.
- Treat `output/{project}/` as the canonical runtime output root for the current production lane.
- Treat `step1_vina_raw/` through `step7_validate/` as a derived interpretation view, not as a replacement for canonical outputs.
- Large raw PyRosetta directories are indexed by manifest and `raw_run_paths.tsv` rather than duplicated into the step layer.
- Start with [../output/README.md](../output/README.md) for output tree navigation.
- Treat `output/egfr_myo1d_vina/` as the current routine project output root for the Vina-centered baseline.
- Treat `output/phase1_ppi/`, `output/phase2_pockets/`, `output/phase3_docking/`, and `output/phase4_perturbation/` as the newer phase-separated output trees.
- Use [../output/phase1_ppi/README.md](../output/phase1_ppi/README.md) for Phase 1 report navigation.
- Use [../output/phase2_pockets/README.md](../output/phase2_pockets/README.md) for Phase 2 report navigation.
- Use [../output/phase3_docking/README.md](../output/phase3_docking/README.md) for Phase 3 report navigation.
- Use [../output/phase4_perturbation/README.md](../output/phase4_perturbation/README.md) for Phase 4 report navigation.
- Some older docs and legacy runs may still mention `vina/`, `ppi/`, or `results/` subdirectories. For the current production step layer, the canonical compact outputs of interest live directly under `output/{project}/`.

## Step Reading Order

Default interpretation order for a completed production run:

1. `output/{project}/step_index.md`
2. `output/{project}/step6_report/project_report.txt`
3. `output/{project}/step5_verdict/valid_sites.csv`
4. `output/{project}/step4_vina_postprocess/vina_pocket_table.csv`
5. `output/{project}/step3_ppi_postprocess/ppi_pyrosetta_residues.csv`

Canonical runtime outputs remain under the project root; the step folders are a derived interpretation view that can be regenerated from canonical outputs.

## Reference Output Layout

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

## Step Folder Map

| Step folder | Inspect first | Canonical source | Notes |
|------|------|------|------|
| `step1_vina_raw/` | `raw_pose_index.csv` | `output/{project}/{receptor_id}/{ligand}_{mode}.pdbqt` | Raw `.pdbqt` files are referenced by path only |
| `step2_ppi_raw/` | `TH1_final_ranking.csv`, `beta_meander_final_ranking.csv` | Target-specific PyRosetta run directories | Heavy raw directories are not mirrored |
| `step3_ppi_postprocess/` | `ppi_pyrosetta_residues.csv` | `output/{project}/ppi_pyrosetta_residues.csv` | Falls back to `output/{project}/ppi/` only when needed |
| `step4_vina_postprocess/` | `vina_pocket_table.csv` | `output/{project}/vina_pose_table.csv` and related canonical CSVs | Pocket-level interpretation layer |
| `step5_verdict/` | `valid_sites.csv` | `output/{project}/valid_sites.csv` and `cross_method_agreement.csv` | Verdict-first review layer |
| `step6_report/` | `project_report.txt` | `output/{project}/project_report.txt` and `combined_residue_evidence.csv` | Narrative summary-first layer |
| `step7_validate/` | `validation_status.json` | `egfr_pipeline.validate` result persisted at run end | Structured validation state |

## Priority Legend

| Label | Meaning |
|------|------|
| `P1` | First file to open for presentation, decision, or top-level interpretation |
| `P2` | Important supporting file; usually the next file to inspect after a `P1` artifact |
| `P3` | Traceability, validation, or provenance file; useful but not first-read |

## 1. Routine Baseline Artifacts

These are the artifacts most relevant to the current Vina-centered operational baseline.

| Area | Artifact | Canonical path | Meaning | Consumed next? | Priority |
|------|------|------|------|------|------|
| Vina poses | `vina_pose_table.csv` | `output/egfr_myo1d_vina/vina_pose_table.csv` | Pose-level table with receptor, ligand, pose rank, affinity, centroid, raw pose file path, pocket ID, and contact residues | Yes -> supports pocket summarization, comparison, report, and validation | `P2` |
| Vina pockets | `vina_pocket_table.csv` | `output/egfr_myo1d_vina/vina_pocket_table.csv` | Pocket-level aggregation of Vina poses within each receptor | Yes -> consumed by verdict, report, and validate | `P1` |
| Ligand-pocket map | `vina_drug_pocket_map.csv` | `output/egfr_myo1d_vina/vina_drug_pocket_map.csv` | Ligand-to-pocket summary showing dominant and multimodal mappings | Yes -> consumed by report and interpretation | `P2` |
| Cross-state pocket comparison | `vina_pocket_comparison.csv` | `output/egfr_myo1d_vina/vina_pocket_comparison.csv` | Cross-receptor pocket comparison table | Yes -> consumed by verdict, report, and validate | `P2` |
| Pocket bootstrap | `vina_pocket_bootstrap.csv` | `output/egfr_myo1d_vina/vina_pocket_bootstrap.csv` | Optional resampling-based support for pocket stability | Optional downstream support only | `P3` |
| PPI residue summary | `ppi_pyrosetta_residues.csv` | `output/egfr_myo1d_vina/ppi_pyrosetta_residues.csv` | Project-root PPI residue extraction table used by the current step layer | Yes -> contributes to verdict and report in the routine baseline | `P2` |
| PPI run summary | `ppi_pyrosetta_summary.csv` | `output/egfr_myo1d_vina/ppi_pyrosetta_summary.csv` | Aggregated summary of PyRosetta-derived PPI results | Supporting only | `P3` |
| Cross-method agreement | `cross_method_agreement.csv` | `output/egfr_myo1d_vina/cross_method_agreement.csv` | Main Vina-vs-PPI overlap and proximity table for the routine baseline | Yes -> used by report and validate; presentation-facing | `P1` |
| Site verdicts | `valid_sites.csv` | `output/egfr_myo1d_vina/valid_sites.csv` | Final rule-based site judgment table from `egfr_pipeline/verdict.py` | Yes -> used by report and validate; often the first baseline judgment table | `P1` |
| Consensus site map | `vina_consensus_sites.csv` | `output/egfr_myo1d_vina/vina_consensus_sites.csv` | Consensus pocket/site grouping across receptor states | No direct next step in baseline; review and presentation support | `P2` |
| Combined residue evidence | `combined_residue_evidence.csv` | `output/egfr_myo1d_vina/combined_residue_evidence.csv` | Residue-level integrated evidence view used for report support | No direct next step; interpretation support | `P2` |
| Top-level report | `project_report.txt` | `output/egfr_myo1d_vina/project_report.txt` | Human-readable summary report for the current baseline | No; presentation/report endpoint | `P1` |

## 2. Phase 1 Structured PPI Artifacts

These artifacts describe the newer structured Phase 1 interface-mapping branch. This branch is scientifically important even though it is not yet the only default operational entrypoint.

| Artifact | Canonical path | Meaning | Consumed next? | Priority |
|------|------|------|------|------|
| Phase 1 interface report | `output/phase1_ppi/phase1_interface_report.md` | Human-readable Phase 1 summary of prepared inputs, hotspot residues, orientation filtering, cross-state robustness, and LightDock support | No machine consumer; review and presentation | `P1` |
| Phase 1 patch handoff | `output/phase1_ppi/phase1_downstream_patch_reference.csv` | Structured receptor-side patch reference file for downstream phases; carries `orientation_validation_status` when orientation filtering is executed | Yes -> primary machine input for Phase 2 patch ingestion | `P1` |
| Cross-state comparison | `output/phase1_ppi/ppi_patch_cross_state_comparison.csv` | Residue-level comparison of Phase 1 patch presence across receptor states | Supporting Phase 1 review; not the direct Phase 2 handoff file | `P2` |
| State robustness | `output/phase1_ppi/ppi_patch_state_robustness.csv` | Robust/moderate/state-specific classification for Phase 1 residues | Supporting Phase 1 review and interpretation | `P2` |
| Interface patch table | `output/phase1_ppi/<state>/ppi_interface_patch_table.csv` | Per-state structured residue occupancy table for the inferred interface patch | No direct phase handoff; source-level support for review | `P2` |
| Cluster summary | `output/phase1_ppi/<state>/ppi_cluster_summary.csv` | Cluster-level summary of PyRosetta interface models | No direct phase handoff; analysis support | `P3` |
| Hotspot residues | `output/phase1_ppi/<state>/ppi_hotspot_residues.csv` | Hotspot residue table per state | Indirectly summarized into downstream patch reference | `P2` |
| Orientation filter log | `output/phase1_ppi/<state>/orientation_filter_log.csv` | PASS/FAIL/AMBIGUOUS orientation filter details for Phase 1 models | No; trace and QA support | `P3` |
| Cross-method convergence | `output/phase1_ppi/<state>/cross_method_convergence.csv` | Per-residue PyRosetta-vs-LightDock convergence table | No direct phase handoff; supports confidence interpretation | `P2` |
| LightDock support table | `output/phase1_ppi/<state>/lightdock/lightdock_interface_support_table.csv` | LightDock-derived residue support table | No; supporting evidence for Phase 1 review | `P3` |
| LightDock metadata | `output/phase1_ppi/<state>/lightdock/lightdock_run_metadata.json` | Run metadata for LightDock setup and execution | No; trace/provenance only | `P3` |

Recommended first-read files for Phase 1:

- `output/phase1_ppi/phase1_interface_report.md`
- `output/phase1_ppi/phase1_downstream_patch_reference.csv`
- `output/phase1_ppi/ppi_patch_state_robustness.csv`

## 3. Phase 2 Pocket Proposal Artifacts

These artifacts represent the planned pocket-proposal layer and its handoff to Phase 3.

| Artifact | Canonical path | Meaning | Consumed next? | Priority |
|------|------|------|------|------|
| Candidate pocket catalog | `output/phase2_pockets/candidate_pockets.csv` | Normalized receptor-local pocket catalog after merge logic | Yes -> supports relationship, druggability, and Phase 3 export | `P1` |
| Raw proposals | `output/phase2_pockets/candidate_pockets_raw.csv` | Tool-native raw pocket proposals before merge | No direct next step; provenance and debugging | `P3` |
| Merge table | `output/phase2_pockets/candidate_pocket_merge_table.csv` | How raw proposals were grouped into normalized pockets | No direct next step; provenance | `P3` |
| Provenance table | `output/phase2_pockets/candidate_pocket_provenance.csv` | Which source proposals contributed to each normalized pocket | No direct next step; provenance | `P3` |
| Patch relationship map | `output/phase2_pockets/pocket_patch_relationship.csv` | Pocket-level classification relative to the Phase 1 patch: orthosteric, rim, allosteric, or low-relevance | Yes -> consumed by Phase 3 export and Phase 4 evidence ingestion | `P1` |
| Relationship metrics | `output/phase2_pockets/pocket_patch_relationship_metrics.csv` | Raw metrics backing the relationship class | No direct next step; supporting QA | `P3` |
| Druggability summary | `output/phase2_pockets/druggability_proposal_summary.csv` | Phase 2 druggability tier/confidence summary for candidate pockets | Yes -> consumed by Phase 3 export and Phase 4 evidence ingestion | `P1` |
| Phase 2 patch reference normalized | `output/phase2_pockets/phase2_patch_reference_normalized.csv` | Internal normalized view of the Phase 1 patch reference | No direct next phase handoff; internal preparation artifact | `P3` |
| State classes | `output/phase2_pockets/candidate_pocket_state_classes.csv` | Cross-state classification for candidate pockets | Yes -> consumed by Phase 3 export and Phase 4 evidence ingestion | `P2` |
| Phase 3 pocket reference | `output/phase2_pockets/phase3_candidate_pocket_reference.csv` | Clean docking-ready pocket reference exported for Phase 3 | Yes -> primary machine input for Phase 3 pocket ingestion | `P1` |
| Phase 2 report | `output/phase2_pockets/phase2_candidate_pocket_report.md` | Human-readable summary of the Phase 2 pocket catalog | No; review and presentation support | `P2` |
| Phase 2 to Phase 3 note | `output/phase2_pockets/phase2_to_phase3_handoff_note.md` | Handoff explanation for downstream Phase 3 use | No machine consumer; handoff narrative | `P2` |

Recommended first-read files for Phase 2:

- `output/phase2_pockets/candidate_pockets.csv`
- `output/phase2_pockets/pocket_patch_relationship.csv`
- `output/phase2_pockets/druggability_proposal_summary.csv`
- `output/phase2_pockets/phase3_candidate_pocket_reference.csv`

## 4. Phase 3 Diverse Docking Artifacts

These artifacts represent the pocket-guided diversity-aware docking path, not the current routine blind-docking baseline.

| Artifact | Canonical path | Meaning | Consumed next? | Priority |
|------|------|------|------|------|
| Phase 3 normalized pocket reference | `output/phase3_docking/phase3_candidate_reference_normalized.csv` | Phase 3 internal normalized pocket reference with budget defaults and docking flags | Yes -> consumed by job construction and Phase 4 export support | `P2` |
| Docking job table | `output/phase3_docking/phase3_docking_job_table.csv` | Generated receptor-pocket-ligand job inventory | Yes -> execution planning only | `P3` |
| Docking box table | `output/phase3_docking/phase3_job_box_table.csv` | Per-pocket docking box definitions | Yes -> execution planning only | `P3` |
| Pocket search status | `output/phase3_docking/pocket_search_status.csv` | Status of each pocket under the saturation/budget policy | Yes -> used during Phase 3 execution and Phase 4 export context | `P1` |
| Budget tracking | `output/phase3_docking/phase3_budget_tracking.csv` | Per-pocket or per-round accounting of search budget usage | No direct next phase requirement; QA and interpretation support | `P2` |
| Round log | `output/phase3_docking/phase3_round_log.csv` | Round-level execution log for diversity-aware docking | No direct next phase requirement; trace support | `P3` |
| Phase 3 pose table | `output/phase3_docking/vina_pose_table.csv` | Phase 3 pose output in a Vina-compatible table shape | Yes -> contributes to Phase 4 docking evidence export | `P2` |
| Diversity report | `output/phase3_docking/phase3_diverse_docking_report.md` | Human-readable summary of the diverse docking run | No; presentation and review support | `P2` |
| Occupancy summary | `output/phase3_docking/phase3_pocket_occupancy_summary.csv` | Pocket occupancy summary across Phase 3 docking outputs | Yes -> consumed by Phase 4 export | `P2` |
| Diversity metrics | `output/phase3_docking/phase3_diversity_metrics.csv` | Diversity-support metrics for pockets and ligands | Yes -> consumed by Phase 4 export | `P2` |
| Phase 4 docking evidence reference | `output/phase3_docking/phase4_docking_evidence_reference.csv` | Clean docking evidence package exported for Phase 4 | Yes -> primary machine input for Phase 4 evidence ingestion | `P1` |
| Phase 3 to Phase 4 note | `output/phase3_docking/phase3_to_phase4_handoff_note.md` | Handoff explanation for Phase 4 consumers | No machine consumer; handoff narrative | `P2` |

Recommended first-read files for Phase 3:

- `output/phase3_docking/pocket_search_status.csv`
- `output/phase3_docking/phase3_diverse_docking_report.md`
- `output/phase3_docking/phase4_docking_evidence_reference.csv`

## 5. Phase 4 Perturbation Artifacts

These artifacts represent the advanced final-ranking stack. They are not the same as the routine `valid_sites.csv` verdict layer.

| Artifact | Canonical path | Meaning | Consumed next? | Priority |
|------|------|------|------|------|
| Perturbation candidate table | `output/phase4_perturbation/perturbation_candidate_table.csv` | Final ranked pocket-ligand candidate table with mechanistic class and score axes | Yes -> feeds review output and final report builders | `P1` |
| Axis score table | `output/phase4_perturbation/perturbation_axis_scores.csv` | Detailed axis-by-axis scoring contributions | Yes -> supports review and interpretation outputs | `P2` |
| Candidate classes | `output/phase4_perturbation/final_candidate_classes.csv` | Mechanistic classification output | Yes -> supports review and report builders | `P2` |
| Integrated Phase 4 report | `output/phase4_perturbation/integrated_phase4_report.md` | Human-readable final Phase 4 narrative report | No; presentation/report endpoint | `P1` |
| Final review table | `output/phase4_perturbation/phase4_final_review_table.csv` | Condensed review-first table for top candidates | No machine consumer; best quick review table | `P1` |
| Expanded evidence table | `output/phase4_perturbation/phase4_expanded_evidence_table.csv` | Full provenance table across Phases 1-4 | No; detailed review and traceability | `P2` |
| Evidence normalized | `output/phase4_perturbation/phase4_evidence_normalized.csv` | Internal normalized multi-phase evidence table | Yes -> supports scoring and classification inside Phase 4 | `P3` |
| Evidence validation | `output/phase4_perturbation/phase4_evidence_validation.md` | Validation and consistency checks for upstream evidence ingestion | No; QA/validation support | `P3` |
| State interpretation | `output/phase4_perturbation/phase4_state_interpretation.csv` | Cross-state interpretation layer for final ranking | Yes -> supports review output and final report | `P2` |
| Presentation shortlist | `output/phase4_perturbation/phase4_presentation_shortlist.csv` | Shortlist intended for presentation or follow-up review | No; presentation support | `P2` |

Recommended first-read files for Phase 4:

- `output/phase4_perturbation/phase4_final_review_table.csv`
- `output/phase4_perturbation/perturbation_candidate_table.csv`
- `output/phase4_perturbation/integrated_phase4_report.md`

## 6. Quick Open Order By Task

Use these shortcuts when a new GPT needs to answer a narrow question quickly.

| If the question is... | Open these artifacts first |
|------|------|
| What does the current routine baseline think are the important sites? | `output/egfr_myo1d_vina/step_index.md`, `output/egfr_myo1d_vina/step5_verdict/valid_sites.csv`, `output/egfr_myo1d_vina/step6_report/project_report.txt` |
| What is the current Phase 1 patch that later phases should trust? | `output/phase1_ppi/phase1_downstream_patch_reference.csv`, `output/phase1_ppi/phase1_interface_report.md` |
| Which pockets matter relative to the Phase 1 patch? | `output/phase2_pockets/pocket_patch_relationship.csv`, `output/phase2_pockets/candidate_pockets.csv`, `output/phase2_pockets/druggability_proposal_summary.csv` |
| What does the Phase 3 experimental docking path actually hand to Phase 4? | `output/phase3_docking/phase4_docking_evidence_reference.csv`, `output/phase3_docking/pocket_search_status.csv` |
| What does the advanced Phase 4 ranking think is most relevant biologically? | `output/phase4_perturbation/phase4_final_review_table.csv`, `output/phase4_perturbation/perturbation_candidate_table.csv`, `output/phase4_perturbation/integrated_phase4_report.md` |

## 7. Current Cautions

- Do not confuse the routine baseline judgment files (`valid_sites.csv`, `project_report.txt`) with the advanced Phase 4 perturbation-ranking outputs.
- Do not treat all files under `output/egfr_myo1d_vina/` as canonical payloads without checking whether they are pointer stubs.
- Treat `orientation_validation_status` as source-dependent: structured `output/phase1_ppi/` runs can carry calibrated classes, while legacy `output/egfr_myo1d_vina/ppi/` postprocess tables often use `not_available`.
- Do not assume that a phase-separated artifact is automatically consumed by the default CLI or production flow just because the file exists.


