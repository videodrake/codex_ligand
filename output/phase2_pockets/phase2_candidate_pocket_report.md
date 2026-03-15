# Phase 2 Review Report: Pocket Proposal and Druggability Mapping

## 1. Executive Summary

- **Receptor states with pocket data:** 1 of 3
  - Present: 3GT8_raw
  - Missing: EGFR_160-185, EGFR_170-200
- **Phase 1 patch hotspot residues:** 7
- **Raw pocket proposals:** 3
- **Merged candidate pockets:** 3
- **Phase 3 docking priorities:** primary=2, secondary=0, exploratory=0, skip=1

> **Note:** Cross-state alignment and state robustness are provisional.
> Re-run Phase 2 after fpocket/P2Rank results for all states are available.

## 2. Pocket Proposal Sources

| Tool | Total Proposals |
|------|----------------|
| fpocket | 3 |

### Per-State Breakdown

| State | Tool | Proposals |
|-------|------|-----------|
| 3GT8_raw | fpocket | 3 |

### Merge Summary

- Raw proposals: 3
- After merge: 3
- Reduction: 0 proposals merged

## 3. Merged Pocket Catalog

| Pocket | State | Sources | N Res | Volume (A3) | Proposal Score |
|--------|-------|---------|-------|-------------|---------------|
| 3GT8_raw_PKT01 | 3GT8_raw | fpocket | 4 | 892.3 | 0.8234 |
| 3GT8_raw_PKT02 | 3GT8_raw | fpocket | 3 | 567.9 | 0.5612 |
| 3GT8_raw_PKT03 | 3GT8_raw | fpocket | 2 | 312.4 | 0.3201 |

### Pocket Lining Residues

| Pocket | Residues |
|--------|----------|
| 3GT8_raw_PKT01 | LEU838;ASP855;ILE857;GLU866 |
| 3GT8_raw_PKT02 | LEU819;ALA822;VAL834 |
| 3GT8_raw_PKT03 | ARG748;GLY750 |

## 4. PPI Patch Relationship Classification

| Class | Count | Pockets |
|-------|-------|---------|
| orthosteric_candidate | 2 | 3GT8_raw_PKT01, 3GT8_raw_PKT02 |
| rim_candidate | 0 | --- |
| allosteric_candidate | 0 | --- |
| low_relevance_candidate | 1 | 3GT8_raw_PKT03 |

### Classification Detail

| Pocket | Class | Hotspot Overlap | Fraction | Centroid Dist (A) | Basis |
|--------|-------|----------------|----------|-------------------|-------|
| 3GT8_raw_PKT01 | orthosteric_candidate | 3 | 0.500 | 132.94 | hotspot_overlap=3 (frac=0.50 >= 0.25) |
| 3GT8_raw_PKT02 | orthosteric_candidate | 3 | 0.500 | 124.62 | hotspot_overlap=3 (frac=0.50 >= 0.25) |
| 3GT8_raw_PKT03 | low_relevance_candidate | 0 | 0.000 | 150.39 | centroid_distant=150.4A > 20.0A, no_overlap |

### Low-Relevance Pockets

1 pocket(s) classified as `low_relevance_candidate`.
These are structurally distant from the PPI patch and are not recommended
for PPI-disruption-oriented docking in Phase 3.

## 5. Druggability Confidence

### Druggability Confidence

| Confidence | Count | Pockets |
|-----------|-------|---------|
| high | 1 | 3GT8_raw_PKT01 |
| medium | 1 | 3GT8_raw_PKT02 |
| low | 1 | 3GT8_raw_PKT03 |

### Overall Druggability Tier

| Tier | Count | Pockets |
|------|-------|---------|
| tier_1 | 1 | 3GT8_raw_PKT01 |
| tier_2 | 1 | 3GT8_raw_PKT02 |
| tier_3 | 1 | 3GT8_raw_PKT03 |

### Multi-Source Support

- Consensus pockets (>=2 tools): 0
- Single-tool pockets: 3
- **Note:** All pockets are currently single-tool (fpocket only).
  Multi-source consensus will improve after P2Rank results are available.

### FTMap Hotspot Support

- FTMap data available: No
- FTMap is the preferred hotspot-support method for PPI-site
  ligandability evidence. Schema is ready for future integration.

## 6. Cross-State Pocket Alignment

| State Class | Count | Pockets |
|-------------|-------|---------|
| state_robust_pocket | 0 | --- |
| state_shifted_pocket | 0 | --- |
| state_specific_pocket | 3 | 3GT8_raw_PKT01, 3GT8_raw_PKT02, 3GT8_raw_PKT03 |
| uncertain_alignment | 0 | --- |

**2 state(s) missing pocket data:** EGFR_160-185, EGFR_170-200

Cross-state classification is provisional. All pockets are currently
`state_specific_pocket` because only one state has data.
Re-run TG 2.5 and TG 2.6 after all states have pocket proposals.

## 7. Phase 3 Docking Recommendation

### Docking Priority Summary

| Priority | Count | Pockets | Action |
|----------|-------|---------|--------|
| primary | 2 | 3GT8_raw_PKT01, 3GT8_raw_PKT02 | Full docking budget |
| secondary | 0 | --- | Reduced budget / diversity round |
| exploratory | 0 | --- | Only if resources allow |
| skip | 1 | 3GT8_raw_PKT03 | Do not dock |

### Comprehensive Pocket Summary

| Pocket | Relationship | Tier | Confidence | State Class | Priority |
|--------|-------------|------|-----------|-------------|----------|
| 3GT8_raw_PKT01 | orthosteric_candidate | tier_1 | high | state_specific_pocket | primary |
| 3GT8_raw_PKT02 | orthosteric_candidate | tier_2 | medium | state_specific_pocket | primary |
| 3GT8_raw_PKT03 | low_relevance_candidate | tier_3 | low | state_specific_pocket | skip |

### Primary Docking Targets — Detail

**3GT8_raw_PKT01** (3GT8_raw)
- Relationship: orthosteric_candidate
- Druggability: high (tier tier_1)
- Hotspot overlap: 3 residues (0.500 of patch)
- Centroid: (44.454, 32.261, 18.584)
- Box: 14.4 x 13.7 x 13.2 A
- Lining residues: LEU838;ASP855;ILE857;GLU866

**3GT8_raw_PKT02** (3GT8_raw)
- Relationship: orthosteric_candidate
- Druggability: medium (tier tier_2)
- Hotspot overlap: 3 residues (0.500 of patch)
- Centroid: (28.490, 45.493, 12.601)
- Box: 11.2 x 11.6 x 12.3 A
- Lining residues: LEU819;ALA822;VAL834

## 8. Limitations and Caveats

1. Only 1/3 receptor states have pocket proposals. Cross-state robustness cannot be assessed yet.
2. All pockets are single-tool (fpocket only). Multi-source consensus awaits P2Rank results.
3. FTMap hotspot support is not yet available. Druggability confidence relies solely on fpocket scores.
4. Druggability tier thresholds are calibrated for fpocket's 0-1 score range. When P2Rank data is added, normalized scores enable cross-tool comparison.

## 9. Phase 2 Output File Inventory

| File | TG | Description |
|------|-----|------------|
| phase2_patch_reference_normalized.csv | 2.0 | Phase 1 patch reference (normalized) |
| phase2_patch_reference_validation.md | 2.0 | Patch reference validation note |
| candidate_pockets_raw.csv | 2.1 | Raw pocket proposals (all tools) |
| candidate_pocket_source_summary.csv | 2.1 | Per-tool proposal counts |
| candidate_pockets.csv | 2.2 | Merged candidate pockets |
| candidate_pocket_merge_table.csv | 2.2 | Pairwise merge decisions |
| candidate_pocket_provenance.csv | 2.2 | Raw-to-merged mapping |
| pocket_patch_relationship.csv | 2.3 | Patch relationship classes |
| pocket_patch_relationship_metrics.csv | 2.3 | Classification raw metrics |
| druggability_proposal_summary.csv | 2.4 | Druggability confidence + tiers |
| candidate_pocket_support_flags.csv | 2.4 | Per-tool support flags |
| candidate_pocket_cross_state_comparison.csv | 2.5 | Pairwise cross-state comparisons |
| candidate_pocket_state_classes.csv | 2.5 | State robustness classes |
| phase3_candidate_pocket_reference.csv | 2.6 | Phase 3 docking reference |
| phase2_to_phase3_handoff_note.md | 2.6 | Phase 3 handoff note |
| phase2_candidate_pocket_report.md | 2.7 | This report |

---

Generated by `egfr_pipeline.phase2.review_report`