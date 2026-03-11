# Integrated Phase 4 Report: Perturbation Relevance Scoring

## 1. Executive Summary

- **Receptor states**: 1 (3GT8_raw)
- **Candidate pockets**: 3
- **Ligands evaluated**: 3 (173940, 97806, VAX-C12_0)
- **Total ranked entries**: 7

### Mechanistic Class Distribution

| Class | Count | Pockets |
|-------|-------|---------|
| ligandable_but_ppi_irrelevant_candidate | 1 | 3GT8_raw_PKT03 |
| orthosteric_disruptor_candidate | 6 | 3GT8_raw_PKT01, 3GT8_raw_PKT02 |

## 2. Why Affinity Alone Is Not the Final Criterion

This pipeline answers a specific biological question: *Which candidate sites are most plausible for disrupting MYO1D attachment to EGFR, and by what mechanism?*

A high-affinity docking score at a site with no connection to the MYO1D interface is chemically interesting but biologically irrelevant to the perturbation goal. Conversely, a moderate-affinity hit at a site that directly overlaps the PPI patch is a stronger perturbation candidate.

The 4-axis scoring framework prevents affinity domination by:

1. **Weighting PPI evidence (A1) and perturbation relevance (A3) at 60% combined** — biological relevance outweighs chemical metrics
2. **Capping affinity-influenced axes** for irrelevant sites — high druggability at a low-relevance site cannot outrank an orthosteric candidate
3. **Preserving mechanistic class labels** — the ranking is interpretable, not a black-box number

## 3. Scoring Framework

| Axis | Weight | Meaning |
|------|--------|---------|
| PPI Interface Confidence | 30% | Strong PPI patch overlap with validated hotspots |
| Druggability Confidence | 25% | High-confidence druggable pocket |
| Perturbation Relevance | 30% | Orthosteric site with ligand support |
| State Robustness / Accessibility | 15% | Robust pocket present in multiple conformational states |

Each axis score is [0, 1]. The perturbation score is their weighted sum, with an affinity cap applied to irrelevant/uncertain sites.

## 4. Ranked Candidates by Mechanistic Class

### 4.1. Orthosteric Disruptors

| Rank | Pocket | Ligand | Score | A1 | A2 | A3 | A4 | Confidence | Overlap |
|------|--------|--------|-------|----|----|----|----|------------|---------|
| 1 | 3GT8_raw_PKT01 | 173940 | 0.6514 | 0.393 | 0.947 | 0.723 | 0.533 | high_provisional | 4 |
| 2 | 3GT8_raw_PKT01 | 97806 | 0.6514 | 0.393 | 0.947 | 0.723 | 0.533 | high_provisional | 4 |
| 3 | 3GT8_raw_PKT01 | VAX-C12_0 | 0.6514 | 0.393 | 0.947 | 0.723 | 0.533 | high_provisional | 4 |
| 4 | 3GT8_raw_PKT02 | 173940 | 0.534 | 0.336 | 0.588 | 0.687 | 0.533 | high_provisional | 3 |
| 5 | 3GT8_raw_PKT02 | 97806 | 0.534 | 0.336 | 0.588 | 0.687 | 0.533 | high_provisional | 3 |
| 6 | 3GT8_raw_PKT02 | VAX-C12_0 | 0.534 | 0.336 | 0.588 | 0.687 | 0.533 | high_provisional | 3 |

### 4.4. Ligandable but PPI-Irrelevant

| Rank | Pocket | Ligand | Score | A1 | A2 | A3 | A4 | Confidence | Overlap |
|------|--------|--------|-------|----|----|----|----|------------|---------|
| 7 | 3GT8_raw_PKT03 | (none) | 0.1493 | 0.164 | 0.236 | 0.045 | 0.283 | high_provisional | 0 |

## 5. State Robustness Assessment

**3GT8_raw_PKT01**
- State class: state_specific_pocket
- Interpretation: state_dependent
- Accessibility: single_state_only
- States matched: 1
- Caveat: Only 1 receptor state has pocket data. Cannot determine if state-specific or under-sampled. NOTE: Single-state data. Cross-state persistence unknown. Re-evaluate when additional receptor states have pocket data.

**3GT8_raw_PKT02**
- State class: state_specific_pocket
- Interpretation: state_dependent
- Accessibility: single_state_only
- States matched: 1
- Caveat: Only 1 receptor state has pocket data. Cannot determine if state-specific or under-sampled. NOTE: Single-state data. Cross-state persistence unknown. Re-evaluate when additional receptor states have pocket data.

**3GT8_raw_PKT03**
- State class: state_specific_pocket
- Interpretation: state_dependent
- Accessibility: single_state_only
- States matched: 1
- Caveat: Only 1 receptor state has pocket data. Cannot determine if state-specific or under-sampled. NOTE: Single-state data. Cross-state persistence unknown. Re-evaluate when additional receptor states have pocket data.

## 6. Evidence Provenance Summary

Each candidate's score traces back to concrete upstream evidence:

| Phase | Data Source | Key Fields |
|-------|-----------|------------|
| Phase 1 | PPI patch reference | hotspot residues, robustness, method agreement |
| Phase 2 | Pocket proposal | relationship class, druggability tier, proposal score |
| Phase 3 | Docking evidence | ligand support, pose count, diversity verdict |
| Phase 4 | Scoring + classification | axis scores, mechanistic class, state interpretation |

Full provenance is available in `phase4_expanded_evidence_table.csv` (7 rows, 42 columns).

## 7. Validation and Caveats

### Data Completeness

- Phase 1 PPI evidence: Complete
- Phase 2 pocket evidence: Complete
- Phase 3 docking evidence: 6/7 candidates (partial)

### Known Limitations

1. **Single receptor state**: Only 3GT8_raw has pocket data. Scores will change when cl38_48 and cl85_100 fpocket/P2Rank results are integrated. State robustness axis (A4) currently has minimal discrimination power.

2. **Pre-execution affinity**: Vina docking jobs are prepared but not yet executed. All ligand support levels are provisional (`pending_*`). Re-run Phase 4 after server-side Vina execution.

3. **Single pocket detection tool**: Only fpocket results available. P2Rank integration will improve druggability axis (A2) through multi-tool consensus.

4. **No FTMap data**: Fragment hotspot mapping not yet available. When integrated, this will strengthen druggability confidence.

5. **Axis weights are initial**: The 30/25/30/15 weighting is a principled starting point but can be tuned after expert review.

### Recommended Next Steps

1. Execute Vina docking on HPC server (`run_phase3_docking.sh`)
2. Re-run Phase 3 TG 3.3-3.7 to populate affinity data
3. Re-run Phase 4 to update scores with real affinity evidence
4. Run fpocket/P2Rank on cl38_48 and cl85_100 receptor states
5. Re-run Phase 2-4 with multi-state data for robust state assessment

## 8. File Inventory

| File | Description |
|------|-------------|
| phase4_evidence_normalized.csv | Merged Phase 1-3 evidence (TG 4.0) |
| phase4_evidence_validation.md | Evidence validation report (TG 4.0) |
| phase4_axis_definition_table.csv | Axis definitions (TG 4.1) |
| phase4_axis_scores.csv | Raw axis scores (TG 4.1) |
| phase4_score_framework.md | Score framework documentation (TG 4.1) |
| final_candidate_classes.csv | Mechanistic classifications (TG 4.2) |
| phase4_mechanistic_classification_note.md | Classification logic (TG 4.2) |
| perturbation_axis_scores.csv | Weighted axis scores + rank (TG 4.3) |
| perturbation_candidate_table.csv | Ranked candidate table (TG 4.3) |
| phase4_ranking_method_note.md | Ranking methodology (TG 4.3) |
| phase4_state_interpretation.csv | State robustness interpretation (TG 4.4) |
| phase4_accessibility_note.md | Accessibility documentation (TG 4.4) |
| phase4_final_review_table.csv | Condensed review table (TG 4.5) |
| phase4_expanded_evidence_table.csv | Full provenance table (TG 4.5) |
| integrated_phase4_report.md | This report (TG 4.6) |

---

Generated by `egfr_pipeline.phase4.final_report`