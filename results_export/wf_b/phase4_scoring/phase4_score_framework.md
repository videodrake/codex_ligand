# Phase 4 Score Framework

## 1. Design Principles

This framework scores candidate sites on **perturbation relevance** rather than simple ligand affinity. The goal is to answer: *Which sites are most plausible for disrupting MYO1D attachment to EGFR, and by what mechanism?*

Key design choices:
- **4 independent axes** prevent any single metric from dominating
- **Raw submetrics preserved** for reviewability
- **Affinity is one input, not the score** — high affinity at an irrelevant site does not outrank moderate affinity at an orthosteric site
- **State-specific sites are flagged, not eliminated**

## 2. Axis Definitions

### A1_ppi_interface: PPI Interface Confidence (weight: 30%)

**Interpretation**: How strong and reliable the receptor-side PPI patch evidence is. Reflects hotspot residue quality, multi-method agreement, and robustness across receptor states.

- Range: 0.0–1.0
- High score means: Strong PPI patch overlap with validated hotspots
- Submetrics: hotspot_overlap_fraction, mean_hotspot_confidence, n_robust_hotspots, n_method_agreement_both

### A2_druggability: Druggability Confidence (weight: 25%)

**Interpretation**: How plausible the pocket is for small-molecule engagement. Reflects fpocket/P2Rank scoring, pocket volume, and multi-tool consensus.

- Range: 0.0–1.0
- High score means: High-confidence druggable pocket
- Submetrics: overall_druggability_tier, druggability_confidence, best_proposal_score

### A3_perturbation_relevance: Perturbation Relevance (weight: 30%)

**Interpretation**: How directly the site can plausibly alter MYO1D attachment. Combines PPI relationship class with docking evidence to assess mechanistic potential.

- Range: 0.0–1.0
- High score means: Orthosteric site with ligand support
- Submetrics: relationship_class, hotspot_overlap_count, ligand_support_strength, pose_support_count

### A4_state_robustness: State Robustness / Accessibility (weight: 15%)

**Interpretation**: Whether the site is persistent across receptor states, conditionally accessible, or strongly state-dependent. State-specific sites are not penalized but flagged.

- Range: 0.0–1.0
- High score means: Robust pocket present in multiple conformational states
- Submetrics: state_class, n_states_matched, diversity_verdict, coverage_fraction

**Total weight**: 100%

## 3. Scoring Rules

Each axis is computed as a weighted sum of its submetrics, each normalized to [0, 1]. Categorical variables are mapped to ordinal scores (e.g., tier_1 = 1.0, tier_2 = 0.6, tier_3 = 0.2). The final perturbation score (TG 4.3) will combine axis scores using the defined weights.

### Categorical Mappings

| Variable | Value | Score |
|----------|-------|-------|
| druggability_tier | tier_1 | 1.0 |
| druggability_tier | tier_2 | 0.6 |
| druggability_tier | tier_3 | 0.2 |
| druggability_confidence | high | 1.0 |
| druggability_confidence | low | 0.2 |
| druggability_confidence | medium | 0.6 |
| relationship_class | allosteric_candidate | 0.5 |
| relationship_class | low_relevance_candidate | 0.1 |
| relationship_class | orthosteric_candidate | 1.0 |
| relationship_class | rim_candidate | 0.7 |
| ligand_support | moderate | 0.6 |
| ligand_support | none | 0.0 |
| ligand_support | pending_moderate | 0.5 |
| ligand_support | pending_strong | 0.8 |
| ligand_support | pending_weak | 0.2 |
| ligand_support | strong | 1.0 |
| ligand_support | weak | 0.3 |
| state_class | robust_pocket | 1.0 |
| state_class | shifted_pocket | 0.7 |
| state_class | state_specific_pocket | 0.4 |

## 4. Current Axis Scores

| Pocket | Ligand | A1 PPI | A2 Drug | A3 Perturb | A4 State |
|--------|--------|--------|---------|------------|----------|
| 3GT8_raw_PKT01 | 173940 | 0.280 | 0.157 | 0.630 | 0.650 |
| 3GT8_raw_PKT02 | 173940 | 0.280 | 0.151 | 0.620 | 0.650 |
| 3GT8_raw_PKT03 | 173940 | 0.253 | 0.149 | 0.475 | 0.650 |
| 3GT8_raw_PKT05 | 173940 | 0.306 | 0.135 | 0.644 | 0.650 |
| 3GT8_raw_PKT07 | 173940 | 0.280 | 0.683 | 0.630 | 0.650 |
| 3GT8_raw_PKT08 | 173940 | 0.306 | 0.123 | 0.644 | 0.567 |
| 3GT8_raw_PKT10 | 173940 | 0.253 | 0.399 | 0.525 | 0.650 |
| 3GT8_raw_PKT11 | 173940 | 0.306 | 0.118 | 0.644 | 0.650 |
| 3GT8_raw_PKT14 | 173940 | 0.280 | 0.115 | 0.630 | 0.650 |
| 3GT8_raw_PKT16 | 173940 | 0.253 | 0.114 | 0.475 | 0.650 |
| 3GT8_raw_PKT20 | 173940 | 0.306 | 0.097 | 0.644 | 0.533 |
| 3GT8_raw_PKT22 | 173940 | 0.253 | 0.024 | 0.475 | 0.567 |
| 3GT8_raw_PKT01 | 97806 | 0.280 | 0.157 | 0.630 | 0.650 |
| 3GT8_raw_PKT02 | 97806 | 0.280 | 0.151 | 0.630 | 0.650 |
| 3GT8_raw_PKT03 | 97806 | 0.253 | 0.149 | 0.475 | 0.650 |
| 3GT8_raw_PKT05 | 97806 | 0.306 | 0.135 | 0.644 | 0.650 |
| 3GT8_raw_PKT07 | 97806 | 0.280 | 0.683 | 0.630 | 0.650 |
| 3GT8_raw_PKT08 | 97806 | 0.306 | 0.123 | 0.564 | 0.567 |
| 3GT8_raw_PKT10 | 97806 | 0.253 | 0.399 | 0.525 | 0.650 |
| 3GT8_raw_PKT11 | 97806 | 0.306 | 0.118 | 0.644 | 0.650 |
| 3GT8_raw_PKT14 | 97806 | 0.280 | 0.115 | 0.630 | 0.650 |
| 3GT8_raw_PKT16 | 97806 | 0.253 | 0.114 | 0.475 | 0.650 |
| 3GT8_raw_PKT20 | 97806 | 0.306 | 0.097 | 0.564 | 0.533 |
| 3GT8_raw_PKT22 | 97806 | 0.253 | 0.024 | 0.475 | 0.567 |
| 3GT8_raw_PKT01 | VAX-C12_0 | 0.280 | 0.157 | 0.630 | 0.650 |
| 3GT8_raw_PKT02 | VAX-C12_0 | 0.280 | 0.151 | 0.630 | 0.650 |
| 3GT8_raw_PKT03 | VAX-C12_0 | 0.253 | 0.149 | 0.475 | 0.650 |
| 3GT8_raw_PKT05 | VAX-C12_0 | 0.306 | 0.135 | 0.644 | 0.650 |
| 3GT8_raw_PKT07 | VAX-C12_0 | 0.280 | 0.683 | 0.630 | 0.650 |
| 3GT8_raw_PKT08 | VAX-C12_0 | 0.306 | 0.123 | 0.564 | 0.567 |
| 3GT8_raw_PKT10 | VAX-C12_0 | 0.253 | 0.399 | 0.525 | 0.650 |
| 3GT8_raw_PKT11 | VAX-C12_0 | 0.306 | 0.118 | 0.644 | 0.650 |
| 3GT8_raw_PKT14 | VAX-C12_0 | 0.280 | 0.115 | 0.630 | 0.650 |
| 3GT8_raw_PKT16 | VAX-C12_0 | 0.253 | 0.114 | 0.475 | 0.650 |
| 3GT8_raw_PKT20 | VAX-C12_0 | 0.306 | 0.097 | 0.644 | 0.533 |
| 3GT8_raw_PKT22 | VAX-C12_0 | 0.253 | 0.024 | 0.475 | 0.567 |
| EGFR_160-185_PKT01 | 173940 | 0.253 | 0.204 | 0.475 | 0.650 |
| EGFR_160-185_PKT02 | 173940 | 0.310 | 0.195 | 0.644 | 0.650 |
| EGFR_160-185_PKT03 | 173940 | 0.253 | 0.180 | 0.475 | 0.650 |
| EGFR_160-185_PKT04 | 173940 | 0.253 | 0.167 | 0.475 | 0.650 |
| EGFR_160-185_PKT05 | 173940 | 0.253 | 0.161 | 0.475 | 0.650 |
| EGFR_160-185_PKT07 | 173940 | 0.281 | 0.153 | 0.620 | 0.650 |
| EGFR_160-185_PKT09 | 173940 | 0.281 | 0.151 | 0.630 | 0.567 |
| EGFR_160-185_PKT11 | 173940 | 0.253 | 0.144 | 0.395 | 0.650 |
| EGFR_160-185_PKT12 | 173940 | 0.281 | 0.143 | 0.630 | 0.650 |
| EGFR_160-185_PKT13 | 173940 | 0.253 | 0.142 | 0.475 | 0.650 |
| EGFR_160-185_PKT14 | 173940 | 0.253 | 0.141 | 0.475 | 0.567 |
| EGFR_160-185_PKT15 | 173940 | 0.281 | 0.141 | 0.630 | 0.650 |
| EGFR_160-185_PKT16 | 173940 | 0.310 | 0.141 | 0.644 | 0.650 |
| EGFR_160-185_PKT19 | 173940 | 0.253 | 0.135 | 0.475 | 0.650 |
| EGFR_160-185_PKT23 | 173940 | 0.281 | 0.129 | 0.630 | 0.650 |
| EGFR_160-185_PKT27 | 173940 | 0.253 | 0.121 | 0.475 | 0.567 |
| EGFR_160-185_PKT29 | 173940 | 0.253 | 0.117 | 0.475 | 0.567 |
| EGFR_160-185_PKT30 | 173940 | 0.253 | 0.114 | 0.475 | 0.650 |
| EGFR_160-185_PKT34 | 173940 | 0.281 | 0.106 | 0.630 | 0.567 |
| EGFR_160-185_PKT37 | 173940 | 0.253 | 0.104 | 0.475 | 0.650 |
| EGFR_160-185_PKT38 | 173940 | 0.253 | 0.104 | 0.475 | 0.650 |
| EGFR_160-185_PKT39 | 173940 | 0.253 | 0.104 | 0.475 | 0.567 |
| EGFR_160-185_PKT40 | 173940 | 0.253 | 0.094 | 0.475 | 0.567 |
| EGFR_160-185_PKT43 | 173940 | 0.253 | 0.088 | 0.475 | 0.567 |
| EGFR_160-185_PKT01 | 97806 | 0.253 | 0.204 | 0.475 | 0.650 |
| EGFR_160-185_PKT02 | 97806 | 0.310 | 0.195 | 0.644 | 0.650 |
| EGFR_160-185_PKT03 | 97806 | 0.253 | 0.180 | 0.475 | 0.650 |
| EGFR_160-185_PKT04 | 97806 | 0.253 | 0.167 | 0.475 | 0.650 |
| EGFR_160-185_PKT05 | 97806 | 0.253 | 0.161 | 0.475 | 0.650 |
| EGFR_160-185_PKT07 | 97806 | 0.281 | 0.153 | 0.630 | 0.650 |
| EGFR_160-185_PKT09 | 97806 | 0.281 | 0.151 | 0.550 | 0.567 |
| EGFR_160-185_PKT11 | 97806 | 0.253 | 0.144 | 0.395 | 0.650 |
| EGFR_160-185_PKT12 | 97806 | 0.281 | 0.143 | 0.630 | 0.650 |
| EGFR_160-185_PKT13 | 97806 | 0.253 | 0.142 | 0.385 | 0.650 |
| EGFR_160-185_PKT14 | 97806 | 0.253 | 0.141 | 0.475 | 0.567 |
| EGFR_160-185_PKT15 | 97806 | 0.281 | 0.141 | 0.630 | 0.650 |
| EGFR_160-185_PKT16 | 97806 | 0.310 | 0.141 | 0.644 | 0.650 |
| EGFR_160-185_PKT19 | 97806 | 0.253 | 0.135 | 0.475 | 0.650 |
| EGFR_160-185_PKT23 | 97806 | 0.281 | 0.129 | 0.630 | 0.650 |
| EGFR_160-185_PKT27 | 97806 | 0.253 | 0.121 | 0.475 | 0.567 |
| EGFR_160-185_PKT29 | 97806 | 0.253 | 0.117 | 0.475 | 0.567 |
| EGFR_160-185_PKT30 | 97806 | 0.253 | 0.114 | 0.475 | 0.650 |
| EGFR_160-185_PKT34 | 97806 | 0.281 | 0.106 | 0.550 | 0.567 |
| EGFR_160-185_PKT37 | 97806 | 0.253 | 0.104 | 0.475 | 0.650 |
| EGFR_160-185_PKT38 | 97806 | 0.253 | 0.104 | 0.475 | 0.650 |
| EGFR_160-185_PKT39 | 97806 | 0.253 | 0.104 | 0.475 | 0.567 |
| EGFR_160-185_PKT40 | 97806 | 0.253 | 0.094 | 0.385 | 0.567 |
| EGFR_160-185_PKT43 | 97806 | 0.253 | 0.088 | 0.475 | 0.567 |
| EGFR_160-185_PKT01 | VAX-C12_0 | 0.253 | 0.204 | 0.475 | 0.650 |
| EGFR_160-185_PKT02 | VAX-C12_0 | 0.310 | 0.195 | 0.644 | 0.650 |
| EGFR_160-185_PKT03 | VAX-C12_0 | 0.253 | 0.180 | 0.475 | 0.650 |
| EGFR_160-185_PKT04 | VAX-C12_0 | 0.253 | 0.167 | 0.475 | 0.650 |
| EGFR_160-185_PKT05 | VAX-C12_0 | 0.253 | 0.161 | 0.475 | 0.650 |
| EGFR_160-185_PKT07 | VAX-C12_0 | 0.281 | 0.153 | 0.630 | 0.650 |
| EGFR_160-185_PKT09 | VAX-C12_0 | 0.281 | 0.151 | 0.630 | 0.567 |
| EGFR_160-185_PKT11 | VAX-C12_0 | 0.253 | 0.144 | 0.395 | 0.650 |
| EGFR_160-185_PKT12 | VAX-C12_0 | 0.281 | 0.143 | 0.630 | 0.650 |
| EGFR_160-185_PKT13 | VAX-C12_0 | 0.253 | 0.142 | 0.475 | 0.650 |
| EGFR_160-185_PKT14 | VAX-C12_0 | 0.253 | 0.141 | 0.475 | 0.567 |
| EGFR_160-185_PKT15 | VAX-C12_0 | 0.281 | 0.141 | 0.630 | 0.650 |
| EGFR_160-185_PKT16 | VAX-C12_0 | 0.310 | 0.141 | 0.644 | 0.650 |
| EGFR_160-185_PKT19 | VAX-C12_0 | 0.253 | 0.135 | 0.475 | 0.650 |
| EGFR_160-185_PKT23 | VAX-C12_0 | 0.281 | 0.129 | 0.630 | 0.650 |
| EGFR_160-185_PKT27 | VAX-C12_0 | 0.253 | 0.121 | 0.475 | 0.567 |
| EGFR_160-185_PKT29 | VAX-C12_0 | 0.253 | 0.117 | 0.475 | 0.567 |
| EGFR_160-185_PKT30 | VAX-C12_0 | 0.253 | 0.114 | 0.475 | 0.650 |
| EGFR_160-185_PKT34 | VAX-C12_0 | 0.281 | 0.106 | 0.630 | 0.567 |
| EGFR_160-185_PKT37 | VAX-C12_0 | 0.253 | 0.104 | 0.475 | 0.650 |
| EGFR_160-185_PKT38 | VAX-C12_0 | 0.253 | 0.104 | 0.475 | 0.650 |
| EGFR_160-185_PKT39 | VAX-C12_0 | 0.253 | 0.104 | 0.475 | 0.567 |
| EGFR_160-185_PKT40 | VAX-C12_0 | 0.253 | 0.094 | 0.475 | 0.567 |
| EGFR_160-185_PKT43 | VAX-C12_0 | 0.253 | 0.088 | 0.475 | 0.567 |
| EGFR_170-200_PKT01 | 173940 | 0.253 | 0.218 | 0.475 | 0.650 |
| EGFR_170-200_PKT04 | 173940 | 0.253 | 0.161 | 0.295 | 0.650 |
| EGFR_170-200_PKT05 | 173940 | 0.253 | 0.157 | 0.475 | 0.650 |
| EGFR_170-200_PKT06 | 173940 | 0.310 | 0.153 | 0.644 | 0.650 |
| EGFR_170-200_PKT10 | 173940 | 0.281 | 0.142 | 0.630 | 0.650 |
| EGFR_170-200_PKT12 | 173940 | 0.281 | 0.139 | 0.620 | 0.650 |
| EGFR_170-200_PKT13 | 173940 | 0.253 | 0.139 | 0.475 | 0.567 |
| EGFR_170-200_PKT15 | 173940 | 0.253 | 0.134 | 0.475 | 0.650 |
| EGFR_170-200_PKT17 | 173940 | 0.339 | 0.133 | 0.659 | 0.650 |
| EGFR_170-200_PKT25 | 173940 | 0.253 | 0.116 | 0.475 | 0.650 |
| EGFR_170-200_PKT26 | 173940 | 0.253 | 0.116 | 0.395 | 0.567 |
| EGFR_170-200_PKT28 | 173940 | 0.310 | 0.115 | 0.564 | 0.650 |
| EGFR_170-200_PKT31 | 173940 | 0.253 | 0.100 | 0.475 | 0.650 |
| EGFR_170-200_PKT33 | 173940 | 0.253 | 0.090 | 0.475 | 0.567 |
| EGFR_170-200_PKT34 | 173940 | 0.253 | 0.645 | 0.525 | 0.650 |
| EGFR_170-200_PKT36 | 173940 | 0.253 | 0.044 | 0.475 | 0.650 |
| EGFR_170-200_PKT01 | 97806 | 0.253 | 0.218 | 0.475 | 0.650 |
| EGFR_170-200_PKT04 | 97806 | 0.253 | 0.161 | 0.295 | 0.650 |
| EGFR_170-200_PKT05 | 97806 | 0.253 | 0.157 | 0.475 | 0.650 |
| EGFR_170-200_PKT06 | 97806 | 0.310 | 0.153 | 0.644 | 0.650 |
| EGFR_170-200_PKT10 | 97806 | 0.281 | 0.142 | 0.630 | 0.650 |
| EGFR_170-200_PKT12 | 97806 | 0.281 | 0.139 | 0.630 | 0.650 |
| EGFR_170-200_PKT13 | 97806 | 0.253 | 0.139 | 0.475 | 0.567 |
| EGFR_170-200_PKT15 | 97806 | 0.253 | 0.134 | 0.475 | 0.650 |
| EGFR_170-200_PKT17 | 97806 | 0.339 | 0.133 | 0.579 | 0.650 |
| EGFR_170-200_PKT25 | 97806 | 0.253 | 0.116 | 0.475 | 0.650 |
| EGFR_170-200_PKT26 | 97806 | 0.253 | 0.116 | 0.395 | 0.567 |
| EGFR_170-200_PKT28 | 97806 | 0.310 | 0.115 | 0.504 | 0.650 |
| EGFR_170-200_PKT31 | 97806 | 0.253 | 0.100 | 0.475 | 0.650 |
| EGFR_170-200_PKT33 | 97806 | 0.253 | 0.090 | 0.475 | 0.567 |
| EGFR_170-200_PKT34 | 97806 | 0.253 | 0.645 | 0.525 | 0.650 |
| EGFR_170-200_PKT36 | 97806 | 0.253 | 0.044 | 0.475 | 0.650 |
| EGFR_170-200_PKT01 | VAX-C12_0 | 0.253 | 0.218 | 0.475 | 0.650 |
| EGFR_170-200_PKT04 | VAX-C12_0 | 0.253 | 0.161 | 0.305 | 0.650 |
| EGFR_170-200_PKT05 | VAX-C12_0 | 0.253 | 0.157 | 0.475 | 0.650 |
| EGFR_170-200_PKT06 | VAX-C12_0 | 0.310 | 0.153 | 0.644 | 0.650 |
| EGFR_170-200_PKT10 | VAX-C12_0 | 0.281 | 0.142 | 0.630 | 0.650 |
| EGFR_170-200_PKT12 | VAX-C12_0 | 0.281 | 0.139 | 0.620 | 0.650 |
| EGFR_170-200_PKT13 | VAX-C12_0 | 0.253 | 0.139 | 0.475 | 0.567 |
| EGFR_170-200_PKT15 | VAX-C12_0 | 0.253 | 0.134 | 0.475 | 0.650 |
| EGFR_170-200_PKT17 | VAX-C12_0 | 0.339 | 0.133 | 0.659 | 0.650 |
| EGFR_170-200_PKT25 | VAX-C12_0 | 0.253 | 0.116 | 0.475 | 0.650 |
| EGFR_170-200_PKT26 | VAX-C12_0 | 0.253 | 0.116 | 0.475 | 0.567 |
| EGFR_170-200_PKT28 | VAX-C12_0 | 0.310 | 0.115 | 0.564 | 0.650 |
| EGFR_170-200_PKT31 | VAX-C12_0 | 0.253 | 0.100 | 0.475 | 0.650 |
| EGFR_170-200_PKT33 | VAX-C12_0 | 0.253 | 0.090 | 0.475 | 0.567 |
| EGFR_170-200_PKT34 | VAX-C12_0 | 0.253 | 0.645 | 0.525 | 0.650 |
| EGFR_170-200_PKT36 | VAX-C12_0 | 0.253 | 0.044 | 0.475 | 0.650 |
| 3GT8_raw_PKT04 | (none) | 0.253 | 0.146 | 0.045 | 0.317 |
| 3GT8_raw_PKT06 | (none) | 0.253 | 0.133 | 0.045 | 0.283 |
| 3GT8_raw_PKT09 | (none) | 0.253 | 0.121 | 0.045 | 0.317 |
| 3GT8_raw_PKT12 | (none) | 0.253 | 0.117 | 0.045 | 0.317 |
| 3GT8_raw_PKT13 | (none) | 0.253 | 0.117 | 0.045 | 0.400 |
| 3GT8_raw_PKT15 | (none) | 0.253 | 0.115 | 0.045 | 0.283 |
| 3GT8_raw_PKT17 | (none) | 0.253 | 0.108 | 0.045 | 0.317 |
| 3GT8_raw_PKT18 | (none) | 0.253 | 0.106 | 0.045 | 0.400 |
| 3GT8_raw_PKT19 | (none) | 0.253 | 0.102 | 0.045 | 0.317 |
| 3GT8_raw_PKT21 | (none) | 0.253 | 0.093 | 0.045 | 0.317 |
| EGFR_160-185_PKT06 | (none) | 0.253 | 0.153 | 0.045 | 0.400 |
| EGFR_160-185_PKT08 | (none) | 0.253 | 0.152 | 0.045 | 0.317 |
| EGFR_160-185_PKT10 | (none) | 0.253 | 0.151 | 0.045 | 0.317 |
| EGFR_160-185_PKT17 | (none) | 0.253 | 0.137 | 0.045 | 0.400 |
| EGFR_160-185_PKT18 | (none) | 0.253 | 0.136 | 0.045 | 0.317 |
| EGFR_160-185_PKT20 | (none) | 0.253 | 0.133 | 0.045 | 0.317 |
| EGFR_160-185_PKT21 | (none) | 0.253 | 0.133 | 0.045 | 0.317 |
| EGFR_160-185_PKT22 | (none) | 0.253 | 0.131 | 0.045 | 0.317 |
| EGFR_160-185_PKT24 | (none) | 0.253 | 0.128 | 0.045 | 0.317 |
| EGFR_160-185_PKT25 | (none) | 0.253 | 0.128 | 0.045 | 0.317 |
| EGFR_160-185_PKT26 | (none) | 0.253 | 0.123 | 0.045 | 0.317 |
| EGFR_160-185_PKT28 | (none) | 0.253 | 0.119 | 0.045 | 0.317 |
| EGFR_160-185_PKT31 | (none) | 0.253 | 0.113 | 0.045 | 0.317 |
| EGFR_160-185_PKT32 | (none) | 0.253 | 0.113 | 0.045 | 0.317 |
| EGFR_160-185_PKT33 | (none) | 0.253 | 0.112 | 0.045 | 0.317 |
| EGFR_160-185_PKT35 | (none) | 0.253 | 0.106 | 0.045 | 0.317 |
| EGFR_160-185_PKT36 | (none) | 0.253 | 0.106 | 0.045 | 0.317 |
| EGFR_160-185_PKT41 | (none) | 0.253 | 0.093 | 0.045 | 0.283 |
| EGFR_160-185_PKT42 | (none) | 0.253 | 0.088 | 0.045 | 0.317 |
| EGFR_160-185_PKT44 | (none) | 0.253 | 0.085 | 0.045 | 0.317 |
| EGFR_170-200_PKT02 | (none) | 0.253 | 0.172 | 0.045 | 0.283 |
| EGFR_170-200_PKT03 | (none) | 0.253 | 0.163 | 0.045 | 0.400 |
| EGFR_170-200_PKT07 | (none) | 0.253 | 0.150 | 0.045 | 0.400 |
| EGFR_170-200_PKT08 | (none) | 0.253 | 0.146 | 0.045 | 0.317 |
| EGFR_170-200_PKT09 | (none) | 0.253 | 0.144 | 0.045 | 0.317 |
| EGFR_170-200_PKT11 | (none) | 0.253 | 0.140 | 0.045 | 0.317 |
| EGFR_170-200_PKT14 | (none) | 0.253 | 0.137 | 0.045 | 0.317 |
| EGFR_170-200_PKT16 | (none) | 0.253 | 0.133 | 0.045 | 0.317 |
| EGFR_170-200_PKT18 | (none) | 0.253 | 0.130 | 0.045 | 0.400 |
| EGFR_170-200_PKT19 | (none) | 0.253 | 0.126 | 0.045 | 0.317 |
| EGFR_170-200_PKT20 | (none) | 0.253 | 0.122 | 0.045 | 0.317 |
| EGFR_170-200_PKT21 | (none) | 0.253 | 0.121 | 0.045 | 0.317 |
| EGFR_170-200_PKT22 | (none) | 0.253 | 0.121 | 0.045 | 0.317 |
| EGFR_170-200_PKT23 | (none) | 0.253 | 0.121 | 0.045 | 0.317 |
| EGFR_170-200_PKT24 | (none) | 0.253 | 0.120 | 0.045 | 0.317 |
| EGFR_170-200_PKT27 | (none) | 0.253 | 0.115 | 0.045 | 0.400 |
| EGFR_170-200_PKT29 | (none) | 0.253 | 0.112 | 0.045 | 0.400 |
| EGFR_170-200_PKT30 | (none) | 0.253 | 0.102 | 0.045 | 0.317 |
| EGFR_170-200_PKT32 | (none) | 0.253 | 0.093 | 0.045 | 0.317 |
| EGFR_170-200_PKT35 | (none) | 0.253 | 0.056 | 0.045 | 0.317 |
| EGFR_170-200_PKT37 | (none) | 0.253 | 0.044 | 0.045 | 0.317 |

## 5. Caveats

1. **Pre-execution affinity**: All ligand support levels are provisional (`pending_*`). Re-score after Vina execution.
2. **Single receptor state**: Only 3GT8_raw has pocket data. State robustness axis has limited discrimination.
3. **Single pocket tool**: Only fpocket. P2Rank integration will improve druggability axis.
4. **Weights are initial**: Axis weights can be tuned after expert review of first-pass results.

---

Generated by `egfr_pipeline.phase4.score_framework`