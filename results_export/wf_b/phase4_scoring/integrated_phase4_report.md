# Integrated Phase 4 Report: Perturbation Relevance Scoring

## 1. Executive Summary

- **Receptor states**: 3 (3GT8_raw, EGFR_160-185, EGFR_170-200)
- **Candidate pockets**: 103
- **Ligands evaluated**: 3 (173940, 97806, VAX-C12_0)
- **Total ranked entries**: 207

### Mechanistic Class Distribution

| Class | Count | Pockets |
|-------|-------|---------|
| allosteric_modulator_candidate | 6 | 3GT8_raw_PKT10, EGFR_170-200_PKT34 |
| interface_rim_modulator_candidate | 63 | 3GT8_raw_PKT01, 3GT8_raw_PKT02, 3GT8_raw_PKT05, 3GT8_raw_PKT07, 3GT8_raw_PKT08, 3GT8_raw_PKT11, 3GT8_raw_PKT14, 3GT8_raw_PKT20, EGFR_160-185_PKT02, EGFR_160-185_PKT07, EGFR_160-185_PKT09, EGFR_160-185_PKT12, EGFR_160-185_PKT15, EGFR_160-185_PKT16, EGFR_160-185_PKT23, EGFR_160-185_PKT34, EGFR_170-200_PKT06, EGFR_170-200_PKT10, EGFR_170-200_PKT12, EGFR_170-200_PKT17, EGFR_170-200_PKT28 |
| ligandable_but_ppi_irrelevant_candidate | 51 | 3GT8_raw_PKT04, 3GT8_raw_PKT06, 3GT8_raw_PKT09, 3GT8_raw_PKT12, 3GT8_raw_PKT13, 3GT8_raw_PKT15, 3GT8_raw_PKT17, 3GT8_raw_PKT18, 3GT8_raw_PKT19, 3GT8_raw_PKT21, EGFR_160-185_PKT06, EGFR_160-185_PKT08, EGFR_160-185_PKT10, EGFR_160-185_PKT17, EGFR_160-185_PKT18, EGFR_160-185_PKT20, EGFR_160-185_PKT21, EGFR_160-185_PKT22, EGFR_160-185_PKT24, EGFR_160-185_PKT25, EGFR_160-185_PKT26, EGFR_160-185_PKT28, EGFR_160-185_PKT31, EGFR_160-185_PKT32, EGFR_160-185_PKT33, EGFR_160-185_PKT35, EGFR_160-185_PKT36, EGFR_160-185_PKT41, EGFR_160-185_PKT42, EGFR_160-185_PKT44, EGFR_170-200_PKT02, EGFR_170-200_PKT03, EGFR_170-200_PKT07, EGFR_170-200_PKT08, EGFR_170-200_PKT09, EGFR_170-200_PKT11, EGFR_170-200_PKT14, EGFR_170-200_PKT16, EGFR_170-200_PKT18, EGFR_170-200_PKT19, EGFR_170-200_PKT20, EGFR_170-200_PKT21, EGFR_170-200_PKT22, EGFR_170-200_PKT23, EGFR_170-200_PKT24, EGFR_170-200_PKT27, EGFR_170-200_PKT29, EGFR_170-200_PKT30, EGFR_170-200_PKT32, EGFR_170-200_PKT35, EGFR_170-200_PKT37 |
| uncertain_mechanism_candidate | 87 | 3GT8_raw_PKT03, 3GT8_raw_PKT16, 3GT8_raw_PKT22, EGFR_160-185_PKT01, EGFR_160-185_PKT03, EGFR_160-185_PKT04, EGFR_160-185_PKT05, EGFR_160-185_PKT11, EGFR_160-185_PKT13, EGFR_160-185_PKT14, EGFR_160-185_PKT19, EGFR_160-185_PKT27, EGFR_160-185_PKT29, EGFR_160-185_PKT30, EGFR_160-185_PKT37, EGFR_160-185_PKT38, EGFR_160-185_PKT39, EGFR_160-185_PKT40, EGFR_160-185_PKT43, EGFR_170-200_PKT01, EGFR_170-200_PKT04, EGFR_170-200_PKT05, EGFR_170-200_PKT13, EGFR_170-200_PKT15, EGFR_170-200_PKT25, EGFR_170-200_PKT26, EGFR_170-200_PKT31, EGFR_170-200_PKT33, EGFR_170-200_PKT36 |

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

### 4.2. Interface Rim Modulators

| Rank | Pocket | Ligand | Score | A1 | A2 | A3 | A4 | Confidence | Overlap |
|------|--------|--------|-------|----|----|----|----|------------|---------|
| 1 | 3GT8_raw_PKT07 | 173940 | 0.5412 | 0.280 | 0.683 | 0.630 | 0.650 | low | 1 |
| 2 | 3GT8_raw_PKT07 | 97806 | 0.5412 | 0.280 | 0.683 | 0.630 | 0.650 | low | 1 |
| 3 | 3GT8_raw_PKT07 | VAX-C12_0 | 0.5412 | 0.280 | 0.683 | 0.630 | 0.650 | low | 1 |
| 7 | EGFR_160-185_PKT02 | 173940 | 0.4327 | 0.310 | 0.195 | 0.644 | 0.650 | medium | 2 |
| 8 | EGFR_160-185_PKT02 | 97806 | 0.4327 | 0.310 | 0.195 | 0.644 | 0.650 | medium | 2 |
| 9 | EGFR_160-185_PKT02 | VAX-C12_0 | 0.4327 | 0.310 | 0.195 | 0.644 | 0.650 | medium | 2 |
| 13 | EGFR_170-200_PKT17 | 173940 | 0.43 | 0.339 | 0.133 | 0.659 | 0.650 | medium | 3 |
| 14 | EGFR_170-200_PKT17 | VAX-C12_0 | 0.43 | 0.339 | 0.133 | 0.659 | 0.650 | medium | 3 |
| 15 | EGFR_170-200_PKT06 | 173940 | 0.4221 | 0.310 | 0.153 | 0.644 | 0.650 | medium | 2 |
| 16 | EGFR_170-200_PKT06 | 97806 | 0.4221 | 0.310 | 0.153 | 0.644 | 0.650 | medium | 2 |
| 17 | EGFR_170-200_PKT06 | VAX-C12_0 | 0.4221 | 0.310 | 0.153 | 0.644 | 0.650 | medium | 2 |
| 18 | EGFR_160-185_PKT16 | 173940 | 0.4191 | 0.310 | 0.141 | 0.644 | 0.650 | medium | 2 |
| 19 | EGFR_160-185_PKT16 | 97806 | 0.4191 | 0.310 | 0.141 | 0.644 | 0.650 | medium | 2 |
| 20 | EGFR_160-185_PKT16 | VAX-C12_0 | 0.4191 | 0.310 | 0.141 | 0.644 | 0.650 | medium | 2 |
| 21 | 3GT8_raw_PKT05 | 173940 | 0.4163 | 0.306 | 0.135 | 0.644 | 0.650 | medium | 2 |
| 22 | 3GT8_raw_PKT05 | 97806 | 0.4163 | 0.306 | 0.135 | 0.644 | 0.650 | medium | 2 |
| 23 | 3GT8_raw_PKT05 | VAX-C12_0 | 0.4163 | 0.306 | 0.135 | 0.644 | 0.650 | medium | 2 |
| 24 | 3GT8_raw_PKT11 | 173940 | 0.4123 | 0.306 | 0.118 | 0.644 | 0.650 | medium | 2 |
| 25 | 3GT8_raw_PKT11 | 97806 | 0.4123 | 0.306 | 0.118 | 0.644 | 0.650 | medium | 2 |
| 26 | 3GT8_raw_PKT11 | VAX-C12_0 | 0.4123 | 0.306 | 0.118 | 0.644 | 0.650 | medium | 2 |
| 27 | 3GT8_raw_PKT01 | 173940 | 0.4095 | 0.280 | 0.157 | 0.630 | 0.650 | low | 1 |
| 28 | 3GT8_raw_PKT01 | 97806 | 0.4095 | 0.280 | 0.157 | 0.630 | 0.650 | low | 1 |
| 29 | 3GT8_raw_PKT01 | VAX-C12_0 | 0.4095 | 0.280 | 0.157 | 0.630 | 0.650 | low | 1 |
| 30 | EGFR_160-185_PKT07 | 97806 | 0.4091 | 0.281 | 0.153 | 0.630 | 0.650 | low | 1 |
| 31 | EGFR_160-185_PKT07 | VAX-C12_0 | 0.4091 | 0.281 | 0.153 | 0.630 | 0.650 | low | 1 |
| 32 | 3GT8_raw_PKT02 | 97806 | 0.408 | 0.280 | 0.151 | 0.630 | 0.650 | low | 1 |
| 33 | 3GT8_raw_PKT02 | VAX-C12_0 | 0.408 | 0.280 | 0.151 | 0.630 | 0.650 | low | 1 |
| 34 | EGFR_160-185_PKT12 | 173940 | 0.4066 | 0.281 | 0.143 | 0.630 | 0.650 | low | 1 |
| 35 | EGFR_160-185_PKT12 | 97806 | 0.4066 | 0.281 | 0.143 | 0.630 | 0.650 | low | 1 |
| 36 | EGFR_160-185_PKT12 | VAX-C12_0 | 0.4066 | 0.281 | 0.143 | 0.630 | 0.650 | low | 1 |
| 37 | EGFR_170-200_PKT10 | 173940 | 0.4063 | 0.281 | 0.142 | 0.630 | 0.650 | low | 1 |
| 38 | EGFR_170-200_PKT10 | 97806 | 0.4063 | 0.281 | 0.142 | 0.630 | 0.650 | low | 1 |
| 39 | EGFR_170-200_PKT10 | VAX-C12_0 | 0.4063 | 0.281 | 0.142 | 0.630 | 0.650 | low | 1 |
| 40 | EGFR_160-185_PKT07 | 173940 | 0.4061 | 0.281 | 0.153 | 0.620 | 0.650 | low | 1 |
| 41 | EGFR_160-185_PKT15 | 173940 | 0.4061 | 0.281 | 0.141 | 0.630 | 0.650 | low | 1 |
| 42 | EGFR_160-185_PKT15 | 97806 | 0.4061 | 0.281 | 0.141 | 0.630 | 0.650 | low | 1 |
| 43 | EGFR_160-185_PKT15 | VAX-C12_0 | 0.4061 | 0.281 | 0.141 | 0.630 | 0.650 | low | 1 |
| 44 | EGFR_170-200_PKT17 | 97806 | 0.406 | 0.339 | 0.133 | 0.579 | 0.650 | medium | 3 |
| 45 | EGFR_170-200_PKT12 | 97806 | 0.4056 | 0.281 | 0.139 | 0.630 | 0.650 | low | 1 |
| 46 | 3GT8_raw_PKT02 | 173940 | 0.405 | 0.280 | 0.151 | 0.620 | 0.650 | low | 1 |
| 47 | EGFR_160-185_PKT23 | 173940 | 0.403 | 0.281 | 0.129 | 0.630 | 0.650 | low | 1 |
| 48 | EGFR_160-185_PKT23 | 97806 | 0.403 | 0.281 | 0.129 | 0.630 | 0.650 | low | 1 |
| 49 | EGFR_160-185_PKT23 | VAX-C12_0 | 0.403 | 0.281 | 0.129 | 0.630 | 0.650 | low | 1 |
| 50 | EGFR_170-200_PKT12 | 173940 | 0.4026 | 0.281 | 0.139 | 0.620 | 0.650 | low | 1 |
| 51 | EGFR_170-200_PKT12 | VAX-C12_0 | 0.4026 | 0.281 | 0.139 | 0.620 | 0.650 | low | 1 |
| 52 | 3GT8_raw_PKT08 | 173940 | 0.4009 | 0.306 | 0.123 | 0.644 | 0.567 | medium | 2 |
| 53 | 3GT8_raw_PKT14 | 173940 | 0.3991 | 0.280 | 0.115 | 0.630 | 0.650 | low | 1 |
| 54 | 3GT8_raw_PKT14 | 97806 | 0.3991 | 0.280 | 0.115 | 0.630 | 0.650 | low | 1 |
| 55 | 3GT8_raw_PKT14 | VAX-C12_0 | 0.3991 | 0.280 | 0.115 | 0.630 | 0.650 | low | 1 |
| 56 | EGFR_160-185_PKT09 | 173940 | 0.3962 | 0.281 | 0.151 | 0.630 | 0.567 | low | 1 |
| 57 | EGFR_160-185_PKT09 | VAX-C12_0 | 0.3962 | 0.281 | 0.151 | 0.630 | 0.567 | low | 1 |
| 58 | 3GT8_raw_PKT20 | 173940 | 0.3893 | 0.306 | 0.097 | 0.644 | 0.533 | medium_provisional | 2 |
| 59 | 3GT8_raw_PKT20 | VAX-C12_0 | 0.3893 | 0.306 | 0.097 | 0.644 | 0.533 | medium_provisional | 2 |
| 60 | EGFR_170-200_PKT28 | 173940 | 0.3885 | 0.310 | 0.115 | 0.564 | 0.650 | medium | 2 |
| 61 | EGFR_170-200_PKT28 | VAX-C12_0 | 0.3885 | 0.310 | 0.115 | 0.564 | 0.650 | medium | 2 |
| 62 | EGFR_160-185_PKT34 | 173940 | 0.3848 | 0.281 | 0.106 | 0.630 | 0.567 | low | 1 |
| 63 | EGFR_160-185_PKT34 | VAX-C12_0 | 0.3848 | 0.281 | 0.106 | 0.630 | 0.567 | low | 1 |
| 64 | 3GT8_raw_PKT08 | 97806 | 0.3769 | 0.306 | 0.123 | 0.564 | 0.567 | medium | 2 |
| 65 | 3GT8_raw_PKT08 | VAX-C12_0 | 0.3769 | 0.306 | 0.123 | 0.564 | 0.567 | medium | 2 |
| 66 | EGFR_160-185_PKT09 | 97806 | 0.3722 | 0.281 | 0.151 | 0.550 | 0.567 | low | 1 |
| 67 | EGFR_170-200_PKT28 | 97806 | 0.3705 | 0.310 | 0.115 | 0.504 | 0.650 | medium | 2 |
| 68 | 3GT8_raw_PKT20 | 97806 | 0.3653 | 0.306 | 0.097 | 0.564 | 0.533 | medium_provisional | 2 |
| 69 | EGFR_160-185_PKT34 | 97806 | 0.3608 | 0.281 | 0.106 | 0.550 | 0.567 | low | 1 |

### 4.3. Allosteric Modulators

| Rank | Pocket | Ligand | Score | A1 | A2 | A3 | A4 | Confidence | Overlap |
|------|--------|--------|-------|----|----|----|----|------------|---------|
| 4 | EGFR_170-200_PKT34 | 173940 | 0.4921 | 0.253 | 0.645 | 0.525 | 0.650 | medium | 0 |
| 5 | EGFR_170-200_PKT34 | 97806 | 0.4921 | 0.253 | 0.645 | 0.525 | 0.650 | medium | 0 |
| 6 | EGFR_170-200_PKT34 | VAX-C12_0 | 0.4921 | 0.253 | 0.645 | 0.525 | 0.650 | medium | 0 |
| 10 | 3GT8_raw_PKT10 | 173940 | 0.4306 | 0.253 | 0.399 | 0.525 | 0.650 | low | 0 |
| 11 | 3GT8_raw_PKT10 | 97806 | 0.4306 | 0.253 | 0.399 | 0.525 | 0.650 | low | 0 |
| 12 | 3GT8_raw_PKT10 | VAX-C12_0 | 0.4306 | 0.253 | 0.399 | 0.525 | 0.650 | low | 0 |

### 4.4. Ligandable but PPI-Irrelevant

| Rank | Pocket | Ligand | Score | A1 | A2 | A3 | A4 | Confidence | Overlap |
|------|--------|--------|-------|----|----|----|----|------------|---------|
| 157 | EGFR_170-200_PKT03 | (none) | 0.1902 | 0.253 | 0.163 | 0.045 | 0.400 | high | 0 |
| 158 | EGFR_160-185_PKT06 | (none) | 0.1878 | 0.253 | 0.153 | 0.045 | 0.400 | high | 0 |
| 159 | EGFR_170-200_PKT07 | (none) | 0.1868 | 0.253 | 0.150 | 0.045 | 0.400 | high | 0 |
| 160 | EGFR_160-185_PKT17 | (none) | 0.1836 | 0.253 | 0.137 | 0.045 | 0.400 | high | 0 |
| 161 | EGFR_170-200_PKT18 | (none) | 0.1818 | 0.253 | 0.130 | 0.045 | 0.400 | high | 0 |
| 162 | 3GT8_raw_PKT13 | (none) | 0.1786 | 0.253 | 0.117 | 0.045 | 0.400 | high | 0 |
| 163 | EGFR_170-200_PKT27 | (none) | 0.1782 | 0.253 | 0.115 | 0.045 | 0.400 | high | 0 |
| 164 | EGFR_170-200_PKT29 | (none) | 0.1773 | 0.253 | 0.112 | 0.045 | 0.400 | high | 0 |
| 165 | 3GT8_raw_PKT18 | (none) | 0.176 | 0.253 | 0.106 | 0.045 | 0.400 | high | 0 |
| 166 | EGFR_160-185_PKT08 | (none) | 0.1749 | 0.253 | 0.152 | 0.045 | 0.317 | high | 0 |
| 167 | EGFR_170-200_PKT02 | (none) | 0.1749 | 0.253 | 0.172 | 0.045 | 0.283 | high_provisional | 0 |
| 168 | EGFR_160-185_PKT10 | (none) | 0.1746 | 0.253 | 0.151 | 0.045 | 0.317 | high | 0 |
| 169 | 3GT8_raw_PKT04 | (none) | 0.1735 | 0.253 | 0.146 | 0.045 | 0.317 | high | 0 |
| 170 | EGFR_170-200_PKT08 | (none) | 0.1733 | 0.253 | 0.146 | 0.045 | 0.317 | high | 0 |
| 171 | EGFR_170-200_PKT09 | (none) | 0.173 | 0.253 | 0.144 | 0.045 | 0.317 | high | 0 |
| 172 | EGFR_170-200_PKT11 | (none) | 0.1718 | 0.253 | 0.140 | 0.045 | 0.317 | high | 0 |
| 173 | EGFR_170-200_PKT14 | (none) | 0.1711 | 0.253 | 0.137 | 0.045 | 0.317 | high | 0 |
| 174 | EGFR_160-185_PKT18 | (none) | 0.1709 | 0.253 | 0.136 | 0.045 | 0.317 | high | 0 |
| 175 | EGFR_160-185_PKT20 | (none) | 0.1703 | 0.253 | 0.133 | 0.045 | 0.317 | high | 0 |
| 176 | EGFR_160-185_PKT21 | (none) | 0.1702 | 0.253 | 0.133 | 0.045 | 0.317 | high | 0 |
| 177 | EGFR_170-200_PKT16 | (none) | 0.1702 | 0.253 | 0.133 | 0.045 | 0.317 | high | 0 |
| 178 | EGFR_160-185_PKT22 | (none) | 0.1697 | 0.253 | 0.131 | 0.045 | 0.317 | high | 0 |
| 179 | EGFR_160-185_PKT24 | (none) | 0.169 | 0.253 | 0.128 | 0.045 | 0.317 | high | 0 |
| 180 | EGFR_160-185_PKT25 | (none) | 0.1689 | 0.253 | 0.128 | 0.045 | 0.317 | high | 0 |
| 181 | EGFR_170-200_PKT19 | (none) | 0.1685 | 0.253 | 0.126 | 0.045 | 0.317 | high | 0 |
| 182 | EGFR_160-185_PKT26 | (none) | 0.1676 | 0.253 | 0.123 | 0.045 | 0.317 | high | 0 |
| 183 | EGFR_170-200_PKT20 | (none) | 0.1673 | 0.253 | 0.122 | 0.045 | 0.317 | high | 0 |
| 184 | EGFR_170-200_PKT21 | (none) | 0.1672 | 0.253 | 0.121 | 0.045 | 0.317 | high | 0 |
| 185 | 3GT8_raw_PKT09 | (none) | 0.1671 | 0.253 | 0.121 | 0.045 | 0.317 | high | 0 |
| 186 | EGFR_170-200_PKT22 | (none) | 0.1671 | 0.253 | 0.121 | 0.045 | 0.317 | high | 0 |
| 187 | EGFR_170-200_PKT23 | (none) | 0.1671 | 0.253 | 0.121 | 0.045 | 0.317 | high | 0 |
| 188 | EGFR_170-200_PKT24 | (none) | 0.167 | 0.253 | 0.120 | 0.045 | 0.317 | high | 0 |
| 189 | EGFR_160-185_PKT28 | (none) | 0.1667 | 0.253 | 0.119 | 0.045 | 0.317 | high | 0 |
| 190 | 3GT8_raw_PKT12 | (none) | 0.1662 | 0.253 | 0.117 | 0.045 | 0.317 | high | 0 |
| 191 | EGFR_160-185_PKT31 | (none) | 0.1652 | 0.253 | 0.113 | 0.045 | 0.317 | high | 0 |
| 192 | EGFR_160-185_PKT32 | (none) | 0.1652 | 0.253 | 0.113 | 0.045 | 0.317 | high | 0 |
| 193 | 3GT8_raw_PKT06 | (none) | 0.165 | 0.253 | 0.133 | 0.045 | 0.283 | high_provisional | 0 |
| 194 | EGFR_160-185_PKT33 | (none) | 0.1649 | 0.253 | 0.112 | 0.045 | 0.317 | high | 0 |
| 195 | 3GT8_raw_PKT17 | (none) | 0.164 | 0.253 | 0.108 | 0.045 | 0.317 | high | 0 |
| 196 | EGFR_160-185_PKT35 | (none) | 0.1634 | 0.253 | 0.106 | 0.045 | 0.317 | high | 0 |
| 197 | EGFR_160-185_PKT36 | (none) | 0.1634 | 0.253 | 0.106 | 0.045 | 0.317 | high | 0 |
| 198 | EGFR_170-200_PKT30 | (none) | 0.1625 | 0.253 | 0.102 | 0.045 | 0.317 | high | 0 |
| 199 | 3GT8_raw_PKT19 | (none) | 0.1623 | 0.253 | 0.102 | 0.045 | 0.317 | high | 0 |
| 200 | 3GT8_raw_PKT15 | (none) | 0.1607 | 0.253 | 0.115 | 0.045 | 0.283 | high_provisional | 0 |
| 201 | 3GT8_raw_PKT21 | (none) | 0.1601 | 0.253 | 0.093 | 0.045 | 0.317 | high | 0 |
| 202 | EGFR_170-200_PKT32 | (none) | 0.1601 | 0.253 | 0.093 | 0.045 | 0.317 | high | 0 |
| 203 | EGFR_160-185_PKT42 | (none) | 0.159 | 0.253 | 0.088 | 0.045 | 0.317 | high | 0 |
| 204 | EGFR_160-185_PKT44 | (none) | 0.1581 | 0.253 | 0.085 | 0.045 | 0.317 | high | 0 |
| 205 | EGFR_160-185_PKT41 | (none) | 0.1551 | 0.253 | 0.093 | 0.045 | 0.283 | high_provisional | 0 |
| 206 | EGFR_170-200_PKT35 | (none) | 0.151 | 0.253 | 0.056 | 0.045 | 0.317 | high | 0 |
| 207 | EGFR_170-200_PKT37 | (none) | 0.1478 | 0.253 | 0.044 | 0.045 | 0.317 | high | 0 |

### 4.5. Uncertain Mechanism

| Rank | Pocket | Ligand | Score | A1 | A2 | A3 | A4 | Confidence | Overlap |
|------|--------|--------|-------|----|----|----|----|------------|---------|
| 70 | 3GT8_raw_PKT03 | 173940 | 0.2668 | 0.253 | 0.149 | 0.475 | 0.650 | low | 0 |
| 71 | 3GT8_raw_PKT03 | 97806 | 0.2668 | 0.253 | 0.149 | 0.475 | 0.650 | low | 0 |
| 72 | 3GT8_raw_PKT03 | VAX-C12_0 | 0.2668 | 0.253 | 0.149 | 0.475 | 0.650 | low | 0 |
| 73 | 3GT8_raw_PKT16 | 173940 | 0.2668 | 0.253 | 0.114 | 0.475 | 0.650 | low | 0 |
| 74 | 3GT8_raw_PKT16 | 97806 | 0.2668 | 0.253 | 0.114 | 0.475 | 0.650 | low | 0 |
| 75 | 3GT8_raw_PKT16 | VAX-C12_0 | 0.2668 | 0.253 | 0.114 | 0.475 | 0.650 | low | 0 |
| 76 | EGFR_160-185_PKT01 | 173940 | 0.2668 | 0.253 | 0.204 | 0.475 | 0.650 | low | 0 |
| 77 | EGFR_160-185_PKT01 | 97806 | 0.2668 | 0.253 | 0.204 | 0.475 | 0.650 | low | 0 |
| 78 | EGFR_160-185_PKT01 | VAX-C12_0 | 0.2668 | 0.253 | 0.204 | 0.475 | 0.650 | low | 0 |
| 79 | EGFR_160-185_PKT03 | 173940 | 0.2668 | 0.253 | 0.180 | 0.475 | 0.650 | low | 0 |
| 80 | EGFR_160-185_PKT03 | 97806 | 0.2668 | 0.253 | 0.180 | 0.475 | 0.650 | low | 0 |
| 81 | EGFR_160-185_PKT03 | VAX-C12_0 | 0.2668 | 0.253 | 0.180 | 0.475 | 0.650 | low | 0 |
| 82 | EGFR_160-185_PKT04 | 173940 | 0.2668 | 0.253 | 0.167 | 0.475 | 0.650 | low | 0 |
| 83 | EGFR_160-185_PKT04 | 97806 | 0.2668 | 0.253 | 0.167 | 0.475 | 0.650 | low | 0 |
| 84 | EGFR_160-185_PKT04 | VAX-C12_0 | 0.2668 | 0.253 | 0.167 | 0.475 | 0.650 | low | 0 |
| 85 | EGFR_160-185_PKT05 | 173940 | 0.2668 | 0.253 | 0.161 | 0.475 | 0.650 | low | 0 |
| 86 | EGFR_160-185_PKT05 | 97806 | 0.2668 | 0.253 | 0.161 | 0.475 | 0.650 | low | 0 |
| 87 | EGFR_160-185_PKT05 | VAX-C12_0 | 0.2668 | 0.253 | 0.161 | 0.475 | 0.650 | low | 0 |
| 88 | EGFR_160-185_PKT11 | 173940 | 0.2668 | 0.253 | 0.144 | 0.395 | 0.650 | low | 0 |
| 89 | EGFR_160-185_PKT11 | 97806 | 0.2668 | 0.253 | 0.144 | 0.395 | 0.650 | low | 0 |
| 90 | EGFR_160-185_PKT11 | VAX-C12_0 | 0.2668 | 0.253 | 0.144 | 0.395 | 0.650 | low | 0 |
| 91 | EGFR_160-185_PKT13 | 173940 | 0.2668 | 0.253 | 0.142 | 0.475 | 0.650 | low | 0 |
| 92 | EGFR_160-185_PKT13 | 97806 | 0.2668 | 0.253 | 0.142 | 0.385 | 0.650 | low | 0 |
| 93 | EGFR_160-185_PKT13 | VAX-C12_0 | 0.2668 | 0.253 | 0.142 | 0.475 | 0.650 | low | 0 |
| 94 | EGFR_160-185_PKT19 | 173940 | 0.2668 | 0.253 | 0.135 | 0.475 | 0.650 | low | 0 |
| 95 | EGFR_160-185_PKT19 | 97806 | 0.2668 | 0.253 | 0.135 | 0.475 | 0.650 | low | 0 |
| 96 | EGFR_160-185_PKT19 | VAX-C12_0 | 0.2668 | 0.253 | 0.135 | 0.475 | 0.650 | low | 0 |
| 97 | EGFR_160-185_PKT30 | 173940 | 0.2668 | 0.253 | 0.114 | 0.475 | 0.650 | low | 0 |
| 98 | EGFR_160-185_PKT30 | 97806 | 0.2668 | 0.253 | 0.114 | 0.475 | 0.650 | low | 0 |
| 99 | EGFR_160-185_PKT30 | VAX-C12_0 | 0.2668 | 0.253 | 0.114 | 0.475 | 0.650 | low | 0 |
| 100 | EGFR_160-185_PKT37 | 173940 | 0.2668 | 0.253 | 0.104 | 0.475 | 0.650 | low | 0 |
| 101 | EGFR_160-185_PKT37 | 97806 | 0.2668 | 0.253 | 0.104 | 0.475 | 0.650 | low | 0 |
| 102 | EGFR_160-185_PKT37 | VAX-C12_0 | 0.2668 | 0.253 | 0.104 | 0.475 | 0.650 | low | 0 |
| 103 | EGFR_160-185_PKT38 | 173940 | 0.2668 | 0.253 | 0.104 | 0.475 | 0.650 | low | 0 |
| 104 | EGFR_160-185_PKT38 | 97806 | 0.2668 | 0.253 | 0.104 | 0.475 | 0.650 | low | 0 |
| 105 | EGFR_160-185_PKT38 | VAX-C12_0 | 0.2668 | 0.253 | 0.104 | 0.475 | 0.650 | low | 0 |
| 106 | EGFR_170-200_PKT01 | 173940 | 0.2668 | 0.253 | 0.218 | 0.475 | 0.650 | low | 0 |
| 107 | EGFR_170-200_PKT01 | 97806 | 0.2668 | 0.253 | 0.218 | 0.475 | 0.650 | low | 0 |
| 108 | EGFR_170-200_PKT01 | VAX-C12_0 | 0.2668 | 0.253 | 0.218 | 0.475 | 0.650 | low | 0 |
| 109 | EGFR_170-200_PKT04 | 173940 | 0.2668 | 0.253 | 0.161 | 0.295 | 0.650 | low | 0 |
| 110 | EGFR_170-200_PKT04 | 97806 | 0.2668 | 0.253 | 0.161 | 0.295 | 0.650 | low | 0 |
| 111 | EGFR_170-200_PKT04 | VAX-C12_0 | 0.2668 | 0.253 | 0.161 | 0.305 | 0.650 | low | 0 |
| 112 | EGFR_170-200_PKT05 | 173940 | 0.2668 | 0.253 | 0.157 | 0.475 | 0.650 | low | 0 |
| 113 | EGFR_170-200_PKT05 | 97806 | 0.2668 | 0.253 | 0.157 | 0.475 | 0.650 | low | 0 |
| 114 | EGFR_170-200_PKT05 | VAX-C12_0 | 0.2668 | 0.253 | 0.157 | 0.475 | 0.650 | low | 0 |
| 115 | EGFR_170-200_PKT15 | 173940 | 0.2668 | 0.253 | 0.134 | 0.475 | 0.650 | low | 0 |
| 116 | EGFR_170-200_PKT15 | 97806 | 0.2668 | 0.253 | 0.134 | 0.475 | 0.650 | low | 0 |
| 117 | EGFR_170-200_PKT15 | VAX-C12_0 | 0.2668 | 0.253 | 0.134 | 0.475 | 0.650 | low | 0 |
| 118 | EGFR_170-200_PKT25 | 173940 | 0.2668 | 0.253 | 0.116 | 0.475 | 0.650 | low | 0 |
| 119 | EGFR_170-200_PKT25 | 97806 | 0.2668 | 0.253 | 0.116 | 0.475 | 0.650 | low | 0 |
| 120 | EGFR_170-200_PKT25 | VAX-C12_0 | 0.2668 | 0.253 | 0.116 | 0.475 | 0.650 | low | 0 |
| 121 | EGFR_170-200_PKT31 | 173940 | 0.2668 | 0.253 | 0.100 | 0.475 | 0.650 | low | 0 |
| 122 | EGFR_170-200_PKT31 | 97806 | 0.2668 | 0.253 | 0.100 | 0.475 | 0.650 | low | 0 |
| 123 | EGFR_170-200_PKT31 | VAX-C12_0 | 0.2668 | 0.253 | 0.100 | 0.475 | 0.650 | low | 0 |
| 124 | EGFR_170-200_PKT36 | 173940 | 0.2668 | 0.253 | 0.044 | 0.475 | 0.650 | low | 0 |
| 125 | EGFR_170-200_PKT36 | 97806 | 0.2668 | 0.253 | 0.044 | 0.475 | 0.650 | low | 0 |
| 126 | EGFR_170-200_PKT36 | VAX-C12_0 | 0.2668 | 0.253 | 0.044 | 0.475 | 0.650 | low | 0 |
| 127 | 3GT8_raw_PKT22 | 173940 | 0.2475 | 0.253 | 0.024 | 0.475 | 0.567 | low | 0 |
| 128 | 3GT8_raw_PKT22 | 97806 | 0.2475 | 0.253 | 0.024 | 0.475 | 0.567 | low | 0 |
| 129 | 3GT8_raw_PKT22 | VAX-C12_0 | 0.2475 | 0.253 | 0.024 | 0.475 | 0.567 | low | 0 |
| 130 | EGFR_160-185_PKT14 | 173940 | 0.2475 | 0.253 | 0.141 | 0.475 | 0.567 | low | 0 |
| 131 | EGFR_160-185_PKT14 | 97806 | 0.2475 | 0.253 | 0.141 | 0.475 | 0.567 | low | 0 |
| 132 | EGFR_160-185_PKT14 | VAX-C12_0 | 0.2475 | 0.253 | 0.141 | 0.475 | 0.567 | low | 0 |
| 133 | EGFR_160-185_PKT27 | 173940 | 0.2475 | 0.253 | 0.121 | 0.475 | 0.567 | low | 0 |
| 134 | EGFR_160-185_PKT27 | 97806 | 0.2475 | 0.253 | 0.121 | 0.475 | 0.567 | low | 0 |
| 135 | EGFR_160-185_PKT27 | VAX-C12_0 | 0.2475 | 0.253 | 0.121 | 0.475 | 0.567 | low | 0 |
| 136 | EGFR_160-185_PKT29 | 173940 | 0.2475 | 0.253 | 0.117 | 0.475 | 0.567 | low | 0 |
| 137 | EGFR_160-185_PKT29 | 97806 | 0.2475 | 0.253 | 0.117 | 0.475 | 0.567 | low | 0 |
| 138 | EGFR_160-185_PKT29 | VAX-C12_0 | 0.2475 | 0.253 | 0.117 | 0.475 | 0.567 | low | 0 |
| 139 | EGFR_160-185_PKT39 | 173940 | 0.2475 | 0.253 | 0.104 | 0.475 | 0.567 | low | 0 |
| 140 | EGFR_160-185_PKT39 | 97806 | 0.2475 | 0.253 | 0.104 | 0.475 | 0.567 | low | 0 |
| 141 | EGFR_160-185_PKT39 | VAX-C12_0 | 0.2475 | 0.253 | 0.104 | 0.475 | 0.567 | low | 0 |
| 142 | EGFR_160-185_PKT40 | 173940 | 0.2475 | 0.253 | 0.094 | 0.475 | 0.567 | low | 0 |
| 143 | EGFR_160-185_PKT40 | 97806 | 0.2475 | 0.253 | 0.094 | 0.385 | 0.567 | low | 0 |
| 144 | EGFR_160-185_PKT40 | VAX-C12_0 | 0.2475 | 0.253 | 0.094 | 0.475 | 0.567 | low | 0 |
| 145 | EGFR_160-185_PKT43 | 173940 | 0.2475 | 0.253 | 0.088 | 0.475 | 0.567 | low | 0 |
| 146 | EGFR_160-185_PKT43 | 97806 | 0.2475 | 0.253 | 0.088 | 0.475 | 0.567 | low | 0 |
| 147 | EGFR_160-185_PKT43 | VAX-C12_0 | 0.2475 | 0.253 | 0.088 | 0.475 | 0.567 | low | 0 |
| 148 | EGFR_170-200_PKT13 | 173940 | 0.2475 | 0.253 | 0.139 | 0.475 | 0.567 | low | 0 |
| 149 | EGFR_170-200_PKT13 | 97806 | 0.2475 | 0.253 | 0.139 | 0.475 | 0.567 | low | 0 |
| 150 | EGFR_170-200_PKT13 | VAX-C12_0 | 0.2475 | 0.253 | 0.139 | 0.475 | 0.567 | low | 0 |
| 151 | EGFR_170-200_PKT26 | 173940 | 0.2475 | 0.253 | 0.116 | 0.395 | 0.567 | low | 0 |
| 152 | EGFR_170-200_PKT26 | 97806 | 0.2475 | 0.253 | 0.116 | 0.395 | 0.567 | low | 0 |
| 153 | EGFR_170-200_PKT26 | VAX-C12_0 | 0.2475 | 0.253 | 0.116 | 0.475 | 0.567 | low | 0 |
| 154 | EGFR_170-200_PKT33 | 173940 | 0.2475 | 0.253 | 0.090 | 0.475 | 0.567 | low | 0 |
| 155 | EGFR_170-200_PKT33 | 97806 | 0.2475 | 0.253 | 0.090 | 0.475 | 0.567 | low | 0 |
| 156 | EGFR_170-200_PKT33 | VAX-C12_0 | 0.2475 | 0.253 | 0.090 | 0.475 | 0.567 | low | 0 |

## 5. State Robustness Assessment

**3GT8_raw_PKT07**
- State class: state_robust_pocket
- Interpretation: unknown
- Accessibility: always_accessible
- States matched: 3

**EGFR_170-200_PKT34**
- State class: state_shifted_pocket
- Interpretation: unknown
- Accessibility: always_accessible
- States matched: 3

**EGFR_160-185_PKT02**
- State class: state_shifted_pocket
- Interpretation: unknown
- Accessibility: always_accessible
- States matched: 3

**3GT8_raw_PKT10**
- State class: state_robust_pocket
- Interpretation: unknown
- Accessibility: always_accessible
- States matched: 3

**EGFR_170-200_PKT17**
- State class: state_shifted_pocket
- Interpretation: unknown
- Accessibility: always_accessible
- States matched: 3

**EGFR_170-200_PKT06**
- State class: state_shifted_pocket
- Interpretation: unknown
- Accessibility: always_accessible
- States matched: 3

**EGFR_160-185_PKT16**
- State class: state_shifted_pocket
- Interpretation: unknown
- Accessibility: always_accessible
- States matched: 3

**3GT8_raw_PKT05**
- State class: state_robust_pocket
- Interpretation: unknown
- Accessibility: always_accessible
- States matched: 3

**3GT8_raw_PKT11**
- State class: state_robust_pocket
- Interpretation: unknown
- Accessibility: always_accessible
- States matched: 3

**3GT8_raw_PKT01**
- State class: state_robust_pocket
- Interpretation: unknown
- Accessibility: always_accessible
- States matched: 3

**EGFR_160-185_PKT07**
- State class: state_shifted_pocket
- Interpretation: unknown
- Accessibility: always_accessible
- States matched: 3

**3GT8_raw_PKT02**
- State class: state_robust_pocket
- Interpretation: unknown
- Accessibility: always_accessible
- States matched: 3

**EGFR_160-185_PKT12**
- State class: state_shifted_pocket
- Interpretation: unknown
- Accessibility: always_accessible
- States matched: 3

**EGFR_170-200_PKT10**
- State class: state_shifted_pocket
- Interpretation: unknown
- Accessibility: always_accessible
- States matched: 3

**EGFR_160-185_PKT15**
- State class: state_shifted_pocket
- Interpretation: unknown
- Accessibility: always_accessible
- States matched: 3

**EGFR_170-200_PKT12**
- State class: state_shifted_pocket
- Interpretation: unknown
- Accessibility: always_accessible
- States matched: 3

**EGFR_160-185_PKT23**
- State class: state_robust_pocket
- Interpretation: unknown
- Accessibility: always_accessible
- States matched: 3

**3GT8_raw_PKT08**
- State class: uncertain_alignment
- Interpretation: unknown
- Accessibility: mostly_accessible
- States matched: 2

**3GT8_raw_PKT14**
- State class: state_robust_pocket
- Interpretation: unknown
- Accessibility: always_accessible
- States matched: 3

**EGFR_160-185_PKT09**
- State class: uncertain_alignment
- Interpretation: unknown
- Accessibility: mostly_accessible
- States matched: 2

**3GT8_raw_PKT20**
- State class: state_specific_pocket
- Interpretation: state_dependent
- Accessibility: single_state_only
- States matched: 1
- Caveat: Only 1 receptor state has pocket data. Cannot determine if state-specific or under-sampled. NOTE: Single-state data. Cross-state persistence unknown. Re-evaluate when additional receptor states have pocket data.

**EGFR_170-200_PKT28**
- State class: state_shifted_pocket
- Interpretation: unknown
- Accessibility: always_accessible
- States matched: 3

**EGFR_160-185_PKT34**
- State class: uncertain_alignment
- Interpretation: unknown
- Accessibility: mostly_accessible
- States matched: 2

**3GT8_raw_PKT03**
- State class: state_robust_pocket
- Interpretation: unknown
- Accessibility: always_accessible
- States matched: 3

**3GT8_raw_PKT16**
- State class: state_robust_pocket
- Interpretation: unknown
- Accessibility: always_accessible
- States matched: 3

**EGFR_160-185_PKT01**
- State class: state_shifted_pocket
- Interpretation: unknown
- Accessibility: always_accessible
- States matched: 3

**EGFR_160-185_PKT03**
- State class: state_shifted_pocket
- Interpretation: unknown
- Accessibility: always_accessible
- States matched: 3

**EGFR_160-185_PKT04**
- State class: state_shifted_pocket
- Interpretation: unknown
- Accessibility: always_accessible
- States matched: 3

**EGFR_160-185_PKT05**
- State class: state_shifted_pocket
- Interpretation: unknown
- Accessibility: always_accessible
- States matched: 3

**EGFR_160-185_PKT11**
- State class: state_shifted_pocket
- Interpretation: unknown
- Accessibility: always_accessible
- States matched: 3

**EGFR_160-185_PKT13**
- State class: state_shifted_pocket
- Interpretation: unknown
- Accessibility: always_accessible
- States matched: 3

**EGFR_160-185_PKT19**
- State class: state_shifted_pocket
- Interpretation: unknown
- Accessibility: always_accessible
- States matched: 3

**EGFR_160-185_PKT30**
- State class: state_shifted_pocket
- Interpretation: unknown
- Accessibility: always_accessible
- States matched: 3

**EGFR_160-185_PKT37**
- State class: state_shifted_pocket
- Interpretation: unknown
- Accessibility: always_accessible
- States matched: 3

**EGFR_160-185_PKT38**
- State class: state_shifted_pocket
- Interpretation: unknown
- Accessibility: always_accessible
- States matched: 3

**EGFR_170-200_PKT01**
- State class: state_shifted_pocket
- Interpretation: unknown
- Accessibility: always_accessible
- States matched: 3

**EGFR_170-200_PKT04**
- State class: state_shifted_pocket
- Interpretation: unknown
- Accessibility: always_accessible
- States matched: 3

**EGFR_170-200_PKT05**
- State class: state_shifted_pocket
- Interpretation: unknown
- Accessibility: always_accessible
- States matched: 3

**EGFR_170-200_PKT15**
- State class: state_shifted_pocket
- Interpretation: unknown
- Accessibility: always_accessible
- States matched: 3

**EGFR_170-200_PKT25**
- State class: state_shifted_pocket
- Interpretation: unknown
- Accessibility: always_accessible
- States matched: 3

**EGFR_170-200_PKT31**
- State class: state_shifted_pocket
- Interpretation: unknown
- Accessibility: always_accessible
- States matched: 3

**EGFR_170-200_PKT36**
- State class: state_shifted_pocket
- Interpretation: unknown
- Accessibility: always_accessible
- States matched: 3

**3GT8_raw_PKT22**
- State class: uncertain_alignment
- Interpretation: unknown
- Accessibility: mostly_accessible
- States matched: 2

**EGFR_160-185_PKT14**
- State class: state_shifted_pocket
- Interpretation: unknown
- Accessibility: mostly_accessible
- States matched: 2

**EGFR_160-185_PKT27**
- State class: state_shifted_pocket
- Interpretation: unknown
- Accessibility: mostly_accessible
- States matched: 2

**EGFR_160-185_PKT29**
- State class: state_shifted_pocket
- Interpretation: unknown
- Accessibility: mostly_accessible
- States matched: 2

**EGFR_160-185_PKT39**
- State class: state_shifted_pocket
- Interpretation: unknown
- Accessibility: mostly_accessible
- States matched: 2

**EGFR_160-185_PKT40**
- State class: state_shifted_pocket
- Interpretation: unknown
- Accessibility: mostly_accessible
- States matched: 2

**EGFR_160-185_PKT43**
- State class: state_shifted_pocket
- Interpretation: unknown
- Accessibility: mostly_accessible
- States matched: 2

**EGFR_170-200_PKT13**
- State class: state_shifted_pocket
- Interpretation: unknown
- Accessibility: mostly_accessible
- States matched: 2

**EGFR_170-200_PKT26**
- State class: state_shifted_pocket
- Interpretation: unknown
- Accessibility: mostly_accessible
- States matched: 2

**EGFR_170-200_PKT33**
- State class: state_shifted_pocket
- Interpretation: unknown
- Accessibility: mostly_accessible
- States matched: 2

**EGFR_170-200_PKT03**
- State class: state_shifted_pocket
- Interpretation: unknown
- Accessibility: always_accessible
- States matched: 3

**EGFR_160-185_PKT06**
- State class: state_shifted_pocket
- Interpretation: unknown
- Accessibility: always_accessible
- States matched: 3

**EGFR_170-200_PKT07**
- State class: state_shifted_pocket
- Interpretation: unknown
- Accessibility: always_accessible
- States matched: 3

**EGFR_160-185_PKT17**
- State class: state_shifted_pocket
- Interpretation: unknown
- Accessibility: always_accessible
- States matched: 3

**EGFR_170-200_PKT18**
- State class: state_shifted_pocket
- Interpretation: unknown
- Accessibility: always_accessible
- States matched: 3

**3GT8_raw_PKT13**
- State class: state_robust_pocket
- Interpretation: unknown
- Accessibility: always_accessible
- States matched: 3

**EGFR_170-200_PKT27**
- State class: state_robust_pocket
- Interpretation: unknown
- Accessibility: always_accessible
- States matched: 3

**EGFR_170-200_PKT29**
- State class: state_shifted_pocket
- Interpretation: unknown
- Accessibility: always_accessible
- States matched: 3

**3GT8_raw_PKT18**
- State class: state_robust_pocket
- Interpretation: unknown
- Accessibility: always_accessible
- States matched: 3

**EGFR_160-185_PKT08**
- State class: state_shifted_pocket
- Interpretation: unknown
- Accessibility: mostly_accessible
- States matched: 2

**EGFR_170-200_PKT02**
- State class: state_specific_pocket
- Interpretation: state_dependent
- Accessibility: single_state_only
- States matched: 1
- Caveat: Only 1 receptor state has pocket data. Cannot determine if state-specific or under-sampled. NOTE: Single-state data. Cross-state persistence unknown. Re-evaluate when additional receptor states have pocket data.

**EGFR_160-185_PKT10**
- State class: state_shifted_pocket
- Interpretation: unknown
- Accessibility: mostly_accessible
- States matched: 2

**3GT8_raw_PKT04**
- State class: uncertain_alignment
- Interpretation: unknown
- Accessibility: mostly_accessible
- States matched: 2

**EGFR_170-200_PKT08**
- State class: state_shifted_pocket
- Interpretation: unknown
- Accessibility: mostly_accessible
- States matched: 2

**EGFR_170-200_PKT09**
- State class: state_shifted_pocket
- Interpretation: unknown
- Accessibility: mostly_accessible
- States matched: 2

**EGFR_170-200_PKT11**
- State class: uncertain_alignment
- Interpretation: unknown
- Accessibility: mostly_accessible
- States matched: 2

**EGFR_170-200_PKT14**
- State class: state_shifted_pocket
- Interpretation: unknown
- Accessibility: mostly_accessible
- States matched: 2

**EGFR_160-185_PKT18**
- State class: uncertain_alignment
- Interpretation: unknown
- Accessibility: mostly_accessible
- States matched: 2

**EGFR_160-185_PKT20**
- State class: state_shifted_pocket
- Interpretation: unknown
- Accessibility: mostly_accessible
- States matched: 2

**EGFR_160-185_PKT21**
- State class: state_shifted_pocket
- Interpretation: unknown
- Accessibility: mostly_accessible
- States matched: 2

**EGFR_170-200_PKT16**
- State class: state_shifted_pocket
- Interpretation: unknown
- Accessibility: mostly_accessible
- States matched: 2

**EGFR_160-185_PKT22**
- State class: state_shifted_pocket
- Interpretation: unknown
- Accessibility: mostly_accessible
- States matched: 2

**EGFR_160-185_PKT24**
- State class: state_shifted_pocket
- Interpretation: unknown
- Accessibility: mostly_accessible
- States matched: 2

**EGFR_160-185_PKT25**
- State class: state_shifted_pocket
- Interpretation: unknown
- Accessibility: mostly_accessible
- States matched: 2

**EGFR_170-200_PKT19**
- State class: uncertain_alignment
- Interpretation: unknown
- Accessibility: mostly_accessible
- States matched: 2

**EGFR_160-185_PKT26**
- State class: state_shifted_pocket
- Interpretation: unknown
- Accessibility: mostly_accessible
- States matched: 2

**EGFR_170-200_PKT20**
- State class: state_shifted_pocket
- Interpretation: unknown
- Accessibility: mostly_accessible
- States matched: 2

**EGFR_170-200_PKT21**
- State class: state_shifted_pocket
- Interpretation: unknown
- Accessibility: mostly_accessible
- States matched: 2

**3GT8_raw_PKT09**
- State class: uncertain_alignment
- Interpretation: unknown
- Accessibility: mostly_accessible
- States matched: 2

**EGFR_170-200_PKT22**
- State class: state_shifted_pocket
- Interpretation: unknown
- Accessibility: mostly_accessible
- States matched: 2

**EGFR_170-200_PKT23**
- State class: state_shifted_pocket
- Interpretation: unknown
- Accessibility: mostly_accessible
- States matched: 2

**EGFR_170-200_PKT24**
- State class: state_shifted_pocket
- Interpretation: unknown
- Accessibility: mostly_accessible
- States matched: 2

**EGFR_160-185_PKT28**
- State class: uncertain_alignment
- Interpretation: unknown
- Accessibility: mostly_accessible
- States matched: 2

**3GT8_raw_PKT12**
- State class: uncertain_alignment
- Interpretation: unknown
- Accessibility: mostly_accessible
- States matched: 2

**EGFR_160-185_PKT31**
- State class: uncertain_alignment
- Interpretation: unknown
- Accessibility: mostly_accessible
- States matched: 2

**EGFR_160-185_PKT32**
- State class: state_shifted_pocket
- Interpretation: unknown
- Accessibility: mostly_accessible
- States matched: 2

**3GT8_raw_PKT06**
- State class: state_specific_pocket
- Interpretation: state_dependent
- Accessibility: single_state_only
- States matched: 1
- Caveat: Only 1 receptor state has pocket data. Cannot determine if state-specific or under-sampled. NOTE: Single-state data. Cross-state persistence unknown. Re-evaluate when additional receptor states have pocket data.

**EGFR_160-185_PKT33**
- State class: state_shifted_pocket
- Interpretation: unknown
- Accessibility: mostly_accessible
- States matched: 2

**3GT8_raw_PKT17**
- State class: uncertain_alignment
- Interpretation: unknown
- Accessibility: mostly_accessible
- States matched: 2

**EGFR_160-185_PKT35**
- State class: state_shifted_pocket
- Interpretation: unknown
- Accessibility: mostly_accessible
- States matched: 2

**EGFR_160-185_PKT36**
- State class: uncertain_alignment
- Interpretation: unknown
- Accessibility: mostly_accessible
- States matched: 2

**EGFR_170-200_PKT30**
- State class: state_shifted_pocket
- Interpretation: unknown
- Accessibility: mostly_accessible
- States matched: 2

**3GT8_raw_PKT19**
- State class: uncertain_alignment
- Interpretation: unknown
- Accessibility: mostly_accessible
- States matched: 2

**3GT8_raw_PKT15**
- State class: state_specific_pocket
- Interpretation: state_dependent
- Accessibility: single_state_only
- States matched: 1
- Caveat: Only 1 receptor state has pocket data. Cannot determine if state-specific or under-sampled. NOTE: Single-state data. Cross-state persistence unknown. Re-evaluate when additional receptor states have pocket data.

**3GT8_raw_PKT21**
- State class: uncertain_alignment
- Interpretation: unknown
- Accessibility: mostly_accessible
- States matched: 2

**EGFR_170-200_PKT32**
- State class: state_shifted_pocket
- Interpretation: unknown
- Accessibility: mostly_accessible
- States matched: 2

**EGFR_160-185_PKT42**
- State class: state_shifted_pocket
- Interpretation: unknown
- Accessibility: mostly_accessible
- States matched: 2

**EGFR_160-185_PKT44**
- State class: state_shifted_pocket
- Interpretation: unknown
- Accessibility: mostly_accessible
- States matched: 2

**EGFR_160-185_PKT41**
- State class: state_specific_pocket
- Interpretation: state_dependent
- Accessibility: single_state_only
- States matched: 1
- Caveat: Only 1 receptor state has pocket data. Cannot determine if state-specific or under-sampled. NOTE: Single-state data. Cross-state persistence unknown. Re-evaluate when additional receptor states have pocket data.

**EGFR_170-200_PKT35**
- State class: state_shifted_pocket
- Interpretation: unknown
- Accessibility: mostly_accessible
- States matched: 2

**EGFR_170-200_PKT37**
- State class: uncertain_alignment
- Interpretation: unknown
- Accessibility: mostly_accessible
- States matched: 2

## 6. Evidence Provenance Summary

Each candidate's score traces back to concrete upstream evidence:

| Phase | Data Source | Key Fields |
|-------|-----------|------------|
| Phase 1 | PPI patch reference | hotspot residues, robustness, method agreement |
| Phase 2 | Pocket proposal | relationship class, druggability tier, proposal score |
| Phase 3 | Docking evidence | ligand support, pose count, diversity verdict |
| Phase 4 | Scoring + classification | axis scores, mechanistic class, state interpretation |

Full provenance is available in `phase4_expanded_evidence_table.csv` (207 rows, 42 columns).

## 7. Validation and Caveats

### Data Completeness

- Phase 1 PPI evidence: Complete
- Phase 2 pocket evidence: Complete
- Phase 3 docking evidence: 156/207 candidates (partial)

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