# Phase 4 Ranking Method Note

## 1. Scoring Formula

The perturbation score is a weighted sum of 4 axis scores:

```
perturbation_score = sum(axis_score_i * weight_i)
```

### Axis Weights

| Axis | Weight | Rationale |
|------|--------|-----------|
| A1_ppi_interface | 30% | Core biological question: does site overlap MYO1D patch? |
| A2_druggability | 25% | Practical requirement: can a small molecule bind here? |
| A3_perturbation_relevance | 30% | Combined biological + chemical evidence |
| A4_state_robustness | 15% | Supporting evidence: is site consistently accessible? |

## 2. Affinity Domination Prevention

For sites classified as `ligandable_but_ppi_irrelevant` or `uncertain_mechanism`, the combined contribution of affinity-influenced axes (A2 + A3) is capped at 35% of the total score. This ensures that a high-affinity hit at a biologically irrelevant site cannot outrank a moderate-affinity hit at an orthosteric site.

## 3. Ranking Results

| Rank | Pocket | Ligand | Class | Score | A1 | A2 | A3 | A4 |
|------|--------|--------|-------|-------|----|----|----|----|
| 1 | 3GT8_raw_PKT07 | 173940 | interface_rim_modulator | 0.5412 | 0.280 | 0.683 | 0.630 | 0.650 |
| 2 | 3GT8_raw_PKT07 | 97806 | interface_rim_modulator | 0.5412 | 0.280 | 0.683 | 0.630 | 0.650 |
| 3 | 3GT8_raw_PKT07 | VAX-C12_0 | interface_rim_modulator | 0.5412 | 0.280 | 0.683 | 0.630 | 0.650 |
| 4 | EGFR_170-200_PKT34 | 173940 | allosteric_modulator | 0.4921 | 0.253 | 0.645 | 0.525 | 0.650 |
| 5 | EGFR_170-200_PKT34 | 97806 | allosteric_modulator | 0.4921 | 0.253 | 0.645 | 0.525 | 0.650 |
| 6 | EGFR_170-200_PKT34 | VAX-C12_0 | allosteric_modulator | 0.4921 | 0.253 | 0.645 | 0.525 | 0.650 |
| 7 | EGFR_160-185_PKT02 | 173940 | interface_rim_modulator | 0.4327 | 0.310 | 0.195 | 0.644 | 0.650 |
| 8 | EGFR_160-185_PKT02 | 97806 | interface_rim_modulator | 0.4327 | 0.310 | 0.195 | 0.644 | 0.650 |
| 9 | EGFR_160-185_PKT02 | VAX-C12_0 | interface_rim_modulator | 0.4327 | 0.310 | 0.195 | 0.644 | 0.650 |
| 10 | 3GT8_raw_PKT10 | 173940 | allosteric_modulator | 0.4306 | 0.253 | 0.399 | 0.525 | 0.650 |
| 11 | 3GT8_raw_PKT10 | 97806 | allosteric_modulator | 0.4306 | 0.253 | 0.399 | 0.525 | 0.650 |
| 12 | 3GT8_raw_PKT10 | VAX-C12_0 | allosteric_modulator | 0.4306 | 0.253 | 0.399 | 0.525 | 0.650 |
| 13 | EGFR_170-200_PKT17 | 173940 | interface_rim_modulator | 0.4300 | 0.339 | 0.133 | 0.659 | 0.650 |
| 14 | EGFR_170-200_PKT17 | VAX-C12_0 | interface_rim_modulator | 0.4300 | 0.339 | 0.133 | 0.659 | 0.650 |
| 15 | EGFR_170-200_PKT06 | 173940 | interface_rim_modulator | 0.4221 | 0.310 | 0.153 | 0.644 | 0.650 |
| 16 | EGFR_170-200_PKT06 | 97806 | interface_rim_modulator | 0.4221 | 0.310 | 0.153 | 0.644 | 0.650 |
| 17 | EGFR_170-200_PKT06 | VAX-C12_0 | interface_rim_modulator | 0.4221 | 0.310 | 0.153 | 0.644 | 0.650 |
| 18 | EGFR_160-185_PKT16 | 173940 | interface_rim_modulator | 0.4191 | 0.310 | 0.141 | 0.644 | 0.650 |
| 19 | EGFR_160-185_PKT16 | 97806 | interface_rim_modulator | 0.4191 | 0.310 | 0.141 | 0.644 | 0.650 |
| 20 | EGFR_160-185_PKT16 | VAX-C12_0 | interface_rim_modulator | 0.4191 | 0.310 | 0.141 | 0.644 | 0.650 |
| 21 | 3GT8_raw_PKT05 | 173940 | interface_rim_modulator | 0.4163 | 0.306 | 0.135 | 0.644 | 0.650 |
| 22 | 3GT8_raw_PKT05 | 97806 | interface_rim_modulator | 0.4163 | 0.306 | 0.135 | 0.644 | 0.650 |
| 23 | 3GT8_raw_PKT05 | VAX-C12_0 | interface_rim_modulator | 0.4163 | 0.306 | 0.135 | 0.644 | 0.650 |
| 24 | 3GT8_raw_PKT11 | 173940 | interface_rim_modulator | 0.4123 | 0.306 | 0.118 | 0.644 | 0.650 |
| 25 | 3GT8_raw_PKT11 | 97806 | interface_rim_modulator | 0.4123 | 0.306 | 0.118 | 0.644 | 0.650 |
| 26 | 3GT8_raw_PKT11 | VAX-C12_0 | interface_rim_modulator | 0.4123 | 0.306 | 0.118 | 0.644 | 0.650 |
| 27 | 3GT8_raw_PKT01 | 173940 | interface_rim_modulator | 0.4095 | 0.280 | 0.157 | 0.630 | 0.650 |
| 28 | 3GT8_raw_PKT01 | 97806 | interface_rim_modulator | 0.4095 | 0.280 | 0.157 | 0.630 | 0.650 |
| 29 | 3GT8_raw_PKT01 | VAX-C12_0 | interface_rim_modulator | 0.4095 | 0.280 | 0.157 | 0.630 | 0.650 |
| 30 | EGFR_160-185_PKT07 | 97806 | interface_rim_modulator | 0.4091 | 0.281 | 0.153 | 0.630 | 0.650 |
| 31 | EGFR_160-185_PKT07 | VAX-C12_0 | interface_rim_modulator | 0.4091 | 0.281 | 0.153 | 0.630 | 0.650 |
| 32 | 3GT8_raw_PKT02 | 97806 | interface_rim_modulator | 0.4080 | 0.280 | 0.151 | 0.630 | 0.650 |
| 33 | 3GT8_raw_PKT02 | VAX-C12_0 | interface_rim_modulator | 0.4080 | 0.280 | 0.151 | 0.630 | 0.650 |
| 34 | EGFR_160-185_PKT12 | 173940 | interface_rim_modulator | 0.4066 | 0.281 | 0.143 | 0.630 | 0.650 |
| 35 | EGFR_160-185_PKT12 | 97806 | interface_rim_modulator | 0.4066 | 0.281 | 0.143 | 0.630 | 0.650 |
| 36 | EGFR_160-185_PKT12 | VAX-C12_0 | interface_rim_modulator | 0.4066 | 0.281 | 0.143 | 0.630 | 0.650 |
| 37 | EGFR_170-200_PKT10 | 173940 | interface_rim_modulator | 0.4063 | 0.281 | 0.142 | 0.630 | 0.650 |
| 38 | EGFR_170-200_PKT10 | 97806 | interface_rim_modulator | 0.4063 | 0.281 | 0.142 | 0.630 | 0.650 |
| 39 | EGFR_170-200_PKT10 | VAX-C12_0 | interface_rim_modulator | 0.4063 | 0.281 | 0.142 | 0.630 | 0.650 |
| 40 | EGFR_160-185_PKT07 | 173940 | interface_rim_modulator | 0.4061 | 0.281 | 0.153 | 0.620 | 0.650 |
| 41 | EGFR_160-185_PKT15 | 173940 | interface_rim_modulator | 0.4061 | 0.281 | 0.141 | 0.630 | 0.650 |
| 42 | EGFR_160-185_PKT15 | 97806 | interface_rim_modulator | 0.4061 | 0.281 | 0.141 | 0.630 | 0.650 |
| 43 | EGFR_160-185_PKT15 | VAX-C12_0 | interface_rim_modulator | 0.4061 | 0.281 | 0.141 | 0.630 | 0.650 |
| 44 | EGFR_170-200_PKT17 | 97806 | interface_rim_modulator | 0.4060 | 0.339 | 0.133 | 0.579 | 0.650 |
| 45 | EGFR_170-200_PKT12 | 97806 | interface_rim_modulator | 0.4056 | 0.281 | 0.139 | 0.630 | 0.650 |
| 46 | 3GT8_raw_PKT02 | 173940 | interface_rim_modulator | 0.4050 | 0.280 | 0.151 | 0.620 | 0.650 |
| 47 | EGFR_160-185_PKT23 | 173940 | interface_rim_modulator | 0.4030 | 0.281 | 0.129 | 0.630 | 0.650 |
| 48 | EGFR_160-185_PKT23 | 97806 | interface_rim_modulator | 0.4030 | 0.281 | 0.129 | 0.630 | 0.650 |
| 49 | EGFR_160-185_PKT23 | VAX-C12_0 | interface_rim_modulator | 0.4030 | 0.281 | 0.129 | 0.630 | 0.650 |
| 50 | EGFR_170-200_PKT12 | 173940 | interface_rim_modulator | 0.4026 | 0.281 | 0.139 | 0.620 | 0.650 |
| 51 | EGFR_170-200_PKT12 | VAX-C12_0 | interface_rim_modulator | 0.4026 | 0.281 | 0.139 | 0.620 | 0.650 |
| 52 | 3GT8_raw_PKT08 | 173940 | interface_rim_modulator | 0.4009 | 0.306 | 0.123 | 0.644 | 0.567 |
| 53 | 3GT8_raw_PKT14 | 173940 | interface_rim_modulator | 0.3991 | 0.280 | 0.115 | 0.630 | 0.650 |
| 54 | 3GT8_raw_PKT14 | 97806 | interface_rim_modulator | 0.3991 | 0.280 | 0.115 | 0.630 | 0.650 |
| 55 | 3GT8_raw_PKT14 | VAX-C12_0 | interface_rim_modulator | 0.3991 | 0.280 | 0.115 | 0.630 | 0.650 |
| 56 | EGFR_160-185_PKT09 | 173940 | interface_rim_modulator | 0.3962 | 0.281 | 0.151 | 0.630 | 0.567 |
| 57 | EGFR_160-185_PKT09 | VAX-C12_0 | interface_rim_modulator | 0.3962 | 0.281 | 0.151 | 0.630 | 0.567 |
| 58 | 3GT8_raw_PKT20 | 173940 | interface_rim_modulator | 0.3893 | 0.306 | 0.097 | 0.644 | 0.533 |
| 59 | 3GT8_raw_PKT20 | VAX-C12_0 | interface_rim_modulator | 0.3893 | 0.306 | 0.097 | 0.644 | 0.533 |
| 60 | EGFR_170-200_PKT28 | 173940 | interface_rim_modulator | 0.3885 | 0.310 | 0.115 | 0.564 | 0.650 |
| 61 | EGFR_170-200_PKT28 | VAX-C12_0 | interface_rim_modulator | 0.3885 | 0.310 | 0.115 | 0.564 | 0.650 |
| 62 | EGFR_160-185_PKT34 | 173940 | interface_rim_modulator | 0.3848 | 0.281 | 0.106 | 0.630 | 0.567 |
| 63 | EGFR_160-185_PKT34 | VAX-C12_0 | interface_rim_modulator | 0.3848 | 0.281 | 0.106 | 0.630 | 0.567 |
| 64 | 3GT8_raw_PKT08 | 97806 | interface_rim_modulator | 0.3769 | 0.306 | 0.123 | 0.564 | 0.567 |
| 65 | 3GT8_raw_PKT08 | VAX-C12_0 | interface_rim_modulator | 0.3769 | 0.306 | 0.123 | 0.564 | 0.567 |
| 66 | EGFR_160-185_PKT09 | 97806 | interface_rim_modulator | 0.3722 | 0.281 | 0.151 | 0.550 | 0.567 |
| 67 | EGFR_170-200_PKT28 | 97806 | interface_rim_modulator | 0.3705 | 0.310 | 0.115 | 0.504 | 0.650 |
| 68 | 3GT8_raw_PKT20 | 97806 | interface_rim_modulator | 0.3653 | 0.306 | 0.097 | 0.564 | 0.533 |
| 69 | EGFR_160-185_PKT34 | 97806 | interface_rim_modulator | 0.3608 | 0.281 | 0.106 | 0.550 | 0.567 |
| 70 | 3GT8_raw_PKT03 | 173940 | uncertain_mechanism | 0.2668 | 0.253 | 0.149 | 0.475 | 0.650 |
| 71 | 3GT8_raw_PKT03 | 97806 | uncertain_mechanism | 0.2668 | 0.253 | 0.149 | 0.475 | 0.650 |
| 72 | 3GT8_raw_PKT03 | VAX-C12_0 | uncertain_mechanism | 0.2668 | 0.253 | 0.149 | 0.475 | 0.650 |
| 73 | 3GT8_raw_PKT16 | 173940 | uncertain_mechanism | 0.2668 | 0.253 | 0.114 | 0.475 | 0.650 |
| 74 | 3GT8_raw_PKT16 | 97806 | uncertain_mechanism | 0.2668 | 0.253 | 0.114 | 0.475 | 0.650 |
| 75 | 3GT8_raw_PKT16 | VAX-C12_0 | uncertain_mechanism | 0.2668 | 0.253 | 0.114 | 0.475 | 0.650 |
| 76 | EGFR_160-185_PKT01 | 173940 | uncertain_mechanism | 0.2668 | 0.253 | 0.204 | 0.475 | 0.650 |
| 77 | EGFR_160-185_PKT01 | 97806 | uncertain_mechanism | 0.2668 | 0.253 | 0.204 | 0.475 | 0.650 |
| 78 | EGFR_160-185_PKT01 | VAX-C12_0 | uncertain_mechanism | 0.2668 | 0.253 | 0.204 | 0.475 | 0.650 |
| 79 | EGFR_160-185_PKT03 | 173940 | uncertain_mechanism | 0.2668 | 0.253 | 0.180 | 0.475 | 0.650 |
| 80 | EGFR_160-185_PKT03 | 97806 | uncertain_mechanism | 0.2668 | 0.253 | 0.180 | 0.475 | 0.650 |
| 81 | EGFR_160-185_PKT03 | VAX-C12_0 | uncertain_mechanism | 0.2668 | 0.253 | 0.180 | 0.475 | 0.650 |
| 82 | EGFR_160-185_PKT04 | 173940 | uncertain_mechanism | 0.2668 | 0.253 | 0.167 | 0.475 | 0.650 |
| 83 | EGFR_160-185_PKT04 | 97806 | uncertain_mechanism | 0.2668 | 0.253 | 0.167 | 0.475 | 0.650 |
| 84 | EGFR_160-185_PKT04 | VAX-C12_0 | uncertain_mechanism | 0.2668 | 0.253 | 0.167 | 0.475 | 0.650 |
| 85 | EGFR_160-185_PKT05 | 173940 | uncertain_mechanism | 0.2668 | 0.253 | 0.161 | 0.475 | 0.650 |
| 86 | EGFR_160-185_PKT05 | 97806 | uncertain_mechanism | 0.2668 | 0.253 | 0.161 | 0.475 | 0.650 |
| 87 | EGFR_160-185_PKT05 | VAX-C12_0 | uncertain_mechanism | 0.2668 | 0.253 | 0.161 | 0.475 | 0.650 |
| 88 | EGFR_160-185_PKT11 | 173940 | uncertain_mechanism | 0.2668 | 0.253 | 0.144 | 0.395 | 0.650 |
| 89 | EGFR_160-185_PKT11 | 97806 | uncertain_mechanism | 0.2668 | 0.253 | 0.144 | 0.395 | 0.650 |
| 90 | EGFR_160-185_PKT11 | VAX-C12_0 | uncertain_mechanism | 0.2668 | 0.253 | 0.144 | 0.395 | 0.650 |
| 91 | EGFR_160-185_PKT13 | 173940 | uncertain_mechanism | 0.2668 | 0.253 | 0.142 | 0.475 | 0.650 |
| 92 | EGFR_160-185_PKT13 | 97806 | uncertain_mechanism | 0.2668 | 0.253 | 0.142 | 0.385 | 0.650 |
| 93 | EGFR_160-185_PKT13 | VAX-C12_0 | uncertain_mechanism | 0.2668 | 0.253 | 0.142 | 0.475 | 0.650 |
| 94 | EGFR_160-185_PKT19 | 173940 | uncertain_mechanism | 0.2668 | 0.253 | 0.135 | 0.475 | 0.650 |
| 95 | EGFR_160-185_PKT19 | 97806 | uncertain_mechanism | 0.2668 | 0.253 | 0.135 | 0.475 | 0.650 |
| 96 | EGFR_160-185_PKT19 | VAX-C12_0 | uncertain_mechanism | 0.2668 | 0.253 | 0.135 | 0.475 | 0.650 |
| 97 | EGFR_160-185_PKT30 | 173940 | uncertain_mechanism | 0.2668 | 0.253 | 0.114 | 0.475 | 0.650 |
| 98 | EGFR_160-185_PKT30 | 97806 | uncertain_mechanism | 0.2668 | 0.253 | 0.114 | 0.475 | 0.650 |
| 99 | EGFR_160-185_PKT30 | VAX-C12_0 | uncertain_mechanism | 0.2668 | 0.253 | 0.114 | 0.475 | 0.650 |
| 100 | EGFR_160-185_PKT37 | 173940 | uncertain_mechanism | 0.2668 | 0.253 | 0.104 | 0.475 | 0.650 |
| 101 | EGFR_160-185_PKT37 | 97806 | uncertain_mechanism | 0.2668 | 0.253 | 0.104 | 0.475 | 0.650 |
| 102 | EGFR_160-185_PKT37 | VAX-C12_0 | uncertain_mechanism | 0.2668 | 0.253 | 0.104 | 0.475 | 0.650 |
| 103 | EGFR_160-185_PKT38 | 173940 | uncertain_mechanism | 0.2668 | 0.253 | 0.104 | 0.475 | 0.650 |
| 104 | EGFR_160-185_PKT38 | 97806 | uncertain_mechanism | 0.2668 | 0.253 | 0.104 | 0.475 | 0.650 |
| 105 | EGFR_160-185_PKT38 | VAX-C12_0 | uncertain_mechanism | 0.2668 | 0.253 | 0.104 | 0.475 | 0.650 |
| 106 | EGFR_170-200_PKT01 | 173940 | uncertain_mechanism | 0.2668 | 0.253 | 0.218 | 0.475 | 0.650 |
| 107 | EGFR_170-200_PKT01 | 97806 | uncertain_mechanism | 0.2668 | 0.253 | 0.218 | 0.475 | 0.650 |
| 108 | EGFR_170-200_PKT01 | VAX-C12_0 | uncertain_mechanism | 0.2668 | 0.253 | 0.218 | 0.475 | 0.650 |
| 109 | EGFR_170-200_PKT04 | 173940 | uncertain_mechanism | 0.2668 | 0.253 | 0.161 | 0.295 | 0.650 |
| 110 | EGFR_170-200_PKT04 | 97806 | uncertain_mechanism | 0.2668 | 0.253 | 0.161 | 0.295 | 0.650 |
| 111 | EGFR_170-200_PKT04 | VAX-C12_0 | uncertain_mechanism | 0.2668 | 0.253 | 0.161 | 0.305 | 0.650 |
| 112 | EGFR_170-200_PKT05 | 173940 | uncertain_mechanism | 0.2668 | 0.253 | 0.157 | 0.475 | 0.650 |
| 113 | EGFR_170-200_PKT05 | 97806 | uncertain_mechanism | 0.2668 | 0.253 | 0.157 | 0.475 | 0.650 |
| 114 | EGFR_170-200_PKT05 | VAX-C12_0 | uncertain_mechanism | 0.2668 | 0.253 | 0.157 | 0.475 | 0.650 |
| 115 | EGFR_170-200_PKT15 | 173940 | uncertain_mechanism | 0.2668 | 0.253 | 0.134 | 0.475 | 0.650 |
| 116 | EGFR_170-200_PKT15 | 97806 | uncertain_mechanism | 0.2668 | 0.253 | 0.134 | 0.475 | 0.650 |
| 117 | EGFR_170-200_PKT15 | VAX-C12_0 | uncertain_mechanism | 0.2668 | 0.253 | 0.134 | 0.475 | 0.650 |
| 118 | EGFR_170-200_PKT25 | 173940 | uncertain_mechanism | 0.2668 | 0.253 | 0.116 | 0.475 | 0.650 |
| 119 | EGFR_170-200_PKT25 | 97806 | uncertain_mechanism | 0.2668 | 0.253 | 0.116 | 0.475 | 0.650 |
| 120 | EGFR_170-200_PKT25 | VAX-C12_0 | uncertain_mechanism | 0.2668 | 0.253 | 0.116 | 0.475 | 0.650 |
| 121 | EGFR_170-200_PKT31 | 173940 | uncertain_mechanism | 0.2668 | 0.253 | 0.100 | 0.475 | 0.650 |
| 122 | EGFR_170-200_PKT31 | 97806 | uncertain_mechanism | 0.2668 | 0.253 | 0.100 | 0.475 | 0.650 |
| 123 | EGFR_170-200_PKT31 | VAX-C12_0 | uncertain_mechanism | 0.2668 | 0.253 | 0.100 | 0.475 | 0.650 |
| 124 | EGFR_170-200_PKT36 | 173940 | uncertain_mechanism | 0.2668 | 0.253 | 0.044 | 0.475 | 0.650 |
| 125 | EGFR_170-200_PKT36 | 97806 | uncertain_mechanism | 0.2668 | 0.253 | 0.044 | 0.475 | 0.650 |
| 126 | EGFR_170-200_PKT36 | VAX-C12_0 | uncertain_mechanism | 0.2668 | 0.253 | 0.044 | 0.475 | 0.650 |
| 127 | 3GT8_raw_PKT22 | 173940 | uncertain_mechanism | 0.2475 | 0.253 | 0.024 | 0.475 | 0.567 |
| 128 | 3GT8_raw_PKT22 | 97806 | uncertain_mechanism | 0.2475 | 0.253 | 0.024 | 0.475 | 0.567 |
| 129 | 3GT8_raw_PKT22 | VAX-C12_0 | uncertain_mechanism | 0.2475 | 0.253 | 0.024 | 0.475 | 0.567 |
| 130 | EGFR_160-185_PKT14 | 173940 | uncertain_mechanism | 0.2475 | 0.253 | 0.141 | 0.475 | 0.567 |
| 131 | EGFR_160-185_PKT14 | 97806 | uncertain_mechanism | 0.2475 | 0.253 | 0.141 | 0.475 | 0.567 |
| 132 | EGFR_160-185_PKT14 | VAX-C12_0 | uncertain_mechanism | 0.2475 | 0.253 | 0.141 | 0.475 | 0.567 |
| 133 | EGFR_160-185_PKT27 | 173940 | uncertain_mechanism | 0.2475 | 0.253 | 0.121 | 0.475 | 0.567 |
| 134 | EGFR_160-185_PKT27 | 97806 | uncertain_mechanism | 0.2475 | 0.253 | 0.121 | 0.475 | 0.567 |
| 135 | EGFR_160-185_PKT27 | VAX-C12_0 | uncertain_mechanism | 0.2475 | 0.253 | 0.121 | 0.475 | 0.567 |
| 136 | EGFR_160-185_PKT29 | 173940 | uncertain_mechanism | 0.2475 | 0.253 | 0.117 | 0.475 | 0.567 |
| 137 | EGFR_160-185_PKT29 | 97806 | uncertain_mechanism | 0.2475 | 0.253 | 0.117 | 0.475 | 0.567 |
| 138 | EGFR_160-185_PKT29 | VAX-C12_0 | uncertain_mechanism | 0.2475 | 0.253 | 0.117 | 0.475 | 0.567 |
| 139 | EGFR_160-185_PKT39 | 173940 | uncertain_mechanism | 0.2475 | 0.253 | 0.104 | 0.475 | 0.567 |
| 140 | EGFR_160-185_PKT39 | 97806 | uncertain_mechanism | 0.2475 | 0.253 | 0.104 | 0.475 | 0.567 |
| 141 | EGFR_160-185_PKT39 | VAX-C12_0 | uncertain_mechanism | 0.2475 | 0.253 | 0.104 | 0.475 | 0.567 |
| 142 | EGFR_160-185_PKT40 | 173940 | uncertain_mechanism | 0.2475 | 0.253 | 0.094 | 0.475 | 0.567 |
| 143 | EGFR_160-185_PKT40 | 97806 | uncertain_mechanism | 0.2475 | 0.253 | 0.094 | 0.385 | 0.567 |
| 144 | EGFR_160-185_PKT40 | VAX-C12_0 | uncertain_mechanism | 0.2475 | 0.253 | 0.094 | 0.475 | 0.567 |
| 145 | EGFR_160-185_PKT43 | 173940 | uncertain_mechanism | 0.2475 | 0.253 | 0.088 | 0.475 | 0.567 |
| 146 | EGFR_160-185_PKT43 | 97806 | uncertain_mechanism | 0.2475 | 0.253 | 0.088 | 0.475 | 0.567 |
| 147 | EGFR_160-185_PKT43 | VAX-C12_0 | uncertain_mechanism | 0.2475 | 0.253 | 0.088 | 0.475 | 0.567 |
| 148 | EGFR_170-200_PKT13 | 173940 | uncertain_mechanism | 0.2475 | 0.253 | 0.139 | 0.475 | 0.567 |
| 149 | EGFR_170-200_PKT13 | 97806 | uncertain_mechanism | 0.2475 | 0.253 | 0.139 | 0.475 | 0.567 |
| 150 | EGFR_170-200_PKT13 | VAX-C12_0 | uncertain_mechanism | 0.2475 | 0.253 | 0.139 | 0.475 | 0.567 |
| 151 | EGFR_170-200_PKT26 | 173940 | uncertain_mechanism | 0.2475 | 0.253 | 0.116 | 0.395 | 0.567 |
| 152 | EGFR_170-200_PKT26 | 97806 | uncertain_mechanism | 0.2475 | 0.253 | 0.116 | 0.395 | 0.567 |
| 153 | EGFR_170-200_PKT26 | VAX-C12_0 | uncertain_mechanism | 0.2475 | 0.253 | 0.116 | 0.475 | 0.567 |
| 154 | EGFR_170-200_PKT33 | 173940 | uncertain_mechanism | 0.2475 | 0.253 | 0.090 | 0.475 | 0.567 |
| 155 | EGFR_170-200_PKT33 | 97806 | uncertain_mechanism | 0.2475 | 0.253 | 0.090 | 0.475 | 0.567 |
| 156 | EGFR_170-200_PKT33 | VAX-C12_0 | uncertain_mechanism | 0.2475 | 0.253 | 0.090 | 0.475 | 0.567 |
| 157 | EGFR_170-200_PKT03 | (none) | ligandable_but_ppi_irrelevant | 0.1902 | 0.253 | 0.163 | 0.045 | 0.400 |
| 158 | EGFR_160-185_PKT06 | (none) | ligandable_but_ppi_irrelevant | 0.1878 | 0.253 | 0.153 | 0.045 | 0.400 |
| 159 | EGFR_170-200_PKT07 | (none) | ligandable_but_ppi_irrelevant | 0.1868 | 0.253 | 0.150 | 0.045 | 0.400 |
| 160 | EGFR_160-185_PKT17 | (none) | ligandable_but_ppi_irrelevant | 0.1836 | 0.253 | 0.137 | 0.045 | 0.400 |
| 161 | EGFR_170-200_PKT18 | (none) | ligandable_but_ppi_irrelevant | 0.1818 | 0.253 | 0.130 | 0.045 | 0.400 |
| 162 | 3GT8_raw_PKT13 | (none) | ligandable_but_ppi_irrelevant | 0.1786 | 0.253 | 0.117 | 0.045 | 0.400 |
| 163 | EGFR_170-200_PKT27 | (none) | ligandable_but_ppi_irrelevant | 0.1782 | 0.253 | 0.115 | 0.045 | 0.400 |
| 164 | EGFR_170-200_PKT29 | (none) | ligandable_but_ppi_irrelevant | 0.1773 | 0.253 | 0.112 | 0.045 | 0.400 |
| 165 | 3GT8_raw_PKT18 | (none) | ligandable_but_ppi_irrelevant | 0.1760 | 0.253 | 0.106 | 0.045 | 0.400 |
| 166 | EGFR_160-185_PKT08 | (none) | ligandable_but_ppi_irrelevant | 0.1749 | 0.253 | 0.152 | 0.045 | 0.317 |
| 167 | EGFR_170-200_PKT02 | (none) | ligandable_but_ppi_irrelevant | 0.1749 | 0.253 | 0.172 | 0.045 | 0.283 |
| 168 | EGFR_160-185_PKT10 | (none) | ligandable_but_ppi_irrelevant | 0.1746 | 0.253 | 0.151 | 0.045 | 0.317 |
| 169 | 3GT8_raw_PKT04 | (none) | ligandable_but_ppi_irrelevant | 0.1735 | 0.253 | 0.146 | 0.045 | 0.317 |
| 170 | EGFR_170-200_PKT08 | (none) | ligandable_but_ppi_irrelevant | 0.1733 | 0.253 | 0.146 | 0.045 | 0.317 |
| 171 | EGFR_170-200_PKT09 | (none) | ligandable_but_ppi_irrelevant | 0.1730 | 0.253 | 0.144 | 0.045 | 0.317 |
| 172 | EGFR_170-200_PKT11 | (none) | ligandable_but_ppi_irrelevant | 0.1718 | 0.253 | 0.140 | 0.045 | 0.317 |
| 173 | EGFR_170-200_PKT14 | (none) | ligandable_but_ppi_irrelevant | 0.1711 | 0.253 | 0.137 | 0.045 | 0.317 |
| 174 | EGFR_160-185_PKT18 | (none) | ligandable_but_ppi_irrelevant | 0.1709 | 0.253 | 0.136 | 0.045 | 0.317 |
| 175 | EGFR_160-185_PKT20 | (none) | ligandable_but_ppi_irrelevant | 0.1703 | 0.253 | 0.133 | 0.045 | 0.317 |
| 176 | EGFR_160-185_PKT21 | (none) | ligandable_but_ppi_irrelevant | 0.1702 | 0.253 | 0.133 | 0.045 | 0.317 |
| 177 | EGFR_170-200_PKT16 | (none) | ligandable_but_ppi_irrelevant | 0.1702 | 0.253 | 0.133 | 0.045 | 0.317 |
| 178 | EGFR_160-185_PKT22 | (none) | ligandable_but_ppi_irrelevant | 0.1697 | 0.253 | 0.131 | 0.045 | 0.317 |
| 179 | EGFR_160-185_PKT24 | (none) | ligandable_but_ppi_irrelevant | 0.1690 | 0.253 | 0.128 | 0.045 | 0.317 |
| 180 | EGFR_160-185_PKT25 | (none) | ligandable_but_ppi_irrelevant | 0.1689 | 0.253 | 0.128 | 0.045 | 0.317 |
| 181 | EGFR_170-200_PKT19 | (none) | ligandable_but_ppi_irrelevant | 0.1685 | 0.253 | 0.126 | 0.045 | 0.317 |
| 182 | EGFR_160-185_PKT26 | (none) | ligandable_but_ppi_irrelevant | 0.1676 | 0.253 | 0.123 | 0.045 | 0.317 |
| 183 | EGFR_170-200_PKT20 | (none) | ligandable_but_ppi_irrelevant | 0.1673 | 0.253 | 0.122 | 0.045 | 0.317 |
| 184 | EGFR_170-200_PKT21 | (none) | ligandable_but_ppi_irrelevant | 0.1672 | 0.253 | 0.121 | 0.045 | 0.317 |
| 185 | 3GT8_raw_PKT09 | (none) | ligandable_but_ppi_irrelevant | 0.1671 | 0.253 | 0.121 | 0.045 | 0.317 |
| 186 | EGFR_170-200_PKT22 | (none) | ligandable_but_ppi_irrelevant | 0.1671 | 0.253 | 0.121 | 0.045 | 0.317 |
| 187 | EGFR_170-200_PKT23 | (none) | ligandable_but_ppi_irrelevant | 0.1671 | 0.253 | 0.121 | 0.045 | 0.317 |
| 188 | EGFR_170-200_PKT24 | (none) | ligandable_but_ppi_irrelevant | 0.1670 | 0.253 | 0.120 | 0.045 | 0.317 |
| 189 | EGFR_160-185_PKT28 | (none) | ligandable_but_ppi_irrelevant | 0.1667 | 0.253 | 0.119 | 0.045 | 0.317 |
| 190 | 3GT8_raw_PKT12 | (none) | ligandable_but_ppi_irrelevant | 0.1662 | 0.253 | 0.117 | 0.045 | 0.317 |
| 191 | EGFR_160-185_PKT31 | (none) | ligandable_but_ppi_irrelevant | 0.1652 | 0.253 | 0.113 | 0.045 | 0.317 |
| 192 | EGFR_160-185_PKT32 | (none) | ligandable_but_ppi_irrelevant | 0.1652 | 0.253 | 0.113 | 0.045 | 0.317 |
| 193 | 3GT8_raw_PKT06 | (none) | ligandable_but_ppi_irrelevant | 0.1650 | 0.253 | 0.133 | 0.045 | 0.283 |
| 194 | EGFR_160-185_PKT33 | (none) | ligandable_but_ppi_irrelevant | 0.1649 | 0.253 | 0.112 | 0.045 | 0.317 |
| 195 | 3GT8_raw_PKT17 | (none) | ligandable_but_ppi_irrelevant | 0.1640 | 0.253 | 0.108 | 0.045 | 0.317 |
| 196 | EGFR_160-185_PKT35 | (none) | ligandable_but_ppi_irrelevant | 0.1634 | 0.253 | 0.106 | 0.045 | 0.317 |
| 197 | EGFR_160-185_PKT36 | (none) | ligandable_but_ppi_irrelevant | 0.1634 | 0.253 | 0.106 | 0.045 | 0.317 |
| 198 | EGFR_170-200_PKT30 | (none) | ligandable_but_ppi_irrelevant | 0.1625 | 0.253 | 0.102 | 0.045 | 0.317 |
| 199 | 3GT8_raw_PKT19 | (none) | ligandable_but_ppi_irrelevant | 0.1623 | 0.253 | 0.102 | 0.045 | 0.317 |
| 200 | 3GT8_raw_PKT15 | (none) | ligandable_but_ppi_irrelevant | 0.1607 | 0.253 | 0.115 | 0.045 | 0.283 |
| 201 | 3GT8_raw_PKT21 | (none) | ligandable_but_ppi_irrelevant | 0.1601 | 0.253 | 0.093 | 0.045 | 0.317 |
| 202 | EGFR_170-200_PKT32 | (none) | ligandable_but_ppi_irrelevant | 0.1601 | 0.253 | 0.093 | 0.045 | 0.317 |
| 203 | EGFR_160-185_PKT42 | (none) | ligandable_but_ppi_irrelevant | 0.1590 | 0.253 | 0.088 | 0.045 | 0.317 |
| 204 | EGFR_160-185_PKT44 | (none) | ligandable_but_ppi_irrelevant | 0.1581 | 0.253 | 0.085 | 0.045 | 0.317 |
| 205 | EGFR_160-185_PKT41 | (none) | ligandable_but_ppi_irrelevant | 0.1551 | 0.253 | 0.093 | 0.045 | 0.283 |
| 206 | EGFR_170-200_PKT35 | (none) | ligandable_but_ppi_irrelevant | 0.1510 | 0.253 | 0.056 | 0.045 | 0.317 |
| 207 | EGFR_170-200_PKT37 | (none) | ligandable_but_ppi_irrelevant | 0.1478 | 0.253 | 0.044 | 0.045 | 0.317 |

### Affinity Cap Applied To

- 3GT8_raw_PKT03: capped score = 0.2668
- 3GT8_raw_PKT03: capped score = 0.2668
- 3GT8_raw_PKT03: capped score = 0.2668
- 3GT8_raw_PKT16: capped score = 0.2668
- 3GT8_raw_PKT16: capped score = 0.2668
- 3GT8_raw_PKT16: capped score = 0.2668
- EGFR_160-185_PKT01: capped score = 0.2668
- EGFR_160-185_PKT01: capped score = 0.2668
- EGFR_160-185_PKT01: capped score = 0.2668
- EGFR_160-185_PKT03: capped score = 0.2668
- EGFR_160-185_PKT03: capped score = 0.2668
- EGFR_160-185_PKT03: capped score = 0.2668
- EGFR_160-185_PKT04: capped score = 0.2668
- EGFR_160-185_PKT04: capped score = 0.2668
- EGFR_160-185_PKT04: capped score = 0.2668
- EGFR_160-185_PKT05: capped score = 0.2668
- EGFR_160-185_PKT05: capped score = 0.2668
- EGFR_160-185_PKT05: capped score = 0.2668
- EGFR_160-185_PKT11: capped score = 0.2668
- EGFR_160-185_PKT11: capped score = 0.2668
- EGFR_160-185_PKT11: capped score = 0.2668
- EGFR_160-185_PKT13: capped score = 0.2668
- EGFR_160-185_PKT13: capped score = 0.2668
- EGFR_160-185_PKT13: capped score = 0.2668
- EGFR_160-185_PKT19: capped score = 0.2668
- EGFR_160-185_PKT19: capped score = 0.2668
- EGFR_160-185_PKT19: capped score = 0.2668
- EGFR_160-185_PKT30: capped score = 0.2668
- EGFR_160-185_PKT30: capped score = 0.2668
- EGFR_160-185_PKT30: capped score = 0.2668
- EGFR_160-185_PKT37: capped score = 0.2668
- EGFR_160-185_PKT37: capped score = 0.2668
- EGFR_160-185_PKT37: capped score = 0.2668
- EGFR_160-185_PKT38: capped score = 0.2668
- EGFR_160-185_PKT38: capped score = 0.2668
- EGFR_160-185_PKT38: capped score = 0.2668
- EGFR_170-200_PKT01: capped score = 0.2668
- EGFR_170-200_PKT01: capped score = 0.2668
- EGFR_170-200_PKT01: capped score = 0.2668
- EGFR_170-200_PKT04: capped score = 0.2668
- EGFR_170-200_PKT04: capped score = 0.2668
- EGFR_170-200_PKT04: capped score = 0.2668
- EGFR_170-200_PKT05: capped score = 0.2668
- EGFR_170-200_PKT05: capped score = 0.2668
- EGFR_170-200_PKT05: capped score = 0.2668
- EGFR_170-200_PKT15: capped score = 0.2668
- EGFR_170-200_PKT15: capped score = 0.2668
- EGFR_170-200_PKT15: capped score = 0.2668
- EGFR_170-200_PKT25: capped score = 0.2668
- EGFR_170-200_PKT25: capped score = 0.2668
- EGFR_170-200_PKT25: capped score = 0.2668
- EGFR_170-200_PKT31: capped score = 0.2668
- EGFR_170-200_PKT31: capped score = 0.2668
- EGFR_170-200_PKT31: capped score = 0.2668
- EGFR_170-200_PKT36: capped score = 0.2668
- EGFR_170-200_PKT36: capped score = 0.2668
- EGFR_170-200_PKT36: capped score = 0.2668
- 3GT8_raw_PKT22: capped score = 0.2475
- 3GT8_raw_PKT22: capped score = 0.2475
- 3GT8_raw_PKT22: capped score = 0.2475
- EGFR_160-185_PKT14: capped score = 0.2475
- EGFR_160-185_PKT14: capped score = 0.2475
- EGFR_160-185_PKT14: capped score = 0.2475
- EGFR_160-185_PKT27: capped score = 0.2475
- EGFR_160-185_PKT27: capped score = 0.2475
- EGFR_160-185_PKT27: capped score = 0.2475
- EGFR_160-185_PKT29: capped score = 0.2475
- EGFR_160-185_PKT29: capped score = 0.2475
- EGFR_160-185_PKT29: capped score = 0.2475
- EGFR_160-185_PKT39: capped score = 0.2475
- EGFR_160-185_PKT39: capped score = 0.2475
- EGFR_160-185_PKT39: capped score = 0.2475
- EGFR_160-185_PKT40: capped score = 0.2475
- EGFR_160-185_PKT40: capped score = 0.2475
- EGFR_160-185_PKT40: capped score = 0.2475
- EGFR_160-185_PKT43: capped score = 0.2475
- EGFR_160-185_PKT43: capped score = 0.2475
- EGFR_160-185_PKT43: capped score = 0.2475
- EGFR_170-200_PKT13: capped score = 0.2475
- EGFR_170-200_PKT13: capped score = 0.2475
- EGFR_170-200_PKT13: capped score = 0.2475
- EGFR_170-200_PKT26: capped score = 0.2475
- EGFR_170-200_PKT26: capped score = 0.2475
- EGFR_170-200_PKT26: capped score = 0.2475
- EGFR_170-200_PKT33: capped score = 0.2475
- EGFR_170-200_PKT33: capped score = 0.2475
- EGFR_170-200_PKT33: capped score = 0.2475
- EGFR_170-200_PKT03: capped score = 0.1902
- EGFR_160-185_PKT06: capped score = 0.1878
- EGFR_170-200_PKT07: capped score = 0.1868
- EGFR_160-185_PKT17: capped score = 0.1836
- EGFR_170-200_PKT18: capped score = 0.1818
- 3GT8_raw_PKT13: capped score = 0.1786
- EGFR_170-200_PKT27: capped score = 0.1782
- EGFR_170-200_PKT29: capped score = 0.1773
- 3GT8_raw_PKT18: capped score = 0.1760
- EGFR_160-185_PKT08: capped score = 0.1749
- EGFR_170-200_PKT02: capped score = 0.1749
- EGFR_160-185_PKT10: capped score = 0.1746
- 3GT8_raw_PKT04: capped score = 0.1735
- EGFR_170-200_PKT08: capped score = 0.1733
- EGFR_170-200_PKT09: capped score = 0.1730
- EGFR_170-200_PKT11: capped score = 0.1718
- EGFR_170-200_PKT14: capped score = 0.1711
- EGFR_160-185_PKT18: capped score = 0.1709
- EGFR_160-185_PKT20: capped score = 0.1703
- EGFR_160-185_PKT21: capped score = 0.1702
- EGFR_170-200_PKT16: capped score = 0.1702
- EGFR_160-185_PKT22: capped score = 0.1697
- EGFR_160-185_PKT24: capped score = 0.1690
- EGFR_160-185_PKT25: capped score = 0.1689
- EGFR_170-200_PKT19: capped score = 0.1685
- EGFR_160-185_PKT26: capped score = 0.1676
- EGFR_170-200_PKT20: capped score = 0.1673
- EGFR_170-200_PKT21: capped score = 0.1672
- 3GT8_raw_PKT09: capped score = 0.1671
- EGFR_170-200_PKT22: capped score = 0.1671
- EGFR_170-200_PKT23: capped score = 0.1671
- EGFR_170-200_PKT24: capped score = 0.1670
- EGFR_160-185_PKT28: capped score = 0.1667
- 3GT8_raw_PKT12: capped score = 0.1662
- EGFR_160-185_PKT31: capped score = 0.1652
- EGFR_160-185_PKT32: capped score = 0.1652
- 3GT8_raw_PKT06: capped score = 0.1650
- EGFR_160-185_PKT33: capped score = 0.1649
- 3GT8_raw_PKT17: capped score = 0.1640
- EGFR_160-185_PKT35: capped score = 0.1634
- EGFR_160-185_PKT36: capped score = 0.1634
- EGFR_170-200_PKT30: capped score = 0.1625
- 3GT8_raw_PKT19: capped score = 0.1623
- 3GT8_raw_PKT15: capped score = 0.1607
- 3GT8_raw_PKT21: capped score = 0.1601
- EGFR_170-200_PKT32: capped score = 0.1601
- EGFR_160-185_PKT42: capped score = 0.1590
- EGFR_160-185_PKT44: capped score = 0.1581
- EGFR_160-185_PKT41: capped score = 0.1551
- EGFR_170-200_PKT35: capped score = 0.1510
- EGFR_170-200_PKT37: capped score = 0.1478

## 4. Provenance

Each ranked candidate preserves references to:
- Phase 1: hotspot overlap count/fraction
- Phase 2: relationship class, druggability tier
- Phase 3: ligand support strength, pose count, best affinity
- Phase 4: mechanistic class, classification basis

These fields are available in `perturbation_candidate_table.csv` for full traceability.

---

Generated by `egfr_pipeline.phase4.perturbation_scoring`