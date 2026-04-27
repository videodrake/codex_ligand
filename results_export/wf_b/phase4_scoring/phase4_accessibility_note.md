# Phase 4 State-Robustness and Accessibility Note

## 1. State Interpretation Framework

| State Class | Interpretation | Confidence Effect | Description |
|-------------|---------------|-------------------|-------------|
| robust_pocket | persistent | strengthen | Pocket is present across all tested receptor states. High confidence i... |
| shifted_pocket | conditionally_accessible | neutral | Pocket geometry shifts between states but remains present. May require... |
| state_specific_pocket | state_dependent | flag | Pocket detected in only one conformational state. Could indicate a cry... |

## 2. Design Principles

- **State-specific sites are NOT auto-penalized**: They may represent cryptic pockets or under-sampled conformations.
- **Robust sites get confidence boost**: Persistence across states strengthens the evidence.
- **Single-state data is flagged**: When only 1 receptor state has pocket data, all pockets are inherently state_specific. This is a data limitation, not a biological conclusion.

## 3. Interpretation Results

| Rank | Pocket | Accessibility | Adjusted Confidence | Caveat |
|------|--------|--------------|--------------------|---------| 
| 1 | 3GT8_raw_PKT07 / 173940 | always_accessible | low |  |
| 2 | 3GT8_raw_PKT07 / 97806 | always_accessible | low |  |
| 3 | 3GT8_raw_PKT07 / VAX-C12_0 | always_accessible | low |  |
| 4 | EGFR_170-200_PKT34 / 173940 | always_accessible | medium |  |
| 5 | EGFR_170-200_PKT34 / 97806 | always_accessible | medium |  |
| 6 | EGFR_170-200_PKT34 / VAX-C12_0 | always_accessible | medium |  |
| 7 | EGFR_160-185_PKT02 / 173940 | always_accessible | medium |  |
| 8 | EGFR_160-185_PKT02 / 97806 | always_accessible | medium |  |
| 9 | EGFR_160-185_PKT02 / VAX-C12_0 | always_accessible | medium |  |
| 10 | 3GT8_raw_PKT10 / 173940 | always_accessible | low |  |
| 11 | 3GT8_raw_PKT10 / 97806 | always_accessible | low |  |
| 12 | 3GT8_raw_PKT10 / VAX-C12_0 | always_accessible | low |  |
| 13 | EGFR_170-200_PKT17 / 173940 | always_accessible | medium |  |
| 14 | EGFR_170-200_PKT17 / VAX-C12_0 | always_accessible | medium |  |
| 15 | EGFR_170-200_PKT06 / 173940 | always_accessible | medium |  |
| 16 | EGFR_170-200_PKT06 / 97806 | always_accessible | medium |  |
| 17 | EGFR_170-200_PKT06 / VAX-C12_0 | always_accessible | medium |  |
| 18 | EGFR_160-185_PKT16 / 173940 | always_accessible | medium |  |
| 19 | EGFR_160-185_PKT16 / 97806 | always_accessible | medium |  |
| 20 | EGFR_160-185_PKT16 / VAX-C12_0 | always_accessible | medium |  |
| 21 | 3GT8_raw_PKT05 / 173940 | always_accessible | medium |  |
| 22 | 3GT8_raw_PKT05 / 97806 | always_accessible | medium |  |
| 23 | 3GT8_raw_PKT05 / VAX-C12_0 | always_accessible | medium |  |
| 24 | 3GT8_raw_PKT11 / 173940 | always_accessible | medium |  |
| 25 | 3GT8_raw_PKT11 / 97806 | always_accessible | medium |  |
| 26 | 3GT8_raw_PKT11 / VAX-C12_0 | always_accessible | medium |  |
| 27 | 3GT8_raw_PKT01 / 173940 | always_accessible | low |  |
| 28 | 3GT8_raw_PKT01 / 97806 | always_accessible | low |  |
| 29 | 3GT8_raw_PKT01 / VAX-C12_0 | always_accessible | low |  |
| 30 | EGFR_160-185_PKT07 / 97806 | always_accessible | low |  |
| 31 | EGFR_160-185_PKT07 / VAX-C12_0 | always_accessible | low |  |
| 32 | 3GT8_raw_PKT02 / 97806 | always_accessible | low |  |
| 33 | 3GT8_raw_PKT02 / VAX-C12_0 | always_accessible | low |  |
| 34 | EGFR_160-185_PKT12 / 173940 | always_accessible | low |  |
| 35 | EGFR_160-185_PKT12 / 97806 | always_accessible | low |  |
| 36 | EGFR_160-185_PKT12 / VAX-C12_0 | always_accessible | low |  |
| 37 | EGFR_170-200_PKT10 / 173940 | always_accessible | low |  |
| 38 | EGFR_170-200_PKT10 / 97806 | always_accessible | low |  |
| 39 | EGFR_170-200_PKT10 / VAX-C12_0 | always_accessible | low |  |
| 40 | EGFR_160-185_PKT07 / 173940 | always_accessible | low |  |
| 41 | EGFR_160-185_PKT15 / 173940 | always_accessible | low |  |
| 42 | EGFR_160-185_PKT15 / 97806 | always_accessible | low |  |
| 43 | EGFR_160-185_PKT15 / VAX-C12_0 | always_accessible | low |  |
| 44 | EGFR_170-200_PKT17 / 97806 | always_accessible | medium |  |
| 45 | EGFR_170-200_PKT12 / 97806 | always_accessible | low |  |
| 46 | 3GT8_raw_PKT02 / 173940 | always_accessible | low |  |
| 47 | EGFR_160-185_PKT23 / 173940 | always_accessible | low |  |
| 48 | EGFR_160-185_PKT23 / 97806 | always_accessible | low |  |
| 49 | EGFR_160-185_PKT23 / VAX-C12_0 | always_accessible | low |  |
| 50 | EGFR_170-200_PKT12 / 173940 | always_accessible | low |  |
| 51 | EGFR_170-200_PKT12 / VAX-C12_0 | always_accessible | low |  |
| 52 | 3GT8_raw_PKT08 / 173940 | mostly_accessible | medium |  |
| 53 | 3GT8_raw_PKT14 / 173940 | always_accessible | low |  |
| 54 | 3GT8_raw_PKT14 / 97806 | always_accessible | low |  |
| 55 | 3GT8_raw_PKT14 / VAX-C12_0 | always_accessible | low |  |
| 56 | EGFR_160-185_PKT09 / 173940 | mostly_accessible | low |  |
| 57 | EGFR_160-185_PKT09 / VAX-C12_0 | mostly_accessible | low |  |
| 58 | 3GT8_raw_PKT20 / 173940 | single_state_only | medium_provisional | Only 1 receptor state has pocket data. Cannot determine if s... |
| 59 | 3GT8_raw_PKT20 / VAX-C12_0 | single_state_only | medium_provisional | Only 1 receptor state has pocket data. Cannot determine if s... |
| 60 | EGFR_170-200_PKT28 / 173940 | always_accessible | medium |  |
| 61 | EGFR_170-200_PKT28 / VAX-C12_0 | always_accessible | medium |  |
| 62 | EGFR_160-185_PKT34 / 173940 | mostly_accessible | low |  |
| 63 | EGFR_160-185_PKT34 / VAX-C12_0 | mostly_accessible | low |  |
| 64 | 3GT8_raw_PKT08 / 97806 | mostly_accessible | medium |  |
| 65 | 3GT8_raw_PKT08 / VAX-C12_0 | mostly_accessible | medium |  |
| 66 | EGFR_160-185_PKT09 / 97806 | mostly_accessible | low |  |
| 67 | EGFR_170-200_PKT28 / 97806 | always_accessible | medium |  |
| 68 | 3GT8_raw_PKT20 / 97806 | single_state_only | medium_provisional | Only 1 receptor state has pocket data. Cannot determine if s... |
| 69 | EGFR_160-185_PKT34 / 97806 | mostly_accessible | low |  |
| 70 | 3GT8_raw_PKT03 / 173940 | always_accessible | low |  |
| 71 | 3GT8_raw_PKT03 / 97806 | always_accessible | low |  |
| 72 | 3GT8_raw_PKT03 / VAX-C12_0 | always_accessible | low |  |
| 73 | 3GT8_raw_PKT16 / 173940 | always_accessible | low |  |
| 74 | 3GT8_raw_PKT16 / 97806 | always_accessible | low |  |
| 75 | 3GT8_raw_PKT16 / VAX-C12_0 | always_accessible | low |  |
| 76 | EGFR_160-185_PKT01 / 173940 | always_accessible | low |  |
| 77 | EGFR_160-185_PKT01 / 97806 | always_accessible | low |  |
| 78 | EGFR_160-185_PKT01 / VAX-C12_0 | always_accessible | low |  |
| 79 | EGFR_160-185_PKT03 / 173940 | always_accessible | low |  |
| 80 | EGFR_160-185_PKT03 / 97806 | always_accessible | low |  |
| 81 | EGFR_160-185_PKT03 / VAX-C12_0 | always_accessible | low |  |
| 82 | EGFR_160-185_PKT04 / 173940 | always_accessible | low |  |
| 83 | EGFR_160-185_PKT04 / 97806 | always_accessible | low |  |
| 84 | EGFR_160-185_PKT04 / VAX-C12_0 | always_accessible | low |  |
| 85 | EGFR_160-185_PKT05 / 173940 | always_accessible | low |  |
| 86 | EGFR_160-185_PKT05 / 97806 | always_accessible | low |  |
| 87 | EGFR_160-185_PKT05 / VAX-C12_0 | always_accessible | low |  |
| 88 | EGFR_160-185_PKT11 / 173940 | always_accessible | low |  |
| 89 | EGFR_160-185_PKT11 / 97806 | always_accessible | low |  |
| 90 | EGFR_160-185_PKT11 / VAX-C12_0 | always_accessible | low |  |
| 91 | EGFR_160-185_PKT13 / 173940 | always_accessible | low |  |
| 92 | EGFR_160-185_PKT13 / 97806 | always_accessible | low |  |
| 93 | EGFR_160-185_PKT13 / VAX-C12_0 | always_accessible | low |  |
| 94 | EGFR_160-185_PKT19 / 173940 | always_accessible | low |  |
| 95 | EGFR_160-185_PKT19 / 97806 | always_accessible | low |  |
| 96 | EGFR_160-185_PKT19 / VAX-C12_0 | always_accessible | low |  |
| 97 | EGFR_160-185_PKT30 / 173940 | always_accessible | low |  |
| 98 | EGFR_160-185_PKT30 / 97806 | always_accessible | low |  |
| 99 | EGFR_160-185_PKT30 / VAX-C12_0 | always_accessible | low |  |
| 100 | EGFR_160-185_PKT37 / 173940 | always_accessible | low |  |
| 101 | EGFR_160-185_PKT37 / 97806 | always_accessible | low |  |
| 102 | EGFR_160-185_PKT37 / VAX-C12_0 | always_accessible | low |  |
| 103 | EGFR_160-185_PKT38 / 173940 | always_accessible | low |  |
| 104 | EGFR_160-185_PKT38 / 97806 | always_accessible | low |  |
| 105 | EGFR_160-185_PKT38 / VAX-C12_0 | always_accessible | low |  |
| 106 | EGFR_170-200_PKT01 / 173940 | always_accessible | low |  |
| 107 | EGFR_170-200_PKT01 / 97806 | always_accessible | low |  |
| 108 | EGFR_170-200_PKT01 / VAX-C12_0 | always_accessible | low |  |
| 109 | EGFR_170-200_PKT04 / 173940 | always_accessible | low |  |
| 110 | EGFR_170-200_PKT04 / 97806 | always_accessible | low |  |
| 111 | EGFR_170-200_PKT04 / VAX-C12_0 | always_accessible | low |  |
| 112 | EGFR_170-200_PKT05 / 173940 | always_accessible | low |  |
| 113 | EGFR_170-200_PKT05 / 97806 | always_accessible | low |  |
| 114 | EGFR_170-200_PKT05 / VAX-C12_0 | always_accessible | low |  |
| 115 | EGFR_170-200_PKT15 / 173940 | always_accessible | low |  |
| 116 | EGFR_170-200_PKT15 / 97806 | always_accessible | low |  |
| 117 | EGFR_170-200_PKT15 / VAX-C12_0 | always_accessible | low |  |
| 118 | EGFR_170-200_PKT25 / 173940 | always_accessible | low |  |
| 119 | EGFR_170-200_PKT25 / 97806 | always_accessible | low |  |
| 120 | EGFR_170-200_PKT25 / VAX-C12_0 | always_accessible | low |  |
| 121 | EGFR_170-200_PKT31 / 173940 | always_accessible | low |  |
| 122 | EGFR_170-200_PKT31 / 97806 | always_accessible | low |  |
| 123 | EGFR_170-200_PKT31 / VAX-C12_0 | always_accessible | low |  |
| 124 | EGFR_170-200_PKT36 / 173940 | always_accessible | low |  |
| 125 | EGFR_170-200_PKT36 / 97806 | always_accessible | low |  |
| 126 | EGFR_170-200_PKT36 / VAX-C12_0 | always_accessible | low |  |
| 127 | 3GT8_raw_PKT22 / 173940 | mostly_accessible | low |  |
| 128 | 3GT8_raw_PKT22 / 97806 | mostly_accessible | low |  |
| 129 | 3GT8_raw_PKT22 / VAX-C12_0 | mostly_accessible | low |  |
| 130 | EGFR_160-185_PKT14 / 173940 | mostly_accessible | low |  |
| 131 | EGFR_160-185_PKT14 / 97806 | mostly_accessible | low |  |
| 132 | EGFR_160-185_PKT14 / VAX-C12_0 | mostly_accessible | low |  |
| 133 | EGFR_160-185_PKT27 / 173940 | mostly_accessible | low |  |
| 134 | EGFR_160-185_PKT27 / 97806 | mostly_accessible | low |  |
| 135 | EGFR_160-185_PKT27 / VAX-C12_0 | mostly_accessible | low |  |
| 136 | EGFR_160-185_PKT29 / 173940 | mostly_accessible | low |  |
| 137 | EGFR_160-185_PKT29 / 97806 | mostly_accessible | low |  |
| 138 | EGFR_160-185_PKT29 / VAX-C12_0 | mostly_accessible | low |  |
| 139 | EGFR_160-185_PKT39 / 173940 | mostly_accessible | low |  |
| 140 | EGFR_160-185_PKT39 / 97806 | mostly_accessible | low |  |
| 141 | EGFR_160-185_PKT39 / VAX-C12_0 | mostly_accessible | low |  |
| 142 | EGFR_160-185_PKT40 / 173940 | mostly_accessible | low |  |
| 143 | EGFR_160-185_PKT40 / 97806 | mostly_accessible | low |  |
| 144 | EGFR_160-185_PKT40 / VAX-C12_0 | mostly_accessible | low |  |
| 145 | EGFR_160-185_PKT43 / 173940 | mostly_accessible | low |  |
| 146 | EGFR_160-185_PKT43 / 97806 | mostly_accessible | low |  |
| 147 | EGFR_160-185_PKT43 / VAX-C12_0 | mostly_accessible | low |  |
| 148 | EGFR_170-200_PKT13 / 173940 | mostly_accessible | low |  |
| 149 | EGFR_170-200_PKT13 / 97806 | mostly_accessible | low |  |
| 150 | EGFR_170-200_PKT13 / VAX-C12_0 | mostly_accessible | low |  |
| 151 | EGFR_170-200_PKT26 / 173940 | mostly_accessible | low |  |
| 152 | EGFR_170-200_PKT26 / 97806 | mostly_accessible | low |  |
| 153 | EGFR_170-200_PKT26 / VAX-C12_0 | mostly_accessible | low |  |
| 154 | EGFR_170-200_PKT33 / 173940 | mostly_accessible | low |  |
| 155 | EGFR_170-200_PKT33 / 97806 | mostly_accessible | low |  |
| 156 | EGFR_170-200_PKT33 / VAX-C12_0 | mostly_accessible | low |  |
| 157 | EGFR_170-200_PKT03 / (none) | always_accessible | high |  |
| 158 | EGFR_160-185_PKT06 / (none) | always_accessible | high |  |
| 159 | EGFR_170-200_PKT07 / (none) | always_accessible | high |  |
| 160 | EGFR_160-185_PKT17 / (none) | always_accessible | high |  |
| 161 | EGFR_170-200_PKT18 / (none) | always_accessible | high |  |
| 162 | 3GT8_raw_PKT13 / (none) | always_accessible | high |  |
| 163 | EGFR_170-200_PKT27 / (none) | always_accessible | high |  |
| 164 | EGFR_170-200_PKT29 / (none) | always_accessible | high |  |
| 165 | 3GT8_raw_PKT18 / (none) | always_accessible | high |  |
| 166 | EGFR_160-185_PKT08 / (none) | mostly_accessible | high |  |
| 167 | EGFR_170-200_PKT02 / (none) | single_state_only | high_provisional | Only 1 receptor state has pocket data. Cannot determine if s... |
| 168 | EGFR_160-185_PKT10 / (none) | mostly_accessible | high |  |
| 169 | 3GT8_raw_PKT04 / (none) | mostly_accessible | high |  |
| 170 | EGFR_170-200_PKT08 / (none) | mostly_accessible | high |  |
| 171 | EGFR_170-200_PKT09 / (none) | mostly_accessible | high |  |
| 172 | EGFR_170-200_PKT11 / (none) | mostly_accessible | high |  |
| 173 | EGFR_170-200_PKT14 / (none) | mostly_accessible | high |  |
| 174 | EGFR_160-185_PKT18 / (none) | mostly_accessible | high |  |
| 175 | EGFR_160-185_PKT20 / (none) | mostly_accessible | high |  |
| 176 | EGFR_160-185_PKT21 / (none) | mostly_accessible | high |  |
| 177 | EGFR_170-200_PKT16 / (none) | mostly_accessible | high |  |
| 178 | EGFR_160-185_PKT22 / (none) | mostly_accessible | high |  |
| 179 | EGFR_160-185_PKT24 / (none) | mostly_accessible | high |  |
| 180 | EGFR_160-185_PKT25 / (none) | mostly_accessible | high |  |
| 181 | EGFR_170-200_PKT19 / (none) | mostly_accessible | high |  |
| 182 | EGFR_160-185_PKT26 / (none) | mostly_accessible | high |  |
| 183 | EGFR_170-200_PKT20 / (none) | mostly_accessible | high |  |
| 184 | EGFR_170-200_PKT21 / (none) | mostly_accessible | high |  |
| 185 | 3GT8_raw_PKT09 / (none) | mostly_accessible | high |  |
| 186 | EGFR_170-200_PKT22 / (none) | mostly_accessible | high |  |
| 187 | EGFR_170-200_PKT23 / (none) | mostly_accessible | high |  |
| 188 | EGFR_170-200_PKT24 / (none) | mostly_accessible | high |  |
| 189 | EGFR_160-185_PKT28 / (none) | mostly_accessible | high |  |
| 190 | 3GT8_raw_PKT12 / (none) | mostly_accessible | high |  |
| 191 | EGFR_160-185_PKT31 / (none) | mostly_accessible | high |  |
| 192 | EGFR_160-185_PKT32 / (none) | mostly_accessible | high |  |
| 193 | 3GT8_raw_PKT06 / (none) | single_state_only | high_provisional | Only 1 receptor state has pocket data. Cannot determine if s... |
| 194 | EGFR_160-185_PKT33 / (none) | mostly_accessible | high |  |
| 195 | 3GT8_raw_PKT17 / (none) | mostly_accessible | high |  |
| 196 | EGFR_160-185_PKT35 / (none) | mostly_accessible | high |  |
| 197 | EGFR_160-185_PKT36 / (none) | mostly_accessible | high |  |
| 198 | EGFR_170-200_PKT30 / (none) | mostly_accessible | high |  |
| 199 | 3GT8_raw_PKT19 / (none) | mostly_accessible | high |  |
| 200 | 3GT8_raw_PKT15 / (none) | single_state_only | high_provisional | Only 1 receptor state has pocket data. Cannot determine if s... |
| 201 | 3GT8_raw_PKT21 / (none) | mostly_accessible | high |  |
| 202 | EGFR_170-200_PKT32 / (none) | mostly_accessible | high |  |
| 203 | EGFR_160-185_PKT42 / (none) | mostly_accessible | high |  |
| 204 | EGFR_160-185_PKT44 / (none) | mostly_accessible | high |  |
| 205 | EGFR_160-185_PKT41 / (none) | single_state_only | high_provisional | Only 1 receptor state has pocket data. Cannot determine if s... |
| 206 | EGFR_170-200_PKT35 / (none) | mostly_accessible | high |  |
| 207 | EGFR_170-200_PKT37 / (none) | mostly_accessible | high |  |

## 4. Confidence Adjustments

- 3GT8_raw_PKT20: medium → medium_provisional (flag)
- 3GT8_raw_PKT20: medium → medium_provisional (flag)
- 3GT8_raw_PKT20: medium → medium_provisional (flag)
- EGFR_170-200_PKT02: high → high_provisional (flag)
- 3GT8_raw_PKT06: high → high_provisional (flag)
- 3GT8_raw_PKT15: high → high_provisional (flag)
- EGFR_160-185_PKT41: high → high_provisional (flag)

## 5. Active Caveats

- Only 1 receptor state has pocket data. Cannot determine if state-specific or under-sampled. NOTE: Single-state data. Cross-state persistence unknown. Re-evaluate when additional receptor states have pocket data.

---

Generated by `egfr_pipeline.phase4.state_interpretation`