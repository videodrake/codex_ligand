# Phase 2 → Phase 3 Handoff Note

## Overview

This document accompanies `phase3_candidate_pocket_reference.csv`,
the primary input for Phase 3 diversity-aware pocket-guided docking.

- Total candidate pockets exported: 103
  - 3GT8_raw: 22 pockets
  - EGFR_160-185: 44 pockets
  - EGFR_170-200: 37 pockets

## Docking Priority Distribution

| Priority | Count | Pockets | Phase 3 Action |
|----------|-------|---------|---------------|
| primary | 2 | 3GT8_raw_PKT07, EGFR_170-200_PKT34 | Full docking budget allocation |
| secondary | 21 | 3GT8_raw_PKT01, 3GT8_raw_PKT02, 3GT8_raw_PKT05, 3GT8_raw_PKT08, 3GT8_raw_PKT10, 3GT8_raw_PKT11, 3GT8_raw_PKT14, 3GT8_raw_PKT20, EGFR_160-185_PKT02, EGFR_160-185_PKT07, EGFR_160-185_PKT09, EGFR_160-185_PKT12, EGFR_160-185_PKT15, EGFR_160-185_PKT16, EGFR_160-185_PKT23, EGFR_160-185_PKT34, EGFR_170-200_PKT06, EGFR_170-200_PKT10, EGFR_170-200_PKT12, EGFR_170-200_PKT17, EGFR_170-200_PKT28 | Reduced budget, include in diversity round |
| exploratory | 29 | 3GT8_raw_PKT03, 3GT8_raw_PKT16, 3GT8_raw_PKT22, EGFR_160-185_PKT01, EGFR_160-185_PKT03, EGFR_160-185_PKT04, EGFR_160-185_PKT05, EGFR_160-185_PKT11, EGFR_160-185_PKT13, EGFR_160-185_PKT14, EGFR_160-185_PKT19, EGFR_160-185_PKT27, EGFR_160-185_PKT29, EGFR_160-185_PKT30, EGFR_160-185_PKT37, EGFR_160-185_PKT38, EGFR_160-185_PKT39, EGFR_160-185_PKT40, EGFR_160-185_PKT43, EGFR_170-200_PKT01, EGFR_170-200_PKT04, EGFR_170-200_PKT05, EGFR_170-200_PKT13, EGFR_170-200_PKT15, EGFR_170-200_PKT25, EGFR_170-200_PKT26, EGFR_170-200_PKT31, EGFR_170-200_PKT33, EGFR_170-200_PKT36 | Minimal budget, only if resources allow |
| skip | 51 | 3GT8_raw_PKT04, 3GT8_raw_PKT06, 3GT8_raw_PKT09, 3GT8_raw_PKT12, 3GT8_raw_PKT13, 3GT8_raw_PKT15, 3GT8_raw_PKT17, 3GT8_raw_PKT18, 3GT8_raw_PKT19, 3GT8_raw_PKT21, EGFR_160-185_PKT06, EGFR_160-185_PKT08, EGFR_160-185_PKT10, EGFR_160-185_PKT17, EGFR_160-185_PKT18, EGFR_160-185_PKT20, EGFR_160-185_PKT21, EGFR_160-185_PKT22, EGFR_160-185_PKT24, EGFR_160-185_PKT25, EGFR_160-185_PKT26, EGFR_160-185_PKT28, EGFR_160-185_PKT31, EGFR_160-185_PKT32, EGFR_160-185_PKT33, EGFR_160-185_PKT35, EGFR_160-185_PKT36, EGFR_160-185_PKT41, EGFR_160-185_PKT42, EGFR_160-185_PKT44, EGFR_170-200_PKT02, EGFR_170-200_PKT03, EGFR_170-200_PKT07, EGFR_170-200_PKT08, EGFR_170-200_PKT09, EGFR_170-200_PKT11, EGFR_170-200_PKT14, EGFR_170-200_PKT16, EGFR_170-200_PKT18, EGFR_170-200_PKT19, EGFR_170-200_PKT20, EGFR_170-200_PKT21, EGFR_170-200_PKT22, EGFR_170-200_PKT23, EGFR_170-200_PKT24, EGFR_170-200_PKT27, EGFR_170-200_PKT29, EGFR_170-200_PKT30, EGFR_170-200_PKT32, EGFR_170-200_PKT35, EGFR_170-200_PKT37 | Not recommended for docking |

## Docking Priority Definitions

| Priority | Criteria |
|----------|----------|
| primary | tier_1, or tier_2 + orthosteric/rim |
| secondary | tier_2 + allosteric, or tier_3 + orthosteric/rim |
| exploratory | tier_2/3 + uncertain class, or tier_3 + allosteric |
| skip | tier_3 + low_relevance (docking budget not justified) |

## PPI Patch Relationship Distribution

| Class | Count |
|-------|-------|
| orthosteric_candidate | 0 |
| rim_candidate | 21 |
| allosteric_candidate | 31 |
| low_relevance_candidate | 51 |

## Druggability Tier Distribution

| Tier | Count |
|------|-------|
| tier_1 | 2 |
| tier_2 | 1 |
| tier_3 | 100 |

## Per-Pocket Reference

| Pocket | State | Relationship | Tier | Priority | Box (Å) |
|--------|-------|-------------|------|----------|---------|
| 3GT8_raw_PKT01 | 3GT8_raw | rim_candidate | tier_3 | secondary | 27.5×25.9×26.7 |
| 3GT8_raw_PKT02 | 3GT8_raw | rim_candidate | tier_3 | secondary | 30.0×19.1×27.8 |
| 3GT8_raw_PKT03 | 3GT8_raw | allosteric_candidate | tier_3 | exploratory | 23.2×23.2×23.1 |
| 3GT8_raw_PKT04 | 3GT8_raw | low_relevance_candidate | tier_3 | skip | 19.2×17.7×18.7 |
| 3GT8_raw_PKT05 | 3GT8_raw | rim_candidate | tier_3 | secondary | 18.3×21.9×27.5 |
| 3GT8_raw_PKT06 | 3GT8_raw | low_relevance_candidate | tier_3 | skip | 18.1×19.2×18.9 |
| 3GT8_raw_PKT07 | 3GT8_raw | rim_candidate | tier_1 | primary | 29.6×30.5×46.2 |
| 3GT8_raw_PKT08 | 3GT8_raw | rim_candidate | tier_3 | secondary | 17.4×20.4×19.6 |
| 3GT8_raw_PKT09 | 3GT8_raw | low_relevance_candidate | tier_3 | skip | 17.7×21.0×21.8 |
| 3GT8_raw_PKT10 | 3GT8_raw | allosteric_candidate | tier_2 | secondary | 25.1×21.1×20.6 |
| 3GT8_raw_PKT11 | 3GT8_raw | rim_candidate | tier_3 | secondary | 24.1×24.1×26.2 |
| 3GT8_raw_PKT12 | 3GT8_raw | low_relevance_candidate | tier_3 | skip | 20.6×18.6×19.0 |
| 3GT8_raw_PKT13 | 3GT8_raw | low_relevance_candidate | tier_3 | skip | 20.6×20.7×17.6 |
| 3GT8_raw_PKT14 | 3GT8_raw | rim_candidate | tier_3 | secondary | 20.7×20.8×20.1 |
| 3GT8_raw_PKT15 | 3GT8_raw | low_relevance_candidate | tier_3 | skip | 16.9×14.1×17.8 |
| 3GT8_raw_PKT16 | 3GT8_raw | allosteric_candidate | tier_3 | exploratory | 21.2×19.7×15.7 |
| 3GT8_raw_PKT17 | 3GT8_raw | low_relevance_candidate | tier_3 | skip | 17.6×17.7×18.0 |
| 3GT8_raw_PKT18 | 3GT8_raw | low_relevance_candidate | tier_3 | skip | 19.7×16.2×20.4 |
| 3GT8_raw_PKT19 | 3GT8_raw | low_relevance_candidate | tier_3 | skip | 19.8×19.3×18.8 |
| 3GT8_raw_PKT20 | 3GT8_raw | rim_candidate | tier_3 | secondary | 19.8×22.0×18.6 |
| 3GT8_raw_PKT21 | 3GT8_raw | low_relevance_candidate | tier_3 | skip | 17.8×17.0×18.7 |
| 3GT8_raw_PKT22 | 3GT8_raw | allosteric_candidate | tier_3 | exploratory | 25.0×21.2×23.9 |
| EGFR_160-185_PKT01 | EGFR_160-185 | allosteric_candidate | tier_3 | exploratory | 21.2×27.4×26.2 |
| EGFR_160-185_PKT02 | EGFR_160-185 | rim_candidate | tier_3 | secondary | 32.3×33.8×30.9 |
| EGFR_160-185_PKT03 | EGFR_160-185 | allosteric_candidate | tier_3 | exploratory | 24.2×28.8×24.8 |
| EGFR_160-185_PKT04 | EGFR_160-185 | allosteric_candidate | tier_3 | exploratory | 24.7×24.8×26.8 |
| EGFR_160-185_PKT05 | EGFR_160-185 | allosteric_candidate | tier_3 | exploratory | 22.1×20.0×23.2 |
| EGFR_160-185_PKT06 | EGFR_160-185 | low_relevance_candidate | tier_3 | skip | 20.3×22.8×15.5 |
| EGFR_160-185_PKT07 | EGFR_160-185 | rim_candidate | tier_3 | secondary | 17.4×23.2×17.9 |
| EGFR_160-185_PKT08 | EGFR_160-185 | low_relevance_candidate | tier_3 | skip | 16.1×17.0×16.4 |
| EGFR_160-185_PKT09 | EGFR_160-185 | rim_candidate | tier_3 | secondary | 19.4×19.9×21.0 |
| EGFR_160-185_PKT10 | EGFR_160-185 | low_relevance_candidate | tier_3 | skip | 23.4×18.1×16.8 |
| EGFR_160-185_PKT11 | EGFR_160-185 | allosteric_candidate | tier_3 | exploratory | 19.1×18.2×17.8 |
| EGFR_160-185_PKT12 | EGFR_160-185 | rim_candidate | tier_3 | secondary | 19.3×22.7×17.8 |
| EGFR_160-185_PKT13 | EGFR_160-185 | allosteric_candidate | tier_3 | exploratory | 28.9×29.2×36.3 |
| EGFR_160-185_PKT14 | EGFR_160-185 | allosteric_candidate | tier_3 | exploratory | 19.7×16.7×21.1 |
| EGFR_160-185_PKT15 | EGFR_160-185 | rim_candidate | tier_3 | secondary | 18.7×24.4×20.5 |
| EGFR_160-185_PKT16 | EGFR_160-185 | rim_candidate | tier_3 | secondary | 20.8×19.0×20.2 |
| EGFR_160-185_PKT17 | EGFR_160-185 | low_relevance_candidate | tier_3 | skip | 20.9×18.8×21.2 |
| EGFR_160-185_PKT18 | EGFR_160-185 | low_relevance_candidate | tier_3 | skip | 18.4×17.8×24.3 |
| EGFR_160-185_PKT19 | EGFR_160-185 | allosteric_candidate | tier_3 | exploratory | 27.4×28.7×30.2 |
| EGFR_160-185_PKT20 | EGFR_160-185 | low_relevance_candidate | tier_3 | skip | 19.8×21.1×17.0 |
| EGFR_160-185_PKT21 | EGFR_160-185 | low_relevance_candidate | tier_3 | skip | 19.5×20.0×17.2 |
| EGFR_160-185_PKT22 | EGFR_160-185 | low_relevance_candidate | tier_3 | skip | 15.5×21.9×23.6 |
| EGFR_160-185_PKT23 | EGFR_160-185 | rim_candidate | tier_3 | secondary | 27.1×35.1×21.7 |
| EGFR_160-185_PKT24 | EGFR_160-185 | low_relevance_candidate | tier_3 | skip | 24.6×26.2×21.6 |
| EGFR_160-185_PKT25 | EGFR_160-185 | low_relevance_candidate | tier_3 | skip | 21.8×18.7×22.1 |
| EGFR_160-185_PKT26 | EGFR_160-185 | low_relevance_candidate | tier_3 | skip | 16.0×19.2×21.8 |
| EGFR_160-185_PKT27 | EGFR_160-185 | allosteric_candidate | tier_3 | exploratory | 30.0×30.7×26.7 |
| EGFR_160-185_PKT28 | EGFR_160-185 | low_relevance_candidate | tier_3 | skip | 23.0×18.6×19.3 |
| EGFR_160-185_PKT29 | EGFR_160-185 | allosteric_candidate | tier_3 | exploratory | 18.8×21.3×22.7 |
| EGFR_160-185_PKT30 | EGFR_160-185 | allosteric_candidate | tier_3 | exploratory | 22.5×23.0×21.7 |
| EGFR_160-185_PKT31 | EGFR_160-185 | low_relevance_candidate | tier_3 | skip | 20.9×21.9×20.4 |
| EGFR_160-185_PKT32 | EGFR_160-185 | low_relevance_candidate | tier_3 | skip | 17.7×17.0×20.7 |
| EGFR_160-185_PKT33 | EGFR_160-185 | low_relevance_candidate | tier_3 | skip | 17.9×19.5×18.6 |
| EGFR_160-185_PKT34 | EGFR_160-185 | rim_candidate | tier_3 | secondary | 15.6×20.9×19.1 |
| EGFR_160-185_PKT35 | EGFR_160-185 | low_relevance_candidate | tier_3 | skip | 18.9×20.1×20.0 |
| EGFR_160-185_PKT36 | EGFR_160-185 | low_relevance_candidate | tier_3 | skip | 17.5×19.1×17.7 |
| EGFR_160-185_PKT37 | EGFR_160-185 | allosteric_candidate | tier_3 | exploratory | 22.3×19.7×24.7 |
| EGFR_160-185_PKT38 | EGFR_160-185 | allosteric_candidate | tier_3 | exploratory | 18.7×21.0×17.0 |
| EGFR_160-185_PKT39 | EGFR_160-185 | allosteric_candidate | tier_3 | exploratory | 20.1×21.6×19.9 |
| EGFR_160-185_PKT40 | EGFR_160-185 | allosteric_candidate | tier_3 | exploratory | 23.3×21.4×21.1 |
| EGFR_160-185_PKT41 | EGFR_160-185 | low_relevance_candidate | tier_3 | skip | 19.5×20.2×16.9 |
| EGFR_160-185_PKT42 | EGFR_160-185 | low_relevance_candidate | tier_3 | skip | 17.4×19.1×15.5 |
| EGFR_160-185_PKT43 | EGFR_160-185 | allosteric_candidate | tier_3 | exploratory | 21.3×20.4×24.0 |
| EGFR_160-185_PKT44 | EGFR_160-185 | low_relevance_candidate | tier_3 | skip | 30.7×31.2×28.2 |
| EGFR_170-200_PKT01 | EGFR_170-200 | allosteric_candidate | tier_3 | exploratory | 22.7×23.3×24.6 |
| EGFR_170-200_PKT02 | EGFR_170-200 | low_relevance_candidate | tier_3 | skip | 18.3×18.2×16.7 |
| EGFR_170-200_PKT03 | EGFR_170-200 | low_relevance_candidate | tier_3 | skip | 19.3×23.0×20.9 |
| EGFR_170-200_PKT04 | EGFR_170-200 | allosteric_candidate | tier_3 | exploratory | 19.1×18.4×18.3 |
| EGFR_170-200_PKT05 | EGFR_170-200 | allosteric_candidate | tier_3 | exploratory | 24.2×18.8×24.1 |
| EGFR_170-200_PKT06 | EGFR_170-200 | rim_candidate | tier_3 | secondary | 24.3×29.5×27.4 |
| EGFR_170-200_PKT07 | EGFR_170-200 | low_relevance_candidate | tier_3 | skip | 17.5×26.3×22.1 |
| EGFR_170-200_PKT08 | EGFR_170-200 | low_relevance_candidate | tier_3 | skip | 23.3×24.0×19.3 |
| EGFR_170-200_PKT09 | EGFR_170-200 | low_relevance_candidate | tier_3 | skip | 33.1×31.4×29.8 |
| EGFR_170-200_PKT10 | EGFR_170-200 | rim_candidate | tier_3 | secondary | 24.3×25.5×21.0 |
| EGFR_170-200_PKT11 | EGFR_170-200 | low_relevance_candidate | tier_3 | skip | 21.1×21.6×16.0 |
| EGFR_170-200_PKT12 | EGFR_170-200 | rim_candidate | tier_3 | secondary | 22.1×21.0×27.0 |
| EGFR_170-200_PKT13 | EGFR_170-200 | allosteric_candidate | tier_3 | exploratory | 21.3×21.8×18.2 |
| EGFR_170-200_PKT14 | EGFR_170-200 | low_relevance_candidate | tier_3 | skip | 18.6×19.4×20.1 |
| EGFR_170-200_PKT15 | EGFR_170-200 | allosteric_candidate | tier_3 | exploratory | 23.0×24.4×21.0 |
| EGFR_170-200_PKT16 | EGFR_170-200 | low_relevance_candidate | tier_3 | skip | 21.7×21.2×15.6 |
| EGFR_170-200_PKT17 | EGFR_170-200 | rim_candidate | tier_3 | secondary | 22.1×20.2×21.2 |
| EGFR_170-200_PKT18 | EGFR_170-200 | low_relevance_candidate | tier_3 | skip | 18.9×19.1×18.8 |
| EGFR_170-200_PKT19 | EGFR_170-200 | low_relevance_candidate | tier_3 | skip | 18.4×19.0×18.1 |
| EGFR_170-200_PKT20 | EGFR_170-200 | low_relevance_candidate | tier_3 | skip | 18.6×19.4×14.0 |
| EGFR_170-200_PKT21 | EGFR_170-200 | low_relevance_candidate | tier_3 | skip | 19.3×20.6×18.7 |
| EGFR_170-200_PKT22 | EGFR_170-200 | low_relevance_candidate | tier_3 | skip | 16.8×19.7×17.9 |
| EGFR_170-200_PKT23 | EGFR_170-200 | low_relevance_candidate | tier_3 | skip | 16.8×21.9×20.9 |
| EGFR_170-200_PKT24 | EGFR_170-200 | low_relevance_candidate | tier_3 | skip | 19.0×18.5×16.0 |
| EGFR_170-200_PKT25 | EGFR_170-200 | allosteric_candidate | tier_3 | exploratory | 26.4×20.4×21.3 |
| EGFR_170-200_PKT26 | EGFR_170-200 | allosteric_candidate | tier_3 | exploratory | 20.4×21.4×17.3 |
| EGFR_170-200_PKT27 | EGFR_170-200 | low_relevance_candidate | tier_3 | skip | 19.0×21.0×19.8 |
| EGFR_170-200_PKT28 | EGFR_170-200 | rim_candidate | tier_3 | secondary | 20.2×21.0×16.5 |
| EGFR_170-200_PKT29 | EGFR_170-200 | low_relevance_candidate | tier_3 | skip | 22.2×24.1×20.0 |
| EGFR_170-200_PKT30 | EGFR_170-200 | low_relevance_candidate | tier_3 | skip | 18.1×16.9×21.3 |
| EGFR_170-200_PKT31 | EGFR_170-200 | allosteric_candidate | tier_3 | exploratory | 30.5×25.6×30.3 |
| EGFR_170-200_PKT32 | EGFR_170-200 | low_relevance_candidate | tier_3 | skip | 24.8×20.5×24.8 |
| EGFR_170-200_PKT33 | EGFR_170-200 | allosteric_candidate | tier_3 | exploratory | 21.5×25.5×23.7 |
| EGFR_170-200_PKT34 | EGFR_170-200 | allosteric_candidate | tier_1 | primary | 35.5×45.7×31.9 |
| EGFR_170-200_PKT35 | EGFR_170-200 | low_relevance_candidate | tier_3 | skip | 20.7×20.7×21.0 |
| EGFR_170-200_PKT36 | EGFR_170-200 | allosteric_candidate | tier_3 | exploratory | 22.8×24.1×25.5 |
| EGFR_170-200_PKT37 | EGFR_170-200 | low_relevance_candidate | tier_3 | skip | 22.9×22.0×20.2 |

## Export Validation

All compatibility checks passed.

## Phase 3 Consumption Instructions

1. Load `phase3_candidate_pocket_reference.csv`
2. Filter by `docking_priority` to allocate docking budget:
   - `primary`: full docking runs with diversity sampling
   - `secondary`: reduced decoy count or single-seed
   - `exploratory`: only if budget allows
   - `skip`: do not dock
3. Use `centroid_x/y/z` and `box_size_x/y/z` for Vina/AutoDock box setup
4. Use `relationship_class` for Phase 4 perturbation relevance context
5. Use `state_class` for cross-state robustness weighting

## Cross-State Note

All three receptor states have pocket data.

---

Generated by `egfr_pipeline.phase2.phase3_export`