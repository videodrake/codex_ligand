# Phase 2 Patch Relationship Classification Note

## Overview

- Patch hotspot residues: 17
- Candidate pockets classified: 103

## Classification Thresholds

| Class | Criteria |
|-------|---------|
| orthosteric | hotspot overlap >= 2 residues AND fraction >= 0.25 |
| rim | hotspot overlap >= 1 residue(s) OR centroid <= 12.0A + overlap |
| allosteric | centroid <= 20.0A, no overlap |
| low_relevance | centroid > 20.0A, no overlap |

## Classification Summary

| Class | Count | Pockets |
|-------|-------|---------|
| orthosteric_candidate | 0 |  |
| rim_candidate | 21 | 3GT8_raw_PKT01, 3GT8_raw_PKT02, 3GT8_raw_PKT05, 3GT8_raw_PKT07, 3GT8_raw_PKT08, 3GT8_raw_PKT11, 3GT8_raw_PKT14, 3GT8_raw_PKT20, EGFR_160-185_PKT02, EGFR_160-185_PKT07, EGFR_160-185_PKT09, EGFR_160-185_PKT12, EGFR_160-185_PKT15, EGFR_160-185_PKT16, EGFR_160-185_PKT23, EGFR_160-185_PKT34, EGFR_170-200_PKT06, EGFR_170-200_PKT10, EGFR_170-200_PKT12, EGFR_170-200_PKT17, EGFR_170-200_PKT28 |
| allosteric_candidate | 31 | 3GT8_raw_PKT03, 3GT8_raw_PKT10, 3GT8_raw_PKT16, 3GT8_raw_PKT22, EGFR_160-185_PKT01, EGFR_160-185_PKT03, EGFR_160-185_PKT04, EGFR_160-185_PKT05, EGFR_160-185_PKT11, EGFR_160-185_PKT13, EGFR_160-185_PKT14, EGFR_160-185_PKT19, EGFR_160-185_PKT27, EGFR_160-185_PKT29, EGFR_160-185_PKT30, EGFR_160-185_PKT37, EGFR_160-185_PKT38, EGFR_160-185_PKT39, EGFR_160-185_PKT40, EGFR_160-185_PKT43, EGFR_170-200_PKT01, EGFR_170-200_PKT04, EGFR_170-200_PKT05, EGFR_170-200_PKT13, EGFR_170-200_PKT15, EGFR_170-200_PKT25, EGFR_170-200_PKT26, EGFR_170-200_PKT31, EGFR_170-200_PKT33, EGFR_170-200_PKT34, EGFR_170-200_PKT36 |
| low_relevance_candidate | 51 | 3GT8_raw_PKT04, 3GT8_raw_PKT06, 3GT8_raw_PKT09, 3GT8_raw_PKT12, 3GT8_raw_PKT13, 3GT8_raw_PKT15, 3GT8_raw_PKT17, 3GT8_raw_PKT18, 3GT8_raw_PKT19, 3GT8_raw_PKT21, EGFR_160-185_PKT06, EGFR_160-185_PKT08, EGFR_160-185_PKT10, EGFR_160-185_PKT17, EGFR_160-185_PKT18, EGFR_160-185_PKT20, EGFR_160-185_PKT21, EGFR_160-185_PKT22, EGFR_160-185_PKT24, EGFR_160-185_PKT25, EGFR_160-185_PKT26, EGFR_160-185_PKT28, EGFR_160-185_PKT31, EGFR_160-185_PKT32, EGFR_160-185_PKT33, EGFR_160-185_PKT35, EGFR_160-185_PKT36, EGFR_160-185_PKT41, EGFR_160-185_PKT42, EGFR_160-185_PKT44, EGFR_170-200_PKT02, EGFR_170-200_PKT03, EGFR_170-200_PKT07, EGFR_170-200_PKT08, EGFR_170-200_PKT09, EGFR_170-200_PKT11, EGFR_170-200_PKT14, EGFR_170-200_PKT16, EGFR_170-200_PKT18, EGFR_170-200_PKT19, EGFR_170-200_PKT20, EGFR_170-200_PKT21, EGFR_170-200_PKT22, EGFR_170-200_PKT23, EGFR_170-200_PKT24, EGFR_170-200_PKT27, EGFR_170-200_PKT29, EGFR_170-200_PKT30, EGFR_170-200_PKT32, EGFR_170-200_PKT35, EGFR_170-200_PKT37 |

## Per-Pocket Detail

| Pocket | Class | Overlap | Fraction | Dist (A) | Basis |
|--------|-------|---------|----------|----------|-------|
| 3GT8_raw_PKT01 | rim_candidate | 1 | 0.067 | 18.32 | partial_hotspot_overlap=1 (frac=0.07) |
| 3GT8_raw_PKT02 | rim_candidate | 1 | 0.067 | 22.30 | partial_hotspot_overlap=1 (frac=0.07) |
| 3GT8_raw_PKT03 | allosteric_candidate | 0 | 0.000 | 6.67 | centroid_nearby=6.7A <= 20.0A, no_hotspot_overlap |
| 3GT8_raw_PKT04 | low_relevance_candidate | 0 | 0.000 | 45.81 | centroid_distant=45.8A > 20.0A, no_overlap |
| 3GT8_raw_PKT05 | rim_candidate | 2 | 0.133 | 18.51 | partial_hotspot_overlap=2 (frac=0.13) |
| 3GT8_raw_PKT06 | low_relevance_candidate | 0 | 0.000 | 58.64 | centroid_distant=58.6A > 20.0A, no_overlap |
| 3GT8_raw_PKT07 | rim_candidate | 1 | 0.067 | 18.73 | partial_hotspot_overlap=1 (frac=0.07) |
| 3GT8_raw_PKT08 | rim_candidate | 2 | 0.133 | 18.74 | partial_hotspot_overlap=2 (frac=0.13) |
| 3GT8_raw_PKT09 | low_relevance_candidate | 0 | 0.000 | 29.96 | centroid_distant=30.0A > 20.0A, no_overlap |
| 3GT8_raw_PKT10 | allosteric_candidate | 0 | 0.000 | 16.80 | centroid_nearby=16.8A <= 20.0A, no_hotspot_overlap |
| 3GT8_raw_PKT11 | rim_candidate | 2 | 0.133 | 18.38 | partial_hotspot_overlap=2 (frac=0.13) |
| 3GT8_raw_PKT12 | low_relevance_candidate | 0 | 0.000 | 24.67 | centroid_distant=24.7A > 20.0A, no_overlap |
| 3GT8_raw_PKT13 | low_relevance_candidate | 0 | 0.000 | 33.00 | centroid_distant=33.0A > 20.0A, no_overlap |
| 3GT8_raw_PKT14 | rim_candidate | 1 | 0.067 | 13.73 | partial_hotspot_overlap=1 (frac=0.07) |
| 3GT8_raw_PKT15 | low_relevance_candidate | 0 | 0.000 | 42.61 | centroid_distant=42.6A > 20.0A, no_overlap |
| 3GT8_raw_PKT16 | allosteric_candidate | 0 | 0.000 | 11.91 | centroid_nearby=11.9A <= 20.0A, no_hotspot_overlap |
| 3GT8_raw_PKT17 | low_relevance_candidate | 0 | 0.000 | 49.36 | centroid_distant=49.4A > 20.0A, no_overlap |
| 3GT8_raw_PKT18 | low_relevance_candidate | 0 | 0.000 | 47.27 | centroid_distant=47.3A > 20.0A, no_overlap |
| 3GT8_raw_PKT19 | low_relevance_candidate | 0 | 0.000 | 41.92 | centroid_distant=41.9A > 20.0A, no_overlap |
| 3GT8_raw_PKT20 | rim_candidate | 2 | 0.133 | 50.98 | partial_hotspot_overlap=2 (frac=0.13) |
| 3GT8_raw_PKT21 | low_relevance_candidate | 0 | 0.000 | 57.93 | centroid_distant=57.9A > 20.0A, no_overlap |
| 3GT8_raw_PKT22 | allosteric_candidate | 0 | 0.000 | 18.12 | centroid_nearby=18.1A <= 20.0A, no_hotspot_overlap |
| EGFR_160-185_PKT01 | allosteric_candidate | 0 | 0.000 | 5.67 | centroid_nearby=5.7A <= 20.0A, no_hotspot_overlap |
| EGFR_160-185_PKT02 | rim_candidate | 2 | 0.143 | 12.96 | partial_hotspot_overlap=2 (frac=0.14) |
| EGFR_160-185_PKT03 | allosteric_candidate | 0 | 0.000 | 13.20 | centroid_nearby=13.2A <= 20.0A, no_hotspot_overlap |
| EGFR_160-185_PKT04 | allosteric_candidate | 0 | 0.000 | 12.09 | centroid_nearby=12.1A <= 20.0A, no_hotspot_overlap |
| EGFR_160-185_PKT05 | allosteric_candidate | 0 | 0.000 | 15.90 | centroid_nearby=15.9A <= 20.0A, no_hotspot_overlap |
| EGFR_160-185_PKT06 | low_relevance_candidate | 0 | 0.000 | 40.55 | centroid_distant=40.5A > 20.0A, no_overlap |
| EGFR_160-185_PKT07 | rim_candidate | 1 | 0.071 | 27.50 | partial_hotspot_overlap=1 (frac=0.07) |
| EGFR_160-185_PKT08 | low_relevance_candidate | 0 | 0.000 | 27.13 | centroid_distant=27.1A > 20.0A, no_overlap |
| EGFR_160-185_PKT09 | rim_candidate | 1 | 0.071 | 23.92 | partial_hotspot_overlap=1 (frac=0.07) |
| EGFR_160-185_PKT10 | low_relevance_candidate | 0 | 0.000 | 32.01 | centroid_distant=32.0A > 20.0A, no_overlap |
| EGFR_160-185_PKT11 | allosteric_candidate | 0 | 0.000 | 14.76 | centroid_nearby=14.8A <= 20.0A, no_hotspot_overlap |
| EGFR_160-185_PKT12 | rim_candidate | 1 | 0.071 | 6.34 | partial_hotspot_overlap=1 (frac=0.07) |
| EGFR_160-185_PKT13 | allosteric_candidate | 0 | 0.000 | 1.78 | centroid_nearby=1.8A <= 20.0A, no_hotspot_overlap |
| EGFR_160-185_PKT14 | allosteric_candidate | 0 | 0.000 | 13.87 | centroid_nearby=13.9A <= 20.0A, no_hotspot_overlap |
| EGFR_160-185_PKT15 | rim_candidate | 1 | 0.071 | 8.11 | partial_hotspot_overlap=1 (frac=0.07) |
| EGFR_160-185_PKT16 | rim_candidate | 2 | 0.143 | 10.30 | partial_hotspot_overlap=2 (frac=0.14) |
| EGFR_160-185_PKT17 | low_relevance_candidate | 0 | 0.000 | 29.86 | centroid_distant=29.9A > 20.0A, no_overlap |
| EGFR_160-185_PKT18 | low_relevance_candidate | 0 | 0.000 | 64.60 | centroid_distant=64.6A > 20.0A, no_overlap |
| EGFR_160-185_PKT19 | allosteric_candidate | 0 | 0.000 | 4.09 | centroid_nearby=4.1A <= 20.0A, no_hotspot_overlap |
| EGFR_160-185_PKT20 | low_relevance_candidate | 0 | 0.000 | 29.26 | centroid_distant=29.3A > 20.0A, no_overlap |
| EGFR_160-185_PKT21 | low_relevance_candidate | 0 | 0.000 | 21.60 | centroid_distant=21.6A > 20.0A, no_overlap |
| EGFR_160-185_PKT22 | low_relevance_candidate | 0 | 0.000 | 37.29 | centroid_distant=37.3A > 20.0A, no_overlap |
| EGFR_160-185_PKT23 | rim_candidate | 1 | 0.071 | 21.74 | partial_hotspot_overlap=1 (frac=0.07) |
| EGFR_160-185_PKT24 | low_relevance_candidate | 0 | 0.000 | 33.36 | centroid_distant=33.4A > 20.0A, no_overlap |
| EGFR_160-185_PKT25 | low_relevance_candidate | 0 | 0.000 | 26.00 | centroid_distant=26.0A > 20.0A, no_overlap |
| EGFR_160-185_PKT26 | low_relevance_candidate | 0 | 0.000 | 65.88 | centroid_distant=65.9A > 20.0A, no_overlap |
| EGFR_160-185_PKT27 | allosteric_candidate | 0 | 0.000 | 11.07 | centroid_nearby=11.1A <= 20.0A, no_hotspot_overlap |
| EGFR_160-185_PKT28 | low_relevance_candidate | 0 | 0.000 | 32.48 | centroid_distant=32.5A > 20.0A, no_overlap |
| EGFR_160-185_PKT29 | allosteric_candidate | 0 | 0.000 | 18.47 | centroid_nearby=18.5A <= 20.0A, no_hotspot_overlap |
| EGFR_160-185_PKT30 | allosteric_candidate | 0 | 0.000 | 1.76 | centroid_nearby=1.8A <= 20.0A, no_hotspot_overlap |
| EGFR_160-185_PKT31 | low_relevance_candidate | 0 | 0.000 | 39.21 | centroid_distant=39.2A > 20.0A, no_overlap |
| EGFR_160-185_PKT32 | low_relevance_candidate | 0 | 0.000 | 35.79 | centroid_distant=35.8A > 20.0A, no_overlap |
| EGFR_160-185_PKT33 | low_relevance_candidate | 0 | 0.000 | 54.28 | centroid_distant=54.3A > 20.0A, no_overlap |
| EGFR_160-185_PKT34 | rim_candidate | 1 | 0.071 | 22.53 | partial_hotspot_overlap=1 (frac=0.07) |
| EGFR_160-185_PKT35 | low_relevance_candidate | 0 | 0.000 | 36.29 | centroid_distant=36.3A > 20.0A, no_overlap |
| EGFR_160-185_PKT36 | low_relevance_candidate | 0 | 0.000 | 40.64 | centroid_distant=40.6A > 20.0A, no_overlap |
| EGFR_160-185_PKT37 | allosteric_candidate | 0 | 0.000 | 1.71 | centroid_nearby=1.7A <= 20.0A, no_hotspot_overlap |
| EGFR_160-185_PKT38 | allosteric_candidate | 0 | 0.000 | 4.72 | centroid_nearby=4.7A <= 20.0A, no_hotspot_overlap |
| EGFR_160-185_PKT39 | allosteric_candidate | 0 | 0.000 | 14.01 | centroid_nearby=14.0A <= 20.0A, no_hotspot_overlap |
| EGFR_160-185_PKT40 | allosteric_candidate | 0 | 0.000 | 19.30 | centroid_nearby=19.3A <= 20.0A, no_hotspot_overlap |
| EGFR_160-185_PKT41 | low_relevance_candidate | 0 | 0.000 | 20.55 | centroid_distant=20.5A > 20.0A, no_overlap |
| EGFR_160-185_PKT42 | low_relevance_candidate | 0 | 0.000 | 76.85 | centroid_distant=76.8A > 20.0A, no_overlap |
| EGFR_160-185_PKT43 | allosteric_candidate | 0 | 0.000 | 17.17 | centroid_nearby=17.2A <= 20.0A, no_hotspot_overlap |
| EGFR_160-185_PKT44 | low_relevance_candidate | 0 | 0.000 | 26.29 | centroid_distant=26.3A > 20.0A, no_overlap |
| EGFR_170-200_PKT01 | allosteric_candidate | 0 | 0.000 | 10.32 | centroid_nearby=10.3A <= 20.0A, no_hotspot_overlap |
| EGFR_170-200_PKT02 | low_relevance_candidate | 0 | 0.000 | 35.83 | centroid_distant=35.8A > 20.0A, no_overlap |
| EGFR_170-200_PKT03 | low_relevance_candidate | 0 | 0.000 | 22.14 | centroid_distant=22.1A > 20.0A, no_overlap |
| EGFR_170-200_PKT04 | allosteric_candidate | 0 | 0.000 | 18.44 | centroid_nearby=18.4A <= 20.0A, no_hotspot_overlap |
| EGFR_170-200_PKT05 | allosteric_candidate | 0 | 0.000 | 2.39 | centroid_nearby=2.4A <= 20.0A, no_hotspot_overlap |
| EGFR_170-200_PKT06 | rim_candidate | 2 | 0.143 | 16.80 | partial_hotspot_overlap=2 (frac=0.14) |
| EGFR_170-200_PKT07 | low_relevance_candidate | 0 | 0.000 | 34.53 | centroid_distant=34.5A > 20.0A, no_overlap |
| EGFR_170-200_PKT08 | low_relevance_candidate | 0 | 0.000 | 35.78 | centroid_distant=35.8A > 20.0A, no_overlap |
| EGFR_170-200_PKT09 | low_relevance_candidate | 0 | 0.000 | 24.04 | centroid_distant=24.0A > 20.0A, no_overlap |
| EGFR_170-200_PKT10 | rim_candidate | 1 | 0.071 | 8.18 | partial_hotspot_overlap=1 (frac=0.07) |
| EGFR_170-200_PKT11 | low_relevance_candidate | 0 | 0.000 | 32.45 | centroid_distant=32.5A > 20.0A, no_overlap |
| EGFR_170-200_PKT12 | rim_candidate | 1 | 0.071 | 14.35 | partial_hotspot_overlap=1 (frac=0.07) |
| EGFR_170-200_PKT13 | allosteric_candidate | 0 | 0.000 | 7.34 | centroid_nearby=7.3A <= 20.0A, no_hotspot_overlap |
| EGFR_170-200_PKT14 | low_relevance_candidate | 0 | 0.000 | 34.07 | centroid_distant=34.1A > 20.0A, no_overlap |
| EGFR_170-200_PKT15 | allosteric_candidate | 0 | 0.000 | 11.35 | centroid_nearby=11.3A <= 20.0A, no_hotspot_overlap |
| EGFR_170-200_PKT16 | low_relevance_candidate | 0 | 0.000 | 36.64 | centroid_distant=36.6A > 20.0A, no_overlap |
| EGFR_170-200_PKT17 | rim_candidate | 3 | 0.214 | 16.75 | partial_hotspot_overlap=3 (frac=0.21) |
| EGFR_170-200_PKT18 | low_relevance_candidate | 0 | 0.000 | 24.26 | centroid_distant=24.3A > 20.0A, no_overlap |
| EGFR_170-200_PKT19 | low_relevance_candidate | 0 | 0.000 | 40.90 | centroid_distant=40.9A > 20.0A, no_overlap |
| EGFR_170-200_PKT20 | low_relevance_candidate | 0 | 0.000 | 29.48 | centroid_distant=29.5A > 20.0A, no_overlap |
| EGFR_170-200_PKT21 | low_relevance_candidate | 0 | 0.000 | 28.02 | centroid_distant=28.0A > 20.0A, no_overlap |
| EGFR_170-200_PKT22 | low_relevance_candidate | 0 | 0.000 | 25.37 | centroid_distant=25.4A > 20.0A, no_overlap |
| EGFR_170-200_PKT23 | low_relevance_candidate | 0 | 0.000 | 54.14 | centroid_distant=54.1A > 20.0A, no_overlap |
| EGFR_170-200_PKT24 | low_relevance_candidate | 0 | 0.000 | 26.65 | centroid_distant=26.7A > 20.0A, no_overlap |
| EGFR_170-200_PKT25 | allosteric_candidate | 0 | 0.000 | 2.43 | centroid_nearby=2.4A <= 20.0A, no_hotspot_overlap |
| EGFR_170-200_PKT26 | allosteric_candidate | 0 | 0.000 | 12.56 | centroid_nearby=12.6A <= 20.0A, no_hotspot_overlap |
| EGFR_170-200_PKT27 | low_relevance_candidate | 0 | 0.000 | 30.54 | centroid_distant=30.5A > 20.0A, no_overlap |
| EGFR_170-200_PKT28 | rim_candidate | 2 | 0.143 | 34.31 | partial_hotspot_overlap=2 (frac=0.14) |
| EGFR_170-200_PKT29 | low_relevance_candidate | 0 | 0.000 | 28.13 | centroid_distant=28.1A > 20.0A, no_overlap |
| EGFR_170-200_PKT30 | low_relevance_candidate | 0 | 0.000 | 34.16 | centroid_distant=34.2A > 20.0A, no_overlap |
| EGFR_170-200_PKT31 | allosteric_candidate | 0 | 0.000 | 0.90 | centroid_nearby=0.9A <= 20.0A, no_hotspot_overlap |
| EGFR_170-200_PKT32 | low_relevance_candidate | 0 | 0.000 | 67.66 | centroid_distant=67.7A > 20.0A, no_overlap |
| EGFR_170-200_PKT33 | allosteric_candidate | 0 | 0.000 | 8.86 | centroid_nearby=8.9A <= 20.0A, no_hotspot_overlap |
| EGFR_170-200_PKT34 | allosteric_candidate | 0 | 0.000 | 9.37 | centroid_nearby=9.4A <= 20.0A, no_hotspot_overlap |
| EGFR_170-200_PKT35 | low_relevance_candidate | 0 | 0.000 | 34.45 | centroid_distant=34.4A > 20.0A, no_overlap |
| EGFR_170-200_PKT36 | allosteric_candidate | 0 | 0.000 | 7.05 | centroid_nearby=7.1A <= 20.0A, no_hotspot_overlap |
| EGFR_170-200_PKT37 | low_relevance_candidate | 0 | 0.000 | 31.07 | centroid_distant=31.1A > 20.0A, no_overlap |

---

Generated by `egfr_pipeline.phase2.patch_relationship`