# Phase 1 Orientation Filter Note

## Task Group 1.2A: Orientation-Aware Filtering

## 1. Problem Statement

The MYO1D beta-meander is a flat, thin beta-sheet ribbon (5 consecutive beta-strands, sheets 8-12). In global blind docking, it can land on the receptor in two orientations:

- Correct: active face of sheets 8/9 contacts the receptor.
- Flipped: back face contacts the receptor.

Contact-count thresholds alone can miss this distinction because backbone and edge-on contacts can satisfy distance cutoffs even when side-chain orientation is wrong.

## 2. Algorithm Summary

1. Collect active-face CA atoms from sheets 8 and 9.
2. Fit the sheet-plane normal via PCA.
3. Orient that normal with CA->CB probe residues (VAL962, VAL964, SER971).
4. Build a local receptor-direction vector from receptor CA atoms near the active face.
5. Compute `orientation_score = dot(active_face_normal, receptor_direction)`.

Classification:

- `score > 0.15`: pass
- `score < -0.15`: fail
- `|score| < 0.15`: ambiguous

## 3. Implementation

File: `egfr_pipeline/phase1/orientation_filter.py`

Core functions:

- `compute_orientation_score()`
- `compute_sheet_plane_normal()`
- `orient_normal_to_active_face()`
- `compute_receptor_interface_centroid()`
- `process_state_orientation()`
- `merge_orientation_into_models()`
- `validate_pilot_structures()`

Primary output file: `orientation_filter_log.csv`

## 4. Integration Points

Pipeline position:

```text
TG 1.1 PyRosetta docking
  -> TG 1.2 interface extraction
  -> TG 1.2A orientation filter
  -> TG 1.3 consensus (pass models only)
```

After orientation filtering, run with `--merge` so `pyrosetta_interface_models.csv` gains `orientation_score` and `orientation_class`.

## 5. Residue Definitions

| Sheet | Residues | Role |
|---|---|---|
| 8 | 961, 962, 963, 964 | Active face |
| 9 | 968, 969, 970, 971, 972 | Active face |
| 10 | ~977-980 | Neutral |
| 11 | ~985-988 | Neutral |
| 12 | ~993-997 | Structural support (role TBD) |

## 6. Commands

```bash
conda activate pyrosetta
python -m egfr_pipeline.phase1.orientation_filter --merge
python -m egfr_pipeline.phase1.orientation_filter --state 3GT8_raw --merge
```

Pilot validation:

```bash
python -m egfr_pipeline.phase1.orientation_filter \
  --pilot_dir /path/to/pilot/final_result/ \
  --pilot_out output/phase1_ppi/orientation_filter_pilot_validation.csv
```

## 7. Limits and Pending Validation

Known limits:

- Requires PyRosetta-accessible CA/CB coordinates.
- Assumes beta-sheet geometry remains stable.
- Sheet 12 functional role is not finalized.

Pending server checks:

- Run on test docking output.
- Validate pilot set classifications.
- Cross-tab class labels vs dG/dSASA/sc.
- Recalibrate ambiguous threshold if needed.
