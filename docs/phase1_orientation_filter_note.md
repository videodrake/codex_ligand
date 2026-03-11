# Phase 1 Orientation Filter Note

## Task Group 1.2A: Orientation-Aware Filtering

---

## 1. Problem Statement

The MYO1D beta-meander is a flat, thin β-sheet ribbon (5 consecutive β-strands, sheets 8–12). In global blind docking, it can land on the receptor in two orientations:

- **Correct:** Active face of sheets 8/9 packs against receptor (functional side chains make direct contacts)
- **Flipped:** Back face contacts receptor (side chains point into solvent)

Simple residue contact counts cannot distinguish these because backbone atoms or edge-on contacts span both faces (sheet thickness ~5–7 Å, contact cutoff 6–10 Å).

---

## 2. Algorithm: Dual-Vector Orientation Test

### Steps

1. **Collect active-face Cα** — 9 residues from sheets 8+9 (VAL961–VAL964, VAL968–LEU972)
2. **PCA plane normal** — Smallest-variance component of SVD on centered Cα coordinates
3. **Multi-probe Cα→Cβ consensus** — Majority vote from VAL962, VAL964, SER971 determines which face is "active"
4. **Local receptor centroid** — Receptor Cα within 10 Å of any active-face Cα (not global centroid)
5. **Dot product** — `active_face_normal · receptor_direction`:
   - `> 0.15` → **PASS** (active face toward receptor)
   - `< -0.15` → **FAIL** (active face away)
   - `|score| < 0.15` → **AMBIGUOUS** (edge-on, needs manual review)

### Why PCA (not cross product)

PCA uses all 9 Cα positions, is robust to sheet curvature, and avoids arbitrary pair selection.

### Why localized receptor centroid

Full kinase domain (~309 residues) geometric center is far from the actual contact area. Using only receptor residues within 10 Å of the active face gives a contact-relevant direction vector.

---

## 3. Implementation

### File

`egfr_pipeline/phase1/orientation_filter.py`

### Core functions

| Function | Purpose |
|----------|---------|
| `compute_orientation_score()` | Full orientation test for one pose → score + class |
| `compute_sheet_plane_normal()` | PCA-based plane normal from Cα coordinates |
| `orient_normal_to_active_face()` | Multi-probe Cα→Cβ consensus disambiguation |
| `compute_receptor_interface_centroid()` | Local receptor centroid near active face |
| `process_state_orientation()` | Batch-process all finals for one receptor state |
| `merge_orientation_into_models()` | Add orientation columns to interface models CSV |
| `validate_pilot_structures()` | Retroactive pilot validation |

### Output schema

`orientation_filter_log.csv`:

| Column | Type | Description |
|--------|------|-------------|
| model_id | str | PDB filename (e.g., Rank01_C01_M01_S-15.20.pdb) |
| receptor_id | str | Receptor state (e.g., 3GT8_raw) |
| seed_index | int | Seed index (0-4) |
| orientation_score | float | Dot product (-1 to +1) |
| orientation_class | str | pass / fail / ambiguous / error codes |
| n_active_face_ca | int | Active-face Cα found (max 9) |
| n_back_face_ca | int | Back-face Cα found |
| n_receptor_contact_ca | int | Receptor Cα within contact_dist |
| sheet_centroid | str | [x,y,z] of sheet 8/9 centroid |
| normal_vector | str | [x,y,z] of oriented active-face normal |
| receptor_direction | str | [x,y,z] toward receptor contact centroid |
| error | str | Error message (empty if successful) |
| source_pdb | str | Relative path to source PDB |

---

## 4. Integration Points

### Pipeline position

```
TG 1.1: PyRosetta docking
    ↓
TG 1.2: Interface residue extraction (all models)
    ↓
★ TG 1.2A: Orientation filter (this module)
    ↓
TG 1.3: Cluster consensus (orientation-validated models only)
```

### Merge into models table

After orientation filtering, run with `--merge` to add `orientation_score` and `orientation_class` columns to `pyrosetta_interface_models.csv`. This enables TG 1.3 to filter on `orientation_class == "pass"` before building consensus.

### Downstream usage

- **TG 1.3** uses only `orientation_class == "pass"` models for consensus
- **TG 1.5** tracks orientation distribution per receptor state
- **TG 1.6** reports orientation statistics in Phase 1 summary

---

## 5. Sheet Residue Definitions

| Sheet | Residues | Role | Experimental evidence |
|-------|----------|------|----------------------|
| 8 | 961, 962, 963, 964 | Active face (primary) | Ala sub. abolished decoy function |
| 9 | 968, 969, 970, 971, 972 | Active face (primary) | Ala sub. abolished decoy function |
| 10 | ~977–980 | Neutral | Ala sub. ≈ wild-type function |
| 11 | ~985–988 | Neutral | Ala sub. ≈ wild-type function |
| 12 | ~993–997 | Structural support | Ala sub. abolished function (role TBD) |

**Note:** Sheets 10–12 boundaries are approximate. They affect supplementary back-face contact counting but NOT the primary orientation score.

---

## 6. Execution Commands

### Server-side orientation filtering (after docking completes)

```bash
conda activate pyrosetta

# All states
python -m egfr_pipeline.phase1.orientation_filter --merge

# Single state
python -m egfr_pipeline.phase1.orientation_filter --state 3GT8_raw --merge
```

### Pilot retroactive validation (server-side only)

```bash
python -m egfr_pipeline.phase1.orientation_filter \
    --pilot_dir /path/to/pilot/final_result/ \
    --pilot_out output/phase1_ppi/orientation_filter_pilot_validation.csv
```

---

## 7. Thresholds and Calibration

| Parameter | Default | Rationale |
|-----------|---------|-----------|
| AMBIGUOUS_BAND | 0.15 | Starting value; calibrate with pilot validation |
| contact_dist | 10.0 Å | Standard for local receptor centroid calculation |
| Probe residues | VAL962, VAL964, SER971 | Span both sheets 8 and 9 |

The ambiguous band (0.15) should be recalibrated after:
1. Pilot retroactive validation (how many pilots pass/fail/ambiguous?)
2. Full docking orientation score distribution (what fraction falls in each class?)
3. Cross-tabulation with energy metrics (do PASS models have better dG/sc on average?)

---

## 8. Limitations

1. **PyRosetta required.** Orientation scoring needs Cα/Cβ coordinates from PDB files via PyRosetta.
2. **Rigid-body assumption.** PCA plane fitting assumes maintained β-sheet geometry (valid for RosettaDock rigid-body docking).
3. **Sheet 12 role unresolved.** Current assumption: structural support, not primary contact face.
4. **Sheets 10–12 boundaries approximate.** Back-face contact counting is informational only.

---

## 9. Validation

### Structural testing (this workspace)

- Module imports and merge logic validated with synthetic data
- 3 models: 1 PASS (score=0.78), 1 FAIL (score=-0.45), 1 AMBIGUOUS (score=0.08)
- Merge correctly adds orientation_score + orientation_class columns to models CSV

### Server-side validation (pending)

- [ ] Run orientation filter on test docking results (1K models)
- [ ] Verify score distribution covers expected range
- [ ] Retroactive pilot validation (5 structures: C02_M01, C02_M03, C04_M01, C04_M02, C07_M03)
- [ ] Cross-tabulate orientation class with dG/dSASA/sc metrics
- [ ] Calibrate ambiguous band threshold if needed
