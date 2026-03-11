# Orientation-Aware Filtering: Algorithm Design Document

## Version: 1.0
## Date: 2026-03
## Status: Implementation ready, pending server-side validation

---

## 1. Problem Statement

The MYO1D beta-meander is a flat, thin β-sheet ribbon (5 consecutive β-strands, sheets 8–12). When used as a docking partner in global blind docking against the EGFR kinase domain, this geometry creates a systematic vulnerability: the structure can land on the receptor in two fundamentally different orientations that are energetically distinguishable but not always separable by simple residue contact counting.

**Correct orientation:** The active face of sheets 8 and 9 (the experimentally validated binding surface) packs against the receptor. The functional side chains of the validated residues make direct contacts with the receptor surface.

**Flipped orientation:** The back face of the beta-meander faces the receptor. Sheet 8/9 residues may still register as "in contact" because backbone atoms or edge-on contacts can form non-specific interactions, but the functional side chains point into solvent rather than toward the receptor.

The pilot docking campaign (C-lobe fragment, 1M decoys) used a contact-count filter (sheet 8+9 residues ≥ 3 in contact) but did not include an orientation check. This means some accepted poses may have been face-flipped artifacts where the correct residues were in proximity but in the wrong geometric configuration.

---

## 2. Why Contact Counting Alone Fails

In a β-sheet, side chains alternate above and below the sheet plane. For a β-meander with 5 consecutive strands, this means:
- One face has the side chains of even-numbered positions (relative to strand register)
- The other face has the side chains of odd-numbered positions

When the sheet is thin (only 5 strands wide, ~15–20 Å across), the distance between the two faces is only about 5–7 Å. A contact cutoff of 8–10 Å (standard in interface analysis) can capture residues from both faces simultaneously. Therefore, a residue can appear "in contact" whether its side chain faces toward or away from the receptor.

Additionally, backbone hydrogen bonds can form from either face, creating non-specific but energetically favorable contacts that inflate interface metrics without reflecting biologically meaningful binding.

---

## 3. Algorithm Design: Dual-Vector Orientation Test

### 3.1 Overview

The filter computes two vectors and their dot product:

1. **Active-face normal vector:** The outward-pointing normal of the sheet 8/9 plane, oriented toward the active face using a Cα→Cβ probe.

2. **Receptor-direction vector:** A vector from the sheet 8/9 centroid toward the receptor interface centroid (receptor residues near the active face).

The dot product of these two vectors gives the orientation score:
- Positive → active face points toward receptor (PASS)
- Negative → active face points away from receptor (FAIL)
- Near zero → edge-on contact (AMBIGUOUS)

### 3.2 Step-by-step

**Step 1: Collect active-face Cα coordinates**

Extract Cα coordinates for all sheet 8 and sheet 9 residues:
- Sheet 8: VAL961, VAL962, ASN963, VAL964
- Sheet 9: VAL968, GLN969, CYS970, SER971, LEU972

These 9 residues define the active face of the beta-meander.

**Step 2: Compute sheet-plane normal via PCA**

Apply principal component analysis (PCA) to the 9 Cα positions. The two largest principal components span the sheet plane; the smallest principal component is the plane normal.

PCA is preferred over a simple cross product of two chosen vectors because:
- It uses all available data points, not just two arbitrarily selected pairs
- It is robust to local deviations in individual residue positions
- It gives a well-defined normal even when the sheet is slightly curved

**Step 3: Orient the normal toward the active face**

The PCA normal has an ambiguous sign (pointing toward one face or the other). To resolve this:
- Compute the Cα→Cβ vector for **multiple** active-face residues (default: VAL962, VAL964, SER971 — spanning both sheets 8 and 9)
- For each probe, check whether the dot product of (normal) · (Cα→Cβ) is positive or negative
- Take a majority vote: if more probes indicate flipping is needed, flip the normal
- This multi-probe consensus is more robust than a single-residue probe because it averages out rotamer-dependent noise and is not vulnerable to any single residue's local geometry being atypical
- Falls back to single-probe (VAL962) if only one probe residue has valid Cα/Cβ coordinates

**Step 4: Compute receptor-direction vector**

Instead of using the entire receptor centroid (which may be far from the contact area), compute a local receptor centroid:
- Find all receptor Cα atoms within 10 Å of any active-face Cα
- Compute the centroid of these receptor contact residues
- The receptor-direction vector goes from the sheet 8/9 centroid toward this receptor contact centroid

This localized approach ensures the orientation test reflects the actual contact geometry, not the global receptor shape.

**Step 5: Compute dot product and classify**

```
orientation_score = dot(active_face_normal, receptor_direction)

if |score| < 0.15:     → ambiguous (edge-on contact)
elif score > 0:         → pass (active face toward receptor)
else:                   → fail (active face away from receptor)
```

### 3.3 Why this approach is scientifically robust

1. **Grounded in β-sheet geometry.** The alternating side-chain pattern of β-sheets is a fundamental structural property. The Cα→Cβ vector reliably indicates which face a side chain belongs to.

2. **PCA-based plane fitting.** Avoids arbitrary vector choices. Works even if the sheet is slightly twisted or curved (common in real beta-meanders).

3. **Localized receptor centroid.** Prevents false results from using a distant global centroid. In a full kinase domain (~280 residues), the geometric center is far from the C-lobe surface where binding occurs.

4. **Configurable thresholds.** The ambiguous band (|score| < 0.15) prevents hard binary decisions for edge cases. These cases are preserved for manual inspection rather than silently discarded.

5. **Residue-independent.** The filter does not depend on specific contact identities. It evaluates the global geometric orientation of the entire active face, making it complementary to (not redundant with) the existing contact-count filter.

---

## 4. Sheet Residue Definitions

### 4.1 Experimentally validated functional sheets

From Ko et al.:
- **8th β-sheet:** Alanine substitution abolished decoy function → functionally essential
- **9th β-sheet:** Alanine substitution abolished decoy function → functionally essential
- **10th β-sheet:** Alanine substitution had roughly the same decoy function as wild-type → not essential for binding
- **11th β-sheet:** Same as 10th → not essential for binding
- **12th β-sheet:** Alanine substitution abolished decoy function → functionally essential

### 4.2 Structural interpretation for filtering

Sheets 8 and 9 define the **active face**: the primary contact surface that must face the receptor.

Sheet 12 is functionally essential but its role is ambiguous:
- It could be direct contact (in which case it should also face the receptor)
- It could be structural support (stabilizing the beta-meander fold so that sheets 8/9 can bind)

Current working assumption: sheet 12 is structural support, not the primary direct-contact face. This assumption may need revision if new data shows otherwise.

For orientation filtering, sheets 8 and 9 are used to define the active-face plane. Sheet 12 residues are monitored but do not define the orientation criterion.

### 4.3 Residue numbering

All residue numbers are MYO1D numbering (no offset):

| Sheet | Residues | Role |
|-------|----------|------|
| 8 | 961, 962, 963, 964 | Active face (primary) |
| 9 | 968, 969, 970, 971, 972 | Active face (primary) |
| 10 | ~977–980 | Neutral (not essential per experiment) |
| 11 | ~985–988 | Neutral (not essential per experiment) |
| 12 | ~993–997 | Structural support (essential but role TBD) |

Sheet 10–12 boundaries are approximate and should be refined with the actual structure.

---

## 5. Validation Strategy

### 5.1 Retroactive pilot validation

Apply the orientation filter to the 5 existing valid pilot structures:
- C02_M01, C02_M03, C04_M01, C04_M02, C07_M03

Expected outcomes:
- Structures that were manually assessed as "biologically reasonable" in PyMOL should PASS
- If any existing valid structure FAILS, this indicates either:
  (a) the filter threshold is too strict, or
  (b) the pilot structure was actually face-flipped (which would be an important finding)

### 5.2 Synthetic face-flip test

Generate a face-flipped version of a known PASS structure by rotating the beta-meander 180° around the strand axis. Apply the filter:
- Original → should PASS
- Flipped → should FAIL

This provides a positive/negative control pair.

### 5.3 Edge case documentation

Structures classified as AMBIGUOUS should be:
- Visually inspected in PyMOL
- Documented with screenshots showing the sheet 8/9 orientation
- Used to calibrate the ambiguous band threshold

### 5.4 Statistical validation on new docking data

After running new full-kinase-domain docking:
- Report the fraction of decoys in each class (pass/fail/ambiguous)
- Cross-tabulate orientation class with energy metrics (dG, I_sc)
- Verify that PASS models have better interface quality metrics on average than FAIL models (expected if the filter captures real binding geometry)

---

## 6. Integration with Pipeline

### 6.1 Position in workflow

```
Task 1.1: PyRosetta docking
    ↓
Task 1.2: Interface residue extraction (all models)
    ↓
★ Task 1.2A: Orientation filter (this algorithm) ← HERE
    ↓
Task 1.3: Cluster consensus (orientation-validated models only)
    ↓
Task 1.5: Multi-state comparison
```

### 6.2 Output schema

`orientation_filter_log.csv`:

| Column | Type | Description |
|--------|------|-------------|
| filename | str | PDB filename |
| orientation_score | float | Dot product (-1 to +1) |
| orientation_class | str | pass / fail / ambiguous / error codes |
| n_active_face_ca | int | Number of active-face Cα found |
| n_back_face_ca | int | Number of back-face Cα found |
| sheet_centroid | str | [x, y, z] of sheet 8/9 centroid |
| normal_vector | str | [x, y, z] of active-face normal |
| receptor_direction | str | [x, y, z] toward receptor contact centroid |
| error | str | Error message if any |

### 6.3 Downstream usage

- Only models with `orientation_class == "pass"` enter consensus building (Task 1.3)
- Models with `orientation_class == "ambiguous"` are flagged for manual review
- Models with `orientation_class == "fail"` are excluded from consensus but preserved in raw data
- The orientation score is available for downstream analysis (e.g., correlation with energy metrics)

---

## 7. Limitations and Caveats

1. **Multi-probe consensus residue availability.** The multi-probe approach (VAL962, VAL964, SER971) requires at least 2 of 3 residues to have valid Cα/Cβ coordinates. If the docked structure has missing atoms in these residues, the consensus degrades to single-probe or fails. Mitigation: these residues are core β-sheet residues unlikely to have missing atoms in properly prepared structures.

2. **Rigid-body assumption.** The filter assumes the beta-meander maintains its β-sheet geometry during docking. If the sheet is distorted (e.g., by Rosetta flexible backbone moves), the PCA-based plane fitting may be less accurate. Mitigation: for standard RosettaDock rigid-body docking, this assumption holds well.

3. **Threshold sensitivity.** The ambiguous band (0.15) is a starting value. It should be calibrated using the retroactive pilot validation and the synthetic flip test. The threshold may need adjustment after examining the actual score distribution.

4. **Sheet 12 role unresolved.** If sheet 12 turns out to be a direct-contact face rather than structural support, the active-face definition would need to be expanded. The current algorithm can be extended by adding sheet 12 residues to `ACTIVE_FACE_RESIDUES`.

---

## 8. Korean Summary (간단 요약)

이 문서는 beta-meander의 face-flip 감지를 위한 orientation filter 알고리즘을 정의한다.

핵심 원리:
- β-sheet의 곁사슬은 sheet 평면 위아래로 교대 배치됨
- Sheet 8/9의 active face가 receptor를 향하는지 반대로 뒤집혔는지를 판별
- PCA로 sheet 평면 법선벡터를 구하고, Cα→Cβ probe로 active face 방향을 결정
- 법선벡터와 receptor 방향 벡터의 dot product가 양수면 PASS, 음수면 FAIL
- 0 근처는 AMBIGUOUS로 분류하여 수동 검토

기존 contact count 필터만으로는 face-flip을 감지할 수 없다. 이 필터는 contact count와 독립적이며 보완적이다.

