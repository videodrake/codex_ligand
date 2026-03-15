# Phase 1 Interface Comparison Report

## Multi-State Receptor-Side Patch Comparison

States compared: 3GT8_raw, EGFR_160-185, EGFR_170-200
Construct type: full_kinase_domain
Evidence: orientation-validated PyRosetta models only

## Per-State Summary

| State | Clusters | Orient-valid models | Hotspot residues (receptor) |
|-------|----------|--------------------|-----------------------------|
| 3GT8_raw | 3 | 1 | 6 |
| EGFR_160-185 | 1 | 4 | 3 |
| EGFR_170-200 | — | — | — |

## Robustness Classification (Receptor-side)

- **Robust** (all 3 states): 0 residues
- **Moderate** (2 states): 3 residues
- **State-specific** (1 state): 4 residues

### Lobe Distribution by Robustness

| Class | N-lobe | C-lobe | Total |
|-------|--------|--------|-------|
| Robust | 0 | 0 | 0 |
| Moderate | 1 | 2 | 3 |
| State-specific | — | — | 4 |

## Partner-side (beta-meander) Robustness

Total partner residues observed: 5
Robust (all states): 0

## Cross-State Numbering Consistency

All states use PDB numbering from 3GT8 crystal structure (699–1007).
MD cluster states (cl38_48, cl85_100) were trimmed to the same range
and chain X was renamed to chain A during TG 1.0 input preparation.

**No numbering mismatches expected.** If residue IDs differ between
states, this indicates genuinely different interface contacts, not
a numbering artifact.

## Interpretation Notes

- **Robust residues** are the strongest candidates for the MYO1D
  attachment patch — they appear regardless of EGFR conformational state.
- **State-specific residues** may indicate conformationally gated
  interaction sites that open only in certain MD-sampled states.
- **N-lobe residues** (< 838) are newly detectable in Phase 1 due to
  the full kinase domain receptor. Their presence/absence is informative.
- All evidence is from **orientation-validated** models only (face-flip
  artifacts excluded).
