# Phase 1 Input Validation Report

Generated: 2026-03-11 02:18:22
Task Group: 1.0 (Receptor and Partner Input Preparation)

## Summary

- **Receptor states prepared:** 3
- **Receptor construct type:** full_kinase_domain
- **Partner construct type:** extended_beta_meander
- **Partner residue range:** 955-1006
- **Numbering system:** PDB (3GT8-consistent)

## Receptor Details

| State | Range | Residues | N-lobe | C-lobe | Gaps |
|-------|-------|----------|--------|--------|------|
| 3GT8_raw | 699-1007 | 309 | 139 | 170 | None |
| 3GT8_cl38_48 | 699-1007 | 309 | 139 | 170 | None |
| 3GT8_cl85_100 | 699-1007 | 309 | 139 | 170 | None |

## Partner Details

- **Construct:** Extended beta-meander
- **Source:** TH1 domain structure
- **Range:** 955-1006 (52 residues)
- **First residue:** SER955
- **VAL962 artifact eliminated:** Yes

### Sheet Annotations

| Sheet | Residues | Role | Complete |
|-------|----------|------|----------|
| sheet_8 | [961, 962, 963, 964] | active_face_primary | Yes |
| sheet_9 | [968, 969, 970, 971, 972] | active_face_primary | Yes |
| sheet_10 | [977, 978, 979, 980] | neutral | Yes |
| sheet_11 | [985, 986, 987, 988] | neutral | Yes |
| sheet_12 | [993, 994, 995, 996, 997] | structural_support | Yes |

## Docking Pairs

| State | Receptor | Partner | Excluded |
|-------|----------|---------|----------|
| 3GT8_raw | A:699-1007 (309 res) | B:955-1006 (52 res) | 43 res |
| 3GT8_cl38_48 | A:699-1007 (309 res) | B:955-1006 (52 res) | 43 res |
| 3GT8_cl85_100 | A:699-1007 (309 res) | B:955-1006 (52 res) | 43 res |

## Validation Checks

```
============================================================
RECEPTOR CROSS-VALIDATION
============================================================
PASS: All states use chain A
PASS: All states share range 699-1007
PASS: All states have 309 residues
PASS: 3GT8_raw has N-lobe (139 res)
PASS: 3GT8_cl38_48 has N-lobe (139 res)
PASS: 3GT8_cl85_100 has N-lobe (139 res)
PASS: All states are full_kinase_domain construct

============================================================
PARTNER VALIDATION
============================================================
PASS: First residue is SER955 (VAL962 artifact eliminated)
PASS: Construct starts at 955 (<= 955)
PASS: sheet_8 complete (4/4 residues)
PASS: sheet_9 complete (5/5 residues)
PASS: sheet_10 complete (4/4 residues)
PASS: sheet_11 complete (4/4 residues)
PASS: sheet_12 complete (5/5 residues)
PASS: All sheets (8-12) are complete
PASS: No gaps in construct
```

## Numbering Reference

- All residue numbers in this report use PDB numbering (3GT8-consistent)
- UniProt (P00533) = PDB + 24 (approximate)
- N-lobe/C-lobe boundary: residue 838 (PDB numbering)
- Membrane-proximal excluded residues: see docking_pair_metadata.csv

## Pilot Data

- Legacy pilot data (C-lobe fragment + truncated beta-meander) registered as historical reference
- See pilot_data_reference.csv for details
- Pilot data is NOT used as validation target for new results
