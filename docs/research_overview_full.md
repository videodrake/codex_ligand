# EGFR-MYO1D Research Overview

## 1. Background

### 1.1 Project Scope

EGFR-MYO1D binding site discovery pipeline using computational docking (Vina blind docking + PyRosetta PPI).

### 1.2 Experimental Facts

Three experimentally confirmed facts constrain the computational analysis:

1. **ATP binding maintained under drug treatment**: Drug treatment abolishes kinase activity while ATP binding is preserved. This means ATP binding site pockets are false positives in this project — the drugs do not compete with ATP. Pockets overlapping >50% with ATP_SITE_RESIDUES (37 residues) are automatically excluded with `exclusion_reason = "ATP_site_experimental"`.

2. **Ko et al. alanine substitution (MYO1D beta-meander)**: Sheets 8 (961-964) and 9 (968-972) constitute the active face directly contacting EGFR. Alanine substitution of these residues abolishes function. Sheets 10/11 substitution has no effect (neutral). Sheet 12 substitution also abolishes function (structural support). PPI hotspots must contain >= 3 active face residues to be considered valid.

3. **Ligand chemical diversity**: The three project ligands (173940, 97806, VAX-C12_0) are chemically diverse (all pairwise Tanimoto < 0.4), validating the use of cross-chemical consensus in Vina scoring.

### 1.3 Vina Result Interpretation

ATP site pockets should be excluded from analysis because:
- ATP binding is maintained under drug treatment (experimental evidence)
- Vina blind docking naturally finds strong poses in the deep, well-defined ATP pocket
- These poses are false positives — the drugs act elsewhere on the kinase surface
- See `region_definitions.py` for the 37 ATP_SITE_RESIDUES definition

### 1.4 Methodology Limitations

This pipeline relies on computational predictions with known limitations.
See [`docs/methodology_limitations.md`](methodology_limitations.md) for details on:
- Rigid-body docking (induced fit not modeled)
- LightDock independence (shared input structures)
- Input structure bias (common blind spots)
- Solvent effects (implicit solvation only)
- Vina scoring bias (hydrophobic overestimation)
