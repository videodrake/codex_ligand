## 2. Materials and Methods

### 2.1 Protein Structure Preparation

**EGFR kinase domain.** The crystal structure of the human EGFR kinase domain in complex with an inhibitor was retrieved from the Protein Data Bank (PDB entry: 3GT8; resolution 2.00 Å). Residues 699–1007 (full kinase domain) were retained. Ligands, water molecules, and non-protein atoms were removed. Missing side chains and loops were rebuilt using the Rosetta FastRelax protocol prior to docking. This structure served as the ground-state receptor conformation (hereafter: 3GT8_raw).

**MYO1D TH1 domain.** No experimental crystal structure is available for the TH1 domain of human Myosin-1D. The three-dimensional model used in all docking calculations was obtained from the AlphaFold Protein Structure Database (UniProt: O94832; AlphaFold model AF-O94832-F1) [12]. Residues 955–1006, corresponding to the beta-meander region of the TH1 domain that encompasses experimentally validated sheets 8 and 9 (the "active face"), were extracted and used as the docking partner [9].

### 2.2 Molecular Dynamics Simulation and Conformational Sampling

A 200 ns all-atom molecular dynamics simulation of the EGFR kinase domain (residues 699–1007) was conducted using GROMACS (version [X.X]) with the CHARMM36 force field [13]. The protein was solvated in a cubic TIP3P water box with a minimum solute-to-box-wall distance of 12 Å. System charge was neutralized by the addition of counter ions (Na⁺/Cl⁻) to a final concentration of 0.15 M. Energy minimization was performed with the steepest descent algorithm until the maximum force fell below 1,000 kJ mol⁻¹ nm⁻¹. The system was equilibrated in the NVT ensemble for 100 ps at [temperature] K, followed by NPT equilibration for 100 ps at [temperature] K and [pressure] bar. Production MD was run for 200 ns under NPT conditions using a 2 fs integration timestep; covalent bonds involving hydrogen were constrained with the LINCS algorithm. Long-range electrostatics were treated by particle-mesh Ewald (PME) summation with a real-space cutoff of 12 Å. Coordinates were saved every 10 ps.

Trajectory frames were clustered using the GROMOS algorithm (RMSD cutoff: 2.0 Å on Cα atoms of residues 699–1007). Two representative cluster centroids were extracted: one from the 38–48 ns interval (EGFR_160-185) and one from the 85–100 ns interval (EGFR_170-200), reflecting distinct kinase domain conformations. Together with 3GT8_raw, these three structures constituted the receptor ensemble for all subsequent calculations.

### 2.3 Global Protein–Protein Docking with PyRosetta

Protein–protein docking was performed using PyRosetta (version 4) [14] implementing the RosettaDock protocol. For each of the three receptor conformations, 10 independent docking runs were executed, each generating 20,000 structural models (total: 600,000 models). Each run was initialized with a distinct random seed to ensure statistical independence of sampling. The MYO1D beta-meander fragment (residues 955–1006) was used as the partner in all calculations. Docking energy was evaluated using the REF2015 score function; binding energy (ΔG_bind) was computed as the difference in total score between the complex and the sum of isolated chain scores after minimization.

### 2.4 Orientation Filter

To enforce biologically plausible docking geometries consistent with the experimentally characterized MYO1D active face, all docking poses were subjected to a dual-vector orientation filter prior to downstream analysis. For each pose, principal component analysis (PCA) was applied to the Cα coordinates of MYO1D sheets 8 and 9 (residues 961–972) to define the sheet plane normal vector (**n**_sheet). A partner-contact vector (**v**_contact) was computed as the unit vector from the centroid of EGFR interface residues to the centroid of MYO1D active-face residues. A pose was classified as orientation-valid if the dot product **n**_sheet · **v**_contact exceeded a threshold of 0.10 (AMBIGUOUS_BAND = 0.10), ensuring that the active face is presented toward the EGFR surface. Poses failing this criterion were excluded from hotspot analysis.

### 2.5 Hotspot Residue Identification

For each docking run, interface residues were identified as those with at least one heavy atom within 5.0 Å of any heavy atom of the opposing chain. Per-residue cluster occupancy was computed as the fraction of orientation-valid models within a cluster in which the residue appeared at the interface. A residue was designated a hotspot if its occupancy reached or exceeded 0.50 in at least one conformational state. Hotspot designations were aggregated across all three receptor states, and residues classified as hotspots in all three states were designated "robust hotspots" of the highest confidence tier.

### 2.6 LightDock Cross-Validation

To provide an independent validation of the PyRosetta interface predictions, global protein–protein docking was independently conducted using LightDock (version 0.9) with the DFIRE2 statistical energy function [15]. The same three receptor conformations and the same MYO1D fragment were used. For each receptor state, 200 swarms of 25 glowworms each were simulated for 200 steps. Top-scoring poses were clustered and their interface residues were extracted. Method agreement between PyRosetta and LightDock was assessed at the residue level: a residue was classified as "both-method" if it appeared in the interface in at least one pose from each method.

### 2.7 Surface Pocket Detection

All-surface pocket detection was performed using fpocket (version 4.0) [16] on each of the three receptor conformations. Raw pocket outputs from the three structures were merged across states using a centroid-distance criterion (overlap threshold: 4.0 Å between pocket centroids), yielding a non-redundant set of candidate pockets. A total of 165 raw pockets were detected and consolidated into 103 merged pocket identities for downstream analysis.

### 2.8 Dual-Workflow Small-Molecule Docking

Small-molecule docking was conducted with AutoDock Vina (version 1.2) [17] using three structurally diverse ligands (pairwise Tanimoto similarity < 0.4) provided for pocket characterization purposes. Receptor and ligand structures were prepared in PDBQT format.

**Workflow A (Blind docking, WF-A):** Each ligand was docked to the entire surface of each receptor conformation using a blind search box encompassing the full kinase domain. Exhaustiveness was set to 32. The best-scoring pose per ligand–pocket combination was retained. Pocket stability was assessed as the fraction of ligand poses within 2.0 Å RMSD of the top-ranked pose across independent runs (pocket_stability metric). Pockets were scored on a five-component Vina quality axis: affinity, convergence, stability, diversity, and ligand consensus. A final site verdict (STRONG / MODERATE / WEAK) was assigned by integrating Vina quality, PPI proximity, and cross-receptor consistency scores.

**Workflow B (PPI-guided focused docking, WF-B):** The PPI patch centroid (Cα-based geometric mean of the 15 highest-occupancy EGFR hotspot residues) was computed for each receptor state. Only fpocket-detected cavities with centroids within 25 Å of this reference point were included in focused docking. Exhaustiveness was set to 128. Docking results were integrated with the PPI interface map via a four-axis perturbation scoring framework (Section 2.9).

### 2.9 Perturbation Scoring Framework

Each candidate pocket was assigned a composite perturbation score (range 0–1) integrating four independent axes, with the following weights:

| Axis | Description | Weight |
|------|-------------|--------|
| A1: PPI Interface | Proximity and residue overlap with hotspot patch | 30% |
| A2: Druggability | Pocket volume, hydrophobicity, shape complementarity | 25% |
| A3: Perturbation Relevance | Ligand docking affinity and convergence in focused runs | 30% |
| A4: State Robustness | Pocket persistence across all three receptor states | 15% |

Pocket–PPI relationships were classified into four categories based on the fraction of hotspot residues in the pocket lining and the centroid-to-patch distance: orthosteric (≥25% hotspot overlap and ≥2 overlapping residues), rim candidate (≥1 overlapping residue), allosteric candidate (centroid ≤20 Å from patch centroid, zero direct overlap), and low-relevance (>20 Å, zero overlap). ATP-binding pockets were identified by comparison against the crystallographic ATP-binding site and excluded from further consideration, consistent with experimental evidence that ATP binding is maintained in the context of MYO1D interaction [9].

### 2.10 Software and Computational Environment

All Python analysis scripts were executed under Python 3.9 with NumPy, SciPy, and pandas. Molecular dynamics simulations and Vina docking runs were executed on a university high-performance computing (HPC) cluster using PBS job scheduling (qsub). Structural visualization was performed in PyMOL [18].
