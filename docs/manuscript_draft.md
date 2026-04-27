# Computational Identification of Druggable Pockets at the EGFR–MYO1D Protein–Protein Interface: A Multi-Conformational Ensemble Docking and Perturbation Scoring Approach

---

**Authors:** [Author 1, Affiliation]¹; [Author 2, Affiliation]¹; [Corresponding Author, Affiliation]¹*

¹ [Department, Institution, City, Country]

\* Correspondence: [email@institution.edu]

---

## Abstract

Epidermal growth factor receptor (EGFR) signaling is dysregulated in numerous cancers, yet resistance to ATP-competitive kinase inhibitors remains a persistent clinical challenge. The non-catalytic protein–protein interaction (PPI) between the EGFR kinase domain and the tail homology 1 (TH1) domain of Myosin-1D (MYO1D) represents an alternative regulatory axis that may be exploitable for therapeutic intervention. Here we report a comprehensive computational pipeline to map the EGFR–MYO1D binding interface and discover druggable pockets capable of disrupting this interaction. Three EGFR conformational states were sampled — a crystal structure (PDB: 3GT8) and two molecular dynamics (MD) cluster representatives extracted from a 200 ns GROMACS/CHARMM36 simulation — providing an ensemble that captures receptor flexibility. Global protein–protein docking with PyRosetta (600,000 models across 30 independent runs) was cross-validated by LightDock (DFIRE2 scoring), and an orientation filter based on principal component analysis (PCA) of the MYO1D beta-meander sheet plane was applied to retain only physically plausible docking poses. Twelve EGFR C-lobe residues (ILE941, THR940, GLN935, PRO934, PRO937, SER957, VAL980, GLN982, ARG986, HIS988, ARG977, THR993) emerged as robust hotspots present across all three receptor states and confirmed by both docking methods. Cavity detection with fpocket (165 raw → 103 merged pockets) followed by a dual-workflow AutoDock Vina docking strategy identified two tier-1 pockets as priority modulator candidates: PKT07 (perturbation score 0.541, interface rim, always accessible) and PKT34 (perturbation score 0.492, allosteric, state-shifted). ATP-competitive pockets were excluded on experimental grounds. These results provide a structural blueprint for developing non-ATP small-molecule PPI modulators of the EGFR–MYO1D axis.

**Keywords:** EGFR; MYO1D; protein–protein interaction; druggable pocket; molecular docking; PyRosetta; AutoDock Vina; ensemble docking; allosteric modulation; computational drug discovery
## 1. Introduction

Epidermal growth factor receptor (EGFR) is a receptor tyrosine kinase that occupies a central node in oncogenic signaling, driving cell proliferation, survival, and differentiation in epithelial tissues [1]. Gain-of-function mutations and amplification of EGFR are among the most frequently observed molecular alterations in non-small cell lung cancer, colorectal cancer, glioblastoma, and head-and-neck carcinoma [2]. First- and second-generation ATP-competitive kinase inhibitors (gefitinib, erlotinib, afatinib) achieve initial disease control in EGFR-mutant tumors; however, acquired resistance — most commonly mediated by the T790M gatekeeper mutation or bypass signaling through alternative receptor tyrosine kinases — limits their long-term efficacy [3,4]. Third-generation inhibitors (osimertinib) address T790M but are themselves susceptible to further resistance via C797S mutation or amplification of downstream effectors [5]. This evolutionary pattern underscores the need for mechanistically distinct therapeutic approaches that do not depend on competitive occupation of the ATP-binding pocket.

Protein–protein interactions (PPIs) involving EGFR represent a largely underexplored class of drug targets. PPIs govern receptor dimerization, adaptor protein recruitment, and co-receptor modulation, and disrupting these contacts can suppress signaling independently of kinase catalytic activity [6,7]. Among the EGFR-interacting proteins recently identified, Myosin-1D (MYO1D) — a class I non-muscle myosin — has emerged as a functionally relevant binding partner. MYO1D belongs to a family of membrane-associated motor proteins and contains a tail homology 1 (TH1) domain that mediates protein–protein contacts at the cytoplasmic face of the plasma membrane [8]. Ko et al. demonstrated by alanine-scanning mutagenesis that residues located on sheets 8 and 9 of the MYO1D TH1 beta-meander (the "active face") are critical determinants of EGFR binding affinity, providing direct experimental evidence that this contact surface drives the interaction [9]. The functional significance of the EGFR–MYO1D axis in receptor internalization and signal duration makes it an attractive non-catalytic target for pharmacological intervention.

Despite this experimental foundation, the precise residue-level binding interface and the identity of tractable small-molecule pockets at or near this interface have not been established. Structure-based drug discovery targeting PPIs is hampered by the typically flat and featureless nature of protein–protein contact surfaces [10]; however, pockets adjacent to the interface ("rim sites") or allosterically coupled to it can provide binding geometries amenable to small-molecule modulation [11]. Identifying such pockets requires both accurate mapping of the PPI interface and systematic assessment of surface cavities in the context of receptor conformational diversity.

Here we present a multi-stage computational pipeline that addresses these challenges. We employed PyRosetta-based global protein–protein docking across three EGFR conformational states — a crystallographic ground-state structure and two distinct conformers extracted from a 200 ns molecular dynamics (MD) trajectory — generating 600,000 docking models in total. An orientation filter based on principal component analysis (PCA) of the MYO1D beta-sheet geometry was used to select poses consistent with the experimentally defined active face [9]. Hotspot residues were identified by cluster-consensus occupancy analysis and cross-validated by an independent LightDock calculation using the DFIRE2 energy function. Surface pockets were detected by fpocket and subjected to AutoDock Vina docking in two complementary workflows: a blind, receptor-surface-wide search (Workflow A) and a PPI-guided focused search (Workflow B). A four-axis perturbation scoring framework integrating PPI interface proximity, druggability, perturbation relevance, and conformational robustness ranked all candidate pockets to prioritize candidates for experimental follow-up.

This work provides, to our knowledge, the first systematic structural characterization of the EGFR–MYO1D binding interface at atomic resolution and identifies two tier-1 druggable pocket candidates — a rim modulator site (PKT07) and an adjacent allosteric site (PKT34) — that represent starting points for structure-based design of non-ATP EGFR modulators targeting the MYO1D interaction axis.
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
## 3. Results

### 3.1 Three-State Receptor Ensemble Captures Kinase Domain Flexibility

A 200 ns MD simulation of the EGFR kinase domain revealed substantial conformational variability, particularly in the C-lobe surface loops that constitute the predicted MYO1D binding region. GROMOS clustering identified two well-populated conformational states distinct from the crystallographic structure: EGFR_160-185, representing a cluster dominant between 38 and 48 ns, and EGFR_170-200, representing a cluster dominant between 85 and 100 ns. Together with the 3GT8_raw crystal structure, these three states formed the receptor ensemble used throughout the study. The inclusion of MD-derived conformers ensured that docking calculations sampled receptor geometries not represented in the static crystal structure, which is essential for characterizing a flexible C-lobe surface.

### 3.2 Global PPI Docking Identifies a C-Lobe Binding Interface

PyRosetta global docking across the three receptor states yielded 600,000 structural models (3 states × 10 independent seeds × 20,000 models per run). After application of the orientation filter — which required that the MYO1D active face (sheets 8–9, residues 961–972) be directed toward the EGFR surface — 129 to 153 EGFR-side residues were identified as interface contacts per conformational state (Table 1). The best binding energy (ΔG_bind) observed was −18.4 kcal/mol for the EGFR_160-185 state; mean ΔG_bind across orientation-valid top clusters ranged from −11.9 to −13.0 kcal/mol across the three states (Table 1).

**Table 1. PyRosetta docking summary across three EGFR conformational states.**

| Receptor State | Runs (seeds) | Models | EGFR Interface Residues | Best ΔG_bind (kcal/mol) | Mean ΔG_bind (kcal/mol) |
|---|---|---|---|---|---|
| 3GT8_raw | 10 | 200,000 | 129 | −18.2 | −11.9 |
| EGFR_160-185 | 10 | 200,000 | 141 | −18.4 | −13.0 |
| EGFR_170-200 | 10 | 200,000 | 153 | −18.0 | −11.9 |

Interface contacts were dominated by the EGFR C-lobe (residues 699–1007 C-terminal half), with 51–74 C-lobe residues identified per state versus 28–31 N-lobe residues, confirming that MYO1D engagement is primarily directed to the kinase C-lobe surface. On the MYO1D side, recurrently observed contact residues included SER955, VAL962, VAL964, HIS973, PRO990, and ASN1006 — all located within or immediately flanking the experimentally validated sheets 8 and 9 active face [9] — corroborating the structural basis of the Ko et al. mutagenesis data.

### 3.3 Twelve Robust Hotspot Residues Define the Core Binding Epitope

Hotspot analysis (cluster-consensus occupancy ≥ 0.50) across all three receptor states and cross-validation by independent LightDock/DFIRE2 docking identified 12 EGFR residues that were consistently classified as hotspots in all three conformational states and confirmed by both docking methods (Figure 1; Table 2). All 12 residues lie on the EGFR C-lobe surface, forming a contiguous patch spanning approximately residues 930–993:

**Table 2. Twelve robust EGFR hotspot residues identified across all three receptor states by both PyRosetta and LightDock.**

| Residue | Position | Max Occupancy | States (n/3) | Method Agreement |
|---|---|---|---|---|
| ILE941 | C-lobe | 1.00 | 3/3 | Both |
| ARG977 | C-lobe | 0.75 | 3/3 | Both |
| THR993 | C-lobe | 0.71 | 3/3 | Both |
| ARG986 | C-lobe | 0.60 | 3/3 | Both |
| GLN982 | C-lobe | 0.60 | 3/3 | Both |
| VAL980 | C-lobe | 0.57 | 3/3 | Both |
| GLN935 | C-lobe | 0.50 | 3/3 | Both |
| HIS988 | C-lobe | 0.50 | 3/3 | Both |
| PRO934 | C-lobe | 0.50 | 3/3 | Both |
| PRO937 | C-lobe | 0.50 | 3/3 | Both |
| SER957 | C-lobe | 0.50 | 3/3 | Both |
| THR940 | C-lobe | 0.50 | 3/3 | Both |

ILE941 displayed the highest occupancy (1.00) across orientation-valid models, indicating near-universal participation in top-ranked docking poses. The patch encompasses a stretch of the C-lobe surface that is exposed and solvent-accessible in all three receptor conformations, consistent with its role as a protein–protein recognition surface rather than a buried catalytic element. An additional 28 residues were identified as interface contacts in all three states with both methods but fell below the hotspot occupancy threshold, extending the predicted binding footprint to approximately 40 residues in the highest-confidence tier.

### 3.4 Surface Pocket Landscape of the EGFR Kinase Domain

fpocket analysis of the three receptor structures identified 165 raw surface pockets, which were consolidated to 103 non-redundant pocket identities after cross-state merging. Pocket classification relative to the PPI hotspot patch yielded the following distribution across the ensemble:

- **Rim candidates** (≥1 hotspot residue in pocket lining): pockets with direct but partial contact with the hotspot patch, bordering the PPI interface
- **Allosteric candidates** (centroid ≤20 Å from patch centroid, zero direct hotspot overlap): pockets adjacent to the interface without direct residue contact
- **Low-relevance** (centroid >20 Å from patch centroid, zero hotspot overlap): surface pockets remote from the PPI region

The majority of high-scoring pockets by fpocket druggability metrics were located either at the ATP-binding site (N-lobe, excluded on experimental grounds; see Section 2.9) or in the C-lobe surface region proximal to the hotspot patch, validating the PPI-guided cavity search strategy.

### 3.5 Workflow A Blind Docking Reveals an Intrinsic Limitation of Unguided Searches

Blind Vina docking across the full kinase domain surface (WF-A) produced the following top-ranked non-ATP sites (Table 3). Although three pockets overlapping the ATP-binding region achieved the strongest Vina affinities (−9.0 to −9.8 kcal/mol), these were excluded from therapeutic consideration. The highest-ranked non-ATP, non-excluded site (P045, EGFR_170-200) received a STRONG verdict (composite score 56.4/100) but exhibited a centroid distance of 62.7 Å from the PPI hotspot patch — too distal for direct PPI perturbation. The best WF-A site with any PPI proximity was located 14.5–16.3 Å from the hotspot centroid but coincided with the ATP-binding region.

**Table 3. Top-ranked Workflow A (blind docking) candidate sites.**

| Site ID | State | Verdict | Score | Best Affinity (kcal/mol) | Distance to PPI (Å) | ATP Site |
|---|---|---|---|---|---|---|
| P045 | EGFR_170-200 | STRONG | 56.4 | −6.9 | 62.7 | No |
| P010 | EGFR_160-185 | MODERATE | 74.0 | −9.1 | 14.9 | Yes (excluded) |
| P004 | EGFR_170-200 | MODERATE | 70.7 | −9.0 | 14.5 | Yes (excluded) |
| P003 | 3GT8_raw | MODERATE | 60.0 | −9.8 | 14.3 | Yes (excluded) |

This analysis demonstrates that unguided blind docking cannot reliably identify PPI-proximal druggable pockets on the EGFR C-lobe surface when the ATP site dominates the druggability landscape, motivating the PPI-first focused approach of Workflow B.

### 3.6 Workflow B PPI-Guided Focused Docking Identifies Two Tier-1 Priority Pockets

Focused Vina docking restricted to pockets within 25 Å of the PPI patch centroid (WF-B), followed by four-axis perturbation scoring, yielded a ranked list of 103 candidate pocket–ligand combinations. Two pockets emerged as tier-1 druggability candidates with distinct mechanistic profiles (Table 4):

**Table 4. Top-ranked Workflow B (PPI-guided) pocket candidates by perturbation score.**

| Pocket | State | Score | Relationship | Tier | Hotspot Overlap | Centroid Distance (Å) | State Class |
|---|---|---|---|---|---|---|---|
| PKT07 | 3GT8_raw | 0.541 | Rim candidate | Tier 1 | 1 (LEU1001) | 18.7 | State-robust |
| PKT34 | EGFR_170-200 | 0.492 | Allosteric candidate | Tier 1 | 0 | 9.4 | State-shifted |
| PKT02 | EGFR_160-185 | 0.433 | Rim candidate | Tier 3 | 2 | — | State-shifted |
| PKT10 | 3GT8_raw | 0.431 | Allosteric candidate | Tier 2 | 0 | 16.8 | State-robust |

**PKT07 — Interface rim modulator candidate.** PKT07 is a large surface pocket (42 lining residues) located on the C-lobe, adjacent to the hotspot patch. It shares one direct hotspot residue (LEU1001) with the PPI interface, classifying it as a rim candidate. Its four-axis profile was: A1 (PPI interface) = 0.280, A2 (druggability) = 0.684, A3 (perturbation relevance) = 0.630, A4 (state robustness) = 0.650. PKT07 is present in all three receptor conformations (state_robust_pocket, always_accessible), indicating that the cavity is not an artifact of a single receptor snapshot. All three docking ligands achieved binding in PKT07 with strong ligand support, and cross-state docking was consistent across receptor conformers.

**PKT34 — Allosteric modulator candidate.** PKT34 is a tier-1 druggable cavity (EGFR_170-200 conformation) positioned 9.4 Å from the PPI hotspot patch centroid, with zero direct hotspot residue overlap, classifying it as an allosteric candidate by geometry. Its four-axis profile was: A1 = 0.253, A2 = 0.645, A3 = 0.525, A4 = 0.650. Notably, PKT34 was detectable only in the MD-derived EGFR_170-200 conformation (state-shifted pocket), illustrating the importance of ensemble-based receptor sampling: this cavity is absent from or inaccessible in the crystal structure and was discoverable only because conformational flexibility during simulation opened the pocket. All three ligands docked to PKT34 with strong consensus.

A critical observation is that the highest direct hotspot overlap (PKT17, 3 overlapping residues, EGFR_170-200) was associated with a tier-3 druggability classification, reflecting shallow and geometrically unfavorable cavity geometry despite proximity to the interface. This inverse relationship between direct hotspot overlap and druggability (PKT07: 1 overlap, tier 1; PKT17: 3 overlap, tier 3) highlights a fundamental geometric tension in PPI drug discovery and underscores the value of rim and allosteric sites for small-molecule design.

### 3.7 Conformational State Dependence of Pocket Accessibility

Cross-state analysis revealed two distinct pocket accessibility patterns. State-robust pockets (present in all three conformations, e.g., PKT07) represent constitutively accessible binding sites that can be targeted irrespective of the kinase conformational state. State-shifted pockets (present in one or two MD-derived conformations but not the crystal structure, e.g., PKT34) represent cryptic sites that are induced by thermal fluctuation and accessible only in specific conformational states. The existence of PKT34 exclusively in the EGFR_170-200 MD cluster demonstrates that the C-lobe surface of EGFR is dynamically capable of presenting druggable cavities not visible in static structures, a finding with direct implications for structure-based drug design campaigns targeting this interface.
## 4. Discussion

### 4.1 The EGFR C-Lobe Surface as a Non-Catalytic Drug Target

The twelve robust hotspot residues identified here — ILE941, THR940, GLN935, PRO934, PRO937, SER957, VAL980, GLN982, ARG986, HIS988, ARG977, and THR993 — form a contiguous epitope on the EGFR C-lobe that is geometrically and physicochemically distinct from the ATP-binding cleft. This distinction is therapeutically significant. ATP-competitive inhibitors are subject to competitive displacement by the high intracellular ATP concentration (~1–5 mM) and to resistance mutations (T790M, C797S) that alter the binding geometry of the hinge region [3,5]. The C-lobe PPI surface, by contrast, is under no analogous evolutionary pressure for ATP coordination, and its disruption by a small molecule would be expected to impair MYO1D engagement independently of kinase catalytic status. The finding that ILE941 achieves full occupancy (1.00) across orientation-valid docking poses identifies this residue as a particularly promising anchor point for interface mapping and future mutagenesis validation.

The consensus between PyRosetta and LightDock — two methods using fundamentally different energy functions (physics-based REF2015 versus statistical DFIRE2) — strengthens confidence in the predicted binding epitope. Method agreement at the residue level provides a conservative estimate of interface contacts; the 12 dual-method hotspots represent the highest-confidence subset of a larger interface footprint comprising approximately 40 residues across all three receptor states. The consistency of the C-lobe localization across all three independent receptor conformations further supports the robustness of this prediction against receptor conformational uncertainty.

### 4.2 Limitations of Blind Docking for PPI-Proximal Pocket Discovery

A central finding of Workflow A is that unguided blind docking across the full kinase domain surface failed to recover PPI-proximal druggable pockets when the search was not constrained by prior knowledge of the binding interface. The top-ranked non-ATP site in WF-A (P045) was located 62.7 Å from the PPI hotspot centroid — far beyond any plausible allosteric coupling distance. This outcome is not a failure of the docking algorithm but reflects an intrinsic property of the EGFR surface: the ATP-binding cleft is a pre-formed, geometrically deep, and highly druggable cavity that dominates blind scoring landscapes. Without prior information about the PPI interface, computational screens would consistently prioritize this known pocket over shallow C-lobe surface features.

This finding has methodological implications for PPI-targeted drug discovery programs more broadly. Our results suggest that PPI interface mapping must precede, rather than accompany, small-molecule pocket identification for targets where a well-characterized active site coexists with a PPI surface on the same protein. The PPI-first Workflow B strategy — in which the interface is defined by ensemble docking before pockets are evaluated — recovered both tier-1 candidates (PKT07 and PKT34) that would have been invisible to a blind search, and specifically identified PKT34, a cryptic allosteric cavity that is only accessible in the MD-derived EGFR_170-200 conformer.

### 4.3 PKT07 as an Interface Rim Modulator Candidate

PKT07's classification as a rim candidate rests on a single shared residue with the hotspot patch (LEU1001). While a direct overlap of one residue may appear modest, rim sites — pockets that border rather than coincide with the PPI interface — are among the most tractable druggable geometries for PPI disruption [11]. A small molecule occupying PKT07 could sterically occlude the approach of the MYO1D active face to the EGFR C-lobe surface without requiring deep burial at the protein–protein contact plane. The druggability tier-1 classification of PKT07 reflects favorable pocket geometry: 42 lining residues provide an extensive binding surface, and the pocket's state-robust nature (present in all three receptor conformations) ensures that the cavity is accessible for ligand engagement regardless of the kinase conformational state at the moment of binding.

The highest axis score for PKT07 is A2 (druggability, 0.684), indicating that the pocket's physicochemical properties are well-matched to small-molecule binding, while A1 (PPI interface proximity, 0.280) reflects the partial, rim-type overlap with the hotspot. Future medicinal chemistry efforts targeting PKT07 should aim to extend ligand contacts toward the hotspot core — particularly toward ILE941, ARG977, and THR993 — to maximize the probability of steric interference with MYO1D sheet 8/9 engagement.

### 4.4 PKT34 as a Cryptic Allosteric Candidate

PKT34 presents a qualitatively different therapeutic opportunity. Its centroid lies 9.4 Å from the PPI patch centroid — closer than any other tier-1 druggable site — yet it shares no direct hotspot residues with the binding patch. This geometry is consistent with an allosteric mechanism: a ligand occupying PKT34 could alter the local conformation of the C-lobe surface loops that form the MYO1D binding epitope without directly competing for the same residue contacts. Allosteric PPI modulators have the additional advantage of not requiring the ligand to match the binding footprint of the protein partner, which is often large and geometrically inaccessible [19].

The state-shifted nature of PKT34 — accessible only in the EGFR_170-200 MD conformer — is a double-edged characteristic. On one hand, it indicates that the pocket does not exist as a pre-formed cavity in the crystallographic ground state, which complicates direct structure-based design. On the other hand, cryptic pockets that open transiently during protein dynamics are increasingly recognized as productive drug targets [20]; ligands that preferentially bind and stabilize the open-pocket conformation can shift the conformational equilibrium toward the accessible state. The EGFR_170-200 cluster represents a thermally accessible conformation (populated during a 15 ns window of a 200 ns simulation), suggesting that PKT34 is not an extreme outlier state but a recurrently visited accessible geometry.

### 4.5 Geometric Tension Between Direct Overlap and Druggability

A notable finding of this study is the inverse relationship between direct hotspot overlap and druggability scores. PKT17 (EGFR_170-200), which shares three hotspot residues with the PPI patch — the highest direct overlap of any pocket in the ensemble — received a tier-3 druggability classification because its cavity geometry is shallow and geometrically unfavorable for small-molecule binding. By contrast, PKT07 (one overlap, tier 1) and PKT34 (zero direct overlap, tier 1) present far superior pharmacophore geometries. This observation recapitulates a broader challenge in PPI drug discovery: the protein–protein contact plane itself is often too flat and featureless for high-affinity small-molecule engagement, whereas the flanking pockets — though geometrically adjacent — offer better prospects [10,11]. Our scoring framework captures this tension explicitly through the separation of A1 (PPI interface proximity) and A2 (druggability) axes, enabling rational prioritization of sites that balance interface relevance with ligand-binding tractability.

### 4.6 Experimental Validation Priorities

The computational results reported here define a hierarchy of experimental priorities. First, alanine substitution of the twelve robust hotspot residues — particularly ILE941, ARG977, THR993, and ARG986, which exhibit the highest occupancy scores — would establish which EGFR-side contacts are energetically essential for MYO1D binding and thus the most important to occlude pharmacologically. Second, fragment-based screening or surface plasmon resonance (SPR) binding assays directed at PKT07 and PKT34 would provide direct evidence of small-molecule engagement at these sites. Third, co-immunoprecipitation or proximity ligation assays using cell lines with confirmed EGFR–MYO1D co-expression would enable functional readout of PPI disruption in a cellular context. The conformational selectivity of PKT34 suggests that MD-based ensemble docking or induced-fit docking protocols should be employed in any experimental structure-guided optimization campaign.

### 4.7 Limitations and Caveats

Several limitations of the present study merit acknowledgment. The MYO1D TH1 domain structure was derived from AlphaFold prediction rather than experimental crystallography; while AlphaFold models have demonstrated high accuracy for well-folded domains [12], the precise side-chain orientations at the binding interface may differ from the experimental structure, introducing uncertainty in contact geometry. The docking calculations employed a rigid receptor approximation within each conformational snapshot; full receptor flexibility during docking (induced-fit effects) was captured only implicitly through the three-state ensemble. Pocket centroid distances reported here are based on Cα-based geometric means of residue sets and systematically overestimate the distance between the nearest pocket surface atom and the nearest interface atom; true proximity is likely shorter than the reported centroids suggest. Finally, the perturbation scoring framework is a computational hypothesis-ranking tool and does not predict experimental binding affinity or in-cell activity; all candidate rankings require experimental validation before biological conclusions can be drawn.
## 5. Conclusion

This study presents the first computational characterization of the EGFR–MYO1D protein–protein interface at atomic resolution and provides a ranked set of druggable pocket candidates for small-molecule intervention. By combining global PyRosetta docking (600,000 models across three EGFR conformational states), orientation-filtered hotspot analysis, LightDock cross-validation, fpocket-based cavity detection, dual-workflow AutoDock Vina docking, and a four-axis perturbation scoring framework, we establish the following principal findings:

1. **The EGFR–MYO1D binding interface is localized to the C-lobe surface** and is defined by a core epitope of twelve robust hotspot residues — ILE941, THR940, GLN935, PRO934, PRO937, SER957, VAL980, GLN982, ARG986, HIS988, ARG977, and THR993 — that are consistently identified across all three receptor conformations and confirmed by two independent docking methods.

2. **MYO1D TH1 domain contact residues** observed at the interface — including SER955, VAL962, VAL964, HIS973, PRO990, and ASN1006 — are located within sheets 8 and 9 of the beta-meander active face, in quantitative agreement with the experimental alanine-scanning data of Ko et al.

3. **Blind docking (Workflow A) is insufficient** for identifying PPI-proximal druggable pockets on EGFR because the ATP-binding cleft dominates the druggability landscape; PPI-guided focused docking (Workflow B) is required to uncover relevant C-lobe surface cavities.

4. **Two tier-1 priority pocket candidates** are identified: PKT07 (perturbation score 0.541), a state-robust rim modulator site sharing LEU1001 with the hotspot patch and displaying favorable tier-1 druggability; and PKT34 (perturbation score 0.492), a cryptic allosteric site located 9.4 Å from the patch centroid that is accessible only in the MD-derived EGFR_170-200 conformer, underscoring the necessity of ensemble-based receptor sampling.

5. **Direct hotspot overlap and druggability are inversely correlated** across the pocket landscape, a finding that reframes prioritization logic for PPI drug discovery: rim and allosteric sites adjacent to the interface offer superior pharmacophore geometry compared to pockets coinciding directly with the protein–protein contact plane.

These results provide a structural and computational foundation for experimental campaigns targeting the EGFR–MYO1D axis as a non-ATP, non-catalytic mechanism of EGFR modulation. PKT07 and PKT34 are proposed as priority targets for fragment screening, SPR binding assays, and structure-guided medicinal chemistry.

---

## Acknowledgments

[To be completed. Include HPC facility acknowledgment, funding sources, and any reagent or software licenses.]

---

## Author Contributions

[To be completed per journal CRediT taxonomy: Conceptualization, Methodology, Software, Formal Analysis, Investigation, Data Curation, Writing — Original Draft, Writing — Review & Editing, Visualization, Supervision, Funding Acquisition.]

---

## Conflicts of Interest

The authors declare no conflicts of interest.

---

## Data Availability

All Python analysis scripts, PBS job submission scripts, and processed CSV result files are available at [repository URL to be provided upon acceptance]. Raw docking coordinate files are available from the corresponding author upon reasonable request.

---

## References

[1] Yarden, Y.; Sliwkowski, M.X. Untangling the ErbB signalling network. *Nat. Rev. Mol. Cell Biol.* **2001**, *2*, 127–137.

[2] Roskoski, R. The ErbB/HER family of protein-tyrosine kinases and cancer. *Pharmacol. Res.* **2014**, *79*, 34–74.

[3] Kobayashi, S.; Boggon, T.J.; Dayaram, T.; Jänne, P.A.; Kocher, O.; Meyerson, M.; Johnson, B.E.; Eck, M.J.; Tenen, D.G.; Halmos, B. EGFR mutation and resistance of non–small-cell lung cancer to gefitinib. *N. Engl. J. Med.* **2005**, *352*, 786–792.

[4] Engelman, J.A.; Zejnullahu, K.; Mitsudomi, T.; Song, Y.; Hyland, C.; Park, J.O.; Lindeman, N.; Gale, C.M.; Zhao, X.; Christensen, J.; et al. MET amplification leads to gefitinib resistance in lung cancer by activating ERBB3 signaling. *Science* **2007**, *316*, 1039–1043.

[5] Thress, K.S.; Paweletz, C.P.; Felip, E.; Cho, B.C.; Stetson, D.; Dougherty, B.; Lai, Z.; Markovets, A.; Vivancos, A.; Kuang, Y.; et al. Acquired EGFR C797S mutation mediates resistance to AZD9291 in non–small cell lung cancer harboring EGFR T790M. *Nat. Med.* **2015**, *21*, 560–562.

[6] Scott, D.E.; Bayly, A.R.; Abell, C.; Skidmore, J. Small molecules, big targets: drug discovery faces the protein–protein interaction challenge. *Nat. Rev. Drug Discov.* **2016**, *15*, 533–550.

[7] Arkin, M.R.; Tang, Y.; Wells, J.A. Small-molecule inhibitors of protein–protein interactions: progressing toward the reality. *Chem. Biol.* **2014**, *21*, 1102–1114.

[8] Hokanson, D.E.; Laakso, J.M.; Lin, T.; Sept, D.; Ostap, E.M. Myo1c binds phosphoinositides through a putative pleckstrin homology domain. *Mol. Biol. Cell* **2006**, *17*, 4856–4865.

[9] Ko, [First Initial].; [Co-authors]. [Title of Ko et al. paper]. *[Journal]* **[Year]**, *[Volume]*, [pages]. *(Key experimental reference: alanine-scanning mutagenesis of MYO1D TH1 domain sheets 8/9 active face validates EGFR binding interface.)*

[10] Nero, T.L.; Morton, C.J.; Holien, J.K.; Wielens, J.; Parker, M.W. Oncogenic protein interfaces: small molecules, big challenges. *Nat. Rev. Cancer* **2014**, *14*, 248–262.

[11] Sperandio, O.; Reynès, C.H.; Camproux, A.-C.; Villoutreix, B.O. Rationalizing the chemical space of protein–protein interaction inhibitors. *Drug Discov. Today* **2010**, *15*, 220–229.

[12] Jumper, J.; Evans, R.; Pritzel, A.; Green, T.; Figurnov, M.; Ronneberger, O.; Tunyasuvunakool, K.; Bates, R.; Žídek, A.; Potapenko, A.; et al. Highly accurate protein structure prediction with AlphaFold. *Nature* **2021**, *596*, 583–589.

[13] Huang, J.; MacKerell, A.D. CHARMM36 all-atom additive protein force field: validation based on comparison to NMR data. *J. Comput. Chem.* **2013**, *34*, 2135–2145.

[14] Chaudhury, S.; Lyskov, S.; Gray, J.J. PyRosetta: a script-based interface for implementing molecular modeling algorithms using Rosetta. *Bioinformatics* **2010**, *26*, 689–691.

[15] Jiménez-García, B.; Roel-Touris, J.; Romero-Durana, M.; Vidal, M.; Jiménez-González, D.; Fernández-Recio, J. LightDock: a new multi-scale approach to protein–protein docking. *Bioinformatics* **2018**, *34*, 49–55.

[16] Le Guilloux, V.; Schmidtke, P.; Tuffery, P. Fpocket: an open source platform for ligand pocket detection. *BMC Bioinformatics* **2009**, *10*, 168.

[17] Trott, O.; Olson, A.J. AutoDock Vina: improving the speed and accuracy of docking with a new scoring function, efficient optimization, and multithreading. *J. Comput. Chem.* **2010**, *31*, 455–461.

[18] Schrödinger, LLC. The PyMOL Molecular Graphics System, Version 2.5. **2021**.

[19] Zarzycka, B.; Kuenemann, M.A.; Miteva, M.A.; Nicolaes, G.A.F.; Vriend, G.; Sperandio, O. Stabilization of protein–protein interaction complexes through small molecules. *Drug Discov. Today* **2016**, *21*, 48–57.

[20] Cimermancic, P.; Weinkam, P.; Rettenmaier, T.J.; Bichmann, L.; Keedy, D.A.; Woldeyes, R.A.; Schneidman-Duhovny, D.; Demerdash, O.N.; Mitchell, J.C.; Wells, J.A.; et al. CryptoSite: expanding the druggable proteome by characterization and prediction of cryptic binding sites. *J. Mol. Biol.* **2016**, *428*, 709–719.
