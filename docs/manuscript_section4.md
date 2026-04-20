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
