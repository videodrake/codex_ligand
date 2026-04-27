# Computational Identification of Druggable Pockets at the EGFR–MYO1D Protein–Protein Interface: A Multi-Conformational Ensemble Docking and Perturbation Scoring Approach

---

**Authors:** [Author 1, Affiliation]¹; [Author 2, Affiliation]¹; [Corresponding Author, Affiliation]¹*

¹ [Department, Institution, City, Country]

\* Correspondence: [email@institution.edu]

---

## Abstract

Epidermal growth factor receptor (EGFR) signaling is dysregulated in numerous cancers, yet resistance to ATP-competitive kinase inhibitors remains a persistent clinical challenge. The non-catalytic protein–protein interaction (PPI) between the EGFR kinase domain and the tail homology 1 (TH1) domain of Myosin-1D (MYO1D) represents an alternative regulatory axis that may be exploitable for therapeutic intervention. Here we report a comprehensive computational pipeline to map the EGFR–MYO1D binding interface and discover druggable pockets capable of disrupting this interaction. Three EGFR conformational states were sampled — a crystal structure (PDB: 3GT8) and two molecular dynamics (MD) cluster representatives extracted from a 200 ns GROMACS/CHARMM36 simulation — providing an ensemble that captures receptor flexibility. Global protein–protein docking with PyRosetta (600,000 models across 30 independent runs) was cross-validated by LightDock (DFIRE2 scoring), and an orientation filter based on principal component analysis (PCA) of the MYO1D beta-meander sheet plane was applied to retain only physically plausible docking poses. Twelve EGFR C-lobe residues (ILE941, THR940, GLN935, PRO934, PRO937, SER957, VAL980, GLN982, ARG986, HIS988, ARG977, THR993) emerged as robust hotspots present across all three receptor states and confirmed by both docking methods. Cavity detection with fpocket (165 raw → 103 merged pockets) followed by a dual-workflow AutoDock Vina docking strategy identified two tier-1 pockets as priority modulator candidates: PKT07 (perturbation score 0.541, interface rim, always accessible) and PKT34 (perturbation score 0.492, allosteric, state-shifted). ATP-competitive pockets were excluded on experimental grounds. These results provide a structural blueprint for developing non-ATP small-molecule PPI modulators of the EGFR–MYO1D axis.

**Keywords:** EGFR; MYO1D; protein–protein interaction; druggable pocket; molecular docking; PyRosetta; AutoDock Vina; ensemble docking; allosteric modulation; computational drug discovery
