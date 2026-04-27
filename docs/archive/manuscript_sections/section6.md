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
