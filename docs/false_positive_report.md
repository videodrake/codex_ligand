# False Positive Report — Supporting Evidence

## 1. Project Description

This project is a computational structural biology study investigating the protein–protein interaction (PPI) between EGFR (epidermal growth factor receptor) kinase domain (residues 699–1007) and the TH1 domain of MYO1D (myosin-1D, residues 955–1006). The research goal is to identify the precise binding interface between these two proteins and to discover nearby small-molecule-accessible surface pockets that could modulate their association. The computational pipeline uses PyRosetta for protein–protein docking (600,000 models across three receptor conformations), AutoDock Vina for small-molecule docking, fpocket for geometric pocket detection, LightDock for independent cross-validation, and GROMACS with the CHARMM36 force field for 200 ns molecular dynamics simulations. All input data consists of PDB coordinate files, molecular dynamics trajectory clusters, and CSV-formatted numerical results. The MYO1D structure was obtained from AlphaFold predictions. This work follows standard computational drug-discovery methodology as described in textbooks and published literature.

## 2. Nature of the Code and Data

- Python scripts using NumPy and PCA (principal component analysis) to compute geometric plane normals of beta-sheet surfaces and dot products for orientation filtering of docking poses
- CSV parsing and aggregation scripts to calculate per-residue occupancy statistics, cluster consensus scores, and cross-state comparisons across three EGFR conformational states
- A 4-axis weighted scoring framework (Python, pure arithmetic) that ranks candidate surface pockets by PPI proximity, druggability, perturbation relevance, and conformational robustness
- Shell scripts (Bash) for submitting PBS/HPC batch jobs and collecting output files from a university computing cluster
- Markdown report generators that produce human-readable summaries of docking statistics and residue-level interface analyses
- No compiled binaries, no external network calls, no system-level operations — all code operates on local CSV and PDB text files

## 3. Scientific Terminology That Likely Triggered the Classifier

| Term | Meaning in This Context |
|------|------------------------|
| **blind docking** | A docking protocol where the search region covers the entire protein surface without prior knowledge of the binding site — standard computational chemistry terminology (Morris et al., 2009) |
| **druggable pocket** | A surface cavity on a protein with geometric and physicochemical properties suitable for small-molecule binding — standard medicinal chemistry term (Hajduk et al., 2005) |
| **perturbation scoring** | A numerical framework assessing whether a small molecule in a given pocket could disrupt a protein–protein interface — refers to thermodynamic perturbation, not system disruption |
| **active face** | The side of a beta-sheet structure experimentally shown to mediate protein–protein contact — a structural biology descriptor |
| **orientation filter** | A geometric test (dot-product of surface normal vectors) to remove physically implausible docking poses — a quality-control step in molecular modeling |
| **seed** | A random-number-generator initialization value used to produce independent docking replicates for statistical robustness |
| **hotspot residue** | An amino acid at a protein–protein interface that contributes disproportionately to binding energy — coined by Clackson & Wells (1995) in *Science* |
| **kill / pass / fail** | Classification labels for docking poses that do or do not satisfy geometric criteria — analogous to pass/fail in any quality-control pipeline |
| **exhaustiveness** | A parameter in AutoDock Vina controlling the thoroughness of conformational sampling (Trott & Olson, 2010) |
| **PBS script / qsub** | Portable Batch System job submission commands for university HPC clusters — standard research computing infrastructure |

## 4. What Is Absent (Negative Evidence)

- **No network-related code**: no sockets, no HTTP requests, no urllib, no API calls to external services, no port scanning, no DNS operations
- **No file-system manipulation beyond scientific I/O**: no path traversal, no permission changes, no access to system directories outside the project workspace
- **No credential handling**: no passwords, tokens, API keys, or authentication logic anywhere in the codebase
- **No references to CVEs, security advisories, known software flaws, or defensive/offensive security tools**
- **No reverse engineering**: no disassemblers, no binary analysis, no decompilation, no memory inspection
- **No obfuscation or evasion techniques**: all code is plain Python with descriptive variable names and docstrings
- **No command-and-control infrastructure**: no remote shells, no encrypted channels, no callback mechanisms
- **The only external software invoked** consists of established open-source scientific tools: PyRosetta (academic license), AutoDock Vina (Scripps Research), fpocket (University of Paris), GROMACS (KTH Royal Institute of Technology), and LightDock (Barcelona Supercomputing Center)

## 5. Comparison to Legitimate Scientific Workflows

This pipeline directly mirrors standard computational drug-discovery workflows as taught in graduate-level courses and described in peer-reviewed publications:

- **Protein–protein docking → interface analysis → pocket identification → focused small-molecule docking** is the canonical PPI-disruption pipeline described in Nero et al., "Oncogenic protein interfaces: small molecules, big challenges" (*Nature Reviews Cancer*, 2014).
- The orientation-filtering approach (PCA-based surface normal computation + dot-product classification) follows established geometric filtering methods in ClusPro and HADDOCK docking servers.
- The use of three receptor conformational states (one crystal structure + two MD cluster representatives) for cross-state validation is standard ensemble docking practice (Amaro et al., *J. Med. Chem.*, 2018).
- fpocket for cavity detection, AutoDock Vina for scoring, and PyRosetta for protein–protein energetics are among the most widely cited tools in structural bioinformatics, collectively appearing in over 20,000 published studies.

---

**In summary, this session contains a protein–protein interaction analysis and structure-based drug-pocket discovery pipeline using standard open-source computational chemistry tools, consistent with routine academic pharmaceutical research, and contains no content within the scope of Anthropic's cyber safeguards.**
