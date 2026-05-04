---
project: EGFR-MYO1D membrane-compatible PPI/pocket/compound discovery
version_date: 2026-04-27
status: FINAL CLEAN PROJECT KNOWLEDGE
language: Korean with English technical terms
recommended_upload: Upload this single Markdown file as the project knowledge source.
---

# EGFR-MYO1D Project Knowledge

## 1. Project definition

본 프로젝트의 목표는 **EGFR와 MYO1D의 상호작용을 교란할 수 있는 약물 후보를 찾는 것**이다. 이를 위해 먼저 세포막 환경과 양립하는 EGFR inactive dimer receptor model을 확보하고, 그 위에서 MYO1D TH1 beta-meander가 인식하는 EGFR-side PPI surface를 정의한 뒤, 해당 표면 주변의 **non-ATP, dimer-accessible, membrane-compatible pocket**에 결합할 수 있는 compound를 선별한다.

프로젝트의 핵심 흐름은 다음과 같다.

```text
Membrane-compatible EGFR inactive dimer model
    -> MYO1D TH1 beta-meander PPI consensus patch
    -> non-ATP / C-lobe / PPI-adjacent druggable pocket
    -> focused compound docking
    -> EGFR-MYO1D PPI-disruptive compound candidates
```

한 문장으로 정리하면 다음과 같다.

> This project aims to identify membrane-compatible, non-ATP EGFR dimer-surface pockets that can perturb the interaction between the MYO1D TH1 beta-meander and the EGFR kinase-domain dimer.

---

## 2. Scientific hypothesis

MYO1D는 EGFR family receptor를 plasma membrane에 유지하는 molecular motor로 작용하며, MYO1D tail domain의 TH1 beta-meander motif가 EGFR family kinase domain과 직접 결합한다. NSCLC 맥락에서는 MYO1D가 EGFR kinase domain의 **C-lobe**, 즉 ATP-binding region과 공간적으로 분리된 부위와 결합하는 것으로 보고되어 있다.

따라서 본 연구의 working hypothesis는 다음과 같다.

```text
MYO1D TH1 beta-meander binds a C-lobe/non-ATP surface of the EGFR kinase-domain dimer.
A compound occupying a nearby EGFR dimer-surface pocket may perturb this PPI
and weaken MYO1D-dependent EGFR membrane retention.
```

이 가설은 기존 ATP-site TKI와 다른 접근을 가진다. 본 프로젝트에서 compound target은 EGFR ATP-binding cleft가 아니라, EGFR-MYO1D PPI와 연결된 **C-lobe / allosteric / dimer-surface pocket**이다.

---

## 3. Literature basis

### 3.1 EGFR dimer regulation

EGFR kinase activation은 activation-loop phosphorylation보다 **kinase-domain asymmetric dimer formation**과 강하게 연결된다. 한쪽 kinase domain의 C-lobe가 activator 역할을 하고, 다른쪽 kinase domain의 N-lobe가 receiver로 활성화되는 allosteric mechanism이 핵심이다.

EGFR intracellular regulation에는 juxtamembrane segment도 포함된다. JM-B는 activator kinase의 C-lobe를 감싸는 JM latch를 형성하며, JM-A와 transmembrane helix는 kinase-domain dimer의 membrane-proximal geometry를 조절한다. Inactive symmetric dimer에서는 proximal C-terminal tail과 AP-2 helix가 JM latch 형성을 차단하는 구조적 요소로 작용할 수 있다.

### 3.2 Membrane-proximal context

EGFR는 single-pass membrane receptor이므로 kinase domain의 배향은 membrane geometry와 함께 해석되어야 한다. TM, JM-A, JM-B, kinase domain, proximal C-terminal tail의 상대 배치는 EGFR inactive/active state를 해석하는 핵심 구조 정보다.

### 3.3 MYO1D-EGFR interaction

MYO1D는 unphosphorylated EGFR family receptor를 plasma membrane에 anchor하는 역할을 한다. MYO1D TH1 domain의 C-terminal beta-meander motif는 EGFR family kinase domain binding에 중요하며, beta-sheet 8, 9, 12가 RTK binding에 필요한 구조 요소로 제시되어 있다.

NSCLC 연구에서는 MYO1D가 EGFR kinase domain의 C-lobe와 상호작용하며, 이 부위는 ATP-binding region을 포함하지 않는 것으로 보고되어 있다. 이 점이 본 프로젝트의 non-ATP PPI-disruptive strategy의 핵심 근거다.

### 3.4 Key references for manuscript citation

| Topic | Reference role |
|---|---|
| EGFR allosteric asymmetric kinase-domain activation | Zhang et al., Cell, 2006 |
| JM latch and inactive symmetric kinase-domain dimer | Jura et al., Cell, 2009 |
| EGFR TM/JM/membrane coupling | Endres et al., Cell, 2013 |
| EGFR structural regulation overview | Kovacs et al., Annual Review of Biochemistry, 2015 |
| EGFR C-terminal tail/AP-2 helix regulation | Gajiwala, Protein Science, 2013 |
| MYO1D anchoring EGFR family before activation | Ko et al., Oncogene, 2019 |
| MYO1D interaction with EGFR C-lobe/non-ATP region in NSCLC | Ko et al., Clinical and Translational Medicine, 2021 |

---

## 4. Confirmed EGFR receptor model

### 4.1 Model identity

현재 프로젝트의 receptor model은 **EGFR symmetric inactive TM-JM-KD dimer**이다. 이 모델은 하나의 실험 구조를 그대로 사용한 것이 아니라, EGFR의 각 구조 요소를 제공하는 3개 구조를 통합하여 만든 composite model이다.

```text
EGFR symmetric inactive TM-JM-KD dimer
= 2M0B TM dimer
+ 2M20 JM-A reference
+ 3GT8 inactive symmetric kinase-domain dimer
+ MODELLER gap filling
+ rotational orientation scan
+ membrane MD validation
```

### 4.2 Structural templates

| Template | Role in model |
|---|---|
| 2M0B | EGFR transmembrane dimer template |
| 2M20 | JM-A orientation reference, including membrane-proximal LRRLL region |
| 3GT8 | EGFR inactive symmetric kinase-domain dimer template |

### 4.3 Model construction

The receptor model was built by aligning the TM/JM and kinase-domain dimer modules into a common dimer coordinate system.

Core construction steps:

1. The TM dimer was positioned according to the membrane-normal axis.
2. JM-A geometry was incorporated by chain-wise alignment to the TM/JM reference.
3. The 3GT8 inactive symmetric kinase-domain dimer was aligned with its dimer symmetry axis in the same coordinate frame.
4. The kinase-domain dimer was positioned relative to the TM-JM module to allow continuous TM-JM-KD connection.
5. Missing JM-B and kinase-domain loop regions were modeled with MODELLER.
6. The working model was converted to a WT-like EGFR sequence state where required for the receptor construct.

Numbering convention used in the project:

| Region/template | Working convention |
|---|---|
| 2M0B | PDB numbering corresponds to UniProt/mature EGFR numbering used for the TM region |
| 3GT8 | PDB numbering is tracked with project-specific offset mapping; mapping file should be stored with the receptor input |
| Dimer protomers | Chain/protomer identity must be retained in all downstream outputs |

**Data required:** store the final receptor PDB, chain mapping table, residue numbering table, and template-to-model alignment record together with the project input files.

---

## 5. EGFR rotational scan and MD validation

### 5.1 Orientation scan

To identify a membrane-compatible TM-KD orientation, five rotational variants of the EGFR TM-JM-KD inactive dimer were generated:

```text
-20°, -10°, 0°, +10°, +20°
```

A TM-KD dihedral angle was defined using four geometric reference points:

| Point | Definition |
|---|---|
| P1 | Chain A TM helix center of mass |
| P2 | Whole dimer midpoint |
| P3 | Kinase-domain dimer midpoint |
| P4 | Chain A kinase-domain N-lobe center of mass |

This dihedral quantifies how the kinase-domain dimer is rotated relative to the TM helix when viewed along the dimer symmetry axis.

### 5.2 MD system

Each rotational variant was embedded in a lipid bilayer and simulated by GROMACS using CHARMM36m.

| Parameter | Setting |
|---|---|
| System | EGFR symmetric inactive TM-JM-KD dimer |
| Force field | CHARMM36m |
| Membrane upper leaflet | POPC 100% |
| Membrane lower leaflet | POPC 70% + POPS 30% |
| Salt | NaCl 0 M |
| Temperature | 298 K |
| Production MD | 200 ns validation run |
| Model variants | -20°, -10°, 0°, +10°, +20° |

### 5.3 Selected receptor orientation

The **+10° model** is the selected membrane-compatible EGFR inactive dimer receptor model.

Key MD comparison:

| Metric | +10° selected model | Reference collapsed orientation (-10°) |
|---|---:|---:|
| TM-KD dihedral | -106.6° | -132.3° |
| RMSD, last 50 ns | 2.8 ± 0.2 Å | 5.0 ± 0.2 Å |
| KD-membrane Z displacement | +0.2 Å | -4.1 Å |
| Inter-chain H-bonds | 13.1 ± 2.9 | 7.0 ± 2.2 |
| Radius of gyration | 28.4 ± 0.2 Å | 28.9 ± 0.1 Å |
| Positive patch Z, chain B | 27.4 Å | 15.4 Å |

Interpretation:

- The +10° orientation preserved stable kinase-domain membrane association over 200 ns.
- The +10° orientation maintained stronger dimer contacts and a compact receptor arrangement.
- The N-lobe positive patch in the +10° orientation was positioned consistently with interaction toward the negatively charged cytoplasmic leaflet.

The +10° receptor model is the current structural basis for downstream EGFR-MYO1D PPI modeling, pocket accessibility interpretation, and compound docking.

**Data required:** final +10° PDB, starting model PDB, final MD frame, trajectory, topology, MD parameter files, and per-variant MD metric table.

---

## 6. Membrane coordinate frame for downstream analysis

The selected +10° EGFR model provides a project-specific membrane coordinate frame. This frame is used to classify PPI patches and pockets by their spatial relationship to the membrane and dimer axis.

Required geometric definitions:

| Term | Project definition |
|---|---|
| Membrane normal | Axis perpendicular to the lipid bilayer in the selected +10° model |
| Cytosolic side | Side where EGFR kinase domains and MYO1D are located |
| Lower/membrane-proximal surface | EGFR surface closer to TM/JM and cytosolic membrane plane |
| Lateral surface | Outward-facing side of the dimer, away from the central protomer-protomer interface |
| Dimer axis | Axis connecting the two EGFR protomer centroids |
| PPI patch centroid | Center of EGFR residues contacted by accepted MYO1D docking poses |
| Pocket centroid | Center of atoms/residues defining a candidate pocket |

Expected geometry output file:

```text
output/geometry/membrane_frame.json
```

Minimum fields:

```json
{
  "receptor_id": "EGFR_plus10_membrane_validated_dimer",
  "membrane_normal": [0.0, 0.0, 1.0],
  "membrane_plane_point": [0.0, 0.0, 0.0],
  "dimer_axis": [1.0, 0.0, 0.0],
  "protomer_a_centroid": [0.0, 0.0, 0.0],
  "protomer_b_centroid": [0.0, 0.0, 0.0],
  "coordinate_convention": "project-specific values must be computed from the selected +10 receptor model"
}
```

**Data required:** computed numeric membrane frame values for the final selected receptor PDB.

---

## 7. MYO1D input definition

### 7.1 Biological target region

The MYO1D region of interest is the C-terminal TH1 beta-meander motif. This region mediates binding to EGFR family kinase domains and contains the beta-sheet elements required for RTK binding.

Core MYO1D binding elements:

| Element | Residue range | Role |
|---|---:|---|
| N-terminal buffer | 955-960 | Structural buffer for beta-meander model |
| Beta-sheet 8 active region | 961-964 | Primary MYO1D-side PPI evidence |
| Beta-sheet 9 active region | 968-972 | Primary MYO1D-side PPI evidence |
| Beta-sheet 12 support region | 993-997 | Structural/support PPI evidence |
| C-terminal cap/noise-monitoring region | 998+ | Structural cap; contact contribution tracked separately |

### 7.2 Working constructs

| Construct | Residues | Use |
|---|---:|---|
| `MYO1D_sheet8_9_12_core_955_1001` | 955-1001 | Primary working construct |
| `MYO1D_ext_beta_meander_955_1006_tail_masked` | 955-1006 | Comparator when extended beta-meander tail is needed structurally |

Residue annotation for docking outputs:

```ini
[MYO1D_BindingEvidence]
primary_active_face = 961-964,968-972
support_region = 993-997
structural_buffer = 955-960
contact_monitoring_cap = 998-1006
key_residue_contact_is_annotation_only = true
```

The key MYO1D residues are used for QC and contact annotation. They are not used as an artificial score bonus in the primary production interpretation.

**Data required:** final MYO1D PDB source, residue numbering validation, secondary-structure annotation for beta-sheet 8/9/12, and prepared docking partner PDB.

---

## 8. EGFR-MYO1D PPI modeling objective

The PPI modeling goal is to define an **EGFR-side consensus patch** contacted by MYO1D TH1 beta-meander in the membrane-compatible dimer context.

Primary PPI output:

```text
output/ppi/ppi_consensus_patch.csv
```

Expected PPI evidence fields:

| Field | Meaning |
|---|---|
| receptor_id | EGFR receptor model or MD representative used |
| seed / replicate | independent docking replicate identifier |
| pose_id | accepted pose identifier |
| egfr_contact_residues | EGFR residues contacted by MYO1D |
| egfr_contact_centroid | 3D centroid of EGFR-side contact residues |
| protomer_id | EGFR protomer contacted by MYO1D |
| myo1d_contact_residues | MYO1D residues contacting EGFR |
| myo1d_active_face_score | sheet 8/9 active-face contact evidence |
| sheet12_support_score | sheet 12 support contact evidence |
| membrane_side_class | lower/lateral/cytosolic-side geometry class |
| ppi_patch_id | consensus patch assignment |

Accepted PPI interpretation is based on repeated EGFR-side patch occurrence across docking replicates and receptor states.

---

## 9. Pocket discovery objective

Pocket discovery is performed after defining the EGFR-MYO1D PPI patch. The target pocket should be spatially and mechanistically connected to the MYO1D-binding surface.

Pocket selection features:

| Feature | Project meaning |
|---|---|
| PPI-adjacent | Pocket overlaps or lies near the EGFR-MYO1D consensus patch/rim |
| Non-ATP | Pocket is outside the canonical ATP-binding cleft |
| Dimer-accessible | Pocket remains accessible in the EGFR dimer assembly |
| Membrane-compatible | Pocket can be approached from the cytosolic/lower/lateral side |
| Druggable | Pocket has suitable size, shape, residue environment, and docking pose convergence |

Expected pocket output:

```text
output/pockets/egfr_myo1d_ppi_adjacent_pockets.csv
```

Expected fields:

| Field | Meaning |
|---|---|
| pocket_id | candidate pocket identifier |
| receptor_id | receptor model/state |
| pocket_centroid | 3D pocket centroid |
| pocket_residues | EGFR residues defining the pocket |
| distance_to_ppi_patch | centroid or residue-distance relationship to PPI patch |
| non_atp_class | ATP-site separation class |
| membrane_accessibility_class | lower/lateral/cytosolic accessibility class |
| dimer_accessibility_class | accessibility in dimer context |
| druggability_score | pocket quality score from fpocket or equivalent method |
| final_priority | prioritized pocket tier |

**Data required:** fpocket or equivalent pocket-detection outputs generated on the selected EGFR dimer receptor states.

---

## 10. Compound docking objective

Compound docking is used to nominate candidate molecules that can occupy the PPI-adjacent EGFR pocket and plausibly perturb MYO1D binding.

Current ligand set carried by the project:

```text
173940
97806
VAX-C12_0
```

Compound interpretation fields:

| Field | Meaning |
|---|---|
| compound_id | compound identifier |
| pocket_id | target pocket |
| receptor_id | receptor state/model |
| vina_affinity | docking affinity |
| pose_cluster_id | docking pose cluster |
| pose_convergence | repeated pose convergence within pocket |
| pocket_contact_residues | residues contacted by compound |
| ppi_overlap_score | overlap or proximity to MYO1D PPI patch |
| membrane_accessibility | compound approach compatibility |
| basic_chemical_filter | basic developability/chemistry screen |
| final_candidate_tier | final candidate class |

Expected compound output:

```text
output/compound/final_candidate_table.csv
```

**Data required:** ligand SDF/PDBQT files, ligand provenance, protonation/tautomer preparation records, docking boxes, Vina parameters, and pose-clustering outputs.

---

## 11. Data package required for a complete run

The following files should exist before a full paper-grade production run.

| Category | Required file/data |
|---|---|
| EGFR receptor | final +10° EGFR TM-JM-KD dimer PDB |
| EGFR model metadata | template alignment, residue numbering, chain/protomer mapping |
| EGFR MD | trajectory, topology, final frame, MD parameters, per-variant metric table |
| Geometry | `membrane_frame.json`, dimer-axis and membrane-plane definitions |
| MYO1D partner | prepared MYO1D 955-1001 or 955-1006-tail-masked PDB |
| MYO1D annotation | sheet 8/9/12 residue annotation and numbering validation |
| PPI docking | docking configs, accepted pose table, EGFR-side consensus patch table |
| Pocket detection | pocket tables for selected receptor states |
| Compound docking | ligand files, docking parameters, pose tables, compound shortlist |
| Manuscript figures | receptor model figure, MD validation figure, PPI patch figure, pocket/compound figure |

---

## 12. Manuscript-ready project summary

본 연구는 MYO1D가 EGFR family kinase domain의 C-lobe/non-ATP region과 상호작용하여 EGFR membrane retention에 기여한다는 문헌 근거를 바탕으로, EGFR-MYO1D PPI를 교란할 수 있는 compound 후보를 탐색한다. EGFR의 intracellular regulation은 kinase-domain dimerization, JM latch, proximal C-terminal tail, and membrane-proximal geometry와 연결되므로, receptor input은 membrane-compatible EGFR inactive dimer model로 정의한다. 이를 위해 2M0B, 2M20, 3GT8을 통합해 EGFR symmetric inactive TM-JM-KD dimer를 구축하고, -20°, -10°, 0°, +10°, +20° rotational variants를 membrane MD로 평가하였다. 그 결과 +10° orientation이 200 ns 동안 안정적인 membrane-compatible kinase-domain positioning과 dimer integrity를 유지하여 downstream EGFR-MYO1D PPI and pocket analysis의 receptor model로 선정되었다. 이후 MYO1D TH1 beta-meander의 sheet 8/9/12 region을 중심으로 EGFR-side PPI consensus patch를 정의하고, 해당 patch 주변의 non-ATP, dimer-accessible, membrane-compatible pocket을 대상으로 compound docking을 수행한다.

---

## 13. Short handoff summary

This project seeks EGFR-MYO1D PPI-disruptive compounds. Use the MD-validated +10° EGFR symmetric inactive TM-JM-KD dimer model built from 2M0B, 2M20, and 3GT8 as the receptor context. Model MYO1D using the TH1 beta-meander region containing beta-sheets 8, 9, and 12. The primary computational output is an EGFR-side MYO1D PPI consensus patch, followed by discovery of a nearby non-ATP, membrane-compatible, dimer-accessible pocket and focused docking of candidate compounds into that pocket.
