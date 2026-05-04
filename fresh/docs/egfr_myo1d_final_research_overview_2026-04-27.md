# EGFR-MYO1D 연구 개요 및 확정된 분석 설계

작성일: 2026-04-27

## 1. 연구의 핵심 질문

이 연구는 EGFR에 MYO1D를 단순히 도킹하는 작업이 아니다. 연구의 핵심 질문은 다음과 같다.

> 세포막에 고정된 비활성 대칭 EGFR 이량체에서 MYO1D TH1 beta-meander가 EGFR kinase domain의 어느 표면을 인식하며, 이 결합을 교란할 수 있는 non-ATP, membrane-compatible, dimer-accessible pocket은 어디인가?

최종 목표는 세 단계로 정리된다.

1. EGFR 비활성 대칭 이량체에서 MYO1D가 인식하는 EGFR-side PPI patch를 정의한다.
2. 그 PPI patch 주변에서 ATP pocket이 아닌 막 근접 lower/lateral pocket을 찾는다.
3. MYO1D activity-associated compound 3종이 해당 pocket에 수렴하는지 분석하여 EGFR-MYO1D PPI 교란 후보 pocket과 compound hypothesis를 만든다.

## 2. 생물학적 배경과 연구 가설

EGFR은 세포막에 존재하는 receptor tyrosine kinase이며, extracellular domain, transmembrane helix, juxtamembrane region, kinase domain, C-terminal tail로 구성된다. EGFR의 활성화는 단순한 ON/OFF 전환이 아니라 transmembrane/JM 구조, kinase-domain dimer geometry, membrane interaction이 함께 관여하는 구조적 전환이다.

이 연구의 출발점은 ligand 결합 전 또는 activation 전후의 membrane-proximal 비활성 상태이다. EGFR은 ligand가 없는 상태에서도 pre-formed dimer 또는 inactive dimer 상태를 취할 수 있으며, 3GT8은 EGFR kinase domain의 비활성 대칭 이량체 구조를 제공한다. 본 연구에서는 이 비활성 대칭 이량체를 MYO1D가 결합할 수 있는 receptor state로 설정한다.

MYO1D는 class I myosin 계열의 actin/membrane-associated molecular motor이다. MYO1D는 EGFR family receptor를 activation 전에 plasma membrane에 유지하는 역할을 하며, EGFR family의 kinase domain과 직접 결합한다. EGFR 쪽 결합 영역은 ATP-binding N-lobe가 아니라 kinase domain C-lobe 쪽으로 정리된다. 따라서 본 연구에서 찾는 pocket은 EGFR ATP pocket이 아니라 MYO1D가 인식하는 C-lobe 주변의 PPI-proximal pocket이다.

연구 가설은 다음과 같다.

> MYO1D TH1 beta-meander는 비활성 대칭 EGFR 이량체의 membrane-proximal C-lobe/lower-lateral surface에 결합하여 EGFR을 plasma membrane-actin axis에 유지하며, 이 결합면 또는 그 주변 pocket에 결합하는 소분자는 EGFR-MYO1D interaction을 교란할 수 있다.

## 3. EGFR receptor model

### 3.1 Receptor state 정의

본 연구에서 사용하는 receptor ensemble은 세 가지 EGFR 구조 상태로 구성된다.

| State | 지위 | 의미 |
|---|---|---|
| EGFR_160-185 | Primary state | 2M0B/2M20/3GT8 기반 TM-JM-kinase membrane model의 stable MD cluster representative |
| EGFR_170-200 | Primary state | 동일한 membrane-validated trajectory에서 선택된 다른 stable MD cluster representative |
| 3GT8_raw | Reference state | 3GT8 crystal-derived inactive symmetric kinase-domain dimer reference |

EGFR_160-185와 EGFR_170-200이 main receptor state이다. 이 두 구조는 단순한 kinase-domain snapshot이 아니라, TM/JM와 kinase domain을 조립한 membrane-embedded EGFR model을 MD로 평가한 뒤 선택한 안정화 상태이다. 3GT8_raw는 crystallographic reference state로 사용한다.

본 연구의 main receptor assumption은 다음과 같다.

- EGFR는 dimer 상태로 분석한다.
- 이 dimer는 비활성 대칭 이량체이다.
- 두 EGFR copy는 계산상 chain A/B로 구분하지만, 해석상 symmetry-equivalent copy로 취급한다.
- 최종 PPI patch와 pocket은 dimer receptor context에서 해석한다.

### 3.2 TM-JM-kinase EGFR model construction

EGFR membrane-compatible receptor model은 다음 구조들을 조합해 만들었다.

| 구성 요소 | 사용 구조 | 역할 |
|---|---|---|
| TM-JM dimer scaffold | 2M0B NMR model 16, residues 634-677 | 비활성 TM/JM 이량체 scaffold |
| JM-A segment | 2M20 JM-A 부분 | 각 chain에 따로 superpose하여 JM geometry 보강 |
| Kinase-domain dimer | 3GT8 inactive symmetric kinase dimer | 비활성 대칭 kinase-domain source |
| Missing linker/gap | MODELLER | JM-B 및 TM/JM-KD 연결부 gap modeling |

2M0B는 NMR ensemble 구조이므로 model 16을 사용하였다. 2M20의 JM-A 부분은 두 chain에 각각 따로 겹쳐 필요한 JM segment만 사용하였다. 3GT8은 inactive symmetric kinase-domain dimer source로 사용하였다.

구조 조립 시 좌표계는 다음처럼 통일하였다.

| 축 | 정의 |
|---|---|
| Z축(+) | C2 대칭축 및 막 법선 방향, extracellular to intracellular |
| X축(+) | EGFR chain A to chain B 방향 |

2M0B TM-JM dimer는 Z축 기준으로 정렬하고, 3GT8 kinase-domain dimer는 X축 기준으로 정렬하였다. 이후 2M0B C-terminal과 3GT8 N-terminal 사이의 Z축 간격을 최대한 가깝게 배치하고, 연결되지 않은 JM-B/gap region은 MODELLER로 모델링하였다.

넘버링은 다음 기준으로 통일하였다.

| Source | Numbering rule |
|---|---|
| 2M0B | PDB residue number = UniProt residue number |
| 3GT8 | PDB residue number + 24 = UniProt residue number |

모델 배향은 Arkhipov 계열의 membrane-proximal EGFR dimer modeling 개념을 참고하여 잡았다. Kinase domain의 양전하 patch, 특히 Lys689, Lys690, Lys692, Lys715가 lower leaflet의 음전하 지질과 접촉할 수 있도록 배향하였다.

### 3.3 Rotation scan과 membrane MD validation

TM/JM 부분과 kinase domain 사이의 최적 배향을 찾기 위해 기준 모델에서 10도 간격으로 회전 변이체를 만들었다.

```text
rotation variants: -20°, -10°, 0°, +10°, +20°
```

각 변이체는 동일하게 MODELLER gap modeling을 거친 뒤 CHARMM-GUI로 membrane system을 만들고 MD를 수행하였다.

최종 MD 조건은 다음과 같다.

| 항목 | 조건 |
|---|---|
| Force field | CHARMM36m |
| Temperature | 298 K |
| Salt | NaCl 0 M |
| Upper leaflet | POPC 100% |
| Lower leaflet | POPC 70%, POPS 30% |
| Simulation system | CHARMM-GUI to GROMACS |

5개 회전 변이체 중 +10° 모델이 membrane-compatible stable model로 선택되었다. 이후 이 stable trajectory에서 EGFR_160-185와 EGFR_170-200 cluster representative를 선택하여 downstream docking receptor state로 사용한다.

### 3.4 Model-derived membrane frame

본 연구에서 사용하는 membrane frame은 TM-JM-kinase EGFR dimer를 공통 좌표계에 정렬하고, membrane MD에서 안정한 배향을 선택한 뒤 그 모델에서 상속한 model-derived membrane frame이다.

Downstream 분석에서 방향성은 다음처럼 해석한다.

| 용어 | 정의 |
|---|---|
| membrane normal | 모델 조립 좌표계의 Z축, extracellular to intracellular 방향 |
| lower side | TM/JM anchor에 가까운 kinase-domain lower surface |
| lateral side | dimer center에서 바깥쪽으로 열린 EGFR copy의 outward-facing surface |
| central interface | 두 EGFR copy 사이에 묻힌 dimer interface region |

최종 pocket 후보는 lower/lateral surface에 있어야 하며, central dimer interface에 깊게 묻힌 pocket이나 ATP cleft가 아니라 MYO1D interaction을 교란할 수 있는 C-lobe 주변 surface pocket이어야 한다.

## 4. MYO1D model and binding region

### 4.1 MYO1D source와 construct

MYO1D partner는 AlphaFold human MYO1D model을 기반으로 한다.

```text
AlphaFold model: AF-O94832-F1-model_v6
Main construct: MYO1D residues 955-1006
Domain context: TH1 C-terminal beta-meander region
```

이 construct는 MYO1D TH1 domain의 C-terminal beta-meander motif를 포함한다. 전체 TH1 domain 대신 EGFR binding에 필요한 beta-meander region을 중심으로 사용한다.

### 4.2 Key beta-sheet region

MYO1D의 EGFR binding 관련 key region은 beta-sheet 8, 9, 12를 포함하는 residue set으로 둔다.

| Region | Residues | 역할 |
|---|---|---|
| beta-sheet 8 | 961-964 | MYO1D beta-meander key region |
| beta-sheet 9 | 968-972 | MYO1D beta-meander key region |
| beta-sheet 12 | 993-997 | MYO1D beta-meander key/support region |

본 분석에서는 beta-sheet 8/9/12를 통합 key region으로 사용한다. Sheet 12가 모든 docking pose에서 EGFR과 직접 접촉해야 한다고 가정하지 않는다. 대신 문헌상 EGFR/RTK binding에 필요한 beta-meander 구조 요소로 포함한다.

### 4.3 MYO1D pose interpretation

MYO1D beta-meander는 얇은 beta-sheet ribbon 구조이므로 단순 contact count만으로 올바른 결합면과 뒤집힌 결합면을 구분하기 어렵다. 따라서 PPI docking 결과는 다음 기준으로 해석한다.

- beta-sheet 8/9/12 key region이 EGFR과 접촉하는지 확인한다.
- MYO1D beta-meander face가 EGFR receptor surface를 향하는지 orientation filter로 확인한다.
- terminal 또는 flexible tail contact가 key region 접촉 없이 단독으로 dominant한 pose는 consensus PPI patch 근거에서 제외한다.
- 단일 best pose가 아니라 여러 seed와 receptor state에서 반복되는 EGFR-side surface patch를 consensus로 본다.

첫 production analysis에서는 key residue에 score bonus를 주지 않는다. Key region은 docking score를 바꾸는 constraint가 아니라, docking 후 PPI pose를 해석하고 QC하는 annotation으로 사용한다.

## 5. Activity-associated compound anchors

본 연구에는 교수님이 제공한 MYO1D activity-associated compound 3종이 있다.

```text
173940
97806
VAX-C12_0
```

이 compound들은 공개된 약물로 다루지 않으며, 논문이나 외부 문서에서는 필요 시 Cpd-A, Cpd-B, Cpd-C처럼 anonymized ID를 사용할 수 있다.

이 compound들의 역할은 단순한 docking probe가 아니라 chemical anchor이다. 즉, MYO1D activity 저하와 관련된 것으로 공유된 compound들이 EGFR dimer의 어느 pocket에 수렴하는지 보고, 그 pocket이 MYO1D PPI patch와 연결되는지를 평가한다.

해석 범위는 다음과 같다.

- 이 compound들은 EGFR-MYO1D PPI disruption hypothesis를 만드는 chemical anchor로 사용한다.
- 이 compound들이 EGFR-side non-ATP C-lobe/lower-lateral pocket에 반복적으로 결합하면, 해당 pocket을 EGFR-MYO1D PPI disruption hypothesis의 chemical anchor pocket으로 우선순위화한다.
- ATP pocket으로 들어가는 pose는 EGFR-MYO1D PPI pocket 후보가 아니라 ATP-site binding pose로 분류한다.

## 6. 최종 분석 workflow

확정된 분석 흐름은 다음 순서로 진행된다.

### Phase 1. Receptor preparation

- EGFR_160-185와 EGFR_170-200을 primary membrane-validated receptor state로 사용한다.
- 3GT8_raw는 crystallographic inactive kinase-domain dimer reference로 사용한다.
- EGFR receptor는 항상 dimer 상태로 유지한다.
- MD-derived PDB에서 두 EGFR copy가 같은 chain ID로 저장되어 있으면, 계산 전에 chain A/B로 명확히 구분한다.
- 두 EGFR copy는 대칭 copy로 해석한다.

### Phase 2. MYO1D PPI docking

- MYO1D 955-1006 beta-meander construct를 EGFR dimer receptor에 docking한다.
- PPI docking은 single best pose가 아니라 seed/state ensemble로 해석한다.
- MYO1D beta-sheet 8/9/12 key region engagement와 orientation을 QC한다.
- Tail-only 또는 face-flipped pose는 consensus patch 계산에서 제외한다.
- EGFR-side contact residues를 모아 state- and symmetry-aware consensus PPI patch를 만든다.

### Phase 3. Dimer/membrane-aware pocket discovery

- Pocket discovery는 EGFR dimer receptor에서 수행한다.
- 후보 pocket은 MYO1D PPI patch와 가까운 C-lobe/lower-lateral surface에 있어야 한다.
- ATP binding pocket은 최종 PPI-disruptive pocket 후보에서 제외한다.
- Central dimer interface에 묻혀 compound 접근성이 낮은 pocket도 제외한다.

### Phase 4. Compound anchor docking

- 173940, 97806, VAX-C12_0을 activity-associated chemical anchors로 사용한다.
- 각 compound가 non-ATP, lower/lateral, PPI-proximal pocket에 수렴하는지 평가한다.
- 복수 compound가 같은 pocket family에 수렴하면 cross-chemical support가 있는 pocket으로 본다.
- Docking affinity 단독으로 최종 후보를 정하지 않고, pocket 위치, PPI proximity, non-ATP status, pose convergence를 함께 평가한다.

### Phase 5. Integrated candidate nomination

최종 후보는 다음 조건을 모두 만족해야 한다.

| 기준 | 내용 |
|---|---|
| Dimer context | EGFR dimer receptor에서 확인된 pocket |
| Membrane context | model-derived membrane frame에서 lower/lateral surface에 위치 |
| PPI relevance | MYO1D consensus PPI patch와 인접하거나 같은 surface network에 있음 |
| Non-ATP | ATP/ANP pocket과 분리됨 |
| Accessibility | central dimer interface에 묻히지 않음 |
| Compound support | activity-associated compound가 해당 pocket에 안정적으로 수렴 |
| Mechanistic interpretation | EGFR-MYO1D PPI disruption hypothesis로 설명 가능 |

## 7. Pocket 해석 기준

본 연구에서 pocket은 단순히 “소분자가 들어간 위치”가 아니라 EGFR-MYO1D PPI를 교란할 수 있는 구조적 위치여야 한다.

### 7.1 ATP exclusion

EGFR ATP-binding pocket은 main target이 아니다. MYO1D는 EGFR kinase domain C-lobe와 결합하고 ATP-binding N-lobe는 이 interaction의 중심이 아니므로, ATP pocket에 결합한 compound는 ATP-site binding pose로 따로 분류한다.

### 7.2 PPI proximity

Pocket은 MYO1D PPI patch와 다음 중 하나의 관계를 가져야 한다.

| Relationship | 의미 |
|---|---|
| Orthosteric | pocket residue가 PPI patch residue와 직접 겹침 |
| Rim | PPI patch 가장자리에 위치하여 MYO1D 접근 경로를 막을 수 있음 |
| Allosteric-near | 직접 겹치지는 않지만 같은 lower/lateral C-lobe surface network에 있음 |

### 7.3 Membrane-compatible lower/lateral geometry

Pocket은 TM/JM anchor에 가까운 kinase-domain lower side 또는 dimer 바깥쪽 lateral side에 있어야 한다. 이 기준은 눈으로 보는 방향성이 아니라, receptor model construction에서 정의한 Z/X coordinate frame과 dimer geometry를 사용해 판정한다.

## 8. 결과 해석 원칙

최종 분석 결과는 다음 범위에서 해석한다.

- MYO1D TH1 beta-meander가 EGFR dimer의 C-lobe/lower-lateral surface를 인식하는 구조적 hypothesis를 세운다.
- 특정 pocket이 non-ATP, membrane-compatible, PPI-proximal pocket으로 계산적으로 지지되는지 평가한다.
- 3개 activity-associated compound가 같은 pocket family에 수렴하는지 확인하여 chemical-anchor support를 부여한다.
- Pocket/compound pair를 EGFR-MYO1D PPI disruption hypothesis의 후보로 제안한다.

## 9. 현재 확정된 연구 자산

| 항목 | 확정 내용 |
|---|---|
| Receptor model type | 비활성 대칭 EGFR TM-JM-kinase dimer |
| Model source structures | 2M0B model 16, 2M20 JM-A, 3GT8 inactive symmetric kinase dimer |
| Primary receptor states | EGFR_160-185, EGFR_170-200 |
| Reference receptor state | 3GT8_raw |
| Membrane frame | 모델 조립 좌표계와 membrane MD stable orientation에서 상속 |
| MYO1D source | AF-O94832-F1-model_v6 |
| MYO1D construct | residues 955-1006 |
| MYO1D key region | beta-sheet 8/9/12, residues 961-964, 968-972, 993-997 |
| Compound anchors | 173940, 97806, VAX-C12_0 |
| Main output | EGFR-MYO1D PPI patch, non-ATP lower/lateral pocket, compound-anchor convergence |

## 10. 최종 산출물

이 연구의 최종 계산 산출물은 다음과 같다.

1. EGFR_160-185, EGFR_170-200, 3GT8_raw receptor state별 dimer-normalized input.
2. TM-JM-kinase model construction provenance와 model-derived membrane frame.
3. MYO1D 955-1006 beta-meander PPI docking ensemble.
4. MYO1D beta-sheet 8/9/12 key region 기반 PPI pose QC table.
5. EGFR-side consensus PPI patch.
6. Non-ATP, lower/lateral, dimer-accessible pocket candidate table.
7. 173940/97806/VAX-C12_0 compound-anchor docking support table.
8. EGFR-MYO1D PPI-disruptive pocket/compound hypothesis shortlist.

## 11. 논문용 요약 문장

본 연구는 2M0B model 16, 2M20 JM-A segment, 3GT8 inactive symmetric kinase-domain dimer를 조립하여 TM-JM-kinase EGFR 비활성 대칭 이량체 모델을 구축하고, 회전 배향 scan과 membrane MD를 통해 안정한 membrane-compatible receptor state를 선택하였다. Stable trajectory에서 얻은 EGFR_160-185와 EGFR_170-200 cluster representative를 primary receptor ensemble로 사용하고, 3GT8_raw는 crystallographic inactive kinase-domain reference로 사용한다. MYO1D는 AF-O94832-F1-model_v6에서 유래한 955-1006 TH1 beta-meander construct로 표현하며, beta-sheet 8/9/12 residue set을 EGFR binding-relevant key region으로 사용한다. Downstream 분석은 EGFR dimer context와 model-derived membrane frame을 유지한 상태에서 MYO1D consensus PPI patch, non-ATP lower/lateral pocket, activity-associated compound anchor convergence를 통합하여 EGFR-MYO1D PPI disruption 후보를 도출한다.

## 12. 주요 근거 구조와 문헌 축

- 2M0B: EGFR TM-JM inactive dimer NMR 구조, model 16 사용.
- 2M20: EGFR TM/JM 구조 중 JM-A segment 보강에 사용.
- 3GT8: EGFR inactive symmetric kinase-domain dimer reference.
- Ko et al.: MYO1D가 EGFR family kinase domain과 결합하여 plasma membrane retention에 관여하며, MYO1D tail domain의 beta-meander motif가 중요하다는 근거.
- MYO1D/EGFR C-lobe 연구: EGFR kinase domain C-lobe 856-979가 MYO1D interaction에 중요하며 ATP-binding region은 중심 결합부위가 아니라는 근거.
- EGFR TM/JM/KD structural literature: TM/JM geometry, membrane coupling, inactive/active dimer transition, kinase-domain dimerization의 구조적 근거.
