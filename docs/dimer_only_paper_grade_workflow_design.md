# EGFR-MYO1D dimer-only 논문급 재실행 설계

Related document:

- `docs/membrane_aware_ppi_pocket_compound_workflow.md` defines the downstream membrane-aware PPI-to-pocket-to-compound workflow. Use that document for pocket gating, dimer-side/lower-side geometry, and compound nomination criteria.

## 1. 결론 먼저

논문에 실을 목적이라면 이번 재실행은 **dimer vs monomer 비교 실험이 아니어야 한다**. 사용자의 조건처럼 EGFR dimer가 생물학적으로 필수라면, monomer는 주 결과에 들어가면 안 된다. 주 실험의 고정 조건은 다음이어야 한다.

- EGFR receptor: 항상 dimer.
- MYO1D partner: 주 분석에서는 하나의 사전 정의된 construct만 사용.
- Receptor state: 3개 conformational state.
- Stochastic replicate: state당 10 seeds 또는 최소 5 seeds.
- 해석 단위: best model이 아니라 seed/state를 넘어 반복되는 receptor-side patch.

따라서 핵심 질문은 이렇게 바뀐다.

> EGFR kinase-domain dimer 표면에서 MYO1D TH1 beta-meander가 반복적으로 인식하는 receptor-side patch는 어디이며, 그 patch 근처에 druggable pocket이 존재하는가?

## 2. 과학적 근거

### EGFR는 dimer로 고정해야 한다

EGFR kinase domain은 단량체 표면만으로 해석하면 생물학적 상태를 잃는다. EGFR kinase activation은 두 kinase domain 사이의 allosteric/asymmetric dimer interaction과 직접 연결되어 있다. 2GS2/2GS6 계열 구조와 Zhang et al. 2006은 EGFR kinase activation이 asymmetric dimer 형성으로 설명됨을 보였다. 3GT8은 EGFR kinase domain의 inactive dimer/juxtamembrane 관련 구조 근거로 쓰이며, RCSB는 3GT8을 inactive EGFR kinase domain + AMP-PNP 구조로 기록한다.

논문 관점에서는 다음 논리가 가장 방어 가능하다.

1. MYO1D는 ligand 결합 전 plasma membrane에서 EGFR family를 붙잡는 축으로 보고되었다.
2. 그 생물학적 맥락의 EGFR은 receptor dimerization/kinase-domain dimer geometry와 분리하기 어렵다.
3. 따라서 docking receptor는 monomer가 아니라 dimer assembly를 유지해야 한다.
4. monomer 결과는 exploratory 또는 negative/technical control로만 둘 수 있고, main conclusion의 evidence로 쓰지 않는다.

### MYO1D는 beta-meander fragment로 설계하되, 너무 짧게 자르면 안 된다

Ko et al. 2019는 MYO1D C-terminal TH1 domain의 beta-meander motif가 EGFR family kinase domain binding에 중요하다고 보고했다. 2021 correction page의 figure title도 beta10/11을 제외한 C-terminal beta-sheets가 RTK-binding site에 중요하다는 방향을 뒷받침한다.

계산 설계에서는 전체 TH1 domain을 그대로 쓰는 것보다 beta-meander를 쓰는 편이 낫다. 이유는 전체 TH1은 너무 커서 비특이적 표면 접촉을 많이 만들고, docking search가 "MYO1D가 어디로든 붙는 문제"가 되기 쉽기 때문이다. 반대로 962-1006처럼 너무 짧게 자르면 VAL962가 artificial N-terminus가 되어 N-terminal charge/backbone freedom artifact를 만들 수 있다. 그래서 주 construct는 `955-1006`이 가장 합리적이다.

## 3. MYO1D construct 설계

### 주 construct

주 분석에는 하나만 사용한다.

- Name: `MYO1D_ext_beta_meander_955_1006`
- Source: AlphaFold human MYO1D model, UniProt O94832, or curated lab PDB.
- Residues: 955-1006.
- Chain: B.
- Required checks:
  - first residue is 955, not 962.
  - residues 961-964 and 968-972 are present.
  - residues 993-997 are present.
  - no missing backbone atoms in 955-1006.
  - no artificial residue renumbering; use author/UniProt-compatible numbering.

### Tail-artifact control for beta-sheet 8/9/12

MYO1D beta sheets 8, 9, and 12 must be included, but the previous long C-terminal tail can create non-specific terminal contacts. The design should therefore separate structural inclusion from binding evidence.

Recommended candidate main construct:

- Name: `MYO1D_sheet8_9_12_core_955_1001`
- Residues: `955-1001`
- Rationale:
  - `955-960` remains as an N-terminal buffer, so `962` is not an artificial terminus.
  - `961-964,968-972` preserve the sheet 8/9 active-face region.
  - `993-997` preserve sheet 12.
  - `998-1001` acts only as a short C-terminal cap/buffer.
  - `1002-1006` is removed from the candidate main construct because it is the likely long-tail artifact zone.

Conservative comparator:

- Name: `MYO1D_ext_beta_meander_955_1006_tail_masked`
- Residues: `955-1006`
- Use only if the `955-1001` cut destabilizes sheet 12 or changes local beta-strand geometry.
- In this comparator, residues `998-1006` must be treated as a non-binding/noise zone.

Recommended pipeline annotation:

```ini
[Constraints]
key_residues_B = 961-964,968-972,993-997
key_residue_bonus_weight = 0.0

[ExperimentalData]
critical_residues_B = 961-964,968-972,993-997
non_binding_residues_B = 998-1006
```

`key_residue_bonus_weight = 0.0` is intentional for the first rerun. It records key-residue contact ratios without changing pose ranking. A positive bonus weight would change model selection and therefore needs explicit approval before any paper-grade production run.

Decision pilot before full production:

1. Prepare both `955-1001` and `955-1006` from the same MYO1D source structure.
2. Confirm that sheet 12 (`993-997`) remains structurally intact after preparation/relaxation in `955-1001`.
3. Run the same small dimer-only pilot for both constructs on the same EGFR dimer state and seed budget.
4. Compare active-face enrichment, sheet-12 support contacts, orientation-pass rate, and tail-dominant false positives.
5. Use `955-1001` for production if it preserves sheet 8/9/12 behavior while reducing terminal-tail contacts. Otherwise use `955-1006` with strict tail masking.

Post-hoc pose classes:

| Class | Definition | Use |
|---|---|---|
| active_8_9_supported_12 | sheet 8/9 contact, sheet 12/support contact, and orientation pass | accepted |
| active_8_9_only | sheet 8/9 contact and orientation pass, but no sheet 12 support | accepted but weaker |
| sheet12_dominant_review | sheet 12 contact dominates without sheet 8/9 support | review separately |
| tail_dominant_artifact | `998+` contacts dominate or terminal contacts occur without sheet 8/9 engagement | reject/quarantine |
| flipped_or_back_face | contacts exist but orientation filter fails | reject |

### Face model

MYO1D beta-meander는 얇은 beta-sheet ribbon이므로 단순 contact count로는 앞면/뒷면을 구분하지 못한다.

주 해석 규칙:

- Primary active-face evidence: sheet 8/9, residues 961-964 and 968-972.
- Structural/support face: sheet 12, residues 993-997.
- Ambiguous/neutral sheets: sheet 10/11은 primary binding evidence로 쓰지 않는다.

중요한 조정:

- sheet 12를 "무관"으로 버리면 안 된다. Ko correction의 표현상 beta12도 기능적으로 중요할 가능성이 있으므로, 논문에서는 sheet 12를 "support/secondary face"로 표현한다.
- 하지만 orientation filter의 primary pass/fail 기준은 sheet 8/9 active-face normal로 둔다. sheet12 contact가 강한 pose는 별도 class로 flag한다.

권장 MYO1D pose classes:

| Class | 정의 | 해석 |
|---|---|---|
| active-face | sheet 8/9 normal이 receptor를 향하고 sheet 8/9 contact가 충분함 | primary accepted |
| support-assisted | active-face 조건을 통과하면서 sheet12도 접촉 | accepted, strong mechanistic interest |
| sheet12-dominant | sheet8/9 약하고 sheet12만 강함 | review-required |
| back-face/flipped | active-face normal이 receptor 반대 | reject |
| ambiguous | dot product near zero | exclude from consensus, report separately |

## 4. EGFR dimer 설계

### 주 receptor unit

각 state마다 docking receptor는 다음 형식으로 준비한다.

- EGFR dimer chains A/B를 하나의 receptor unit으로 유지.
- PyRosetta docking 편의를 위해 EGFR dimer를 chain A 하나로 병합 가능.
- original EGFR chain A residues: 699-1007.
- original EGFR chain B residues: 699-1007 또는 구조상 존재 범위.
- 병합 시 original chain B residue number는 +1000 offset.
- MYO1D는 chain B로 배치.

즉 docking PDB의 chain 구성은 다음처럼 된다.

| Docking chain | 의미 |
|---|---|
| A:699-1007 | EGFR protomer 1 |
| A:1699-2007 | EGFR protomer 2, +1000 offset |
| B:955-1006 | MYO1D beta-meander |

결과 해석 시에는 mapping CSV로 protomer identity를 복원한다.

### 3개 receptor state 문제

현재 local input 기준으로 `3GT8_raw.pdb`는 chain A/B가 있어 dimer source로 바로 쓸 수 있다. 반면 `EGFR_160-185.pdb`, `EGFR_170-200.pdb`는 chain X monomer만 있다. 이 상태에서 논문급 dimer-only 설계를 하려면 MD states의 dimer를 반드시 새로 정의해야 한다.

추천 방식:

1. 3GT8 dimer를 template으로 둔다.
2. MD cluster monomer를 3GT8 chain A에 superpose하여 receiver/monomer-A 위치를 만든다.
3. 같은 MD cluster monomer 사본을 3GT8 chain B에 superpose하여 activator/monomer-B 위치를 만든다.
4. 두 protomer의 clash, interface geometry, RMSD, buried surface를 기록한다.
5. 이 dimer가 modeled dimer임을 Methods와 Limitations에 명시한다.

대체 방식:

- MD trajectory 원본에 dimer simulation이 있으면 거기서 cluster representative dimer를 직접 추출한다. 이쪽이 가장 강하다.
- 없으면 위 template-superposition dimer를 쓰되, 논문에서는 "state-specific EGFR kinase conformers were embedded into a common 3GT8-derived dimer geometry"라고 투명하게 쓴다.

중단 조건:

- MD-state dimer source/modeling 방법을 기록할 수 없으면 3-state dimer conclusion을 쓰면 안 된다.

## 5. 실험 설계

### 변수 정의

고정:

- EGFR dimer receptor.
- MYO1D 955-1006.
- docking engine, score function, filters.
- residue numbering.

변수:

- EGFR conformational state: `3GT8_raw`, `EGFR_160-185`, `EGFR_170-200`.
- stochastic seed.

주 실험:

| Factor | Level |
|---|---|
| EGFR state | 3 |
| EGFR oligomer | dimer only |
| MYO1D construct | 955-1006 only |
| seeds | 10 per state, or minimum 5 |
| models/seed | 20,000 |

총량:

- 3 x 10 x 20,000 = 600,000 models.
- 만약 5 seeds만 사용하면 300,000 models이고, 논문에는 "five independent seeds"로 명시한다.

### 왜 MYO1D construct 여러 개를 주 실험에 섞지 않는가

논문 주 실험에서는 955-1006 하나만 써야 한다. 960/962 시작 construct를 섞으면 다시 "partner range"가 confounder가 된다. 다만 supplementary control로는 가능하다.

Supplementary controls:

- `960-1006` or `962-1006`: artifact demonstration only.
- sheet8/9 alanine mimic: computational negative control, 가능하면 나중에.
- randomized/back-face orientation: orientation filter validation only.

## 6. 검증 구조

### 입력 구조 검증

논문급으로는 input validation table이 필요하다.

필수 컬럼:

- receptor_state
- dimer_source
- protomer_A_range
- protomer_B_range
- protomer_B_offset_range
- MYO1D_range
- n_missing_backbone
- n_clashes_after_preparation
- dimer_interface_RMSD_to_template
- dimer_BSA
- mapping_csv
- input_pdb_sha256

### Docking output 검증

모든 seed에 대해 다음을 기록한다.

- total generated models.
- final ranked models.
- filter pass counts.
- mean/median dG_separated.
- mean dSASA.
- shape complementarity.
- orientation pass/fail/ambiguous counts.
- seed_complete marker.

### MYO1D orientation 검증

orientation filter는 단순 후처리가 아니라 논문 Methods의 핵심이다.

- active-face residues: 961-964 and 968-972.
- normal vector: PCA of active-face C-alpha plane.
- direction vector: active-face centroid to receptor interface centroid.
- pass: dot product > +0.10.
- fail: dot product < -0.10.
- ambiguous: absolute value <= 0.10.

ambiguous는 버리지 않았다고 쓰기보다 "excluded from consensus but reported"라고 써야 한다.

### Consensus 기준

결론은 single best pose가 아니라 consensus로만 낸다.

Hotspot acceptance:

- orientation-pass models only.
- within-cluster occupancy >= 0.50.
- appears in at least 2 independent seeds within a state.
- robust hotspot: appears in at least 2 of 3 receptor states.
- pan-state hotspot: appears in all 3 states.

MYO1D-side sanity check:

- accepted clusters should enrich sheet 8/9 contacts.
- sheet12-only clusters are not primary evidence.
- sheet8/9 enrichment should be reported as consistency with Ko et al., not as independent experimental proof.

## 7. 교차검증

LightDock은 필수는 아니지만 논문 설득력에 중요하다.

규칙:

- 같은 EGFR dimer input과 같은 MYO1D 955-1006 input 사용.
- scoring은 LightDock/DFIRE2 or fastdfire.
- PyRosetta와 점수를 합산하지 않는다.
- residue-level receptor patch overlap과 centroid agreement만 비교한다.

Convergence classes:

- PyRosetta + LightDock both support: high confidence.
- PyRosetta only: docking-engine-specific, review-required.
- LightDock only: secondary candidate.
- neither robust: reject.

## 8. Downstream pocket 설계

PPI patch가 확정된 뒤에만 pocket analysis로 넘어간다.

1. fpocket/P2Rank로 EGFR dimer 표면 포켓 탐지.
2. 포켓은 protomer identity를 유지해서 annotate한다.
3. ATP site는 strong candidate에서 제외.
4. pocket-PPI relationship:
   - interface-overlap.
   - rim.
   - allosteric-near.
   - allosteric-far.
   - irrelevant.
5. focused Vina는 dimer-derived pocket box를 사용한다.

논문용으로는 "MYO1D를 직접 막는 포켓"과 "MYO1D patch에 allosterically coupled된 포켓"을 구분해야 한다.

## 9. Methods에 쓸 수 있는 문장 골격

> EGFR kinase-domain dimers were used throughout the protein-protein docking workflow because EGFR kinase activation and intracellular receptor organization are governed by kinase-domain dimer geometry. For each receptor conformational state, a dimeric receptor model was prepared and represented as a single receptor chain for docking, with the second protomer residue numbers offset by +1000 to preserve protomer identity during downstream analysis.

> The MYO1D partner was restricted to residues 955-1006 of the C-terminal TH1 beta-meander. This construct retains the experimentally implicated beta-sheet region while avoiding an artificial N-terminal boundary at VAL962 observed in shorter constructs.

> Docking poses were not interpreted at the single-model level. Instead, receptor-side interface residues were aggregated across independent stochastic seeds, filtered by an active-face orientation criterion, and summarized as state-specific and cross-state consensus patches.

## 10. References

- Zhang X. et al. An allosteric mechanism for activation of the kinase domain of epidermal growth factor receptor. Cell. 2006. DOI: https://doi.org/10.1016/j.cell.2006.05.013
- RCSB PDB 2GS6, active EGFR kinase domain complex: https://www.rcsb.org/structure/2GS6
- Jura N. et al. Mechanism for activation of the EGF receptor catalytic domain by the juxtamembrane segment. Cell. 2009. DOI: https://doi.org/10.1016/j.cell.2009.04.025
- RCSB PDB 3GT8, inactive EGFR kinase domain with AMP-PNP: https://www.rcsb.org/structure/3GT8
- Ko Y.S. et al. MYO1D binds with kinase domain of the EGFR family to anchor them to plasma membrane before their activation and contributes carcinogenesis. Oncogene. 2019. DOI: https://doi.org/10.1038/s41388-019-0954-8
- Correction to Ko Y.S. et al. Oncogene. 2021. DOI: https://doi.org/10.1038/s41388-021-01675-y
- NCBI Gene MYO1D, conserved Myosin_TH1 domain annotation: https://www.ncbi.nlm.nih.gov/gene/4642

## 11. 이 설계가 이전 계획과 다른 점

이전의 seed 0-4 dimer, seed 5-9 monomer 설계는 "dimer vs monomer"라는 질문에는 맞지만, 사용자의 목표인 "EGFR는 dimer로 반드시 사용"에는 맞지 않는다. 논문 주장의 중심이 EGFR-MYO1D 생물학이라면 receptor oligomer state를 변수로 만들면 안 된다. dimer는 고정 조건이어야 하고, 재현성은 seed와 receptor conformational ensemble로 평가해야 한다.

따라서 실행 계획은 다음 방향으로 수정되어야 한다.

1. monomer main run 제거.
2. 3개 state 모두 dimer input 생성.
3. 30개 config 모두 dimer input 사용.
4. seed는 dimer replicate로만 해석.
5. MYO1D는 955-1006 하나를 주 construct로 고정.
6. shorter MYO1D constructs는 supplementary artifact/control로만 사용.
