# GPT Pro full-reasoning handoff: EGFR-MYO1D PPI, membrane-aware pocket discovery, and compound nomination

작성일: 2026-04-27  
작성 목적: 사용자가 GPT Pro의 "full reasoning" 또는 고강도 추론 모드로 전체 연구 설계를 다시 검증하고 완성하기 위해, 현재 프로젝트 상황과 사용자의 핵심 요청을 손실 없이 전달한다.

Important framing:

- 사용자는 기존 결과를 복구하거나 그대로 재검증하려는 것이 아니다.
- 이번 연구는 새로 수행할 fresh run이다.
- 기존 Workflow A/B 결과, PKT07/PKT34 후보, manuscript draft의 결론은 참고용 reasoning context일 뿐이다.
- GPT Pro는 과거 output을 찾거나 보존하는 데 시간을 쓰기보다, 기존 문서에서 배운 문제점과 제약을 이용해 새 dimer-only membrane-aware workflow를 설계해야 한다.

---

## 1. 사용자의 최종 목표

사용자의 최종 목표는 단순한 docking 결과 생성이 아니다. 최종 연구 목표는 다음과 같다.

> 세포막 근처에서 dimer 상태로 존재하는 EGFR kinase domain과 MYO1D TH1 beta-meander 사이의 PPI를 원자 수준에서 분석하고, 이 PPI를 교란할 수 있는 실제 druggable pocket을 찾은 뒤, 그 pocket에 docking 가능한 compound 후보를 도출하는 것.

즉 최종 산출물은 아래 3단계를 모두 만족해야 한다.

1. **EGFR-MYO1D PPI interface 정의**
   - MYO1D가 EGFR kinase-domain dimer의 어느 표면을 반복적으로 인식하는지 찾는다.
   - 단일 best pose가 아니라 seed/state를 넘어서 재현되는 consensus patch를 찾아야 한다.

2. **PPI 교란 가능 pocket 정의**
   - pocket은 ATP pocket이 아니어야 한다.
   - pocket은 EGFR-MYO1D interface 또는 그 바로 주변 rim/allosteric 위치에 있어야 한다.
   - pocket은 EGFR dimer와 membrane-proximal geometry를 고려했을 때 접근 가능해야 한다.

3. **Compound 후보 도출**
   - pocket이 실제 소분자 결합에 적합한지 검증한다.
   - probe ligand가 아니라 실제 후보 compound 또는 fragment library 기반 shortlist를 제시할 수 있어야 한다.
   - docking score만으로 결론내리지 않고, pose convergence, pocket relevance, dimer/membrane accessibility, ADMET/basic chemical filters를 함께 본다.

---

## 2. 사용자의 핵심 요청 정리

이번 대화에서 사용자가 명시한 핵심 조건은 다음과 같다.

### 2.1 EGFR 조건

- EGFR는 **무조건 dimer**를 사용해야 한다.
- monomer 결과는 main conclusion에 들어가면 안 된다.
- EGFR는 세포막에 붙어 있는 receptor tyrosine kinase이므로, kinase domain의 공간적 배치가 membrane-proximal geometry와 맞아야 한다.
- pocket은 임의의 표면 pocket이 아니라, EGFR dimer의 **측면/lateral side** 또는 **아래쪽/lower, membrane-proximal side**에 있어야 한다.
- central dimer interface 내부에 묻힌 pocket이나, membrane 방향에서 접근 불가능한 pocket은 compound target으로 부적절하다.

### 2.2 MYO1D 조건

- MYO1D의 핵심 결합부위인 **beta-sheet 8, 9, 12**는 반드시 포함되어야 한다.
- 이전 설계에서는 뒤쪽 C-terminal tail이 너무 길어 terminal-tail contact noise가 발생할 수 있었다.
- 따라서 핵심 beta-sheet는 유지하되, tail이 단백질과 비특이적으로 결합해 false positive를 만드는 문제를 막아야 한다.
- 사용자는 이전에 설계한 "말단이 단백질과 결합해서 노이즈가 발생하지 않도록 marking"하는 개념을 언급했다.
- 이 marking은 pipeline에서 `critical_residues_B`, `non_binding_residues_B`, `key_residues_B` 같은 residue annotation으로 구현 가능하다.

### 2.3 연구 품질 조건

- 단순한 실행 계획이 아니라 **논문에 기재 가능한 수준**의 과학적, 논리적 검증이 필요하다.
- 변수 confounding을 피해야 한다.
- dimer/monomer, MYO1D fragment length, receptor state, seed를 동시에 바꾸면 안 된다.
- 최종 결론은 single best model이 아니라 ensemble/consensus evidence에 기반해야 한다.
- pocket 또는 compound 후보는 "점수가 좋아 보인다"가 아니라 "EGFR-MYO1D PPI를 교란할 수 있는 생물물리적 위치와 접근성"을 가져야 한다.

---

## 3. 현재 프로젝트 상황

### 3.1 작업 디렉토리

현재 로컬 작업 디렉토리:

```text
Z:\!%data1\eunae\Sync_pororo\protein_ligand_hwang\last_result\codex_ligand
```

현재 날짜:

```text
2026-04-27
```

사용자 시간대:

```text
Asia/Seoul
```

### 3.2 Git 상태

사용자가 "깃허브 메인에 들어있을텐데 깃풀하고 확인해봐"라고 요청하여 `origin/main`에서 `git pull --ff-only origin main`을 수행했고, `docs/rerun_agent_prompt.md`가 추가된 것을 확인했다.

현재 worktree에는 사용자가 만들었거나 이전 작업에서 생긴 변경이 존재한다. 되돌리면 안 된다.

관찰된 dirty/untracked 파일:

```text
M .claude/hooks/csv-schema-guard.py
M .claude/hooks/pre-commit.sh
M scripts/deploy_fresh_n2000.sh
M scripts/extract_all_results.sh
?? .claude/settings.local.json
?? docs/dimer_only_paper_grade_workflow_design.md
?? docs/membrane_aware_ppi_pocket_compound_workflow.md
?? docs/rerun_workflow_design.md
?? docs/gpt_pro_full_reasoning_context_and_prompt.md
```

주의:

- 위 변경 중 `.claude/*`, `scripts/*`는 현재 요청과 직접 관련 없는 기존 변경이다. 되돌리지 않는다.
- 새로 작성한 문서는 현재 연구 설계 정리를 위한 것이다.

### 3.3 이미 읽고 확인한 주요 문서

#### `docs/rerun_agent_prompt.md`

핵심 내용:

- 기존 seed 0-4와 seed 5-9 실험은 confounding 문제가 있었다.
- seed 0-4는 dimer + MYO1D 960-1006, seed 5-9는 monomer + MYO1D 955-1006였던 것으로 정리되어 있었다.
- 이 경우 receptor construct와 partner range가 동시에 바뀌므로 결과 차이를 해석할 수 없다.
- 원래 rerun prompt는 "모든 seed에서 partner를 955-1006으로 통일하고 dimer/monomer만 변수화"하는 방향을 제시했다.
- 하지만 사용자가 이후 **EGFR는 dimer를 무조건 사용해야 한다**고 명시했으므로, 지금 최종 설계에서는 dimer-vs-monomer 비교가 main experiment가 아니다.

#### `docs/CONTEXT.md`

핵심 내용:

- 이전 workflow A/B가 완료된 이력이 있다.
- MYO1D side active-face top contacts로 VAL962, VAL964, CYS970, SER971 등이 관찰되었다.
- sheet12는 낮은 occupancy의 support signal로 관찰되었다.
- `_capped` PDB 관련 주의가 있다. `_capped` 파일이 실제로 다른 구조였던 문제가 있었고, 3GT8_raw seed5 `_capped` 결과는 폐기 후 normal PDB로 rerun한 이력이 있다.
- `input/` symlink 문제로 `TH1 domain.pdb`가 missing될 수 있다. HPC 환경에서는 원본 `input/` 연결이 필요할 수 있다.

#### `docs/archive/design_intent.md`

매우 중요한 설계 의도:

- 전체 TH1 domain을 docking partner로 쓰면 비특이적 표면 접촉 noise가 커졌다.
- 962-1006처럼 너무 짧게 자르면 VAL962가 artificial N-terminus가 되어 terminal charge/backbone freedom artifact가 생길 수 있었다.
- 그래서 955-1006으로 N-terminal buffer를 확장한 설계가 나왔다.
- MYO1D beta-meander sheet는 얇은 ribbon 형태이므로 contact count만으로 active face/back face를 구분할 수 없다.
- orientation filter가 필수이며, sheet8/9 active-face normal이 receptor 방향을 향해야 한다.

### 3.4 현재 코드 구조에서 확인한 점

#### PPI input/config generation 관련

- `egfr_pipeline/phase1/generate_configs.py`
  - 현재 모든 seed를 monomer full_kinase_domain, partner 955-1006으로 생성하는 상태가 관찰되었다.
  - dimer-only production 설계와 맞지 않을 수 있으므로 수정 또는 별도 dimer-only config generator가 필요하다.

- `run_production.py`
  - `_build_ppi_targets`에서 모든 seed에 `construct_type="full_kinase_domain"`와 `mapping_csv=""`를 넣는 구조가 확인되었다.
  - dimer/protomer mapping metadata가 충분히 반영되지 않을 수 있다.

- `egfr_pipeline/ppi/prepare_dimer_pdb.py`
  - EGFR dimer chain A+B를 chain A 하나로 병합하고, 원래 chain B residue에 +1000 offset을 붙이는 기능이 있다.
  - MYO1D는 chain B로 붙인다.
  - mapping CSV와 info JSON을 출력할 수 있다.

#### receptor input 상태

확인된 local receptor PDB chain 상태:

- `3GT8_raw.pdb`
  - chain A: residues 699-1007
  - chain B: residues 701-1007
  - dimer source로 바로 사용할 수 있음

- `EGFR_160-185.pdb`
  - chain X only
  - residues 634-1014
  - 현재 monomer만 있음

- `EGFR_170-200.pdb`
  - chain X only
  - residues 634-1014
  - 현재 monomer만 있음

해석:

- 3GT8_raw는 dimer source가 있다.
- MD-derived states는 현재 monomer-only로 보이며, dimer-only 논문 결론에 사용하려면 template-superposition dimer 또는 true dimer trajectory source가 필요하다.
- 이 부분은 최종 연구 설계의 핵심 취약점이다.

#### PPI scoring/annotation 관련

`egfr_pipeline/pyrosetta_docking/pipeline_manager.py`에서 확인된 기능:

- `[Constraints] key_residues_B`
- `[Constraints] key_residue_weights`
- `[Constraints] key_residue_bonus_weight`
- `[ExperimentalData] critical_residues_B`
- `[ExperimentalData] non_binding_residues_B`

중요한 해석:

- `critical_residues_B`와 `non_binding_residues_B`는 validation/report 쪽에서 sensitivity/specificity/false-positive 성격의 해석 지표로 쓰인다.
- 이것만으로 docking search를 hard-block하지는 않는다.
- `key_residue_bonus_weight > 0`이면 adjusted score에 영향을 주어 ranking이 바뀔 수 있다.
- 따라서 paper-grade rerun의 첫 번째 설계에서는 `key_residue_bonus_weight = 0.0`으로 두고, key-contact ratio를 기록만 하는 것이 가장 방어 가능하다.
- scoring weight를 바꾸려면 사용자의 명시적 승인이 필요하다.

#### pocket/compound pipeline 관련

이미 존재하는 advanced pipeline:

```text
PPI docking -> PPI patch -> fpocket/P2Rank -> patch relationship -> focused Vina -> perturbation score
```

관련 모듈:

- `egfr_pipeline/phase2/pocket_proposal.py`
  - fpocket/P2Rank pocket proposal setup/parser.

- `egfr_pipeline/phase2/patch_relationship.py`
  - pocket과 PPI patch의 관계를 orthosteric/rim/allosteric/low_relevance로 분류.
  - 현재 기준:
    - orthosteric: hotspot overlap >= 2 and fraction >= 0.25
    - rim: hotspot overlap >= 1
    - allosteric: centroid distance <= 20 A and no overlap
    - low relevance: >20 A and no overlap

- `egfr_pipeline/phase2/druggability_confidence.py`
  - fpocket/P2Rank druggability score 및 tier 부여.

- `egfr_pipeline/phase2/phase3_export.py`
  - Phase 3 focused docking candidate pocket reference를 export.

- `egfr_pipeline/phase3/*`
  - focused Vina docking, diversity-aware docking, pose attribution, phase4 export.

현재 부족한 점:

- dimer/membrane-aware pocket geometry filter가 없다.
- "lower/lateral side"와 "membrane-proximal accessibility"가 pocket gating에 포함되어 있지 않다.
- pocket proposal이 monomer receptor 기준으로 수행되었을 가능성이 있으며, dimer receptor에서 재실행되어야 한다.

---

## 4. 이미 작성/추가한 설계 문서

### 4.1 `docs/dimer_only_paper_grade_workflow_design.md`

역할:

- EGFR를 dimer-only로 고정한 논문급 rerun 설계를 정리한다.

주요 내용:

- EGFR receptor는 항상 dimer.
- monomer는 main conclusion에서 제외.
- MYO1D construct는 beta-sheet 8/9/12를 포함해야 한다.
- `955-1001`을 candidate main construct로 제안.
- `955-1006`은 conservative tail-masked comparator로 제안.
- `key_residue_bonus_weight=0.0`으로 시작할 것을 권장.
- sheet8/9 active face와 sheet12 support를 구분.
- dimer receptor state별 준비와 mapping CSV 필요성을 정리.

### 4.2 `docs/membrane_aware_ppi_pocket_compound_workflow.md`

역할:

- PPI 분석에서 pocket discovery와 compound nomination까지 이어지는 최종 workflow를 새로 설계한다.

주요 내용:

- 최종 objective를 PPI-to-pocket-to-compound로 재정의.
- EGFR dimer와 membrane geometry를 필수 조건으로 둠.
- pocket은 lower/lateral membrane-compatible side에 있어야 함.
- virtual membrane frame 정의 필요성을 제시.
- `3GOP`, `2M20`, `3GT8`를 reference로 사용하자고 제안.
- dimer receptor 기반 pocket discovery, ATP exclusion, PPI proximity, membrane-side gate, dimer-accessibility gate, state support를 hard gate로 정의.
- 최종 compound 후보는 focused docking, pose convergence, PPI perturbation mechanism, ADMET/basic chemistry filter를 통과해야 함.

---

## 5. 과학적 근거 요약

### 5.1 EGFR dimer와 membrane geometry가 필요한 이유

EGFR는 receptor tyrosine kinase이고, kinase activation은 asymmetric kinase-domain dimer와 juxtamembrane/TM geometry에 강하게 연결되어 있다. EGFR kinase domain만 떼어낸 monomer 표면은 계산상 접근 가능해 보여도, 실제 세포막 근처의 dimer context에서는 접근 불가능하거나 생물학적으로 무의미할 수 있다.

설계상 중요한 구조 reference:

- `3GOP`: EGFR juxtamembrane and kinase domains. EGFR JM/kinase asymmetric dimer reference로 중요.
- `2M20`: EGFR transmembrane-juxtamembrane segment in bicelles. membrane/TM-JM orientation reference로 중요.
- `3GT8`: inactive EGFR kinase-domain structure/dimer context reference.

핵심 설계 원칙:

- EGFR dimer는 receptor unit으로 유지한다.
- pocket은 dimer context에서 다시 찾아야 한다.
- pocket은 dimer interface 내부나 membrane-occluded region에 있으면 안 된다.
- lower/lateral side에 존재해야 compound가 세포질 쪽에서 접근 가능하고 MYO1D interaction을 perturb할 가능성이 있다.

### 5.2 MYO1D construct 설계 이유

Ko et al. 2019는 MYO1D TH1 domain의 beta-meander motif가 EGFR family kinase domain binding에 중요하다고 보고했다. 사용자의 조건상 beta-sheet 8, 9, 12는 반드시 포함되어야 한다.

하지만 계산 설계에서는 양쪽 artifact를 피해야 한다.

- 너무 긴 TH1: 비특이적 표면 접촉 noise.
- 962부터 시작하는 짧은 fragment: VAL962 artificial N-terminus artifact.
- 955-1006: N-terminal artifact는 줄였지만 C-terminal tail 998-1006이 terminal-tail false contact를 만들 가능성.

따라서 제안된 설계:

- candidate main: `955-1001`
  - 955-960: N-terminal buffer
  - 961-964, 968-972: beta-sheet 8/9 active-face
  - 993-997: beta-sheet 12 support
  - 998-1001: short C-terminal cap/buffer
  - 1002-1006: main candidate에서 제거

- conservative comparator: `955-1006_tail_masked`
  - sheet12 안정성이 `955-1001`에서 무너지면 사용
  - 998-1006을 non-binding/noise zone으로 annotation

권장 annotation:

```ini
[Constraints]
key_residues_B = 961-964,968-972,993-997
key_residue_bonus_weight = 0.0

[ExperimentalData]
critical_residues_B = 961-964,968-972,993-997
non_binding_residues_B = 998-1006
```

해석:

- `critical_residues_B`는 "이 residue가 결합 증거로 중요하다"는 marking이다.
- `non_binding_residues_B`는 "tail contact는 false positive/noise로 봐야 한다"는 marking이다.
- `key_residue_bonus_weight = 0.0`은 ranking을 바꾸지 않기 위한 선택이다.

### 5.3 Pocket discovery의 핵심 논리

PPI-targeted pocket discovery에서는 가장 깊고 점수가 좋은 pocket이 꼭 좋은 target은 아니다. EGFR에서는 ATP cleft가 매우 강한 false attractor가 될 수 있다.

따라서 pocket selection은 다음 순서여야 한다.

1. ATP pocket 제외.
2. dimer receptor에서 존재 확인.
3. MYO1D PPI patch와 직접 overlap, rim, 또는 가까운 allosteric relationship.
4. lower/lateral membrane-compatible location.
5. dimer interface 내부에 묻히지 않음.
6. fpocket/P2Rank 등 다중 tool 또는 state recurrence로 support.
7. compound docking이 pocket 내부에서 수렴.
8. compound pose가 MYO1D approach path 또는 PPI-coupled residues를 perturb할 수 있음.

---

## 6. 최적 워크플로우 제안

### Phase 0. Input and biological-frame preparation

목표:

- dimer receptor와 virtual membrane coordinate frame을 준비한다.

필수 작업:

1. 모든 receptor state에 대해 EGFR dimer를 준비한다.
2. 3GT8_raw는 dimer chain A/B를 사용할 수 있다.
3. EGFR_160-185, EGFR_170-200은 현재 monomer-only이므로 다음 중 하나를 선택해야 한다.
   - true dimer trajectory representative가 있으면 그것을 사용.
   - 없으면 3GT8 또는 3GOP dimer template에 state-specific monomer를 superpose하여 modeled dimer를 생성.
4. dimer construction method, RMSD, clash, dimer interface geometry, protomer mapping을 모두 기록한다.
5. `3GOP`와 `2M20` reference를 이용해 virtual membrane frame을 정의한다.

출력 제안:

```text
output/phase0_geometry/membrane_frame.json
output/phase0_geometry/dimer_geometry_qc.csv
```

중요 QC:

- dimer interface clash.
- protomer identity.
- dimer interface RMSD to template.
- membrane normal vector.
- lower/lateral coordinate convention.

### Phase 1. Dimer-only MYO1D PPI docking

목표:

- MYO1D가 dimeric EGFR의 lower/lateral membrane-compatible surface를 반복적으로 인식하는지 찾는다.

권장 실험:

- receptor: EGFR dimer only.
- MYO1D: 먼저 `955-1001` vs `955-1006_tail_masked` small pilot.
- pilot에서 sheet12 구조 보존과 tail false contact 감소를 비교.
- main production에는 하나의 MYO1D construct만 사용.
- receptor state x seed 반복.

필터:

- MYO1D active-face orientation pass.
- sheet8/9 contact present.
- sheet12 support tracked.
- tail-dominant pose reject/quarantine.
- membrane-compatible MYO1D placement.
- EGFR contact centroid lower/lateral.
- central dimer-interface buried contact reject.
- ATP-site overlap reject.

출력 제안:

```text
output/phase1_ppi/ppi_pose_membrane_qc.csv
output/phase1_ppi/ppi_consensus_patch_dimer_aware.csv
```

### Phase 2. Dimer/membrane-aware pocket proposal

목표:

- dimer receptor에서 실제 접근 가능한 pocket을 찾는다.

필수:

- fpocket/P2Rank는 dimer receptor input에 대해 실행.
- pocket residue IDs는 protomer-aware하게 기록.
- pocket centroid를 membrane frame에 투영.

Hard gate:

| Gate | 내용 |
|---|---|
| G0 | dimer receptor에서 발견된 pocket |
| G1 | ATP site 아님 |
| G2 | dimer-aware PPI patch와 orthosteric/rim/allosteric 관계 |
| G3 | lower/lateral membrane-compatible side |
| G4 | central dimer interface 내부에 묻히지 않음 |
| G5 | state robust 또는 reproducible cryptic pocket |

출력 제안:

```text
output/phase2_pockets/pocket_membrane_geometry.csv
output/phase2_pockets/pocket_dimer_accessibility.csv
output/phase2_pockets/phase3_candidate_pocket_reference_membrane_aware.csv
```

### Phase 3. Focused compound docking

목표:

- surviving pocket에만 compound docking을 수행한다.

중요:

- 기존 3개 ligand는 probe로만 사용.
- 실제 compound discovery에는 fragment library 또는 PPI-oriented library가 필요하다.
- compound library는 PAINS/reactive filter, basic physicochemical filter, diversity selection을 거쳐야 한다.

권장 docking:

- surviving pocket별 focused box.
- multiple receptor states.
- multiple seeds.
- pose clustering.
- ligand diversity tracking.

### Phase 4. Compound mechanism and stability validation

목표:

- docking score가 아니라 PPI perturbation mechanism을 가진 compound shortlist를 만든다.

필터:

- ATP site binding 아님.
- pocket pose가 lower/lateral dimer side에서 가능.
- compound가 MYO1D approach path를 sterically block하거나 PPI-coupled residues와 상호작용.
- pose convergence가 있음.
- local relaxation 또는 short MD에서 pocket occupancy가 유지.
- basic ADMET/PAINS/reactive flags 통과.

출력:

```text
output/final_compounds/membrane_aware_compound_shortlist.csv
```

---

## 7. 논문급 결론을 위한 판정 기준

### 7.1 Final pocket claim 기준

어떤 pocket을 EGFR-MYO1D PPI modulator pocket으로 주장하려면 모두 만족해야 한다.

1. EGFR dimer receptor에서 발견되었다.
2. ATP pocket이 아니다.
3. orientation-filtered MYO1D PPI patch와 연결된다.
4. lower/lateral membrane-compatible side에 있다.
5. dimer geometry에서 compound 접근이 가능하다.
6. state robust 또는 reproducible cryptic pocket이다.
7. focused compound docking에서 pose convergence가 있다.
8. 제안 mechanism이 generic kinase inhibition이 아니라 PPI perturbation이다.

### 7.2 Final compound claim 기준

어떤 compound를 후보로 제안하려면 모두 만족해야 한다.

1. accepted pocket에 반복적으로 docking된다.
2. ATP-site binding으로 설명되지 않는다.
3. membrane/dimer geometry와 충돌하지 않는다.
4. PPI patch와 연결된 residues 또는 approach path를 perturb한다.
5. chemically plausible하다.
6. local relaxation 또는 short MD에서 pose가 유지된다.

---

## 8. 현재 설계에서 가장 큰 리스크

### Risk 1. MD-derived states가 monomer-only일 가능성

문제:

- 3GT8_raw는 dimer지만 EGFR_160-185/170-200은 local input 기준 monomer-only였다.
- dimer-only 논문 결론을 위해서는 이 두 state도 dimer로 만들어야 한다.

해결:

- true dimer MD trajectory가 있으면 그 representative를 사용.
- 없으면 3GT8/3GOP template-superposition dimer를 만들고, "modeled dimer"로 명시.
- Methods/Limitations에 반드시 기록.

### Risk 2. MYO1D tail artifact

문제:

- `955-1006` tail이 EGFR에 비특이적으로 붙으면 false positive hotspot이 생길 수 있다.

해결:

- `955-1001` candidate와 `955-1006_tail_masked`를 small pilot 비교.
- sheet12 preservation을 확인한 뒤 main construct를 하나로 결정.
- tail contact fraction을 계산해 tail-dominant pose reject/quarantine.

### Risk 3. ATP pocket dominance

문제:

- blind docking은 ATP pocket을 너무 잘 찾는다.
- 이는 EGFR-MYO1D PPI modulator discovery에는 false attractor이다.

해결:

- ATP pocket은 final STRONG 후보에서 제외.
- PPI-first focused docking을 중심으로 한다.

### Risk 4. Pocket이 실제 membrane context에서 접근 불가능할 수 있음

문제:

- kinase-only monomer surface에서 보이는 pocket이 실제 dimer/membrane geometry에서는 접근 불가능할 수 있다.

해결:

- virtual membrane frame을 만들고 lower/lateral gate를 추가.
- dimer interface buried pocket을 제외.

### Risk 5. 기존 manuscript draft의 후보를 새 연구의 출발점으로 오해할 수 있음

문제:

- 기존 `PKT07`, `PKT34` style conclusions는 historical result이다.
- 사용자는 어차피 새로 수행할 계획이므로, 이 후보들을 반드시 복구하거나 재현할 필요는 없다.

해결:

- 이 후보들은 "이전 workflow가 어떤 종류의 pocket을 뽑았는지 보여주는 참고 사례"로만 사용한다.
- 새 workflow의 pocket 후보는 fresh dimer/membrane-aware run에서 새로 도출한다.

---

## 9. 구현 관점의 권장 변경

### 새 모듈 제안

```text
egfr_pipeline/phase0/membrane_frame.py
egfr_pipeline/phase1/membrane_pose_filter.py
egfr_pipeline/phase2/membrane_geometry.py
egfr_pipeline/final_compounds/compound_shortlist.py
```

### 수정 대상 제안

```text
egfr_pipeline/phase1/generate_configs.py
egfr_pipeline/phase1/launch_docking.py
run_production.py
egfr_pipeline/phase2/patch_relationship.py
egfr_pipeline/phase2/phase3_export.py
egfr_pipeline/phase3/run_diverse_docking.py
```

### 수정 시 원칙

- `egfr_pipeline/paths.py`는 되도록 수정하지 않는다.
- 기존 CSV column은 rename/delete하지 않는다.
- 필요한 정보는 새 column으로 추가한다.
- scoring weight는 승인 없이 변경하지 않는다.
- 기존 ligand 3종은 유지하되, 실제 compound discovery에는 별도 library를 추가하는 방향으로 설계한다.

---

## 10. GPT Pro가 풀 추론으로 검토해야 할 질문

다음 질문에 대해 깊게 검토해야 한다.

1. EGFR kinase-domain dimer와 membrane frame을 가장 방어 가능하게 정의하는 방법은 무엇인가?
2. 3GT8/3GOP/2M20를 어떻게 조합해야 lower/lateral membrane-proximal side를 정량화할 수 있는가?
3. MD-derived monomer state를 dimer로 embed하는 것이 논문에서 허용 가능한가? 허용된다면 어떤 QC가 필요한가?
4. MYO1D construct는 `955-1001`이 좋은가, 아니면 `955-1006_tail_masked`가 더 방어 가능한가?
5. beta-sheet 8/9/12를 모두 포함하면서 terminal artifact를 최소화하는 최적 residue range는 무엇인가?
6. MYO1D active-face orientation filter와 membrane-compatible pose filter를 어떻게 결합해야 하는가?
7. EGFR-side PPI patch를 consensus로 정의할 때 seed/state/protomer를 어떻게 통합해야 하는가?
8. dimer-aware pocket discovery에서 protomer A/B pocket을 같은 pocket으로 merge해야 하는가, 아니면 별도 pocket으로 유지해야 하는가?
9. lower/lateral pocket gate의 수학적 정의는 어떻게 해야 하는가?
10. compound docking 전에 pocket을 hard gate로 거를 것과 soft score로 반영할 것을 어떻게 나눌 것인가?
11. 기존 `PKT07/PKT34` 후보는 새 설계에서 살아남을 가능성이 있는가? 살아남으려면 어떤 증거가 필요할까?
12. 실제 compound shortlist를 논문급으로 만들기 위해 최소한 어떤 library, filter, docking replicate, rescoring이 필요한가?

---

## 11. 참고 문헌 및 구조

설계 근거로 확인한 자료:

- Ko et al., 2019, MYO1D binds with kinase domain of the EGFR family to anchor them to plasma membrane before activation: https://www.nature.com/articles/s41388-019-0954-8
- PubMed entry for Ko et al. 2019: https://pubmed.ncbi.nlm.nih.gov/31420606/
- Ko et al., 2021 correction: https://www.nature.com/articles/s41388-021-01675-y
- RCSB 3GOP, EGFR juxtamembrane and kinase domains: https://www.rcsb.org/structure/3GOP
- RCSB 2M20, EGFR TM-JM segment in bicelles: https://www.rcsb.org/structure/2m20
- RCSB 3GT8, inactive EGFR kinase-domain structure: https://www.rcsb.org/structure/3GT8
- Jura et al., 2009, mechanism for activation of EGFR catalytic domain by juxtamembrane segment: https://doi.org/10.1016/j.cell.2009.04.025
- Zhang/Kuriyan EGFR asymmetric kinase dimer activation context, 2GS6 reference: https://www.rcsb.org/structure/2GS6

---

## 12. 그대로 사용할 GPT Pro 프롬프트

아래 프롬프트를 GPT Pro full-reasoning 모드에 그대로 넣는다.

```text
나는 EGFR-MYO1D PPI를 논문급으로 분석하고, 그 PPI를 교란할 수 있는 실제 druggable pocket과 compound 후보를 도출하려는 연구를 하고 있다. 단순 docking 실행 계획이 아니라, 논문 Methods/Results/Limitations에 기재 가능한 수준의 과학적이고 논리적인 전체 workflow를 완성해줘.

반드시 지켜야 하는 조건:

1. EGFR는 무조건 dimer로 사용해야 한다. monomer 결과는 main conclusion에 넣지 않는다.
2. EGFR는 세포막에 붙어 있는 receptor tyrosine kinase이므로, pocket은 dimer의 lower/lateral, membrane-proximal side에 있어야 한다.
3. pocket은 central dimer interface 내부에 묻혀 있거나 membrane geometry상 접근 불가능하면 안 된다.
4. ATP pocket은 최종 PPI-modulator pocket 후보에서 제외해야 한다.
5. MYO1D construct에는 핵심 결합 부위인 beta-sheet 8, 9, 12가 반드시 포함되어야 한다.
6. 이전에는 MYO1D C-terminal tail이 너무 길어 단백질에 비특이적으로 붙는 noise가 있었을 가능성이 있다. terminal-tail artifact를 줄이는 설계를 포함해줘.
7. MYO1D active-face는 sheet 8/9 residues 961-964, 968-972이고, sheet12 support region은 993-997로 본다.
8. `955-1001` construct와 `955-1006_tail_masked` comparator 중 어떤 설계가 더 논문급으로 방어 가능한지 검토해줘.
9. scoring weight를 임의로 바꾸는 설계는 피하고, 바꿔야 한다면 왜 필요한지와 승인/민감도 분석이 필요한 이유를 명시해줘.
10. 최종 목표는 EGFR-MYO1D PPI patch를 찾고, 그 주변의 실제 druggable pocket을 찾은 뒤, 그곳에 docking 가능한 compound 후보를 도출하는 것이다.

현재 프로젝트 상황:

- repo: `Z:\!%data1\eunae\Sync_pororo\protein_ligand_hwang\last_result\codex_ligand`
- 기존 pipeline에는 PyRosetta PPI docking, LightDock cross-validation, fpocket/P2Rank pocket proposal, AutoDock Vina focused docking, perturbation scoring framework가 있다.
- 기존 output/result 후보는 복구 대상이 아니다. 이번 목적은 기존 문서와 코드에서 배운 문제점을 반영해 fresh workflow를 설계하는 것이다.
- `egfr_pipeline/ppi/prepare_dimer_pdb.py`는 EGFR dimer chain A/B를 chain A 하나로 병합하고 원래 chain B residue에 +1000 offset을 붙이는 기능이 있다.
- `3GT8_raw.pdb`는 local 기준 chain A/B dimer가 있다.
- `EGFR_160-185.pdb`, `EGFR_170-200.pdb`는 local 기준 chain X monomer-only로 보인다. dimer-only conclusion을 위해서는 true dimer source를 찾거나 3GT8/3GOP template-superposition dimer를 만들어야 할 수 있다.
- 기존 design intent에 따르면 whole TH1은 noise가 크고, 962-1006은 VAL962 artificial N-terminus artifact가 가능해서 955-1006으로 확장했었다.
- 그러나 955-1006은 C-terminal 998-1006 tail이 noise를 만들 수 있어 `955-1001` 또는 tail-masked comparator가 필요하다.
- pipeline에는 `[Constraints] key_residues_B`, `key_residue_bonus_weight`, `[ExperimentalData] critical_residues_B`, `non_binding_residues_B`가 있다. `key_residue_bonus_weight=0.0`이면 ranking을 바꾸지 않고 key-contact ratio를 기록할 수 있다.
- 현재 advanced pipeline은 PPI docking -> PPI patch -> fpocket/P2Rank -> patch relationship -> focused Vina -> perturbation score 구조지만, 아직 dimer/membrane-aware pocket geometry gate가 충분하지 않다.

내가 원하는 출력:

1. 전체 workflow를 Phase별로 다시 설계해줘.
2. 각 Phase의 input, output, QC metric, fail condition, 논문에 쓸 수 있는 Methods 문장을 제시해줘.
3. EGFR dimer와 membrane frame을 어떻게 정의할지 구체적으로 제안해줘. 3GOP, 2M20, 3GT8를 어떻게 사용할지 설명해줘.
4. MYO1D construct 최종 선택 전략을 제시해줘. beta-sheet 8/9/12 보존, tail artifact 방지, orientation filter, tail marking을 포함해줘.
5. PPI patch consensus를 seed/state/protomer across-analysis로 어떻게 정의할지 제안해줘.
6. pocket discovery에서 hard gate와 soft score를 나눠줘. 특히 ATP exclusion, PPI proximity, lower/lateral membrane geometry, dimer accessibility, state robustness를 포함해줘.
7. compound docking과 compound shortlist 도출 전략을 제시해줘. 기존 3개 ligand를 probe로만 볼지, fragment/PPI-oriented library를 어떻게 도입할지 설명해줘.
8. 기존 PKT07/PKT34 같은 후보는 참고 사례로만 보고, 새 fresh run에서는 어떤 기준으로 새로운 pocket 후보를 도출해야 하는지 제시해줘.
9. 논문용 Results/Discussion에서 어떤 결론은 말할 수 있고 어떤 결론은 아직 말하면 안 되는지 구분해줘.
10. 마지막에 이 계획을 구현하기 위한 repo-level task list를 파일/모듈 단위로 제시해줘.

참고 구조/문헌:

- Ko et al. 2019 EGFR-MYO1D: https://www.nature.com/articles/s41388-019-0954-8
- Ko et al. 2021 correction: https://www.nature.com/articles/s41388-021-01675-y
- RCSB 3GOP EGFR juxtamembrane and kinase domains: https://www.rcsb.org/structure/3GOP
- RCSB 2M20 EGFR TM-JM segment in bicelles: https://www.rcsb.org/structure/2m20
- RCSB 3GT8 inactive EGFR kinase-domain structure: https://www.rcsb.org/structure/3GT8
- Jura et al. 2009 EGFR juxtamembrane activation mechanism: https://doi.org/10.1016/j.cell.2009.04.025

답변은 한국어로 작성하고, 과학적 근거와 계산 파이프라인 설계를 분리해서 매우 체계적으로 설명해줘. 단순한 아이디어가 아니라, 실제 논문과 repo 구현으로 이어질 수 있는 수준으로 작성해줘.
```
