# Final current-state summary: EGFR-MYO1D dimer/membrane-aware rerun

작성일: 2026-04-27  
목적: `docs/` 전체와 현재 workspace 상태를 대조하여, 다음 GPT Pro full-reasoning 또는 실제 구현 단계가 혼동하지 않도록 현재 상태, 충돌점, 확정 원칙, 다음 작업을 한 번에 정리한다.

Important framing update:

- The next study is a fresh rerun, not an attempt to preserve or reproduce the previous result set.
- Historical Workflow A/B outputs are useful only as reasoning context: they reveal prior design assumptions, false positives, runtime pitfalls, and useful module boundaries.
- Existing `PKT07/PKT34`-style candidates do not need to be recovered before planning the new workflow.
- Missing old `output/workflow_a` or `output/workflow_b` payloads are not blockers for the fresh design.
- The real blockers are fresh-run inputs and methods: dimer receptor construction, MYO1D construct selection, membrane frame definition, pocket gating, and compound-library strategy.

---

## 1. 최종 연구 목표

최종 목표는 EGFR-MYO1D docking 자체가 아니라 다음 연쇄를 논문급으로 완성하는 것이다.

```text
EGFR dimer membrane-proximal context
    -> MYO1D TH1 beta-meander PPI patch
    -> dimer/membrane-compatible druggable pocket
    -> focused compound docking
    -> PPI-disruptive compound shortlist
```

따라서 최종 결론은 다음 조건을 모두 만족해야 한다.

- EGFR는 main analysis에서 반드시 dimer이다.
- MYO1D는 beta-sheet 8, 9, 12를 포함해야 한다.
- PPI patch는 MYO1D active-face orientation과 tail-noise control을 통과해야 한다.
- pocket은 ATP pocket이 아니어야 한다.
- pocket은 EGFR dimer의 lower/lateral, membrane-proximal, compound-accessible 위치여야 한다.
- compound는 단순 Vina affinity가 아니라 PPI perturbation mechanism을 가져야 한다.

---

## 2. 이번 문서 검토에서 확인한 문서 그룹

### 2.1 새로 작성된 설계 문서

- `docs/dimer_only_paper_grade_workflow_design.md`
  - EGFR dimer-only 논문급 rerun 설계.
  - MYO1D `955-1001` candidate와 `955-1006_tail_masked` comparator를 제안.
  - scoring weight를 바꾸지 않고 `key_residue_bonus_weight = 0.0`으로 key-contact를 기록하는 방향.

- `docs/membrane_aware_ppi_pocket_compound_workflow.md`
  - PPI -> pocket -> compound까지 이어지는 membrane-aware workflow.
  - virtual membrane frame, lower/lateral gate, dimer-accessibility gate를 제안.
  - `3GOP`, `2M20`, `3GT8`를 geometry reference로 제안.

- `docs/gpt_pro_full_reasoning_context_and_prompt.md`
  - GPT Pro full-reasoning에 그대로 전달할 handoff 문서.
  - 사용자의 요구, 현재 repo 상태, scientific constraints, full prompt 포함.

### 2.2 기존 active 문서

- `docs/rerun_agent_prompt.md`
  - 원래 rerun의 출발점.
  - 과거 실험에서 seed 0-4와 seed 5-9 사이에 receptor construct와 MYO1D range가 동시에 바뀌어 confounding이 생긴 문제를 지적.
  - 단, 이 문서는 dimer-vs-monomer 비교 설계를 포함한다. 현재 사용자 조건에서는 monomer comparison은 main analysis에서 제외되어야 한다.

- `docs/CONTEXT.md`
  - 과거 HPC 실행 결과와 문제 해결 기록의 가장 중요한 상태 문서.
  - Workflow A/B 완료, 600K PPI, fpocket 165 -> 103, PKT07/PKT34 후보 등을 기록.
  - 하지만 이 결과는 현재 새 dimer/membrane-aware 기준으로는 final claim이 아니라 재검증 대상이다.

- `docs/data_inventory.md`
  - 입력/출력 파일의 기대 구조를 설명.
  - 현재 local snapshot과 불일치가 있다. 문서에는 raw PPI partner와 ligands가 active처럼 기술되어 있지만, 현재 filesystem에는 존재하지 않는다.

- `docs/architecture.md`, `docs/runbook.md`, `docs/workflow_comparison_guide.md`
  - Workflow A와 Workflow B의 구조를 설명.
  - Workflow B는 PPI-first -> pocket -> focused docking -> perturbation scoring 구조를 이미 갖고 있다.
  - 그러나 membrane/dimer-aware gate는 아직 충분히 반영되어 있지 않다.

- `docs/manuscript_draft.md`
  - 기존 결과 기반 논문 초안.
  - PKT07/PKT34를 tier-1 후보로 제시.
  - 현재 새 조건에서는 이 초안의 결론을 그대로 쓰면 안 된다. dimer/membrane-aware rerun 이후 수정해야 한다.

- `docs/methodology_limitations.md`, `docs/false_positive_report.md`, `docs/phase1_notes.md`
  - orientation filter, false-positive control, method limitation 관련 보조 문서.
  - tail artifact, ATP pocket dominance, single best pose 회피라는 현재 설계 방향과 일관된다.

### 2.3 archive 문서에서 유지해야 할 핵심 설계 의도

- `docs/archive/design_intent.md`
  - 전체 TH1 domain은 비특이적 noise가 커서 부적절.
  - `962-1006`은 VAL962 artificial N-terminus artifact 위험.
  - `955-1006`은 N-terminal artifact를 줄이기 위해 도입.
  - 그러나 현재 사용자의 문제 제기에 따라 C-terminal tail artifact는 별도로 관리해야 한다.
  - MYO1D beta-meander는 얇은 sheet 구조이므로 contact count만으로 active/back face를 구분할 수 없고 orientation filter가 필수.

---

## 3. 현재 local filesystem 상태

### 3.1 현재 존재하는 input 파일

현재 `input/` 아래에서 확인된 실제 파일은 다음뿐이다.

```text
input/PPI/phase1/docking_pair_metadata.csv
input/PPI/phase1/partner_metadata.csv
input/PPI/phase1/phase1_input_validation_report.md
input/PPI/phase1/pilot_data_reference.csv
input/PPI/phase1/receptor_metadata.csv
input/receptors/3GT8_raw.pdb
input/receptors/EGFR_160-185.pdb
input/receptors/EGFR_170-200.pdb
```

### 3.2 현재 missing인 중요 파일

다음 파일들은 문서 또는 metadata에서 참조되지만 현재 local snapshot에는 없다.

```text
input/PPI/TH1 domain.pdb
input/PPI/beta_meander.pdb
input/PPI/phase1/partner_extended_beta_meander.pdb
input/PPI/phase1/receptor_3GT8_raw.pdb
input/PPI/phase1/docking_3GT8_raw_ext_beta_meander.pdb
input/ligands/173940_ligand.sdf
input/ligands/97806_ligand.sdf
input/ligands/VAX-C12_0_ligand.sdf
output/workflow_a
output/workflow_b
```

해석:

- 현재 로컬 clone은 full payload가 아니다.
- HPC 또는 원본 repo에 있는 `input/`과 `output/workflow_a` symlink/payload가 필요할 수 있다.
- 이 상태에서 즉시 rerun하거나 결과 검증을 수행하면 파일 없음으로 실패할 가능성이 높다.

### 3.3 receptor PDB 직접 확인 결과

현재 local `input/receptors/*.pdb`의 chain/range:

| File | Chain | Residue count | Range | Atom count |
|---|---:|---:|---:|---:|
| `3GT8_raw.pdb` | A | 309 | 699-1007 | 2485 |
| `3GT8_raw.pdb` | B | 307 | 701-1007 | 2470 |
| `EGFR_160-185.pdb` | X | 381 | 634-1014 | 12226 |
| `EGFR_170-200.pdb` | X | 381 | 634-1014 | 12226 |

해석:

- `3GT8_raw.pdb`는 local 기준 dimer chain A/B가 존재한다.
- `EGFR_160-185.pdb`와 `EGFR_170-200.pdb`는 local 기준 monomer-only chain X이다.
- dimer-only 논문 결론을 위해서는 MD states를 true dimer source로 교체하거나, 3GT8/3GOP template-superposition 방식으로 modeled dimer를 생성해야 한다.

---

## 4. 현재 config/code 상태

### 4.1 config 파일 존재 상태

`config/phase1/`에는 production INI가 3 states x 10 seeds = 30개 존재한다.

```text
3GT8_raw_seed0-9
EGFR_160-185_seed0-9
EGFR_170-200_seed0-9
```

### 4.2 하지만 현재 INI는 dimer-only 설계가 아니다

샘플 `config/phase1/phase1_prod_3GT8_raw_seed0.ini` 기준:

```text
# Construct: full_kinase_domain (NOT legacy dimer)
# Receptor is monomer (chain A only, ~309 residues)
# Partner is extended beta-meander (955-1006, not 960-1006)
input_pdb_name = output/workflow_a/phase2_ppi_docking/runtime_inputs/docking_3GT8_raw_ext_beta_meander.pdb
key_residues_b =
critical_residues_b =
non_binding_residues_b =
```

해석:

- 현재 production config는 10 seed로 확장되어 있지만 monomer baseline이다.
- 사용자의 최신 조건인 "EGFR는 dimer를 무조건 사용"과 맞지 않는다.
- MYO1D key/critical/non-binding residue annotation도 아직 비어 있다.
- input_pdb가 `output/workflow_a/.../runtime_inputs`를 가리키지만 현재 local에는 해당 output tree가 없다.

### 4.3 코드 상태

확인된 핵심 코드 상태:

- `run_production.py`
  - 주석과 target builder가 monomer-based Phase 1을 전제로 한다.
  - `PRODUCTION_N_SEEDS = 10`.
  - `_build_ppi_targets`에서 `mapping_csv=""`, `construct_type="full_kinase_domain"`로 들어간다.
  - chain restoration은 monomer target으로 판단해 skip한다.

- `egfr_pipeline/phase1/generate_configs.py`
  - monomer setup을 명시.
  - partner는 `955-1006`.
  - dimer +1000 offset residue handling은 현재 production config에 들어가지 않는다.

- `egfr_pipeline/phase1/launch_docking.py`
  - run metadata에 `partner_id = extended_beta_meander_955_1006`, `construct_type = full_kinase_domain`를 기록한다.
  - seed argument가 없을 경우 `range(5)`를 쓰는 경로가 있어, 10-seed 설계와 불일치 가능성이 있다.

- `egfr_pipeline/ppi/prepare_dimer_pdb.py`
  - dimer chain A/B를 chain A 하나로 병합하고 원래 chain B residue에 +1000 offset을 부여하는 기능이 있다.
  - dimer-only rerun에서 재사용할 핵심 유틸리티다.

### 4.4 Phase 2-4 코드 구조

Workflow B의 뼈대는 이미 있다.

- Phase 2:
  - `pocket_proposal.py`
  - `pocket_merge.py`
  - `patch_relationship.py`
  - `druggability_confidence.py`
  - `cross_state_alignment.py`
  - `phase3_export.py`

- Phase 3:
  - `job_construction.py`
  - `run_diverse_docking.py`
  - `pose_attribution.py`
  - `phase4_export.py`

- Phase 4:
  - `evidence_ingestion.py`
  - `score_framework.py`
  - `perturbation_scoring.py`
  - `mechanistic_classification.py`
  - `final_report.py`

부족한 점:

- dimer/membrane frame 정의 모듈이 없다.
- PPI pose의 membrane compatibility filter가 없다.
- pocket의 lower/lateral, dimer-interface accessibility gate가 없다.
- compound shortlist를 membrane/dimer/PPI mechanism으로 최종 통합하는 layer가 없다.

---

## 5. 기존 결과의 의미와 현재 재해석

### 5.1 문서상 완료된 기존 결과

`docs/CONTEXT.md`에 따르면 HPC `codex_ligand2`에서 과거 다음이 완료되었다.

- Workflow A Phase 1-7 완료.
- PPI 30 seeds, 3 states x 10 seeds = 600K models 완료.
- ATP site STRONG 차단 구현.
- Workflow B Phase 1-4 완료.
- fpocket 165 raw pockets -> 103 merged.
- P2Rank은 Java 11 requirement 때문에 미실행, fpocket 단독.
- focused Vina 168 jobs, 324 PDBQT output.
- Phase 4에서 207 candidates scored, 103 shortlisted.
- top candidates:
  - `3GT8_raw_PKT07`: rim, score 0.541, PPI 18.7 A, tier_1, state_robust.
  - `EGFR_170-200_PKT34`: allosteric, score 0.492, PPI 9.4 A, tier_1.

### 5.2 현재 새 기준에서의 해석

위 결과는 valuable historical evidence지만 final claim은 아니다.

이유:

- 현재 code/config는 monomer baseline으로 되어 있다.
- 현재 local output payload가 없다.
- pocket이 dimer receptor에서 발견된 것인지, monomer receptor에서 발견된 것인지 새 기준으로 확인해야 한다.
- pocket이 lower/lateral membrane-proximal side에 있는지 아직 정량화되어 있지 않다.
- dimer central interface 또는 membrane-occluded 위치인지 확인되어 있지 않다.
- `PKT07/PKT34`는 유지 후보가 아니라 "dimer/membrane-aware 재검증 대상"이다.

---

## 6. 확정된 설계 원칙

### 6.1 EGFR

- main analysis는 dimer-only.
- monomer는 technical/historical comparison으로만 사용.
- MD states가 monomer-only이면 dimer embedding 또는 true dimer source 확보가 필요.
- protomer mapping과 +1000 offset mapping은 반드시 기록.

### 6.2 MYO1D

- beta-sheet 8/9/12는 반드시 포함.
- active face:

```text
961-964, 968-972
```

- sheet12 support:

```text
993-997
```

- candidate main:

```text
955-1001
```

- conservative comparator:

```text
955-1006_tail_masked
```

- recommended annotation:

```ini
[Constraints]
key_residues_B = 961-964,968-972,993-997
key_residue_bonus_weight = 0.0

[ExperimentalData]
critical_residues_B = 961-964,968-972,993-997
non_binding_residues_B = 998-1006
```

### 6.3 Pocket

최종 pocket 후보 hard gate:

1. dimer receptor에서 발견.
2. ATP pocket 아님.
3. dimer-aware PPI patch와 orthosteric/rim/allosteric 관계.
4. lower/lateral membrane-compatible side.
5. central dimer interface 내부에 묻히지 않음.
6. state robust 또는 reproducible cryptic pocket.
7. focused docking pose convergence 있음.

### 6.4 Compound

- 기존 3개 ligand는 probe로만 사용.
- 실제 compound discovery에는 fragment/PPI-oriented library가 필요.
- 최종 shortlist는 affinity가 아니라 mechanism-aware ranking이어야 한다.

---

## 7. 가장 중요한 충돌점

### 충돌 1. 문서상 완료 결과 vs 현재 로컬 snapshot

문서에는 Workflow A/B 완료 결과가 있지만, 현재 local filesystem에는 `output/workflow_a`, `output/workflow_b`가 없다.

판단:

- 로컬은 lightweight/current-code snapshot이다.
- 과거 결과는 HPC 또는 원본 repo symlink에 있을 가능성이 있다.
- 재검증 전까지 결과 수치만으로 final manuscript를 확정하면 안 된다.

### 충돌 2. 사용자 최신 조건 vs 현재 config/code

사용자 최신 조건:

```text
EGFR는 dimer를 무조건 사용
```

현재 config/code:

```text
monomer full_kinase_domain baseline
```

판단:

- dimer-only rerun을 위해 config generator, run_production target builder, runtime input preparation을 다시 맞춰야 한다.

### 충돌 3. MYO1D tail control 설계 vs 현재 config

현재 config:

```text
key_residues_b =
critical_residues_b =
non_binding_residues_b =
```

새 설계:

```text
critical/key: 961-964,968-972,993-997
non-binding/tail: 998-1006
```

판단:

- 아직 tail-noise marking은 구현/config 반영되지 않았다.

### 충돌 4. 기존 pocket ranking vs membrane-aware criteria

기존 ranking은 PPI proximity, druggability, focused docking, state robustness 중심이다.

새 기준은 여기에 다음을 추가한다.

- dimer receptor origin.
- membrane lower/lateral geometry.
- dimer central-interface accessibility.
- compound approach feasibility.

판단:

- 기존 `PKT07/PKT34`는 재검증 전까지 "우선 검토 후보"이지 "최종 후보"가 아니다.

---

## 8. 바로 다음 작업 우선순위

### Priority 0. fresh-run input source 확정

기존 `output/workflow_a` 또는 `output/workflow_b` 복구는 필수 작업이 아니다. 새 연구는 fresh run이므로, 먼저 새 실행에 필요한 원천 입력을 확정해야 한다.

필수 fresh inputs:

```text
EGFR receptor states
MYO1D TH1/beta-meander source structure
compound or fragment library
reference structures for dimer/membrane geometry
```

현재 local에 `input/PPI/TH1 domain.pdb`, `input/ligands/*.sdf` 등은 없으므로, 새 run에서는 다음 중 하나를 선택해야 한다.

- 원본/HPC에서 필요한 input만 가져온다.
- AlphaFold/UniProt/RCSB 등에서 fresh source를 다시 준비한다.
- compound 후보는 기존 3개 probe ligand에 의존하지 않고 새 library를 정의한다.

### Priority 1. dimer input preparation 설계 확정

- `3GT8_raw`: existing dimer chain A/B 사용 가능.
- `EGFR_160-185`, `EGFR_170-200`: true dimer source가 있는지 확인.
- 없으면 3GT8/3GOP template-superposition dimer 생성.
- dimer QC table 작성.

### Priority 2. membrane frame 구현

신규 모듈 제안:

```text
egfr_pipeline/phase0/membrane_frame.py
```

출력:

```text
output/phase0_geometry/membrane_frame.json
output/phase0_geometry/dimer_geometry_qc.csv
```

### Priority 3. MYO1D construct pilot

비교:

```text
955-1001
955-1006_tail_masked
```

판정:

- sheet12 구조 보존.
- active-face contact 유지.
- tail-dominant false positive 감소.
- orientation pass 유지.

### Priority 4. dimer-only config generator 수정

필수 반영:

- dimer runtime input PDB.
- mapping CSV.
- `construct_type = dimer_offset` 또는 새 명명.
- `key_residues_b`, `critical_residues_b`, `non_binding_residues_b`.
- seed 0-9 동기화.

### Priority 5. membrane-aware pocket gate 추가

신규/수정 모듈:

```text
egfr_pipeline/phase1/membrane_pose_filter.py
egfr_pipeline/phase2/membrane_geometry.py
egfr_pipeline/phase2/phase3_export.py
```

### Priority 6. compound discovery funnel 확장

- existing 3 ligands = probe.
- 실제 compound 후보 = fragment/PPI-oriented library 필요.
- PAINS/reactive/basic property filter.
- focused docking + pose clustering + rescoring.

---

## 9. GPT Pro에 넘길 때 한 줄 결론

현재 프로젝트는 기존 결과를 복구해서 확정하는 단계가 아니라, 기존 문서에서 배운 문제점과 설계 의도를 바탕으로 **fresh dimer-only, membrane-aware, MYO1D-tail-controlled PPI-to-pocket-to-compound workflow**를 새로 설계해야 하는 단계다. 과거 결과는 참고용이며, 새 결론은 새 dimer/membrane-aware run에서만 낸다.
