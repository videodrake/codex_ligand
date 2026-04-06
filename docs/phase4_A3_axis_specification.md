# Phase 4 A3 Axis Specification — Perturbation Relevance

Workflow B Phase 4 스코어링 체계의 A3 축(Perturbation Relevance) 기술 명세.

## 1. 축 정의

| 항목 | 값 |
|------|-----|
| Axis ID | `A3_perturbation_relevance` |
| 가중치 | 0.30 (전체 perturbation score의 30%) |
| 범위 | 0.0 - 1.0 |
| 의미 | 포켓이 MYO1D 결합을 교란할 수 있는 기계적 가능성 |
| 핵심 함수 | `score_framework.py:224` `compute_a3_perturbation_relevance()` |

A1(PPI Interface Confidence, 30%) + A3 = **60%** → MYO1D와 무관한 고친화도 포켓이 상위에 올라오는 것을 방지한다.

## 2. 스코어링 공식

```
A3 = 0.45 * relationship_score
   + 0.25 * hotspot_overlap_norm
   + 0.20 * ligand_support_score
   + 0.10 * pose_support_norm
```

### 2.1 Relationship Class (45%)

PPI 패치와의 기계적 근접도. Phase 2 pocket-patch relationship에서 결정된다.

| 분류 | 점수 | 정의 |
|------|------|------|
| `orthosteric_candidate` | 1.0 | MYO1D 부착 패치와 직접 중첩 |
| `rim_candidate` | 0.7 | PPI 계면 가장자리; 간접적 결합 약화 가능 |
| `allosteric_candidate` | 0.5 | 공간적으로 원격이나 kinase domain 내; 구조 변화 유도 가능 |
| `low_relevance_candidate` | 0.1 | 약물 결합 가능하나 MYO1D 부착과 기계적 연결 없음 |

**분류 기준** (`mechanistic_classification.py:41-98`):

| 기계적 분류 | 조건 |
|------------|------|
| orthosteric_disruptor | relationship=orthosteric AND hotspot_overlap >= 2 |
| interface_rim_modulator | relationship=rim AND hotspot_overlap >= 1 |
| allosteric_modulator | relationship=allosteric AND druggability_tier in (tier_1, tier_2) |
| ligandable_but_ppi_irrelevant | relationship=low_relevance AND hotspot_overlap = 0 |
| uncertain_mechanism | 위 조건 모두 불일치 시 기본값 |

**Edge case**: orthosteric이지만 overlap=1인 경우 → `uncertain_mechanism_candidate`로 분류. 신뢰도 있는 orthosteric 판정에는 overlap >= 2가 필요하다.

### 2.2 Hotspot Overlap (25%)

Phase 1 PPI hotspot 잔기와 포켓 잔기의 중첩 정도.

```python
overlap_norm = min(hotspot_overlap_count / n_hotspot_residues, 1.0)
```

- `n_hotspot_residues`: Phase 1에서 `is_hotspot_any_state = True`인 잔기 수 (최소 1로 처리)
- 예: hotspot 7개 중 4개 중첩 → `overlap_norm = 4/7 = 0.571`

### 2.3 Ligand Support Strength (20%)

Phase 3 도킹 결과에 기반한 리간드 지지 수준.

| 분류 | 점수 | 조건 (실행 후) | 조건 (실행 전) |
|------|------|---------------|---------------|
| `strong` | 1.0 | >= 5 poses AND best affinity <= -6.0 kcal/mol | -- |
| `moderate` | 0.6 | >= 2 poses AND best affinity <= -5.0 kcal/mol | -- |
| `weak` | 0.3 | > 0 poses (위 기준 미달) | -- |
| `none` | 0.0 | 0 poses | -- |
| `pending_strong` | 0.8 | -- | >= 5 poses (친화도 미확인) |
| `pending_moderate` | 0.5 | -- | >= 2 poses (친화도 미확인) |
| `pending_weak` | 0.2 | -- | < 2 poses (친화도 미확인) |

분류 로직: `phase3/phase4_export.py:187` `_classify_support()`

**주의**: `pending_*` 상태는 Vina 실행 전 잠정 평가이다. 실행 후 재점수가 필수이다.

### 2.4 Pose Support Count (10%)

Phase 3에서 해당 포켓에 귀속된 총 포즈 수.

```python
pose_norm = min(pose_support_count / 10.0, 1.0)
```

10개 이상은 1.0으로 cap된다.

## 3. 친화도 지배 방지 (Affinity Cap)

생물학적 무관 포켓이 높은 친화도로만 상위에 오르는 것을 방지한다 (`perturbation_scoring.py:133-173`).

**대상**: `ligandable_but_ppi_irrelevant_candidate`, `uncertain_mechanism_candidate`

```python
affinity_cap_fraction = 0.35
max_affinity = (0.35 / 0.65) * biological_component  # A1 + A4
capped_affinity = min(actual_affinity_component, max_affinity)  # A2 + A3
final_score = biological_component + capped_affinity
```

효과: 생물학적 관련성 없는 고친화도 포켓은 최종 점수의 35%를 초과할 수 없다.

## 4. 기계적 분류 신뢰도

Orthosteric 분류의 신뢰도 등급 (`mechanistic_classification.py:215-226`):

| 신뢰도 | 조건 |
|--------|------|
| high | overlap_fraction >= 0.40 AND (strong OR pending_strong OR moderate OR pending_moderate) |
| medium | overlap_fraction >= 0.25 OR moderate support |
| low | 위 조건 불일치 OR 도킹 데이터 없음 |

## 5. 데이터 흐름

```
Phase 1 (PPI 분석)
  └→ phase1_downstream_patch_reference.csv
       └→ n_hotspot_residues, is_hotspot_any_state

Phase 2 (Pocket Analysis)
  └→ pocket_patch_relationship.csv
       └→ hotspot_overlap_count, hotspot_overlap_fraction, relationship_class

Phase 3 (Focused Docking)
  └→ phase4_docking_evidence_reference.csv
       └→ pose_support_count, ligand_support_strength

Phase 4 (Evidence Ingestion)
  └→ phase4_evidence_normalized.csv (병합된 증거)
       └→ compute_a3_perturbation_relevance() → A3 점수
```

핵심 함수:
- `evidence_ingestion.py:265` `load_phase3_evidence()`
- `evidence_ingestion.py:368` `build_normalized_evidence()`
- `mechanistic_classification.py:129` `classify_candidate()`
- `score_framework.py:224` `compute_a3_perturbation_relevance()`

## 6. 출력 컬럼

`phase4_axis_scores.csv`에 기록되는 A3 관련 컬럼:

| 컬럼 | 타입 | 설명 |
|------|------|------|
| `A3_perturbation_relevance` | float | 최종 A3 점수 (소수점 4자리) |
| `relationship_class` | str | 기계적 분류 (categorical) |
| `hotspot_overlap_count` | int | 중첩 hotspot 잔기 수 |
| `hotspot_overlap_fraction` | float | 정규화된 중첩 비율 |
| `ligand_support_strength` | str | 리간드 지지 수준 (categorical) |
| `pose_support_count` | int | 귀속 포즈 수 |

## 7. Edge Cases

| 상황 | 처리 |
|------|------|
| hotspot 잔기 0개 | `n_hotspot_residues`를 1로 처리 (0 나누기 방지) |
| orthosteric + overlap=1 | `uncertain_mechanism_candidate`로 분류 |
| 도킹 증거 없음 (pending_*) | A3 정상 계산, 신뢰도 = low |
| 포즈 0개 | `ligand_support_strength = "none"`, 기여도 0.0 |
| Phase 2 포켓이 Phase 3에 미진입 | 별도 행으로 처리, `ligand_id=""`, support="none" |
| `pose_support_count > 10` | 1.0으로 cap |
| 상태 특이적 포켓 | A4에서 flag, A3에서는 페널티 없음 (데이터 한계) |

## 8. A3와 다른 축의 관계

| 축 | 가중치 | A3와의 관계 |
|----|--------|------------|
| A1 PPI Interface | 30% | A1은 PPI 패치 품질, A3는 도킹 증거 결합. 함께 60% |
| A2 Druggability | 25% | A2는 실현 가능성, A3는 기계적 관련성. 무관 포켓은 A2+A3 cap 적용 |
| A4 State Robustness | 15% | A4가 상태 견고성 판단. A3의 상태 특이성은 A4에 위임 |

## 9. 변경 시 주의사항

- A3 하위 지표 가중치(0.45/0.25/0.20/0.10) 변경은 **사람 승인 필수** (절대 규칙 #8)
- `RELATIONSHIP_SCORES` 또는 `SUPPORT_SCORES` 매핑 변경 시 `mechanistic_classification.py`의 분류 로직과 일관성 확인 필요
- 리간드 지지 임계값(-6.0/-5.0 kcal/mol) 변경 시 `phase3/phase4_export.py`와 `score_framework.py` 동시 수정
