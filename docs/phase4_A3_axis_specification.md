# Phase 4 A3 축 (Perturbation Relevance) 스펙

> AC-4.5: Phase 4의 A3 축 계산 로직을 코드에서 추출하여 문서화

## 1. 입력 데이터

A3는 Phase 2와 Phase 3의 증거를 결합한다:

| 입력 | 출처 | 설명 |
|------|------|------|
| `relationship_class` | Phase 2 `pocket_patch_relationship.csv` | 포켓-패치 공간 관계 (orthosteric/rim/allosteric/low_relevance) |
| `hotspot_overlap_count` | Phase 2 `pocket_patch_relationship.csv` | 포켓과 PPI hotspot 잔기의 직접 겹침 수 |
| `n_hotspot_residues` | Phase 1 `phase1_downstream_patch_reference.csv` | 전체 hotspot 잔기 수 (정규화 분모) |
| `ligand_support_strength` | Phase 3 `phase4_docking_evidence_reference.csv` | 도킹 증거 강도 (strong/moderate/weak/none/pending_*) |
| `pose_support_count` | Phase 3 `phase4_docking_evidence_reference.csv` | 지지 포즈 수 (정수) |

입력 경로: `phase4_evidence_normalized.csv` (TG 4.0 evidence_ingestion이 위 소스를 병합)

## 2. 계산 로직

**소스:** `egfr_pipeline/phase4/score_framework.py:221` — `compute_a3_perturbation_relevance()`

```
A3 = 0.45 × relationship_score
   + 0.25 × overlap_norm
   + 0.20 × support_score
   + 0.10 × pose_norm
```

### 서브메트릭 변환

#### 2.1 relationship_score (45%)
범주형 → 수치 매핑:

| relationship_class | 점수 | 근거 |
|---|---|---|
| orthosteric_candidate | 1.0 | PPI 패치 직접 겹침 |
| rim_candidate | 0.7 | PPI 인터페이스 경계 |
| allosteric_candidate | 0.5 | 공간적으로 인접하나 직접 겹침 없음 |
| low_relevance_candidate | 0.1 | PPI와 무관 |

#### 2.2 overlap_norm (25%)
정규화: `min(hotspot_overlap_count / n_hotspot_residues, 1.0)`
- 0개 겹침 → 0.0, 전체 겹침 → 1.0

#### 2.3 support_score (20%)
범주형 → 수치 매핑:

| ligand_support_strength | 점수 | 비고 |
|---|---|---|
| strong | 1.0 | Phase 3 도킹 완료, 강한 증거 |
| pending_strong | 0.8 | 도킹 전, 구조적 예측 기반 |
| moderate | 0.6 | |
| pending_moderate | 0.5 | |
| weak | 0.3 | |
| pending_weak | 0.2 | |
| none | 0.0 | 도킹 증거 없음 |

#### 2.4 pose_norm (10%)
정규화: `min(pose_support_count / 10.0, 1.0)`
- 10개 이상 포즈 → 1.0으로 캡

## 3. 가중치

Phase 4 전체 스코어링에서 A3의 위치:

| 축 | 가중치 | 역할 |
|---|---|---|
| A1 PPI Interface Confidence | 30% | PPI 패치 증거 품질 |
| A2 Druggability Confidence | 25% | 소분자 결합 가능성 |
| **A3 Perturbation Relevance** | **30%** | **MYO1D 교란 메커니즘** |
| A4 State Robustness | 15% | 구조 상태 간 일관성 |

**A1 + A3 합산 = 60%** — MYO1D 교란 관련 축이 전체 점수의 과반.

최종 점수: `perturbation_score = Σ(axis_score_i × weight_i)`

## 4. 출력 범위

| 항목 | 범위 | 설명 |
|------|------|------|
| A3 축 점수 | 0.0 – 1.0 | 서브메트릭 가중합, clamp 적용 |
| A3 가중 기여 | 0.0 – 0.30 | A3 × 0.30 (전체 점수 대비) |
| perturbation_score | 0.0 – 1.0 | 4축 가중합 |

### Affinity Domination 방지 (TG 4.3)

`ligandable_but_ppi_irrelevant` 또는 `uncertain_mechanism` 사이트에 대해:
- A2 + A3 기여분이 전체 점수의 **35%**를 초과하지 못하도록 캡 적용
- 수식: `max_affinity = (0.35 / 0.65) × bio_component` (bio = A1 + A4 기여)
- 생물학적 무관 사이트에서 높은 affinity만으로 순위가 올라가는 것을 방지

## 5. A3 해석 가이드

| A3 점수 | 해석 |
|---------|------|
| ≥ 0.80 | Orthosteric + 강한 도킹 증거 → 직접 PPI 차단 후보 |
| 0.50–0.79 | Rim/allosteric + 중간 증거 → 간접 교란 가능성 |
| 0.20–0.49 | 약한 공간적 관련성 또는 증거 부족 |
| < 0.20 | PPI 교란 메커니즘 불명확 |

## 6. 코드 참조

| 파일 | 함수/상수 | 역할 |
|------|-----------|------|
| `phase4/score_framework.py:221` | `compute_a3_perturbation_relevance()` | A3 계산 |
| `phase4/score_framework.py:203-208` | `RELATIONSHIP_SCORES` | 관계 클래스 매핑 |
| `phase4/score_framework.py:210-218` | `SUPPORT_SCORES` | 리간드 지지 매핑 |
| `phase4/perturbation_scoring.py:39-44` | `DEFAULT_WEIGHTS` | 축 가중치 |
| `phase4/perturbation_scoring.py:106` | `compute_perturbation_score()` | 4축 합산 |
| `phase4/perturbation_scoring.py:130` | `apply_affinity_cap()` | Affinity 캡 |
| `phase4/mechanistic_classification.py:38` | `MECHANISTIC_CLASSES` | 기계적 분류 |
| `phase4/evidence_ingestion.py` | `run_evidence_ingestion()` | 증거 수집 |
