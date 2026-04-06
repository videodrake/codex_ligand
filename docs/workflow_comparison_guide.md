# Workflow Comparison Guide

Workflow A (Blind)와 Workflow B (PPI-First)의 구조, 스코어링, 사용 시나리오를 비교한다.

## 1. 개요

| 항목 | Workflow A (Blind) | Workflow B (PPI-First) |
|------|-------------------|----------------------|
| 전략 | 넓은 탐색 → 증거 통합 | PPI 기반 정밀 탐색 |
| Phase 수 | 7 (도킹 4 + 판정 3) | 4 |
| 스코어링 | 3축 적응형 (100점) | 4축 가중형 (0.0-1.0) |
| 의존 구조 | Phase 1/2 병렬 실행 가능 | 엄격한 순차 실행 |
| 선행 조건 | 없음 | Workflow A Phase 2 (PPI 도킹) 완료 |

## 2. Phase 대응 관계

| Phase | Workflow A | 디렉토리 | Workflow B | 디렉토리 |
|-------|-----------|----------|-----------|----------|
| 1 | Vina Blind 도킹 | `vina/` | PPI 분석 (TG 1.0-1.6) | `phase1/` |
| 2 | PPI Global Blind | `ppi/`, `pyrosetta_docking/` | Pocket Analysis (TG 2.0-2.7) | `phase2/` |
| 3 | PPI Postprocess | `ppi/` | Focused Vina (TG 3.0-3.6) | `phase3/` |
| 4 | Vina Postprocess | `vina/` | Perturbation Scoring (TG 4.0-4.6) | `phase4/` |
| 5 | Verdict | `verdict.py` | -- | |
| 6 | Report | `report.py` | -- | |
| 7 | Validate | `validate.py` | -- | |

**공유 엔진**: `egfr_pipeline/pyrosetta_docking/` — Workflow A Phase 2와 Workflow B Phase 1이 동일한 PyRosetta 도킹 엔진을 사용한다.

## 3. Workflow A 상세

### 의존 그래프

```
Phase 1 (Vina Blind) ──┐
                        ├→ Phase 4 (Vina Post) ──┐
Phase 2 (PPI Blind)  ──┤                         ├→ Phase 5 (Verdict) → Phase 6 (Report) → Phase 7 (Validate)
                        └→ Phase 3 (PPI Post)  ──┘
```

Phase 1과 2는 **병렬 실행 가능**하다. 둘 사이에 데이터 의존이 없다.

### Phase별 입출력

**Phase 1: Vina Blind 도킹**
- 70 A+ blind box로 EGFR 전체 표면 탐색
- 3 수용체 상태 x 3 리간드, exhaustiveness=384
- 출력: raw poses, 친화도, centroid 좌표

**Phase 2: PPI Global Blind**
- PyRosetta 글로벌 도킹 (3 상태 x 5 seeds x 20K models = 300K)
- Orientation filter로 Ko et al. active-face 방향성 검증
- 출력: `final_ranking.csv`, interface residue tables

**Phase 3: PPI Postprocess**
- Chain 복원 (merged dimer → 원래 chain ID)
- PPI 잔기 추출 → project-level CSV
- 출력: `ppi_pyrosetta_residues.csv`, `ppi_interface_patch_table.csv`

**Phase 4: Vina Postprocess**
- Pose 파싱 → contact 분석 → DBSCAN 클러스터링 → cross-receptor 비교
- 출력: `vina_pocket_table.csv`, `vina_pocket_comparison.csv`

**Phase 5: Verdict (3축 스코어링)**
- 아래 "스코어링 비교" 섹션 참조
- 출력: `valid_sites.csv`, `cross_method_agreement.csv`

**Phase 6-7: Report / Validate**
- 사람이 읽을 수 있는 보고서 생성 + 회귀 검증

### 3축 스코어링 체계

| 축 | 이름 | 배점 (PPI 있음) | 배점 (PPI 없음) |
|----|------|----------------|----------------|
| A1 | Vina Quality | 50 | 60 |
| A2 | PPI Spatial | 20 | 0 |
| A3 | Cross-Receptor | 30 | 40 |
| **합계** | | **100** | **100** |

- **STRONG**: >= 55점, **MODERATE**: 25-54점, **WEAK**: < 25점
- **적응형**: PPI 데이터 부재 시 A2를 제거하고 A1/A3에 가중치 재배분
- ATP site가 STRONG 판정을 받으면 안 된다 (절대 규칙 #2)

## 4. Workflow B 상세

### 의존 그래프

```
Phase 1 → Phase 2 → Phase 3 → Phase 4
```

**엄격한 순차 실행**. 각 Phase는 이전 Phase의 핸드오프 CSV를 소비한다.

### 핸드오프 CSV 계약

| From | To | CSV | 생성 TG |
|------|----|-----|---------|
| Phase 1 | 2 | `phase1_downstream_patch_reference.csv` | TG 1.6 |
| Phase 2 | 3 | `phase3_candidate_pocket_reference.csv` | TG 2.6 |
| Phase 3 | 4 | `phase4_docking_evidence_reference.csv` | TG 3.6 |

각 Phase는 입력 CSV를 `_validate_adv_handoff()`로 검증한 후 진행한다.

### Phase별 입출력

**Phase 1: PPI 분석**
- Full kinase domain (699-1007) + extended beta-meander (955-1006)로 입력 준비
- PyRosetta 도킹 + LightDock 보조 검증
- Orientation filter + consensus scoring
- **중단 조건**: Ko et al. sheet 8/9 활성면 잔기 3개 미만 → 리뷰 필요

**Phase 2: Pocket Analysis**
- Phase 1 패치 참조를 수신하여 후보 포켓 제안 (fpocket)
- 중복 포켓 병합, druggability 등급 분류 (TIER-1 ~ TIER-4)
- 교차 상태 포켓 정렬

**Phase 3: Focused Vina**
- Phase 2 포켓별 targeted box로 집중 도킹
- Budget-aware 실행: 포켓별 라운드 기반 + 포화도 감지
- Pose attribution: 각 포즈를 원래 포켓에 추적

**Phase 4: Perturbation Scoring (4축)**
- 아래 "스코어링 비교" 섹션 참조
- 기계적 분류: orthosteric (직접) / rim (가장자리) / allosteric (원격)

### 4축 스코어링 체계

| 축 | 이름 | 가중치 | 주요 하위 지표 |
|----|------|--------|---------------|
| A1 | PPI Interface Confidence | 30% | hotspot overlap, confidence, method agreement |
| A2 | Druggability | 25% | 전체 tier, pocket volume, fpocket/P2Rank consensus |
| A3 | Perturbation Relevance | 30% | relationship class, ligand support, pose count |
| A4 | State Robustness | 15% | state class, n_states_matched, coverage |

- A1 + A3 = 60% → MYO1D와 무관한 포켓이 상위에 오르는 것을 방지
- 축 가중치 변경은 사람 승인 필수 (절대 규칙 #8)

## 5. 스코어링 체계 비교

| 기준 | Workflow A (3축) | Workflow B (4축) |
|------|-----------------|-----------------|
| 척도 | 0-100 절대 점수 | 0.0-1.0 가중 합산 |
| PPI 핵심 비중 | 20% (A2) | 60% (A1+A3) |
| Vina 비중 | 50% (A1) | A3에 간접 포함 |
| 교차 상태 | 독립 축 30% (A3) | 15% (A4) |
| 적응형 | PPI 부재 시 재배분 | 적응 없음 (PPI 필수) |
| 기계적 분류 | 없음 | orthosteric/rim/allosteric |
| 결과 직접 비교 | -- | **불가** (척도와 축 구성이 다름) |

## 6. 사용 시나리오

### Workflow A 적합
- 초기 탐색 단계: EGFR 전체 표면에서 후보 포켓을 넓게 탐색
- PPI 데이터가 예비적이거나 불확실할 때
- 빠른 baseline 결과가 필요할 때 (Phase 1/2 병렬 가능)

### Workflow B 적합
- Workflow A Phase 2 (PPI 도킹)가 완료된 후
- PPI 패치 기반으로 포켓을 정밀하게 좁히고 싶을 때
- 기계적 분류(orthosteric vs allosteric)가 필요할 때
- Budget-aware 집중 도킹으로 컴퓨팅 자원을 효율적으로 사용하고 싶을 때

### 주의사항
- **같은 Phase 번호가 다른 모듈을 가리킨다** (절대 규칙 #3). "Phase 2 수정" 요청 시 반드시 Workflow A/B를 확인한다.
- 두 Workflow의 점수를 직접 비교할 수 없다. 척도, 축 구성, 탐색 범위가 모두 다르다.
- Workflow B는 Workflow A를 대체하는 것이 아니라 **보완**한다.

## 7. 디렉토리-워크플로우 매핑

```
egfr_pipeline/vina/              → Workflow A (Phase 1 + 4)
egfr_pipeline/ppi/               → Workflow A (Phase 2 + 3)
egfr_pipeline/phase1/            → Workflow B Phase 1
egfr_pipeline/phase2/            → Workflow B Phase 2
egfr_pipeline/phase3/            → Workflow B Phase 3
egfr_pipeline/phase4/            → Workflow B Phase 4
egfr_pipeline/pyrosetta_docking/ → 공유 (WA Phase 2 + WB Phase 1)
verdict.py, report.py, validate.py → Workflow A Phase 5-7
```
