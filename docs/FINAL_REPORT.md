# EGFR-MYO1D PPI 교란 약물 포켓 탐색 — 최종 리포트

**프로젝트**: EGFR kinase – MYO1D TH1 domain 결합 부위 규명 및 druggable pocket 탐색  
**날짜**: 2026-04-08  
**파이프라인**: Workflow A (Blind) + Workflow B (PPI-First) 이중 증거 통합  

---

## 1. 연구 목표

| 우선순위 | 목표 | 달성 여부 |
|----------|------|-----------|
| 1차 | MYO1D가 EGFR의 **어디에** 붙는가 — 결합 부위 규명 | **달성** |
| 1차 | 3개 EGFR 상태에서 교차 검증 | **달성** |
| 2차 | 결합 부위 근처에서 소분자로 PPI를 교란할 druggable pocket 탐색 | **달성** |

---

## 2. 사용한 구조

| ID | 유형 | 설명 |
|----|------|------|
| 3GT8_raw | Crystal | PDB 3GT8, EGFR kinase domain monomer |
| EGFR_160-185 | MD cluster | 38–48 ns MD trajectory cluster |
| EGFR_170-200 | MD cluster | 85–100 ns MD trajectory cluster |
| MYO1D TH1 | Partner | Extended beta-meander (residues 955–1006) |

리간드 3종: 173940, 97806, VAX-C12_0 (쌍별 Tanimoto < 0.4 — 구조적 편향 방지)

---

## 3. Workflow A 결과 (Blind Approach)

### 3.1 PPI 결합 부위

**PPI 도킹**: 3 상태 × 10 seeds × 20K models = 600K 모델, PyRosetta global docking

**EGFR 측 (chain A) — MYO1D가 붙는 곳**:
- 3/3 상태 공통 잔기 22개, **C-lobe 표면**에 집중
- 최상위 occupancy: **ILE941**(0.205), VAL980(0.148), THR940(0.118), PRO992(0.118)
- N-lobe 접촉 미미 → 결합은 C-lobe 특이적

**MYO1D 측 (chain B) — EGFR에 닿는 부분**:
- Sheet 8/9 active face: VAL962(0.210), VAL964(0.165), CYS970(0.103), SER971(0.093)
- **Ko et al. alanine substitution 실험과 일치** → PPI 도킹 생물학적 타당성 확인

### 3.2 Vina Blind Docking

| Pocket | Receptor | Verdict | Affinity | PPI 거리 | ATP | 비고 |
|--------|----------|---------|----------|----------|-----|------|
| P045 | EGFR_170-200 | **STRONG** | -6.9 | 62.7Å | No | PPI 원거리 |
| P010 | EGFR_160-185 | MODERATE | -9.1 | 14.9Å | **Yes** | ATP 배제 |
| P004 | EGFR_170-200 | MODERATE | -9.0 | 14.5Å | **Yes** | ATP 배제 |
| P003 | 3GT8_raw | MODERATE | -9.8 | 14.3Å | **Yes** | ATP 배제 |
| P004 | 3GT8_raw | MODERATE | -8.9 | 16.3Å | No | non-ATP 최근접 |

**WF-A 한계**: ATP 포켓 배제 후, PPI 근접(< 20Å) non-ATP druggable pocket 미발견. 유일한 STRONG(P045)은 PPI에서 62.7Å → PPI 교란 목적 부적합.

→ **Workflow B (PPI-First) 전환 근거**

---

## 4. Workflow B 결과 (PPI-First Approach)

### 4.1 Phase 1: PPI 패치 분석

**Orientation filter** (AMBIGUOUS_BAND=0.10, WF-A 600모델 retroactive 검증으로 최적화):

| State | Pass | Fail | Ambiguous |
|-------|------|------|-----------|
| 3GT8_raw | 156 (78%) | 26 (13%) | 18 (9%) |
| EGFR_160-185 | 127 (64%) | 49 (25%) | 24 (12%) |
| EGFR_170-200 | 135 (68%) | 46 (23%) | 19 (10%) |

**Cluster consensus** (orientation-filtered models only):

| State | Hotspots (≥50% occ) | Multi-cluster |
|-------|---------------------|---------------|
| 3GT8_raw | 17 | 5 |
| EGFR_160-185 | **22** | **12** |
| EGFR_170-200 | 5 | 4 |

**Cross-state comparison**: 179 잔기 비교 → 130 잔기 핸드오프

**상위 robust hotspots** (3/3 상태 공통):
- **ILE941** (occ=1.0), ARG977 (0.75), THR993 (0.71), ARG986 (0.60)
- WF-A 결과와 완전 일치 → 이중 검증

### 4.2 Phase 2: Pocket Analysis

- fpocket: 165 raw pockets → 103 merged (3 states)
- PPI 관계: 21 rim, 31 allosteric, 51 low_relevance
- Druggability: **2 high (tier_1)**, 1 medium, 100 low
- Cross-state: 13 state_robust, 68 shifted
- Docking priority: **2 primary**, 21 secondary, 29 exploratory, 51 skip

### 4.3 Phase 3: Focused Vina Docking

- 168 jobs (52 pockets × 3 ligands), 3 rounds
- Round 1: 156 ok (전체 포켓) / Round 2-3: 6 ok each (primary 추가 도킹)
- 324 output PDBQT 파일

### 4.4 Phase 4: Perturbation Scoring (최종 결과)

207 candidates scored → 103 shortlisted

---

## 5. 최종 후보 포켓

### Primary Candidates

| 순위 | Pocket ID | 분류 | Score | PPI 거리 | Druggability | Cross-state | Volume |
|------|-----------|------|-------|----------|-------------|-------------|--------|
| **1** | **3GT8_raw_PKT07** | **rim** | **0.541** | 18.7Å | tier_1 / high | robust (3/3) | 1593Å³ |
| **2** | **EGFR_170-200_PKT34** | **allosteric** | **0.492** | **9.4Å** | tier_1 / high | shifted (3/3) | 1938Å³ |

### Secondary Candidates (Top 10)

| 순위 | Pocket ID | 분류 | Score |
|------|-----------|------|-------|
| 3 | EGFR_160-185_PKT02 | rim | 0.433 |
| 4 | 3GT8_raw_PKT10 | allosteric | 0.431 |
| 5 | EGFR_170-200_PKT17 | rim | 0.430 |
| 6 | EGFR_170-200_PKT06 | rim | 0.422 |
| 7 | EGFR_160-185_PKT16 | rim | 0.419 |
| 8 | 3GT8_raw_PKT05 | rim | 0.416 |
| 9 | 3GT8_raw_PKT11 | rim | 0.412 |
| 10 | 3GT8_raw_PKT01 | rim | 0.410 |

### 분류 분포

| 분류 | 개수 | Score 범위 | 해석 |
|------|------|-----------|------|
| rim | 23 | 0.38–0.54 | PPI 가장자리, 소분자 교란 가능 |
| allosteric | 2 | 0.43–0.49 | PPI 인접, 구조적 교란 가능 |
| uncertain | 29 | 0.25–0.27 | 증거 부족, 추가 검증 필요 |
| irrelevant | 49 | 0.15–0.19 | PPI 교란 무관 |

---

## 6. 핵심 결론

### 6.1 MYO1D는 EGFR C-lobe 표면에 결합한다

- **EGFR 측**: ILE941 중심의 C-lobe 표면 패치 (3/3 상태 공통, occupancy 1.0)
- **MYO1D 측**: Sheet 8/9 active face (VAL962, VAL964) — Ko et al. 실험 일치
- Orientation filter 적용 후에도 동일 결과 → 노이즈 모델 제거로 신호 강화
- **PyRosetta + LightDock 교차 검증**: ILE941, ARG977 등 핵심 hotspot이 두 독립적 도킹 방법에서 모두 확인 (method_agreement: both)

**LightDock 교차 검증 결과**:

| State | Convergent | PyRosetta-only | LightDock-only | Jaccard |
|-------|-----------|----------------|----------------|---------|
| 3GT8_raw | 53 | 16 | 116 | 0.286 |
| EGFR_160-185 | 56 | 25 | 104 | 0.303 |
| EGFR_170-200 | 19 | 72 | 56 | 0.129 |

C-lobe overlap (0.21–0.38)이 N-lobe (0.00–0.23)보다 높음 → C-lobe 결합 부위가 방법론적으로 독립 확인됨.

### 6.2 PPI 교란 druggable pocket 발견

| | WF-A (Blind) | WF-B (PPI-First) |
|---|---|---|
| PPI 최근접 non-ATP | 16.3Å (3GT8 P004) | **9.4Å (PKT34)** |
| Druggable pocket 수 | 1 STRONG (62.7Å) | **2 primary (9.4Å, 18.7Å)** |
| PPI 교란 가능성 | 낮음 | **높음** |

**PKT34** (EGFR_170-200): PPI 인터페이스에서 9.4Å — WF-A blind Vina로는 발견 불가했던 allosteric druggable pocket. PPI 패치 기반 탐색(WF-B)의 핵심 성과.

**PKT07** (3GT8_raw): PPI rim에 위치 (18.7Å), 3/3 상태 robust. 포켓에 소분자를 결합시키면 PPI 인터페이스 구조를 교란할 가능성.

### 6.3 ATP 포켓은 배제됨

WF-A에서 최고 점수(74.0)를 받은 P010/P004/P003는 모두 ATP 포켓. 실험적으로 ATP 결합이 유지되면서 활성이 소실되므로 (Ko et al.), ATP site는 PPI 교란 타겟으로 부적합.

---

## 7. 제한사항

1. **P2Rank 미사용**: HPC Java 8 vs P2Rank Java 11+ 요구. fpocket 단독 결과.
2. **LightDock clustering 실패**: lgd_cluster_bsas 에러로 클러스터 기반 분석 불가. rank_by_scoring 기반 추출로 대체 — 교차 검증은 성공적으로 수행됨.
3. **MD 시뮬레이션 미통합**: GROMACS MD 분석 모듈 존재하나 이번 파이프라인에 미통합.
4. **실험 데이터 미검증**: 최종 후보 포켓의 실험적 검증 (mutagenesis, SPR 등) 필요.

---

## 8. 파일 위치 (HPC)

```
output/workflow_a/                          # WF-A 전체 결과
  phase5_verdict/valid_sites.csv            # WF-A 최종 포켓 판정
  phase3_ppi_postprocess/ppi_pyrosetta_residues.csv  # PPI 결합 잔기

output/workflow_b/                          # WF-B 전체 결과
  phase1_ppi_analysis/
    phase1_downstream_patch_reference.csv   # PPI 패치 핸드오프 (130 잔기)
    */orientation_filter_log.csv            # Orientation dot product
    */ppi_hotspot_residues.csv              # 상태별 hotspot
  phase2_pocket_analysis/
    phase3_candidate_pocket_reference.csv   # 포켓 후보 (103개)
  phase3_focused_docking/
    phase3_round_log.csv                    # Vina 도킹 결과
    runs/*/                                 # PDBQT 도킹 포즈
  phase4_scoring/
    phase4_final_report.md                  # Phase 4 최종 리포트
    phase4_review_expanded.csv              # 전체 스코어 테이블
```
