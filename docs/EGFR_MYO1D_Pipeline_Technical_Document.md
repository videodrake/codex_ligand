# EGFR-MYO1D 다중 방법론 결합 부위 탐색 파이프라인

## 기술 문서 (Technical Documentation)

> **버전**: 2.0 (2026-03-10)
> **목적**: EGFR kinase C-lobe 표면의 MYO1D 결합 후보 부위를 AutoDock Vina 리간드 도킹과
> PyRosetta 단백질-단백질 도킹의 교차 검증으로 탐색한다.

---

## 목차

1. [연구 배경 및 전략](#1-연구-배경-및-전략)
2. [시스템 구성](#2-시스템-구성)
3. [Phase 1: AutoDock Vina Blind Docking](#3-phase-1-autodock-vina-blind-docking)
4. [Phase 2: PyRosetta Global Blind PPI Docking](#4-phase-2-pyrosetta-global-blind-ppi-docking)
5. [Phase 3: PPI 후처리 — 체인 복원 및 잔기 추출](#5-phase-3-ppi-후처리)
6. [Phase 4: Vina 후처리 파이프라인](#6-phase-4-vina-후처리-파이프라인)
7. [Phase 5: 3축 Verdict 통합 평가](#7-phase-5-3축-verdict-통합-평가)
8. [Phase 6–7: 리포트 생성 및 검증](#8-phase-6-7-리포트-생성-및-검증)
9. [통계적 신뢰성 분석](#9-통계적-신뢰성-분석)
10. [결과 해석 가이드](#10-결과-해석-가이드)
11. [파라미터 요약 및 근거](#11-파라미터-요약-및-근거)
12. [제한사항 및 향후 과제](#12-제한사항-및-향후-과제)
13. [출력 파일 명세](#13-출력-파일-명세)

---

## 1. 연구 배경 및 전략

### 1.1 생물학적 맥락

EGFR (Epidermal Growth Factor Receptor)의 kinase domain C-lobe는 단백질-단백질
상호작용(PPI)의 핵심 매개 영역이다. MYO1D (Myosin-1D)와의 상호작용 부위를 규명하는
것이 본 연구의 목표이며, 이를 위해 **두 가지 독립적인 계산적 접근**을 교차 검증한다:

1. **소분자 프로브 도킹 (AutoDock Vina)**: 3종의 소분자 리간드를 C-lobe 표면에
   blind docking하여 "druggable pocket"을 탐색한다. 소분자가 결합하는 표면 포켓은
   단백질-단백질 상호작용 인터페이스와 공간적으로 겹칠 가능성이 높다.

2. **단백질-단백질 도킹 (PyRosetta)**: MYO1D 파트너 도메인(TH1, β-meander)을
   EGFR C-lobe에 global blind docking하여 직접적인 PPI 결합 인터페이스를 예측한다.

두 방법의 결과를 **공간적 근접성(centroid distance)**과 **잔기 수준 중첩(Jaccard
similarity)**으로 비교하여, 독립적 증거가 수렴하는 부위를 높은 신뢰도로 보고한다.

### 1.2 수용체 상태 다양성

단일 정적 구조의 한계를 극복하기 위해 3종의 수용체 상태를 사용한다:

| 수용체 ID | 출처 | 특성 |
|-----------|------|------|
| **3GT8_raw** | X-ray 결정구조 (PDB: 3GT8) | 기준 상태, 잔기 699–1007 |
| **3GT8_cl38_48** | MD 시뮬레이션 클러스터 | 동적 상태 1, 잔기 634–1014 |
| **3GT8_cl85_100** | MD 시뮬레이션 클러스터 | 동적 상태 2, 잔기 634–1014 |

MD 클러스터는 분자동역학 시뮬레이션 궤적에서 RMSD 기반 클러스터링으로 추출한 대표
구조로, 결정구조에서 포착되지 않는 C-lobe의 **동적 구조 변이(conformational
heterogeneity)**를 반영한다. 동일 포켓이 여러 수용체 상태에서 발견되면 구조적으로
안정한 결합 부위임을 시사한다.

### 1.3 리간드 선택

| 리간드 ID | 역할 |
|-----------|------|
| **173940** | 소분자 프로브 1 |
| **97806** | 소분자 프로브 2 |
| **VAX-C12_0** | 소분자 프로브 3 |

3종의 화학적으로 다양한 리간드를 사용하는 이유: 단일 리간드로는 특정 화학적 친화성에
편향된 결과를 얻을 수 있으나, **다수 리간드가 동일 포켓에 결합하면 해당 포켓의
"druggability"가 화학적으로 일반화된다** (cross-chemical consensus).

---

## 2. 시스템 구성

### 2.1 실행 환경

- **하드웨어**: Linux HPC (node05), 32 CPU cores
- **소프트웨어**: Python 3.x, AutoDock Vina 1.2, PyRosetta (ref2015 force field)
- **네트워크**: 차단 환경 (설치된 라이브러리만 사용)

### 2.2 전체 파이프라인 흐름도

```
┌─────────────────────────────────────────────────────────────┐
│  입력: 3 수용체 PDB/PDBQT + 3 리간드 SDF/PDBQT            │
│        2 PPI 타겟 (TH1, β-meander) PDB                    │
└───────────────┬─────────────────────────────┬───────────────┘
                │                             │
     ┌──────────▼──────────┐       ┌──────────▼──────────┐
     │ Phase 1: Vina Dock  │       │ Phase 2: PyRosetta  │
     │ 3×3×100 = 900 포즈  │       │   PPI Blind Dock    │
     │ ~30–45분            │       │ 50K×2 = 100K 모델   │
     └──────────┬──────────┘       │ ~36–48시간          │
                │                  └──────────┬──────────┘
                │                             │
                │                  ┌──────────▼──────────┐
                │                  │ Phase 3: PPI 후처리  │
                │                  │  체인 복원 + 잔기추출 │
                │                  └──────────┬──────────┘
                │                             │
     ┌──────────▼──────────┐                  │
     │ Phase 4: Vina 후처리 │                  │
     │  파싱→접촉→클러스터   │                  │
     │  →요약→비교→부트스트랩 │                  │
     └──────────┬──────────┘                  │
                │                             │
     ┌──────────▼─────────────────────────────▼───────────┐
     │ Phase 5: 3축 Verdict — Vina + PPI + Cross 통합     │
     └──────────┬─────────────────────────────────────────┘
                │
     ┌──────────▼──────────┐     ┌─────────────────────┐
     │ Phase 6: Report     │────▶│ Phase 7: Validation  │
     │  통합 리포트 생성    │     │  스키마/일관성 검증   │
     └─────────────────────┘     └─────────────────────┘
```

---

## 3. Phase 1: AutoDock Vina Blind Docking

### 3.1 방법론적 기초

AutoDock Vina는 **경험적 점수 함수(empirical scoring function)**에 기반한 분자
도킹 프로그램으로, Monte Carlo 시뮬레이션과 Broyden–Fletcher–Goldfarb–Shanno
(BFGS) 국소 최적화를 결합한다. Blind docking은 수용체 표면 전체를 탐색 영역으로
설정하여, 사전 지식 없이 결합 부위를 발견하는 비편향적(unbiased) 접근이다.

### 3.2 탐색 영역 설정

```
탐색 상자: 수용체 전체를 포함하는 직육면체
  최소 크기: 70 Å (각 축)
  여유 공간(padding): 5.0 Å
  모드: blind (수용체 경계 + padding으로 자동 계산)
```

Blind docking에서는 특정 결합 부위를 지정하지 않고 수용체 표면 전체를 탐색한다.
이를 위해 수용체의 모든 원자를 포함하는 최소 경계 상자에 5Å 여유를 추가하되,
최소 70Å를 보장한다. EGFR kinase C-lobe (~300 잔기)의 경우 약 70–80Å 크기의
탐색 상자가 설정된다.

### 3.3 점수 함수 (Vina Scoring Function)

AutoDock Vina 1.2의 점수 함수는 다음 항의 가중합이다:

| 에너지 항 | 물리적 의미 |
|-----------|------------|
| **가우시안 인력** | Lennard-Jones 유형 인력 (원자 간 최적 거리 ~3.5–4.0Å) |
| **반발력** | 원자 중첩에 대한 기하급수적 페널티 |
| **수소결합** | donor–acceptor 기하 구조 기반, 거리/각도 의존 |
| **소수성 상호작용** | 탄소-탄소 접촉에 대한 에너지 보상 |
| **회전 엔트로피 손실** | 0.595 kcal/mol × 회전 가능 결합 수 |

최종 결합 에너지(affinity)는 **kcal/mol** 단위이며, 음수값이 클수록 강한 결합을
의미한다.

### 3.4 샘플링 파라미터

| 파라미터 | 프로덕션 값 | 의미 |
|---------|------------|------|
| `exhaustiveness` | 384 | Monte Carlo 탐색 강도. 높을수록 탐색 공간을 치밀하게 샘플링 |
| `n_poses` | 100 | 수용체-리간드 쌍당 출력 포즈 수 |

**총 포즈 수**: 3 수용체 × 3 리간드 × 100 포즈 = **900 포즈**

`exhaustiveness=384`는 blind docking에서 높은 수준이다. Vina의 기본값은 8이며,
표적 도킹에서는 32–64가 일반적이나, C-lobe 전체 표면(~300잔기)을 blind로 탐색할
때는 높은 exhaustiveness가 희소 포켓의 누락을 방지한다.

### 3.5 출력 데이터

각 수용체-리간드 쌍에 대해:
```
output/{project}/{receptor_id}/{ligand}_blind.pdbqt
```
PDBQT 파일에 100개의 MODEL이 포함되며, 각 모델의 REMARK에 결합 에너지(kcal/mol),
RMSD lower bound, RMSD upper bound가 기록된다.

---

## 4. Phase 2: PyRosetta Global Blind PPI Docking

### 4.1 방법론적 기초

PyRosetta는 Rosetta 분자 모델링 플랫폼의 Python 인터페이스로, 물리 기반 에너지
함수(physics-based energy function)를 사용한다. PPI docking은 두 단백질 체인의
상대적 배향을 최적화하여 결합 인터페이스를 예측한다.

Global blind docking은 회전·병진의 전체 6차원 자유도를 탐색하여, 알려진 결합
부위 정보 없이 모든 가능한 결합 배향을 샘플링한다.

### 4.2 에너지 함수: ref2015

Rosetta ref2015 (Reference Energy Function 2015)는 다음 에너지 항의 가중합이다:

| 에너지 항 | 기호 | 물리적 의미 |
|-----------|------|------------|
| **Lennard-Jones 인력** | `fa_atr` | 원자 간 van der Waals 인력 (r⁻⁶) |
| **Lennard-Jones 반발** | `fa_rep` | 원자 중첩 페널티 (r⁻¹²) |
| **용매화** | `fa_sol` | 극성/비극성 원자의 매장 에너지 (Lazaridis-Karplus 모델) |
| **정전기** | `fa_elec` | 거리 의존 Coulomb 상호작용 |
| **수소결합** | `hbond_*` | donor–acceptor 기하 의존 (bb-bb, bb-sc, sc-sc) |
| **주쇄 비틀림** | `rama_prepro` | Ramachandran 선호 각도 |
| **로타머 에너지** | `fa_dun` | Dunbrack 로타머 라이브러리 기반 |
| **참조 에너지** | `ref` | 아미노산별 보정 상수 |

단위는 **Rosetta Energy Unit (REU)**이며, kcal/mol과 대략적으로 대응하나 정확한
변환 계수는 없다.

### 4.3 7단계 워크플로우

#### Step 1: Relax (구조 이완)

```
프로토콜: FastRelax (ref2015)
목적: 입력 PDB의 입체 충돌(steric clash) 해소, 에너지 최소화
캐싱: relaxed_cache/ 에 저장 → 반복 실행 시 재사용
```

FastRelax는 주쇄(backbone)와 곁사슬(sidechain)을 동시에 이완하되, Rosetta의
`fa_rep` 가중치를 점진적으로 증가시키는 "ramp-up" 전략을 사용한다. 이는 구조가
국소 최소(local minimum)에 빠지는 것을 방지하며, 도킹 이전에 에너지적으로 타당한
출발점을 제공한다.

#### Step 2: Global Docking (전역 도킹)

```
1단계: RigidBodyPerturbMover (360°, 100Å)
  → 완전 무작위 배향: 회전 360°, 병진 100Å
  → 표면 전체를 편향 없이 탐색

2단계: DockingSlideIntoContact
  → 두 체인을 접촉 거리까지 이동 (van der Waals 접촉)
  → ★ 필수 단계: 이 없이 DockMCMProtocol 실행 시 모든 dG가 0.0

3단계: DockMCMProtocol
  → Monte Carlo + Minimization 반복
  → 리지드바디 섭동 + 곁사슬 리패킹 + 에너지 최소화
  → Metropolis 수용 기준으로 결합 배향 최적화
```

**Early Rejection (선택적)**: `enable_early_rejection=True` 시, SlideIntoContact
후 DockMCMProtocol 이전에 Chain B가 금지 영역(excluded_residues_A)과 접촉하는지
검사한다. 접촉 시 해당 모델을 즉시 폐기하여 ~80–90%의 연산을 절약한다.

**샘플링 규모**: 50,000 모델 × 2 타겟 = 100,000 독립 도킹 시행

#### Step 3: 2-Pass 필터링 (v2.0)

도킹 생존자를 **2단계 필터링**으로 선별한다:

**Pass 1 (전체 대상, 저비용 메트릭)**:

| 메트릭 | 계산 비용 | 물리적 의미 |
|--------|----------|------------|
| `dG_separated` | ~1초 | 인터페이스 결합 에너지 = E(복합체) − E(A) − E(B) |
| `dSASA` | ~1초 | 결합 시 매장되는 용매 접근 가능 표면적 (Å²) |
| `sc_value` | ~1초 | 형상 상보성 (Shape Complementarity, 0–1) |
| `total_score` | ~0초 | Rosetta 전체 에너지 (필터링용) |

**Stage 1 필터**: 금지영역 접촉 제거 → dG > 0 제거 → total_score 백분위 컷오프
(상위 10%)

**Mini Refinement**: Stage 1 생존자에 인터페이스 곁사슬 리패킹 적용

```
방법: IncludeCurrent + ExtraRotamersGeneric (ex1 + ex2 chi 확장)
라운드: 3회 반복
목적: 리지드바디 도킹 후 불량한 곁사슬 패킹 보완
      → Stage 2 통과율 2–5% 향상
```

**Pass 2 (Stage 1 생존자만, 고비용 메트릭)**:

| 메트릭 | 계산 비용 | 물리적 의미 |
|--------|----------|------------|
| `packstat` | ~5초 | 원자 패킹 밀도 (0–1, Rosetta 고유) |
| `delta_unsatHbonds` | ~5초 | 결합 시 미충족 수소결합 수 |
| `nres_int` | ~3초 | 인터페이스 잔기 수 |
| `hbonds_int` | ~3초 | 인터페이스 수소결합 수 |

**Stage 2 필터**: dSASA ≥ 500Å², dG_density ≤ −1.0, sc ≥ 0.50,
packstat ≥ 0.55, delta_unsatHbonds ≤ 8, nres_int ≥ 10, hbonds_int ≥ 0

**Graduated Fallback**: 엄격한 필터로 생존자가 부족하면 단계적으로 완화:

```
Level 0: 모든 필터 적용 (충분한 생존자)
Level 1: 고비용 메트릭 해제 (dSASA + sc + dG_density만 유지)
Level 2: sc + dG_density 해제, dSASA 임계값 50% 완화
Level 3: Stage 2 전체 해제 → Stage 1 생존자에서 dG 상위 N개
```

#### Step 4: L_RMSD 기반 Greedy 클러스터링

필터 통과 구조들을 **결합 부위(binding site)**별로 그룹화한다.

```
알고리즘: Greedy 클러스터링 (dG 순 정렬)
거리 메트릭: Ligand RMSD (Chain B의 Cα 원자 RMSD)
```

**CoM Pre-gate** (속도 최적화):
```
CoM 임계값 = max(cluster_threshold × 8, 30Å)
Chain B 질량 중심(Center of Mass) 간 거리가 이 임계값을 초과하면
L_RMSD 계산을 건너뛴다 → O(N²) RMSD 계산 회피
```

**클러스터링 루프**:
```
for candidate in sorted_by_dG:
    for leader in existing_leaders:
        if CoM_dist(candidate, leader) < CoM_threshold:
            if L_RMSD(candidate, leader) < cluster_threshold:
                → 해당 클러스터에 배정
                break
    if 미배정:
        if len(leaders) < cluster_top_n:
            → 새 클러스터 생성
        else:
            → Shadow Archive에 기록 (dropped_candidates.csv)
```

**적응적 확장 (Adaptive Escalation)**:
```
potential_extra_clusters (cluster_top_n 초과 후보 수)가
max(cluster_top_n × 0.2, 3) 이상이면:
  → cluster_top_n을 1.5배 확장 (1회 한정)
  → 추가 클러스터 생성 허용
```

**Shadow Archive**: 클러스터 용량 초과로 탈락한 모든 후보의 메타데이터를
`dropped_candidates.csv`에 기록한다 (PDB 파일 없이 좌표/에너지/거리만).
이는 탈락 후보의 사후 분석과 데이터 손실 추적을 가능하게 한다.

#### Step 5: 다양성 기반 최종 선택

```
알고리즘: Round-Robin 다양성 선택
  1. 클러스터를 에너지순 정렬
  2. 각 라운드에서 클러스터당 1개씩 선택
  3. 모든 클러스터를 순환하며 목표 수(save_top_n)까지 채움
  4. L_RMSD < 1.0Å 중복 제거
```

**Archive Ranking**: save_top_n보다 넓은 범위(archive_top_n, 기본 3배)의 모델
메타데이터를 `archive_ranking.csv`에 보존한다 (PDB 파일 없이).

#### Step 6: 국소 정밀화 (Refinement)

선택된 대표 구조에 DockMCMProtocol을 재적용한다:
```
병진 섭동: 0.1Å (global의 1/1000)
회전 섭동: 1.0° (global의 1/360)
→ 결합 배향의 미세 조정
```

#### Step 7: 시각화 및 검증 리포트

**PyMOL 스크립트 자동 생성**:
- `1_OVERVIEW_Clusters.pml`: 전체 결합 부위 분포 (클러스터별 색상)
- `2_DETAIL_C##.pml`: 개별 부위 포즈 수렴도
- `view_results.pml`: 최종 모델 B-factor 컬러링

**10-Point 품질 검증 (Validation Report)**:

| 검사 | 항목 | 판정 기준 |
|------|------|----------|
| C1 | 파이프라인 실행 성공률 | ≥90% PASS |
| C2 | 결합 에너지 | best dG < −10 REU: Excellent |
| C3 | 에너지 퍼널 | dG gap > 10 REU + P_near > 0.1: Excellent |
| C4 | 인터페이스 크기 | mean dSASA > 500Å²: PASS |
| C5 | 형상 상보성 | mean sc > 0.50: PASS |
| C6 | 부위 수렴도 | 지배 부위 > 50%: Excellent |
| C7 | 샘플링 규모 | ≥10,000 모델: PASS |
| C8 | dG 밀도 | mean dG/dSASA×100 < −1.0: PASS |
| C9 | 인터페이스 잔기 수 | mean nres > 15: PASS |
| C10 | 수소결합 | mean hbonds ≥ 1: PASS |

### 4.4 스코어링 메트릭 해석 기준

| 메트릭 | 우수 | 양호 | 경계 | 불량 |
|--------|------|------|------|------|
| dG_separated (REU) | < −10 | −10 ~ −3 | −3 ~ 0 | > 0 |
| dSASA (Å²) | > 1000 | 500–1000 | 300–500 | < 300 |
| sc_value (0–1) | > 0.70 | 0.50–0.70 | 0.35–0.50 | < 0.35 |
| packstat (0–1) | > 0.55 | 0.40–0.55 | — | < 0.40 |
| dG_density | < −1.5 | −1.0 ~ −1.5 | −0.5 ~ −1.0 | > −0.5 |
| delta_unsatHbonds | < 5 | 5–10 | — | > 10 |
| nres_int | > 20 | 15–20 | 8–15 | < 8 |
| hbonds_int | ≥ 3 | 1–3 | — | 0 |

### 4.5 제약 시스템 (Constraints)

**금지 영역 (Hard Filter)**:
```
excluded_residues_A: Chain A의 막면/다이머 인터페이스 잔기
  → 이 영역과 접촉하는 포즈는 물리적으로 불가능한 결합 배향
  → 즉시 제거 (Early Rejection 또는 필터링)
```

**핵심 잔기 (Soft Bonus)**:
```
key_residues_B: 실험 데이터 기반 Chain B 핵심 잔기
  → key_contact_ratio = 핵심 잔기 접촉 비율
  → adjusted_dG = dG − bonus_weight × key_contact_ratio
  → 에너지 순위에 가산점 부여 (정보적, 강제 아님)
```

---

## 5. Phase 3: PPI 후처리

### 5.1 체인 복원 (Chain Restoration)

PyRosetta 도킹은 두 체인을 단일 번호 체계로 병합하여 작동한다 (A: 699–1007,
B: 1708–2015). 결과 구조에서 원래 체인 ID와 잔기 번호를 복원해야 한다.

```
입력: 병합된 PDB (모든 잔기가 단일 체인)
매핑: mapping.csv (merged_resnum → original_chain + original_resnum)
출력: 복원된 PDB (Chain A + Chain B 분리)
```

복원 대상: `final_result/*.pdb`, `cluster_results/*.pdb`, `final_ranking.csv`,
`cluster_summary.csv`

### 5.2 PPI 잔기 추출

복원된 `final_ranking.csv`에서 인터페이스 잔기 정보를 추출한다:

```
per residue:
  receptor_id: 수용체 ID
  residue_id: 잔기 이름+번호 (예: ALA702)
  occupancy: final ranking 모델들 중 이 잔기가 인터페이스에 나타나는 빈도
  mean_interface_delta_e: 해당 잔기의 평균 인터페이스 에너지 기여
```

`occupancy`는 **해당 잔기가 PPI 인터페이스의 핵심 구성요소인지**를 나타내는 핵심
지표이다. Occupancy > 0.5이면 최종 모델의 절반 이상에서 인터페이스에 참여하므로
안정적인 접촉 잔기로 해석한다.

---

## 6. Phase 4: Vina 후처리 파이프라인

### 6.1 포즈 테이블 구축 (Parse)

Vina 출력 PDBQT 파일을 파싱하여 구조화된 테이블로 변환한다:

```
per pose:
  receptor_id, ligand_id, pose_rank
  affinity (kcal/mol)
  centroid_x, centroid_y, centroid_z (리간드 원자 평균 좌표)
  raw_pose_file (PDBQT 경로)
```

**중심점(centroid)**: 리간드의 모든 비수소 원자 좌표의 산술 평균. 결합 포켓 내부의
위치를 나타내며, 포즈 간 공간적 유사성 비교의 기초가 된다.

### 6.2 접촉 잔기 추출 (Contacts)

각 도킹 포즈에 대해 수용체의 접촉 잔기를 식별한다:

```
알고리즘: 원자쌍 거리 스캔
  for each receptor heavy atom:
      for each ligand heavy atom:
          if distance ≤ 4.0Å:
              record residue (min_distance로 중복 제거)
```

**4.0Å 컷오프의 근거**: van der Waals 반지름의 합(~3.4Å)에 약간의 여유를 더한
값으로, 비결합 상호작용(van der Waals 접촉, 수소결합)이 일어나는 범위이다.
너무 넓으면(>5Å) 용매 매개 상호작용이 포함되고, 너무 좁으면(<3Å) 실제 접촉을
놓친다.

### 6.3 반복적 중심점 기반 클러스터링 (Clustering)

**이 단계는 데이터 손실 방지를 위해 전면 재설계되었다.**

#### 6.3.1 알고리즘 개요

k-means 스타일의 반복적 클러스터링으로, 기존 단일 패스 탐욕적(greedy) 알고리즘의
**중심점 표류(centroid drift)**와 **입력 순서 의존성** 문제를 해결한다.

```
Phase A: 시드 생성 — 친화력 순서로 중심점 스캔,
         기존 시드와 cutoff(7Å) 이내이면 건너뛰고 아니면 새 시드 생성
         ★ 시드 위치는 생성 시 고정 (이동하지 않음)

Phase B: 초기 배정 — 각 포즈를 가장 가까운 시드에 배정

Phase C: 반복 정밀화 (최대 15회)
         1. 시드 위치를 배정된 포즈들의 중심점 평균으로 재계산
         2. 모든 포즈를 새 시드 위치에 재배정
         3. 배정이 변하지 않으면 수렴 → 종료

Phase D: 소포켓 흡수 — min_pocket_size(3) 미만의 포켓을
         가장 가까운 대형 포켓에 병합

Phase E: 포켓 ID 부여 — P001, P002, ... (출현 순서)
```

#### 6.3.2 기존 알고리즘의 문제점과 해결

| 문제 | 원인 | 해결 |
|------|------|------|
| **중심점 표류** | 시드에 포즈 추가 시 즉시 중심점 이동 | 시드 생성 단계에서 위치 고정 |
| **입력 순서 의존** | 순서에 따라 시드 생성/배정 결과 변동 | 반복 수렴 → 동일 결과 보장 |
| **포켓 분절** | 좁은 cutoff(4Å)로 인접 포즈 분리 | cutoff 7Å + 잔기 기반 병합 |
| **노이즈 포켓** | 단일 포즈가 독립 포켓 형성 | min_pocket_size=3 흡수 |

#### 6.3.3 잔기 기반 포켓 병합

클러스터링 후, 접촉 잔기 중첩도에 기반한 포켓 병합을 수행한다:

```
알고리즘: Union-Find (경로 압축)

두 포켓 A, B에 대해:
  residues_A = 포켓 A 전체 포즈의 접촉 잔기 합집합
  residues_B = 포켓 B 전체 포즈의 접촉 잔기 합집합

  Jaccard = |A ∩ B| / |A ∪ B|
  Overlap = |A ∩ B| / min(|A|, |B|)

  병합 조건 (OR):
    Jaccard ≥ 0.25  (동일 잔기 집합)
    Overlap ≥ 0.4   (부분집합 관계)
    centroid 거리 ≤ 8.0Å  (잔기 데이터 부족 시 공간적 근접성 대체)
```

**왜 Jaccard와 Overlap을 OR로 결합하는가?**

Jaccard는 대칭적이므로 크기가 비슷한 포켓 간 비교에 적합하다. 그러나 큰 포켓의
부분이 작은 포켓으로 분절된 경우, Jaccard는 분모(합집합)가 커서 낮게 나오지만
Overlap coefficient는 작은 포켓 기준으로 높게 나온다. OR 조건으로 두 경우 모두
포착한다.

**centroid fallback의 역할**: 접촉 잔기 추출이 불완전한 경우(포즈가 표면에
불완전하게 결합, 또는 수소 원자 누락 등), 잔기 기반 기준만으로는 분절을 해소할 수
없다. 8.0Å 이내 중심점이면 동일 포켓으로 판정하는 안전망을 제공한다.

**Union-Find 전이적 닫힘**: A-B 병합 + B-C 병합 → A-C도 자동 병합된다.

### 6.4 포켓 요약 (Summarize)

각 포켓에 대해 집계 통계를 계산한다:

| 통계 | 계산 | 해석 |
|------|------|------|
| **centroid** | 포즈 중심점 평균 | 포켓의 3차원 위치 |
| **centroid_spread** | √(Σd²/N) | 포즈 수렴도 (<3Å: 수렴, >5Å: 분산) |
| **n_pose** | 포즈 수 | 포켓의 샘플링 빈도 |
| **n_ligand** | 결합 리간드 종류 수 | 화학적 일반성 (3=높음) |
| **best_affinity** | min(affinity) | 최적 결합 에너지 |
| **mean_affinity** | mean(affinity) | 평균 결합 에너지 |
| **top_residues** | 빈도 상위 5 잔기 | 핵심 접촉 잔기 |
| **union_contact_residues** | 모든 포즈 접촉 잔기 합집합 | 포켓 전체 풋프린트 |

**Drug-Pocket Map** (다중 결합 모드 감지):

```
per (receptor, ligand):
  dominant_pocket: 최다 포즈를 포함한 포켓
  dominant_pocket_fraction: 해당 포켓의 포즈 비율
  is_multimodal_binding: fraction < 1.0이면 다중 모드
```

해석: `fraction > 0.8` → 단일 결합 모드 (monodal),
`fraction 0.5–0.8` → 두 부위 경쟁 (bimodal),
`fraction < 0.5` → 다수 부위 분산 (promiscuous)

### 6.5 교차 수용체 포켓 비교 (Compare)

3종 수용체 상태에서 발견된 포켓들의 대응 관계를 정량화한다:

| 비교 메트릭 | 계산 | 해석 |
|------------|------|------|
| **centroid_distance** | 유클리드 거리 | < 8Å: 동일 포켓 |
| **residue_jaccard** | \|A∩B\| / \|A∪B\| | > 0.3: 유의미한 잔기 중첩 |
| **residue_overlap_coeff** | \|A∩B\| / min(\|A\|,\|B\|) | > 0.5: 부분집합 관계 |
| **shared_ligands** | 공유 리간드 수 | > 0: 화학적 교차 검증 |

**Same-Patch Candidate** 판정:
```
centroid_dist < 8Å AND (Jaccard ≥ 0.3 OR Overlap ≥ 0.5)
→ "동일 결합 부위가 다른 수용체 상태에서도 보존됨"
```

이는 구조적 안정성의 강력한 증거이다: 결정구조와 MD 클러스터에서 동일 포켓이
발견되면, 해당 포켓은 열역학적으로 접근 가능한(thermodynamically accessible)
결합 부위일 가능성이 높다.

### 6.6 부트스트랩 안정성 분석 (Bootstrap)

포켓 배정의 **통계적 강건성(robustness)**을 정량화한다.

#### 방법

```
for i in 1..200:
    1. 원본 포즈의 80%를 복원 추출(with replacement)로 리샘플링
    2. 리샘플링 데이터로 전체 클러스터링 재실행
    3. 레플리카 포켓을 참조 포켓에 centroid 근접성으로 매칭
    4. 매칭된 포켓의 통계 기록
```

#### 출력 메트릭

| 메트릭 | 의미 | 해석 |
|--------|------|------|
| `pocket_exists_frac` | 200회 중 해당 포켓 재현 비율 | > 0.8: 안정적 |
| `centroid_std_A` | 중심점 위치 표준편차 (Å) | < 2Å: 위치 안정 |
| `affinity_mean/std` | 결합 에너지 평균/표준편차 | std 작을수록 안정 |
| `affinity_iqr` | 사분위 범위 | 이상치에 강건한 산포 지표 |
| `n_pose_mean/std` | 포즈 수 평균/표준편차 | 포켓 크기 안정성 |

**복원 추출을 사용하는 이유**: 일부 포즈가 중복 출현하고 일부는 누락되는 효과를
통해, "약간 다른 도킹 시행을 했다면 결과가 어떻게 달라졌을까?"라는 인식론적
불확실성(epistemic uncertainty)을 추정한다.

---

## 7. Phase 5: 3축 Verdict 통합 평가

### 7.1 설계 철학

Verdict 시스템은 **증거 강도 분류(evidence classification)**이며, 타당성 판정
(validity judgment)이 아니다. STRONG으로 분류된 포켓이라도 반드시 시각적 검토가
필요하고, WEAK 포켓도 생물학적으로 유의미할 수 있다.

### 7.2 3축 점수 체계

#### 축 1: Vina 품질 (직접적 증거, 최대 50점)

EGFR C-lobe 표면 포켓에 대해 보정된 점수:

```
결합 에너지 (0–20점):
  < −8.0 kcal/mol → 20점 (표면 포켓에서는 매우 강한 결합)
  < −6.5 kcal/mol → 15점 (전형적인 표면 포켓)
  < −5.0 kcal/mol →  8점 (경계값)
  ≥ −5.0 kcal/mol →  0점

수렴도 (0–15점):
  n_pose ≥ 8  → 15점 (blind docking의 ~40%가 한 포켓에 수렴)
  n_pose ≥ 3  → 10점
  n_pose ≥ 1  →  5점

다중 리간드 합의 (0–15점):
  3종 리간드 결합 → 15점 (화학적으로 일반화된 포켓)
  2종 리간드 결합 →  8점
  1종 리간드 결합 →  2점
```

**보정 맥락**: ATP 활성 부위(−9 ~ −12 kcal/mol)가 아닌 C-lobe **표면 포켓**(−5
~ −8 kcal/mol)이 주 탐색 대상이므로, 임계값이 일반적인 활성 부위 도킹보다 높게
설정되어 있다.

#### 축 2: PPI 공간적 근접성 (간접적 증거, 최대 20점)

"이 약물 결합 포켓이 MYO1D 결합 인터페이스 근처인가?"

```
PPI 중심점까지 거리:
  < 8Å  → 20점 (직접 인터페이스, 경쟁적 결합 가능)
  < 15Å → 15점 (알로스테릭 범위)
  < 25Å → 8점  (동일 도메인)
  ≥ 25Å → 0점

잔기 중첩 보너스:
  공유 잔기 > 0 → +5점
```

**중심점 의미론적 차이에 대한 주의**:
- Vina 중심점 = 리간드 원자 평균 (결합 포켓 **내부**)
- PPI 중심점 = 수용체 Cα 원자 평균 (단백질 **표면**)
- 체계적 오프셋: Vina 중심점이 표면에서 ~3–5Å 더 깊음
- 임계값(8/15/25Å)은 이 오프셋을 반영하여 보정됨

#### 축 3: 교차 수용체 일관성 (구조적 증거, 최대 30점)

```
Same-Patch Candidate가 2+ 수용체에서 발견 → 30점
Same-Patch Candidate가 1 수용체에서 발견 → 15점
없음 → 0점

부트스트랩 안정성 보너스:
  centroid_std < 2.0Å → +5점
```

#### 적응적 가중치

PPI 데이터가 없는 수용체에 대해 페널티를 주지 않기 위한 정규화:

```
PPI 데이터 있음: 총점 = 축1(50) + 축2(20) + 축3(30) = 100
PPI 데이터 없음: 총점 = 축1(60) + 축2(0)  + 축3(40) = 100
```

### 7.3 증거 분류

| 총점 | 분류 | 해석 |
|------|------|------|
| ≥ 55 | **STRONG** | 복수 축에서 수렴하는 강한 증거 |
| 30–54 | **MODERATE** | 부분적 증거, 추가 검증 권장 |
| < 30 | **WEAK** | 탐색적 결과, 시각적 검토 필요 |

### 7.4 합의 부위 그룹화 (Consensus Site)

동일 포켓이 여러 수용체에서 발견되면 **합의 부위 ID** (CS_001, CS_002, ...)를
부여한다:

```
per consensus site:
  n_receptors: 포함 수용체 수
  receptor_list: 수용체 ID 목록
  centroid: 수용체 간 중심점 평균
  best_affinity: 전체 최적 에너지
  consensus_residues: 잔기 합집합
```

---

## 8. Phase 6–7: 리포트 생성 및 검증

### 8.1 통합 리포트 구조

| 섹션 | 내용 |
|------|------|
| **1. 수용체별 포켓 요약** | 포켓 수, 리간드-포켓 매핑, 다중 결합 모드 |
| **2. 교차 수용체 비교** | Same-Patch Candidate 테이블, 공유 잔기, 부트스트랩 CI |
| **3. PPI 보조 증거** | PyRosetta 인터페이스 잔기, occupancy 테이블 |
| **4. Verdict 판정** | STRONG/MODERATE/WEAK 분포, 포켓별 점수 상세 |
| **5. 핵심 관찰** | 통계 요약, 주요 발견 |

### 8.2 출력 검증 (Validation)

8단계 자동 검증 시스템:

```
1. 출력 파일 존재 확인 (필수/선택 구분)
2. CSV 스키마 회귀 검사 (필드 누락 = FAIL, 추가 = WARN)
3. 잔기 번호 일관성 (수용체 간 번호 중첩, 아미노산 정합성)
4. 알려진 변이 확인 (결정화 아티팩트 등)
5. 오프셋 감지 (수용체 간 번호 체계 이동)
6. 교차 테이블 일관성 (receptor_id, ligand_id 정합)
7. 파일 참조 추적 (raw_pose_file 경로 유효성)
8. 배포 준비 확인 (필수 문서/모듈 존재)
```

종료 코드: 0 (통과), 1 (경고), 2 (실패)

---

## 9. 통계적 신뢰성 분석

### 9.1 부트스트랩 분석의 해석

200회 부트스트랩 레플리카에서:

- **pocket_exists_frac > 0.8**: 데이터의 80% 이상에서 재현되는 강건한 포켓
- **centroid_std_A < 2.0Å**: 포켓 위치가 데이터 변동에 안정적
- **affinity_iqr < 1.0 kcal/mol**: 에너지 추정이 안정적

200회 레플리카에서 pocket_exists_frac의 95% 신뢰구간 폭은 ±~7%이다
(이항 분포 근사).

### 9.2 파라미터 민감도 분석 (Sweep)

pocket_cutoff를 2.0–12.0Å 범위에서 0.5Å 단위로 변화시키며 포켓 수와 크기 분포의
변화를 관찰한다:

```
낮은 cutoff (2–4Å): 과도한 분절 → 다수의 작은 포켓
중간 cutoff (6–8Å): 균형적 → 의미 있는 포켓 그룹화
높은 cutoff (10–12Å): 과도한 병합 → 구분되어야 할 포켓이 합쳐짐
```

**Elbow point** (포켓 수 대 cutoff의 변곡점)에서 최적 cutoff를 선택한다.
EGFR C-lobe에서는 전형적으로 6–8Å가 적절하다.

### 9.3 교차 수용체 일관성의 통계적 의미

3종 독립 수용체 상태에서 동일 포켓이 발견될 확률:
- **무작위 기대**: 표면적 ~8000Å²에서 15Å 반경 영역이 겹칠 확률은 ~5–10%
- **Same-Patch Candidate**: 유의미한 일치 (p < 0.05 수준)

---

## 10. 결과 해석 가이드

### 10.1 권장 분석 순서

```
1. Validation Report 확인
   → C1–C10 PASS/FAIL 항목 검토
   → 핵심 실패(C2, C3, C4)가 있으면 결과 해석에 주의

2. valid_sites.csv 검토
   → STRONG 포켓부터 확인
   → confidence_score와 각 축 점수 분포 파악

3. 교차 수용체 비교
   → same_patch_candidate 포켓 우선 검토
   → 여러 수용체에서 보존된 포켓 = 높은 구조적 신뢰도

4. PyMOL 시각적 검증 (필수)
   → 1_OVERVIEW_Clusters.pml로 전체 분포 확인
   → 에너지 좋은 클러스터부터 2_DETAIL_C##.pml 검토
   → 생물학적 타당성 판단 (막면 접촉 여부, 포켓 깊이 등)

5. 정량적 메트릭 비교
   → final_ranking.csv에서 dG + dSASA + sc + packstat 종합 판단
   → cluster_summary.csv에서 포켓별 population 비교
```

### 10.2 흔한 질문

**Q: Vina에서 STRONG인데 PPI 근접성이 낮으면?**
A: 해당 포켓은 druggable하지만 MYO1D 인터페이스와 무관할 수 있다. 다른 PPI
파트너의 결합 부위이거나, allosteric 조절 부위일 가능성을 고려한다.

**Q: PPI에서는 높은 occupancy인데 Vina 포켓이 없으면?**
A: PPI 인터페이스가 소분자로 druggable하지 않은 평탄한 표면일 수 있다. 이는
PPI 저해제 개발의 일반적인 한계이다.

**Q: L_RMSD가 높은데 괜찮은가?**
A: Global blind docking에서 **높은 L_RMSD는 정상**이다. 표면 전체를 탐색하므로
서로 먼 위치에 부위가 발견되는 것이 자연스럽다. 낮은 L_RMSD는 대부분의 포즈가
한 곳에 수렴한 것으로 높은 신뢰도를 의미한다.

---

## 11. 파라미터 요약 및 근거

### 11.1 Vina 도킹

| 파라미터 | 값 | 근거 |
|---------|-----|------|
| exhaustiveness | 384 | ~300잔기 blind docking에 필요한 탐색 밀도 |
| n_poses | 100 | 900 total → min_pocket_size=3 기준 최대 300 포켓 커버 |
| padding | 5.0Å | 수용체 경계 외 탐색 여유 |
| min_box | 70.0Å | C-lobe 전체 포함 최소 크기 |

### 11.2 클러스터링

| 파라미터 | 값 | 근거 |
|---------|-----|------|
| pocket_cutoff | 7.0Å | 초기 분리 보수적, merge_centroid_fallback이 보완 |
| merge_jaccard | 0.25 | sparse 접촉 데이터에서도 동일 포켓 감지 |
| merge_overlap | 0.4 | 큰 포켓의 부분 분절 흡수 |
| merge_centroid_fallback | 8.0Å | 포켓 직경 ~15Å의 절반, 잔기 데이터 부족 시 안전망 |
| cluster_max_iterations | 15 | 900포즈 규모 수렴 보장 |
| min_pocket_size | 3 | 3개 독립 도킹에서 재현 = 노이즈 아님 |

**설계 원칙**:
```
pocket_cutoff(7Å) < merge_centroid_fallback(8Å) < comparison_cutoff(12Å)
  초기분리               같은포켓병합              교차수용체비교
```

### 11.3 부트스트랩

| 파라미터 | 값 | 근거 |
|---------|-----|------|
| n_replicates | 200 | 95% CI 폭 ±7% (이항 분포 기준) |
| sample_fraction | 0.8 | 분산 추정과 샘플 대표성의 균형 |
| seed | 42 | 재현성 보장 (seed + i로 레플리카별 분리) |

### 11.4 교차 비교

| 파라미터 | 값 | 근거 |
|---------|-----|------|
| comparison_centroid_cutoff | 12.0Å | 15Å에서 축소: false positive 감소 |
| same_patch: centroid < 8Å | 8.0Å | 동일 포켓 직접 겹침 범위 |
| same_patch: Jaccard ≥ 0.3 | 0.3 | 잔기 30% 공유 = 유의미 |

---

## 12. 제한사항 및 향후 과제

### 12.1 방법론적 제한

1. **점수 함수의 근사**: Vina와 Rosetta 모두 경험적/반경험적 에너지 함수를
   사용하며, 엔트로피 효과와 물 매개 상호작용을 완전히 포착하지 못한다.

2. **리지드 수용체 가정**: Vina blind docking은 수용체를 강체로 처리한다.
   유도 적합(induced fit)은 포착되지 않으며, MD 클러스터 사용이 부분적 보완이다.

3. **PPI 도킹의 샘플링 한계**: 50K 모델은 6차원 구성 공간의 미세한 부분만
   탐색한다. 에너지 최소 구배가 얕은 결합 부위는 놓칠 수 있다.

4. **중심점 기반 비교의 한계**: Vina(리간드 원자)와 PPI(수용체 Cα)의 중심점은
   물리적 의미가 다르며, 3–5Å의 체계적 오프셋이 존재한다.

5. **소분자 프로브의 대표성**: 3종 리간드가 모든 druggable pocket을 탐지한다고
   보장할 수 없다. 화학 공간의 작은 부분만 탐색한다.

### 12.2 향후 과제

- **AlphaFold-Multimer** 통합: 3축에 4번째 독립 증거 추가
- **앙상블 도킹**: MD 궤적의 다중 프레임에 대한 Vina 도킹
- **자유 에너지 계산**: MM-PBSA/GBSA를 통한 결합 에너지 정밀화
- **실험적 검증**: 교차 결합(cross-linking) 질량분석 또는 HDX-MS

---

## 13. 출력 파일 명세

### 13.1 Vina 파이프라인 출력

| 파일 | 설명 | 주요 필드 |
|------|------|----------|
| `vina_pose_table.csv` | 전체 포즈 테이블 | receptor_id, ligand_id, affinity, centroid_xyz, pocket_id, contact_residues |
| `vina_pocket_table.csv` | 포켓별 요약 | pocket_id, n_pose, n_ligand, best/mean_affinity, centroid_spread, union_contact_residues |
| `vina_drug_pocket_map.csv` | 리간드-포켓 매핑 | ligand_id, dominant_pocket, fraction, is_multimodal |
| `vina_pocket_comparison.csv` | 교차 수용체 비교 | receptor_a/b, pocket_a/b, centroid_dist, jaccard, same_patch_candidate |
| `vina_pocket_bootstrap.csv` | 부트스트랩 통계 | pocket_id, pocket_exists_frac, centroid_std_A, affinity_mean/std |

### 13.2 PyRosetta PPI 출력

| 파일 | 설명 | 주요 필드 |
|------|------|----------|
| `final_ranking.csv` | 최종 모델 순위 | Rank, Cluster_ID, dG, dSASA, sc, packstat, Binding_Residues |
| `cluster_summary.csv` | 클러스터별 요약 | Cluster_ID, Population, mean_dG, mean_dSASA |
| `dropped_candidates.csv` | 탈락 후보 기록 | center_xyz, dG, dSASA, drop_reason |
| `archive_ranking.csv` | 확장 메타데이터 | top N×3 모델의 전체 메트릭 |
| `validation_report.txt` | 10-point 품질 검증 | C1–C10 PASS/FAIL |

### 13.3 통합 출력

| 파일 | 설명 | 주요 필드 |
|------|------|----------|
| `valid_sites.csv` | Verdict 판정 | pocket_id, verdict, confidence_score, 3축 점수 |
| `cross_method_agreement.csv` | Vina↔PPI 비교 | spatial_dist_A, jaccard, shared_residues |
| `ppi_pyrosetta_residues.csv` | PPI 잔기 추출 | residue_id, occupancy, mean_delta_e |
| `combined_residue_evidence.csv` | 다중 증거 잔기 | evidence_sources (vina/ppi/afm) |
| `project_report.txt` | 통합 리포트 | 5개 섹션 텍스트 + 테이블 |

---

*본 문서는 파이프라인의 모든 계산적 방법론, 알고리즘, 파라미터, 출력을 기술한다.
논문 작성 시 Methods 섹션의 기초 자료로 사용할 수 있다.*
