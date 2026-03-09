# PyRosetta PPI Docking Pipeline - Manual

## 1. Overview

PyRosetta 기반 **단백질-단백질 상호작용(PPI) Global Blind Docking** 파이프라인.
입력 PDB 구조에 대해 구조 이완, 전역 도킹, 다중 단계 필터링, Full Scoring, L_RMSD 클러스터링을 수행하여 최적의 결합 포즈를 예측한다.

- **버전**: v2.0 (필터링 파이프라인 개편) / v1.0 하위 호환
- **실행 환경**: Linux HPC Cluster (PBS/qsub), 32 CPU cores
- **Python**: 3.x + PyRosetta (네트워크 차단 환경, 구버전 호환성 확보)
- **병렬 처리**: `multiprocessing.Pool` 기반 CPU 병렬

### v2.0 주요 변경사항

| 항목 | v1.0 | v2.0 |
|------|------|------|
| 필터링 단계 | 3단계 (dSASA → dG 백분위 → sc) | 2-Stage 다중 메트릭 (Stage 1 Coarse + Stage 2 Quality) |
| 스코어링 Pass | 1-Pass (Fast만) | 2-Pass (Fast + Intermediate) |
| 필터 메트릭 수 | 3개 | 8개 (total_score, dG, dSASA, sc, dG_density, packstat, unsatHb, nres_int, hbonds_int) |
| Fallback | 3단계 | 4단계 (Level 0~3) |
| Validation Report | 7개 품질 체크 | 10개 품질 체크 (C1~C10) |
| 하위 호환 | - | `[FilterStage1]` 섹션 없으면 자동으로 v1.0 경로 |

---

## 2. File Structure

```
Pyrosetta_PPI/
├── pipeline_manager.py        # 파이프라인 오케스트레이터 (~1660줄)
├── docking.py                 # 워커: 도킹/이완/정밀화 (~250줄)
├── analysis.py                # 워커: 스코어링/분석 (~590줄)
├── common.py                  # 유틸리티: PyRosetta 초기화, Pose 변환 (~150줄)
├── config.ini                 # 기본 설정 (1K 모델, v2.0 필터 포함)
├── config_20k.ini             # 중간 설정 (20K 모델)
├── config_100k.ini            # 프로덕션 설정 (100K 모델, v2.0 필터 포함)
├── run_v1.pbs                 # PBS 스크립트 v1.0 (input_PDB/ 전체 순차 처리)
├── run_v2_test.pbs            # PBS 스크립트 v2.0 (단일 PDB 테스트)
├── check_filter_v2_compat.py  # v2.0 메트릭 호환성 점검 (25항목)
├── check_improvements_compat.py # v1.0 개선사항 호환성 점검
├── input_PDB/                 # 입력 PDB 파일 (2-chain 구조)
│   ├── EGFR_TH1.pdb
│   ├── EGFR_beta.pdb
│   ├── C-lobe_TH1.pdb
│   └── C-lobe_beta.pdb
├── relaxed_cache/             # FastRelax 결과 캐시 (전역, PDB별 재사용)
├── <PDB_NAME>/                # PDB별 출력 디렉토리 (자동 생성)
│   ├── filter_passed/         # 필터 통과 구조
│   ├── cluster_results/       # 클러스터 대표 구조
│   ├── final_result/          # 최종 랭킹 구조 + CSV + PyMOL 스크립트
│   ├── final_ranking.csv      # 최종 랭킹 테이블
│   ├── energy_funnel.png      # 에너지 Funnel Plot
│   └── 1_OVERVIEW_Clusters.pml
├── pipeline.log               # v1.0 진행 로그
├── pipeline_v2_test.log       # v2.0 테스트 진행 로그
└── worker_debug.log           # 워커 프로세스 디버그 로그
```

---

## 3. Quick Start

### 3.1 서버 배포 전 호환성 점검

```bash
conda activate pyrosetta

# v2.0 필터 메트릭 호환성 (25항목 점검)
python check_filter_v2_compat.py

# v1.0 개선사항 호환성
python check_improvements_compat.py
```

점검 결과가 `[OK]` 또는 `[WARN]`이면 실행 가능. `[FAIL]`이 있으면 해당 기능을 비활성화하거나 코드 수정 필요.

### 3.2 PBS 배치 실행

```bash
# v2.0 테스트 (단일 PDB, 기본: C-lobe_beta.pdb)
qsub run_v2_test.pbs

# v2.0 테스트 (PDB/설정 지정)
qsub -v CONFIG_FILE=config_100k.ini,INPUT_PDB=input_PDB/EGFR_TH1.pdb run_v2_test.pbs

# v1.0 (input_PDB/ 내 모든 PDB 순차 처리)
qsub run_v1.pbs
qsub -v CONFIG_FILE=config_100k.ini run_v1.pbs
```

### 3.3 직접 실행

```bash
conda activate pyrosetta

# 단일 PDB
python pipeline_manager.py config.ini input_PDB/EGFR_TH1.pdb

# config만 지정 (config.ini의 input_pdb_name 사용)
python pipeline_manager.py config.ini
```

### 3.4 로그 모니터링

```bash
tail -f pipeline_v2_test.log   # v2.0 테스트 진행 상황
tail -f pipeline.log           # v1.0 파이프라인 진행 상황
tail -f worker_debug.log       # 워커 프로세스 디버그
```

---

## 4. Pipeline Steps

### Step 1: Relax (구조 이완)

입력 PDB에 `FastRelax` (ref2015 scorefxn) 적용.
결과는 `relaxed_cache/`에 캐싱되어 동일 PDB 재실행 시 생략된다.

### Step 2: Global Blind Docking

각 워커 프로세스에서 다음을 수행:

1. **FoldTree Setup** - `setup_foldtree(pose, "A_B", ...)` 로 도킹 점프 설정
2. **Randomization** - `RigidBodyPerturbMover(360°, 100Å)` 로 Chain B를 완전 랜덤 위치로 이동
3. **Slide Into Contact** - `DockingSlideIntoContact` 로 두 체인을 접촉 거리까지 이동
4. **Early Rejection** (선택) - 금지구역 접촉 시 `DockMCMProtocol` **이전에** 거부 (연산 절약)
5. **High-Res Docking** - `DockMCMProtocol` (ref2015) 로 고해상도 도킹 수행

> **핵심 주의사항**: `DockingSlideIntoContact`이 없으면 체인이 100Å 떨어진 상태에서
> `DockMCMProtocol`이 접촉을 회복하지 못해 모든 에너지가 0.0으로 나온다 (V1.0 핵심 버그).

### Step 2.5: Scoring & Multi-Stage Filtering

이 단계는 config 파일에 `[FilterStage1]` 섹션의 존재 여부에 따라 자동 분기된다.

---

#### v2.0 경로 (`[FilterStage1]` 섹션 존재 시)

**Pass 1: Fast Scoring** (모든 decoy 대상)

`InterfaceAnalyzerMover`를 통해 각 구조의 기본 메트릭을 계산:
- `dG_separated`: 인터페이스 결합 에너지 (REU)
- `dSASA`: 인터페이스 표면적 (Å²)
- `sc_value`: Shape Complementarity (0~1)
- `total_score`: 전체 Rosetta 에너지 (scorefxn(pose) 반환값)
- Chain B 중심 좌표 (center_x/y/z)
- Constraint 메트릭 (excluded_contacts, key_contact_ratio)

비용: 기존 v1.0 Fast Scoring과 거의 동일 (~0% 오버헤드, total_score 캡처 1줄 추가)

**Stage 1: Coarse Energy Filter**

구조적으로 불량한 모델을 빠르게 제거:

```
1a. Excluded Residue Hard Filter (기존 유지)
    - 금지구역 접촉 수 > max_excluded_contacts → 제거

1b. Repulsive Interface Removal
    - dG_separated > max_dG_separated (기본: 0.0) → 즉시 제거
    - 반발적 인터페이스(양의 에너지)는 물리적으로 의미 없음

1c. Total Score Percentile Filter
    - total_score 기준 상위 total_score_percentile% (기본: 10%) 유지
    - 전체 에너지가 비정상적으로 높은 구조 제거
```

**Stage 2: Interface Quality Filter (Cheap Metrics)**

Pass 1에서 이미 계산된 메트릭으로 인터페이스 품질 평가:

```
2a. dSASA >= min_dSASA (기본: 800.0 Å²)
    - 생물학적으로 유의미한 인터페이스 최소 면적
    - 문헌 근거: 단백질-단백질 인터페이스 평균 ~1600 Å²

2b. dG_density <= max_dG_density (기본: -1.5)
    - dG_density = dG_separated / dSASA × 100
    - 에너지 밀도: 면적 대비 결합 에너지의 효율성
    - 넓지만 약한 인터페이스 vs 좁지만 강한 인터페이스 구분

2c. sc_value >= min_sc_value (기본: 0.65)
    - 기하학적 형태 상보성
    - 문헌 근거: 단백질-단백질 인터페이스 평균 0.70 (Lawrence & Colman 1993)
```

**Pass 2: Intermediate Scoring** (Stage 1 생존자만 대상)

비싼 메트릭을 Stage 1 생존자(전체의 ~5-10%)에만 계산하여 효율성 확보:

```python
# analysis.py - run_intermediate_scoring_task()
# IAM 설정: set_compute_packstat(True) → 비용이 높은 패킹 분석 활성화

출력 메트릭:
- packstat:          원자 패킹 밀도 (0~1)
- delta_unsatHbonds: 결합 시 새로 미충족된 수소결합 수
- nres_int:          인터페이스 잔기 수
- hbonds_int:        인터페이스 수소결합 수
```

> `enable_expensive_metrics = false`로 설정하면 Pass 2를 건너뛴다.

**Stage 2: Interface Quality Filter (Expensive Metrics)**

```
2d. packstat >= min_packstat (기본: 0.65)
    - 인터페이스 패킹 품질
    - 0.0으로 설정 시 비활성화

2e. delta_unsatHbonds <= max_delta_unsatHbonds (기본: 5)
    - 결합으로 인해 새로 미충족된 수소결합 수
    - 99로 설정 시 비활성화
    - 문헌 근거: 좋은 인터페이스는 기존 H-bond를 보존함

2f. nres_int >= min_nres_int (기본: 15)
    - 인터페이스에 참여하는 잔기 수
    - 0으로 설정 시 비활성화

2g. hbonds_int >= min_hbonds_int (기본: 1)
    - 인터페이스 수소결합 수
    - 0으로 설정 시 비활성화
    - 문헌 근거: 자연 PPI 인터페이스는 대부분 최소 1개 이상 H-bond 보유
```

**Graduated Fallback v2.0**

Stage 2 필터 후 생존자가 `min_survivors` 미만이면 단계적으로 완화:

| Level | 동작 | 유지되는 필터 |
|-------|------|--------------|
| 0 | 모든 필터 통과 (정상) | 전부 |
| 1 | 비싼 메트릭만 해제 | dSASA + sc + dG_density |
| 2 | sc + dG_density 해제, dSASA 50% 완화 | dSASA (완화) |
| 3 | 모든 Stage 2 해제 | 없음 (Stage 1 생존자에서 dG 상위 N개) |

---

#### v1.0 경로 (`[FilterStage1]` 섹션 없을 때)

기존 3단계 필터링 그대로 동작:

**필터 1 - dSASA**:
```
dSASA >= min_dSASA (기본: 100.0 Å²)
```

**필터 1.5 - 금지구역**:
```
excluded_contacts <= max_excluded_contacts
```

**필터 2 - dG Percentile**:
```
dG_cutoff = np.percentile(scores, filter_percentile)
```
- 분포 독립적: 비대칭/skewed 에너지 분포에서도 안정적으로 작동
- 기본값: 하위 15% (filter_percentile=15)

**필터 3 - sc_value**:
```
sc_value >= min_sc_value (기본: 0.3)
```

**v1.0 Fallback**: 생존자 < `min_survivors`이면:
1. sc_value 필터 완화
2. dG 필터 완화
3. 모든 필터 무시, dG 상위 N개

---

### Step 3: Full Scoring

필터 생존자 전체(~50-200개)에 대해 `run_scoring_task`를 호출하여 모든 메트릭을 확보한다.
기존 파이프라인에서 클러스터링 이후에 수행하던 Final Scoring을 **클러스터링 이전으로 이동**한 것이다.

계산 메트릭:
- dG_separated, dSASA, sc_value, packstat (모든 버전)
- delta_unsatHbonds, nres_int, hbonds_int (v2.0에서 추가)
- dG_density (full scoring의 dG/dSASA 기반 재계산)
- L_RMSD, B-factor injection, per-residue CSV, binding residues

이점:
- **클러스터링 입력에 모든 메트릭 포함** → cluster_summary.csv에 15+ 컬럼 기록
- **Refinement(~40% 연산) 제거** → 유력 포켓 탐색 우선 전략
- 필터 생존자 ~50-200개에만 적용하므로 연산 부담 제한적

> `pool.imap` (ordered) + `zip`으로 입출력 1:1 매핑을 유지한다.

### Step 4: L_RMSD Greedy Clustering

에너지 순 정렬된 후보군에 대해 Greedy Leader-Follower 알고리즘 적용:

1. 첫 번째 구조 = Cluster 1의 Leader
2. 다음 구조를 순서대로, 모든 기존 Leader와 L_RMSD를 계산
3. L_RMSD <= threshold (4.0Å) 이면 해당 클러스터에 배정
   - **Closest-match**: 가장 가까운 Leader에 배정 (첫 매칭 아님)
   - 배정 시, 기존 멤버들과 **Pairwise L_RMSD >= member_diversity_dist (2.0Å)** 확인
   - 조건 충족 시만 멤버 추가 (최대 members_per_cluster개)
4. 모든 Leader와 L_RMSD > threshold 이면 새 클러스터 생성 (최대 cluster_top_n개)

**L_RMSD 계산 방법**:
```
1. calpha_superimpose_pose(mobile, reference)  # Chain A 기준 정렬
2. Chain B의 CA 원자 RMSD 측정                  # Ligand의 위치 차이
```

**CoM Pre-filter**: Chain B 중심 좌표(Center of Mass) 거리로 먼 후보를 사전 제외하여 O(N²) RMSD 계산 회피.

> L_RMSD(Ligand RMSD)는 Rosetta 도킹 논문의 표준 클러스터링 메트릭이다
> (Gray 2003, Chaudhury 2011).

**cluster_summary.csv** (15+ 컬럼, Step 3의 Full Scoring 결과 포함):
```
Cluster_ID, Member_ID, dG_separated, Total_Score, dSASA, dG_density, sc_value,
packstat, delta_unsatHbonds, nres_int, hbonds_int, L_RMSD, Binding_Residues,
Population, File_PDB [, key_contact_ratio]
```

### Step 5: Diversity-Aware Selection & Save

1. **Round-robin Selection**: 클러스터 간 균등하게 구조 선택 (다양성 보장)
2. **L_RMSD Dedup**: 이미 선택된 구조와 L_RMSD < 2.0Å인 구조 제외
3. **Key Residue Bonus**: `adjusted_dG = dG - bonus_weight × key_contact_ratio`

저장:
- `final_result/Rank01_C01_M01_S-12.34.pdb`
- `final_result/Rank01_C01_M01_Energies.csv` (per-residue 에너지)
- `final_ranking.csv` (전체 랭킹 테이블)

**final_ranking.csv 컬럼** (v2.0):
```
Rank, File, dG_separated, Total_Score, adjusted_dG, dSASA, sc_value, packstat,
dG_density, delta_unsatHbonds, nres_int, hbonds_int,
L_RMSD, Cluster, Parent, excl_contacts, key_contact_ratio, Binding_Residues
```

### Step 6: Visualization & Validation Report

파이프라인 완료 시 다음 시각화/분석 파일이 자동 생성된다:

| 파일 | 위치 | 설명 |
|------|------|------|
| `1_OVERVIEW_Clusters.pml` | `root_dir/` | 전체 결합 사이트 비교. chain A=회색, chain B=클러스터별 색상. PyMOL 그룹으로 사이트별 on/off 가능 |
| `2_DETAIL_C##.pml` | `final_result/` | 개별 사이트 상세 뷰. B-factor 컬러링 (blue=안정, red=불안정) |
| `view_results.pml` | `final_result/` | 최종 랭킹 모델. B-factor 컬러링 |
| `energy_funnel.png` | `root_dir/` | L_RMSD vs dG 산점도. 퍼널 형태면 수렴, 평탄하면 비수렴 |
| `docking_validation_report.txt` | `final_result/` | PPI 사이트 탐색 리포트 (10개 품질 체크 C1~C10) |

**PML 특징**:
- **상대경로** 사용: 폴더째 로컬 복사 시 PyMOL에서 바로 열림
- **chain A 자동 정렬**: 모든 구조가 수용체 기준으로 정렬됨
- **`python`/`python end` 블록**: PyMOL의 줄단위 실행 특성에 맞게 Python 코드를 블록으로 감쌈

**실행 시 자동 정리**: `filter_passed/`, `cluster_results/`, `final_result/`는 파이프라인 시작 시 `shutil.rmtree`로 삭제 후 재생성. `relaxed_cache/`는 보존.

---

## 7. 결과 해석 가이드

### 7.1 핵심 개념: PPI 사이트 탐색

이 파이프라인은 **단일 최적 포즈 예측이 아닌, 후보 결합 사이트 탐색**이 목적이다.

- 각 클러스터 = 하나의 **후보 결합 포켓**
- 수렴도가 낮음 ≠ 실패 → 여러 후보 사이트가 발견된 것
- 높은 L_RMSD = 정상 → 표면 전체를 탐색하므로 서로 먼 사이트가 발견됨
- **생물학적 타당성은 사용자가 PyMOL에서 최종 판단**해야 함

### 7.2 결과 확인 순서

#### Step 1: Validation Report 확인

`docking_validation_report.txt`의 품질 체크 확인:

| 항목 | 의미 | 주의사항 |
|------|------|----------|
| C1_파이프라인 | 성공률 (조기거부 제외) | 너무 낮으면 입력 PDB 확인 |
| C2_결합에너지 | 최저 dG | < -10이면 PASS |
| C3_에너지퍼널 | L_RMSD vs dG 상관 | 음의 상관 = 좋은 퍼널 |
| C4_인터페이스크기 | dSASA 중앙값 | > 800 Å² 권장 |
| C5_형상상보성 | sc_value 중앙값 | > 0.65 권장 |
| C6_사이트탐색 | 사이트 다양성 | **항상 PASS** (다양한 사이트 = 정보적) |
| C7_샘플링규모 | 모델 수 | 10K+ 권장 |
| C8_dG밀도 | dG_density 중앙값 | < -1.5 권장 |
| C9_인터페이스잔기 | nres_int 중앙값 | > 15 권장 |
| C10_수소결합 | hbonds_int 중앙값 | ≥ 1 필수 |

#### Step 2: 정량적 메트릭 비교

`final_ranking.csv`와 `cluster_summary.csv`에서:
- 파일명 `Rank01_C01_M01_S-18.32.pdb`의 `S-18.32` = dG -18.32 REU
- dG 기준: < -10 의미있는 결합, -10~-5 약한 결합, > -5 무의미
- 단일 메트릭이 아닌 **dG + dSASA + sc + packstat 종합 판단** 필요

#### Step 3: PyMOL 시각적 분석 (가장 중요)

```bash
# 서버에서 폴더째 로컬 복사
scp -r eunae@node04:/path/to/C-lobe_beta/ ~/Desktop/C-lobe_beta/
cd ~/Desktop/C-lobe_beta/

# 1) 전체 사이트 분포 확인
pymol 1_OVERVIEW_Clusters.pml
# → chain A=회색, chain B=클러스터별 색상
# → PyMOL 패널에서 C01, C02... 그룹 on/off

# 2) 개별 사이트 상세 확인 (에너지 좋은 순서)
pymol final_result/2_DETAIL_C01.pml
# → 한 사이트의 포즈들이 같은 포켓에 모이면 수렴된 사이트

# 3) 전장 원본 PDB와 겹치기 (수동)
# PyMOL 내에서:
#   load /path/to/full_length.pdb
#   align full_length and chain A, C01_M01_S-18 and chain A
```

### 7.3 결과 파일별 용도

| 파일 | 용도 | 핵심 확인 사항 |
|------|------|---------------|
| `docking_validation_report.txt` | 자동 품질 평가 | C1~C10 PASS/FAIL, 전체 신뢰도 판정 |
| `final_ranking.csv` | 정량적 메트릭 | dG, dSASA, sc, packstat 수치 비교 |
| `cluster_summary.csv` | 사이트별 비교 | Population, 에너지, 메트릭 분포 |
| `1_OVERVIEW_Clusters.pml` | 전체 사이트 분포 | 어느 표면에 사이트가 모이는지 |
| `2_DETAIL_C##.pml` | 개별 사이트 수렴도 | 포즈들이 같은 포켓에 모이는지 |
| `view_results.pml` | 최종 구조 | B-factor 기반 잔기별 에너지 분포 |
| `energy_funnel.png` | 에너지 랜드스케이프 | 퍼널(수렴) vs 평탄(비수렴) |

---

## 8. Configuration Reference

### [Path]

| Parameter | Default | Description |
|-----------|---------|-------------|
| `input_pdb_name` | - | 기본 입력 PDB 경로 (CLI 인자로 override 가능) |

### [System]

| Parameter | Default | Description |
|-----------|---------|-------------|
| `n_cpus` | 32 | 사용 CPU 코어 수 (실제 코어 수와 min 적용) |

### [Docking]

| Parameter | Default | Description |
|-----------|---------|-------------|
| `total_global_models` | 1000 | 전역 도킹 샘플링 수. 프로덕션: 10K~100K 권장 |

### [Filter] (v1.0 레거시, v2.0에서도 일부 참조)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `filter_percentile` | 15 | dG 백분위 컷오프 (하위 N%) |
| `min_dSASA` | 100.0 | 최소 인터페이스 dSASA (Å²) |
| `contact_distance` | 6.0 | 결합 잔기 탐지 거리 (Å). 6.0 = PPI 표준 |
| `min_sc_value` | 0.3 | 최소 Shape Complementarity |
| `min_survivors` | 50 | 최소 생존 구조 수 (fallback 트리거) |

### [FilterStage1] (v2.0 신규 - 선택 사항)

이 섹션이 **존재하면** v2.0 필터링 경로가 자동 활성화됨. 섹션이 **없으면** v1.0 경로 사용.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `total_score_percentile` | 10 | total_score 상위 N% 유지 (Coarse filter) |
| `max_dG_separated` | 0.0 | dG > 이 값인 구조 즉시 제거 (반발적 인터페이스) |
| `enabled` | true | Stage 1 활성화 여부 (false면 Stage 1 건너뜀) |

### [FilterStage2] (v2.0 신규 - 선택 사항)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `min_dSASA` | 800.0 | 최소 인터페이스 면적 (Å²). 문헌 권장: 800+ |
| `max_dG_density` | -1.5 | 최대 에너지 밀도 (dG/dSASA×100). 낮을수록 좋음 |
| `min_sc_value` | 0.65 | 최소 Shape Complementarity. 문헌 권장: 0.65+ |
| `min_packstat` | 0.65 | 최소 패킹 품질. 0.0 = 비활성 |
| `max_delta_unsatHbonds` | 5 | 최대 미충족 수소결합 수. 99 = 비활성 |
| `min_nres_int` | 15 | 최소 인터페이스 잔기 수. 0 = 비활성 |
| `min_hbonds_int` | 1 | 최소 인터페이스 수소결합 수. 0 = 비활성 |
| `enable_expensive_metrics` | true | Pass 2 (packstat, unsatHb, nres, hbonds) 활성화 |
| `min_survivors` | 50 | Fallback 트리거 최소 생존자 수 |

### [Cluster]

| Parameter | Default | Description |
|-----------|---------|-------------|
| `threshold` | 4.0 | L_RMSD 클러스터링 임계값 (Å) |
| `cluster_top_n` | 20 | 최대 클러스터 수 |
| `members_per_cluster` | 5 | 클러스터당 최대 멤버 수 |
| `member_diversity_dist` | 2.0 | 멤버 간 최소 L_RMSD (Å) |

### [Refinement]

| Parameter | Default | Description |
|-----------|---------|-------------|
| `refine_per_struct` | 10 | 대표 구조당 정밀화 반복 횟수 |
| `perturb_trans` | 0.1 | 정밀화 병진 섭동 크기 (Å) |
| `perturb_rot` | 1.0 | 정밀화 회전 섭동 크기 (°) |

### [Constraints] (선택 사항)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `excluded_residues_A` | (비어있으면 비활성) | Chain A 금지 잔기. 예: `700-807,840-854` |
| `key_residues_B` | (비어있으면 비활성) | Chain B 핵심 잔기 (soft bonus) |
| `enable_early_rejection` | false | Global Docking 중 DockMCMProtocol 이전 조기 거부 |
| `exclusion_contact_dist` | 10.0 | 금지구역 접촉 판정 거리 (Å) |
| `max_excluded_contacts` | 0 | 허용 금지잔기 접촉 수 (0=전면 거부) |
| `key_residue_bonus_weight` | 0.0 | 핵심잔기 접촉 보너스 가중치 (0.0=비활성) |

### [Output]

| Parameter | Default | Description |
|-----------|---------|-------------|
| `save_filter_max` | 200 | 필터 단계 저장 구조 수 |
| `save_top_n` | 20 | 최종 결과 저장 수 |

### 설정 파일 비교

| 항목 | config.ini | config_20k.ini | config_100k.ini |
|------|-----------|---------------|----------------|
| `total_global_models` | 1,000 | 20,000 | 100,000 |
| `[FilterStage1/2]` | 포함 (v2.0) | 미포함 (v1.0) | 포함 (v2.0) |
| 용도 | 테스트/디버깅 | 논문용 | 대규모 프로덕션 |

---

## 9. Module Reference

### pipeline_manager.py (Orchestrator, ~1660줄)

파이프라인 전체 흐름을 제어하는 메인 모듈.

| Method | Description |
|--------|-------------|
| `__init__(config_file, override_input_pdb)` | 설정 로드, 디렉토리 생성 |
| `_setup_logging()` | 로거 2개 설정 (progress + internal) |
| `_validate_config()` | 설정값 파싱 (v1.0 + v2.0 파라미터), 유효성 검사 |
| `_compute_l_rmsd(mobile, ref)` | L_RMSD 계산 (클러스터링용) |
| `execute()` | 6단계 파이프라인 실행 (v2.0/v1.0 자동 분기 포함) |
| `generate_validation_report()` | 자동 품질 평가 (v2.0: C1~C10, v1.0: C1~C7) |
| `generate_pymol_script()` | PyMOL .pml 생성 |
| `plot_energy_funnel()` | L_RMSD vs dG 산점도 생성 |
| `generate_cluster_overview()` | 1_OVERVIEW_Clusters.pml |
| `generate_per_cluster_views()` | 2_DETAIL_C*.pml |

**v2.0 관련 주요 속성**:
```python
self.filter_v2_enabled     # bool: [FilterStage1] 섹션 존재 여부
self.s1_enabled            # bool: Stage 1 활성화 여부
self.s1_total_score_pct    # int: total_score 백분위 (기본: 10)
self.s1_max_dG             # float: 최대 dG (기본: 0.0)
self.s2_min_dSASA          # float: 최소 dSASA (기본: 800.0)
self.s2_max_dG_density     # float: 최대 dG_density (기본: -1.5)
self.s2_min_sc             # float: 최소 sc_value (기본: 0.65)
self.s2_min_packstat       # float: 최소 packstat (기본: 0.65)
self.s2_max_unsatHb        # int: 최대 unsatHbonds (기본: 5)
self.s2_min_nres           # int: 최소 nres_int (기본: 15)
self.s2_min_hbonds_int     # int: 최소 hbonds_int (기본: 1)
self.s2_enable_expensive   # bool: Pass 2 활성화 (기본: True)
self.s2_min_survivors      # int: Fallback 트리거 (기본: 50)
```

### docking.py (Worker: Sampling, ~250줄)

실제 Rosetta 프로토콜을 수행하는 워커 함수 모음.

| Function | Input | Output |
|----------|-------|--------|
| `run_relax_task(pdb_path)` | PDB 파일 경로 | `{status, pdb_data}` |
| `run_global_docking_task(args)` | `(idx, pdb_str, [excl_res, excl_dist, max_contacts])` | `{status, id, pdb_data, center_x/y/z}` 또는 `{status:"rejected"}` |
| `run_refinement_task(args)` | `(idx, pdb_str, parent, trans, rot, sub_idx)` | `{status, id, parent, pdb_data}` |
| `_check_excluded_contact(pose, excl, dist)` | Pose, 금지잔기 set, 거리 | 접촉 수 (int) |

Global Docking 프로토콜:
```
setup_foldtree → RigidBodyPerturbMover(360°,100Å) → DockingSlideIntoContact
→ [Early Rejection: _check_excluded_contact] → DockMCMProtocol
```

### analysis.py (Worker: Scoring, ~590줄)

인터페이스 에너지 및 구조 분석.

| Function | Description |
|----------|-------------|
| `run_fast_scoring_task(args)` | **Pass 1**: 빠른 스코어링 (dG, dSASA, sc, total_score, center, constraint metrics) |
| `run_intermediate_scoring_task(args)` | **Pass 2 (v2.0 신규)**: 비싼 메트릭 (packstat, delta_unsatHbonds, nres_int, hbonds_int) |
| `run_scoring_task(args)` | **Final**: 전체 스코어링 (위 전부 + L_RMSD, B-factor, per-residue CSV, binding residues) |
| `inject_bfactors(pose)` | Residue 에너지를 B-factor에 주입 (모든 원자 순회) |
| `get_residue_energy_csv_string(pose, sfxn)` | Per-residue 에너지 CSV 생성 (fa_atr, fa_rep, fa_sol, fa_elec) |
| `analyze_interface_contacts(pose, dist)` | NeighborhoodResidueSelector로 Chain A 인터페이스 잔기 탐지 |
| `calculate_l_rmsd(mobile, ref)` | L_RMSD (Chain A 정렬 후 Chain B CA RMSD) |
| `count_excluded_contacts(pose, excl, dist)` | 금지구역 접촉 수 (스코어링 단계용) |
| `calculate_key_contact_ratio(pose, key_B, dist)` | 핵심잔기 접촉 비율 (soft metric) |
| `get_chain_b_center(pose)` | Chain B 기하학적 중심 좌표 |

**run_fast_scoring_task 출력** (v2.0):
```python
{
    "status": "success",
    "dG_separated": float,     # 인터페이스 결합 에너지
    "total_score": float,      # 전체 Rosetta 에너지 (v2.0 추가)
    "dSASA": float,            # 인터페이스 표면적
    "sc_value": float,         # Shape Complementarity
    "center_x/y/z": float,    # Chain B 중심 좌표
    "excluded_contacts": int,  # 금지구역 접촉 수
    "key_contact_ratio": float # 핵심잔기 접촉 비율
}
```

**run_intermediate_scoring_task 출력** (v2.0 신규):
```python
{
    "status": "success",
    "packstat": float,             # 패킹 밀도 (0~1), 실패 시 0.0
    "delta_unsatHbonds": int,      # 미충족 수소결합, 실패 시 99
    "nres_int": int,               # 인터페이스 잔기 수, 실패 시 0
    "hbonds_int": int              # 인터페이스 수소결합 수, 실패 시 0
}
```

### common.py (Utility, ~150줄)

| Function | Description |
|----------|-------------|
| `init_rosetta()` | PyRosetta 초기화 (멱등, PID 기반 랜덤 시드, 배너 억제) |
| `setup_worker_logging()` | 워커 프로세스 로깅 핸들러 설정 (멱등) |
| `string_to_pose(pdb_string)` | PDB 문자열 → Pose 변환 (실패 시 None) |
| `pose_to_string(pose)` | Pose → PDB 문자열 변환 |
| `parse_residue_ranges(range_str)` | `"700-807,840-854"` → `{700,701,...,807,840,...,854}` |

---

## 10. Scoring Metrics Reference

### 기본 메트릭 (모든 버전)

| Metric | Excellent | Good | Marginal | Poor | 설명 |
|--------|-----------|------|----------|------|------|
| dG_separated (REU) | < -15 | -15 ~ -5 | -5 ~ 0 | > 0 | 인터페이스 결합 에너지. 낮을수록 강한 결합 |
| dSASA (Å²) | > 1200 | 800~1200 | 500~800 | < 500 | 결합 시 묻히는 표면적 |
| sc_value | > 0.70 | 0.65~0.70 | 0.50~0.65 | < 0.50 | Shape Complementarity |
| packstat | > 0.70 | 0.65~0.70 | 0.55~0.65 | < 0.55 | 원자 패킹 밀도 |
| L_RMSD (Å) | 낮을수록 | - | - | - | Relaxed 대비 Ligand RMSD |

### v2.0 추가 메트릭

| Metric | Excellent | Good | Marginal | Poor | 설명 |
|--------|-----------|------|----------|------|------|
| total_score (REU) | 낮을수록 | - | - | - | 전체 Rosetta 에너지 |
| dG_density | < -2.0 | -2.0 ~ -1.5 | -1.5 ~ -1.0 | > -1.0 | dG/dSASA×100. 면적 대비 에너지 효율 |
| delta_unsatHbonds | 0~2 | 3~5 | 6~10 | > 10 | 결합 시 새로 미충족된 수소결합 |
| nres_int | > 25 | 15~25 | 10~15 | < 10 | 인터페이스 잔기 수 |
| hbonds_int | > 5 | 2~5 | 1 | 0 | 인터페이스 수소결합 수 |

### Constraint 메트릭

| Metric | 설명 |
|--------|------|
| excluded_contacts | 금지구역 잔기 접촉 수. 0이어야 통과 |
| key_contact_ratio | 핵심 잔기 접촉 비율 (0~1). 높을수록 좋음 |
| adjusted_dG | dG - bonus_weight × key_contact_ratio. 핵심잔기 보너스 적용 에너지 |

---

## 11. Filtering Reliability & Scientific Basis

본 섹션은 필터링 전략의 과학적 신뢰성을 논문 수준으로 기술한다.

### 8.1 Overview

Global blind docking은 단백질 표면 전체를 탐색하므로, 수천~수만 개의 decoy 중 물리화학적으로 무의미한 포즈가 대다수를 차지한다. 본 파이프라인은 8개의 독립적 인터페이스 품질 메트릭을 계층적으로 적용하는 **2-Pass, 2-Stage 필터링** 전략을 사용하여, 계산 비용을 최소화하면서 신뢰도 높은 후보만을 클러스터링 단계로 전달한다.

### 8.2 2-Pass Design Rationale

필터링 메트릭은 계산 비용에 따라 두 그룹으로 나뉜다:

| Pass | 대상 | 메트릭 | 계산 비용 | PyRosetta API |
|------|------|--------|----------|---------------|
| **Pass 1** (Fast) | 전체 decoy (N) | dG_separated, dSASA, sc_value, total_score | Low | `InterfaceAnalyzerMover` (packstat=**OFF**) |
| **Pass 2** (Expensive) | Stage 1 생존자만 (~5-10% N) | packstat, delta_unsatHbonds, nres_int, hbonds_int | High | `InterfaceAnalyzerMover` (packstat=**ON**) |

이 설계의 핵심은 비용이 큰 메트릭(packstat 등)을 전체 decoy가 아닌 1차 생존자에만 적용하여, 100K decoy 기준 expensive metric 연산량을 **90-95% 절감**하는 것이다.

### 8.3 Stage 1: Coarse Energy Filter

Stage 1은 에너지적으로 비합리적인 decoy를 조기 제거한다:

1. **Excluded zone hard filter**: 생물학적으로 불가능한 결합 부위(막면, 다이머 인터페이스 등)에 접촉하는 포즈를 제거
2. **dG_separated ≤ 0 REU**: 양의 결합 에너지(반발적 인터페이스)를 가진 decoy 제거
3. **total_score percentile (상위 10%)**: Rosetta 전체 에너지 기준 구조적으로 불량한 decoy 제거

Stage 1은 일반적으로 전체 decoy의 5-10%만 생존시키며, 이 단계까지의 계산 비용은 사실상 0이다 (Pass 1에서 이미 확보된 값 사용).

### 8.4 Stage 2: Interface Quality Filter (Literature-Based Thresholds)

Stage 2는 인터페이스의 물리화학적 품질을 8개 독립 메트릭으로 평가한다. 각 임계값은 Rosetta docking 문헌에 기반한다.

#### Cheap Metrics (Pass 1에서 이미 계산)

| Metric | Threshold | Biophysical Meaning | Literature Basis |
|--------|-----------|---------------------|-----------------|
| **dSASA** ≥ 800 Å² | 결합 시 용매로부터 묻히는 표면적 | 천연 PPI 인터페이스 평균 ~1,600 Å² (Lo Conte et al., 1999). 800 Å²는 생물학적으로 유의미한 인터페이스의 하한 |
| **dG_density** ≤ −1.5 | dG/dSASA × 100; 단위 면적당 결합 에너지 효율 | 넓지만 에너지적으로 비효율적인 인터페이스 배제. Rosetta docking 벤치마크에서 네이티브 복합체 평균 −1.5~−2.5 범위 (Chaudhury et al., 2011) |
| **sc_value** ≥ 0.65 | Shape Complementarity (0-1); 두 표면의 기하학적 맞물림 정도 | Lawrence & Colman (1993)의 원 정의에 따르면 천연 PPI 평균 0.70, 항체-항원 0.64-0.68 |

#### Expensive Metrics (Pass 2에서 계산)

| Metric | Threshold | Biophysical Meaning | Literature Basis |
|--------|-----------|---------------------|-----------------|
| **packstat** ≥ 0.65 | 인터페이스 원자 패킹 밀도 (0-1) | Sheffler & Baker (2009)의 RosettaHoles 알고리즘 기반. 잘 패킹된 인터페이스는 >0.65 |
| **delta_unsatHbonds** ≤ 5 | 결합 시 새로 발생하는 미충족 수소결합 수 | 미충족 극성 원자는 결합의 열역학적 페널티를 의미. 천연 복합체에서 일반적으로 <5 (Kortemme et al., 2003) |
| **nres_int** ≥ 15 | 인터페이스에 참여하는 잔기 수 | 안정적 PPI는 최소 15-20개 잔기가 인터페이스를 형성 (Chakrabarti & Janin, 2002) |
| **hbonds_int** ≥ 1 | 인터페이스 간 수소결합 수 | 생물학적 PPI는 거의 예외 없이 최소 1개 이상의 체인 간 수소결합을 가짐 (Xu et al., 1997) |

### 8.5 Conservative Default Strategy

PyRosetta 버전 간 API 호환성 문제를 처리하기 위해, 모든 메트릭 추출은 **3-tier fallback** 전략을 사용한다:

```
1st: getPoseExtraScore(pose, key)    — 가장 안정적
2nd: IAM method (hasattr guard)      — 버전 의존적
3rd: Conservative default            — 안전 기본값
```

기본값은 모두 필터 **탈락 방향**으로 설정된다:
- packstat → 0.0 (임계값 0.65 미달)
- delta_unsatHbonds → 99 (임계값 5 초과)
- nres_int → 0 (임계값 15 미달)
- hbonds_int → 0 (임계값 1 미달)

이는 메트릭 추출 실패가 거짓 양성(false positive)으로 이어지지 않음을 보장한다.

### 8.6 Graduated Fallback System

엄격한 다중 필터는 과소 생존(under-survival) 위험이 있다. 이를 방지하기 위해 4단계 점진적 완화(Graduated Fallback)를 적용한다:

| Level | 조건 | 완화 내용 | 유지되는 필터 |
|-------|------|-----------|--------------|
| **0** (정상) | 생존자 ≥ min_survivors | 완화 없음 | 전체 8개 메트릭 |
| **1** | Stage 2 Full 부족 | Expensive 메트릭만 해제 | dSASA + sc + dG_density |
| **2** | Stage 2 Cheap 부족 | sc, dG_density 추가 해제, dSASA 50% 완화 | dSASA only (400 Å²) |
| **3** | Stage 1 부족 | 모든 Stage 2 해제 | dG 상위 N개만 선택 |

이 시스템은 항상 최소 `min_survivors`(기본 50)개의 후보를 보장하면서, 가능한 한 가장 엄격한 수준의 필터를 적용한다. Fallback 발동 여부와 레벨은 Validation Report에 명시적으로 기록된다.

### 8.7 Filtering Reliability의 근거

본 필터링 전략의 신뢰성은 다음 네 가지 설계 원칙에 기반한다:

1. **다중 독립 메트릭**: 단일 메트릭의 한계를 8개 직교적 메트릭의 교차 검증으로 보완. 에너지(dG), 기하학(sc, packstat), 크기(dSASA, nres_int), 극성 상보성(unsatHbonds, hbonds_int)의 네 차원을 모두 평가
2. **문헌 기반 임계값**: 모든 임계값은 천연 PPI 구조 통계 및 Rosetta docking 벤치마크에서 도출된 값을 사용
3. **보수적 실패 모드**: 메트릭 추출 실패 시 항상 탈락 방향의 기본값을 사용하여 false positive 방지
4. **Graduated Fallback**: 과도한 필터링으로 인한 정보 손실을 방지하면서, 데이터가 허용하는 최대한 엄격한 기준을 적용

### 8.8 References

- Chakrabarti, P. & Janin, J. (2002) Dissecting protein-protein recognition sites. *Proteins* 47, 334-343
- Chaudhury, S. et al. (2011) Benchmarking and analysis of protein docking performance in Rosetta v3.2. *PLoS ONE* 6, e22477
- Gray, J.J. et al. (2003) Protein-protein docking with simultaneous optimization of rigid-body displacement and side-chain conformations. *J. Mol. Biol.* 331, 281-299
- Kortemme, T. et al. (2003) Computational alanine scanning of protein-protein interfaces. *Sci. STKE* 2003, pl2
- Lawrence, M.C. & Colman, P.M. (1993) Shape complementarity at protein/protein interfaces. *J. Mol. Biol.* 234, 946-950
- Lo Conte, L. et al. (1999) The atomic structure of protein-protein recognition sites. *J. Mol. Biol.* 285, 2177-2198
- Sheffler, W. & Baker, D. (2009) RosettaHoles: rapid assessment of protein core packing for structure prediction. *Protein Sci.* 18, 229-239
- Xu, D. et al. (1997) Hydrogen bonds and salt bridges across protein-protein interfaces. *Protein Eng.* 10, 999-1012

---

## 12. Validation Report

`docking_validation_report.txt`는 **PPI 사이트 탐색 관점**으로 자동 생성되는 품질 평가 리포트이다.

### v2.0 품질 체크 항목 (C1~C10)

| Check | 이름 | PASS | WARNING | FAIL |
|-------|------|------|---------|------|
| C1 | C1_파이프라인 | 유효 성공률 ≥ 기준 | marginal | 낮음 |
| C2 | C2_결합에너지 | dG < -10 | -10 ~ -5 | > -5 |
| C3 | C3_에너지퍼널 | 음의 상관 | 약한 상관 | 양의 상관 |
| C4 | C4_인터페이스크기 | dSASA > 800 | 500~800 | < 500 |
| C5 | C5_형상상보성 | sc > 0.65 | 0.50~0.65 | < 0.50 |
| C6 | C6_사이트탐색 | **항상 PASS** | - | - |
| C7 | C7_샘플링규모 | ≥ 10K | 1K~10K | < 1K |
| C8 | C8_dG밀도 | < -1.5 | -1.5 ~ -1.0 | > -1.0 |
| C9 | C9_인터페이스잔기 | > 15 | 10~15 | < 10 |
| C10 | C10_수소결합 | ≥ 1 | - | 0 |

> - v1.0은 C1~C7만 포함.
> - **C1**: 조기거부(early rejection) 제외한 유효 성공률 사용.
> - **C6**: PPI 사이트 탐색에서 다양한 사이트 발견은 정상이므로 항상 PASS.

### 리포트 구성

```
=== PPI 결합 사이트 탐색 리포트 ===
  Global Blind Docking Pipeline (필터 v2.0)

1. 실행 요약 (입력, 샘플링, 필터 통계)
2. 에너지 분포 (dG, dSASA, sc, packstat 등)
3. 새 메트릭 분포 (v2.0: dG_density, unsatHb, nres, hbonds)
4. 결합 사이트 탐색 분석
   - 발견 사이트 수, 분포, Boltzmann 확률
   - 사이트 탐색 결과 판정 (지배적/유력/N개 후보)
   - L_RMSD 분포 (높은 L_RMSD = 정상)
   - 결합 핫스팟 잔기
5. 품질 체크 (C1~C10)
6. 종합 판정 (높은/보통/낮은 신뢰도 / 신뢰 불가)
7. 추천사항 + PyMOL 후속 분석 가이드
```

---

## 13. Logging System

| File | Content | Writer |
|------|---------|--------|
| `pipeline.log` / `pipeline_v2_test.log` | 파이프라인 진행 상황 (단계별 시간, 필터 통계, 클러스터 수) | PBS가 stdout 캡처 |
| `worker_debug.log` | 워커 프로세스 디버그 로그 (PID, 에러 traceback) | 파일 핸들러 직접 기록 |

- 메인 프로세스: `progress_logger` (stdout) + `internal_logger` (file)
- 워커 프로세스: `internal_logger` (file only, `setup_worker_logging()`으로 핸들러 설정)

v2.0 필터링 로그 예시:
```
[Step 2.5] v2.0 Multi-Stage Filtering
  Stage 1: Excluded filter: 1000 → 850
  Stage 1: dG > 0.0 removal: 850 → 720
  Stage 1: Total score top 10%: 720 → 72
  Stage 2 (cheap): dSASA >= 800: 72 → 65
  Stage 2 (cheap): dG_density <= -1.5: 65 → 58
  Stage 2 (cheap): sc >= 0.65: 58 → 52
  Pass 2: Intermediate scoring on 52 structures...
  Stage 2 (expensive): packstat >= 0.65: 52 → 48
  Stage 2 (expensive): unsatHb <= 5: 48 → 45
  Stage 2 (expensive): nres_int >= 15: 45 → 42
  Stage 2 (expensive): hbonds_int >= 1: 42 → 40
  Graduated Fallback: Level 0 (all filters passed, 40 survivors)
```

---

## 14. PyRosetta Compatibility

서버 환경은 네트워크가 차단되어 패키지 업데이트가 불가능하다.
아래 호환성 대응이 코드에 반영되어 있다.

### 메트릭 추출 Fallback 체인

각 메트릭은 3단계 fallback으로 추출:

| Metric | 1st: getPoseExtraScore | 2nd: IAM method | 3rd: Safe Default |
|--------|----------------------|-----------------|-------------------|
| sc_value | `"sc_value"` | `iam.get_sc_value()` | 0.0 |
| packstat | `"packstat"` | `iam.get_packstat()` → `iam.get_interface_packstat()` | 0.0 |
| delta_unsatHbonds | `"delta_unsatHbonds"` | `iam.get_interface_delta_hbond_unsat()` | 99 |
| nres_int | `"nres_int"` | `iam.get_num_interface_residues()` | 0 |
| hbonds_int | `"hbonds_int"` | `iam.get_interface_hbonds()` | 0 |

> **안전 기본값 설계**: 모든 기본값은 필터에서 **탈락** 방향으로 설정 (보수적).
> Graduated Fallback이 안전망 역할을 하여 추출 실패가 파이프라인 중단을 유발하지 않음.

### 기타 호환성 대응

| Issue | Solution | Location |
|-------|----------|----------|
| `residue_total_energies` 미지원 | `hasattr` 체크 후 조건부 호출 | analysis.py |
| `pyrosetta.init()` 중복 호출 | `RuntimeError` catch + "already initialized" 확인 | common.py |
| `-constant_seed` 문제 | `-jran {PID}` 로 프로세스별 고유 시드 | common.py |
| PyRosetta 배너 출력 | stdout/stderr 임시 redirect | common.py |

### 호환성 점검 스크립트

`check_filter_v2_compat.py`는 25개 항목을 점검:

| 번호 | 점검 항목 | 중요도 |
|------|-----------|--------|
| 1 | delta_unsatHbonds (getPoseExtraScore) | 필수 |
| 2 | delta_unsatHbonds (IAM method) | fallback |
| 3 | nres_int (getPoseExtraScore) | 필수 |
| 3.5 | hbonds_int (getPoseExtraScore) | 필수 |
| 4 | nres_int (IAM method) | fallback |
| 5 | packstat (getPoseExtraScore) | 필수 |
| 6 | packstat (get_packstat method) | fallback |
| 6.5 | packstat (get_interface_packstat method) | fallback |
| 7 | scorefxn(pose) 반환값 타입 | 필수 |
| 8 | dG_density 계산 검증 | 필수 |
| 9 | packstat 성능 벤치마크 | 정보 |
| 10 | sklearn 가용성 | 선택 |
| 11 | scipy 가용성 | 선택 |

서버 실행 결과: **25 OK / 3 WARNING / 0 FAIL** (2025.02 기준)

---

## 15. Multiprocessing Patterns

| Step | Pool Method | 순서 보장 | 이유 |
|------|-------------|-----------|------|
| Step 2 (Global Docking) | `pool.imap_unordered` | 불필요 | 최대 처리량, 결과를 리스트에 append |
| Step 2.5 Pass 1 (Fast Scoring) | `pool.imap` + enumerate | **필수** | 입출력 1:1 매핑 (docking_results[i]와 대응) |
| Step 2.5 Pass 2 (Intermediate) | `pool.imap` + enumerate | **필수** | Stage 1 생존자 인덱스와 대응 |
| Step 3 (Full Scoring) | `pool.imap` + enumerate | **필수** | 필터 생존자 인덱스와 대응 |

> **절대 주의**: `pool.imap`을 `pool.imap_unordered`로 바꾸면 인덱스가 깨져서
> 잘못된 구조에 잘못된 스코어가 배정된다.

### 메모리 관리

대규모 실행 시 PDB 문자열이 수만 개 메모리에 올라감:
```python
docking_results[i] = None  # 사용 완료 즉시 해제
del variable                # 참조 제거
gc.collect()               # 가비지 컬렉션 강제 실행
```

---

## 16. Estimated Runtime (32 cores)

| Config | Docking | Pass 1 Scoring | Pass 2 Scoring | Full Scoring | Clustering | Total |
|--------|---------|----------------|----------------|--------------|------------|-------|
| 1K 모델 | ~1h | ~20min | ~5min | ~10min | ~5min | ~2h |
| 20K 모델 | ~12-24h | ~4h | ~30min | ~20min | ~30min | ~18-30h |
| 100K 모델 | ~60-80h | ~15h | ~2h | ~1h | ~1h | ~80-100h |

- Pass 2 오버헤드: Stage 1 생존자만 대상 (~5-10%) → 전체 대비 +15~25%
- Full Scoring: 필터 생존자 ~50-200개 대상 → 소요 시간 제한적
- Clustering bottleneck: `string_to_pose` (~160ms/op) + `calpha_superimpose_pose` (~3ms/op)
- Refinement 제거로 기존 대비 ~40% 연산 절감

---

## 17. Development History

| Commit | Description |
|--------|-------------|
| **v2.0 필터링 개편** | |
| `b39ca4c` | hbonds_int 메트릭 추가 + 임계값 문헌 권장값으로 강화 |
| `d6aa898` | get_interface_packstat 3단계 fallback 수정 |
| `a47955b` | 호환성 체크 test pose 체인 분리 수정 |
| `ba598fd` | **v2.0 다중 단계 필터링 파이프라인 대개편** |
| `c2e816a` | v1.0-stable 백업 태그 (필터링 개편 전) |
| **v1.0** | |
| `ae91d63` | Tier 1-3 종합 개선 (Constraints, Early Rejection, Selection, Dedup) |
| `10c2576` | Validation Report (자동 품질 평가) |
| `169e170` | Per-residue 에너지 분석 + Summary CSV |
| `aa9a4a4` | L_RMSD 클러스터링 + 다중 기준 필터링 구현 |
| `6c05508` | Contact residue cutoff 4.0Å → 6.0Å |
| `4f15e33` | SlideIntoContact 누락 버그 수정 (핵심 V1.0 fix) |

---

## 18. Troubleshooting

### 모든 dG가 0.0으로 나옴
→ `DockingSlideIntoContact` 누락. docking.py에서 `DockMCMProtocol` 이전에 반드시 호출.

### 필터 후 생존자가 0개
→ Graduated Fallback이 자동 작동하여 최소 `min_survivors`개 확보. 로그에 Fallback Level 표시.
→ 임계값이 너무 엄격할 수 있음. `[FilterStage2]` 값 완화 고려.

### Pass 2 메트릭이 모두 기본값 (packstat=0.0, unsatHb=99 등)
→ `getPoseExtraScore`와 IAM method 모두 실패. `check_filter_v2_compat.py` 실행하여 PyRosetta 호환성 확인.
→ `enable_expensive_metrics = false`로 설정하면 Pass 2 건너뛰기 가능.

### v2.0 필터링이 작동하지 않음 (v1.0 경로로 실행됨)
→ `config.ini`에 `[FilterStage1]` 섹션이 있는지 확인. 이 섹션의 존재가 v2.0 활성화 트리거.

### 메모리 부족 (100K 모델)
→ `docking_results[i] = None` 패턴이 정상 작동하는지 확인.
→ `save_filter_max` 값을 줄여 디스크 I/O 감소.

### PBS 작업 시간 초과
→ `run_v2_test.pbs`의 walltime 확인 (기본 12시간).
→ 100K 모델은 `walltime=168:00:00` (7일) 이상 필요.

### 호환성 점검 FAIL
→ `check_filter_v2_compat.py`에서 `[FAIL]` 항목 확인.
→ getPoseExtraScore가 실패하면 IAM method fallback 사용.
→ 모든 방법이 실패하면 해당 메트릭의 필터를 비활성화 (예: `min_packstat = 0.0`).
