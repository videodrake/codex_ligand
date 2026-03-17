# EGFR-MYO1D Docking Pipeline — 전체 아키텍처 상세 보고서

> 이 보고서는 `egfr_pipeline/` 내 모든 모듈의 실제 소스코드를 분석하여 작성되었습니다.
> 각 Phase와 Task Group이 수행하는 구체적인 작업, 입출력, 알고리즘을 기술합니다.

---

## 목차

1. [프로젝트 개요](#1-프로젝트-개요)
2. [두 가지 실행 워크플로우](#2-두-가지-실행-워크플로우)
3. [Vina 모듈 — AutoDock Vina Blind Docking](#3-vina-모듈)
4. [PyRosetta Docking Core — PPI Global Docking 엔진](#4-pyrosetta-docking-core)
5. [Phase 1 — PPI-First Interface Mapping](#5-phase-1)
6. [Phase 2 — Pocket Analysis & Druggability Assessment](#6-phase-2)
7. [Phase 3 — Diversity-Aware Focused Docking](#7-phase-3)
8. [Phase 4 — Perturbation Relevance Scoring](#8-phase-4)
9. [Cross-Cutting 유틸리티](#9-cross-cutting-유틸리티)
10. [Step View / 출력 조직화](#10-step-view)
11. [전체 데이터 흐름 요약](#11-전체-데이터-흐름)

---

## 1. 프로젝트 개요

**목적**: EGFR 키나아제 도메인과 MYO1D beta-meander 간의 단백질-단백질 상호작용(PPI)을 교란할 수 있는 약물 결합 부위를 탐색하는 통합 파이프라인.

**핵심 전략**:
- 소분자 도킹(Vina)과 단백질 도킹(PyRosetta PPI)을 수행
- 두 증거원의 공간적 일치를 분석하여 PPI 교란 가능성이 높은 사이트를 우선순위 지정
- 3가지 EGFR 구조 상태(3GT8_raw, EGFR_160-185, EGFR_170-200)에 걸친 교차 상태 로버스트니스 분석

**실행 환경**: Linux HPC (PBS/qsub), 32 CPU cores, 네트워크 차단 환경. **모든 도킹/연산은 반드시 qsub를 통해 HPC 서버에서 실행**한다.

**진입점**: `main.py` — 대화형 CLI, `run_production.py` — 프로덕션 자동화 (lane 기반 PBS 제출)

> **중요**: 이 파이프라인에는 **두 가지 실행 워크플로우**가 있다. 다음 섹션에서 상세히 설명한다.

---

## 2. 두 가지 실행 워크플로우

이 파이프라인은 **두 가지 독립적인 워크플로우**로 구성된다. 보고서의 모든 모듈은 이 두 워크플로우 중 하나 또는 양쪽에 소속된다.

### Workflow A: Standard Production (자동화 완료)

> 진입점: `run_production.py` + lane별 PBS 스크립트
> 상태: **자동화 완료**, `qsub` 한 번으로 전체 실행 가능

Vina blind docking과 PPI blind docking을 **독립적으로** 수행한 뒤, 마지막에 결과를 통합하는 워크플로우.

```
Phase 1: Vina Blind Docking        (소분자 → EGFR 표면 전체)
Phase 2: PPI Global Blind Docking  (MYO1D → EGFR, 3 states × 5 seeds, 300K models)
Phase 3: PPI Postprocess           (인터페이스 잔기 추출 + 스코어 표준화)
Phase 4: Vina Postprocess          (파싱 → 접촉 → 클러스터링 → 요약 → 교차비교 → 부트스트랩)
Phase 5: Site Verdict              (Vina 포켓 × PPI 증거 3축 통합 판정)
Phase 6: Report                    (종합 보고서)
Phase 7: Validate                  (출력 검증)
```

**핵심 특징**: Phase 1(Vina)과 Phase 2(PPI)는 서로의 결과를 사용하지 않는다. 독립 증거원이며, Phase 5(Verdict)에서 처음으로 병합된다.

**`run_production.py` 전체 lane 목록**: `vina-cpu`, `ppi`, `ppi-post`, `vina-post`, `finalize`, `status`, `vina-gpu`, `phase3-gpu`, `adv-phase1`, `adv-phase2`, `adv-phase3-setup`, `adv-phase3-execute`, `adv-phase3-post`, `adv-phase4` (총 14개)

**서버 제출 방법**:
```bash
# 방법 1: 올인원 (순차 실행)
PRECHECK=$(qsub config/run_pre_qsub_checks.pbs)
qsub -W depend=afterok:${PRECHECK} config/run_production.pbs

# 방법 2: Lane별 병렬 제출 (권장 — Vina와 PPI 동시 실행)
PRECHECK=$(qsub config/run_pre_qsub_checks.pbs)
VINA=$(qsub -W depend=afterok:${PRECHECK} config/run_vina_cpu.pbs)
PPI_JOBS=""
for STATE in 3GT8_raw EGFR_160-185 EGFR_170-200; do
  for SEED in 0 1 2 3 4; do
    JOB=$(qsub -W depend=afterok:${PRECHECK} -v STATE=${STATE},SEED=${SEED} config/run_ppi_state_seed.pbs)
    PPI_JOBS="${PPI_JOBS}:${JOB}"
  done
done
VINA_POST=$(qsub -W depend=afterok:${VINA} config/run_vina_postprocess.pbs)
PPI_POST=$(qsub -W depend=afterok${PPI_JOBS} config/run_ppi_postprocess.pbs)
qsub -W depend=afterok:${VINA_POST}:${PPI_POST} config/run_finalize.pbs
```

### Workflow B: Advanced PPI-First Pipeline (자동화 완료)

> 진입점: `qsub config/run_advanced_pipeline.pbs` (전체 자동) 또는 `run_production.py --lane adv-*` (개별)
> 상태: **자동화 완료** — `run_production.py`에 `adv-*` lane 6개 통합, PBS 스크립트 완비

PPI 도킹 결과를 기반으로 포켓을 좁혀가며 순차적으로 분석하는 워크플로우. **각 Phase의 출력이 다음 Phase의 입력**이 되는 실질적 데이터 의존성이 있다.

```
Phase 1: PPI 도킹 + 인터페이스 분석 (TG 1.0~1.6)
    ↓ phase1_downstream_patch_reference.csv
Phase 2: Pocket Analysis (fpocket/P2Rank + PPI 관계 분석, TG 2.0~2.7)
    ↓ phase3_candidate_pocket_reference.csv
Phase 3: Focused Vina Docking (포켓 타겟, budget-aware, TG 3.0~3.6)
    ↓ phase4_docking_evidence_reference.csv
Phase 4: Perturbation Relevance Scoring (4축 통합, TG 4.0~4.6)
```

**전제 조건**: Workflow A의 Phase 2(PPI docking)가 완료되어야 시작 가능.

**서버 제출 방법**:
```bash
# 전체 자동 (PBS 의존성 체인으로 Phase 1~4 순차 제출)
qsub config/run_advanced_pipeline.pbs

# Phase 2부터 시작 (Phase 1 분석 이미 완료 시)
qsub -v ADV_FROM=2 config/run_advanced_pipeline.pbs

# Round 수 변경 (기본 3)
qsub -v ADV_ROUNDS=5 config/run_advanced_pipeline.pbs

# 개별 Phase 수동 제출
qsub config/run_adv_phase1.pbs
qsub -W depend=afterok:<job> config/run_adv_phase2.pbs
qsub -W depend=afterok:<job> config/run_adv_phase3_setup.pbs
qsub -v ROUND=0 -W depend=afterok:<job> config/run_adv_phase3_execute.pbs
qsub -W depend=afterok:<job> config/run_adv_phase3_post.pbs
qsub -W depend=afterok:<job> config/run_adv_phase4.pbs
```

**Lane별 PBS 스크립트**:

| Lane | PBS 스크립트 | CPU | 역할 |
|------|-------------|-----|------|
| `adv-phase1` | `run_adv_phase1.pbs` | 4 | PPI 분석 (TG 1.1.5~1.6) |
| `adv-phase2` | `run_adv_phase2.pbs` | 16 | Pocket cascade (TG 2.0~2.7) |
| `adv-phase3-setup` | `run_adv_phase3_setup.pbs` | 2 | 도킹 잡 계획 (TG 3.0~3.2) |
| `adv-phase3-execute` | `run_adv_phase3_execute.pbs` | 16 | Vina focused 실행 (라운드별) |
| `adv-phase3-post` | `run_adv_phase3_post.pbs` | 4 | 도킹 후 분석 (TG 3.4~3.6) |
| `adv-phase4` | `run_adv_phase4.pbs` | 2 | 통합 스코어링 (TG 4.0~4.6) |

### 워크플로우 비교

| | Workflow A (Standard) | Workflow B (Advanced PPI-First) |
|---|---|---|
| **전략** | Blind 탐색 → 나중에 비교 | PPI 결과로 포켓 좁히기 → Focused 도킹 |
| **Vina↔PPI 관계** | 독립 (병렬 가능) | 순차 의존 (PPI → Pocket → Vina) |
| **자동화** | `run_production.py` 완전 자동 | `run_production.py --lane adv-*` 완전 자동 |
| **PBS** | lane 스크립트 5개 + status/vina-gpu/phase3-gpu | `run_advanced_pipeline.pbs` + lane 스크립트 6개 |
| **사용 모듈** | Vina모듈 + PyRosetta Core + Verdict/Report/Validate | Phase 1~4 전체 + PyRosetta Core |

### 모듈-워크플로우 소속 매핑

| 보고서 섹션 | 모듈 | Workflow A | Workflow B |
|---|---|---|---|
| 섹션 3: Vina 모듈 | `egfr_pipeline/vina/` | Phase 1 (blind) + Phase 4 (postprocess) | Phase 3 (focused) |
| 섹션 4: PyRosetta Core | `egfr_pipeline/pyrosetta_docking/` | Phase 2 (엔진) | Phase 1 (엔진) |
| 섹션 5: Phase 1 (PPI-First) | `egfr_pipeline/phase1/` | Phase 2~3 (도킹+후처리) | Phase 1 전체 (TG 1.0~1.6) |
| 섹션 6: Phase 2 (Pocket) | `egfr_pipeline/phase2/` | — | Phase 2 (TG 2.0~2.7) |
| 섹션 7: Phase 3 (Focused) | `egfr_pipeline/phase3/` | — | Phase 3 (TG 3.0~3.6) |
| 섹션 8: Phase 4 (Scoring) | `egfr_pipeline/phase4/` | — | Phase 4 (TG 4.0~4.6) |
| 섹션 9: Cross-Cutting | `verdict.py`, `report.py`, `validate.py` | Phase 5~7 (finalize) | — |

---

## 3. Vina 모듈 `[A: Phase 1,4 | B: Phase 3]`

> 경로: `egfr_pipeline/vina/`
>
> Workflow A에서는 blind docking(Phase 1) + postprocess(Phase 4)로 사용. Workflow B에서는 focused docking(Phase 3)으로 사용.

### 3.1 Vina 도킹 실행 (`vina_executor.py`)

**역할**: AutoDock Vina를 이용한 blind/focused 소분자 도킹 실행

**주요 기능**:
- **Blind 모드**: 수용체 전체 표면을 탐색하는 글로벌 도킹. 수용체 PDB에서 원자 좌표의 min/max를 계산하여 자동으로 검색 박스를 설정
- **Focused 모드**: 사전 정의된 영역(예: C-lobe pocket) 또는 사용자 지정 좌표 중심으로 탐색
- **Region 프리셋**: `REGION_PRESETS` 딕셔너리에 사전 정의된 포켓 좌표 (예: clobe)
- **시드 재현성**: `derive_docking_seed(base_seed, receptor_id, ligand_id)` — receptor_id + ligand_id를 해시하여 결정론적 시드 생성
- **병렬 실행**: `ThreadPoolExecutor`를 이용한 다중 리간드 동시 도킹
- **SMILES 지원**: SMILES 문자열 입력 시 OpenBabel을 통해 3D 구조 생성 후 PDBQT로 변환
- **Config 저장/로드**: YAML 우선, JSON fallback. 재실행을 위해 모든 파라미터를 기록

**입력**: 수용체 PDB/PDBQT, 리간드 PDBQT/SDF/SMILES
**출력**: `{receptor_id}/{ligand_name}_{mode}.pdbqt` (포즈 파일)

---

### 3.2 포즈 파싱 (`parse_poses.py`)

**역할**: Vina 출력 PDBQT 파일에서 개별 포즈를 추출하고 구조화된 테이블로 변환

**알고리즘**:
1. PDBQT 파일을 `MODEL`/`ENDMDL` 블록 단위로 분리
2. 각 포즈에서:
   - `REMARK VINA RESULT:` 줄에서 affinity, rmsd_lb, rmsd_ub 파싱
   - 비수소 원자(heavy atom)의 좌표에서 centroid(무게중심) 계산
3. 메타데이터(도킹 모드, exhaustiveness, 시드 등)와 결합하여 행 생성

**출력**:
- `vina_pose_table.csv` — 전체 포즈 테이블 (21개 컬럼: receptor_id, ligand_id, pose_rank, affinity, centroid_xyz, ...)
- `vina_postprocess_coverage.csv` — 요청 대비 실제 파싱된 포즈 수 커버리지

---

### 3.3 접촉 잔기 추출 (`pose_contacts.py`)

**역할**: 각 도킹 포즈와 수용체 간 접촉 잔기를 거리 기반으로 계산

**알고리즘**:
1. 수용체 PDB에서 비수소 원자 좌표 파싱 (chain, residue_id, coord)
2. 포즈의 비수소 원자 좌표 파싱
3. 모든 수용체-리간드 원자 쌍 거리 계산
4. cutoff(기본 4.0 A) 이내이면 접촉 잔기로 분류
5. 각 접촉 잔기의 **최소 거리**를 기록

**출력**:
- `vina_pose_table.csv` 업데이트 (contact_residues, n_contact_residues, contact_distances 컬럼 추가)
- `vina_contact_distances.csv` — 잔기별 상세 거리 (long-form)

---

### 3.4 포켓 클러스터링 (`pocket_cluster.py`)

**역할**: 도킹 포즈들을 3D 공간 위치 기반으로 포켓 단위로 그룹화

**알고리즘** (반복 centroid 기반 클러스터링 — k-means 변형):
1. **시드 단계**: 포즈를 affinity 순으로 스캔. 기존 시드에서 cutoff(기본 8 A) 이내면 해당 시드에 할당, 아니면 새 시드 생성
2. **할당 단계**: 모든 포즈를 가장 가까운 시드에 할당
3. **업데이트 단계**: 시드 위치를 할당된 멤버의 평균 centroid로 재계산
4. **반복**: 2-3단계를 수렴 또는 max_iterations(기본 10)까지 반복
5. **소규모 흡수**: min_pocket_size(기본 2) 미만의 포켓은 가장 가까운 큰 포켓에 흡수

**후처리 잔기 기반 병합** (Union-Find):
- 동일 수용체 내 포켓 쌍의 접촉 잔기 Jaccard 유사도 또는 overlap coefficient가 임계값 이상이면 병합
- Centroid 거리가 6 A 이내면 잔기 데이터 없이도 병합 (centroid fallback)
- 전이적 폐쇄(Union-Find)로 A-B, B-C 병합 시 A-C도 병합

**포켓 캡** (선택 사항):
- 포켓당 최대 포즈 수 제한 (`max_per_pocket`)
- 리간드 간 라운드-로빈 선택으로 다양성 보장

**출력**:
- `vina_pose_table.csv` 업데이트 (pocket_id 컬럼)
- `vina_clustering_merge_log.csv` — 병합 결정 로그
- `vina_clustering_parameters.json` — 재현성을 위한 파라미터 기록
- `vina_pocket_cap_report.csv` — 캡 적용 시 통계

---

### 3.5 포켓 요약 (`pocket_summary.py`)

**역할**: 포켓 단위 통계 집계 및 리간드-포켓 매핑

**산출물** 3개 테이블:

1. **vina_pocket_table.csv** — 포켓별 요약 (16개 컬럼):
   - centroid 좌표, 포즈/리간드 수
   - best/mean affinity
   - 접촉 잔기 합집합 + 상위 5개 잔기
   - centroid_spread_A (포즈 centroid의 RMSD — 포켓 내 수렴도)
   - affinity_std, affinity_iqr
   - dominant_ligand_fraction, ligand_pose_entropy (리간드 다양성 지표)

2. **vina_drug_pocket_map.csv** — 리간드별 포켓 할당 (10개 컬럼):
   - 지배적 포켓 ID, 포즈 비율, best affinity, best pose rank
   - 대안 포켓 목록, multimodal binding 여부

3. **vina_pocket_residue_occupancy.csv** — 잔기별 점유율 (6개 컬럼):
   - 잔기가 포켓 내 포즈에서 관찰된 빈도
   - occupancy >= 0.5이면 is_hotspot = True

---

### 3.6 교차 수용체 비교 (`cross_receptor.py`)

**역할**: 서로 다른 수용체 상태의 포켓을 쌍별 비교

**비교 메트릭**:
- **Centroid 거리**: 3D 유클리드 거리
- **잔기 Jaccard 유사도**: 접촉 잔기 집합 간 Jaccard 지수
- **잔기 Overlap coefficient**: |교집합| / min(|A|, |B|)
- **공유 리간드**: 두 포켓 모두에서 도킹된 리간드
- **Same-patch 후보 플래그**: centroid 8 A 이내 AND (Jaccard >= 0.3 OR overlap >= 0.5)
- **Bootstrap CI**: 부트스트랩 centroid 표준편차로부터 95% 신뢰구간 계산

**출력**: `vina_pocket_comparison.csv` (29개 컬럼)

---

### 3.7 부트스트랩 안정성 분석 (`pocket_stability.py`)

**역할**: 포켓 클러스터링의 통계적 안정성을 부트스트랩 리샘플링으로 평가

**알고리즘**:
1. 원본 데이터로 레퍼런스 포켓 생성
2. N회(기본 100) 리샘플링: 80% 포즈를 복원추출하여 재클러스터링
3. 각 리플리케이트 포켓을 레퍼런스에 centroid 근접도로 매칭
4. 포켓별 안정성 지표 계산:
   - **pocket_exists_frac**: 리플리케이트에서 포켓이 재현된 비율
   - **centroid_std_A**: centroid 위치의 표준편차
   - **affinity_mean/std/iqr**: 결합 에너지 분포 안정성
   - **n_pose_mean/std**: 포즈 수 안정성

**출력**: `vina_pocket_bootstrap.csv`

---

### 3.8 Phase 3 브릿지 (`phase3_bridge.py`)

**역할**: Vina blind docking 결과를 Phase 3 focused docking의 입력 형식으로 변환

**분류 기준** (best affinity 기반):
- <= -7.0 kcal/mol -> **primary** (필수 도킹)
- <= -5.0 -> **secondary** (권장 도킹)
- <= -3.0 -> **exploratory** (탐색적)
- \> -3.0 -> **skip**

**박스 크기 자동 결정**: centroid_spread * 2.5 + 10 A padding, 18-30 A 범위로 클램핑

**출력**: `phase2_pockets/phase3_candidate_pocket_reference.csv` — Phase 3 TG 3.0 입력

---

## 4. PyRosetta Docking Core `[A: Phase 2 엔진 | B: Phase 1 엔진]`

> 경로: `egfr_pipeline/pyrosetta_docking/`

**역할**: PyRosetta 기반 단백질-단백질 글로벌 블라인드 도킹의 핵심 엔진

### 4.1 pipeline_manager.py (~1660줄)

7단계 파이프라인을 순차 실행하는 메인 오케스트레이터:

1. **Relax**: FastRelax (ref2015 score function)로 입력 구조를 에너지 최소화. 결과는 `relaxed_cache/`에 캐싱하여 재실행 시 건너뜀
2. **Global Docking**: `RigidBodyPerturbMover(360deg, 100A)` -> `DockingSlideIntoContact` -> `DockMCMProtocol`로 무작위 위치/방향에서 도킹 탐색
3. **Fast Scoring & Filtering**: v2.0 2-Pass 필터링 시스템
   - Pass 1: dG, dSASA, sc, total_score (전체 decoy 대상, 비용 ~0)
   - Stage 1: 에너지 기반 조기 탈락 (dG > 0 제거, total_score 백분위)
   - Mini Refinement: Stage 1 생존자에 인터페이스 사이드체인 리패킹
   - Stage 2 Cheap: dSASA >= 500, dG_density <= -1.0, sc >= 0.50
   - Stage 2 Expensive: packstat >= 0.55, delta_unsatHbonds <= 8, nres_int >= 10
   - Graduated Fallback (Level 0-3)
4. **L_RMSD Greedy Clustering**: CoM(Center of Mass) pre-filter + closest-match 클러스터링, 멤버 다양성 보장
5. **Refinement**: `DockMCMProtocol` 재적용으로 미세 조정
6. **Final Scoring & Selection**: Round-robin 다양성 선택 + L_RMSD 중복 제거
7. **Visualization & Report**: PyMOL 스크립트, Energy Funnel Plot, Validation Report (10개 품질 체크)

### 4.2 movers.py

Relax, Global Docking, Refinement의 구체적 PyRosetta Mover 구성

### 4.3 scoring.py

Scoring, RMSD, InterfaceAnalyzer, Per-residue 에너지 분석 워커들. PyRosetta 구버전 호환성을 위한 try-except/hasattr 가드 포함

### 4.4 pyrosetta_init.py

PyRosetta 초기화 유틸리티: PID 기반 랜덤 시드, `-mute all`, stdout/stderr 리다이렉트, Pose<->String 직렬화

### 4.5 logging_config.py

`pipeline`/`pipeline.worker` 계층적 로거 구성

### 4.6 run_metadata.py

실행 메타데이터 JSON 생성: config 경로, 모델 수, 필터 버전, 완료 상태 등

---

## 5. Phase 1 — PPI-First Interface Mapping `[A: Phase 2~3 | B: Phase 1]`

> 경로: `egfr_pipeline/phase1/`
>
> Workflow A에서는 PPI 도킹(Phase 2) + 후처리(Phase 3)로 사용. Workflow B에서는 전체 TG 1.0~1.6 순차 실행.
> **목적**: 3가지 EGFR 구조 상태에서 MYO1D beta-meander와의 PPI 인터페이스를 매핑하여, 약물이 결합하면 PPI를 교란할 수 있는 수용체 표면 패치를 정의

### TG 1.0: 입력 준비 (`prepare_inputs.py`)

- EGFR monomer 수용체 PDB + MYO1D beta-meander partner PDB 준비
- 3개 수용체 상태: 3GT8_raw (결정 구조), EGFR_160-185 (MD 클러스터), EGFR_170-200 (MD 클러스터)
- 메타데이터 CSV 생성: receptor_metadata.csv, partner_metadata.csv, docking_pair_metadata.csv
- 입력 검증 리포트 생성

### TG 1.1: Config 생성 + PyRosetta 도킹 실행 (`generate_configs.py`, `launch_docking.py`)

`generate_configs.py`는 3 states × 5 seeds에 대한 `.ini` 설정 파일을 자동 생성하고, `launch_docking.py`가 해당 config로 실제 도킹을 실행한다.

**역할**: PyRosetta PPI 도킹을 표준화된 방식으로 실행

**주요 작업**:
1. Phase 1 입력 검증 (TG 1.0 완료 확인)
2. `pyrosetta_run_metadata.json` 생성:
   - receptor_id, partner_id (extended_beta_meander_955_1006)
   - construct_type (full_kinase_domain — 이전 pilot의 C-lobe fragment와 차별)
   - filter_version (v1.0/v2.0), mini_refinement 여부
   - excluded_residues_A, enable_early_rejection
3. 표준화된 출력 경로: `output/phase1_ppi/{state}/{prod|test}_seed{n}/`
4. Multi-seed 프로덕션 실행: 3 states x 5 seeds = 15회 독립 실행 (300K+ models)
5. 스크래치 디렉토리 지원: HPC 로컬 디스크 사용 후 결과 동기화
6. `PipelineManager`를 내부적으로 호출하여 7단계 파이프라인 실행

**실행 모드**:
- `--dry-run`: 설정 검증 + 메타데이터 생성만 (도킹 없음)
- `--test`: 1K models/state
- `--production`: 20K models/state/seed x 5 seeds

### TG 1.1.5: 스코어 표준화 (`standardize_scores.py`)

**역할**: 여러 시드/실행에서 나온 final_ranking.csv를 하나의 표준 형식으로 통합

**출력**: `pyrosetta_decoy_scores.csv` (25개 컬럼) — decoy_id, receptor_id, partner_id, construct_type, seed_index, 전체 스코어링 메트릭

### TG 1.2: 인터페이스 잔기 추출 (`extract_interface.py`)

**역할**: PyRosetta 도킹 결과에서 수용체-파트너 인터페이스 잔기를 구조화된 테이블로 추출

**핵심 처리**:
1. `final_ranking.csv`의 `Binding_Residues_A` / `Binding_Residues_B` 컬럼 파싱
2. 잔기 번호 기반 N-lobe / C-lobe 분류 (경계: 잔기 838)
3. `InterfaceEnergies.csv`에서 잔기별 에너지(delta_e: fa_atr, fa_rep, fa_sol, fa_elec) 로드
4. 잔기 이름 정규화 (HSD->HIS, CYX->CYS 등)

**출력**:
- `pyrosetta_interface_residue_table.csv` — long-form: 모델 x 잔기 (17개 컬럼)
- `pyrosetta_interface_models.csv` — 모델별 요약 (18개 컬럼, N-lobe/C-lobe 잔기 수 포함)

### TG 1.2A: Orientation 필터 (`orientation_filter.py`)

**역할**: MYO1D beta-meander의 active face (sheets 8+9)가 수용체 방향을 향하는지 검증

**배경**: 리지드바디 도킹은 파트너의 방향을 제어하지 않으므로, beta-meander가 뒤집힌(flipped) 아티팩트 포즈가 발생할 수 있음. 이 필터가 생물학적으로 올바른 배향만 선별.

**Dual-Vector Orientation Test 알고리즘**:
1. Active face (sheets 8+9, 잔기 961-972)의 Ca 좌표 수집
2. PCA로 sheet-plane normal 계산 (최소 분산 방향)
3. Multi-probe Ca->Cb consensus로 normal을 active face 방향으로 정렬
   - 3개 probe 잔기 (VAL962, VAL964, SER971)의 Ca->Cb 벡터와 투표
4. Sheet centroid에서 local receptor centroid 방향 벡터 계산
5. Dot product -> orientation 분류:
   - \> +0.15: **PASS** (active face -> receptor)
   - < -0.15: **FAIL** (active face -> away)
   - |dot| < 0.15: **AMBIGUOUS** (edge-on)

**출력**: `orientation_filter_log.csv` + `pyrosetta_interface_models.csv`에 orientation_score/class 병합

### TG 1.3: 클러스터 합의 (`cluster_consensus.py`)

**역할**: 여러 도킹 모델에서 반복적으로 나타나는 인터페이스 패치(hotspot 잔기)를 식별

**핵심 분석**:
1. 모델을 클러스터별로 그룹화
2. Orientation-validated 모델만 합의에 참여 (가능한 경우)
3. 클러스터 내 잔기 점유율(occupancy) 계산: 해당 잔기가 orientation-valid 모델 중 몇 %에서 관찰되는가
4. Hotspot 판정: occupancy >= 50% (기본값)
5. 수용체 PDB에서 CA 좌표를 로드하여 클러스터 centroid와 centroid_spread 계산

**출력 3개 테이블**:
- `ppi_cluster_summary.csv` — 클러스터별 요약 (30개 컬럼): 멤버 수, 평균 메트릭, N-lobe/C-lobe 분포, hotspot 잔기
- `ppi_hotspot_residues.csv` — 잔기별 점유율 (14개 컬럼)
- `ppi_interface_patch_table.csv` — 잔기별 전체 통계 (17개 컬럼): 몇 개 클러스터에 출현, 글로벌 점유율, CA 좌표

### TG 1.4: LightDock 2차 검증 (`lightdock_validation.py`)

**역할**: PyRosetta와 독립적인 두 번째 도킹 방법(LightDock)으로 인터페이스 결과를 교차 검증

**서브모듈**:
- **1.4.1 Setup**: LightDock 실행 스크립트 생성 (서버 측)
- **1.4.2 Extract**: Top-ranked swarm 포즈에서 인터페이스 잔기 추출 (CA-CA < 10 A) + orientation filter 적용
- **1.4.3 Convergence**: PyRosetta vs LightDock 잔기 교집합/차집합 분석

**핵심 개념**: LightDock은 swarm 기반 도킹(DFIRE2 scoring). PyRosetta와 완전히 다른 알고리즘이므로 양쪽에서 모두 나타나는 잔기는 방법 독립적 증거가 됨.

### TG 1.5: Multi-State 인터페이스 비교 (`compare_states.py`)

**역할**: 3개 수용체 상태에 걸쳐 인터페이스 잔기의 로버스트니스(강건성) 분류

**로버스트니스 분류**:
- **Robust**: 3개 상태 모두에서 관찰 -> MYO1D 결합 패치의 최강 후보
- **Moderate**: 2개 상태에서 관찰
- **State-specific**: 1개 상태에서만 관찰 -> conformationally gated 상호작용 가능성

**출력**:
- `ppi_patch_cross_state_comparison.csv` — 잔기별 다중 상태 비교
- `ppi_patch_state_robustness.csv` — 로버스트니스 분류
- `phase1_interface_comparison_report.md` — Markdown 리포트

### TG 1.6: Phase 1 리뷰 리포트 (`review_report.py`)

**역할**: Phase 1 전체를 요약하고 Phase 2로의 핸드오프 파일을 생성

**출력**:
- `phase1_interface_report.md` — 8섹션 종합 리포트:
  1. 입력 요약 (레거시 vs Phase 1 차이)
  2. 상태별 PyRosetta 증거
  3. Orientation 필터링 요약
  4. Hotspot 잔기 상세
  5. Cross-state 로버스트니스
  6. LightDock 수렴 분석
  7. 증거 계층 (Primary: PyRosetta, Secondary: LightDock, Auxiliary: AlphaFold-Multimer)
  8. Phase 2 핸드오프 내용

- `phase1_downstream_patch_reference.csv` — **Phase 2 핸드오프 파일** (17개 컬럼):
  - 각 수용체측 잔기의 evidence_source, robustness_class, method_agreement, confidence(high/medium/low)

### TG 1.7: 파일럿 데이터 비교 (`pilot_comparison.py`)

**역할**: 이전 파일럿 실행(C-lobe fragment 기반)과 현재 Phase 1(full kinase domain) 결과를 비교

- 파일럿 데이터가 존재하는 경우에만 실행 (선택적)
- 잔기 수준 겹침, 에너지 분포 비교, construct_type 차이 분석

---

## 6. Phase 2 — Pocket Analysis & Druggability Assessment `[B only]`

> 경로: `egfr_pipeline/phase2/`
>
> **Workflow B 전용**. `run_production.py --lane adv-phase2` 또는 `qsub config/run_adv_phase2.pbs`로 실행.
>
> **목적**: Phase 1에서 정의한 PPI 패치와 독립적인 포켓 탐지 도구(fpocket, P2Rank)의 결과를 통합하여, 약물 결합 가능성이 있는 포켓을 식별하고 PPI 패치와의 관계를 파악
>
> **선행 조건**: `phase1_downstream_patch_reference.csv` (Phase 1 TG 1.6 출력)

### TG 2.0: Patch Reference 수신 (`patch_ingestion.py`)

- Phase 1의 `phase1_downstream_patch_reference.csv` 로드
- 스키마 검증 (17개 필수 컬럼 확인)
- Phase 2 내부 형식으로 정규화
- 검증 리포트 생성

### TG 2.1: Pocket 제안 (`pocket_proposal.py`)

**역할**: fpocket/P2Rank을 이용한 리간드 결합 가능 포켓 탐지

**작업**:
1. **Setup 생성**: 서버 측 실행용 셸 스크립트 + 메타데이터 JSON 생성
2. **Output 파싱**: fpocket (알파 구 기반) 및 P2Rank (기계학습 기반) 출력을 도구 무관(tool-agnostic) 스키마로 변환
3. **포켓 메타데이터 추출**: centroid, 점수, 부피, 잔기 목록, 도킹 박스 크기

**출력**: `candidate_pockets_raw.csv` — 도구 무관 포켓 후보 목록

### TG 2.2: Pocket 병합 (`pocket_merge.py`)

**역할**: 서로 다른 도구에서 제안된 중복 포켓을 병합

**병합 기준**:
- Centroid 거리 <= 6.0 A
- 잔기 Jaccard >= 0.30

**출력**:
- `candidate_pockets.csv` — 병합된 최종 포켓 목록
- `candidate_pocket_merge_table.csv` — 병합 결정 로그
- `candidate_pocket_provenance.csv` — raw -> merged 매핑 (출처 추적)

### TG 2.3: Patch-Pocket 관계 분석 (`patch_relationship.py`)

**역할**: Phase 1 PPI 패치와 Phase 2 포켓 후보 간의 공간적/잔기적 관계 분류

**관계 클래스**:
- **orthosteric_candidate**: 포켓이 PPI 패치와 직접 겹침 -> 리간드가 PPI를 직접 차단
- **rim_candidate**: 포켓이 PPI 인터페이스 가장자리 -> 간접적 교란 가능
- **allosteric_candidate**: 공간적으로 떨어져 있지만 잠재적 알로스테릭 효과
- **no_relationship**: PPI와 관련 없는 포켓

### TG 2.4: Druggability 신뢰도 (`druggability_confidence.py`)

**역할**: 포켓의 약물 결합 가능성을 정량화

**작업**:
1. **점수 정규화**: 도구별 raw score 보존 + 정규화 점수 추가
2. **Multi-source 지원**: fpocket + P2Rank 양쪽에서 모두 탐지된 포켓은 consensus 표시
3. **Druggability 등급**: high (>= 0.50) / medium (>= 0.25) / low
4. **FTMap 확장 스키마**: 핫스팟 분석 데이터를 추가할 수 있는 확장 가능 구조

### TG 2.5: Cross-State 포켓 정렬 (`cross_state_alignment.py`)

**역할**: 3개 수용체 상태에 걸쳐 포켓이 보존/이동/상태특이적인지 분류

**정렬 기준**:
- Centroid <= 8 A + 잔기 Jaccard >= 0.20 -> **state_robust_pocket**
- Centroid 8-15 A -> **state_shifted_pocket**
- 1개 상태에서만 발견 -> **state_specific_pocket**
- 모호한 경우 -> **uncertain_alignment**

### TG 2.6: Phase 3 Export (`phase3_export.py`)

**역할**: Phase 3 도킹에 필요한 깨끗한 포켓 카탈로그 내보내기

**출력**: `phase3_candidate_pocket_reference.csv` — Phase 3이 소비하는 단일 파일:
- pocket_id, centroid_xyz, box_size_xyz
- relationship_class, druggability tier
- docking_priority (primary/secondary/exploratory/skip)
- state_class, n_residues, residue_ids

### TG 2.7: Phase 2 리뷰 리포트 (`review_report.py`)

Phase 2 전체 결과를 요약한 Markdown 리포트 생성

### Cascade Runner (`rerun_cascade.py`)

**역할**: TG 2.0 -> 2.7을 순차 자동 실행

- `--parse-only`: fpocket/P2Rank 서버 실행 완료 후 분석만 수행
- `--from-tg 2.5`: 특정 TG부터 재실행
- `--ftmap-dir`: FTMap 데이터 경로 지정

---

## 7. Phase 3 — Diversity-Aware Focused Docking `[B only]`

> 경로: `egfr_pipeline/phase3/`
>
> **Workflow B 전용**. `run_production.py --lane adv-phase3-*` 또는 `qsub config/run_adv_phase3_*.pbs`로 실행.
> Cascade runner: `egfr_pipeline/phase3/rerun_cascade.py` (setup / execute / post 3모드)
>
> **목적**: Phase 2에서 식별된 포켓에 대해 리간드별 focused Vina 도킹을 수행하되, 검색 노력을 포켓 간에 균등 분배하여 지배적 포켓에 편중되지 않도록 함
>
> **선행 조건**: `phase3_candidate_pocket_reference.csv` (Phase 2 TG 2.6 출력)

### TG 3.0: Phase 2 Reference 수신 (`pocket_reference_ingestion.py`)

- `phase3_candidate_pocket_reference.csv` 로드 및 검증
- 필수 필드 확인: pocket_id, centroid, docking_priority 등
- 내부 정규화 버전 생성

### TG 3.1: 도킹 잡 구성 (`job_construction.py`)

**역할**: 포켓-리간드-시드 조합별 Vina 도킹 잡 매트릭스 생성

**출력**:
- `phase3_docking_job_table.csv` — 전체 잡 목록 (receptor x pocket x ligand x seed)
- `phase3_job_box_table.csv` — 포켓별 도킹 박스 메타데이터

### TG 3.2: Search Budget Policy (`budget_policy.py`)

**역할**: 포켓별 검색 예산을 정의하여 이미 충분히 탐색된 포켓이 추가 자원을 소모하지 않도록 함

**예산 파라미터**:
- `max_rounds`: 포켓당 최대 도킹 라운드 (기본 3)
- `saturation_pose_threshold`: 허용 가능한 포즈가 이 수에 도달하면 포화 (기본 10)
- `saturation_affinity_window_kcal`: best affinity에서 이 범위 내 포즈가 "허용 가능" (기본 2.0 kcal/mol)

**포켓 상태**:
- `open` -> `saturated` (충분한 포즈 확보)
- `open` -> `exhausted` (라운드 소진)
- `skipped` (priority=skip)

### TG 3.3: 도킹 실행 (`run_diverse_docking.py`)

**역할**: Budget-aware, pocket-guided Vina 실행

**특징**:
- Round 기반 실행: 각 라운드에서 open 상태의 포켓만 도킹
- 포화 감지 시 자동 중단 -> 예산을 다른 포켓으로 재할당
- Workspace(개발) vs Server(프로덕션) 규칙 자동 구분
- Setup 모드: 서버 실행용 스크립트만 생성
- Dry-run 모드: 검증만 수행

### TG 3.4: Pose Attribution (`pose_attribution.py`)

**역할**: Phase 3 도킹 결과에 포켓 출처 정보를 추가

**확장 컬럼**: candidate_pocket_id, round_id, seed, docking_mode, job_id — 기존 포즈 테이블과 하위 호환

### TG 3.5: Diversity 검증 (`diversity_validation.py`)

**역할**: 검색 노력이 포켓 간에 실제로 균등 분배되었는지 검증

**메트릭**:
- 포켓 점유율 요약 (포즈 수, 시드 수, best/mean affinity)
- 다양성 메트릭 (concentration ratio, 비탐색 포켓 수)
- Blind vs Diverse 기준선 비교 (선택 사항)

### TG 3.6: Phase 4 Export (`phase4_export.py`)

**역할**: Phase 4에 필요한 도킹 증거 패키지 내보내기

**출력**: `phase4_docking_evidence_reference.csv`:
- 포켓별 리간드 지원 강도 (pose_support_count, best/mean affinity)
- 포켓 컨텍스트 (relationship_class, druggability_tier)
- 다양성 메타데이터 (검색 상태, 검색 완전성)

---

## 8. Phase 4 — Perturbation Relevance Scoring `[B only]`

> 경로: `egfr_pipeline/phase4/`
>
> **Workflow B 전용**. `run_production.py --lane adv-phase4` 또는 `qsub config/run_adv_phase4.pbs`로 실행.
> Cascade runner: `egfr_pipeline/phase4/rerun_cascade.py`
>
> **목적**: Phase 1-3의 모든 증거를 통합하여, 각 포켓-리간드 조합이 EGFR-MYO1D PPI를 실제로 교란할 가능성을 다축(multi-axis) 스코어링으로 평가
>
> **선행 조건**: `phase4_docking_evidence_reference.csv` (Phase 3 TG 3.6 출력)

### TG 4.0: 증거 수집 (`evidence_ingestion.py`)

- Phase 1 패치 레퍼런스 + hotspot 요약 로드
- Phase 2 포켓 관계 + druggability + state class 로드
- Phase 3 도킹 증거 + 리간드 지원 + diversity 데이터 로드
- Cross-phase 일관성 검증

### TG 4.1: Score Framework 정의 (`score_framework.py`)

**4개 핵심 스코어링 축**:

| 축 | 이름 | 가중치 | 설명 |
|---|------|--------|------|
| A1 | PPI Interface Confidence | 0.30 | PPI 패치 증거의 강도/신뢰도 (hotspot overlap, multi-method agreement, state robustness) |
| A2 | Druggability | 0.25 | 포켓의 약물 결합 가능성 (fpocket/P2Rank 점수, multi-source consensus, 부피) |
| A3 | Perturbation Relevance | 0.30 | 리간드가 PPI를 실제로 교란할 가능성 (Vina affinity, pocket-PPI 공간 근접도, 접촉 잔기 overlap) |
| A4 | State Robustness | 0.15 | 구조 상태에 걸친 포켓 안정성 (robust/shifted/state-specific) |

### TG 4.2: 기계적 분류 (`mechanistic_classification.py`)

**기계적 클래스**:
- **orthosteric_disruptor_candidate**: PPI 패치와 직접 겹침. 리간드가 PPI를 입체적으로 차단
- **interface_rim_modulator_candidate**: PPI 인터페이스 가장자리. 주변 접촉을 간접 약화
- **allosteric_modulator_candidate**: PPI 패치에서 떨어져 있지만 알로스테릭 효과 가능
- **insufficient_evidence**: 증거 부족으로 분류 불가

**분류 신뢰도**: high / medium / low / speculative — 증거 강도에 따라

### TG 4.3: Perturbation Scoring (`perturbation_scoring.py`)

**역할**: 4축 점수를 가중 합산하여 최종 순위 생성

**Affinity Cap**: affinity 유래 신호의 최대 기여를 35%로 제한 -> 무관한 사이트의 높은 affinity가 순위를 지배하는 것 방지

### TG 4.4: State Interpretation (`state_interpretation.py`)

**역할**: Cross-state 증거를 해석하여 포켓 접근성 분류

**해석 라벨**:
- `persistent`: 모든 상태에서 존재 -> 높은 접근성 신뢰도
- `conditionally_accessible`: 상태 간 이동 -> 상태 인식 전략 필요
- `state_dependent`: 1개 상태에서만 -> cryptic site 또는 불완전 샘플링 가능

### TG 4.5: Review-First 출력 (`review_report.py`)

**출력**:
- `phase4_final_review_table.csv` — 포켓 x 리간드 1행 축약 리뷰 (20개 컬럼)
- `phase4_expanded_evidence_table.csv` — 전체 provenance (upstream 추적 가능)

### TG 4.6: 최종 보고서 (`final_report.py`)

**역할**: 기계적 클래스별 최종 순위 후보를 제시하는 통합 Markdown 리포트

**내용**:
1. Executive Summary (포켓/리간드/수용체 수, 클래스별 분포)
2. 기계적 클래스별 상세 순위
3. 왜 affinity만으로는 최종 기준이 될 수 없는지 설명
4. 검증/주의 사항

### Presentation Summary (`presentation_summary.py`)

슬라이드용 축약 요약 생성

---

## 9. Cross-Cutting 유틸리티 `[A: Phase 5~7]`

### 9.1 Site Verdict (`verdict.py`)

**역할**: Vina 포켓을 3개 독립 축으로 증거 강도 분류

**3개 축**:
- **Axis 1 — Vina Quality**: 결합 에너지, 포즈 수렴, 다중 리간드 합의
- **Axis 2 — PPI Spatial Proximity**: Vina 포켓 centroid와 PPI 인터페이스 centroid 간 3D 거리
- **Axis 3 — Cross-Receptor Consistency**: 동일 포켓이 여러 수용체 상태에서 발견되는지

**분류**: STRONG / MODERATE / WEAK (증거 강도 요약, 타당성 판단이 아닌 조사 우선순위 가이드)

**출력**:
- `cross_method_agreement.csv` — Vina <-> PPI 공간/잔기 분석
- `valid_sites.csv` — 증거 분류

### 9.2 Report (`report.py`)

프로젝트 수준 종합 리포트 생성: Vina 포켓 요약, cross-receptor 비교 하이라이트, 리간드-포켓 매핑, PPI 잔기 증거 (보조)

### 9.3 Validate (`validate.py`)

파이프라인 출력 검증:
1. 핵심 출력 파일 존재 확인
2. Receptor/ligand ID 일관성 검증
3. CSV 필드 구조 스키마 회귀 검증
4. 잔기 번호 일관성
5. 핸드오프 문서 존재 확인

---

## 10. Step View / 출력 조직화

### 10.1 Step View (`step_view.py`)

**역할**: 프로덕션 출력을 단계별로 조직화된 뷰로 제공

**7개 Step**:

| Step | 이름 | 내용 |
|------|------|------|
| 1 | vina_raw | Raw blind docking 포즈 인벤토리 |
| 2 | ppi_docking | PPI 도킹 결과 |
| 3 | ppi_evidence | PPI 인터페이스 증거 (잔기 테이블, 요약) |
| 4 | vina_analysis | Vina 포켓 분석 (포켓 테이블, 포즈 테이블, 약물-포켓 맵) |
| 5 | site_verdict | 사이트 판정 (유효 사이트, cross-method agreement) |
| 6 | report | 종합 보고서 + 잔기 증거 |
| 7 | validation | 검증 상태 (JSON + 텍스트 요약) |

### 10.2 Output Organizer (`output_organizer.py`)

`steps/` 디렉토리 아래에 심볼릭 링크 기반 조직화된 뷰를 생성. `STEP_INDEX.txt` 인덱스 파일 포함.

---

## 11. 전체 데이터 흐름

### Workflow A: Standard Production

> `run_production.py` 자동화. Vina와 PPI는 **독립 증거원** — Phase 5에서 처음 병합.

```
입력: EGFR PDB (3 states) + MYO1D beta-meander + 리간드 라이브러리
         |                              |
         v                              v
  [Phase 1] Vina Blind           [Phase 2] PPI Global Blind
  Docking (소분자)                Docking (PyRosetta, 3×5 seeds)
         |                              |
         v                              v
  [Phase 4] Vina Postprocess     [Phase 3] PPI Postprocess
  (파싱→클러스터→포켓 요약)       (인터페이스 잔기 추출)
         |                              |
         +-------------+----------------+
                        |
                        v
              [Phase 5] Site Verdict
              (Vina 포켓 × PPI 증거 3축 통합)
                        |
                        v
              [Phase 6] Report → [Phase 7] Validate
```

**서버 제출**: `qsub config/run_production.pbs` (순차) 또는 lane별 PBS 병렬 제출 (섹션 2 참고)

### Workflow B: Advanced PPI-First Pipeline

> 순차 의존. 각 Phase의 출력 CSV가 다음 Phase의 **필수 입력**.

```
[Phase 2 PPI 도킹 완료 (Workflow A에서)]
         |
         v
  [Phase 1 분석] TG 1.0~1.6
  인터페이스 매핑 + 패치 정의
         |
         v  phase1_downstream_patch_reference.csv
  [Phase 2 Pocket Analysis] TG 2.0~2.7
  fpocket/P2Rank + PPI 패치 관계 분석
         |
         v  phase3_candidate_pocket_reference.csv
  [Phase 3 Focused Docking] TG 3.0~3.6
  Budget-aware Vina (포켓 타겟)
         |
         v  phase4_docking_evidence_reference.csv
  [Phase 4 Perturbation Scoring] TG 4.0~4.6
  4축 통합 스코어링 → 최종 순위
```

**서버 제출**: `qsub config/run_advanced_pipeline.pbs` — 내부에서 PBS 의존성 체인으로 자동 제출.

**안전성 가드**: 각 `_adv_*` lane 함수는 `_validate_adv_handoff()`를 통해 이전 Phase의 핸드오프 파일 존재를 사전 검증한다 (defense-in-depth). Phase 3 cascade는 모드별 사전 조건(setup → job table, post → round log)을 추가 검증하며, Vina 가용성 가드가 silent all-skip을 방지한다.

### 핵심 핸드오프 파일 (Workflow B)

| 구간 | 파일 | 생성 TG |
|------|------|---------|
| Phase 1 → Phase 2 | `phase1_downstream_patch_reference.csv` | TG 1.6 |
| Phase 2 → Phase 3 | `phase3_candidate_pocket_reference.csv` | TG 2.6 |
| Phase 3 → Phase 4 | `phase4_docking_evidence_reference.csv` | TG 3.6 |
| Vina blind → Phase 3 (보조) | `phase3_candidate_pocket_reference.csv` | via `phase3_bridge.py` |

---

*이 보고서는 2026-03-17 기준 소스코드를 분석하여 생성되었습니다. 워크플로우 구분 반영: 2026-03-17. 견고성 가드 추가: 2026-03-17.*
