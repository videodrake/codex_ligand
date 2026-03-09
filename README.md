# PyRosetta PPI Global Blind Docking Pipeline

PyRosetta 기반 **단백질-단백질 상호작용(PPI) Global Blind Docking** 자동화 파이프라인.
2-chain PDB 구조를 입력받아 수천~수만 개의 도킹 모델을 생성하고, 다중 기준 필터링 및 L_RMSD 클러스터링을 거쳐 **후보 결합 사이트를 탐색**한다. 단일 최적 포즈 예측이 아닌, 단백질 표면 전체에서 가능한 결합 포켓을 발견하는 것이 목적이다.

---

## 1. 프로젝트 상태

- **버전**: v2.0 (필터링 파이프라인 개편) / v1.0 하위 호환
- **타겟**: EGFR 도메인 - TH1/beta 상호작용 연구
- **환경**: Linux HPC Cluster (PBS/qsub), 32 CPU cores, **네트워크 차단**
- **v2.0 상태**: 구현 완료, 1K 모델 테스트 진행 중

---

## 2. 파일 구조

### 핵심 코드 (4개 파일)

| 파일 | 역할 | 크기 |
|------|------|------|
| `pipeline_manager.py` | 메인 오케스트레이터. 7단계 파이프라인 제어, v2.0/v1.0 자동 분기 | ~1660줄 |
| `docking.py` | 워커: FastRelax, Global Docking, Refinement 수행 | ~250줄 |
| `analysis.py` | 워커: Fast/Intermediate/Full 스코어링, RMSD, 잔기 분석 | ~590줄 |
| `common.py` | 유틸리티: PyRosetta 초기화, Pose<->String 변환 | ~150줄 |

### 설정 파일

| 파일 | total_global_models | 용도 |
|------|---------------------|------|
| `config.ini` | 1,000 | 테스트/디버깅 (~3-4시간) |
| `config_20k.ini` | 20,000 | 논문용 (~24-36시간) |
| `config_100k.ini` | 100,000 | 대규모 프로덕션 (~60-120시간) |

### 기타

| 파일/디렉토리 | 설명 |
|---------------|------|
| `run_v1.pbs` | PBS 배치 스크립트 v1.0 (`input_PDB/*.pdb` 순차 처리) |
| `run_v2_test.pbs` | PBS 배치 스크립트 v2.0 (단일 PDB 테스트, 기본: C-lobe_beta) |
| `input_PDB/` | 입력 PDB 파일 (EGFR_TH1, EGFR_beta, C-lobe_TH1, C-lobe_beta) |
| `relaxed_cache/` | FastRelax 결과 캐시 (전역, PDB별 재사용) |
| `CLAUDE.md` | Claude Code 자동 컨텍스트 (세션 시작 시 자동 로드) |
| `MANUAL.md` | 상세 매뉴얼 (알고리즘, 설정, 모듈 레퍼런스) |
| `check_filter_v2_compat.py` | v2.0 필터 메트릭 호환성 점검 (25 항목) |
| `check_improvements_compat.py` | v1.0 개선사항 호환성 점검 |

---

## 3. 파이프라인 아키텍처

```
입력 PDB (2-chain: A_B)
        |
   [Step 1] FastRelax (ref2015) ──── relaxed_cache/ 캐싱
        |
   [Step 2] Global Blind Docking (N천~N만 모델)
        |     RigidBodyPerturb(360deg,100A)
        |     -> DockingSlideIntoContact
        |     -> DockMCMProtocol
        |     (+ Early Rejection: 금지구역 접촉 검사)
        |
   [Step 2.5] Scoring & Multi-Stage Filtering (v2.0 / v1.0 자동 분기)
        |
        |  ┌─ v2.0 ([FilterStage1] 섹션 존재 시) ──────────────┐
        |  │ Pass 1: Fast Scoring (dG, dSASA, sc, total_score)  │
        |  │ Stage 1: Coarse Energy Filter                      │
        |  │   - 제외구역 hard filter                            │
        |  │   - dG > 0 제거 (반발적 인터페이스)                  │
        |  │   - total_score 백분위 필터 (상위 N%)                │
        |  │ Mini Refinement: 인터페이스 리패킹 (다중 라운드)    │
        |  │   (IncludeCurrent + ExtraRotamers ex1/ex2)          │
        |  │ Stage 2 Cheap: dSASA, dG_density, sc_value         │
        |  │ Pass 2: Intermediate Scoring (Stage 1 생존자만)     │
        |  │ Stage 2 Expensive: packstat, unsatHb, nres, hbonds │
        |  │ Graduated Fallback v2.0 (Level 0~3)                │
        |  └────────────────────────────────────────────────────┘
        |  ┌─ v1.0 (레거시, [FilterStage1] 없을 때) ───────────┐
        |  │ dSASA → 제외구역 → dG 백분위 → sc (기존 fallback)  │
        |  └────────────────────────────────────────────────────┘
        |
   [Step 3] L_RMSD Greedy Clustering
        |     CoM pre-filter -> Closest-match assignment
        |     Member diversity 보장 (pairwise L_RMSD >= 2.0A)
        |     최대 cluster_top_n개 클러스터, 클러스터당 members_per_cluster명
        |
   [Step 4] Targeted Refinement (DockMCMProtocol)
        |     각 대표 구조 x refine_per_struct회
        |
   [Step 5] Final Full Scoring
        |     dG, dSASA, sc, packstat, L_RMSD, binding residues, per-residue CSV
        |
   [Step 6] Diversity-Aware Selection
        |     Round-robin (클러스터 균등 선택)
        |     + L_RMSD 중복 제거
        |     + Key residue bonus (adjusted_dG)
        |
   [Step 7] Visualization & Validation Report
              PyMOL 스크립트, Energy Funnel Plot, 자동 품질 평가
```

---

## 4. 실행 방법

### 4.1 직접 실행

```bash
conda activate pyrosetta

# 특정 PDB 지정
python pipeline_manager.py config.ini input_PDB/EGFR_TH1.pdb

# config의 기본 PDB 사용
python pipeline_manager.py config.ini
```

### 4.2 PBS 배치 실행

```bash
# v1.0: input_PDB/ 내 모든 PDB 순차 처리
qsub run_v1.pbs
qsub -v CONFIG_FILE=config_100k.ini run_v1.pbs

# v2.0 테스트: 단일 PDB (기본: C-lobe_beta.pdb)
qsub run_v2_test.pbs
qsub -v CONFIG_FILE=config_100k.ini,INPUT_PDB=input_PDB/EGFR_TH1.pdb run_v2_test.pbs
```

### 4.3 로그 모니터링

```bash
tail -f pipeline.log             # v1.0 파이프라인 진행 상황
tail -f pipeline_v2_test.log     # v2.0 테스트 진행 상황
tail -f worker_debug.log         # 워커 프로세스 디버그
```

### 4.4 호환성 점검 (서버 배포 전)

```bash
python check_filter_v2_compat.py      # v2.0 메트릭 25항목 점검
python check_improvements_compat.py   # v1.0 개선사항 점검
```

---

## 5. 설정 가이드 (config.ini)

### [Docking]

| 파라미터 | 기본값 | 설명 |
|----------|--------|------|
| `total_global_models` | 1000 | 전역 도킹 샘플 수. 프로덕션: 10K~100K 권장 |

### [Filter]

| 파라미터 | 기본값 | 설명 |
|----------|--------|------|
| `filter_percentile` | 15 | dG 하위 N% 통과 (백분위 기반, 분포 독립적) |
| `min_dSASA` | 100.0 | 최소 인터페이스 표면적 (A^2). 비접촉 포즈 제거 |
| `contact_distance` | 6.0 | 결합 잔기 탐지 거리 (A). 6.0 = PPI 표준 |
| `min_sc_value` | 0.3 | 최소 Shape Complementarity |
| `min_survivors` | 50 | 최소 생존 구조 수 (fallback 트리거) |

### [Cluster]

| 파라미터 | 기본값 | 설명 |
|----------|--------|------|
| `threshold` | 4.0 | L_RMSD 클러스터링 임계값 (A) |
| `cluster_top_n` | 20 | 최대 클러스터 수 |
| `members_per_cluster` | 5 | 클러스터당 최대 멤버 |
| `member_diversity_dist` | 2.0 | 멤버 간 최소 L_RMSD (A) |

### [Refinement]

| 파라미터 | 기본값 | 설명 |
|----------|--------|------|
| `refine_per_struct` | 10 | 대표 구조당 정밀화 반복 |
| `perturb_trans` | 0.1 | 병진 섭동 (A) |
| `perturb_rot` | 1.0 | 회전 섭동 (deg) |

### [Constraints] (선택 사항)

| 파라미터 | 기본값 | 설명 |
|----------|--------|------|
| `excluded_residues_A` | (비어있으면 비활성) | Chain A 금지 잔기. 예: `700-807,840-854` |
| `key_residues_B` | (비어있으면 비활성) | Chain B 핵심 잔기 (soft bonus) |
| `enable_early_rejection` | false | Global Docking 중 조기 거부 활성화 |
| `exclusion_contact_dist` | 10.0 | 금지구역 접촉 판정 거리 (A) |
| `max_excluded_contacts` | 0 | 허용 금지잔기 접촉 수 (0=전면 거부) |
| `key_residue_bonus_weight` | 0.0 | 핵심잔기 접촉 보너스 가중치 |

### [MiniRefinement] (v2.0 선택 사항)

Stage 1과 Stage 2 사이에 인터페이스 사이드체인 리패킹을 수행하여 리지드바디 도킹 후 불량한 패킹을 개선합니다. IncludeCurrent(현재 로타머 보존) + ExtraRotamersGeneric(ex1+ex2 chi 확장) + 다중 라운드 리패킹으로 Stage 2 통과율을 높입니다.

| 파라미터 | 기본값 | 설명 |
|----------|--------|------|
| `enabled` | false | Mini Refinement 활성화 여부 |
| `mode` | pack_only | `pack_only`=리패킹만, `pack_min`=리패킹+최소화 |
| `n_rounds` | 3 | 리패킹 반복 횟수 (많을수록 더 나은 패킹, ~선형 비용) |

**성능 영향**: 100K 모델 → Stage 1 생존 ~1,500개 × 3라운드 ≈ 15분 (16코어). 전체 파이프라인 대비 무시 가능.

### [FilterStage1] (v2.0 신규 - 선택 사항)

이 섹션이 존재하면 v2.0 필터링 경로가 자동 활성화됨.

| 파라미터 | 기본값 | 설명 |
|----------|--------|------|
| `total_score_percentile` | 10 | total_score 상위 N% 유지 (Coarse filter) |
| `max_dG_separated` | 0.0 | dG > 이 값인 구조 즉시 제거 (반발적 인터페이스) |
| `enabled` | true | Stage 1 활성화 여부 |

### [FilterStage2] (v2.0 신규 - 선택 사항)

| 파라미터 | 기본값 | 설명 |
|----------|--------|------|
| `min_dSASA` | 500.0 | 최소 인터페이스 면적 (A^2). Mini Refinement 보완으로 완화 |
| `max_dG_density` | -1.0 | 최대 에너지 밀도 (dG/dSASA×100). 낮을수록 좋음 |
| `min_sc_value` | 0.50 | 최소 Shape Complementarity. Mini Refinement 후 기준 |
| `min_packstat` | 0.55 | 최소 패킹 품질. 0.0 = 비활성 |
| `max_delta_unsatHbonds` | 8 | 최대 미충족 수소결합 수. 99 = 비활성 |
| `min_nres_int` | 10 | 최소 인터페이스 잔기 수. 0 = 비활성 |
| `min_hbonds_int` | 0 | 최소 인터페이스 수소결합 수. 0 = 비활성 |
| `enable_expensive_metrics` | true | Pass 2 (packstat, unsatHb, nres, hbonds) 활성화 |
| `min_survivors` | 50 | Fallback 트리거 최소 생존자 수 |

### [Output]

| 파라미터 | 기본값 | 설명 |
|----------|--------|------|
| `save_filter_max` | 200 | 필터 단계 최대 저장 수 |
| `save_top_n` | 20 | 최종 결과 저장 수 |

---

## 6. 출력 결과

### 디렉토리 구조

```
<PDB_NAME>/
├── filter_passed/
│   ├── F0001_S-15.23.pdb             # 필터 통과 1위
│   └── ...
├── cluster_results/
│   ├── C01_M01_S-15.23.pdb           # 사이트 1 Leader
│   ├── C01_M02_S-14.10.pdb           # 사이트 1 Member 2
│   ├── cluster_summary.csv            # 사이트별 요약
│   └── ...
├── final_result/
│   ├── Rank01_C01_M01_S-18.45.pdb    # 최종 1위 (S뒤 숫자 = dG)
│   ├── Rank01_C01_M01_Energies.csv   # Per-residue 에너지
│   ├── view_results.pml               # 최종 랭킹 모델 (B-factor 컬러링)
│   ├── 2_DETAIL_C01.pml              # 사이트별 상세 뷰
│   └── docking_validation_report.txt  # PPI 사이트 탐색 리포트
├── final_ranking.csv                  # 종합 랭킹 (모든 메트릭 포함)
├── energy_funnel.png                  # L_RMSD vs dG 산점도
└── 1_OVERVIEW_Clusters.pml            # 전체 사이트 비교 (클러스터별 색상)
```

> **PML 파일 특징**: 상대경로 사용 (폴더째 로컬 복사 시 바로 작동), chain A 기준 자동 정렬, OVERVIEW는 클러스터별 색상 구분 (chain A=회색, chain B=클러스터별 색상) + PyMOL 그룹 지원.

> **실행 시 자동 정리**: `filter_passed/`, `cluster_results/`, `final_result/`는 파이프라인 시작 시 자동 삭제 후 재생성됨. `relaxed_cache/`는 보존.

### 스코어링 메트릭

| Metric | Excellent | Good | Marginal | Poor |
|--------|-----------|------|----------|------|
| dG_separated (REU) | < -15 | -15 ~ -5 | -5 ~ 0 | > 0 |
| dSASA (A^2) | > 1200 | 800~1200 | 500~800 | < 500 |
| sc_value | > 0.70 | 0.65~0.70 | 0.50~0.65 | < 0.50 |
| packstat | > 0.70 | 0.65~0.70 | 0.55~0.65 | < 0.55 |
| dG_density (dG/dSASA×100) | < -2.0 | -2.0 ~ -1.5 | -1.5 ~ -1.0 | > -1.0 |
| delta_unsatHbonds | 0~2 | 3~5 | 6~10 | > 10 |
| nres_int | > 25 | 15~25 | 10~15 | < 10 |
| hbonds_int | > 5 | 2~5 | 1 | 0 |
| Energy Gap (REU) | > 15 | 8~15 | 3~8 | < 3 |

---

## 7. 결과 해석 가이드

### 7.1 핵심 개념

이 파이프라인은 **PPI 결합 사이트 탐색**이 목적이다. 각 클러스터는 하나의 후보 결합 포켓을 나타내며:
- 수렴도가 낮음 ≠ 실패 (여러 사이트가 발견된 것)
- 높은 L_RMSD = 정상 (표면 전체 탐색이므로)
- 에너지가 좋은 사이트부터 생물학적 타당성을 개별 평가

### 7.2 결과 확인 순서

**Step 1: Validation Report** (`docking_validation_report.txt`)
- 10개 품질 체크(C1~C10)의 PASS/WARNING/FAIL 확인
- 핵심 항목: C2_결합에너지, C3_에너지퍼널, C4_인터페이스크기
- C6_사이트탐색은 항상 PASS (다양한 사이트 발견 = 정보적)
- 전체 신뢰도 판정: 높은/보통/낮은/신뢰불가

**Step 2: 정량적 메트릭 비교** (`final_ranking.csv`, `cluster_summary.csv`)
- 파일명 `S-18.32` = dG -18.32 REU (낮을수록 강한 결합)
- dG < -10이면 의미 있는 결합, > -5이면 거의 무의미
- dSASA, sc_value, packstat 등 복합 판단 필요

**Step 3: PyMOL 시각적 분석** (가장 중요)

```bash
# 서버에서 폴더째 로컬 복사
scp -r eunae@node04:/path/to/C-lobe_beta/ ~/Desktop/C-lobe_beta/
cd ~/Desktop/C-lobe_beta/

# 전체 사이트 분포 확인 (chain A=회색, chain B=클러스터별 색상)
pymol 1_OVERVIEW_Clusters.pml

# 개별 사이트 상세 확인 (에너지 좋은 순서대로)
pymol final_result/2_DETAIL_C01.pml

# 최종 랭킹 모델 (B-factor 컬러링)
pymol final_result/view_results.pml

# 전장 원본 PDB 겹치기 (PyMOL 내에서):
#   load /path/to/full_length.pdb
#   align full_length and chain A, C01_M01_S-18 and chain A
```

### 7.3 결과 파일별 용도

| 파일 | 용도 | 핵심 확인 사항 |
|------|------|---------------|
| `docking_validation_report.txt` | 자동 품질 평가 | C1~C10 PASS/FAIL, 전체 신뢰도 |
| `final_ranking.csv` | 정량적 메트릭 비교 | dG, dSASA, sc, packstat 수치 |
| `cluster_summary.csv` | 사이트별 비교 | Population, 평균 에너지 |
| `1_OVERVIEW_Clusters.pml` | 전체 사이트 분포 | 어느 표면에 사이트가 모이는지 |
| `2_DETAIL_C##.pml` | 개별 사이트 수렴도 | 포즈들이 같은 포켓에 모이는지 |
| `view_results.pml` | 최종 구조 확인 | B-factor 기반 에너지 분포 |
| `energy_funnel.png` | 에너지 랜드스케이프 | 퍼널(수렴) vs 평탄(비수렴) |

---

## 8. 모듈 상세

### pipeline_manager.py - PipelineManager 클래스

파이프라인 전체 흐름을 제어하는 싱글 클래스.

```python
PipelineManager(config_file, override_input_pdb=None)
├── __init__()         # 설정 로드, 디렉토리 생성
├── _setup_logging()   # progress_logger (stdout) + internal_logger (file)
├── _validate_config() # 설정값 파싱 및 유효성 검사
├── _compute_l_rmsd()  # 클러스터링용 L_RMSD (clone -> superimpose -> measure)
├── execute()          # 메인 파이프라인 (Step 1~7)
│   ├── Step 1: Relax (캐시 확인 -> FastRelax)
│   ├── Step 2: Global Docking (multiprocessing.Pool)
│   ├── Step 2.5: Fast Scoring & Filtering (4단계 필터 + fallback)
│   ├── Step 3: L_RMSD Greedy Clustering (CoM pre-filter)
│   ├── Step 4: Refinement (DockMCMProtocol)
│   ├── Step 5: Final Scoring (full metrics)
│   ├── Step 6: Selection (round-robin + dedup)
│   └── Step 7: Visualization & Report
├── generate_validation_report()  # 자동 품질 평가 (v2.0: 10개 기준 C1~C10, v1.0: 7개)
├── generate_pymol_script()       # view_results.pml
├── plot_energy_funnel()          # energy_funnel.png
├── generate_cluster_overview()   # 1_OVERVIEW_Clusters.pml
└── generate_per_cluster_views()  # 2_DETAIL_C*.pml
```

### docking.py - 워커 함수

| 함수 | 입력 | 출력 | 설명 |
|------|------|------|------|
| `run_relax_task(pdb_path)` | PDB 파일 경로 | `{status, pdb_data}` | FastRelax |
| `run_global_docking_task(args)` | `(idx, pdb_str, [excl, dist, max])` | `{status, id, pdb_data, Center_*}` 또는 `{status:"rejected"}` | 글로벌 도킹 + 조기 거부 |
| `run_refinement_task(args)` | `(idx, pdb_str, parent, t, r, sub)` | `{status, id, parent, pdb_data}` | DockMCMProtocol 재적용 |
| `_check_excluded_contact(pose, excl, dist)` | Pose, 금지잔기 set, 거리 | 접촉 수 (int) | 금지구역 CA-CA 거리 체크 |

### analysis.py - 스코어링 함수

| 함수 | 설명 |
|------|------|
| `run_fast_scoring_task(args)` | Pass 1: 빠른 스코어링 (dG, dSASA, sc, **total_score**, center, constraint metrics) |
| `run_intermediate_scoring_task(args)` | **Pass 2 (v2.0 신규)**: 비싼 메트릭 (packstat, delta_unsatHbonds, nres_int, hbonds_int). Stage 1 생존자만 대상 |
| `run_scoring_task(args)` | 전체 스코어링 (위 전부 + L_RMSD, B-factor, per-residue CSV, binding residues) |
| `inject_bfactors(pose)` | Residue 에너지를 B-factor에 주입 (모든 원자 순회) |
| `get_residue_energy_csv_string(pose, sfxn)` | Per-residue 에너지 CSV (fa_atr, fa_rep, fa_sol, fa_elec) |
| `analyze_interface_contacts(pose, dist)` | NeighborhoodResidueSelector로 Chain A 접촉 잔기 탐지 |
| `calculate_l_rmsd(mobile, ref)` | Chain A 정렬 후 Chain B CA RMSD |
| `count_excluded_contacts(pose, excl, dist)` | 금지구역 접촉 수 (스코어링 단계용) |
| `calculate_key_contact_ratio(pose, key_B, dist)` | 핵심잔기 접촉 비율 (soft metric) |

### common.py - 유틸리티

| 함수 | 설명 |
|------|------|
| `init_rosetta()` | PyRosetta 초기화 (PID 시드, 멱등, 배너 억제) |
| `setup_worker_logging()` | 워커 프로세스 로거 설정 (멱등) |
| `string_to_pose(pdb_string)` | PDB 문자열 -> Pose (실패 시 None) |
| `pose_to_string(pose)` | Pose -> PDB 문자열 |
| `parse_residue_ranges(range_str)` | `"700-807,840-854"` -> `{700,701,...,807,840,...,854}` |

---

## 9. 핵심 설계 결정 및 주의사항

### PyRosetta 버전 호환성

서버 환경은 네트워크 차단으로 패키지 업데이트 불가. 구버전 호환 방어 코드가 곳곳에 존재:

| 이슈 | 대응 | 위치 |
|------|------|------|
| `get_sc_value()` 미지원 | `getPoseExtraScore` 우선 → fallback `iam.get_sc_value()` | analysis.py |
| `get_packstat()` 미지원 | `getPoseExtraScore` 우선 → `get_packstat()` → `get_interface_packstat()` (3단계) | analysis.py |
| `get_interface_delta_hbond_unsat()` 미지원 | `getPoseExtraScore` 우선 → `hasattr` 체크 → 기본값 99 | analysis.py |
| `get_num_interface_residues()` 미지원 | `getPoseExtraScore` 우선 → `hasattr` 체크 → 기본값 0 | analysis.py |
| `get_interface_hbonds()` 미지원 | `getPoseExtraScore("hbonds_int")` 우선 → `hasattr` 체크 → 기본값 0 | analysis.py |
| `residue_total_energies` 미지원 | `hasattr` 체크 후 조건부 호출 | analysis.py |
| `pyrosetta.init()` 중복 | RuntimeError catch | common.py |
| PyRosetta 배너 | stdout/stderr redirect | common.py |

### Multiprocessing 패턴

- **Step 2 (Global Docking)**: `pool.imap_unordered` - 순서 무관, 최대 처리량
- **Step 2.5 Pass 1 (Fast Scoring)**: `pool.imap` + enumerate - 입출력 1:1 매핑 필수
- **Step 2.5 Pass 2 (Intermediate Scoring, v2.0)**: `pool.imap` + enumerate - 동일 패턴
- **Step 4 (Refinement)**: `pool.imap_unordered` - 순서 무관
- **Step 5 (Final Scoring)**: `pool.imap` + zip - 입출력 1:1 매핑 필수

### 메모리 관리 (대규모 실행)

- `docking_results[i] = None` + `del` + `gc.collect()` 패턴
- PDB 문자열이 수만 개 메모리에 올라가므로, 사용 완료 즉시 해제 필수

### 클러스터링 최적화

- **CoM pre-filter**: Chain B 중심 좌표 거리로 먼 후보 사전 제외 (O(N^2) RMSD 계산 회피)
- **Closest-match**: 첫 번째 매칭이 아닌, 가장 가까운 리더에 배정 (클러스터 품질 향상)

### Graduated Fallback (필터링)

**v2.0** (4단계):
- Level 0: 모든 필터 통과 (충분한 생존자)
- Level 1: 비싼 메트릭만 해제 (dSASA + sc + dG_density 유지)
- Level 2: sc + dG_density 해제 (dSASA만, 임계값 50% 완화)
- Level 3: 모든 Stage 2 해제 → Stage 1 생존자에서 dG 상위 N개

**v1.0** (3단계):
1. sc_value 필터 완화
2. dG 필터 완화
3. 모든 필터 무시, dG 상위 N개

---

## 10. 로깅 시스템

| 로거 | 대상 | 출력 | 용도 |
|------|------|------|------|
| `pipeline_progress` | 메인 프로세스 | stdout (PBS가 pipeline.log로 캡처) | 단계별 진행, 필터 통계, 타이밍 |
| `python_internal` | 메인 + 워커 | `worker_debug.log` (파일) | 디버그, 에러 traceback, PID |

---

## 11. 예상 실행 시간 (32 cores 기준)

| 설정 | Docking | Scoring+Filter | Clustering | Refinement | 총 |
|------|---------|----------------|------------|------------|-----|
| 1K 모델 | ~1h | ~20min | ~5min | ~2h | ~3-4h |
| 20K 모델 | ~12-24h | ~4h | ~30min | ~4h | ~24-36h |
| 100K 모델 | ~60-80h | ~15h | ~1h | ~4h | ~60-120h |

---

## 12. 개발 이력 (주요 커밋)

| 커밋 | 내용 |
|------|------|
| **v2.0 필터링 개편** | |
| `b39ca4c` | hbonds_int 메트릭 추가 + 임계값 문헌 권장값으로 강화 |
| `d6aa898` | get_interface_packstat 3단계 fallback 수정 |
| `a47955b` | 호환성 체크 test pose 체인 분리 수정 |
| `ba598fd` | **v2.0 다중 단계 필터링 파이프라인 대개편** (Stage 1+2, 2-Pass, Graduated Fallback v2.0) |
| `c2e816a` | v1.0-stable 백업 태그 (필터링 개편 전) |
| **v1.0 기반** | |
| `ae91d63` | Tier 1-3 종합 개선 (Constraints, Early Rejection, Round-robin Selection, Dedup) |
| `10c2576` | Validation Report (자동 품질 평가 7개 기준) |
| `169e170` | Per-residue 에너지 분석 + Summary CSV |
| `aa9a4a4` | L_RMSD 클러스터링 + 다중 기준 필터링 구현 |
| `6c05508` | Contact residue cutoff 4.0A -> 6.0A 수정 |
| `4f15e33` | Global Docking 실패 + SlideIntoContact 누락 버그 수정 (핵심 V1.0 수정) |

---

## 13. 추가 참고

- **상세 매뉴얼**: `MANUAL.md` 참조 (알고리즘 상세, L_RMSD 수식, 설정 테이블)
- **호환성 점검**: 새 서버 배포 전 `python check_filter_v2_compat.py` (v2.0) 및 `python check_improvements_compat.py` (v1.0) 실행
- **참고 논문**: Gray (2003), Chaudhury (2011) - Rosetta docking protocol
