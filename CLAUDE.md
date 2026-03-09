# CLAUDE.md - PyRosetta PPI Docking Pipeline

> 이 파일은 Claude Code가 세션 시작 시 자동으로 읽는 프로젝트 컨텍스트 파일입니다.
> 새로운 세션이 시작될 때마다 이 파일을 통해 프로젝트를 빠르게 파악할 수 있습니다.

## 프로젝트 요약

PyRosetta 기반 **단백질-단백질 상호작용(PPI) Global Blind Docking 파이프라인**.
2-chain PDB를 입력받아 Relax -> Global Docking -> Filtering -> Clustering -> Refinement -> Final Scoring을 수행하여 최적 결합 포즈를 예측한다.

- **상태**: v2.0 (필터링 파이프라인 개편, 1K 테스트 진행 중) / v1.0 하위 호환
- **실행 환경**: Linux HPC (PBS/qsub), 32 CPU cores, 네트워크 차단
- **언어**: Python 3.x + PyRosetta (구버전 호환성 필수)

## 파일 구조

```
pipeline_manager.py        # 메인 오케스트레이터 (~1660줄) - 7단계 파이프라인 실행
docking.py                 # 워커: Relax, Global Docking, Refinement (253줄)
analysis.py                # 워커: Scoring, RMSD, Interface 분석 (~590줄)
common.py                  # 유틸리티: PyRosetta 초기화, Pose<->String 변환 (150줄)
config.ini                 # 기본 설정 (1K 모델, 테스트용) - v2.0 필터 포함
config_100k.ini            # 프로덕션 설정 (100K 모델) - v2.0 필터 포함
config_20k.ini             # 중간 설정 (20K 모델)
run_v1.pbs                 # PBS 배치 (v1.0, input_PDB/ 순차처리)
run_v2_test.pbs            # PBS 배치 (v2.0 테스트, 단일 PDB)
check_filter_v2_compat.py  # v2.0 메트릭 호환성 점검 (서버 배포 전)
check_improvements_compat.py # v1.0 개선사항 호환성 점검
input_PDB/                 # 입력 PDB 파일 (EGFR_TH1, EGFR_beta, C-lobe_TH1, C-lobe_beta)
```

## 파이프라인 7단계 흐름

1. **Relax** - FastRelax (ref2015), 결과는 `relaxed_cache/`에 캐싱
2. **Global Docking** - RigidBodyPerturbMover(360deg,100A) -> SlideIntoContact -> DockMCMProtocol
3. **Fast Scoring & Filtering** (v2.0 / v1.0 자동 분기)
   - **v2.0** (`[FilterStage1]` 섹션 존재 시):
     - Pass 1: Fast Scoring (dG, dSASA, sc, total_score)
     - Stage 1: Coarse Energy (제외구역 → dG>0 제거 → total_score 백분위)
     - **Mini Refinement**: Stage 1 생존자에 인터페이스 리패킹 (IncludeCurrent + ExtraRotamers ex1/ex2, 다중 라운드). 리지드바디 도킹 후 불량한 인터페이스 패킹을 개선하여 Stage 2 통과율을 높임.
     - Stage 2 Cheap: dSASA ≥ 500, dG_density ≤ -1.0, sc ≥ 0.50
     - Stage 2 Expensive (Pass 2): packstat ≥ 0.55, delta_unsatHbonds ≤ 8, nres_int ≥ 10, hbonds_int ≥ 0
     - Graduated Fallback v2.0 (Level 0~3)
   - **v1.0** (레거시): dSASA → 제외구역 → dG 백분위 → sc (기존 fallback)
4. **L_RMSD Greedy Clustering** - CoM pre-filter + closest-match, member diversity 보장
5. **Refinement** - DockMCMProtocol 재적용 (미세 조정)
6. **Final Scoring & Selection** - Round-robin 다양성 선택 + L_RMSD 중복 제거
7. **Visualization & Report** - PyMOL 스크립트, Energy Funnel Plot, Validation Report

## 절대 주의사항 (코드 수정 시)

1. **PyRosetta 버전 호환성**: `analysis.py`의 `try-except`/`hasattr` 체크 절대 제거 금지. 서버 PyRosetta가 구버전일 수 있음.
2. **Multiprocessing 순서**: `pool.imap` + `zip` 으로 입출력 1:1 매핑 유지. `imap_unordered`로 바꾸면 인덱스 깨짐.
3. **네트워크 차단**: `pip install`/`conda update` 불가. 현재 설치된 라이브러리로만 동작해야 함.
4. **DockingSlideIntoContact 필수**: 이것 없이 DockMCMProtocol 호출하면 모든 dG가 0.0 (역사적 V1.0 핵심 버그).
5. **FoldTree Setup**: `setup_foldtree(pose, "A_B", movable_jumps)` 역직렬화 후 반드시 재설정 필요.
6. **stdout/stderr 리다이렉트**: `common.py`에서 PyRosetta 배너 억제를 위해 사용. 구조 변경 시 주의.

## Constraints 시스템

- `excluded_residues_A`: Chain A 금지 구역 (막면/다이머 인터페이스). Hard filter.
- `key_residues_B`: Chain B 핵심 잔기 (실험 데이터 기반). Soft bonus (adjusted_dG 계산).
- `enable_early_rejection`: Global Docking 단계에서 DockMCMProtocol **이전에** 금지구역 접촉 검사 (연산 절약).

## 실행 방법

```bash
# 직접 실행
conda activate pyrosetta
python pipeline_manager.py config.ini input_PDB/EGFR_TH1.pdb

# PBS 배치 - v1.0 (input_PDB/ 전체 순차 처리)
qsub run_v1.pbs
qsub -v CONFIG_FILE=config_100k.ini run_v1.pbs

# PBS 배치 - v2.0 테스트 (단일 PDB, 기본: C-lobe_beta)
qsub run_v2_test.pbs
qsub -v INPUT_PDB=input_PDB/EGFR_TH1.pdb run_v2_test.pbs
qsub -v CONFIG_FILE=config_100k.ini,INPUT_PDB=input_PDB/C-lobe_TH1.pdb run_v2_test.pbs

# PBS 배치 - v2.0 전체 PDB 처리
qsub -v INPUT_PDB=ALL run_v2_test.pbs

# 호환성 점검 (서버 배포 전)
python check_filter_v2_compat.py
```

## 출력 디렉토리 구조

```
<PDB_NAME>/
  filter_passed/         # 필터 통과 구조 (F0001_S-15.23.pdb)
  cluster_results/       # 클러스터 대표 (C01_M01_S-15.23.pdb) + cluster_summary.csv
  final_result/          # 최종 랭킹 + PyMOL 스크립트 + 리포트
    Rank01_C01_M01_S-18.45.pdb   # 최종 구조 (S뒤 숫자 = dG)
    Rank01_C01_M01_Energies.csv  # Per-residue 에너지
    view_results.pml             # 최종 랭킹 모델 (B-factor 컬러링)
    2_DETAIL_C01.pml             # 사이트별 상세 뷰
    docking_validation_report.txt # PPI 사이트 탐색 리포트
  final_ranking.csv      # 종합 랭킹 테이블 (모든 메트릭 포함)
  energy_funnel.png      # L_RMSD vs dG 산점도
  1_OVERVIEW_Clusters.pml # 전체 결합 사이트 비교 (클러스터별 색상 구분)
```

> **PML 파일은 상대경로** 사용 → 폴더째 로컬로 복사하면 PyMOL에서 바로 열림.
> 모든 PML에 chain A 기준 자동 정렬 포함.

## 결과 해석 가이드

### 파이프라인의 목적

이 파이프라인은 **단일 최적 포즈 예측이 아닌, PPI 결합 사이트 탐색**이 목적이다.
각 클러스터는 하나의 **후보 결합 사이트(포켓)**를 나타내며, 수렴도가 낮다고 해서 실패가 아니라 여러 후보 사이트가 발견된 것이다.

### 결과 해석 순서

1. **Validation Report 확인** (`docking_validation_report.txt`)
   - 10개 품질 체크(C1~C10)의 PASS/WARNING/FAIL 확인
   - C6_사이트탐색은 항상 PASS (다양한 사이트 = 정보적)
   - 핵심 실패 항목(C2_결합에너지, C3_에너지퍼널, C4_인터페이스크기)에 주목

2. **final_ranking.csv 분석**
   - 파일명의 `S-18.32` = dG_separated (-18.32 REU)
   - dG < -10이면 의미 있는 결합, > -5이면 거의 무의미
   - 여러 메트릭을 종합 판단: dG + dSASA + sc_value + packstat

3. **PyMOL 시각적 분석** (핵심 단계)
   - `1_OVERVIEW_Clusters.pml`: 전체 사이트 분포 확인 (클러스터별 색상 구분)
   - `2_DETAIL_C01.pml` ~ `C20.pml`: 개별 사이트의 포즈 수렴도 확인
   - 에너지가 좋은 상위 클러스터부터 순서대로 확인
   - 각 사이트가 생물학적으로 타당한 위치인지 판단

4. **cluster_summary.csv 비교**
   - 사이트별 에너지/메트릭 수치 비교
   - Population(멤버 수)이 많은 사이트 = 더 많이 수렴된 사이트

### PyMOL 시각화 워크플로우

```bash
# 서버에서 폴더째 로컬로 복사
scp -r eunae@node04:/path/to/C-lobe_beta/ ~/Desktop/C-lobe_beta/

# 전체 사이트 분포 확인 (chain A=회색, chain B=클러스터별 색상)
cd ~/Desktop/C-lobe_beta/
pymol 1_OVERVIEW_Clusters.pml

# 개별 사이트 상세 확인 (에너지 좋은 순서대로)
pymol final_result/2_DETAIL_C01.pml
pymol final_result/2_DETAIL_C02.pml

# 최종 랭킹 모델 (B-factor 컬러링)
pymol final_result/view_results.pml

# 전장 원본 PDB와 겹치기 (수동)
# PyMOL 내에서:
#   load /path/to/full_length.pdb
#   align full_length and chain A, C01_M01_S-18 and chain A
```

### 결과 파일별 용도 요약

| 파일 | 용도 | 핵심 확인 사항 |
|------|------|---------------|
| `docking_validation_report.txt` | 자동 품질 평가 | C1~C10 PASS/FAIL, 전체 신뢰도 판정 |
| `final_ranking.csv` | 정량적 메트릭 비교 | dG, dSASA, sc, packstat 수치 |
| `cluster_summary.csv` | 사이트별 비교 | Population, 평균 에너지 |
| `1_OVERVIEW_Clusters.pml` | 전체 사이트 분포 | 어느 표면에 사이트가 모이는지 |
| `2_DETAIL_C##.pml` | 개별 사이트 수렴도 | 포즈들이 같은 포켓에 모이는지 |
| `energy_funnel.png` | 에너지 랜드스케이프 | 퍼널 형태(수렴) vs 평탄(비수렴) |

### L_RMSD 해석 주의

Global Blind Docking에서 **높은 L_RMSD는 정상**이다. 표면 전체를 탐색하므로 서로 먼 위치에 사이트가 발견되는 것이 자연스럽다. L_RMSD가 낮다면 대부분의 포즈가 한 곳에 수렴한 것이므로 오히려 해당 사이트의 신뢰도가 높다는 의미이다.

## 스코어링 메트릭 해석

| Metric | 좋은 값 | 설명 |
|--------|---------|------|
| dG_separated | < -10 REU | 인터페이스 결합 에너지 (낮을수록 강한 결합) |
| dSASA | > 800 A^2 | 결합 시 묻히는 표면적 (클수록 넓은 인터페이스) |
| sc_value | > 0.65 | Shape Complementarity (높을수록 좋은 기하적 맞물림) |
| packstat | > 0.65 | 원자 패킹 밀도 (높을수록 조밀) |
| dG_density | < -1.5 | dG/dSASA×100 (에너지 밀도, 낮을수록 효율적) |
| delta_unsatHbonds | < 5 | 결합 시 미충족 수소결합 수 (적을수록 양호) |
| nres_int | > 15 | 인터페이스 잔기 수 (클수록 넓은 접촉면) |
| hbonds_int | ≥ 1 | 인터페이스 수소결합 수 (많을수록 강한 상호작용) |
| L_RMSD | 낮을수록 | Relaxed 대비 Ligand RMSD |

## v2.0 필터링 아키텍처

### 자동 분기
- `config.ini`에 `[FilterStage1]` 섹션이 **있으면** → v2.0 경로
- **없으면** → v1.0 경로 (완전 하위 호환)

### 2-Pass 설계
- **Pass 1**: 모든 decoy에 `run_fast_scoring_task()` (dG, dSASA, sc, total_score) — 비용 ~0
- **Mini Refinement** (Stage 1 ~ Stage 2 사이): Stage 1 생존자에 인터페이스 사이드체인 리패킹
  - `[MiniRefinement]` 섹션: `enabled`, `mode` (pack_only/pack_min), `n_rounds` (기본 3)
  - IncludeCurrent (현재 로타머 보존) + ExtraRotamersGeneric (ex1+ex2 chi 확장) + 다중 라운드
  - 리파인 전후 dG/dSASA/sc 비교 → Validation Report 섹션 1.5에 보고
  - 리지드바디 도킹 후 불량한 패킹을 보완하여 Stage 2 통과율 향상
- **Pass 2**: Stage 1 생존자(~5-10%)만 `run_intermediate_scoring_task()` (packstat, unsatHb, nres_int, hbonds_int) — 비싼 메트릭

### Graduated Fallback v2.0
- Level 0: 모든 필터 통과 (충분한 생존자)
- Level 1: 비싼 메트릭만 해제 (dSASA + sc + dG_density 유지)
- Level 2: sc + dG_density 해제 (dSASA만, 임계값 50% 완화)
- Level 3: 모든 Stage 2 해제 → Stage 1 생존자에서 dG 상위 N개

### 새 메트릭 안전 기본값 (추출 실패 시)
- packstat → 0.0, delta_unsatHbonds → 99, nres_int → 0, hbonds_int → 0
- 모두 필터에서 **탈락** 방향 (보수적) → Fallback이 안전망 역할

## 기술적 디테일 (참고용)

- PyRosetta 초기화: PID 기반 랜덤 시드, `-mute all`, 멱등성 보장
- 클러스터링: CoM(Center of Mass) pre-filter로 O(N^2) RMSD 계산 회피
- 메모리 관리: `docking_results[i] = None` + `gc.collect()` 로 대규모 실행 시 메모리 절약
- Boltzmann-weighted cluster probability: Validation Report에 포함 (kT=1.0 REU)
- v2.0 Validation Report: 10개 품질 체크 (C1~C10) + 새 메트릭 분포 섹션
