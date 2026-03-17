# CLAUDE.md - PyRosetta PPI Docking Pipeline

> 이 파일은 Claude Code가 세션 시작 시 자동으로 읽는 프로젝트 컨텍스트 파일입니다.
> 새로운 세션이 시작될 때마다 이 파일을 통해 프로젝트를 빠르게 파악할 수 있습니다.

## 프로젝트 요약

EGFR-MYO1D 상호작용 결합 사이트 탐색을 위한 **이중 워크플로우 파이프라인**.
AutoDock Vina(소분자 blind docking)와 PyRosetta(PPI global blind docking) 결과를 통합하여 후보 결합 포켓을 식별한다.

- **상태**: Workflow A 자동화 완료 (3 states × 5 seeds, 300K PPI 모델), Workflow B 모듈 구현 완료 (통합 자동화 미완)
- **실행 환경**: Linux HPC (PBS/qsub), 32 CPU cores, 네트워크 차단
- **언어**: Python 3.x + PyRosetta (구버전 호환성 필수)

## 두 가지 워크플로우

### Workflow A: Standard Production (자동화 완료)

Vina blind docking과 PPI blind docking을 **독립적으로** 수행 후 결과 통합.

```
Phase 1: Vina blind docking       (vina/)
Phase 2: PPI global blind docking (pyrosetta_docking/, phase1/)
Phase 3: PPI 후처리               (ppi/)
Phase 4: Vina 후처리              (vina/ postprocess chain)
Phase 5: 사이트 판정              (verdict.py)
Phase 6: 리포트                   (report.py)
Phase 7: 검증                     (validate.py)
```

진입점: `run_production.py` (Phase 1→7 자동 순차), `main.py` (대화형 CLI)

### Workflow B: Advanced PPI-First Pipeline (모듈 구현 완료, 통합 자동화 미완)

PPI 결과를 기반으로 포켓을 **순차적으로 좁혀가며** 분석. 각 Phase 출력이 다음 Phase 입력.

```
[전제: Workflow A Phase 2 완료]
  → Adv Phase 2: 포켓 분석 & 제안       (phase2/, TG 2.0→2.7)
  → Adv Phase 3: 포켓 유도 focused 도킹 (phase3/, TG 3.0→3.7)
  → Adv Phase 4: 교란 분석 & 최종 스코어 (phase4/, TG 4.0→4.6)
```

## 파일 구조

```
main.py                          # 대화형 CLI 진입점 (932줄)
run_production.py                # Workflow A PBS 배치 오케스트레이터 (1609줄)

egfr_pipeline/
  config.py                      # YAML/JSON 설정 로드
  runtime.py                     # 런타임 자원 해석, 환경 감지
  verdict.py                     # 3축 증거 통합 사이트 판정 (1791줄)
  report.py                      # Vina+PPI 통합 리포트
  validate.py                    # 출력 품질 검증
  output_organizer.py            # steps/ 디렉토리 정리
  step_view.py                   # Step 인덱스 & 심볼릭 링크 뷰
  residue_utils.py               # 잔기 파싱 유틸리티

  pyrosetta_docking/             # ── PPI Global Blind Docking 엔진 ──
    pipeline_manager.py          # 코어 오케스트레이터 (4079줄) - 7단계 PPI 파이프라인
    scoring.py                   # dG, dSASA, sc, packstat, per-residue 에너지 (1057줄)
    run_metadata.py              # 실행 메타데이터, 스코어 CSV 빌더 (448줄)
    movers.py                    # FastRelax, RigidBodyPerturb, DockMCM, MiniRefine (414줄)
    pyrosetta_init.py            # PyRosetta 초기화, Pose↔String 변환 (230줄)
    logging_config.py            # 중앙 로깅 설정 (127줄)

  vina/                          # ── AutoDock Vina 도킹 & 분석 ──
    vina_executor.py             # Vina CLI 오케스트레이터 (2813줄)
    pocket_cluster.py            # k-means + Union-Find 포켓 클러스터링 (600줄)
    pocket_stability.py          # Bootstrap 안정성 분석 (405줄)
    cross_receptor.py            # 교차 수용체 포켓 비교 (282줄)
    parse_poses.py               # PDBQT → pose_table.csv (259줄)
    phase3_bridge.py             # Tier 1 포켓 → Phase 3 후보 (253줄)
    pocket_summary.py            # 포켓 통계 요약 (226줄)
    param_sweep.py               # 파라미터 민감도 스윕 (223줄)
    pose_contacts.py             # 접촉 잔기 추출 (167줄)

  ppi/                           # ── PPI 준비 & 후처리 ──
    pyrosetta_extract.py         # PPI 잔기 추출 (691줄)
    prepare_dimer_pdb.py         # Dimer PDB 준비, chain 번호 매기기 (410줄)
    postprocess_ppi.py           # 도킹 후 chain 복원, 잔기 추출 (318줄)
    afm_extract.py               # AlphaFold-Multimer PPI 추출 (195줄)
    pbs_submit.py                # PBS qsub 래퍼 (70줄)

  phase1/                        # ── Workflow A Phase 2: PPI 인터페이스 매핑 ──
    prepare_inputs.py            # 전체 키나아제 도메인 PDB 준비 (1063줄)
    lightdock_validation.py      # LightDock 교차 검증 (1311줄)
    cluster_consensus.py         # 클러스터 인터페이스 합의 (817줄)
    orientation_filter.py        # 방향 인식 필터링 (754줄)
    review_report.py             # Phase 1 리뷰 리포트 (661줄)
    extract_interface.py         # 수용체측 인터페이스 잔기 추출 (544줄)
    compare_states.py            # 다중 상태 패치 비교 (530줄)
    pilot_comparison.py          # 파일럿 데이터 비교 (378줄)
    launch_docking.py            # 도킹 런처 (351줄)
    standardize_scores.py        # 스코어 표준화 (257줄)
    generate_configs.py          # .ini 설정 자동 생성 (239줄)

  phase2/                        # ── Workflow B Adv Phase 2: 포켓 분석 ──
    pocket_proposal.py           # fpocket/P2Rank 기반 포켓 제안 (807줄)
    patch_ingestion.py           # Phase 1 패치 참조 수집 (734줄)
    cross_state_alignment.py     # 교차 상태 포켓 정렬 (583줄)
    review_report.py             # Phase 2 리뷰 (551줄)
    patch_relationship.py        # 패치 연결성/중첩 분류 (550줄)
    druggability_confidence.py   # 약물성 신뢰도 스코어링 (520줄)
    phase3_export.py             # Phase 3 후보 테이블 내보내기 (486줄)
    pocket_merge.py              # 후보 포켓 정규화/병합 (440줄)
    rerun_cascade.py             # Cascade 러너 TG 2.0→2.7 (172줄)

  phase3/                        # ── Workflow B Adv Phase 3: 포켓 유도 도킹 ──
    run_diverse_docking.py       # 다양성 인식 Vina 실행 (718줄)
    budget_policy.py             # 검색 예산 정책 (562줄)
    pocket_reference_ingestion.py # Phase 2 포켓 참조 수집 (507줄)
    pose_attribution.py          # 포즈 파싱 & 포켓 귀속 (495줄)
    phase4_export.py             # Phase 4 내보내기 (474줄)
    job_construction.py          # 도킹 job 구성 (413줄)
    review_report.py             # Phase 3 리뷰 (389줄)
    diversity_validation.py      # 포켓 점유율 & 다양성 검증 (384줄)
    rerun_cascade.py             # Cascade 러너 TG 3.0→3.7 (284줄)

  phase4/                        # ── Workflow B Adv Phase 4: 교란 분석 ──
    evidence_ingestion.py        # 다중 Phase 증거 수집 (703줄)
    score_framework.py           # 4축 스코어 프레임워크 (549줄)
    state_interpretation.py      # 상태 강건성 해석 (465줄)
    mechanistic_classification.py # 메커니즘 분류 (445줄)
    perturbation_scoring.py      # 교란 관련성 스코어링 (441줄)
    final_report.py              # 최종 리포트 (381줄)
    presentation_summary.py      # 프레젠테이션 요약 (373줄)
    review_report.py             # 리뷰 리포트 (266줄)
    rerun_cascade.py             # Cascade 러너 TG 4.0→4.6 (145줄)

  md/                            # ── MD 분석 (선택적) ──
    gromacs_analysis.py          # GROMACS 궤적 분석 (1273줄)
    ligand_contacts.py           # 단백질-리간드 접촉 분석 (507줄)

config/
  example-project.yaml           # 메인 프로젝트 설정 (Vina + PPI 통합)
  phase1/*.ini                   # Phase 1 PyRosetta 설정 (3 test + 15 production)
  run_production.pbs             # Workflow A 프로덕션 PBS 진입점
  run_advanced_pipeline.pbs      # Workflow B 전체 PBS 오케스트레이터
  run_adv_phase*.pbs             # Workflow B 개별 lane PBS (6개)
  run_lightdock.pbs              # Phase 1 LightDock 검증 PBS
  run_lightdock_test.pbs         # Phase 1 LightDock 테스트 PBS
  run_pre_qsub_checks.pbs       # 사전 검증 PBS 진입점
input/PPI/phase1/                # Phase 1 입력 PDB (monomer receptor + partner)
input/PPI/prepared/              # 레거시 dimer 준비 PDB
```

## PPI 도킹 7단계 흐름 (pipeline_manager.py 내부)

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

0. **실행 환경 대원칙**: 파이프라인의 모든 무거운 연산(Vina, PyRosetta PPI, LightDock 등)은 **반드시 qsub를 통해 HPC 서버에서 실행**한다. 이 Codespace(개발 환경)에서는 코드 수정, 테스트, 세부 조정만 수행한다. **절대로 개발 환경에서 도킹을 직접 실행하거나 실행을 제안하지 말 것.** 도킹 실행이 필요한 경우 항상 PBS 스크립트(`qsub`)를 통한 서버 제출 방법을 안내해야 한다.
1. **PyRosetta 버전 호환성**: `scoring.py`의 `try-except`/`hasattr` 체크 절대 제거 금지. 서버 PyRosetta가 구버전일 수 있음.
2. **Multiprocessing 순서**: `pool.imap` + `zip` 으로 입출력 1:1 매핑 유지. `imap_unordered`로 바꾸면 인덱스 깨짐.
3. **네트워크 차단**: `pip install`/`conda update` 불가. 현재 설치된 라이브러리로만 동작해야 함.
4. **DockingSlideIntoContact 필수**: 이것 없이 DockMCMProtocol 호출하면 모든 dG가 0.0 (역사적 V1.0 핵심 버그).
5. **FoldTree Setup**: `setup_foldtree(pose, "A_B", movable_jumps)` 역직렬화 후 반드시 재설정 필요.
6. **stdout/stderr 리다이렉트**: `pyrosetta_init.py`에서 PyRosetta 배너 억제를 위해 사용. 구조 변경 시 주의.

## Constraints 시스템

- `excluded_residues_A`: Chain A 금지 구역 (막면/다이머 인터페이스). Hard filter.
- `key_residues_B`: Chain B 핵심 잔기 (실험 데이터 기반). Soft bonus (adjusted_dG 계산).
- `enable_early_rejection`: Global Docking 단계에서 DockMCMProtocol **이전에** 금지구역 접촉 검사 (연산 절약).

## 실행 방법

> **대원칙**: 모든 도킹/연산은 반드시 qsub를 통해 HPC 서버에서 실행한다.

### Workflow A: Standard Production (Vina blind + PPI blind → 통합)

```bash
# 올인원 (precheck → production 순차)
PRECHECK_JOB=$(qsub config/run_pre_qsub_checks.pbs)
qsub -W depend=afterok:${PRECHECK_JOB} config/run_production.pbs

# 또는 lane별 병렬 제출 (Vina와 PPI 동시 실행)
qsub config/run_vina_cpu.pbs
qsub -v STATE=3GT8_raw,SEED=0 config/run_ppi_state_seed.pbs
# ... (3 states × 5 seeds)
qsub -W depend=afterok:<vina> config/run_vina_postprocess.pbs
qsub -W depend=afterok:<ppi_jobs> config/run_ppi_postprocess.pbs
qsub -W depend=afterok:<vina_post>:<ppi_post> config/run_finalize.pbs
```

### Workflow B: Advanced PPI-First Pipeline (PPI → 포켓 분석 → Focused 도킹 → 스코어링)

```bash
# 전체 자동 (PBS 의존성 체인으로 Phase 1~4 순차 제출)
qsub config/run_advanced_pipeline.pbs

# Phase 2부터 시작
qsub -v ADV_FROM=2 config/run_advanced_pipeline.pbs

# 개별 Phase 수동 제출
qsub config/run_adv_phase1.pbs                                      # PPI 분석
qsub -W depend=afterok:<job> config/run_adv_phase2.pbs              # Pocket cascade
qsub -W depend=afterok:<job> config/run_adv_phase3_setup.pbs        # 도킹 계획
qsub -v ROUND=0 -W depend=afterok:<job> config/run_adv_phase3_execute.pbs  # Vina 실행
qsub -W depend=afterok:<job> config/run_adv_phase3_post.pbs         # 분석
qsub -W depend=afterok:<job> config/run_adv_phase4.pbs              # 통합 스코어링
```

### 기타

```bash
# Phase 1 LightDock 검증
qsub config/run_lightdock.pbs                          # 전체 state
qsub config/run_lightdock_test.pbs                     # 테스트

# 로컬 검증/보고서 (개발 환경에서 가능)
python main.py -c config/example-project.yaml validate
python main.py -c config/example-project.yaml organize
```

## Step별 출력 구조 (steps/)

프로덕션 실행 완료 후 `output/{project}/steps/` 아래에 자동 생성됨.
각 디렉토리에는 해당 스텝의 모든 결과 파일에 대한 심볼릭 링크가 포함됨.

```
output/{project}/steps/
├── STEP_INDEX.txt           # 전체 인덱스
├── 01_vina_docking/         # Step 1: Vina 원시 포즈 (.pdbqt)
│   ├── 3GT8_raw/            →  {project}/3GT8_raw/
│   ├── EGFR_160-185/        →  {project}/EGFR_160-185/
│   └── EGFR_170-200/        →  {project}/EGFR_170-200/
├── 02_ppi_docking/          # Step 2: PPI 도킹 결과 (전체)
│   ├── 3GT8_raw/            →  phase1_ppi/3GT8_raw/prod_seed0/
│   ├── EGFR_160-185/        →  phase1_ppi/EGFR_160-185/prod_seed0/
│   └── EGFR_170-200/        →  phase1_ppi/EGFR_170-200/prod_seed0/
├── 03_ppi_evidence/         # Step 3: PPI 인터페이스 증거
│   ├── ppi_pyrosetta_residues.csv
│   ├── ppi_pyrosetta_summary.csv
│   └── ...
├── 04_vina_analysis/        # Step 4: Vina 포켓 분석
│   ├── vina_pocket_table.csv
│   ├── vina_pose_table.csv
│   ├── vina_drug_pocket_map.csv
│   └── ...
├── 05_site_verdict/         # Step 5: 사이트 판정
│   ├── valid_sites.csv
│   └── cross_method_agreement.csv
├── 06_report/               # Step 6: 보고서
│   ├── project_report.txt
│   └── combined_residue_evidence.csv
└── 07_validation/           # Step 7: 검증
    ├── validation_status.json
    └── validation_summary.txt
```

> 심볼릭 링크이므로 원본 수정 없이 공간 절약.
> 서버에서 로컬로 복사 시: `scp -r steps/ local_dest/` 또는 `cp -rL steps/ dest/` (링크 해제)

## 출력 디렉토리 구조 (PPI 도킹 개별 결과)

```
<PDB_NAME>/
  seed_complete.json     # 시드 완료 마커 (재실행 시 스킵 판단, --force로 초기화)
  scored_all_models.csv  # 전체 스코어링 모델 메트릭 + filter_status (재도킹 없이 재분석 가능)
  scored_stage2_models.csv # Stage 2 대상 모델의 비싼 메트릭 (packstat, unsatHb 등)
  filter_thresholds.csv  # 필터 임계값 기록 (재현성 보장)
  all_scored_summary.csv # 에너지 퍼널용 데이터 (dG, dSASA, sc, total_score + L_RMSD)
  logs/                  # 파이프라인 로그 (pipeline.log + workers.log)
  filter_passed/         # 필터 통과 구조 (F0001_S-15.23.pdb)
  cluster_results/       # 클러스터 대표 + 전체 멤버십
    cluster_summary.csv          # 대표 모델 메트릭 (centroid, dSASA 분해, 에너지 분포 통계 포함)
    cluster_membership.csv       # 전체 model→cluster 매핑 (비대표 포함)
    dropped_candidates.csv       # 중복제거로 제외된 후보
  final_result/          # 최종 랭킹 + PyMOL 스크립트 + 리포트
    Rank01_C01_M01_S-18.45.pdb   # 최종 구조 (S뒤 숫자 = dG)
    Rank01_C01_M01_Energies.csv  # Per-residue 에너지 (9항: fa_atr/rep/sol/elec + hbond 4항 + fa_dun)
    Rank01_C01_M01_InterfaceEnergies.csv  # 인터페이스 잔기 ΔE (9항 분해)
    Rank01_C01_M01_ContactPairs.csv      # 잔기 쌍별 최소 거리 (A-B 접촉 네트워크)
    view_results.pml             # 최종 랭킹 모델 (B-factor 컬러링)
    2_DETAIL_C01.pml             # 사이트별 상세 뷰
    docking_validation_report.txt # PPI 사이트 탐색 리포트 (I_RMSD 분포 포함)
  final_ranking.csv      # 종합 랭킹 (dSASA_polar/hydrophobic, I_RMSD, centroid 포함)
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
| dSASA_polar | — | dSASA 중 극성 기여분 (A^2) |
| dSASA_hydrophobic | — | dSASA 중 소수성 기여분 (A^2) |
| sc_value | > 0.65 | Shape Complementarity (높을수록 좋은 기하적 맞물림) |
| packstat | > 0.65 | 원자 패킹 밀도 (높을수록 조밀) |
| dG_density | < -1.5 | dG/dSASA×100 (에너지 밀도, 낮을수록 효율적) |
| delta_unsatHbonds | < 5 | 결합 시 미충족 수소결합 수 (적을수록 양호) |
| nres_int | > 15 | 인터페이스 잔기 수 (클수록 넓은 접촉면) |
| hbonds_int | ≥ 1 | 인터페이스 수소결합 수 (많을수록 강한 상호작용) |
| I_RMSD | 낮을수록 | 인터페이스 잔기 기준 RMSD (결합 부위 구조 변화) |
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
- v2.0 Validation Report: 10개 품질 체크 (C1~C10) + 새 메트릭 분포 섹션 + I_RMSD 분포
- Per-residue 에너지 CSV: 9개 가중 에너지 항 (fa_atr, fa_rep, fa_sol, fa_elec, hbond_sr_bb, hbond_lr_bb, hbond_bb_sc, hbond_sc, fa_dun)
- ContactPairs CSV: 잔기 쌍별 최소 거리 (CB-CB 기본, GLY는 CA, 실패 시 heavy-atom 최소)
- scored_all_models.csv: 전체 모델의 Pass 1 메트릭 + filter_status (재도킹 없이 재분석 가능)
- filter_thresholds.csv: 사용된 필터 임계값 + 입출력 카운트 (재현성 보장)
- Per-seed 완료 마커: `seed_complete.json`으로 시드별 완료 추적. 재실행 시 완료된 시드 스킵, `--force`로 마커 초기화. 기존 `final_ranking.csv`만 있는 시드는 자동 backfill (하위호환)
- PPI_TARGETS 동적 생성: `RECEPTOR_STATES` × `PRODUCTION_N_SEEDS`로 `_build_ppi_targets()` 생성 (run_production.py). `--status` 시 시드별 상태 그리드 출력
