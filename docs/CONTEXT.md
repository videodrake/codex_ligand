# 프로젝트 컨텍스트

## 현재 작업 상태
- 워크플로우: **Workflow A + B 모두 완료**
- 현재 작업: 결과 분석 및 해석 단계
- 다음 작업: **WF-A + WF-B 결과 종합, 최종 druggable pocket 후보 확정, 논문/보고서 준비**

## 다음 세션에서 해야 할 일 (에이전트 필독)

### 현재 완료된 것
- Workflow A Phase 1~7 전체 실행 완료 (HPC codex_ligand2)
- PPI 30 seeds (3 상태 × 10 seeds = 600K 모델) 완료
- Vina blind docking (3 수용체 × 3 리간드) 완료
- partner(MYO1D, chain B) 인터페이스 잔기 추출 버그 수정 → Phase 3~7 재실행 완료
- ATP site STRONG 차단 구현 (절대 규칙 #2)
- **Workflow B Phase 1~4 전체 완료** (HPC codex_ligand2)
  - Phase 1: PPI 패치 분석 (orientation filter + cluster consensus + cross-state)
  - Phase 2: Pocket Analysis (fpocket 165→103 merged, 2 primary candidates)
  - Phase 3: Focused Vina (168 jobs, 3 rounds, 156+6+6 ok)
  - Phase 4: Perturbation Scoring (207 candidates scored, 103 shortlisted)

### 핵심 발견 사항

**Vina 결과**: STRONG 포켓 상위 3개(P003, P010, P004)가 모두 ATP 포켓 (overlap 16~18/20). 배제 후 유일한 STRONG = EGFR_170-200 P045 (allosteric 후보, PPI에서 62.7Å 원거리).

**PPI 결과 — EGFR 측 (chain A)**: 3/3 상태에서 반복 출현하는 잔기 22개 확인. 고 occupancy 순: ILE941(avg 0.205), VAL980(0.148), THR940(0.118), PRO992(0.118), PRO937(0.090), GLN982(0.102). C-lobe 표면에 집중.

**PPI 결과 — MYO1D 측 (chain B)**: Sheet 8/9 active face 잔기가 최상위 접촉 — VAL962(avg 0.210), VAL964(0.165), CYS970(0.103), SER971(0.093). Ko et al. 실험 결과와 일치. Sheet 12 잔기는 낮은 occupancy → 구조적 지지 판정 재확인.

### 해결된 미완료 투두 [2026-04-07]

1. **AMBIGUOUS_BAND 검증 — 완료, 0.15→0.10 축소**:
   - WF-A PPI 600개 모델에 retroactive dot product 분석 실시
   - 결과: pass 65.2%, fail 17.5%, ambiguous 17.3% (0.15 기준)
   - 17.3% ambiguous는 과도 → 0.10으로 축소 (ambiguous 10.2%)
   - orientation_filter.py + lightdock_validation.py 양쪽 수정

2. **pending_* 리간드 재점수 — Workflow A N/A**:
   - Workflow A에서 Vina 실행 완료 상태, 모든 사이트에 실제 affinity 값 존재
   - pending_* 개념은 Workflow B Phase 3-4 전용

3. **Centroid 거리 임계값 — 현재 유지**:
   - Cross-receptor 2쌍: P010↔P004(4.52Å), P023↔P045(5.11Å)
   - 6Å/8Å/15Å 임계값 체계 적정, Bootstrap CI 상한이 8Å 근처이나 변경 불필요

4. **Phase 4 축 가중치 — 현재 유지, WF-B 후 재검토 (승인 완료)**:
   - 유일한 STRONG P045가 PPI 0점 (62.67Å 원거리) — 가중치 변경으로 해결 불가
   - PPI 근접 non-ATP 포켓이 없는 구조적 한계
   - Workflow B (PPI-First)로 진행하여 PPI 패치 기반 포켓 탐색이 더 생산적

### Workflow B Phase 1 결과 [2026-04-07]

**Orientation filter** (AMBIGUOUS_BAND=0.10 적용):
- 3GT8_raw: 156/200 pass (78%), 26 fail, 18 ambiguous
- EGFR_160-185: 127/200 pass (64%), 49 fail, 24 ambiguous
- EGFR_170-200: 135/200 pass (68%), 46 fail, 19 ambiguous

**Cluster consensus** (orientation-filtered):
- 3GT8_raw: 17 hotspots, 5 multi-cluster
- EGFR_160-185: 22 hotspots, 12 multi-cluster ← 가장 풍부
- EGFR_170-200: 5 hotspots, 4 multi-cluster

**Cross-state comparison**: 179개 잔기 비교, 130개 핸드오프

**상위 robust hotspots** (3/3 상태 공통):
- ILE941 (occ=1.0), ARG977 (0.75), THR993 (0.71), ARG986 (0.60)
- 모두 C-lobe 표면, WF-A 결과와 일치

**Phase 2 핸드오프**: `phase1_downstream_patch_reference.csv` (130 잔기)

### Workflow B Phase 2 결과 [2026-04-08]

- fpocket: 165 raw pockets → 103 merged
- PPI 관계 분류: 21 rim, 31 allosteric, 51 low_relevance
- Druggability: 2 high (tier_1), 1 medium (tier_2), 100 low
- Cross-state: 13 state_robust, 68 shifted, 17 uncertain, 5 specific
- Phase 3 우선순위: **2 primary**, 21 secondary, 29 exploratory, 51 skip
- P2Rank 미실행 (Java 11 필요, HPC는 Java 8). fpocket 단독 결과.

**Primary candidates (Phase 2 선정)**:
1. **3GT8_raw_PKT07**: rim_candidate, PPI 18.73Å, tier_1, state_robust(3/3), 42 residues, 1593Å³
2. **EGFR_170-200_PKT34**: allosteric_candidate, **PPI 9.37Å**, tier_1, state_shifted(3/3), 62 residues, 1938Å³

### Workflow B Phase 3 결과 [2026-04-08]

- 168 focused Vina docking jobs (52 pockets × 3 ligands, 3 rounds)
- Round 1: 156 ok, 0 error / Round 2-3: 6 ok each (primary 추가 도킹)
- 324 output PDBQT 파일 생성

### Workflow B Phase 4 최종 결과 [2026-04-08]

207 candidates scored (103 pockets × ~2 seeds), 103 shortlisted.

**최종 순위 Top 10**:

| 순위 | Pocket | 분류 | Score | 비고 |
|------|--------|------|-------|------|
| **1** | **3GT8_raw_PKT07** | **rim** | **0.541** | PPI 18.7Å, tier_1, state_robust |
| **2** | **EGFR_170-200_PKT34** | **allosteric** | **0.492** | **PPI 9.4Å**, tier_1 |
| 3 | EGFR_160-185_PKT02 | rim | 0.433 | |
| 4 | 3GT8_raw_PKT10 | allosteric | 0.431 | |
| 5 | EGFR_170-200_PKT17 | rim | 0.430 | |
| 6 | EGFR_170-200_PKT06 | rim | 0.422 | |
| 7 | EGFR_160-185_PKT16 | rim | 0.419 | |
| 8 | 3GT8_raw_PKT05 | rim | 0.416 | |
| 9 | 3GT8_raw_PKT11 | rim | 0.412 | |
| 10 | 3GT8_raw_PKT01 | rim | 0.410 | |

**분류 분포**: rim 23개, allosteric 2개, uncertain 29개, irrelevant 49개

**핵심 결론**:
- Phase 2 primary 2개가 최종 1, 2위 — 파이프라인 일관성 검증됨
- PKT34 (PPI 9.4Å, allosteric)는 WF-A blind Vina로 발견 불가했던 PPI 근접 druggable pocket
- PKT07 (rim)은 PPI 가장자리에서 소분자 교란 가능성 있는 포켓
- 상위 23개 rim 포켓이 PPI 교란 약물 설계의 출발점

### 다음 세션에서 진행할 작업

1. **WF-A + WF-B 결과 종합**:
   - WF-A의 valid_sites.csv와 WF-B Phase 4 결과 교차 비교
   - PKT07/PKT34의 WF-A Vina blind 포켓과의 공간적 일치 확인
   - 최종 druggable pocket 후보 확정 (2~5개)

2. **구조 기반 심층 분석**:
   - PyMOL 시각화: PKT07/PKT34 + PPI hotspot 잔기 매핑
   - 포켓 내부 잔기의 druggability 세부 분석
   - 리간드 결합 포즈 분석 (focused Vina 출력)

3. **교차 상태 결합 부위 심층 분석**:
   - EGFR 측 3/3 상태 공통 잔기 3D 시각화 (PyMOL)
   - 공간적 패치 연속성 확인

### HPC 환경 정보
- 최신 코드: `/work4/hwang/onepack/my_second_project/codex_ligand2`
- 심볼릭 링크: `output/workflow_a` → codex_ligand(원본), `input` → codex_ligand(원본)
- 링크가 끊어질 수 있음 (git pull 시). 확인: `ls output/workflow_a/ input/PPI/`

## 작업 로그
- [2026-04-06] Phase 0-1: 레거시 파일 삭제 완료 — 삭제 3개, 미존재 12개
- [2026-04-06] Phase 0-2: 레거시 참조 정리 — README.md 2개 행 삭제, docs/README.md 2개 행 삭제, tests/test_nightly_review.py 삭제, test_e2e_group7.py 미존재 확인
- [2026-04-06] Phase 0-3: docs/ 스캔 완료 — 활성 7개, 아카이브 35개, paths.py DEPRECATED 없음, input/PPI/prepared/ 활성 참조 4개소 확인
- [2026-04-06] Phase 0-4: CLAUDE.md 재작성 — 설계서 섹션 3 기반 68줄, 3 File System 내용 전부 제거
- [2026-04-06] Phase 1-1: CLAUDE.md 재작성 완료
- [2026-04-06] Phase 1-2: CONTEXT.md 재작성 완료
- [2026-04-06] Phase 2-1: 스킬 생성 — phase-dependencies, bug-history (2/7)
- [2026-04-06] Phase 2-2: 스킬 생성 — ppi-analysis, vina-docking, hpc-operations, scoring-system, testing (7/7 완료)
- [2026-04-06] Phase 3-1: 에이전트 3개 생성 — pipeline-dev, reviewer, science-qa (권한 경계 원칙 포함)
- [2026-04-06] Phase 3-2: 훅 2개 생성 — pre-commit.sh, csv-schema-guard.py
- [2026-04-06] Phase 4: 하네스 구축 완료 — README.md 역할별 재구성, 참조 무결성 검증 통과, 스킬 7 + 에이전트 3 + 훅 2 확인
- [2026-04-07] CLAUDE.md에 "⚠️ 이 환경과 HPC는 완전히 분리되어 있다" 섹션 추가 — 절대 규칙 바로 위, output/ 접근 시도 원천 차단 목적
- [2026-04-07] Workflow A 미완료 투두 4건 해결: AMBIGUOUS_BAND 축소(0.15→0.10), pending_* N/A, centroid 임계값 유지, 축 가중치 유지
- [2026-04-07] Workflow B Phase 1 완료 — _adv_phase1() 구현 + 경로 버그 6건 수정 + HPC 실행 성공
- [2026-04-07] CONTEXT.md 아카이브 규칙 추가 — 작업 로그 50줄 초과 시 docs/archive/로 이동
- [2026-04-08] Workflow B Phase 2 완료 — fpocket 165→103 merged, primary 2개 선정
- [2026-04-08] Workflow B Phase 3 완료 — Focused Vina 168 jobs 성공 (경로/round log 버그 3건 수정 후)
- [2026-04-08] Workflow B Phase 4 완료 — 207 candidates scored, PKT07(rim, 0.541) + PKT34(allosteric, 0.492) 최종 1,2위

## 최근 결정 사항
- [2026-04-08] Workflow B 전체 완료. 최종 후보: PKT07 (rim, 0.541) + PKT34 (allosteric, 0.492)
- [2026-04-08] P2Rank 미사용 결정 — HPC Java 8, P2Rank 2.4.2는 Java 11+ 필요. fpocket 단독으로 충분한 결과 획득
- [2026-04-07] AMBIGUOUS_BAND 0.15→0.10 축소 (사용자 승인). 근거: WF-A 600모델 retroactive 분석에서 17.3%→10.2% ambiguous 감소
- [2026-04-07] Phase 4 축 가중치 현재 유지, Workflow B 완료 후 재검토 (사용자 승인)
- [2026-04-06] 하네스 엔지니어링 적용 시작. 설계서: harness_engineering_design.md
- [2026-04-06] PRODUCTION_N_SEEDS 5→10 확장. generate_configs.py와 run_production.py 양쪽 동기화 필수
- [2026-04-06] HPC 환경을 codex_ligand2로 이전. output/workflow_a와 input/은 codex_ligand(원본)에서 심볼릭 링크로 연결

## 발견된 이슈 (미해결)
- (없음)

## 해결된 이슈
- [2026-04-08] WF-B Phase 2~4 경로/실행 버그 5건:
  1. Phase 2 receptor PDB 경로 — 3파일 input/PPI/phase1/ → input/receptors/, 파일명 receptor_{state}.pdb → {state}.pdb
  2. Phase 2 patch_relationship chain A 하드코딩 — 자동 감지로 변경 (MD cluster chain X 지원)
  3. Phase 2 parse_only=True — fpocket/P2Rank 실행 가능하도록 False로 변경
  4. Phase 3 TG 3.0 phase2_dir 미전달 — rerun_cascade에서 paths.wb_phase2_pocket_analysis() 명시 전달
  5. Phase 3 job_construction PDBQT 우선순위 — resolve_receptor_path()가 PDB를 PDBQT보다 먼저 반환. PDBQT 우선으로 변경. LIGAND_EXTENSIONS도 .pdbqt > .sdf 순으로 변경
  6. Phase 3 run_diverse_docking --execute가 round_log.csv에 결과 미기록 — CSV append 로직 추가
- [2026-04-07] WF-B Phase 1 경로 버그 6건:
  1. orientation_filter metadata 탐색 — docking_* 하위 디렉토리 탐색 추가
  2. orientation_filter runs_base 누락 — wa_phase2_ppi_docking() 명시 전달
  3. orientation_filter source_pdb relative_to 충돌 — resolve() 적용
  4. cluster_consensus receptor PDB — input/PPI/phase1/ → input/receptors/ + chain 자동 감지
  5. compare_across_states CSV 미쓰기 — _write_csv() + generate_comparison_report() 호출 추가
  6. merge_orientation_into_models() 호출 누락 — process_state_orientation 후 병합 추가
- [2026-04-06] input/PPI/prepared/ 활성 참조 4개소 정리 완료

## HPC 배포 주의사항
- codex_ligand2 (최신 repo)에서 실행 시 input/과 output/workflow_a를 codex_ligand(원본)에서 심볼릭 링크해야 함. `mkdir -p output && ln -s ...`
- `_capped` PDB는 일반 `ext_beta_meander` PDB와 **다른 구조** (Chain B 원자 수 상이: 5362 vs 5653). 혼용하면 결과 비일관. 3GT8_raw seed5 `_capped` 결과 폐기 후 일반 PDB로 재실행
- PRODUCTION_N_SEEDS 변경 시 generate_configs.py와 run_production.py **양쪽 모두** 수정해야 함. 한쪽만 바꾸면 config INI 부재로 RuntimeError 발생
- seed별 PBS 제출: `qsub -v STATE=$STATE,SEED=$SEED config/run_ppi_state_seed.pbs` (node 지정: `-l nodes=node06:ppn=16`)
- Workflow B Phase 2~3는 setup→execute 2단계: setup이 스크립트를 생성하고, execute는 `bash <script>` 또는 qsub로 별도 실행해야 함
- fpocket은 conda-forge (`conda install -c conda-forge fpocket`), P2Rank은 Java 11+ 필요 ($HOME/tools/ 수동 설치)
- Focused Vina 도킹은 계산 집약적 — ppn=32로 PBS 제출 권장 (ppn=4 대비 ~8x 빠름)

## 실패 패턴 (반복 방지)
- [2026-04-08] Phase 3 도킹 스크립트가 Vina를 실행했지만 결과를 round_log.csv에 안 씀 → post 단계 진입 불가. 원인: --execute 모드가 결과를 화면 출력만 하고 CSV 저장 안 함. 교훈: **실행 결과는 반드시 파일에 persist해야 다음 단계가 읽을 수 있다**. 화면 출력 ≠ 저장
- [2026-04-08] Phase 3 job_construction이 receptor PDB를 PDBQT보다 먼저 선택 → Vina 168 jobs 전부 실패. 원인: resolve_receptor_path()가 .pdb 존재 시 .pdbqt를 확인 안 함. 교훈: **Vina는 PDBQT만 받는다**. 입력 파일 형식 우선순위를 소비자(Vina) 관점에서 설정
- [2026-04-08] Phase 2 receptor PDB 경로가 3파일에서 동일 패턴으로 잘못됨 (input/PPI/phase1/). 교훈: **레거시 경로 상수가 여러 파일에 복제되어 있으면, 한 곳 고치면 끝이 아니다**. grep으로 전체 검색 필수
- [2026-04-08] Phase 3 TG 3.0이 phase2_dir를 기본값(레거시 경로)으로 사용. 원인: rerun_cascade의 _run_tg30()이 output_dir만 전달하고 phase2_dir를 안 넘김. 교훈: **cascade runner의 thin wrapper가 모든 매개변수를 전달하는지 확인**
- [2026-04-08] fpocket 실행 스크립트를 생성만 하고 PBS에서 자동 실행 안 됨 (Phase 2, Phase 3 동일 패턴). 교훈: **setup→execute 2단계 패턴에서, execute가 자동인지 수동인지 반드시 확인**. 코드가 스크립트를 "생성"하는 것과 "실행"하는 것은 다르다
- [2026-04-08] PBS 스크립트 코어 수(ppn=4)가 도킹 성능에 직접 영향. 교훈: **계산 집약적 작업의 PBS 리소스를 사전에 확인**
- [2026-04-07] _adv_phase1() 구현 시 모듈 함수가 데이터만 반환하고 CSV를 쓰지 않는 패턴을 놓침. 원인: CLI main()에서만 쓰기를 하는 구조를 확인 안 함. 교훈: 새 lane 함수 작성 시 각 모듈의 main()이 하는 일을 확인
- [2026-04-07] HPC 디렉토리 구조(docking_{state}_ext_beta_meander/ 중간 디렉토리)를 코드에서 가정한 구조와 불일치. 교훈: 새 스크립트 작성 시 실제 HPC 디렉토리 구조를 먼저 확인 (find 명령)
- [2026-04-06] CSV 컬럼 추가 시 하위 전파 판단 누락 + DoD 미확인 커밋
- [2026-04-06] BUG-001 관련 변경 시 버그 번호 명시적 인용 누락 + 디렉토리 매핑 누락(pyrosetta_docking/)
- [2026-04-06] 스코어링 가중치 변경을 사람 승인 없이 코드 수정 + 커밋. 원인: 권한 경계 원칙이 scoring-system 스킬에 명시되지 않았음. 해결: pipeline-dev, scoring-system, CLAUDE.md 3곳에 승인 필요 명시
- [2026-04-06] DoD smoke test를 '영향 없음' 판단으로 생략. 원인: DoD 문구가 조건부 생략을 허용하는 것처럼 읽힘. 해결: 무조건 실행으로 문구 강화
- [2026-04-06] PRODUCTION_N_SEEDS=5 상태에서 seed 5~9 실행 시도 → config INI 부재로 "completed without valid outputs" RuntimeError. 원인: generate_configs.py와 run_production.py의 N_SEEDS 불일치. 해결: 양쪽 10으로 통일 + seed 5~9 INI 생성
- [2026-04-06] codex_ligand2 클론 후 input/ 미링크 상태에서 PBS 제출 → FileNotFoundError (TH1 domain.pdb). 원인: git repo에 대용량 PDB 미포함. 해결: 원본 repo에서 심볼릭 링크
- [2026-04-06] 3GT8_raw seed5에 _capped PDB(5653 atoms) 사용 → 일반 PDB(5362 atoms)와 비일관. 원인: 이전 코드에서 다른 입력 PDB 이름 사용. 해결: seed5 결과 폐기 후 일반 PDB로 재실행

## 아카이브 규칙
- 작업 로그가 50줄을 넘으면, 가장 오래된 항목부터 `docs/archive/context_YYYY_MM.md`로 이동
- CONTEXT.md에는 최근 2주치 로그만 유지
- 아카이브 시 "최근 결정 사항"과 "실패 패턴"은 이동하지 않는다 (항상 현행 유지)
- "발견된 이슈"는 해결된 것만 아카이브, 미해결은 유지
- 사용자가 "로그 정리해줘"라고 요청할 때 이 규칙대로 실행
