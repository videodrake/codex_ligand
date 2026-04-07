# 프로젝트 컨텍스트

## 현재 작업 상태
- 워크플로우: Workflow A 완료 (Phase 1~7 전체 실행됨, 30 seeds PPI + Vina)
- 현재 작업: Workflow B Phase 1 완료, Phase 2 준비 단계
- 다음 작업: **Workflow B Phase 2 실행 (Pocket Analysis, TG 2.0~2.7)**

## 다음 세션에서 해야 할 일 (에이전트 필독)

### 현재 완료된 것
- Workflow A Phase 1~7 전체 실행 완료 (HPC codex_ligand2)
- PPI 30 seeds (3 상태 × 10 seeds = 600K 모델) 완료
- Vina blind docking (3 수용체 × 3 리간드) 완료
- partner(MYO1D, chain B) 인터페이스 잔기 추출 버그 수정 → Phase 3~7 재실행 완료
- ATP site STRONG 차단 구현 (절대 규칙 #2)
- **Workflow B Phase 1 완료** (TG 1.1.5~1.7 전체 실행, HPC codex_ligand2)

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

### 다음 세션에서 진행할 작업

1. **Workflow B Phase 2 (Pocket Analysis, TG 2.0~2.7)**:
   - 입력: `phase1_downstream_patch_reference.csv`
   - PPI 패치 기반 druggable pocket 탐색
   - 실행: `qsub config/run_adv_phase2.pbs`

2. **이후**: Phase 3 (Focused Vina) → Phase 4 (Perturbation Scoring)

3. **교차 상태 결합 부위 심층 분석** (병행 가능):
   - EGFR 측 3/3 상태 공통 잔기 3D 시각화 (PyMOL)
   - 공간적 패치 연속성 확인, occupancy + deltaE 결합 분석

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
- [2026-04-07] Workflow B Phase 1 완료 — _adv_phase1() 구현 + 경로 버그 4건 수정 + HPC 실행 성공
- [2026-04-07] CONTEXT.md 아카이브 규칙 추가 — 작업 로그 50줄 초과 시 docs/archive/로 이동

## 최근 결정 사항
- [2026-04-07] AMBIGUOUS_BAND 0.15→0.10 축소 (사용자 승인). 근거: WF-A 600모델 retroactive 분석에서 17.3%→10.2% ambiguous 감소
- [2026-04-07] Phase 4 축 가중치 현재 유지, Workflow B 완료 후 재검토 (사용자 승인)
- [2026-04-06] 하네스 엔지니어링 적용 시작. 설계서: harness_engineering_design.md
- [2026-04-06] PRODUCTION_N_SEEDS 5→10 확장. generate_configs.py와 run_production.py 양쪽 동기화 필수
- [2026-04-06] HPC 환경을 codex_ligand2로 이전. output/workflow_a와 input/은 codex_ligand(원본)에서 심볼릭 링크로 연결

## 발견된 이슈 (미해결)
- (없음)

## 해결된 이슈
- [2026-04-07] WF-B Phase 1 경로 버그 4건:
  1. orientation_filter metadata 탐색 — docking_* 하위 디렉토리 탐색 추가
  2. orientation_filter runs_base 누락 — wa_phase2_ppi_docking() 명시 전달
  3. orientation_filter source_pdb relative_to 충돌 — resolve() 적용
  4. cluster_consensus receptor PDB — input/PPI/phase1/ → input/receptors/ + chain 자동 감지
  5. compare_across_states CSV 미쓰기 — _write_csv() + generate_comparison_report() 호출 추가
  6. merge_orientation_into_models() 호출 누락 — process_state_orientation 후 병합 추가
- [2026-04-06] input/PPI/prepared/ 활성 참조 4개소 정리 완료 — register_pilot_data() 제거(dead code), main.py/postprocess_ppi.py/test 경로 업데이트, data_inventory.md 반영

## HPC 배포 주의사항
- codex_ligand2 (최신 repo)에서 실행 시 input/과 output/workflow_a를 codex_ligand(원본)에서 심볼릭 링크해야 함. `mkdir -p output && ln -s ...`
- `_capped` PDB는 일반 `ext_beta_meander` PDB와 **다른 구조** (Chain B 원자 수 상이: 5362 vs 5653). 혼용하면 결과 비일관. 3GT8_raw seed5 `_capped` 결과 폐기 후 일반 PDB로 재실행
- PRODUCTION_N_SEEDS 변경 시 generate_configs.py와 run_production.py **양쪽 모두** 수정해야 함. 한쪽만 바꾸면 config INI 부재로 RuntimeError 발생
- seed별 PBS 제출: `qsub -v STATE=$STATE,SEED=$SEED config/run_ppi_state_seed.pbs` (node 지정: `-l nodes=node06:ppn=16`)

## 실패 패턴 (반복 방지)
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
