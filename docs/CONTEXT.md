# 프로젝트 컨텍스트

## 현재 작업 상태
- 워크플로우: Workflow A 완료 (Phase 1~7 전체 실행됨, 30 seeds PPI + Vina)
- 현재 작업: PPI 결합 부위 분석 단계
- 다음 작업: **교차 상태 결합 부위 심층 분석 → Workflow B 진행 여부 판단**

## 다음 세션에서 해야 할 일 (에이전트 필독)

### 현재 완료된 것
- Workflow A Phase 1~7 전체 실행 완료 (HPC codex_ligand2)
- PPI 30 seeds (3 상태 × 10 seeds = 600K 모델) 완료
- Vina blind docking (3 수용체 × 3 리간드) 완료
- partner(MYO1D, chain B) 인터페이스 잔기 추출 버그 수정 → Phase 3~7 재실행 완료
- ATP site STRONG 차단 구현 (절대 규칙 #2)

### 핵심 발견 사항

**Vina 결과**: STRONG 포켓 상위 3개(P003, P010, P004)가 모두 ATP 포켓 (overlap 16~18/20). 배제 후 유일한 STRONG = EGFR_170-200 P045 (allosteric 후보, PPI에서 62.7Å 원거리).

**PPI 결과 — EGFR 측 (chain A)**: 3/3 상태에서 반복 출현하는 잔기 22개 확인. 고 occupancy 순: ILE941(avg 0.205), VAL980(0.148), THR940(0.118), PRO992(0.118), PRO937(0.090), GLN982(0.102). C-lobe 표면에 집중.

**PPI 결과 — MYO1D 측 (chain B)**: Sheet 8/9 active face 잔기가 최상위 접촉 — VAL962(avg 0.210), VAL964(0.165), CYS970(0.103), SER971(0.093). Ko et al. 실험 결과와 일치. Sheet 12 잔기는 낮은 occupancy → 구조적 지지 판정 재확인.

### 다음 세션에서 진행할 분석

1. **교차 상태 결합 부위 심층 분석**:
   - EGFR 측 3/3 상태 공통 잔기를 3D 구조에서 시각화 (PyMOL)
   - 공간적으로 연속된 패치(patch)인지, 분산된 잔기인지 확인
   - occupancy + deltaE 결합 에너지를 함께 분석

2. **Workflow B 진행 여부 판단**:
   - PPI 패치 기반 포켓 탐색 (Workflow B Phase 1~4)
   - Ko et al. sheet 8/9 잔기 3개 이상 확인됨 → Workflow B 진행 가능

3. **투두리스트 미완료 항목** (결과 기반 분석):
   - Orientation filter AMBIGUOUS_BAND 검증 (dot product 분포)
   - pending_* 리간드 지지 수준 재점수
   - Centroid 거리 편향 임계값 리뷰
   - Phase 4 축 가중치 조정 (사람 승인 필수)

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

## 최근 결정 사항
- [2026-04-06] 하네스 엔지니어링 적용 시작. 설계서: harness_engineering_design.md
- [2026-04-06] PRODUCTION_N_SEEDS 5→10 확장. generate_configs.py와 run_production.py 양쪽 동기화 필수
- [2026-04-06] HPC 환경을 codex_ligand2로 이전. output/workflow_a와 input/은 codex_ligand(원본)에서 심볼릭 링크로 연결

## 발견된 이슈 (미해결)
- (없음)

## 해결된 이슈
- [2026-04-06] input/PPI/prepared/ 활성 참조 4개소 정리 완료 — register_pilot_data() 제거(dead code), main.py/postprocess_ppi.py/test 경로 업데이트, data_inventory.md 반영

## HPC 배포 주의사항
- codex_ligand2 (최신 repo)에서 실행 시 input/과 output/workflow_a를 codex_ligand(원본)에서 심볼릭 링크해야 함. `mkdir -p output && ln -s ...`
- `_capped` PDB는 일반 `ext_beta_meander` PDB와 **다른 구조** (Chain B 원자 수 상이: 5362 vs 5653). 혼용하면 결과 비일관. 3GT8_raw seed5 `_capped` 결과 폐기 후 일반 PDB로 재실행
- PRODUCTION_N_SEEDS 변경 시 generate_configs.py와 run_production.py **양쪽 모두** 수정해야 함. 한쪽만 바꾸면 config INI 부재로 RuntimeError 발생
- seed별 PBS 제출: `qsub -v STATE=$STATE,SEED=$SEED config/run_ppi_state_seed.pbs` (node 지정: `-l nodes=node06:ppn=16`)

## 실패 패턴 (반복 방지)
- [2026-04-06] CSV 컬럼 추가 시 하위 전파 판단 누락 + DoD 미확인 커밋
- [2026-04-06] BUG-001 관련 변경 시 버그 번호 명시적 인용 누락 + 디렉토리 매핑 누락(pyrosetta_docking/)
- [2026-04-06] 스코어링 가중치 변경을 사람 승인 없이 코드 수정 + 커밋. 원인: 권한 경계 원칙이 scoring-system 스킬에 명시되지 않았음. 해결: pipeline-dev, scoring-system, CLAUDE.md 3곳에 승인 필요 명시
- [2026-04-06] DoD smoke test를 '영향 없음' 판단으로 생략. 원인: DoD 문구가 조건부 생략을 허용하는 것처럼 읽힘. 해결: 무조건 실행으로 문구 강화
- [2026-04-06] PRODUCTION_N_SEEDS=5 상태에서 seed 5~9 실행 시도 → config INI 부재로 "completed without valid outputs" RuntimeError. 원인: generate_configs.py와 run_production.py의 N_SEEDS 불일치. 해결: 양쪽 10으로 통일 + seed 5~9 INI 생성
- [2026-04-06] codex_ligand2 클론 후 input/ 미링크 상태에서 PBS 제출 → FileNotFoundError (TH1 domain.pdb). 원인: git repo에 대용량 PDB 미포함. 해결: 원본 repo에서 심볼릭 링크
- [2026-04-06] 3GT8_raw seed5에 _capped PDB(5653 atoms) 사용 → 일반 PDB(5362 atoms)와 비일관. 원인: 이전 코드에서 다른 입력 PDB 이름 사용. 해결: seed5 결과 폐기 후 일반 PDB로 재실행
