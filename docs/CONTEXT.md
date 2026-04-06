# 프로젝트 컨텍스트

## 현재 작업 상태
- 워크플로우: 양쪽 (하네스는 A/B 공통 인프라)
- 현재 작업: 하네스 구축 완료, 일상 유지보수 모드 전환
- 다음 작업: 일상 파이프라인 개발/유지보수

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
