name: hpc-operations
description: PBS/qsub/서버 관련 작업 시 로딩. 트리거 — qsub, PBS, 서버, HPC, 프로덕션 실행 언급 시. 비트리거 — 로컬에서 도는 단위 테스트, 코드 로직만 바꿀 때는 이 스킬이 아님.

## 대원칙
모든 도킹/연산은 반드시 qsub로 HPC 서버에서 실행한다.
에이전트는 PBS 스크립트를 생성/수정만 하고, 직접 실행하지 않는다.

## 실행 방법
- Workflow A: qsub config/run_production.pbs (또는 lane별 병렬)
- Workflow B: qsub config/run_advanced_pipeline.pbs
- 개별 lane: run_production.py --lane {lane_name}

## Lane 목록 (14개)
vina-cpu, ppi, ppi-post, vina-post, finalize, status,
vina-gpu, phase3-gpu,
adv-phase1, adv-phase2, adv-phase3-setup, adv-phase3-execute,
adv-phase3-post, adv-phase4

## 모드 구분
- dry-run: 검증만 수행 (실제 도킹 안 함)
- setup: 서버 실행용 스크립트만 생성
- execute: 실제 도킹 실행

## Phase 3 cascade (Workflow B)
rerun_cascade.py가 setup → execute → post 3모드를 순차 관리
모드별 사전조건: setup → job table, post → round log

## 상세 참조
- docs/runbook.md — 실행 가이드
- docs/environment_setup.md — 환경 설정
- docs/pre_qsub_test_line.md — 사전 제출 테스트 절차
- config/README.md — Config 파일 의미 (YAML, INI, PBS 설명)
