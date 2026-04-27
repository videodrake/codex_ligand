# 문서 인덱스

## 핵심 문서

| 문서 | 설명 |
|------|------|
| [PROJECT_USAGE_OVERVIEW.md](PROJECT_USAGE_OVERVIEW.md) | 프로젝트 사용 개요 — 각 도구의 역할과 워크플로우 |
| [architecture.md](architecture.md) | 전체 아키텍처 상세 보고서 — 모든 모듈/알고리즘/입출력 |
| [runbook.md](runbook.md) | 실행 가이드 — qsub 명령, 실행 순서, 결과 확인 |
| [environment_setup.md](environment_setup.md) | 환경 설정 — conda, PyRosetta, 서버 설정 |
| [harness_design.md](harness_design.md) | Claude Code 하네스 엔지니어링 설계 |
| [harness_execution.md](harness_execution.md) | 하네스 실행 가이드 |

## 참고 문서

| 문서 | 설명 |
|------|------|
| [manual_vina.md](manual_vina.md) | AutoDock Vina 도킹 상세 매뉴얼 |
| [manual_pyrosetta.md](manual_pyrosetta.md) | PyRosetta PPI 도킹 상세 매뉴얼 |
| [phase1_notes.md](phase1_notes.md) | Phase 1 참고 노트 (실행, 샘플링, 필터, LightDock, 핸드오프) |
| [data_inventory.md](data_inventory.md) | 입출력 데이터 인벤토리 |
| [../config/README.md](../config/README.md) | Config 파일 의미 (YAML, INI, PBS) |

## 개발/테스트

| 문서 | 설명 |
|------|------|
| [pre_qsub_test_line.md](pre_qsub_test_line.md) | 사전 제출 테스트 절차 |
| [test_suite_triage.md](test_suite_triage.md) | 테스트 분류 가이드 |

## 아카이브

`archive/` 디렉토리에는 다음이 포함됨:
- 구 버전 문서 (구 architecture.md, data_flow_guide.md 등 → 새 architecture.md로 대체)
- Phase PRD/Task 기획 문서 (prd_phase_*.md, tasks_phase_*.md)
- 구 버전 Phase 1 노트 (phase1_*.md → phase1_notes.md로 병합)
- 히스토리 문서 (GEMINI.md, handoff docs 등)
- Manuscript 섹션별 초안 (manuscript_sections/section1~6.md, 병합본은 docs/manuscript_draft.md)
- 구 설계의도/투두/매핑 문서 (design_intent.md, todo.md, document_mapping_plan.md)
