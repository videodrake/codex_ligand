# 프로젝트 컨텍스트

## 현재 작업 상태
- 워크플로우: 양쪽 (하네스는 A/B 공통 인프라)
- 현재 작업: 하네스 엔지니어링 구축 Phase 2 스킬 구축 진행 중
- 다음 작업: 나머지 스킬 5개 생성

## 작업 로그
- [2026-04-06] Phase 0-1: 레거시 파일 삭제 완료 — 삭제 3개, 미존재 12개
- [2026-04-06] Phase 0-2: 레거시 참조 정리 — README.md 2개 행 삭제, docs/README.md 2개 행 삭제, tests/test_nightly_review.py 삭제, test_e2e_group7.py 미존재 확인
- [2026-04-06] Phase 0-3: docs/ 스캔 완료 — 활성 7개, 아카이브 35개, paths.py DEPRECATED 없음, input/PPI/prepared/ 활성 참조 4개소 확인
- [2026-04-06] Phase 0-4: CLAUDE.md 재작성 — 설계서 섹션 3 기반 68줄, 3 File System 내용 전부 제거
- [2026-04-06] Phase 1-1: CLAUDE.md 재작성 완료
- [2026-04-06] Phase 1-2: CONTEXT.md 재작성 완료
- [2026-04-06] Phase 2-1: 스킬 생성 — phase-dependencies, bug-history (2/7)

## 최근 결정 사항
- [2026-04-06] 하네스 엔지니어링 적용 시작. 설계서: harness_engineering_design.md

## 발견된 이슈 (미해결)
- input/PPI/prepared/ 가 활성 코드 4개소에서 참조됨 — 데이터/코드 동시 정리 필요 (Phase 0-3에서 발견)
- 참조 문서 4개 미생성: methodology_limitations.md, workflow_comparison_guide.md, phase4_A3_axis_specification.md, 투두리스트.md

## 실패 패턴 (반복 방지)
(아직 없음)
