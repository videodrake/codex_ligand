name: testing
description: 테스트 작성/실행/검증 관련 작업 시 로딩. 트리거 — 테스트, 검증, validate, 회귀, smoke test 언급 시. 비트리거 — 기능 구현 자체만 할 때 (구현 후 테스트는 별도 단계로).

## validate.py (Workflow A 출력 검증)
4개 그룹, 8개 함수:
- 그룹 1: 파일 존재 + ID 일관성 + 추적성 + coverage
- 그룹 2: CSV 스키마 회귀 검사
- 그룹 3: 잔기 번호 일관성 + 알려진 변이 확인
- 그룹 4: 핸드오프 준비 확인
종료 코드: 0(통과) / 1(경고) / 2(실패)

## pytest 마커
- smoke: 빠른 기본 검증 (도킹 없이 실행 가능)
- full: PyRosetta/Vina 필요한 통합 테스트

## mock 전략
- PyRosetta 없는 환경: mock pose 객체 사용
- Vina 없는 환경: 사전 생성된 결과 파일로 후처리 테스트

## paths.py 수정 후
반드시 pytest tests/ -m smoke --tb=short 전체 실행

## 상세 참조
- docs/test_suite_triage.md — 테스트 분류 가이드
- docs/pre_qsub_test_line.md — 사전 제출 테스트 절차
