name: ppi-analysis
description: PPI/PyRosetta 관련 작업 시 로딩. 트리거 — PPI, PyRosetta, interface, patch, orientation filter, MYO1D 언급 시. 비트리거 — Vina 소분자 도킹만 다룰 때, MD 분석 시에는 이 스킬이 아님.

## 핵심 정보
- PyRosetta PPI 도킹: 3 states × 5 seeds = 15 세트, 각 20K 모델 = 300K 총
- orientation_filter: sheet 8/9의 active-face normal vs receptor 방향 dot product
  - 양수 = pass (active face가 receptor 향함)
  - 음수 = fail (뒤집힘)
  - consensus 계산에는 pass 모델만 사용

## 실험 데이터 매핑
- Ko et al. alanine substitution: sheet 8/9 잔기 = active face
- PPI hotspot에 이 잔기가 3개 미만이면 Workflow B 중단 조건

## LightDock 검증
- PyRosetta와 독립적인 교차 검증 수단
- LightDock 결과가 PyRosetta와 일치하면 신뢰도 상승
- 단, LightDock 자체의 한계 있음

## 위험한 코드 변경
- DockingSlideIntoContact 누락 → 모든 dG가 0.0 (V1.0 역사적 버그)
- FoldTree: 역직렬화 후 반드시 setup_foldtree 재설정 필요
- excluded_residues_A: 막면/다이머 인터페이스 금지 구역, hard filter
- key_residues_B: 실험 데이터 기반, soft bonus (adjusted_dG)
- enable_early_rejection: DockMCMProtocol 전 금지구역 접촉 검사 (연산 절약)

## 이 스킬을 쓰지 말아야 할 때
- Verdict/스코어링 점수 체계만 바꿀 때 → scoring-system 스킬
- Phase 간 CSV 핸드오프 문제 → phase-dependencies 스킬
- PBS 스크립트만 수정할 때 → hpc-operations 스킬

## 상세 참조
- docs/manual_pyrosetta.md — PyRosetta PPI 도킹 상세 매뉴얼
- docs/phase1_notes.md — Phase 1 참고 노트 (실행, 샘플링, 필터, LightDock)
- PIPELINE_ARCHITECTURE_REPORT.md — Phase 1 PPI, LightDock 섹션
