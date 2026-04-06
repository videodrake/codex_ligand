name: vina-docking
description: Vina 도킹 관련 작업 시 로딩. 트리거 — Vina, 도킹, 리간드, PDBQT, exhaustiveness, affinity 언급 시. 비트리거 — PPI/PyRosetta 작업, 스코어링 로직만 변경할 때는 이 스킬이 아님.

## 핵심 정보
- vina_executor.py: prepare_receptor → prepare_ligand → run_vina
- exhaustiveness=384 (기본값 8의 48배, blind docking 70Å+ box 대응)
- 3종 리간드: 173940, 97806, VAX-C12_0 (쌍별 Tanimoto < 0.4)
- 에너지 단위: 항상 kcal/mol

## Workflow별 역할
- Workflow A: Phase 1 (blind, 전체 표면) + Phase 4 (postprocess)
- Workflow B: Phase 3 (focused, 포켓별 집중, budget-aware)

## PDBQT 변환 fallback 순서
Meeko → ADFR → MGLTools → OpenBabel (코드: vina_executor.py)

## 주의사항
- Workflow A blind box는 EGFR 전체를 감싸는 70Å+
- Workflow B focused box는 Phase 2에서 정의된 포켓별 좌표 사용
- 두 모드의 결과를 직접 비교하면 안 됨 (탐색 범위가 다름)

## 이 스킬을 쓰지 말아야 할 때
- Verdict/스코어링 점수 체계만 바꿀 때 → scoring-system 스킬
- Phase 간 CSV 핸드오프 문제 → phase-dependencies 스킬
- PBS 스크립트만 수정할 때 → hpc-operations 스킬

## 상세 참조
- docs/manual_vina.md — AutoDock Vina 상세 매뉴얼
- PIPELINE_ARCHITECTURE_REPORT.md — Vina 모듈 섹션
