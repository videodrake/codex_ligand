# EGFR-MYO1D Docking Pipeline

EGFR-MYO1D PPI 교란 약물 포켓 탐색 파이프라인.
Vina(소분자) + PyRosetta(PPI) 이중 증거 통합, 3개 EGFR 구조 상태 교차 비교.

## 절대 규칙

1. **HPC 전용 실행**: 도킹/연산 코드를 직접 실행하지 않는다. PBS 스크립트를 생성하고 qsub 명령을 안내한다.
2. **ATP 포켓 배제**: is_atp_site=True인 포켓을 STRONG 판정하면 안 된다. ATP 결합 유지 + 활성 소실이 실험적 사실이다.
3. **워크플로우 확인 의무**: "Phase N 수정" 요청 시 반드시 Workflow A/B 중 어느 쪽인지 확인한다. 같은 번호가 다른 모듈을 가리킨다. 확인 없이 작업을 시작하지 않는다.
4. **디렉토리-워크플로우 매핑**: 코드 수정 시 이 매핑을 따른다:
   - `egfr_pipeline/vina/` → Workflow A (Phase 1 + 4)
   - `egfr_pipeline/ppi/` → Workflow A (Phase 2 + 3)
   - `egfr_pipeline/phase1/` → Workflow B Phase 1
   - `egfr_pipeline/phase2/` → Workflow B Phase 2
   - `egfr_pipeline/phase3/` → Workflow B Phase 3
   - `egfr_pipeline/phase4/` → Workflow B Phase 4
   - `egfr_pipeline/pyrosetta_docking/` → Workflow A Phase 2 + Workflow B Phase 1 (PPI 도킹 엔진, 양쪽에서 공유)
   - `verdict.py`, `report.py`, `validate.py` → Workflow A Phase 5~7
5. **CSV 스키마 보존**: 기존 CSV 출력의 컬럼명/타입을 변경할 때는 하위 Phase의 ingestion 코드와 validate.py를 함께 수정한다.
6. **paths.py 보호**: `egfr_pipeline/paths.py`를 수정하면 전체 Phase의 경로 해석이 바뀐다. 수정 후 반드시 모든 Phase의 smoke test를 실행한다.
7. **잔기 번호 체계**: PDB 잔기 번호(author numbering)를 사용한다. 내부 인덱스(0-based)와 혼동하면 전체 분석이 틀어진다.
8. **스코어링 가중치 변경 승인**: verdict.py의 축별 점수 배분(vina_max, ppi_max, cross_max) 및 score_framework.py의 A1~A4 가중치 변경은 과학적 판단이므로 반드시 사람 승인 후 진행한다. 승인 없이 코드를 수정하지 않는다.

## 워크플로우 구분

| Phase | Workflow A (Blind) | Workflow B (PPI-First) |
|-------|-------------------|----------------------|
| 1 | Vina Blind (`vina/`) | PPI 분석 (`phase1/`, TG 1.0~1.6) |
| 2 | PPI Global Blind (`ppi/`) | Pocket Analysis (`phase2/`, TG 2.0~2.7) |
| 3 | PPI Postprocess (`ppi/`) | Focused Vina (`phase3/`, TG 3.0~3.6) |
| 4 | Vina Postprocess (`vina/`) | Perturbation Scoring (`phase4/`, TG 4.0~4.6) |
| 5~7 | Verdict → Report → Validate | — |

## Definition of Done

코드 변경을 완료로 간주하려면:
- `pytest -m smoke` 전체 통과
- 변경된 Phase의 `validate.py` 검증 통과
- CSV 스키마 변경 시 하위 ingestion 코드 동시 수정 확인
- `paths.py` 변경 시 전체 smoke test 통과
- 이 체크리스트를 확인하기 전에 커밋하지 않는다.

## 스킬 (.claude/skills/)

작업 시 관련 스킬을 먼저 읽는다:
- **vina-docking** — Vina 도킹 파라미터, PDBQT 변환, blind/focused 구분
- **ppi-analysis** — PyRosetta PPI, orientation filter, 실험 데이터 매핑
- **hpc-operations** — PBS/qsub 스크립트 생성, lane 목록, 모드 구분
- **phase-dependencies** — Phase 간 핸드오프 CSV 계약, 의존 그래프
- **scoring-system** — Verdict 3축(A)/4축(B) 체계, 판정 원칙
- **testing** — validate.py 4그룹, pytest 마커, mock 전략
- **bug-history** — PyRosetta/Vina 코드 수정 전 반드시 확인 (역사적 버그 4건)

## 독립 모듈

- `egfr_pipeline/md/` — GROMACS MD 분석 (선택적, MDAnalysis 필요). 핵심 파이프라인과 독립.

## 참조 문서

| 역할 | 문서 |
|------|------|
| 구조 | `PIPELINE_ARCHITECTURE_REPORT.md`, `docs/data_inventory.md` |
| 실행 | `docs/runbook.md` |
| 추적 | `docs/CONTEXT.md` |

## 실험적 근거

1. **ATP 결합 유지**: EGFR kinase의 ATP 포켓은 MYO1D 결합 후에도 유지됨 → ATP site를 교란 포켓으로 판정 금지
2. **Ko et al. alanine substitution**: MYO1D TH1 도메인 sheet 8/9 잔기가 active face → PPI hotspot에 이 잔기 3개 미만이면 Workflow B 중단 조건
3. **리간드 다양성**: 3종 리간드(173940, 97806, VAX-C12_0)는 쌍별 Tanimoto < 0.4 → 구조적 편향 방지 목적, 임의 교체 금지
