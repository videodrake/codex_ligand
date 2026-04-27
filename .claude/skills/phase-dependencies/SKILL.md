name: phase-dependencies
description: Phase 간 전환, 핸드오프, 의존 관계 작업 시 로딩. 트리거 — Phase 간 전환, 핸드오프, 다음 Phase, 입력 참조 언급 시. 비트리거 — 단일 Phase 내부 로직만 바꿀 때, 스코어링 점수 체계만 변경할 때.

## Workflow A 의존 그래프 (독립 → 병합)
Phase 1(Vina Blind) ──────┐
                           ├→ Phase 5(Verdict) → Phase 6(Report) → Phase 7(Validate)
Phase 2(PPI Blind) ───────┘
  └→ Phase 3(PPI Post) ───┘
Phase 4(Vina Post) ────────┘

핵심: Phase 1과 Phase 2는 서로의 결과를 사용하지 않는다.
Phase 5(Verdict)에서 처음으로 병합.

## Workflow B 순차 의존 그래프
Phase 1 → Phase 2 → Phase 3 → Phase 4
전제 조건: Workflow A의 PPI 도킹(Phase 2)이 완료되어야 시작 가능

## 핸드오프 CSV 파일 (Workflow B)
| 구간 | 파일 | 생성 TG |
|------|------|---------|
| Phase 1 → 2 | phase1_downstream_patch_reference.csv | TG 1.6 |
| Phase 2 → 3 | phase3_candidate_pocket_reference.csv | TG 2.6 |
| Phase 3 → 4 | phase4_docking_evidence_reference.csv | TG 3.6 |

## 기존 안전장치 (코드 내장)
- _validate_adv_handoff(): 각 Phase 시작 전 핸드오프 파일 존재 사전 검증
- Vina 가용성 가드: silent all-skip 방지
- Phase 3 cascade 모드별 사전조건 검증

→ 이 안전장치들을 건드리지 않는다. 추가 검증이 필요하면 기존 패턴을 따른다.

## CSV 컬럼 추가 시 전파 판단 규칙
CSV 컬럼을 추가할 때는 '하위 코드가 깨지는가'뿐 아니라
'새 컬럼이 핸드오프 CSV까지 전파되어야 하는가'도 확인한다.
특히 SUMMARY_COLUMNS → EXPORT_COLUMNS 전파 여부를 매번 판단한다.

## 상세 참조
- docs/architecture.md — 전체 아키텍처, 데이터 흐름
- docs/data_inventory.md — 입출력 인벤토리
