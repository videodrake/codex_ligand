# 기존 문서 → 하네스 매핑 및 체계화 계획

> CLAUDE.md가 "짧은 지도"로서 docs/를 가리키려면, docs/ 자체가 역할별로 정리되어 있어야 한다.
> 이 계획은 기존 문서들을 하네스 구조에 매핑하고, Claude Code가 실행할 체계화 작업을 정의한다.

---

## 1. 현재 문서 전체 인벤토리

README.md와 테스트 코드에서 확인된 문서들:

```
루트 레벨:
  README.md
  CLAUDE.md                          ← 재작성 대상
  CLAUDE_org.md                      ← 폐기
  PIPELINE_ARCHITECTURE_REPORT.md
  설계의도.md (workflow_a_design_rationale)
  투두리스트.md

docs/:
  PROJECT_USAGE_OVERVIEW.md
  runbook.md
  environment_setup.md
  data_inventory.md
  manual_vina.md
  manual_pyrosetta.md
  manual_execution.md
  phase1_notes.md
  phase4_A3_axis_specification.md
  methodology_limitations.md
  workflow_comparison_guide.md
  workflow_comparison_design.md
  output_path_guide.md
  output_artifact_map.md
  module_separation_analysis.md
  before_after_comparison.md
  ac_coverage_checklist.md
  research_overview_full.md
  pre_qsub_test_line.md
  test_suite_triage.md
  architecture.md
  current_pipeline_status.md
  AI_START_HERE.md
  CONTEXT.md                         ← 내용 재작성
  prd.md                             ← 폐기
  tasks.md                           ← 폐기
  nightly_review_automation.md       ← 폐기
  nightly_incremental_improvement_automation.md  ← 폐기
  archive/                           ← Claude Code가 스캔

config/:
  README.md
```

---

## 2. 역할별 매핑

하네스에서 문서는 5가지 역할 중 하나를 담당한다:
- **ARCHITECTURE** — 시스템이 어떻게 생겼는지 (구조, 설계 의도, 데이터 흐름)
- **OPERATIONS** — 시스템을 어떻게 실행하는지 (실행 가이드, 환경, 명령어)
- **SCIENCE** — 과학적 맥락과 제약 (방법론, 한계, 실험 데이터)
- **QUALITY** — 품질과 검증 (테스트, 출력 검증, 커버리지)
- **TRACKING** — 진행 상태와 미완료 항목 (투두, 컨텍스트, 비교)

### 매핑 테이블

| 문서 | 역할 | 스킬 연결 | CLAUDE.md 참조 |
|------|------|-----------|---------------|
| **ARCHITECTURE** | | | |
| PIPELINE_ARCHITECTURE_REPORT.md | 전체 아키텍처 | phase-dependencies | ✅ 직접 참조 |
| 설계의도.md | 설계 판단 근거 | 전체 | — (스킬에서 참조) |
| docs/architecture.md | 아키텍처 요약 | phase-dependencies | — |
| docs/module_separation_analysis.md | 모듈 분리 분석 | phase-dependencies | — |
| docs/output_path_guide.md | 출력 경로 가이드 | phase-dependencies | — |
| docs/output_artifact_map.md | 출력 아티팩트 맵 | testing | — |
| docs/data_inventory.md | 입출력 인벤토리 | phase-dependencies | ✅ 직접 참조 |
| docs/PROJECT_USAGE_OVERVIEW.md | 프로젝트 사용 개요 | — | — (README에서 참조) |
| **OPERATIONS** | | | |
| docs/runbook.md | 실행 가이드 | hpc-operations | ✅ 직접 참조 |
| docs/environment_setup.md | 환경 설정 | hpc-operations | — |
| docs/manual_execution.md | 수동 실행 명령 | hpc-operations | — |
| docs/manual_vina.md | Vina 매뉴얼 | vina-docking | — (스킬에서 참조) |
| docs/manual_pyrosetta.md | PyRosetta 매뉴얼 | ppi-analysis | — (스킬에서 참조) |
| docs/phase1_notes.md | Phase 1 실행 노트 | ppi-analysis | — |
| docs/pre_qsub_test_line.md | 사전 제출 테스트 | hpc-operations, testing | — |
| config/README.md | Config 파일 의미 | hpc-operations | — |
| docs/current_pipeline_status.md | 현재 파이프라인 상태 | — | — |
| docs/AI_START_HERE.md | 온보딩 순서 | — | — (역할이 CLAUDE.md로 이전) |
| **SCIENCE** | | | |
| docs/methodology_limitations.md | 방법론 한계 | ppi-analysis, vina-docking | ✅ 직접 참조 |
| docs/workflow_comparison_guide.md | A↔B 불일치 해석 | scoring-system | ✅ 직접 참조 |
| docs/workflow_comparison_design.md | 비교 설계 | scoring-system | — |
| docs/phase4_A3_axis_specification.md | A3 축 계산 로직 | scoring-system | ✅ 직접 참조 |
| docs/research_overview_full.md | 연구 개요 전체 | — | — |
| docs/before_after_comparison.md | 전후 비교 | — | — |
| **QUALITY** | | | |
| docs/test_suite_triage.md | 테스트 분류 | testing | — (스킬에서 참조) |
| docs/ac_coverage_checklist.md | AC 커버리지 | testing | — |
| **TRACKING** | | | |
| 투두리스트.md | 미구현 항목 | — | ✅ 직접 참조 |
| docs/CONTEXT.md | 세션 간 메모리 | — | ✅ 직접 참조 |

---

## 3. CLAUDE.md 참조 문서 설계

CLAUDE.md의 "참조 문서" 섹션은 모든 문서를 나열하지 않는다. 에이전트가 작업 유형별로 어디를 찾아가야 하는지 안내하는 "지도"만 둔다.

```markdown
## 참조 문서 (docs/ = source of truth)

구조를 이해하려면:
- PIPELINE_ARCHITECTURE_REPORT.md — 전체 아키텍처, 모듈별 입출력, 데이터 흐름
- docs/data_inventory.md — 입출력 파일 인벤토리

실행하려면:
- docs/runbook.md — 실행 순서, qsub 명령, 결과 확인

과학적 맥락:
- docs/methodology_limitations.md — 방법론 한계 5개 섹션
- docs/workflow_comparison_guide.md — Workflow A↔B 불일치 해석
- docs/phase4_A3_axis_specification.md — Phase 4 A3 축 계산 로직

진행 추적:
- 투두리스트.md — 미구현 항목
- docs/CONTEXT.md — 세션 간 메모리

상세 매뉴얼은 관련 스킬 파일이 안내한다.
```

---

## 4. 스킬 → 문서 연결

각 스킬의 SKILL.md에 "상세 참조" 섹션을 두어, 해당 도메인의 심층 문서를 가리킨다. 에이전트가 스킬을 로딩했을 때 필요하면 문서를 더 읽을 수 있도록.

| 스킬 | 참조할 문서 |
|------|------------|
| vina-docking | docs/manual_vina.md, docs/methodology_limitations.md (Vina scoring 섹션) |
| ppi-analysis | docs/manual_pyrosetta.md, docs/phase1_notes.md, docs/methodology_limitations.md (rigid-body, LightDock 섹션) |
| hpc-operations | docs/runbook.md, docs/environment_setup.md, docs/manual_execution.md, docs/pre_qsub_test_line.md, config/README.md |
| phase-dependencies | PIPELINE_ARCHITECTURE_REPORT.md, docs/data_inventory.md, docs/output_path_guide.md, docs/architecture.md |
| scoring-system | docs/phase4_A3_axis_specification.md, docs/workflow_comparison_guide.md, docs/workflow_comparison_design.md |
| testing | docs/test_suite_triage.md, docs/ac_coverage_checklist.md, docs/output_artifact_map.md |
| bug-history | 설계의도.md (PyRosetta 절대 주의사항 섹션) |

---

## 5. 정리가 필요한 문서들

### 5.1 역할이 중복되는 문서

| 문서 A | 문서 B | 문제 | 처리 |
|--------|--------|------|------|
| docs/AI_START_HERE.md | CLAUDE.md (재작성) | 온보딩 순서 역할이 겹침 | AI_START_HERE.md를 폐기하거나, CLAUDE.md가 이를 대체한다고 명시 |
| docs/architecture.md | PIPELINE_ARCHITECTURE_REPORT.md | 아키텍처 요약 vs 상세 | 유지 — architecture.md가 축약본 역할이면 그대로 두되, 불일치 여부 확인 필요 |
| docs/manual_execution.md | docs/runbook.md | 수동 실행 vs 실행 가이드 | 유지 — 역할이 다름 (manual은 개별 명령, runbook은 전체 순서) |

### 5.2 폐기 대상 (이전 섹션에서 이미 확정)

- docs/prd.md
- docs/tasks.md
- docs/nightly_review_automation.md
- docs/nightly_incremental_improvement_automation.md
- CLAUDE_org.md

### 5.3 Claude Code가 확인해야 할 것

- docs/AI_START_HERE.md — 내용이 CLAUDE.md 재작성과 겹치는지 확인, 겹치면 폐기 또는 CLAUDE.md로 병합
- docs/current_pipeline_status.md — 내용이 최신인지 확인, 오래됐으면 업데이트 또는 투두리스트.md에 흡수
- docs/before_after_comparison.md — 일회성 비교 문서인지, 지속적으로 쓰이는지 확인
- docs/archive/ — 전체 스캔, 현재 코드와 맞지 않는 문서 보고

---

## 6. Claude Code 실행 계획

### Step 1: 폐기 + 참조 정리 (설계서 섹션 11)
이미 정의된 삭제 목록 실행

### Step 2: 문서 역할 확인
```
Claude Code에게:
"docs/ 폴더의 모든 .md 파일을 읽고, 각 문서의 첫 3줄과 마지막 수정일을 보고해줘.
특히 다음을 확인:
1. docs/AI_START_HERE.md가 CLAUDE.md와 역할이 겹치는지
2. docs/current_pipeline_status.md가 최신인지
3. docs/before_after_comparison.md가 일회성인지
4. docs/archive/ 안에 뭐가 있는지"
```

### Step 3: CLAUDE.md 재작성
설계서 섹션 3 기반으로 재작성. 참조 문서 섹션은 위 섹션 3의 "지도" 형태로.

### Step 4: 스킬 파일에 참조 문서 추가
각 스킬의 SKILL.md 끝에 "상세 참조" 섹션 추가 (섹션 4 테이블 기반).

### Step 5: README.md 문서 안내 테이블 업데이트
- 폐기된 문서 행 삭제 (nightly_review 2개)
- 역할별 그룹핑으로 테이블 재구성 (현재는 평탄한 나열)

### Step 6: 문서 drift 초기 점검
```
Claude Code에게:
"docs/ 내 모든 .md 파일에서 다음을 grep:
1. 'nightly_review' — 폐기된 시스템 참조
2. '3 File System' 또는 'stage1.md' 또는 'stage2.md' — 폐기된 시스템 참조
3. 'prd.md' 또는 'tasks.md' — 폐기된 파일 참조
4. 'CLAUDE_org' — 폐기된 파일 참조
발견되면 해당 참조를 제거하거나 수정"
```

---

## 7. 최종 상태 (목표)

```
CLAUDE.md                    ← 짧은 지도 (60줄)
  ├→ 절대 규칙 7개
  ├→ 워크플로우 매핑표
  ├→ Definition of Done
  ├→ 스킬 목록 (7개)
  └→ 참조 문서 (역할별 5그룹, 핵심 7개만)
       │
       ├→ 구조: PIPELINE_ARCHITECTURE_REPORT.md, data_inventory.md
       ├→ 실행: runbook.md
       ├→ 과학: methodology_limitations.md, workflow_comparison_guide.md, phase4_A3_axis_specification.md
       └→ 추적: 투두리스트.md, CONTEXT.md

.claude/skills/
  각 SKILL.md
  └→ "상세 참조" 섹션에서 해당 도메인 문서를 가리킴
       예: vina-docking → manual_vina.md, methodology_limitations.md
       예: hpc-operations → runbook.md, environment_setup.md, manual_execution.md, config/README.md

README.md
  └→ 문서 안내 테이블 (역할별 그룹핑)

docs/
  └→ source of truth (변경 없음, 정리만)
```

에이전트가 작업할 때의 흐름:
1. CLAUDE.md를 읽는다 (자동) → 절대 규칙 + 워크플로우 확인
2. 작업 유형에 따라 스킬이 로딩된다 → 도메인 지식 + 위험한 변경 안내
3. 더 깊은 정보가 필요하면 스킬의 "상세 참조"가 가리키는 docs/를 읽는다
4. docs/CONTEXT.md에서 이전 세션 맥락을 파악한다
