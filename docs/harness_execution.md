# 하네스 구축 실행 가이드

이 문서를 Claude Code에서 열고 단계별로 실행한다.
함께 참조할 문서: `harness_design.md`

---

## 실행 순서

```
Phase 0: 레거시 정리              ✅ 완료
Phase 1: CLAUDE.md + CONTEXT.md   ✅ 완료
Phase 2: 스킬 구축               ← 현재
Phase 3: 에이전트 + 훅
Phase 4: 문서 체계화 + 최종 점검
```

---

## Phase 0: 레거시 정리 ✅

완료 사항:
- 폐기 파일 삭제 (nightly_review 3개 삭제, 나머지 12개는 이미 미존재)
- 깨진 참조 수정 (README.md, docs/README.md, test_nightly_review.py)
- 문서 상태 스캔 완료 — 활성 docs/ 11개, archive/ 35개

---

## Phase 1: CLAUDE.md + CONTEXT.md ✅

완료 사항:
- CLAUDE.md 재작성 (472줄 → 68줄)
- docs/CONTEXT.md 재작성 (5개 섹션 구조)

---

## Phase 2: 스킬 구축

### 현재 프로젝트의 활성 문서 현황

스킬의 "상세 참조"는 실제 존재하는 문서만 가리켜야 한다.

```
활성 문서 (docs/):
  PROJECT_USAGE_OVERVIEW.md, runbook.md, environment_setup.md,
  data_inventory.md, manual_vina.md, manual_pyrosetta.md,
  phase1_notes.md, pre_qsub_test_line.md, test_suite_triage.md,
  CONTEXT.md, README.md

루트 레벨:
  CLAUDE.md, architecture.md, archive/design_intent.md

config/:
  README.md

미존재 (archive에도 없음):
  methodology_limitations.md, workflow_comparison_guide.md,
  phase4_A3_axis_specification.md, archive/todo.md,
  output_path_guide.md, workflow_comparison_design.md,
  ac_coverage_checklist.md, module_separation_analysis.md,
  manual_execution.md
```

### 2-1 시작 프롬프트 (핵심 스킬 2개)

```
하네스 스킬 구축을 시작한다.
harness_design.md의 섹션 4를 기반으로 스킬 파일을 생성해줘.

먼저 .claude/skills/ 디렉토리 구조를 만들고,
가장 중요한 2개부터 시작:

1. .claude/skills/phase-dependencies/SKILL.md
   — 설계서 섹션 4.4 내용 기반
   — 끝에 "상세 참조" 섹션 추가:
     - architecture.md (전체 아키텍처, 데이터 흐름)
     - docs/data_inventory.md (입출력 인벤토리)

2. .claude/skills/bug-history/SKILL.md
   — 설계서 섹션 4.7 내용 기반
   — 끝에 "상세 참조" 섹션 추가:
     - archive/design_intent.md (PyRosetta 절대 주의사항 + 설계 판단 근거)

완료 후 docs/CONTEXT.md 작업 로그에 기록
```

### 2-2 시작 프롬프트 (나머지 스킬 5개)

```
나머지 5개 스킬을 생성해줘.
harness_design.md 섹션 4의 각 스킬 내용을 기반으로.
"상세 참조"에는 실제 존재하는 문서만 넣는다.

3. .claude/skills/ppi-analysis/SKILL.md (섹션 4.2)
   상세 참조:
   - docs/manual_pyrosetta.md
   - docs/phase1_notes.md
   - architecture.md (Phase 1 PPI, LightDock 섹션)

4. .claude/skills/vina-docking/SKILL.md (섹션 4.1)
   상세 참조:
   - docs/manual_vina.md
   - architecture.md (Vina 모듈 섹션)

5. .claude/skills/hpc-operations/SKILL.md (섹션 4.3)
   상세 참조:
   - docs/runbook.md
   - docs/environment_setup.md
   - docs/pre_qsub_test_line.md
   - config/README.md

6. .claude/skills/scoring-system/SKILL.md (섹션 4.5)
   상세 참조:
   - architecture.md (Phase 4 Scoring + Verdict 섹션)
   - archive/design_intent.md (3축/4축 설계 의도)

7. .claude/skills/testing/SKILL.md (섹션 4.6)
   상세 참조:
   - docs/test_suite_triage.md
   - docs/pre_qsub_test_line.md

완료 후 docs/CONTEXT.md 작업 로그에 기록
```

---

## Phase 3: 에이전트 + 훅

### 3-1 시작 프롬프트

```
에이전트 정의와 훅을 생성해줘.
harness_design.md 섹션 5, 7을 기반으로.

에이전트:
1. .claude/agents/pipeline-dev.md (섹션 5.1)
2. .claude/agents/reviewer.md (섹션 5.2)
3. .claude/agents/science-qa.md (섹션 5.3)

훅:
4. .claude/hooks/pre-commit.sh (섹션 7.1)
5. .claude/hooks/csv-schema-guard.py (섹션 7.2)

에이전트 파일 앞에 섹션 5.0의 권한 경계 원칙을 포함해줘:
- "에이전트 자율 가능" 목록 (버그 수정, 테스트 추가, 문서 업데이트 등)
- "반드시 사람 승인 필요" 목록 (CSV 스키마 변경, 스코어링 변경, 워크플로우 구조 변경 등)

완료 후 docs/CONTEXT.md 작업 로그에 기록
```

---

## Phase 4: 문서 체계화 + 최종 점검

### 4-1 시작 프롬프트

```
하네스 구축 마무리 작업이다.

1. README.md 문서 안내 테이블을 역할별로 재구성해줘.
   현재 존재하는 문서만 포함한다:

   구조 이해:
   - architecture.md — 전체 아키텍처
   - docs/PROJECT_USAGE_OVERVIEW.md — 프로젝트 개요
   - docs/data_inventory.md — 입출력 인벤토리

   실행:
   - docs/runbook.md — 실행 가이드
   - docs/environment_setup.md — 환경 설정
   - docs/manual_vina.md — Vina 매뉴얼
   - docs/manual_pyrosetta.md — PyRosetta 매뉴얼

   테스트/검증:
   - docs/test_suite_triage.md — 테스트 분류
   - docs/pre_qsub_test_line.md — 사전 제출 테스트

   설정:
   - config/README.md — Config 파일 의미

   AI 에이전트:
   - CLAUDE.md — Claude Code 컨텍스트

   폐기된 문서(nightly_review 등)가 남아있으면 제거

2. 최종 점검:
   - 각 스킬의 "상세 참조"에 나온 문서가 실제로 존재하는지 확인
   - CLAUDE.md의 참조 문서가 실제로 존재하는지 확인
   - 폐기된 시스템 참조가 아직 남아있는지 프로젝트 전체 grep:
     "nightly_review", "3 File System", "CLAUDE_org",
     "docs/prd.md", "docs/tasks.md"

3. 하네스 구조 요약 보고:
   - .claude/ 디렉토리 트리를 보여줘
   - CLAUDE.md 줄 수 확인
   - 스킬 7개, 에이전트 3개, 훅 2개 존재 확인

완료 후 docs/CONTEXT.md:
- 작업 로그에 "Phase 4: 하네스 구축 완료" 기록
- 현재 작업 상태를 "하네스 구축 완료, 일상 유지보수 모드 전환"으로 업데이트
```

---

## 작업 로깅 규칙

모든 Phase에 공통 적용.

```
docs/CONTEXT.md 업데이트 규칙:
- 세션 시작 시: "현재 작업 상태" 섹션을 읽고 맥락을 파악한다
- 작업 완료 시: "작업 로그"에 날짜 + 한 줄 요약을 추가한다
- 중요한 결정 시: "최근 결정 사항"에 이유와 함께 기록한다
- 이슈 발견 시: "발견된 이슈"에 추가한다
- 실패 시: "실패 패턴"에 유형 + 원인 + 해결 방법을 기록한다
```
