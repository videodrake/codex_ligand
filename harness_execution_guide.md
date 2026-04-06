# 하네스 구축 실행 가이드

이 문서를 Claude Code에서 열고 단계별로 실행한다.
함께 참조할 문서: `harness_engineering_design.md`, `document_mapping_plan.md`

---

## 실행 순서

```
Phase 0: 레거시 정리
  0-1. 폐기 파일 삭제
  0-2. 깨진 참조 수정
  0-3. 문서 상태 스캔 + 보고
      ↓
Phase 1: CLAUDE.md 재작성
  1-1. CLAUDE.md 재작성 (짧은 지도 + 절대 규칙 + 참조 문서)
  1-2. docs/CONTEXT.md 재작성 (세션 메모리 + 작업 로그)
      ↓
Phase 2: 스킬 구축
  2-1. skill-phase-dependencies
  2-2. skill-bug-history
  2-3. skill-ppi-analysis
  2-4. skill-vina-docking
  2-5. skill-hpc-operations
  2-6. skill-scoring-system
  2-7. skill-testing
      ↓
Phase 3: 에이전트 + 훅
  3-1. pipeline-dev 에이전트
  3-2. reviewer 에이전트
  3-3. science-qa 에이전트
  3-4. pre-commit.sh 훅
  3-5. csv-schema-guard.py 훅
      ↓
Phase 4: 문서 체계화
  4-1. README.md 문서 테이블 업데이트
  4-2. 스킬 → 문서 참조 연결 확인
  4-3. 문서 drift 최종 점검
```

---

## Phase 0: 레거시 정리

### 0-1 시작 프롬프트

```
이 프로젝트에 하네스 엔지니어링을 적용하기 위한 레거시 정리를 시작한다.
아래 파일들을 삭제해줘. 삭제 전에 각 파일이 실제로 존재하는지 확인하고,
존재하는 것만 삭제한 뒤 결과를 보고해줘.

삭제 대상:
- templates/stage1.md
- templates/stage2.md
- templates/ (폴더, 비어있으면)
- projects/ (폴더 전체)
- docs/prd.md
- docs/tasks.md
- .claude/commands/recover.md
- .claude/commands/execute.md
- .claude/commands/review.md
- .claude/commands/test.md
- .claude/commands/ (폴더, 비어있으면)
- scripts/nightly_review.py
- docs/nightly_review_automation.md
- docs/nightly_incremental_improvement_automation.md
- CLAUDE_org.md

삭제 완료 후 docs/CONTEXT.md에 아래 형식으로 기록해줘:

## 작업 로그
- [오늘 날짜] Phase 0-1: 레거시 파일 삭제 완료 — 삭제 N개, 미존재 N개
```

### 0-2 시작 프롬프트

```
Phase 0-1에서 삭제한 파일들을 참조하던 곳을 수정해줘.

1. README.md 문서 안내 테이블에서 nightly_review 관련 2개 행 삭제
2. tests/test_e2e_group7.py의 TestDocumentExistence에서
   "docs/prd.md", "docs/tasks.md" 항목 제거
3. tests/test_nightly_review.py — 파일 전체가 폐기된 nightly_review 테스트이면 삭제
4. 프로젝트 전체에서 아래 키워드를 grep해서 깨진 참조를 찾아줘:
   - "nightly_review"
   - "stage1.md" 또는 "stage2.md"
   - "3 File System"
   - "CLAUDE_org"
   - "docs/prd.md" 또는 "docs/tasks.md"
   발견되면 해당 참조를 제거하거나 맥락에 맞게 수정

완료 후 docs/CONTEXT.md 작업 로그에 기록
```

### 0-3 시작 프롬프트

```
docs/ 폴더의 현재 상태를 스캔해줘.

1. docs/ 아래 모든 .md 파일 목록을 보여줘
2. docs/archive/ 폴더가 있으면 안에 뭐가 있는지 보여줘
3. 다음 문서들의 상태를 확인해줘:
   - docs/AI_START_HERE.md — 존재 여부, 있으면 첫 5줄 보여줘 (CLAUDE.md와 역할 중복 가능)
   - docs/current_pipeline_status.md — 존재 여부, 있으면 마지막 수정 내용이 언제 기준인지
   - docs/before_after_comparison.md — 일회성 비교인지 지속 문서인지
4. input/PPI/prepared/ — 다른 코드에서 참조하는지 grep 확인
5. paths.py의 DEPRECATED 함수 — 다른 모듈에서 호출하는지 grep 확인

결과를 보고하고, 폐기/유지/보류 판단을 제안해줘.
완료 후 docs/CONTEXT.md 작업 로그에 기록
```

---

## Phase 1: CLAUDE.md + CONTEXT.md

### 1-1 시작 프롬프트

```
CLAUDE.md를 재작성해줘. harness_engineering_design.md의 섹션 3을 기반으로 하되,
다음 원칙을 따라:

1. 60줄 이내 — 백과사전이 아니라 짧은 지도
2. 절대 규칙 7개 (설계서 섹션 3.2 그대로)
3. 워크플로우 A/B Phase 번호 매핑표 + 디렉토리 매핑
4. Definition of Done 4개 조건
5. 스킬 목록 7개 (이름 + 한 줄 설명)
6. 독립 모듈: egfr_pipeline/md/
7. 참조 문서를 역할별로 그룹핑:
   - 구조: PIPELINE_ARCHITECTURE_REPORT.md, data_inventory.md
   - 실행: runbook.md
   - 과학: methodology_limitations.md, workflow_comparison_guide.md, phase4_A3_axis_specification.md
   - 추적: 투두리스트.md, docs/CONTEXT.md
8. 실험적 근거 3개 (ATP, Ko et al., 리간드 다양성) — 기존 내용 유지

기존 CLAUDE.md의 3 File System 관련 내용은 전부 제거.
완료 후 docs/CONTEXT.md 작업 로그에 기록
```

### 1-2 시작 프롬프트

```
docs/CONTEXT.md를 아래 구조로 재작성해줘.
기존 내용은 전부 교체한다.

---
# 프로젝트 컨텍스트

## 현재 작업 상태
- 워크플로우: (A / B / 양쪽)
- 현재 작업: 하네스 엔지니어링 구축 Phase 1 완료
- 다음 작업: Phase 2 스킬 구축

## 작업 로그
(Phase 0에서 기록한 내용 유지)
- [오늘 날짜] Phase 1-1: CLAUDE.md 재작성 완료
- [오늘 날짜] Phase 1-2: CONTEXT.md 재작성 완료

## 최근 결정 사항
- [오늘 날짜] 하네스 엔지니어링 적용 시작. 설계서: harness_engineering_design.md

## 발견된 이슈 (미해결)
(Phase 0-3에서 발견된 사항이 있으면 여기에)

## 실패 패턴 (반복 방지)
(아직 없음)
---

이 파일은 앞으로 모든 작업 세션에서 업데이트한다.
작업 시작 시 읽고, 작업 완료 시 작업 로그에 기록하는 것이 규칙이다.
```

---

## Phase 2: 스킬 구축

### 2-1 시작 프롬프트 (첫 스킬)

```
하네스 스킬 구축을 시작한다.
harness_engineering_design.md의 섹션 4를 기반으로 스킬 파일을 생성해줘.

먼저 .claude/skills/ 디렉토리 구조를 만들고,
가장 중요한 2개부터 시작:

1. .claude/skills/phase-dependencies/SKILL.md
   — 설계서 섹션 4.4 내용 기반
   — 끝에 "상세 참조" 섹션 추가:
     PIPELINE_ARCHITECTURE_REPORT.md, docs/data_inventory.md,
     docs/output_path_guide.md, docs/architecture.md

2. .claude/skills/bug-history/SKILL.md
   — 설계서 섹션 4.7 내용 기반
   — 끝에 "상세 참조" 섹션 추가:
     설계의도.md (PyRosetta 절대 주의사항 섹션)

완료 후 docs/CONTEXT.md 작업 로그에 기록
```

### 2-2 시작 프롬프트 (나머지 스킬)

```
나머지 5개 스킬을 생성해줘.
harness_engineering_design.md 섹션 4의 각 스킬 내용을 기반으로.

3. .claude/skills/ppi-analysis/SKILL.md (섹션 4.2)
   상세 참조: docs/manual_pyrosetta.md, docs/phase1_notes.md,
   docs/methodology_limitations.md

4. .claude/skills/vina-docking/SKILL.md (섹션 4.1)
   상세 참조: docs/manual_vina.md, docs/methodology_limitations.md

5. .claude/skills/hpc-operations/SKILL.md (섹션 4.3)
   상세 참조: docs/runbook.md, docs/environment_setup.md,
   docs/manual_execution.md, docs/pre_qsub_test_line.md, config/README.md

6. .claude/skills/scoring-system/SKILL.md (섹션 4.5)
   상세 참조: docs/phase4_A3_axis_specification.md,
   docs/workflow_comparison_guide.md, docs/workflow_comparison_design.md

7. .claude/skills/testing/SKILL.md (섹션 4.6)
   상세 참조: docs/test_suite_triage.md, docs/ac_coverage_checklist.md,
   docs/output_artifact_map.md

완료 후 docs/CONTEXT.md 작업 로그에 기록
```

---

## Phase 3: 에이전트 + 훅

### 3-1 시작 프롬프트

```
에이전트 정의와 훅을 생성해줘.
harness_engineering_design.md 섹션 5, 7을 기반으로.

에이전트:
1. .claude/agents/pipeline-dev.md (섹션 5.1)
2. .claude/agents/reviewer.md (섹션 5.2)
3. .claude/agents/science-qa.md (섹션 5.3)

훅:
4. .claude/hooks/pre-commit.sh (섹션 7.1)
5. .claude/hooks/csv-schema-guard.py (섹션 7.2)

에이전트 파일 앞에 섹션 5.0의 권한 경계 원칙
("에이전트 자율 가능" vs "반드시 사람 승인 필요")을 포함해줘.

완료 후 docs/CONTEXT.md 작업 로그에 기록
```

---

## Phase 4: 문서 체계화

### 4-1 시작 프롬프트

```
문서 체계화 마무리 작업이다.

1. README.md 문서 안내 테이블을 역할별로 그룹핑해서 재구성해줘:
   - 구조 이해 (PIPELINE_ARCHITECTURE_REPORT, PROJECT_USAGE_OVERVIEW, data_inventory)
   - 실행 (runbook, environment_setup, manual_vina, manual_pyrosetta)
   - 과학적 맥락 (methodology_limitations, workflow_comparison_guide, phase4_A3)
   - 테스트/검증 (test_suite_triage, pre_qsub_test_line)
   - 설정 (config/README.md)
   - AI 에이전트 (CLAUDE.md)
   폐기된 문서(nightly_review 등)는 이미 제거되었으므로 빠져 있어야 함

2. 프로젝트 전체에서 문서 drift 최종 점검:
   - 폐기된 시스템 참조가 남아있는지 grep
   - 각 스킬의 "상세 참조"에 나온 문서가 실제로 존재하는지 확인

완료 후 docs/CONTEXT.md 작업 로그에 기록하고,
현재 작업 상태를 "하네스 구축 완료"로 업데이트
```

---

## 작업 로깅 규칙

모든 Phase에 공통 적용. CLAUDE.md에 이 규칙을 포함한다.

```
docs/CONTEXT.md 업데이트 규칙:
- 세션 시작 시: "현재 작업 상태" 섹션을 읽고 맥락을 파악한다
- 작업 완료 시: "작업 로그"에 날짜 + 한 줄 요약을 추가한다
- 중요한 결정 시: "최근 결정 사항"에 이유와 함께 기록한다
- 이슈 발견 시: "발견된 이슈"에 추가한다
- 실패 시: "실패 패턴"에 유형 + 원인 + 해결 방법을 기록한다
```
