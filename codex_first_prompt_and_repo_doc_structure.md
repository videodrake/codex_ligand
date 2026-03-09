## Context Summary

- Project: EGFR–MYO1D structural analysis pipeline for standardized docking, pocket, and residue-level comparison across three receptor states
- Current Use: Immediate handoff package for starting Codex-assisted implementation
- Document Purpose: Provide (1) the first Codex execution prompt and (2) the recommended GitHub repository documentation structure
- Key Constraints: Reuse existing GitHub codebase; practical CPU ceiling is 16 usable cores on a shared 32-core server; start with Task Groups 0–2 only; preserve traceability and avoid full rewrites
- Primary Audience: Project owner, future collaborators, and Codex working inside the repository

---

# Codex 첫 실행 프롬프트 + GitHub 문서 구조 가이드

이 문서는 두 가지 목적을 가진다.

1. **Codex에게 처음 넘길 실행 프롬프트를 제공하는 것**
2. **GitHub 저장소 안에서 문서들을 어떻게 배치해야 하는지 가이드를 제공하는 것**

이 문서는 실제 구현 전에 Codex가 길을 잃지 않도록 하기 위한 **실행 준비 문서**다. 핵심 원칙은 다음과 같다.

- 기존 GitHub 코드베이스를 버리지 않고 재사용한다.
- 현재 연구의 핵심 우선순위는 **Vina 중심 표준화**다.
- receptor는 정확히 세 개의 상태를 비교한다.
  1. 3GT8 raw
  2. MD cluster representative 38–48
  3. MD cluster representative 85–100
- 서버는 32코어지만, 이 프로젝트에서는 **실사용 가능 16코어**를 기준으로 병렬 실행을 설계한다.
- 기존 보고서의 residue/site 해석은 **후순위 참고자료**이며, 새 계산 결과가 더 높은 우선순위를 가진다.
- 처음부터 전체 파이프라인을 다 고치지 않고, **Task Group 0\~2**만 먼저 시작한다.

---

# 1. Codex 첫 실행 프롬프트

아래 프롬프트는 Codex에게 그대로 전달할 수 있다. 이 프롬프트의 목적은 저장소 구조를 먼저 이해하게 하고, 가장 낮은 리스크의 리팩터링부터 시작하게 만드는 것이다.

## Codex Prompt

```text
You are working inside a GitHub repository for an internal research pipeline project.

Before changing code, read these documents in this order:

1. docs/brief-egfr-myo1d-pipeline.md
2. docs/prd-egfr-myo1d-pipeline.md
3. docs/tasks-egfr-myo1d-pipeline.md
4. CODEX_HANDOFF_EGFR_MYO1D_PIPELINE.md

This repository is for an EGFR–MYO1D structural analysis pipeline.
The main near-term goal is to standardize the Vina-centered workflow first.

Important project rules:
- Reuse the existing GitHub codebase instead of rewriting everything from scratch.
- The project compares exactly three receptor states:
  1) 3GT8 raw
  2) MD cluster representative 38–48
  3) MD cluster representative 85–100
- The practical server constraint is: 32 cores exist, but only 16 cores should be treated as safely usable for this project.
- Parallel execution must therefore be configurable, with 16 as the practical normal upper bound.
- New computational outputs should be treated as higher-priority evidence than old legacy residue/site labels.
- Do not hard-code old site names from earlier reports as truth.
- Keep the workflow file-based, resumable, and easy to hand off later.

Your first assignment is NOT to implement the whole pipeline.
Your first assignment is to work only on:
- Task Group 0: Project Setup and Repository Baseline
- Task Group 1: Structured Input and Run Management
- Task Group 2: Parallel Batch Docking Execution

What I want from you in this first pass:

1. Inspect the current repository and identify:
   - the current Vina-related scripts
   - current config handling, if any
   - current output folders and naming logic
   - likely entry points for batch execution
   - obvious dead code, duplicate code, or risky code paths

2. Propose a minimal refactor plan for Task Groups 0–2 only.
   - Do not propose a full rewrite.
   - Prefer wrapping and refactoring existing code.
   - Show which files should be modified, which new files should be added, and why.

3. Implement only the smallest safe first step needed to support:
   - structured receptor/ligand/config input
   - configurable max_workers for batch docking
   - stable receptor/ligand-specific output placement
   - logging that makes failed jobs visible

4. Keep output structure consistent between sequential and parallel execution.

5. After you inspect the repo, give me:
   - a repository map
   - a file-by-file change plan
   - any assumptions you need me to confirm
   - then begin with the lowest-risk implementation step

Constraints:
- Do not start Task Group 3 or later yet.
- Do not invent web UI, database layers, authentication, or deployment infrastructure.
- Do not over-automate scientific interpretation.
- Do not collapse raw evidence into hidden logic.
- Keep the system Codex-friendly and document any new conventions clearly.

If there is ambiguity, prefer making the smallest safe refactor that preserves traceability and future extension.
```

---

## 추가 제어 문장 (선택 사항)

Codex가 너무 많은 파일을 한 번에 바꾸려 하거나 대규모 재작성으로 흐를 가능성이 있다면, 아래 문장을 프롬프트 마지막에 추가한다.

```text
Do not modify more than the minimum number of files needed in the first pass.
Prefer a narrow, inspectable patch over a broad architectural rewrite.
```

이 문장은 첫 구현 라운드를 더 안전하게 만들기 위한 제어문이다.

---

# 2. 이 프롬프트를 이렇게 설계한 이유

이 프롬프트는 단순한 기능 요청이 아니라, Codex가 잘못된 방향으로 들어가지 않도록 경계를 미리 설정하는 문서다.

## 2.1 문서를 먼저 읽게 하는 이유

Codex가 코드만 보고 수정에 들어가면, 프로젝트의 실제 우선순위를 오해할 가능성이 높다. 현재 프로젝트는 일반 제품 개발이 아니라, **연구용 계산 파이프라인 표준화**가 목적이다. 따라서 아래 문서를 먼저 읽게 해야 한다.

- Brief: 프로젝트 범위와 목표
- PRD: 기능 요구사항과 성공 기준
- Tasks: 구현 순서와 분할 단위
- Handoff 문서: 저장소와 분석 흐름의 구체적 맥락

## 2.2 Task Group 0\~2만 먼저 보게 하는 이유

처음부터 parser, pocket comparison, report, PyRosetta, AFM까지 다 건드리면 저장소를 이해하기 전에 범위가 터질 가능성이 높다. 따라서 첫 라운드에서는 아래까지만 다루게 한다.

- 저장소 점검 및 baseline
- 입력 구조화
- 16코어 병렬 batch 실행

이 세 단계가 안정화되어야 이후 parsing과 pocket summary가 안전해진다.

## 2.3 full rewrite를 금지하는 이유

이 저장소는 이미 스크립트가 있고, 이를 기반으로 점진적 개선을 해야 한다. 따라서 “처음부터 새로 만들기”가 아니라 **existing code wrapping + targeted refactoring**이 맞다.

## 2.4 16코어 제약을 먼저 못 박는 이유

공유 서버 환경에서는 병렬화 자체보다 **안전한 병렬화**가 중요하다. 따라서 Codex가 무제한 worker나 공격적인 병렬화를 제안하지 않도록, “32코어 서버지만 이 프로젝트는 16코어만 실사용 가능”이라는 제약을 명시해야 한다.

---

# 3. Codex가 첫 응답에서 ideally 내놓아야 할 내용

첫 번째 Codex 응답은 곧바로 대규모 코드 수정이 아니라, 아래 내용을 포함하는 것이 이상적이다.

1. **Repository map**

   - 현재 어떤 폴더와 스크립트가 있는지
   - Vina 실행 entry point가 어디인지
   - config 파일이나 유사 구조가 있는지
   - output 저장 방식이 어떤지

2. **File-by-file change plan**

   - 어떤 기존 파일을 수정할지
   - 어떤 새 파일을 만들지
   - 왜 그 파일이 필요한지

3. **Assumptions to confirm**

   - receptor 입력 형식
   - ligand 입력 형식
   - output 디렉터리 규칙
   - 기존 코드 중 재사용 가능한 부분

4. **Smallest safe first patch**

   - 입력 구조화 또는 max\_workers 설정 추가
   - output path isolation
   - 실패 로그 노출

즉, 첫 응답은 “전부 구현 완료”가 아니라 **저장소 이해 + 최소 위험 첫 수정안**이어야 한다.

---

# 4. GitHub 저장소에 넣을 문서 구조

아래 구조를 권장한다.

```text
repo-root/
├── README.md
├── CLAUDE.md
├── CODEX_HANDOFF_EGFR_MYO1D_PIPELINE.md
├── docs/
│   ├── brief-egfr-myo1d-pipeline.md
│   ├── prd-outline-egfr-myo1d-pipeline.md
│   ├── prd-egfr-myo1d-pipeline.md
│   ├── tasks-outline-egfr-myo1d-pipeline.md
│   ├── tasks-egfr-myo1d-pipeline.md
│   ├── project-context.md
│   ├── repository-map.md
│   └── runbook.md
├── config/
│   └── example-project.yaml
├── scripts/
├── receptors/
├── ligands/
├── outputs/
│   ├── raw/
│   ├── parsed/
│   ├── reports/
│   └── logs/
└── tests/
```

이 구조의 핵심은 다음과 같다.

- **repo root**는 입구 역할
- \*\*docs/\*\*는 설계/운영 문서 모음
- \*\*config/\*\*는 실행 설정 예시
- \*\*outputs/\*\*는 raw/parsed/report/log를 분리
- \*\*tests/\*\*는 회귀 방지 및 validation용

---

# 5. 각 문서의 역할

## 5.1 `README.md`

저장소에 처음 들어온 사람이 가장 먼저 읽는 문서다.

### 역할

- 프로젝트 한 줄 설명
- 이 저장소가 무엇을 위한 것인지 설명
- 현재 최우선이 Vina 중심 표준화라는 점 설명
- 어떤 문서를 어떤 순서로 읽어야 하는지 안내

### 반드시 포함할 내용

- 프로젝트 개요
- 빠른 시작 경로
- 문서 읽기 순서
- 현재 범위가 연구용 MVP라는 점

### 추천 첫 문장 예시

```md
# EGFR–MYO1D Pipeline

Research pipeline for standardized docking, pocket clustering, residue extraction, and cross-receptor comparison across three EGFR receptor states.

Start here:
1. docs/brief-egfr-myo1d-pipeline.md
2. docs/prd-egfr-myo1d-pipeline.md
3. docs/tasks-egfr-myo1d-pipeline.md
4. CODEX_HANDOFF_EGFR_MYO1D_PIPELINE.md
```

---

## 5.2 `CLAUDE.md`

이 문서는 프로젝트 운영 규칙, 단계 규칙, 3-file system 같은 **프로세스 규칙**을 담는다.

### 역할

- 문서 생성 순서
- 승인 방식
- 단계 건너뛰기 금지 규칙
- 문서 작성 규칙

즉, “무엇을 만들 것인가”보다 “어떤 순서로 만들 것인가”를 정의하는 문서다.

---

## 5.3 `CODEX_HANDOFF_EGFR_MYO1D_PIPELINE.md`

이 문서는 개발자/코딩 에이전트용 **핵심 인수인계 문서**다.

### 역할

- 연구 목적 설명
- receptor 3개 정의
- Vina/PyRosetta/AFM의 역할 정리
- output schema 개요
- Vina 우선 개발 원칙
- 16코어 병렬 실행 제약
- 구현 우선순위와 acceptance criteria

즉, Codex가 “왜 이 저장소를 바꾸는지” 이해하는 문서다.

---

## 5.4 `docs/brief-egfr-myo1d-pipeline.md`

Phase 0 Project Brief다.

### 역할

- 프로젝트 한 줄 정의
- 대상 사용자 정의
- 핵심 기능 5개 정의
- 비범위 정의
- 성공 기준 정의

이 문서는 프로젝트 범위를 짧고 단단하게 고정하는 역할을 한다.

---

## 5.5 `docs/prd-outline-egfr-myo1d-pipeline.md`

PRD 작성 전 outline이다.

### 역할

- full PRD 작성 전 섹션 구조를 보존
- 문서 설계 히스토리 유지

필수는 아니지만, 3-file system을 충실히 따를 경우 보존하는 것이 좋다.

---

## 5.6 `docs/prd-egfr-myo1d-pipeline.md`

이 문서는 **핵심 제품 요구 문서**다.

### 역할

- 기능 정의
- 사용자 스토리
- acceptance criteria
- non-goals
- technical considerations
- open questions

Codex가 실제 구현을 시작하기 전에 반드시 읽어야 한다.

---

## 5.7 `docs/tasks-outline-egfr-myo1d-pipeline.md`

Task 구조 초안이다.

### 역할

- Task Group 순서 정의
- 왜 그 순서인지 설명
- full task list 작성 전 구조 보여주기

---

## 5.8 `docs/tasks-egfr-myo1d-pipeline.md`

이 문서는 실제 구현용 핵심 문서다.

### 역할

- Task Group 0\~8 정의
- objective
- subtasks
- test tasks
- dependencies
- deliverables

즉, Codex에게 가장 자주 참조하게 될 구현 문서다.

---

## 5.9 `docs/project-context.md`

짧은 운영 맥락 메모다.

### 역할

- receptor 3개가 무엇인지
- 왜 이 3개를 쓰는지
- 기존 보고서는 후순위 참고라는 점
- 새 결과 우선 원칙
- wet 실험은 현재 범위 밖이라는 점

이 문서는 repo만 보고도 맥락을 유지하게 해준다.

---

## 5.10 `docs/repository-map.md`

저장소 구조 안내 문서다.

### 역할

- 어떤 폴더에 어떤 코드가 있는지
- Vina entry point 위치
- parser 위치
- report generator 위치
- legacy script 위치
- 조심해서 수정해야 할 파일

이 문서는 Codex 첫 응답을 바탕으로 생성해도 좋다.

---

## 5.11 `docs/runbook.md`

실행 설명서다.

### 역할

- 설정 파일 예시
- 배치 실행 예시
- 병렬 worker 설정법
- output 생성 위치
- 로그 위치
- 오류 발생 시 어디를 봐야 하는지

README보다 더 실무적인 실행 문서다.

---

# 6. 문서 배치 전략

가장 좋은 방식은 다음과 같다.

## repo root에 둘 문서

- `README.md`
- `CLAUDE.md`
- `CODEX_HANDOFF_EGFR_MYO1D_PIPELINE.md`

이 문서들은 저장소 진입점 역할을 한다.

## `docs/`에 둘 문서

- `brief-egfr-myo1d-pipeline.md`
- `prd-outline-egfr-myo1d-pipeline.md`
- `prd-egfr-myo1d-pipeline.md`
- `tasks-outline-egfr-myo1d-pipeline.md`
- `tasks-egfr-myo1d-pipeline.md`
- `project-context.md`
- `repository-map.md`
- `runbook.md`

이 문서들은 설계와 운영의 장기 기록 역할을 한다.

---

# 7. 지금 당장 repo에 반영할 권장 순서

## 1단계

repo root에 다음 문서를 둔다.

- `README.md`
- `CLAUDE.md`
- `CODEX_HANDOFF_EGFR_MYO1D_PIPELINE.md`

## 2단계

`docs/` 폴더를 만들고 다음 문서를 넣는다.

- `brief-egfr-myo1d-pipeline.md`
- `prd-outline-egfr-myo1d-pipeline.md`
- `prd-egfr-myo1d-pipeline.md`
- `tasks-outline-egfr-myo1d-pipeline.md`
- `tasks-egfr-myo1d-pipeline.md`

## 3단계

Codex에게는 처음에 전체 구현이 아니라, **Task Group 0\~2만** 시작하라고 한다.

즉, 초기 구현 범위는 다음으로 제한한다.

- 저장소 점검 및 baseline
- 입력 구조화
- 16코어 병렬 batch docking 실행

---

# 8. 최종 요약

이 문서의 핵심은 단순하다.

1. **Codex에게는 처음부터 전부 맡기지 않는다.**
2. **먼저 문서를 읽게 하고, Task Group 0\~2만 보게 한다.**
3. **repo root와 docs 구조를 명확히 나눠서 문서를 배치한다.**
4. **16코어 병렬 실행은 정식 요구사항으로 다룬다.**
5. **새 계산 결과가 기존 보고서보다 더 높은 해석 우선순위를 가진다.**

가장 좋은 시작 문장은 이걸로 요약된다.

> GitHub repo 루트에 handoff 문서를 두고, docs 폴더에 brief/PRD/tasks 문서를 정리한 뒤, Codex에게 Task Group 0\~2만 먼저 진행하라고 지시한다.

---

# Korean Translation (간단 요약)

이 문서는 두 가지를 위한 것이다.

- Codex에게 처음 줄 실행 프롬프트 제공
- GitHub 저장소 안에 문서를 어떻게 놓을지 가이드 제공

핵심 원칙은 다음과 같다.

- 기존 코드 재사용
- Vina 중심 우선 개발
- receptor 3개 고정
- 16코어 병렬 실행 지원
- 기존 보고서는 참고, 새 계산 결과 우선
- 첫 구현은 Task Group 0\~2까지만

추천 repo 구조는 다음과 같다.

- 루트: `README.md`, `CLAUDE.md`, `CODEX_HANDOFF_EGFR_MYO1D_PIPELINE.md`
- docs: brief, PRD, tasks, project context, repository map, runbook

Codex 첫 프롬프트는 문서를 먼저 읽고,

- 저장소 구조 파악
- Vina 관련 코드 파악
- config/output 구조 파악
- 최소 리팩터링 계획 제안
- 가장 낮은 리스크 첫 수정부터 시작 하도록 설계해야 한다.

