## Context Summary
- Project: EGFR–MYO1D structural analysis pipeline for standardized docking, pocket, and residue-level comparison across three receptor states
- Current Phase: Phase 1 (Full PRD)
- Status: Full PRD drafted from approved brief and PRD outline
- Phase Inputs: `brief-egfr-myo1d-pipeline.md`, `prd-outline-egfr-myo1d-pipeline.md`
- Key Constraints: Existing GitHub codebase must be reused; practical parallel ceiling is 16 usable CPU cores on a shared 32-core server; new computational outputs take precedence over legacy report residue/site labels; no low-level implementation details in this phase
- Output Rule: English first, Korean translation after the English PRD

---

# Product Requirements Document (PRD)
## EGFR–MYO1D Pipeline

## 1. Introduction & Overview

This project is a research-focused computational pipeline for the EGFR–MYO1D study. Its purpose is to standardize the way docking, pocket, and residue-level outputs are generated, compared, and reviewed across three receptor structural states:

1. 3GT8 raw structure
2. MD cluster representative from frames 38–48
3. MD cluster representative from frames 85–100

The current workflow depends on multiple disconnected scripts and a large amount of manual interpretation. This creates several recurring problems:

- raw docking outputs are difficult to compare across receptor states,
- pose-to-pocket assignments are not stored in a reusable format,
- residue-level contact summaries are not extracted consistently,
- supporting PPI outputs from PyRosetta and AlphaFold-Multimer are difficult to compare with ligand-based pocket outputs,
- and coding agents such as Codex cannot easily improve the repository without repeated scientific re-explanation.

The goal of this product is not to replace scientific judgment. Instead, it should convert a scattered research workflow into a reproducible evidence-producing system. The system must preserve raw metrics, standardize outputs, and support interpretation across receptor states without hard-coding old assumptions from previous reports.

This is an internal research MVP, not a public software product.

---

## 2. User Personas

### Persona 1: Primary Researcher
**Role:** Project owner and day-to-day user

**Profile:**
- Runs docking and structural analysis repeatedly
- Compares receptor states and ligand behaviors
- Needs residue-level outputs that can be reused across analysis sessions
- Wants results to be easy to hand off to coding agents or collaborators

**Goals:**
- Run the same analysis across multiple receptor states consistently
- Identify pockets and residue patches without manually re-parsing raw outputs
- Preserve enough detail to inspect, challenge, and reinterpret results later
- Improve or refactor the pipeline incrementally using repository-based tooling

**Pain Points:**
- Results are spread across scripts, folders, and ad hoc notes
- Different runs are hard to compare cleanly
- Output formats are not standardized
- Manual parsing is slow and error-prone
- Scientific context has to be repeatedly re-explained to new tooling sessions

### Persona 2: Collaborator / Reviewer
**Role:** Supervising professor, collaborator, or technically literate reviewer

**Profile:**
- Does not need to run every step manually
- Needs interpretable summaries of pockets, residues, and receptor-state differences
- May review the project intermittently rather than continuously

**Goals:**
- Understand what changed between receptor states
- Review which ligands map to which pockets
- See residue-level evidence without opening raw docking files
- Evaluate whether the analysis logic is scientifically defensible

**Pain Points:**
- Raw outputs are too fragmented
- Scientific conclusions are harder to audit when intermediate data is missing
- Legacy residue/site labels may bias interpretation unless new outputs are clearly separated

---

## 3. Feature Requirements

Each feature below is intended to be independently scorable, buildable, and testable in Phase 2.

---

### Feature 1: Structured Input and Run Management
**Priority:** Must-Have  
**Zone:** 🟢 Green

#### Description
The system must provide a consistent way to define receptor structures, ligands, shared run parameters, and project metadata so that all downstream results remain traceable and comparable.

This feature exists because the pipeline compares multiple receptor states and multiple ligands. If receptor identity, residue numbering, or run settings are not recorded consistently, then downstream pocket comparisons become unreliable.

The system should treat the following receptor states as first-class named inputs:
- 3GT8 raw
- MD cluster representative 38–48
- MD cluster representative 85–100

The system must preserve receptor identifiers, ligand identifiers, and run metadata in a form that later features can read directly.

#### User Story
As the primary researcher, I want to define receptors, ligands, and shared run settings once in a structured format, so that all later outputs can be traced back to the exact structural state and run context that produced them.

#### Acceptance Criteria
- [ ] The system accepts the three receptor states as explicitly named project inputs.
- [ ] The system accepts multiple ligands through a structured project-level input format.
- [ ] Receptor identity, ligand identity, and run metadata are preserved in downstream outputs.
- [ ] The system records enough metadata to distinguish results from different receptor states without ambiguity.
- [ ] The system prevents silent mixing of outputs between receptor states.
- [ ] The system supports later session recovery by storing project context in files rather than relying on memory alone.

#### Dependencies
None. This feature is the prerequisite for all later features.

---

### Feature 2: Parallel Batch Docking Execution
**Priority:** Must-Have  
**Zone:** 🟢 Green

#### Description
The system must execute receptor × ligand docking jobs as a configurable batch workflow and support practical parallel execution on the available shared server environment.

The working environment is a 32-core server where only 16 CPU cores are assumed to be safely available for this pipeline in routine use. The product must therefore support user-defined worker limits and allow the pipeline to run efficiently within that practical ceiling.

The goal of this feature is not simply speed for its own sake. The goal is to make repeated research runs practically usable while preserving traceability, deterministic output structure, and clear failure reporting.

#### User Story
As the primary researcher, I want the docking matrix to run in parallel using up to the available safe CPU budget, so that repeated multi-receptor multi-ligand analyses finish faster without becoming harder to inspect or debug.

#### Acceptance Criteria
- [ ] The system can execute docking jobs across the receptor × ligand matrix as a batch.
- [ ] The user can set a maximum worker count for parallel execution.
- [ ] The pipeline can operate within a practical 16-core ceiling on the shared server.
- [ ] Parallel execution does not change the logical structure of output files.
- [ ] Failed jobs are visible and traceable instead of silently disappearing.
- [ ] Results remain attributable to the exact receptor and ligand combination that generated them.
- [ ] Sequential and parallel runs produce logically equivalent result structures.

#### Dependencies
Depends on Feature 1 because batch execution requires structured input definitions.

---

### Feature 3: Pose Parsing, Contact Extraction, and Pocket Clustering
**Priority:** Must-Have  
**Zone:** 🟢 Green

#### Description
The system must convert raw Vina docking outputs into standardized pose-level data, extract receptor contact residues for each pose, compute pose centroids, and assign poses into reproducible pocket groupings using a user-defined centroid rule.

This feature is the central evidence-conversion layer of the product. Raw docking files are not sufficient for scientific comparison unless they are transformed into reusable intermediate data structures.

The system should preserve both interpreted and raw information. It should not only assign a pose to a pocket, but also keep the raw centroid coordinates, pose ranking, affinity values, and residue contact evidence needed for later review.

This feature must not assume that all pockets from all receptor states are directly comparable. It only standardizes pose-level and receptor-local pocket structure.

#### User Story
As the primary researcher, I want raw docking outputs converted into reusable pose-level and pocket-level datasets, so that I can inspect residue contacts and pocket assignments without manually re-parsing docking files every time.

#### Acceptance Criteria
- [ ] The system parses raw docking outputs into a standardized pose-level dataset.
- [ ] Each pose record includes enough information to identify the receptor, ligand, pose ranking, and core pose metrics.
- [ ] The system extracts receptor contact residues for each pose in a reproducible way.
- [ ] The system computes pose centroids for later pocket grouping.
- [ ] The system groups poses into receptor-local pockets using the configured centroid rule.
- [ ] The system preserves raw metrics needed for later manual inspection and reinterpretation.
- [ ] Pocket assignment is reproducible from the same inputs and settings.

#### Dependencies
Depends on Feature 2 because raw docking outputs must exist first.

---

### Feature 4: Pocket Summary and Cross-Receptor Comparison
**Priority:** Must-Have  
**Zone:** 🟡 Yellow

#### Description
The system must summarize pocket-level evidence for each receptor state and compare pockets across receptor states using location, residue overlap, and ligand assignment evidence.

This feature must help the user answer a key scientific question: whether pockets observed in different receptor states represent the same structural patch, a partial overlap, or distinct pockets. The system must support this judgment without pretending to settle it automatically.

Pocket comparison is therefore evidence-driven rather than conclusion-driven. The product should generate structured metrics that support manual interpretation, not force a single black-box verdict.

#### User Story
As the primary researcher or reviewer, I want pocket summaries and cross-receptor comparison metrics, so that I can judge whether pockets recur, shift, overlap partially, or disappear across receptor states.

#### Acceptance Criteria
- [ ] The system generates a pocket summary for each receptor state.
- [ ] Each pocket summary includes counts, representative information, and residue-level summary information.
- [ ] The system compares pockets across receptor states using raw comparison metrics such as location and residue-overlap evidence.
- [ ] The output makes it possible to distinguish same-patch candidates from partial overlaps and clearly separate pockets.
- [ ] The system does not hide the raw evidence behind a single opaque classification.
- [ ] The comparison output can be reviewed without opening raw docking files.

#### Dependencies
Depends on Feature 3 because cross-receptor comparison requires standardized pocket outputs.

---

### Feature 5: Integrated Research Reporting with Supporting PPI Outputs
**Priority:** Should-Have  
**Zone:** 🟡 Yellow

#### Description
The system must generate readable reports that combine Vina-derived pocket outputs with supporting residue-level outputs from PyRosetta global docking and AlphaFold-Multimer.

These supporting modules are not the primary pocket-definition engine. Their role is to provide auxiliary residue-side structural context that can be reviewed alongside ligand-derived pocket evidence.

The reports must support three use cases:
- self-review by the primary researcher,
- intermittent review by collaborators,
- and repository improvement by coding agents who need compact but faithful context.

#### User Story
As the primary researcher, I want the main pocket outputs and supporting PPI residue summaries compiled into readable reports, so that I can review the project, hand it off, and resume later without reconstructing the scientific context from scratch.

#### Acceptance Criteria
- [ ] The system generates readable summaries of Vina-derived pocket outputs.
- [ ] The system can include supporting PyRosetta residue-side summaries.
- [ ] The system can include supporting AlphaFold-Multimer residue-side summaries.
- [ ] Supporting PPI outputs are clearly labeled as auxiliary evidence, not hard-coded truth.
- [ ] Reports are usable by both the primary researcher and technically literate collaborators.
- [ ] Reports preserve enough context for a future coding session to continue work from files rather than from memory alone.

#### Dependencies
Depends primarily on Features 3 and 4, and secondarily on whatever standardized PyRosetta/AFM outputs are available.

---

## 4. User Flows

### Flow 1: Run a Standard Multi-Receptor Docking Batch
1. The researcher prepares or selects the three receptor states.
2. The researcher provides a ligand set and project-level run settings.
3. The system validates the structured inputs.
4. The system runs the receptor × ligand docking batch using the allowed parallel worker budget.
5. Raw outputs are stored with receptor and ligand identity preserved.
6. The researcher confirms that the run completed and that failed jobs, if any, are visible.

### Flow 2: Convert Raw Docking Outputs into Reviewable Pocket Data
1. The researcher starts from completed docking outputs.
2. The system parses pose-level results.
3. The system extracts contact residues for each pose.
4. The system computes centroids and assigns poses into receptor-local pockets.
5. The system generates pose-level and pocket-level structured outputs.
6. The researcher reviews pocket summaries without reopening raw docking files.

### Flow 3: Compare Pockets Across Receptor States
1. The researcher selects completed pocket outputs for the three receptor states.
2. The system calculates pocket comparison evidence across receptor states.
3. The system produces overlap/location evidence and comparison summaries.
4. The researcher reviews whether pockets appear recurrent, shifted, partially overlapping, or distinct.
5. The researcher uses the results to prioritize further manual interpretation.

### Flow 4: Review Supporting PPI Context
1. The researcher supplies or points to standardized PyRosetta and/or AFM outputs.
2. The system reads receptor-side residue summaries from those modules.
3. The system incorporates them into an integrated report alongside Vina-based pocket summaries.
4. The researcher or collaborator reviews the report as an interpretation aid.

---

## 5. Non-Goals (Out of Scope)

The MVP explicitly does **not** include the following:

- A web application or interactive dashboard
- Production deployment or cloud scaling architecture
- User authentication or account systems
- Billing, subscription, or payment handling
- Automated wet-lab integration workflows
- Fully automatic scientific decision-making that replaces human interpretation
- Treating old report residue/site labels as fixed ground truth
- Broad general-purpose bioinformatics platform behavior outside the EGFR–MYO1D use case

These items are excluded because the current objective is a focused research MVP for structured computational output generation and comparison.

---

## 6. Technical Considerations

This section records product-level technical constraints without specifying detailed implementation architecture.

### Required Constraints
- The existing GitHub codebase must be reused and improved rather than replaced from scratch.
- Python is the default implementation language unless explicitly changed later with approval.
- The pipeline must remain file-based and session-recoverable.
- Receptor identity and residue numbering consistency must be preserved wherever possible.
- The pipeline must support routine use on a shared server where only 16 CPU cores are assumed safely available for this work.
- Parallel execution must be configurable rather than fixed.
- New computational outputs must take precedence over older manually labeled report interpretations.
- Outputs must be readable by both humans and downstream tooling.

### Deferred to Later Phases
The following are intentionally deferred to later planning phases:
- exact project directory structure,
- exact CSV/JSON schema details,
- exact library selection,
- exact clustering algorithm choice,
- exact logging framework,
- exact retry/failure handling implementation,
- and exact report rendering format.

Those choices belong in the task and implementation phases, not in the PRD itself.

---

## 7. MVP Success Criteria

The MVP is successful when the following are all true:

- [ ] The pipeline accepts the three defined receptor states and a ligand set through a structured project input.
- [ ] The pipeline can run docking jobs in batch within a configurable parallel worker budget and operate practically within a 16-core shared-server ceiling.
- [ ] The pipeline produces reusable pose-level outputs rather than requiring raw docking files to be re-read manually.
- [ ] The pipeline produces receptor-local pocket summaries with residue-level evidence.
- [ ] The pipeline produces cross-receptor pocket comparison outputs that support manual interpretation.
- [ ] The outputs are sufficiently standardized that future sessions can resume from files and project docs rather than from memory alone.
- [ ] Supporting PyRosetta and AFM residue-side outputs can be attached as auxiliary evidence in integrated review documents.
- [ ] The resulting repository state is clear enough for iterative improvement by a coding agent such as Codex.

---

## 8. Open Questions

### Q1. How much of the current GitHub codebase is reusable without major restructuring?
**Impact:** High

Why it matters: This affects implementation scope, migration strategy, and the boundary between refactor and rewrite.

### Q2. Should pocket comparison include an automatic same-patch candidate heuristic in the MVP, or only raw comparison metrics?
**Impact:** High

Why it matters: This changes how much interpretation logic the system performs versus how much it leaves to the researcher.

### Q3. How much standardized detail is currently available from the existing PyRosetta and AFM scripts?
**Impact:** Medium

Why it matters: This determines whether Feature 5 is mostly report assembly or also requires substantial parser work.

### Q4. Will all receptor states preserve residue numbering cleanly enough for direct comparison, or is a normalization layer needed?
**Impact:** High

Why it matters: Residue mismatch could undermine pocket and residue comparison validity.

### Q5. What is the minimum report form that is useful enough for review without introducing unnecessary complexity?
**Impact:** Medium

Why it matters: This affects the boundary between a lightweight research report and a prematurely overbuilt reporting system.

---

# Korean Translation (한글 번역)

## 1. 소개 및 개요

이 프로젝트는 EGFR–MYO1D 연구를 위한 연구용 계산 파이프라인이다. 목적은 세 개의 receptor 구조 상태에 대해 생성되는 docking, pocket, residue 수준 결과를 표준화하고 비교 가능하게 만드는 것이다.

비교 대상 receptor 상태는 다음 세 가지다.

1. 3GT8 원본 구조
2. MD 38–48 프레임 구간 클러스터 대표 구조
3. MD 85–100 프레임 구간 클러스터 대표 구조

현재 워크플로우는 여러 개의 분리된 스크립트와 많은 수동 해석에 의존하고 있다. 그래서 다음과 같은 문제가 반복된다.

- raw docking output을 receptor 상태 간 비교하기 어렵다.
- pose-to-pocket assignment가 재사용 가능한 형식으로 저장되지 않는다.
- residue 수준 contact 요약이 일관되게 추출되지 않는다.
- PyRosetta와 AlphaFold-Multimer의 보조 PPI 결과를 ligand 기반 pocket 결과와 함께 보기가 어렵다.
- Codex 같은 코딩 에이전트가 저장소를 개선하려 해도 과학적 맥락을 반복해서 다시 설명해야 한다.

이 제품의 목표는 과학적 판단을 대체하는 것이 아니다. 대신 산발적인 연구 워크플로우를 재현 가능하고 증거를 생산하는 시스템으로 바꾸는 것이다. 시스템은 raw metric을 보존하고, 출력을 표준화하며, 이전 보고서의 가정을 하드코딩하지 않은 채 receptor 상태 간 해석을 지원해야 한다.

이것은 공개 소프트웨어 제품이 아니라 내부 연구용 MVP다.

---

## 2. 사용자 페르소나

### 페르소나 1: 주 연구자
**역할:** 프로젝트 소유자이자 일상 사용자

**프로필:**
- docking과 구조 분석을 반복 수행한다.
- receptor 상태와 ligand 거동을 비교한다.
- 여러 세션에 걸쳐 재사용 가능한 residue-level output이 필요하다.
- 결과를 코딩 에이전트나 협업자에게 쉽게 넘기고 싶다.

**목표:**
- 여러 receptor 상태에서 같은 분석을 일관되게 수행하고 싶다.
- raw output을 매번 손으로 다시 열지 않고 pocket과 residue patch를 보고 싶다.
- 나중에 다시 검토하고 재해석할 수 있을 만큼 충분한 디테일을 보존하고 싶다.
- 저장소 기반 도구를 사용해 파이프라인을 점진적으로 개선하고 싶다.

**불편점:**
- 결과가 스크립트, 폴더, 임시 메모에 흩어져 있다.
- 다른 실행 결과를 깔끔하게 비교하기 어렵다.
- 출력 형식이 표준화되어 있지 않다.
- 수동 파싱이 느리고 오류가 많다.
- 새 도구 세션마다 과학적 맥락을 반복 설명해야 한다.

### 페르소나 2: 협업자 / 검토자
**역할:** 지도교수, 공동연구자, 또는 기술 이해도가 있는 리뷰어

**프로필:**
- 모든 단계를 직접 실행할 필요는 없다.
- pocket, residue, receptor 상태 차이에 대한 해석 가능한 요약이 필요하다.
- 간헐적으로 프로젝트를 검토할 가능성이 높다.

**목표:**
- receptor 상태에 따라 무엇이 달라졌는지 이해하고 싶다.
- 어떤 ligand가 어떤 pocket에 들어가는지 보고 싶다.
- raw docking 파일을 열지 않고 residue 수준 증거를 보고 싶다.
- 분석 논리가 과학적으로 방어 가능한지 평가하고 싶다.

**불편점:**
- raw output이 너무 파편화되어 있다.
- 중간 데이터가 없으면 과학적 결론을 검토하기 어렵다.
- 기존 residue/site 라벨이 분리되어 있지 않으면 해석이 편향될 수 있다.

---

## 3. 기능 요구사항

각 기능은 Phase 2에서 독립적으로 점수화, 구현, 테스트될 수 있어야 한다.

### Feature 1: Structured Input and Run Management
**우선순위:** Must-Have  
**Zone:** 🟢 Green

#### 설명
시스템은 receptor 구조, ligand, 공통 실행 파라미터, 프로젝트 메타데이터를 일관되게 정의할 수 있어야 하며, 이를 통해 downstream 결과가 추적 가능하고 비교 가능해야 한다.

이 기능이 필요한 이유는 파이프라인이 여러 receptor 상태와 여러 ligand를 비교하기 때문이다. receptor identity, residue numbering, run setting이 일관되게 기록되지 않으면 이후 pocket 비교는 신뢰하기 어렵다.

시스템은 다음 receptor 상태를 1급 입력으로 다뤄야 한다.
- 3GT8 원본
- MD cluster representative 38–48
- MD cluster representative 85–100

#### 사용자 스토리
주 연구자로서, receptor, ligand, 공통 실행 설정을 한 번 구조화해서 정의하고 싶다. 그래야 이후 모든 output이 정확히 어떤 구조 상태와 실행 맥락에서 나왔는지 추적할 수 있다.

#### 수락 기준
- [ ] 시스템이 세 receptor 상태를 명시적으로 이름 붙은 프로젝트 입력으로 받는다.
- [ ] 시스템이 여러 ligand를 구조화된 프로젝트 입력 형식으로 받는다.
- [ ] receptor identity, ligand identity, run metadata가 downstream output에 보존된다.
- [ ] 결과가 어떤 receptor 상태에서 나왔는지 모호하지 않도록 충분한 메타데이터가 기록된다.
- [ ] receptor 상태 간 output이 조용히 섞이지 않도록 방지한다.
- [ ] 프로젝트 컨텍스트가 파일로 저장되어 다음 세션 복구를 지원한다.

#### 의존성
없음. 모든 후속 기능의 전제다.

---

### Feature 2: Parallel Batch Docking Execution
**우선순위:** Must-Have  
**Zone:** 🟢 Green

#### 설명
시스템은 receptor × ligand docking job을 batch workflow로 실행할 수 있어야 하며, 공유 서버 환경에서 실질적으로 사용 가능한 병렬 실행을 지원해야 한다.

작업 환경은 32코어 서버지만, 이 파이프라인에서는 통상 16코어만 안전하게 사용할 수 있다고 가정한다. 따라서 제품은 사용자 정의 worker 제한을 지원해야 하며, 이 실질적 상한선 안에서 효율적으로 실행되어야 한다.

이 기능의 목표는 단순 속도 향상이 아니라, 반복 연구 실행을 현실적으로 가능하게 하면서도 추적 가능성과 출력 구조, 실패 보고를 유지하는 것이다.

#### 사용자 스토리
주 연구자로서, 여러 receptor와 ligand 조합의 docking을 병렬로 돌리고 싶다. 그래야 반복 분석을 더 빨리 끝내면서도 결과 추적성과 디버깅 가능성을 유지할 수 있다.

#### 수락 기준
- [ ] 시스템이 receptor × ligand 조합 전체를 batch로 실행할 수 있다.
- [ ] 사용자가 병렬 worker 수 상한을 지정할 수 있다.
- [ ] 파이프라인이 공유 서버에서 실질적 16코어 범위 안에서 운영 가능하다.
- [ ] 병렬 실행이 결과 파일의 논리 구조를 바꾸지 않는다.
- [ ] 실패한 job이 숨겨지지 않고 추적 가능하다.
- [ ] 모든 결과가 정확한 receptor × ligand 조합에 귀속된다.
- [ ] 순차 실행과 병렬 실행의 결과 구조가 논리적으로 동일하다.

#### 의존성
Feature 1에 의존한다.

---

### Feature 3: Pose Parsing, Contact Extraction, and Pocket Clustering
**우선순위:** Must-Have  
**Zone:** 🟢 Green

#### 설명
시스템은 raw Vina docking output을 표준 pose-level 데이터로 바꾸고, 각 pose의 receptor contact residue를 추출하고, pose centroid를 계산하고, 사용자 정의 기준으로 재현 가능한 pocket grouping을 수행해야 한다.

이 기능은 제품의 핵심 증거 변환 계층이다. raw docking 파일은 재사용 가능한 중간 데이터 구조로 변환되지 않으면 과학적 비교에 충분하지 않다.

시스템은 해석된 정보와 raw 정보를 모두 보존해야 한다. 즉, pose를 pocket에 배정하는 것뿐 아니라 centroid 좌표, pose rank, affinity, residue contact 증거도 함께 남겨야 한다.

#### 사용자 스토리
주 연구자로서, raw docking output이 재사용 가능한 pose-level / pocket-level 데이터셋으로 자동 변환되길 원한다. 그래야 매번 raw 파일을 다시 파싱하지 않고 residue contact와 pocket assignment를 볼 수 있다.

#### 수락 기준
- [ ] 시스템이 raw docking output을 표준 pose-level 데이터셋으로 파싱한다.
- [ ] 각 pose에 receptor, ligand, pose ranking, 핵심 pose metric이 포함된다.
- [ ] 시스템이 각 pose의 receptor contact residue를 재현 가능하게 추출한다.
- [ ] 시스템이 centroid를 계산해 이후 pocket grouping에 사용한다.
- [ ] 시스템이 설정된 centroid 규칙으로 pose를 receptor-local pocket으로 grouping한다.
- [ ] 이후 수동 검토와 재해석에 필요한 raw metric을 보존한다.
- [ ] 같은 입력과 설정에서 pocket assignment가 재현된다.

#### 의존성
Feature 2에 의존한다.

---

### Feature 4: Pocket Summary and Cross-Receptor Comparison
**우선순위:** Must-Have  
**Zone:** 🟡 Yellow

#### 설명
시스템은 각 receptor 상태의 pocket-level 증거를 요약하고, 위치, residue overlap, ligand assignment를 근거로 receptor 상태 간 pocket 비교를 수행해야 한다.

이 기능은 중요한 과학 질문을 다룬다. 서로 다른 receptor 상태에서 보이는 pocket이 같은 구조 patch인지, 부분 overlap인지, 완전히 다른 pocket인지 평가할 수 있게 해야 한다. 하지만 시스템이 이를 자동으로 확정 판정하는 것은 아니다.

즉, pocket 비교는 결론 강요형이 아니라 **증거 제공형**이어야 한다.

#### 사용자 스토리
주 연구자나 검토자로서, pocket summary와 receptor 간 비교 metric을 보고 싶다. 그래야 pocket이 반복되는지, 이동하는지, 부분 overlap인지, 사라지는지를 판단할 수 있다.

#### 수락 기준
- [ ] 시스템이 각 receptor 상태에 대한 pocket summary를 생성한다.
- [ ] 각 summary에는 count, representative 정보, residue-level summary가 포함된다.
- [ ] 시스템이 receptor 상태 간 pocket을 raw 비교 metric으로 비교한다.
- [ ] 출력만으로 same-patch candidate, partial overlap, distinct pocket을 구분할 근거를 제공한다.
- [ ] 시스템이 단일 불투명 분류 뒤에 raw evidence를 숨기지 않는다.
- [ ] raw docking 파일을 열지 않고도 비교 결과를 검토할 수 있다.

#### 의존성
Feature 3에 의존한다.

---

### Feature 5: Integrated Research Reporting with Supporting PPI Outputs
**우선순위:** Should-Have  
**Zone:** 🟡 Yellow

#### 설명
시스템은 Vina 기반 pocket output과 함께 PyRosetta global docking, AlphaFold-Multimer의 residue-level 보조 output을 읽기 쉬운 report로 정리해야 한다.

이 보조 모듈들은 pocket 정의의 주된 엔진이 아니다. 역할은 ligand-derived pocket evidence와 나란히 볼 수 있는 receptor-side structural context를 제공하는 것이다.

이 report는 세 가지 용도를 지원해야 한다.
- 주 연구자의 자기 검토
- 협업자의 간헐적 리뷰
- 코딩 에이전트가 저장소를 개선할 수 있도록 하는 컨텍스트 전달

#### 사용자 스토리
주 연구자로서, pocket 결과와 보조 PPI residue summary가 한 문서로 정리되길 원한다. 그래야 프로젝트를 검토하고, 인계하고, 나중에 다시 시작할 때 과학적 맥락을 다시 처음부터 만들지 않아도 된다.

#### 수락 기준
- [ ] 시스템이 Vina 기반 pocket output을 읽기 쉬운 형태로 요약한다.
- [ ] 시스템이 PyRosetta residue-side summary를 보조 정보로 포함할 수 있다.
- [ ] 시스템이 AFM residue-side summary를 보조 정보로 포함할 수 있다.
- [ ] 보조 PPI output이 고정 진실이 아니라 auxiliary evidence로 명확히 표기된다.
- [ ] report가 연구자와 기술 이해도가 있는 협업자 모두에게 유용하다.
- [ ] future coding session이 memory가 아니라 file context로 이어질 수 있을 만큼 충분한 맥락을 담는다.

#### 의존성
주로 Features 3, 4에 의존하고, 부차적으로 표준화된 PyRosetta/AFM output 가용성에 의존한다.

---

## 4. 사용자 흐름

### Flow 1: 표준 multi-receptor docking batch 실행
1. 연구자가 세 receptor 상태를 준비하거나 선택한다.
2. ligand 세트와 프로젝트 수준 실행 설정을 제공한다.
3. 시스템이 구조화된 입력을 검증한다.
4. 시스템이 허용된 병렬 worker 범위 안에서 receptor × ligand batch를 실행한다.
5. raw output이 receptor와 ligand identity를 보존한 채 저장된다.
6. 연구자는 실행 완료 여부와 실패 job 유무를 확인한다.

### Flow 2: raw docking output을 검토 가능한 pocket 데이터로 변환
1. 연구자가 완료된 docking output을 준비한다.
2. 시스템이 pose-level 결과를 파싱한다.
3. 시스템이 각 pose의 contact residue를 추출한다.
4. 시스템이 centroid를 계산하고 receptor-local pocket에 배정한다.
5. 시스템이 pose-level 및 pocket-level 구조화 output을 생성한다.
6. 연구자는 raw docking 파일을 다시 열지 않고 pocket summary를 검토한다.

### Flow 3: receptor 상태 간 pocket 비교
1. 연구자가 세 receptor 상태의 pocket output을 준비한다.
2. 시스템이 receptor 간 pocket 비교 근거를 계산한다.
3. 시스템이 overlap/location evidence와 비교 summary를 생성한다.
4. 연구자는 pocket이 반복, 이동, 부분 overlap, 독립 pocket인지 검토한다.
5. 연구자는 결과를 바탕으로 후속 수동 해석 우선순위를 정한다.

### Flow 4: 보조 PPI 컨텍스트 검토
1. 연구자가 표준화된 PyRosetta 및/또는 AFM output을 제공한다.
2. 시스템이 receptor-side residue summary를 읽는다.
3. 시스템이 이를 Vina 기반 pocket summary와 함께 integrated report에 포함한다.
4. 연구자나 협업자가 이를 해석 보조 자료로 검토한다.

---

## 5. 비목표 (Out of Scope)

이번 MVP는 다음을 포함하지 않는다.

- 웹 애플리케이션 또는 인터랙티브 대시보드
- 프로덕션 배포나 클라우드 확장 아키텍처
- 사용자 인증/계정 시스템
- 결제/구독/과금 처리
- wet-lab 자동 통합 워크플로우
- 인간 해석을 완전히 대체하는 자동 과학 판단
- 기존 보고서 residue/site 라벨을 고정 진실처럼 취급하는 기능
- EGFR–MYO1D 범위를 벗어난 범용 바이오인포매틱스 플랫폼 역할

이는 현재 목표가 구조화된 계산 output 생성과 비교에 집중된 연구용 MVP이기 때문이다.

---

## 6. 기술적 고려사항

이 섹션은 제품 수준 제약만 기록하며, 상세 구현 구조는 아직 확정하지 않는다.

### 필수 제약
- 기존 GitHub 코드베이스를 버리지 않고 재사용해야 한다.
- Python을 기본 구현 언어로 사용한다.
- 파이프라인은 파일 기반이며 세션 복구 가능해야 한다.
- receptor identity와 residue numbering 일관성을 가능한 한 유지해야 한다.
- 공유 서버에서 실사용 가능 16코어 기준 병렬 실행을 지원해야 한다.
- 병렬 실행은 고정값이 아니라 설정 가능해야 한다.
- 새 계산 output이 기존 수동 라벨보다 더 높은 우선순위를 가져야 한다.
- 출력은 사람과 downstream tooling 모두가 읽을 수 있어야 한다.

### 후속 단계로 미루는 항목
다음 항목은 일부러 이후 단계로 미룬다.
- 정확한 프로젝트 디렉터리 구조
- 정확한 CSV/JSON schema
- 정확한 라이브러리 선택
- 정확한 clustering 알고리즘
- 정확한 logging framework
- 정확한 retry/failure handling 구현
- 정확한 report rendering format

이 선택들은 PRD가 아니라 task/implementation 단계의 범위다.

---

## 7. MVP 성공 기준

다음이 모두 충족되면 MVP가 성공한 것으로 본다.

- [ ] 파이프라인이 세 receptor 상태와 ligand 세트를 구조화된 프로젝트 입력으로 받을 수 있다.
- [ ] 파이프라인이 설정 가능한 worker budget 안에서 batch docking을 수행하고 실질적 16코어 범위 안에서 운영 가능하다.
- [ ] 파이프라인이 raw docking 파일을 사람이 다시 읽지 않아도 되는 pose-level output을 생성한다.
- [ ] 파이프라인이 receptor-local pocket summary와 residue-level evidence를 생성한다.
- [ ] 파이프라인이 receptor 상태 간 pocket 비교 output을 생성해 수동 해석을 지원한다.
- [ ] 출력이 충분히 표준화되어 향후 세션이 memory가 아니라 file/document로 이어질 수 있다.
- [ ] PyRosetta와 AFM residue-side output을 auxiliary evidence로 붙일 수 있다.
- [ ] 결과 저장소 상태가 Codex 같은 코딩 에이전트가 반복 개선할 수 있을 정도로 명확하다.

---

## 8. 열린 질문

### Q1. 현재 GitHub 코드베이스 중 얼마나 많은 부분을 큰 구조 변경 없이 재사용할 수 있는가?
**영향도:** 높음

중요한 이유: 리팩터링 범위와 migration 전략, 그리고 refactor와 rewrite의 경계를 결정한다.

### Q2. pocket 비교에서 MVP 단계부터 same-patch candidate heuristic을 자동 제공할 것인가, 아니면 raw metric만 제공할 것인가?
**영향도:** 높음

중요한 이유: 시스템이 어느 정도까지 해석을 수행할지 결정한다.

### Q3. 기존 PyRosetta/AFM 스크립트에서 얼마나 표준화된 residue-level detail을 바로 뽑을 수 있는가?
**영향도:** 중간

중요한 이유: Feature 5가 단순 report 조립인지, 추가 parser 개발까지 필요한지 결정된다.

### Q4. 세 receptor 상태에서 residue numbering이 direct comparison이 가능할 만큼 잘 유지되는가, 아니면 normalization layer가 필요한가?
**영향도:** 높음

중요한 이유: numbering mismatch는 pocket/residue 비교 타당성을 무너뜨릴 수 있다.

### Q5. 리뷰에 충분히 유용하면서도 과도하게 무겁지 않은 최소 report 형태는 무엇인가?
**영향도:** 중간

중요한 이유: lightweight research report와 과도한 reporting system 사이 경계를 정해야 한다.

