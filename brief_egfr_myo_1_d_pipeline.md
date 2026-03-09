## Context Summary
- Project: EGFR–MYO1D structural analysis pipeline for standardized docking, pocket, and residue-level comparison across receptor states
- Current Phase: Phase 0 (Project Brief)
- Status: Brief drafted and ready for user review before Phase 1 PRD generation
- Key Constraints: Use existing GitHub codebase; support up to 16 usable CPU cores for parallel execution; prioritize new computational outputs over legacy residue/site interpretations
- Tech Stack: Python (default, detailed stack deferred to Phase 2)

---

# Project Brief

## One-Line Description
A research-focused computational pipeline that standardizes Vina, PyRosetta, and AlphaFold-Multimer outputs into comparable pose-, pocket-, and residue-level datasets for EGFR–MYO1D analysis across three receptor structural states.

## Target User
**Primary user:** the researcher maintaining and extending the EGFR–MYO1D project.

**Secondary users:** a small group of collaborators, such as a supervising professor or research partner, who need to inspect results, compare receptor-state-dependent pockets, and review residue-level summaries without manually reconstructing raw outputs.

### Why this user needs the tool
The current workflow relies on multiple separate scripts and manually interpreted outputs. That makes it hard to:
- compare results across receptor states,
- track which ligands converge on which pockets,
- extract residue-level contact summaries consistently,
- reuse outputs in later analysis sessions,
- and hand off the project cleanly to coding agents such as Codex.

This project exists to turn a scattered computational workflow into a repeatable, inspectable, and codex-friendly research pipeline.

## Core Features (MVP Scope) — MAX 5

1. **Structured Input and Run Management**
   - The system must accept three receptor structures (3GT8 raw, MD cluster representative 38–48, MD cluster representative 85–100), multiple ligand inputs, and shared run settings in a consistent project-level configuration format.
   - The system must preserve receptor identity, residue numbering consistency, and run metadata so downstream comparisons remain valid.

2. **Parallel Batch Docking Execution for Available CPU Capacity**
   - The system must run docking jobs across the receptor × ligand matrix using configurable parallel execution.
   - The practical operating target is a server with 32 CPU cores where only 16 cores are considered safely available for this project.
   - The user must be able to define a maximum worker count, with 16 as the expected routine upper operating limit.
   - Batch execution must improve throughput without changing output structure or making results harder to trace.

3. **Pose Parsing, Contact Extraction, and Pocket Clustering**
   - The system must parse Vina outputs into a standardized pose table, extract receptor contact residues for each pose, compute pose centroids, and cluster poses into pockets using the user-defined centroid rule.
   - Pocket assignment must remain reproducible and must preserve raw metrics needed for later review.
   - The output must be easy to inspect both manually and programmatically.

4. **Pocket Summary and Cross-Receptor Comparison**
   - The system must summarize pockets for each receptor state, including pose counts, ligand counts, representative poses, and residue frequency summaries.
   - It must also compare pockets across receptor states using location and residue-overlap metrics, so the user can judge whether pockets are likely to represent the same patch, a partial overlap, or a distinct pocket.
   - The goal is not to force automatic conclusions, but to provide structured evidence for manual scientific interpretation.

5. **Integrated Research Reporting with Supporting PPI Outputs**
   - The system must generate readable reports that combine Vina-derived pocket summaries with supporting residue-level outputs from PyRosetta global docking and AlphaFold-Multimer.
   - These supporting modules are not the primary source of truth for pocket definition, but they should provide auxiliary receptor-side residue information that can be reviewed alongside docking-derived pocket data.
   - The final outputs must support review by the researcher, collaborators, and coding agents working on the repository.

> ⚠️ Each core feature will be treated as an independent unit in Phase 2 (Task Breakdown).
> They will be planned, built, and tested separately to avoid hidden coupling.

## Explicitly Out of Scope
- **Web application or graphical dashboard** — Reason: the MVP is a research pipeline, not a user-facing product.
- **Authentication, user accounts, or permissions** — Reason: this is not a multi-tenant platform.
- **Payments, billing, or subscription systems** — Reason: completely unrelated to the MVP.
- **Production deployment, CI/CD, cloud infrastructure, or scaling architecture** — Reason: the immediate goal is a working research MVP on an existing server.
- **Automated wet-lab data integration** — Reason: current scope is computational pipeline standardization only.
- **Overwriting scientific judgment with fully automatic residue conclusions** — Reason: the system should structure evidence, not pretend to replace manual interpretation.
- **Treating legacy report residues/sites as fixed truth** — Reason: newly generated computational outputs should carry higher interpretive weight than legacy labels.

## Success Criteria
The MVP is considered successful when all of the following are true:

- [ ] The researcher can run the pipeline on the three defined receptor states using a shared configuration file.
- [ ] The system can execute docking jobs in parallel with a configurable worker limit and operate cleanly within a practical 16-core usage ceiling.
- [ ] The system produces standardized outputs for Vina results, including a pose-level table, pocket-level table, ligand-to-pocket mapping table, and cross-receptor pocket comparison table.
- [ ] Each output preserves enough raw detail to support later manual scientific interpretation rather than only providing black-box rankings.
- [ ] The receptor-state-specific pocket landscape can be reviewed without manually re-parsing raw docking files.
- [ ] Supporting PyRosetta and AlphaFold-Multimer outputs can be captured in residue-level summary form and inspected alongside the pocket results.
- [ ] The pipeline output is clean enough to hand off to Codex for iterative repository improvement without re-explaining the scientific context from scratch.
- [ ] The workflow is reproducible enough that a future session can resume from saved files and project documents rather than from chat memory alone.

## Key Constraints
- **Existing GitHub codebase must be reused rather than replaced from scratch.**
- **Python is the default implementation language unless explicitly changed later with approval.** fileciteturn8file1L67-L78
- **Phase-based development must be respected:** Brief → PRD → Tasks → Execution, without skipping steps. fileciteturn8file1L10-L28
- **All documents should be written in English first, then followed by a Korean translation section.** fileciteturn8file1L91-L100
- **The document set should remain self-contained and suitable for session recovery.** fileciteturn8file8L30-L48
- **The system must preserve residue numbering consistency across receptor states wherever possible.**
- **Parallel execution must be configurable and safe for shared-server use.**
- **Legacy interpretations from older reports must remain reference material only, not hard-coded ground truth.**

## Definition of “Done” for Phase 0
Phase 0 is complete when:
- the project is clearly defined in one sentence,
- the target user is explicit,
- the MVP core features are capped at five and independently describable,
- out-of-scope boundaries are explicit,
- success criteria are concrete,
- and the user approves this brief as the basis for PRD generation.

---

## 한글 번역 (Korean Translation)

## Context Summary
- 프로젝트: EGFR–MYO1D 연구를 위한 구조 상태 기반 docking, pocket, residue 비교 파이프라인
- 현재 단계: Phase 0 (Project Brief)
- 상태: 브리프 초안 작성 완료, Phase 1 PRD 작성 전 사용자 검토 대기
- 핵심 제약: 기존 GitHub 코드베이스 재사용, 병렬 실행 시 실사용 가능 16코어 기준, 기존 보고서보다 새 계산 결과 우선
- 기술 스택: Python 기본 사용 (세부 스택은 Phase 2에서 확정)

---

# 프로젝트 브리프

## 한 줄 설명
EGFR–MYO1D 연구를 위해 3개의 receptor 구조 상태에서 얻어지는 Vina, PyRosetta, AlphaFold-Multimer 결과를 pose, pocket, residue 수준에서 서로 비교 가능한 표준 데이터셋으로 변환하는 연구용 계산 파이프라인이다.

## 대상 사용자
**주 사용자:** EGFR–MYO1D 프로젝트를 직접 유지·확장하는 연구자 본인.

**보조 사용자:** 지도교수나 공동연구자처럼 receptor 상태별 pocket 비교, residue 요약 검토, 리간드-포켓 매핑 확인이 필요한 소규모 협업자.

### 왜 이 도구가 필요한가
현재 워크플로우는 여러 개의 분리된 스크립트와 수동 해석에 의존하고 있어서 다음이 어렵다.
- receptor 상태 간 결과 비교
- 어떤 약물이 어떤 pocket으로 수렴하는지 추적
- residue 수준 contact 요약의 일관된 추출
- 이후 세션에서 결과 재사용
- Codex 같은 코딩 에이전트에게 프로젝트를 명확히 인계

이 프로젝트의 목적은 산발적인 계산 워크플로우를 반복 가능하고, 검토 가능하며, 코덱스 친화적인 연구 파이프라인으로 바꾸는 것이다.

## 핵심 기능 (MVP 범위) — 최대 5개

1. **구조화된 입력 및 실행 관리**
   - 시스템은 세 개의 receptor 구조(3GT8 원본, MD 클러스터 대표 38–48, MD 클러스터 대표 85–100), 여러 ligand 입력, 공통 실행 설정을 일관된 프로젝트 설정 형식으로 받아야 한다.
   - receptor ID, residue numbering, 실행 메타데이터가 downstream 비교에 문제 없도록 보존되어야 한다.

2. **가용 CPU 범위를 활용한 병렬 batch docking 실행**
   - 시스템은 receptor × ligand 조합 전체를 대상으로 설정 가능한 병렬 실행을 지원해야 한다.
   - 실제 운영 환경은 32코어 서버지만, 이 프로젝트에서는 통상 16코어만 안전하게 사용 가능한 것으로 본다.
   - 사용자는 최대 worker 수를 지정할 수 있어야 하며, 실질적인 일상 상한은 16코어를 기준으로 한다.
   - 병렬 실행은 처리량을 높이되, 출력 구조를 흐트러뜨리거나 추적 가능성을 떨어뜨리면 안 된다.

3. **pose 파싱, contact 추출, pocket clustering**
   - 시스템은 Vina 출력에서 pose 테이블을 표준화해서 만들고, 각 pose의 receptor contact residue를 추출하고, pose centroid를 계산하고, 사용자 정의 기준으로 pose를 pocket으로 clustering해야 한다.
   - pocket assignment는 재현 가능해야 하며, 이후 검토를 위한 raw metric도 함께 남겨야 한다.
   - 출력은 사람이 직접 보기에도 쉽고, 프로그램이 후처리하기에도 쉬워야 한다.

4. **pocket 요약 및 receptor 간 비교**
   - 시스템은 각 receptor 상태별 pocket을 요약해야 하며, pose 수, ligand 수, 대표 pose, residue frequency 요약을 포함해야 한다.
   - 또한 receptor 상태 간 pocket을 위치와 residue overlap 기준으로 비교해서, 같은 patch인지, 부분 overlap인지, 완전히 다른 pocket인지 사용자가 판단할 수 있게 해야 한다.
   - 목표는 자동 결론 강제가 아니라, 과학적 해석을 위한 구조화된 증거 제공이다.

5. **보조 PPI 결과를 포함한 통합 연구 리포트**
   - 시스템은 Vina 기반 pocket 요약과 함께 PyRosetta global docking 및 AlphaFold-Multimer에서 얻은 residue-level 보조 정보를 읽기 쉬운 리포트로 제공해야 한다.
   - 이 보조 모듈들은 pocket 정의의 주된 진실 소스는 아니지만, docking-derived pocket 데이터와 나란히 검토할 수 있는 receptor-side residue 참고자료를 제공해야 한다.
   - 최종 출력은 연구자, 협업자, 그리고 저장소를 개선할 코딩 에이전트가 함께 검토할 수 있어야 한다.

> ⚠️ 각 핵심 기능은 Phase 2에서 서로 독립된 단위로 다뤄진다.
> 숨은 결합을 막기 위해 별도로 계획, 구현, 테스트된다.

## 명시적 비범위
- **웹 애플리케이션 또는 그래픽 대시보드** — 이유: 현재 MVP는 사용자용 제품이 아니라 연구 파이프라인이다.
- **인증, 계정, 권한 관리** — 이유: 멀티테넌트 플랫폼이 아니다.
- **결제, 청구, 구독 시스템** — 이유: MVP와 전혀 관련 없다.
- **프로덕션 배포, CI/CD, 클라우드 인프라, 확장 아키텍처** — 이유: 당장의 목표는 기존 서버에서 작동하는 연구용 MVP다.
- **wet-lab 데이터 자동 통합** — 이유: 현재 범위는 계산 파이프라인 표준화에 한정된다.
- **과학적 해석을 자동 결론으로 대체하는 기능** — 이유: 시스템은 증거를 구조화해야지 해석을 강제로 대신하면 안 된다.
- **기존 보고서 residue/site를 고정 진실로 취급하는 기능** — 이유: 새 계산 결과가 더 높은 해석 우선순위를 가져야 한다.

## 성공 기준
다음 조건이 모두 충족되면 MVP는 성공으로 본다.

- [ ] 연구자가 하나의 설정 파일로 세 개의 receptor 상태를 대상으로 파이프라인을 실행할 수 있다.
- [ ] 시스템이 worker 수를 설정 가능하게 두고, 실사용 16코어 범위 안에서 병렬 batch 실행을 수행할 수 있다.
- [ ] 시스템이 Vina 결과를 표준 출력으로 생성하며, 최소한 pose-level table, pocket-level table, ligand-to-pocket mapping table, cross-receptor pocket comparison table을 만든다.
- [ ] 각 출력이 단순 black-box ranking이 아니라 이후 수동 과학 해석에 필요한 raw detail을 충분히 보존한다.
- [ ] receptor 상태별 pocket landscape를 raw docking 파일을 다시 손으로 파싱하지 않고도 검토할 수 있다.
- [ ] PyRosetta와 AlphaFold-Multimer의 보조 결과도 residue-level summary 형태로 저장되어 pocket 결과와 함께 검토할 수 있다.
- [ ] 결과물이 Codex에게 추가 개선을 맡길 수 있을 정도로 정돈되어 있어, 과학적 맥락을 매번 처음부터 다시 설명하지 않아도 된다.
- [ ] 이후 세션에서 대화 기억이 아니라 저장된 문서와 파일만으로 작업을 이어갈 수 있을 만큼 재현 가능하다.

## 핵심 제약조건
- **기존 GitHub 코드베이스를 버리지 않고 재사용해야 한다.**
- **Python을 기본 구현 언어로 사용한다.** fileciteturn8file1L67-L78
- **개발 단계는 Brief → PRD → Tasks → Execution 순서를 반드시 지켜야 한다.** fileciteturn8file1L10-L28
- **모든 문서는 영어 본문 먼저, 이후 한글 번역 섹션을 붙여야 한다.** fileciteturn8file1L91-L100
- **문서 세트는 self-contained 해야 하며, 새 세션에서도 복구 가능해야 한다.** fileciteturn8file8L30-L48
- **가능한 한 receptor 간 residue numbering 일관성을 유지해야 한다.**
- **병렬 실행은 공유 서버 환경에서 안전하게 설정 가능해야 한다.**
- **기존 보고서 해석은 참고자료일 뿐, 하드코딩된 진실이 되어서는 안 된다.**

## Phase 0 완료 조건
다음이 충족되면 Phase 0이 완료된 것으로 본다.
- 프로젝트가 한 문장으로 명확하게 정의되었다.
- 대상 사용자가 분명하다.
- MVP 핵심 기능이 5개 이하이고 각각 독립적으로 설명 가능하다.
- 비범위가 명시되어 있다.
- 성공 기준이 구체적이다.
- 사용자가 이 브리프를 PRD 작성의 기준 문서로 승인한다.

