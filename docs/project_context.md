# Project Context
## EGFR–MYO1D Pipeline

---

## 1. Why this project exists

This project exists to standardize the computational workflow used in the EGFR–MYO1D study.

The current analysis environment contains multiple scripts and partial tools for:
- AutoDock Vina ligand docking,
- PyRosetta global docking,
- and local AlphaFold-Multimer runs.

However, those tools have not yet been fully unified into a single reproducible research pipeline.
The main problem is not the absence of computation, but the lack of **stable, comparable, and reviewable outputs**.

This repository is meant to solve that problem.

---

## 2. Current scientific frame

The current scientific focus is the **EGFR kinase domain C-lobe** in the context of a broader **EGFR–MYO1D interaction study**.

The immediate computational goal is to support analysis of:
- receptor-state-specific ligand docking behavior,
- pocket formation and clustering,
- residue-level contact evidence,
- and cross-state comparison of pockets and residue patches.

This repository is therefore not a generic docking project. It is a **state-comparison research pipeline**.

---

## 3. Fixed receptor states

The current receptor ensemble is explicitly limited to these three structures:

1. **3GT8 raw structure**
2. **MD cluster representative from frames 38–48**
3. **MD cluster representative from frames 85–100**

All major pipeline logic should assume that these are the core comparison states unless the project is intentionally expanded later.

These receptor states must remain clearly separated in all outputs.
Residue numbering consistency across them is a high-priority requirement.

---

## 4. Current computational priority

The current near-term priority is **not** to solve every structure question at once.
The priority is to build a reliable evidence layer around the Vina-centered workflow first.

That means the implementation order is intentionally biased toward:

1. Vina batch execution
2. Pose parsing
3. Contact residue extraction
4. Pocket clustering
5. Pocket summary generation
6. Cross-receptor pocket comparison

PyRosetta and AlphaFold-Multimer are important, but at this stage they are treated as **supporting structural evidence modules**, not the main pocket-definition engine.

---

## 5. Interpretation rule: new outputs outrank old labels

A very important project rule is this:

> Legacy residue/site interpretations from older reports are not fixed truth.

Older report residue labels, pocket labels, or named sites may still be useful as historical reference, but they must remain **secondary**.

The current repository should prioritize:
- newly generated pose data,
- newly generated pocket assignments,
- receptor-state-specific overlap evidence,
- and residue-level metrics extracted directly from current runs.

In practical terms, this means:
- do not hard-code old site names into logic,
- do not assume old labels are automatically correct,
- and do not let legacy naming override current structured output.

---

## 6. Pocket interpretation philosophy

Pocket analysis in this project is evidence-driven, not conclusion-driven.

This means the system should help answer questions like:
- does a pocket recur across receptor states?
- does it shift?
- does it partially overlap with another pocket?
- is it clearly distinct?

But the system should **not** pretend that every pocket relationship can be fully decided automatically.

In particular:
- if pockets occupy different locations, global “common residue” logic may be meaningless,
- therefore pocket comparison should preserve raw metrics such as location and residue overlap,
- and same-patch interpretation should remain conditional rather than automatic truth.

---

## 7. Current operational constraint

The main execution environment is a server with **32 CPU cores**, but for practical shared use, the project should assume only **16 cores are safely usable** during routine execution.

This constraint is not optional. It should be treated as a real operating boundary.

That means:
- batch docking must support configurable parallel execution,
- 16 workers should be treated as the practical normal upper bound,
- and performance optimization should not come at the cost of traceability or clean output structure.

---

## 8. Current implementation philosophy

This project should be improved by **refactoring and standardizing the existing GitHub codebase**, not by discarding everything and starting over.

The preferred development style is:
- inspect the current repo first,
- identify the real entry points,
- preserve what works,
- wrap and standardize where needed,
- and expand in small safe steps.

Broad rewrites are discouraged unless absolutely necessary.

---

## 9. What the repository should eventually provide

When the repository matures, it should be able to provide at least the following classes of outputs:

- receptor metadata
- ligand metadata
- Vina pose-level parsed output
- Vina pocket-level summary output
- ligand-to-pocket mapping output
- cross-receptor pocket comparison output
- PyRosetta receptor-side residue summaries
- AlphaFold-Multimer receptor-side residue summaries
- markdown summary reports
- manual-review helper exports where useful

This is the long-term data organization goal.

---

## 10. What the repository is NOT trying to be

This is not:
- a public web application,
- a cloud deployment project,
- a user account system,
- a commercial SaaS tool,
- or an automated wet-lab integration platform.

It is a focused research pipeline for structured computational output generation and comparison.

---

## 11. Recommended first implementation scope

The first implementation pass should remain limited to:

- **Task Group 0: Project Setup and Repository Baseline**
- **Task Group 1: Structured Input and Run Management**
- **Task Group 2: Parallel Batch Docking Execution**

Only after those are stable should the project move into:
- pose parsing,
- contact extraction,
- pocket clustering,
- and cross-receptor comparison.

This protects the repository from premature complexity.

---

## 12. Practical reading order for future contributors

Anyone entering this repository should read documents in this order:

1. `README.md`
2. `docs/brief-egfr-myo1d-pipeline.md`
3. `docs/prd-egfr-myo1d-pipeline.md`
4. `docs/tasks-egfr-myo1d-pipeline.md`
5. `CODEX_HANDOFF_EGFR_MYO1D_PIPELINE.md`

This reading order ensures that the contributor understands:
- the project scope,
- the product requirements,
- the implementation order,
- and the deeper structural-analysis context.

---

## 13. Korean summary (간단 요약)

이 프로젝트는 EGFR–MYO1D 연구를 위한 계산 파이프라인 표준화 프로젝트다.

핵심 요점은 다음과 같다.
- receptor는 3GT8 raw / 38–48 cluster rep / 85–100 cluster rep 세 개다.
- 현재 최우선은 Vina 중심 evidence layer 구축이다.
- PyRosetta와 AFM은 보조 구조 증거 모듈이다.
- 기존 보고서 residue/site 해석은 후순위 참고자료다.
- 새 계산 결과가 더 높은 해석 우선순위를 가진다.
- 서버는 32코어지만 실사용 16코어 기준으로 병렬 실행을 설계해야 한다.
- 전체 재작성보다 기존 GitHub 코드의 점진적 리팩터링이 원칙이다.

