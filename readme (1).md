# EGFR–MYO1D Pipeline

Research pipeline for standardized docking, pocket clustering, residue extraction, and cross-receptor comparison across three EGFR receptor states.

---

## What this repository is for

This repository supports an internal research workflow for the EGFR–MYO1D project.

The main goal is to convert scattered structural-analysis scripts into a reproducible, comparable, and reviewable pipeline that can:

- run ligand docking across multiple EGFR receptor states,
- standardize raw docking outputs into reusable pose- and pocket-level datasets,
- compare receptor-state-dependent pockets and residue patterns,
- preserve outputs in forms that are readable by both humans and downstream tools,
- and support iterative repository improvement by coding agents such as Codex.

This is a **research pipeline**, not a public software product.

---

## Current scientific focus

The current computational focus is the **EGFR kinase domain C-lobe** in the context of **EGFR–MYO1D analysis**.

The pipeline is designed to help compare:
- receptor-state-specific pocket behavior,
- ligand-to-pocket assignments,
- residue-level contact evidence,
- and supporting receptor-side structural signals from PPI-oriented tools.

The current receptor ensemble is fixed to the following three states:

1. **3GT8 raw structure**
2. **MD cluster representative from frames 38–48**
3. **MD cluster representative from frames 85–100**

---

## Current implementation priority

The near-term implementation priority is:

1. **Vina-centered workflow standardization**
2. **Pose parsing and contact extraction**
3. **Pocket clustering and pocket summary generation**
4. **Cross-receptor pocket comparison**
5. **Supporting PyRosetta / AlphaFold-Multimer output standardization**

In other words, the repository is currently centered on making the **Vina-based evidence layer** reliable first.

---

## Practical compute constraint

The working server has **32 CPU cores**, but for this project the practical safe assumption is that only **16 cores are available for routine use**.

Because of that:
- parallel execution must be configurable,
- 16 workers should be treated as the practical normal upper bound,
- and performance improvements must not break traceability or output clarity.

---

## Important interpretation rule

Legacy residue/site interpretations from older reports are **reference material only**.

Newly generated computational outputs should be treated as higher-priority evidence than legacy labels.
The repository should therefore avoid hard-coding older site names or treating older report conclusions as fixed truth.

---

## Start here

If you are new to this repository, read these documents in this order:

1. `docs/brief-egfr-myo1d-pipeline.md`
2. `docs/prd-egfr-myo1d-pipeline.md`
3. `docs/tasks-egfr-myo1d-pipeline.md`
4. `CODEX_HANDOFF_EGFR_MYO1D_PIPELINE.md`

If you are a coding agent or developer, do **not** begin with a full rewrite.
Start with repository inspection and the first implementation batch only.

---

## Recommended first implementation scope

The recommended first implementation scope is limited to:

- **Task Group 0: Project Setup and Repository Baseline**
- **Task Group 1: Structured Input and Run Management**
- **Task Group 2: Parallel Batch Docking Execution**

This means the first implementation pass should focus on:
- understanding the current repository,
- standardizing input/config handling,
- and making Vina batch execution safe, configurable, and traceable.

---

## Recommended repository documents

### Root-level documents
- `README.md`
- `CLAUDE.md`
- `CODEX_HANDOFF_EGFR_MYO1D_PIPELINE.md`

### Documentation directory
- `docs/brief-egfr-myo1d-pipeline.md`
- `docs/prd-outline-egfr-myo1d-pipeline.md`
- `docs/prd-egfr-myo1d-pipeline.md`
- `docs/tasks-outline-egfr-myo1d-pipeline.md`
- `docs/tasks-egfr-myo1d-pipeline.md`
- `docs/project-context.md`
- `docs/repository-map.md`
- `docs/runbook.md`

---

## What this repository should eventually produce

At minimum, the standardized workflow should produce outputs such as:

- pose-level parsed tables,
- pocket-level summary tables,
- ligand-to-pocket mapping tables,
- receptor-to-receptor pocket comparison tables,
- supporting PyRosetta residue summaries,
- supporting AlphaFold-Multimer residue summaries,
- and readable markdown reports.

The exact schema belongs to the implementation layer, but the repository should be organized around these classes of outputs.

---

## Development philosophy

This repository should evolve according to the following principles:

- reuse existing code where possible,
- prefer narrow and inspectable refactors over broad rewrites,
- keep raw evidence visible,
- separate primary docking evidence from auxiliary PPI evidence,
- preserve residue numbering consistency whenever possible,
- and keep the project resumable from files and docs rather than from chat memory.

---

## Current status

This repository has already been documented through:
- a Phase 0 project brief,
- a full PRD,
- a full task breakdown,
- and a Codex handoff specification.

The next practical step is implementation, beginning with the first Vina-centered task batch.

---

## Korean summary (간단 요약)

이 저장소는 EGFR–MYO1D 연구를 위한 계산 파이프라인이다.

현재 핵심은 다음과 같다.
- receptor 3개 상태 비교
- Vina 중심 표준화
- pose / pocket / residue 수준 출력 생성
- receptor 상태 간 pocket 비교
- PyRosetta / AFM 보조 output 정리
- 16코어 병렬 실행 지원

먼저 읽을 문서:
1. `docs/brief-egfr-myo1d-pipeline.md`
2. `docs/prd-egfr-myo1d-pipeline.md`
3. `docs/tasks-egfr-myo1d-pipeline.md`
4. `CODEX_HANDOFF_EGFR_MYO1D_PIPELINE.md`

처음 구현은 Task Group 0~2까지만 시작하는 것이 권장된다.

