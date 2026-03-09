# Runbook
## EGFR–MYO1D Pipeline

---

## 1. Purpose of this runbook

This runbook explains how the EGFR–MYO1D pipeline is expected to be run, inspected, and debugged once the repository refactor is underway.

It is not a low-level implementation spec. Instead, it serves as an **operator-facing execution guide** for:
- the project owner,
- collaborators reviewing results,
- and coding agents improving execution flow.

The runbook defines the intended execution sequence, input expectations, output locations, and practical operating rules.

---

## 2. Project scope reminder

This repository is a research pipeline for comparing structural-analysis outputs across three EGFR receptor states in the EGFR–MYO1D project.

The current receptor ensemble is fixed to:

1. **3GT8 raw structure**
2. **MD cluster representative from frames 38–48**
3. **MD cluster representative from frames 85–100**

The current implementation priority is Vina-centered:

1. batch docking execution
2. pose parsing
3. contact residue extraction
4. pocket clustering
5. pocket summary generation
6. cross-receptor pocket comparison
7. supporting PyRosetta / AlphaFold-Multimer summaries

---

## 3. Practical execution constraint

The main server has **32 CPU cores**, but the practical safe assumption for this project is that only **16 cores are available for routine use**.

Because of this:
- parallel execution must be configurable,
- 16 workers should be treated as the normal upper operating bound,
- and execution speed must never come at the cost of traceability or result clarity.

If there is uncertainty, use fewer workers rather than more.

---

## 4. Recommended document reading order before running anything

Before executing the pipeline, read the following documents in this order:

1. `README.md`
2. `docs/project-context.md`
3. `docs/brief-egfr-myo1d-pipeline.md`
4. `docs/prd-egfr-myo1d-pipeline.md`
5. `docs/tasks-egfr-myo1d-pipeline.md`
6. `CODEX_HANDOFF_EGFR_MYO1D_PIPELINE.md`

This ensures that the operator understands:
- the project goal,
- the three receptor states,
- the current implementation priority,
- and the difference between primary Vina evidence and auxiliary PPI evidence.

---

## 5. Expected repository structure

The intended repository structure is:

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
├── receptors/
├── ligands/
├── scripts/
├── outputs/
│   ├── raw/
│   ├── parsed/
│   ├── reports/
│   └── logs/
└── tests/
```

This runbook assumes the repository gradually converges toward this structure.

---

## 6. Input expectations

The intended input model is project-level and file-based.

### 6.1 Receptors
The pipeline should accept exactly three named receptor states during the current phase:

- `3GT8_raw`
- `3GT8_cl38_48`
- `3GT8_cl85_100`

Each receptor definition should preserve:
- receptor ID
- source type
- PDB path
- PDBQT path
- chain information if applicable
- optional preparation notes

### 6.2 Ligands
The pipeline should accept a ligand list through config, with each ligand preserving at least:
- ligand ID
- file path
- optional scientific annotation

### 6.3 Shared run settings
The pipeline should be runnable from a single project config file that stores:
- receptor definitions
- ligand definitions
- docking settings
- contact cutoff
- pocket clustering cutoff
- worker count
- output root

---

## 7. Example config skeleton

The exact implementation may vary, but the intended config model looks like this:

```yaml
project_name: egfr_myo1d_pipeline
output_root: ./outputs
max_workers: 16
receptors:
  - id: 3GT8_raw
    pdb: receptors/3GT8_raw/3GT8_raw.pdb
    pdbqt: receptors/3GT8_raw/3GT8_raw.pdbqt
    chain: A
    source_type: raw_3GT8
  - id: 3GT8_cl38_48
    pdb: receptors/cluster_38_48/rep_38_48.pdb
    pdbqt: receptors/cluster_38_48/rep_38_48.pdbqt
    chain: A
    source_type: md_cluster_38_48
  - id: 3GT8_cl85_100
    pdb: receptors/cluster_85_100/rep_85_100.pdb
    pdbqt: receptors/cluster_85_100/rep_85_100.pdbqt
    chain: A
    source_type: md_cluster_85_100
ligands:
  - id: drugA
    pdbqt: ligands/drugA.pdbqt
    annotation_myo1d_inhibition: strong
  - id: drugB
    pdbqt: ligands/drugB.pdbqt
    annotation_myo1d_inhibition: medium
vina:
  center: [0.0, 0.0, 0.0]
  size: [30.0, 30.0, 30.0]
  exhaustiveness: 32
  num_modes: 20
  energy_range: 6
contacts:
  cutoff: 4.0
pocket_clustering:
  centroid_cutoff: 4.0
```

The exact schema may evolve, but this is the intended operational model.

---

## 8. Recommended execution order

### Phase A: Repository and input sanity check
Do this before running any heavy computation.

1. Confirm that the repository contains the expected documents.
2. Confirm that receptor file paths exist.
3. Confirm that ligand file paths exist.
4. Confirm that config syntax is valid.
5. Confirm that receptor IDs are exactly what the pipeline expects.
6. Confirm that worker count is appropriate for the current shared-server load.

### Phase B: Run Vina batch execution
This is the current highest-priority execution layer.

Intended command shape:

```bash
python scripts/run_vina_batch.py --config config/example-project.yaml
```

or, if a unified runner exists later:

```bash
python scripts/run_pipeline.py --config config/example-project.yaml --steps vina
```

### Phase C: Parse pose-level outputs
After docking outputs are produced, parse them into a reusable table.

Intended command shape:

```bash
python scripts/parse_vina_results.py --config config/example-project.yaml
```

### Phase D: Extract contact residues
After pose parsing, extract receptor contact residues.

Intended command shape:

```bash
python scripts/extract_contacts.py --config config/example-project.yaml
```

### Phase E: Cluster pockets and summarize them
After pose-level data exists, generate pocket-level summaries.

Intended command shape:

```bash
python scripts/cluster_pockets.py --config config/example-project.yaml
python scripts/summarize_pockets.py --config config/example-project.yaml
```

### Phase F: Compare pockets across receptor states
After pocket tables exist, generate cross-receptor comparison outputs.

Intended command shape:

```bash
python scripts/compare_pockets.py --config config/example-project.yaml
```

### Phase G: Build reports
Once core parsed outputs exist, generate readable summaries.

Intended command shape:

```bash
python scripts/build_vina_report.py --config config/example-project.yaml
```

### Phase H: Optional supporting PPI modules
Run only after the Vina-centered layer is stable.

Possible later commands:

```bash
python scripts/run_pyrosetta_global_docking.py --config config/example-project.yaml
python scripts/parse_pyrosetta_scores.py --config config/example-project.yaml
python scripts/run_afm_batch.py --config config/example-project.yaml
python scripts/parse_afm_outputs.py --config config/example-project.yaml
```

---

## 9. Recommended first implementation scope

Until the repository is stabilized, the recommended first implementation scope is limited to:

- **Task Group 0: Project Setup and Repository Baseline**
- **Task Group 1: Structured Input and Run Management**
- **Task Group 2: Parallel Batch Docking Execution**

This means the first real execution milestone is:

- config-driven Vina batch execution
- configurable `max_workers`
- safe 16-core operation
- stable receptor/ligand-specific output placement
- visible failure logging

Do not expand into pocket comparison, integrated reporting, or PPI standardization until this layer is stable.

---

## 10. Expected output locations

### 10.1 Raw outputs
Raw execution outputs should go under:

```text
outputs/raw/
```

This includes:
- raw Vina output files
- raw logs for individual jobs
- optional intermediate outputs from supporting tools

### 10.2 Parsed outputs
Structured machine-readable outputs should go under:

```text
outputs/parsed/
```

Expected files include:
- `vina_pose_table.csv`
- `vina_pocket_table.csv`
- `vina_drug_pocket_map.csv`
- `vina_pocket_overlap_table.csv`
- and later, standardized PyRosetta / AFM summary files

### 10.3 Reports
Human-readable markdown summaries should go under:

```text
outputs/reports/
```

Expected files include:
- `vina_report.md`
- later `pyrosetta_report.md`
- later `afm_report.md`
- later `integrated_report.md`

### 10.4 Logs
Execution logs should go under:

```text
outputs/logs/
```

This directory should preserve:
- batch-level logs
- per-job logs when applicable
- failure visibility

---

## 11. What to check after each run

After any meaningful run, check the following.

### 11.1 Basic completion checks
- Did the script exit cleanly?
- Are output directories present?
- Are expected files created?
- Are receptor IDs correct in outputs?
- Are ligand IDs correct in outputs?

### 11.2 Parallel run checks
- Was `max_workers` honored?
- Were outputs separated correctly by receptor and ligand?
- Did any job fail silently?
- Did logs clearly identify failures?

### 11.3 Parsed output checks
- Does every pose row include receptor and ligand identity?
- Are centroid coordinates present?
- Are contact residues present where expected?
- Are pocket assignments reproducible when rerun?

### 11.4 Comparison checks
- Does the pocket comparison table exist?
- Does it preserve raw metrics rather than only labels?
- Can a human review same-patch vs distinct-pocket possibilities without reopening raw docking files?

---

## 12. Failure handling expectations

The pipeline should fail loudly and specifically, not silently.

### Examples of failures that must be surfaced clearly
- missing receptor file
- missing ligand file
- malformed config
- malformed or empty docking output
- residue numbering inconsistency warning
- output path collision
- failed docking subprocess

When a failure happens, the operator should be able to answer:
- what step failed,
- which receptor/ligand pair failed,
- where the relevant logs are,
- and whether the rest of the run continued or stopped.

---

## 13. Interpretation guardrails during operation

While running and reviewing outputs, the operator should remember the following rules.

### 13.1 New outputs are primary
Newly generated structured outputs carry more interpretive weight than older manually labeled report residues or site names.

### 13.2 Pocket identity is conditional
Pocket comparison should preserve evidence rather than force absolute conclusions.
A pocket in one receptor state may shift, overlap partially, or disappear in another state.

### 13.3 Primary vs auxiliary evidence must remain separate
- Vina-centered pocket outputs are the current primary evidence layer.
- PyRosetta and AlphaFold-Multimer outputs are supporting evidence.

The run workflow and reports should preserve this distinction.

---

## 14. Minimal operator checklist before handing work to Codex

Before asking Codex to modify code, make sure the repo contains at least:

- `README.md`
- `docs/project-context.md`
- `docs/brief-egfr-myo1d-pipeline.md`
- `docs/prd-egfr-myo1d-pipeline.md`
- `docs/tasks-egfr-myo1d-pipeline.md`
- `CODEX_HANDOFF_EGFR_MYO1D_PIPELINE.md`
- `docs/runbook.md`

And make sure Codex is explicitly told to start only with:
- Task Group 0
- Task Group 1
- Task Group 2

---

## 15. Recommended Codex start condition

The cleanest starting condition is:

- repository documents are placed correctly,
- receptor and ligand paths are known,
- a draft config file exists,
- and Codex is instructed to first inspect the repo and propose a narrow refactor plan.

This reduces the risk of premature full rewrites.

---

## 16. Korean summary (간단 요약)

이 runbook은 이 저장소를 실제로 어떻게 돌릴지를 정리한 운영 문서다.

핵심 내용은 다음과 같다.
- receptor는 3GT8 raw / 38–48 cluster rep / 85–100 cluster rep 세 개다.
- 현재 최우선은 Vina 중심 실행과 표준화다.
- 서버는 32코어지만 16코어만 실사용 기준으로 병렬 실행을 설정한다.
- 먼저 해야 할 것은 config 기반 batch 실행, pose parsing, contact extraction, pocket clustering의 실행 흐름을 고정하는 것이다.
- raw / parsed / reports / logs 출력 위치를 분리해야 한다.
- 실패는 조용히 넘어가면 안 되고, 어느 단계에서 무엇이 실패했는지 보여야 한다.
- Codex에게 작업을 넘기기 전에는 최소 문서 세트를 repo에 두어야 한다.

