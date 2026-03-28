# PyRosetta Manual (Current Repository)

Last updated: 2026-03-13

This manual describes the **current** PyRosetta Phase 1 execution paths in this repository. It replaces legacy references to standalone files such as `pipeline_manager.py`, `run_v1.pbs`, and `run_v2_test.pbs` that are not part of the active codebase.

## 1. Scope

Use this document for:

- running the PyRosetta branch from the unified CLI (`main.py`)
- running dedicated PPI PBS jobs for test/production lanes
- checking where Phase 1 outputs are written

For overall operator sequence and stop/go policy, use `docs/runbook.md`.

## 2. Active Entry Points

### 2.1 Unified CLI (recommended)

```bash
python main.py -c config/example-project.yaml pyrosetta
python main.py -c config/example-project.yaml ppi-postprocess
```

- `pyrosetta`: runs the current Phase 1 PyRosetta branch
- `ppi-postprocess`: rebuilds or re-extracts current PPI postprocess outputs

### 2.2 PBS wrappers

```bash
qsub config/run_lightdock.pbs                          # Phase 1 LightDock 전체 state
qsub config/run_lightdock_test.pbs                     # Phase 1 LightDock 테스트
```

- `run_lightdock.pbs`: Phase 1 LightDock secondary validation
- `run_lightdock_test.pbs`: LightDock test lane (3GT8_raw)

Precheck/production wrappers for the full baseline flow:

- `config/run_pre_qsub_checks.pbs`
- `config/run_production.pbs`

## 3. Active Configuration Files

Project-level config:

- `config/example-project.yaml`

Phase 1 state/seed configs:

- `config/phase1/phase1_test_3GT8_raw.ini`
- `config/phase1/phase1_test_EGFR_160-185.ini`
- `config/phase1/phase1_test_EGFR_170-200.ini`
- `config/phase1/phase1_prod_3GT8_*_seed*.ini`

Phase 1 LightDock PBS:

- `config/run_lightdock.pbs`
- `config/run_lightdock_test.pbs`

## 4. Output Locations

Primary Phase 1 outputs are written under:

- `output/workflow_a/phase2_ppi_docking/`

Key artifacts include:

- `output/workflow_a/phase2_ppi_docking/phase1_interface_report.md`
- `output/workflow_a/phase2_ppi_docking/phase1_interface_comparison_report.md`
- `output/workflow_a/phase2_ppi_docking/phase1_pilot_comparison_note.md`
- `output/workflow_a/phase2_ppi_docking/phase1_downstream_patch_reference.csv`
- `output/workflow_a/phase2_ppi_docking/<state>/ppi_interface_patch_table.csv`
- `output/workflow_a/phase2_ppi_docking/<state>/ppi_cluster_summary.csv`
- `output/workflow_a/phase2_ppi_docking/<state>/ppi_hotspot_residues.csv`
- `output/workflow_a/phase2_ppi_docking/<state>/orientation_filter_log.csv` (state-dependent)
- `output/workflow_a/phase2_ppi_docking/<state>/lightdock/lightdock_interface_support_table.csv` (when present)

`ppi-postprocess` 체인 복원 산출물은 기본적으로 다음 경로에 기록됩니다.

- `output/workflow_a/phase3_ppi_postprocess/restored_runs/{receptor_id}/{partner_name}/`

Legacy 위치인 `{docking_dir}/restored/`는 기본 저장 경로가 아니며, 기존 legacy 디렉토리가 이미 있으면
이동 안내용 `MOVED_TO.txt`가 남을 수 있습니다.

## 5. Minimal Run Patterns

### 5.1 Local/interactive Phase 1 run

```bash
python main.py -c config/example-project.yaml pyrosetta
python main.py -c config/example-project.yaml ppi-postprocess
```

### 5.2 Cluster submission

```bash
qsub config/run_lightdock.pbs                          # LightDock 전체 state
qsub config/run_lightdock_test.pbs                     # LightDock 테스트
```

### 5.3 Production baseline with precheck guard

```bash
PRECHECK_JOB=$(qsub config/run_pre_qsub_checks.pbs)
qsub -W depend=afterok:${PRECHECK_JOB} config/run_production.pbs
```

## 6. Validation Checklist

After a Phase 1 run:

1. Confirm `output/workflow_a/phase2_ppi_docking/` exists and is newly updated.
2. Confirm `phase1_interface_report.md` is present.
3. Confirm per-state patch tables exist for expected receptor states.
4. If LightDock validation is expected, confirm lightdock support tables exist.
5. Confirm the downstream handoff file exists: `phase1_downstream_patch_reference.csv`.

## 7. Known Legacy References (do not use)

The following names are from older layouts and are not active entry points in this repository:

- `pipeline_manager.py`
- `movers.py` / `scoring.py` / `pyrosetta_init.py` (as a standalone top-level pipeline set)
- `run_v1.pbs`
- `run_v2_test.pbs`
- standalone `config_10k.ini` / `config_100k.ini` workflows

If another document still references those paths as active commands, treat that section as historical only.

## 8. Related Documents

- `docs/manual_execution.md`
- `docs/runbook.md`
- `docs/first_time_environment_setup.md`
- `docs/phase1_pyrosetta_execution_note.md`
- `config/README.md`
