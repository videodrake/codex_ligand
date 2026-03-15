# Execution Manual

Last updated: 2026-03-12

This document is the command reference for the current repository. Use it when you need the exact shell entry points, PBS submission paths, or the quickest way to run a specific lane. For operator sequencing, stop/go rules, and interpretation boundaries, use [runbook.md](runbook.md). For config field meaning, use [../config/README.md](../config/README.md).

## Baseline Assumptions For Command Use

These command examples assume the current baseline:

- the fixed receptor set is `3GT8_raw`, `EGFR_160-185`, and `EGFR_170-200`
- routine ligand work is Vina-centered
- Phase 1 primary evidence is PyRosetta
- Phase 1 secondary validation is LightDock
- AFM is legacy optional and inactive unless explicitly re-enabled
- production and pre-qsub lanes reuse the shared `pyrosetta` conda environment
- `max_workers = 16` is the routine safe default

## Environment Entry

Recommended server-side baseline:

```bash
cd ~/codex_ligand
conda activate pyrosetta
```

Important note:

- do not assume a separate test-only conda environment is the current baseline

## Config Entry Point

Use the project config explicitly unless you are intentionally targeting another config file.

```bash
config/example-project.yaml
```

Recommended pattern:

```bash
python main.py -c config/example-project.yaml <command>
```

Important: `-c/--config` is a top-level `main.py` option and must be placed before the subcommand.

## Interactive Entry Surface

Open the interactive entry surface:

```bash
python main.py
```

Use this when you want to inspect the available command surface from the current CLI.

## Main CLI Commands

These are the primary current commands exposed by `main.py`.

| Command | Use when | Primary output area |
|------|------|------|
| `python main.py -c config/example-project.yaml vina` | Run docking for the configured receptors and ligands | project raw docking output root |
| `python main.py -c config/example-project.yaml postprocess` | Parse poses, summarize pockets, and compare receptor states | `output/egfr_myo1d_vina/vina/` |
| `python main.py -c config/example-project.yaml pyrosetta` | Run the PyRosetta Phase 1 branch | PyRosetta result directories plus current PPI exports |
| `python main.py -c config/example-project.yaml ppi-postprocess` | Rebuild or extract current PPI postprocess outputs | PPI export areas referenced by the routine baseline |
| `python main.py -c config/example-project.yaml verdict` | Generate the routine site-judgment layer | `output/egfr_myo1d_vina/results/valid_sites.csv` |
| `python main.py -c config/example-project.yaml report` | Generate the routine text and combined evidence report | `output/egfr_myo1d_vina/results/project_report.txt` |
| `python main.py -c config/example-project.yaml validate` | Run validation on the current output state | validation outputs and checks |
| `python main.py -c config/example-project.yaml full` | Run the default integrated CLI path | current routine baseline output tree |

Important interpretation:

- `full` is still aligned to the routine Vina-centered baseline
- it should not be read as proof that the scientific Phase 1 -> 2 -> 3 -> 4 plan is the single default end-to-end path

## Routine Command Sequences

### Minimal Routine Baseline Sequence

Use this when you want the standard current ligand-facing evidence path.

```bash
python main.py -c config/example-project.yaml vina
python main.py -c config/example-project.yaml postprocess
python main.py -c config/example-project.yaml verdict
python main.py -c config/example-project.yaml report
python main.py -c config/example-project.yaml validate
```

Primary review files after this lane:

- `output/egfr_myo1d_vina/results/valid_sites.csv`
- `output/egfr_myo1d_vina/results/cross_method_agreement.csv`
- `output/egfr_myo1d_vina/results/project_report.txt`

### Phase 1-Focused Sequence

Use this when receptor-side patch evidence is the main target.

```bash
python main.py -c config/example-project.yaml pyrosetta
python main.py -c config/example-project.yaml ppi-postprocess
```

Primary review files after this lane:

- `output/phase1_ppi/phase1_downstream_patch_reference.csv`
- `output/phase1_ppi/phase1_interface_report.md`

### Fast Single-Command Path

Use this when the default integrated CLI path is appropriate.

```bash
python main.py -c config/example-project.yaml full
```

Review the same routine baseline result files listed above after completion.

## PBS Submission Paths

These are the current server-oriented submission entry points.

### Precheck Lane

```bash
qsub config/run_pre_qsub_checks.pbs
```

Expected checkpoint:

```text
output/pre_qsub_status/last_pass.json
```

### Production Lane

```bash
qsub config/run_production.pbs
```

### Safe Chained Submission

Use this when you want production to wait for precheck success.

```bash
PRECHECK_JOB=$(qsub config/run_pre_qsub_checks.pbs)
qsub -W depend=afterok:${PRECHECK_JOB} config/run_production.pbs
```

### Phase 1 LightDock PBS Paths

Use these for Phase 1 LightDock secondary validation.

```bash
qsub config/run_lightdock.pbs                          # 전체 state
qsub config/run_lightdock_test.pbs                     # 테스트
```

## Output Orientation By Command

Use this table when you need to jump from a command to the first output location to inspect.

| Command family | Look here first |
|------|------|
| `vina`, `postprocess` | `output/egfr_myo1d_vina/vina/` |
| `verdict`, `report`, `validate`, `full` | `output/egfr_myo1d_vina/results/` |
| `pyrosetta`, `ppi-postprocess` | `output/phase1_ppi/` and current registered PPI result directories |
| pre-qsub PBS | `output/pre_qsub_status/` |
| advanced scientific phases | `output/phase2_pockets/`, `output/phase3_docking/`, `output/phase4_perturbation/` only when those lanes are explicitly in scope |

## Common Command Mistakes

- forgetting `-c config/example-project.yaml` and accidentally relying on an unintended config path
- assuming `full` means the entire scientific 4-phase plan now runs by default
- running AFM-dependent work even though AFM is not active in the baseline config
- treating machine core count as permission to exceed the routine safe worker bound of 16
- reviewing pointer stub files in `output/egfr_myo1d_vina/` instead of the actual payload directories

## Use These Docs Next

- [runbook.md](runbook.md): run order, checkpoints, and escalation rules
- [../config/README.md](../config/README.md): config semantics and field meaning
- [output_artifact_map.md](output_artifact_map.md): what the major artifacts mean
- [data_inventory.md](data_inventory.md): where the inputs and outputs physically live


