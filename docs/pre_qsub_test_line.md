# Pre-Qsub Test Line

This repository is intended to run on a higher-spec Linux server through `qsub`, but the code should pass a lightweight local validation step before any large run is submitted.

## Goal

Catch the most common pre-submission failures early:

- Python syntax errors
- broken CLI entrypoints
- phase 2 to 4 regression failures
- validation contract breakage
- configuration loading issues

This test line is intentionally lighter than a full production run. It does not launch heavy Vina, PyRosetta production docking, or MD production workflows.

## 1. Verify the existing `pyrosetta` environment

```bash
cd ~/codex_ligand
bash scripts/setup_test_env.sh
```

This script does not create or modify an environment anymore.
It assumes you already use the existing `pyrosetta` conda environment and only checks that the required test packages are installed.

If the package check fails, install the missing packages manually in `pyrosetta` first.
See [server_environment_setup.md](/Users/admin/Desktop/hwang/codex/codex_ligand/docs/server_environment_setup.md).

## 2. Run the pre-qsub checks

```bash
cd ~/codex_ligand
bash scripts/run_pre_qsub_checks.sh
```

If you prefer your normal PBS workflow, use:

```bash
cd ~/codex_ligand
qsub config/run_pre_qsub_checks.pbs
```

On success, the PBS lane now writes a status marker to:

```bash
output/pre_qsub_status/last_pass.json
```

The script runs:

1. `compileall` for `main.py`, `egfr_pipeline/`, and `tests/`
2. UTF-8 CLI smoke checks for:
   - `python main.py --help`
   - `python main.py validate --help`
3. `pytest` on:
   - `tests/test_phase2.py`
   - `tests/test_phase3.py`
   - `tests/test_phase4.py`
   - `tests/test_smoke_cli.py`
   - `tests/test_validation_smoke.py`

## 3. What this means

If this line passes, the repository is in a much safer state for:

- cloning onto the compute server
- preparing production config values
- and submitting heavier `qsub` jobs

If this line fails, fix the code first and do not spend cluster resources yet.

## 4. Suggested server workflow

```bash
git clone <repo-url>
cd codex_ligand
conda activate pyrosetta
bash scripts/setup_test_env.sh
bash scripts/run_pre_qsub_checks.sh
# if green, then proceed to qsub-based runs
```

## 5. Qsub-first workflow

If you want the whole preflight to run through the scheduler instead of directly through `bash`, use:

```bash
cd ~/codex_ligand
qsub config/run_pre_qsub_checks.pbs
```

After it passes, submit the real production job separately:

```bash
qsub config/run_production.pbs
```

`config/run_production.pbs` now checks for the success marker by default.
If the marker is missing, it exits early and prints the recommended dependency
submission command instead of starting the heavy production run.

If you want to chain them safely, use PBS job dependency:

```bash
PRECHECK_JOB=$(qsub config/run_pre_qsub_checks.pbs)
qsub -W depend=afterok:${PRECHECK_JOB} config/run_production.pbs
```

This is the best fit if you usually work entirely through `qsub`. It prevents the production run from starting unless the precheck job exits successfully.

If you intentionally need to bypass the guard, do it explicitly:

```bash
qsub -v SKIP_PRECHECK_GUARD=1 config/run_production.pbs
```

## 6. Current server baseline

The repository no longer treats `.venv-tests` or a separate test-only conda env as the official server test path.
The current baseline is:

- shared production environment: `pyrosetta`
- manual package installation by the user
- package reference list: [requirements-test.txt](/Users/admin/Desktop/hwang/codex/codex_ligand/requirements-test.txt)
