# Server Environment Setup

This repository now assumes that the user manages the server conda environment manually.

The official server baseline is:

- conda environment name: `pyrosetta`
- production and pre-qsub checks use the same environment
- repository scripts do not create a separate test environment

## Required packages for pre-qsub checks

The pre-qsub lane expects these packages to be available in `pyrosetta`:

- `pytest`
- `pyyaml`
- `numpy`
- `pandas`
- `scipy`
- `matplotlib`

The package reference list is also stored in [requirements-test.txt](/Users/admin/Desktop/hwang/codex/codex_ligand/requirements-test.txt).

## Recommended manual install

Activate the environment first:

```bash
source ~/.bashrc
conda activate pyrosetta
```

Then install the required packages manually in that environment using the package source available on your server.

Example:

```bash
conda install pytest pyyaml numpy pandas scipy matplotlib
```

If your server is offline, use your site's normal internal mirror, preloaded channel, or administrator-approved package path.

## How repository scripts behave

- `scripts/setup_test_env.sh`
  - does not create an environment
  - only verifies that `pyrosetta` has the required packages

- `scripts/run_pre_qsub_checks.sh`
  - activates `pyrosetta`
  - runs the lightweight pre-qsub validation lane

- `config/run_pre_qsub_checks.pbs`
  - runs the same validation lane through PBS
  - writes `output/pre_qsub_status/last_pass.json`
  - records `status: failed` if the environment check or test run fails
