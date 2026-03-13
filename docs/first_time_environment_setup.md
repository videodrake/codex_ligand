# First-Time Environment Setup

This document is the installation and environment guide for someone using this
repository for the first time.

The goal is simple:

- know which environment this project expects
- know which tools are mandatory versus optional
- know which external binaries must be available on `PATH`
- verify the setup before submitting heavier jobs

This repository does **not** create or manage your production environment for
you. The current baseline assumes that you manage the shared server environment
manually.

## 1. Current Baseline

The current runtime baseline is:

- shared conda environment name: `pyrosetta`
- production and pre-qsub validation use the same environment
- LightDock is the active Phase 1 secondary validation path
- AlphaFold-Multimer is not part of the routine baseline

If you are starting from zero, assume you should prepare **one** working
environment called `pyrosetta` and install the tools you actually need into
that environment or make them available on the server `PATH`.

## 2. What Every First-Time User Needs

Before thinking about Vina, PyRosetta, or LightDock, make sure the following
baseline exists.

### 2.1 Conda environment

Use a shared environment named `pyrosetta`:

```bash
source ~/.bashrc
conda activate pyrosetta
```

This repository currently assumes that:

- `scripts/setup_test_env.sh`
- `scripts/run_pre_qsub_checks.sh`
- `config/run_pre_qsub_checks.pbs`
- `config/run_production.pbs`

all run inside that environment.

### 2.2 Baseline Python packages

These packages are required even if you are only doing the lightweight
validation lane:

- `pytest`
- `pyyaml`
- `numpy`
- `pandas`
- `scipy`
- `matplotlib`

These are the packages checked by:

- `scripts/setup_test_env.sh`
- `docs/pre_qsub_test_line.md`
- `docs/server_environment_setup.md`

## 3. Tool Installation by Pipeline Stage

The easiest way to avoid confusion is to separate the project into
installation profiles.

### 3.1 Profile A: Pre-qsub validation only

If you only want to verify that the repository is healthy before a large
submission, you need:

- the shared `pyrosetta` conda environment
- the baseline Python packages listed above

This profile is enough for:

- `bash scripts/setup_test_env.sh`
- `bash scripts/run_pre_qsub_checks.sh`
- `qsub config/run_pre_qsub_checks.pbs`

It is **not** enough for full production docking.

### 3.2 Profile B: Current active baseline

If you want to run the current active project workflow, install or provide the
following.

#### A. PyRosetta

Required for:

- Phase 1 receptor-side PPI mapping
- `qsub config/run_ppi_test.pbs`
- `qsub config/run_ppi_prod.pbs`
- `qsub config/run_production.pbs` when Phase 2/3 PPI steps are involved

Notes:

- PyRosetta is not installed by this repository
- you must install it yourself in the `pyrosetta` environment
- license and installation method depend on your local lab/server policy

#### B. AutoDock Vina Python API

Required for:

- `python main.py vina`
- Vina docking inside `run_production.py`

Expected Python package:

- `vina`

#### C. Ligand/receptor PDBQT preparation tools

Vina uses `.pdbqt` files. The repository can reuse prebuilt `.pdbqt` files, but
if they are missing it will try to prepare them automatically.

Current fallback order in `egfr_pipeline/vina/dock.py` is:

For receptor preparation:

- `prepare_receptor` (ADFR)
- `prepare_receptor4.py` via `pythonsh` (MGLTools)
- `obabel`

For ligand preparation:

- `mk_prepare_ligand.py` (Meeko)
- `prepare_ligand` (ADFR)
- `prepare_ligand4.py` via `pythonsh` (MGLTools)
- `obabel`

For practical use, install at least:

- one receptor-preparation tool
- one ligand-preparation tool

The most robust combinations are:

- Meeko + OpenBabel
- ADFR suite + OpenBabel

If your input `.pdbqt` files are already prepared and valid, these tools are
less critical for routine reruns.

#### D. RDKit

Helpful for:

- ligand file handling
- SDF-oriented workflows
- some local molecule preparation paths

Package:

- `rdkit`

This is strongly recommended for the current repository, even if a subset of
the workflow may still run without it.

#### E. LightDock

Required for:

- active Phase 1 secondary validation
- `egfr_pipeline.phase1.lightdock_validation`
- generation and later execution of `run_lightdock_<state>.sh`

The repository currently expects LightDock command-line tools such as:

- `lightdock3_setup.py`
- `lightdock3.py`
- `lgd_generate_conformations.py`
- `lgd_cluster_bsas.py`

The implementation note in
`egfr_pipeline/phase1/lightdock_validation.py` currently documents two common
install routes:

- `pip install lightdock3`
- `conda install -c bioconda lightdock`

Important practical rule:

- it is not enough that the Python package exists
- the LightDock command-line executables must also be callable on `PATH`

### 3.3 Profile C: Extended research tools

These are not required for the minimal pre-qsub lane, but they matter for
extended or later-stage work.

#### A. Pocket proposal tools

Phase 2 pocket proposal modules are designed around:

- `fpocket`
- `P2Rank`
- optionally later hotspot-style evidence such as FTMap

Current state:

- the code supports parsing and setup generation
- these tools may be run server-side rather than directly from the default CLI

If you plan to use the full Phase 2 branch, install:

- `fpocket`
- `P2Rank`

#### B. MD stack

MD is not the first-trust onboarding surface and is not automatically executed
by the current routine production path, but it is part of the broader research
pipeline.

For MD-related work, you likely need:

- `gromacs` / `gmx`
- `MDAnalysis`
- optionally `gmx_MMPBSA` or your site-standard equivalent

## 4. Minimal Installation Matrix

| Component | Pre-qsub only | Current active baseline | Extended research |
|---|---|---|---|
| `pytest`, `pyyaml`, `numpy`, `pandas`, `scipy`, `matplotlib` | required | required | required |
| `vina` Python package | not required | required | required |
| PDBQT prep tools (`mk_prepare_ligand.py`, `prepare_receptor`, `obabel`, etc.) | not required | recommended, often required | recommended |
| `rdkit` | optional | recommended | recommended |
| PyRosetta | not required | required | required |
| LightDock | not required | required for active Phase 1 validation | required if using Phase 1 convergence |
| `fpocket` / `P2Rank` | not required | optional | required for Phase 2 pocket proposal work |
| `gromacs` / `MDAnalysis` | not required | optional | required for MD stability work |

## 5. First-Time Verification Checklist

Once you believe the environment is ready, check it in this order.

### 5.1 Activate the shared environment

```bash
source ~/.bashrc
conda activate pyrosetta
```

### 5.2 Check baseline repository health

```bash
cd ~/codex_ligand
bash scripts/setup_test_env.sh
bash scripts/run_pre_qsub_checks.sh
```

Or through PBS:

```bash
qsub config/run_pre_qsub_checks.pbs
```

### 5.3 Check the main CLI

```bash
python main.py --help
python main.py validate --help
```

### 5.4 Check LightDock availability

At minimum, these commands should resolve without “command not found”:

```bash
which lightdock3_setup.py
which lightdock3.py
which lgd_generate_conformations.py
which lgd_cluster_bsas.py
```

You can also confirm the repository module is reachable:

```bash
python -m egfr_pipeline.phase1.lightdock_validation --help
```

### 5.5 Check PyRosetta availability

If your site policy allows a quick import test:

```bash
python -c "import pyrosetta; print('PyRosetta import OK')"
```

### 5.6 Check Vina stack

If you plan to run ligand docking:

```bash
python -c "import vina; print('vina import OK')"
which obabel
```

If you rely on Meeko or ADFR preparation:

```bash
which mk_prepare_ligand.py
which prepare_receptor
which prepare_ligand
```

## 6. What This Repository Does Not Do For You

This repository currently does **not**:

- create a conda environment automatically
- install PyRosetta automatically
- install LightDock automatically
- install server-side external tools automatically
- choose your package mirror or offline package source

You must manage those pieces according to your lab/server policy.

## 7. Recommended Read Order After Environment Setup

After the environment is ready, read these next:

1. `docs/AI_START_HERE.md`
2. `docs/current_pipeline_status.md`
3. `docs/runbook.md`
4. `docs/data_flow_guide.md`
5. `docs/phase1_lightdock_validation_note.md`
6. `docs/server_environment_setup.md`
7. `docs/pre_qsub_test_line.md`

## 8. Common First-Time Failure Modes

Most early failures fall into one of these groups:

- `pyrosetta` environment is active, but required Python packages are missing
- LightDock package is installed, but its command-line tools are not on `PATH`
- Vina input conversion fails because no PDBQT preparation tool is available
- users assume AFM is still part of the routine baseline
- users assume MD is part of the default production path when it is currently a downstream gate

## 9. If You Only Want One Practical Rule

For a first successful setup, make sure this is true:

- shared conda env `pyrosetta` is active
- baseline Python packages are installed
- PyRosetta works
- `vina` works
- at least one PDBQT-preparation route works
- LightDock executables are on `PATH`
- pre-qsub checks pass before heavy submission

