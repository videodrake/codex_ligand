# HPC Tool Environment Status

This document records the verified HPC runtime state for the fresh EGFR-MYO1D
PPI-surface pocket workflow. Use it when implementing Milestone 3 adapters and
job scripts.

## Verified Date

- Checked on HPC node: `node04`
- Repository branch: `codex/m2-ppi-input-generation-clean-v2`
- Main conda env: `pyrosetta`
- Additional isolated envs created:
  - `ppi_surface`
  - `p2rank_java11`
  - `pesto`
  - `pocketminer`

## Core `pyrosetta` Env

The `pyrosetta` env remains the primary workflow environment.

Verified OK:

- PyRosetta import
- RDKit import, version `2025.09.2`
- BioPython import, version `1.79`
- numpy import, version `1.23.5`
- pandas import, version `2.3.3`
- AutoDock Vina: `/usr/local/anaconda/3/2023.09/bin/vina`
- Vina version: `AutoDock Vina 52ec525-mod`
- fpocket: `/home/eunae/.conda/envs/pyrosetta/bin/fpocket`
- fpocket version: `4.0`
- mdpocket: `/home/eunae/.conda/envs/pyrosetta/bin/mdpocket`
- Open Babel: `/usr/local/anaconda/3/2023.09/bin/obabel`
- Open Babel version: `3.1.0`
- PBS tools:
  - `/opt/pbs/bin/qsub`
  - `/opt/pbs/bin/qstat`
  - `/opt/pbs/bin/qdel`
- GROMACS: `/home/eunae/.conda/envs/pyrosetta/bin.AVX2_256/gmx`
- GROMACS version: `2025.4-conda_forge`

Use this env for:

- M1/M2/M3 workflow CLI orchestration
- PyRosetta PPI adapters
- Vina ligand docking adapters
- fpocket/mdpocket adapters
- Open Babel ligand conversion
- GROMACS/MSMD lightweight checks, where needed

## `ppi_surface` Env

This env isolates pyKVFinder and its Python dependencies.

Verified OK:

- Python: `3.9.25`
- pyKVFinder: `0.6.16`
- PyYAML: `6.0.3`

Install commands used:

```bash
conda create -n ppi_surface python=3.9 -y
conda activate ppi_surface
python -m pip install pyKVFinder
python -m pip install PyYAML
```

Smoke test:

```bash
python - <<'PY'
import pyKVFinder
print("pyKVFinder import OK")
print(getattr(pyKVFinder, "__version__", "version_unknown"))
PY
```

Use this env for:

- pyKVFinder cavity detection only
- Run through subprocess/PBS adapter from the main workflow

Do not assume the `ppi_surface` env has PyRosetta, Vina, fpocket, RDKit, or
other core tools.

## `p2rank_java11` Env

This env isolates Java 11 for P2Rank. Do not modify system Java.

System Java observed in `pyrosetta`:

```text
/bin/java
openjdk version "1.8.0_352"
```

P2Rank with system Java fails:

```text
UnsupportedClassVersionError
class file version 55.0; Java runtime recognizes up to 52.0
```

Isolated env verified OK:

- Java: OpenJDK `11.0.27`
- P2Rank: `2.4.2`
- P2Rank path: `$HOME/tools/p2rank_2.4.2/prank`

Install/setup commands used:

```bash
conda create -n p2rank_java11 python=3.9 openjdk=11 -y
conda activate p2rank_java11
export PATH="$HOME/tools/p2rank_2.4.2:$PATH"
java -version
prank help
```

Use this env for:

- Optional P2Rank pocket prediction only
- Run through subprocess/PBS adapter with explicit conda activation

Do not change `/bin/java`, global module state, or shared lab Java settings.

## Optional AI / Server Tools

Detailed optional AI runtime notes are in
`fresh/docs/optional_ai_tool_runtime_status.md`.

Verified locally after the first core preflight:

- PeSTo:
  - env: `/home/eunae/.conda/envs/pesto`
  - repo: `/home/eunae/tools/PeSTo`
  - model: `/home/eunae/tools/PeSTo/model/save/i_v4_1_2021-09-07_11-21/model_ckpt.pt`
  - status: dependency import OK; CPU one-PDB prediction OK
- PocketMiner:
  - env: `/home/eunae/.conda/envs/pocketminer`
  - repo: `/home/eunae/tools/gvp_pocketminer`
  - branch/commit: `pocket_pred` / `187062d`
  - checkpoint basename: `/home/eunae/tools/gvp_pocketminer/models/pocketminer`
  - status: source imports OK; checkpoint restore OK; CPU one-PDB prediction OK
  - required env var: `export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"`
- MaSIF:
  - runtime: rootless `podman 4.2.0`
  - image: `docker.io/pablogainza/masif:latest`
  - image id: `6b3c808b7bf7fabdadfee4c6dc2a48c4761b4a118d94983131f34e5a76754a12`
  - status: container pull OK; basic container smoke OK; prediction smoke pending

Not installed locally:

- InDeep

External-server only:

- PASSer

These tools must remain optional. The M2/M3 core workflow must pass without
them.

## PBS Note

The qsub smoke job completed but stderr showed:

```text
GMXRC: line 10: shell: unbound variable
```

This is likely a `set -u` interaction with GROMACS shell initialization. Future
PBS scripts that source GMXRC should wrap that section with `set +u` / `set -u`
or avoid sourcing GMXRC when the conda `gmx` executable is already on `PATH`.

## Adapter Policy

Implement tool adapters with explicit environment boundaries:

- `pyrosetta` env: core workflow, PyRosetta, Vina, fpocket, mdpocket, Open Babel
- `ppi_surface` env: pyKVFinder only
- `p2rank_java11` env: P2Rank only
- `pesto` env: PeSTo only
- `pocketminer` env: PocketMiner only
- rootless `podman`: MaSIF container only

Adapters should record:

- activated env name
- executable or Python module path
- command line
- stdout/stderr log path
- return code
- output manifest path

No adapter should silently fall back from an isolated env to the system tool.
