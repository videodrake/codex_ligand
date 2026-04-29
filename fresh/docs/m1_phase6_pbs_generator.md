# M1 Phase 6 — PBS Generator

Closes M1 §23 #6 (qsub smoke can be generated) per `milestone1_foundation_codex_handoff_v0_5.md` §10.

Actual `qsub` execution is HPC-pending (user-side step). This phase validates concrete PBS file generation.

## What it does

`fresh/src/egfr_myo1d/hpc/pbs.py` renders a concrete PBS job file under `runs/<run_id>/scripts/<job_name>.pbs` with:

- Absolute stdout/stderr paths (no PBS variable expansion)
- `set -euo pipefail`
- `cd <ABSOLUTE_REPO_ROOT>`
- `source <conda_sh> && conda activate pyrosetta`
- Locale and `PYTHONPATH=$REPO_ROOT/fresh/src` exports
- All five thread limits (`OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS`, `NUMEXPR_NUM_THREADS`, `VECLIB_MAXIMUM_THREADS`) set to 1
- Mode-specific command body:
  - `smoke_env`: `cli preflight --mode smoke_env --profile hpc_strict` + `cli cleanup --mode test --dry-run false`
  - `smoke_input`: `cli prepare-inputs --mode smoke_input --profile hpc_strict --input-root fresh/data/raw` + `cli cleanup --mode test --dry-run false`
  - `mini` / `scaling` / `production`: placeholder for user-supplied M2/M3 commands

Mode → ppn / walltime mapping is read from `fresh/configs/hpc.yaml` (no hardcoding). `smoke_env` and `smoke_input` both map to the `smoke` key in hpc.yaml.

The module **never calls qsub**. A regression test (`test_no_auto_qsub_call_in_module_or_scripts`) enforces this for both `hpc/pbs.py` and the wrapper scripts. PBS files use Unix line endings (`\n` only).

## CLI

```bash
python -m egfr_myo1d.cli prepare-pbs \
    --run-id RUN \
    --job-name NAME \
    --mode smoke_env|smoke_input|mini|scaling|production \
    [--node node04|node05|node06] \
    [--ppn N] \
    [--walltime HH:MM:SS] \
    [--output-path PATH] \
    [--input-root fresh/data/raw] \
    [--profile codex_dev|hpc_strict]
```

Wrapper scripts:

```bash
bash fresh/scripts/submit_smoke_env.sh   [<run_id>] [<node>]
bash fresh/scripts/submit_smoke_input.sh [<run_id>] [<node>]
```

Both scripts emit a PBS file under `fresh/runs/<run_id>/scripts/` and print the `qsub <path>` instruction for the user to run manually on HPC. Neither auto-calls qsub.

## Module additions

```text
fresh/src/egfr_myo1d/hpc/__init__.py    new
fresh/src/egfr_myo1d/hpc/pbs.py         new
```

Public API:

```python
THREAD_ENV_KEYS = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                   "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")
DEFAULT_CONDA_SH = "/usr/local/anaconda/3/2023.09/etc/profile.d/conda.sh"
DEFAULT_CONDA_ENV = "pyrosetta"
WALLTIME_PATTERN = re.compile(r"^\d+:\d{2}:\d{2}$")

PBSFile (dataclass)
load_hpc_gate(ctx) -> dict
render_pbs_content(...) -> str             # pure function used by tests
derive_smoke_env_command_lines(run_id) -> list[str]
derive_smoke_input_command_lines(run_id, input_root) -> list[str]
derive_command_lines_for_mode(mode, run_id, input_root) -> list[str]
generate_pbs(ctx, job_name, mode, node=None, ppn=None, walltime=None,
             output_path=None, profile="codex_dev", input_root="fresh/data/raw",
             command_lines=None) -> PBSFile
```

## Severity

| Status | Conditions |
| --- | --- |
| `PASS` | PBS file generated cleanly; no overrides |
| `WARN` | `--ppn` or `--walltime` overrides hpc.yaml default (recorded but accepted) |
| `FAIL` | unknown node, unknown mode, malformed walltime, malformed job name, `--output-path` outside run_dir |

(The dataclass `PBSFile.status` reflects per-call result. The CLI handler maps WARN → exit 0 with a message; FAIL is raised as `ValueError` and caught by the CLI top-level error handler returning exit 2.)

## Outputs

```text
fresh/runs/<run_id>/scripts/<job_name>.pbs
fresh/runs/<run_id>/logs/job_status.jsonl     (appended; status=GENERATED)
fresh/runs/<run_id>/logs/phase_status.jsonl   (appended; phase=prepare-pbs)
fresh/runs/<run_id>/logs/master.log           (appended)
```

## PBS file structure (handoff §10.1)

```bash
#!/bin/bash
#PBS -q workq
#PBS -N <job_name>
#PBS -l nodes=<node>:ppn=<n>
#PBS -l walltime=<HH:MM:SS>
#PBS -V
#PBS -o /ABS/REPO/fresh/runs/<run_id>/logs/jobs/<job_name>.stdout
#PBS -e /ABS/REPO/fresh/runs/<run_id>/logs/jobs/<job_name>.stderr

set -euo pipefail
cd /ABS/REPO

source /usr/local/anaconda/3/2023.09/etc/profile.d/conda.sh
conda activate pyrosetta

export LANG=C.UTF-8
export LC_ALL=C.UTF-8
export PYTHONIOENCODING=utf-8
export REPO_ROOT="/ABS/REPO"
export PYTHONPATH="$REPO_ROOT/fresh/src:${PYTHONPATH:-}"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

<command_lines for the mode>
```

## What is intentionally not in this phase

- Actual `qsub` execution (HPC-pending; user-side `bash fresh/scripts/submit_smoke_*.sh` then `qsub <generated.pbs>`)
- M2/M3 command bodies (mini/scaling/production modes emit placeholder commands; the user fills in real commands later)
- ligand manifest (Phase 7)
- prepare-inputs orchestrator (Phase 8)
- Tasks 4-9 schema realignment (Phase 9)

## HPC-pending validation (user-side)

```bash
# On HPC after placing real input files (if smoke_input):
bash fresh/scripts/submit_smoke_env.sh
qsub fresh/runs/<run_id>/scripts/<job>.pbs
qstat
# After completion
python -m egfr_myo1d.cli status --run-id <run_id>
```
