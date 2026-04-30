#!/usr/bin/env bash
# Emit a PBS file that discovers/smoke-tests PPI-surface tool availability on HPC.
# This script does not call qsub.
#
# Usage:
#   bash fresh/scripts/submit_tool_preflight.sh [<run_id>] [<node>] [discover|smoke]
set -euo pipefail

REPO_ROOT="${REPO_ROOT:-$(pwd)}"
RUN_ID="${1:-tool_preflight_$(date +%Y%m%d_%H%M%S)}"
NODE="${2:-node04}"
MODE="${3:-smoke}"
JOB_NAME="tool_preflight_${RUN_ID}"

export PYTHONPATH="$REPO_ROOT/fresh/src:${PYTHONPATH:-}"

python -m egfr_myo1d.cli init-run --mode smoke_env --run-id "$RUN_ID" >/dev/null

PBS_DIR="$REPO_ROOT/fresh/runs/$RUN_ID/scripts"
mkdir -p "$PBS_DIR"
PBS_FILE="$PBS_DIR/${JOB_NAME}.pbs"

cat > "$PBS_FILE" <<PBS
#!/bin/bash
#PBS -q workq
#PBS -N $JOB_NAME
#PBS -l nodes=$NODE:ppn=4
#PBS -l walltime=00:20:00
#PBS -V
#PBS -o $REPO_ROOT/fresh/runs/$RUN_ID/logs/jobs/$JOB_NAME.stdout
#PBS -e $REPO_ROOT/fresh/runs/$RUN_ID/logs/jobs/$JOB_NAME.stderr

set -euo pipefail
cd "$REPO_ROOT"

source /usr/local/anaconda/3/2023.09/etc/profile.d/conda.sh
conda activate pyrosetta

export LANG=C.UTF-8
export LC_ALL=C.UTF-8
export PYTHONIOENCODING=utf-8
export PYTHONPATH="$REPO_ROOT/fresh/src:\${PYTHONPATH:-}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

python -m egfr_myo1d.cli tool-preflight --run-id "$RUN_ID" --mode "$MODE" --profile hpc_strict
python -m egfr_myo1d.cli status --run-id "$RUN_ID"
PBS

chmod +x "$PBS_FILE"

cat <<MSG
Tool-preflight PBS file generated:

    $PBS_FILE

Submit on HPC from repo root:

    qsub $PBS_FILE

After completion, inspect:

    fresh/runs/$RUN_ID/manifest/tool_status.json
    fresh/runs/$RUN_ID/reports/tool_installation_report.md
MSG
