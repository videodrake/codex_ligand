#!/usr/bin/env bash
# M1 Phase 6: emit a smoke-env PBS file under fresh/runs/<run_id>/scripts/
# and print the qsub command for the user to run manually on HPC.
#
# Does NOT auto-call qsub.
#
# Usage:
#   bash fresh/scripts/submit_smoke_env.sh [<run_id>] [<node>]
set -euo pipefail

REPO_ROOT="${REPO_ROOT:-$(pwd)}"
RUN_ID="${1:-smoke_env_$(date +%Y%m%d_%H%M%S)}"
NODE="${2:-node04}"
JOB_NAME="smoke_env_${RUN_ID}"

export PYTHONPATH="$REPO_ROOT/fresh/src:${PYTHONPATH:-}"

# 1. Create or reuse the run directory
python -m egfr_myo1d.cli init-run --mode smoke_env --run-id "$RUN_ID" >/dev/null

# 2. Generate the concrete PBS file under runs/<run_id>/scripts/
PBS_FILE="$REPO_ROOT/fresh/runs/$RUN_ID/scripts/${JOB_NAME}.pbs"
python -m egfr_myo1d.cli prepare-pbs \
    --run-id "$RUN_ID" \
    --job-name "$JOB_NAME" \
    --mode smoke_env \
    --node "$NODE"

cat <<MSG
Smoke-env PBS file generated:

    $PBS_FILE

To submit on HPC, run from repo root:

    qsub $PBS_FILE

This script does NOT auto-call qsub. Inspect the generated file before submission.
MSG
