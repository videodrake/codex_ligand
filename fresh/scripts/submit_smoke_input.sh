#!/usr/bin/env bash
# M1 Phase 6: emit a smoke-input PBS file under fresh/runs/<run_id>/scripts/
# and print the qsub command for the user to run manually on HPC.
#
# Smoke-input mode runs `prepare-inputs` (Phase 8 orchestrator) which expects
# real EGFR/MYO1D inputs under fresh/data/raw/. If real files are absent the
# CLI reports `missing_required_inputs` cleanly without crashing.
#
# Does NOT auto-call qsub.
#
# Usage:
#   bash fresh/scripts/submit_smoke_input.sh [<run_id>] [<node>]
set -euo pipefail

REPO_ROOT="${REPO_ROOT:-$(pwd)}"
RUN_ID="${1:-smoke_input_$(date +%Y%m%d_%H%M%S)}"
NODE="${2:-node04}"
JOB_NAME="smoke_input_${RUN_ID}"

export PYTHONPATH="$REPO_ROOT/fresh/src:${PYTHONPATH:-}"

python -m egfr_myo1d.cli init-run --mode smoke_input --run-id "$RUN_ID" >/dev/null

PBS_FILE="$REPO_ROOT/fresh/runs/$RUN_ID/scripts/${JOB_NAME}.pbs"
python -m egfr_myo1d.cli prepare-pbs \
    --run-id "$RUN_ID" \
    --job-name "$JOB_NAME" \
    --mode smoke_input \
    --node "$NODE" \
    --input-root fresh/data/raw

cat <<MSG
Smoke-input PBS file generated:

    $PBS_FILE

To submit on HPC, run from repo root:

    qsub $PBS_FILE

This script does NOT auto-call qsub. Inspect the generated file before submission.
Place real EGFR/MYO1D/ligand input PDBs under fresh/data/raw/ before running on HPC.
MSG
