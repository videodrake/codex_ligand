#!/usr/bin/env bash
# Run PPI-surface tool discovery/smoke checks in the current shell.
# This script does not install packages or submit scheduler jobs.
#
# Usage:
#   bash fresh/scripts/preflight_tools.sh [<run_id>] [discover|smoke] [codex_dev|hpc_strict]
set -euo pipefail

REPO_ROOT="${REPO_ROOT:-$(pwd)}"
RUN_ID="${1:-tool_preflight_$(date +%Y%m%d_%H%M%S)}"
MODE="${2:-discover}"
PROFILE="${3:-codex_dev}"

export PYTHONPATH="$REPO_ROOT/fresh/src:${PYTHONPATH:-}"

python -m egfr_myo1d.cli init-run --mode smoke_env --run-id "$RUN_ID" >/dev/null
python -m egfr_myo1d.cli tool-preflight --run-id "$RUN_ID" --mode "$MODE" --profile "$PROFILE"

cat <<MSG

Tool preflight finished.

Inspect:
  fresh/runs/$RUN_ID/manifest/tool_status.json
  fresh/runs/$RUN_ID/reports/tool_installation_report.md
MSG
