#!/usr/bin/env bash
# Internal helper script.
# Canonical server submission entrypoint: qsub config/run_pre_qsub_checks.pbs
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

if ! command -v conda >/dev/null 2>&1; then
  echo "conda command not found."
  echo "Load your conda setup first, then rerun this script."
  exit 1
fi

eval "$(conda shell.bash hook)"
bash "${ROOT_DIR}/scripts/setup_test_env.sh"
conda activate pyrosetta

export PYTHONUTF8=1

echo "[1/4] Syntax compile check"
python -m compileall "${ROOT_DIR}/main.py" "${ROOT_DIR}/egfr_pipeline" "${ROOT_DIR}/tests"

echo "[2/4] CLI smoke check"
python -X utf8 "${ROOT_DIR}/main.py" --help > /dev/null
python -X utf8 "${ROOT_DIR}/main.py" validate --help > /dev/null

# AC-6.2: Workflow B external tool availability check
if [ "${CHECK_WORKFLOW_B:-0}" = "1" ]; then
  echo "[2.5/4] Workflow B tool availability check"
  TOOLS_OK=1

  # fpocket
  if command -v fpocket >/dev/null 2>&1; then
    echo "  fpocket: $(fpocket --version 2>&1 | head -1 || echo 'available')"
  else
    echo "  WARNING: fpocket not found. Phase 2 pocket proposal will fail."
    TOOLS_OK=0
  fi

  # P2Rank
  if command -v prank >/dev/null 2>&1 || [ -x "${PRANK_HOME:-}/prank" ]; then
    echo "  P2Rank: available"
  else
    echo "  WARNING: P2Rank (prank) not found. Phase 2 will use fpocket only."
  fi

  # LightDock
  if command -v lightdock3 >/dev/null 2>&1 || command -v lightdock3_setup.py >/dev/null 2>&1; then
    echo "  LightDock: available"
  else
    echo "  WARNING: LightDock not found. Phase 1 cross-method validation will be skipped."
  fi

  if [ "$TOOLS_OK" -eq 0 ]; then
    echo "  CRITICAL: fpocket missing — Workflow B cannot proceed without it."
    exit 1
  fi
fi

echo "[3/4] Pytest baseline suite"
python -m pytest \
  -m "unit and not slow" \
  "${ROOT_DIR}/tests/unit" \
  -q

echo "[4/4] Pre-qsub checks passed"
