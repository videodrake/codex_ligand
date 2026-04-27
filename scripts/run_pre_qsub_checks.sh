#!/usr/bin/env bash
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

echo "[3/4] Pytest smoke suite"
python -m pytest \
  -m "smoke and not reporting" \
  --strict-markers \
  "${ROOT_DIR}/tests" \
  -q

echo "[4/4] Pre-qsub checks passed"
