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

echo "[3/4] Pytest phase and smoke suite"
python -m pytest \
  -m "not reporting" \
  "${ROOT_DIR}/tests/test_phase2.py" \
  "${ROOT_DIR}/tests/test_phase3.py" \
  "${ROOT_DIR}/tests/test_phase4.py" \
  "${ROOT_DIR}/tests/test_cluster_consensus.py" \
  "${ROOT_DIR}/tests/test_compare_states.py" \
  "${ROOT_DIR}/tests/test_phase1_smoke.py" \
  "${ROOT_DIR}/tests/test_precheck_guard.py" \
  "${ROOT_DIR}/tests/test_review_report.py" \
  "${ROOT_DIR}/tests/test_lightdock_validation.py" \
  "${ROOT_DIR}/tests/test_pyrosetta_extract.py" \
  "${ROOT_DIR}/tests/test_pyrosetta_metadata.py" \
  "${ROOT_DIR}/tests/test_postprocess_ppi.py" \
  "${ROOT_DIR}/tests/test_smoke_cli.py" \
  "${ROOT_DIR}/tests/test_validation_smoke.py" \
  -q

echo "[4/4] Pre-qsub checks passed"
