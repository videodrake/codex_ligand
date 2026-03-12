#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
VENV_DIR="${ROOT_DIR}/.venv-tests"

python3 -m venv "${VENV_DIR}"
source "${VENV_DIR}/bin/activate"

python -m pip install --upgrade pip setuptools wheel
python -m pip install -r "${ROOT_DIR}/requirements-test.txt"

echo
echo "Test environment is ready."
echo "Activate with:"
echo "  source ${VENV_DIR}/bin/activate"
echo "Run pre-qsub checks with:"
echo "  ${ROOT_DIR}/scripts/run_pre_qsub_checks.sh"
