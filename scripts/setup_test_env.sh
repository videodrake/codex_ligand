#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TEST_CONDA_ENV="${TEST_CONDA_ENV:-codex-tests}"
TEST_PYTHON_VERSION="${TEST_PYTHON_VERSION:-3.11}"

if ! command -v conda >/dev/null 2>&1; then
  echo "conda command not found."
  echo "Load your conda setup first, then rerun this script."
  exit 1
fi

eval "$(conda shell.bash hook)"

if conda env list | awk '{print $1}' | grep -Fxq "${TEST_CONDA_ENV}"; then
  echo "Updating existing conda test environment: ${TEST_CONDA_ENV}"
else
  echo "Creating conda test environment: ${TEST_CONDA_ENV}"
  conda create -y -n "${TEST_CONDA_ENV}" "python=${TEST_PYTHON_VERSION}"
fi

mapfile -t TEST_PACKAGES < "${ROOT_DIR}/requirements-test.txt"
conda install -y -n "${TEST_CONDA_ENV}" "${TEST_PACKAGES[@]}"

echo
echo "Test environment is ready."
echo "Activate with:"
echo "  conda activate ${TEST_CONDA_ENV}"
echo "Run pre-qsub checks with:"
echo "  ${ROOT_DIR}/scripts/run_pre_qsub_checks.sh"
