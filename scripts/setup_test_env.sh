#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONDA_ENV_NAME="pyrosetta"

if ! command -v conda >/dev/null 2>&1; then
  echo "conda command not found."
  echo "Load your conda setup first, then rerun this script."
  exit 1
fi

eval "$(conda shell.bash hook)"
conda activate "${CONDA_ENV_NAME}"

python - <<'PY'
import importlib
import sys

required = {
    "pytest": "pytest",
    "pyyaml": "yaml",
    "numpy": "numpy",
    "pandas": "pandas",
    "scipy": "scipy",
    "matplotlib": "matplotlib",
}

missing = []
for package_name, module_name in required.items():
    try:
        importlib.import_module(module_name)
    except Exception:
        missing.append(package_name)

if missing:
    print("Missing required packages in the 'pyrosetta' conda environment:")
    for name in missing:
        print(f"  - {name}")
    print("")
    print("Install them manually before running pre-qsub checks.")
    print("See: docs/server_environment_setup.md")
    sys.exit(1)
PY

echo
echo "Pre-qsub environment check passed."
echo "Active conda environment: ${CONDA_ENV_NAME}"
echo "Run pre-qsub checks with:"
echo "  ${ROOT_DIR}/scripts/run_pre_qsub_checks.sh"
