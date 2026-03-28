#!/usr/bin/env bash
# Internal helper script.
# Canonical server submission entrypoint: qsub config/run_precision_smoke.pbs
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BASE_CONFIG="${ROOT_DIR}/config/example-project.yaml"
SMOKE_CONFIG="${ROOT_DIR}/config/smoke-fast.yaml"
BACKUP_CONFIG="$(mktemp)"

cleanup() {
  cp "${BACKUP_CONFIG}" "${BASE_CONFIG}" || true
  rm -f "${BACKUP_CONFIG}" || true
}
trap cleanup EXIT

cp "${BASE_CONFIG}" "${BACKUP_CONFIG}"
cp "${SMOKE_CONFIG}" "${BASE_CONFIG}"

echo "[precision-smoke] Using smoke profile: ${SMOKE_CONFIG}"

echo "[1/5] lane=vina-cpu"
python "${ROOT_DIR}/run_production.py" --lane vina-cpu --allocated-cpus 1

echo "[2/5] lane=ppi (single state/seed)"
python "${ROOT_DIR}/run_production.py" --lane ppi --state 3GT8_raw --seed 0 --allocated-cpus 1

echo "[3/5] lane=ppi-post"
python "${ROOT_DIR}/run_production.py" --lane ppi-post --allocated-cpus 1

echo "[4/5] lane=vina-post"
python "${ROOT_DIR}/run_production.py" --lane vina-post --allocated-cpus 1

echo "[5/5] lane=finalize"
python "${ROOT_DIR}/run_production.py" --lane finalize --allocated-cpus 1

echo "[precision-smoke] Completed. Outputs under ${ROOT_DIR}/output_smoke"
