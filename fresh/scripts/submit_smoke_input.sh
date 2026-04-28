#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="${REPO_ROOT:-$(pwd)}"

cat <<'MSG'
Milestone 1 Task 1 placeholder: this script does not call qsub.

Future smoke-input PBS jobs should include:

#PBS -q workq
#PBS -l nodes=node04:ppn=4
#PBS -l walltime=02:00:00

export PYTHONPATH="$REPO_ROOT/fresh/src:${PYTHONPATH:-}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

source /usr/local/anaconda/3/2023.09/etc/profile.d/conda.sh
conda activate pyrosetta

Real input preparation is intentionally not implemented in Task 1.
No qsub command is run by this placeholder.
MSG

export PYTHONPATH="$REPO_ROOT/fresh/src:${PYTHONPATH:-}"
python -m egfr_myo1d.cli version
