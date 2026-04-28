#!/usr/bin/env python3
"""Milestone 1 placeholder for future PBS generation."""

import argparse
import textwrap


PBS_PREVIEW = """\
#!/bin/bash
#PBS -q workq
#PBS -l nodes=node04:ppn=4
#PBS -l walltime=02:00:00

export REPO_ROOT="${REPO_ROOT:-$(pwd)}"
export PYTHONPATH="$REPO_ROOT/fresh/src:${PYTHONPATH:-}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

source /usr/local/anaconda/3/2023.09/etc/profile.d/conda.sh
conda activate pyrosetta

python -m egfr_myo1d.cli version
"""


def build_parser():
    parser = argparse.ArgumentParser(
        description="Preview the future fresh workflow PBS skeleton without writing files."
    )
    parser.add_argument(
        "--preview",
        action="store_true",
        help="Print a non-submitting PBS preview snippet.",
    )
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.preview:
        print(PBS_PREVIEW)
    else:
        print(
            textwrap.dedent(
                """\
                Milestone 1 Task 1 placeholder: no PBS files are generated yet.
                Use --preview to display the required non-submitting PBS principles.
                """
            ).strip()
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
