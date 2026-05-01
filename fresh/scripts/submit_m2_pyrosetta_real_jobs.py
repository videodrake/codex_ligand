#!/usr/bin/env python
"""Generate M2.2 real PyRosetta PBS scripts for HPC.

By default this script only writes PBS files and manifests. It calls qsub only
when both --submit and --i-understand-this-submits-hpc-jobs are provided.
"""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "fresh" / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from egfr_myo1d.core.logging_utils import initialize_logs  # noqa: E402
from egfr_myo1d.core.run_context import RunContext  # noqa: E402
from egfr_myo1d.ppi.pyrosetta_real_jobs import plan_pyrosetta_real_pbs  # noqa: E402


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Generate M2.2 real PyRosetta chunked PBS scripts; qsub is disabled unless explicitly confirmed."
    )
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--profile", default="production", choices=["smoke", "mini", "scaling", "production"])
    parser.add_argument("--job-manifest", default=None)
    parser.add_argument("--nodes", default="node04,node05,node06")
    parser.add_argument("--ppn", type=int, default=None)
    parser.add_argument("--walltime", default=None)
    parser.add_argument("--queue", default=None)
    parser.add_argument("--conda-env", default=None)
    parser.add_argument("--conda-sh", default=None)
    parser.add_argument("--python-executable", default=None)
    parser.add_argument(
        "--models-per-chunk",
        type=int,
        default=None,
        help="Models executed by each isolated subprocess chunk. Default auto-splits each state/seed into about ppn chunks.",
    )
    parser.add_argument("--submit", action="store_true", default=False)
    parser.add_argument("--i-understand-this-submits-hpc-jobs", action="store_true", default=False)
    return parser


def _submit_pbs_files(pbs_manifest: Path) -> int:
    with pbs_manifest.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    failures = 0
    for row in rows:
        if row.get("pbs_generation_status") != "PASS" or not row.get("qsub_command"):
            continue
        pbs_file = REPO_ROOT / row["pbs_file"]
        proc = subprocess.run(["qsub", str(pbs_file)], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False)
        if proc.returncode != 0:
            failures += 1
            sys.stderr.write(proc.stderr)
    return failures


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    ctx = RunContext.for_existing(args.run_id, repo_root=REPO_ROOT)
    initialize_logs(ctx)
    nodes = [item.strip() for item in args.nodes.split(",") if item.strip()]
    plan = plan_pyrosetta_real_pbs(
        ctx,
        profile=args.profile,
        job_manifest=Path(args.job_manifest) if args.job_manifest else None,
        nodes=nodes,
        ppn=args.ppn,
        walltime=args.walltime,
        queue=args.queue,
        conda_env=args.conda_env,
        conda_sh=args.conda_sh,
        python_executable=args.python_executable,
        models_per_chunk=args.models_per_chunk,
    )
    print(
        "m2 PyRosetta real PBS plan {0}: chunks={1} models={2} pbs_files={3}".format(
            plan.status,
            plan.counts.get("planned_chunks", 0),
            plan.counts.get("planned_models", 0),
            plan.counts.get("pbs_files", 0),
        )
    )
    print("chunk_manifest={0}".format(plan.chunk_manifest_csv))
    print("pbs_manifest={0}".format(plan.pbs_manifest_csv))
    if not args.submit:
        print("No qsub submission requested. To submit on HPC, pass both --submit and --i-understand-this-submits-hpc-jobs.")
        return 0 if plan.status != "FAIL" else 1
    if not args.i_understand_this_submits_hpc_jobs:
        print("Refusing qsub: --submit requires --i-understand-this-submits-hpc-jobs.")
        return 2
    if plan.status == "FAIL":
        print("Refusing qsub: generated plan failed QC.")
        return 2
    failures = _submit_pbs_files(plan.pbs_manifest_csv)
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
