#!/usr/bin/env python
"""Safe M3-T5 helper for generating compound docking PBS plans.

By default this script only generates/prints the plan. It calls qsub only when
both --submit and --i-understand-this-submits-hpc-jobs are provided.
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

from egfr_myo1d.compound.vina_jobs import run_m3_vina_jobs  # noqa: E402
from egfr_myo1d.core.logging_utils import initialize_logs  # noqa: E402
from egfr_myo1d.core.run_context import RunContext  # noqa: E402


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Generate M3 compound docking PBS scripts; qsub is disabled unless explicitly confirmed.")
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--m2-run-id", required=True)
    parser.add_argument("--mode", default="generate", choices=["dry-run", "generate"])
    parser.add_argument("--profile", default="mini", choices=["smoke", "mini", "scaling", "production"])
    parser.add_argument("--compound-public-id", choices=["Cpd-A", "Cpd-B", "Cpd-C"], default=None)
    parser.add_argument("--state-id", default=None)
    parser.add_argument("--pocket-family-id", default=None)
    parser.add_argument("--box-id", default=None)
    parser.add_argument("--protomer-id", default=None)
    parser.add_argument("--vina-executable", default="vina")
    parser.add_argument("--nodes", default="node04,node05,node06")
    parser.add_argument("--ppn", type=int, default=None)
    parser.add_argument("--walltime", default=None)
    parser.add_argument("--queue", default=None)
    parser.add_argument("--conda-env", default=None)
    parser.add_argument("--conda-sh", default=None)
    parser.add_argument("--python-executable", default=None)
    parser.add_argument("--allow-independent", action="store_true")
    parser.add_argument("--allow-production-plan-without-mini", action="store_true")
    parser.add_argument("--submit", action="store_true", default=False, help="Submit generated PBS files only with confirmation flag.")
    parser.add_argument("--i-understand-this-submits-hpc-jobs", action="store_true", default=False)
    return parser


def _submit_pbs_files(pbs_manifest: Path) -> int:
    with pbs_manifest.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    failures = 0
    for row in rows:
        if row.get("pbs_generation_status") != "PASS" or not row.get("qsub_command") or int(row.get("assigned_job_count") or 0) <= 0:
            continue
        pbs_file = REPO_ROOT / row["pbs_file"]
        proc = subprocess.run(["qsub", str(pbs_file)], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False)
        if proc.returncode != 0:
            failures += 1
    return failures


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    ctx = RunContext.create(args.run_id, "smoke_input", repo_root=REPO_ROOT)
    initialize_logs(ctx)
    result = run_m3_vina_jobs(
        ctx,
        m2_run_id=args.m2_run_id,
        mode=args.mode,
        profile=args.profile,
        compound_public_id=args.compound_public_id,
        state_id=args.state_id,
        pocket_family_id=args.pocket_family_id,
        box_id=args.box_id,
        protomer_id=args.protomer_id,
        vina_executable=args.vina_executable,
        nodes=[item.strip() for item in args.nodes.split(",") if item.strip()],
        ppn=args.ppn,
        walltime=args.walltime,
        queue=args.queue,
        conda_env=args.conda_env,
        conda_sh=args.conda_sh,
        python_executable=args.python_executable,
        allow_independent=args.allow_independent,
        allow_production_plan_without_mini=args.allow_production_plan_without_mini,
    )
    print("m3 compound docking PBS plan {0}: jobs={1} pbs_files={2} qsub_submission_allowed={3}".format(result.status, result.counts.get("planned_vina_commands", 0), result.counts.get("pbs_files", 0), result.qsub_submission_allowed))
    print("pbs_manifest={0}".format(result.pbs_manifest_csv))
    if not args.submit:
        print("No qsub submission requested. To submit on HPC, pass both --submit and --i-understand-this-submits-hpc-jobs.")
        return 0 if result.status != "FAIL" else 1
    if not args.i_understand_this_submits_hpc_jobs:
        print("Refusing qsub: --submit requires --i-understand-this-submits-hpc-jobs.")
        return 2
    if not result.qsub_submission_allowed:
        print("Refusing qsub: generated plan did not open qsub_submission_allowed.")
        return 2
    failures = _submit_pbs_files(result.pbs_manifest_csv)
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
