#!/usr/bin/env python3
"""Phase 3 Cascade Runner — sequential execution of TG 3.0 through 3.6.

Phase 3 has three execution modes because TG 3.3 (Vina docking) requires
the Vina binary on the compute node and runs in multiple rounds:

  setup   : TG 3.0 → 3.2 + docking script generation (no Vina needed)
  execute : TG 3.3 round-by-round Vina execution (Vina required)
  post    : TG 3.4 → 3.6 post-docking analysis (no Vina needed)

Usage:
  # Full cascade (setup + execute all rounds + post) — requires Vina
  python -m egfr_pipeline.phase3.rerun_cascade

  # Setup only (planning, no Vina needed)
  python -m egfr_pipeline.phase3.rerun_cascade --mode setup

  # Execute a specific round
  python -m egfr_pipeline.phase3.rerun_cascade --mode execute --round 0

  # Post-docking analysis only
  python -m egfr_pipeline.phase3.rerun_cascade --mode post

  # Resume from a specific TG
  python -m egfr_pipeline.phase3.rerun_cascade --from-tg 3.4
"""

import argparse
import json
from pathlib import Path
from typing import Optional

from egfr_pipeline import paths
from egfr_pipeline.cascade_runner import run_cascade

_CFG = {"output_root": str(paths.REPO_ROOT / "output")}


def run_phase3_cascade(
    mode: str = "full",
    from_tg: str = "3.0",
    round_id: Optional[int] = None,
    allocated_cpus: Optional[int] = None,
) -> None:
    """Run Phase 3 TG 3.0 -> 3.6 sequentially.

    Parameters
    ----------
    mode : str
        "full" (all TGs), "setup" (3.0-3.2 + script gen), "execute" (3.3 per round),
        "post" (3.4-3.6).
    from_tg : str
        Start from this TG (e.g. "3.4" to skip earlier steps).
    round_id : int, optional
        When mode="execute", run only this round.
    allocated_cpus : int, optional
        Scheduler-allocated CPUs for worker cap.
    """
    output_dir = paths.wb_phase3_focused_docking(_CFG)

    # Handoff file validation
    handoff = paths.wb_phase2_pocket_analysis(_CFG) / "phase3_candidate_pocket_reference.csv"
    if not handoff.exists():
        raise FileNotFoundError(
            f"Phase 2 handoff file not found: {handoff}\n"
            "Phase 2 (Pocket Analysis) must complete before Phase 3.\n"
            "Run: qsub config/run_phase2_cascade.pbs"
        )

    setup_steps = [
        ("3.0", "Pocket Reference Ingestion", lambda: _run_tg30(output_dir)),
        ("3.1", "Job Construction", lambda: _run_tg31(output_dir)),
        ("3.2", "Budget Policy", lambda: _run_tg32(output_dir)),
    ]

    setup_gen_step = (
        "3.2s", "Docking Script Generation",
        lambda: _run_tg33_setup(output_dir, allocated_cpus),
    )

    post_steps = [
        ("3.4", "Pose Attribution", lambda: _run_tg34(output_dir)),
        ("3.5", "Diversity Validation", lambda: _run_tg35(output_dir)),
        ("3.6", "Phase 4 Export", lambda: _run_tg36(output_dir)),
    ]

    review_step = ("3.7", "Review Report", lambda: _run_tg37(output_dir))

    # Build step list based on mode
    if mode == "setup":
        steps = setup_steps + [setup_gen_step]
    elif mode == "execute":
        # setup 완료 검증
        job_table = output_dir / "phase3_docking_job_table.csv"
        if not job_table.exists():
            raise FileNotFoundError(
                f"Phase 3 job table 없음: {job_table}\n"
                "먼저 --mode setup 을 실행하세요."
            )
        steps = [
            ("3.3", f"Docking Execution (round {round_id})",
             lambda: _run_tg33_execute(output_dir, round_id, allocated_cpus)),
        ]
    elif mode == "post":
        # execute 완료 검증 (round log에 데이터 행이 있어야 함)
        round_log = output_dir / "phase3_round_log.csv"
        if not _has_data_rows(round_log):
            raise FileNotFoundError(
                f"Phase 3 round log 없음 또는 비어있음: {round_log}\n"
                "최소 1회의 execute round가 완료되어야 합니다.\n"
                "--mode execute --round 0 을 먼저 실행하세요."
            )
        steps = post_steps + [review_step]
    else:  # full
        steps = setup_steps + [setup_gen_step]
        # For full mode, read round count and execute all rounds
        steps.append(
            ("3.3", "Docking Execution (all rounds)",
             lambda: _run_tg33_all_rounds(output_dir, allocated_cpus)),
        )
        steps += post_steps + [review_step]

    header_lines = [
        f"Output:   {output_dir}",
        f"Mode:     {mode}",
    ]
    if round_id is not None:
        header_lines.append(f"Round:    {round_id}")

    run_cascade(
        phase_label="Phase 3 Cascade Runner",
        steps=steps,
        from_tg=from_tg,
        header_lines=header_lines,
    )


# ---------------------------------------------------------------------------
# Per-TG wrappers
# ---------------------------------------------------------------------------

def _run_tg30(output_dir: Path) -> None:
    from egfr_pipeline.phase3.pocket_reference_ingestion import run_pocket_reference_ingestion
    run_pocket_reference_ingestion(output_dir=output_dir)


def _run_tg31(output_dir: Path) -> None:
    from egfr_pipeline.phase3.job_construction import run_job_construction
    run_job_construction(output_dir=output_dir)


def _run_tg32(output_dir: Path) -> None:
    from egfr_pipeline.phase3.budget_policy import run_budget_policy
    run_budget_policy(output_dir=output_dir)


def _run_tg33_setup(output_dir: Path, allocated_cpus: Optional[int]) -> None:
    from egfr_pipeline.phase3.run_diverse_docking import run_setup
    kwargs = {"output_dir": output_dir}
    if allocated_cpus:
        kwargs["allocated_cpus"] = allocated_cpus
    run_setup(**kwargs)


def _run_tg33_execute(
    output_dir: Path,
    round_id: Optional[int],
    allocated_cpus: Optional[int],
) -> None:
    from egfr_pipeline.phase3.run_diverse_docking import execute_round
    if round_id is None:
        raise ValueError("--round is required when mode=execute")
    jobs = _load_round_jobs(output_dir, round_id)
    kwargs = {"round_id": round_id}
    if allocated_cpus:
        kwargs["allocated_cpus"] = allocated_cpus
    execute_round(jobs, **kwargs)


def _run_tg33_all_rounds(output_dir: Path, allocated_cpus: Optional[int]) -> None:
    """Execute all rounds sequentially (for full mode on server)."""
    round_count = _get_round_count(output_dir)
    print(f"  Executing {round_count} rounds sequentially...")
    for r in range(round_count):
        print(f"\n  --- Round {r}/{round_count - 1} ---")
        _run_tg33_execute(output_dir, r, allocated_cpus)


def _run_tg34(output_dir: Path) -> None:
    from egfr_pipeline.phase3.pose_attribution import run_pose_attribution
    run_pose_attribution(output_dir=output_dir)


def _run_tg35(output_dir: Path) -> None:
    from egfr_pipeline.phase3.diversity_validation import run_diversity_validation
    run_diversity_validation(output_dir=output_dir)


def _run_tg36(output_dir: Path) -> None:
    from egfr_pipeline.phase3.phase4_export import run_phase4_export
    run_phase4_export(output_dir=output_dir)


def _run_tg37(output_dir: Path) -> None:
    from egfr_pipeline.phase3.review_report import run_review_report
    run_review_report(output_dir=output_dir)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _has_data_rows(path: Path) -> bool:
    """CSV 파일이 존재하고 데이터 행이 1개 이상인지 확인."""
    if not path.exists():
        return False
    import csv
    try:
        with open(path, newline="", encoding="utf-8") as f:
            reader = csv.reader(f)
            next(reader, None)  # header
            return next(reader, None) is not None
    except Exception:
        return False


def _get_round_count(output_dir: Path) -> int:
    """Read round count from budget policy output, default 3."""
    budget_file = output_dir / "phase3_budget_policy.json"
    if budget_file.exists():
        try:
            data = json.loads(budget_file.read_text())
            return int(data.get("max_rounds", 3))
        except (json.JSONDecodeError, KeyError, ValueError):
            pass
    return 3


def _load_round_jobs(output_dir: Path, round_id: int) -> list:
    """Load job list for a specific round from the job table."""
    import csv
    job_table = output_dir / "phase3_docking_job_table.csv"
    if not job_table.exists():
        raise FileNotFoundError(
            f"Job table not found: {job_table}\n"
            "Run Phase 3 setup first (--mode setup)."
        )
    jobs = []
    with open(job_table) as f:
        reader = csv.DictReader(f)
        for row in reader:
            if int(row.get("round_id", 0)) == round_id:
                jobs.append(row)
    if not jobs:
        print(f"  [WARN] No jobs found for round {round_id}")
    return jobs


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Phase 3 Cascade Runner (TG 3.0 -> 3.6)",
    )
    parser.add_argument(
        "--mode", choices=["full", "setup", "execute", "post"], default="full",
        help="Execution mode: full, setup (3.0-3.2), execute (3.3 per round), post (3.4-3.6)",
    )
    parser.add_argument(
        "--from-tg", default="3.0",
        help="Start from this TG (e.g. '3.4' to skip earlier steps)",
    )
    parser.add_argument(
        "--round", type=int, default=None, dest="round_id",
        help="Round ID for execute mode",
    )
    parser.add_argument(
        "--allocated-cpus", type=int, default=None,
        help="Scheduler-allocated CPUs for worker cap",
    )
    args = parser.parse_args()

    run_phase3_cascade(
        mode=args.mode,
        from_tg=args.from_tg,
        round_id=args.round_id,
        allocated_cpus=args.allocated_cpus,
    )


if __name__ == "__main__":
    main()
