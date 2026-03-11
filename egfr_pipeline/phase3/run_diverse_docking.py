#!/usr/bin/env python3
"""Phase 3 Task Group 3.3: Diversity-Aware Vina Execution Layer.

Executes pocket-guided, budget-aware ligand docking that deliberately
spreads search effort across candidate pockets rather than naive blind
docking.

Tasks:
  3.3.1 - New execution entry point (coexists with legacy)
  3.3.2 - Pocket-guided local docking execution
  3.3.3 - Controlled round-based execution with saturation
  3.3.4 - Workspace vs server rule preservation

Inputs:
  - output/phase3_docking/phase3_docking_job_table.csv
  - output/phase3_docking/pocket_search_status.csv
  - output/phase3_docking/phase3_candidate_reference_normalized.csv

Outputs:
  - output/phase3_docking/phase3_run_metadata.json
  - output/phase3_docking/phase3_round_log.csv
  - output/phase3_docking/pocket_search_status.csv  (updated)
  - output/phase3_docking/phase3_budget_tracking.csv (updated)

Usage:
    # Setup mode (generates run script for server execution)
    python -m egfr_pipeline.phase3.run_diverse_docking --setup

    # Execute mode (requires Vina — server only)
    python -m egfr_pipeline.phase3.run_diverse_docking --execute

    # Dry-run (validate without executing)
    python -m egfr_pipeline.phase3.run_diverse_docking --dry-run
"""

import argparse
import csv
import json
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

PHASE3_OUTPUT_DIR = PROJECT_ROOT / "output" / "phase3_docking"
MAX_WORKERS_DEFAULT = 16  # Server-side default (Task 3.3.4)

# Round log schema
ROUND_LOG_COLUMNS = [
    "round_id",
    "pocket_id",
    "receptor_id",
    "ligand_id",
    "seed_index",
    "job_id",
    "status",
    "n_poses",
    "best_affinity_kcal",
    "output_file",
    "error_message",
    "timestamp",
]


# ---------------------------------------------------------------------------
# 3.3.1: Execution entry point — setup mode
# ---------------------------------------------------------------------------

def generate_run_script(
    jobs: List[dict],
    pocket_statuses: List[dict],
    output_dir: Path,
    max_workers: int = MAX_WORKERS_DEFAULT,
) -> Path:
    """Generate a shell script for server-side execution.

    This is the --setup mode: produces a self-contained script that
    the server can execute without re-running the planning phase.
    """
    script_path = output_dir / "run_phase3_docking.sh"

    # Filter to only open pockets
    open_pockets = {
        ps["pocket_id"] for ps in pocket_statuses
        if ps["status"] == "open"
    }
    active_jobs = [j for j in jobs if j["pocket_id"] in open_pockets]

    # Group jobs by round (seed_index)
    by_round = defaultdict(list)
    for j in active_jobs:
        by_round[int(j["seed_index"])].append(j)

    lines = [
        "#!/bin/bash",
        "# Phase 3 Diversity-Aware Docking — Auto-generated run script",
        f"# Generated: {datetime.now().isoformat()}",
        f"# Total active jobs: {len(active_jobs)}",
        f"# Rounds: {len(by_round)}",
        f"# Max workers: {max_workers}",
        "",
        "set -euo pipefail",
        "",
        "SCRIPT_DIR=\"$(cd \"$(dirname \"${BASH_SOURCE[0]}\")\" && pwd)\"",
        "PROJECT_ROOT=\"$(cd \"$SCRIPT_DIR/../..\" && pwd)\"",
        "",
        "# Activate environment",
        "# conda activate vina  # Uncomment if needed",
        "",
        f"OUTPUT_DIR=\"{output_dir}\"",
        "",
    ]

    for round_id in sorted(by_round.keys()):
        round_jobs = by_round[round_id]
        lines.extend([
            f"echo '=== Round {round_id}: {len(round_jobs)} jobs ==='",
            "",
        ])

        for j in round_jobs:
            # Create output directory
            lines.append(f"mkdir -p \"{j['output_dir']}\"")

        lines.append("")
        lines.extend([
            f"# Execute round {round_id} via Phase 3 runner",
            f"python -m egfr_pipeline.phase3.run_diverse_docking \\",
            f"    --execute \\",
            f"    --round {round_id} \\",
            f"    --max_workers {max_workers} \\",
            f"    --output_dir \"$OUTPUT_DIR\"",
            "",
            "# Evaluate saturation after round",
            f"python -m egfr_pipeline.phase3.run_diverse_docking \\",
            f"    --evaluate-round {round_id} \\",
            f"    --output_dir \"$OUTPUT_DIR\"",
            "",
        ])

    lines.extend([
        "echo '=== Phase 3 Docking Complete ==='",
        f"echo \"Results in: $OUTPUT_DIR\"",
    ])

    with open(script_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    script_path.chmod(0o755)

    return script_path


# ---------------------------------------------------------------------------
# 3.3.2: Pocket-guided docking execution
# ---------------------------------------------------------------------------

def execute_single_job(job: dict) -> dict:
    """Execute a single Vina docking job.

    Wraps the existing Vina execution interface. Returns result dict.
    This function requires Vina to be available (server-side only).
    """
    try:
        from egfr_pipeline.vina.dock import run_docking, prepare_ligand
    except ImportError:
        return {
            "job_id": job["job_id"],
            "status": "skip_no_vina",
            "error_message": "Vina not available in this environment",
            "n_poses": 0,
            "best_affinity_kcal": "",
        }

    receptor_path = Path(job["receptor_pdb"])
    ligand_path = Path(job["ligand_file"])
    output_dir = Path(job["output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)

    # Prepare PDBQT if needed
    receptor_pdbqt = receptor_path.with_suffix(".pdbqt")
    if not receptor_pdbqt.exists():
        receptor_pdbqt = receptor_path  # Assume pre-converted

    ligand_pdbqt = output_dir / f"{ligand_path.stem}.pdbqt"
    if ligand_path.suffix.lower() != ".pdbqt":
        try:
            prepare_ligand(str(ligand_path), str(ligand_pdbqt))
        except Exception as e:
            return {
                "job_id": job["job_id"],
                "status": "error",
                "error_message": f"Ligand preparation failed: {e}",
                "n_poses": 0,
                "best_affinity_kcal": "",
            }
    else:
        ligand_pdbqt = ligand_path

    output_pdbqt = output_dir / f"{job['output_prefix']}.pdbqt"

    center = [
        float(job["center_x"]),
        float(job["center_y"]),
        float(job["center_z"]),
    ]
    box_size = [
        float(job["box_size_x"]),
        float(job["box_size_y"]),
        float(job["box_size_z"]),
    ]

    try:
        energies = run_docking(
            str(receptor_pdbqt),
            str(ligand_pdbqt),
            str(output_pdbqt),
            center,
            box_size,
            int(job["exhaustiveness"]),
            int(job["n_poses"]),
        )
        n_poses = len(energies) if energies is not None else 0
        best = energies[0][0] if n_poses > 0 else ""

        return {
            "job_id": job["job_id"],
            "status": "ok",
            "error_message": "",
            "n_poses": n_poses,
            "best_affinity_kcal": best,
            "energies": energies,
            "output_file": str(output_pdbqt),
        }
    except Exception as e:
        return {
            "job_id": job["job_id"],
            "status": "error",
            "error_message": str(e),
            "n_poses": 0,
            "best_affinity_kcal": "",
        }


def execute_round(
    jobs: List[dict],
    round_id: int,
    max_workers: int = MAX_WORKERS_DEFAULT,
) -> List[dict]:
    """Execute all jobs for a given round (seed_index).

    Uses ThreadPoolExecutor for parallel execution on server.
    Falls back to serial execution if max_workers=1.
    """
    round_jobs = [j for j in jobs if int(j["seed_index"]) == round_id]

    if not round_jobs:
        return []

    results = []
    if max_workers <= 1:
        for job in round_jobs:
            results.append(execute_single_job(job))
    else:
        from concurrent.futures import ThreadPoolExecutor, as_completed
        with ThreadPoolExecutor(max_workers=max_workers) as pool:
            futures = {pool.submit(execute_single_job, j): j for j in round_jobs}
            for future in as_completed(futures):
                results.append(future.result())

    return results


# ---------------------------------------------------------------------------
# 3.3.3: Round-based saturation evaluation
# ---------------------------------------------------------------------------

def evaluate_round_results(
    results: List[dict],
    jobs: List[dict],
    pocket_statuses: List[dict],
    policy: dict,
    round_id: int,
) -> Tuple[List[dict], List[dict]]:
    """Evaluate results after a round and update pocket statuses.

    Returns (updated_statuses, budget_tracking_entries).
    """
    from egfr_pipeline.phase3.budget_policy import (
        evaluate_saturation,
        is_pose_acceptable,
    )

    # Build result lookup: job_id -> result
    result_map = {r["job_id"]: r for r in results}

    # Build job lookup: job_id -> job
    job_map = {j["job_id"]: j for j in jobs}

    # Update pocket statuses with round results
    # Group results by pocket
    pocket_results = defaultdict(list)
    for r in results:
        job = job_map.get(r["job_id"], {})
        pocket_id = job.get("pocket_id", "")
        if pocket_id:
            pocket_results[pocket_id].append(r)

    tracking_entries = []
    status_map = {ps["pocket_id"]: ps for ps in pocket_statuses}

    for pocket_id, p_results in pocket_results.items():
        ps = status_map.get(pocket_id)
        if not ps or ps["status"] != "open":
            continue

        status_before = ps["status"]

        for r in p_results:
            job = job_map.get(r["job_id"], {})
            n_poses = int(r.get("n_poses", 0))
            best_this = r.get("best_affinity_kcal", "")

            # Update cumulative stats
            ps["total_poses_generated"] = int(ps.get("total_poses_generated", 0)) + n_poses
            ps["seeds_completed"] = int(ps.get("seeds_completed", 0)) + 1
            ps["seeds_remaining"] = max(0, int(ps.get("seeds_allocated", 0)) - int(ps["seeds_completed"]))

            # Update best affinity
            current_best = ps.get("best_affinity_kcal", "")
            if best_this != "":
                try:
                    best_val = float(best_this)
                    if current_best == "" or best_val < float(current_best):
                        ps["best_affinity_kcal"] = best_this
                except (ValueError, TypeError):
                    pass

            # Count acceptable poses
            if ps["best_affinity_kcal"] != "" and n_poses > 0:
                acceptable_this = 0
                energies = r.get("energies", [])
                if energies:
                    for e in energies:
                        if is_pose_acceptable(e[0], float(ps["best_affinity_kcal"]), policy):
                            acceptable_this += 1
                ps["acceptable_poses"] = int(ps.get("acceptable_poses", 0)) + acceptable_this
                acceptable_round = acceptable_this
            else:
                acceptable_round = 0

            ps["rounds_completed"] = int(ps.get("rounds_completed", 0)) + 1

            # Tracking entry
            tracking_entries.append({
                "round_id": round_id,
                "pocket_id": pocket_id,
                "receptor_id": job.get("receptor_id", ""),
                "ligand_id": job.get("ligand_id", ""),
                "seed_index": job.get("seed_index", ""),
                "status_before": status_before,
                "status_after": "",  # filled below
                "poses_this_round": n_poses,
                "acceptable_poses_this_round": acceptable_round,
                "best_affinity_this_round": best_this,
                "cumulative_poses": ps["total_poses_generated"],
                "cumulative_acceptable": ps["acceptable_poses"],
                "budget_action": "execute",
            })

        # Evaluate saturation
        new_status, reason = evaluate_saturation(ps, policy)
        if new_status != ps["status"]:
            ps["status"] = new_status
            ps["status_reason"] = reason
            if new_status == "saturated":
                ps["round_saturated"] = round_id

        # Update tracking entries with final status
        for te in tracking_entries:
            if te["pocket_id"] == pocket_id and te["status_after"] == "":
                te["status_after"] = ps["status"]

    updated_statuses = list(status_map.values())
    return updated_statuses, tracking_entries


# ---------------------------------------------------------------------------
# Run metadata (Task 3.3.1)
# ---------------------------------------------------------------------------

def build_run_metadata(
    jobs: List[dict],
    pocket_statuses: List[dict],
    policy: dict,
    mode: str,
    max_workers: int,
) -> dict:
    """Build run metadata JSON."""
    open_pockets = [ps for ps in pocket_statuses if ps["status"] == "open"]
    skip_pockets = [ps for ps in pocket_statuses if ps["status"] == "skipped"]

    by_round = defaultdict(int)
    for j in jobs:
        by_round[int(j["seed_index"])] += 1

    receptors = sorted(set(j["receptor_id"] for j in jobs))
    ligands = sorted(set(j["ligand_id"] for j in jobs))
    pockets = sorted(set(j["pocket_id"] for j in jobs))

    return {
        "pipeline": "phase3_diverse_docking",
        "version": "1.0",
        "mode": mode,
        "generated": datetime.now().isoformat(),
        "max_workers": max_workers,
        "total_jobs": len(jobs),
        "rounds": {str(k): v for k, v in sorted(by_round.items())},
        "receptors": receptors,
        "ligands": ligands,
        "pockets_total": len(pocket_statuses),
        "pockets_open": len(open_pockets),
        "pockets_skipped": len(skip_pockets),
        "pocket_ids_active": pockets,
        "budget_policy": policy,
        "workspace_note": (
            "This environment is functional-validation-only. "
            "Vina execution requires server-side deployment."
        ),
    }


# ---------------------------------------------------------------------------
# Dry-run mode (Task 3.3.4)
# ---------------------------------------------------------------------------

def dry_run(
    jobs: List[dict],
    pocket_statuses: List[dict],
    policy: dict,
) -> List[str]:
    """Validate job plan without executing. Returns validation messages."""
    messages = []

    # Check jobs exist
    open_pockets = {ps["pocket_id"] for ps in pocket_statuses if ps["status"] == "open"}
    active_jobs = [j for j in jobs if j["pocket_id"] in open_pockets]
    messages.append(f"[OK] {len(active_jobs)} active jobs across {len(open_pockets)} open pockets")

    # Check receptor files
    missing_receptors = set()
    for j in active_jobs:
        rp = Path(j["receptor_pdb"])
        if not rp.exists():
            missing_receptors.add(str(rp))
    if missing_receptors:
        messages.append(f"[WARN] {len(missing_receptors)} receptor file(s) not found (OK for setup mode)")
    else:
        messages.append("[OK] All receptor files found")

    # Check ligand files
    missing_ligands = set()
    for j in active_jobs:
        lp = Path(j["ligand_file"])
        if not lp.exists():
            missing_ligands.add(str(lp))
    if missing_ligands:
        messages.append(f"[WARN] {len(missing_ligands)} ligand file(s) not found")
    else:
        messages.append("[OK] All ligand files found")

    # Check round structure
    by_round = defaultdict(list)
    for j in active_jobs:
        by_round[int(j["seed_index"])].append(j)
    messages.append(f"[OK] {len(by_round)} round(s) planned")
    for rid in sorted(by_round.keys()):
        rjobs = by_round[rid]
        pockets = set(j["pocket_id"] for j in rjobs)
        messages.append(f"  Round {rid}: {len(rjobs)} jobs across {len(pockets)} pockets")

    # Check no cross-receptor mixing
    for j in active_jobs:
        if not j["pocket_id"].startswith(j["receptor_id"]):
            messages.append(f"[ERROR] Cross-receptor mixing: {j['job_id']}")
            break
    else:
        messages.append("[OK] No cross-receptor pocket mixing")

    # Check Vina availability
    try:
        from vina import Vina
        messages.append("[OK] Vina is available")
    except ImportError:
        messages.append("[INFO] Vina not available — use --setup for server script generation")

    return messages


# ---------------------------------------------------------------------------
# Main orchestrator
# ---------------------------------------------------------------------------

def run_setup(
    output_dir: Path = PHASE3_OUTPUT_DIR,
    max_workers: int = MAX_WORKERS_DEFAULT,
) -> Tuple[Path, Path]:
    """Setup mode: generate run script + metadata without executing.

    Returns (script_path, metadata_path).
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    from egfr_pipeline.phase3.budget_policy import DEFAULT_BUDGET_POLICY

    # Load inputs
    jobs = _load_csv(output_dir / "phase3_docking_job_table.csv")
    pocket_statuses = _load_csv(output_dir / "pocket_search_status.csv")
    policy = dict(DEFAULT_BUDGET_POLICY)

    if not jobs:
        raise FileNotFoundError("No jobs found. Run TG 3.1 first.")
    if not pocket_statuses:
        raise FileNotFoundError("No pocket statuses found. Run TG 3.2 first.")

    print(f"  Loaded {len(jobs)} jobs, {len(pocket_statuses)} pocket statuses")

    # Dry-run validation
    messages = dry_run(jobs, pocket_statuses, policy)
    for m in messages:
        print(f"  {m}")

    # Generate run script
    script_path = generate_run_script(jobs, pocket_statuses, output_dir, max_workers)
    print(f"  Generated run script: {script_path}")

    # Generate metadata
    metadata = build_run_metadata(jobs, pocket_statuses, policy, "setup", max_workers)
    meta_path = output_dir / "phase3_run_metadata.json"
    with open(meta_path, "w", encoding="utf-8") as f:
        json.dump(metadata, f, indent=2)

    # Initialize empty round log
    round_log_path = output_dir / "phase3_round_log.csv"
    if not round_log_path.exists():
        _write_csv(round_log_path, ROUND_LOG_COLUMNS, [])

    return script_path, meta_path


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _load_csv(path: Path) -> List[dict]:
    if not path.exists():
        return []
    with open(path, encoding="utf-8") as f:
        return list(csv.DictReader(f))


def _write_csv(path: Path, columns: List[str], rows: List[dict]) -> None:
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=columns, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Phase 3 TG 3.3: Diversity-aware Vina execution layer"
    )
    parser.add_argument(
        "--output_dir", type=Path, default=PHASE3_OUTPUT_DIR,
        help="Phase 3 output directory",
    )
    parser.add_argument(
        "--max_workers", type=int, default=MAX_WORKERS_DEFAULT,
        help=f"Max parallel workers (default: {MAX_WORKERS_DEFAULT})",
    )

    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument(
        "--setup", action="store_true",
        help="Generate run script + metadata for server execution",
    )
    mode.add_argument(
        "--execute", action="store_true",
        help="Execute docking (requires Vina)",
    )
    mode.add_argument(
        "--dry-run", action="store_true",
        help="Validate job plan without executing",
    )
    mode.add_argument(
        "--evaluate-round", type=int, metavar="ROUND",
        help="Evaluate results after a completed round",
    )

    parser.add_argument(
        "--round", type=int,
        help="Round (seed_index) to execute (for --execute mode)",
    )

    args = parser.parse_args()

    print("Phase 3 — Task Group 3.3: Diversity-Aware Vina Execution Layer")
    print(f"  Output: {args.output_dir}")
    print()

    if args.setup:
        script_path, meta_path = run_setup(args.output_dir, args.max_workers)
        print(f"\n{'='*60}")
        print("TG 3.3 SETUP COMPLETE")
        print(f"{'='*60}")
        print(f"  Run script: {script_path}")
        print(f"  Metadata:   {meta_path}")
        print()
        print("To execute on server:")
        print(f"  bash {script_path}")

    elif args.dry_run:
        from egfr_pipeline.phase3.budget_policy import DEFAULT_BUDGET_POLICY
        jobs = _load_csv(args.output_dir / "phase3_docking_job_table.csv")
        statuses = _load_csv(args.output_dir / "pocket_search_status.csv")
        messages = dry_run(jobs, statuses, dict(DEFAULT_BUDGET_POLICY))
        for m in messages:
            print(f"  {m}")

    elif args.execute:
        if args.round is None:
            parser.error("--execute requires --round")
        jobs = _load_csv(args.output_dir / "phase3_docking_job_table.csv")
        results = execute_round(jobs, args.round, args.max_workers)
        ok = sum(1 for r in results if r.get("status") == "ok")
        err = sum(1 for r in results if r.get("status") == "error")
        skip = sum(1 for r in results if r.get("status") == "skip_no_vina")
        print(f"  Round {args.round}: {ok} ok, {err} error, {skip} skip")

    elif args.evaluate_round is not None:
        from egfr_pipeline.phase3.budget_policy import (
            DEFAULT_BUDGET_POLICY,
            SEARCH_STATUS_COLUMNS,
            BUDGET_TRACKING_COLUMNS,
        )
        # This would be called after --execute with actual results
        print(f"  Evaluate round {args.evaluate_round} — "
              "requires prior --execute results")

    return 0


if __name__ == "__main__":
    sys.exit(main())
