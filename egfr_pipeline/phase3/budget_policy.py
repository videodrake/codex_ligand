#!/usr/bin/env python3
"""Phase 3 Task Group 3.2: Search Budget Policy and Saturation Rules.

Defines and implements an explicit search-budget policy so that already
well-sampled pockets stop consuming disproportionate docking effort.

Tasks:
  3.2.1 - Budget parameter definition
  3.2.2 - Saturation rule definition
  3.2.3 - Budget reallocation policy
  3.2.4 - Saturation transparency

Inputs:
  - output/phase3_docking/phase3_docking_job_table.csv  (from TG 3.1)
  - output/phase3_docking/phase3_candidate_reference_normalized.csv  (from TG 3.0)

Outputs:
  - output/phase3_docking/phase3_budget_policy.md
  - output/phase3_docking/pocket_search_status.csv
  - output/phase3_docking/phase3_budget_tracking.csv

Usage:
    python -m egfr_pipeline.phase3.budget_policy [--output_dir ...]
"""

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

PHASE3_OUTPUT_DIR = PROJECT_ROOT / "output" / "phase3_docking"

# Pocket status labels (Task 3.2.2)
STATUS_OPEN = "open"
STATUS_SATURATED = "saturated"
STATUS_SKIPPED = "skipped"
STATUS_EXHAUSTED = "exhausted"

VALID_STATUSES = {STATUS_OPEN, STATUS_SATURATED, STATUS_SKIPPED, STATUS_EXHAUSTED}

# ---------------------------------------------------------------------------
# 3.2.1: Budget parameter definition
# ---------------------------------------------------------------------------

# Default budget policy — configurable per deployment
DEFAULT_BUDGET_POLICY = {
    # Max rounds of docking per pocket (each round = 1 seed)
    "max_rounds": 3,
    # Per-pocket saturation: stop when this many acceptable poses accumulated
    "saturation_pose_threshold": 10,
    # Affinity window: poses within this many kcal/mol of the best are "acceptable"
    "saturation_affinity_window_kcal": 2.0,
    # Minimum poses before saturation can trigger (avoid premature saturation)
    "min_poses_before_saturation": 5,
    # Per-pocket max total poses across all rounds
    "max_poses_per_pocket": 60,
    # Per-pocket max seeds
    "max_seeds_per_pocket": 3,
    # Budget reallocation: when a pocket saturates, redistribute its remaining
    # seeds to unsaturated pockets of the same or lower priority tier
    "enable_reallocation": True,
    # Priority order for reallocation (higher priority gets budget first)
    "priority_order": ["primary", "secondary", "exploratory"],
}

BUDGET_POLICY_COLUMNS = [
    "parameter",
    "value",
    "description",
]

# Per-pocket search status schema (Task 3.2.2)
SEARCH_STATUS_COLUMNS = [
    "pocket_id",
    "receptor_id",
    "docking_priority",
    "status",
    "status_reason",
    # Budget consumed
    "seeds_allocated",
    "seeds_completed",
    "seeds_remaining",
    "total_poses_generated",
    "acceptable_poses",
    "best_affinity_kcal",
    # Budget limits
    "max_seeds",
    "max_poses",
    "saturation_threshold",
    # Round tracking
    "rounds_completed",
    "round_saturated",
]

# Budget tracking per round (Task 3.2.4)
BUDGET_TRACKING_COLUMNS = [
    "round_id",
    "pocket_id",
    "receptor_id",
    "ligand_id",
    "seed_index",
    "status_before",
    "status_after",
    "poses_this_round",
    "acceptable_poses_this_round",
    "best_affinity_this_round",
    "cumulative_poses",
    "cumulative_acceptable",
    "budget_action",
]


# ---------------------------------------------------------------------------
# 3.2.2: Saturation rule definition
# ---------------------------------------------------------------------------

def evaluate_saturation(
    pocket_state: dict,
    policy: dict,
) -> Tuple[str, str]:
    """Evaluate whether a pocket should transition to saturated/exhausted.

    Args:
        pocket_state: Current pocket tracking state with keys:
            - acceptable_poses: int
            - total_poses_generated: int
            - seeds_completed: int
            - seeds_allocated: int
        policy: Budget policy dict.

    Returns:
        (new_status, reason)
    """
    acceptable = pocket_state.get("acceptable_poses", 0)
    total_poses = pocket_state.get("total_poses_generated", 0)
    seeds_completed = pocket_state.get("seeds_completed", 0)
    seeds_allocated = pocket_state.get("seeds_allocated", 0)

    # Check saturation: enough acceptable poses
    sat_threshold = policy.get("saturation_pose_threshold", 10)
    min_before = policy.get("min_poses_before_saturation", 5)
    if (total_poses >= min_before and acceptable >= sat_threshold):
        return STATUS_SATURATED, (
            f"accumulated {acceptable} acceptable poses "
            f"(threshold={sat_threshold}, min_poses={min_before})"
        )

    # Check exhausted: all allocated seeds used up
    max_poses = policy.get("max_poses_per_pocket", 60)
    if total_poses >= max_poses:
        return STATUS_EXHAUSTED, (
            f"reached max_poses_per_pocket={max_poses} "
            f"({total_poses} total)"
        )

    if seeds_completed >= seeds_allocated:
        return STATUS_EXHAUSTED, (
            f"all {seeds_allocated} allocated seeds completed"
        )

    return STATUS_OPEN, ""


def is_pose_acceptable(
    affinity: float,
    best_affinity: float,
    policy: dict,
) -> bool:
    """Check if a pose falls within the saturation affinity window.

    Args:
        affinity: This pose's affinity (kcal/mol, negative = better).
        best_affinity: Best affinity seen so far for this pocket.
        policy: Budget policy dict.
    """
    window = policy.get("saturation_affinity_window_kcal", 2.0)
    # Both values are negative; acceptable if within window of best
    return affinity <= (best_affinity + window)


# ---------------------------------------------------------------------------
# 3.2.3: Budget reallocation policy
# ---------------------------------------------------------------------------

def compute_reallocation(
    pocket_statuses: List[dict],
    policy: dict,
) -> List[dict]:
    """Compute budget reallocation from saturated to unsaturated pockets.

    When a pocket saturates, its remaining seed budget is redistributed
    to open pockets by priority order (primary > secondary > exploratory).

    Returns list of reallocation actions:
        [{"from_pocket", "to_pocket", "seeds_reallocated", "reason"}]
    """
    if not policy.get("enable_reallocation", True):
        return []

    priority_order = policy.get("priority_order",
                                ["primary", "secondary", "exploratory"])

    # Collect surplus from saturated pockets
    surplus_seeds = 0
    donors = []
    for ps in pocket_statuses:
        if ps["status"] in (STATUS_SATURATED, STATUS_EXHAUSTED):
            remaining = ps["seeds_remaining"]
            if remaining > 0:
                surplus_seeds += remaining
                donors.append(ps["pocket_id"])

    if surplus_seeds == 0:
        return []

    # Distribute to open pockets by priority
    actions = []
    for priority in priority_order:
        open_pockets = [
            ps for ps in pocket_statuses
            if ps["status"] == STATUS_OPEN
            and ps["docking_priority"] == priority
        ]
        for ps in open_pockets:
            if surplus_seeds <= 0:
                break
            max_extra = ps["max_seeds"] - ps["seeds_allocated"]
            give = min(surplus_seeds, max(max_extra, 0))
            if give > 0:
                actions.append({
                    "from_pocket": ",".join(donors),
                    "to_pocket": ps["pocket_id"],
                    "seeds_reallocated": give,
                    "reason": f"reallocation from saturated pocket(s)",
                })
                surplus_seeds -= give

    return actions


# ---------------------------------------------------------------------------
# Initialize pocket search status from job table + normalized reference
# ---------------------------------------------------------------------------

def initialize_pocket_status(
    normalized_pockets: List[dict],
    jobs: List[dict],
    policy: dict,
) -> List[dict]:
    """Create initial pocket_search_status from normalized pockets and job table.

    Before any docking has run, all dockable pockets start as 'open'
    and non-dockable as 'skipped'.
    """
    # Count jobs per pocket
    jobs_per_pocket = defaultdict(int)
    seeds_per_pocket = defaultdict(set)
    for j in jobs:
        pid = j["pocket_id"]
        jobs_per_pocket[pid] += 1
        seeds_per_pocket[pid].add(int(j.get("seed_index", 0)))

    statuses = []
    for p in normalized_pockets:
        pocket_id = p["pocket_id"]
        dockable = p.get("dockable", "False") == "True"
        priority = p.get("docking_priority", "skip")

        max_seeds = int(p.get("recommended_seeds", 0))
        max_poses = policy.get("max_poses_per_pocket", 60)
        sat_threshold = policy.get("saturation_pose_threshold", 10)

        if not dockable:
            status = STATUS_SKIPPED
            reason = p.get("skip_reason", "not dockable")
        else:
            status = STATUS_OPEN
            reason = ""

        seeds_allocated = len(seeds_per_pocket.get(pocket_id, set()))

        statuses.append({
            "pocket_id": pocket_id,
            "receptor_id": p.get("receptor_id", ""),
            "docking_priority": priority,
            "status": status,
            "status_reason": reason,
            "seeds_allocated": seeds_allocated,
            "seeds_completed": 0,
            "seeds_remaining": seeds_allocated,
            "total_poses_generated": 0,
            "acceptable_poses": 0,
            "best_affinity_kcal": "",
            "max_seeds": max_seeds,
            "max_poses": max_poses,
            "saturation_threshold": sat_threshold,
            "rounds_completed": 0,
            "round_saturated": "",
        })

    return statuses


# ---------------------------------------------------------------------------
# Budget tracking initialization (empty — filled during execution)
# ---------------------------------------------------------------------------

def initialize_budget_tracking() -> List[dict]:
    """Return empty budget tracking table (populated during TG 3.3 execution)."""
    return []


# ---------------------------------------------------------------------------
# Policy document generation (Task 3.2.4)
# ---------------------------------------------------------------------------

def build_policy_document(policy: dict, pocket_statuses: List[dict]) -> List[str]:
    """Generate the budget policy markdown document."""
    lines = [
        "# Phase 3 Search Budget Policy",
        "",
        "## Policy Parameters",
        "",
        "| Parameter | Value | Description |",
        "|-----------|-------|-------------|",
    ]

    param_descriptions = {
        "max_rounds": "Maximum docking rounds per pocket",
        "saturation_pose_threshold": "Acceptable poses needed to saturate a pocket",
        "saturation_affinity_window_kcal": "Affinity window (kcal/mol) from best pose",
        "min_poses_before_saturation": "Minimum poses before saturation can trigger",
        "max_poses_per_pocket": "Hard cap on total poses per pocket",
        "max_seeds_per_pocket": "Maximum seeds per pocket",
        "enable_reallocation": "Redistribute budget from saturated pockets",
    }

    for key, desc in param_descriptions.items():
        val = policy.get(key, "N/A")
        lines.append(f"| {key} | {val} | {desc} |")
    lines.append("")

    # Saturation rules
    lines.extend([
        "## Saturation Rules",
        "",
        "A pocket transitions from `open` to `saturated` when:",
        "",
        f"1. At least **{policy.get('min_poses_before_saturation', 5)}** "
        f"total poses have been generated, AND",
        f"2. At least **{policy.get('saturation_pose_threshold', 10)}** "
        f"poses fall within **{policy.get('saturation_affinity_window_kcal', 2.0)} "
        f"kcal/mol** of the best affinity.",
        "",
        "A pocket transitions to `exhausted` when all allocated seeds are "
        "completed or max_poses_per_pocket is reached.",
        "",
    ])

    # Status labels
    lines.extend([
        "## Pocket Status Labels",
        "",
        "| Status | Meaning |",
        "|--------|---------|",
        "| open | Eligible for more docking rounds |",
        "| saturated | Enough acceptable poses; stop allocating budget |",
        "| skipped | Not dockable (tier_3 + low_relevance) |",
        "| exhausted | All allocated seeds used; no more budget available |",
        "",
    ])

    # Reallocation
    lines.extend([
        "## Budget Reallocation",
        "",
    ])
    if policy.get("enable_reallocation", True):
        priority_order = policy.get("priority_order",
                                    ["primary", "secondary", "exploratory"])
        lines.extend([
            "When a pocket saturates, its remaining seed budget is "
            "redistributed to open pockets in priority order:",
            "",
            f"  {' > '.join(priority_order)}",
            "",
            "This ensures under-sampled pockets receive additional search "
            "effort from over-sampled ones.",
            "",
        ])
    else:
        lines.extend([
            "Budget reallocation is **disabled**. Each pocket uses only its "
            "originally allocated budget.",
            "",
        ])

    # Current pocket status summary
    lines.extend([
        "## Initial Pocket Status",
        "",
        "| Pocket | Receptor | Priority | Status | Seeds | Max Poses |",
        "|--------|----------|----------|--------|-------|-----------|",
    ])
    for ps in pocket_statuses:
        lines.append(
            f"| {ps['pocket_id']} | {ps['receptor_id']} | "
            f"{ps['docking_priority']} | {ps['status']} | "
            f"{ps['seeds_allocated']} | {ps['max_poses']} |"
        )
    lines.append("")

    # Budget summary
    open_pockets = [ps for ps in pocket_statuses if ps["status"] == STATUS_OPEN]
    skip_pockets = [ps for ps in pocket_statuses if ps["status"] == STATUS_SKIPPED]
    total_seeds = sum(ps["seeds_allocated"] for ps in open_pockets)
    lines.extend([
        "## Budget Summary",
        "",
        f"- Open pockets: {len(open_pockets)}",
        f"- Skipped pockets: {len(skip_pockets)}",
        f"- Total seeds allocated: {total_seeds}",
        f"- Estimated total jobs (all ligands): "
        f"see phase3_docking_job_table.csv",
        "",
        "---",
        "",
        "Generated by `egfr_pipeline.phase3.budget_policy`",
    ])

    return lines


# ---------------------------------------------------------------------------
# Main orchestrator
# ---------------------------------------------------------------------------

def run_budget_policy(
    output_dir: Path = PHASE3_OUTPUT_DIR,
    policy: Optional[dict] = None,
) -> Tuple[Path, Path, Path]:
    """Run full TG 3.2 pipeline.

    Returns (policy_path, status_path, tracking_path).
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    if policy is None:
        policy = dict(DEFAULT_BUDGET_POLICY)

    # Load TG 3.0 normalized pockets
    norm_path = output_dir / "phase3_candidate_reference_normalized.csv"
    if not norm_path.exists():
        raise FileNotFoundError(
            f"Normalized pocket reference not found: {norm_path}\n"
            "Run TG 3.0 first."
        )
    normalized = _load_csv(norm_path)
    print(f"  Loaded {len(normalized)} pockets from normalized reference")

    # Load TG 3.1 job table
    job_path = output_dir / "phase3_docking_job_table.csv"
    if not job_path.exists():
        raise FileNotFoundError(
            f"Job table not found: {job_path}\n"
            "Run TG 3.1 first."
        )
    jobs = _load_csv(job_path)
    print(f"  Loaded {len(jobs)} docking jobs")

    # Task 3.2.1-3.2.2: Initialize pocket search status
    pocket_statuses = initialize_pocket_status(normalized, jobs, policy)
    open_count = sum(1 for ps in pocket_statuses if ps["status"] == STATUS_OPEN)
    skip_count = sum(1 for ps in pocket_statuses if ps["status"] == STATUS_SKIPPED)
    print(f"  Pocket status: {open_count} open, {skip_count} skipped")

    # Task 3.2.3: Check reallocation (at init, no saturation yet)
    realloc = compute_reallocation(pocket_statuses, policy)
    if realloc:
        print(f"  Reallocation actions: {len(realloc)}")
    else:
        print("  No reallocation needed (no saturated pockets yet)")

    # Task 3.2.4: Initialize empty budget tracking
    tracking = initialize_budget_tracking()

    # Write outputs
    status_path = output_dir / "pocket_search_status.csv"
    _write_csv(status_path, SEARCH_STATUS_COLUMNS, pocket_statuses)

    tracking_path = output_dir / "phase3_budget_tracking.csv"
    _write_csv(tracking_path, BUDGET_TRACKING_COLUMNS, tracking)

    # Generate policy document
    policy_lines = build_policy_document(policy, pocket_statuses)
    policy_path = output_dir / "phase3_budget_policy.md"
    with open(policy_path, "w", encoding="utf-8") as f:
        f.write("\n".join(policy_lines))

    return policy_path, status_path, tracking_path


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
        description="Phase 3 TG 3.2: Search budget policy and saturation rules"
    )
    parser.add_argument(
        "--output_dir", type=Path, default=PHASE3_OUTPUT_DIR,
        help="Phase 3 output directory",
    )
    args = parser.parse_args()

    print("Phase 3 — Task Group 3.2: Search Budget Policy and Saturation Rules")
    print(f"  Phase 3 output: {args.output_dir}")
    print()

    policy_path, status_path, tracking_path = run_budget_policy(args.output_dir)

    print(f"\n{'='*60}")
    print("TG 3.2 COMPLETE")
    print(f"{'='*60}")
    print(f"  Policy:   {policy_path}")
    print(f"  Status:   {status_path}")
    print(f"  Tracking: {tracking_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
