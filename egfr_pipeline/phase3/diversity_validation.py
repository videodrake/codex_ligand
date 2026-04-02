#!/usr/bin/env python3
"""Phase 3 Task Group 3.5: Pocket Occupancy and Diversity Validation.

Validates that the diversity-aware workflow distributes search effort
across candidate pockets rather than over-concentrating on dominant ones.

Tasks:
  3.5.1 - Pocket occupancy summary (per receptor, ligand, pocket)
  3.5.2 - Diversity metrics (concentration ratio, budget fraction)
  3.5.3 - Naive blind baseline comparison (optional)

Inputs:
  - output/phase3_docking/vina_pose_table.csv
  - output/phase3_docking/pocket_search_status.csv
  - output/phase3_docking/phase3_candidate_reference_normalized.csv

Outputs:
  - output/phase3_docking/phase3_pocket_occupancy_summary.csv
  - output/phase3_docking/phase3_diversity_metrics.csv
  - output/phase3_docking/phase3_blind_vs_diverse_comparison.csv (optional)

Usage:
    python -m egfr_pipeline.phase3.diversity_validation [--output_dir ...]
"""

import argparse
import csv
import math
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from egfr_pipeline.csv_utils import load_csv

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

PHASE3_OUTPUT_DIR = PROJECT_ROOT / "output" / "phase3_docking"

OCCUPANCY_COLUMNS = [
    "receptor_id",
    "ligand_id",
    "candidate_pocket_id",
    "docking_priority",
    "n_poses",
    "n_seeds",
    "best_affinity_kcal",
    "mean_affinity_kcal",
    "pocket_status",
    "fraction_of_receptor_poses",
]

DIVERSITY_COLUMNS = [
    "receptor_id",
    "ligand_id",
    "n_pockets_total",
    "n_pockets_sampled",
    "n_pockets_unsampled",
    "n_pockets_saturated",
    "coverage_fraction",
    "top1_concentration",
    "top2_concentration",
    "shannon_entropy",
    "max_possible_entropy",
    "normalized_entropy",
    "budget_top1_fraction",
    "diversity_verdict",
]

COMPARISON_COLUMNS = [
    "receptor_id",
    "ligand_id",
    "metric",
    "blind_value",
    "diverse_value",
    "improvement",
]


# ---------------------------------------------------------------------------
# 3.5.1: Pocket occupancy summary
# ---------------------------------------------------------------------------

def compute_occupancy(
    pose_rows: List[dict],
    pocket_statuses: List[dict],
    normalized_pockets: List[dict],
) -> List[dict]:
    """Compute per-receptor, per-ligand, per-pocket occupancy summary."""
    # Build status lookup
    status_map = {ps["pocket_id"]: ps for ps in pocket_statuses}
    priority_map = {p["pocket_id"]: p.get("docking_priority", "")
                    for p in normalized_pockets}

    # Group poses by (receptor, ligand, pocket)
    groups = defaultdict(list)
    for row in pose_rows:
        key = (
            row.get("receptor_id", ""),
            row.get("ligand_id", ""),
            row.get("candidate_pocket_id", ""),
        )
        groups[key].append(row)

    # Total poses per (receptor, ligand) for fraction calculation
    receptor_ligand_totals = defaultdict(int)
    for row in pose_rows:
        receptor_ligand_totals[(row.get("receptor_id", ""),
                                row.get("ligand_id", ""))] += 1

    # Build all possible (receptor, ligand, pocket) combinations
    receptors = sorted(set(r.get("receptor_id", "") for r in pose_rows))
    ligands = sorted(set(r.get("ligand_id", "") for r in pose_rows))
    pockets_by_receptor = defaultdict(set)
    for p in normalized_pockets:
        pockets_by_receptor[p["receptor_id"]].add(p["pocket_id"])

    results = []
    for receptor_id in receptors:
        for ligand_id in ligands:
            total = receptor_ligand_totals.get((receptor_id, ligand_id), 0)
            for pocket_id in sorted(pockets_by_receptor.get(receptor_id, set())):
                key = (receptor_id, ligand_id, pocket_id)
                rows = groups.get(key, [])

                affinities = []
                seeds = set()
                for r in rows:
                    aff = r.get("affinity", "")
                    if aff != "" and aff is not None:
                        try:
                            affinities.append(float(aff))
                        except (ValueError, TypeError):
                            pass
                    seed = r.get("seed", "")
                    if seed:
                        seeds.add(seed)

                best = min(affinities) if affinities else ""
                mean = sum(affinities) / len(affinities) if affinities else ""
                fraction = len(rows) / total if total > 0 else 0.0

                ps = status_map.get(pocket_id, {})

                results.append({
                    "receptor_id": receptor_id,
                    "ligand_id": ligand_id,
                    "candidate_pocket_id": pocket_id,
                    "docking_priority": priority_map.get(pocket_id, ""),
                    "n_poses": len(rows),
                    "n_seeds": len(seeds),
                    "best_affinity_kcal": f"{best:.2f}" if isinstance(best, float) else "",
                    "mean_affinity_kcal": f"{mean:.2f}" if isinstance(mean, float) else "",
                    "pocket_status": ps.get("status", ""),
                    "fraction_of_receptor_poses": f"{fraction:.3f}",
                })

    return results


# ---------------------------------------------------------------------------
# 3.5.2: Diversity metrics
# ---------------------------------------------------------------------------

def compute_diversity_metrics(
    occupancy: List[dict],
    normalized_pockets: List[dict],
) -> List[dict]:
    """Compute diversity metrics per (receptor, ligand) pair."""
    # Group occupancy by (receptor, ligand)
    groups = defaultdict(list)
    for row in occupancy:
        key = (row["receptor_id"], row["ligand_id"])
        groups[key].append(row)

    # Count dockable pockets per receptor
    dockable_by_receptor = defaultdict(int)
    for p in normalized_pockets:
        if p.get("dockable", "False") == "True":
            dockable_by_receptor[p["receptor_id"]] += 1

    results = []
    for (receptor_id, ligand_id), occ_rows in sorted(groups.items()):
        n_total = dockable_by_receptor.get(receptor_id, 0)
        pose_counts = [int(r["n_poses"]) for r in occ_rows]
        sampled = [r for r in occ_rows if int(r["n_poses"]) > 0
                   and r["docking_priority"] != "skip"]
        saturated = [r for r in occ_rows if r["pocket_status"] == "saturated"]

        n_sampled = len(sampled)
        n_unsampled = n_total - n_sampled
        total_poses = sum(pose_counts)

        # Coverage
        coverage = n_sampled / n_total if n_total > 0 else 0.0

        # Concentration ratios (fraction of poses in top-N pockets)
        sorted_counts = sorted(pose_counts, reverse=True)
        top1 = sorted_counts[0] / total_poses if total_poses > 0 else 0.0
        top2 = (sum(sorted_counts[:2]) / total_poses
                if total_poses > 0 and len(sorted_counts) >= 2 else top1)

        # Shannon entropy of pose distribution
        entropy = 0.0
        for c in pose_counts:
            if c > 0 and total_poses > 0:
                p = c / total_poses
                entropy -= p * math.log2(p)
        max_entropy = math.log2(n_total) if n_total > 1 else 0.0
        norm_entropy = entropy / max_entropy if max_entropy > 0 else 0.0

        # Budget fraction for top pocket
        budget_top1 = top1  # Same as concentration for now

        # Verdict
        if n_total <= 1:
            verdict = "single_pocket"
        elif coverage >= 0.8 and norm_entropy >= 0.7:
            verdict = "well_distributed"
        elif coverage >= 0.5 and norm_entropy >= 0.4:
            verdict = "moderately_distributed"
        elif top1 >= 0.8:
            verdict = "over_concentrated"
        else:
            verdict = "under_sampled"

        results.append({
            "receptor_id": receptor_id,
            "ligand_id": ligand_id,
            "n_pockets_total": n_total,
            "n_pockets_sampled": n_sampled,
            "n_pockets_unsampled": n_unsampled,
            "n_pockets_saturated": len(saturated),
            "coverage_fraction": f"{coverage:.3f}",
            "top1_concentration": f"{top1:.3f}",
            "top2_concentration": f"{top2:.3f}",
            "shannon_entropy": f"{entropy:.3f}",
            "max_possible_entropy": f"{max_entropy:.3f}",
            "normalized_entropy": f"{norm_entropy:.3f}",
            "budget_top1_fraction": f"{budget_top1:.3f}",
            "diversity_verdict": verdict,
        })

    return results


# ---------------------------------------------------------------------------
# 3.5.3: Naive blind baseline comparison (optional)
# ---------------------------------------------------------------------------

def compare_with_blind_baseline(
    diversity_metrics: List[dict],
    blind_pose_table_path: Optional[Path] = None,
) -> List[dict]:
    """Compare diversity-aware metrics against naive blind baseline.

    If no blind baseline exists, generates a comparison note explaining
    that no baseline is available.
    """
    comparisons = []

    if blind_pose_table_path and blind_pose_table_path.exists():
        # Load blind baseline and compute metrics
        # (Future: implement when blind baseline data is available)
        pass

    # Always generate a placeholder comparison showing the diverse metrics
    for dm in diversity_metrics:
        for metric_name in ["coverage_fraction", "top1_concentration",
                            "normalized_entropy"]:
            comparisons.append({
                "receptor_id": dm["receptor_id"],
                "ligand_id": dm["ligand_id"],
                "metric": metric_name,
                "blind_value": "",  # No baseline yet
                "diverse_value": dm[metric_name],
                "improvement": "",  # Cannot compute without baseline
            })

    return comparisons


# ---------------------------------------------------------------------------
# Main orchestrator
# ---------------------------------------------------------------------------

def run_diversity_validation(
    output_dir: Path = PHASE3_OUTPUT_DIR,
) -> Tuple[Path, Path, Path]:
    """Run full TG 3.5 pipeline.

    Returns (occupancy_path, diversity_path, comparison_path).
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load inputs
    pose_rows = load_csv(output_dir / "vina_pose_table.csv")
    pocket_statuses = load_csv(output_dir / "pocket_search_status.csv")
    normalized = load_csv(output_dir / "phase3_candidate_reference_normalized.csv")

    print(f"  Loaded {len(pose_rows)} pose rows, "
          f"{len(pocket_statuses)} pocket statuses, "
          f"{len(normalized)} pockets")

    # 3.5.1: Occupancy
    occupancy = compute_occupancy(pose_rows, pocket_statuses, normalized)
    print(f"  Computed occupancy: {len(occupancy)} entries")

    # 3.5.2: Diversity metrics
    diversity = compute_diversity_metrics(occupancy, normalized)
    for dm in diversity:
        print(f"    {dm['receptor_id']}/{dm['ligand_id']}: "
              f"coverage={dm['coverage_fraction']}, "
              f"entropy={dm['normalized_entropy']}, "
              f"verdict={dm['diversity_verdict']}")

    # 3.5.3: Blind baseline comparison
    comparison = compare_with_blind_baseline(diversity)

    # Write outputs
    occ_path = output_dir / "phase3_pocket_occupancy_summary.csv"
    _write_csv(occ_path, OCCUPANCY_COLUMNS, occupancy)

    div_path = output_dir / "phase3_diversity_metrics.csv"
    _write_csv(div_path, DIVERSITY_COLUMNS, diversity)

    cmp_path = output_dir / "phase3_blind_vs_diverse_comparison.csv"
    _write_csv(cmp_path, COMPARISON_COLUMNS, comparison)

    return occ_path, div_path, cmp_path


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

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
        description="Phase 3 TG 3.5: Pocket occupancy and diversity validation"
    )
    parser.add_argument(
        "--output_dir", type=Path, default=PHASE3_OUTPUT_DIR,
        help="Phase 3 output directory",
    )
    args = parser.parse_args()

    print("Phase 3 — Task Group 3.5: Pocket Occupancy and Diversity Validation")
    print(f"  Output: {args.output_dir}")
    print()

    occ_path, div_path, cmp_path = run_diversity_validation(args.output_dir)

    print(f"\n{'='*60}")
    print("TG 3.5 COMPLETE")
    print(f"{'='*60}")
    print(f"  Occupancy:  {occ_path}")
    print(f"  Diversity:  {div_path}")
    print(f"  Comparison: {cmp_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
