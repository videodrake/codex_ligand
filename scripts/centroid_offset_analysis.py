#!/usr/bin/env python3
"""Centroid offset correction simulation (AC-4.2, D-2.2).

Simulates how pocket_depth-based centroid correction affects PPI spatial
proximity scoring. For each pocket with known depth, computes:

    corrected_distance = raw_distance - alpha × pocket_depth

at alpha values [0.5, 0.6, 0.7, 0.8, 0.9, 1.0] and re-evaluates the
ppi_spatial_pts assignment at each alpha.

Inputs:
  - valid_sites.csv (for spatial_dist_to_ppi + pocket_depth_A)
  - OR synthetic data (--synthetic)

Outputs:
  - centroid_offset_correction.csv
  - Summary printed to stdout

Usage:
    python scripts/centroid_offset_analysis.py --synthetic
    python scripts/centroid_offset_analysis.py --input output/.../valid_sites.csv
"""

import argparse
import csv
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

PROJECT_ROOT = Path(__file__).resolve().parent.parent

# PPI spatial scoring thresholds (from verdict.py)
PPI_DIST_ADJACENT = 8.0
PPI_DIST_NEAR = 15.0
PPI_DIST_MODERATE = 25.0

ALPHA_VALUES = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0]


def classify_spatial(dist: float) -> Tuple[str, float]:
    """Classify PPI distance into spatial category and points.

    Returns (category, spatial_pts).
    """
    if dist < PPI_DIST_ADJACENT:
        return "adjacent", 10.0
    elif dist < PPI_DIST_NEAR:
        return "near", 6.0
    elif dist < PPI_DIST_MODERATE:
        return "moderate", 2.0
    else:
        return "distant", 0.0


def simulate_offset_correction(
    pockets: List[dict],
    alphas: List[float] = None,
) -> Tuple[List[dict], List[str]]:
    """Simulate centroid offset correction at multiple alpha values.

    Each pocket dict must have: pocket_id, raw_dist, pocket_depth.

    Returns (result_rows, warnings).
    """
    if alphas is None:
        alphas = list(ALPHA_VALUES)

    results = []
    change_count = 0

    for pocket in pockets:
        raw_dist = pocket["raw_dist"]
        depth = pocket["pocket_depth"]
        raw_cat, raw_pts = classify_spatial(raw_dist)

        row = {
            "pocket_id": pocket["pocket_id"],
            "receptor_id": pocket.get("receptor_id", ""),
            "raw_dist_A": round(raw_dist, 2),
            "pocket_depth_A": round(depth, 2),
            "raw_category": raw_cat,
            "raw_spatial_pts": raw_pts,
        }

        any_change = False
        for alpha in alphas:
            corrected = max(raw_dist - alpha * depth, 0.0)
            cat, pts = classify_spatial(corrected)
            row[f"alpha_{alpha:.1f}_dist"] = round(corrected, 2)
            row[f"alpha_{alpha:.1f}_cat"] = cat
            row[f"alpha_{alpha:.1f}_pts"] = pts
            if cat != raw_cat:
                any_change = True

        row["verdict_change"] = any_change
        if any_change:
            change_count += 1
        results.append(row)

    warnings = []
    if change_count >= 5:
        warnings.append(
            f"EC-4.2: 보정이 {change_count}개 포켓의 판정을 변경합니다. "
            f"적용 전 개별 포켓의 생물학적 의미를 확인해야 합니다."
        )

    return results, warnings


def generate_synthetic_offset_data() -> List[dict]:
    """Generate synthetic pockets with depth and distance for testing."""
    return [
        # Deep pocket far from PPI → correction moves to near
        {"pocket_id": "P_deep_far", "receptor_id": "R1",
         "raw_dist": 18.0, "pocket_depth": 8.0},
        # Shallow pocket near PPI → correction minimal
        {"pocket_id": "P_shallow_near", "receptor_id": "R1",
         "raw_dist": 10.0, "pocket_depth": 2.0},
        # Deep pocket moderate → correction may push to near
        {"pocket_id": "P_deep_moderate", "receptor_id": "R1",
         "raw_dist": 16.0, "pocket_depth": 6.0},
        # Surface pocket adjacent → stays adjacent
        {"pocket_id": "P_surface_adj", "receptor_id": "R1",
         "raw_dist": 6.0, "pocket_depth": 1.0},
        # Deep pocket at boundary 15Å → correction pushes to adjacent
        {"pocket_id": "P_boundary_deep", "receptor_id": "R1",
         "raw_dist": 14.5, "pocket_depth": 7.0},
        # Distant pocket with moderate depth
        {"pocket_id": "P_distant", "receptor_id": "R1",
         "raw_dist": 30.0, "pocket_depth": 5.0},
    ]


def main():
    parser = argparse.ArgumentParser(
        description="Centroid offset correction simulation (AC-4.2)"
    )
    parser.add_argument("--input", type=Path, default=None)
    parser.add_argument("--synthetic", action="store_true")
    parser.add_argument("--output_dir", type=Path,
                        default=PROJECT_ROOT / "output")
    args = parser.parse_args()

    print("Centroid Offset Correction Simulation (AC-4.2)")
    print()

    if args.synthetic or args.input is None:
        pockets = generate_synthetic_offset_data()
        print(f"  Using {len(pockets)} synthetic pockets")
    else:
        if not args.input.exists():
            print(f"  ERROR: {args.input} not found. Use --synthetic.")
            return 1
        with open(args.input, encoding="utf-8") as f:
            raw_rows = list(csv.DictReader(f))
        pockets = []
        for r in raw_rows:
            try:
                raw_dist = float(r.get("spatial_dist_to_ppi", ""))
                depth = float(r.get("pocket_depth_A", ""))
            except (ValueError, TypeError):
                continue
            pockets.append({
                "pocket_id": r.get("pocket_id", ""),
                "receptor_id": r.get("receptor_id", ""),
                "raw_dist": raw_dist,
                "pocket_depth": depth,
            })
        print(f"  Loaded {len(pockets)} pockets with depth data")

    results, warnings = simulate_offset_correction(pockets)

    # Print summary
    for r in results:
        change = "CHANGED" if r["verdict_change"] else "same"
        print(f"  {r['pocket_id']:20s} raw={r['raw_dist_A']:5.1f}Å depth={r['pocket_depth_A']:4.1f}Å "
              f"→ α=0.7: {r['alpha_0.7_dist']:5.1f}Å ({r['alpha_0.7_cat']}) [{change}]")

    for w in warnings:
        print(f"  WARNING: {w}")

    # Write CSV
    args.output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = args.output_dir / "centroid_offset_correction.csv"
    if results:
        cols = list(results[0].keys())
        with open(csv_path, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=cols)
            w.writeheader()
            w.writerows(results)
        print(f"\n  CSV: {csv_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
