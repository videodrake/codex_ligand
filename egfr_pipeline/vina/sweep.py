#!/usr/bin/env python3
"""Pocket cutoff sensitivity sweep.

Re-clusters an existing vina_pose_table.csv across a range of pocket_cutoff
values and reports how pocket counts and compositions change.

Usage:
    python -m egfr_pipeline.vina.sweep config/example-project.yaml
    python -m egfr_pipeline.vina.sweep config/example-project.yaml --min 2.0 --max 8.0 --step 0.5
"""
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from egfr_pipeline.config import load_config, project_root_from_config
from egfr_pipeline.vina.cluster import (
    assign_pockets,
    load_pose_table,
    merge_pockets_by_residue,
)


def _count_pockets(rows: List[dict]) -> Dict[str, int]:
    """Return {receptor_id: n_pockets}."""
    pockets: Dict[str, set] = defaultdict(set)
    for row in rows:
        pid = row.get("pocket_id", "")
        if pid:
            pockets[row["receptor_id"]].add(pid)
    return {k: len(v) for k, v in sorted(pockets.items())}


def _pocket_sizes(rows: List[dict]) -> Dict[str, Dict[str, int]]:
    """Return {receptor_id: {pocket_id: n_poses}}."""
    result: Dict[str, Dict[str, int]] = defaultdict(lambda: Counter())
    for row in rows:
        pid = row.get("pocket_id", "")
        if pid:
            result[row["receptor_id"]][pid] += 1
    return dict(result)


def sweep(
    pose_table_path: str,
    cutoff_min: float = 2.0,
    cutoff_max: float = 12.0,
    cutoff_step: float = 0.5,
    merge_by_residue: bool = True,
    merge_jaccard: float = 0.3,
    merge_overlap: float = 0.5,
    merge_centroid_fallback: float = 6.0,
) -> List[dict]:
    """Run pocket clustering at each cutoff and collect statistics.

    Returns list of dicts with fields:
      cutoff, receptor_id, n_pockets, pocket_sizes, largest_pocket, smallest_pocket
    """
    import copy

    base_rows = load_pose_table(Path(pose_table_path))
    results: List[dict] = []
    cutoff = cutoff_min

    while cutoff <= cutoff_max + 1e-9:
        # Deep copy to avoid mutation
        rows = copy.deepcopy(base_rows)
        # Strip existing pocket_id so assign_pockets starts fresh
        for row in rows:
            row.pop("pocket_id", None)

        clustered = assign_pockets(rows, cutoff)
        if merge_by_residue:
            clustered = merge_pockets_by_residue(
                clustered,
                jaccard_threshold=merge_jaccard,
                overlap_threshold=merge_overlap,
                centroid_fallback_cutoff=merge_centroid_fallback,
            )

        sizes_by_receptor = _pocket_sizes(clustered)

        for receptor_id in sorted(sizes_by_receptor.keys()):
            sizes = sizes_by_receptor[receptor_id]
            size_values = sorted(sizes.values(), reverse=True)
            results.append({
                "cutoff": round(cutoff, 2),
                "receptor_id": receptor_id,
                "n_pockets": len(sizes),
                "pocket_sizes": ";".join(f"{pid}:{cnt}" for pid, cnt in
                                         sorted(sizes.items(), key=lambda x: (-x[1], x[0]))),
                "largest_pocket": size_values[0] if size_values else 0,
                "smallest_pocket": size_values[-1] if size_values else 0,
            })

        cutoff += cutoff_step

    return results


SWEEP_FIELDS = [
    "cutoff",
    "receptor_id",
    "n_pockets",
    "largest_pocket",
    "smallest_pocket",
    "pocket_sizes",
]


def sweep_from_config(
    config_path: str,
    cutoff_min: float = 2.0,
    cutoff_max: float = 8.0,
    cutoff_step: float = 0.5,
    output_path: Optional[str] = None,
) -> Path:
    """Run sweep and save results CSV."""
    config = load_config(config_path)
    project_root = project_root_from_config(config)
    pose_table = project_root / "vina_pose_table.csv"

    if not pose_table.exists():
        raise FileNotFoundError(f"Pose table not found: {pose_table}")

    pp = config.get("postprocess", {})
    merge = pp.get("merge_by_residue", True)
    merge_j = pp.get("merge_jaccard", 0.3)
    merge_oc = pp.get("merge_overlap", 0.5)
    merge_cf = pp.get("merge_centroid_fallback", 6.0)

    results = sweep(
        str(pose_table),
        cutoff_min=cutoff_min,
        cutoff_max=cutoff_max,
        cutoff_step=cutoff_step,
        merge_by_residue=merge,
        merge_jaccard=merge_j,
        merge_overlap=merge_oc,
        merge_centroid_fallback=merge_cf,
    )

    out_csv = Path(output_path) if output_path else project_root / "vina_pocket_cutoff_sweep.csv"
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=SWEEP_FIELDS)
        writer.writeheader()
        writer.writerows(results)

    return out_csv


def print_sweep_summary(results: List[dict]):
    """Print a human-readable sweep summary to stdout."""
    if not results:
        print("No results.")
        return

    receptors = sorted({r["receptor_id"] for r in results})
    cutoffs = sorted({r["cutoff"] for r in results})

    print(f"\n{'cutoff':>8}", end="")
    for rec in receptors:
        print(f"  {rec:>16}", end="")
    print()
    print("-" * (8 + 18 * len(receptors)))

    lookup = {(r["cutoff"], r["receptor_id"]): r for r in results}
    for c in cutoffs:
        print(f"{c:8.1f}", end="")
        for rec in receptors:
            r = lookup.get((c, rec))
            if r:
                print(f"  {r['n_pockets']:>4} pockets", end="")
                print(f" ({r['largest_pocket']:>3}~{r['smallest_pocket']:<3})", end="")
            else:
                print(f"  {'N/A':>16}", end="")
        print()


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def main():
    import argparse

    parser = argparse.ArgumentParser(description="Pocket cutoff sensitivity sweep")
    parser.add_argument("config", help="Project config YAML path")
    parser.add_argument("--min", type=float, default=2.0, help="Min cutoff (default: 2.0)")
    parser.add_argument("--max", type=float, default=8.0, help="Max cutoff (default: 8.0)")
    parser.add_argument("--step", type=float, default=0.5, help="Step size (default: 0.5)")
    parser.add_argument("-o", "--output", help="Output CSV path")

    args = parser.parse_args()

    print(f"Sweeping pocket_cutoff: {args.min} → {args.max} (step {args.step})")
    out = sweep_from_config(
        args.config,
        cutoff_min=args.min,
        cutoff_max=args.max,
        cutoff_step=args.step,
        output_path=args.output,
    )

    # Also print summary
    config = load_config(args.config)
    project_root = project_root_from_config(config)
    with open(out, newline="", encoding="utf-8") as f:
        results = list(csv.DictReader(f))
    # Convert types back
    for r in results:
        r["cutoff"] = float(r["cutoff"])
        r["n_pockets"] = int(r["n_pockets"])
        r["largest_pocket"] = int(r["largest_pocket"])
        r["smallest_pocket"] = int(r["smallest_pocket"])
    print_sweep_summary(results)
    print(f"\nSaved: {out}")


if __name__ == "__main__":
    main()
