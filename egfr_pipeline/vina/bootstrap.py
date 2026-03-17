#!/usr/bin/env python3
"""Bootstrap resampling for pocket clustering stability analysis.

Re-samples Vina poses per receptor, re-clusters, and measures how
consistently each pocket is recovered across replicates.

Usage:
    python -m egfr_pipeline.vina.bootstrap config/example-project.yaml
    python -m egfr_pipeline.vina.bootstrap config/example-project.yaml --n 200 --seed 123
"""
import copy
import csv
import math
import random
import statistics
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

from egfr_pipeline.config import load_config, project_root_from_config
from egfr_pipeline.vina.cluster import (
    assign_pockets,
    euclidean_distance,
    load_pose_table,
    merge_pockets_by_residue,
)
from egfr_pipeline.vina.summarize import summarize_pose_rows


# ---------------------------------------------------------------------------
# A1: Core bootstrap function
# ---------------------------------------------------------------------------

def bootstrap_once(
    rows: List[dict],
    pocket_cutoff: float,
    rng: random.Random,
    sample_fraction: float = 0.8,
    max_iterations: int = 10,
    min_pocket_size: int = 2,
    merge_by_residue: bool = True,
    merge_jaccard: float = 0.3,
    merge_overlap: float = 0.5,
    merge_centroid_fallback: float = 6.0,
) -> List[dict]:
    """Single bootstrap replicate: subsample poses, re-cluster.

    Samples ``sample_fraction`` of poses per receptor (with replacement),
    strips existing pocket_id, and re-clusters via assign_pockets().

    Returns clustered rows (same schema as vina_pose_table.csv rows).
    """
    # Group by receptor
    by_receptor: Dict[str, List[dict]] = defaultdict(list)
    for row in rows:
        by_receptor[row["receptor_id"]].append(row)

    sampled: List[dict] = []
    for receptor_id, rec_rows in by_receptor.items():
        k = max(1, int(len(rec_rows) * sample_fraction))
        chosen = rng.choices(rec_rows, k=k)
        for row in chosen:
            r = copy.copy(row)
            r.pop("pocket_id", None)
            sampled.append(r)

    clustered = assign_pockets(
        sampled,
        pocket_cutoff,
        max_iterations=max_iterations,
        min_pocket_size=min_pocket_size,
    )
    if merge_by_residue:
        clustered, _merge_log = merge_pockets_by_residue(
            clustered,
            jaccard_threshold=merge_jaccard,
            overlap_threshold=merge_overlap,
            centroid_fallback_cutoff=merge_centroid_fallback,
        )
    return clustered


# ---------------------------------------------------------------------------
# A2: Multi-replicate runner
# ---------------------------------------------------------------------------

def run_bootstrap_replicates(
    pose_table_path: str,
    n_replicates: int = 100,
    pocket_cutoff: float = 8.0,
    sample_fraction: float = 0.8,
    seed: int = 42,
    max_iterations: int = 10,
    min_pocket_size: int = 2,
    merge_by_residue: bool = True,
    merge_jaccard: float = 0.3,
    merge_overlap: float = 0.5,
    merge_centroid_fallback: float = 6.0,
) -> Tuple[List[dict], List[List[dict]]]:
    """Run N bootstrap replicates.

    Returns:
        reference_pockets: pocket summary rows from the full (non-resampled) data
        replicate_summaries: list of N pocket summary row lists
    """
    base_rows = load_pose_table(Path(pose_table_path))

    # Reference: cluster full data
    ref_rows = copy.deepcopy(base_rows)
    for r in ref_rows:
        r.pop("pocket_id", None)
    ref_clustered = assign_pockets(
        ref_rows,
        pocket_cutoff,
        max_iterations=max_iterations,
        min_pocket_size=min_pocket_size,
    )
    if merge_by_residue:
        ref_clustered, _merge_log = merge_pockets_by_residue(
            ref_clustered,
            jaccard_threshold=merge_jaccard,
            overlap_threshold=merge_overlap,
            centroid_fallback_cutoff=merge_centroid_fallback,
        )
    ref_pockets, _, _ = summarize_pose_rows(ref_clustered)

    # Replicates
    replicate_summaries: List[List[dict]] = []
    for i in range(n_replicates):
        rng = random.Random(seed + i)
        rep_rows = bootstrap_once(
            base_rows, pocket_cutoff, rng,
            sample_fraction=sample_fraction,
            max_iterations=max_iterations,
            min_pocket_size=min_pocket_size,
            merge_by_residue=merge_by_residue,
            merge_jaccard=merge_jaccard,
            merge_overlap=merge_overlap,
            merge_centroid_fallback=merge_centroid_fallback,
        )
        rep_pockets, _, _ = summarize_pose_rows(rep_rows)
        replicate_summaries.append(rep_pockets)

    return ref_pockets, replicate_summaries


# ---------------------------------------------------------------------------
# A3: Pocket matching across replicates
# ---------------------------------------------------------------------------

def _pocket_centroid(pocket: dict) -> List[float]:
    return [float(pocket["centroid_x"]),
            float(pocket["centroid_y"]),
            float(pocket["centroid_z"])]


def match_pockets_to_reference(
    reference_pockets: List[dict],
    replicate_pockets: List[dict],
    centroid_cutoff: float = 8.0,
) -> Dict[Tuple[str, str], Optional[dict]]:
    """Match replicate pockets to reference pockets by centroid proximity.

    Returns {(receptor_id, ref_pocket_id): matched_replicate_pocket_or_None}.
    """
    # Group replicate pockets by receptor
    rep_by_receptor: Dict[str, List[dict]] = defaultdict(list)
    for p in replicate_pockets:
        rep_by_receptor[p["receptor_id"]].append(p)

    result: Dict[Tuple[str, str], Optional[dict]] = {}
    used: Set[Tuple[str, str]] = set()

    for ref in reference_pockets:
        rec_id = ref["receptor_id"]
        ref_c = _pocket_centroid(ref)
        key = (rec_id, ref["pocket_id"])

        best_match = None
        best_dist = centroid_cutoff + 1.0

        for rep in rep_by_receptor.get(rec_id, []):
            rep_key = (rec_id, rep["pocket_id"])
            if rep_key in used:
                continue
            dist = euclidean_distance(ref_c, _pocket_centroid(rep))
            if dist <= centroid_cutoff and dist < best_dist:
                best_match = rep
                best_dist = dist

        if best_match is not None:
            used.add((rec_id, best_match["pocket_id"]))

        result[key] = best_match

    return result


# ---------------------------------------------------------------------------
# A4: Aggregate bootstrap statistics
# ---------------------------------------------------------------------------

BOOTSTRAP_FIELDS = [
    "receptor_id",
    "pocket_id",
    "pocket_exists_frac",
    "centroid_std_A",
    "affinity_mean",
    "affinity_std",
    "affinity_iqr",
    "n_pose_mean",
    "n_pose_std",
    "n_replicates",
    "sample_fraction",
    "stability_scope",
]


def compute_bootstrap_stats(
    reference_pockets: List[dict],
    replicate_summaries: List[List[dict]],
    centroid_cutoff: float = 8.0,
) -> List[dict]:
    """Compute per-pocket stability statistics across bootstrap replicates."""
    n_rep = len(replicate_summaries)
    results: List[dict] = []

    for ref in reference_pockets:
        key = (ref["receptor_id"], ref["pocket_id"])
        ref_c = _pocket_centroid(ref)

        matched_centroids: List[List[float]] = []
        matched_affinities: List[float] = []
        matched_nposes: List[int] = []

        for rep_pockets in replicate_summaries:
            mapping = match_pockets_to_reference(
                [ref], rep_pockets, centroid_cutoff,
            )
            m = mapping.get(key)
            if m is not None:
                matched_centroids.append(_pocket_centroid(m))
                aff = m.get("best_affinity", "")
                if aff not in ("", None):
                    matched_affinities.append(float(aff))
                np_val = m.get("n_pose", "")
                if np_val not in ("", None):
                    matched_nposes.append(int(np_val))

        exists_frac = len(matched_centroids) / n_rep if n_rep > 0 else 0.0

        # Centroid std
        if len(matched_centroids) >= 2:
            dists = [euclidean_distance(c, ref_c) for c in matched_centroids]
            centroid_std = statistics.stdev(dists)
        else:
            centroid_std = 0.0

        # Affinity stats
        aff_mean = statistics.mean(matched_affinities) if matched_affinities else 0.0
        aff_std = statistics.stdev(matched_affinities) if len(matched_affinities) >= 2 else 0.0
        if len(matched_affinities) >= 4:
            s = sorted(matched_affinities)
            q1 = s[len(s) // 4]
            q3 = s[3 * len(s) // 4]
            aff_iqr = q3 - q1
        else:
            aff_iqr = 0.0

        # Pose count stats
        np_mean = statistics.mean(matched_nposes) if matched_nposes else 0.0
        np_std = statistics.stdev(matched_nposes) if len(matched_nposes) >= 2 else 0.0

        results.append({
            "receptor_id": ref["receptor_id"],
            "pocket_id": ref["pocket_id"],
            "pocket_exists_frac": round(exists_frac, 4),
            "centroid_std_A": round(centroid_std, 4),
            "affinity_mean": round(aff_mean, 4),
            "affinity_std": round(aff_std, 4),
            "affinity_iqr": round(aff_iqr, 4),
            "n_pose_mean": round(np_mean, 2),
            "n_pose_std": round(np_std, 2),
            "n_replicates": n_rep,
        })

    return results


# ---------------------------------------------------------------------------
# A5: Config integration and CLI entry point
# ---------------------------------------------------------------------------

def bootstrap_from_config(
    config_path: str,
    n_replicates: Optional[int] = None,
    seed: Optional[int] = None,
    output_path: Optional[str] = None,
) -> Path:
    """Run bootstrap analysis from config, save vina_pocket_bootstrap.csv."""
    config = load_config(config_path)
    project_root = project_root_from_config(config)
    pose_table = project_root / "vina_pose_table.csv"

    if not pose_table.exists():
        raise FileNotFoundError(f"Pose table not found: {pose_table}")

    pp = config.get("postprocess", {})
    bs = config.get("bootstrap", {})

    n = n_replicates or bs.get("n_replicates", 100)
    s = seed or bs.get("seed", 42)
    cutoff = pp.get("pocket_cutoff", 8.0)
    fraction = bs.get("sample_fraction", 0.8)
    max_iterations = pp.get("cluster_max_iterations", 10)
    min_pocket_size = pp.get("min_pocket_size", 2)
    merge = pp.get("merge_by_residue", True)
    merge_j = pp.get("merge_jaccard", 0.3)
    merge_oc = pp.get("merge_overlap", 0.5)
    merge_cf = pp.get("merge_centroid_fallback", 6.0)

    ref_pockets, rep_summaries = run_bootstrap_replicates(
        str(pose_table),
        n_replicates=n,
        pocket_cutoff=cutoff,
        sample_fraction=fraction,
        seed=s,
        max_iterations=max_iterations,
        min_pocket_size=min_pocket_size,
        merge_by_residue=merge,
        merge_jaccard=merge_j,
        merge_overlap=merge_oc,
        merge_centroid_fallback=merge_cf,
    )

    stats = compute_bootstrap_stats(ref_pockets, rep_summaries, centroid_cutoff=cutoff)
    for row in stats:
        row["sample_fraction"] = round(float(fraction), 4)
        row["stability_scope"] = "pose_resampling"

    out_csv = Path(output_path) if output_path else project_root / "vina_pocket_bootstrap.csv"
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=BOOTSTRAP_FIELDS)
        writer.writeheader()
        writer.writerows(stats)

    return out_csv


def print_bootstrap_summary(stats: List[dict]):
    """Print human-readable bootstrap summary."""
    if not stats:
        print("No results.")
        return

    print(f"\n{'receptor':<20} {'pocket':<8} {'stability':>10} "
          f"{'centroid_std':>12} {'aff_mean':>10} {'aff_std':>9}")
    print("-" * 75)
    for s in stats:
        print(f"{s['receptor_id']:<20} {s['pocket_id']:<8} "
              f"{s['pocket_exists_frac']:>10.2%} "
              f"{s['centroid_std_A']:>10.2f} A "
              f"{s['affinity_mean']:>10.2f} "
              f"{s['affinity_std']:>8.2f}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    import argparse

    parser = argparse.ArgumentParser(description="Bootstrap pocket stability analysis")
    parser.add_argument("config", help="Project config YAML path")
    parser.add_argument("--n", type=int, default=None, help="Number of replicates (default: 100)")
    parser.add_argument("--seed", type=int, default=None, help="Random seed (default: 42)")
    parser.add_argument("-o", "--output", help="Output CSV path")

    args = parser.parse_args()

    n = args.n or 100
    print(f"Bootstrap pocket stability: {n} replicates, seed={args.seed or 42}")
    out = bootstrap_from_config(
        args.config,
        n_replicates=args.n,
        seed=args.seed,
        output_path=args.output,
    )

    # Print summary
    with open(out, newline="", encoding="utf-8") as f:
        stats = list(csv.DictReader(f))
    for s in stats:
        s["pocket_exists_frac"] = float(s["pocket_exists_frac"])
        s["centroid_std_A"] = float(s["centroid_std_A"])
        s["affinity_mean"] = float(s["affinity_mean"])
        s["affinity_std"] = float(s["affinity_std"])
    print_bootstrap_summary(stats)
    print(f"\nSaved: {out}")


if __name__ == "__main__":
    main()
