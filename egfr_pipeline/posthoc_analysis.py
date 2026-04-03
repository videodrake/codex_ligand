# -*- coding: utf-8 -*-
"""Post-hoc analysis tools for validating PPI docking results.

Operates on existing output files (no re-docking required):
  1. Filter threshold sensitivity analysis
  2. Cross-seed convergence analysis
  3. Entropy-corrected binding energy

Usage:
  # Threshold sensitivity (single seed)
  python -m egfr_pipeline.posthoc_analysis sensitivity \
      --scored output/workflow_a/phase2_ppi_docking/3GT8_raw/prod_seed0/scored_all_models.csv

  # Cross-seed convergence (all seeds for one state)
  python -m egfr_pipeline.posthoc_analysis convergence \
      --state_dir output/workflow_a/phase2_ppi_docking/3GT8_raw/

  # Entropy correction (single seed)
  python -m egfr_pipeline.posthoc_analysis entropy \
      --ranking output/workflow_a/phase2_ppi_docking/3GT8_raw/prod_seed0/final_ranking.csv
"""

import argparse
import csv
import logging
import math
import os
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

logger = logging.getLogger("pipeline.posthoc")

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Rigid-body translational+rotational entropy loss upon PPI binding (REU).
# Estimated from Finkelstein & Janin (1989): ~2-3 kcal/mol for typical PPI.
# For a 52-residue partner (MYO1D beta-meander), conformational restriction
# adds ~0.5-1.5 kcal/mol additional penalty.  Total ~2.5 REU conservative.
RIGID_BODY_ENTROPY_PENALTY_REU = 2.5


# ---------------------------------------------------------------------------
# 1. Filter Threshold Sensitivity Analysis
# ---------------------------------------------------------------------------

SENSITIVITY_FIELDS = [
    "parameter", "original_value", "variation_pct", "varied_value",
    "n_survivors", "n_survivors_original", "delta_n",
    "survivor_overlap_frac", "top20_overlap_frac",
    "mean_dG_survivors", "median_dG_survivors",
]


def _load_scored_models(csv_path: str) -> List[Dict[str, Any]]:
    """Load scored_all_models.csv with safe numeric parsing."""
    rows = []
    with open(csv_path, "r", newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            parsed = {}
            for k, v in row.items():
                try:
                    parsed[k] = float(v)
                except (ValueError, TypeError):
                    parsed[k] = v
            rows.append(parsed)
    return rows


def _apply_filter(
    models: List[Dict],
    dG_max: float,
    dSASA_min: float,
    sc_min: float,
    total_score_pct: float,
) -> List[Dict]:
    """Apply Stage 1+2 cheap filters to models."""
    # Stage 1a: dG cutoff
    passed = [m for m in models if m.get("dG_separated", 0) <= dG_max]
    if not passed:
        return []

    # Stage 1b: total_score percentile
    scores = sorted([m.get("total_score", 0) for m in passed])
    n = len(scores)
    idx = max(0, int(n * total_score_pct / 100.0) - 1)
    cutoff = scores[idx]
    passed = [m for m in passed if m.get("total_score", 0) <= cutoff]

    # Stage 2 cheap: dSASA + sc
    passed = [m for m in passed
              if m.get("dSASA", 0) >= dSASA_min
              and m.get("sc_value", 0) >= sc_min]

    return passed


def _model_ids(models: List[Dict]) -> set:
    return {m.get("id") or m.get("model_id") or i for i, m in enumerate(models)}


def run_sensitivity_analysis(
    scored_csv: str,
    output_path: str,
    variations: Optional[List[int]] = None,
) -> List[Dict]:
    """Vary each filter threshold by ±N% and measure survivor stability.

    Args:
        scored_csv: Path to scored_all_models.csv
        output_path: Output CSV path
        variations: Percentage variations to test (default: [-20, -10, +10, +20])
    """
    if variations is None:
        variations = [-20, -10, 10, 20]

    models = _load_scored_models(scored_csv)
    if not models:
        logger.error(f"No models in {scored_csv}")
        return []

    logger.info(f"Loaded {len(models)} models from {scored_csv}")

    # Original thresholds (from production config defaults)
    orig = {
        "dG_max": 0.0,
        "dSASA_min": 500.0,
        "sc_min": 0.50,
        "total_score_pct": 10.0,
    }

    # Baseline survivors
    baseline = _apply_filter(models, **orig)
    baseline_ids = _model_ids(baseline)
    baseline_top20 = _model_ids(sorted(baseline, key=lambda m: m.get("dG_separated", 0))[:20])
    n_baseline = len(baseline)

    logger.info(f"Baseline: {n_baseline} survivors")

    results = []

    for param_name, orig_val in orig.items():
        for var_pct in variations:
            # Apply variation
            if param_name in ("dG_max",):
                # dG_max is 0.0; vary as absolute offset
                varied_val = orig_val + var_pct * 0.1  # ±2, ±1 REU
            elif param_name == "total_score_pct":
                varied_val = max(1, min(50, orig_val + var_pct * orig_val / 100))
            else:
                varied_val = orig_val * (1 + var_pct / 100.0)

            # Run filter with varied threshold
            test_params = dict(orig)
            test_params[param_name] = varied_val
            test_survivors = _apply_filter(models, **test_params)
            test_ids = _model_ids(test_survivors)

            # Overlap analysis
            overlap = len(baseline_ids & test_ids)
            overlap_frac = overlap / max(1, len(baseline_ids | test_ids))

            test_top20 = _model_ids(sorted(test_survivors, key=lambda m: m.get("dG_separated", 0))[:20])
            top20_overlap = len(baseline_top20 & test_top20)
            top20_frac = top20_overlap / max(1, len(baseline_top20))

            # Stats
            dGs = [m.get("dG_separated", 0) for m in test_survivors]
            mean_dG = sum(dGs) / len(dGs) if dGs else 0
            sorted_dGs = sorted(dGs)
            median_dG = sorted_dGs[len(sorted_dGs) // 2] if sorted_dGs else 0

            results.append({
                "parameter": param_name,
                "original_value": round(orig_val, 3),
                "variation_pct": var_pct,
                "varied_value": round(varied_val, 3),
                "n_survivors": len(test_survivors),
                "n_survivors_original": n_baseline,
                "delta_n": len(test_survivors) - n_baseline,
                "survivor_overlap_frac": round(overlap_frac, 4),
                "top20_overlap_frac": round(top20_frac, 4),
                "mean_dG_survivors": round(mean_dG, 2),
                "median_dG_survivors": round(median_dG, 2),
            })

    # Write results
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=SENSITIVITY_FIELDS)
        writer.writeheader()
        writer.writerows(results)

    _print_sensitivity_summary(results, n_baseline)
    logger.info(f"Sensitivity analysis: {out}")
    return results


def _print_sensitivity_summary(results: List[Dict], n_baseline: int) -> None:
    print(f"\n{'='*70}")
    print(f"  Filter Threshold Sensitivity Analysis")
    print(f"  Baseline survivors: {n_baseline}")
    print(f"{'='*70}")

    # Group by parameter
    by_param: Dict[str, List[Dict]] = defaultdict(list)
    for r in results:
        by_param[r["parameter"]].append(r)

    for param, rows in by_param.items():
        print(f"\n  {param} (original={rows[0]['original_value']}):")
        for r in rows:
            stability = "STABLE" if r["top20_overlap_frac"] >= 0.80 else \
                        "MODERATE" if r["top20_overlap_frac"] >= 0.50 else "UNSTABLE"
            print(f"    {r['variation_pct']:+d}% → {r['varied_value']:.3f}  "
                  f"survivors={r['n_survivors']:>5d} ({r['delta_n']:+d})  "
                  f"top20_overlap={r['top20_overlap_frac']:.0%}  [{stability}]")

    # Overall verdict
    top20_fracs = [r["top20_overlap_frac"] for r in results]
    min_overlap = min(top20_fracs)
    mean_overlap = sum(top20_fracs) / len(top20_fracs)
    print(f"\n  Overall: min_top20_overlap={min_overlap:.0%}, "
          f"mean_top20_overlap={mean_overlap:.0%}")
    if mean_overlap >= 0.75:
        print(f"  → Filtering is ROBUST to ±20% threshold variation")
    elif mean_overlap >= 0.50:
        print(f"  → Filtering is MODERATELY SENSITIVE — consider tightening criteria")
    else:
        print(f"  → WARNING: Filtering is HIGHLY SENSITIVE — results may not be reproducible")
    print(f"{'='*70}\n")


# ---------------------------------------------------------------------------
# 2. Cross-Seed Convergence Analysis
# ---------------------------------------------------------------------------

CONVERGENCE_FIELDS = [
    "state", "n_seeds", "n_clusters_per_seed_mean", "n_clusters_per_seed_std",
    "n_consensus_sites", "consensus_site_list",
    "centroid_agreement_mean_A", "centroid_agreement_max_A",
    "dG_best_mean", "dG_best_std", "dG_best_range",
    "convergence_verdict",
]


def _load_cluster_summaries(state_dir: str) -> Dict[str, List[Dict]]:
    """Load cluster_summary.csv from each seed directory.

    Returns {seed_name: [cluster_rows]}.
    """
    root = Path(state_dir)
    seed_data = {}
    for seed_dir in sorted(root.glob("prod_seed*")):
        csv_path = seed_dir / "cluster_results" / "cluster_summary.csv"
        if not csv_path.exists():
            continue
        rows = []
        with open(csv_path, "r", newline="", encoding="utf-8") as f:
            for row in csv.DictReader(f):
                parsed = {}
                for k, v in row.items():
                    try:
                        parsed[k] = float(v)
                    except (ValueError, TypeError):
                        parsed[k] = v
                rows.append(parsed)
        seed_data[seed_dir.name] = rows
    return seed_data


def _cluster_centroids(clusters: List[Dict]) -> List[Tuple[float, float, float, float]]:
    """Extract (x, y, z, dG) from cluster rows."""
    centroids = []
    for c in clusters:
        try:
            x = float(c.get("centroid_x") or c.get("center_x") or 0)
            y = float(c.get("centroid_y") or c.get("center_y") or 0)
            z = float(c.get("centroid_z") or c.get("center_z") or 0)
            dG = float(c.get("best_dG") or c.get("dG_separated") or 0)
            centroids.append((x, y, z, dG))
        except (ValueError, TypeError):
            continue
    return centroids


def _euclidean(a: Tuple[float, ...], b: Tuple[float, ...]) -> float:
    return math.sqrt(sum((ai - bi) ** 2 for ai, bi in zip(a[:3], b[:3])))


def run_convergence_analysis(
    state_dir: str,
    output_path: str,
    match_distance: float = 8.0,
) -> Dict[str, Any]:
    """Analyse cross-seed cluster convergence.

    Two clusters from different seeds are considered "matching" if their
    centroids are within match_distance Angstroms.

    Args:
        state_dir: Directory containing prod_seed*/ subdirs.
        output_path: Output CSV path.
        match_distance: Maximum centroid distance for cluster matching (Å).
    """
    seed_data = _load_cluster_summaries(state_dir)
    n_seeds = len(seed_data)
    state_name = Path(state_dir).name

    if n_seeds < 2:
        logger.warning(f"Need ≥2 seeds for convergence analysis, found {n_seeds}")
        return {}

    logger.info(f"Loaded cluster data from {n_seeds} seeds for {state_name}")

    # Extract centroids per seed
    seed_centroids: Dict[str, List[Tuple[float, float, float, float]]] = {}
    for seed_name, clusters in seed_data.items():
        seed_centroids[seed_name] = _cluster_centroids(clusters)

    # Find consensus sites: centroids that appear in ≥2 seeds
    all_centroids = []
    for seed_name, centroids in seed_centroids.items():
        for c in centroids:
            all_centroids.append((seed_name, c))

    # Greedy matching across seeds
    consensus_groups: List[List[Tuple[str, Tuple]]] = []
    used = set()

    for i, (s1, c1) in enumerate(all_centroids):
        if i in used:
            continue
        group = [(s1, c1)]
        used.add(i)
        for j, (s2, c2) in enumerate(all_centroids):
            if j in used or s2 == s1:
                continue
            # Check if already have this seed in group
            if any(gs == s2 for gs, _ in group):
                continue
            if _euclidean(c1, c2) <= match_distance:
                group.append((s2, c2))
                used.add(j)
        if len(group) >= 2:
            consensus_groups.append(group)

    # Compute statistics
    n_clusters_list = [len(v) for v in seed_centroids.values()]
    n_clusters_mean = sum(n_clusters_list) / len(n_clusters_list)
    n_clusters_std = math.sqrt(
        sum((x - n_clusters_mean) ** 2 for x in n_clusters_list) / max(1, len(n_clusters_list) - 1)
    ) if len(n_clusters_list) > 1 else 0

    # Best dG per seed
    best_dGs = []
    for centroids in seed_centroids.values():
        if centroids:
            best_dGs.append(min(c[3] for c in centroids))

    dG_mean = sum(best_dGs) / len(best_dGs) if best_dGs else 0
    dG_std = math.sqrt(
        sum((x - dG_mean) ** 2 for x in best_dGs) / max(1, len(best_dGs) - 1)
    ) if len(best_dGs) > 1 else 0
    dG_range = max(best_dGs) - min(best_dGs) if best_dGs else 0

    # Consensus site centroid agreement
    centroid_dists = []
    consensus_labels = []
    for gi, group in enumerate(consensus_groups):
        coords = [c for _, c in group]
        # Mean centroid
        mean_x = sum(c[0] for c in coords) / len(coords)
        mean_y = sum(c[1] for c in coords) / len(coords)
        mean_z = sum(c[2] for c in coords) / len(coords)
        for _, c in group:
            d = math.sqrt((c[0] - mean_x)**2 + (c[1] - mean_y)**2 + (c[2] - mean_z)**2)
            centroid_dists.append(d)
        seeds_in_group = len(group)
        consensus_labels.append(f"CS{gi+1}({seeds_in_group}/{n_seeds}seeds)")

    cent_mean = sum(centroid_dists) / len(centroid_dists) if centroid_dists else 0
    cent_max = max(centroid_dists) if centroid_dists else 0

    # Verdict
    if len(consensus_groups) >= 2 and cent_mean < 5.0 and dG_std < 3.0:
        verdict = "CONVERGED"
    elif len(consensus_groups) >= 1 and cent_mean < 8.0:
        verdict = "PARTIALLY_CONVERGED"
    else:
        verdict = "NOT_CONVERGED"

    result = {
        "state": state_name,
        "n_seeds": n_seeds,
        "n_clusters_per_seed_mean": round(n_clusters_mean, 1),
        "n_clusters_per_seed_std": round(n_clusters_std, 1),
        "n_consensus_sites": len(consensus_groups),
        "consensus_site_list": "; ".join(consensus_labels),
        "centroid_agreement_mean_A": round(cent_mean, 2),
        "centroid_agreement_max_A": round(cent_max, 2),
        "dG_best_mean": round(dG_mean, 2),
        "dG_best_std": round(dG_std, 2),
        "dG_best_range": round(dG_range, 2),
        "convergence_verdict": verdict,
    }

    # Write
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=CONVERGENCE_FIELDS)
        writer.writeheader()
        writer.writerow(result)

    _print_convergence_summary(result)
    return result


def _print_convergence_summary(r: Dict) -> None:
    print(f"\n{'='*70}")
    print(f"  Cross-Seed Convergence Analysis: {r['state']}")
    print(f"{'='*70}")
    print(f"  Seeds analysed:          {r['n_seeds']}")
    print(f"  Clusters per seed:       {r['n_clusters_per_seed_mean']:.1f} ± {r['n_clusters_per_seed_std']:.1f}")
    print(f"  Consensus sites (≥2 seeds): {r['n_consensus_sites']}")
    if r['consensus_site_list']:
        print(f"    {r['consensus_site_list']}")
    print(f"  Centroid agreement:      {r['centroid_agreement_mean_A']:.1f} Å mean, "
          f"{r['centroid_agreement_max_A']:.1f} Å max")
    print(f"  Best dG across seeds:    {r['dG_best_mean']:.1f} ± {r['dG_best_std']:.1f} REU "
          f"(range {r['dG_best_range']:.1f})")
    print(f"  Verdict:                 {r['convergence_verdict']}")
    if r["convergence_verdict"] == "CONVERGED":
        print(f"  → Seeds converge to consistent binding sites (good sampling)")
    elif r["convergence_verdict"] == "PARTIALLY_CONVERGED":
        print(f"  → Some consensus sites found, but sampling may be insufficient")
    else:
        print(f"  → WARNING: Seeds do not converge — consider increasing decoys or seeds")
    print(f"{'='*70}\n")


# ---------------------------------------------------------------------------
# 3. Entropy-Corrected Binding Energy
# ---------------------------------------------------------------------------

def apply_entropy_correction(
    ranking_csv: str,
    output_path: str,
    penalty_reu: float = RIGID_BODY_ENTROPY_PENALTY_REU,
) -> List[Dict]:
    """Add entropy-corrected dG column to final_ranking.csv.

    dG_corrected = dG_separated + penalty_reu

    The penalty accounts for rigid-body translational+rotational entropy
    loss upon PPI complex formation, which ref2015 does not include.

    Args:
        ranking_csv: Path to final_ranking.csv
        output_path: Output CSV path (same schema + dG_corrected column)
        penalty_reu: Entropy penalty in REU (positive value, added to dG)
    """
    rows = []
    fieldnames = None
    with open(ranking_csv, "r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        fieldnames = list(reader.fieldnames or [])
        for row in reader:
            rows.append(row)

    if not rows:
        logger.warning(f"No rows in {ranking_csv}")
        return []

    # Add corrected columns
    new_fields = []
    if "dG_corrected" not in fieldnames:
        new_fields.append("dG_corrected")
    if "dG_density_corrected" not in fieldnames:
        new_fields.append("dG_density_corrected")
    if "rank_corrected" not in fieldnames:
        new_fields.append("rank_corrected")
    fieldnames.extend(new_fields)

    for row in rows:
        try:
            dG = float(row.get("dG_separated", 0))
        except (ValueError, TypeError):
            dG = 0
        dG_corr = dG + penalty_reu
        row["dG_corrected"] = round(dG_corr, 3)

        try:
            dSASA = float(row.get("dSASA", 0))
        except (ValueError, TypeError):
            dSASA = 0
        row["dG_density_corrected"] = round(
            dG_corr / dSASA * 100 if abs(dSASA) > 1e-6 else 0, 4
        )

    # Re-rank by corrected dG
    sortable = [(i, float(r.get("dG_corrected", 0))) for i, r in enumerate(rows)]
    sortable.sort(key=lambda x: x[1])
    for new_rank, (orig_idx, _) in enumerate(sortable, 1):
        rows[orig_idx]["rank_corrected"] = new_rank

    # Write
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    # Check for rank changes
    rank_changes = 0
    for row in rows:
        orig_rank = row.get("rank") or row.get("Rank")
        if orig_rank and str(orig_rank) != str(row.get("rank_corrected", "")):
            rank_changes += 1

    print(f"\n  Entropy correction applied: +{penalty_reu} REU")
    print(f"  Rank changes: {rank_changes}/{len(rows)} models")
    print(f"  Output: {out}\n")

    return rows


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Post-hoc analysis of PPI docking results (no re-docking)"
    )
    subparsers = parser.add_subparsers(dest="command")

    # Sensitivity
    sens = subparsers.add_parser("sensitivity",
        help="Filter threshold sensitivity analysis")
    sens.add_argument("--scored", required=True,
        help="scored_all_models.csv path")
    sens.add_argument("--output", default="",
        help="Output CSV (default: alongside input)")

    # Convergence
    conv = subparsers.add_parser("convergence",
        help="Cross-seed convergence analysis")
    conv.add_argument("--state_dir", required=True,
        help="State directory containing prod_seed*/ subdirs")
    conv.add_argument("--output", default="",
        help="Output CSV (default: in state_dir)")
    conv.add_argument("--match_distance", type=float, default=8.0,
        help="Max centroid distance for cluster matching (Å)")

    # Entropy
    ent = subparsers.add_parser("entropy",
        help="Apply entropy correction to final_ranking.csv")
    ent.add_argument("--ranking", required=True,
        help="final_ranking.csv path")
    ent.add_argument("--output", default="",
        help="Output CSV (default: alongside input)")
    ent.add_argument("--penalty", type=float, default=RIGID_BODY_ENTROPY_PENALTY_REU,
        help=f"Entropy penalty in REU (default: {RIGID_BODY_ENTROPY_PENALTY_REU})")

    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    if args.command == "sensitivity":
        output = args.output or str(
            Path(args.scored).parent / "sensitivity_analysis.csv")
        run_sensitivity_analysis(args.scored, output)

    elif args.command == "convergence":
        output = args.output or str(
            Path(args.state_dir) / "convergence_analysis.csv")
        run_convergence_analysis(args.state_dir, output, args.match_distance)

    elif args.command == "entropy":
        output = args.output or str(
            Path(args.ranking).parent / "final_ranking_entropy_corrected.csv")
        apply_entropy_correction(args.ranking, output, args.penalty)

    else:
        parser.print_help()


if __name__ == "__main__":
    main()
