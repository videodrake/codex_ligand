#!/usr/bin/env python3
"""Analyze C-lobe surface pocket affinity distribution and evaluate thresholds.

AC-2.4: Assess whether current affinity thresholds (-8.0/-6.5 kcal/mol)
provide sufficient discriminative power for C-lobe surface pockets.

Usage:
    python scripts/analyze_affinity_distribution.py <pocket_table_csv>
    python scripts/analyze_affinity_distribution.py --config <config.yaml>
"""

import argparse
import csv
import statistics
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple


ANALYSIS_FIELDS = [
    "metric",
    "value",
]


def load_pocket_table(path: Path) -> List[dict]:
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def analyze_affinity_distribution(
    pocket_rows: List[dict],
    thresholds: Tuple[float, float] = (-8.0, -6.5),
) -> Dict[str, object]:
    """Analyze affinity distribution for non-ATP-site pockets.

    Returns analysis dict with percentiles, bracket counts, and
    discriminative verdict ("sufficient" or "insufficient").
    """
    affinities: List[float] = []
    for row in pocket_rows:
        is_atp = row.get("is_atp_site", "")
        if str(is_atp).lower() in ("true", "1", "yes"):
            continue
        aff_str = row.get("best_affinity", "")
        if aff_str in ("", None):
            continue
        try:
            affinities.append(float(aff_str))
        except (TypeError, ValueError):
            continue

    if not affinities:
        return {"n_pockets": 0, "discriminative_power": "no_data"}

    affinities.sort()
    n = len(affinities)

    def _percentile(data: List[float], p: float) -> float:
        k = (len(data) - 1) * (p / 100.0)
        f = int(k)
        c = min(f + 1, len(data) - 1)
        return data[f] + (k - f) * (data[c] - data[f])

    p25 = _percentile(affinities, 25)
    p50 = _percentile(affinities, 50)
    p75 = _percentile(affinities, 75)
    p90 = _percentile(affinities, 90)

    # Bracket analysis: strong (<= -8.0), good (-8.0 to -6.5], weak (> -6.5)
    t_strong, t_good = sorted(thresholds)
    bracket_strong = sum(1 for a in affinities if a <= t_strong)
    bracket_good = sum(1 for a in affinities if t_strong < a <= t_good)
    bracket_weak = sum(1 for a in affinities if a > t_good)

    pct_strong = round(bracket_strong / n * 100, 1)
    pct_good = round(bracket_good / n * 100, 1)
    pct_weak = round(bracket_weak / n * 100, 1)

    max_pct = max(pct_strong, pct_good, pct_weak)
    dominant = (
        "strong" if pct_strong == max_pct
        else "good" if pct_good == max_pct
        else "weak"
    )
    discriminative = "sufficient" if max_pct < 70 else "insufficient"

    recommendation = "maintain_current_thresholds"
    if discriminative == "insufficient":
        if dominant == "weak":
            recommendation = (
                "consider_adding_intermediate_threshold_around_"
                f"{round((t_good + p50) / 2, 1)}"
            )
        elif dominant == "strong":
            recommendation = (
                "consider_splitting_strong_bracket_around_"
                f"{round((t_strong + p50) / 2, 1)}"
            )
        else:
            recommendation = "consider_adjusting_threshold_boundaries"

    return {
        "n_pockets": n,
        "min_affinity": round(min(affinities), 4),
        "max_affinity": round(max(affinities), 4),
        "mean_affinity": round(statistics.mean(affinities), 4),
        "std_affinity": round(statistics.stdev(affinities), 4) if n >= 2 else 0.0,
        "p25": round(p25, 4),
        "p50_median": round(p50, 4),
        "p75": round(p75, 4),
        "p90": round(p90, 4),
        "threshold_strong": t_strong,
        "threshold_good": t_good,
        "n_bracket_strong": bracket_strong,
        "n_bracket_good": bracket_good,
        "n_bracket_weak": bracket_weak,
        "pct_bracket_strong": pct_strong,
        "pct_bracket_good": pct_good,
        "pct_bracket_weak": pct_weak,
        "max_bracket_pct": max_pct,
        "dominant_bracket": dominant,
        "discriminative_power": discriminative,
        "recommendation": recommendation,
    }


def write_analysis_csv(path: Path, analysis: Dict[str, object]) -> Path:
    """Write analysis results as metric/value pairs."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=ANALYSIS_FIELDS)
        writer.writeheader()
        for metric, value in analysis.items():
            writer.writerow({"metric": metric, "value": value})
    return path


def format_summary(analysis: Dict[str, object]) -> str:
    """Format a human-readable summary of the analysis."""
    if analysis.get("n_pockets", 0) == 0:
        return "No non-ATP-site pockets found for analysis."

    lines = [
        "C-lobe Surface Pocket Affinity Distribution Analysis",
        "=" * 55,
        f"  Pockets analyzed (is_atp_site=False): {analysis['n_pockets']}",
        f"  Range: {analysis['min_affinity']} to {analysis['max_affinity']} kcal/mol",
        f"  Mean:  {analysis['mean_affinity']} +/- {analysis['std_affinity']}",
        "",
        "  Percentiles:",
        f"    25th: {analysis['p25']} kcal/mol",
        f"    50th: {analysis['p50_median']} kcal/mol (median)",
        f"    75th: {analysis['p75']} kcal/mol",
        f"    90th: {analysis['p90']} kcal/mol",
        "",
        f"  Current thresholds: {analysis['threshold_strong']} / {analysis['threshold_good']}",
        f"  Bracket distribution:",
        f"    Strong (<= {analysis['threshold_strong']}): "
        f"{analysis['n_bracket_strong']} ({analysis['pct_bracket_strong']}%)",
        f"    Good   ({analysis['threshold_strong']} to {analysis['threshold_good']}]: "
        f"{analysis['n_bracket_good']} ({analysis['pct_bracket_good']}%)",
        f"    Weak   (> {analysis['threshold_good']}): "
        f"{analysis['n_bracket_weak']} ({analysis['pct_bracket_weak']}%)",
        "",
        f"  Discriminative power: {analysis['discriminative_power'].upper()}",
        f"  Dominant bracket: {analysis['dominant_bracket']} ({analysis['max_bracket_pct']}%)",
        f"  Recommendation: {analysis['recommendation']}",
    ]
    return "\n".join(lines)


def main(argv: Optional[List[str]] = None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "pocket_table",
        nargs="?",
        help="Path to vina_pocket_table.csv",
    )
    parser.add_argument(
        "--config",
        help="Path to project YAML config (alternative to direct CSV path)",
    )
    parser.add_argument(
        "--output",
        help="Output CSV path (default: alongside input)",
    )
    args = parser.parse_args(argv)

    if args.config:
        sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
        from egfr_pipeline.config import load_config
        from egfr_pipeline import paths

        config = load_config(args.config)
        pocket_table_path = paths.wa_phase4_vina_postprocess(config) / "vina_pocket_table.csv"
    elif args.pocket_table:
        pocket_table_path = Path(args.pocket_table)
    else:
        parser.error("Provide either a pocket_table path or --config")

    pocket_rows = load_pocket_table(pocket_table_path)
    analysis = analyze_affinity_distribution(pocket_rows)

    out_path = (
        Path(args.output)
        if args.output
        else pocket_table_path.parent / "affinity_distribution_analysis.csv"
    )
    write_analysis_csv(out_path, analysis)

    print(format_summary(analysis))
    print(f"\nResults written to: {out_path}")


if __name__ == "__main__":
    main()
