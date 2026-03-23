#!/usr/bin/env python3
"""Verdict weight sensitivity simulation (AC-4.1).

Simulates 6 weight combinations on existing verdict data (valid_sites.csv)
or synthetic pocket data. For each combination, re-normalizes the 3-axis
scores and recomputes STRONG/MODERATE/WEAK verdicts.

Outputs:
  - verdict_weight_sensitivity.csv  (pocket × combination matrix)
  - verdict_weight_sensitivity_report.md  (analysis + recommendation)

Weight combinations:
  baseline:   50/20/30  (current)
  ppi_boost:  40/30/30  (PPI 강화)
  ppi_max:    35/35/30  (Vina-PPI 동등)
  vina_weak:  30/40/30  (PPI 우선)
  equal:      33/33/34  (균등)
  cross_boost: 40/20/40  (Cross 강화)

Usage:
    python scripts/verdict_weight_sensitivity.py [--input valid_sites.csv]
    python scripts/verdict_weight_sensitivity.py --synthetic  # no CSV needed

EC-4.1: If all combinations yield identical STRONG set → "가중치에 robust" 판정.
"""

import argparse
import csv
import sys
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

PROJECT_ROOT = Path(__file__).resolve().parent.parent

# ---------------------------------------------------------------------------
# Weight combinations
# ---------------------------------------------------------------------------

WEIGHT_COMBINATIONS = {
    "baseline":    {"vina": 50, "ppi": 20, "cross": 30},
    "ppi_boost":   {"vina": 40, "ppi": 30, "cross": 30},
    "ppi_max":     {"vina": 35, "ppi": 35, "cross": 30},
    "vina_weak":   {"vina": 30, "ppi": 40, "cross": 30},
    "equal":       {"vina": 33, "ppi": 33, "cross": 34},
    "cross_boost": {"vina": 40, "ppi": 20, "cross": 40},
}

# Verdict thresholds (from verdict.py)
STRONG_MIN = 55
MODERATE_MIN = 30

# Raw sub-score maximums (from verdict.py score_pocket)
VINA_RAW_MAX = 60.0   # affinity(20) + convergence(15) + stability(10) + diversity(15)
PPI_RAW_MAX = 20.0     # spatial(10) + overlap(4+2) + repro(4) → capped at 20
CROSS_RAW_MAX = 30.0   # coverage(20) + support(10) → capped at 30

# Output columns
SENSITIVITY_COLUMNS = [
    "pocket_id",
    "receptor_id",
    "has_ppi_data",
    "vina_raw",
    "ppi_raw",
    "cross_raw",
]
# Dynamic columns added per combination: {combo}_total, {combo}_verdict


# ---------------------------------------------------------------------------
# Core simulation logic
# ---------------------------------------------------------------------------

def extract_raw_scores(row: dict) -> dict:
    """Extract raw axis scores from a valid_sites.csv row.

    Raw scores are the pre-normalization sub-totals from score_pocket().
    We reconstruct them from the stored component fields.
    """
    vina_raw = (
        _safe_float(row.get("vina_affinity_pts"))
        + _safe_float(row.get("vina_convergence_pts"))
        + _safe_float(row.get("vina_stability_pts"))
        + _safe_float(row.get("vina_diversity_pts"))
    )
    ppi_raw = (
        _safe_float(row.get("ppi_spatial_pts"))
        + _safe_float(row.get("ppi_overlap_pts"))
        + _safe_float(row.get("ppi_reproducibility_pts"))
    )
    cross_raw = _safe_float(row.get("cross_receptor_pts"))
    has_ppi = row.get("ppi_data_available", "").lower() in ("yes", "true", "1")

    return {
        "pocket_id": row.get("pocket_id", ""),
        "receptor_id": row.get("receptor_id", ""),
        "has_ppi_data": has_ppi,
        "vina_raw": round(vina_raw, 2),
        "ppi_raw": round(min(ppi_raw, PPI_RAW_MAX), 2),
        "cross_raw": round(min(cross_raw, CROSS_RAW_MAX), 2),
        "baseline_verdict": row.get("verdict", ""),
        "baseline_total": _safe_float(row.get("confidence_score")),
    }


def rescore(
    pocket: dict,
    weights: Dict[str, int],
) -> Tuple[float, str]:
    """Re-normalize raw scores with given weight allocation.

    Returns (total_score, verdict).
    """
    has_ppi = pocket["has_ppi_data"]

    if has_ppi:
        vina_max = float(weights["vina"])
        ppi_max = float(weights["ppi"])
        cross_max = float(weights["cross"])
    else:
        # Adaptive: redistribute PPI weight to Vina and Cross proportionally
        vina_base = weights["vina"]
        cross_base = weights["cross"]
        ppi_redistribute = weights["ppi"]
        ratio = vina_base / (vina_base + cross_base) if (vina_base + cross_base) > 0 else 0.5
        vina_max = float(vina_base + ppi_redistribute * ratio)
        ppi_max = 0.0
        cross_max = float(cross_base + ppi_redistribute * (1 - ratio))

    # Normalize each axis: raw / raw_max * allocated_max
    vina_score = min(pocket["vina_raw"], VINA_RAW_MAX) / VINA_RAW_MAX * vina_max
    ppi_score = min(pocket["ppi_raw"], PPI_RAW_MAX) / PPI_RAW_MAX * ppi_max if ppi_max > 0 else 0.0
    cross_score = min(pocket["cross_raw"], CROSS_RAW_MAX) / CROSS_RAW_MAX * cross_max

    total = round(vina_score + ppi_score + cross_score, 2)

    if total >= STRONG_MIN:
        verdict = "STRONG"
    elif total >= MODERATE_MIN:
        verdict = "MODERATE"
    else:
        verdict = "WEAK"

    return total, verdict


def run_sensitivity(
    pockets: List[dict],
    combinations: Dict[str, Dict[str, int]] = None,
) -> Tuple[List[dict], List[str]]:
    """Run weight sensitivity simulation.

    Args:
        pockets: List of extracted raw score dicts.
        combinations: Weight combos to test.

    Returns:
        (result_rows, warnings)
    """
    if combinations is None:
        combinations = WEIGHT_COMBINATIONS

    combo_names = list(combinations.keys())
    results = []

    for pocket in pockets:
        row = {
            "pocket_id": pocket["pocket_id"],
            "receptor_id": pocket["receptor_id"],
            "has_ppi_data": str(pocket["has_ppi_data"]),
            "vina_raw": pocket["vina_raw"],
            "ppi_raw": pocket["ppi_raw"],
            "cross_raw": pocket["cross_raw"],
        }

        for combo_name in combo_names:
            total, verdict = rescore(pocket, combinations[combo_name])
            row[f"{combo_name}_total"] = total
            row[f"{combo_name}_verdict"] = verdict

        results.append(row)

    # Analyze verdict changes
    warnings = []
    baseline_strong = {
        r["pocket_id"] for r in results
        if r.get("baseline_verdict") == "STRONG"
    }

    # Track changes across combinations
    all_strong_sets = {}
    for combo in combo_names:
        strong_set = {
            r["pocket_id"] for r in results
            if r.get(f"{combo}_verdict") == "STRONG"
        }
        all_strong_sets[combo] = strong_set

    # EC-4.1: identical STRONG set across all combinations
    strong_values = list(all_strong_sets.values())
    if strong_values and all(s == strong_values[0] for s in strong_values):
        warnings.append(
            "가중치에 robust: 모든 조합에서 STRONG 포켓 집합이 동일합니다. "
            "현재 배분 유지를 권장합니다."
        )

    # PPI promotion: pockets that become STRONG with higher PPI weight
    baseline_strongs = all_strong_sets.get("baseline", set())
    for combo in ["ppi_boost", "ppi_max", "vina_weak"]:
        promoted = all_strong_sets.get(combo, set()) - baseline_strongs
        if promoted:
            warnings.append(
                f"{combo}: PPI 가중치 증가로 {len(promoted)}개 포켓 STRONG 승격 — "
                f"{', '.join(sorted(promoted))}"
            )

    return results, warnings


# ---------------------------------------------------------------------------
# Report generation
# ---------------------------------------------------------------------------

def build_sensitivity_report(
    results: List[dict],
    warnings: List[str],
    combinations: Dict[str, Dict[str, int]] = None,
) -> List[str]:
    """Build verdict_weight_sensitivity_report.md."""
    if combinations is None:
        combinations = WEIGHT_COMBINATIONS

    combo_names = list(combinations.keys())

    lines = [
        "# Verdict Weight Sensitivity Report",
        "",
        "## 1. 가중치 조합",
        "",
        "| 조합 | Vina | PPI | Cross | 의도 |",
        "|------|------|-----|-------|------|",
    ]
    intents = {
        "baseline": "현재 설정",
        "ppi_boost": "PPI 강화",
        "ppi_max": "Vina-PPI 동등",
        "vina_weak": "PPI 우선",
        "equal": "균등 배분",
        "cross_boost": "Cross 강화",
    }
    for name, w in combinations.items():
        lines.append(
            f"| {name} | {w['vina']} | {w['ppi']} | {w['cross']} | "
            f"{intents.get(name, '')} |"
        )
    lines.append("")

    # Summary table
    lines.extend([
        "## 2. 포켓별 판정 비교",
        "",
        "| pocket_id | receptor_id | " + " | ".join(combo_names) + " |",
        "|-----------|-------------|" + "|".join(["-----"] * len(combo_names)) + "|",
    ])
    for r in results:
        verdicts = [f"{r[f'{c}_verdict']}({r[f'{c}_total']:.0f})" for c in combo_names]
        lines.append(
            f"| {r['pocket_id']} | {r['receptor_id']} | " + " | ".join(verdicts) + " |"
        )
    lines.append("")

    # Verdict changes
    lines.extend([
        "## 3. 판정 변화 분석",
        "",
    ])

    changed_pockets = []
    for r in results:
        verdicts = {r[f"{c}_verdict"] for c in combo_names}
        if len(verdicts) > 1:
            changed_pockets.append(r)

    if changed_pockets:
        lines.append(f"**{len(changed_pockets)}개 포켓**에서 조합 간 판정이 변합니다:")
        lines.append("")
        for r in changed_pockets:
            changes = []
            for c in combo_names:
                changes.append(f"{c}={r[f'{c}_verdict']}")
            lines.append(f"- **{r['pocket_id']}** ({r['receptor_id']}): {', '.join(changes)}")
        lines.append("")
    else:
        lines.append("판정 변화 없음 — 모든 포켓이 가중치 변화에 robust합니다.")
        lines.append("")

    # Stability classification
    lines.extend([
        "## 4. 안정성 분류",
        "",
        "| 분류 | 정의 | 포켓 수 |",
        "|------|------|---------|",
    ])
    stable = [r for r in results if len({r[f"{c}_verdict"] for c in combo_names}) == 1]
    sensitive = [r for r in results if len({r[f"{c}_verdict"] for c in combo_names}) > 1]
    lines.append(f"| 안정 | 모든 조합에서 동일 판정 | {len(stable)} |")
    lines.append(f"| 민감 | 조합 간 판정 변화 | {len(sensitive)} |")
    lines.append("")

    # Warnings
    if warnings:
        lines.extend(["## 5. 경고 및 권장", ""])
        for w in warnings:
            lines.append(f"- {w}")
        lines.append("")

    lines.extend([
        "---",
        "",
        "Generated by `scripts/verdict_weight_sensitivity.py`",
    ])

    return lines


# ---------------------------------------------------------------------------
# Synthetic data for offline simulation
# ---------------------------------------------------------------------------

def generate_synthetic_pockets() -> List[dict]:
    """Generate synthetic pocket data for testing without valid_sites.csv.

    Creates a representative set of pockets spanning the scoring space.
    """
    return [
        # High Vina, high PPI, high cross → STRONG everywhere
        {"pocket_id": "P_strong_all", "receptor_id": "R1", "has_ppi_data": True,
         "vina_raw": 50.0, "ppi_raw": 18.0, "cross_raw": 28.0},
        # High Vina only → borderline
        {"pocket_id": "P_vina_only", "receptor_id": "R1", "has_ppi_data": True,
         "vina_raw": 55.0, "ppi_raw": 0.0, "cross_raw": 0.0},
        # PPI strong, Vina moderate → sensitive to PPI weight
        {"pocket_id": "P_ppi_sensitive", "receptor_id": "R1", "has_ppi_data": True,
         "vina_raw": 30.0, "ppi_raw": 18.0, "cross_raw": 15.0},
        # Cross strong, others weak
        {"pocket_id": "P_cross_only", "receptor_id": "R1", "has_ppi_data": True,
         "vina_raw": 20.0, "ppi_raw": 5.0, "cross_raw": 28.0},
        # Balanced moderate
        {"pocket_id": "P_balanced_mod", "receptor_id": "R1", "has_ppi_data": True,
         "vina_raw": 35.0, "ppi_raw": 10.0, "cross_raw": 15.0},
        # No PPI data
        {"pocket_id": "P_no_ppi", "receptor_id": "R2", "has_ppi_data": False,
         "vina_raw": 45.0, "ppi_raw": 0.0, "cross_raw": 25.0},
        # Weak all around
        {"pocket_id": "P_weak", "receptor_id": "R1", "has_ppi_data": True,
         "vina_raw": 10.0, "ppi_raw": 2.0, "cross_raw": 5.0},
        # STRONG boundary (55 at baseline)
        {"pocket_id": "P_boundary", "receptor_id": "R1", "has_ppi_data": True,
         "vina_raw": 38.0, "ppi_raw": 12.0, "cross_raw": 18.0},
    ]


# ---------------------------------------------------------------------------
# AC-4.3: Cross-receptor coverage differentiation simulation
# ---------------------------------------------------------------------------

CROSS_COVERAGE_SCHEMES = {
    "current":  {"3of3": 20.0, "2of3": 20.0, "1of3": 10.0},   # before: 2/3 and 3/3 same
    "scheme_A": {"3of3": 20.0, "2of3": 14.0, "1of3": 0.0},    # default in verdict.py now
    "scheme_B": {"3of3": 22.0, "2of3": 12.0, "1of3": 0.0},    # stronger 3/3 premium
}


def simulate_cross_differentiation(
    pockets: List[dict],
) -> Tuple[List[dict], List[str]]:
    """Simulate cross-receptor coverage differentiation (AC-4.3).

    Each pocket needs: n_cross, cross_coverage_pts, cross_support_pts.
    Re-computes total with different coverage schemes.

    Returns (result_rows, warnings).
    """
    scheme_names = list(CROSS_COVERAGE_SCHEMES.keys())
    results = []

    for pocket in pockets:
        n_cross = pocket.get("n_cross", 0)
        support_pts = pocket.get("cross_support_pts", 0.0)
        has_ppi = pocket.get("has_ppi_data", True)
        cross_max = 30.0 if has_ppi else 40.0

        # Base scores (non-cross)
        vina_max = 50.0 if has_ppi else 60.0
        ppi_max = 20.0 if has_ppi else 0.0
        vina_score = min(pocket["vina_raw"], VINA_RAW_MAX) / VINA_RAW_MAX * vina_max
        ppi_score = min(pocket["ppi_raw"], PPI_RAW_MAX) / PPI_RAW_MAX * ppi_max if ppi_max > 0 else 0.0

        row = {
            "pocket_id": pocket["pocket_id"],
            "receptor_id": pocket["receptor_id"],
            "n_cross": n_cross,
            "n_states": n_cross + 1,
        }

        for scheme_name in scheme_names:
            scheme = CROSS_COVERAGE_SCHEMES[scheme_name]
            if n_cross >= 2:
                coverage = scheme["3of3"]
            elif n_cross >= 1:
                coverage = scheme["2of3"]
            else:
                coverage = scheme["1of3"]

            cross_score = min(coverage + support_pts, cross_max)
            total = round(vina_score + ppi_score + cross_score, 2)
            verdict = "STRONG" if total >= STRONG_MIN else ("MODERATE" if total >= MODERATE_MIN else "WEAK")

            row[f"{scheme_name}_total"] = total
            row[f"{scheme_name}_verdict"] = verdict

        results.append(row)

    # Find pockets where verdict differs between schemes
    warnings = []
    changed = []
    for r in results:
        verdicts = {r[f"{s}_verdict"] for s in scheme_names}
        if len(verdicts) > 1:
            changed.append(r)

    if changed:
        for r in changed:
            diffs = ", ".join(f"{s}={r[f'{s}_verdict']}({r[f'{s}_total']:.0f})" for s in scheme_names)
            warnings.append(
                f"Cross 차등 영향: {r['pocket_id']} ({r['n_states']}/3 states) — {diffs}"
            )

    return results, warnings


def generate_synthetic_cross_pockets() -> List[dict]:
    """Generate synthetic pockets with n_cross for cross differentiation test."""
    return [
        # 3/3 states, high support → should be STRONG in all schemes
        {"pocket_id": "P_3of3_strong", "receptor_id": "R1", "has_ppi_data": True,
         "vina_raw": 40.0, "ppi_raw": 15.0, "cross_raw": 28.0,
         "n_cross": 2, "cross_coverage_pts": 20.0, "cross_support_pts": 8.0},
        # 2/3 states, good support → sensitive to differentiation
        {"pocket_id": "P_2of3_border", "receptor_id": "R1", "has_ppi_data": True,
         "vina_raw": 35.0, "ppi_raw": 10.0, "cross_raw": 26.0,
         "n_cross": 1, "cross_coverage_pts": 20.0, "cross_support_pts": 6.0},
        # 2/3 states, weak support → may drop with differentiation
        {"pocket_id": "P_2of3_weak", "receptor_id": "R1", "has_ppi_data": True,
         "vina_raw": 30.0, "ppi_raw": 8.0, "cross_raw": 13.0,
         "n_cross": 1, "cross_coverage_pts": 10.0, "cross_support_pts": 3.0},
        # 1/3 state only
        {"pocket_id": "P_1of3_only", "receptor_id": "R1", "has_ppi_data": True,
         "vina_raw": 40.0, "ppi_raw": 12.0, "cross_raw": 10.0,
         "n_cross": 0, "cross_coverage_pts": 10.0, "cross_support_pts": 0.0},
        # 3/3 vs 2/3 tied at current scheme → should differentiate
        {"pocket_id": "P_tied_3of3", "receptor_id": "R1", "has_ppi_data": True,
         "vina_raw": 32.0, "ppi_raw": 10.0, "cross_raw": 26.0,
         "n_cross": 2, "cross_coverage_pts": 20.0, "cross_support_pts": 6.0},
        {"pocket_id": "P_tied_2of3", "receptor_id": "R1", "has_ppi_data": True,
         "vina_raw": 32.0, "ppi_raw": 10.0, "cross_raw": 26.0,
         "n_cross": 1, "cross_coverage_pts": 20.0, "cross_support_pts": 6.0},
    ]


# ---------------------------------------------------------------------------
# AC-4.4: PPI spatial distance distribution analysis
# ---------------------------------------------------------------------------

PPI_DIST_THRESHOLDS = [8.0, 15.0, 25.0]
PPI_DIST_LABELS = {
    "adjacent": (0.0, 8.0),
    "near": (8.0, 15.0),
    "moderate": (15.0, 25.0),
    "distant": (25.0, float("inf")),
}
CENTROID_OFFSET_RANGE = (3.0, 5.0)  # systematic Vina centroid offset (Å)


def analyze_ppi_distance_distribution(
    distances: List[float],
    thresholds: List[float] = None,
    offset: float = 0.0,
) -> Tuple[List[dict], List[str]]:
    """Analyze PPI distance distribution against current thresholds.

    Args:
        distances: List of spatial_dist_to_ppi values (Å).
        thresholds: Bin edges [8, 15, 25] (default).
        offset: Centroid offset correction to apply (Å).

    Returns:
        (distribution_rows, warnings)
    """
    if thresholds is None:
        thresholds = list(PPI_DIST_THRESHOLDS)

    if not distances:
        return [], ["WARNING: No PPI distance data available"]

    # Apply offset correction
    corrected = [max(d - offset, 0.0) for d in distances]
    n_total = len(corrected)

    # Bin into categories
    bins = thresholds + [float("inf")]
    bin_labels = ["adjacent", "near", "moderate", "distant"]
    bin_counts = [0] * len(bin_labels)

    for d in corrected:
        for i, upper in enumerate(bins):
            lower = bins[i - 1] if i > 0 else 0.0
            if lower <= d < upper:
                bin_counts[i] += 1
                break

    distribution_rows = []
    for i, label in enumerate(bin_labels):
        lower = bins[i - 1] if i > 0 else 0.0
        upper = bins[i] if i < len(bins) else float("inf")
        count = bin_counts[i]
        frac = count / n_total if n_total > 0 else 0.0

        distribution_rows.append({
            "bin": label,
            "range_lower_A": round(lower, 1),
            "range_upper_A": round(upper, 1) if upper != float("inf") else "inf",
            "n_pockets": count,
            "fraction": round(frac, 4),
            "ppi_spatial_pts": _bin_to_pts(label),
            "offset_applied_A": round(offset, 1),
        })

    # Warnings
    warnings = []

    # AC-4.4: 70% concentration check
    for row in distribution_rows:
        if row["fraction"] > 0.70:
            warnings.append(
                f"포켓 {row['fraction']:.0%}가 '{row['bin']}' 구간에 집중 (> 70%). "
                f"구간 세분화를 권장합니다."
            )

    # Percentile summary
    sorted_d = sorted(corrected)
    percentiles = {}
    for p in [10, 25, 50, 75, 90]:
        idx = int(len(sorted_d) * p / 100)
        idx = min(idx, len(sorted_d) - 1)
        percentiles[f"p{p}"] = round(sorted_d[idx], 2)

    warnings.append(
        f"거리 분포 percentiles (offset={offset:.1f}Å): "
        + ", ".join(f"{k}={v}Å" for k, v in percentiles.items())
    )

    # Effective range with offset
    if offset > 0:
        effective = [f"{t - offset:.1f}–{t:.1f}Å" for t in thresholds]
        warnings.append(
            f"오프셋 {offset:.1f}Å 감안 실질 범위: " + ", ".join(
                f"{label}: {eff}" for label, eff in zip(bin_labels[:3], effective)
            )
        )

    return distribution_rows, warnings


def _bin_to_pts(label: str) -> float:
    """Map bin label to current PPI spatial points."""
    return {"adjacent": 10.0, "near": 6.0, "moderate": 2.0, "distant": 0.0}.get(label, 0.0)


def generate_synthetic_distances(n: int = 30) -> List[float]:
    """Generate synthetic PPI distances for testing."""
    import random
    random.seed(42)
    # Realistic distribution: most pockets are moderate-to-distant
    distances = []
    for _ in range(n):
        r = random.random()
        if r < 0.10:       # 10% adjacent
            distances.append(random.uniform(3.0, 8.0))
        elif r < 0.25:     # 15% near
            distances.append(random.uniform(8.0, 15.0))
        elif r < 0.55:     # 30% moderate
            distances.append(random.uniform(15.0, 25.0))
        else:               # 45% distant
            distances.append(random.uniform(25.0, 50.0))
    return [round(d, 2) for d in distances]


def build_distance_report(
    distribution_rows: List[dict],
    warnings: List[str],
) -> List[str]:
    """Build PPI distance distribution section for the report."""
    lines = [
        "## PPI 거리 분포 분석 (AC-4.4)",
        "",
        "### 현재 임계값 분할",
        "",
        "| 구간 | 범위 (Å) | 포켓 수 | 비율 | spatial_pts |",
        "|------|----------|---------|------|-------------|",
    ]
    for r in distribution_rows:
        upper = r["range_upper_A"] if r["range_upper_A"] != "inf" else "∞"
        lines.append(
            f"| {r['bin']} | {r['range_lower_A']}–{upper} | "
            f"{r['n_pockets']} | {r['fraction']:.0%} | {r['ppi_spatial_pts']} |"
        )
    lines.append("")

    for w in warnings:
        lines.append(f"- {w}")
    lines.append("")

    return lines


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _safe_float(val, default=0.0):
    try:
        return float(val) if val not in (None, "", "NA", "None") else default
    except (ValueError, TypeError):
        return default


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Verdict weight sensitivity simulation (AC-4.1)"
    )
    parser.add_argument(
        "--input", type=Path, default=None,
        help="Path to valid_sites.csv (if available)"
    )
    parser.add_argument(
        "--synthetic", action="store_true",
        help="Use synthetic pocket data (no valid_sites.csv needed)"
    )
    parser.add_argument(
        "--output_dir", type=Path,
        default=PROJECT_ROOT / "output",
        help="Output directory for results"
    )
    args = parser.parse_args()

    print("Verdict Weight Sensitivity Simulation (AC-4.1)")
    print()

    # Load data
    if args.synthetic or args.input is None:
        if args.input and args.input.exists():
            # Load real data
            with open(args.input, encoding="utf-8") as f:
                raw_rows = list(csv.DictReader(f))
            pockets = [extract_raw_scores(r) for r in raw_rows]
            print(f"  Loaded {len(pockets)} pockets from {args.input}")
        else:
            pockets = generate_synthetic_pockets()
            print(f"  Using {len(pockets)} synthetic pockets")
    else:
        if not args.input.exists():
            print(f"  ERROR: {args.input} not found. Use --synthetic.")
            return 1
        with open(args.input, encoding="utf-8") as f:
            raw_rows = list(csv.DictReader(f))
        pockets = [extract_raw_scores(r) for r in raw_rows]
        print(f"  Loaded {len(pockets)} pockets from {args.input}")

    # Run simulation
    results, warnings = run_sensitivity(pockets)

    for w in warnings:
        print(f"  WARNING: {w}")

    # Write CSV
    args.output_dir.mkdir(parents=True, exist_ok=True)
    combo_names = list(WEIGHT_COMBINATIONS.keys())
    all_columns = list(SENSITIVITY_COLUMNS)
    for c in combo_names:
        all_columns.extend([f"{c}_total", f"{c}_verdict"])

    csv_path = args.output_dir / "verdict_weight_sensitivity.csv"
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=all_columns, extrasaction="ignore")
        w.writeheader()
        w.writerows(results)
    print(f"  CSV: {csv_path}")

    # Write report
    report_lines = build_sensitivity_report(results, warnings)
    report_path = args.output_dir / "verdict_weight_sensitivity_report.md"
    report_path.write_text("\n".join(report_lines), encoding="utf-8")
    print(f"  Report: {report_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
