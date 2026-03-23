#!/usr/bin/env python3
"""Workflow A↔B comparison module (AC-5.1).

Matches pockets between Workflow A (Verdict 3-axis) and Workflow B
(Phase 4 Perturbation 4-axis) using centroid distance + residue Jaccard,
then classifies into Consensus / A-only / B-only / Conflict.

Inputs:
  Workflow A: valid_sites.csv + vina_pocket_table.csv (for centroids/residues)
  Workflow B: phase4_final_review_table.csv + candidate_pockets.csv (for centroids/residues)

Outputs:
  workflow_comparison.csv
  workflow_comparison_report.md

Usage:
    python -m egfr_pipeline.workflow_comparison --config config.yaml
    python -m egfr_pipeline.workflow_comparison --synthetic
"""

import argparse
import csv
import math
import sys
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

PROJECT_ROOT = Path(__file__).resolve().parent.parent

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

CENTROID_MATCH_CUTOFF = 8.0   # Å
JACCARD_MATCH_CUTOFF = 0.3

# Workflow B "top" threshold
B_TOP_SCORE_CUTOFF = 0.5     # perturbation_score > 0.5 → "top"

COMPARISON_COLUMNS = [
    "pocket_id_a",
    "receptor_id_a",
    "pocket_id_b",
    "receptor_id_b",
    "comparison_category",    # Consensus / A-only / B-only / Conflict
    "verdict_a",              # STRONG / MODERATE / WEAK
    "score_a",
    "class_b",                # mechanistic_class from Phase 4
    "score_b",                # perturbation_score
    "rank_b",
    "centroid_dist_A",
    "residue_jaccard",
    "allosteric_flag",        # True if A-only + allosteric_candidate
    "bias_flag",              # True if B-only (blind docking bias?)
    "is_atp_site",
]


# ---------------------------------------------------------------------------
# Pocket data loading
# ---------------------------------------------------------------------------

def load_workflow_a(
    verdict_path: Path,
    pocket_table_path: Path,
) -> List[dict]:
    """Load Workflow A pockets with centroid and residues."""
    verdicts = _load_csv(verdict_path)
    pocket_table = _load_csv(pocket_table_path)

    # Index pocket table by (receptor_id, pocket_id)
    pt_index = {}
    for row in pocket_table:
        key = (row.get("receptor_id", ""), row.get("pocket_id", ""))
        pt_index[key] = row

    pockets = []
    for v in verdicts:
        # Skip ATP site pockets
        if str(v.get("is_atp_site", "")).lower() in ("true", "1"):
            continue

        key = (v.get("receptor_id", ""), v.get("pocket_id", ""))
        pt = pt_index.get(key, {})

        pockets.append({
            "pocket_id": v.get("pocket_id", ""),
            "receptor_id": v.get("receptor_id", ""),
            "verdict": v.get("verdict", ""),
            "score": _safe_float(v.get("confidence_score")),
            "vina_quality_score": _safe_float(v.get("vina_quality_score")),
            "ppi_proximity_score": _safe_float(v.get("ppi_proximity_score")),
            "allosteric_candidate": str(v.get("allosteric_candidate", "")).lower() in ("true", "1"),
            "is_atp_site": False,
            "centroid": (
                _safe_float(pt.get("centroid_x")),
                _safe_float(pt.get("centroid_y")),
                _safe_float(pt.get("centroid_z")),
            ),
            "residues": _parse_residues(pt.get("union_contact_residues", "")),
        })

    return pockets


def load_workflow_b(
    review_table_path: Path,
    candidate_pockets_path: Path,
) -> List[dict]:
    """Load Workflow B pockets with centroid and residues."""
    reviews = _load_csv(review_table_path)
    candidates = _load_csv(candidate_pockets_path)

    # Index candidates by (receptor_id, pocket_id)
    cp_index = {}
    for row in candidates:
        key = (row.get("receptor_id", ""), row.get("pocket_id", ""))
        cp_index[key] = row

    pockets = []
    for r in reviews:
        key = (r.get("receptor_id", ""), r.get("pocket_id", ""))
        cp = cp_index.get(key, {})

        pockets.append({
            "pocket_id": r.get("pocket_id", ""),
            "receptor_id": r.get("receptor_id", ""),
            "mechanistic_class": r.get("mechanistic_class", ""),
            "perturbation_score": _safe_float(r.get("perturbation_score")),
            "rank": _safe_int(r.get("rank")),
            "centroid": (
                _safe_float(cp.get("centroid_x")),
                _safe_float(cp.get("centroid_y")),
                _safe_float(cp.get("centroid_z")),
            ),
            "residues": _parse_residues(cp.get("residue_ids", "")),
        })

    return pockets


# ---------------------------------------------------------------------------
# Matching
# ---------------------------------------------------------------------------

def compute_centroid_distance(a: Tuple[float, ...], b: Tuple[float, ...]) -> float:
    """Euclidean distance between two 3D points."""
    return math.sqrt(sum((ai - bi) ** 2 for ai, bi in zip(a, b)))


def compute_jaccard(set_a: Set[str], set_b: Set[str]) -> float:
    """Jaccard similarity between two residue sets."""
    if not set_a and not set_b:
        return 1.0
    union = len(set_a | set_b)
    return len(set_a & set_b) / union if union > 0 else 0.0


def match_pockets(
    a_pockets: List[dict],
    b_pockets: List[dict],
    centroid_cutoff: float = CENTROID_MATCH_CUTOFF,
    jaccard_cutoff: float = JACCARD_MATCH_CUTOFF,
) -> Tuple[List[Tuple[dict, dict, float, float]], List[dict], List[dict]]:
    """Match A↔B pockets by centroid distance + Jaccard.

    Returns:
        (matched_pairs, unmatched_a, unmatched_b)
        where each pair is (a_pocket, b_pocket, centroid_dist, jaccard)
    """
    # Compute all pairwise scores
    candidates = []
    for a in a_pockets:
        for b in b_pockets:
            dist = compute_centroid_distance(a["centroid"], b["centroid"])
            if dist > centroid_cutoff:
                continue
            jacc = compute_jaccard(a["residues"], b["residues"])
            if jacc < jaccard_cutoff:
                continue
            candidates.append((a, b, dist, jacc))

    # Greedy matching: closest first
    candidates.sort(key=lambda x: x[2])
    matched_a = set()
    matched_b = set()
    matched_pairs = []

    for a, b, dist, jacc in candidates:
        a_key = (a["receptor_id"], a["pocket_id"])
        b_key = (b["receptor_id"], b["pocket_id"])
        if a_key in matched_a or b_key in matched_b:
            continue
        matched_pairs.append((a, b, dist, jacc))
        matched_a.add(a_key)
        matched_b.add(b_key)

    unmatched_a = [a for a in a_pockets
                   if (a["receptor_id"], a["pocket_id"]) not in matched_a]
    unmatched_b = [b for b in b_pockets
                   if (b["receptor_id"], b["pocket_id"]) not in matched_b]

    return matched_pairs, unmatched_a, unmatched_b


# ---------------------------------------------------------------------------
# Classification
# ---------------------------------------------------------------------------

def _is_b_top(b_pocket: dict, n_b_total: int) -> bool:
    """Check if B pocket is in the top tier."""
    if b_pocket["perturbation_score"] > B_TOP_SCORE_CUTOFF:
        return True
    if n_b_total > 0 and b_pocket["rank"] <= math.ceil(n_b_total / 2):
        return True
    return False


def _is_b_irrelevant(b_pocket: dict) -> bool:
    mc = b_pocket.get("mechanistic_class", "")
    return mc in (
        "ligandable_but_ppi_irrelevant_candidate",
        "uncertain_mechanism_candidate",
        "low_relevance_candidate",
        "",
    )


def classify_comparison(
    matched_pairs: List[Tuple[dict, dict, float, float]],
    unmatched_a: List[dict],
    unmatched_b: List[dict],
    n_b_total: int,
) -> List[dict]:
    """Classify each pocket into Consensus/A-only/B-only/Conflict."""
    results = []

    for a, b, dist, jacc in matched_pairs:
        a_strong = a["verdict"] in ("STRONG", "MODERATE")
        b_top = _is_b_top(b, n_b_total)
        b_irr = _is_b_irrelevant(b)

        if a_strong and b_top:
            category = "Consensus"
        elif a_strong and b_irr:
            category = "Conflict"
        elif not a_strong and b_top:
            category = "Conflict"
        elif a_strong and not b_top:
            category = "A-only"
        else:
            category = "Consensus" if not a_strong and not b_top else "B-only"

        results.append({
            "pocket_id_a": a["pocket_id"],
            "receptor_id_a": a["receptor_id"],
            "pocket_id_b": b["pocket_id"],
            "receptor_id_b": b["receptor_id"],
            "comparison_category": category,
            "verdict_a": a["verdict"],
            "score_a": a["score"],
            "class_b": b.get("mechanistic_class", ""),
            "score_b": b["perturbation_score"],
            "rank_b": b["rank"],
            "centroid_dist_A": round(dist, 2),
            "residue_jaccard": round(jacc, 4),
            "allosteric_flag": a.get("allosteric_candidate", False) and category == "A-only",
            "bias_flag": False,
            "is_atp_site": False,
        })

    # Unmatched A → A-only
    for a in unmatched_a:
        results.append({
            "pocket_id_a": a["pocket_id"],
            "receptor_id_a": a["receptor_id"],
            "pocket_id_b": "",
            "receptor_id_b": "",
            "comparison_category": "A-only",
            "verdict_a": a["verdict"],
            "score_a": a["score"],
            "class_b": "",
            "score_b": "",
            "rank_b": "",
            "centroid_dist_A": "",
            "residue_jaccard": "",
            "allosteric_flag": a.get("allosteric_candidate", False),
            "bias_flag": False,
            "is_atp_site": False,
        })

    # Unmatched B → B-only
    for b in unmatched_b:
        results.append({
            "pocket_id_a": "",
            "receptor_id_a": "",
            "pocket_id_b": b["pocket_id"],
            "receptor_id_b": b["receptor_id"],
            "comparison_category": "B-only",
            "verdict_a": "",
            "score_a": "",
            "class_b": b.get("mechanistic_class", ""),
            "score_b": b["perturbation_score"],
            "rank_b": b["rank"],
            "centroid_dist_A": "",
            "residue_jaccard": "",
            "allosteric_flag": False,
            "bias_flag": True,
            "is_atp_site": False,
        })

    return results


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

def build_comparison_report(results: List[dict]) -> List[str]:
    """Build workflow_comparison_report.md."""
    lines = [
        "# Workflow A↔B Comparison Report",
        "",
        "## Summary",
        "",
    ]

    cats = {}
    for r in results:
        cat = r["comparison_category"]
        cats.setdefault(cat, []).append(r)

    lines.append("| Category | Count | Description |")
    lines.append("|----------|-------|-------------|")
    descs = {
        "Consensus": "Both workflows agree — strongest computational evidence",
        "A-only": "Verdict STRONG/MODERATE but B unmatched/irrelevant — allosteric?",
        "B-only": "B top-ranked but A weak/unmatched — blind docking bias?",
        "Conflict": "Workflows disagree — manual review needed",
    }
    for cat in ["Consensus", "A-only", "B-only", "Conflict"]:
        n = len(cats.get(cat, []))
        lines.append(f"| {cat} | {n} | {descs.get(cat, '')} |")
    lines.append("")

    # Consensus detail
    consensus = cats.get("Consensus", [])
    if consensus:
        lines.extend(["## Consensus Pockets", ""])
        for r in consensus:
            lines.append(
                f"- **{r['pocket_id_a']}** ({r['receptor_id_a']}): "
                f"A={r['verdict_a']}({r['score_a']}), "
                f"B={r['class_b']}({r['score_b']}), "
                f"dist={r['centroid_dist_A']}Å, J={r['residue_jaccard']}"
            )
        lines.append("")

    # Allosteric flags
    allosteric = [r for r in results if r.get("allosteric_flag")]
    if allosteric:
        lines.extend(["## Allosteric Candidates (A-only)", ""])
        for r in allosteric:
            lines.append(f"- **{r['pocket_id_a']}** ({r['receptor_id_a']}): A={r['verdict_a']}({r['score_a']})")
        lines.append("")

    # Conflict detail
    conflicts = cats.get("Conflict", [])
    if conflicts:
        lines.extend(["## Conflict Pockets (Manual Review Needed)", ""])
        for r in conflicts:
            lines.append(
                f"- **{r['pocket_id_a'] or r['pocket_id_b']}**: "
                f"A={r['verdict_a'] or 'N/A'}, B={r['class_b'] or 'N/A'}"
            )
        lines.append("")

    lines.extend(["---", "", "Generated by `egfr_pipeline.workflow_comparison`"])
    return lines


# ---------------------------------------------------------------------------
# Synthetic data
# ---------------------------------------------------------------------------

def generate_synthetic_data() -> Tuple[List[dict], List[dict]]:
    """Generate synthetic A and B pockets for testing."""
    a_pockets = [
        {"pocket_id": "PA1", "receptor_id": "R1", "verdict": "STRONG", "score": 72.0,
         "vina_quality_score": 40.0, "ppi_proximity_score": 12.0,
         "allosteric_candidate": False, "is_atp_site": False,
         "centroid": (10.0, 20.0, 30.0), "residues": {"ALA700", "ALA701", "ALA702"}},
        {"pocket_id": "PA2", "receptor_id": "R1", "verdict": "MODERATE", "score": 45.0,
         "vina_quality_score": 38.0, "ppi_proximity_score": 2.0,
         "allosteric_candidate": True, "is_atp_site": False,
         "centroid": (50.0, 60.0, 70.0), "residues": {"LEU800", "LEU801"}},
        {"pocket_id": "PA3", "receptor_id": "R1", "verdict": "WEAK", "score": 20.0,
         "vina_quality_score": 15.0, "ppi_proximity_score": 0.0,
         "allosteric_candidate": False, "is_atp_site": False,
         "centroid": (80.0, 90.0, 10.0), "residues": {"GLY900"}},
    ]
    b_pockets = [
        {"pocket_id": "PB1", "receptor_id": "R1",
         "mechanistic_class": "orthosteric_disruptor_candidate",
         "perturbation_score": 0.85, "rank": 1,
         "centroid": (10.5, 20.2, 30.1), "residues": {"ALA700", "ALA701", "ALA703"}},
        {"pocket_id": "PB2", "receptor_id": "R1",
         "mechanistic_class": "ligandable_but_ppi_irrelevant_candidate",
         "perturbation_score": 0.20, "rank": 3,
         "centroid": (50.3, 60.1, 70.2), "residues": {"LEU800", "LEU802"}},
        {"pocket_id": "PB3", "receptor_id": "R1",
         "mechanistic_class": "interface_rim_modulator_candidate",
         "perturbation_score": 0.65, "rank": 2,
         "centroid": (30.0, 40.0, 50.0), "residues": {"VAL850", "VAL851"}},
    ]
    return a_pockets, b_pockets


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _safe_float(val, default=0.0):
    try:
        return float(val) if val not in (None, "", "NA") else default
    except (ValueError, TypeError):
        return default


def _safe_int(val, default=0):
    try:
        return int(val) if val not in (None, "", "NA") else default
    except (ValueError, TypeError):
        return default


def _parse_residues(val: str) -> Set[str]:
    if not val:
        return set()
    return {r.strip() for r in val.split(";") if r.strip()}


def _load_csv(path: Path) -> List[dict]:
    if not path.exists():
        return []
    with open(path, encoding="utf-8") as f:
        return list(csv.DictReader(f))


def _write_csv(path: Path, rows: List[dict], columns: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=columns, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def run_comparison(
    a_pockets: List[dict],
    b_pockets: List[dict],
) -> Tuple[List[dict], List[str]]:
    """Run full comparison pipeline.

    Returns (comparison_rows, report_lines).
    """
    matched, unmatched_a, unmatched_b = match_pockets(a_pockets, b_pockets)
    results = classify_comparison(matched, unmatched_a, unmatched_b, len(b_pockets))
    report = build_comparison_report(results)
    return results, report


def main():
    parser = argparse.ArgumentParser(
        description="Workflow A↔B comparison (AC-5.1)"
    )
    parser.add_argument("--synthetic", action="store_true",
                        help="Use synthetic data for testing")
    parser.add_argument("--output_dir", type=Path,
                        default=PROJECT_ROOT / "output",
                        help="Output directory")
    args = parser.parse_args()

    print("Workflow A↔B Comparison (AC-5.1)")

    if args.synthetic:
        a_pockets, b_pockets = generate_synthetic_data()
        print(f"  Synthetic: {len(a_pockets)} A pockets, {len(b_pockets)} B pockets")
    else:
        print("  ERROR: Provide --synthetic or implement file loading path.")
        return 1

    results, report = run_comparison(a_pockets, b_pockets)

    # Write outputs
    args.output_dir.mkdir(parents=True, exist_ok=True)
    _write_csv(args.output_dir / "workflow_comparison.csv", results, COMPARISON_COLUMNS)
    (args.output_dir / "workflow_comparison_report.md").write_text(
        "\n".join(report), encoding="utf-8"
    )

    # Print summary
    cats = {}
    for r in results:
        cats.setdefault(r["comparison_category"], []).append(r)
    for cat in ["Consensus", "A-only", "B-only", "Conflict"]:
        print(f"  {cat}: {len(cats.get(cat, []))}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
