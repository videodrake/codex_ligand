#!/usr/bin/env python3
"""Phase 2 Task Group 2.4: Druggability Confidence Layer.

Separates simple geometric pocket presence from stronger evidence of
ligandability/druggability by normalizing proposal scores, flagging
multi-source support, and providing an extensible schema for future
hotspot-based evidence (e.g., FTMap).

Tasks:
  2.4.1 - Proposal-score normalization (raw preserved, normalized added)
  2.4.2 - Multi-source support scoring (single-tool vs consensus)
  2.4.3 - Optional hotspot/druggability support (FTMap-ready schema)

Inputs:
  - output/phase2_pockets/candidate_pockets.csv
  - output/phase2_pockets/candidate_pocket_provenance.csv
  - output/phase2_pockets/pocket_patch_relationship.csv

Outputs:
  - output/phase2_pockets/druggability_proposal_summary.csv
  - output/phase2_pockets/candidate_pocket_support_flags.csv
  - output/phase2_pockets/phase2_druggability_confidence_note.md

Usage:
    python -m egfr_pipeline.phase2.druggability_confidence [--output_dir ...]
"""

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from egfr_pipeline.csv_utils import load_csv
from egfr_pipeline.parsing_utils import safe_float, safe_int

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

PHASE2_OUTPUT_DIR = PROJECT_ROOT / "output" / "phase2_pockets"

# Druggability confidence thresholds (based on fpocket druggability score)
# These are calibrated for fpocket's 0–1 druggability_score range.
# P2Rank's probability field is also 0–1 and mapped to druggability_score
# in pocket_proposal.py, so these thresholds apply to both tools.
# When additional tools are integrated, raw scores are preserved and
# the normalized_druggability_score provides cross-tool comparability.
HIGH_DRUGGABILITY_THRESHOLD = 0.50     # ≥ 0.50 → high
MEDIUM_DRUGGABILITY_THRESHOLD = 0.25   # ≥ 0.25 → medium, else low

# Multi-source consensus: a pocket supported by ≥2 independent tools
CONSENSUS_MIN_SOURCES = 2

# Known proposal tools (for support flag columns)
KNOWN_TOOLS = ["fpocket", "p2rank"]

# ---------------------------------------------------------------------------
# Output schemas
# ---------------------------------------------------------------------------

SUMMARY_COLUMNS = [
    "receptor_id",
    "pocket_id",
    "relationship_class",
    # Raw scores (Task 2.4.1: preserved)
    "best_proposal_score",
    "best_druggability_score",
    # Normalized scores (Task 2.4.1: added)
    "normalized_proposal_score",
    "normalized_druggability_score",
    "druggability_confidence",        # high / medium / low
    # Multi-source support (Task 2.4.2)
    "n_sources",
    "sources",
    "multi_source_support",           # True / False
    # Pocket geometry
    "mean_volume_A3",
    "n_residues",
    # Hotspot support (Task 2.4.3: extensible, empty until FTMap available)
    "ftmap_hotspot_count",
    "ftmap_overlap_residues",
    "hotspot_support_available",
    # Composite tier
    "overall_druggability_tier",      # tier_1 / tier_2 / tier_3
]

SUPPORT_FLAGS_COLUMNS = [
    "receptor_id",
    "pocket_id",
    # Per-tool support flags
    "fpocket_support",
    "fpocket_score",
    "fpocket_druggability",
    "p2rank_support",
    "p2rank_score",
    "p2rank_probability",
    # Aggregate
    "n_tools_supporting",
    "consensus_pocket",               # True if ≥ CONSENSUS_MIN_SOURCES
    # Future tools (Task 2.4.3)
    "ftmap_support",
    "ftmap_hotspot_count",
]


# ---------------------------------------------------------------------------
# Score normalization (Task 2.4.1)
# ---------------------------------------------------------------------------

def normalize_scores_within_state(
    pockets: List[dict],
    score_field: str,
) -> Dict[str, float]:
    """Min-max normalize a score field within each receptor state.

    Returns {pocket_id: normalized_score} mapping.
    """
    by_state: Dict[str, List[Tuple[str, float]]] = defaultdict(list)
    for p in pockets:
        val = safe_float(p.get(score_field, ""))
        if val is not None:
            by_state[p["receptor_id"]].append((p["pocket_id"], val))

    result = {}
    for state, entries in by_state.items():
        vals = [v for _, v in entries]
        vmin, vmax = min(vals), max(vals)
        span = vmax - vmin
        for pid, v in entries:
            if span > 0:
                result[pid] = (v - vmin) / span
            else:
                # All same value — assign 0.5 (neutral)
                result[pid] = 0.5

    return result


def assign_druggability_confidence(
    raw_drug_score: Optional[float],
    normalized_drug_score: Optional[float],
) -> str:
    """Assign druggability confidence label.

    Primary signal: raw druggability score (tool-native scale).
    Fallback: normalized score if raw is unavailable.
    """
    score = raw_drug_score if raw_drug_score is not None else normalized_drug_score
    if score is None:
        return "low"
    if score >= HIGH_DRUGGABILITY_THRESHOLD:
        return "high"
    if score >= MEDIUM_DRUGGABILITY_THRESHOLD:
        return "medium"
    return "low"


# ---------------------------------------------------------------------------
# Multi-source support (Task 2.4.2)
# ---------------------------------------------------------------------------

def build_per_tool_support(
    provenance: List[dict],
) -> Dict[str, dict]:
    """Build per-tool support flags for each merged pocket.

    Returns {pocket_id: {tool: {score, druggability, ...}}}.
    """
    support: Dict[str, Dict[str, dict]] = defaultdict(dict)

    for row in provenance:
        pid = row.get("merged_pocket_id", "")
        tool = row.get("proposal_source", "").lower().strip()
        if not pid or not tool:
            continue

        support[pid][tool] = {
            "score": row.get("proposal_score", ""),
            "druggability": row.get("druggability_score", ""),
        }

    return support


# ---------------------------------------------------------------------------
# Overall druggability tier (composite)
# ---------------------------------------------------------------------------

def assign_overall_tier(
    druggability_confidence: str,
    multi_source: bool,
    relationship_class: str,
) -> str:
    """Assign an overall druggability tier integrating multiple signals.

    tier_1: High-confidence druggable pocket near PPI patch
    tier_2: Medium confidence or single-evidence support
    tier_3: Low confidence or structurally distant
    """
    # PPI-relevant classes get a boost
    ppi_relevant = relationship_class in (
        "orthosteric_candidate", "rim_candidate", "allosteric_candidate"
    )

    if druggability_confidence == "high" and ppi_relevant:
        return "tier_1"
    if druggability_confidence == "high" and not ppi_relevant:
        return "tier_2"
    if druggability_confidence == "medium" and (ppi_relevant or multi_source):
        return "tier_2"
    if druggability_confidence == "medium":
        return "tier_3"
    # low confidence
    if multi_source and ppi_relevant:
        return "tier_2"
    return "tier_3"


# ---------------------------------------------------------------------------
# Report generation
# ---------------------------------------------------------------------------

def build_confidence_note(
    summary_rows: List[dict],
    support_rows: List[dict],
) -> List[str]:
    """Build the druggability confidence markdown note."""
    lines = [
        "# Phase 2 Druggability Confidence Note",
        "",
        "## Overview",
        "",
        "This note documents TG 2.4 outputs: score normalization, multi-source",
        "support assessment, and overall druggability tier assignment.",
        "",
        f"- Candidate pockets assessed: {len(summary_rows)}",
        "",
    ]

    # Thresholds
    lines.extend([
        "## Druggability Confidence Thresholds",
        "",
        "| Confidence | Criterion |",
        "|-----------|-----------|",
        f"| high | raw druggability score >= {HIGH_DRUGGABILITY_THRESHOLD} |",
        f"| medium | raw druggability score >= {MEDIUM_DRUGGABILITY_THRESHOLD} |",
        "| low | below medium threshold or no score |",
        "",
    ])

    # Tier definition
    lines.extend([
        "## Overall Druggability Tier Definition",
        "",
        "| Tier | Description |",
        "|------|------------|",
        "| tier_1 | High druggability + PPI-relevant (orthosteric/rim/allosteric) |",
        "| tier_2 | High but distant, or medium + PPI-relevant/multi-source, or low + multi-source + PPI-relevant |",
        "| tier_3 | Low confidence or distant with weak evidence |",
        "",
    ])

    # Summary table
    if summary_rows:
        lines.extend([
            "## Per-Pocket Summary",
            "",
            "| Pocket | Drugg. Conf. | Tier | Multi-Source | Relationship | Raw Score |",
            "|--------|-------------|------|-------------|-------------|-----------|",
        ])
        for r in summary_rows:
            lines.append(
                f"| {r['pocket_id']} | {r['druggability_confidence']} | "
                f"{r['overall_druggability_tier']} | {r['multi_source_support']} | "
                f"{r['relationship_class']} | {r['best_druggability_score']} |"
            )
        lines.append("")

    # Tier distribution
    tier_counts: Dict[str, int] = defaultdict(int)
    for r in summary_rows:
        tier_counts[r["overall_druggability_tier"]] += 1
    lines.extend([
        "## Tier Distribution",
        "",
        "| Tier | Count |",
        "|------|-------|",
    ])
    for tier in ["tier_1", "tier_2", "tier_3"]:
        lines.append(f"| {tier} | {tier_counts.get(tier, 0)} |")
    lines.append("")

    # FTMap note
    lines.extend([
        "## FTMap Hotspot Support",
        "",
        "FTMap columns are included in the output schema but currently empty.",
        "When FTMap data becomes available, re-run this module with `--ftmap_dir`",
        "to populate `ftmap_hotspot_count` and `ftmap_overlap_residues`.",
        "FTMap is the preferred hotspot-support method for PPI-site ligandability",
        "evidence (per PRD Phase 2, Task 2.4.3).",
        "",
        "---",
        "",
        "Generated by `egfr_pipeline.phase2.druggability_confidence`",
    ])

    return lines


# ---------------------------------------------------------------------------
# Main orchestrator
# ---------------------------------------------------------------------------

def run_druggability_confidence(
    output_dir: Path,
    ftmap_dir: Optional[Path] = None,
) -> Tuple[Path, Path, Path]:
    """Run full TG 2.4 pipeline.

    Returns (summary_csv, support_csv, note_path).
    """
    # Load inputs
    pockets = load_csv(output_dir / "candidate_pockets.csv")
    provenance = load_csv(output_dir / "candidate_pocket_provenance.csv")
    relationships = load_csv(output_dir / "pocket_patch_relationship.csv")

    if not pockets:
        raise FileNotFoundError("No candidate pockets found. Run TG 2.2 first.")

    print(f"  Loaded {len(pockets)} merged pockets, "
          f"{len(provenance)} provenance rows, "
          f"{len(relationships)} relationship rows")

    # Build lookup maps
    rel_map = {r["pocket_id"]: r for r in relationships}
    tool_support = build_per_tool_support(provenance)

    # Task 2.4.1: Score normalization
    norm_proposal = normalize_scores_within_state(pockets, "best_proposal_score")
    norm_druggability = normalize_scores_within_state(pockets, "best_druggability_score")

    print(f"  Normalized scores for {len(norm_proposal)} pockets")

    # Task 2.4.3: FTMap support (placeholder — populated when data available)
    ftmap_data: Dict[str, dict] = {}
    if ftmap_dir and ftmap_dir.exists():
        print(f"  FTMap directory provided: {ftmap_dir}")
        # Future: parse FTMap output here
        # ftmap_data = parse_ftmap_output(ftmap_dir)
    else:
        print("  FTMap data not available (placeholder columns will be empty)")

    # Build summary rows
    summary_rows = []
    support_rows = []

    for p in pockets:
        pid = p["pocket_id"]
        state = p["receptor_id"]

        # Raw scores
        raw_proposal = safe_float(p.get("best_proposal_score", ""))
        raw_druggability = safe_float(p.get("best_druggability_score", ""))

        # Normalized
        norm_p = norm_proposal.get(pid)
        norm_d = norm_druggability.get(pid)

        # Druggability confidence
        confidence = assign_druggability_confidence(raw_druggability, norm_d)

        # Multi-source
        n_sources = safe_int(p.get("n_sources", "1")) or 1
        sources = p.get("sources", "")
        multi_source = n_sources >= CONSENSUS_MIN_SOURCES

        # Relationship class
        rel = rel_map.get(pid, {})
        rel_class = rel.get("relationship_class", "")

        # Overall tier
        tier = assign_overall_tier(confidence, multi_source, rel_class)

        # FTMap (placeholder)
        ftmap_info = ftmap_data.get(pid, {})

        summary_rows.append({
            "receptor_id": state,
            "pocket_id": pid,
            "relationship_class": rel_class,
            "best_proposal_score": p.get("best_proposal_score", ""),
            "best_druggability_score": p.get("best_druggability_score", ""),
            "normalized_proposal_score": f"{norm_p:.4f}" if norm_p is not None else "",
            "normalized_druggability_score": f"{norm_d:.4f}" if norm_d is not None else "",
            "druggability_confidence": confidence,
            "n_sources": n_sources,
            "sources": sources,
            "multi_source_support": "true" if multi_source else "false",
            "mean_volume_A3": p.get("mean_volume_A3", ""),
            "n_residues": p.get("n_residues", ""),
            "ftmap_hotspot_count": ftmap_info.get("hotspot_count", ""),
            "ftmap_overlap_residues": ftmap_info.get("overlap_residues", ""),
            "hotspot_support_available": "true" if ftmap_info else "false",
            "overall_druggability_tier": tier,
        })

        # Per-tool support flags
        tools = tool_support.get(pid, {})
        support_rows.append({
            "receptor_id": state,
            "pocket_id": pid,
            "fpocket_support": "true" if "fpocket" in tools else "false",
            "fpocket_score": tools.get("fpocket", {}).get("score", ""),
            "fpocket_druggability": tools.get("fpocket", {}).get("druggability", ""),
            "p2rank_support": "true" if "p2rank" in tools else "false",
            "p2rank_score": tools.get("p2rank", {}).get("score", ""),
            "p2rank_probability": tools.get("p2rank", {}).get("druggability", ""),
            "n_tools_supporting": len(tools),
            "consensus_pocket": "true" if len(tools) >= CONSENSUS_MIN_SOURCES else "false",
            "ftmap_support": "true" if ftmap_info else "false",
            "ftmap_hotspot_count": ftmap_info.get("hotspot_count", ""),
        })

    # Print summary
    conf_counts: Dict[str, int] = defaultdict(int)
    tier_counts: Dict[str, int] = defaultdict(int)
    for r in summary_rows:
        conf_counts[r["druggability_confidence"]] += 1
        tier_counts[r["overall_druggability_tier"]] += 1

    print(f"  Druggability confidence: {dict(conf_counts)}")
    print(f"  Overall tiers: {dict(tier_counts)}")

    # Write outputs
    summary_path = output_dir / "druggability_proposal_summary.csv"
    _write_csv(summary_path, SUMMARY_COLUMNS, summary_rows)

    support_path = output_dir / "candidate_pocket_support_flags.csv"
    _write_csv(support_path, SUPPORT_FLAGS_COLUMNS, support_rows)

    # Write note
    note_lines = build_confidence_note(summary_rows, support_rows)
    note_path = output_dir / "phase2_druggability_confidence_note.md"
    with open(note_path, "w", encoding="utf-8") as f:
        f.write("\n".join(note_lines))

    return summary_path, support_path, note_path


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
        description="Phase 2 TG 2.4: Druggability confidence layer"
    )
    parser.add_argument(
        "--output_dir", type=Path, default=PHASE2_OUTPUT_DIR,
        help="Phase 2 output directory",
    )
    parser.add_argument(
        "--ftmap_dir", type=Path, default=None,
        help="Optional FTMap output directory for hotspot support",
    )
    args = parser.parse_args()

    print("Phase 2 — Task Group 2.4: Druggability Confidence Layer")
    print(f"  Output: {args.output_dir}")
    print()

    summary_path, support_path, note_path = run_druggability_confidence(
        args.output_dir, args.ftmap_dir
    )

    print(f"\n{'='*60}")
    print("TG 2.4 COMPLETE")
    print(f"{'='*60}")
    print(f"  Summary:       {summary_path}")
    print(f"  Support flags: {support_path}")
    print(f"  Note:          {note_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
