#!/usr/bin/env python3
"""Phase 4 Task Group 4.4: State-Robustness and Accessibility Interpretation.

Interprets cross-state evidence to distinguish robust, shifted, and
state-specific pockets. State-specificity is preserved as potentially
meaningful rather than auto-penalized.

Tasks:
  4.4.1 - Cross-state support interpretation
  4.4.2 - Accessibility-aware scoring influence
  4.4.3 - Cryptic/allosteric caution handling

Inputs:
  - output/phase4_perturbation/perturbation_candidate_table.csv
  - output/phase4_perturbation/phase4_evidence_normalized.csv

Outputs:
  - output/phase4_perturbation/phase4_state_interpretation.csv
  - output/phase4_perturbation/phase4_accessibility_note.md

Usage:
    python -m egfr_pipeline.phase4.state_interpretation [--output_dir ...]
"""

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
PHASE4_OUTPUT_DIR = PROJECT_ROOT / "output" / "phase4_perturbation"

# ---------------------------------------------------------------------------
# State interpretation labels
# ---------------------------------------------------------------------------

STATE_LABELS = {
    "robust_pocket": {
        "interpretation": "persistent",
        "description": (
            "Pocket is present across all tested receptor states. "
            "High confidence in accessibility regardless of conformational state."
        ),
        "confidence_modifier": "strengthen",
    },
    "shifted_pocket": {
        "interpretation": "conditionally_accessible",
        "description": (
            "Pocket geometry shifts between states but remains present. "
            "May require state-aware targeting strategy."
        ),
        "confidence_modifier": "neutral",
    },
    "state_specific_pocket": {
        "interpretation": "state_dependent",
        "description": (
            "Pocket detected in only one conformational state. "
            "Could indicate a cryptic site or incomplete sampling. "
            "Not penalized but flagged for cautious interpretation."
        ),
        "confidence_modifier": "flag",
    },
}

# Output schema
INTERPRETATION_COLUMNS = [
    "rank",
    "pocket_id",
    "receptor_id",
    "ligand_id",
    "mechanistic_class",
    "perturbation_score",
    # State interpretation
    "state_class",
    "n_states_matched",
    "state_interpretation",
    "accessibility_label",
    "confidence_modifier",
    "state_caveat",
    # Adjusted confidence
    "classification_confidence",
    "adjusted_confidence",
    # Context
    "A4_state_robustness",
    "diversity_verdict",
    "coverage_fraction",
]


# ---------------------------------------------------------------------------
# 4.4.1: Cross-state support interpretation
# ---------------------------------------------------------------------------

def interpret_state(row: dict) -> dict:
    """Interpret state evidence for a single candidate.

    Returns dict with state interpretation fields.
    """
    state_class = row.get("state_class", "")
    n_matched = int(_safe_float(row.get("n_states_matched", "1")))

    label_info = STATE_LABELS.get(state_class, {
        "interpretation": "unknown",
        "description": "State class not recognized",
        "confidence_modifier": "flag",
    })

    # Accessibility label
    if n_matched >= 3:
        accessibility = "always_accessible"
    elif n_matched == 2:
        accessibility = "mostly_accessible"
    else:
        accessibility = "single_state_only"

    # State caveat
    caveat = ""
    if state_class == "state_specific_pocket":
        if n_matched == 1:
            caveat = (
                "Only 1 receptor state has pocket data. "
                "Cannot determine if state-specific or under-sampled."
            )
    elif state_class == "shifted_pocket":
        caveat = "Pocket geometry varies between states; box coordinates may need per-state adjustment."

    return {
        "state_interpretation": label_info["interpretation"],
        "accessibility_label": accessibility,
        "confidence_modifier": label_info["confidence_modifier"],
        "state_caveat": caveat,
    }


# ---------------------------------------------------------------------------
# 4.4.2: Accessibility-aware confidence adjustment
# ---------------------------------------------------------------------------

def adjust_confidence(
    base_confidence: str,
    modifier: str,
    state_class: str,
    n_states: int,
) -> str:
    """Adjust classification confidence based on state evidence.

    Rules:
    - robust_pocket (3 states): can upgrade medium→high
    - shifted_pocket (2 states): no change
    - state_specific (1 state): no downgrade, but add qualifier
    """
    conf_order = ["low", "medium", "high"]

    if modifier == "strengthen" and n_states >= 3:
        idx = conf_order.index(base_confidence) if base_confidence in conf_order else 0
        return conf_order[min(idx + 1, 2)]

    if modifier == "flag" and n_states == 1:
        # Don't downgrade, but mark as provisional
        return f"{base_confidence}_provisional"

    return base_confidence


# ---------------------------------------------------------------------------
# 4.4.3: Cryptic/allosteric caution
# ---------------------------------------------------------------------------

def check_cryptic_caution(row: dict) -> str:
    """Flag potential cryptic sites that deserve cautious treatment."""
    state_class = row.get("state_class", "")
    mech_class = row.get("mechanistic_class", "")
    n_matched = int(_safe_float(row.get("n_states_matched", "1")))

    if (state_class == "state_specific_pocket" and
            mech_class == "allosteric_modulator_candidate"):
        return (
            "CAUTION: State-specific allosteric site — may be cryptic pocket "
            "only accessible in specific conformational states. "
            "Retain as candidate but interpret with care."
        )

    if n_matched == 1:
        return (
            "NOTE: Single-state data. Cross-state persistence unknown. "
            "Re-evaluate when additional receptor states have pocket data."
        )

    return ""


# ---------------------------------------------------------------------------
# Process all candidates
# ---------------------------------------------------------------------------

def interpret_all(candidates: List[dict], evidence: List[dict]) -> List[dict]:
    """Interpret state evidence for all candidates."""
    # Build evidence lookup for additional context
    evidence_map = {}
    for e in evidence:
        key = (e["pocket_id"], e.get("ligand_id", ""))
        evidence_map[key] = e

    results = []
    for row in candidates:
        key = (row["pocket_id"], row.get("ligand_id", ""))
        ev = evidence_map.get(key, {})

        # Merge evidence fields not in candidate table
        merged = dict(row)
        for field in ("state_class", "n_states_matched", "diversity_verdict",
                      "coverage_fraction", "A4_state_robustness"):
            if not merged.get(field):
                merged[field] = ev.get(field, "")

        state_info = interpret_state(merged)
        cryptic_note = check_cryptic_caution(merged)
        if cryptic_note and state_info["state_caveat"]:
            state_info["state_caveat"] += " " + cryptic_note
        elif cryptic_note:
            state_info["state_caveat"] = cryptic_note

        base_conf = merged.get("classification_confidence", "low")
        adjusted = adjust_confidence(
            base_conf,
            state_info["confidence_modifier"],
            merged.get("state_class", ""),
            int(_safe_float(merged.get("n_states_matched", "1"))),
        )

        result = {
            "rank": merged.get("rank", ""),
            "pocket_id": merged["pocket_id"],
            "receptor_id": merged["receptor_id"],
            "ligand_id": merged.get("ligand_id", ""),
            "mechanistic_class": merged.get("mechanistic_class", ""),
            "perturbation_score": merged.get("perturbation_score", ""),
            "state_class": merged.get("state_class", ""),
            "n_states_matched": merged.get("n_states_matched", ""),
            "state_interpretation": state_info["state_interpretation"],
            "accessibility_label": state_info["accessibility_label"],
            "confidence_modifier": state_info["confidence_modifier"],
            "state_caveat": state_info["state_caveat"],
            "classification_confidence": base_conf,
            "adjusted_confidence": adjusted,
            "A4_state_robustness": merged.get("A4_state_robustness", ""),
            "diversity_verdict": merged.get("diversity_verdict", ""),
            "coverage_fraction": merged.get("coverage_fraction", ""),
        }
        results.append(result)

    return results


# ---------------------------------------------------------------------------
# Accessibility note
# ---------------------------------------------------------------------------

def build_accessibility_note(interpreted: List[dict]) -> List[str]:
    """Build the accessibility interpretation documentation."""
    lines = [
        "# Phase 4 State-Robustness and Accessibility Note",
        "",
        "## 1. State Interpretation Framework",
        "",
        "| State Class | Interpretation | Confidence Effect | Description |",
        "|-------------|---------------|-------------------|-------------|",
    ]
    for sc, info in STATE_LABELS.items():
        lines.append(
            f"| {sc} | {info['interpretation']} | "
            f"{info['confidence_modifier']} | "
            f"{info['description'][:70]}... |"
        )
    lines.append("")

    lines.extend([
        "## 2. Design Principles",
        "",
        "- **State-specific sites are NOT auto-penalized**: They may represent "
        "cryptic pockets or under-sampled conformations.",
        "- **Robust sites get confidence boost**: Persistence across states "
        "strengthens the evidence.",
        "- **Single-state data is flagged**: When only 1 receptor state has "
        "pocket data, all pockets are inherently state_specific. This is a "
        "data limitation, not a biological conclusion.",
        "",
    ])

    # Current results
    lines.extend([
        "## 3. Interpretation Results",
        "",
        "| Rank | Pocket | Accessibility | Adjusted Confidence | Caveat |",
        "|------|--------|--------------|--------------------|---------| ",
    ])
    seen = set()
    for r in interpreted:
        key = (r["pocket_id"], r.get("ligand_id", ""))
        if key in seen:
            continue
        seen.add(key)
        lig = r.get("ligand_id", "") or "(none)"
        caveat_short = r["state_caveat"][:60] + "..." if len(r["state_caveat"]) > 60 else r["state_caveat"]
        lines.append(
            f"| {r['rank']} | {r['pocket_id']} / {lig} | "
            f"{r['accessibility_label']} | "
            f"{r['adjusted_confidence']} | {caveat_short} |"
        )
    lines.append("")

    # Confidence adjustments
    adjusted = [r for r in interpreted
                if r["adjusted_confidence"] != r["classification_confidence"]]
    if adjusted:
        lines.extend([
            "## 4. Confidence Adjustments",
            "",
        ])
        for r in adjusted:
            lines.append(
                f"- {r['pocket_id']}: {r['classification_confidence']} → "
                f"{r['adjusted_confidence']} ({r['confidence_modifier']})")
        lines.append("")
    else:
        lines.extend([
            "## 4. Confidence Adjustments",
            "",
            "No confidence adjustments applied in current data. "
            "This is expected when only 1 receptor state has pocket data.",
            "",
        ])

    # Caveats
    caveats = [r for r in interpreted if r["state_caveat"]]
    if caveats:
        lines.extend([
            "## 5. Active Caveats",
            "",
        ])
        seen_caveats = set()
        for r in caveats:
            if r["state_caveat"] not in seen_caveats:
                seen_caveats.add(r["state_caveat"])
                lines.append(f"- {r['state_caveat']}")
        lines.append("")

    lines.extend([
        "---",
        "",
        "Generated by `egfr_pipeline.phase4.state_interpretation`",
    ])

    return lines


# ---------------------------------------------------------------------------
# Main orchestrator
# ---------------------------------------------------------------------------

def run_state_interpretation(
    output_dir: Path = PHASE4_OUTPUT_DIR,
) -> Tuple[Path, Path]:
    """Run full TG 4.4 pipeline.

    Returns (interpretation_path, note_path).
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load candidate table from TG 4.3
    cand_path = output_dir / "perturbation_candidate_table.csv"
    if not cand_path.exists():
        raise FileNotFoundError(
            f"Candidate table not found: {cand_path}\n"
            "Run TG 4.3 (perturbation_scoring) first."
        )
    candidates = _load_csv(cand_path)

    # Load normalized evidence for additional context
    ev_path = output_dir / "phase4_evidence_normalized.csv"
    evidence = _load_csv(ev_path) if ev_path.exists() else []

    print(f"  Loaded {len(candidates)} candidates, "
          f"{len(evidence)} evidence rows")

    # Interpret
    interpreted = interpret_all(candidates, evidence)

    # Summary
    by_access = defaultdict(int)
    for r in interpreted:
        by_access[r["accessibility_label"]] += 1
    for label, count in sorted(by_access.items()):
        print(f"    {label}: {count}")

    # Write outputs
    interp_path = output_dir / "phase4_state_interpretation.csv"
    _write_csv(interp_path, INTERPRETATION_COLUMNS, interpreted)

    note_lines = build_accessibility_note(interpreted)
    note_path = output_dir / "phase4_accessibility_note.md"
    with open(note_path, "w", encoding="utf-8") as f:
        f.write("\n".join(note_lines))

    return interp_path, note_path


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _safe_float(val) -> float:
    try:
        return float(val)
    except (ValueError, TypeError):
        return 0.0


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
        description="Phase 4 TG 4.4: State-Robustness Interpretation"
    )
    parser.add_argument(
        "--output_dir", type=Path, default=PHASE4_OUTPUT_DIR,
        help="Phase 4 output directory",
    )
    args = parser.parse_args()

    print("Phase 4 — Task Group 4.4: State-Robustness Interpretation")
    print(f"  Output: {args.output_dir}")
    print()

    interp_path, note_path = run_state_interpretation(args.output_dir)

    print(f"\n{'='*60}")
    print("TG 4.4 COMPLETE")
    print(f"{'='*60}")
    print(f"  State interpretation: {interp_path}")
    print(f"  Accessibility note:   {note_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
