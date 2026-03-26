#!/usr/bin/env python3
"""Phase 1 Task 1.7: Pilot data comparison layer.

Compares new full-kinase-domain docking results against existing C-lobe fragment
pilot data. This comparison is INFORMATIONAL ONLY — pilot data must not constrain
interpretation of new results.

Critical framing rule:
  The new full-kinase-domain + extended-beta-meander docking is a fresh start with
  a fundamentally different system. Pilot data (C-lobe fragment 45 res × truncated
  beta-meander 962–1006) carries known systematic artifacts:
    - N-lobe absence (45 res fragment vs ~309 res monomer)
    - VAL962 N-terminal artifact (truncation at 962, not 955)
    - No orientation filtering

Tasks:
  1.7.1 - Register pilot data as historical reference (done in TG 1.0)
  1.7.2 - Optional methodological comparison
  1.7.3 - N-terminal artifact resolution check (VAL962)

Usage:
    python -m egfr_pipeline.phase1.pilot_comparison \
        [--output_dir output/workflow_b/phase1_ppi_analysis]
"""

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

from egfr_pipeline import paths

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

DEFAULT_PATH_CONFIG = {"output_root": str(PROJECT_ROOT / "output")}
PHASE1_OUTPUT_DIR = paths.wb_phase1_ppi_analysis(DEFAULT_PATH_CONFIG)
PHASE1_INPUT_DIR = paths.wa_phase2_runtime_inputs(DEFAULT_PATH_CONFIG)
RECEPTOR_STATES = ["3GT8_raw", "EGFR_160-185", "EGFR_170-200"]

# Known pilot artifacts
PILOT_ARTIFACTS = {
    "N-lobe_absence": (
        "Pilot used C-lobe fragment (~45 res, 962-1006 region only). "
        "N-lobe residues (699-837) could not be detected."
    ),
    "VAL962_terminal": (
        "Pilot partner started at VAL962 (truncation artifact). The extended "
        "beta-meander starts at SER955, eliminating the terminal bias."
    ),
    "no_orientation_filter": (
        "Pilot did not apply orientation-aware filtering. Face-flipped poses "
        "may have contaminated pilot consensus."
    ),
    "dimer_receptor": (
        "Pilot used EGFR dimer (both chains merged, ~600 res with +1000 offset). "
        "Phase 1 uses monomer (chain A only, ~309 res)."
    ),
}

# Pilot clusters of interest (from historical reference)
PILOT_VALID_CLUSTERS = ["C02", "C04", "C07"]
PILOT_VALID_MODELS = ["C02_M01", "C02_M03", "C04_M01", "C04_M02", "C07_M03"]

# VAL962 — the residue to check for N-terminal artifact resolution
ARTIFACT_RESIDUE = "VAL962"
ARTIFACT_RESNUM = 962


# ---------------------------------------------------------------------------
# Comparison logic
# ---------------------------------------------------------------------------

def _is_historical_reference_entry(row: dict) -> bool:
    """Return True only for explicit non-runtime historical reference rows."""
    status_ok = row.get("status") == "historical_reference_only"
    historical_only_ok = str(row.get("historical_only", "")).strip().lower() == "true"
    has_namespace_ref = bool(
        row.get("historical_reference.results_dir") or row.get("legacy_reference.results_dir")
    )
    return status_ok and historical_only_ok and has_namespace_ref


def generate_pilot_comparison(
    output_base: Path,
    states: List[str] = None,
) -> Path:
    """Generate pilot comparison note.

    Reads:
      - pilot_data_reference.csv (from TG 1.0)
      - ppi_interface_patch_table.csv (per state, from TG 1.3)
      - ppi_hotspot_residues.csv (per state, from TG 1.3)

    Returns path to comparison note.
    """
    if states is None:
        states = RECEPTOR_STATES

    # Load pilot reference
    pilot_ref_path = PHASE1_INPUT_DIR / "pilot_data_reference.csv"
    legacy_pilot_ref_path = PROJECT_ROOT / "input" / "PPI" / "phase1" / "pilot_data_reference.csv"
    pilot_entries = []
    if pilot_ref_path.exists():
        with open(pilot_ref_path, encoding="utf-8") as f:
            pilot_entries = [row for row in csv.DictReader(f) if _is_historical_reference_entry(row)]
    elif legacy_pilot_ref_path.exists():
        with open(legacy_pilot_ref_path, encoding="utf-8") as f:
            pilot_entries = [row for row in csv.DictReader(f) if _is_historical_reference_entry(row)]

    # Load new results per state
    state_data = {}
    for state in states:
        state_dir = output_base / state
        patches = _load_csv(state_dir / "ppi_interface_patch_table.csv")
        hotspots = _load_csv(state_dir / "ppi_hotspot_residues.csv")
        if patches or hotspots:
            state_data[state] = {"patches": patches, "hotspots": hotspots}

    # --- VAL962 artifact check (Task 1.7.3) ---
    val962_status = check_val962_artifact(state_data)

    # --- Build comparison note ---
    note_path = output_base / "phase1_pilot_comparison_note.md"
    lines = _build_comparison_note(pilot_entries, state_data, val962_status)

    note_path.parent.mkdir(parents=True, exist_ok=True)
    with open(note_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))

    return note_path


def check_val962_artifact(
    state_data: Dict[str, dict],
) -> Dict[str, dict]:
    """Check whether VAL962 appears as a dominant contact in new results.

    Returns per-state VAL962 status.
    """
    result = {}
    for state, data in state_data.items():
        patches = data.get("patches", [])
        hotspots = data.get("hotspots", [])

        # Find VAL962 in patch table (partner chain B)
        val962_patch = None
        for p in patches:
            if p.get("chain") == "B" and p.get("residue_id") == ARTIFACT_RESIDUE:
                val962_patch = p
                break

        # Find VAL962 in hotspot table
        val962_hotspot = False
        val962_max_occ = 0.0
        for h in hotspots:
            if (h.get("chain") == "B" and h.get("residue_id") == ARTIFACT_RESIDUE
                    and h.get("is_hotspot") in ("True", "true", True)):
                val962_hotspot = True
                occ = _safe_float(h.get("occupancy", ""))
                if occ and occ > val962_max_occ:
                    val962_max_occ = occ

        if val962_patch:
            max_occ = _safe_float(val962_patch.get("max_occupancy", "")) or 0
            n_clusters = int(val962_patch.get("n_clusters_present", 0) or 0)
        else:
            max_occ = 0
            n_clusters = 0

        # Determine status
        if val962_patch is None:
            status = "absent"
            interpretation = "VAL962 not detected in any cluster → artifact resolved"
        elif val962_hotspot:
            status = "dominant_contact"
            interpretation = (
                f"VAL962 is a hotspot (occ={val962_max_occ:.2f}) in {n_clusters} cluster(s). "
                "This could be independent confirmation of genuine binding or residual artifact. "
                "Compare with other partner residues for context."
            )
        elif max_occ > 0.3:
            status = "moderate_contact"
            interpretation = (
                f"VAL962 present (occ={max_occ:.2f}) but not dominant. "
                "Artifact partially resolved — not a primary contact."
            )
        else:
            status = "minor_contact"
            interpretation = (
                f"VAL962 present at low occupancy ({max_occ:.2f}). "
                "Artifact effectively resolved — background-level contact."
            )

        result[state] = {
            "status": status,
            "max_occupancy": max_occ,
            "n_clusters": n_clusters,
            "is_hotspot": val962_hotspot,
            "interpretation": interpretation,
        }

    return result


def _build_comparison_note(
    pilot_entries: List[dict],
    state_data: Dict[str, dict],
    val962_status: Dict[str, dict],
) -> List[str]:
    """Build the pilot comparison markdown note."""
    lines = [
        "# Phase 1 Pilot Data Comparison Note",
        "",
        "## Framing",
        "",
        "This comparison is **informational only**. New full-kinase-domain results",
        "are interpreted independently of pilot expectations.",
        "",
        "The pilot campaign used a fundamentally different system:",
        "",
    ]

    # Known artifacts table
    lines.extend([
        "### Known Pilot Artifacts",
        "",
        "| Artifact | Description |",
        "|----------|-------------|",
    ])
    for name, desc in PILOT_ARTIFACTS.items():
        lines.append(f"| {name} | {desc} |")
    lines.append("")

    # Pilot registration
    lines.extend([
        "## Pilot Data Registration (Task 1.7.1)",
        "",
    ])
    if pilot_entries:
        lines.extend([
            "| Label | Construct | Partner range | Artifacts |",
            "|-------|-----------|---------------|-----------|",
        ])
        for e in pilot_entries:
            lines.append(
                f"| {e.get('label', '')} | {e.get('construct_type', '')} | "
                f"{e.get('partner_range', '')} | {e.get('known_artifacts', '')} |"
            )
        lines.append("")
        lines.append(f"Status: All pilot data labeled as `{pilot_entries[0].get('status', 'historical_reference_only')}`.")
        lines.append("Pilot data is NOT integrated into new consensus or patch definitions.")
    else:
        lines.append("No pilot reference data found.")
    lines.append("")

    # VAL962 artifact resolution (Task 1.7.3)
    lines.extend([
        "## VAL962 N-Terminal Artifact Resolution (Task 1.7.3)",
        "",
        "In the pilot campaign, VAL962 was the first residue of the truncated",
        "beta-meander (962–1006). Terminal residues in docking can accumulate",
        "artifactual contacts due to exposed backbone and increased flexibility.",
        "",
        "The extended beta-meander (955–1006) eliminates this artifact by",
        "placing VAL962 in the middle of the chain, 7 residues from the N-terminus.",
        "",
    ])

    if val962_status:
        lines.extend([
            "### Per-State VAL962 Status",
            "",
            "| State | Status | Max occupancy | Clusters | Hotspot | Interpretation |",
            "|-------|--------|--------------|----------|---------|---------------|",
        ])
        for state in RECEPTOR_STATES:
            v = val962_status.get(state)
            if v:
                lines.append(
                    f"| {state} | {v['status']} | {v['max_occupancy']:.2f} | "
                    f"{v['n_clusters']} | {'Yes' if v['is_hotspot'] else 'No'} | "
                    f"{v['interpretation']} |"
                )
            else:
                lines.append(f"| {state} | no data | — | — | — | — |")
        lines.append("")

        # Overall assessment
        statuses = [v["status"] for v in val962_status.values()]
        if all(s == "absent" for s in statuses):
            assessment = ("VAL962 is absent from all states → N-terminal artifact "
                          "fully resolved by extended beta-meander construct.")
        elif all(s in ("absent", "minor_contact") for s in statuses):
            assessment = ("VAL962 appears at background levels only → artifact "
                          "effectively resolved.")
        elif any(s == "dominant_contact" for s in statuses):
            assessment = ("VAL962 remains a dominant contact in some states. This "
                          "may indicate genuine binding (independent confirmation) "
                          "or residual proximity effects. Compare with neighboring "
                          "partner residues for context.")
        else:
            assessment = "Mixed results — see per-state details above."

        lines.extend([
            "### Overall Assessment",
            "",
            assessment,
            "",
        ])
    else:
        lines.append("No new docking data available for VAL962 analysis.")
        lines.append("")

    # Methodological comparison (Task 1.7.2) — brief, not constraining
    if state_data:
        lines.extend([
            "## Methodological Notes (Task 1.7.2)",
            "",
            "### System Differences",
            "",
            "| Parameter | Pilot | Phase 1 |",
            "|-----------|-------|---------|",
            "| Receptor | C-lobe fragment (~45 res) | Full kinase domain (~309 res) |",
            "| Partner | Truncated beta-meander (962–1006) | Extended beta-meander (955–1006) |",
            "| Receptor construct | EGFR dimer (merged chains) | EGFR monomer (chain A) |",
            "| Orientation filter | None | Mandatory (dual-vector PCA) |",
            "| N-lobe coverage | Impossible | Available (699–837) |",
            "",
            "### Expected Differences in Results",
            "",
            "- **New N-lobe contacts**: Expected, since the pilot could not detect them",
            "- **Different C-lobe patch shape**: Expected, since N-lobe presence changes the",
            "  steric landscape and may redirect partner binding",
            "- **Different cluster ranking**: Expected, since the energy landscape is",
            "  fundamentally different with a larger receptor",
            "- **Fewer face-flipped poses**: Expected, due to mandatory orientation filtering",
            "",
            "These differences are **expected improvements**, not failures.",
            "",
        ])

    lines.extend([
        "---",
        "",
        "Generated by `egfr_pipeline.phase1.pilot_comparison`",
    ])

    return lines


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _load_csv(path: Path) -> List[dict]:
    if not path.exists():
        return []
    with open(path, encoding="utf-8") as f:
        return list(csv.DictReader(f))


def _safe_float(val: str) -> Optional[float]:
    try:
        return float(val)
    except (ValueError, TypeError):
        return None


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Phase 1 TG 1.7: Pilot data comparison layer"
    )
    parser.add_argument("--output_dir", type=Path,
                        default=PHASE1_OUTPUT_DIR,
                        help="Workflow B Phase 1 analysis output directory")
    args = parser.parse_args()

    print("Phase 1 — Task 1.7: Pilot Data Comparison Layer")
    note_path = generate_pilot_comparison(args.output_dir)
    print(f"\nComparison note: {note_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
