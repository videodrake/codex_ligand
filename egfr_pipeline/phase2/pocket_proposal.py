#!/usr/bin/env python3
"""Phase 2 Task Group 2.1: Multi-Tool Candidate Pocket Proposal.

Enumerates candidate ligandable pockets for each receptor state using fpocket
and/or P2Rank. Since these tools may not be available in the current environment,
this module provides:

  1. Setup generation: shell scripts + metadata for server-side execution
  2. Output parsers: fpocket and P2Rank output → tool-agnostic pocket schema
  3. Source summary: per-tool statistics

Tasks:
  2.1.1 - Pocket proposal source integration (fpocket + P2Rank parsers)
  2.1.2 - Receptor-local pocket proposal (per-state execution)
  2.1.3 - Pocket metadata extraction (centroid, score, box size)

Inputs:
  - input/PPI/phase1/receptor_*.pdb (3 receptor states)
  - fpocket/P2Rank output directories (after server-side execution)

Outputs:
  - output/phase2_pockets/pocket_proposal_run.sh
  - output/phase2_pockets/pocket_proposal_metadata.json
  - output/phase2_pockets/candidate_pockets_raw.csv
  - output/phase2_pockets/candidate_pocket_source_summary.csv
  - output/phase2_pockets/phase2_pocket_proposal_note.md

Usage:
    # Generate run scripts (before server execution)
    python -m egfr_pipeline.phase2.pocket_proposal --setup

    # Parse results (after server execution)
    python -m egfr_pipeline.phase2.pocket_proposal --parse

    # Both
    python -m egfr_pipeline.phase2.pocket_proposal --setup --parse
"""

import argparse
import csv
import json
import os
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

PHASE1_INPUT_DIR = PROJECT_ROOT / "input" / "PPI" / "phase1"
PHASE2_OUTPUT_DIR = PROJECT_ROOT / "output" / "phase2_pockets"

RECEPTOR_STATES = ["3GT8_raw", "3GT8_cl38_48", "3GT8_cl85_100"]
RECEPTOR_PDB_PATTERN = "receptor_{state}.pdb"

# Supported pocket proposal tools
SUPPORTED_TOOLS = ["fpocket", "p2rank"]

# Tool-agnostic candidate pocket schema
CANDIDATE_POCKET_COLUMNS = [
    "receptor_id",
    "candidate_pocket_id",       # e.g., 3GT8_raw_fpocket_P01
    "proposal_source",           # fpocket / p2rank
    "pocket_rank",               # tool-native rank (1-based)
    "proposal_score",            # tool-native score
    "druggability_score",        # tool-native druggability (if available)
    "centroid_x",
    "centroid_y",
    "centroid_z",
    "volume_A3",                 # pocket volume in Å³
    "n_alpha_spheres",           # fpocket: alpha sphere count; p2rank: pocket points
    "residue_ids",               # semicolon-delimited residue IDs (e.g., LEU838;ASP855)
    "residue_nums",              # semicolon-delimited residue numbers
    "n_residues",                # number of lining residues
    "box_center_x",              # suggested docking box center
    "box_center_y",
    "box_center_z",
    "box_size_x",                # suggested box dimensions
    "box_size_y",
    "box_size_z",
]

SOURCE_SUMMARY_COLUMNS = [
    "receptor_id",
    "proposal_source",
    "n_pockets_found",
    "best_score",
    "best_druggability",
    "mean_volume_A3",
    "output_dir",
    "status",                    # found / not_found / parse_error
]

# fpocket output structure
FPOCKET_OUT_SUFFIX = "_out"      # e.g., receptor_3GT8_raw_out/
FPOCKET_INFO_FILE = "*_info.txt"
FPOCKET_POCKET_PDB_PATTERN = "pocket{n}_atm.pdb"

# P2Rank output structure
P2RANK_PREDICTIONS_CSV = "predictions.csv"
P2RANK_RESIDUES_CSV = "residues.csv"

# Box padding around pocket centroid (Å)
DEFAULT_BOX_PADDING = 5.0


# ---------------------------------------------------------------------------
# Setup generation (Task 2.1.2)
# ---------------------------------------------------------------------------

def generate_pocket_proposal_setup(
    output_dir: Path,
    receptor_dir: Path = PHASE1_INPUT_DIR,
) -> Tuple[Path, Path]:
    """Generate run scripts and metadata for server-side pocket proposal.

    Returns (script_path, metadata_path).
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # Locate receptor PDBs
    receptors = {}
    for state in RECEPTOR_STATES:
        pdb_name = RECEPTOR_PDB_PATTERN.format(state=state)
        pdb_path = receptor_dir / pdb_name
        if pdb_path.exists():
            receptors[state] = pdb_path
        else:
            print(f"  WARNING: receptor PDB not found: {pdb_path}")

    if not receptors:
        raise FileNotFoundError(
            f"No receptor PDBs found in {receptor_dir}. "
            "Run Phase 1 TG 1.0 first."
        )

    # Build run script
    script_lines = _build_run_script(receptors, output_dir)
    script_path = output_dir / "pocket_proposal_run.sh"
    with open(script_path, "w", encoding="utf-8") as f:
        f.write("\n".join(script_lines))
    os.chmod(script_path, 0o755)

    # Build metadata
    metadata = {
        "phase": "phase2",
        "task_group": "2.1",
        "description": "Multi-tool candidate pocket proposal",
        "tools": SUPPORTED_TOOLS,
        "receptor_states": list(receptors.keys()),
        "receptor_pdbs": {s: str(p) for s, p in receptors.items()},
        "output_dir": str(output_dir),
        "expected_outputs": {
            "fpocket": {s: str(output_dir / s / f"receptor_{s}{FPOCKET_OUT_SUFFIX}")
                        for s in receptors},
            "p2rank": {s: str(output_dir / s / "p2rank_output")
                       for s in receptors},
        },
    }
    meta_path = output_dir / "pocket_proposal_metadata.json"
    with open(meta_path, "w", encoding="utf-8") as f:
        json.dump(metadata, f, indent=2)

    return script_path, meta_path


def _build_run_script(
    receptors: Dict[str, Path],
    output_dir: Path,
) -> List[str]:
    """Build shell script for fpocket + P2Rank execution."""
    lines = [
        "#!/bin/bash",
        "# Phase 2 TG 2.1: Pocket Proposal Execution Script",
        "# Generated by egfr_pipeline.phase2.pocket_proposal",
        "#",
        "# Run this script on the production server where fpocket and/or P2Rank",
        "# are installed. Both tools are optional — the parser handles missing output.",
        "#",
        "# Usage:",
        "#   bash pocket_proposal_run.sh",
        "#   bash pocket_proposal_run.sh fpocket    # fpocket only",
        "#   bash pocket_proposal_run.sh p2rank     # P2Rank only",
        "",
        "set -euo pipefail",
        "",
        f'OUTPUT_BASE="{output_dir}"',
        'TOOL="${1:-all}"',
        "",
        "echo '=== Phase 2 TG 2.1: Pocket Proposal ==='",
        "echo",
        "",
    ]

    # fpocket section
    lines.extend([
        '# --- fpocket ---',
        'if [[ "$TOOL" == "all" || "$TOOL" == "fpocket" ]]; then',
        '  echo "--- Running fpocket ---"',
        '  if ! command -v fpocket &>/dev/null; then',
        '    echo "WARNING: fpocket not found in PATH. Skipping."',
        '  else',
    ])

    for state, pdb_path in receptors.items():
        state_dir = output_dir / state
        lines.extend([
            f'    echo "  State: {state}"',
            f'    mkdir -p "{state_dir}"',
            f'    cp "{pdb_path}" "{state_dir}/receptor_{state}.pdb"',
            f'    cd "{state_dir}"',
            f'    fpocket -f "receptor_{state}.pdb"',
            f'    cd "$OUTPUT_BASE"',
            f'    echo "    Done: {state}"',
        ])

    lines.extend([
        '  fi',
        '  echo',
        'fi',
        '',
    ])

    # P2Rank section
    lines.extend([
        '# --- P2Rank ---',
        'if [[ "$TOOL" == "all" || "$TOOL" == "p2rank" ]]; then',
        '  echo "--- Running P2Rank ---"',
        '  if ! command -v prank &>/dev/null; then',
        '    echo "WARNING: prank (P2Rank) not found in PATH. Skipping."',
        '  else',
    ])

    for state, pdb_path in receptors.items():
        state_dir = output_dir / state
        p2rank_out = state_dir / "p2rank_output"
        lines.extend([
            f'    echo "  State: {state}"',
            f'    mkdir -p "{state_dir}"',
            f'    prank predict -f "{pdb_path}" -o "{p2rank_out}"',
            f'    echo "    Done: {state}"',
        ])

    lines.extend([
        '  fi',
        '  echo',
        'fi',
        '',
        'echo "=== Pocket proposal complete ==="',
        'echo "Run parser: python -m egfr_pipeline.phase2.pocket_proposal --parse"',
    ])

    return lines


# ---------------------------------------------------------------------------
# fpocket parser (Task 2.1.1 + 2.1.3)
# ---------------------------------------------------------------------------

def parse_fpocket_output(
    state: str,
    state_dir: Path,
) -> List[dict]:
    """Parse fpocket output for one receptor state.

    fpocket creates: receptor_{state}_out/
      - receptor_{state}_out.pdb  (all pockets)
      - receptor_{state}_info.txt  (pocket summary)
      - pockets/
          pocket0001_atm.pdb ... pocketNNNN_atm.pdb

    Returns list of candidate pocket dicts.
    """
    fpocket_dir = state_dir / f"receptor_{state}{FPOCKET_OUT_SUFFIX}"
    if not fpocket_dir.exists():
        return []

    # Parse info file for pocket-level stats
    info_files = list(fpocket_dir.glob("*_info.txt"))
    if not info_files:
        return []

    pocket_stats = _parse_fpocket_info(info_files[0])
    if not pocket_stats:
        return []

    # Parse pocket PDB files for residues and centroids
    pockets_dir = fpocket_dir / "pockets"
    results = []

    for i, stats in enumerate(pocket_stats, start=1):
        pocket_pdb = pockets_dir / f"pocket{i}_atm.pdb"
        residues = []
        coords = []

        if pocket_pdb.exists():
            residues, coords = _parse_pocket_pdb(pocket_pdb)

        # Compute centroid from alpha sphere centers or pocket atoms
        if coords:
            cx = sum(c[0] for c in coords) / len(coords)
            cy = sum(c[1] for c in coords) / len(coords)
            cz = sum(c[2] for c in coords) / len(coords)
        else:
            cx, cy, cz = 0.0, 0.0, 0.0

        # Compute bounding box
        if coords:
            box_size_x = max(c[0] for c in coords) - min(c[0] for c in coords) + 2 * DEFAULT_BOX_PADDING
            box_size_y = max(c[1] for c in coords) - min(c[1] for c in coords) + 2 * DEFAULT_BOX_PADDING
            box_size_z = max(c[2] for c in coords) - min(c[2] for c in coords) + 2 * DEFAULT_BOX_PADDING
        else:
            box_size_x = box_size_y = box_size_z = 0.0

        # Residue IDs
        res_ids = sorted(set(residues), key=lambda r: _resnum_from_id(r))
        res_nums = [str(_resnum_from_id(r)) for r in res_ids]

        pocket_id = f"{state}_fpocket_P{i:02d}"

        results.append({
            "receptor_id": state,
            "candidate_pocket_id": pocket_id,
            "proposal_source": "fpocket",
            "pocket_rank": i,
            "proposal_score": stats.get("score", ""),
            "druggability_score": stats.get("druggability_score", ""),
            "centroid_x": f"{cx:.3f}",
            "centroid_y": f"{cy:.3f}",
            "centroid_z": f"{cz:.3f}",
            "volume_A3": stats.get("volume", ""),
            "n_alpha_spheres": stats.get("n_alpha_spheres", ""),
            "residue_ids": ";".join(res_ids),
            "residue_nums": ";".join(res_nums),
            "n_residues": len(res_ids),
            "box_center_x": f"{cx:.3f}",
            "box_center_y": f"{cy:.3f}",
            "box_center_z": f"{cz:.3f}",
            "box_size_x": f"{box_size_x:.1f}",
            "box_size_y": f"{box_size_y:.1f}",
            "box_size_z": f"{box_size_z:.1f}",
        })

    return results


def _parse_fpocket_info(info_path: Path) -> List[dict]:
    """Parse fpocket *_info.txt file.

    Format (per pocket block):
        Pocket 1 :
            Score :                0.4521
            Druggability Score :   0.234
            Number of Alpha Spheres :    35
            Total SASA :           123.45
            Polar SASA :            45.67
            Apolar SASA :           77.78
            Volume :               456.78
            ...
    """
    pockets = []
    current = None

    with open(info_path, encoding="utf-8") as f:
        for line in f:
            line = line.rstrip()

            # New pocket header
            m = re.match(r"^Pocket\s+(\d+)\s*:", line)
            if m:
                if current is not None:
                    pockets.append(current)
                current = {"pocket_num": int(m.group(1))}
                continue

            if current is None:
                continue

            # Key-value pairs
            if ":" in line:
                key, _, val = line.partition(":")
                key = key.strip().lower().replace(" ", "_")
                val = val.strip()

                if key == "score":
                    current["score"] = _safe_float(val) or ""
                elif key == "druggability_score":
                    current["druggability_score"] = _safe_float(val) or ""
                elif key == "number_of_alpha_spheres":
                    current["n_alpha_spheres"] = _safe_int(val) or ""
                elif key == "volume":
                    current["volume"] = _safe_float(val) or ""

    if current is not None:
        pockets.append(current)

    return pockets


def _parse_pocket_pdb(pdb_path: Path) -> Tuple[List[str], List[Tuple[float, float, float]]]:
    """Parse a fpocket pocket PDB for lining residues and atom coordinates.

    Returns (residue_ids, coordinates).
    """
    residues = []
    coords = []

    with open(pdb_path, encoding="utf-8") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                resname = line[17:20].strip()
                resnum_str = line[22:26].strip()
                resnum = _safe_int(resnum_str)

                if resname and resnum is not None:
                    residues.append(f"{resname}{resnum}")

                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append((x, y, z))
                except (ValueError, IndexError):
                    pass

    return residues, coords


# ---------------------------------------------------------------------------
# P2Rank parser (Task 2.1.1 + 2.1.3)
# ---------------------------------------------------------------------------

def parse_p2rank_output(
    state: str,
    state_dir: Path,
) -> List[dict]:
    """Parse P2Rank output for one receptor state.

    P2Rank creates: p2rank_output/
      - receptor_{state}.pdb_predictions.csv
      - receptor_{state}.pdb_residues.csv

    predictions.csv columns:
      name, rank, score, probability, sas_points, surf_atoms,
      center_x, center_y, center_z, residue_ids

    Returns list of candidate pocket dicts.
    """
    p2rank_dir = state_dir / "p2rank_output"
    if not p2rank_dir.exists():
        return []

    # Find predictions CSV (P2Rank may use various naming)
    pred_files = list(p2rank_dir.glob("*predictions.csv"))
    if not pred_files:
        return []

    results = []

    with open(pred_files[0], encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rank = _safe_int(row.get("rank", "")) or 0
            score = _safe_float(row.get("score", "")) or 0.0
            prob = _safe_float(row.get("probability", "")) or 0.0

            cx = _safe_float(row.get("center_x", "")) or 0.0
            cy = _safe_float(row.get("center_y", "")) or 0.0
            cz = _safe_float(row.get("center_z", "")) or 0.0

            # P2Rank residue_ids format: "ALA_123 LEU_456 ..."
            raw_residues = row.get("residue_ids", "").strip()
            res_ids = []
            res_nums = []
            if raw_residues:
                for token in raw_residues.split():
                    token = token.strip()
                    if "_" in token:
                        parts = token.split("_")
                        resname = parts[0]
                        resnum = parts[1] if len(parts) > 1 else ""
                        res_ids.append(f"{resname}{resnum}")
                        res_nums.append(resnum)
                    else:
                        res_ids.append(token)

            pocket_id = f"{state}_p2rank_P{rank:02d}"

            # Estimate box size from SAS points count (heuristic)
            sas_points = _safe_int(row.get("sas_points", "")) or 0
            # Rough estimate: box_dim ≈ 2 * (sas_points)^(1/3) + padding
            estimated_dim = 2.0 * (sas_points ** (1.0/3.0)) + 2 * DEFAULT_BOX_PADDING if sas_points > 0 else 0.0

            results.append({
                "receptor_id": state,
                "candidate_pocket_id": pocket_id,
                "proposal_source": "p2rank",
                "pocket_rank": rank,
                "proposal_score": f"{score:.4f}" if score else "",
                "druggability_score": f"{prob:.4f}" if prob else "",
                "centroid_x": f"{cx:.3f}",
                "centroid_y": f"{cy:.3f}",
                "centroid_z": f"{cz:.3f}",
                "volume_A3": "",  # P2Rank doesn't report volume directly
                "n_alpha_spheres": sas_points,
                "residue_ids": ";".join(res_ids),
                "residue_nums": ";".join(res_nums),
                "n_residues": len(res_ids),
                "box_center_x": f"{cx:.3f}",
                "box_center_y": f"{cy:.3f}",
                "box_center_z": f"{cz:.3f}",
                "box_size_x": f"{estimated_dim:.1f}",
                "box_size_y": f"{estimated_dim:.1f}",
                "box_size_z": f"{estimated_dim:.1f}",
            })

    return results


# ---------------------------------------------------------------------------
# Aggregation and summary (Task 2.1.2 + 2.1.3)
# ---------------------------------------------------------------------------

def parse_all_pocket_proposals(
    output_dir: Path,
) -> Tuple[List[dict], List[dict]]:
    """Parse all pocket proposal outputs across states and tools.

    Returns (all_pockets, source_summary).
    """
    all_pockets = []
    source_entries = []

    for state in RECEPTOR_STATES:
        state_dir = output_dir / state

        for tool in SUPPORTED_TOOLS:
            if tool == "fpocket":
                pockets = parse_fpocket_output(state, state_dir)
                tool_dir = state_dir / f"receptor_{state}{FPOCKET_OUT_SUFFIX}"
            elif tool == "p2rank":
                pockets = parse_p2rank_output(state, state_dir)
                tool_dir = state_dir / "p2rank_output"
            else:
                pockets = []
                tool_dir = state_dir / tool

            if pockets:
                all_pockets.extend(pockets)
                scores = [_safe_float(p.get("proposal_score", "")) for p in pockets]
                scores = [s for s in scores if s is not None]
                drug_scores = [_safe_float(p.get("druggability_score", "")) for p in pockets]
                drug_scores = [s for s in drug_scores if s is not None]
                volumes = [_safe_float(p.get("volume_A3", "")) for p in pockets]
                volumes = [v for v in volumes if v is not None]

                source_entries.append({
                    "receptor_id": state,
                    "proposal_source": tool,
                    "n_pockets_found": len(pockets),
                    "best_score": f"{max(scores):.4f}" if scores else "",
                    "best_druggability": f"{max(drug_scores):.4f}" if drug_scores else "",
                    "mean_volume_A3": f"{sum(volumes)/len(volumes):.1f}" if volumes else "",
                    "output_dir": str(tool_dir),
                    "status": "found",
                })
            else:
                source_entries.append({
                    "receptor_id": state,
                    "proposal_source": tool,
                    "n_pockets_found": 0,
                    "best_score": "",
                    "best_druggability": "",
                    "mean_volume_A3": "",
                    "output_dir": str(tool_dir),
                    "status": "not_found",
                })

    return all_pockets, source_entries


def build_proposal_note(
    all_pockets: List[dict],
    source_summary: List[dict],
) -> List[str]:
    """Build the pocket proposal markdown note."""
    lines = [
        "# Phase 2 Pocket Proposal Note",
        "",
        "## Overview",
        "",
        f"Total candidate pockets found: **{len(all_pockets)}**",
        "",
    ]

    # Source summary table
    lines.extend([
        "## Source Summary",
        "",
        "| Receptor | Tool | Pockets | Best Score | Best Druggability | Status |",
        "|----------|------|---------|-----------|-------------------|--------|",
    ])
    for s in source_summary:
        lines.append(
            f"| {s['receptor_id']} | {s['proposal_source']} | "
            f"{s['n_pockets_found']} | {s['best_score']} | "
            f"{s['best_druggability']} | {s['status']} |"
        )
    lines.append("")

    # Per-state pocket summary
    if all_pockets:
        by_state = defaultdict(list)
        for p in all_pockets:
            by_state[p["receptor_id"]].append(p)

        for state in RECEPTOR_STATES:
            pockets = by_state.get(state, [])
            if not pockets:
                lines.append(f"### {state}: No pockets found")
                lines.append("")
                continue

            lines.extend([
                f"### {state}",
                "",
                "| Pocket ID | Source | Rank | Score | Drugg. | Centroid | N res |",
                "|-----------|--------|------|-------|--------|---------|-------|",
            ])
            for p in sorted(pockets, key=lambda x: int(x.get("pocket_rank", 99))):
                centroid = f"({p['centroid_x']}, {p['centroid_y']}, {p['centroid_z']})"
                lines.append(
                    f"| {p['candidate_pocket_id']} | {p['proposal_source']} | "
                    f"{p['pocket_rank']} | {p['proposal_score']} | "
                    f"{p['druggability_score']} | {centroid} | {p['n_residues']} |"
                )
            lines.append("")
    else:
        lines.extend([
            "## No Pockets Found",
            "",
            "Neither fpocket nor P2Rank output was found.",
            "Run `pocket_proposal_run.sh` on the production server first,",
            "then re-run with `--parse`.",
            "",
        ])

    lines.extend([
        "---",
        "",
        "Generated by `egfr_pipeline.phase2.pocket_proposal`",
    ])

    return lines


# ---------------------------------------------------------------------------
# Main orchestrator
# ---------------------------------------------------------------------------

def run_pocket_proposal(
    do_setup: bool = False,
    do_parse: bool = False,
    output_dir: Path = PHASE2_OUTPUT_DIR,
    receptor_dir: Path = PHASE1_INPUT_DIR,
) -> None:
    """Run TG 2.1 pipeline."""
    output_dir.mkdir(parents=True, exist_ok=True)

    if do_setup:
        print("--- Setup: Generating run scripts ---")
        script_path, meta_path = generate_pocket_proposal_setup(
            output_dir, receptor_dir
        )
        print(f"  Run script: {script_path}")
        print(f"  Metadata:   {meta_path}")
        print()

    if do_parse:
        print("--- Parse: Reading pocket proposal outputs ---")
        all_pockets, source_summary = parse_all_pocket_proposals(output_dir)

        n_found = sum(1 for s in source_summary if s["status"] == "found")
        n_total = len(source_summary)
        print(f"  Sources with data: {n_found}/{n_total}")
        print(f"  Total pockets: {len(all_pockets)}")

        # Write candidate_pockets_raw.csv
        raw_path = output_dir / "candidate_pockets_raw.csv"
        with open(raw_path, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=CANDIDATE_POCKET_COLUMNS,
                               extrasaction="ignore")
            w.writeheader()
            w.writerows(all_pockets)
        print(f"  Raw pockets: {raw_path}")

        # Write source summary
        summary_path = output_dir / "candidate_pocket_source_summary.csv"
        with open(summary_path, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=SOURCE_SUMMARY_COLUMNS,
                               extrasaction="ignore")
            w.writeheader()
            w.writerows(source_summary)
        print(f"  Source summary: {summary_path}")

        # Write proposal note
        note_lines = build_proposal_note(all_pockets, source_summary)
        note_path = output_dir / "phase2_pocket_proposal_note.md"
        with open(note_path, "w", encoding="utf-8") as f:
            f.write("\n".join(note_lines))
        print(f"  Proposal note: {note_path}")
        print()


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _safe_float(val) -> Optional[float]:
    try:
        return float(val)
    except (ValueError, TypeError):
        return None


def _safe_int(val) -> Optional[int]:
    try:
        return int(float(val))
    except (ValueError, TypeError):
        return None


def _resnum_from_id(residue_id: str) -> int:
    """Extract residue number from ID like 'ASP855' → 855."""
    m = re.search(r"(\d+)$", residue_id)
    return int(m.group(1)) if m else 0


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Phase 2 TG 2.1: Multi-tool candidate pocket proposal"
    )
    parser.add_argument(
        "--setup", action="store_true",
        help="Generate run scripts for server-side execution",
    )
    parser.add_argument(
        "--parse", action="store_true",
        help="Parse fpocket/P2Rank outputs into candidate pockets",
    )
    parser.add_argument(
        "--output_dir", type=Path, default=PHASE2_OUTPUT_DIR,
        help="Phase 2 output directory",
    )
    parser.add_argument(
        "--receptor_dir", type=Path, default=PHASE1_INPUT_DIR,
        help="Directory containing receptor PDB files",
    )
    args = parser.parse_args()

    if not args.setup and not args.parse:
        print("Specify --setup and/or --parse. Use --help for details.")
        return 1

    print("Phase 2 — Task Group 2.1: Pocket Proposal")
    print(f"  Output: {args.output_dir}")
    print()

    run_pocket_proposal(
        do_setup=args.setup,
        do_parse=args.parse,
        output_dir=args.output_dir,
        receptor_dir=args.receptor_dir,
    )

    print(f"{'='*60}")
    print("TG 2.1 COMPLETE")
    print(f"{'='*60}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
