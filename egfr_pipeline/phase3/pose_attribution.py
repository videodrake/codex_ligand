#!/usr/bin/env python3
"""Phase 3 Task Group 3.4: Pose Parsing Compatibility and Pocket Attribution.

Ensures Phase 3 docking outputs remain compatible with the existing pose
parsing flow while preserving pocket-guided provenance.

Tasks:
  3.4.1 - Pose provenance extension (pocket_id, round_id, seed, docking_mode)
  3.4.2 - Compatibility with existing parse_poses / extract_contacts
  3.4.3 - Pocket-attribution-aware pose table for Phase 4

Inputs:
  - output/phase3_docking/runs/{receptor}/{pocket}/{ligand}/*.pdbqt
  - output/phase3_docking/phase3_docking_job_table.csv
  - output/phase3_docking/phase3_round_log.csv

Outputs:
  - output/phase3_docking/vina_pose_table.csv  (extended with provenance)
  - output/phase3_docking/phase3_pose_provenance_note.md

Usage:
    python -m egfr_pipeline.phase3.pose_attribution [--output_dir ...]
"""

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from egfr_pipeline.csv_utils import load_csv

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

PHASE3_OUTPUT_DIR = PROJECT_ROOT / "output" / "phase3_docking"

# Extended pose table schema — backward-compatible superset of legacy schema
# Legacy fields (13): receptor_id, ligand_id, pose_rank, affinity, rmsd_lb,
#   rmsd_ub, centroid_x, centroid_y, centroid_z, raw_pose_file, pocket_id,
#   contact_residues, n_contact_residues
# Phase 3 additions (5): candidate_pocket_id, round_id, seed, docking_mode,
#   job_id
POSE_TABLE_COLUMNS = [
    # Legacy fields (backward compatible)
    "receptor_id",
    "ligand_id",
    "pose_rank",
    "affinity",
    "rmsd_lb",
    "rmsd_ub",
    "centroid_x",
    "centroid_y",
    "centroid_z",
    "raw_pose_file",
    "pocket_id",             # Legacy clustering pocket (filled later by cluster)
    "contact_residues",      # Filled by contacts enrichment
    "n_contact_residues",    # Filled by contacts enrichment
    # Phase 3 provenance fields (Task 3.4.1)
    "candidate_pocket_id",   # Phase 2/3 pocket identity
    "round_id",              # Docking round (seed_index)
    "seed",                  # Seed index
    "docking_mode",          # "focused" (Phase 3) vs "blind" (legacy)
    "job_id",                # Full job ID for traceability
]


# ---------------------------------------------------------------------------
# 3.4.1: Pose provenance — parse Phase 3 PDBQT outputs
# ---------------------------------------------------------------------------

def parse_pose_blocks(pdbqt_path: Path) -> List[dict]:
    """Parse PDBQT multi-model file into pose dicts.

    Reuses logic from egfr_pipeline.vina.parse_poses.parse_pose_blocks
    but kept standalone to avoid import dependency issues.
    """
    poses = []
    current_coords = []
    current_affinity = None
    current_rmsd_lb = None
    current_rmsd_ub = None

    if not pdbqt_path.exists():
        return poses

    with open(pdbqt_path, encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if line.startswith("REMARK VINA RESULT:"):
                parts = line.split()
                if len(parts) >= 6:
                    try:
                        current_affinity = float(parts[3])
                        current_rmsd_lb = float(parts[4])
                        current_rmsd_ub = float(parts[5])
                    except (ValueError, TypeError):
                        continue
            if line.startswith(("ATOM", "HETATM")):
                atom_name = line[12:16].strip()
                if atom_name.upper().startswith("H"):
                    continue
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                except ValueError:
                    continue
                current_coords.append([x, y, z])
            if line.startswith("ENDMDL"):
                if current_coords:
                    centroid = [
                        sum(c[0] for c in current_coords) / len(current_coords),
                        sum(c[1] for c in current_coords) / len(current_coords),
                        sum(c[2] for c in current_coords) / len(current_coords),
                    ]
                    poses.append({
                        "pose_rank": len(poses) + 1,
                        "affinity": current_affinity,
                        "rmsd_lb": current_rmsd_lb,
                        "rmsd_ub": current_rmsd_ub,
                        "centroid_x": round(centroid[0], 4),
                        "centroid_y": round(centroid[1], 4),
                        "centroid_z": round(centroid[2], 4),
                    })
                current_coords = []
                current_affinity = None
                current_rmsd_lb = None
                current_rmsd_ub = None

    # Handle last block without ENDMDL
    if current_coords:
        centroid = [
            sum(c[0] for c in current_coords) / len(current_coords),
            sum(c[1] for c in current_coords) / len(current_coords),
            sum(c[2] for c in current_coords) / len(current_coords),
        ]
        poses.append({
            "pose_rank": len(poses) + 1,
            "affinity": current_affinity,
            "rmsd_lb": current_rmsd_lb,
            "rmsd_ub": current_rmsd_ub,
            "centroid_x": round(centroid[0], 4),
            "centroid_y": round(centroid[1], 4),
            "centroid_z": round(centroid[2], 4),
        })

    return poses


def build_pose_table_from_jobs(
    jobs: List[dict],
    runs_dir: Path,
) -> List[dict]:
    """Build pose table by parsing PDBQT outputs from Phase 3 docking jobs.

    For each job, look for output PDBQT file and parse its poses.
    Each pose gets full provenance metadata from the job.
    """
    rows = []

    for job in jobs:
        # Only process dockable jobs
        if job.get("docking_priority") == "skip":
            continue

        # Find output PDBQT
        output_dir = Path(job["output_dir"])
        output_prefix = job["output_prefix"]
        pdbqt_path = output_dir / f"{output_prefix}.pdbqt"

        if not pdbqt_path.exists():
            continue

        poses = parse_pose_blocks(pdbqt_path)
        for pose in poses:
            rows.append({
                # Legacy fields
                "receptor_id": job["receptor_id"],
                "ligand_id": job["ligand_id"],
                "pose_rank": pose["pose_rank"],
                "affinity": pose["affinity"],
                "rmsd_lb": pose["rmsd_lb"],
                "rmsd_ub": pose["rmsd_ub"],
                "centroid_x": pose["centroid_x"],
                "centroid_y": pose["centroid_y"],
                "centroid_z": pose["centroid_z"],
                "raw_pose_file": str(pdbqt_path),
                "pocket_id": "",               # Filled by clustering later
                "contact_residues": "",         # Filled by contacts later
                "n_contact_residues": 0,        # Filled by contacts later
                # Phase 3 provenance
                "candidate_pocket_id": job["pocket_id"],
                "round_id": job["seed_index"],
                "seed": job["seed_index"],
                "docking_mode": job.get("docking_mode", "focused"),
                "job_id": job["job_id"],
            })

    return rows


# ---------------------------------------------------------------------------
# 3.4.2: Backward-compatible pose table generation (no PDBQT yet)
# ---------------------------------------------------------------------------

def build_pose_table_from_job_metadata(
    jobs: List[dict],
) -> List[dict]:
    """Build a placeholder pose table from job metadata when no PDBQT
    output files exist yet (pre-execution / setup mode).

    Creates one row per job with provenance fields populated and
    affinity/pose fields empty. This enables downstream tools to
    validate the schema before server-side execution.
    """
    rows = []
    for job in jobs:
        if job.get("docking_priority") == "skip":
            continue

        rows.append({
            "receptor_id": job["receptor_id"],
            "ligand_id": job["ligand_id"],
            "pose_rank": "",
            "affinity": "",
            "rmsd_lb": "",
            "rmsd_ub": "",
            "centroid_x": "",
            "centroid_y": "",
            "centroid_z": "",
            "raw_pose_file": "",
            "pocket_id": "",
            "contact_residues": "",
            "n_contact_residues": 0,
            "candidate_pocket_id": job["pocket_id"],
            "round_id": job["seed_index"],
            "seed": job["seed_index"],
            "docking_mode": job.get("docking_mode", "focused"),
            "job_id": job["job_id"],
        })

    return rows


# ---------------------------------------------------------------------------
# 3.4.3: Provenance validation
# ---------------------------------------------------------------------------

def validate_pose_provenance(rows: List[dict]) -> List[str]:
    """Validate that all pose rows retain Phase 3 provenance fields."""
    messages = []

    if not rows:
        messages.append("[WARN] Empty pose table")
        return messages

    messages.append(f"[OK] {len(rows)} pose rows")

    # Check provenance fields
    provenance_fields = ["candidate_pocket_id", "round_id", "seed",
                         "docking_mode", "job_id"]
    missing_provenance = []
    for field in provenance_fields:
        if field not in rows[0]:
            missing_provenance.append(field)

    if missing_provenance:
        messages.append(f"[ERROR] Missing provenance fields: {missing_provenance}")
    else:
        messages.append("[OK] All provenance fields present")

    # Check no empty candidate_pocket_id
    empty_pocket = sum(1 for r in rows if not r.get("candidate_pocket_id"))
    if empty_pocket:
        messages.append(f"[WARN] {empty_pocket} rows with empty candidate_pocket_id")
    else:
        messages.append("[OK] All rows have candidate_pocket_id")

    # Check backward compatibility
    legacy_fields = ["receptor_id", "ligand_id", "pose_rank", "affinity",
                     "pocket_id", "contact_residues", "n_contact_residues"]
    missing_legacy = [f for f in legacy_fields if f not in rows[0]]
    if missing_legacy:
        messages.append(f"[WARN] Missing legacy fields: {missing_legacy}")
    else:
        messages.append("[OK] Backward-compatible with legacy pose table schema")

    # Pocket distribution
    by_pocket = defaultdict(int)
    for r in rows:
        by_pocket[r.get("candidate_pocket_id", "")] += 1
    messages.append(f"[OK] Poses by pocket: {dict(sorted(by_pocket.items()))}")

    return messages


# ---------------------------------------------------------------------------
# Provenance note generation
# ---------------------------------------------------------------------------

def build_provenance_note(
    rows: List[dict],
    jobs: List[dict],
    validation: List[str],
) -> List[str]:
    """Generate the Phase 3 pose provenance markdown note."""
    lines = [
        "# Phase 3 Pose Provenance Note",
        "",
        "## Schema Extension",
        "",
        "Phase 3 extends the legacy `vina_pose_table.csv` schema with 5 provenance fields:",
        "",
        "| Field | Type | Description |",
        "|-------|------|-------------|",
        "| candidate_pocket_id | str | Phase 2/3 pocket identity (e.g. 3GT8_raw_PKT01) |",
        "| round_id | int | Docking round (maps to seed_index) |",
        "| seed | int | Seed index for this docking run |",
        "| docking_mode | str | 'focused' (Phase 3) or 'blind' (legacy) |",
        "| job_id | str | Full job ID for traceability |",
        "",
        "These fields are **appended** to the existing 13-column schema.",
        "Legacy tools that read only the first 13 columns remain unaffected.",
        "",
        "## Compatibility",
        "",
        "- `pocket_id` (column 11): Reserved for clustering-assigned pocket ID (legacy flow)",
        "- `candidate_pocket_id` (column 14): Phase 2/3 pre-defined pocket",
        "- Both can coexist — `pocket_id` reflects pose-level clustering,",
        "  `candidate_pocket_id` reflects the intended search target",
        "",
        "## Phase 4 Readiness",
        "",
        "Each pose can be traced to:",
        "- Receptor state (`receptor_id`)",
        "- Target pocket (`candidate_pocket_id`)",
        "- Ligand (`ligand_id`)",
        "- Search round (`round_id`, `seed`)",
        "- Docking mode (`docking_mode`)",
        "",
        "This is sufficient for Phase 4 perturbation scoring which requires",
        "pocket-attributed ligand support evidence.",
        "",
    ]

    # Summary stats
    if rows:
        by_receptor = defaultdict(int)
        by_pocket = defaultdict(int)
        by_ligand = defaultdict(int)
        for r in rows:
            by_receptor[r.get("receptor_id", "")] += 1
            by_pocket[r.get("candidate_pocket_id", "")] += 1
            by_ligand[r.get("ligand_id", "")] += 1

        lines.extend([
            "## Current Pose Table Summary",
            "",
            f"- Total rows: {len(rows)}",
            f"- Receptors: {len(by_receptor)}",
            f"- Pockets: {len(by_pocket)}",
            f"- Ligands: {len(by_ligand)}",
            "",
        ])

        if any(r.get("affinity") for r in rows):
            lines.append("Pose table contains docking results (post-execution).")
        else:
            lines.append(
                "Pose table contains job metadata only (pre-execution). "
                "Re-run after server-side docking to populate affinity values."
            )
        lines.append("")

    # Validation results
    lines.extend([
        "## Validation",
        "",
    ])
    for v in validation:
        lines.append(f"- {v}")
    lines.append("")

    lines.extend([
        "---",
        "",
        "Generated by `egfr_pipeline.phase3.pose_attribution`",
    ])

    return lines


# ---------------------------------------------------------------------------
# Main orchestrator
# ---------------------------------------------------------------------------

def run_pose_attribution(
    output_dir: Path = PHASE3_OUTPUT_DIR,
) -> Tuple[Path, Path]:
    """Run full TG 3.4 pipeline.

    Returns (pose_table_path, provenance_note_path).
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load job table
    job_path = output_dir / "phase3_docking_job_table.csv"
    if not job_path.exists():
        raise FileNotFoundError(
            f"Job table not found: {job_path}\nRun TG 3.1 first."
        )
    jobs = load_csv(job_path)
    print(f"  Loaded {len(jobs)} jobs")

    # Try to parse actual PDBQT outputs first
    runs_dir = output_dir / "runs"
    rows = build_pose_table_from_jobs(jobs, runs_dir)

    if rows:
        print(f"  Parsed {len(rows)} poses from PDBQT output files")
    else:
        # Pre-execution: build metadata-only table
        rows = build_pose_table_from_job_metadata(jobs)
        print(f"  No PDBQT outputs found — built {len(rows)} metadata-only rows")
        print("  (Re-run after server-side docking for full pose table)")

    # Validate
    validation = validate_pose_provenance(rows)
    for v in validation:
        print(f"  {v}")

    # Write pose table
    pose_path = output_dir / "vina_pose_table.csv"
    _write_csv(pose_path, POSE_TABLE_COLUMNS, rows)

    # Write provenance note
    note_lines = build_provenance_note(rows, jobs, validation)
    note_path = output_dir / "phase3_pose_provenance_note.md"
    with open(note_path, "w", encoding="utf-8") as f:
        f.write("\n".join(note_lines))

    return pose_path, note_path


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
        description="Phase 3 TG 3.4: Pose parsing compatibility and pocket attribution"
    )
    parser.add_argument(
        "--output_dir", type=Path, default=PHASE3_OUTPUT_DIR,
        help="Phase 3 output directory",
    )
    args = parser.parse_args()

    print("Phase 3 — Task Group 3.4: Pose Parsing and Pocket Attribution")
    print(f"  Output: {args.output_dir}")
    print()

    pose_path, note_path = run_pose_attribution(args.output_dir)

    print(f"\n{'='*60}")
    print("TG 3.4 COMPLETE")
    print(f"{'='*60}")
    print(f"  Pose table: {pose_path}")
    print(f"  Provenance: {note_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
