#!/usr/bin/env python3
"""Phase 3 Task Group 3.1: Pocket-Guided Docking Job Construction.

Builds receptor- and pocket-specific docking jobs from the Phase 3
normalized pocket reference and a ligand inventory.

Tasks:
  3.1.1 - Receptor-local job generation (one job set per receptor state)
  3.1.2 - Pocket-local box generation (per-pocket box metadata)
  3.1.3 - Ligand dispatch matrix (receptor x pocket x ligand x seed)

Inputs:
  - output/phase3_docking/phase3_candidate_reference_normalized.csv
  - input/ligands/*.sdf  (auto-detected)

Outputs:
  - output/phase3_docking/phase3_docking_job_table.csv
  - output/phase3_docking/phase3_job_box_table.csv

Usage:
    python -m egfr_pipeline.phase3.job_construction [--output_dir ...]
"""

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

from egfr_pipeline.csv_utils import load_csv

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

PHASE3_OUTPUT_DIR = PROJECT_ROOT / "output" / "phase3_docking"
LIGAND_INPUT_DIR = PROJECT_ROOT / "input" / "ligands"
RECEPTOR_INPUT_DIR = PROJECT_ROOT / "input" / "receptors"

# Supported ligand input formats (SDF preferred, PDBQT if pre-converted)
LIGAND_EXTENSIONS = (".pdbqt", ".sdf", ".mol2")  # Prefer PDBQT (Vina-ready)

# Output schemas
JOB_TABLE_COLUMNS = [
    "job_id",
    "receptor_id",
    "pocket_id",
    "ligand_id",
    "seed_index",
    # Docking geometry
    "center_x",
    "center_y",
    "center_z",
    "box_size_x",
    "box_size_y",
    "box_size_z",
    # Docking parameters
    "exhaustiveness",
    "n_poses",
    "docking_mode",
    "cpu_per_job",
    "engine",
    "gpu_required",
    # Paths
    "receptor_pdb",
    "ligand_file",
    "output_dir",
    "output_prefix",
    # Provenance
    "docking_priority",
    "relationship_class",
    "overall_druggability_tier",
]

BOX_TABLE_COLUMNS = [
    "pocket_id",
    "receptor_id",
    "center_x",
    "center_y",
    "center_z",
    "box_size_x",
    "box_size_y",
    "box_size_z",
    "box_source",
    "docking_priority",
    "relationship_class",
    "dockable",
]


# ---------------------------------------------------------------------------
# 3.1.1: Receptor-local job generation
# ---------------------------------------------------------------------------

def discover_ligands(ligand_dir: Path) -> List[Dict[str, str]]:
    """Auto-detect ligand files from the input directory.

    Returns list of dicts with 'ligand_id' and 'ligand_file' keys,
    sorted by ligand_id for determinism.
    """
    ligands = []
    if not ligand_dir.exists():
        return ligands

    seen_ids = set()
    for ext in LIGAND_EXTENSIONS:
        for path in sorted(ligand_dir.glob(f"*{ext}")):
            # Extract ligand ID: strip _ligand suffix and extension
            stem = path.stem
            if stem.endswith("_ligand"):
                lig_id = stem[: -len("_ligand")]
            else:
                lig_id = stem
            # Prefer first-found (SDF > PDBQT > MOL2 by extension order)
            if lig_id not in seen_ids:
                seen_ids.add(lig_id)
                ligands.append({
                    "ligand_id": lig_id,
                    "ligand_file": str(path),
                })

    return sorted(ligands, key=lambda x: x["ligand_id"])


def resolve_receptor_path(receptor_id: str, receptor_dir: Path) -> str:
    """Find the receptor file for a given receptor_id. Prefer PDBQT for Vina."""
    # Prefer PDBQT (Vina-ready)
    pdbqt_path = receptor_dir / f"{receptor_id}_receptor.pdbqt"
    if pdbqt_path.exists():
        return str(pdbqt_path)
    pdbqt_plain = receptor_dir / f"{receptor_id}.pdbqt"
    if pdbqt_plain.exists():
        return str(pdbqt_plain)
    # Fallback: PDB (will need conversion)
    pdb_path = receptor_dir / f"{receptor_id}.pdb"
    return str(pdb_path)


def build_job_id(receptor_id: str, pocket_id: str, ligand_id: str,
                 seed_index: int) -> str:
    """Generate a stable, deterministic job ID.

    Format: {receptor_id}__{pocket_id}__{ligand_id}__seed{N}
    Double-underscore separators for unambiguous parsing.
    """
    return f"{receptor_id}__{pocket_id}__{ligand_id}__seed{seed_index}"


def build_output_prefix(receptor_id: str, pocket_id: str, ligand_id: str,
                        seed_index: int) -> str:
    """Output file prefix for docking results."""
    return f"{pocket_id}_{ligand_id}_seed{seed_index}"


def generate_jobs(
    normalized_pockets: List[dict],
    ligands: List[Dict[str, str]],
    receptor_dir: Path,
    output_root: Path,
    *,
    cpu_per_job: int = 4,
    engine: str = "cpu-vina",
    gpu_required: bool = False,
) -> List[dict]:
    """Generate the full job table: receptor x pocket x ligand x seed.

    Only dockable pockets are included. Receptor-local separation is enforced
    by construction (pockets carry receptor_id, never mixed).
    """
    jobs = []

    for pocket in normalized_pockets:
        # Skip non-dockable pockets
        if pocket.get("dockable", "False") != "True":
            continue

        receptor_id = pocket["receptor_id"]
        pocket_id = pocket["pocket_id"]
        n_seeds = int(pocket.get("recommended_seeds", 1))
        exhaustiveness = int(pocket.get("recommended_exhaustiveness", 8))
        n_poses = int(pocket.get("recommended_n_poses", 5))
        receptor_pdb = resolve_receptor_path(receptor_id, receptor_dir)

        for ligand in ligands:
            ligand_id = ligand["ligand_id"]
            ligand_file = ligand["ligand_file"]

            for seed_idx in range(1, n_seeds + 1):
                job_id = build_job_id(receptor_id, pocket_id, ligand_id, seed_idx)
                out_prefix = build_output_prefix(
                    receptor_id, pocket_id, ligand_id, seed_idx
                )
                # Receptor-local output directory
                job_output_dir = str(
                    output_root / receptor_id / pocket_id / ligand_id
                )

                jobs.append({
                    "job_id": job_id,
                    "receptor_id": receptor_id,
                    "pocket_id": pocket_id,
                    "ligand_id": ligand_id,
                    "seed_index": seed_idx,
                    "center_x": pocket.get("centroid_x", ""),
                    "center_y": pocket.get("centroid_y", ""),
                    "center_z": pocket.get("centroid_z", ""),
                    "box_size_x": pocket.get("box_size_x", ""),
                    "box_size_y": pocket.get("box_size_y", ""),
                    "box_size_z": pocket.get("box_size_z", ""),
                    "exhaustiveness": exhaustiveness,
                    "n_poses": n_poses,
                    "docking_mode": "focused",
                    "cpu_per_job": int(cpu_per_job),
                    "engine": engine,
                    "gpu_required": bool(gpu_required),
                    "receptor_pdb": receptor_pdb,
                    "ligand_file": ligand_file,
                    "output_dir": job_output_dir,
                    "output_prefix": out_prefix,
                    "docking_priority": pocket.get("docking_priority", ""),
                    "relationship_class": pocket.get("relationship_class", ""),
                    "overall_druggability_tier": pocket.get(
                        "overall_druggability_tier", ""
                    ),
                })

    return jobs


# ---------------------------------------------------------------------------
# 3.1.2: Pocket-local box generation
# ---------------------------------------------------------------------------

def build_box_table(normalized_pockets: List[dict]) -> List[dict]:
    """Extract per-pocket box metadata.

    Preserves all pockets (including skip) for reference.
    """
    boxes = []
    for pocket in normalized_pockets:
        boxes.append({
            "pocket_id": pocket["pocket_id"],
            "receptor_id": pocket["receptor_id"],
            "center_x": pocket.get("centroid_x", ""),
            "center_y": pocket.get("centroid_y", ""),
            "center_z": pocket.get("centroid_z", ""),
            "box_size_x": pocket.get("box_size_x", ""),
            "box_size_y": pocket.get("box_size_y", ""),
            "box_size_z": pocket.get("box_size_z", ""),
            "box_source": "phase2_fpocket",
            "docking_priority": pocket.get("docking_priority", ""),
            "relationship_class": pocket.get("relationship_class", ""),
            "dockable": pocket.get("dockable", "False"),
        })
    return boxes


# ---------------------------------------------------------------------------
# 3.1.3: Ligand dispatch matrix summary
# ---------------------------------------------------------------------------

def summarize_dispatch(jobs: List[dict]) -> Dict[str, dict]:
    """Summarize job dispatch by receptor state.

    Returns {receptor_id: {n_pockets, n_ligands, n_seeds_total, n_jobs}}.
    """
    by_receptor: Dict[str, dict] = {}
    receptor_jobs: Dict[str, List[dict]] = defaultdict(list)

    for j in jobs:
        receptor_jobs[j["receptor_id"]].append(j)

    for receptor_id, rjobs in sorted(receptor_jobs.items()):
        pockets = set(j["pocket_id"] for j in rjobs)
        ligands = set(j["ligand_id"] for j in rjobs)
        by_receptor[receptor_id] = {
            "n_pockets": len(pockets),
            "n_ligands": len(ligands),
            "n_jobs": len(rjobs),
        }

    return by_receptor


# ---------------------------------------------------------------------------
# Main orchestrator
# ---------------------------------------------------------------------------

def run_job_construction(
    output_dir: Path = PHASE3_OUTPUT_DIR,
    ligand_dir: Path = LIGAND_INPUT_DIR,
    receptor_dir: Path = RECEPTOR_INPUT_DIR,
) -> Tuple[Path, Path]:
    """Run full TG 3.1 pipeline.

    Returns (job_table_path, box_table_path).
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load normalized pocket reference from TG 3.0
    norm_path = output_dir / "phase3_candidate_reference_normalized.csv"
    if not norm_path.exists():
        raise FileNotFoundError(
            f"Normalized pocket reference not found: {norm_path}\n"
            "Run TG 3.0 (pocket_reference_ingestion) first."
        )

    normalized = load_csv(norm_path)
    dockable = [n for n in normalized if n.get("dockable") == "True"]
    print(f"  Loaded {len(normalized)} pockets ({len(dockable)} dockable)")

    # Discover ligands
    ligands = discover_ligands(ligand_dir)
    if not ligands:
        raise FileNotFoundError(
            f"No ligand files found in {ligand_dir}\n"
            f"Expected formats: {LIGAND_EXTENSIONS}"
        )
    print(f"  Discovered {len(ligands)} ligand(s): "
          f"{', '.join(l['ligand_id'] for l in ligands)}")

    # Task 3.1.1: Generate job table
    docking_output_root = output_dir / "runs"
    jobs = generate_jobs(normalized, ligands, receptor_dir, docking_output_root)
    print(f"  Generated {len(jobs)} docking jobs")

    # Task 3.1.2: Build box table
    box_table = build_box_table(normalized)
    print(f"  Built box table: {len(box_table)} pocket entries")

    # Task 3.1.3: Dispatch summary
    dispatch = summarize_dispatch(jobs)
    for receptor_id, stats in dispatch.items():
        print(f"    {receptor_id}: {stats['n_pockets']} pockets x "
              f"{stats['n_ligands']} ligands = {stats['n_jobs']} jobs")

    # Validate: no cross-receptor pocket mixing
    for j in jobs:
        assert j["pocket_id"].startswith(j["receptor_id"]), (
            f"Cross-receptor mixing detected: job {j['job_id']} "
            f"pocket {j['pocket_id']} != receptor {j['receptor_id']}"
        )

    # Write outputs
    job_path = output_dir / "phase3_docking_job_table.csv"
    _write_csv(job_path, JOB_TABLE_COLUMNS, jobs)

    box_path = output_dir / "phase3_job_box_table.csv"
    _write_csv(box_path, BOX_TABLE_COLUMNS, box_table)

    return job_path, box_path


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
        description="Phase 3 TG 3.1: Pocket-guided docking job construction"
    )
    parser.add_argument(
        "--output_dir", type=Path, default=PHASE3_OUTPUT_DIR,
        help="Phase 3 output directory",
    )
    parser.add_argument(
        "--ligand_dir", type=Path, default=LIGAND_INPUT_DIR,
        help="Directory containing ligand files (SDF/PDBQT)",
    )
    parser.add_argument(
        "--receptor_dir", type=Path, default=RECEPTOR_INPUT_DIR,
        help="Directory containing receptor PDB files",
    )
    args = parser.parse_args()

    print("Phase 3 — Task Group 3.1: Pocket-Guided Docking Job Construction")
    print(f"  Phase 3 output: {args.output_dir}")
    print(f"  Ligand input:   {args.ligand_dir}")
    print(f"  Receptor input: {args.receptor_dir}")
    print()

    job_path, box_path = run_job_construction(
        args.output_dir, args.ligand_dir, args.receptor_dir
    )

    print(f"\n{'='*60}")
    print("TG 3.1 COMPLETE")
    print(f"{'='*60}")
    print(f"  Job table: {job_path}")
    print(f"  Box table: {box_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
