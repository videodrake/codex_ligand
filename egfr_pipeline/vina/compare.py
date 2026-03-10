#!/usr/bin/env python3
"""Cross-receptor pocket comparison.

Compares pockets across receptor states using raw metrics:
- Centroid distance (Angstrom)
- Residue overlap (Jaccard similarity)
- Shared ligand evidence

Outputs a comparison table preserving raw metrics for human interpretation.
"""
import csv
import math
from collections import defaultdict
from itertools import combinations
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

from egfr_pipeline.config import load_config, project_root_from_config
from egfr_pipeline.residue_utils import normalize_residue_id, parse_residue_set


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def load_pocket_table(path: Path) -> List[dict]:
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def load_drug_pocket_map(path: Path) -> List[dict]:
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: List[dict], fieldnames: List[str]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return path


# ---------------------------------------------------------------------------
# Comparison metrics
# ---------------------------------------------------------------------------

def centroid_distance(a: dict, b: dict) -> float:
    return math.sqrt(
        (float(a["centroid_x"]) - float(b["centroid_x"])) ** 2
        + (float(a["centroid_y"]) - float(b["centroid_y"])) ** 2
        + (float(a["centroid_z"]) - float(b["centroid_z"])) ** 2
    )


def jaccard_similarity(set_a: Set[str], set_b: Set[str]) -> float:
    if not set_a and not set_b:
        return 0.0
    union = set_a | set_b
    if not union:
        return 0.0
    return len(set_a & set_b) / len(union)


def overlap_coefficient(set_a: Set[str], set_b: Set[str]) -> float:
    """Overlap coefficient = |intersection| / min(|A|, |B|).

    More lenient than Jaccard when pocket sizes differ greatly.
    """
    if not set_a or not set_b:
        return 0.0
    return len(set_a & set_b) / min(len(set_a), len(set_b))


# ---------------------------------------------------------------------------
# Ligand sharing
# ---------------------------------------------------------------------------

def build_pocket_ligands(
    drug_map_rows: List[dict],
) -> Dict[Tuple[str, str], Set[str]]:
    """Map (receptor_id, pocket_id) -> set of ligand_ids that use this pocket.

    A ligand "uses" a pocket if it is the dominant pocket OR an alternative.
    """
    result: Dict[Tuple[str, str], Set[str]] = defaultdict(set)
    for row in drug_map_rows:
        receptor_id = row["receptor_id"]
        ligand_id = row["ligand_id"]
        dominant = row["dominant_pocket_id"]
        result[(receptor_id, dominant)].add(ligand_id)
        alternatives = row.get("alternative_pockets", "")
        if alternatives:
            for pocket_id in alternatives.split(";"):
                pocket_id = pocket_id.strip()
                if pocket_id:
                    result[(receptor_id, pocket_id)].add(ligand_id)
    return dict(result)


# ---------------------------------------------------------------------------
# Core comparison
# ---------------------------------------------------------------------------

COMPARISON_FIELDS = [
    "receptor_a",
    "pocket_a",
    "receptor_b",
    "pocket_b",
    "centroid_dist",
    "residue_jaccard",
    "residue_overlap_coeff",
    "shared_residues",
    "n_shared_residues",
    "residues_only_a",
    "residues_only_b",
    "n_residues_a",
    "n_residues_b",
    "shared_ligands",
    "n_shared_ligands",
    "n_ligands_a",
    "n_ligands_b",
    "affinity_a",
    "affinity_b",
    "n_pose_a",
    "n_pose_b",
    "same_patch_candidate",
]


def compare_all_pockets(
    pocket_rows: List[dict],
    drug_map_rows: List[dict],
    centroid_cutoff: float = 15.0,
    keep_chain: bool = False,
) -> List[dict]:
    """Compare every pocket pair across different receptor states.

    Args:
        pocket_rows: rows from vina_pocket_table.csv
        drug_map_rows: rows from vina_drug_pocket_map.csv
        centroid_cutoff: only emit pairs within this centroid distance (A).
            Set to 0 to emit all pairs.
        keep_chain: preserve chain prefix in residue identifiers.
    """
    # Group pockets by receptor
    by_receptor: Dict[str, List[dict]] = defaultdict(list)
    for row in pocket_rows:
        by_receptor[row["receptor_id"]].append(row)

    pocket_ligands = build_pocket_ligands(drug_map_rows)
    receptor_ids = sorted(by_receptor.keys())

    results: List[dict] = []

    for rec_a, rec_b in combinations(receptor_ids, 2):
        for pocket_a in by_receptor[rec_a]:
            for pocket_b in by_receptor[rec_b]:
                dist = centroid_distance(pocket_a, pocket_b)
                if centroid_cutoff > 0 and dist > centroid_cutoff:
                    continue

                res_a = parse_residue_set(pocket_a.get("union_contact_residues", ""), keep_chain=keep_chain)
                res_b = parse_residue_set(pocket_b.get("union_contact_residues", ""), keep_chain=keep_chain)
                shared = sorted(res_a & res_b)
                only_a = sorted(res_a - res_b)
                only_b = sorted(res_b - res_a)

                lig_a = pocket_ligands.get((rec_a, pocket_a["pocket_id"]), set())
                lig_b = pocket_ligands.get((rec_b, pocket_b["pocket_id"]), set())
                shared_lig = sorted(lig_a & lig_b)

                j = jaccard_similarity(res_a, res_b)
                oc = overlap_coefficient(res_a, res_b)

                # Auxiliary same-patch candidate flag:
                # centroid < 8 A AND (residue jaccard >= 0.3 OR overlap_coeff >= 0.5)
                # This is a soft label for human review, not a definitive classification.
                is_candidate = (
                    dist < 8.0
                    and (j >= 0.3 or oc >= 0.5)
                )

                results.append({
                    "receptor_a": rec_a,
                    "pocket_a": pocket_a["pocket_id"],
                    "receptor_b": rec_b,
                    "pocket_b": pocket_b["pocket_id"],
                    "centroid_dist": round(dist, 2),
                    "residue_jaccard": round(j, 4),
                    "residue_overlap_coeff": round(oc, 4),
                    "shared_residues": ";".join(shared),
                    "n_shared_residues": len(shared),
                    "residues_only_a": ";".join(only_a),
                    "residues_only_b": ";".join(only_b),
                    "n_residues_a": len(res_a),
                    "n_residues_b": len(res_b),
                    "shared_ligands": ";".join(shared_lig),
                    "n_shared_ligands": len(shared_lig),
                    "n_ligands_a": len(lig_a),
                    "n_ligands_b": len(lig_b),
                    "affinity_a": pocket_a.get("best_affinity", ""),
                    "affinity_b": pocket_b.get("best_affinity", ""),
                    "n_pose_a": pocket_a.get("n_pose", ""),
                    "n_pose_b": pocket_b.get("n_pose", ""),
                    "same_patch_candidate": is_candidate,
                })

    # Sort by centroid distance (closest pairs first)
    results.sort(key=lambda r: r["centroid_dist"])
    return results


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def compare_from_config(
    config_path: str,
    pocket_table_path: Optional[str] = None,
    drug_map_path: Optional[str] = None,
    centroid_cutoff: float = 15.0,
    output_path: Optional[str] = None,
) -> Path:
    config = load_config(config_path)
    project_root = project_root_from_config(config)

    pocket_csv = Path(pocket_table_path) if pocket_table_path else project_root / "vina_pocket_table.csv"
    drug_csv = Path(drug_map_path) if drug_map_path else project_root / "vina_drug_pocket_map.csv"
    out_csv = Path(output_path) if output_path else project_root / "vina_pocket_comparison.csv"

    pocket_rows = load_pocket_table(pocket_csv)
    drug_map_rows = load_drug_pocket_map(drug_csv)

    keep_chain = config.get("postprocess", {}).get("keep_chain", False)
    comparison = compare_all_pockets(pocket_rows, drug_map_rows, centroid_cutoff, keep_chain=keep_chain)
    return write_csv(out_csv, comparison, COMPARISON_FIELDS)
