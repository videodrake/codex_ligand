#!/usr/bin/env python3
import csv
import statistics
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from egfr_pipeline.config import load_config, project_root_from_config


def load_pose_table(path: Path) -> List[dict]:
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def split_contact_residues(value: str) -> List[str]:
    if not value:
        return []
    return [item for item in value.split(";") if item]


def summarize_pose_rows(rows: List[dict]) -> Tuple[List[dict], List[dict]]:
    pocket_groups: Dict[Tuple[str, str], List[dict]] = defaultdict(list)
    ligand_groups: Dict[Tuple[str, str], List[dict]] = defaultdict(list)

    for row in rows:
        if not row.get("pocket_id"):
            continue
        pocket_groups[(row["receptor_id"], row["pocket_id"])].append(row)
        ligand_groups[(row["receptor_id"], row["ligand_id"])].append(row)

    pocket_rows: List[dict] = []
    for (receptor_id, pocket_id), group in sorted(pocket_groups.items()):
        centroids_x = [float(item["centroid_x"]) for item in group]
        centroids_y = [float(item["centroid_y"]) for item in group]
        centroids_z = [float(item["centroid_z"]) for item in group]
        affinities = [float(item["affinity"]) for item in group if item["affinity"] not in ("", None)]
        ligands = sorted({item["ligand_id"] for item in group})
        residue_counter: Counter = Counter()
        for item in group:
            residue_counter.update(split_contact_residues(item.get("contact_residues", "")))
        union_contact_residues = sorted(residue_counter)
        top_residues = [residue for residue, _ in sorted(residue_counter.items(), key=lambda item: (-item[1], item[0]))[:5]]
        pocket_rows.append({
            "receptor_id": receptor_id,
            "pocket_id": pocket_id,
            "centroid_x": round(sum(centroids_x) / len(centroids_x), 4),
            "centroid_y": round(sum(centroids_y) / len(centroids_y), 4),
            "centroid_z": round(sum(centroids_z) / len(centroids_z), 4),
            "n_pose": len(group),
            "n_ligand": len(ligands),
            "best_affinity": round(min(affinities), 4) if affinities else "",
            "mean_affinity": round(sum(affinities) / len(affinities), 4) if affinities else "",
            "union_contact_residues": ";".join(union_contact_residues),
            "top_residues": ";".join(top_residues),
        })

    drug_map_rows: List[dict] = []
    for (receptor_id, ligand_id), group in sorted(ligand_groups.items()):
        pocket_counter: Counter = Counter(item["pocket_id"] for item in group)
        best_by_pocket: Dict[str, float] = {}
        best_pose_rank_by_pocket: Dict[str, int] = {}
        residues_by_pocket: Dict[str, str] = {}
        for item in group:
            pocket_id = item["pocket_id"]
            affinity = float(item["affinity"]) if item["affinity"] not in ("", None) else float("inf")
            pose_rank = int(item["pose_rank"])
            if pocket_id not in best_by_pocket or affinity < best_by_pocket[pocket_id]:
                best_by_pocket[pocket_id] = affinity
                best_pose_rank_by_pocket[pocket_id] = pose_rank
                residues_by_pocket[pocket_id] = item.get("contact_residues", "")

        dominant_pocket_id = sorted(
            pocket_counter,
            key=lambda pocket_id: (-pocket_counter[pocket_id], best_by_pocket.get(pocket_id, float("inf")), pocket_id),
        )[0]
        alternative_pockets = sorted(pocket_id for pocket_id in pocket_counter if pocket_id != dominant_pocket_id)
        dominant_pose_count = pocket_counter[dominant_pocket_id]
        drug_map_rows.append({
            "receptor_id": receptor_id,
            "ligand_id": ligand_id,
            "dominant_pocket_id": dominant_pocket_id,
            "dominant_pocket_pose_count": dominant_pose_count,
            "dominant_pocket_fraction": round(dominant_pose_count / len(group), 4),
            "best_affinity": round(best_by_pocket[dominant_pocket_id], 4),
            "best_pose_rank": best_pose_rank_by_pocket[dominant_pocket_id],
            "top_pose_residues": residues_by_pocket.get(dominant_pocket_id, ""),
            "alternative_pockets": ";".join(alternative_pockets),
            "is_multimodal_binding": len(pocket_counter) > 1,
        })

    return pocket_rows, drug_map_rows


def write_csv(path: Path, rows: List[dict], fieldnames: List[str]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return path


def summarize_from_config(config_path: str, pose_table_path: Optional[str] = None) -> Tuple[Path, Path]:
    config = load_config(config_path)
    project_root = project_root_from_config(config)
    pose_table = Path(pose_table_path) if pose_table_path else project_root / "vina_pose_table.csv"
    rows = load_pose_table(pose_table)

    # Cross-receptor comparison is not implemented yet; residue numbering consistency
    # should be verified on the real server before Task Group 5 uses these summaries.
    pocket_rows, drug_map_rows = summarize_pose_rows(rows)
    pocket_csv = write_csv(
        project_root / "vina_pocket_table.csv",
        pocket_rows,
        [
            "receptor_id",
            "pocket_id",
            "centroid_x",
            "centroid_y",
            "centroid_z",
            "n_pose",
            "n_ligand",
            "best_affinity",
            "mean_affinity",
            "union_contact_residues",
            "top_residues",
        ],
    )
    drug_csv = write_csv(
        project_root / "vina_drug_pocket_map.csv",
        drug_map_rows,
        [
            "receptor_id",
            "ligand_id",
            "dominant_pocket_id",
            "dominant_pocket_pose_count",
            "dominant_pocket_fraction",
            "best_affinity",
            "best_pose_rank",
            "top_pose_residues",
            "alternative_pockets",
            "is_multimodal_binding",
        ],
    )
    return pocket_csv, drug_csv
