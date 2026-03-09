#!/usr/bin/env python3
import argparse
import csv
import json
import math
from pathlib import Path
from typing import Dict, List, Optional

from parse_vina_results import load_config, project_root_from_config


def load_pose_table(path: Path) -> List[dict]:
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_pose_table(path: Path, rows: List[dict]) -> Path:
    if not rows:
        return path
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    return path


def euclidean_distance(a: List[float], b: List[float]) -> float:
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


def row_centroid(row: dict) -> List[float]:
    return [float(row["centroid_x"]), float(row["centroid_y"]), float(row["centroid_z"])]


def sort_pose_rows(rows: List[dict]) -> List[dict]:
    def sort_key(row: dict):
        affinity = float(row["affinity"]) if row["affinity"] not in ("", None) else float("inf")
        return (
            row["receptor_id"],
            affinity,
            row["ligand_id"],
            int(row["pose_rank"]),
            row["raw_pose_file"],
        )
    return sorted(rows, key=sort_key)


def assign_pockets(rows: List[dict], cutoff: float) -> List[dict]:
    sorted_rows = sort_pose_rows(rows)
    pocket_state: Dict[str, List[dict]] = {}

    for row in sorted_rows:
        receptor_id = row["receptor_id"]
        centroid = row_centroid(row)
        receptor_pockets = pocket_state.setdefault(receptor_id, [])

        selected = None
        min_distance = None
        for pocket in receptor_pockets:
            dist = euclidean_distance(centroid, pocket["center"])
            if dist <= cutoff and (min_distance is None or dist < min_distance):
                selected = pocket
                min_distance = dist

        if selected is None:
            pocket_index = len(receptor_pockets) + 1
            selected = {
                "id": f"P{pocket_index:03d}",
                "center": centroid[:],
                "count": 0,
            }
            receptor_pockets.append(selected)

        selected["count"] += 1
        count = selected["count"]
        selected["center"] = [
            ((selected["center"][axis] * (count - 1)) + centroid[axis]) / count
            for axis in range(3)
        ]
        row["pocket_id"] = selected["id"]

    return sorted_rows


def cluster_pose_table(config_path: str, pose_table_path: Optional[str] = None, cutoff: float = 4.0) -> Path:
    config = load_config(config_path)
    project_root = project_root_from_config(config)
    target = Path(pose_table_path) if pose_table_path else project_root / "vina_pose_table.csv"
    rows = load_pose_table(target)
    clustered_rows = assign_pockets(rows, cutoff)
    return write_pose_table(target, clustered_rows)


def parse_args():
    parser = argparse.ArgumentParser(description="Assign receptor-level pocket ids from vina_pose_table.csv.")
    parser.add_argument("--config", required=True, help="Project YAML/JSON config path")
    parser.add_argument("--pose-table", default=None, help="Pose table CSV path (default: <project_root>/vina_pose_table.csv)")
    parser.add_argument("--cutoff", type=float, default=4.0, help="Pocket centroid clustering cutoff in Angstrom (default: 4.0)")
    return parser.parse_args()


def main():
    args = parse_args()
    output_csv = cluster_pose_table(args.config, args.pose_table, args.cutoff)
    print(f"[OK] Updated pose table with pocket assignments: {output_csv}")


if __name__ == "__main__":
    main()
