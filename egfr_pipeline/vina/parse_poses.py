#!/usr/bin/env python3
import csv
from pathlib import Path
from typing import Dict, Iterable, List, Optional

from egfr_pipeline.config import load_config, project_root_from_config


def parse_pose_blocks(pdbqt_path: Path) -> List[dict]:
    poses: List[dict] = []
    current_lines: List[str] = []
    current_coords: List[List[float]] = []
    current_affinity: Optional[float] = None
    current_rmsd_lb: Optional[float] = None
    current_rmsd_ub: Optional[float] = None

    with open(pdbqt_path, encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            current_lines.append(line)
            if line.startswith("REMARK VINA RESULT:"):
                parts = line.split()
                if len(parts) >= 6:
                    current_affinity = float(parts[3])
                    current_rmsd_lb = float(parts[4])
                    current_rmsd_ub = float(parts[5])
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
                        sum(coord[0] for coord in current_coords) / len(current_coords),
                        sum(coord[1] for coord in current_coords) / len(current_coords),
                        sum(coord[2] for coord in current_coords) / len(current_coords),
                    ]
                    poses.append({
                        "pose_rank": len(poses) + 1,
                        "affinity": current_affinity,
                        "rmsd_lb": current_rmsd_lb,
                        "rmsd_ub": current_rmsd_ub,
                        "centroid_x": centroid[0],
                        "centroid_y": centroid[1],
                        "centroid_z": centroid[2],
                        "lines": current_lines[:],
                    })
                current_lines = []
                current_coords = []
                current_affinity = None
                current_rmsd_lb = None
                current_rmsd_ub = None

    if current_lines and current_coords:
        centroid = [
            sum(coord[0] for coord in current_coords) / len(current_coords),
            sum(coord[1] for coord in current_coords) / len(current_coords),
            sum(coord[2] for coord in current_coords) / len(current_coords),
        ]
        poses.append({
            "pose_rank": len(poses) + 1,
            "affinity": current_affinity,
            "rmsd_lb": current_rmsd_lb,
            "rmsd_ub": current_rmsd_ub,
            "centroid_x": centroid[0],
            "centroid_y": centroid[1],
            "centroid_z": centroid[2],
            "lines": current_lines[:],
        })

    return poses


def iter_pose_rows(config: dict) -> Iterable[dict]:
    project_root = project_root_from_config(config)
    mode = config.get("mode", config.get("vina", {}).get("mode", "blind"))

    for receptor in config.get("receptors", []):
        receptor_id = receptor["id"]
        receptor_dir = project_root / receptor_id
        for ligand in config.get("ligands", []):
            ligand_id = ligand["id"]
            ligand_name = Path(ligand["pdbqt"]).stem.replace("_ligand", "")
            raw_pose_file = receptor_dir / f"{ligand_name}_{mode}.pdbqt"
            if not raw_pose_file.exists():
                continue
            for pose in parse_pose_blocks(raw_pose_file):
                yield {
                    "receptor_id": receptor_id,
                    "ligand_id": ligand_id,
                    "pose_rank": pose["pose_rank"],
                    "affinity": pose["affinity"],
                    "rmsd_lb": pose["rmsd_lb"],
                    "rmsd_ub": pose["rmsd_ub"],
                    "centroid_x": round(pose["centroid_x"], 4),
                    "centroid_y": round(pose["centroid_y"], 4),
                    "centroid_z": round(pose["centroid_z"], 4),
                    "raw_pose_file": str(raw_pose_file),
                    "pocket_id": "",
                    "contact_residues": "",
                    "n_contact_residues": 0,
                }


def write_pose_table(rows: List[dict], output_csv: Path) -> Path:
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
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
        "pocket_id",
        "contact_residues",
        "n_contact_residues",
    ]
    with open(output_csv, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return output_csv


def build_pose_table_from_config(config_path: str, output_csv: Optional[str] = None) -> Path:
    config = load_config(config_path)
    project_root = project_root_from_config(config)
    target = Path(output_csv) if output_csv else project_root / "vina_pose_table.csv"
    rows = list(iter_pose_rows(config))
    return write_pose_table(rows, target)
