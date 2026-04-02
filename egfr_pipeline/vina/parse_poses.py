#!/usr/bin/env python3
import csv
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

from egfr_pipeline.config import load_config
from egfr_pipeline import paths
from egfr_pipeline.parsing_utils import safe_int
from egfr_pipeline.vina.vina_executor import derive_docking_seed


POSE_TABLE_FIELDS = [
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
    "docking_mode",
    "exhaustiveness",
    "requested_n_poses",
    "energy_range",
    "cpu_per_job",
    "vina_seed",
    "pose_source_status",
    "pocket_id",
    "contact_residues",
    "n_contact_residues",
    "contact_distances",
]

POSE_COVERAGE_FIELDS = [
    "receptor_id",
    "ligand_id",
    "raw_pose_file",
    "status",
    "parsed_n_poses",
    "requested_n_poses",
    "coverage_fraction",
    "docking_mode",
    "exhaustiveness",
    "energy_range",
    "cpu_per_job",
    "base_seed",
    "vina_seed",
]


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




def _pose_job_metadata(config: dict) -> Dict[str, object]:
    vina_cfg = config.get("vina", {}) if isinstance(config.get("vina"), dict) else {}
    return {
        "mode": config.get("mode", vina_cfg.get("mode", "blind")),
        "exhaustiveness": vina_cfg.get("exhaustiveness", ""),
        "requested_n_poses": vina_cfg.get("n_poses", ""),
        "energy_range": vina_cfg.get("energy_range", ""),
        "cpu_per_job": vina_cfg.get("cpu", ""),
        "base_seed": vina_cfg.get("seed", ""),
    }


def collect_pose_rows(config: dict) -> Tuple[List[dict], List[dict]]:
    vina_docking_root = paths.wa_phase1_vina_docking(config)
    metadata = _pose_job_metadata(config)
    mode = str(metadata["mode"])
    requested_n_poses = safe_int(metadata["requested_n_poses"])
    base_seed = safe_int(metadata["base_seed"])
    rows: List[dict] = []
    coverage_rows: List[dict] = []

    for receptor in config.get("receptors", []):
        receptor_id = receptor["id"]
        receptor_dir = vina_docking_root / receptor_id
        for ligand in config.get("ligands", []):
            ligand_id = ligand["id"]
            ligand_name = Path(ligand["pdbqt"]).stem.replace("_ligand", "")
            raw_pose_file = receptor_dir / f"{ligand_name}_{mode}.pdbqt"
            poses: List[dict] = []
            status = "missing"
            if raw_pose_file.exists():
                poses = parse_pose_blocks(raw_pose_file)
                status = "parsed" if poses else "empty"

            vina_seed = (
                derive_docking_seed(base_seed, receptor_id, ligand_id)
                if base_seed is not None
                else ""
            )
            parsed_n_poses = len(poses)
            coverage_fraction = ""
            if requested_n_poses and requested_n_poses > 0:
                coverage_fraction = round(parsed_n_poses / requested_n_poses, 4)

            coverage_rows.append(
                {
                    "receptor_id": receptor_id,
                    "ligand_id": ligand_id,
                    "raw_pose_file": str(raw_pose_file),
                    "status": status,
                    "parsed_n_poses": parsed_n_poses,
                    "requested_n_poses": metadata["requested_n_poses"],
                    "coverage_fraction": coverage_fraction,
                    "docking_mode": mode,
                    "exhaustiveness": metadata["exhaustiveness"],
                    "energy_range": metadata["energy_range"],
                    "cpu_per_job": metadata["cpu_per_job"],
                    "base_seed": metadata["base_seed"],
                    "vina_seed": vina_seed,
                }
            )

            for pose in poses:
                rows.append(
                    {
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
                    "docking_mode": mode,
                    "exhaustiveness": metadata["exhaustiveness"],
                    "requested_n_poses": metadata["requested_n_poses"],
                    "energy_range": metadata["energy_range"],
                    "cpu_per_job": metadata["cpu_per_job"],
                    "vina_seed": vina_seed,
                    "pose_source_status": status,
                    "pocket_id": "",
                    "contact_residues": "",
                    "n_contact_residues": 0,
                    "contact_distances": "",
                    }
                )

    return rows, coverage_rows


def iter_pose_rows(config: dict) -> Iterable[dict]:
    rows, _coverage_rows = collect_pose_rows(config)
    yield from rows


def write_pose_table(rows: List[dict], output_csv: Path) -> Path:
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with open(output_csv, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=POSE_TABLE_FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    return output_csv


def write_pose_coverage_table(rows: List[dict], output_csv: Path) -> Path:
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with open(output_csv, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=POSE_COVERAGE_FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    return output_csv


def build_pose_table_from_config(
    config_path: str,
    output_csv: Optional[str] = None,
    coverage_csv: Optional[str] = None,
) -> Path:
    config = load_config(config_path)
    postprocess_root = paths.wa_phase4_vina_postprocess(config)
    target = Path(output_csv) if output_csv else postprocess_root / "vina_pose_table.csv"
    coverage_target = (
        Path(coverage_csv) if coverage_csv else postprocess_root / "vina_postprocess_coverage.csv"
    )
    rows, coverage_rows = collect_pose_rows(config)
    write_pose_coverage_table(coverage_rows, coverage_target)
    return write_pose_table(rows, target)
