#!/usr/bin/env python3
import argparse
import csv
import json
import math
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from parse_vina_results import load_config, parse_pose_blocks, project_root_from_config


def parse_receptor_atoms(pdb_path: Path) -> List[dict]:
    atoms: List[dict] = []
    with open(pdb_path, encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            atom_name = line[12:16].strip()
            element = line[76:78].strip() or atom_name[:1]
            if element.upper() == "H" or atom_name.upper().startswith("H"):
                continue
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue
            chain_id = line[21].strip() or "?"
            resname = line[17:20].strip()
            resid = line[22:26].strip()
            residue_id = f"{chain_id}:{resname}{resid}"
            atoms.append({
                "residue_id": residue_id,
                "coord": (x, y, z),
            })
    return atoms


def parse_pose_atoms(pose_lines: List[str]) -> List[Tuple[float, float, float]]:
    atoms: List[Tuple[float, float, float]] = []
    for line in pose_lines:
        if not line.startswith(("ATOM", "HETATM")):
            continue
        atom_name = line[12:16].strip()
        element = line[76:78].strip() or atom_name[:1]
        if element.upper() == "H" or atom_name.upper().startswith("H"):
            continue
        try:
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
        except ValueError:
            continue
        atoms.append((x, y, z))
    return atoms


def compute_contacts(
    receptor_atoms: List[dict],
    ligand_atoms: List[Tuple[float, float, float]],
    cutoff: float,
) -> List[str]:
    cutoff_sq = cutoff * cutoff
    min_dist_by_residue: Dict[str, float] = {}
    for receptor_atom in receptor_atoms:
        rx, ry, rz = receptor_atom["coord"]
        for lx, ly, lz in ligand_atoms:
            dist_sq = (rx - lx) ** 2 + (ry - ly) ** 2 + (rz - lz) ** 2
            if dist_sq > cutoff_sq:
                continue
            residue_id = receptor_atom["residue_id"]
            dist = math.sqrt(dist_sq)
            current = min_dist_by_residue.get(residue_id)
            if current is None or dist < current:
                min_dist_by_residue[residue_id] = dist
    return sorted(min_dist_by_residue)


def load_pose_table(path: Path) -> List[dict]:
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_pose_table(path: Path, rows: List[dict]) -> None:
    if not rows:
        return
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def build_receptor_pdb_map(config: dict) -> Dict[str, Path]:
    return {item["id"]: Path(item["pdb"]) for item in config.get("receptors", []) if item.get("pdb")}


def enrich_pose_table_with_contacts(config_path: str, pose_table_path: Optional[str] = None, cutoff: float = 4.0) -> Path:
    config = load_config(config_path)
    project_root = project_root_from_config(config)
    pose_table = Path(pose_table_path) if pose_table_path else project_root / "vina_pose_table.csv"
    rows = load_pose_table(pose_table)
    receptor_pdb_map = build_receptor_pdb_map(config)
    receptor_cache: Dict[str, List[dict]] = {}
    pose_cache: Dict[str, List[dict]] = {}

    for row in rows:
        receptor_id = row["receptor_id"]
        receptor_pdb = receptor_pdb_map.get(receptor_id)
        if receptor_pdb is None or not receptor_pdb.exists():
            raise FileNotFoundError(f"Missing receptor PDB for {receptor_id}")

        if receptor_id not in receptor_cache:
            receptor_cache[receptor_id] = parse_receptor_atoms(receptor_pdb)

        raw_pose_file = row["raw_pose_file"]
        if raw_pose_file not in pose_cache:
            pose_cache[raw_pose_file] = parse_pose_blocks(Path(raw_pose_file))

        pose_rank = int(row["pose_rank"])
        pose = pose_cache[raw_pose_file][pose_rank - 1]
        contacts = compute_contacts(receptor_cache[receptor_id], parse_pose_atoms(pose["lines"]), cutoff)
        row["contact_residues"] = ";".join(contacts)
        row["n_contact_residues"] = len(contacts)

    write_pose_table(pose_table, rows)
    return pose_table


def parse_args():
    parser = argparse.ArgumentParser(description="Extract receptor contact residues for Vina poses.")
    parser.add_argument("--config", required=True, help="Project YAML/JSON config path")
    parser.add_argument("--pose-table", default=None, help="Pose table CSV path (default: <project_root>/vina_pose_table.csv)")
    parser.add_argument("--cutoff", type=float, default=4.0, help="Heavy-atom contact cutoff in Angstrom (default: 4.0)")
    return parser.parse_args()


def main():
    args = parse_args()
    output_csv = enrich_pose_table_with_contacts(args.config, args.pose_table, args.cutoff)
    print(f"[OK] Updated pose table with contacts: {output_csv}")


if __name__ == "__main__":
    main()
