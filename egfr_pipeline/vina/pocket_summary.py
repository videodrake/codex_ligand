#!/usr/bin/env python3
import csv
import math
import statistics
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from egfr_pipeline.config import load_config
from egfr_pipeline import paths
from egfr_pipeline.region_definitions import (
    ATP_SITE_RESIDUES,
    KO_SHEET_8_9,
    KO_SHEET_10_11,
    KO_SHEET_12,
)
from egfr_pipeline.residue_utils import extract_resnum


def load_pose_table(path: Path) -> List[dict]:
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def split_contact_residues(value: str) -> List[str]:
    if not value:
        return []
    return [item for item in value.split(";") if item]


def compute_pocket_depth(
    centroid: Tuple[float, float, float],
    contact_resnums: set,
    receptor_pdb: Path,
    chain: str = "A",
) -> Optional[float]:
    """Compute pocket depth: centroid → nearest non-contact surface Cα (AC-4.2).

    Pocket depth measures how "buried" the centroid is below the receptor
    surface. Contact residues are excluded so that the distance reflects
    the distance to the actual surface, not to the pocket lining.

    Returns distance in Å, or None if PDB not found or no surface atoms.
    """
    if not receptor_pdb.exists():
        return None

    min_dist = float("inf")
    with open(receptor_pdb, encoding="utf-8") as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            atom_name = line[12:16].strip()
            if atom_name != "CA":
                continue
            pdb_chain = line[21].strip()
            if pdb_chain and pdb_chain != chain:
                continue
            try:
                resnum = int(line[22:26].strip())
            except (ValueError, IndexError):
                continue
            # Skip contact residues — we want the surface, not pocket lining
            if resnum in contact_resnums:
                continue
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except (ValueError, IndexError):
                continue
            dist = math.sqrt(
                (x - centroid[0]) ** 2
                + (y - centroid[1]) ** 2
                + (z - centroid[2]) ** 2
            )
            if dist < min_dist:
                min_dist = dist

    return round(min_dist, 2) if min_dist < float("inf") else None


def summarize_pose_rows(rows: List[dict]) -> Tuple[List[dict], List[dict], List[dict]]:
    # Exclude capped poses when pocket cap was applied
    if rows and "cap_status" in rows[0]:
        rows = [r for r in rows if r.get("cap_status") != "capped"]

    pocket_groups: Dict[Tuple[str, str], List[dict]] = defaultdict(list)
    ligand_groups: Dict[Tuple[str, str], List[dict]] = defaultdict(list)

    n_skipped = 0
    for row in rows:
        if not row.get("pocket_id"):
            n_skipped += 1
            continue
        pocket_groups[(row["receptor_id"], row["pocket_id"])].append(row)
        ligand_groups[(row["receptor_id"], row["ligand_id"])].append(row)
    if n_skipped:
        import logging
        logging.getLogger(__name__).warning(
            "summarize: %d pose(s) without pocket_id excluded from summary",
            n_skipped,
        )

    pocket_rows: List[dict] = []
    residue_occupancy_rows: List[dict] = []
    for (receptor_id, pocket_id), group in sorted(pocket_groups.items()):
        centroids_x = [float(item["centroid_x"]) for item in group]
        centroids_y = [float(item["centroid_y"]) for item in group]
        centroids_z = [float(item["centroid_z"]) for item in group]
        affinities = [float(item["affinity"]) for item in group if item["affinity"] not in ("", None)]
        ligands = sorted({item["ligand_id"] for item in group})
        ligand_counts: Counter = Counter(item["ligand_id"] for item in group)
        residue_counter: Counter = Counter()
        for item in group:
            residue_counter.update(split_contact_residues(item.get("contact_residues", "")))
        union_contact_residues = sorted(residue_counter)
        top_residues = [residue for residue, _ in sorted(residue_counter.items(), key=lambda item: (-item[1], item[0]))[:5]]

        # Centroid spread (RMSD of pose centroids from pocket centroid)
        mean_x = sum(centroids_x) / len(centroids_x)
        mean_y = sum(centroids_y) / len(centroids_y)
        mean_z = sum(centroids_z) / len(centroids_z)
        if len(centroids_x) >= 2:
            spread = math.sqrt(sum(
                (x - mean_x) ** 2 + (y - mean_y) ** 2 + (z - mean_z) ** 2
                for x, y, z in zip(centroids_x, centroids_y, centroids_z)
            ) / len(centroids_x))
        else:
            spread = 0.0

        # Affinity uncertainty
        aff_std = round(statistics.stdev(affinities), 4) if len(affinities) >= 2 else 0.0
        if len(affinities) >= 4:
            s_aff = sorted(affinities)
            aff_iqr = round(s_aff[3 * len(s_aff) // 4] - s_aff[len(s_aff) // 4], 4)
        else:
            aff_iqr = 0.0

        dominant_ligand_fraction = 0.0
        ligand_pose_entropy = 0.0
        if ligand_counts:
            dominant_ligand_fraction = max(ligand_counts.values()) / len(group)
            for count in ligand_counts.values():
                fraction = count / len(group)
                if fraction > 0:
                    ligand_pose_entropy -= fraction * math.log(fraction)

        # ATP site false positive detection (AC-1.2)
        resnums = [extract_resnum(r) for r in union_contact_residues]
        resnums = [n for n in resnums if n is not None]
        atp_overlap = sum(1 for n in resnums if n in ATP_SITE_RESIDUES)
        atp_frac = atp_overlap / len(resnums) if resnums else 0.0
        is_atp_site = atp_frac > 0.5
        is_atp_site_borderline = 0.4 <= atp_frac <= 0.6

        pocket_rows.append({
            "receptor_id": receptor_id,
            "pocket_id": pocket_id,
            "centroid_x": round(mean_x, 4),
            "centroid_y": round(mean_y, 4),
            "centroid_z": round(mean_z, 4),
            "n_pose": len(group),
            "n_ligand": len(ligands),
            "best_affinity": round(min(affinities), 4) if affinities else "",
            "mean_affinity": round(sum(affinities) / len(affinities), 4) if affinities else "",
            "union_contact_residues": ";".join(union_contact_residues),
            "top_residues": ";".join(top_residues),
            "centroid_spread_A": round(spread, 4),
            "affinity_std": aff_std,
            "affinity_iqr": aff_iqr,
            "dominant_ligand_fraction": round(dominant_ligand_fraction, 4),
            "ligand_pose_entropy": round(ligand_pose_entropy, 4),
            "is_atp_site": is_atp_site,
            "is_atp_site_borderline": is_atp_site_borderline,
            # Ko et al. sheet contact counts (AC-1.5)
            "contacts_sheet_8_9": sum(1 for n in resnums if n in KO_SHEET_8_9),
            "contacts_sheet_10_11": sum(1 for n in resnums if n in KO_SHEET_10_11),
            "contacts_sheet_12": sum(1 for n in resnums if n in KO_SHEET_12),
            "pocket_depth_A": "",  # Computed in summarize_from_config when PDB available
        })

        # Build full per-residue occupancy rows for this pocket
        n_poses = len(group)
        for residue_id, n_occ in sorted(residue_counter.items(), key=lambda item: (-item[1], item[0])):
            occ = round(n_occ / n_poses, 4) if n_poses > 0 else 0.0
            residue_occupancy_rows.append({
                "receptor_id": receptor_id,
                "pocket_id": pocket_id,
                "residue_id": residue_id,
                "n_occurrences": n_occ,
                "occupancy": occ,
                "is_hotspot": occ >= 0.5,
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

    return pocket_rows, drug_map_rows, residue_occupancy_rows


def write_csv(path: Path, rows: List[dict], fieldnames: List[str]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return path


def summarize_from_config(config_path: str, pose_table_path: Optional[str] = None) -> Tuple[Path, Path, Path]:
    config = load_config(config_path)
    postprocess_root = paths.wa_phase4_vina_postprocess(config)
    pose_table = Path(pose_table_path) if pose_table_path else postprocess_root / "vina_pose_table.csv"
    rows = load_pose_table(pose_table)

    # Cross-receptor comparison is not implemented yet; residue numbering consistency
    # should be verified on the real server before Task Group 5 uses these summaries.
    pocket_rows, drug_map_rows, residue_occupancy_rows = summarize_pose_rows(rows)

    # AC-4.2: Compute pocket_depth_A from receptor PDBs when available
    receptors = config.get("receptors", [])
    receptor_pdbs = {r["id"]: Path(r.get("pdb", "")) for r in receptors}
    for pocket in pocket_rows:
        rid = pocket["receptor_id"]
        pdb_path = receptor_pdbs.get(rid)
        if pdb_path and pdb_path.exists():
            centroid = (
                float(pocket["centroid_x"]),
                float(pocket["centroid_y"]),
                float(pocket["centroid_z"]),
            )
            contact_str = pocket.get("union_contact_residues", "")
            contact_resnums = set()
            for res_id in contact_str.split(";"):
                try:
                    rn = extract_resnum(res_id.strip())
                    if rn is not None:
                        contact_resnums.add(rn)
                except (ValueError, TypeError):
                    pass
            depth = compute_pocket_depth(centroid, contact_resnums, pdb_path)
            pocket["pocket_depth_A"] = depth if depth is not None else ""

    pocket_csv = write_csv(
        postprocess_root / "vina_pocket_table.csv",
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
            "centroid_spread_A",
            "affinity_std",
            "affinity_iqr",
            "dominant_ligand_fraction",
            "ligand_pose_entropy",
            "is_atp_site",
            "is_atp_site_borderline",
            "contacts_sheet_8_9",
            "contacts_sheet_10_11",
            "contacts_sheet_12",
            "pocket_depth_A",
        ],
    )
    drug_csv = write_csv(
        postprocess_root / "vina_drug_pocket_map.csv",
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
    occupancy_csv = write_csv(
        postprocess_root / "vina_pocket_residue_occupancy.csv",
        residue_occupancy_rows,
        [
            "receptor_id",
            "pocket_id",
            "residue_id",
            "n_occurrences",
            "occupancy",
            "is_hotspot",
        ],
    )
    return pocket_csv, drug_csv, occupancy_csv
