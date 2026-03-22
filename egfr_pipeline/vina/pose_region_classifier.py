"""Pose-level region classification for Vina blind docking bias analysis.

Classifies each Vina pose into one of four EGFR kinase domain regions
based on the majority of its contact residues, then aggregates pose counts,
fractions, and mean affinities per (receptor_id, ligand_id, region).

AC-2.1: vina_pose_distribution_by_region.csv output.
EC-2.1: WARNING when C-lobe surface < 10% or 0%.
"""

import csv
import logging
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from egfr_pipeline.config import load_config
from egfr_pipeline import paths
from egfr_pipeline.region_definitions import get_region
from egfr_pipeline.residue_utils import extract_resnum
from egfr_pipeline.vina.pocket_summary import split_contact_residues

logger = logging.getLogger(__name__)

REGION_ORDER = (
    "n_lobe",
    "atp_site",
    "c_lobe_surface",
    "c_lobe_core",
    "mixed",
    "unknown",
)

DISTRIBUTION_FIELDS = [
    "receptor_id",
    "ligand_id",
    "region",
    "n_poses",
    "fraction",
    "mean_affinity",
]


def _classify_single_pose(contact_residues: str) -> str:
    """Determine the primary region for a single pose.

    Returns the region containing >50% of contact residues,
    ``"mixed"`` if no majority, or ``"unknown"`` if no classifiable residues.
    """
    tokens = split_contact_residues(contact_residues)
    if not tokens:
        return "unknown"

    region_counts: Counter = Counter()
    for token in tokens:
        resnum = extract_resnum(token)
        if resnum is None:
            continue
        region = get_region(resnum)
        if region is not None:
            region_counts[region] += 1

    total = sum(region_counts.values())
    if total == 0:
        return "unknown"

    best_region, best_count = region_counts.most_common(1)[0]
    if best_count > total * 0.5:
        return best_region
    return "mixed"


def classify_poses_by_region(
    rows: List[dict],
) -> Tuple[List[dict], List[str]]:
    """Classify poses by region and aggregate per (receptor_id, ligand_id).

    Args:
        rows: Pose table rows with ``receptor_id``, ``ligand_id``,
              ``contact_residues``, and ``affinity`` fields.

    Returns:
        ``(distribution_rows, warnings)`` where *distribution_rows* has one
        entry per (receptor_id, ligand_id, region) combination and *warnings*
        contains any C-lobe surface coverage alerts.
    """
    if rows and "cap_status" in rows[0]:
        rows = [r for r in rows if r.get("cap_status") != "capped"]

    groups: Dict[Tuple[str, str, str], dict] = defaultdict(
        lambda: {"count": 0, "affinities": []}
    )
    totals: Counter = Counter()

    for row in rows:
        receptor_id = row.get("receptor_id", "")
        ligand_id = row.get("ligand_id", "")
        region = _classify_single_pose(row.get("contact_residues", ""))

        groups[(receptor_id, ligand_id, region)]["count"] += 1
        totals[(receptor_id, ligand_id)] += 1

        aff_str = row.get("affinity", "")
        if aff_str not in ("", None):
            try:
                groups[(receptor_id, ligand_id, region)]["affinities"].append(
                    float(aff_str)
                )
            except (TypeError, ValueError):
                pass

    distribution_rows: List[dict] = []
    warnings: List[str] = []

    for receptor_id, ligand_id in sorted(totals):
        total = totals[(receptor_id, ligand_id)]

        for region in REGION_ORDER:
            data = groups.get((receptor_id, ligand_id, region))
            if data is None or data["count"] == 0:
                continue

            n = data["count"]
            affs = data["affinities"]
            distribution_rows.append({
                "receptor_id": receptor_id,
                "ligand_id": ligand_id,
                "region": region,
                "n_poses": n,
                "fraction": round(n / total, 4) if total else 0.0,
                "mean_affinity": round(sum(affs) / len(affs), 4) if affs else "",
            })

        # C-lobe surface check (AC-2.1 threshold 10%, EC-2.1 zero case)
        c_data = groups.get((receptor_id, ligand_id, "c_lobe_surface"))
        c_count = c_data["count"] if c_data else 0
        c_frac = c_count / total if total else 0.0

        if c_count == 0:
            warnings.append(
                f"WARNING [{receptor_id}/{ligand_id}]: No poses reached C-lobe "
                f"surface. Rely on Workflow B focused docking for this region."
            )
        elif c_frac < 0.10:
            warnings.append(
                f"WARNING [{receptor_id}/{ligand_id}]: C-lobe surface poses "
                f"= {c_count} ({c_frac:.1%} of {total}), below 10% threshold."
            )

    return distribution_rows, warnings


def classify_from_config(
    config_path: str, pose_table_path: Optional[str] = None
) -> Path:
    """Run region classification and write CSV to the postprocess directory.

    Args:
        config_path: Path to the project YAML config.
        pose_table_path: Override path to vina_pose_table.csv.

    Returns:
        Path to the generated vina_pose_distribution_by_region.csv.
    """
    config = load_config(config_path)
    postprocess_root = paths.wa_phase4_vina_postprocess(config)
    pose_table = (
        Path(pose_table_path)
        if pose_table_path
        else postprocess_root / "vina_pose_table.csv"
    )

    with open(pose_table, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))

    dist_rows, warnings = classify_poses_by_region(rows)
    for w in warnings:
        logger.warning(w)

    out_path = postprocess_root / "vina_pose_distribution_by_region.csv"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=DISTRIBUTION_FIELDS)
        writer.writeheader()
        writer.writerows(dist_rows)

    logger.info(
        "Pose region distribution: %d rows written to %s", len(dist_rows), out_path
    )
    return out_path
