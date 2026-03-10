"""Multi-evidence pocket evaluation: evidence strength classification.

Classifies each Vina pocket by evidence strength (STRONG / MODERATE / WEAK)
based on three independent axes. This is NOT a validity judgment — it is
an evidence summary to guide the researcher's manual inspection priority.

  Axis 1 — Vina Quality (direct evidence)
    Drug binding affinity, pose convergence, multi-ligand consensus.
    Calibrated for EGFR C-lobe surface pockets (-6.5 ~ -8.0 kcal/mol range).

  Axis 2 — PPI Spatial Proximity (indirect evidence, when available)
    3D distance between Vina pocket centroid and PPI interface centroid.
    Answers: "Is this drug pocket near the MYO1D binding surface?"
    Drug pockets and PPI interfaces are fundamentally different binding
    modes — spatial proximity is more meaningful than residue overlap.

  Axis 3 — Cross-Receptor Consistency (structural evidence)
    Same pocket found across multiple receptor conformational states
    (3GT8_raw, MD cl38_48, MD cl85_100). Indicates structural stability.

Adaptive scoring: when PPI data is absent, scoring redistributes weights
so pockets are never penalized for missing data.

Outputs:
  - cross_method_agreement.csv  (Vina ↔ PPI spatial + residue analysis)
  - valid_sites.csv             (evidence classification per pocket)
"""
import csv
import math
import re
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

from egfr_pipeline.config import load_config, project_root_from_config
from egfr_pipeline.residue_utils import (
    normalize_residue_id,
    extract_resnum,
    parse_residue_set,
)

# ---------------------------------------------------------------------------
# CSV field definitions
# ---------------------------------------------------------------------------

AGREEMENT_FIELDS = [
    "receptor_id",
    "pocket_id",
    "n_vina_residues",
    "n_ppi_residues",
    "n_shared_residues",
    "jaccard",
    "overlap_coeff",
    "shared_residue_list",
    "ppi_mean_occupancy_of_shared",
    "spatial_dist_A",
    "spatial_proximity",
    "vina_best_affinity_kcal",
    "ppi_best_dg_REU",
    "agreement_level",
]

VERDICT_FIELDS = [
    "receptor_id",
    "pocket_id",
    "verdict",
    "confidence_score",
    "vina_quality_score",
    "ppi_proximity_score",
    "cross_receptor_score",
    "ppi_data_available",
    "best_affinity",
    "n_pose",
    "n_ligand",
    "spatial_dist_to_ppi",
    "n_shared_with_ppi",
    "cross_receptor_matches",
    "reasons",
]

# ---------------------------------------------------------------------------
# Default scoring thresholds (overridable via config.yaml verdict section)
# ---------------------------------------------------------------------------

DEFAULT_THRESHOLDS = {
    # --- Axis 1: Vina quality ---
    # With PPI data: max 50 pts; Without PPI: max 60 pts (adaptive)
    #
    # EGFR C-lobe context:
    #   ATP binding site (active site): -9 ~ -12 kcal/mol (not our target)
    #   C-lobe surface pockets: -5 ~ -8 kcal/mol (shallow, our main target)
    #   Allosteric sites: -6 ~ -9 kcal/mol
    # Thresholds calibrated for C-lobe surface, not active site.
    "affinity_good": -6.5,       # surface pocket hit
    "affinity_great": -8.0,      # strong surface binding (rare on C-lobe)
    "n_pose_good": 3,            # 30% convergence in blind docking
    "n_pose_great": 8,           # high convergence
    "n_ligand_good": 2,          # 2/3 ligands → cross-chemical consensus
    "n_ligand_all": 3,           # all 3 ligands

    # --- Axis 2: PPI spatial proximity ---
    # Max 20 pts (only when PPI data available)
    #
    # EGFR kinase domain is ~40 Å across.
    # MYO1D binds C-lobe surface (not ATP site).
    # "Adjacent" = drug pocket directly at PPI interface (competition/overlap)
    # "Near" = drug pocket on same face, potential allosteric crosstalk
    # "Moderate" = same domain but different face
    "ppi_dist_adjacent": 8.0,    # < 8 Å — direct PPI interface overlap
    "ppi_dist_near": 15.0,       # < 15 Å — same face, allosteric range
    "ppi_dist_moderate": 25.0,   # < 25 Å — same domain (~half of 40 Å width)
    "ppi_residue_bonus": True,   # extra credit for shared residues

    # --- Axis 3: Cross-receptor consistency ---
    # With PPI: max 30 pts; Without PPI: max 40 pts (adaptive)
    #
    # 3 receptor states: 3GT8_raw, MD cl38_48, MD cl85_100
    # MD clusters represent conformational sampling.
    # Same pocket in 2+ states = structurally stable site (not artifact).
    "cross_receptor_centroid_cutoff": 8.0,

    # --- Evidence level thresholds (always out of 100) ---
    # NOT a validity judgment — an evidence strength classification.
    # Researcher must still inspect STRONG pockets visually.
    "valid_min": 55,             # STRONG evidence
    "uncertain_min": 30,         # MODERATE evidence
    #                            # < 30 = WEAK evidence
}


def _get_thresholds(config: dict) -> dict:
    thresholds = dict(DEFAULT_THRESHOLDS)
    user = config.get("verdict", {})
    for key in DEFAULT_THRESHOLDS:
        if key in user:
            thresholds[key] = user[key]
    return thresholds


# ---------------------------------------------------------------------------
# CSV I/O helpers
# ---------------------------------------------------------------------------

def _load_csv(path: Path) -> List[dict]:
    if not path.exists():
        return []
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def _write_csv(path: Path, rows: List[dict], fieldnames: List[str]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return path


# ---------------------------------------------------------------------------
# PDB coordinate extraction (for PPI interface centroid)
# ---------------------------------------------------------------------------

def _parse_ca_coords(pdb_path: Path) -> Dict[int, Tuple[float, float, float]]:
    """Extract CA atom coordinates from a PDB, keyed by residue number.

    Only reads chain A CA atoms. Returns {resnum: (x, y, z)}.
    """
    coords: Dict[int, Tuple[float, float, float]] = {}
    try:
        with open(pdb_path, encoding="utf-8", errors="ignore") as f:
            for line in f:
                if not line.startswith("ATOM"):
                    continue
                atom_name = line[12:16].strip()
                chain = line[21]
                if atom_name != "CA" or chain not in ("A", "X"):
                    continue
                try:
                    resnum = int(line[22:26].strip())
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords[resnum] = (x, y, z)
                except ValueError:
                    continue
    except (OSError, IOError):
        pass
    return coords


def _compute_centroid_and_spread(
    coords: Dict[int, Tuple[float, float, float]],
    resnums: Set[int],
) -> Tuple[Optional[Tuple[float, float, float]], float]:
    """Compute centroid and spatial spread (RMSD from centroid) of CA atoms.

    Returns (centroid, spread_A). High spread (>15 Å) indicates the PPI
    interface is dispersed across multiple patches, making the centroid
    a poor representative — spatial proximity to this centroid should be
    interpreted with caution.
    """
    found = [coords[r] for r in resnums if r in coords]
    if not found:
        return None, 0.0
    n = len(found)
    cx = sum(c[0] for c in found) / n
    cy = sum(c[1] for c in found) / n
    cz = sum(c[2] for c in found) / n
    centroid = (cx, cy, cz)
    # RMSD from centroid (spread)
    spread = math.sqrt(
        sum((c[0]-cx)**2 + (c[1]-cy)**2 + (c[2]-cz)**2 for c in found) / n
    )
    return centroid, round(spread, 2)


def _euclidean_dist(
    a: Tuple[float, float, float],
    b: Tuple[float, float, float],
) -> float:
    return math.sqrt((a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2)


# ---------------------------------------------------------------------------
# Evidence loading
# ---------------------------------------------------------------------------

def load_all_evidence(project_root: Path) -> dict:
    return {
        "pocket_table": _load_csv(project_root / "vina_pocket_table.csv"),
        "drug_pocket_map": _load_csv(project_root / "vina_drug_pocket_map.csv"),
        "pocket_comparison": _load_csv(project_root / "vina_pocket_comparison.csv"),
        "ppi_residues": _load_csv(project_root / "ppi_pyrosetta_residues.csv"),
        "ppi_summary": _load_csv(project_root / "ppi_pyrosetta_summary.csv"),
    }


# ---------------------------------------------------------------------------
# PPI interface centroid computation
# ---------------------------------------------------------------------------

def _build_ppi_interface_centroids(
    ppi_residues: List[dict],
    config: dict,
) -> Tuple[Dict[str, Tuple[float, float, float]], Dict[str, float]]:
    """Compute PPI interface centroid per receptor from CA coordinates.

    Uses receptor PDB files from config to look up coordinates.
    Returns {receptor_id: (cx, cy, cz)}.

    IMPORTANT — Coordinate system assumption:
    This function reads CA coords from the receptor PDB listed in config.
    The resulting centroid is compared to Vina pocket centroids, which come
    from the PDBQT used by Vina. Both must share the same coordinate frame.

    For crystal structures (3GT8_raw), PDB→PDBQT conversion preserves coords.
    For MD cluster representatives (cl38_48, cl85_100), the PDB may have a
    different origin if extracted from a trajectory without re-alignment.
    If PPI is ever run on MD clusters, ensure the PDB used here is the same
    structure (or aligned to the same reference) as the Vina PDBQT input.
    """
    # Group PPI residue numbers by receptor
    ppi_resnums: Dict[str, Set[int]] = defaultdict(set)
    for row in ppi_residues:
        rid = row["receptor_id"]
        resnum = _safe_int(row.get("residue_num"), 0)
        if resnum > 0:
            ppi_resnums[rid].add(resnum)

    if not ppi_resnums:
        return {}

    # Build receptor PDB path index
    receptor_pdbs: Dict[str, Path] = {}
    for rec in config.get("receptors", []):
        pdb_path = rec.get("pdb", "")
        if pdb_path:
            receptor_pdbs[rec["id"]] = Path(pdb_path)

    centroids: Dict[str, Tuple[float, float, float]] = {}
    spreads: Dict[str, float] = {}
    for rid, resnums in ppi_resnums.items():
        pdb_path = receptor_pdbs.get(rid)
        if not pdb_path or not pdb_path.exists():
            continue
        ca_coords = _parse_ca_coords(pdb_path)
        centroid, spread = _compute_centroid_and_spread(ca_coords, resnums)
        if centroid:
            centroids[rid] = centroid
            spreads[rid] = spread

    return centroids, spreads


# ---------------------------------------------------------------------------
# PPI residue index
# ---------------------------------------------------------------------------

def _build_ppi_residue_index(
    ppi_residues: List[dict],
) -> Dict[str, Dict[str, dict]]:
    """Index PPI residues by receptor_id -> normalized_residue_id -> row."""
    index: Dict[str, Dict[str, dict]] = defaultdict(dict)
    for row in ppi_residues:
        rid = row["receptor_id"]
        res_id = normalize_residue_id(row["residue_id"])
        index[rid][res_id] = row
    return dict(index)


# ---------------------------------------------------------------------------
# Offset residue detection
# ---------------------------------------------------------------------------

def check_ppi_residue_offsets(ppi_residues: List[dict]) -> List[str]:
    """Detect PPI residues with suspiciously high numbers (>1700).

    These likely indicate unrestored chain B offset residues.
    Returns list of warning messages.
    """
    warnings = []
    for row in ppi_residues:
        resnum = _safe_int(row.get("residue_num"), 0)
        if resnum > 1700:
            warnings.append(
                f"  Offset residue detected: {row.get('receptor_id')} "
                f"{row.get('residue_id')} (num={resnum}) — "
                f"run PPI postprocess to restore chain numbering"
            )
    return warnings


# ---------------------------------------------------------------------------
# Step 1: Cross-method agreement (spatial + residue)
# ---------------------------------------------------------------------------

def compute_cross_method_agreement(
    pocket_rows: List[dict],
    ppi_residues: List[dict],
    ppi_centroids: Dict[str, Tuple[float, float, float]],
    thresholds: Optional[dict] = None,
) -> List[dict]:
    """Compute Vina ↔ PPI agreement for each pocket.

    Includes both residue overlap (informational) and spatial proximity
    (distance between pocket centroid and PPI interface centroid).

    Note on centroid semantics: Vina pocket centroids are ligand-atom
    averages (inside binding cavity), while PPI centroids are receptor
    CA-atom averages (protein surface). The Vina centroid sits ~3-5 Å
    deeper than the pocket entrance, so true surface-to-surface distance
    may be shorter than the computed centroid-to-centroid distance.
    Thresholds account for this systematic offset.
    """
    T = thresholds or DEFAULT_THRESHOLDS
    ppi_index = _build_ppi_residue_index(ppi_residues)
    # Collect all receptors that have PPI data
    ppi_receptors = set(ppi_index.keys()) | set(ppi_centroids.keys())

    results: List[dict] = []

    for pocket in pocket_rows:
        receptor_id = pocket["receptor_id"]
        pocket_id = pocket["pocket_id"]

        # Vina residues
        vina_residues = parse_residue_set(pocket.get("union_contact_residues", ""))

        # PPI residues for this receptor
        ppi_res_map = ppi_index.get(receptor_id, {})
        ppi_residue_set = set(ppi_res_map.keys())

        # Residue overlap (informational — weak signal for drug vs protein)
        shared = vina_residues & ppi_residue_set
        n_vina = len(vina_residues)
        n_ppi = len(ppi_residue_set)
        n_shared = len(shared)

        union_size = len(vina_residues | ppi_residue_set)
        jaccard = n_shared / union_size if union_size > 0 else 0.0
        min_size = min(n_vina, n_ppi)
        overlap = n_shared / min_size if min_size > 0 else 0.0

        # Mean occupancy of shared residues
        mean_occ = 0.0
        if shared:
            occs = []
            for res_id in shared:
                occ_val = ppi_res_map.get(res_id, {}).get("occupancy", "0")
                try:
                    occs.append(float(occ_val))
                except (ValueError, TypeError):
                    pass
            mean_occ = sum(occs) / len(occs) if occs else 0.0

        # Spatial proximity (primary PPI signal)
        spatial_dist = None
        spatial_proximity = "no_data"
        ppi_centroid = ppi_centroids.get(receptor_id)
        if ppi_centroid:
            try:
                pocket_centroid = (
                    float(pocket.get("centroid_x", 0)),
                    float(pocket.get("centroid_y", 0)),
                    float(pocket.get("centroid_z", 0)),
                )
                spatial_dist = round(_euclidean_dist(pocket_centroid, ppi_centroid), 2)
                if spatial_dist < T["ppi_dist_adjacent"]:
                    spatial_proximity = "adjacent"
                elif spatial_dist < T["ppi_dist_near"]:
                    spatial_proximity = "near"
                elif spatial_dist < T["ppi_dist_moderate"]:
                    spatial_proximity = "moderate"
                else:
                    spatial_proximity = "distant"
            except (ValueError, TypeError):
                pass

        # Agreement level (combining spatial + residue)
        if spatial_proximity in ("adjacent",) and n_shared > 0:
            agreement_level = "strong"
        elif spatial_proximity in ("adjacent", "near"):
            agreement_level = "moderate"
        elif n_shared > 0:
            agreement_level = "weak"
        elif receptor_id in ppi_receptors:
            agreement_level = "none"
        else:
            agreement_level = "no_data"

        results.append({
            "receptor_id": receptor_id,
            "pocket_id": pocket_id,
            "n_vina_residues": n_vina,
            "n_ppi_residues": n_ppi,
            "n_shared_residues": n_shared,
            "jaccard": round(jaccard, 4),
            "overlap_coeff": round(overlap, 4),
            "shared_residue_list": ";".join(sorted(shared)),
            "ppi_mean_occupancy_of_shared": round(mean_occ, 4),
            "spatial_dist_A": spatial_dist if spatial_dist is not None else "",
            "spatial_proximity": spatial_proximity,
            "vina_best_affinity_kcal": pocket.get("best_affinity", ""),
            "ppi_best_dg_REU": "",
            "agreement_level": agreement_level,
        })

    return results


# ---------------------------------------------------------------------------
# Step 2: Cross-receptor consistency
# ---------------------------------------------------------------------------

def compute_cross_receptor_consistency(
    comparison_rows: List[dict],
) -> Dict[Tuple[str, str], List[str]]:
    """For each pocket, find which other receptors have a same_patch_candidate."""
    matches: Dict[Tuple[str, str], Set[str]] = defaultdict(set)
    for row in comparison_rows:
        is_candidate = str(row.get("same_patch_candidate", "")).lower()
        if is_candidate not in ("true", "1", "yes"):
            continue
        rec_a, rec_b = row["receptor_a"], row["receptor_b"]
        pkt_a, pkt_b = row["pocket_a"], row["pocket_b"]
        matches[(rec_a, pkt_a)].add(rec_b)
        matches[(rec_b, pkt_b)].add(rec_a)
    return {k: sorted(v) for k, v in matches.items()}


# ---------------------------------------------------------------------------
# Membrane-proximal surface check
# ---------------------------------------------------------------------------

# EGFR C-lobe membrane-facing residues (from structure inspection).
# Vina has no such constraint, so drug pockets can appear on membrane face.
# Tag these as informational warning — not scored, just flagged.
_MEMBRANE_RESIDUES = set()
for _rng in ["709-720", "724-731", "736-739", "783-785", "799-805",
             "871-873", "917-921"]:
    _parts = _rng.split("-")
    _MEMBRANE_RESIDUES.update(range(int(_parts[0]), int(_parts[1]) + 1))


def _check_membrane_overlap(pocket: dict) -> Optional[str]:
    """Check if pocket contact residues overlap with membrane-proximal surface.

    Returns a warning tag like 'membrane_face(3res)' or None.
    """
    union_raw = pocket.get("union_contact_residues", "")
    if not union_raw:
        return None
    membrane_count = 0
    for res in union_raw.split(";"):
        res = res.strip()
        if not res:
            continue
        # Extract residue number
        num_match = re.search(r'(\d+)$', res)
        if num_match:
            resnum = int(num_match.group(1))
            if resnum in _MEMBRANE_RESIDUES:
                membrane_count += 1
    if membrane_count > 0:
        return f"membrane_face({membrane_count}res)"
    return None


# ---------------------------------------------------------------------------
# Step 3: Per-pocket scoring (adaptive)
# ---------------------------------------------------------------------------

def score_pocket(
    pocket: dict,
    ppi_agreement: Optional[dict],
    cross_receptor_matches: List[str],
    thresholds: dict,
    has_ppi_data: bool,
) -> Tuple[float, str, List[str], float, float, float]:
    """Score a single pocket with adaptive weighting.

    Scoring adapts to available evidence:
      With PPI data:    Vina(50) + PPI_proximity(20) + Cross_receptor(30) = 100
      Without PPI data: Vina(60) + Cross_receptor(40) = 100

    This ensures pockets are never penalized for missing PPI data.
    """
    reasons: List[str] = []
    T = thresholds

    # ---- Axis 1: Vina Quality ----
    affinity = _safe_float(pocket.get("best_affinity"), 0.0)
    n_pose = _safe_int(pocket.get("n_pose"), 0)
    n_ligand = _safe_int(pocket.get("n_ligand"), 0)

    vina_raw = 0.0

    # Affinity: graduated scoring (0 / 10 / 20)
    if affinity <= T["affinity_great"]:
        vina_raw += 20
        reasons.append(f"affinity={affinity:.1f}")
    elif affinity <= T["affinity_good"]:
        vina_raw += 10
        reasons.append(f"affinity={affinity:.1f}")

    # Pose convergence: graduated (0 / 8 / 15)
    if n_pose >= T["n_pose_great"]:
        vina_raw += 15
        reasons.append(f"n_pose={n_pose}")
    elif n_pose >= T["n_pose_good"]:
        vina_raw += 8
        reasons.append(f"n_pose={n_pose}")

    # Multi-ligand consensus: graduated (0 / 10 / 20)
    if n_ligand >= T["n_ligand_all"]:
        vina_raw += 20
        reasons.append(f"n_ligand={n_ligand}")
    elif n_ligand >= T["n_ligand_good"]:
        vina_raw += 10
        reasons.append(f"n_ligand={n_ligand}")

    # Membrane-face tag (informational, not scored)
    membrane_tag = _check_membrane_overlap(pocket)
    if membrane_tag:
        reasons.append(membrane_tag)

    # Normalize to adaptive max
    # Raw max = affinity(20) + convergence(15) + consensus(20) = 55
    VINA_RAW_MAX = 20.0 + 15.0 + 20.0  # keep in sync with components above
    vina_max = 50.0 if has_ppi_data else 60.0
    vina_score = min(vina_raw, VINA_RAW_MAX) / VINA_RAW_MAX * vina_max

    # ---- Axis 2: PPI Spatial Proximity ----
    # Note: spatial_dist is centroid-to-centroid. Vina centroid sits inside
    # the pocket cavity (~3-5 Å below surface), so actual pocket-entrance-to-
    # PPI-surface distance is shorter. Thresholds are set conservatively to
    # account for this systematic overestimate.
    ppi_score = 0.0
    ppi_max = 20.0 if has_ppi_data else 0.0

    if has_ppi_data and ppi_agreement:
        spatial = ppi_agreement.get("spatial_proximity", "no_data")
        n_shared = _safe_int(ppi_agreement.get("n_shared_residues"), 0)
        spatial_dist = ppi_agreement.get("spatial_dist_A", "")

        if spatial == "adjacent":
            ppi_score += 15.0
            reasons.append(f"ppi_adjacent({spatial_dist}A)")
        elif spatial == "near":
            ppi_score += 10.0
            reasons.append(f"ppi_near({spatial_dist}A)")
        elif spatial == "moderate":
            ppi_score += 4.0
            reasons.append(f"ppi_moderate({spatial_dist}A)")
        elif spatial == "no_data":
            # PPI residues exist but no PDB for centroid → use residue overlap
            if n_shared > 0:
                ppi_score += 8.0
                reasons.append(f"ppi_shared={n_shared}res")
        else:
            reasons.append(f"ppi_distant({spatial_dist}A)")

        # Residue overlap bonus (weak but informative)
        if T.get("ppi_residue_bonus") and n_shared > 0:
            occ = _safe_float(ppi_agreement.get("ppi_mean_occupancy_of_shared"), 0)
            if occ >= 0.5:
                ppi_score += 5.0
                reasons.append(f"shared_highocc={n_shared}")
            else:
                ppi_score += 2.0

        ppi_score = min(ppi_score, ppi_max)
    elif not has_ppi_data:
        reasons.append("no_ppi_data")

    # ---- Axis 3: Cross-Receptor Consistency ----
    cross_max = 30.0 if has_ppi_data else 40.0
    n_cross = len(cross_receptor_matches)
    cross_raw = 0.0

    if n_cross >= 2:
        cross_raw = cross_max
        reasons.append(f"cross={n_cross + 1}/3")
    elif n_cross >= 1:
        cross_raw = cross_max * 0.5
        reasons.append(f"cross=2/3")

    cross_score = min(cross_raw, cross_max)

    # ---- Total (always out of 100) ----
    total = vina_score + ppi_score + cross_score

    if total >= T["valid_min"]:
        verdict = "STRONG"
    elif total >= T["uncertain_min"]:
        verdict = "MODERATE"
    else:
        verdict = "WEAK"

    return total, verdict, reasons, vina_score, ppi_score, cross_score


# ---------------------------------------------------------------------------
# Step 4: Main entry point
# ---------------------------------------------------------------------------

def generate_verdict(
    config_path: str,
    output_dir: Optional[str] = None,
) -> Tuple[Path, Path]:
    """Run full verdict pipeline.

    Returns: (cross_method_agreement_csv_path, valid_sites_csv_path)
    """
    config = load_config(config_path)
    project_root = project_root_from_config(config)
    out_dir = Path(output_dir) if output_dir else project_root
    thresholds = _get_thresholds(config)

    # Load evidence
    evidence = load_all_evidence(project_root)
    pocket_rows = evidence["pocket_table"]
    if not pocket_rows:
        print("[verdict] No pocket table found — skipping.")
        return Path(), Path()

    ppi_residues = evidence["ppi_residues"]

    # Check for offset residues (unrestored chain B)
    offset_warnings = check_ppi_residue_offsets(ppi_residues)
    if offset_warnings:
        print("\n[verdict] WARNING: Unrestored chain B offset residues detected!")
        for w in offset_warnings[:5]:
            print(w)
        if len(offset_warnings) > 5:
            print(f"  ... and {len(offset_warnings) - 5} more")
        print("  Run PPI Postprocess (option 8) to fix before verdict.\n")

    # Determine which receptors have PPI data
    ppi_receptor_ids = {r["receptor_id"] for r in ppi_residues}
    has_any_ppi = len(ppi_receptor_ids) > 0

    # Compute PPI interface centroids from receptor PDBs
    ppi_centroids, ppi_spreads = _build_ppi_interface_centroids(ppi_residues, config)

    if has_any_ppi:
        print(f"[verdict] PPI data available for: {', '.join(sorted(ppi_receptor_ids))}")
        if ppi_centroids:
            for rid, c in sorted(ppi_centroids.items()):
                spread = ppi_spreads.get(rid, 0)
                spread_warn = " (DISPERSED — centroid unreliable)" if spread > 15.0 else ""
                print(f"  PPI interface centroid ({rid}): "
                      f"({c[0]:.1f}, {c[1]:.1f}, {c[2]:.1f})  "
                      f"spread={spread:.1f}A{spread_warn}")
        else:
            print("  (No receptor PDBs found for centroid — using residue overlap only)")

    # Step 1: Cross-method agreement
    agreement_rows = compute_cross_method_agreement(
        pocket_rows, ppi_residues, ppi_centroids, thresholds,
    )

    # Enrich with PPI summary best_dg (receptor-level, not pocket-specific;
    # same value for all pockets of a receptor — use for reference only)
    ppi_summary_map = {s.get("receptor_id", ""): s for s in evidence["ppi_summary"]}
    for row in agreement_rows:
        summary = ppi_summary_map.get(row["receptor_id"])
        if summary:
            row["ppi_best_dg_REU"] = summary.get("best_dg", "")

    agreement_csv = _write_csv(
        out_dir / "cross_method_agreement.csv",
        agreement_rows,
        AGREEMENT_FIELDS,
    )

    # Step 2: Cross-receptor consistency
    cross_receptor = compute_cross_receptor_consistency(
        evidence["pocket_comparison"]
    )

    # Step 3: Build lookup indices
    agreement_index = {
        (r["receptor_id"], r["pocket_id"]): r for r in agreement_rows
    }
    pocket_ligand_count = _count_pocket_ligands(evidence["drug_pocket_map"])

    # Step 4: Score each pocket
    verdict_rows: List[dict] = []
    for pocket in pocket_rows:
        rid = pocket["receptor_id"]
        pid = pocket["pocket_id"]
        key = (rid, pid)

        if not pocket.get("n_ligand"):
            pocket["n_ligand"] = pocket_ligand_count.get(key, 0)

        ppi_agr = agreement_index.get(key)
        cross_matches = cross_receptor.get(key, [])

        # Per-pocket PPI availability check
        pocket_has_ppi = rid in ppi_receptor_ids

        total, verdict, reasons, v_score, p_score, c_score = score_pocket(
            pocket, ppi_agr, cross_matches, thresholds, pocket_has_ppi,
        )

        spatial_dist = ""
        if ppi_agr:
            spatial_dist = ppi_agr.get("spatial_dist_A", "")

        verdict_rows.append({
            "receptor_id": rid,
            "pocket_id": pid,
            "verdict": verdict,
            "confidence_score": round(total, 1),
            "vina_quality_score": round(v_score, 1),
            "ppi_proximity_score": round(p_score, 1),
            "cross_receptor_score": round(c_score, 1),
            "ppi_data_available": "yes" if pocket_has_ppi else "no",
            "best_affinity": pocket.get("best_affinity", ""),
            "n_pose": pocket.get("n_pose", ""),
            "n_ligand": pocket.get("n_ligand", ""),
            "spatial_dist_to_ppi": spatial_dist,
            "n_shared_with_ppi": ppi_agr.get("n_shared_residues", 0) if ppi_agr else 0,
            "cross_receptor_matches": ";".join(cross_matches),
            "reasons": "; ".join(reasons),
        })

    # Sort: VALID first, then by score descending
    verdict_order = {"STRONG": 0, "MODERATE": 1, "WEAK": 2}
    verdict_rows.sort(
        key=lambda r: (verdict_order.get(r["verdict"], 9), -r["confidence_score"])
    )

    verdict_csv = _write_csv(
        out_dir / "valid_sites.csv",
        verdict_rows,
        VERDICT_FIELDS,
    )

    _print_summary(verdict_rows, ppi_receptor_ids)
    return agreement_csv, verdict_csv


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _count_pocket_ligands(
    drug_map_rows: List[dict],
) -> Dict[Tuple[str, str], int]:
    counter: Dict[Tuple[str, str], Set[str]] = defaultdict(set)
    for row in drug_map_rows:
        rid = row["receptor_id"]
        lid = row["ligand_id"]
        dominant = row.get("dominant_pocket_id", "")
        if dominant:
            counter[(rid, dominant)].add(lid)
        alternatives = row.get("alternative_pockets", "")
        if alternatives:
            for pid in alternatives.split(";"):
                pid = pid.strip()
                if pid:
                    counter[(rid, pid)].add(lid)
    return {k: len(v) for k, v in counter.items()}


def _safe_float(val, default: float = 0.0) -> float:
    try:
        return float(val)
    except (TypeError, ValueError):
        return default


def _safe_int(val, default: int = 0) -> int:
    try:
        return int(val)
    except (TypeError, ValueError):
        try:
            return int(float(val))
        except (TypeError, ValueError):
            return default


def _print_summary(
    verdict_rows: List[dict],
    ppi_receptor_ids: Set[str],
) -> None:
    counts = defaultdict(int)
    for row in verdict_rows:
        counts[row["verdict"]] += 1

    total = len(verdict_rows)
    n_receptors = len(set(r["receptor_id"] for r in verdict_rows))
    n_ppi = len(ppi_receptor_ids)

    print(f"\n{'='*65}")
    print(f"  Site Verdict Summary")
    print(f"  {total} pockets across {n_receptors} receptors  "
          f"(PPI data: {n_ppi}/{n_receptors} receptors)")
    print(f"{'='*65}")
    print(f"  STRONG:   {counts.get('STRONG', 0)}")
    print(f"  MODERATE: {counts.get('MODERATE', 0)}")
    print(f"  WEAK:     {counts.get('WEAK', 0)}")

    if n_ppi < n_receptors:
        print(f"\n  Note: {n_receptors - n_ppi} receptor(s) scored without PPI data")
        print(f"        (weights redistributed — no penalty)")

    print(f"{'='*65}")
    for row in verdict_rows[:5]:
        ppi_flag = " [+PPI]" if row.get("ppi_data_available") == "yes" else ""
        print(
            f"  {row['receptor_id']:20s} {row['pocket_id']:6s} "
            f"{row['verdict']:10s} {row['confidence_score']:5.1f}/100{ppi_flag}  "
            f"{row['reasons']}"
        )
    if total > 5:
        print(f"  ... and {total - 5} more")
    print()
