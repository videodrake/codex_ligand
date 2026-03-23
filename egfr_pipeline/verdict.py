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

from egfr_pipeline.config import load_config
from egfr_pipeline import paths
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
    "closest_ppi_partner",
    "n_ppi_partners_near",
    "ppi_frac_runs_supporting",
    "ppi_best_frac_runs_supporting",
    "ppi_best_interface_delta_e",
    "ppi_shared_partner_count",
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
    "vina_affinity_pts",
    "vina_convergence_pts",
    "vina_stability_pts",
    "vina_diversity_pts",
    "vina_consensus_pts",
    "ppi_spatial_pts",
    "ppi_overlap_pts",
    "ppi_reproducibility_pts",
    "cross_receptor_pts",
    "cross_receptor_support_pts",
    "score_denominator",
    "ppi_data_available",
    "best_affinity",
    "n_pose",
    "n_ligand",
    "dominant_ligand_fraction",
    "ligand_pose_entropy",
    "spatial_dist_to_ppi",
    "closest_ppi_partner",
    "n_ppi_partners_near",
    "n_shared_with_ppi",
    "ppi_frac_runs_supporting",
    "ppi_best_interface_delta_e",
    "cross_receptor_matches",
    "cross_receptor_support",
    "consensus_site_id",
    "exp_sensitivity",
    "exp_specificity",
    "exp_enrichment",
    "exp_rank_impact",
    "pocket_stability",
    "bootstrap_confidence",
    "evidence_profile",
    "reason_tags",
    "decision_trace",
    "reasons",
    "is_atp_site",
    "exclusion_reason",
    "allosteric_candidate",   # AC-5.3: Vina ≥ 35 AND PPI ≤ 5
]

CONSENSUS_FIELDS = [
    "consensus_site_id",
    "n_receptors",
    "receptor_list",
    "pocket_list",
    "centroid_x",
    "centroid_y",
    "centroid_z",
    "best_affinity",
    "total_n_ligand",
    "total_n_pose",
    "consensus_residues",
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
    "stability_min": 0.40,       # weakly recurring cluster
    "stability_good": 0.60,      # bootstrap cluster recovery
    "stability_great": 0.80,     # robust cluster recovery
    "n_ligand_good": 2,          # 2/3 ligands → cross-chemical consensus
    "n_ligand_all": 3,           # all 3 ligands
    "dominant_ligand_good": 0.70,  # lower = less single-ligand domination
    "dominant_ligand_warn": 0.85,
    "ligand_entropy_good": 0.90,   # broader ligand distribution
    "ligand_entropy_ok": 0.50,
    # Sub-score max point allocation (sensitivity analysis 용).
    # 합계(VINA_RAW_MAX)가 정규화 분모가 됨.
    "vina_affinity_max": 20.0,
    "vina_convergence_max": 15.0,
    "vina_stability_max": 10.0,
    "vina_diversity_max": 15.0,

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
    "ppi_repro_good": 0.50,      # mean frac_runs_supporting of shared residues
    "ppi_repro_great": 0.75,

    # --- Axis 3: Cross-receptor consistency ---
    # With PPI: max 30 pts; Without PPI: max 40 pts (adaptive)
    #
    # 3 receptor states: 3GT8_raw, MD cl38_48, MD cl85_100
    # MD clusters represent conformational sampling.
    # Same pocket in 2+ states = structurally stable site (not artifact).
    "cross_receptor_centroid_cutoff": 8.0,
    "cross_support_good": 0.50,
    "cross_support_great": 0.75,
    # Coverage differentiation (AC-4.3): 3/3 vs 2/3 vs 1/3 states
    # n_cross = len(cross_receptor_matches): 0=1/3, 1=2/3, 2=3/3
    # Total cross_receptor_pts = coverage + support(max 10):
    #   3/3: 20 + [0-10] = 20-30 (PRD target: 30)
    #   2/3: 14 + [0-10] = 14-24 (PRD target: 22-24)
    #   1/3:  0 + [0-10] =  0-10 (PRD baseline: 10-15)
    # 1/3 coverage set to 0 (not PRD's 10): single-state pockets should earn
    # cross points only through support quality, not mere existence.
    "cross_coverage_3of3": 20.0,     # pocket found in all 3 states
    "cross_coverage_2of3": 14.0,     # pocket found in 2 of 3 states
    "cross_coverage_1of3": 0.0,      # pocket found in only 1 state (coverage via support only)

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


def _bootstrap_confidence(pocket_exists_frac: object) -> str:
    """Map bootstrap pocket_exists_frac to a categorical confidence level.

    AC-2.3: high (>= 0.80), medium (0.50-0.80), low (< 0.50), not_assessed.
    """
    if pocket_exists_frac in ("", None):
        return "not_assessed"
    try:
        frac = float(pocket_exists_frac)
    except (TypeError, ValueError):
        return "not_assessed"
    if frac >= 0.80:
        return "high"
    if frac >= 0.50:
        return "medium"
    return "low"


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

def load_all_evidence(vina_dir: Path, ppi_dir: Path) -> dict:
    return {
        "pocket_table": _load_csv(vina_dir / "vina_pocket_table.csv"),
        "drug_pocket_map": _load_csv(vina_dir / "vina_drug_pocket_map.csv"),
        "pocket_comparison": _load_csv(vina_dir / "vina_pocket_comparison.csv"),
        "ppi_residues": _load_csv(ppi_dir / "ppi_pyrosetta_residues.csv"),
        "ppi_summary": _load_csv(ppi_dir / "ppi_pyrosetta_summary.csv"),
        "afm_residues": _load_csv(ppi_dir / "ppi_afm_residues.csv"),
        "afm_summary": _load_csv(ppi_dir / "ppi_afm_summary.csv"),
    }


# ---------------------------------------------------------------------------
# AlphaFold-Multimer field adaptation
# ---------------------------------------------------------------------------

def _adapt_afm_to_ppi_format(afm_rows: List[dict]) -> List[dict]:
    """Convert AFM residue rows to PPI-compatible format for merge.

    AFM provides min_ca_distance but no occupancy/frequency.
    Maps CA-CA distance to synthetic occupancy via sigmoid:
      4 A -> 0.95, 8 A -> 0.50, 12 A -> 0.15
    """
    adapted = []
    for row in afm_rows:
        dist = _safe_float(row.get("min_ca_distance"), 99.0)
        # Sigmoid: closer distance = higher synthetic occupancy
        synthetic_occ = 1.0 / (1.0 + math.exp((dist - 8.0) / 2.0))
        adapted.append({
            "receptor_id": row.get("receptor_id", ""),
            "source": f"afm:{row.get('receptor_id', '')}",
            "residue_id": row.get("residue_id", ""),
            "residue_num": row.get("residue_num", ""),
            "occupancy": round(synthetic_occ, 3),
            "frequency": 1.0,
            "min_ca_distance": dist,
        })
    return adapted


# ---------------------------------------------------------------------------
# PPI interface centroid computation
# ---------------------------------------------------------------------------


def _ppi_partner_name(row: dict) -> str:
    partner = str(row.get("partner_id", "")).strip()
    if partner:
        return partner
    source = str(row.get("source", "pyrosetta_ppi")).strip()
    if source.startswith("pyrosetta_ppi:"):
        return source.split(":", 1)[1]
    return "default"


def _ppi_row_priority(row: dict) -> Tuple[float, float, float]:
    return (
        _safe_float(row.get("frac_runs_supporting"), 0.0),
        _safe_float(row.get("occupancy"), 0.0),
        _safe_float(row.get("n_runs_supporting"), 0.0),
    )


def _best_summary_by_receptor(ppi_summary: List[dict]) -> Dict[str, dict]:
    best: Dict[str, dict] = {}
    for row in ppi_summary:
        rid = row.get("receptor_id", "")
        if not rid:
            continue
        current = best.get(rid)
        current_best = _safe_float(current.get("best_dg"), float("inf")) if current else float("inf")
        candidate_best = _safe_float(row.get("best_dg"), float("inf"))
        if current is None or candidate_best < current_best:
            best[rid] = row
    return best


def _enrich_agreement_with_ppi_reproducibility(
    agreement_rows: List[dict],
    ppi_residues: List[dict],
) -> None:
    """Attach reproducibility summaries for shared PPI residues.

    Step 3 now preserves run-level reproducibility in ``frac_runs_supporting``.
    Step 5 should not throw that signal away. This helper projects the most
    relevant residue-level reproducibility metrics onto each pocket-level
    agreement row.
    """
    residue_index: Dict[Tuple[str, str], dict] = {}
    partner_labels: Dict[Tuple[str, str], Set[str]] = defaultdict(set)
    for row in ppi_residues:
        key = (row.get("receptor_id", ""), row.get("residue_id", ""))
        if not key[0] or not key[1]:
            continue
        partner_labels[key].add(_ppi_partner_name(row))
        existing = residue_index.get(key)
        if existing is None or _ppi_row_priority(row) > _ppi_row_priority(existing):
            residue_index[key] = row

    for row in agreement_rows:
        rid = row.get("receptor_id", "")
        shared_ids = [
            item.strip()
            for item in str(row.get("shared_residue_list", "")).split(";")
            if item.strip()
        ]
        matched = [
            residue_index[(rid, res_id)]
            for res_id in shared_ids
            if (rid, res_id) in residue_index
        ]
        if not matched:
            row["ppi_frac_runs_supporting"] = ""
            row["ppi_best_frac_runs_supporting"] = ""
            row["ppi_best_interface_delta_e"] = ""
            row["ppi_shared_partner_count"] = ""
            continue

        mean_frac = sum(
            _safe_float(item.get("frac_runs_supporting"), 0.0)
            for item in matched
        ) / len(matched)
        best_frac = max(_safe_float(item.get("frac_runs_supporting"), 0.0) for item in matched)
        delta_values = [
            _safe_float(item.get("best_interface_delta_e"), 0.0)
            for item in matched
            if item.get("best_interface_delta_e") not in ("", None)
        ]
        partner_count = len(
            {
                partner
                for item in matched
                for partner in partner_labels.get(
                    (item.get("receptor_id", ""), item.get("residue_id", "")),
                    set(),
                )
            }
        )
        row["ppi_frac_runs_supporting"] = round(mean_frac, 4)
        row["ppi_best_frac_runs_supporting"] = round(best_frac, 4)
        row["ppi_best_interface_delta_e"] = (
            round(min(delta_values), 4) if delta_values else ""
        )
        row["ppi_shared_partner_count"] = partner_count

def _merge_multi_partner_residues(
    ppi_residues: List[dict],
) -> List[dict]:
    """Merge PPI residues from multiple partners into a unified set.

    When multiple partners (e.g., beta_meander + TH1) provide residue data
    for the same receptor, keeps the most reproducible entry per
    (receptor_id, residue_id) pair. Annotates merged source.

    Returns merged list (safe to use as drop-in replacement).
    """
    # Group by (receptor_id, residue_id), keep best reproducibility/occupancy.
    best: Dict[Tuple[str, str], dict] = {}
    sources: Dict[Tuple[str, str], Set[str]] = defaultdict(set)
    partners: Dict[Tuple[str, str], Set[str]] = defaultdict(set)

    for row in ppi_residues:
        key = (row["receptor_id"], row.get("residue_id", ""))
        source = row.get("source", "pyrosetta_ppi")
        sources[key].add(source)
        partners[key].add(_ppi_partner_name(row))
        existing = best.get(key)
        if existing is None or _ppi_row_priority(row) > _ppi_row_priority(existing):
            best[key] = dict(row)

    # Annotate merged source/partner label
    for key, row in best.items():
        src_set = sources[key]
        partner_names = sorted(partners[key])
        if len(partner_names) > 1:
            merged_partner = "+".join(partner_names)
            row["source"] = f"pyrosetta_ppi:merged({merged_partner})"
            row["partner_id"] = f"merged({merged_partner})"
        elif partner_names:
            row["partner_id"] = partner_names[0]
            if len(src_set) > 1:
                row["source"] = f"pyrosetta_ppi:{partner_names[0]}"

    return list(best.values())


def _build_ppi_partner_centroids(
    ppi_residues: List[dict],
    config: dict,
) -> Dict[str, List[Tuple[str, Optional[Tuple[float, float, float]], float]]]:
    """Compute per-partner PPI interface centroids for each receptor.

    Returns {receptor_id: [(partner_name, centroid, spread), ...]}.
    Each PPI partner gets its own centroid, enabling min-distance matching.
    """
    # Group PPI residue numbers by (receptor, partner)
    partner_resnums: Dict[Tuple[str, str], Set[int]] = defaultdict(set)
    for row in ppi_residues:
        rid = row["receptor_id"]
        partner = _ppi_partner_name(row)
        resnum = _safe_int(row.get("residue_num"), 0)
        if resnum > 0:
            partner_resnums[(rid, partner)].add(resnum)

    if not partner_resnums:
        return {}

    # Build receptor PDB path index
    receptor_pdbs: Dict[str, Path] = {}
    for rec in config.get("receptors", []):
        pdb_path = rec.get("pdb", "")
        if pdb_path:
            receptor_pdbs[rec["id"]] = Path(pdb_path)

    result: Dict[str, List[Tuple[str, Optional[Tuple[float, float, float]], float]]] = defaultdict(list)
    for (rid, partner), resnums in partner_resnums.items():
        pdb_path = receptor_pdbs.get(rid)
        if not pdb_path or not pdb_path.exists():
            result[rid].append((partner, None, 0.0))
            continue
        ca_coords = _parse_ca_coords(pdb_path)
        centroid, spread = _compute_centroid_and_spread(ca_coords, resnums)
        result[rid].append((partner, centroid, spread))

    return dict(result)


def _build_ppi_interface_centroids(
    ppi_residues: List[dict],
    config: dict,
) -> Tuple[Dict[str, Tuple[float, float, float]], Dict[str, float]]:
    """Compute merged PPI interface centroid per receptor (backward compatible).

    Uses merged residues from all partners. Returns {receptor_id: centroid}.
    """
    merged = _merge_multi_partner_residues(ppi_residues)
    ppi_resnums: Dict[str, Set[int]] = defaultdict(set)
    for row in merged:
        rid = row["receptor_id"]
        resnum = _safe_int(row.get("residue_num"), 0)
        if resnum > 0:
            ppi_resnums[rid].add(resnum)

    if not ppi_resnums:
        return {}, {}

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
        existing = index[rid].get(res_id)
        if existing is None or _ppi_row_priority(row) > _ppi_row_priority(existing):
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
    ppi_partner_centroids: Dict[str, List[Tuple[str, Optional[Tuple[float, float, float]], float]]],
    thresholds: Optional[dict] = None,
) -> List[dict]:
    """Compute Vina ↔ PPI agreement for each pocket.

    Uses per-partner centroids: for each pocket, finds the closest PPI
    partner centroid and reports the minimum distance. This is biologically
    meaningful: "Is this drug pocket near ANY known PPI interface?"

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
    ppi_receptors = set(ppi_index.keys())
    for rid in ppi_partner_centroids:
        ppi_receptors.add(rid)

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

        # Spatial proximity: find closest PPI partner centroid
        spatial_dist = None
        spatial_proximity = "no_data"
        closest_partner = ""
        n_partners_near = 0

        partner_entries = ppi_partner_centroids.get(receptor_id, [])
        if partner_entries:
            try:
                pocket_centroid = (
                    float(pocket.get("centroid_x", 0)),
                    float(pocket.get("centroid_y", 0)),
                    float(pocket.get("centroid_z", 0)),
                )
                best_dist = float("inf")
                for pname, pcentroid, _ in partner_entries:
                    if pcentroid is None:
                        continue
                    d = _euclidean_dist(pocket_centroid, pcentroid)
                    # Count partners within "moderate" distance
                    if d < T["ppi_dist_moderate"]:
                        n_partners_near += 1
                    if d < best_dist:
                        best_dist = d
                        closest_partner = pname
                spatial_dist = round(best_dist, 2) if best_dist < float("inf") else None

                if spatial_dist is not None:
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
            "closest_ppi_partner": closest_partner,
            "n_ppi_partners_near": n_partners_near,
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


def compute_cross_receptor_support(
    comparison_rows: List[dict],
) -> Dict[Tuple[str, str], dict]:
    """Summarize graded cross-state support for each pocket.

    ``same_patch_candidate`` is still the hard gate for consensus, but the
    downstream verdict should distinguish weak candidate matches from stronger
    ones using centroid, residue, and ligand agreement signals.
    """
    raw: Dict[Tuple[str, str], List[dict]] = defaultdict(list)
    for row in comparison_rows:
        is_candidate = str(row.get("same_patch_candidate", "")).lower()
        if is_candidate not in ("true", "1", "yes"):
            continue
        rec_a, rec_b = row["receptor_a"], row["receptor_b"]
        pkt_a, pkt_b = row["pocket_a"], row["pocket_b"]
        centroid_dist = _safe_float(row.get("centroid_dist"), 999.0)
        residue_jaccard = _safe_float(
            row.get("residue_jaccard", row.get("jaccard")),
            0.0,
        )
        overlap_coeff = _safe_float(
            row.get("residue_overlap_coeff", row.get("overlap_coeff")),
            0.0,
        )
        n_shared_ligands = _safe_int(row.get("n_shared_ligands"), 0)

        pair_strength = 0.40
        if centroid_dist <= 6.0:
            pair_strength += 0.20
        elif centroid_dist <= 8.0:
            pair_strength += 0.10

        if residue_jaccard >= 0.50:
            pair_strength += 0.20
        elif residue_jaccard >= 0.30:
            pair_strength += 0.10

        if overlap_coeff >= 0.60:
            pair_strength += 0.10
        elif overlap_coeff >= 0.40:
            pair_strength += 0.05

        if n_shared_ligands >= 2:
            pair_strength += 0.10
        elif n_shared_ligands >= 1:
            pair_strength += 0.05

        pair_strength = min(pair_strength, 1.0)
        entry = {
            "other_receptor": rec_b,
            "centroid_dist": centroid_dist,
            "residue_jaccard": residue_jaccard,
            "overlap_coeff": overlap_coeff,
            "n_shared_ligands": n_shared_ligands,
            "pair_strength": round(pair_strength, 4),
        }
        raw[(rec_a, pkt_a)].append(entry)
        raw[(rec_b, pkt_b)].append(
            dict(entry, other_receptor=rec_a)
        )

    summary: Dict[Tuple[str, str], dict] = {}
    for key, items in raw.items():
        if not items:
            continue
        strength = sum(item["pair_strength"] for item in items)
        support_frac = min(strength / 2.0, 1.0)
        summary[key] = {
            "match_count": len(items),
            "support_frac": round(support_frac, 4),
            "best_centroid_dist": round(min(item["centroid_dist"] for item in items), 4),
            "mean_residue_jaccard": round(
                sum(item["residue_jaccard"] for item in items) / len(items),
                4,
            ),
            "max_shared_ligands": max(item["n_shared_ligands"] for item in items),
            "matched_receptors": sorted(
                {item["other_receptor"] for item in items if item["other_receptor"]}
            ),
        }
    return summary


# ---------------------------------------------------------------------------
# Step 2.5: Consensus site identification
# ---------------------------------------------------------------------------

def identify_consensus_sites(
    comparison_rows: List[dict],
    pocket_rows: List[dict],
) -> Tuple[List[dict], Dict[Tuple[str, str], str]]:
    """Group pockets into consensus sites via transitive closure.

    A consensus site = group of pockets from different receptors that are
    same_patch_candidates. Pockets appearing in 2+ receptors form a
    consensus site (structurally stable across conformational states).

    Returns:
      - List of consensus site rows (for CSV)
      - Mapping {(receptor_id, pocket_id): consensus_site_id}
    """
    # Build adjacency from same_patch_candidate pairs
    edges: List[Tuple[Tuple[str, str], Tuple[str, str]]] = []
    for row in comparison_rows:
        is_candidate = str(row.get("same_patch_candidate", "")).lower()
        if is_candidate not in ("true", "1", "yes"):
            continue
        a = (row["receptor_a"], row["pocket_a"])
        b = (row["receptor_b"], row["pocket_b"])
        edges.append((a, b))

    if not edges:
        return [], {}

    # Union-Find for transitive closure
    parent: Dict[Tuple[str, str], Tuple[str, str]] = {}

    def find(x):
        while parent.get(x, x) != x:
            parent[x] = parent.get(parent[x], parent[x])
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for a, b in edges:
        parent.setdefault(a, a)
        parent.setdefault(b, b)
        union(a, b)

    # Group by root
    groups: Dict[Tuple[str, str], List[Tuple[str, str]]] = defaultdict(list)
    for node in parent:
        groups[find(node)].append(node)

    # Filter: consensus = 2+ different receptors
    pocket_index = {
        (p["receptor_id"], p["pocket_id"]): p for p in pocket_rows
    }

    consensus_rows: List[dict] = []
    pocket_to_cs: Dict[Tuple[str, str], str] = {}
    cs_id = 0

    for members in sorted(groups.values(), key=lambda m: len(set(r for r, _ in m)), reverse=True):
        receptors = set(r for r, _ in members)
        if len(receptors) < 2:
            continue
        cs_id += 1
        cs_name = f"CS{cs_id:03d}"

        # Aggregate metrics
        all_residues: Set[str] = set()
        best_aff = 0.0
        total_ligand = 0
        total_pose = 0
        centroids_x, centroids_y, centroids_z = [], [], []

        for key in members:
            pocket_to_cs[key] = cs_name
            p = pocket_index.get(key, {})
            residues = parse_residue_set(p.get("union_contact_residues", ""))
            all_residues |= residues
            aff = _safe_float(p.get("best_affinity"), 0.0)
            if aff < best_aff:
                best_aff = aff
            total_ligand += _safe_int(p.get("n_ligand"), 0)
            total_pose += _safe_int(p.get("n_pose"), 0)
            cx = _safe_float(p.get("centroid_x"), None)
            cy = _safe_float(p.get("centroid_y"), None)
            cz = _safe_float(p.get("centroid_z"), None)
            if cx is not None and cy is not None and cz is not None:
                centroids_x.append(cx)
                centroids_y.append(cy)
                centroids_z.append(cz)

        avg_cx = sum(centroids_x) / len(centroids_x) if centroids_x else ""
        avg_cy = sum(centroids_y) / len(centroids_y) if centroids_y else ""
        avg_cz = sum(centroids_z) / len(centroids_z) if centroids_z else ""

        consensus_rows.append({
            "consensus_site_id": cs_name,
            "n_receptors": len(receptors),
            "receptor_list": ";".join(sorted(receptors)),
            "pocket_list": ";".join(f"{r}:{p}" for r, p in sorted(members)),
            "centroid_x": round(avg_cx, 2) if isinstance(avg_cx, float) else "",
            "centroid_y": round(avg_cy, 2) if isinstance(avg_cy, float) else "",
            "centroid_z": round(avg_cz, 2) if isinstance(avg_cz, float) else "",
            "best_affinity": round(best_aff, 2) if best_aff < 0 else "",
            "total_n_ligand": total_ligand,
            "total_n_pose": total_pose,
            "consensus_residues": ";".join(sorted(all_residues)),
        })

    return consensus_rows, pocket_to_cs


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
# Experimental priors (known binding/non-binding residues)
# ---------------------------------------------------------------------------

def _parse_residue_ranges(raw: str) -> Set[int]:
    """Parse comma-separated residue numbers and ranges into a set of ints.

    Accepts: "744, 752, 831-859" → {744, 752, 831, 832, ..., 859}
    """
    result: Set[int] = set()
    if not raw:
        return result
    for token in str(raw).split(","):
        token = token.strip()
        if not token:
            continue
        if "-" in token:
            parts = token.split("-", 1)
            try:
                lo, hi = int(parts[0].strip()), int(parts[1].strip())
                result.update(range(lo, hi + 1))
            except ValueError:
                continue
        else:
            try:
                result.add(int(token))
            except ValueError:
                continue
    return result


def _pocket_resnums(pocket: dict) -> Set[int]:
    """Extract residue numbers from a pocket's union_contact_residues."""
    resnums: Set[int] = set()
    raw = pocket.get("union_contact_residues", "")
    for res in raw.split(";"):
        res = res.strip()
        if not res:
            continue
        m = re.search(r'(\d+)$', res)
        if m:
            resnums.add(int(m.group(1)))
    return resnums


def compute_experimental_correlation(
    pocket_rows: List[dict],
    known_binding: Set[int],
    known_non_binding: Set[int],
) -> Dict[Tuple[str, str], dict]:
    """Compute per-pocket correlation with experimental residue data.

    For each pocket, computes:
      - sensitivity: fraction of known binding residues contacted
      - specificity: fraction of known non-binding residues NOT contacted
      - enrichment: observed/expected ratio of binding residue hits
      - exp_hit_count: number of known binding residues in pocket contacts
      - exp_false_pos: number of known non-binding residues in pocket contacts

    Returns {(receptor_id, pocket_id): stats_dict}.
    Does NOT change the 100-point scoring — adds informational tags only.
    """
    if not known_binding and not known_non_binding:
        return {}

    n_known_binding = len(known_binding)
    n_known_non_binding = len(known_non_binding)
    total_known = n_known_binding + n_known_non_binding

    results: Dict[Tuple[str, str], dict] = {}

    for pocket in pocket_rows:
        rid = pocket["receptor_id"]
        pid = pocket["pocket_id"]
        pocket_nums = _pocket_resnums(pocket)
        n_pocket = len(pocket_nums)

        # Hits: known binding residues found in pocket contacts
        hits = pocket_nums & known_binding
        n_hits = len(hits)

        # False positives: known non-binding residues in pocket contacts
        false_pos = pocket_nums & known_non_binding
        n_false_pos = len(false_pos)

        # Sensitivity: how many known binders does this pocket capture?
        sensitivity = n_hits / n_known_binding if n_known_binding > 0 else 0.0

        # Specificity: how many known non-binders does this pocket avoid?
        specificity = (
            (n_known_non_binding - n_false_pos) / n_known_non_binding
            if n_known_non_binding > 0 else 1.0
        )

        # Enrichment: observed / expected ratio
        # Expected = n_pocket * (n_known_binding / total_residue_pool)
        # Use total_known as a proxy for the testable pool
        if total_known > 0 and n_pocket > 0:
            expected = n_pocket * (n_known_binding / total_known)
            enrichment = n_hits / expected if expected > 0 else 0.0
        else:
            enrichment = 0.0

        results[(rid, pid)] = {
            "exp_sensitivity": round(sensitivity, 4),
            "exp_specificity": round(specificity, 4),
            "exp_enrichment": round(enrichment, 2),
            "exp_hit_count": n_hits,
            "exp_false_pos": n_false_pos,
            "exp_hit_residues": ";".join(str(r) for r in sorted(hits)),
        }

    return results


# ---------------------------------------------------------------------------
# Step 3: Per-pocket scoring (adaptive)
# ---------------------------------------------------------------------------

def score_pocket(
    pocket: dict,
    ppi_agreement: Optional[dict],
    cross_receptor_matches: List[str],
    cross_support: Optional[dict],
    thresholds: dict,
    has_ppi_data: bool,
    exp_correlation: Optional[dict] = None,
) -> Tuple[float, str, List[str], float, float, float, dict]:
    """Score a single pocket with adaptive weighting.

    Scoring adapts to available evidence:
      With PPI data:    Vina(50) + PPI_proximity(20) + Cross_receptor(30) = 100
      Without PPI data: Vina(60) + Cross_receptor(40) = 100

    This ensures pockets are never penalized for missing PPI data.

    exp_correlation: optional dict from compute_experimental_correlation().
      Does NOT change the 100-point total — adds informational reason tags only.

    Returns:
      (total, verdict, reasons, vina_score, ppi_score, cross_score, raw_components)
      where raw_components is a dict of pre-normalization point values.
    """
    reasons: List[str] = []
    reason_tags: List[str] = []
    T = thresholds

    # ---- Axis 1: Vina Quality ----
    affinity = _safe_float(pocket.get("best_affinity"), 0.0)
    n_pose = _safe_int(pocket.get("n_pose"), 0)
    n_ligand = _safe_int(pocket.get("n_ligand"), 0)
    pocket_stability = _safe_float(pocket.get("pocket_exists_frac"), 0.0)
    dominant_ligand_fraction = _safe_float(
        pocket.get("dominant_ligand_fraction"),
        1.0 if n_ligand <= 1 else 0.0,
    )
    ligand_pose_entropy = _safe_float(pocket.get("ligand_pose_entropy"), 0.0)

    vina_raw = 0.0
    affinity_pts = 0.0
    convergence_pts = 0.0
    stability_pts = 0.0
    diversity_pts = 0.0

    # Affinity: graduated scoring (0 / 10 / 20)
    if affinity <= T["affinity_great"]:
        affinity_pts = 20.0
        vina_raw += 20
        reasons.append(f"affinity={affinity:.1f}")
        reason_tags.append("affinity_strong")
    elif affinity <= T["affinity_good"]:
        affinity_pts = 10.0
        vina_raw += 10
        reasons.append(f"affinity={affinity:.1f}")
        reason_tags.append("affinity_good")

    # Pose convergence is informative, but must be discounted when the pocket
    # is unstable under pose resampling or dominated by a single ligand.
    base_convergence = 0.0
    if n_pose >= T["n_pose_great"]:
        base_convergence = 15.0
    elif n_pose >= T["n_pose_good"]:
        base_convergence = 8.0

    convergence_factor = 1.0
    if pocket_stability > 0 and pocket_stability < T["stability_good"]:
        convergence_factor *= 0.5
    if n_ligand <= 1 and dominant_ligand_fraction >= T["dominant_ligand_warn"]:
        convergence_factor *= 0.75
    convergence_pts = round(base_convergence * convergence_factor, 2)
    if convergence_pts > 0:
        vina_raw += convergence_pts
        reasons.append(f"n_pose={n_pose}")
        reason_tags.append("pose_convergence")

    # Stability from pose-resampling bootstrap remains distinct from raw pose
    # abundance so that blind-docking sampling bias does not masquerade as
    # reproducibility.
    if pocket_stability >= T["stability_great"]:
        stability_pts = 10.0
        reasons.append(f"stability={pocket_stability:.2f}")
        reason_tags.append("vina_stable")
    elif pocket_stability >= T["stability_good"]:
        stability_pts = 6.0
        reasons.append(f"stability={pocket_stability:.2f}")
        reason_tags.append("vina_stable")
    elif pocket_stability >= T["stability_min"]:
        stability_pts = 3.0
        reasons.append(f"stability={pocket_stability:.2f}")
        reason_tags.append("vina_borderline_stable")
    elif pocket_stability > 0:
        reasons.append(f"low_stability={pocket_stability:.2f}")
        reason_tags.append("vina_low_stability")
    vina_raw += stability_pts

    # Multi-ligand support is scored separately from raw pose abundance.
    if n_ligand >= T["n_ligand_all"]:
        diversity_pts += 8.0
    elif n_ligand >= T["n_ligand_good"]:
        diversity_pts += 4.0

    if n_ligand >= T["n_ligand_good"]:
        if dominant_ligand_fraction <= T["dominant_ligand_good"]:
            diversity_pts += 4.0
        elif dominant_ligand_fraction <= T["dominant_ligand_warn"]:
            diversity_pts += 2.0
        else:
            reasons.append(f"ligand_dominance={dominant_ligand_fraction:.2f}")
            reason_tags.append("ligand_dominated")

        if ligand_pose_entropy >= T["ligand_entropy_good"]:
            diversity_pts += 3.0
        elif ligand_pose_entropy >= T["ligand_entropy_ok"]:
            diversity_pts += 1.5

    diversity_pts = min(diversity_pts, 15.0)
    if diversity_pts > 0:
        vina_raw += diversity_pts
        reasons.append(f"n_ligand={n_ligand}")
        reason_tags.append("multi_ligand")

    # Membrane-face tag (informational, not scored)
    membrane_tag = _check_membrane_overlap(pocket)
    if membrane_tag:
        reasons.append(membrane_tag)
        reason_tags.append("membrane_face")

    # Normalize to adaptive max — config에서 조정 가능 (sensitivity analysis)
    _AFFINITY_MAX = T["vina_affinity_max"]
    _CONVERGENCE_MAX = T["vina_convergence_max"]
    _STABILITY_MAX = T["vina_stability_max"]
    _DIVERSITY_MAX = T["vina_diversity_max"]
    VINA_RAW_MAX = _AFFINITY_MAX + _CONVERGENCE_MAX + _STABILITY_MAX + _DIVERSITY_MAX
    vina_max = 50.0 if has_ppi_data else 60.0
    vina_score = min(vina_raw, VINA_RAW_MAX) / VINA_RAW_MAX * vina_max

    # ---- Axis 2: PPI Spatial Proximity ----
    # Note: spatial_dist is centroid-to-centroid. Vina centroid sits inside
    # the pocket cavity (~3-5 Å below surface), so actual pocket-entrance-to-
    # PPI-surface distance is shorter. Thresholds are set conservatively to
    # account for this systematic overestimate.
    ppi_score = 0.0
    ppi_max = 20.0 if has_ppi_data else 0.0
    spatial_pts = 0.0
    overlap_pts = 0.0
    reproducibility_pts = 0.0

    if has_ppi_data and ppi_agreement:
        spatial = ppi_agreement.get("spatial_proximity", "no_data")
        n_shared = _safe_int(ppi_agreement.get("n_shared_residues"), 0)
        spatial_dist = ppi_agreement.get("spatial_dist_A", "")
        ppi_frac_runs = _safe_float(ppi_agreement.get("ppi_frac_runs_supporting"), 0.0)
        ppi_best_delta_e = _safe_float(
            ppi_agreement.get("ppi_best_interface_delta_e"),
            0.0,
        )

        if spatial == "adjacent":
            spatial_pts = 10.0
            reasons.append(f"ppi_adjacent({spatial_dist}A)")
            reason_tags.append("ppi_adjacent")
        elif spatial == "near":
            spatial_pts = 6.0
            reasons.append(f"ppi_near({spatial_dist}A)")
            reason_tags.append("ppi_near")
        elif spatial == "moderate":
            spatial_pts = 2.0
            reasons.append(f"ppi_moderate({spatial_dist}A)")
        elif spatial == "no_data":
            # PPI residues exist but no PDB for centroid → use residue overlap
            if n_shared > 0:
                spatial_pts = 5.0
                reasons.append(f"ppi_shared={n_shared}res")
                reason_tags.append("ppi_shared_only")
        else:
            reasons.append(f"ppi_distant({spatial_dist}A)")
            reason_tags.append("ppi_distant")

        # Residue overlap bonus (weak but informative)
        if T.get("ppi_residue_bonus") and n_shared > 0:
            occ = _safe_float(ppi_agreement.get("ppi_mean_occupancy_of_shared"), 0)
            if occ >= 0.5:
                overlap_pts = 4.0
                reasons.append(f"shared_highocc={n_shared}")
                reason_tags.append("ppi_shared_highocc")
            else:
                overlap_pts = 2.0
                reason_tags.append("ppi_shared")

        # Multi-partner corroboration bonus: near both beta_meander AND TH1
        n_partners_near = _safe_int(ppi_agreement.get("n_ppi_partners_near"), 0)
        if n_partners_near >= 2:
            spatial_pts += 2.0
            reasons.append(f"multi_ppi={n_partners_near}partners")
            reason_tags.append("ppi_multi_partner")

        if ppi_frac_runs >= T["ppi_repro_great"]:
            reproducibility_pts = 4.0
            reasons.append(f"ppi_repro={ppi_frac_runs:.2f}")
            reason_tags.append("ppi_multi_seed_strong")
        elif ppi_frac_runs >= T["ppi_repro_good"]:
            reproducibility_pts = 2.0
            reasons.append(f"ppi_repro={ppi_frac_runs:.2f}")
            reason_tags.append("ppi_multi_seed")
        elif ppi_frac_runs > 0:
            reasons.append(f"ppi_repro={ppi_frac_runs:.2f}")

        if ppi_best_delta_e < 0:
            reasons.append(f"ppi_dE={ppi_best_delta_e:.1f}")

        ppi_score = min(spatial_pts + overlap_pts + reproducibility_pts, ppi_max)
    elif not has_ppi_data:
        reasons.append("no_ppi_data")
        reason_tags.append("no_ppi_data")

    # ---- Axis 3: Cross-Receptor Consistency ----
    # 2-tier 스코어링:
    #   (1) Coverage (n_cross): 몇 개 receptor state에서 동일 포켓이 발견되었는가
    #       AC-4.3 차등: 3/3 → 20 pts, 2/3 → 14 pts, 1/3 → 0 pts (configurable)
    #   (2) Support quality (support_frac): 매칭 품질 (centroid 근접도 + 잔기
    #       Jaccard/overlap + 공유 리간드 수로 산출된 pair_strength의 합, 0-1 정규화)
    #       - great(≥0.75) → 10 pts, good(≥0.50) → 6 pts, >0 → 3 pts
    #       (PPI 데이터 없으면 각 2배)
    cross_max = 30.0 if has_ppi_data else 40.0
    n_cross = len(cross_receptor_matches)
    cross_coverage_pts = 0.0
    cross_support_pts = 0.0
    support_frac = _safe_float(
        (cross_support or {}).get("support_frac"),
        0.0,
    )

    # AC-4.3: Differentiated coverage scoring (3/3 > 2/3 > 1/3)
    if n_cross >= 2:
        cross_coverage_pts = T["cross_coverage_3of3"]
        reasons.append(f"cross={n_cross + 1}/3")
        reason_tags.append("cross_state_recurrent")
    elif n_cross >= 1:
        cross_coverage_pts = T["cross_coverage_2of3"]
        reasons.append("cross=2/3")
        reason_tags.append("cross_state_partial")

    if support_frac >= T["cross_support_great"]:
        cross_support_pts = 10.0 if has_ppi_data else 20.0
        reasons.append(f"cross_support={support_frac:.2f}")
    elif support_frac >= T["cross_support_good"]:
        cross_support_pts = 6.0 if has_ppi_data else 12.0
        reasons.append(f"cross_support={support_frac:.2f}")
    elif support_frac > 0:
        cross_support_pts = 3.0 if has_ppi_data else 6.0
        reasons.append(f"cross_support={support_frac:.2f}")

    cross_score = min(cross_coverage_pts + cross_support_pts, cross_max)

    # ---- Experimental priors (informational only, no score impact) ----
    if exp_correlation:
        sens = exp_correlation.get("exp_sensitivity", 0)
        n_hits = exp_correlation.get("exp_hit_count", 0)
        n_fp = exp_correlation.get("exp_false_pos", 0)
        enrichment = exp_correlation.get("exp_enrichment", 0)
        if n_hits > 0:
            reasons.append(f"exp_hit={n_hits}res(sens={sens:.0%})")
            reason_tags.append("experimental_hit")
        if enrichment > 2.0:
            reasons.append(f"exp_enriched({enrichment:.1f}x)")
        if n_fp > 0:
            reasons.append(f"exp_fp={n_fp}res")

    # ---- Total (always out of 100) ----
    total = vina_score + ppi_score + cross_score

    if total >= T["valid_min"]:
        verdict = "STRONG"
    elif total >= T["uncertain_min"]:
        verdict = "MODERATE"
    else:
        verdict = "WEAK"

    evidence_profile: List[str] = ["exploratory"]
    if pocket_stability >= T["stability_good"]:
        evidence_profile.append("stable")
    if n_ligand >= T["n_ligand_good"]:
        evidence_profile.append("multi_ligand")
    if spatial_pts >= 6.0 or overlap_pts > 0 or reproducibility_pts > 0:
        evidence_profile.append("ppi_supported")
    if n_cross > 0:
        evidence_profile.append("cross_state")
    if support_frac >= T["cross_support_good"]:
        evidence_profile.append("recurrent")
    evidence_profile = sorted(set(evidence_profile))
    decision_trace = (
        "VINA[aff={aff:.1f},conv={conv:.1f},stab={stab:.1f},div={div:.1f}] "
        "PPI[spatial={spatial:.1f},overlap={overlap:.1f},repro={repro:.1f}] "
        "CROSS[count={count},support={support:.1f}]"
    ).format(
        aff=affinity_pts,
        conv=convergence_pts,
        stab=stability_pts,
        div=diversity_pts,
        spatial=spatial_pts,
        overlap=overlap_pts,
        repro=reproducibility_pts,
        count=cross_coverage_pts,
        support=cross_support_pts,
    )

    # Raw component breakdown (pre-normalization points)
    score_denominator = vina_max + ppi_max + cross_max
    raw_components = {
        "vina_affinity_pts": affinity_pts,
        "vina_convergence_pts": convergence_pts,
        "vina_stability_pts": stability_pts,
        "vina_diversity_pts": diversity_pts,
        # consensus_pts는 diversity_pts의 legacy alias (동일 값).
        # Multi-ligand support 점수를 그대로 반영하며, 별도 consensus 로직 없음.
        "vina_consensus_pts": diversity_pts,
        "ppi_spatial_pts": spatial_pts,
        "ppi_overlap_pts": overlap_pts,
        "ppi_reproducibility_pts": reproducibility_pts,
        "cross_receptor_pts": cross_coverage_pts + cross_support_pts,
        "cross_receptor_support_pts": cross_support_pts,
        "score_denominator": score_denominator,
        "dominant_ligand_fraction": dominant_ligand_fraction,
        "ligand_pose_entropy": ligand_pose_entropy,
        "ppi_frac_runs_supporting": (
            _safe_float(ppi_agreement.get("ppi_frac_runs_supporting"), 0.0)
            if ppi_agreement else 0.0
        ),
        "ppi_best_interface_delta_e": (
            _safe_float(ppi_agreement.get("ppi_best_interface_delta_e"), 0.0)
            if ppi_agreement else 0.0
        ),
        "cross_receptor_support": support_frac,
        "evidence_profile": "+".join(evidence_profile),
        "reason_tags": ";".join(sorted(set(reason_tags))),
        "decision_trace": decision_trace,
    }

    return total, verdict, reasons, vina_score, ppi_score, cross_score, raw_components


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
    vina_post = paths.wa_phase4_vina_postprocess(config)
    ppi_post = paths.wa_phase3_ppi_postprocess(config)
    out_dir = Path(output_dir) if output_dir else paths.wa_phase5_verdict(config)
    thresholds = _get_thresholds(config)

    # Load evidence
    evidence = load_all_evidence(vina_post, ppi_post)
    pocket_rows = evidence["pocket_table"]
    if not pocket_rows:
        print("[verdict] No pocket table found — skipping.")
        return Path(), Path()

    ppi_residues = evidence["ppi_residues"]
    afm_residues = evidence.get("afm_residues", [])

    # Merge AFM data into PPI residues (adapted to common format)
    if afm_residues:
        adapted_afm = _adapt_afm_to_ppi_format(afm_residues)
        ppi_residues = ppi_residues + adapted_afm
        n_afm_receptors = len({r.get("receptor_id") for r in afm_residues})
        print(f"[verdict] AFM data merged: {len(afm_residues)} residues "
              f"from {n_afm_receptors} receptor(s)")

    # Check for offset residues (unrestored chain B)
    offset_warnings = check_ppi_residue_offsets(ppi_residues)
    if offset_warnings:
        print("\n[verdict] WARNING: Unrestored chain B offset residues detected!")
        for w in offset_warnings[:5]:
            print(w)
        if len(offset_warnings) > 5:
            print(f"  ... and {len(offset_warnings) - 5} more")
        print("  Run PPI Postprocess (option 8) to fix before verdict.\n")

    # Merge multi-partner PPI residues
    ppi_residues_merged = _merge_multi_partner_residues(ppi_residues)

    # Determine which receptors have PPI data
    ppi_receptor_ids = {r["receptor_id"] for r in ppi_residues_merged}
    has_any_ppi = len(ppi_receptor_ids) > 0

    # Compute per-partner PPI centroids (for min-distance matching)
    ppi_partner_centroids = _build_ppi_partner_centroids(ppi_residues, config)
    # Also compute merged centroids for backward-compatible display
    ppi_centroids, ppi_spreads = _build_ppi_interface_centroids(ppi_residues, config)

    if has_any_ppi:
        print(f"[verdict] PPI data available for: {', '.join(sorted(ppi_receptor_ids))}")
        for rid, partners in sorted(ppi_partner_centroids.items()):
            for pname, centroid, spread in partners:
                if centroid:
                    spread_warn = " (DISPERSED)" if spread > 15.0 else ""
                    print(f"  PPI centroid ({rid}/{pname}): "
                          f"({centroid[0]:.1f}, {centroid[1]:.1f}, {centroid[2]:.1f})  "
                          f"spread={spread:.1f}A{spread_warn}")
                else:
                    print(f"  PPI centroid ({rid}/{pname}): no PDB coords")
        if not ppi_partner_centroids:
            print("  (No receptor PDBs found for centroid — using residue overlap only)")

    # Step 1: Cross-method agreement (uses per-partner centroids)
    agreement_rows = compute_cross_method_agreement(
        pocket_rows, ppi_residues_merged, ppi_partner_centroids, thresholds,
    )
    _enrich_agreement_with_ppi_reproducibility(agreement_rows, ppi_residues_merged)

    # Enrich with PPI summary best_dg (receptor-level, not pocket-specific;
    # same value for all pockets of a receptor — use for reference only)
    ppi_summary_map = _best_summary_by_receptor(evidence["ppi_summary"])
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
    cross_support = compute_cross_receptor_support(
        evidence["pocket_comparison"]
    )

    # Step 2.5: Consensus site identification
    consensus_rows, pocket_to_cs = identify_consensus_sites(
        evidence["pocket_comparison"], pocket_rows,
    )
    if consensus_rows:
        consensus_csv = _write_csv(
            out_dir / "vina_consensus_sites.csv",
            consensus_rows,
            CONSENSUS_FIELDS,
        )
        print(f"[verdict] {len(consensus_rows)} consensus sites identified → {consensus_csv}")
    else:
        print("[verdict] No consensus sites found (need same_patch_candidate across receptors)")

    # Step 2.7: Experimental priors (known binding/non-binding residues)
    exp_config = config.get("experimental") or {}
    known_binding_raw = exp_config.get("known_binding_residues", "")
    known_non_binding_raw = exp_config.get("known_non_binding_residues", "")
    known_binding = _parse_residue_ranges(
        ",".join(str(r) for r in known_binding_raw)
        if isinstance(known_binding_raw, list) else str(known_binding_raw)
    )
    known_non_binding = _parse_residue_ranges(
        ",".join(str(r) for r in known_non_binding_raw)
        if isinstance(known_non_binding_raw, list) else str(known_non_binding_raw)
    )

    exp_correlations: Dict[Tuple[str, str], dict] = {}
    if known_binding or known_non_binding:
        exp_source = exp_config.get("source", "unknown")
        print(f"[verdict] Experimental priors: {len(known_binding)} binding, "
              f"{len(known_non_binding)} non-binding residues (source: {exp_source})")
        exp_correlations = compute_experimental_correlation(
            pocket_rows, known_binding, known_non_binding,
        )

    # Step 2.8: Bootstrap stability data (optional)
    bootstrap_path = vina_post / "vina_pocket_bootstrap.csv"
    bootstrap_index: Dict[Tuple[str, str], dict] = {}
    if bootstrap_path.exists():
        bootstrap_rows = _load_csv(bootstrap_path)
        for brow in bootstrap_rows:
            bkey = (brow["receptor_id"], brow["pocket_id"])
            bootstrap_index[bkey] = brow
        print(f"[verdict] Bootstrap data loaded: {len(bootstrap_rows)} pocket entries")

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
        if key in bootstrap_index:
            pocket["pocket_exists_frac"] = bootstrap_index[key].get("pocket_exists_frac", "")

        ppi_agr = agreement_index.get(key)
        cross_matches = cross_receptor.get(key, [])

        # Per-pocket PPI availability check
        pocket_has_ppi = rid in ppi_receptor_ids

        exp_corr = exp_correlations.get(key)
        total, verdict, reasons, v_score, p_score, c_score, raw_comp = score_pocket(
            pocket, ppi_agr, cross_matches, cross_support.get(key), thresholds, pocket_has_ppi,
            exp_correlation=exp_corr,
        )

        spatial_dist = ""
        closest_partner = ""
        n_partners_near = 0
        if ppi_agr:
            spatial_dist = ppi_agr.get("spatial_dist_A", "")
            closest_partner = ppi_agr.get("closest_ppi_partner", "")
            n_partners_near = _safe_int(ppi_agr.get("n_ppi_partners_near"), 0)

        cs_id = pocket_to_cs.get(key, "")

        verdict_rows.append({
            "receptor_id": rid,
            "pocket_id": pid,
            "verdict": verdict,
            "confidence_score": round(total, 1),
            "vina_quality_score": round(v_score, 1),
            "ppi_proximity_score": round(p_score, 1),
            "cross_receptor_score": round(c_score, 1),
            "vina_affinity_pts": raw_comp["vina_affinity_pts"],
            "vina_convergence_pts": raw_comp["vina_convergence_pts"],
            "vina_stability_pts": raw_comp["vina_stability_pts"],
            "vina_diversity_pts": raw_comp["vina_diversity_pts"],
            "vina_consensus_pts": raw_comp["vina_consensus_pts"],
            "ppi_spatial_pts": raw_comp["ppi_spatial_pts"],
            "ppi_overlap_pts": raw_comp["ppi_overlap_pts"],
            "ppi_reproducibility_pts": raw_comp["ppi_reproducibility_pts"],
            "cross_receptor_pts": raw_comp["cross_receptor_pts"],
            "cross_receptor_support_pts": raw_comp["cross_receptor_support_pts"],
            "score_denominator": raw_comp["score_denominator"],
            "ppi_data_available": "yes" if pocket_has_ppi else "no",
            "best_affinity": pocket.get("best_affinity", ""),
            "n_pose": pocket.get("n_pose", ""),
            "n_ligand": pocket.get("n_ligand", ""),
            "dominant_ligand_fraction": (
                raw_comp["dominant_ligand_fraction"]
                if raw_comp["dominant_ligand_fraction"] not in (0.0, 1.0) or _safe_int(pocket.get("n_ligand"), 0) > 1
                else pocket.get("dominant_ligand_fraction", "")
            ),
            "ligand_pose_entropy": raw_comp["ligand_pose_entropy"],
            "spatial_dist_to_ppi": spatial_dist,
            "closest_ppi_partner": closest_partner,
            "n_ppi_partners_near": n_partners_near,
            "n_shared_with_ppi": ppi_agr.get("n_shared_residues", 0) if ppi_agr else 0,
            "ppi_frac_runs_supporting": raw_comp["ppi_frac_runs_supporting"] or "",
            "ppi_best_interface_delta_e": raw_comp["ppi_best_interface_delta_e"] or "",
            "cross_receptor_matches": ";".join(cross_matches),
            "cross_receptor_support": raw_comp["cross_receptor_support"] or "",
            "consensus_site_id": cs_id,
            "exp_sensitivity": exp_corr["exp_sensitivity"] if exp_corr else "",
            "exp_specificity": exp_corr["exp_specificity"] if exp_corr else "",
            "exp_enrichment": exp_corr["exp_enrichment"] if exp_corr else "",
            "exp_rank_impact": _compute_exp_rank_impact(exp_corr) if exp_corr else "",
            "pocket_stability": bootstrap_index[key]["pocket_exists_frac"]
                if key in bootstrap_index else "",
            "bootstrap_confidence": _bootstrap_confidence(
                bootstrap_index.get(key, {}).get("pocket_exists_frac", "")
            ),
            "evidence_profile": raw_comp["evidence_profile"],
            "reason_tags": raw_comp["reason_tags"],
            "decision_trace": raw_comp["decision_trace"],
            "reasons": "; ".join(reasons),
            "is_atp_site": pocket.get("is_atp_site", False),
            "exclusion_reason": "ATP_site_experimental"
                if pocket.get("is_atp_site") in (True, "True", "true")
                else "",
            # AC-5.3: Allosteric candidate — high Vina quality but PPI-distant
            "allosteric_candidate": (
                round(v_score, 1) >= 35.0 and round(p_score, 1) <= 5.0
            ),
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

    # AC-5.5: Write disclaimer metadata file (EC-5.3: separate file, not CSV comment)
    disclaimer_path = out_dir / "valid_sites_disclaimer.md"
    disclaimer_path.write_text(
        "# Disclaimer — valid_sites.csv\n\n"
        "This file contains computational predictions generated by automated\n"
        "molecular docking (AutoDock Vina) and PPI analysis (PyRosetta/LightDock)\n"
        "pipelines. All binding site predictions, affinity estimates, and verdict\n"
        "classifications (STRONG/MODERATE/WEAK) are hypotheses requiring\n"
        "experimental validation.\n\n"
        "**Do not interpret these results as confirmed biological effects.**\n\n"
        "Known limitations:\n"
        "- Rigid-body docking (induced fit not modeled)\n"
        "- Vina scoring function bias (hydrophobic overestimation)\n"
        "- Shared input structures across methods (common blind spots)\n"
        "- Implicit solvation only (water-mediated contacts missed)\n\n"
        "See `docs/methodology_limitations.md` for full details.\n",
        encoding="utf-8",
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


def _compute_exp_rank_impact(exp_corr: Optional[dict]) -> str:
    """Qualitative rank impact tag from experimental correlation.

    Helps the researcher see at a glance whether experimental data
    supports, contradicts, or is neutral for this pocket.
    """
    if not exp_corr:
        return ""
    sens = exp_corr.get("exp_sensitivity", 0)
    enrichment = exp_corr.get("exp_enrichment", 0)
    n_fp = exp_corr.get("exp_false_pos", 0)

    if sens >= 0.3 and enrichment >= 2.0 and n_fp == 0:
        return "supports"
    elif sens >= 0.1 and enrichment >= 1.5:
        return "consistent"
    elif n_fp > 0 and sens == 0:
        return "contradicts"
    else:
        return "neutral"


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
