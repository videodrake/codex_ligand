#!/usr/bin/env python3
"""Phase 1 Task 1.2A: Orientation-aware filtering for beta-meander docking poses.

Determines whether the active face (sheets 8 & 9) of the MYO1D beta-meander
is oriented toward the EGFR receptor surface (correct) or flipped away (artifact).

Algorithm — Dual-Vector Orientation Test:
  1. Collect Cα coordinates for sheets 8+9 (active face)
  2. Compute sheet-plane normal via PCA (smallest variance direction)
  3. Orient normal toward active face using multi-probe Cα→Cβ consensus
  4. Compute receptor-direction vector from sheet centroid to local receptor centroid
  5. Dot product determines orientation: positive=PASS, negative=FAIL, near-zero=AMBIGUOUS

Integration with Phase 1 pipeline:
  - Reads raw PyRosetta run outputs from
    output/workflow_a/phase2_ppi_docking/{state}/{run}/docking_*/final_result/
  - Produces orientation_filter_log.csv per state under
    output/workflow_b/phase1_ppi_analysis/{state}/
  - Merges orientation columns into pyrosetta_interface_models.csv
  - Retroactive pilot validation (server-side only, requires pilot PDB files)

Dependencies:
  - PyRosetta 4 (for pose loading and residue access)
  - NumPy (for linear algebra)

Usage:
    python -m egfr_pipeline.phase1.orientation_filter [--state 3GT8_raw]
    python -m egfr_pipeline.phase1.orientation_filter --pilot_dir /path/to/pilot/pdbs
"""

import argparse
import csv
import json
import os
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
from egfr_pipeline import paths

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent

# ---------------------------------------------------------------------------
# Constants: beta-meander sheet definitions (MYO1D residue numbering)
# ---------------------------------------------------------------------------

# Sheet boundaries — residue numbers in MYO1D (no offset)
SHEET_8_RESIDUES = [961, 962, 963, 964]
SHEET_9_RESIDUES = [968, 969, 970, 971, 972]
SHEET_10_RESIDUES = [977, 978, 979, 980]       # Approximate — verify from structure
SHEET_11_RESIDUES = [985, 986, 987, 988]       # Approximate — verify from structure
SHEET_12_RESIDUES = [993, 994, 995, 996, 997]  # Approximate — verify from structure

# Active face = sheets 8 + 9 (experimentally validated: Ko et al.)
ACTIVE_FACE_RESIDUES = SHEET_8_RESIDUES + SHEET_9_RESIDUES

# Back face = sheets 10, 11, 12 (structural scaffold)
BACK_FACE_RESIDUES = SHEET_10_RESIDUES + SHEET_11_RESIDUES + SHEET_12_RESIDUES

# Multi-probe consensus: span both sheets 8 and 9
PROBE_RESIDUES = [962, 964, 971]  # VAL962, VAL964 (sheet 8), SER971 (sheet 9)
PROBE_RESIDUE_FALLBACK = 962      # Single-probe fallback

# Orientation thresholds
AMBIGUOUS_BAND = 0.15  # |dot product| < this → ambiguous (edge-on)

# Phase 1 output structure
DEFAULT_PATH_CONFIG = {"output_root": str(PROJECT_ROOT / "output")}
PHASE1_RUNS_DIR = paths.wa_phase2_ppi_docking(DEFAULT_PATH_CONFIG)
PHASE1_OUTPUT_DIR = paths.wb_phase1_ppi_analysis(DEFAULT_PATH_CONFIG)
RECEPTOR_STATES = ["3GT8_raw", "EGFR_160-185", "EGFR_170-200"]

# Output schema
ORIENTATION_LOG_COLUMNS = [
    "model_id",
    "receptor_id",
    "seed_index",
    "orientation_score",
    "orientation_class",    # pass / fail / ambiguous / error codes
    "n_active_face_ca",
    "n_back_face_ca",
    "n_receptor_contact_ca",
    "sheet_centroid",       # [x, y, z]
    "normal_vector",        # [x, y, z]
    "receptor_direction",   # [x, y, z]
    "error",
    "source_pdb",
]


# ---------------------------------------------------------------------------
# Geometry utilities (adapted from docs/orientation_filter.py)
# ---------------------------------------------------------------------------

def get_ca_coord(pose, chain: str, resnum: int) -> Optional[np.ndarray]:
    """Extract Cα coordinate for a given chain and residue number."""
    for i in range(1, pose.total_residue() + 1):
        pdb_info = pose.pdb_info()
        if pdb_info.chain(i) == chain and pdb_info.number(i) == resnum:
            ca_idx = pose.residue(i).atom_index("CA")
            xyz = pose.residue(i).xyz(ca_idx)
            return np.array([xyz.x, xyz.y, xyz.z])
    return None


def get_cb_coord(pose, chain: str, resnum: int) -> Optional[np.ndarray]:
    """Extract Cβ coordinate (or Cα for GLY) for a given chain and residue."""
    for i in range(1, pose.total_residue() + 1):
        pdb_info = pose.pdb_info()
        if pdb_info.chain(i) == chain and pdb_info.number(i) == resnum:
            res = pose.residue(i)
            if res.has("CB"):
                cb_idx = res.atom_index("CB")
                xyz = res.xyz(cb_idx)
                return np.array([xyz.x, xyz.y, xyz.z])
            else:
                ca_idx = res.atom_index("CA")
                xyz = res.xyz(ca_idx)
                return np.array([xyz.x, xyz.y, xyz.z])
    return None


def collect_ca_coords(pose, chain: str, resnums: List[int]) -> np.ndarray:
    """Collect Cα coordinates for a list of residue numbers. Skip missing."""
    coords = []
    for rn in resnums:
        c = get_ca_coord(pose, chain, rn)
        if c is not None:
            coords.append(c)
    return np.array(coords) if coords else np.array([]).reshape(0, 3)


def compute_sheet_plane_normal(coords: np.ndarray) -> np.ndarray:
    """Compute sheet-plane normal via PCA (smallest variance component)."""
    if len(coords) < 3:
        raise ValueError(f"Need ≥3 Cα coordinates to define a plane, got {len(coords)}")

    centroid = coords.mean(axis=0)
    centered = coords - centroid
    _, s, vt = np.linalg.svd(centered, full_matrices=False)
    normal = vt[-1]
    return normal / np.linalg.norm(normal)


def orient_normal_to_active_face(
    normal: np.ndarray,
    pose,
    partner_chain: str,
    probe_resnums: List[int] = None,
) -> np.ndarray:
    """Orient the normal toward the active face using multi-probe Cα→Cβ consensus."""
    if probe_resnums is None:
        probe_resnums = PROBE_RESIDUES

    votes = 0
    n_valid = 0

    for rn in probe_resnums:
        ca = get_ca_coord(pose, partner_chain, rn)
        cb = get_cb_coord(pose, partner_chain, rn)
        if ca is None or cb is None:
            continue
        ca_cb = cb - ca
        norm_cacb = np.linalg.norm(ca_cb)
        if norm_cacb < 1e-6:
            continue
        ca_cb = ca_cb / norm_cacb
        dot = np.dot(normal, ca_cb)
        votes += 1 if dot > 0 else -1
        n_valid += 1

    if n_valid == 0:
        # Last resort: single-probe fallback
        ca = get_ca_coord(pose, partner_chain, PROBE_RESIDUE_FALLBACK)
        cb = get_cb_coord(pose, partner_chain, PROBE_RESIDUE_FALLBACK)
        if ca is None or cb is None:
            raise ValueError(f"No valid probe residues found for chain {partner_chain}")
        ca_cb = (cb - ca)
        ca_cb = ca_cb / np.linalg.norm(ca_cb)
        if np.dot(normal, ca_cb) < 0:
            normal = -normal
        return normal

    if votes < 0:
        normal = -normal
    return normal


def compute_receptor_interface_centroid(
    pose,
    receptor_chain: str,
    partner_chain: str,
    contact_dist: float = 10.0,
    active_face_resnums: List[int] = None,
) -> Tuple[np.ndarray, int]:
    """Compute centroid of receptor residues near active-face residues.

    Returns (centroid, n_receptor_contact_ca).
    """
    if active_face_resnums is None:
        active_face_resnums = ACTIVE_FACE_RESIDUES

    partner_cas = collect_ca_coords(pose, partner_chain, active_face_resnums)
    if len(partner_cas) == 0:
        raise ValueError("No active-face Cα coordinates found")

    receptor_contact_cas = []
    for i in range(1, pose.total_residue() + 1):
        pdb_info = pose.pdb_info()
        if pdb_info.chain(i) != receptor_chain:
            continue
        res = pose.residue(i)
        if not res.has("CA"):
            continue
        ca_idx = res.atom_index("CA")
        xyz = res.xyz(ca_idx)
        r_ca = np.array([xyz.x, xyz.y, xyz.z])

        dists = np.linalg.norm(partner_cas - r_ca, axis=1)
        if dists.min() < contact_dist:
            receptor_contact_cas.append(r_ca)

    if len(receptor_contact_cas) == 0:
        # Fallback: entire receptor centroid
        all_receptor_cas = []
        for i in range(1, pose.total_residue() + 1):
            if pose.pdb_info().chain(i) != receptor_chain:
                continue
            res = pose.residue(i)
            if not res.has("CA"):
                continue
            ca_idx = res.atom_index("CA")
            xyz = res.xyz(ca_idx)
            all_receptor_cas.append(np.array([xyz.x, xyz.y, xyz.z]))
        receptor_contact_cas = all_receptor_cas

    return np.mean(receptor_contact_cas, axis=0), len(receptor_contact_cas)


# ---------------------------------------------------------------------------
# Core orientation scoring
# ---------------------------------------------------------------------------

def compute_orientation_score(
    pose,
    partner_chain: str = "B",
    receptor_chain: str = "A",
    contact_dist: float = 10.0,
    active_face_residues: Optional[List[int]] = None,
) -> Dict:
    """Compute orientation score for a single docked pose.

    Args:
        active_face_residues: Override the default ACTIVE_FACE_RESIDUES
            (sheets 8+9).  Pass sheets 8+9+12 for sensitivity analysis
            (AC-3.4).

    Returns a dictionary with orientation_score, orientation_class, and
    geometric details.
    """
    if active_face_residues is None:
        active_face_residues = ACTIVE_FACE_RESIDUES
    # 1. Collect active-face Cα coordinates
    active_cas = collect_ca_coords(pose, partner_chain, active_face_residues)
    n_active = len(active_cas)

    if n_active < 3:
        return {
            "orientation_score": 0.0,
            "orientation_class": "insufficient_data",
            "n_active_face_ca": n_active,
            "n_back_face_ca": 0,
            "n_receptor_contact_ca": 0,
            "error": f"Only {n_active} active-face Cα found (need ≥3)",
        }

    # 2. Sheet centroid
    sheet_centroid = active_cas.mean(axis=0)

    # 3. Sheet-plane normal via PCA
    try:
        normal = compute_sheet_plane_normal(active_cas)
    except ValueError as e:
        return {
            "orientation_score": 0.0,
            "orientation_class": "computation_error",
            "n_active_face_ca": n_active,
            "error": str(e),
        }

    # 4. Orient normal toward active face
    try:
        normal = orient_normal_to_active_face(normal, pose, partner_chain)
    except ValueError as e:
        return {
            "orientation_score": 0.0,
            "orientation_class": "probe_error",
            "n_active_face_ca": n_active,
            "error": str(e),
        }

    # 5. Receptor-direction vector
    try:
        receptor_centroid, n_receptor_contact = compute_receptor_interface_centroid(
            pose, receptor_chain, partner_chain, contact_dist
        )
    except ValueError as e:
        return {
            "orientation_score": 0.0,
            "orientation_class": "receptor_error",
            "n_active_face_ca": n_active,
            "error": str(e),
        }

    receptor_dir = receptor_centroid - sheet_centroid
    receptor_dir_norm = np.linalg.norm(receptor_dir)
    if receptor_dir_norm < 1e-6:
        return {
            "orientation_score": 0.0,
            "orientation_class": "degenerate_geometry",
            "n_active_face_ca": n_active,
            "error": "Sheet centroid coincides with receptor centroid",
        }
    receptor_dir = receptor_dir / receptor_dir_norm

    # 6. Dot product
    dot = float(np.dot(normal, receptor_dir))

    # 7. Classify
    if abs(dot) < AMBIGUOUS_BAND:
        orientation_class = "ambiguous"
    elif dot > 0:
        orientation_class = "pass"
    else:
        orientation_class = "fail"

    # 8. Back-face Cα count (supplementary)
    back_cas = collect_ca_coords(pose, partner_chain, BACK_FACE_RESIDUES)

    return {
        "orientation_score": round(dot, 4),
        "orientation_class": orientation_class,
        "n_active_face_ca": n_active,
        "n_back_face_ca": len(back_cas),
        "n_receptor_contact_ca": n_receptor_contact,
        "sheet_centroid": sheet_centroid.tolist(),
        "normal_vector": normal.tolist(),
        "receptor_direction": receptor_dir.tolist(),
    }


# ---------------------------------------------------------------------------
# Phase 1 pipeline integration
# ---------------------------------------------------------------------------

def find_final_pdbs(run_dir: Path) -> List[Path]:
    """Find all final-result PDB files in a pipeline run directory."""
    pdbs = []
    for subdir in run_dir.iterdir():
        if not subdir.is_dir():
            continue
        result_dir = subdir / "final_result"
        if result_dir.exists():
            pdbs.extend(sorted(result_dir.glob("Rank*.pdb")))
        # Also check for PDBs directly in the subdir (legacy layout)
        pdbs.extend(sorted(subdir.glob("Rank*.pdb")))
    # Deduplicate by resolved path
    seen = set()
    unique = []
    for p in pdbs:
        rp = p.resolve()
        if rp not in seen:
            seen.add(rp)
            unique.append(p)
    return unique


def process_state_orientation(
    state_name: str,
    output_base: Path,
    contact_dist: float = 10.0,
    runs_base: Optional[Path] = None,
    active_face_residues: Optional[List[int]] = None,
) -> Optional[Path]:
    """Apply orientation filter to all final models for a receptor state.

    Returns path to orientation_filter_log.csv, or None if no data.
    """
    import pyrosetta
    pyrosetta.init("-mute all", silent=True)

    source_base = runs_base or output_base
    source_state_dir = source_base / state_name
    if not source_state_dir.exists():
        print(f"  State directory not found: {source_state_dir}")
        return None
    state_dir = output_base / state_name
    state_dir.mkdir(parents=True, exist_ok=True)

    all_results = []

    for run_dir in sorted(source_state_dir.iterdir()):
        if not run_dir.is_dir():
            continue
        meta_path = run_dir / "pyrosetta_run_metadata.json"
        if not meta_path.exists():
            continue

        with open(meta_path, encoding="utf-8") as f:
            metadata = json.load(f)

        if metadata.get("dry_run"):
            continue

        receptor_id = metadata.get("receptor_id", state_name)
        seed_index = metadata.get("seed_index", 0)

        pdbs = find_final_pdbs(run_dir)
        if not pdbs:
            print(f"    [SKIP] {run_dir.name}: no final PDB files")
            continue

        n_pass = n_fail = n_ambig = n_error = 0
        for pdb_path in pdbs:
            model_id = pdb_path.name
            try:
                pose = pyrosetta.pose_from_pdb(str(pdb_path))
                result = compute_orientation_score(
                    pose,
                    partner_chain="B",
                    receptor_chain="A",
                    contact_dist=contact_dist,
                    active_face_residues=active_face_residues,
                )
            except Exception as e:
                result = {
                    "orientation_score": 0.0,
                    "orientation_class": "load_error",
                    "n_active_face_ca": 0,
                    "n_back_face_ca": 0,
                    "n_receptor_contact_ca": 0,
                    "error": str(e),
                }

            result["model_id"] = model_id
            result["receptor_id"] = receptor_id
            result["seed_index"] = seed_index
            result["source_pdb"] = str(pdb_path.relative_to(PROJECT_ROOT))

            # Format array fields as strings
            for key in ("sheet_centroid", "normal_vector", "receptor_direction"):
                if key in result and isinstance(result[key], list):
                    result[key] = f"[{','.join(f'{v:.3f}' for v in result[key])}]"

            all_results.append(result)

            oc = result.get("orientation_class", "error")
            if oc == "pass":
                n_pass += 1
            elif oc == "fail":
                n_fail += 1
            elif oc == "ambiguous":
                n_ambig += 1
            else:
                n_error += 1

        print(f"    [OK] {run_dir.name}: {len(pdbs)} models → "
              f"P:{n_pass} F:{n_fail} A:{n_ambig} E:{n_error}")

    if not all_results:
        print(f"  No orientation data for {state_name}")
        return None

    # Write orientation log
    log_path = state_dir / "orientation_filter_log.csv"
    _write_csv(log_path, all_results, ORIENTATION_LOG_COLUMNS)

    # Summary
    classes = [r.get("orientation_class", "error") for r in all_results]
    n_total = len(classes)
    n_pass = classes.count("pass")
    n_fail = classes.count("fail")
    n_ambig = classes.count("ambiguous")
    n_other = n_total - n_pass - n_fail - n_ambig

    print(f"\n  Orientation filter summary for {state_name}:")
    print(f"    Total:     {n_total}")
    print(f"    PASS:      {n_pass:>6d} ({100*n_pass/n_total:.1f}%)")
    print(f"    FAIL:      {n_fail:>6d} ({100*n_fail/n_total:.1f}%)")
    print(f"    AMBIGUOUS: {n_ambig:>6d} ({100*n_ambig/n_total:.1f}%)")
    if n_other:
        print(f"    OTHER:     {n_other:>6d}")
    print(f"    Log: {log_path}")

    return log_path


# ---------------------------------------------------------------------------
# Merge orientation results into interface models table
# ---------------------------------------------------------------------------

def merge_orientation_into_models(
    state_name: str,
    output_base: Path,
) -> bool:
    """Add orientation_score and orientation_class columns to
    pyrosetta_interface_models.csv by joining on model_id.

    Returns True if merge was performed.
    """
    state_dir = output_base / state_name
    log_path = state_dir / "orientation_filter_log.csv"
    models_path = state_dir / "pyrosetta_interface_models.csv"

    if not log_path.exists() or not models_path.exists():
        print(f"  [SKIP] merge for {state_name}: missing files")
        return False

    # Load orientation results keyed by (model_id, seed_index), with legacy
    # fallback indexes for globally unique model ids.
    orient_by_model_seed = {}
    orient_by_model_only = {}
    orient_by_basename_only = {}
    with open(log_path, encoding="utf-8") as f:
        for row in csv.DictReader(f):
            payload = {
                "orientation_score": row.get("orientation_score", ""),
                "orientation_class": row.get("orientation_class", ""),
            }
            model_id = row.get("model_id", "")
            seed_index = _normalize_seed_index(row.get("seed_index", ""))
            basename = Path(model_id).name if model_id else ""

            orient_by_model_seed[(model_id, seed_index)] = payload
            if basename and basename != model_id:
                orient_by_model_seed[(basename, seed_index)] = payload

            _mark_unique_orientation(orient_by_model_only, model_id, payload)
            _mark_unique_orientation(orient_by_basename_only, basename, payload)

    # Read existing models table
    with open(models_path, encoding="utf-8") as f:
        reader = csv.DictReader(f)
        existing_cols = list(reader.fieldnames)
        model_rows = list(reader)

    # Add orientation columns if not already present
    new_cols = list(existing_cols)
    for col in ("orientation_score", "orientation_class"):
        if col not in new_cols:
            new_cols.append(col)

    # Merge
    n_matched = 0
    for row in model_rows:
        mid = row.get("model_id", "")
        seed_index = _normalize_seed_index(row.get("seed_index", ""))
        orient = _lookup_orientation_payload(
            mid,
            seed_index,
            orient_by_model_seed,
            orient_by_model_only,
            orient_by_basename_only,
        )
        row["orientation_score"] = orient.get("orientation_score", "")
        row["orientation_class"] = orient.get("orientation_class", "")
        if orient:
            n_matched += 1

    # Write back
    with open(models_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=new_cols, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(model_rows)

    print(f"  Merged orientation for {state_name}: "
          f"{n_matched}/{len(model_rows)} models matched")
    return True


# ---------------------------------------------------------------------------
# Threshold sensitivity sweep (AC-3.1)
# ---------------------------------------------------------------------------

SWEEP_THRESHOLDS = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30]

SENSITIVITY_COLUMNS = [
    "threshold",
    "n_total",
    "n_pass",
    "n_fail",
    "n_ambiguous",
    "pass_rate",
    "fail_rate",
    "ambiguous_rate",
    "n_hotspot_residues",
    "hotspot_residues",          # semicolon-separated residue IDs
    "n_sheet89_in_hotspot",      # sheet 8/9 residues in hotspot
    "sheet89_stable",            # True if all core sheet 8/9 residues remain
    "hotspot_robustness",        # Jaccard with next-lower threshold
]


def _reclassify_orientation(
    score: float, threshold: float,
) -> str:
    """Re-classify an orientation score at a given threshold."""
    if abs(score) < threshold:
        return "ambiguous"
    return "pass" if score > 0 else "fail"


def _compute_occupancy_hotspots(
    pass_model_ids: set,
    residues: List[dict],
    hotspot_threshold: float = 0.50,
) -> set:
    """Compute hotspot residue IDs from PASS models' interface residues.

    Returns set of (chain, residue_id) tuples that exceed occupancy threshold.
    """
    if not pass_model_ids:
        return set()

    residue_counts: Dict[Tuple[str, str], int] = {}
    for r in residues:
        if r.get("model_id", "") not in pass_model_ids:
            continue
        key = (r.get("chain", ""), r.get("residue_id", ""))
        residue_counts[key] = residue_counts.get(key, 0) + 1

    n_pass = len(pass_model_ids)
    return {
        key for key, count in residue_counts.items()
        if count / n_pass >= hotspot_threshold
    }


def sweep_threshold_sensitivity(
    orientation_log: List[dict],
    residues: List[dict],
    thresholds: Optional[List[float]] = None,
    hotspot_threshold: float = 0.50,
) -> Tuple[List[dict], List[str]]:
    """Run threshold sensitivity sweep on existing orientation scores.

    Re-classifies models at each threshold and computes hotspot sets.
    Does NOT require PyRosetta — works entirely on pre-computed scores.

    Args:
        orientation_log: Rows from orientation_filter_log.csv
        residues: Rows from pyrosetta_interface_residue_table.csv
        thresholds: List of thresholds to sweep (default: SWEEP_THRESHOLDS)
        hotspot_threshold: Occupancy fraction for hotspot designation

    Returns:
        (sensitivity_rows, warnings) for orientation_threshold_sensitivity.csv
    """
    if thresholds is None:
        thresholds = list(SWEEP_THRESHOLDS)

    # Filter to models with valid numeric orientation scores
    scored_models = []
    for row in orientation_log:
        try:
            score = float(row.get("orientation_score", ""))
        except (TypeError, ValueError):
            continue
        scored_models.append((row.get("model_id", ""), score))

    n_total = len(scored_models)
    if n_total == 0:
        return [], ["WARNING: No scored models found in orientation log"]

    from egfr_pipeline.region_definitions import KO_SHEET_8_9

    prev_hotspot_ids: Optional[set] = None
    sensitivity_rows: List[dict] = []
    warnings: List[str] = []

    for threshold in sorted(thresholds):
        # Re-classify
        pass_ids = set()
        n_pass = n_fail = n_ambig = 0
        for model_id, score in scored_models:
            cls = _reclassify_orientation(score, threshold)
            if cls == "pass":
                n_pass += 1
                pass_ids.add(model_id)
            elif cls == "fail":
                n_fail += 1
            else:
                n_ambig += 1

        # Compute hotspot residues for this threshold
        hotspot_keys = _compute_occupancy_hotspots(pass_ids, residues, hotspot_threshold)
        hotspot_ids = sorted({rid for _, rid in hotspot_keys})

        # Check sheet 8/9 stability
        hotspot_resnums = set()
        for rid in hotspot_ids:
            try:
                hotspot_resnums.add(int("".join(c for c in rid if c.isdigit())))
            except ValueError:
                pass
        n_sheet89 = len(hotspot_resnums & KO_SHEET_8_9)
        # Core sheet 8/9: at least the probe residues (962, 964, 971)
        core_sheet89 = {962, 964, 971}
        sheet89_stable = core_sheet89.issubset(hotspot_resnums)

        # Hotspot robustness vs previous threshold (Jaccard)
        current_set = set(hotspot_ids)
        if prev_hotspot_ids is not None and (current_set or prev_hotspot_ids):
            intersection = len(current_set & prev_hotspot_ids)
            union = len(current_set | prev_hotspot_ids)
            robustness = round(intersection / union, 4) if union > 0 else 1.0
        else:
            robustness = ""

        sensitivity_rows.append({
            "threshold": threshold,
            "n_total": n_total,
            "n_pass": n_pass,
            "n_fail": n_fail,
            "n_ambiguous": n_ambig,
            "pass_rate": round(n_pass / n_total, 4) if n_total else 0.0,
            "fail_rate": round(n_fail / n_total, 4) if n_total else 0.0,
            "ambiguous_rate": round(n_ambig / n_total, 4) if n_total else 0.0,
            "n_hotspot_residues": len(hotspot_ids),
            "hotspot_residues": ";".join(hotspot_ids),
            "n_sheet89_in_hotspot": n_sheet89,
            "sheet89_stable": sheet89_stable,
            "hotspot_robustness": robustness,
        })

        prev_hotspot_ids = current_set

    # EC-3.1: all thresholds same pass rate → warning
    pass_rates = [r["pass_rate"] for r in sensitivity_rows]
    if len(set(pass_rates)) == 1:
        warnings.append(
            "Threshold에 무관한 결과, 현재 값 유지 (all thresholds yield identical pass rate)"
        )

    # C-1.3: Generate recommendation
    if sensitivity_rows:
        # Find threshold with best balance: sheet89 stable + reasonable ambiguous rate
        best = None
        for row in sensitivity_rows:
            if row["sheet89_stable"] and row["ambiguous_rate"] < 0.20:
                if best is None or row["ambiguous_rate"] < best["ambiguous_rate"]:
                    best = row
        if best:
            warnings.append(
                f"Recommended threshold: {best['threshold']} "
                f"(pass_rate={best['pass_rate']}, ambiguous={best['ambiguous_rate']}, "
                f"sheet89_stable=True)"
            )

    return sensitivity_rows, warnings


def sweep_from_state_dir(
    state_dir: Path,
    output_path: Optional[Path] = None,
    thresholds: Optional[List[float]] = None,
) -> Optional[Path]:
    """Run threshold sweep for a single state directory.

    Reads orientation_filter_log.csv and pyrosetta_interface_residue_table.csv
    from state_dir.

    Returns path to output CSV, or None if insufficient data.
    """
    from egfr_pipeline.phase1.cluster_consensus import load_residue_table

    log_path = state_dir / "orientation_filter_log.csv"
    if not log_path.exists():
        print(f"  [SKIP] No orientation log in {state_dir}")
        return None

    with open(log_path, encoding="utf-8") as f:
        orientation_log = list(csv.DictReader(f))

    residues = load_residue_table(state_dir)
    if not residues:
        print(f"  [SKIP] No residue table in {state_dir}")
        return None

    rows, warnings = sweep_threshold_sensitivity(orientation_log, residues, thresholds)
    for w in warnings:
        print(f"  {w}")

    if not rows:
        return None

    out = output_path or (state_dir / "orientation_threshold_sensitivity.csv")
    _write_csv(out, rows, SENSITIVITY_COLUMNS)
    print(f"  Threshold sensitivity: {out}")
    return out


# ---------------------------------------------------------------------------
# Ambiguous model report (AC-3.2)
# ---------------------------------------------------------------------------

AMBIGUOUS_REPORT_COLUMNS = [
    "receptor_id",
    "seed_index",
    "n_total",
    "n_ambiguous",
    "ambiguous_rate",
    "mean_dG_pass",
    "mean_dG_ambiguous",
    "dG_difference",
    "unique_ambiguous_residues",
    "n_unique_ambiguous_residues",
]


def compute_ambiguous_report(
    orientation_log: List[dict],
    models: List[dict],
    residues: List[dict],
) -> Tuple[List[dict], List[str]]:
    """Generate ambiguous model report per (receptor_id, seed_index).

    Args:
        orientation_log: Rows from orientation_filter_log.csv
        models: Rows from pyrosetta_interface_models.csv (for dG_separated)
        residues: Rows from pyrosetta_interface_residue_table.csv

    Returns:
        (report_rows, warnings)
    """
    # Build model_id → dG_separated index from models table
    dg_index: Dict[str, float] = {}
    for m in models:
        mid = m.get("model_id", "")
        try:
            dg_index[mid] = float(m.get("dG_separated", ""))
        except (TypeError, ValueError):
            pass

    # Build model_id → set of interface residue IDs
    residues_by_model: Dict[str, set] = {}
    for r in residues:
        mid = r.get("model_id", "")
        residues_by_model.setdefault(mid, set()).add(r.get("residue_id", ""))

    # Group orientation log by (receptor_id, seed_index)
    from collections import defaultdict
    groups: Dict[Tuple[str, str], List[dict]] = defaultdict(list)
    for row in orientation_log:
        key = (row.get("receptor_id", ""), row.get("seed_index", ""))
        groups[key].append(row)

    report_rows: List[dict] = []
    warnings: List[str] = []

    for (receptor_id, seed_index), rows in sorted(groups.items()):
        n_total = len(rows)
        pass_models = [r for r in rows if r.get("orientation_class") == "pass"]
        ambig_models = [r for r in rows if r.get("orientation_class") == "ambiguous"]
        n_ambig = len(ambig_models)

        if n_total == 0:
            continue

        ambig_rate = round(n_ambig / n_total, 4)

        # Mean dG for pass vs ambiguous
        pass_dgs = [dg_index[m["model_id"]] for m in pass_models if m["model_id"] in dg_index]
        ambig_dgs = [dg_index[m["model_id"]] for m in ambig_models if m["model_id"] in dg_index]

        mean_dg_pass = round(sum(pass_dgs) / len(pass_dgs), 4) if pass_dgs else ""
        mean_dg_ambig = round(sum(ambig_dgs) / len(ambig_dgs), 4) if ambig_dgs else ""

        dg_diff = ""
        if isinstance(mean_dg_pass, float) and isinstance(mean_dg_ambig, float):
            dg_diff = round(mean_dg_ambig - mean_dg_pass, 4)

        # Unique residues: in ambiguous models but NOT in any pass model
        pass_residue_union: set = set()
        for m in pass_models:
            pass_residue_union |= residues_by_model.get(m["model_id"], set())

        ambig_residue_union: set = set()
        for m in ambig_models:
            ambig_residue_union |= residues_by_model.get(m["model_id"], set())

        unique_ambig = sorted(ambig_residue_union - pass_residue_union)

        report_rows.append({
            "receptor_id": receptor_id,
            "seed_index": seed_index,
            "n_total": n_total,
            "n_ambiguous": n_ambig,
            "ambiguous_rate": ambig_rate,
            "mean_dG_pass": mean_dg_pass,
            "mean_dG_ambiguous": mean_dg_ambig,
            "dG_difference": dg_diff,
            "unique_ambiguous_residues": ";".join(unique_ambig),
            "n_unique_ambiguous_residues": len(unique_ambig),
        })

    # Summary warnings
    all_unique = set()
    for r in report_rows:
        res_str = r.get("unique_ambiguous_residues", "")
        if res_str:
            all_unique.update(res_str.split(";"))

    if len(all_unique) >= 3:
        warnings.append(
            f"Ambiguous-only unique residues ≥ 3 ({len(all_unique)}): "
            f"추가 조사 권장. Residues: {';'.join(sorted(all_unique)[:10])}"
        )

    return report_rows, warnings


def ambiguous_report_from_state_dir(
    state_dir: Path,
    output_path: Optional[Path] = None,
) -> Optional[Path]:
    """Generate ambiguous report for a single state.

    Returns path to output CSV, or None if insufficient data.
    """
    from egfr_pipeline.phase1.cluster_consensus import (
        load_models_table,
        load_residue_table,
    )

    log_path = state_dir / "orientation_filter_log.csv"
    if not log_path.exists():
        return None

    with open(log_path, encoding="utf-8") as f:
        orientation_log = list(csv.DictReader(f))

    models = load_models_table(state_dir)
    residues = load_residue_table(state_dir)

    rows, warnings = compute_ambiguous_report(orientation_log, models, residues)
    for w in warnings:
        print(f"  {w}")

    if not rows:
        return None

    out = output_path or (state_dir / "orientation_ambiguous_report.csv")
    _write_csv(out, rows, AMBIGUOUS_REPORT_COLUMNS)
    print(f"  Ambiguous report: {out}")
    return out


# ---------------------------------------------------------------------------
# Sheet 12 sensitivity analysis (AC-3.4)
# ---------------------------------------------------------------------------

SHEET12_SENSITIVITY_COLUMNS = [
    "setting",
    "active_face_sheets",
    "n_total",
    "n_pass",
    "n_fail",
    "n_ambiguous",
    "pass_rate",
    "n_hotspot_residues",
    "hotspot_residues",
    "hotspot_overlap",
]

# Sheet 8+9+12 combined for sensitivity variant
ACTIVE_FACE_WITH_SHEET12 = ACTIVE_FACE_RESIDUES + SHEET_12_RESIDUES


def compare_sheet12_sensitivity(
    log_a: List[dict],
    log_b: List[dict],
    residues: List[dict],
    hotspot_threshold: float = 0.50,
) -> List[dict]:
    """Compare orientation results from two active face settings.

    Args:
        log_a: Orientation log from setting A (sheet 8+9 only)
        log_b: Orientation log from setting B (sheet 8+9+12)
        residues: Interface residue table (shared — same models)
        hotspot_threshold: Occupancy threshold for hotspot designation

    Returns:
        2-row list for sheet12_sensitivity.csv with hotspot_overlap.
    """
    def _summarize(log: List[dict], setting: str, sheets: str) -> dict:
        n_total = len(log)
        pass_ids = set()
        n_pass = n_fail = n_ambig = 0
        for row in log:
            cls = row.get("orientation_class", "")
            if cls == "pass":
                n_pass += 1
                pass_ids.add(row.get("model_id", ""))
            elif cls == "fail":
                n_fail += 1
            elif cls == "ambiguous":
                n_ambig += 1

        hotspot_keys = _compute_occupancy_hotspots(pass_ids, residues, hotspot_threshold)
        hotspot_ids = sorted({rid for _, rid in hotspot_keys})

        return {
            "setting": setting,
            "active_face_sheets": sheets,
            "n_total": n_total,
            "n_pass": n_pass,
            "n_fail": n_fail,
            "n_ambiguous": n_ambig,
            "pass_rate": round(n_pass / n_total, 4) if n_total else 0.0,
            "n_hotspot_residues": len(hotspot_ids),
            "hotspot_residues": ";".join(hotspot_ids),
            "_hotspot_set": set(hotspot_ids),
        }

    row_a = _summarize(log_a, "A", "sheet_8+9")
    row_b = _summarize(log_b, "B", "sheet_8+9+12")

    # Compute hotspot overlap (Jaccard)
    set_a = row_a.pop("_hotspot_set")
    set_b = row_b.pop("_hotspot_set")

    if set_a or set_b:
        intersection = len(set_a & set_b)
        union = len(set_a | set_b)
        overlap = round(intersection / union, 4) if union > 0 else 0.0
    else:
        overlap = ""

    row_a["hotspot_overlap"] = overlap
    row_b["hotspot_overlap"] = overlap

    return [row_a, row_b]


def sheet12_sensitivity_from_logs(
    log_a_path: Path,
    log_b_path: Path,
    residue_table_path: Path,
    output_path: Path,
) -> Path:
    """Compare two orientation logs and write sheet12_sensitivity.csv.

    Args:
        log_a_path: Orientation log from setting A (sheet 8+9)
        log_b_path: Orientation log from setting B (sheet 8+9+12)
        residue_table_path: pyrosetta_interface_residue_table.csv
        output_path: Output CSV path

    Returns:
        Path to output CSV.
    """
    with open(log_a_path, encoding="utf-8") as f:
        log_a = list(csv.DictReader(f))
    with open(log_b_path, encoding="utf-8") as f:
        log_b = list(csv.DictReader(f))
    with open(residue_table_path, encoding="utf-8") as f:
        residues = list(csv.DictReader(f))

    rows = compare_sheet12_sensitivity(log_a, log_b, residues)
    _write_csv(output_path, rows, SHEET12_SENSITIVITY_COLUMNS)
    print(f"  Sheet 12 sensitivity: {output_path}")

    # Interpretation
    if rows and rows[0]["hotspot_overlap"] != "":
        overlap = float(rows[0]["hotspot_overlap"])
        if overlap >= 0.80:
            print(f"  → Overlap {overlap:.0%} ≥ 80%: sheet 12 inclusion has minimal effect. Maintain current (8+9).")
        elif overlap >= 0.50:
            print(f"  → Overlap {overlap:.0%} 50-80%: sheet 12 has moderate effect. Review recommended.")
        else:
            print(f"  → Overlap {overlap:.0%} < 50%: sheet 12 has major effect. Report both settings.")

    return output_path


# ---------------------------------------------------------------------------
# Pilot retroactive validation
# ---------------------------------------------------------------------------

def validate_pilot_structures(
    pilot_dir: Path,
    out_csv: Path,
    contact_dist: float = 10.0,
) -> None:
    """Apply orientation filter to pilot structures for calibration.

    Expected pilot structures: C02_M01, C02_M03, C04_M01, C04_M02, C07_M03
    """
    import pyrosetta
    pyrosetta.init("-mute all", silent=True)

    pdb_files = sorted(pilot_dir.glob("*.pdb"))
    if not pdb_files:
        print(f"No PDB files found in {pilot_dir}")
        return

    print(f"Pilot validation: {len(pdb_files)} structures in {pilot_dir}")
    results = []

    for pdb_path in pdb_files:
        fname = pdb_path.name
        try:
            pose = pyrosetta.pose_from_pdb(str(pdb_path))
            result = compute_orientation_score(
                pose,
                partner_chain="B",
                receptor_chain="A",
                contact_dist=contact_dist,
            )
        except Exception as e:
            result = {
                "orientation_score": 0.0,
                "orientation_class": "load_error",
                "error": str(e),
            }

        result["filename"] = fname
        result["source"] = "pilot_retroactive"

        for key in ("sheet_centroid", "normal_vector", "receptor_direction"):
            if key in result and isinstance(result[key], list):
                result[key] = f"[{','.join(f'{v:.3f}' for v in result[key])}]"

        results.append(result)
        oc = result.get("orientation_class", "error")
        print(f"  {fname}: orientation_score={result.get('orientation_score', 0):.4f} → {oc}")

    # Write pilot validation CSV
    pilot_cols = [
        "filename", "source", "orientation_score", "orientation_class",
        "n_active_face_ca", "n_back_face_ca", "n_receptor_contact_ca",
        "sheet_centroid", "normal_vector", "receptor_direction", "error",
    ]
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=pilot_cols, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(results)

    # Summary
    n_pass = sum(1 for r in results if r.get("orientation_class") == "pass")
    n_fail = sum(1 for r in results if r.get("orientation_class") == "fail")
    n_ambig = sum(1 for r in results if r.get("orientation_class") == "ambiguous")
    print(f"\nPilot validation summary: P:{n_pass} F:{n_fail} A:{n_ambig} / {len(results)} total")
    print(f"Output: {out_csv}")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _write_csv(path: Path, rows: List[dict], columns: List[str]) -> None:
    """Write rows to CSV with specified columns."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _normalize_seed_index(seed_index: object) -> str:
    """Normalize seed index for stable cross-table joins."""
    return str(seed_index).strip() if seed_index is not None else ""


def _mark_unique_orientation(index: Dict[str, Optional[dict]], key: str, payload: dict) -> None:
    """Track a payload only when the key is globally unique."""
    if not key:
        return
    if key in index:
        index[key] = None
        return
    index[key] = payload


def _lookup_orientation_payload(
    model_id: str,
    seed_index: str,
    orient_by_model_seed: Dict[tuple, dict],
    orient_by_model_only: Dict[str, Optional[dict]],
    orient_by_basename_only: Dict[str, Optional[dict]],
) -> dict:
    """Resolve orientation payload with seed-aware and legacy-compatible fallbacks."""
    basename = Path(model_id).name if model_id else ""

    for candidate in ((model_id, seed_index), (basename, seed_index)):
        payload = orient_by_model_seed.get(candidate)
        if payload:
            return payload

    for index, key in ((orient_by_model_only, model_id), (orient_by_basename_only, basename)):
        payload = index.get(key)
        if payload:
            return payload

    return {}


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Phase 1 TG 1.2A: Orientation-aware filtering for beta-meander docking"
    )
    parser.add_argument("--state", choices=RECEPTOR_STATES,
                        help="Process a single receptor state")
    parser.add_argument("--runs_dir", type=Path,
                        default=PHASE1_RUNS_DIR,
                        help="Workflow A Phase 2 PyRosetta run directory")
    parser.add_argument("--output_dir", type=Path,
                        default=PHASE1_OUTPUT_DIR,
                        help="Workflow B Phase 1 analysis output directory")
    parser.add_argument("--contact_dist", type=float, default=10.0,
                        help="Distance cutoff for receptor contact centroid (default: 10.0 Å)")
    parser.add_argument("--merge", action="store_true",
                        help="Merge orientation columns into interface models table")
    parser.add_argument("--sweep", action="store_true",
                        help="Run threshold sensitivity sweep (AC-3.1)")
    parser.add_argument("--sheet12_log", type=Path,
                        help="Compare current log with sheet 8+9+12 log (AC-3.4)")
    parser.add_argument("--active_face", choices=["8+9", "8+9+12"],
                        help="Override active face residues for orientation scoring")
    parser.add_argument("--pilot_dir", type=Path,
                        help="Directory with pilot PDB files for retroactive validation")
    parser.add_argument("--pilot_out", type=Path,
                        help="Output CSV for pilot validation results")
    args = parser.parse_args()

    # Threshold sensitivity sweep mode (AC-3.1)
    if args.sweep:
        print("Phase 1 — Threshold Sensitivity Sweep (AC-3.1)")
        states = [args.state] if args.state else RECEPTOR_STATES
        for state in states:
            state_dir = args.output_dir / state
            print(f"\nSweeping {state}:")
            sweep_from_state_dir(state_dir)
        return 0

    # Sheet 12 sensitivity comparison mode (AC-3.4)
    if args.sheet12_log:
        print("Phase 1 — Sheet 12 Sensitivity Comparison (AC-3.4)")
        states = [args.state] if args.state else RECEPTOR_STATES
        for state in states:
            state_dir = args.output_dir / state
            log_a = state_dir / "orientation_filter_log.csv"
            residue_table = state_dir / "pyrosetta_interface_residue_table.csv"
            if not log_a.exists():
                print(f"  [SKIP] {state}: no orientation log")
                continue
            if not residue_table.exists():
                print(f"  [SKIP] {state}: no residue table")
                continue
            out = state_dir / "sheet12_sensitivity.csv"
            sheet12_sensitivity_from_logs(log_a, args.sheet12_log, residue_table, out)
        return 0

    # Pilot validation mode
    if args.pilot_dir:
        pilot_out = args.pilot_out or (
            args.output_dir / "orientation_filter_pilot_validation.csv"
        )
        validate_pilot_structures(args.pilot_dir, pilot_out, args.contact_dist)
        return 0

    print("Phase 1 — Task 1.2A: Orientation-Aware Filtering")
    print(f"Runs base: {args.runs_dir}")
    print(f"Output base: {args.output_dir}")
    print(f"Contact distance: {args.contact_dist} Å")
    print(f"Ambiguous band: |dot product| < {AMBIGUOUS_BAND}")
    print(f"Active face residues: sheets 8 ({SHEET_8_RESIDUES}) + 9 ({SHEET_9_RESIDUES})")
    print(f"Probe residues: {PROBE_RESIDUES}")

    states = [args.state] if args.state else RECEPTOR_STATES
    results = []

    for state in states:
        print(f"\nProcessing {state}:")
        active_override = None
        if getattr(args, "active_face", None) == "8+9+12":
            active_override = ACTIVE_FACE_WITH_SHEET12
        log_path = process_state_orientation(
            state,
            args.output_dir,
            args.contact_dist,
            runs_base=args.runs_dir,
            active_face_residues=active_override,
        )
        if log_path:
            results.append((state, log_path))

    if not results:
        print("\nNo docking results found. Run docking first.")
        return 0

    # Merge orientation into models table
    if args.merge:
        print(f"\n{'='*60}")
        print("Merging orientation results into interface models tables")
        for state, _ in results:
            merge_orientation_into_models(state, args.output_dir)

    print(f"\n{'='*60}")
    print("ORIENTATION FILTERING COMPLETE")
    print(f"{'='*60}")
    for state, log_path in results:
        print(f"  {state}: {log_path.name}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
