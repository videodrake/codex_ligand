#!/usr/bin/env python3
"""
orientation_filter.py
=====================
Beta-meander face-flip detection for EGFR–MYO1D PPI docking poses.

Purpose:
    Determines whether the active face (sheets 8 & 9) of the MYO1D
    beta-meander is oriented toward the EGFR receptor surface (correct)
    or flipped away from it (artifact).

Scientific rationale:
    A beta-meander is a flat, thin β-sheet ribbon. In global docking,
    the structure can land on the receptor in two orientations:
    - Active-face-down (correct): functional side chains of sheets 8/9
      pack against receptor surface.
    - Active-face-up  (flipped): backbone / back-face of the sheet
      faces the receptor; side chains point into solvent.
    Simple residue contact counts cannot distinguish these because
    backbone atoms or edge-on contacts can occur in both orientations.

Algorithm — Dual-Vector Orientation Test:
    1. Define the beta-sheet plane from Cα coordinates of sheets 8+9.
    2. Compute the sheet-plane normal via cross product of two
       in-plane vectors (strand direction × inter-strand direction).
    3. Choose the normal direction that points toward the active face
       using a Cα→Cβ probe on a canonical active-face residue.
    4. Compute a receptor-direction vector from the sheet centroid
       toward the receptor interface centroid.
    5. The dot product of (active-face normal) · (receptor direction)
       determines orientation:
         > 0 → active face toward receptor (PASS)
         < 0 → active face away from receptor (FAIL)
         ≈ 0 → edge-on / ambiguous (AMBIGUOUS)

    This is more robust than a single Cα→Cβ vector because:
    - A single residue's Cβ may be rotamer-dependent.
    - The sheet-plane normal is a global geometric property of the
      entire active-face β-sheet, making it resilient to local noise.

Dependencies:
    - PyRosetta 4 (for pose loading and residue access)
    - NumPy (for linear algebra)

Usage:
    python orientation_filter.py \\
        --pdb_dir ./docking_output \\
        --out_csv ./orientation_filter_log.csv \\
        --partner_chain B \\
        --receptor_chain A

Author: EGFR-MYO1D Pipeline
Date: 2026-03
"""

import argparse
import csv
import glob
import os
import sys
from typing import Dict, List, Optional, Tuple

import numpy as np

# ---------------------------------------------------------------------------
# Constants: beta-meander sheet definitions (MYO1D residue numbering)
# ---------------------------------------------------------------------------

# Sheet boundaries — residue numbers in MYO1D (no offset)
SHEET_8_RESIDUES = [961, 962, 963, 964]
SHEET_9_RESIDUES = [968, 969, 970, 971, 972]
SHEET_10_RESIDUES = [977, 978, 979, 980]       # MUST VERIFY from actual structure / Ko et al. Table S2
SHEET_11_RESIDUES = [985, 986, 987, 988]       # MUST VERIFY from actual structure / Ko et al. Table S2
SHEET_12_RESIDUES = [993, 994, 995, 996, 997]  # MUST VERIFY from actual structure / Ko et al. Table S2
# WARNING: These back-face residue definitions are approximate pending
# verification against the actual MYO1D beta-meander structure.
# They affect back-face contact counting but NOT the primary orientation
# score (which uses only sheets 8+9). Verify before using back-face
# contact ratios for any gating decision.

# Active face is defined by sheets 8 + 9 (experimentally validated)
ACTIVE_FACE_RESIDUES = SHEET_8_RESIDUES + SHEET_9_RESIDUES

# Back face is primarily sheets 10, 11, and structural scaffold of 12
BACK_FACE_RESIDUES = SHEET_10_RESIDUES + SHEET_11_RESIDUES + SHEET_12_RESIDUES

# Canonical probe residues for Cα→Cβ direction disambiguation
# Multiple probes provide consensus-based disambiguation, reducing
# sensitivity to any single residue's local rotamer state.
# All should be on the active face (sheets 8/9).
PROBE_RESIDUES = [962, 964, 971]  # VAL962 (sheet 8), VAL964 (sheet 8), SER971 (sheet 9)
PROBE_RESIDUE = 962  # Legacy single-probe fallback

# Orientation thresholds
PASS_THRESHOLD = 0.0      # dot product > 0 → active face toward receptor
AMBIGUOUS_BAND = 0.15     # |dot product| < this → ambiguous (edge-on)


# ---------------------------------------------------------------------------
# Geometry utilities
# ---------------------------------------------------------------------------

def get_ca_coord(pose, chain: str, resnum: int) -> Optional[np.ndarray]:
    """Extract Cα coordinate for a given chain and residue number."""
    for i in range(1, pose.total_residue() + 1):
        pdb_info = pose.pdb_info()
        if (pdb_info.chain(i) == chain and
                pdb_info.number(i) == resnum):
            ca_idx = pose.residue(i).atom_index("CA")
            xyz = pose.residue(i).xyz(ca_idx)
            return np.array([xyz.x, xyz.y, xyz.z])
    return None


def get_cb_coord(pose, chain: str, resnum: int) -> Optional[np.ndarray]:
    """Extract Cβ coordinate (or Cα for GLY) for a given chain and residue."""
    for i in range(1, pose.total_residue() + 1):
        pdb_info = pose.pdb_info()
        if (pdb_info.chain(i) == chain and
                pdb_info.number(i) == resnum):
            res = pose.residue(i)
            if res.has("CB"):
                cb_idx = res.atom_index("CB")
                xyz = res.xyz(cb_idx)
                return np.array([xyz.x, xyz.y, xyz.z])
            else:
                # GLY: fall back to CA
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
    return np.array(coords)


def compute_sheet_plane_normal(coords: np.ndarray) -> np.ndarray:
    """
    Compute the sheet-plane normal from a set of Cα coordinates using
    PCA (principal component analysis).

    The sheet plane is defined by the two largest principal components.
    The normal is the third (smallest variance) component.

    This is more robust than a cross product of two manually chosen
    vectors because it uses all available Cα positions and is not
    sensitive to the choice of specific residue pairs.
    """
    if len(coords) < 3:
        raise ValueError("Need at least 3 Cα coordinates to define a plane")

    centroid = coords.mean(axis=0)
    centered = coords - centroid

    # SVD: columns of V^T are principal components
    _, s, vt = np.linalg.svd(centered, full_matrices=False)
    # The last row of V^T is the direction of least variance = plane normal
    normal = vt[-1]
    return normal / np.linalg.norm(normal)


def orient_normal_to_active_face(
    normal: np.ndarray,
    pose,
    partner_chain: str,
    probe_resnums: List[int] = None
) -> np.ndarray:
    """
    Ensure the normal vector points toward the active face (side chains
    of sheets 8/9) rather than the back face.

    Method: compute the Cα→Cβ vector for multiple active-face residues
    and take a consensus vote. If the majority of probes indicate the
    normal should be flipped, flip it. This is more robust than relying
    on a single probe residue, which may have an atypical rotamer.

    Falls back to single-probe if only one probe is available.
    """
    if probe_resnums is None:
        probe_resnums = PROBE_RESIDUES

    votes = 0  # positive = keep normal, negative = flip
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
        # Last resort: try legacy single probe
        ca = get_ca_coord(pose, partner_chain, PROBE_RESIDUE)
        cb = get_cb_coord(pose, partner_chain, PROBE_RESIDUE)
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
    active_face_resnums: List[int] = None
) -> np.ndarray:
    """
    Compute the centroid of receptor residues that are near the
    active-face residues of the partner.

    This gives a more targeted receptor-direction vector than using
    the entire receptor centroid, because the full receptor is large
    and its geometric center may be far from the actual contact area.
    """
    if active_face_resnums is None:
        active_face_resnums = ACTIVE_FACE_RESIDUES

    # Collect partner active-face Cα positions
    partner_cas = collect_ca_coords(pose, partner_chain, active_face_resnums)
    if len(partner_cas) == 0:
        raise ValueError("No active-face Cα coordinates found")

    # Find receptor residues within contact_dist of any active-face Cα
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

        # Check distance to any partner active-face Cα
        dists = np.linalg.norm(partner_cas - r_ca, axis=1)
        if dists.min() < contact_dist:
            receptor_contact_cas.append(r_ca)

    if len(receptor_contact_cas) == 0:
        # Fallback: use entire receptor centroid
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

    return np.mean(receptor_contact_cas, axis=0)


# ---------------------------------------------------------------------------
# Core orientation scoring
# ---------------------------------------------------------------------------

def compute_orientation_score(
    pose,
    partner_chain: str = "B",
    receptor_chain: str = "A",
    probe_resnum: int = PROBE_RESIDUE,
    contact_dist: float = 10.0,
) -> Dict:
    """
    Compute the full orientation score for a single docked pose.

    Returns a dictionary with:
      - orientation_score: float (dot product, -1 to +1)
      - orientation_class: str (pass / fail / ambiguous)
      - sheet_centroid: array
      - normal_vector: array
      - receptor_direction: array
      - n_active_face_ca: int (how many active-face Cα were found)
      - n_receptor_contact_ca: int
    """
    # 1. Collect active-face Cα coordinates
    active_cas = collect_ca_coords(pose, partner_chain, ACTIVE_FACE_RESIDUES)
    n_active = len(active_cas)

    if n_active < 3:
        return {
            "orientation_score": 0.0,
            "orientation_class": "insufficient_data",
            "n_active_face_ca": n_active,
            "n_receptor_contact_ca": 0,
            "error": f"Only {n_active} active-face Cα found (need ≥ 3)"
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
            "error": str(e)
        }

    # 4. Orient normal toward active face using multi-probe consensus
    try:
        normal = orient_normal_to_active_face(
            normal, pose, partner_chain, PROBE_RESIDUES
        )
    except ValueError as e:
        return {
            "orientation_score": 0.0,
            "orientation_class": "probe_error",
            "n_active_face_ca": n_active,
            "error": str(e)
        }

    # 5. Receptor-direction vector
    try:
        receptor_centroid = compute_receptor_interface_centroid(
            pose, receptor_chain, partner_chain, contact_dist
        )
    except ValueError as e:
        return {
            "orientation_score": 0.0,
            "orientation_class": "receptor_error",
            "n_active_face_ca": n_active,
            "error": str(e)
        }

    receptor_dir = receptor_centroid - sheet_centroid
    receptor_dir_norm = np.linalg.norm(receptor_dir)
    if receptor_dir_norm < 1e-6:
        return {
            "orientation_score": 0.0,
            "orientation_class": "degenerate_geometry",
            "n_active_face_ca": n_active,
            "error": "Sheet centroid coincides with receptor centroid"
        }
    receptor_dir = receptor_dir / receptor_dir_norm

    # 6. Dot product = orientation score
    dot = float(np.dot(normal, receptor_dir))

    # 7. Classify
    if abs(dot) < AMBIGUOUS_BAND:
        orientation_class = "ambiguous"
    elif dot > 0:
        orientation_class = "pass"
    else:
        orientation_class = "fail"

    # 8. Supplementary: count active vs back face contacts
    back_cas = collect_ca_coords(pose, partner_chain, BACK_FACE_RESIDUES)
    # (contact counting would require per-residue distance checks,
    #  but the primary filter is the geometric orientation score above)

    return {
        "orientation_score": round(dot, 4),
        "orientation_class": orientation_class,
        "n_active_face_ca": n_active,
        "n_back_face_ca": len(back_cas),
        "sheet_centroid": sheet_centroid.tolist(),
        "normal_vector": normal.tolist(),
        "receptor_direction": receptor_dir.tolist(),
    }


# ---------------------------------------------------------------------------
# Batch processing
# ---------------------------------------------------------------------------

def process_directory(
    pdb_dir: str,
    out_csv: str,
    partner_chain: str = "B",
    receptor_chain: str = "A",
    probe_resnum: int = PROBE_RESIDUE,
    contact_dist: float = 10.0,
) -> None:
    """Process all PDB files in a directory and write orientation log."""

    import pyrosetta
    pyrosetta.init("-mute all", silent=True)

    pdb_files = sorted(glob.glob(os.path.join(pdb_dir, "*.pdb")))
    if not pdb_files:
        print(f"No PDB files found in {pdb_dir}")
        sys.exit(1)

    print(f"Processing {len(pdb_files)} PDB files...")

    results = []
    for pdb_path in pdb_files:
        fname = os.path.basename(pdb_path)
        try:
            pose = pyrosetta.pose_from_pdb(pdb_path)
            score = compute_orientation_score(
                pose, partner_chain, receptor_chain,
                probe_resnum, contact_dist
            )
            score["filename"] = fname
            results.append(score)
        except Exception as e:
            results.append({
                "filename": fname,
                "orientation_score": 0.0,
                "orientation_class": "load_error",
                "error": str(e)
            })

        # Progress indicator
        if len(results) % 100 == 0:
            print(f"  ... processed {len(results)} / {len(pdb_files)}")

    # Write CSV
    fieldnames = [
        "filename", "orientation_score", "orientation_class",
        "n_active_face_ca", "n_back_face_ca",
        "sheet_centroid", "normal_vector", "receptor_direction",
        "error"
    ]
    with open(out_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for r in results:
            writer.writerow(r)

    # Summary
    classes = [r.get("orientation_class", "unknown") for r in results]
    n_pass = classes.count("pass")
    n_fail = classes.count("fail")
    n_ambig = classes.count("ambiguous")
    n_other = len(classes) - n_pass - n_fail - n_ambig

    print(f"\nOrientation filter summary:")
    print(f"  PASS:      {n_pass:>6d} ({100*n_pass/len(classes):.1f}%)")
    print(f"  FAIL:      {n_fail:>6d} ({100*n_fail/len(classes):.1f}%)")
    print(f"  AMBIGUOUS: {n_ambig:>6d} ({100*n_ambig/len(classes):.1f}%)")
    print(f"  OTHER:     {n_other:>6d}")
    print(f"\nResults written to: {out_csv}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Beta-meander orientation filter for EGFR-MYO1D docking"
    )
    parser.add_argument("--pdb_dir", required=True,
                        help="Directory containing docked PDB files")
    parser.add_argument("--out_csv", required=True,
                        help="Output CSV path for orientation log")
    parser.add_argument("--partner_chain", default="B",
                        help="Chain ID of beta-meander partner (default: B)")
    parser.add_argument("--receptor_chain", default="A",
                        help="Chain ID of EGFR receptor (default: A)")
    parser.add_argument("--probe_residue", type=int, default=PROBE_RESIDUE,
                        help=f"Probe residue for Cα→Cβ disambiguation "
                             f"(default: {PROBE_RESIDUE})")
    parser.add_argument("--contact_dist", type=float, default=10.0,
                        help="Distance cutoff for receptor contact centroid "
                             "(default: 10.0 Å)")
    args = parser.parse_args()

    process_directory(
        args.pdb_dir, args.out_csv,
        args.partner_chain, args.receptor_chain,
        args.probe_residue, args.contact_dist
    )


if __name__ == "__main__":
    main()
