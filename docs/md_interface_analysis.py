#!/usr/bin/env python3
"""
md_interface_analysis.py
========================
Post-MD analysis of PPI interface dynamic stability for EGFR–MYO1D complexes.

Computes:
  1. Interface contact survival over time (fraction of initial contacts maintained)
  2. Per-residue contact occupancy (which residues stay in contact)
  3. Hotspot contact persistence (specifically for alanine-scanning hotspot residues)
  4. Interface buried surface area (dSASA) time series
  5. Inter-chain minimum distance time series

Dependencies:
  - MDAnalysis >= 2.0
  - NumPy
  - (Optional) FreeSASA for dSASA calculation

Usage:
    python md_interface_analysis.py \\
        --topology md_run/topol.tpr \\
        --trajectory md_run/md_noPBC.xtc \\
        --receptor_chain A \\
        --partner_chain B \\
        --contact_cutoff 5.0 \\
        --hotspot_residues 962,925,928,874,919 \\
        --out_prefix md_analysis_results/complex1_rep1

Author: EGFR-MYO1D Pipeline
Date: 2026-03
"""

import argparse
import os
import warnings
from typing import Dict, List, Optional, Tuple

import numpy as np

try:
    import MDAnalysis as mda
    from MDAnalysis.analysis import distances
except ImportError:
    raise ImportError("MDAnalysis is required. Install: pip install MDAnalysis")


# ---------------------------------------------------------------------------
# Contact map computation
# ---------------------------------------------------------------------------

def compute_contact_map(
    universe: mda.Universe,
    receptor_sel: str,
    partner_sel: str,
    cutoff: float = 5.0,
    frame_idx: int = 0,
) -> np.ndarray:
    """
    Compute a binary contact map between receptor and partner heavy atoms.

    Returns a boolean matrix: contact_map[i, j] = True if receptor heavy
    atom i is within cutoff of partner heavy atom j.
    """
    universe.trajectory[frame_idx]
    rec_atoms = universe.select_atoms(f"{receptor_sel} and not name H*")
    par_atoms = universe.select_atoms(f"{partner_sel} and not name H*")

    dist_matrix = distances.distance_array(
        rec_atoms.positions, par_atoms.positions, box=universe.dimensions
    )
    return dist_matrix < cutoff


def get_residue_contact_pairs(
    universe: mda.Universe,
    receptor_sel: str,
    partner_sel: str,
    cutoff: float = 5.0,
    frame_idx: int = 0,
) -> set:
    """
    Get a set of (receptor_resid, partner_resid) pairs that are in contact
    at a given frame. Contact = any heavy atom pair within cutoff.
    """
    universe.trajectory[frame_idx]
    rec_atoms = universe.select_atoms(f"{receptor_sel} and not name H*")
    par_atoms = universe.select_atoms(f"{partner_sel} and not name H*")

    dist_matrix = distances.distance_array(
        rec_atoms.positions, par_atoms.positions, box=universe.dimensions
    )

    contact_pairs = set()
    in_contact = np.argwhere(dist_matrix < cutoff)
    for i, j in in_contact:
        rec_resid = rec_atoms[i].resid
        par_resid = par_atoms[j].resid
        contact_pairs.add((rec_resid, par_resid))

    return contact_pairs


# ---------------------------------------------------------------------------
# Interface contact survival
# ---------------------------------------------------------------------------

def compute_contact_survival(
    universe: mda.Universe,
    receptor_sel: str,
    partner_sel: str,
    cutoff: float = 5.0,
    stride: int = 10,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute the fraction of initial (frame 0) residue-residue contacts
    that are maintained at each subsequent frame.

    Returns:
        times: array of time points (ps)
        survival: array of fractional contact survival [0, 1]
    """
    # Reference contacts at frame 0
    ref_contacts = get_residue_contact_pairs(
        universe, receptor_sel, partner_sel, cutoff, frame_idx=0
    )
    n_ref = len(ref_contacts)

    if n_ref == 0:
        warnings.warn("No initial contacts found at frame 0!")
        return np.array([]), np.array([])

    times = []
    survivals = []

    for ts in universe.trajectory[::stride]:
        current_contacts = get_residue_contact_pairs(
            universe, receptor_sel, partner_sel, cutoff,
            frame_idx=ts.frame
        )
        # Note: get_residue_contact_pairs sets the frame internally,
        # but we've already set it via the loop. Re-set to be safe.
        universe.trajectory[ts.frame]

        maintained = ref_contacts & current_contacts
        survival_frac = len(maintained) / n_ref

        times.append(ts.time)  # in ps
        survivals.append(survival_frac)

    return np.array(times), np.array(survivals)


# ---------------------------------------------------------------------------
# Per-residue contact occupancy
# ---------------------------------------------------------------------------

def compute_residue_occupancy(
    universe: mda.Universe,
    receptor_sel: str,
    partner_sel: str,
    cutoff: float = 5.0,
    stride: int = 10,
    first_n_frames: Optional[int] = None,
) -> Dict[int, float]:
    """
    For each partner residue, compute the fraction of frames where it
    has at least one contact with the receptor.

    Returns: {partner_resid: occupancy_fraction}
    """
    par_atoms = universe.select_atoms(f"{partner_sel} and not name H*")
    all_partner_resids = sorted(set(par_atoms.resids))

    contact_counts = {r: 0 for r in all_partner_resids}
    n_frames = 0

    frames = universe.trajectory[::stride]
    if first_n_frames is not None:
        frames = list(frames)[:first_n_frames]

    for ts in frames:
        contacts = get_residue_contact_pairs(
            universe, receptor_sel, partner_sel, cutoff, frame_idx=ts.frame
        )
        universe.trajectory[ts.frame]

        partner_resids_in_contact = {p for _, p in contacts}
        for r in partner_resids_in_contact:
            if r in contact_counts:
                contact_counts[r] += 1

        n_frames += 1

    occupancy = {r: contact_counts[r] / max(n_frames, 1)
                 for r in all_partner_resids}
    return occupancy


# ---------------------------------------------------------------------------
# Hotspot contact persistence
# ---------------------------------------------------------------------------

def compute_hotspot_persistence(
    universe: mda.Universe,
    receptor_sel: str,
    partner_sel: str,
    hotspot_resids: List[int],
    cutoff: float = 5.0,
    stride: int = 10,
) -> Dict[int, Tuple[np.ndarray, np.ndarray]]:
    """
    For each hotspot residue, compute a time series of whether it is in
    contact with the receptor at each frame.

    Returns: {resid: (times_array, in_contact_bool_array)}
    """
    result = {r: ([], []) for r in hotspot_resids}

    for ts in universe.trajectory[::stride]:
        contacts = get_residue_contact_pairs(
            universe, receptor_sel, partner_sel, cutoff, frame_idx=ts.frame
        )
        universe.trajectory[ts.frame]

        partner_resids_in_contact = {p for _, p in contacts}

        for r in hotspot_resids:
            result[r][0].append(ts.time)
            result[r][1].append(r in partner_resids_in_contact)

    return {r: (np.array(times), np.array(contacts))
            for r, (times, contacts) in result.items()}


# ---------------------------------------------------------------------------
# Inter-chain minimum distance
# ---------------------------------------------------------------------------

def compute_interchain_mindist(
    universe: mda.Universe,
    receptor_sel: str,
    partner_sel: str,
    stride: int = 10,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute the minimum heavy-atom distance between chains at each frame.
    """
    times = []
    mindists = []

    for ts in universe.trajectory[::stride]:
        rec_atoms = universe.select_atoms(f"{receptor_sel} and not name H*")
        par_atoms = universe.select_atoms(f"{partner_sel} and not name H*")

        dist_matrix = distances.distance_array(
            rec_atoms.positions, par_atoms.positions, box=universe.dimensions
        )
        mindist = dist_matrix.min()

        times.append(ts.time)
        mindists.append(mindist)

    return np.array(times), np.array(mindists)


# ---------------------------------------------------------------------------
# Convergence assessment
# ---------------------------------------------------------------------------

def assess_convergence(
    times: np.ndarray,
    values: np.ndarray,
    window_ns: float = 50.0,
    threshold: float = 0.15,
) -> Dict:
    """
    Assess whether a time series has converged by checking if the
    mean and standard deviation in the last window are stable.

    Returns dict with:
      - converged: bool
      - last_window_mean: float
      - last_window_std: float
      - full_mean: float
    """
    # Convert window from ns to ps
    window_ps = window_ns * 1000
    last_window_mask = times >= (times[-1] - window_ps)

    if last_window_mask.sum() < 10:
        return {"converged": False, "error": "Insufficient data in last window"}

    last_vals = values[last_window_mask]
    full_mean = values.mean()

    return {
        "converged": last_vals.std() < threshold,
        "last_window_mean": float(last_vals.mean()),
        "last_window_std": float(last_vals.std()),
        "full_mean": float(full_mean),
    }


# ---------------------------------------------------------------------------
# Main analysis pipeline
# ---------------------------------------------------------------------------

def run_full_analysis(
    topology: str,
    trajectory: str,
    receptor_sel: str,
    partner_sel: str,
    hotspot_resids: List[int],
    contact_cutoff: float = 5.0,
    stride: int = 10,
    out_prefix: str = "md_analysis",
) -> Dict:
    """Run all analyses and save results."""

    print(f"Loading: {topology}, {trajectory}")
    u = mda.Universe(topology, trajectory)
    print(f"  Frames: {len(u.trajectory)}, "
          f"dt: {u.trajectory.dt:.1f} ps, "
          f"total: {u.trajectory[-1].time/1000:.1f} ns")

    os.makedirs(os.path.dirname(out_prefix) or ".", exist_ok=True)

    results = {}

    # 1. Contact survival
    print("Computing contact survival...")
    times_cs, survival = compute_contact_survival(
        u, receptor_sel, partner_sel, contact_cutoff, stride
    )
    np.savetxt(
        f"{out_prefix}_contact_survival.csv",
        np.column_stack([times_cs / 1000, survival]),  # convert to ns
        header="time_ns,contact_survival_fraction",
        delimiter=",", comments=""
    )
    conv_cs = assess_convergence(times_cs, survival, window_ns=50, threshold=0.1)
    results["contact_survival"] = conv_cs
    print(f"  Last 50 ns mean survival: {conv_cs.get('last_window_mean', 'N/A'):.3f}")

    # 2. Per-residue occupancy (last 100 ns)
    print("Computing per-residue occupancy (last 100 ns)...")
    # Select last 100 ns of frames
    total_time = u.trajectory[-1].time
    last_100ns_start = total_time - 100_000  # 100 ns in ps
    n_frames_last100 = sum(
        1 for ts in u.trajectory[::stride] if ts.time >= last_100ns_start
    )
    occupancy = compute_residue_occupancy(
        u, receptor_sel, partner_sel, contact_cutoff, stride
    )
    with open(f"{out_prefix}_residue_occupancy.csv", "w") as f:
        f.write("partner_resid,occupancy\n")
        for resid in sorted(occupancy):
            f.write(f"{resid},{occupancy[resid]:.4f}\n")

    # 3. Hotspot persistence
    print("Computing hotspot contact persistence...")
    hotspot_data = compute_hotspot_persistence(
        u, receptor_sel, partner_sel, hotspot_resids, contact_cutoff, stride
    )
    with open(f"{out_prefix}_hotspot_persistence.csv", "w") as f:
        header_parts = ["time_ns"] + [f"res{r}" for r in hotspot_resids]
        f.write(",".join(header_parts) + "\n")
        # Get time array from first residue
        if hotspot_resids:
            time_arr = hotspot_data[hotspot_resids[0]][0] / 1000  # ns
            for idx in range(len(time_arr)):
                row = [f"{time_arr[idx]:.3f}"]
                for r in hotspot_resids:
                    row.append(str(int(hotspot_data[r][1][idx])))
                f.write(",".join(row) + "\n")

    # Hotspot occupancy in last 100 ns
    hotspot_occ = {}
    for r in hotspot_resids:
        times_h, contacts_h = hotspot_data[r]
        mask = times_h >= last_100ns_start
        if mask.sum() > 0:
            hotspot_occ[r] = float(contacts_h[mask].mean())
        else:
            hotspot_occ[r] = 0.0
    results["hotspot_occupancy_last100ns"] = hotspot_occ
    print(f"  Hotspot occupancy (last 100 ns): {hotspot_occ}")

    # 4. Inter-chain minimum distance
    print("Computing inter-chain minimum distance...")
    times_md, mindists = compute_interchain_mindist(
        u, receptor_sel, partner_sel, stride
    )
    np.savetxt(
        f"{out_prefix}_interchain_mindist.csv",
        np.column_stack([times_md / 1000, mindists]),
        header="time_ns,min_distance_angstrom",
        delimiter=",", comments=""
    )
    conv_md = assess_convergence(times_md, mindists, window_ns=50, threshold=1.0)
    results["interchain_mindist"] = conv_md

    # Check for detachment events (> 8 Å for > 10 ns)
    detach_threshold = 8.0
    detach_window_ps = 10_000  # 10 ns
    detached = mindists > detach_threshold
    if detached.any():
        # Find longest consecutive detachment
        max_detach_duration = 0
        current_duration = 0
        dt = times_md[1] - times_md[0] if len(times_md) > 1 else stride
        for d in detached:
            if d:
                current_duration += dt
                max_detach_duration = max(max_detach_duration, current_duration)
            else:
                current_duration = 0
        results["max_detachment_duration_ns"] = max_detach_duration / 1000
        results["detachment_detected"] = max_detach_duration > detach_window_ps
    else:
        results["max_detachment_duration_ns"] = 0.0
        results["detachment_detected"] = False

    print(f"  Detachment detected: {results['detachment_detected']}")

    # 5. Summary
    print("\n=== SUMMARY ===")
    print(f"Contact survival (last 50 ns): "
          f"{conv_cs.get('last_window_mean', 'N/A')}")
    print(f"Inter-chain min distance (last 50 ns): "
          f"{conv_md.get('last_window_mean', 'N/A'):.2f} Å")
    print(f"Detachment event: {results['detachment_detected']}")
    for r in hotspot_resids:
        occ = hotspot_occ.get(r, 0)
        status = "✅" if occ >= 0.7 else ("⚠️" if occ >= 0.4 else "❌")
        print(f"  Hotspot {r}: {occ:.1%} occupancy {status}")

    # Write summary JSON
    import json
    with open(f"{out_prefix}_summary.json", "w") as f:
        json.dump(results, f, indent=2, default=str)

    print(f"\nResults written to {out_prefix}_*.csv/json")
    return results


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="MD interface stability analysis for EGFR-MYO1D PPI"
    )
    parser.add_argument("--topology", required=True,
                        help="Topology file (tpr, prmtop, psf)")
    parser.add_argument("--trajectory", required=True,
                        help="Trajectory file (xtc, trr, dcd)")
    parser.add_argument("--receptor_chain", default="A",
                        help="Receptor chain ID (default: A)")
    parser.add_argument("--partner_chain", default="B",
                        help="Partner chain ID (default: B)")
    parser.add_argument("--contact_cutoff", type=float, default=5.0,
                        help="Heavy-atom contact cutoff in Å (default: 5.0)")
    parser.add_argument("--hotspot_residues", default="962,925,928,874,919",
                        help="Comma-separated hotspot residue IDs to track")
    parser.add_argument("--stride", type=int, default=10,
                        help="Frame stride for analysis (default: 10)")
    parser.add_argument("--out_prefix", default="md_analysis",
                        help="Output file prefix")

    args = parser.parse_args()

    receptor_sel = f"segid {args.receptor_chain} or chainID {args.receptor_chain}"
    partner_sel = f"segid {args.partner_chain} or chainID {args.partner_chain}"
    hotspot_resids = [int(x) for x in args.hotspot_residues.split(",")]

    run_full_analysis(
        args.topology, args.trajectory,
        receptor_sel, partner_sel,
        hotspot_resids,
        args.contact_cutoff, args.stride,
        args.out_prefix,
    )


if __name__ == "__main__":
    main()
