#!/usr/bin/env python3
"""Phase 1 Task 1.3: Cluster-level interface consensus.

Aggregates receptor-side interface evidence across PyRosetta clusters to
identify stable interface patches and hotspot residues. Only orientation-
validated models contribute to consensus.

Reads:
  - pyrosetta_interface_models.csv    (per-model summary with orientation columns)
  - pyrosetta_interface_residue_table.csv  (per-residue long-form)

Produces:
  - ppi_cluster_summary.csv       (cluster-level summary)
  - ppi_hotspot_residues.csv      (hotspot residues per cluster)
  - ppi_interface_patch_table.csv (per-residue occupancy across clusters)

Usage:
    python -m egfr_pipeline.phase1.cluster_consensus [--state 3GT8_raw]
    python -m egfr_pipeline.phase1.cluster_consensus --hotspot_threshold 0.4
"""

import argparse
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

PHASE1_OUTPUT_DIR = PROJECT_ROOT / "output" / "phase1_ppi"
RECEPTOR_STATES = ["3GT8_raw", "3GT8_cl38_48", "3GT8_cl85_100"]

# Hotspot threshold: a residue must appear in ≥ this fraction of
# orientation-validated models within a cluster to be a hotspot.
DEFAULT_HOTSPOT_THRESHOLD = 0.50

# N-lobe / C-lobe boundary
NLOBE_CLOBE_BOUNDARY = 838

# Output schemas
CLUSTER_SUMMARY_COLUMNS = [
    "receptor_id",
    "construct_type",
    "orientation_validation_status",
    "cluster_id",
    "n_members_total",
    "n_members_orientation_valid",
    "orientation_pass_rate",
    "representative_model",
    "representative_dG",
    "mean_dG",
    "mean_dSASA",
    "mean_sc_value",
    "mean_packstat",
    "mean_nres_int",
    "n_unique_receptor_residues",
    "n_unique_partner_residues",
    "n_nlobe_hotspots",
    "n_clobe_hotspots",
    "receptor_lobe_distribution",   # "N-lobe:X%,C-lobe:Y%"
    "receptor_hotspot_residues",    # Semicolon-separated
    "partner_hotspot_residues",     # Semicolon-separated
]

HOTSPOT_COLUMNS = [
    "receptor_id",
    "construct_type",
    "orientation_validation_status",
    "cluster_id",
    "chain",                # A (receptor) or B (partner)
    "residue_id",           # e.g., "LEU819"
    "residue_num",
    "residue_name",
    "lobe_label",           # N-lobe / C-lobe / partner
    "occupancy",            # fraction of orientation-valid models in cluster
    "n_models_with_residue",
    "n_orientation_valid_models",
    "is_hotspot",           # True if occupancy >= threshold
    "mean_delta_e_total",   # Mean per-residue energy (if available)
]

PATCH_TABLE_COLUMNS = [
    "receptor_id",
    "construct_type",
    "orientation_validation_status",
    "chain",
    "residue_id",
    "residue_num",
    "residue_name",
    "lobe_label",
    "n_clusters_present",       # How many clusters contain this residue
    "clusters_present",         # Semicolon-separated cluster IDs
    "max_occupancy",            # Highest occupancy across clusters
    "mean_occupancy",           # Mean occupancy across clusters where present
    "is_multi_cluster_hotspot", # Present as hotspot in ≥2 clusters
    "global_model_count",       # Total orientation-valid models with this residue
    "global_model_fraction",    # Fraction of all orientation-valid models
]


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_models_table(state_dir: Path) -> List[dict]:
    """Load pyrosetta_interface_models.csv."""
    path = state_dir / "pyrosetta_interface_models.csv"
    if not path.exists():
        return []
    with open(path, encoding="utf-8") as f:
        return list(csv.DictReader(f))


def load_residue_table(state_dir: Path) -> List[dict]:
    """Load pyrosetta_interface_residue_table.csv."""
    path = state_dir / "pyrosetta_interface_residue_table.csv"
    if not path.exists():
        return []
    with open(path, encoding="utf-8") as f:
        return list(csv.DictReader(f))


# ---------------------------------------------------------------------------
# Consensus computation
# ---------------------------------------------------------------------------

def compute_cluster_consensus(
    models: List[dict],
    residues: List[dict],
    receptor_id: str,
    hotspot_threshold: float = DEFAULT_HOTSPOT_THRESHOLD,
) -> Tuple[List[dict], List[dict], List[dict]]:
    """Compute cluster consensus from models and residue tables.

    Returns (cluster_summaries, hotspot_rows, patch_rows).
    """
    construct_type = _derive_construct_type(models, residues)
    orientation_filter_applied = _has_orientation_filter(models)
    orientation_validation_status = (
        "orientation_validated" if orientation_filter_applied else "not_available"
    )

    # --- Group models by cluster ---
    cluster_models = defaultdict(list)
    for m in models:
        cid = m.get("cluster_id", "")
        if cid:
            cluster_models[cid].append(m)

    # --- Index residues by model_id ---
    residues_by_model = defaultdict(list)
    for r in residues:
        residues_by_model[r["model_id"]].append(r)

    # --- Build per-cluster analysis ---
    cluster_summaries = []
    all_hotspot_rows = []
    # Track global residue presence for patch table
    residue_cluster_info = defaultdict(lambda: {
        "clusters": set(), "occupancies": {}, "global_count": 0,
        "chain": "", "residue_num": "", "residue_name": "", "lobe_label": "",
    })

    # Count total orientation-valid models for global fraction
    total_orient_valid = _count_consensus_models(models, orientation_filter_applied)

    for cid in sorted(cluster_models.keys()):
        members = cluster_models[cid]
        n_total = len(members)

        # Filter to orientation-validated
        valid_members = _select_consensus_members(members, orientation_filter_applied)
        n_valid = len(valid_members)
        pass_rate = n_valid / n_total if n_total > 0 else 0.0

        # Representative model: best dG among orientation-valid models
        representative = None
        representative_dG = None
        if valid_members:
            for m in valid_members:
                dg = _safe_float(m.get("dG_separated", ""))
                if dg is not None:
                    if representative_dG is None or dg < representative_dG:
                        representative = m.get("model_id", "")
                        representative_dG = dg
            if representative is None:
                representative = valid_members[0].get("model_id", "")

        # Mean metrics (orientation-valid only)
        mean_dG = _mean_field(valid_members, "dG_separated")
        mean_dSASA = _mean_field(valid_members, "dSASA")
        mean_sc = _mean_field(valid_members, "sc_value")
        mean_packstat = _mean_field(valid_members, "packstat")
        mean_nres = _mean_field(valid_members, "nres_int")

        # --- Per-residue occupancy within cluster (orientation-valid only) ---
        valid_model_ids = {m["model_id"] for m in valid_members}
        residue_counts = Counter()  # residue_id → count of valid models
        residue_meta = {}           # residue_id → metadata

        for mid in valid_model_ids:
            seen_in_model = set()
            for r in residues_by_model.get(mid, []):
                rid = r.get("residue_id", "")
                if rid and rid not in seen_in_model:
                    seen_in_model.add(rid)
                    key = (r.get("chain", ""), rid)
                    residue_counts[key] += 1
                    if key not in residue_meta:
                        residue_meta[key] = {
                            "chain": r.get("chain", ""),
                            "residue_id": rid,
                            "residue_num": r.get("residue_num", ""),
                            "residue_name": r.get("residue_name", ""),
                            "lobe_label": r.get("lobe_label", ""),
                        }

        # Compute per-residue energy means for orientation-valid models
        residue_energies = defaultdict(list)
        for mid in valid_model_ids:
            for r in residues_by_model.get(mid, []):
                key = (r.get("chain", ""), r.get("residue_id", ""))
                de = _safe_float(r.get("delta_e_total", ""))
                if de is not None:
                    residue_energies[key].append(de)

        # Build hotspot rows
        receptor_hotspots = []
        partner_hotspots = []
        n_nlobe_hotspot = 0
        n_clobe_hotspot = 0

        unique_receptor = set()
        unique_partner = set()

        for key, count in sorted(residue_counts.items(), key=lambda x: -x[1]):
            chain, rid = key
            meta = residue_meta.get(key, {})
            occupancy = count / n_valid if n_valid > 0 else 0.0
            is_hotspot = occupancy >= hotspot_threshold

            if chain == "A":
                unique_receptor.add(rid)
            elif chain == "B":
                unique_partner.add(rid)

            energies = residue_energies.get(key, [])
            mean_de = sum(energies) / len(energies) if energies else None

            hotspot_row = {
                "receptor_id": receptor_id,
                "construct_type": construct_type,
                "orientation_validation_status": orientation_validation_status,
                "cluster_id": cid,
                "chain": chain,
                "residue_id": rid,
                "residue_num": meta.get("residue_num", ""),
                "residue_name": meta.get("residue_name", ""),
                "lobe_label": meta.get("lobe_label", ""),
                "occupancy": round(occupancy, 4),
                "n_models_with_residue": count,
                "n_orientation_valid_models": n_valid,
                "is_hotspot": is_hotspot,
                "mean_delta_e_total": round(mean_de, 3) if mean_de is not None else "",
            }
            all_hotspot_rows.append(hotspot_row)

            if is_hotspot:
                if chain == "A":
                    receptor_hotspots.append(rid)
                    lobe = meta.get("lobe_label", "")
                    if lobe == "N-lobe":
                        n_nlobe_hotspot += 1
                    elif lobe == "C-lobe":
                        n_clobe_hotspot += 1
                elif chain == "B":
                    partner_hotspots.append(rid)

            # Track for global patch table
            rinfo = residue_cluster_info[(chain, rid)]
            rinfo["clusters"].add(cid)
            rinfo["occupancies"][cid] = occupancy
            rinfo["global_count"] += count
            rinfo["chain"] = chain
            rinfo["residue_num"] = meta.get("residue_num", "")
            rinfo["residue_name"] = meta.get("residue_name", "")
            rinfo["lobe_label"] = meta.get("lobe_label", "")

        # Lobe distribution
        n_nlobe_res = sum(1 for (ch, _), m in residue_meta.items()
                         if ch == "A" and m.get("lobe_label") == "N-lobe")
        n_clobe_res = sum(1 for (ch, _), m in residue_meta.items()
                         if ch == "A" and m.get("lobe_label") == "C-lobe")
        total_receptor_res = n_nlobe_res + n_clobe_res
        if total_receptor_res > 0:
            nlobe_pct = round(100 * n_nlobe_res / total_receptor_res, 1)
            clobe_pct = round(100 * n_clobe_res / total_receptor_res, 1)
            lobe_dist = f"N-lobe:{nlobe_pct}%,C-lobe:{clobe_pct}%"
        else:
            lobe_dist = ""

        cluster_summaries.append({
            "receptor_id": receptor_id,
            "construct_type": construct_type,
            "orientation_validation_status": orientation_validation_status,
            "cluster_id": cid,
            "n_members_total": n_total,
            "n_members_orientation_valid": n_valid,
            "orientation_pass_rate": round(pass_rate, 4),
            "representative_model": representative or "",
            "representative_dG": round(representative_dG, 2) if representative_dG is not None else "",
            "mean_dG": round(mean_dG, 2) if mean_dG is not None else "",
            "mean_dSASA": round(mean_dSASA, 1) if mean_dSASA is not None else "",
            "mean_sc_value": round(mean_sc, 4) if mean_sc is not None else "",
            "mean_packstat": round(mean_packstat, 4) if mean_packstat is not None else "",
            "mean_nres_int": round(mean_nres, 1) if mean_nres is not None else "",
            "n_unique_receptor_residues": len(unique_receptor),
            "n_unique_partner_residues": len(unique_partner),
            "n_nlobe_hotspots": n_nlobe_hotspot,
            "n_clobe_hotspots": n_clobe_hotspot,
            "receptor_lobe_distribution": lobe_dist,
            "receptor_hotspot_residues": ";".join(receptor_hotspots),
            "partner_hotspot_residues": ";".join(partner_hotspots),
        })

    # --- Build global patch table ---
    patch_rows = []
    for (chain, rid), info in sorted(residue_cluster_info.items()):
        n_clusters = len(info["clusters"])
        occupancies = list(info["occupancies"].values())
        max_occ = max(occupancies) if occupancies else 0
        mean_occ = sum(occupancies) / len(occupancies) if occupancies else 0

        # Is hotspot in ≥2 clusters?
        n_hotspot_clusters = sum(1 for o in occupancies if o >= hotspot_threshold)
        is_multi = n_hotspot_clusters >= 2

        global_frac = (info["global_count"] / total_orient_valid
                       if total_orient_valid > 0 else 0)

        patch_rows.append({
            "receptor_id": receptor_id,
            "construct_type": construct_type,
            "orientation_validation_status": orientation_validation_status,
            "chain": chain,
            "residue_id": rid,
            "residue_num": info["residue_num"],
            "residue_name": info["residue_name"],
            "lobe_label": info["lobe_label"],
            "n_clusters_present": n_clusters,
            "clusters_present": ";".join(sorted(info["clusters"])),
            "max_occupancy": round(max_occ, 4),
            "mean_occupancy": round(mean_occ, 4),
            "is_multi_cluster_hotspot": is_multi,
            "global_model_count": info["global_count"],
            "global_model_fraction": round(global_frac, 4),
        })

    # Sort patch table: receptor first, then by global_model_fraction descending
    patch_rows.sort(key=lambda r: (
        0 if r["chain"] == "A" else 1,
        -r.get("global_model_fraction", 0),
    ))

    return cluster_summaries, all_hotspot_rows, patch_rows


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _safe_float(val: str) -> Optional[float]:
    """Convert string to float, returning None on failure."""
    try:
        return float(val)
    except (ValueError, TypeError):
        return None


def _has_orientation_filter(models: List[dict]) -> bool:
    """Return True when non-empty orientation labels are available."""
    if not models:
        return False
    if "orientation_class" not in models[0]:
        return False
    return any((m.get("orientation_class", "") or "").strip() for m in models)


def _select_consensus_members(
    members: List[dict],
    orientation_filter_applied: bool,
) -> List[dict]:
    """Select models that contribute to consensus.

    If orientation filtering is available, only pass models contribute.
    Otherwise, all members contribute as a backward-compatible fallback.
    """
    if not orientation_filter_applied:
        return list(members)
    return [m for m in members if m.get("orientation_class", "") == "pass"]


def _count_consensus_models(
    models: List[dict],
    orientation_filter_applied: bool,
) -> int:
    """Count the models contributing to occupancy denominators."""
    if not orientation_filter_applied:
        return len(models)
    return sum(1 for m in models if m.get("orientation_class", "") == "pass")


def _derive_construct_type(models: List[dict], residues: List[dict]) -> str:
    """Derive construct_type from upstream tables without hiding mismatches."""
    values = {
        row.get("construct_type", "").strip()
        for row in [*models, *residues]
        if row.get("construct_type", "").strip()
    }
    if not values:
        return ""
    if len(values) == 1:
        return next(iter(values))
    return "mixed_construct_type"


def _mean_field(rows: List[dict], field: str) -> Optional[float]:
    """Compute mean of a numeric field across rows."""
    values = []
    for r in rows:
        v = _safe_float(r.get(field, ""))
        if v is not None:
            values.append(v)
    return sum(values) / len(values) if values else None


def _write_csv(path: Path, rows: List[dict], columns: List[str]) -> None:
    """Write rows to CSV with specified columns."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


# ---------------------------------------------------------------------------
# State-level processing
# ---------------------------------------------------------------------------

def process_state(
    state_name: str,
    output_base: Path,
    hotspot_threshold: float = DEFAULT_HOTSPOT_THRESHOLD,
) -> bool:
    """Process cluster consensus for one receptor state.

    Returns True if outputs were generated.
    """
    state_dir = output_base / state_name

    models = load_models_table(state_dir)
    residues = load_residue_table(state_dir)

    if not models:
        print(f"  No models table for {state_name}")
        return False
    if not residues:
        print(f"  No residue table for {state_name}")
        return False

    orientation_filter_applied = _has_orientation_filter(models)

    # Check orientation labels exist
    if not orientation_filter_applied:
        if "orientation_class" not in models[0]:
            print(f"  WARNING: orientation_class not in models table for {state_name}")
        else:
            print(
                f"  WARNING: orientation_class labels missing/empty in models table for {state_name}"
            )
        print(f"  Run orientation filter with --merge first.")
        print(f"  Proceeding without orientation filtering (all models used).")

    summaries, hotspots, patches = compute_cluster_consensus(
        models, residues, state_name, hotspot_threshold
    )

    if not summaries:
        print(f"  No clusters found for {state_name}")
        return False

    # Write outputs
    summary_path = state_dir / "ppi_cluster_summary.csv"
    hotspot_path = state_dir / "ppi_hotspot_residues.csv"
    patch_path = state_dir / "ppi_interface_patch_table.csv"

    _write_csv(summary_path, summaries, CLUSTER_SUMMARY_COLUMNS)
    _write_csv(hotspot_path, hotspots, HOTSPOT_COLUMNS)
    _write_csv(patch_path, patches, PATCH_TABLE_COLUMNS)

    # Print summary
    n_clusters = len(summaries)
    total_models = sum(s.get("n_members_total", 0) for s in summaries)
    total_valid = sum(s.get("n_members_orientation_valid", 0) for s in summaries)
    n_hotspot_res = sum(1 for h in hotspots if h.get("is_hotspot") and h.get("chain") == "A")
    n_multi_patch = sum(1 for p in patches if p.get("is_multi_cluster_hotspot"))

    print(f"\n  Cluster consensus for {state_name}:")
    print(f"    Clusters: {n_clusters}")
    print(f"    Models: {total_models} total, {total_valid} orientation-valid")
    print(f"    Receptor hotspot residues: {n_hotspot_res} (≥{hotspot_threshold*100:.0f}% occupancy)")
    print(f"    Multi-cluster hotspots: {n_multi_patch}")
    print(f"    Outputs:")
    print(f"      {summary_path.name}")
    print(f"      {hotspot_path.name}")
    print(f"      {patch_path.name}")

    # Per-cluster detail
    for s in summaries:
        cid = s["cluster_id"]
        nv = s["n_members_orientation_valid"]
        nt = s["n_members_total"]
        mdg = s.get("mean_dG", "")
        lobe = s.get("receptor_lobe_distribution", "")
        rhs = s.get("receptor_hotspot_residues", "")
        n_rh = len(rhs.split(";")) if rhs else 0
        print(f"    {cid}: {nv}/{nt} valid, mean_dG={mdg}, lobe=[{lobe}], "
              f"receptor_hotspots={n_rh}")

    return True


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Phase 1 TG 1.3: Cluster-level interface consensus"
    )
    parser.add_argument("--state", choices=RECEPTOR_STATES,
                        help="Process a single receptor state")
    parser.add_argument("--output_dir", type=Path,
                        default=PHASE1_OUTPUT_DIR,
                        help="Phase 1 output base directory")
    parser.add_argument("--hotspot_threshold", type=float,
                        default=DEFAULT_HOTSPOT_THRESHOLD,
                        help=f"Minimum occupancy for hotspot designation "
                             f"(default: {DEFAULT_HOTSPOT_THRESHOLD})")
    args = parser.parse_args()

    print("Phase 1 — Task 1.3: Cluster-Level Interface Consensus")
    print(f"Output base: {args.output_dir}")
    print(f"Hotspot threshold: {args.hotspot_threshold*100:.0f}% occupancy")
    print(f"N-lobe/C-lobe boundary: residue {NLOBE_CLOBE_BOUNDARY}")

    states = [args.state] if args.state else RECEPTOR_STATES
    results = []

    for state in states:
        print(f"\nProcessing {state}:")
        success = process_state(state, args.output_dir, args.hotspot_threshold)
        if success:
            results.append(state)

    if not results:
        print("\nNo consensus data generated. Run extraction + orientation filter first.")
        return 0

    print(f"\n{'='*60}")
    print("CLUSTER CONSENSUS COMPLETE")
    print(f"{'='*60}")
    for state in results:
        print(f"  {state}: ppi_cluster_summary.csv, ppi_hotspot_residues.csv, ppi_interface_patch_table.csv")

    return 0


if __name__ == "__main__":
    sys.exit(main())
