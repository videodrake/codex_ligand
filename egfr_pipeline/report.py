"""Reporting and manual review exports.

Generates a project-level summary report combining:
- Vina pocket summaries (per-receptor)
- Cross-receptor pocket comparison highlights
- Ligand-to-pocket mappings
- Auxiliary PPI residue evidence (when available)

Output: plain-text report + CSV summary tables.
PPI sections are explicitly marked as auxiliary evidence.
"""
import csv
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from egfr_pipeline.config import load_config, project_root_from_config


# ---------------------------------------------------------------------------
# CSV loaders
# ---------------------------------------------------------------------------

def load_csv(path: Path) -> List[dict]:
    if not path.exists():
        return []
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


# ---------------------------------------------------------------------------
# Report sections
# ---------------------------------------------------------------------------

def section_header(title: str, level: int = 1) -> str:
    marker = "=" * 70 if level == 1 else "-" * 50
    return f"\n{marker}\n{title}\n{marker}\n"


def format_receptor_pocket_section(
    pocket_rows: List[dict],
    drug_map_rows: List[dict],
) -> str:
    """Per-receptor pocket summary."""
    by_receptor: Dict[str, List[dict]] = defaultdict(list)
    for row in pocket_rows:
        by_receptor[row["receptor_id"]].append(row)

    drug_by_receptor: Dict[str, List[dict]] = defaultdict(list)
    for row in drug_map_rows:
        drug_by_receptor[row["receptor_id"]].append(row)

    lines = []
    for receptor_id in sorted(by_receptor):
        pockets = by_receptor[receptor_id]
        drugs = drug_by_receptor.get(receptor_id, [])

        lines.append(section_header(f"Receptor: {receptor_id}", 2))
        lines.append(f"  Pockets found: {len(pockets)}")
        lines.append(f"  Ligands docked: {len(drugs)}")
        lines.append("")

        # Pocket table
        lines.append(f"  {'Pocket':<8} {'Poses':>6} {'Ligs':>5} {'Best aff':>9} {'Mean aff':>9} {'Top residues'}")
        lines.append(f"  {'------':<8} {'-----':>6} {'----':>5} {'--------':>9} {'--------':>9} {'------------'}")
        for p in sorted(pockets, key=lambda r: float(r.get("best_affinity", 0) or 0)):
            lines.append(
                f"  {p['pocket_id']:<8} {p.get('n_pose',''):>6} {p.get('n_ligand',''):>5} "
                f"{p.get('best_affinity',''):>9} {p.get('mean_affinity',''):>9} "
                f"{p.get('top_residues','')[:40]}"
            )
        lines.append("")

        # Ligand-pocket mapping
        if drugs:
            lines.append("  Ligand -> Dominant Pocket:")
            for d in drugs:
                multi = " [multimodal]" if d.get("is_multimodal_binding", "").lower() == "true" else ""
                lines.append(
                    f"    {d['ligand_id']:<15} -> {d['dominant_pocket_id']} "
                    f"(affinity {d.get('best_affinity','')}, "
                    f"fraction {d.get('dominant_pocket_fraction','')}){multi}"
                )
            lines.append("")

    return "\n".join(lines)


def format_comparison_highlights(
    comparison_rows: List[dict],
    max_pairs: int = 20,
) -> str:
    """Cross-receptor pocket comparison highlights."""
    if not comparison_rows:
        return "  No cross-receptor comparison data available.\n"

    # Filter to same_patch_candidates first, then top by distance
    candidates = [r for r in comparison_rows if r.get("same_patch_candidate", "").lower() == "true"]
    if not candidates:
        # Show top closest pairs instead
        candidates = sorted(comparison_rows, key=lambda r: float(r.get("centroid_dist", 999)))[:max_pairs]
        header_note = "  (No same-patch candidates found. Showing closest pocket pairs.)\n"
    else:
        candidates = sorted(candidates, key=lambda r: float(r.get("centroid_dist", 999)))[:max_pairs]
        header_note = f"  Same-patch candidates: {len(candidates)} pairs\n"

    lines = [header_note]
    lines.append(
        f"  {'Rec_A':<16} {'Pkt_A':<6} {'Rec_B':<16} {'Pkt_B':<6} "
        f"{'Dist':>6} {'Jaccard':>8} {'Overlap':>8} {'Shared':>7} {'Ligs':>5} {'Candidate'}"
    )
    lines.append(
        f"  {'----':<16} {'-----':<6} {'----':<16} {'-----':<6} "
        f"{'----':>6} {'-------':>8} {'-------':>8} {'------':>7} {'----':>5} {'---------'}"
    )
    for r in candidates:
        lines.append(
            f"  {r['receptor_a']:<16} {r['pocket_a']:<6} {r['receptor_b']:<16} {r['pocket_b']:<6} "
            f"{r.get('centroid_dist',''):>6} {r.get('residue_jaccard',''):>8} "
            f"{r.get('residue_overlap_coeff',''):>8} {r.get('n_shared_residues',''):>7} "
            f"{r.get('n_shared_ligands',''):>5} {r.get('same_patch_candidate','')}"
        )
    lines.append("")

    # Summarize shared residues across all candidates
    all_shared = set()
    for r in candidates:
        shared = r.get("shared_residues", "")
        if shared:
            all_shared.update(shared.split(";"))
    if all_shared:
        lines.append(f"  Residues shared across candidate pairs: {', '.join(sorted(all_shared))}")
        lines.append("")

    return "\n".join(lines)


def format_ppi_section(
    pyrosetta_summary: List[dict],
    pyrosetta_residues: List[dict],
    afm_residues: List[dict],
) -> str:
    """Auxiliary PPI evidence section."""
    lines = []
    lines.append("  ** These are AUXILIARY evidence sources, not pocket definitions. **")
    lines.append("  ** Use alongside Vina pocket data for cross-validation only. **")
    lines.append("")

    if pyrosetta_summary:
        lines.append("  PyRosetta Global Docking Summary:")
        for s in pyrosetta_summary:
            lines.append(f"    Receptor: {s.get('receptor_id','')}")
            lines.append(f"      Models: {s.get('n_final_models','')}, Clusters: {s.get('n_clusters','')}")
            lines.append(f"      Interface residues: {s.get('n_interface_residues','')}")
            lines.append(f"      Best dG: {s.get('best_dg','')} REU, Mean dG: {s.get('mean_dg','')} REU")
            lines.append(f"      Top residues: {s.get('top_residues','')}")
            lines.append("")

    if pyrosetta_residues:
        # Show top occupancy residues per receptor
        by_receptor: Dict[str, List[dict]] = defaultdict(list)
        for r in pyrosetta_residues:
            by_receptor[r["receptor_id"]].append(r)

        lines.append("  PyRosetta Top Interface Residues (by occupancy):")
        for rec_id in sorted(by_receptor):
            residues = sorted(by_receptor[rec_id], key=lambda r: -float(r.get("occupancy", 0)))[:10]
            lines.append(f"    {rec_id}:")
            for r in residues:
                e_str = f", deltaE={r['mean_interface_delta_e']}" if r.get("mean_interface_delta_e") else ""
                lines.append(f"      {r['residue_id']:<12} occupancy={r.get('occupancy','')}{e_str}")
            lines.append("")

    if afm_residues:
        by_receptor_afm: Dict[str, List[dict]] = defaultdict(list)
        for r in afm_residues:
            by_receptor_afm[r["receptor_id"]].append(r)

        lines.append("  AlphaFold-Multimer Interface Residues:")
        for rec_id in sorted(by_receptor_afm):
            residues = sorted(by_receptor_afm[rec_id], key=lambda r: float(r.get("min_ca_distance", 999)))[:10]
            lines.append(f"    {rec_id}: {len(by_receptor_afm[rec_id])} contact residues")
            for r in residues:
                lines.append(f"      {r['residue_id']:<12} CA-CA dist={r.get('min_ca_distance','')} A")
            lines.append("")

    if not pyrosetta_summary and not pyrosetta_residues and not afm_residues:
        lines.append("  No PPI auxiliary data available.")
        lines.append("  Configure ppi.pyrosetta_result_dirs and/or ppi.afm_models in project config.")
        lines.append("")

    return "\n".join(lines)


def format_combined_residue_table(
    pocket_rows: List[dict],
    pyrosetta_residues: List[dict],
) -> List[dict]:
    """Build a combined residue evidence table: Vina pocket residues + PPI occupancy.

    This allows direct comparison of which residues are highlighted by
    both Vina docking and PPI analysis.
    """
    # Collect Vina pocket residues per receptor
    vina_residues: Dict[str, Dict[str, dict]] = defaultdict(dict)
    for pocket in pocket_rows:
        receptor_id = pocket["receptor_id"]
        pocket_id = pocket["pocket_id"]
        union = pocket.get("union_contact_residues", "")
        if not union:
            continue
        for res in union.split(";"):
            res = res.strip()
            if not res:
                continue
            # Normalize: strip chain prefix for comparison
            if ":" in res:
                norm = res.split(":", 1)[1]
            else:
                norm = res
            if norm not in vina_residues[receptor_id]:
                vina_residues[receptor_id][norm] = {
                    "receptor_id": receptor_id,
                    "residue_id": norm,
                    "vina_pockets": [],
                }
            vina_residues[receptor_id][norm]["vina_pockets"].append(pocket_id)

    # Merge PPI data
    ppi_by_receptor: Dict[str, Dict[str, dict]] = defaultdict(dict)
    for r in pyrosetta_residues:
        ppi_by_receptor[r["receptor_id"]][r["residue_id"]] = r

    # Build combined rows
    all_receptors = sorted(set(list(vina_residues.keys()) + list(ppi_by_receptor.keys())))
    combined = []
    for rec in all_receptors:
        all_res = sorted(
            set(list(vina_residues.get(rec, {}).keys()) + list(ppi_by_receptor.get(rec, {}).keys()))
        )
        for res in all_res:
            v = vina_residues.get(rec, {}).get(res, {})
            p = ppi_by_receptor.get(rec, {}).get(res, {})
            combined.append({
                "receptor_id": rec,
                "residue_id": res,
                "vina_pockets": ";".join(v.get("vina_pockets", [])),
                "n_vina_pockets": len(v.get("vina_pockets", [])),
                "ppi_occupancy": p.get("occupancy", ""),
                "ppi_frequency": p.get("frequency_final_ranking", ""),
                "ppi_delta_e": p.get("mean_interface_delta_e", ""),
                "evidence_sources": "+".join(
                    s for s in [
                        "vina" if v else "",
                        "ppi" if p else "",
                    ] if s
                ),
            })

    return combined


COMBINED_FIELDS = [
    "receptor_id",
    "residue_id",
    "vina_pockets",
    "n_vina_pockets",
    "ppi_occupancy",
    "ppi_frequency",
    "ppi_delta_e",
    "evidence_sources",
]


# ---------------------------------------------------------------------------
# Main report generator
# ---------------------------------------------------------------------------

def generate_report(
    config_path: str,
    output_dir: Optional[str] = None,
) -> Tuple[Path, Path]:
    config = load_config(config_path)
    project_root = project_root_from_config(config)
    out_root = Path(output_dir) if output_dir else project_root

    # Load all available data
    pocket_rows = load_csv(project_root / "vina_pocket_table.csv")
    drug_map_rows = load_csv(project_root / "vina_drug_pocket_map.csv")
    comparison_rows = load_csv(project_root / "vina_pocket_comparison.csv")
    pyrosetta_summary = load_csv(project_root / "ppi_pyrosetta_summary.csv")
    pyrosetta_residues = load_csv(project_root / "ppi_pyrosetta_residues.csv")
    afm_residues = load_csv(project_root / "ppi_afm_residues.csv")

    # --- Text report ---
    report_lines = []
    report_lines.append(section_header("EGFR-MYO1D Docking Analysis Report"))
    report_lines.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}")
    report_lines.append(f"Config: {config_path}")
    report_lines.append(f"Project: {config.get('project_name', 'N/A')}")
    report_lines.append("")

    # Data availability
    report_lines.append("Data sources loaded:")
    report_lines.append(f"  Vina pocket table:     {'YES' if pocket_rows else 'NO'} ({len(pocket_rows)} pockets)")
    report_lines.append(f"  Vina drug-pocket map:  {'YES' if drug_map_rows else 'NO'} ({len(drug_map_rows)} entries)")
    report_lines.append(f"  Pocket comparison:     {'YES' if comparison_rows else 'NO'} ({len(comparison_rows)} pairs)")
    report_lines.append(f"  PyRosetta PPI:         {'YES' if pyrosetta_summary else 'NO'}")
    report_lines.append(f"  AlphaFold-Multimer:    {'YES' if afm_residues else 'NO'}")
    report_lines.append("")

    # Section 1: Receptor pocket summaries
    report_lines.append(section_header("1. Receptor-Level Pocket Summary"))
    if pocket_rows:
        report_lines.append(format_receptor_pocket_section(pocket_rows, drug_map_rows))
    else:
        report_lines.append("  No Vina pocket data available.\n")

    # Section 2: Cross-receptor comparison
    report_lines.append(section_header("2. Cross-Receptor Pocket Comparison"))
    report_lines.append(format_comparison_highlights(comparison_rows))

    # Section 3: PPI auxiliary evidence
    report_lines.append(section_header("3. Auxiliary PPI Evidence"))
    report_lines.append(format_ppi_section(pyrosetta_summary, pyrosetta_residues, afm_residues))

    # Section 4: Key observations
    report_lines.append(section_header("4. Key Observations"))
    n_candidates = sum(1 for r in comparison_rows if r.get("same_patch_candidate", "").lower() == "true")
    n_receptors = len(set(r.get("receptor_id", "") for r in pocket_rows))
    n_total_pockets = len(pocket_rows)
    report_lines.append(f"  Total receptor states analyzed: {n_receptors}")
    report_lines.append(f"  Total pockets identified: {n_total_pockets}")
    report_lines.append(f"  Cross-receptor same-patch candidates: {n_candidates}")
    if pyrosetta_residues:
        high_occ = [r for r in pyrosetta_residues if float(r.get("occupancy", 0)) >= 0.5]
        report_lines.append(f"  PPI high-occupancy residues (>=0.5): {len(high_occ)}")
    report_lines.append("")
    report_lines.append("  NOTE: This report presents raw evidence for human interpretation.")
    report_lines.append("  Biological significance should be assessed by the researcher.")
    report_lines.append("")

    report_text = "\n".join(report_lines)
    report_path = out_root / "project_report.txt"
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(report_text, encoding="utf-8")

    # --- Combined residue evidence table ---
    combined = format_combined_residue_table(pocket_rows, pyrosetta_residues)
    combined_csv = out_root / "combined_residue_evidence.csv"
    if combined:
        with open(combined_csv, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=COMBINED_FIELDS)
            writer.writeheader()
            writer.writerows(combined)

    return report_path, combined_csv
