"""Organize production outputs into step-based directories.

Creates a clean browsing structure under ``output/{project}/steps/`` with
symbolic links pointing back to the canonical output locations.  This avoids
data duplication while providing a single, intuitive entry point for analysis.

Usage::

    python main.py organize
    python main.py -c config/example-project.yaml organize
"""

from __future__ import annotations

import os
import shutil
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from egfr_pipeline.config import load_config, project_root_from_config

# ---------------------------------------------------------------------------
# Step metadata: (folder_name, display_title, description)
# ---------------------------------------------------------------------------

STEPS: Dict[int, Tuple[str, str, str]] = {
    1: (
        "01_vina_docking",
        "Vina Blind Docking",
        "Raw Vina docking poses (.pdbqt) for each receptor x ligand pair.",
    ),
    2: (
        "02_ppi_docking",
        "PPI Global Blind Docking",
        (
            "PyRosetta PPI docking results per receptor state:\n"
            "  ranked poses, clusters, energy funnels, validation reports, PyMOL scripts.\n"
            "\n"
            "  Key files per state:\n"
            "    final_ranking.csv            Overall ranked models\n"
            "    scored_all_models.csv        All models with metrics + filter_status\n"
            "    scored_stage2_models.csv     Stage 2 expensive metrics (packstat, unsatHb)\n"
            "    energy_funnel.png            L_RMSD vs dG scatter plot\n"
            "    energy_funnel_best.png       Best-cluster energy funnel\n"
            "    pyrosetta_run_metadata.json  Run parameters and timing\n"
            "    1_OVERVIEW_Clusters.pml      PyMOL: all binding sites overview\n"
            "    final_result/                Final structures + per-residue energies\n"
            "      dashboard.html            Interactive results dashboard\n"
            "    cluster_results/             Cluster summaries + membership\n"
            "    logs/                        Pipeline & worker logs\n"
            "    docking_validation_report    (inside final_result/)"
        ),
    ),
    3: (
        "03_ppi_evidence",
        "PPI Interface Evidence",
        "Interface residue evidence tables extracted from PPI docking results.",
    ),
    4: (
        "04_vina_analysis",
        "Vina Pocket Analysis",
        (
            "Pocket clustering, drug-pocket mapping, cross-receptor comparison,\n"
            "  bootstrap stability analysis."
        ),
    ),
    5: (
        "05_site_verdict",
        "Site Verdict",
        "Multi-method evidence integration and site prioritization.",
    ),
    6: (
        "06_report",
        "Report",
        "Narrative project report and combined residue evidence.",
    ),
    7: (
        "07_validation",
        "Validation",
        "Output quality checks and recovery guidance.",
    ),
    8: (
        "08_pocket_analysis",
        "Pocket Analysis (Phase 2)",
        "Druggability assessment, candidate pockets, cross-state alignment.",
    ),
    9: (
        "09_diverse_docking",
        "Diverse Docking (Phase 3)",
        "Expanded docking campaigns with diversity metrics and pocket search status.",
    ),
    10: (
        "10_perturbation",
        "Perturbation Analysis (Phase 4)",
        "Perturbation evidence, normalized scores, final candidate classes.",
    ),
}


# ---------------------------------------------------------------------------
# Symlink helpers
# ---------------------------------------------------------------------------


def _safe_symlink(target: Path, link: Path) -> bool:
    """Create a relative symlink if *target* exists.  Returns True on success."""
    target_abs = target.resolve() if target.exists() else target
    if not target_abs.exists():
        return False
    link.parent.mkdir(parents=True, exist_ok=True)
    if link.exists() or link.is_symlink():
        link.unlink()
    # Compute relative path from link location to target
    link_parent_abs = link.parent.resolve()
    rel = os.path.relpath(target_abs, link_parent_abs)
    link.symlink_to(rel)
    return True


def _link_dir(src: Path, link: Path) -> bool:
    """Symlink an entire directory.  Returns True on success."""
    src_abs = src.resolve() if src.exists() else src
    if not src_abs.is_dir():
        return False
    link.parent.mkdir(parents=True, exist_ok=True)
    if link.exists() or link.is_symlink():
        if link.is_symlink():
            link.unlink()
        else:
            shutil.rmtree(link)
    link_parent_abs = link.parent.resolve()
    rel = os.path.relpath(src_abs, link_parent_abs)
    link.symlink_to(rel)
    return True


# ---------------------------------------------------------------------------
# Per-step organizers
# ---------------------------------------------------------------------------


def _step1_vina(
    project_root: Path,
    config: dict,
    step_dir: Path,
) -> List[Tuple[str, str]]:
    """Step 1 — symlink each receptor's Vina output directory."""
    mode = str(
        config.get("mode", (config.get("vina") or {}).get("mode", "blind"))
    )
    linked: List[Tuple[str, str]] = []
    for receptor in config.get("receptors", []):
        rid = str(receptor.get("id", ""))
        if not rid:
            continue
        rec_dir = project_root / rid
        if _link_dir(rec_dir, step_dir / rid):
            n = len(list(rec_dir.glob(f"*_{mode}.pdbqt")))
            linked.append((f"{rid}/", f"{n} pose file(s)"))
    return linked


def _discover_ppi_dirs(repo_root: Path) -> Dict[str, Path]:
    """Auto-discover PPI run directories under ``output/phase1_ppi/``."""
    ppi_base = repo_root / "output" / "phase1_ppi"
    result: Dict[str, Path] = {}
    if not ppi_base.is_dir():
        return result
    for state_dir in sorted(ppi_base.iterdir()):
        if not state_dir.is_dir() or state_dir.name.startswith("."):
            continue
        for run_dir in sorted(state_dir.iterdir()):
            if not run_dir.is_dir():
                continue
            has_ranking = (
                (run_dir / "final_result" / "final_ranking.csv").exists()
                or (run_dir / "final_ranking.csv").exists()
            )
            has_csv = any(run_dir.glob("*.csv"))
            if has_ranking or has_csv:
                # Use state name as key; prefer prod over test,
                # lowest seed number among prod dirs.
                key = state_dir.name
                existing = result.get(key)
                is_prod = "prod" in run_dir.name
                existing_is_prod = existing is not None and "prod" in existing.name
                if existing is None:
                    result[key] = run_dir
                elif is_prod and not existing_is_prod:
                    result[key] = run_dir  # prod beats test
                elif is_prod and existing_is_prod:
                    if run_dir.name < existing.name:  # lexicographic: seed0 < seed1
                        result[key] = run_dir
    return result


def _step2_ppi(
    repo_root: Path,
    step_dir: Path,
    ppi_dirs: Optional[Dict[str, Path]] = None,
) -> List[Tuple[str, str]]:
    """Step 2 — symlink each receptor state's PPI run directory."""
    if ppi_dirs is None:
        ppi_dirs = _discover_ppi_dirs(repo_root)

    linked: List[Tuple[str, str]] = []
    for state_name, run_dir in sorted(ppi_dirs.items()):
        if _link_dir(run_dir, step_dir / state_name):
            # Summarise what's inside
            items = []
            if (run_dir / "final_ranking.csv").exists():
                items.append("final_ranking.csv")
            if (run_dir / "final_result").is_dir():
                items.append("final_result/")
            if (run_dir / "cluster_results").is_dir():
                items.append("cluster_results/")
            if (run_dir / "energy_funnel.png").exists():
                items.append("energy_funnel.png")
            summary = ", ".join(items[:4])
            if len(items) > 4:
                summary += f" +{len(items) - 4} more"
            linked.append((f"{state_name}/", summary or "PPI run directory"))

    # Also link LightDock cross-method summaries at top level
    ppi_base = repo_root / "output" / "phase1_ppi"
    if ppi_base.is_dir():
        for pattern, desc in [
            ("*_interface_summary.csv", "LightDock interface summary"),
            ("*_convergence.csv", "PyRosetta-LightDock convergence"),
        ]:
            for csv_file in sorted(ppi_base.glob(pattern)):
                if _safe_symlink(csv_file, step_dir / csv_file.name):
                    linked.append((csv_file.name, desc))

        # Phase 1 review files
        for name, desc in [
            ("phase1_review_report.txt", "Phase 1 comprehensive review"),
            ("phase1_summary.csv", "Phase 1 statistics summary"),
        ]:
            src = ppi_base / name
            if _safe_symlink(src, step_dir / name):
                linked.append((name, desc))

    return linked


def _step3_ppi_evidence(
    project_root: Path,
    repo_root: Path,
    step_dir: Path,
) -> List[Tuple[str, str]]:
    """Step 3 — PPI interface evidence tables."""
    linked: List[Tuple[str, str]] = []
    for name, desc in [
        ("ppi_pyrosetta_residues.csv", "Interface residues per receptor"),
        ("ppi_pyrosetta_summary.csv", "Partner-aware summaries"),
        ("ppi_pyrosetta_residue_long.csv", "Provenance-preserving long-form table"),
        ("ppi_pyrosetta_model_table.csv", "Model-level PPI table"),
        ("ppi_afm_residues.csv", "AlphaFold-Multimer contact residues"),
        ("ppi_afm_summary.csv", "AlphaFold-Multimer summary"),
    ]:
        if _safe_symlink(project_root / name, step_dir / name):
            linked.append((name, desc))

    # Phase 1 interface report
    report = repo_root / "output" / "phase1_ppi" / "phase1_interface_report.md"
    if _safe_symlink(report, step_dir / "phase1_interface_report.md"):
        linked.append(("phase1_interface_report.md", "Phase 1 interface report"))

    # Phase 1 cross-state comparison files
    ppi_base = repo_root / "output" / "phase1_ppi"
    for name, desc in [
        ("ppi_patch_cross_state_comparison.csv", "Cross-state patch comparison"),
        ("ppi_patch_state_robustness.csv", "Patch robustness across states"),
        ("phase1_interface_comparison_report.md", "Cross-state interface comparison report"),
    ]:
        src = ppi_base / name
        if _safe_symlink(src, step_dir / name):
            linked.append((name, desc))

    return linked


def _step4_vina_analysis(
    project_root: Path,
    step_dir: Path,
) -> List[Tuple[str, str]]:
    """Step 4 — Vina pocket clustering and analysis."""
    linked: List[Tuple[str, str]] = []
    for name, desc in [
        ("vina_pose_table.csv", "All poses with pocket assignments"),
        ("vina_pocket_table.csv", "Pocket summaries (centroid, affinity, contacts)"),
        ("vina_drug_pocket_map.csv", "Ligand-pocket relationships"),
        ("vina_pocket_comparison.csv", "Cross-receptor pocket comparison"),
        ("vina_pocket_bootstrap.csv", "Bootstrap stability analysis"),
        ("vina_postprocess_coverage.csv", "Parse coverage metrics"),
        ("vina_contact_distances.csv", "Contact residue distances"),
        ("vina_pocket_residue_occupancy.csv", "Per-residue pocket occupancy"),
        ("vina_clustering_merge_log.csv", "Cluster merge decisions"),
        ("vina_clustering_parameters.json", "Clustering parameters (reproducibility)"),
    ]:
        if _safe_symlink(project_root / name, step_dir / name):
            linked.append((name, desc))
    return linked


def _step5_verdict(
    project_root: Path,
    step_dir: Path,
) -> List[Tuple[str, str]]:
    """Step 5 — site verdict."""
    linked: List[Tuple[str, str]] = []
    for name, desc in [
        ("valid_sites.csv", "Prioritized candidate sites"),
        ("cross_method_agreement.csv", "Vina vs PPI agreement scores"),
        ("vina_consensus_sites.csv", "Vina-only consensus sites"),
    ]:
        if _safe_symlink(project_root / name, step_dir / name):
            linked.append((name, desc))
    return linked


def _step6_report(
    project_root: Path,
    step_dir: Path,
) -> List[Tuple[str, str]]:
    """Step 6 — report."""
    linked: List[Tuple[str, str]] = []
    for name, desc in [
        ("project_report.txt", "Comprehensive project report"),
        ("combined_residue_evidence.csv", "Integrated residue findings"),
    ]:
        if _safe_symlink(project_root / name, step_dir / name):
            linked.append((name, desc))
    return linked


def _step7_validation(
    project_root: Path,
    step_dir: Path,
) -> List[Tuple[str, str]]:
    """Step 7 — validation results."""
    linked: List[Tuple[str, str]] = []

    # Pull from existing step7_validate dir (includes subdirs)
    existing = project_root / "step7_validate"
    if existing.is_dir():
        if _link_dir(existing, step_dir / "step7_validate"):
            for f in sorted(existing.rglob("*")):
                if f.is_file():
                    rel = f.relative_to(existing)
                    linked.append((f"step7_validate/{rel}", ""))
    else:
        # Fallback: link any validation files from project root
        for name, desc in [
            ("validation_report.txt", "Validation checklist"),
            ("validation_summary.csv", "Validation summary"),
        ]:
            if _safe_symlink(project_root / name, step_dir / name):
                linked.append((name, desc))

    # Root metadata files produced by output_steps.py
    for name, desc in [
        ("run_status.json", "Run status metadata"),
        ("run_overview.md", "Run overview (Markdown)"),
        ("run_overview.html", "Run overview (HTML)"),
        ("current_run_manifest.json", "Current run manifest"),
        ("step_index.md", "Step index (from output_steps)"),
        ("report_digest.md", "Report digest"),
    ]:
        if _safe_symlink(project_root / name, step_dir / name):
            linked.append((name, desc))

    return linked


def _step_optional_phase(
    repo_root: Path,
    phase_dir_name: str,
    step_dir: Path,
) -> List[Tuple[str, str]]:
    """Generic organizer for optional Phase 2/3/4 analysis directories."""
    linked: List[Tuple[str, str]] = []
    src = repo_root / "output" / phase_dir_name
    if not src.is_dir():
        return linked
    if _link_dir(src, step_dir / phase_dir_name):
        for f in sorted(src.rglob("*")):
            if f.is_file():
                rel = f.relative_to(src)
                linked.append((f"{phase_dir_name}/{rel}", ""))
    return linked


# ---------------------------------------------------------------------------
# README and index generators
# ---------------------------------------------------------------------------


def _write_readme(
    step_dir: Path,
    step_num: int,
    items: List[Tuple[str, str]],
) -> None:
    """Write a README.txt describing step contents."""
    _, title, desc = STEPS[step_num]
    header = f"Step {step_num}: {title}"
    lines = [
        header,
        "=" * len(header),
        "",
        desc,
        "",
    ]
    if items:
        lines.append("Contents:")
        lines.append("-" * 60)
        for name, note in items:
            if note:
                lines.append(f"  {name:<45} {note}")
            else:
                lines.append(f"  {name}")
    else:
        lines.append("(No output files found for this step yet.)")
    lines.extend([
        "",
        "Note: Files are symbolic links to the original output locations.",
        "      Use 'ls -la' to see where each link points.",
        "",
    ])
    (step_dir / "README.txt").write_text("\n".join(lines), encoding="utf-8")


def _write_index(
    steps_root: Path,
    all_items: Dict[int, List[Tuple[str, str]]],
) -> Path:
    """Write STEP_INDEX.txt master overview."""
    lines = [
        "Pipeline Output Index",
        "=" * 40,
        "",
        "Each numbered directory contains symbolic links to all output files",
        "for that pipeline step.  Open a step directory to browse results.",
        "",
        "Steps 1-7: Production pipeline (Vina + PPI docking -> analysis -> report)",
        "Steps 8-10: Extended analysis phases (present only when executed)",
        "",
        "These step views mirror the structure described in CLAUDE.md",
        "(output/{project}/steps/) and complement the legacy step7_validate/ directory.",
        "",
    ]
    total = 0
    for step_num in sorted(STEPS):
        folder, title, _ = STEPS[step_num]
        items = all_items.get(step_num, [])
        n = len(items)
        total += n
        status = f"{n} item(s)" if n > 0 else "(no output yet)"
        lines.append(f"  {folder}/")
        lines.append(f"    Step {step_num}: {title}  [{status}]")
        lines.append("")

    lines.extend([
        "-" * 40,
        f"Total: {total} items across {len(STEPS)} steps",
        "",
        "Tip: Use 'cp -rL steps/ dest/' to make a self-contained copy",
        "     with all symlinks resolved to real files.",
        "",
    ])

    path = steps_root / "STEP_INDEX.txt"
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------


def organize_outputs(
    config_path: str,
    *,
    ppi_dirs: Optional[Dict[str, Path]] = None,
    repo_root: Optional[Path] = None,
    clean: bool = True,
) -> Path:
    """Create step-based output view under ``{project_root}/steps/``.

    Parameters
    ----------
    config_path : str
        Path to the project YAML config.
    ppi_dirs : dict, optional
        Map of state name -> run directory for PPI results.
        Auto-discovered from ``output/phase1_ppi/`` if not provided.
    repo_root : Path, optional
        Repository root.  Auto-detected from *config_path* if not provided.
    clean : bool
        Remove existing ``steps/`` directory before recreating.

    Returns
    -------
    Path
        Path to the created ``STEP_INDEX.txt``.
    """
    config = load_config(config_path)
    project_root = project_root_from_config(config)

    if repo_root is None:
        cfg = Path(config_path).resolve()
        repo_root = cfg.parents[1] if cfg.parent.name == "config" else cfg.parent

    # Resolve to absolute for reliable relative symlinks
    project_root = project_root.resolve()
    repo_root = repo_root.resolve()

    steps_root = project_root / "steps"

    if clean and steps_root.exists():
        shutil.rmtree(steps_root)
    steps_root.mkdir(parents=True, exist_ok=True)

    print(f"\n  Project root: {project_root}")
    print(f"  Organizing outputs -> {steps_root}\n")

    all_items: Dict[int, List[Tuple[str, str]]] = {}

    organizers = [
        (1, lambda sd: _step1_vina(project_root, config, sd)),
        (2, lambda sd: _step2_ppi(repo_root, sd, ppi_dirs)),
        (3, lambda sd: _step3_ppi_evidence(project_root, repo_root, sd)),
        (4, lambda sd: _step4_vina_analysis(project_root, sd)),
        (5, lambda sd: _step5_verdict(project_root, sd)),
        (6, lambda sd: _step6_report(project_root, sd)),
        (7, lambda sd: _step7_validation(project_root, sd)),
        (8, lambda sd: _step_optional_phase(repo_root, "phase2_pockets", sd)),
        (9, lambda sd: _step_optional_phase(repo_root, "phase3_docking", sd)),
        (10, lambda sd: _step_optional_phase(repo_root, "phase4_perturbation", sd)),
    ]

    for step_num, organizer_fn in organizers:
        folder_name, title, _ = STEPS[step_num]
        step_dir = steps_root / folder_name
        step_dir.mkdir(exist_ok=True)

        items = organizer_fn(step_dir)
        all_items[step_num] = items
        _write_readme(step_dir, step_num, items)

        n = len(items)
        if n > 0:
            print(f"  Step {step_num} ({title}): {n} item(s)")
        else:
            print(f"  Step {step_num} ({title}): (no output yet)")

    index_path = _write_index(steps_root, all_items)

    total = sum(len(v) for v in all_items.values())
    print(f"\n  Done: {total} items organized into {len(STEPS)} steps")
    print(f"  Index: {index_path}")

    return index_path
