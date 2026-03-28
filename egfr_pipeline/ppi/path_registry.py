"""PPI path resolution helpers with explicit provenance tags.

This module centralizes path resolution policies used by step-view builders.
Each resolver returns a ``ResolvedPath`` object that captures whether the
path came from canonical layout rules, legacy fallback rules, or is missing.
"""

from __future__ import annotations

import configparser
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Union

from egfr_pipeline.pyrosetta_docking.run_metadata import build_output_root_name


@dataclass(frozen=True)
class ResolvedPath:
    """Path resolution result with provenance."""

    path: Optional[Path]
    source_type: str  # canonical | legacy | missing


def resolve_target_run_dir(target: dict, repo_root: Path) -> ResolvedPath:
    """Resolve a Step 2 PPI target run directory.

    Resolution order:
    1) explicit ``docking_dir``
    2) derive from ``config_ini`` + ``input_pdb``
    """
    direct_dir = target.get("docking_dir")
    if direct_dir:
        path = Path(str(direct_dir))
        resolved = path if path.is_absolute() else repo_root / path
        return ResolvedPath(path=resolved, source_type="canonical")

    config_ini = target.get("config_ini")
    input_pdb_hint = target.get("input_pdb")
    if not config_ini or not input_pdb_hint:
        return ResolvedPath(path=None, source_type="missing")

    config_path = repo_root / str(config_ini)
    if not config_path.exists():
        return ResolvedPath(path=None, source_type="missing")

    parser = configparser.ConfigParser()
    parser.read(config_path, encoding="utf-8")

    input_pdb_rel = parser.get("Path", "input_pdb_name", fallback=str(input_pdb_hint))
    input_pdb = repo_root / input_pdb_rel
    if not input_pdb.exists():
        alt_input = repo_root / str(input_pdb_hint)
        if alt_input.exists():
            input_pdb = alt_input
        else:
            return ResolvedPath(path=None, source_type="missing")

    root_name = build_output_root_name(parser, str(input_pdb), input_pdb.stem)
    return ResolvedPath(path=repo_root / root_name, source_type="canonical")


def resolve_ranking_path(run_dir: Path) -> ResolvedPath:
    """Resolve the ranking CSV from a run directory."""
    preferred = run_dir / "final_result" / "final_ranking.csv"
    if preferred.exists():
        return ResolvedPath(path=preferred, source_type="canonical")

    fallback = run_dir / "final_ranking.csv"
    if fallback.exists():
        return ResolvedPath(path=fallback, source_type="legacy")

    return ResolvedPath(path=preferred, source_type="missing")


def resolve_metadata_path(
    run_dir: Path,
    explicit_path: Optional[Union[Path, str]] = None,
) -> ResolvedPath:
    """Resolve run metadata JSON with canonical/legacy source type."""
    if explicit_path:
        path = Path(explicit_path)
        if path.exists():
            return ResolvedPath(path=path, source_type="canonical")
        if not path.is_absolute():
            candidate = run_dir / path
            if candidate.exists():
                return ResolvedPath(path=candidate, source_type="canonical")

    for candidate, source in (
        (run_dir / "pyrosetta_run_metadata.json", "canonical"),
        (run_dir / "final_result" / "pyrosetta_run_metadata.json", "canonical"),
        (run_dir / "restored" / "pyrosetta_run_metadata.json", "legacy"),
    ):
        if candidate.exists():
            return ResolvedPath(path=candidate, source_type=source)

    return ResolvedPath(path=None, source_type="missing")


def resolve_phase1_interface_report(
    repo_root: Path,
    *,
    allow_legacy_paths: bool = False,
) -> ResolvedPath:
    """Resolve phase1 interface report with optional legacy fallback."""
    new_path = repo_root / "output" / "workflow_b" / "phase1_ppi_analysis" / "phase1_interface_report.md"
    if new_path.exists():
        return ResolvedPath(path=new_path, source_type="canonical")

    legacy = repo_root / "output" / "phase1_ppi" / "phase1_interface_report.md"
    if allow_legacy_paths and legacy.exists():
        return ResolvedPath(path=legacy, source_type="legacy")

    return ResolvedPath(path=None, source_type="missing")
