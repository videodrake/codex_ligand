#!/usr/bin/env python3
"""Remove production outputs so the full pipeline can be rerun from scratch.

This script is intentionally conservative about what it removes:

- current project output root from ``config/example-project.yaml``
- phase-separated downstream output directories
- PyRosetta PPI docking directories derived from the production INI files
- legacy ``EGFR_dimer_*`` docking directories at repository root

It does **not** remove:

- input data
- config files
- pre-qsub success marker under ``output/pre_qsub_status/``
"""

from __future__ import annotations

import argparse
import configparser
import shutil
import sys
from pathlib import Path
from typing import Iterable, List, Set

import yaml


REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from egfr_pipeline.pyrosetta_docking.metadata import build_output_root_name

PROJECT_CONFIG = REPO_ROOT / "config" / "example-project.yaml"
PPI_CONFIGS = (
    REPO_ROOT / "config" / "ppi_prod_TH1.ini",
    REPO_ROOT / "config" / "ppi_prod_beta_meander.ini",
)
PHASE_OUTPUT_DIRS = (
    REPO_ROOT / "output" / "phase1_ppi",
    REPO_ROOT / "output" / "phase2_pockets",
    REPO_ROOT / "output" / "phase3_docking",
    REPO_ROOT / "output" / "phase4_perturbation",
)


def _project_output_dir() -> Path:
    payload = yaml.safe_load(PROJECT_CONFIG.read_text(encoding="utf-8"))
    output_root_raw = payload.get("output_root", "./output")
    output_root = Path(output_root_raw)
    if not output_root.is_absolute():
        output_root = (REPO_ROOT / output_root).resolve()
    project_name = payload.get("project_name", "").strip()
    return output_root / project_name if project_name else output_root


def _ppi_output_dirs() -> Set[Path]:
    targets: Set[Path] = set()

    for config_path in PPI_CONFIGS:
        if not config_path.exists():
            continue

        parser = configparser.ConfigParser()
        parser.read(config_path, encoding="utf-8")
        input_pdb_rel = parser.get("Path", "input_pdb_name", fallback="")
        if not input_pdb_rel:
            continue

        input_pdb = REPO_ROOT / input_pdb_rel
        stem = input_pdb.stem
        targets.add(REPO_ROOT / stem)
        tagged = build_output_root_name(parser, str(input_pdb), stem)
        targets.add(REPO_ROOT / tagged)

    for path in REPO_ROOT.glob("EGFR_dimer_*"):
        if path.is_dir():
            targets.add(path)

    return targets


def collect_cleanup_targets() -> List[Path]:
    targets: List[Path] = []
    seen: Set[Path] = set()

    def add(path: Path) -> None:
        resolved = path.resolve()
        if resolved not in seen:
            seen.add(resolved)
            targets.append(resolved)

    add(_project_output_dir())
    for path in PHASE_OUTPUT_DIRS:
        add(path)
    for path in sorted(_ppi_output_dirs()):
        add(path)

    return targets


def _remove_path(path: Path, execute: bool) -> str:
    if not path.exists():
        return f"[SKIP] {path} (not found)"

    if not execute:
        return f"[DRY-RUN] would remove {path}"

    if path.is_dir():
        shutil.rmtree(path)
        return f"[REMOVED] {path}"

    path.unlink()
    return f"[REMOVED] {path}"


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Reset production outputs for a clean full rerun."
    )
    parser.add_argument(
        "--execute",
        action="store_true",
        help="Actually remove the targets. Without this flag, only print a dry-run plan.",
    )
    args = parser.parse_args()

    targets = collect_cleanup_targets()

    print("Production output reset plan")
    print(f"Repository root: {REPO_ROOT}")
    print(f"Execute removals: {args.execute}")
    print("")
    for path in targets:
        print(_remove_path(path, execute=args.execute))

    if not args.execute:
        print("")
        print("Dry-run only. Re-run with --execute to remove the listed targets.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
