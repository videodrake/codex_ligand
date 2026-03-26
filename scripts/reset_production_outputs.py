#!/usr/bin/env python3
"""Remove production outputs so the full pipeline can be rerun from scratch.

This script is intentionally conservative about what it removes:

- active workflow outputs (``workflow_a/workflow_b/precheck``)
- optional legacy outputs only when ``--include-legacy`` is explicitly enabled

It does **not** remove:

- input data
- config files
- pre-qsub success marker under ``output/pre_qsub_status/``
"""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path
from typing import List, Set

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

ACTIVE_WORKFLOW_OUTPUT_DIRS = (
    REPO_ROOT / "output" / "workflow_a",
    REPO_ROOT / "output" / "workflow_b",
    REPO_ROOT / "output" / "precheck",
)
LEGACY_OPTIONAL_OUTPUT_DIRS = (
    REPO_ROOT / "output" / "phase1_ppi",
    REPO_ROOT / "output" / "egfr_myo1d_vina",
)
LEGACY_OPTIONAL_GLOBS = (
    "output/step*",
    "step*",
)


def collect_cleanup_targets(*, include_legacy: bool) -> tuple[List[Path], List[Path]]:
    active_targets: List[Path] = []
    legacy_targets: List[Path] = []
    seen: Set[Path] = set()

    def add(path: Path, *, group: List[Path]) -> None:
        resolved = path.resolve()
        if resolved not in seen:
            seen.add(resolved)
            group.append(resolved)

    for path in ACTIVE_WORKFLOW_OUTPUT_DIRS:
        add(path, group=active_targets)

    if include_legacy:
        for path in LEGACY_OPTIONAL_OUTPUT_DIRS:
            add(path, group=legacy_targets)
        for pattern in LEGACY_OPTIONAL_GLOBS:
            for path in sorted(REPO_ROOT.glob(pattern)):
                add(path, group=legacy_targets)

    return active_targets, legacy_targets


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
    parser.add_argument(
        "--include-legacy",
        action="store_true",
        help="Also clean optional legacy outputs (phase1_ppi, egfr_myo1d_vina, and step* artifacts).",
    )
    args = parser.parse_args()

    active_targets, legacy_targets = collect_cleanup_targets(
        include_legacy=args.include_legacy
    )

    print("Production output reset plan")
    print(f"Repository root: {REPO_ROOT}")
    print(f"Execute removals: {args.execute}")
    print(f"Include legacy targets: {args.include_legacy}")
    print("")
    print("=== Active workflow targets (always) ===")
    for path in active_targets:
        print(_remove_path(path, execute=args.execute))

    print("")
    if legacy_targets:
        print("=== Legacy optional targets (--include-legacy) ===")
        for path in legacy_targets:
            print(_remove_path(path, execute=args.execute))
    else:
        print("=== Legacy optional targets (--include-legacy) ===")
        print("[SKIP] legacy cleanup disabled")

    if not args.execute:
        print("")
        print("Dry-run only. Re-run with --execute to remove the listed targets.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())