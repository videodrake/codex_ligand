#!/usr/bin/env python3
"""Phase 2 Cascade Runner — sequential re-execution of TG 2.0 through 2.7.

Usage:
  # Full cascade (setup + parse + all TGs)
  python -m egfr_pipeline.phase2.rerun_cascade

  # Parse-only (fpocket already done on server)
  python -m egfr_pipeline.phase2.rerun_cascade --parse-only

  # Resume from a specific TG
  python -m egfr_pipeline.phase2.rerun_cascade --parse-only --from-tg 2.5

  # With FTMap data
  python -m egfr_pipeline.phase2.rerun_cascade --parse-only --ftmap-dir /path/to/ftmap
"""

import argparse
import sys
import time
from pathlib import Path
from typing import Optional

from egfr_pipeline import paths

_CFG = {"output_root": str(paths.REPO_ROOT / "output")}
PHASE1_INPUT_DIR = paths.REPO_ROOT / "input" / "PPI" / "phase1"


def run_phase2_cascade(
    parse_only: bool = False,
    from_tg: str = "2.0",
    ftmap_dir: Optional[Path] = None,
) -> None:
    """Run Phase 2 TG 2.0 -> 2.7 sequentially.

    Parameters
    ----------
    parse_only : bool
        If True, skip --setup in TG 2.1 (assume fpocket/P2Rank already ran).
    from_tg : str
        Start from this TG (e.g. "2.3" to skip 2.0-2.2).
    ftmap_dir : Path, optional
        FTMap output directory for TG 2.4 druggability enrichment.
    """
    output_dir = paths.wb_phase2_pocket_analysis(_CFG)
    receptor_dir = PHASE1_INPUT_DIR
    phase1_dir = paths.wb_phase1_ppi_analysis(_CFG)

    steps = [
        ("2.0", "Patch Ingestion", lambda: _run_tg20(phase1_dir, output_dir)),
        ("2.1", "Pocket Proposal", lambda: _run_tg21(parse_only, output_dir, receptor_dir)),
        ("2.2", "Pocket Merge", lambda: _run_tg22(output_dir)),
        ("2.3", "Patch Relationship", lambda: _run_tg23(output_dir, receptor_dir)),
        ("2.4", "Druggability Confidence", lambda: _run_tg24(output_dir, ftmap_dir)),
        ("2.5", "Cross-State Alignment", lambda: _run_tg25(output_dir)),
        ("2.6", "Phase 3 Export", lambda: _run_tg26(output_dir)),
        ("2.7", "Review Report", lambda: _run_tg27(output_dir)),
    ]

    print("=" * 60)
    print("  Phase 2 Cascade Runner")
    print(f"  Output:     {output_dir}")
    print(f"  Parse-only: {parse_only}")
    print(f"  From TG:    {from_tg}")
    print(f"  FTMap dir:  {ftmap_dir or '(none)'}")
    print("=" * 60)

    t0 = time.time()
    for tg_id, label, func in steps:
        if tg_id < from_tg:
            print(f"\n--- TG {tg_id} {label} --- SKIPPED (before --from-tg {from_tg})")
            continue
        print(f"\n{'=' * 60}")
        print(f"  TG {tg_id}: {label}")
        print(f"{'=' * 60}")
        ts = time.time()
        try:
            func()
        except Exception as e:
            print(f"\n  [ERROR] TG {tg_id} failed: {e}")
            print("  Cascade stopped.")
            sys.exit(1)
        elapsed = time.time() - ts
        print(f"  TG {tg_id} done ({elapsed:.1f}s)")

    total = time.time() - t0
    print(f"\n{'=' * 60}")
    print(f"  Phase 2 Cascade COMPLETE ({total:.1f}s)")
    print(f"{'=' * 60}")


# ---------------------------------------------------------------------------
# Per-TG wrappers (thin delegation to existing run_* functions)
# ---------------------------------------------------------------------------

def _run_tg20(phase1_dir: Path, output_dir: Path) -> None:
    from egfr_pipeline.phase2.patch_ingestion import run_patch_ingestion
    run_patch_ingestion(phase1_dir, output_dir)


def _run_tg21(parse_only: bool, output_dir: Path, receptor_dir: Path) -> None:
    from egfr_pipeline.phase2.pocket_proposal import run_pocket_proposal
    run_pocket_proposal(
        do_setup=not parse_only,
        do_parse=True,
        output_dir=output_dir,
        receptor_dir=receptor_dir,
    )


def _run_tg22(output_dir: Path) -> None:
    from egfr_pipeline.phase2.pocket_merge import run_pocket_merge
    run_pocket_merge(output_dir)


def _run_tg23(output_dir: Path, receptor_dir: Path) -> None:
    from egfr_pipeline.phase2.patch_relationship import run_patch_relationship
    run_patch_relationship(output_dir, receptor_dir)


def _run_tg24(output_dir: Path, ftmap_dir: Optional[Path]) -> None:
    from egfr_pipeline.phase2.druggability_confidence import run_druggability_confidence
    run_druggability_confidence(output_dir, ftmap_dir)


def _run_tg25(output_dir: Path) -> None:
    from egfr_pipeline.phase2.cross_state_alignment import run_cross_state_alignment
    run_cross_state_alignment(output_dir)


def _run_tg26(output_dir: Path) -> None:
    from egfr_pipeline.phase2.phase3_export import run_phase3_export
    run_phase3_export(output_dir)


def _run_tg27(output_dir: Path) -> None:
    from egfr_pipeline.phase2.review_report import run_review_report
    run_review_report(output_dir)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Phase 2 Cascade Runner (TG 2.0 -> 2.7)",
    )
    parser.add_argument(
        "--parse-only", action="store_true",
        help="Skip --setup in TG 2.1 (fpocket/P2Rank already executed)",
    )
    parser.add_argument(
        "--from-tg", default="2.0",
        help="Start from this TG (e.g. '2.5' to skip earlier steps)",
    )
    parser.add_argument(
        "--ftmap-dir", type=Path, default=None,
        help="FTMap output directory for druggability enrichment",
    )
    args = parser.parse_args()

    run_phase2_cascade(
        parse_only=args.parse_only,
        from_tg=args.from_tg,
        ftmap_dir=args.ftmap_dir,
    )


if __name__ == "__main__":
    main()
