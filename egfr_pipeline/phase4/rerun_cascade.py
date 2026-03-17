#!/usr/bin/env python3
"""Phase 4 Cascade Runner — sequential execution of TG 4.0 through 4.6.

Usage:
  # Full cascade
  python -m egfr_pipeline.phase4.rerun_cascade

  # Resume from a specific TG
  python -m egfr_pipeline.phase4.rerun_cascade --from-tg 4.3
"""

import argparse
import sys
import time
from pathlib import Path

from egfr_pipeline import paths

_CFG = {"output_root": str(paths.REPO_ROOT / "output")}


def run_phase4_cascade(
    from_tg: str = "4.0",
) -> None:
    """Run Phase 4 TG 4.0 -> 4.6 sequentially.

    Parameters
    ----------
    from_tg : str
        Start from this TG (e.g. "4.3" to skip earlier steps).
    """
    output_dir = paths.wb_phase4_scoring(_CFG)

    # Handoff file validation
    handoff = paths.wb_phase3_focused_docking(_CFG) / "phase4_docking_evidence_reference.csv"
    if not handoff.exists():
        raise FileNotFoundError(
            f"Phase 3 handoff file not found: {handoff}\n"
            "Phase 3 (Focused Docking) must complete before Phase 4.\n"
            "Run Phase 3 post-docking analysis first."
        )

    steps = [
        ("4.0", "Evidence Ingestion", lambda: _run_tg40(output_dir)),
        ("4.1", "Score Framework", lambda: _run_tg41(output_dir)),
        ("4.2", "Mechanistic Classification", lambda: _run_tg42(output_dir)),
        ("4.3", "Perturbation Scoring", lambda: _run_tg43(output_dir)),
        ("4.4", "State Interpretation", lambda: _run_tg44(output_dir)),
        ("4.5", "Review Output", lambda: _run_tg45(output_dir)),
        ("4.6", "Final Report + Presentation", lambda: _run_tg46(output_dir)),
    ]

    print("=" * 60)
    print("  Phase 4 Cascade Runner")
    print(f"  Output:   {output_dir}")
    print(f"  From TG:  {from_tg}")
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
    print(f"  Phase 4 Cascade COMPLETE ({total:.1f}s)")
    print(f"{'=' * 60}")


# ---------------------------------------------------------------------------
# Per-TG wrappers
# ---------------------------------------------------------------------------

def _run_tg40(output_dir: Path) -> None:
    from egfr_pipeline.phase4.evidence_ingestion import run_evidence_ingestion
    run_evidence_ingestion(output_dir=output_dir)


def _run_tg41(output_dir: Path) -> None:
    from egfr_pipeline.phase4.score_framework import run_score_framework
    run_score_framework(output_dir=output_dir)


def _run_tg42(output_dir: Path) -> None:
    from egfr_pipeline.phase4.mechanistic_classification import run_mechanistic_classification
    run_mechanistic_classification(output_dir=output_dir)


def _run_tg43(output_dir: Path) -> None:
    from egfr_pipeline.phase4.perturbation_scoring import run_perturbation_scoring
    run_perturbation_scoring(output_dir=output_dir)


def _run_tg44(output_dir: Path) -> None:
    from egfr_pipeline.phase4.state_interpretation import run_state_interpretation
    run_state_interpretation(output_dir=output_dir)


def _run_tg45(output_dir: Path) -> None:
    from egfr_pipeline.phase4.review_report import run_review_output
    run_review_output(output_dir=output_dir)


def _run_tg46(output_dir: Path) -> None:
    from egfr_pipeline.phase4.final_report import run_final_report
    run_final_report(output_dir=output_dir)

    from egfr_pipeline.phase4.presentation_summary import run_presentation_summary
    run_presentation_summary(output_dir=output_dir)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Phase 4 Cascade Runner (TG 4.0 -> 4.6)",
    )
    parser.add_argument(
        "--from-tg", default="4.0",
        help="Start from this TG (e.g. '4.3' to skip earlier steps)",
    )
    args = parser.parse_args()

    run_phase4_cascade(from_tg=args.from_tg)


if __name__ == "__main__":
    main()
