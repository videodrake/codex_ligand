# -*- coding: utf-8 -*-
"""Generic cascade runner for sequential TG execution.

Extracts the shared orchestration pattern from phase2/3/4 rerun_cascade.py:
- TG version comparison (handles numeric and suffixed IDs like "3.2s")
- Timed execution with error handling
- Skip logic for --from-tg resume
"""
import re
import sys
import time
from typing import Callable, Dict, List, Optional, Sequence, Tuple

Step = Tuple[str, str, Callable[[], None]]  # (tg_id, label, func)


def parse_tg(tg: str) -> Tuple[int, int, str]:
    """Parse TG ID like '3.2' or '3.2s' into a comparable tuple."""
    m = re.match(r"(\d+)\.(\d+)(.*)", tg)
    if m:
        return (int(m.group(1)), int(m.group(2)), m.group(3))
    return (0, 0, tg)


def run_cascade(
    phase_label: str,
    steps: Sequence[Step],
    from_tg: str,
    header_lines: Optional[Sequence[str]] = None,
) -> None:
    """Execute a sequence of TG steps with timing and error handling.

    Parameters
    ----------
    phase_label : str
        Display name (e.g. "Phase 2 Cascade Runner").
    steps : sequence of (tg_id, label, callable)
        Steps to execute in order.
    from_tg : str
        Skip steps whose TG ID is earlier than this.
    header_lines : sequence of str, optional
        Extra lines to print in the banner (e.g. output dir, mode).
    """
    print("=" * 60)
    print(f"  {phase_label}")
    for line in (header_lines or []):
        print(f"  {line}")
    print(f"  From TG:  {from_tg}")
    print("=" * 60)

    from_key = parse_tg(from_tg)
    t0 = time.time()

    for tg_id, label, func in steps:
        if parse_tg(tg_id) < from_key:
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
    print(f"  {phase_label} COMPLETE ({total:.1f}s)")
    print(f"{'=' * 60}")
