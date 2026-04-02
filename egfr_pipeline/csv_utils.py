# -*- coding: utf-8 -*-
"""Shared CSV I/O utilities used across all pipeline phases.

Consolidates the repeated _load_csv/_write_csv patterns found in 50+ files.
"""
import csv
from pathlib import Path
from typing import Dict, List, Optional, Sequence


def load_csv(path: Path, encoding: str = "utf-8") -> List[Dict[str, str]]:
    """Load a CSV file as a list of dicts. Returns [] if file does not exist."""
    if not path.exists():
        return []
    with open(path, encoding=encoding) as f:
        return list(csv.DictReader(f))


def write_csv(
    path: Path,
    rows: Sequence[Dict[str, object]],
    fieldnames: Sequence[str],
    encoding: str = "utf-8",
) -> None:
    """Write *rows* as a CSV file with the given *fieldnames*.

    Creates parent directories if they don't exist.
    Extra keys in rows are silently ignored.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding=encoding) as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
