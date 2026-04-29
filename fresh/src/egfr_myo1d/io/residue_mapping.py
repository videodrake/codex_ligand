"""Residue mapping CSV writer/reader for M1 receptor normalization.

Schema per ``milestone1_foundation_codex_handoff_v0_5.md`` §14.3.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any, Iterable


MAPPING_CSV_COLUMNS = [
    "state",
    "source_file",
    "protomer_id",
    "source_chain",
    "source_resseq",
    "source_icode",
    "source_resname",
    "runtime_chain",
    "runtime_resseq",
    "atom_count",
    "role",
]


def write_residue_mapping(path: Path, rows: Iterable[dict[str, Any]], ctx=None) -> None:
    """Write a residue mapping CSV.

    If ``ctx`` is provided, the path is checked against ``ctx.run_dir`` to
    refuse writes outside the run directory.
    """
    target = Path(path)
    if ctx is not None:
        target = ctx.require_within_run_dir(target)
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=MAPPING_CSV_COLUMNS, extrasaction="ignore"
        )
        writer.writeheader()
        for row in rows:
            normalized = {col: row.get(col, "") for col in MAPPING_CSV_COLUMNS}
            # Empty insertion code → "_" placeholder for CSV readability
            if not normalized["source_icode"]:
                normalized["source_icode"] = "_"
            writer.writerow(normalized)


def read_residue_mapping(path: Path) -> list[dict[str, Any]]:
    """Read a residue mapping CSV. Restores empty insertion code from ``_``."""
    rows: list[dict[str, Any]] = []
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            if row.get("source_icode") == "_":
                row["source_icode"] = ""
            rows.append(row)
    return rows
