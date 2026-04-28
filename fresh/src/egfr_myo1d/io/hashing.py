"""Hashing and file description helpers."""

from __future__ import annotations

import hashlib
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


def sha256_file(path: Path) -> str:
    hasher = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            hasher.update(chunk)
    return hasher.hexdigest()


def describe_file(path: Path) -> dict[str, Any]:
    candidate = Path(path)
    if not candidate.exists():
        return {
            "present": False,
            "sha256": None,
            "size_bytes": None,
            "mtime": None,
            "notes": "missing",
        }

    stat = candidate.stat()
    return {
        "present": True,
        "sha256": sha256_file(candidate),
        "size_bytes": stat.st_size,
        "mtime": datetime.fromtimestamp(stat.st_mtime, timezone.utc)
        .replace(microsecond=0)
        .isoformat(),
        "notes": "",
    }
