"""Confidentiality helpers for M3 compound-stage public outputs."""

from __future__ import annotations

import csv
import fnmatch
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from egfr_myo1d.core.run_context import RunContext, ensure_within


PUBLIC_COMPOUND_IDS = ["Cpd-A", "Cpd-B", "Cpd-C"]
PRIVATE_PATH_REDACTED = "PRIVATE_PATH_REDACTED"

RECOMMENDED_GITIGNORE_PATTERNS = [
    "fresh/data/raw/ligands/*",
    "fresh/data/private/*",
    "fresh/data/private/**",
    "fresh/runs/*/phase3_compounds/ligands/*",
    "fresh/runs/*/phase3_compounds/ligands/**",
    "fresh/runs/*/phase3_compounds/prepared_ligands/*",
    "fresh/runs/*/phase3_compounds/prepared_ligands/**",
    "*.pdbqt.tmp",
    "*.vina.tmp",
]


@dataclass(frozen=True)
class PrivateMapEntry:
    public_id: str
    internal_id: str
    notes_present: bool


def load_private_map(path: Path) -> tuple[dict[str, PrivateMapEntry], list[str]]:
    """Read private ID mapping without exposing internal IDs to callers' outputs."""
    warnings: list[str] = []
    if not path.is_file():
        warnings.append("private mapping not found")
        return {}, warnings

    entries: dict[str, PrivateMapEntry] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"public_id", "internal_id", "notes"}
        if not required.issubset(set(reader.fieldnames or [])):
            warnings.append("private mapping schema is missing required columns")
            return {}, warnings
        for row in reader:
            public_id = (row.get("public_id") or "").strip()
            internal_id = (row.get("internal_id") or "").strip()
            if public_id in PUBLIC_COMPOUND_IDS and internal_id:
                entries[public_id] = PrivateMapEntry(
                    public_id=public_id,
                    internal_id=internal_id,
                    notes_present=bool((row.get("notes") or "").strip()),
                )
    return entries, warnings


def update_gitignore(repo_root: Path) -> tuple[bool, list[str]]:
    """Append missing sensitive-data ignore rules without removing existing lines."""
    gitignore = repo_root / ".gitignore"
    existing = []
    if gitignore.is_file():
        existing = gitignore.read_text(encoding="utf-8").splitlines()
    missing = [pattern for pattern in RECOMMENDED_GITIGNORE_PATTERNS if pattern not in existing]
    if not missing:
        return False, RECOMMENDED_GITIGNORE_PATTERNS
    gitignore.parent.mkdir(parents=True, exist_ok=True)
    with gitignore.open("a", encoding="utf-8") as handle:
        if existing and existing[-1].strip():
            handle.write("\n")
        handle.write("\n# M3 confidential ligand and docking scratch protection\n")
        for pattern in missing:
            handle.write(pattern + "\n")
    return True, RECOMMENDED_GITIGNORE_PATTERNS


def git_ls_files(repo_root: Path, path: Path) -> list[str]:
    try:
        rel = Path(path).resolve().relative_to(Path(repo_root).resolve()).as_posix()
    except ValueError:
        return []
    try:
        proc = subprocess.run(
            ["git", "ls-files", "--", rel],
            cwd=str(repo_root),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )
    except OSError:
        return []
    if proc.returncode != 0:
        return []
    return [line for line in proc.stdout.splitlines() if line.strip()]


def git_check_ignored(repo_root: Path, path: Path) -> bool:
    try:
        rel = Path(path).resolve().relative_to(Path(repo_root).resolve()).as_posix()
    except ValueError:
        return False
    try:
        proc = subprocess.run(
            ["git", "check-ignore", "-q", "--", rel],
            cwd=str(repo_root),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )
    except OSError:
        proc = None
    if proc is not None and proc.returncode == 0:
        return True
    gitignore = repo_root / ".gitignore"
    if not gitignore.is_file():
        return False
    for line in gitignore.read_text(encoding="utf-8").splitlines():
        pattern = line.strip()
        if not pattern or pattern.startswith("#") or pattern.startswith("!"):
            continue
        if fnmatch.fnmatch(rel, pattern) or fnmatch.fnmatch(rel, pattern.rstrip("/") + "/*"):
            return True
    return False


def _token_pattern(token: str) -> re.Pattern[str]:
    escaped = re.escape(token)
    return re.compile(r"(?<![A-Za-z0-9_.-])" + escaped + r"(?![A-Za-z0-9_.-])")


def token_present(text: str, token: str) -> bool:
    if not token:
        return False
    if token.isdigit():
        return bool(_token_pattern(token).search(text))
    return token in text


def public_output_scan_paths(ctx: RunContext) -> list[Path]:
    roots = [
        ctx.run_dir / "phase3_compounds" / "manifests",
        ctx.run_dir / "phase3_compounds" / "qc",
        ctx.run_dir / "phase3_compounds" / "reports",
        ctx.logs_dir / "jobs",
    ]
    files: list[Path] = []
    for root in roots:
        if root.is_dir():
            for path in sorted(root.rglob("*")):
                if path.is_file():
                    files.append(ctx.require_within_run_dir(path))
    for path in [
        ctx.logs_dir / "phase3_compounds.log",
        ctx.logs_dir / "master.log",
        ctx.logs_dir / "phase_status.jsonl",
        ctx.logs_dir / "job_status.jsonl",
    ]:
        if path.is_file():
            files.append(ctx.require_within_run_dir(path))
    return files


def scan_internal_id_leaks(
    ctx: RunContext,
    entries: Iterable[PrivateMapEntry],
) -> list[dict[str, str]]:
    leaks: list[dict[str, str]] = []
    scan_paths = public_output_scan_paths(ctx)
    for entry in entries:
        for path in scan_paths:
            try:
                text = path.read_text(encoding="utf-8", errors="ignore")
            except OSError:
                continue
            if token_present(text, entry.internal_id):
                leaks.append(
                    {
                        "check_id": "INTERNAL_ID_LEAK",
                        "category": "confidentiality",
                        "scope": "public_outputs_and_logs",
                        "status": "FAIL",
                        "severity": "BLOCKER",
                        "public_id": entry.public_id,
                        "observed_file": ctx.relative_to_repo(path),
                        "details": "leaked_token_label=INTERNAL_ID_FOR_{0}".format(entry.public_id),
                        "recommended_fix": "Remove the private token from public outputs/logs and regenerate M3-T1.",
                    }
                )
    return leaks


def assert_fresh_input_path(ctx: RunContext, path: Path) -> Path:
    return ensure_within(path, ctx.fresh_root)
