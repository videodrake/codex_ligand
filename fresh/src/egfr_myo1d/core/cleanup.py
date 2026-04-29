"""M1 Phase 1: cleanup manager.

Implements the cleanup policy from milestone1_foundation_codex_handoff_v0_5.md §18:

- mode=test:       deletes intermediate/scratch files inside ctx.run_dir
- mode=production: defaults to dry-run unless explicitly overridden via --dry-run false
- always preserves: manifest/, logs/, qc/, reports/ (all files within)
- refuses to delete anything outside ctx.run_dir
- always writes manifest/cleanup_report.json
- appends phase status
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from egfr_myo1d.core.logging_utils import append_phase_status
from egfr_myo1d.core.manifest import now_iso, write_json
from egfr_myo1d.core.run_context import RunContext, RunContextError


# Per handoff §18.2 deletion candidates
DELETION_FILE_PATTERNS = (
    "*.tmp",
    "*.pdbqt.tmp",
    "*.vina.tmp",
    "*.silent",
)
DELETION_DIR_BASENAMES = ("scratch", "tmp")

# Per handoff §18.3 must-NEVER-delete preservation list
PRESERVED_DIR_BASENAMES = ("manifest", "logs", "qc", "reports")


@dataclass
class CleanupReport:
    run_id: str
    run_dir: str
    mode: str
    dry_run: bool
    profile: str
    timestamp: str
    deleted_files: list[dict[str, Any]] = field(default_factory=list)
    preserved_files: list[dict[str, Any]] = field(default_factory=list)
    errors: list[dict[str, Any]] = field(default_factory=list)
    deleted_count: int = 0
    preserved_count: int = 0
    status: str = "PASS"

    def to_json(self) -> dict[str, Any]:
        return {
            "run_id": self.run_id,
            "run_dir": self.run_dir,
            "mode": self.mode,
            "dry_run": self.dry_run,
            "profile": self.profile,
            "timestamp": self.timestamp,
            "deleted_files": self.deleted_files,
            "preserved_files": self.preserved_files,
            "deleted_count": self.deleted_count,
            "preserved_count": self.preserved_count,
            "errors": self.errors,
            "status": self.status,
        }


def resolve_dry_run_default(mode: str) -> bool:
    """Return the default dry_run value per mode (handoff §18.4)."""
    if mode == "production":
        return True
    return False


def _is_under(child: Path, parent: Path) -> bool:
    try:
        child.resolve().relative_to(parent.resolve())
        return True
    except ValueError:
        return False


def _safe_rel(path: Path, run_dir: Path) -> str:
    try:
        return path.resolve().relative_to(run_dir.resolve()).as_posix()
    except ValueError:
        return str(path)


def _enumerate_deletion_candidates(run_dir: Path) -> list[Path]:
    """Collect file paths eligible for deletion under run_dir per §18.2.

    Two source sets:
    1. Any file (recursive) inside scratch/ or tmp/ directories.
    2. Files anywhere under run_dir (excluding scratch/, tmp/, and preserved
       dirs) whose name matches DELETION_FILE_PATTERNS.

    Returns deduplicated list in deterministic order.
    """
    seen: set[Path] = set()
    candidates: list[Path] = []

    for dir_name in DELETION_DIR_BASENAMES:
        target = run_dir / dir_name
        if not target.is_dir():
            continue
        for entry in sorted(target.rglob("*")):
            if entry.is_dir():
                continue
            if entry.is_symlink() or entry.is_file():
                resolved = entry
                if resolved not in seen:
                    seen.add(resolved)
                    candidates.append(resolved)

    for entry in sorted(run_dir.rglob("*")):
        if entry.is_dir():
            continue
        if not (entry.is_file() or entry.is_symlink()):
            continue
        try:
            rel = entry.relative_to(run_dir)
        except ValueError:
            continue
        first_part = rel.parts[0] if rel.parts else ""
        if first_part in PRESERVED_DIR_BASENAMES:
            continue
        if first_part in DELETION_DIR_BASENAMES:
            # Already collected above
            continue
        for pattern in DELETION_FILE_PATTERNS:
            if entry.match(pattern):
                if entry not in seen:
                    seen.add(entry)
                    candidates.append(entry)
                break

    return candidates


def _enumerate_preserved_files(run_dir: Path) -> list[Path]:
    """Collect every file under PRESERVED_DIR_BASENAMES for the report."""
    preserved: list[Path] = []
    for dir_name in PRESERVED_DIR_BASENAMES:
        target = run_dir / dir_name
        if not target.is_dir():
            continue
        for entry in sorted(target.rglob("*")):
            if entry.is_file():
                preserved.append(entry)
    return preserved


def run_cleanup(
    ctx: RunContext,
    mode: str,
    dry_run: bool | None = None,
    profile: str = "codex_dev",
) -> CleanupReport:
    """Execute cleanup per handoff §18 policy.

    Parameters
    ----------
    ctx
        RunContext for an existing run directory. Cleanup operates only
        inside ``ctx.run_dir``.
    mode
        One of ``"test"`` or ``"production"``.
    dry_run
        If ``None``, defaults to ``False`` for ``test`` and ``True`` for
        ``production`` (handoff §18.4). Otherwise honors the explicit value.
    profile
        Recorded in the report only; ``codex_dev`` or ``hpc_strict``.

    Returns
    -------
    CleanupReport
        Populated report. Also written to ``manifest/cleanup_report.json``.

    Raises
    ------
    ValueError
        If ``mode`` is not ``"test"`` or ``"production"``.
    RunContextError
        If ``run_dir`` does not exist or a deletion candidate would escape it.
    """
    if mode not in {"test", "production"}:
        raise ValueError(
            "mode must be 'test' or 'production'; got: {0}".format(mode)
        )

    if dry_run is None:
        dry_run = resolve_dry_run_default(mode)

    run_dir = ctx.run_dir
    if not run_dir.is_dir():
        raise RunContextError("run_dir does not exist: {0}".format(run_dir))

    candidates = _enumerate_deletion_candidates(run_dir)
    preserved_files = _enumerate_preserved_files(run_dir)

    deleted_records: list[dict[str, Any]] = []
    errors: list[dict[str, Any]] = []
    deleted_count = 0

    for candidate in candidates:
        # Defensive double-check: refuse outside run_dir even if collector slips
        if not _is_under(candidate, run_dir):
            raise RunContextError(
                "cleanup refuses to delete outside run_dir: {0}".format(candidate)
            )

        try:
            size_bytes = candidate.stat().st_size if candidate.exists() else 0
        except OSError:
            size_bytes = 0

        record: dict[str, Any] = {
            "path": _safe_rel(candidate, run_dir),
            "size_bytes": size_bytes,
            "deleted": False,
        }

        if dry_run:
            deleted_records.append(record)
            continue

        try:
            candidate.unlink()
            record["deleted"] = True
            deleted_count += 1
        except OSError as exc:
            errors.append(
                {
                    "path": _safe_rel(candidate, run_dir),
                    "error": str(exc),
                }
            )
            record["error"] = str(exc)
        deleted_records.append(record)

    preserved_records = [
        {"path": _safe_rel(p, run_dir), "reason": "preservation_list"}
        for p in preserved_files
    ]

    status = "WARN" if errors else "PASS"

    report = CleanupReport(
        run_id=ctx.run_id,
        run_dir=str(run_dir),
        mode=mode,
        dry_run=dry_run,
        profile=profile,
        timestamp=now_iso(),
        deleted_files=deleted_records,
        preserved_files=preserved_records,
        deleted_count=deleted_count if not dry_run else len(deleted_records),
        preserved_count=len(preserved_records),
        errors=errors,
        status=status,
    )

    report_path = ctx.manifest_dir / "cleanup_report.json"
    write_json(report_path, report.to_json(), ctx)

    append_phase_status(
        ctx,
        phase="cleanup",
        status=status,
        message=(
            "cleanup mode={0} dry_run={1} candidates={2} deleted={3} "
            "preserved={4} errors={5}".format(
                mode,
                dry_run,
                len(deleted_records),
                report.deleted_count,
                report.preserved_count,
                len(errors),
            )
        ),
        details={
            "mode": mode,
            "dry_run": dry_run,
            "candidate_count": len(deleted_records),
            "deleted_count": report.deleted_count,
            "preserved_count": report.preserved_count,
            "errors_count": len(errors),
            "profile": profile,
            "cleanup_report": _safe_rel(report_path, run_dir),
        },
    )

    return report
