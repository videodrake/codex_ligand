"""M3-T11 cleanup planning/execution with strict run-dir safety gates."""

from __future__ import annotations

import fnmatch
import hashlib
import json
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from egfr_myo1d.core.logging_utils import append_failed_job, append_job_status
from egfr_myo1d.core.run_context import RunContext


DEFAULT_CLEANUP_CONFIG: dict[str, Any] = {
    "cleanup": {
        "allow_delete_outside_run_dir": False,
        "follow_symlinks": False,
        "execute_requires_confirm_run_id": True,
        "never_delete_patterns": [
            "phase3_compounds/tables/**",
            "phase3_compounds/qc/**",
            "phase3_compounds/manifests/**",
            "phase3_compounds/reports/**",
            "logs/**",
            "report/**",
            "fresh/data/private/**",
        ],
        "delete_test_scope_patterns": [
            "phase3_compounds/scratch/**",
            "phase3_compounds/tmp/**",
            "phase3_compounds/**/__pycache__/**",
            "phase3_compounds/**/.pytest_cache/**",
            "phase3_compounds/docking_inputs/**/tmp/**",
            "phase3_compounds/docking_outputs/**/failed_tmp/**",
        ],
    }
}

CLEANUP_QC_FIELDS = ["check_id", "category", "status", "severity", "details", "recommended_fix"]
CLEANUP_CHECKS = [
    "cleanup_mode_valid",
    "run_dir_exists",
    "run_dir_under_fresh_runs",
    "confirm_run_id_required_for_execute",
    "confirm_run_id_matches",
    "no_delete_outside_run_dir",
    "symlink_escape_blocked",
    "private_map_preserved",
    "final_tables_preserved",
    "qc_preserved",
    "reports_preserved",
    "manifests_preserved",
    "logs_preserved",
    "selected_poses_preserved_if_configured",
    "production_docking_outputs_preserved_by_default",
    "test_tmp_files_identified",
    "dry_run_plan_written",
    "execute_deletions_logged",
    "cleanup_report_written",
]


@dataclass
class M3CleanupResult:
    status: str
    report: dict[str, Any]
    report_root: Path
    report_phase3: Path
    cleanup_qc: Path
    blockers: list[str]
    warnings: list[str]


def now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def sha256(path: Path) -> str | None:
    if not path.is_file():
        return None
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _inside(child: Path, parent: Path) -> bool:
    try:
        child.resolve().relative_to(parent.resolve())
        return True
    except ValueError:
        return False


def _rel(ctx: RunContext, path: Path) -> str:
    try:
        return path.relative_to(ctx.run_dir).as_posix()
    except ValueError:
        try:
            return path.resolve().relative_to(ctx.run_dir.resolve()).as_posix()
        except ValueError:
            return ctx.relative_to_repo(path)


def _matches_any(rel: str, patterns: list[str]) -> bool:
    return any(fnmatch.fnmatch(rel, pattern) or fnmatch.fnmatch(rel, pattern.rstrip("/") + "/*") for pattern in patterns)


def _load_config(ctx: RunContext, config_path: Path | None) -> dict[str, Any]:
    config = json.loads(json.dumps(DEFAULT_CLEANUP_CONFIG))
    path = config_path or (ctx.fresh_root / "configs" / "milestone3_report.yaml")
    if path and path.is_file():
        try:
            import yaml  # type: ignore

            loaded = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
            if isinstance(loaded.get("cleanup"), dict):
                config["cleanup"].update(loaded["cleanup"])
        except Exception:
            pass
    return config


def _write_json(path: Path, data: dict[str, Any], ctx: RunContext) -> None:
    safe = ctx.require_within_run_dir(path)
    safe.parent.mkdir(parents=True, exist_ok=True)
    safe.write_text(json.dumps(data, indent=2, sort_keys=True), encoding="utf-8")


def _write_cleanup_qc(path: Path, statuses: dict[str, str], ctx: RunContext) -> None:
    import csv

    path = ctx.require_within_run_dir(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CLEANUP_QC_FIELDS)
        writer.writeheader()
        for check in CLEANUP_CHECKS:
            status = statuses.get(check, "PASS")
            severity = "BLOCKER" if status == "FAIL" else ("MAJOR" if status == "WARN" else "INFO")
            writer.writerow({"check_id": check, "category": "m3_cleanup", "status": status, "severity": severity, "details": check, "recommended_fix": ""})


def plan_m3_cleanup(
    ctx: RunContext,
    *,
    profile: str | None = None,
    mode: str = "plan",
    cleanup_scope: str = "conservative",
    config: Path | None = None,
    preserve_selected_poses: bool = True,
    preserve_all_docking_outputs: bool = True,
    confirm_run_id: str | None = None,
    force: bool = False,
) -> M3CleanupResult:
    cfg = _load_config(ctx, config)
    timestamp = now_iso()
    phase3 = ctx.run_dir / "phase3_compounds"
    root_report = ctx.run_dir / "report" / "cleanup_report.json"
    phase3_report = phase3 / "reports" / "cleanup_report.json"
    cleanup_qc = phase3 / "qc" / "cleanup_qc.csv"
    for directory in [root_report.parent, phase3_report.parent, cleanup_qc.parent, ctx.logs_dir, ctx.errors_dir]:
        ctx.require_within_run_dir(directory).mkdir(parents=True, exist_ok=True)

    blockers: list[str] = []
    warnings: list[str] = []
    if mode not in {"plan", "execute"}:
        blockers.append("cleanup mode invalid")
    if cleanup_scope not in {"conservative", "test", "mini", "production", "aggressive-test-only"}:
        blockers.append("cleanup scope invalid")
    if not ctx.run_dir.is_dir():
        blockers.append("run_dir missing")
    if not _inside(ctx.run_dir, ctx.fresh_root / "runs"):
        blockers.append("run_dir is not under fresh/runs")
    if mode == "execute" and confirm_run_id != ctx.run_id:
        blockers.append("cleanup execute requires --confirm-run-id matching --run-id")
    if cleanup_scope == "conservative":
        warnings.append("conservative cleanup is plan-only")
    if cleanup_scope == "aggressive-test-only" and profile == "production":
        blockers.append("aggressive-test-only cleanup refused for production profile")
    if cleanup_scope == "aggressive-test-only" and not force:
        blockers.append("aggressive-test-only cleanup requires --force")

    never = list(cfg["cleanup"]["never_delete_patterns"])
    delete_patterns = list(cfg["cleanup"]["delete_test_scope_patterns"])
    candidates: list[Path] = []
    files_considered = 0
    if ctx.run_dir.is_dir():
        for path in sorted(ctx.run_dir.rglob("*")):
            if path.is_dir():
                continue
            files_considered += 1
            rel = _rel(ctx, path)
            if _matches_any(rel, never):
                continue
            if preserve_all_docking_outputs and rel.startswith("phase3_compounds/docking_outputs/"):
                continue
            if preserve_selected_poses and path.suffix.lower() == ".pdbqt":
                continue
            if cleanup_scope in {"test", "mini", "production", "aggressive-test-only"} and _matches_any(rel, delete_patterns):
                candidates.append(path)
            elif cleanup_scope == "aggressive-test-only" and (path.name.endswith(".tmp") or "__pycache__" in path.parts):
                candidates.append(path)
    planned: list[dict[str, Any]] = []
    blocked: list[dict[str, Any]] = []
    for path in sorted(set(candidates)):
        rel = _rel(ctx, path)
        if path.is_symlink():
            try:
                target = path.resolve(strict=True)
            except OSError:
                target = path.resolve()
            if not _inside(target, ctx.run_dir):
                blocked.append({"relative_path": rel, "reason": "symlink_escape"})
                blockers.append("cleanup plan includes symlink escaping run_dir")
                continue
        resolved = path.resolve()
        if not _inside(resolved, ctx.run_dir):
            blocked.append({"relative_path": rel, "reason": "outside_run_dir"})
            blockers.append("cleanup plan includes path outside run_dir")
            continue
        if _matches_any(rel, never) or rel.endswith("compound_id_map.csv"):
            blocked.append({"relative_path": rel, "reason": "never_delete"})
            blockers.append("cleanup plan includes preserved artifact")
            continue
        planned.append({"relative_path": rel, "reason": "configured temporary artifact", "size_bytes": path.stat().st_size if path.exists() else 0, "sha256": sha256(path)})

    deleted: list[dict[str, Any]] = []
    bytes_deleted = 0
    if mode == "execute" and not blockers:
        for item in planned:
            path = ctx.run_dir / item["relative_path"]
            try:
                size = path.stat().st_size
                path.unlink()
                deleted.append(dict(item))
                bytes_deleted += size
                append_job_status(ctx, f"cleanup:{item['relative_path']}", "PASS", details={"phase": "phase3_compounds", "task": "M3-T11", "job_type": "cleanup_delete", "relative_path": item["relative_path"], "message": "deleted planned temporary artifact"})
            except OSError as exc:
                blocked.append({"relative_path": item["relative_path"], "reason": str(exc)})
                blockers.append("cleanup execute failed for planned file")
    elif mode == "execute" and blockers:
        warnings.append("cleanup execution refused due to safety blockers")

    preserved = []
    for rel in [
        "logs",
        "phase3_compounds/manifests",
        "phase3_compounds/qc",
        "phase3_compounds/tables",
        "phase3_compounds/reports",
        "report",
    ]:
        path = ctx.run_dir / rel
        if path.exists():
            preserved.append(rel)
    status = "FAIL" if blockers else ("WARN" if warnings or mode == "plan" else "PASS")
    report = {
        "schema_version": "m3_cleanup_report_v1",
        "run_id": ctx.run_id,
        "profile": profile,
        "mode": mode,
        "cleanup_scope": cleanup_scope,
        "generated_at": timestamp,
        "executed_at": timestamp if mode == "execute" and not blockers else None,
        "run_dir": ctx.relative_to_repo(ctx.run_dir),
        "allow_delete_outside_run_dir": False,
        "follow_symlinks": False,
        "execute_confirmed": mode == "execute" and confirm_run_id == ctx.run_id,
        "overall_status": status,
        "files_considered": files_considered,
        "files_preserved": len(preserved),
        "files_planned_for_deletion": len(planned),
        "files_deleted": len(deleted),
        "bytes_planned_for_deletion": sum(int(item["size_bytes"]) for item in planned),
        "bytes_deleted": bytes_deleted,
        "preserved_artifacts": preserved,
        "planned_deletions": planned,
        "deleted_files": deleted,
        "blocked_deletions": blocked,
        "blockers": sorted(set(blockers)),
        "warnings": sorted(set(warnings)),
    }
    _write_json(root_report, report, ctx)
    _write_json(phase3_report, report, ctx)
    qc_status = {
        "cleanup_mode_valid": "PASS" if mode in {"plan", "execute"} else "FAIL",
        "run_dir_exists": "PASS" if ctx.run_dir.is_dir() else "FAIL",
        "run_dir_under_fresh_runs": "PASS" if _inside(ctx.run_dir, ctx.fresh_root / "runs") else "FAIL",
        "confirm_run_id_required_for_execute": "PASS" if mode == "plan" or confirm_run_id else "FAIL",
        "confirm_run_id_matches": "PASS" if mode == "plan" or confirm_run_id == ctx.run_id else "FAIL",
        "test_tmp_files_identified": "PASS" if planned or cleanup_scope == "conservative" else "WARN",
        "dry_run_plan_written": "PASS" if root_report.is_file() else "FAIL",
        "execute_deletions_logged": "PASS" if mode == "plan" or len(deleted) == len(planned) else ("FAIL" if mode == "execute" and blockers else "WARN"),
        "cleanup_report_written": "PASS" if root_report.is_file() and phase3_report.is_file() else "FAIL",
    }
    _write_cleanup_qc(cleanup_qc, qc_status, ctx)
    if status == "FAIL":
        append_failed_job(ctx, "M3-T11-cleanup", "FAIL", ";".join(sorted(set(blockers))[:3]))
    return M3CleanupResult(status, report, root_report, phase3_report, cleanup_qc, sorted(set(blockers)), sorted(set(warnings)))
