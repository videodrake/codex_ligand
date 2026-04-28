"""Centralized logging utilities for fresh runs."""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from egfr_myo1d.core.run_context import RunContext


FAILED_JOBS_HEADER = "timestamp,job_name,status,message\n"


def now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def initialize_logs(ctx: RunContext) -> None:
    ctx.create_directories()
    for path in [
        ctx.logs_dir / "master.log",
        ctx.logs_dir / "phase_status.jsonl",
        ctx.logs_dir / "job_status.jsonl",
        ctx.errors_dir / "error_summary.txt",
    ]:
        ctx.require_within_run_dir(path)
        path.touch(exist_ok=True)

    failed_jobs = ctx.errors_dir / "failed_jobs.csv"
    ctx.require_within_run_dir(failed_jobs)
    if not failed_jobs.exists() or failed_jobs.stat().st_size == 0:
        failed_jobs.write_text(FAILED_JOBS_HEADER, encoding="utf-8")

    log_master(ctx, "INFO", "centralized logs initialized")


def log_master(ctx: RunContext, level: str, message: str) -> None:
    path = ctx.require_within_run_dir(ctx.logs_dir / "master.log")
    line = f"{now_iso()} [{level.upper()}] {message}\n"
    with path.open("a", encoding="utf-8") as handle:
        handle.write(line)


def append_jsonl(path: Path, record: dict[str, Any], ctx: RunContext) -> None:
    safe_path = ctx.require_within_run_dir(path)
    with safe_path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(record, sort_keys=True) + "\n")


def append_phase_status(
    ctx: RunContext,
    phase: str,
    status: str,
    message: str,
    details: dict[str, Any] | None = None,
) -> dict[str, Any]:
    record = {
        "timestamp": now_iso(),
        "phase": phase,
        "status": status,
        "run_id": ctx.run_id,
        "message": message,
        "details": details or {},
    }
    append_jsonl(ctx.logs_dir / "phase_status.jsonl", record, ctx)
    log_master(ctx, status, f"{phase}: {message}")
    if status in {"WARN", "FAIL"}:
        append_error_summary(ctx, f"{phase} {status}: {message}")
    return record


def append_job_status(
    ctx: RunContext,
    job_name: str,
    status: str,
    node: str | None = None,
    ppn: int | None = None,
    stdout: str | None = None,
    stderr: str | None = None,
    details: dict[str, Any] | None = None,
) -> dict[str, Any]:
    record = {
        "timestamp": now_iso(),
        "job_name": job_name,
        "status": status,
        "run_id": ctx.run_id,
        "node": node,
        "ppn": ppn,
        "stdout": stdout,
        "stderr": stderr,
        "details": details or {},
    }
    append_jsonl(ctx.logs_dir / "job_status.jsonl", record, ctx)
    if status == "FAIL":
        append_failed_job(ctx, job_name, status, details.get("message", "") if details else "")
    return record


def append_error_summary(ctx: RunContext, message: str) -> None:
    path = ctx.require_within_run_dir(ctx.errors_dir / "error_summary.txt")
    with path.open("a", encoding="utf-8") as handle:
        handle.write(f"{now_iso()} {message}\n")


def append_failed_job(ctx: RunContext, job_name: str, status: str, message: str) -> None:
    path = ctx.require_within_run_dir(ctx.errors_dir / "failed_jobs.csv")
    with path.open("a", encoding="utf-8") as handle:
        safe_message = str(message).replace("\n", " ").replace(",", ";")
        handle.write(f"{now_iso()},{job_name},{status},{safe_message}\n")
