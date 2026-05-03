"""PBS chunk runner contract for M3-T5 planned Vina jobs.

The runner executes only rows assigned to a requested chunk from a previously
written vina_job_manifest.csv. It does not parse affinities, poses, clusters,
scores, rankings, or candidate claims.
"""

from __future__ import annotations

import csv
import hashlib
import json
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from egfr_myo1d.core.logging_utils import append_job_status, now_iso
from egfr_myo1d.core.run_context import RunContext


@dataclass
class VinaChunkRunResult:
    status: str
    chunk_id: str
    rows_seen: int
    rows_run: int
    rows_failed: int
    dry_run_runner: bool


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _resolve_repo_path(ctx: RunContext, value: str) -> Path:
    path = Path(value)
    if not path.is_absolute():
        path = ctx.repo_root / path
    return path.resolve()


def _inside(path: Path, parent: Path) -> bool:
    try:
        path.resolve().relative_to(parent.resolve())
        return True
    except ValueError:
        return False


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _argv_value(argv: list[str], flag: str) -> str | None:
    try:
        index = argv.index(flag)
    except ValueError:
        return None
    if index + 1 >= len(argv):
        return None
    return argv[index + 1]


def _same_path(left: Path, right: Path) -> bool:
    return left.resolve() == right.resolve()


def _validate_row_paths(ctx: RunContext, row: dict[str, str]) -> tuple[bool, str]:
    if row.get("allowed_for_execution") != "true":
        return False, "row is not marked allowed_for_execution=true"
    if row.get("submission_status") != "SUBMIT_READY":
        return False, "row submission_status is not SUBMIT_READY"
    try:
        argv = json.loads(row.get("vina_argv_json", ""))
    except json.JSONDecodeError:
        return False, "vina_argv_json is invalid"
    if not isinstance(argv, list) or not argv or not all(isinstance(item, str) and item for item in argv):
        return False, "vina_argv_json must be a non-empty string argv list"
    expected_executable = row.get("vina_executable", "")
    if expected_executable and argv[0] != expected_executable:
        return False, "argv executable does not match manifest vina_executable"
    ligand_path = _resolve_repo_path(ctx, row.get("ligand_pdbqt_file", ""))
    receptor_path = _resolve_repo_path(ctx, row.get("receptor_pdbqt_file", ""))
    output_path = _resolve_repo_path(ctx, row.get("output_pdbqt_file", ""))
    log_path = _resolve_repo_path(ctx, row.get("vina_log_file", ""))
    stdout_path = _resolve_repo_path(ctx, row.get("stdout_file", ""))
    stderr_path = _resolve_repo_path(ctx, row.get("stderr_file", ""))
    if not _inside(ligand_path, ctx.run_dir / "phase3_compounds"):
        return False, "ligand_pdbqt_file outside phase3_compounds"
    if not _inside(receptor_path, ctx.run_dir / "phase3_compounds" / "receptor_pdbqt"):
        return False, "receptor_pdbqt_file outside receptor_pdbqt"
    if not _inside(output_path, ctx.run_dir / "phase3_compounds" / "docking_outputs" / "focused_pocket_first"):
        return False, "output_pdbqt_file outside focused_pocket_first run directory"
    if not _inside(log_path, ctx.run_dir / "phase3_compounds" / "docking_outputs" / "focused_pocket_first"):
        return False, "vina_log_file outside focused_pocket_first run directory"
    if not _inside(stdout_path, ctx.jobs_log_dir) or not _inside(stderr_path, ctx.jobs_log_dir):
        return False, "stdout/stderr outside logs/jobs"
    for flag, expected in [
        ("--ligand", ligand_path),
        ("--receptor", receptor_path),
        ("--out", output_path),
    ]:
        value = _argv_value(argv, flag)
        if value is None:
            return False, "argv missing {0}".format(flag)
        if not _same_path(_resolve_repo_path(ctx, value), expected):
            return False, "argv {0} path does not match manifest".format(flag)
    if row.get("ligand_pdbqt_sha256") and (not ligand_path.is_file() or _sha256(ligand_path) != row["ligand_pdbqt_sha256"]):
        return False, "ligand PDBQT hash mismatch"
    if row.get("receptor_pdbqt_sha256") and (not receptor_path.is_file() or _sha256(receptor_path) != row["receptor_pdbqt_sha256"]):
        return False, "receptor PDBQT hash mismatch"
    return True, ""


def _run_one(ctx: RunContext, row: dict[str, str]) -> tuple[str, int | None, str | None]:
    stdout_path = _resolve_repo_path(ctx, row["stdout_file"])
    stderr_path = _resolve_repo_path(ctx, row["stderr_file"])
    output_path = _resolve_repo_path(ctx, row["output_pdbqt_file"])
    vina_log = _resolve_repo_path(ctx, row["vina_log_file"])
    stdout_path.parent.mkdir(parents=True, exist_ok=True)
    stderr_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    vina_log.parent.mkdir(parents=True, exist_ok=True)
    argv = json.loads(row["vina_argv_json"])
    timeout = int(row.get("timeout_seconds") or 1800)
    started = time.monotonic()
    try:
        proc = subprocess.run(argv, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False, timeout=timeout)
        stdout_path.write_text(proc.stdout or "", encoding="utf-8")
        stderr_path.write_text(proc.stderr or "", encoding="utf-8")
        if proc.returncode == 0 and not vina_log.exists() and stdout_path.is_file():
            vina_log.write_text(stdout_path.read_text(encoding="utf-8", errors="ignore"), encoding="utf-8")
        status = "PASS" if proc.returncode == 0 else "FAIL"
        append_job_status(
            ctx,
            row["job_id"],
            status,
            node=row.get("node"),
            ppn=int(row.get("ppn") or 0) or None,
            stdout=ctx.relative_to_repo(stdout_path),
            stderr=ctx.relative_to_repo(stderr_path),
            details={
                "phase": "phase3_compounds",
                "task": "M3-T5",
                "job_type": "vina_planned_row",
                "chunk_id": row["chunk_id"],
                "output_pdbqt_file": ctx.relative_to_repo(output_path),
                "vina_log_file": ctx.relative_to_repo(vina_log),
                "exit_code": proc.returncode,
                "duration_seconds": round(time.monotonic() - started, 3),
            },
        )
        return status, proc.returncode, None
    except subprocess.TimeoutExpired as exc:
        stdout_path.write_text(exc.stdout or "", encoding="utf-8")
        stderr_path.write_text(exc.stderr or "Vina chunk row timed out\n", encoding="utf-8")
        append_job_status(
            ctx,
            row["job_id"],
            "FAIL",
            node=row.get("node"),
            ppn=int(row.get("ppn") or 0) or None,
            stdout=ctx.relative_to_repo(stdout_path),
            stderr=ctx.relative_to_repo(stderr_path),
            details={"phase": "phase3_compounds", "task": "M3-T5", "job_type": "vina_planned_row", "chunk_id": row["chunk_id"], "timeout": True},
        )
        return "FAIL", None, "timeout"
    except OSError as exc:
        stderr_path.write_text("Vina row invocation failed: {0}\n".format(exc.__class__.__name__), encoding="utf-8")
        append_job_status(
            ctx,
            row["job_id"],
            "FAIL",
            node=row.get("node"),
            ppn=int(row.get("ppn") or 0) or None,
            stdout=ctx.relative_to_repo(stdout_path),
            stderr=ctx.relative_to_repo(stderr_path),
            details={"phase": "phase3_compounds", "task": "M3-T5", "job_type": "vina_planned_row", "chunk_id": row["chunk_id"], "message": "OSError"},
        )
        return "FAIL", None, "oserror"


def run_vina_chunk(
    ctx: RunContext,
    *,
    job_manifest: Path,
    chunk_id: str,
    max_workers: int = 1,
    dry_run_runner: bool = False,
) -> VinaChunkRunResult:
    manifest = job_manifest if job_manifest.is_absolute() else ctx.repo_root / job_manifest
    manifest = ctx.require_within_run_dir(manifest.resolve())
    rows = _read_csv(manifest)
    selected = [row for row in rows if row.get("chunk_id") == chunk_id]
    if not selected:
        raise ValueError("no rows assigned to requested chunk_id={0}".format(chunk_id))
    failures = 0
    for row in selected:
        ok, reason = _validate_row_paths(ctx, row)
        if not ok:
            raise ValueError("unsafe planned path for job_id={0}: {1}".format(row.get("job_id"), reason))
    worker_limit = max(1, int(max_workers or 1))
    # The contract permits pools, but sequential execution preserves strict
    # provenance and is enough for tests; PBS controls chunk-level parallelism.
    if dry_run_runner:
        for row in selected:
            append_job_status(
                ctx,
                row["job_id"],
                "PLANNED",
                node=row.get("node"),
                ppn=int(row.get("ppn") or 0) or None,
                stdout=row.get("stdout_file"),
                stderr=row.get("stderr_file"),
                details={
                    "phase": "phase3_compounds",
                    "task": "M3-T5",
                    "job_type": "vina_planned_row",
                    "chunk_id": chunk_id,
                    "dry_run_runner": True,
                    "max_workers": worker_limit,
                    "timestamp": now_iso(),
                },
            )
        return VinaChunkRunResult("PASS", chunk_id, len(selected), 0, 0, True)
    for row in selected:
        status, _, _ = _run_one(ctx, row)
        if status != "PASS":
            failures += 1
    return VinaChunkRunResult("FAIL" if failures else "PASS", chunk_id, len(selected), len(selected), failures, False)
