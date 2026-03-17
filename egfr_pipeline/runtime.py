"""Runtime resource resolution, lane manifests, and scratch helpers."""

from __future__ import annotations

import json
import os
import shutil
import socket
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Mapping, Optional


THREAD_ENV_VARS = (
    "OMP_NUM_THREADS",
    "MKL_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
)


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def _safe_int(value: object) -> Optional[int]:
    if value in (None, "", False):
        return None
    try:
        parsed = int(str(value).strip())
    except (TypeError, ValueError):
        return None
    return parsed if parsed > 0 else None


def resolve_visible_gpu_ids(
    environ: Optional[Mapping[str, str]] = None,
) -> list[str]:
    env = environ or os.environ
    raw = str(env.get("CUDA_VISIBLE_DEVICES", "")).strip()
    if not raw or raw in {"-1", "none", "None"}:
        return []
    return [item.strip() for item in raw.split(",") if item.strip()]


@dataclass(frozen=True)
class RuntimeResources:
    requested_cpus: Optional[int]
    allocated_cpus: int
    thread_limit: Optional[int]
    effective_cpus: int
    visible_gpu_ids: list[str]
    pbs_job_id: str
    host: str
    source: str


def resolve_runtime_resources(
    *,
    requested_cpus: Optional[int] = None,
    allocated_cpus: Optional[int] = None,
    environ: Optional[Mapping[str, str]] = None,
) -> RuntimeResources:
    """Resolve safe CPU/GPU runtime limits from scheduler and thread env."""
    env = environ or os.environ
    source = "default"

    allocated = _safe_int(allocated_cpus)
    if allocated is not None:
        source = "override"
    else:
        allocated = _safe_int(env.get("PBS_NP"))
        if allocated is not None:
            source = "pbs_np"
        else:
            thread_candidates = [_safe_int(env.get(name)) for name in THREAD_ENV_VARS]
            thread_candidates = [value for value in thread_candidates if value]
            if thread_candidates:
                allocated = min(thread_candidates)
                source = "thread_env"
            else:
                allocated = max(1, os.cpu_count() or 1)

    requested = _safe_int(requested_cpus)
    thread_values = [_safe_int(env.get(name)) for name in THREAD_ENV_VARS]
    thread_values = [value for value in thread_values if value]
    thread_limit = min(thread_values) if thread_values else None

    effective = allocated
    if requested is not None:
        effective = min(effective, requested)
    if thread_limit is not None:
        effective = min(effective, thread_limit)
    effective = max(1, effective)

    return RuntimeResources(
        requested_cpus=requested,
        allocated_cpus=max(1, allocated),
        thread_limit=thread_limit,
        effective_cpus=effective,
        visible_gpu_ids=resolve_visible_gpu_ids(env),
        pbs_job_id=str(env.get("PBS_JOBID", "")),
        host=str(env.get("HOSTNAME") or socket.gethostname()),
        source=source,
    )


def cap_worker_count(
    *,
    requested_workers: int,
    n_jobs: int,
    cpu_per_job: int,
    available_cpus: int,
) -> int:
    workers = min(max(1, int(requested_workers)), max(1, int(n_jobs)))
    if cpu_per_job <= 0:
        return 1
    budget = max(1, int(available_cpus)) // max(1, int(cpu_per_job))
    return max(1, min(workers, budget))


def lane_identity(
    lane: str,
    *,
    state: Optional[str] = None,
    seed: Optional[int] = None,
    engine: Optional[str] = None,
) -> str:
    parts = [lane]
    if state:
        parts.append(state)
    if seed is not None:
        parts.append(f"seed{seed}")
    if engine:
        parts.append(engine)
    return "__".join(parts)


def runtime_root(project_root: Path) -> Path:
    return Path(project_root) / "runtime"


def lane_runtime_dir(
    project_root: Path,
    lane: str,
    *,
    state: Optional[str] = None,
    seed: Optional[int] = None,
    engine: Optional[str] = None,
) -> Path:
    return runtime_root(project_root) / "lanes" / lane_identity(
        lane,
        state=state,
        seed=seed,
        engine=engine,
    )


def lane_marker_path(
    project_root: Path,
    lane: str,
    *,
    state: Optional[str] = None,
    seed: Optional[int] = None,
    engine: Optional[str] = None,
) -> Path:
    return lane_runtime_dir(
        project_root,
        lane,
        state=state,
        seed=seed,
        engine=engine,
    ) / "lane_complete.json"


def lane_manifest_path(
    project_root: Path,
    lane: str,
    *,
    state: Optional[str] = None,
    seed: Optional[int] = None,
    engine: Optional[str] = None,
) -> Path:
    return lane_runtime_dir(
        project_root,
        lane,
        state=state,
        seed=seed,
        engine=engine,
    ) / "runtime_manifest.json"


def lane_is_complete(
    project_root: Path,
    lane: str,
    *,
    state: Optional[str] = None,
    seed: Optional[int] = None,
    engine: Optional[str] = None,
) -> bool:
    marker = lane_marker_path(
        project_root,
        lane,
        state=state,
        seed=seed,
        engine=engine,
    )
    if not marker.exists():
        return False
    try:
        payload = json.loads(marker.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return False
    return payload.get("status") == "completed"


def copy_config_snapshot(config_path: Path, target_dir: Path) -> Optional[Path]:
    try:
        target_dir.mkdir(parents=True, exist_ok=True)
        destination = target_dir / config_path.name
        shutil.copy2(config_path, destination)
        return destination
    except OSError:
        return None


def write_lane_manifest(
    project_root: Path,
    lane: str,
    *,
    resources: RuntimeResources,
    state: Optional[str] = None,
    seed: Optional[int] = None,
    engine: Optional[str] = None,
    config_path: Optional[Path] = None,
    extra: Optional[dict] = None,
) -> Path:
    runtime_dir = lane_runtime_dir(
        project_root,
        lane,
        state=state,
        seed=seed,
        engine=engine,
    )
    runtime_dir.mkdir(parents=True, exist_ok=True)
    snapshot_path = copy_config_snapshot(config_path, runtime_dir) if config_path else None
    payload = {
        "lane": lane,
        "state": state or "",
        "seed": seed,
        "engine": engine or "",
        "started_at": utc_now_iso(),
        "resources": asdict(resources),
        "source_config": str(config_path) if config_path else "",
        "config_snapshot": str(snapshot_path) if snapshot_path else "",
    }
    if extra:
        payload.update(extra)
    manifest_path = runtime_dir / "runtime_manifest.json"
    manifest_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return manifest_path


def finalize_lane_manifest(
    project_root: Path,
    lane: str,
    *,
    status: str,
    state: Optional[str] = None,
    seed: Optional[int] = None,
    engine: Optional[str] = None,
    extra: Optional[dict] = None,
) -> Path:
    manifest_path = lane_manifest_path(
        project_root,
        lane,
        state=state,
        seed=seed,
        engine=engine,
    )
    payload: dict = {}
    if manifest_path.exists():
        try:
            payload = json.loads(manifest_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            payload = {}

    payload["status"] = status
    payload["completed_at"] = utc_now_iso()
    if extra:
        payload.update(extra)

    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return manifest_path


def write_lane_completion(
    project_root: Path,
    lane: str,
    *,
    status: str,
    state: Optional[str] = None,
    seed: Optional[int] = None,
    engine: Optional[str] = None,
    extra: Optional[dict] = None,
) -> Path:
    marker = lane_marker_path(
        project_root,
        lane,
        state=state,
        seed=seed,
        engine=engine,
    )
    marker.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "lane": lane,
        "state": state or "",
        "seed": seed,
        "engine": engine or "",
        "status": status,
        "completed_at": utc_now_iso(),
    }
    if extra:
        payload.update(extra)
    marker.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    finalize_lane_manifest(
        project_root,
        lane,
        status=status,
        state=state,
        seed=seed,
        engine=engine,
        extra=extra,
    )
    return marker


def resolve_scratch_base(
    *,
    scratch_dir: Optional[Path] = None,
    environ: Optional[Mapping[str, str]] = None,
) -> Optional[Path]:
    env = environ or os.environ
    if scratch_dir:
        return Path(scratch_dir)
    tmpdir = env.get("TMPDIR") or env.get("TMP")
    if tmpdir:
        return Path(tmpdir)
    return None


def prepare_scratch_workspace(
    project_root: Path,
    lane: str,
    *,
    state: Optional[str] = None,
    seed: Optional[int] = None,
    engine: Optional[str] = None,
    scratch_dir: Optional[Path] = None,
    environ: Optional[Mapping[str, str]] = None,
) -> Optional[Path]:
    base = resolve_scratch_base(scratch_dir=scratch_dir, environ=environ)
    if base is None:
        return None
    scratch = base / "ppi_vina_runtime" / lane_identity(
        lane,
        state=state,
        seed=seed,
        engine=engine,
    )
    scratch.mkdir(parents=True, exist_ok=True)
    runtime_dir = lane_runtime_dir(
        project_root,
        lane,
        state=state,
        seed=seed,
        engine=engine,
    )
    runtime_dir.mkdir(parents=True, exist_ok=True)
    return scratch


def sync_tree(src: Path, dst: Path) -> None:
    src = Path(src)
    dst = Path(dst)
    if not src.exists():
        return
    dst.mkdir(parents=True, exist_ok=True)
    for item in src.iterdir():
        target = dst / item.name
        if item.is_dir():
            sync_tree(item, target)
        else:
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(item, target)


def reset_lane_completion(
    project_root: Path,
    lane: str,
    *,
    state: Optional[str] = None,
    seed: Optional[int] = None,
    engine: Optional[str] = None,
) -> None:
    marker = lane_marker_path(
        project_root,
        lane,
        state=state,
        seed=seed,
        engine=engine,
    )
    if marker.exists():
        marker.unlink()
