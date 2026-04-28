"""Run directory context for fresh workflow commands."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
from pathlib import Path


class RunContextError(ValueError):
    """Raised when a run context would be unsafe or invalid."""


def utc_timestamp_for_id() -> str:
    return datetime.utcnow().strftime("%Y%m%d_%H%M%S")


def generate_run_id(mode: str) -> str:
    prefixes = {
        "smoke_env": "smoke_env",
        "smoke_input": "smoke_input",
        "mini": "mini",
        "production": "prod",
        "prod": "prod",
    }
    prefix = prefixes.get(mode, mode or "run")
    return f"{prefix}_{utc_timestamp_for_id()}"


def validate_run_id(run_id: str) -> str:
    if not run_id:
        raise RunContextError("run_id must not be empty")
    if ".." in run_id or "/" in run_id or "\\" in run_id:
        raise RunContextError("run_id must not contain '..', '/', or '\\'")
    allowed = set("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_.-")
    bad = sorted({char for char in run_id if char not in allowed})
    if bad:
        raise RunContextError(f"run_id contains unsupported characters: {''.join(bad)}")
    return run_id


def discover_repo_root(repo_root: Path | None = None) -> Path:
    if repo_root is not None:
        root = Path(repo_root).resolve()
        if not (root / "fresh").is_dir():
            raise RunContextError(f"repo_root does not contain fresh/: {root}")
        return root

    cwd = Path.cwd().resolve()
    for candidate in [cwd, *cwd.parents]:
        if (candidate / "fresh" / "configs" / "fresh_run.yaml").is_file():
            return candidate

    return Path(__file__).resolve().parents[4]


def ensure_within(child: Path, parent: Path) -> Path:
    child_resolved = Path(child).resolve()
    parent_resolved = Path(parent).resolve()
    try:
        child_resolved.relative_to(parent_resolved)
    except ValueError as exc:
        raise RunContextError(f"path escapes allowed directory: {child_resolved}") from exc
    return child_resolved


@dataclass(frozen=True)
class RunContext:
    repo_root: Path
    fresh_root: Path
    run_id: str
    run_dir: Path
    manifest_dir: Path
    logs_dir: Path
    jobs_log_dir: Path
    errors_dir: Path
    qc_dir: Path
    reports_dir: Path
    scratch_dir: Path
    tmp_dir: Path

    @classmethod
    def create(
        cls,
        run_id: str | None,
        mode: str,
        repo_root: Path | None = None,
    ) -> "RunContext":
        root = discover_repo_root(repo_root)
        fresh_root = (root / "fresh").resolve()
        final_run_id = validate_run_id(run_id or generate_run_id(mode))
        runs_root = ensure_within(fresh_root / "runs", fresh_root)
        run_dir = ensure_within(runs_root / final_run_id, runs_root)
        ctx = cls._from_paths(root, fresh_root, final_run_id, run_dir)
        ctx.create_directories()
        return ctx

    @classmethod
    def for_existing(
        cls,
        run_id: str,
        repo_root: Path | None = None,
    ) -> "RunContext":
        root = discover_repo_root(repo_root)
        fresh_root = (root / "fresh").resolve()
        final_run_id = validate_run_id(run_id)
        runs_root = ensure_within(fresh_root / "runs", fresh_root)
        run_dir = ensure_within(runs_root / final_run_id, runs_root)
        if not run_dir.is_dir():
            raise RunContextError(f"run directory does not exist: {run_dir}")
        return cls._from_paths(root, fresh_root, final_run_id, run_dir)

    @classmethod
    def _from_paths(
        cls,
        repo_root: Path,
        fresh_root: Path,
        run_id: str,
        run_dir: Path,
    ) -> "RunContext":
        return cls(
            repo_root=repo_root,
            fresh_root=fresh_root,
            run_id=run_id,
            run_dir=run_dir,
            manifest_dir=run_dir / "manifest",
            logs_dir=run_dir / "logs",
            jobs_log_dir=run_dir / "logs" / "jobs",
            errors_dir=run_dir / "logs" / "errors",
            qc_dir=run_dir / "qc",
            reports_dir=run_dir / "reports",
            scratch_dir=run_dir / "scratch",
            tmp_dir=run_dir / "tmp",
        )

    def create_directories(self) -> None:
        for path in [
            self.run_dir,
            self.manifest_dir,
            self.logs_dir,
            self.jobs_log_dir,
            self.errors_dir,
            self.qc_dir,
            self.reports_dir,
            self.scratch_dir,
            self.tmp_dir,
        ]:
            self.require_within_run_dir(path)
            path.mkdir(parents=True, exist_ok=True)

    def require_within_run_dir(self, path: Path) -> Path:
        return ensure_within(path, self.run_dir)

    def relative_to_repo(self, path: Path) -> str:
        try:
            return Path(path).resolve().relative_to(self.repo_root).as_posix()
        except ValueError:
            return str(Path(path).resolve())
