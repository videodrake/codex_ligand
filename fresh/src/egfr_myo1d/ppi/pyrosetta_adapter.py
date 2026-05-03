"""M2.2 PyRosetta adapter dry-run harness.

The adapter consumes M2.1 PyRosetta input specs and emits job manifests,
commands, and QC reports. It intentionally does not import PyRosetta or run
docking/relaxation in this phase.
"""

from __future__ import annotations

import csv
import importlib.util
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

from egfr_myo1d.core.logging_utils import append_job_status, append_phase_status
from egfr_myo1d.core.manifest import now_iso, write_json
from egfr_myo1d.core.run_context import RunContext, RunContextError, ensure_within
from egfr_myo1d.io.hashing import describe_file


M2_2_PHASE = "m2.2-pyrosetta-adapter-harness"
HARNESS_ROOT = "phase1_ppi/pyrosetta_adapter"
NON_GOALS_CONFIRMED = [
    "no_pyrosetta_import",
    "no_pyrosetta_docking",
    "no_pyrosetta_relaxation",
    "no_lightdock_execution",
    "no_vina_execution",
    "no_fpocket_or_p2rank_runtime",
    "no_compound_docking_or_scoring",
    "no_candidate_nomination",
    "no_qsub_submission",
    "no_cleanup_deletion",
]


JOB_FIELDS = [
    "job_name",
    "state_id",
    "method",
    "seed",
    "models_per_seed",
    "input_spec_json",
    "receptor_pdb",
    "partner_pdb",
    "mapping_csv",
    "output_dir",
    "stdout_path",
    "stderr_path",
    "execution_allowed",
    "status",
    "command",
]


@dataclass
class PyRosettaHarnessJob:
    job_name: str
    state_id: str
    method: str
    seed: int
    models_per_seed: int
    input_spec_json: str
    receptor_pdb: str
    partner_pdb: str
    mapping_csv: str
    output_dir: str
    stdout_path: str
    stderr_path: str
    execution_allowed: bool
    status: str
    command: str
    warnings: list[str] = field(default_factory=list)


@dataclass
class PyRosettaHarnessReport:
    run_id: str
    mode: str
    profile: str
    status: str
    jobs: list[PyRosettaHarnessJob] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    blockers: list[str] = field(default_factory=list)
    manifest_path: Path | None = None
    job_manifest_csv: Path | None = None
    job_manifest_jsonl: Path | None = None
    launch_script_path: Path | None = None
    qc_csv_path: Path | None = None
    summary_report_path: Path | None = None


def _mode_settings(mode: str) -> tuple[list[int], int, int | None]:
    if mode == "smoke_input" or mode == "smoke_env":
        return [0], 5, 1
    if mode == "mini":
        return [0, 1], 20, None
    if mode == "production":
        return list(range(10)), 100, None
    raise ValueError("unsupported M2.2 mode: {0}".format(mode))


def _repo_path(ctx: RunContext, path: Path) -> str:
    return ctx.relative_to_repo(path)


def _resolve_repo_or_abs(ctx: RunContext, text: str) -> Path:
    path = Path(text)
    if path.is_absolute():
        return path.resolve()
    return (ctx.repo_root / path).resolve()


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _load_m2_1_manifest(ctx: RunContext, manifest_path: Path | None) -> tuple[dict[str, Any] | None, Path]:
    path = manifest_path or (ctx.manifest_dir / "m2_1_ppi_input_manifest.json")
    path = path if path.is_absolute() else (ctx.repo_root / path)
    path = path.resolve()
    ensure_within(path, ctx.run_dir)
    if not path.is_file():
        return None, path
    return _load_json(path), path


def _pyrosetta_available() -> bool:
    return importlib.util.find_spec("pyrosetta") is not None


def _write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str], ctx: RunContext) -> None:
    safe = ctx.require_within_run_dir(path)
    safe.parent.mkdir(parents=True, exist_ok=True)
    with safe.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            encoded = {}
            for key, value in row.items():
                if isinstance(value, (list, dict)):
                    encoded[key] = json.dumps(value, sort_keys=True)
                else:
                    encoded[key] = value
            writer.writerow(encoded)


def _write_jsonl(path: Path, rows: Iterable[dict[str, Any]], ctx: RunContext) -> None:
    safe = ctx.require_within_run_dir(path)
    safe.parent.mkdir(parents=True, exist_ok=True)
    with safe.open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(row, sort_keys=True) + "\n")


def _job_dict(job: PyRosettaHarnessJob) -> dict[str, Any]:
    return {
        "job_name": job.job_name,
        "state_id": job.state_id,
        "method": job.method,
        "seed": job.seed,
        "models_per_seed": job.models_per_seed,
        "input_spec_json": job.input_spec_json,
        "receptor_pdb": job.receptor_pdb,
        "partner_pdb": job.partner_pdb,
        "mapping_csv": job.mapping_csv,
        "output_dir": job.output_dir,
        "stdout_path": job.stdout_path,
        "stderr_path": job.stderr_path,
        "execution_allowed": job.execution_allowed,
        "status": job.status,
        "command": job.command,
        "warnings": job.warnings,
    }


def _job_rows(jobs: list[PyRosettaHarnessJob]) -> list[dict[str, Any]]:
    return [_job_dict(job) for job in jobs]


def _select_packs(m2_1_manifest: dict[str, Any], states: Iterable[str] | None, max_packs: int | None) -> list[dict[str, Any]]:
    requested = set(states or [])
    packs = []
    for pack in m2_1_manifest.get("packs", []):
        if pack.get("status") not in ("PASS", "WARN"):
            continue
        if requested and pack.get("state_id") not in requested:
            continue
        if "pyrosetta_global_ppi" not in pack.get("method_specs", {}):
            continue
        packs.append(pack)
    if max_packs is not None:
        return packs[:max_packs]
    return packs


def _build_job(
    ctx: RunContext,
    pack: dict[str, Any],
    seed: int,
    models_per_seed: int,
    harness_dir: Path,
) -> PyRosettaHarnessJob:
    state_id = pack["state_id"]
    spec_path = ensure_within(
        _resolve_repo_or_abs(ctx, pack["method_specs"]["pyrosetta_global_ppi"]),
        ctx.run_dir,
    )
    spec = _load_json(spec_path)
    receptor_path = ensure_within(_resolve_repo_or_abs(ctx, spec["receptor"]["pack_path"]), ctx.run_dir)
    partner_path = ensure_within(_resolve_repo_or_abs(ctx, spec["partner"]["primary_pack_path"]), ctx.run_dir)
    mapping_path = ensure_within(_resolve_repo_or_abs(ctx, spec["receptor"]["mapping_csv"]), ctx.run_dir)
    job_name = "m2_2_pyrosetta_{0}_seed{1:03d}".format(state_id, seed)
    output_dir = ctx.require_within_run_dir(harness_dir / "outputs" / job_name)
    stdout_path = ctx.require_within_run_dir(ctx.jobs_log_dir / "{0}.stdout".format(job_name))
    stderr_path = ctx.require_within_run_dir(ctx.jobs_log_dir / "{0}.stderr".format(job_name))
    command = (
        "python -m egfr_myo1d.ppi.run_ppi_job "
        "--run-id {run_id} --job-name {job_name} --dry-run true"
    ).format(run_id=ctx.run_id, job_name=job_name)
    warnings = []
    for label, path in {
        "input_spec_json": spec_path,
        "receptor_pdb": receptor_path,
        "partner_pdb": partner_path,
        "mapping_csv": mapping_path,
    }.items():
        if not path.is_file():
            warnings.append("missing_{0}:{1}".format(label, path))
    return PyRosettaHarnessJob(
        job_name=job_name,
        state_id=state_id,
        method="pyrosetta_global_ppi",
        seed=seed,
        models_per_seed=models_per_seed,
        input_spec_json=_repo_path(ctx, spec_path),
        receptor_pdb=_repo_path(ctx, receptor_path),
        partner_pdb=_repo_path(ctx, partner_path),
        mapping_csv=_repo_path(ctx, mapping_path),
        output_dir=_repo_path(ctx, output_dir),
        stdout_path=_repo_path(ctx, stdout_path),
        stderr_path=_repo_path(ctx, stderr_path),
        execution_allowed=False,
        status="READY_DRY_RUN" if not warnings else "WARN",
        command=command,
        warnings=warnings,
    )


def _write_launch_script(ctx: RunContext, path: Path, jobs: list[PyRosettaHarnessJob]) -> Path:
    safe = ctx.require_within_run_dir(path)
    safe.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "#!/usr/bin/env bash",
        "set -euo pipefail",
        'export PYTHONPATH="$(pwd)/fresh/src:${PYTHONPATH:-}"',
        "",
        "# Dry-run commands only. These validate harness wiring without PyRosetta execution.",
    ]
    for job in jobs:
        lines.append(job.command)
    safe.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return safe


def _write_summary(ctx: RunContext, report: PyRosettaHarnessReport) -> Path:
    target = ctx.require_within_run_dir(ctx.reports_dir / "m2_2_pyrosetta_adapter.md")
    lines = [
        "# M2.2 PyRosetta Adapter Harness",
        "",
        "Status: {0}".format(report.status),
        "",
        "This phase writes adapter manifests and dry-run commands only. PyRosetta is not imported and no docking or relaxation is executed.",
        "",
        "Jobs: {0}".format(len(report.jobs)),
        "",
        "## Non-goals Confirmed",
    ]
    for item in NON_GOALS_CONFIRMED:
        lines.append("- {0}".format(item))
    target.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return target


def generate_pyrosetta_harness(
    ctx: RunContext,
    mode: str = "smoke_input",
    profile: str = "codex_dev",
    states: Iterable[str] | None = None,
    m2_1_manifest_path: Path | None = None,
    models_per_seed_override: int | None = None,
    seed_start_override: int | None = None,
    seed_count_override: int | None = None,
) -> PyRosettaHarnessReport:
    """Generate dry-run PyRosetta adapter manifests from M2.1 specs."""
    ctx.create_directories()
    harness_dir = ctx.require_within_run_dir(ctx.run_dir / HARNESS_ROOT)
    harness_dir.mkdir(parents=True, exist_ok=True)
    seeds, models_per_seed, max_packs = _mode_settings(mode)
    if models_per_seed_override is not None:
        if models_per_seed_override <= 0:
            raise ValueError("models_per_seed_override must be positive")
        models_per_seed = int(models_per_seed_override)
    if seed_start_override is not None or seed_count_override is not None:
        seed_count = int(seed_count_override) if seed_count_override is not None else len(seeds)
        seed_start = int(seed_start_override) if seed_start_override is not None else 0
        if seed_start < 0:
            raise ValueError("seed_start_override must be non-negative")
        if seed_count <= 0:
            raise ValueError("seed_count_override must be positive")
        seeds = list(range(seed_start, seed_start + seed_count))
    m2_1_manifest, manifest_path = _load_m2_1_manifest(ctx, m2_1_manifest_path)

    warnings = []
    blockers = []
    jobs: list[PyRosettaHarnessJob] = []
    if m2_1_manifest is None:
        blockers.append("missing_m2_1_manifest:{0}".format(manifest_path))
    else:
        packs = _select_packs(m2_1_manifest, states, max_packs)
        if not packs:
            blockers.append("no_ready_m2_1_pyrosetta_packs")
        for pack in packs:
            for seed in seeds:
                jobs.append(_build_job(ctx, pack, seed, models_per_seed, harness_dir))

    if not _pyrosetta_available():
        warnings.append("pyrosetta_module_not_available_adapter_manifest_only")
    for job in jobs:
        warnings.extend("{0}:{1}".format(job.job_name, warning) for warning in job.warnings)

    if blockers:
        status = "FAIL"
    elif warnings or any(job.status == "WARN" for job in jobs):
        status = "PASS_WITH_WARNINGS"
    else:
        status = "PASS"

    report = PyRosettaHarnessReport(
        run_id=ctx.run_id,
        mode=mode,
        profile=profile,
        status=status,
        jobs=jobs,
        warnings=warnings,
        blockers=blockers,
    )

    job_csv = harness_dir / "pyrosetta_job_manifest.csv"
    job_jsonl = harness_dir / "pyrosetta_job_manifest.jsonl"
    launch = harness_dir / "launch_pyrosetta_dry_runs.sh"
    _write_csv(job_csv, _job_rows(jobs), JOB_FIELDS, ctx)
    _write_jsonl(job_jsonl, _job_rows(jobs), ctx)
    _write_launch_script(ctx, launch, jobs)
    qc_path = ctx.qc_dir / "m2_2_pyrosetta_adapter_qc.csv"
    _write_csv(
        qc_path,
        [
            {
                "job_name": job.job_name,
                "state_id": job.state_id,
                "seed": job.seed,
                "models_per_seed": job.models_per_seed,
                "execution_allowed": job.execution_allowed,
                "status": job.status,
                "warnings": job.warnings,
            }
            for job in jobs
        ],
        [
            "job_name",
            "state_id",
            "seed",
            "models_per_seed",
            "execution_allowed",
            "status",
            "warnings",
        ],
        ctx,
    )
    summary = _write_summary(ctx, report)

    payload = {
        "run_id": ctx.run_id,
        "created_at": now_iso(),
        "phase": M2_2_PHASE,
        "mode": mode,
        "profile": profile,
        "status": status,
        "seed_start": min(seeds) if seeds else None,
        "seed_count": len(seeds),
        "seeds": seeds,
        "source_m2_1_manifest": _repo_path(ctx, manifest_path),
        "pyrosetta_available": _pyrosetta_available(),
        "execution_allowed": False,
        "pyrosetta_imported": False,
        "docking_executed": False,
        "relaxation_executed": False,
        "jobs": [_job_dict(job) for job in jobs],
        "warnings": warnings,
        "blockers": blockers,
        "outputs": {
            "job_manifest_csv": _repo_path(ctx, job_csv),
            "job_manifest_jsonl": _repo_path(ctx, job_jsonl),
            "launch_script": _repo_path(ctx, launch),
            "qc_csv": _repo_path(ctx, qc_path),
            "summary_report": _repo_path(ctx, summary),
        },
        "non_goals_confirmed": NON_GOALS_CONFIRMED,
    }
    manifest = ctx.manifest_dir / "m2_2_pyrosetta_adapter_manifest.json"
    write_json(manifest, payload, ctx)

    report.manifest_path = manifest
    report.job_manifest_csv = job_csv
    report.job_manifest_jsonl = job_jsonl
    report.launch_script_path = launch
    report.qc_csv_path = qc_path
    report.summary_report_path = summary

    phase_status = "FAIL" if status == "FAIL" else (
        "WARN" if status == "PASS_WITH_WARNINGS" else "PASS"
    )
    append_phase_status(
        ctx,
        M2_2_PHASE,
        phase_status,
        "M2.2 PyRosetta adapter harness {0} with {1} dry-run job(s)".format(
            status, len(jobs)
        ),
        {
            "manifest": _repo_path(ctx, manifest),
            "execution_allowed": False,
            "pyrosetta_imported": False,
            "docking_executed": False,
        },
    )
    append_job_status(
        ctx,
        "m2_2_pyrosetta_adapter_harness_local",
        phase_status,
        details={
            "message": "adapter dry-run manifest generation only",
            "job_count": len(jobs),
            "execution_allowed": False,
        },
    )
    return report


def load_job_by_name(ctx: RunContext, job_name: str) -> dict[str, Any]:
    manifest = ctx.run_dir / HARNESS_ROOT / "pyrosetta_job_manifest.jsonl"
    if not manifest.is_file():
        raise RunContextError("PyRosetta job manifest not found: {0}".format(manifest))
    with manifest.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            row = json.loads(line)
            if row.get("job_name") == job_name:
                return row
    raise ValueError("job not found in PyRosetta manifest: {0}".format(job_name))
