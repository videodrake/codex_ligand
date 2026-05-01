"""PBS planner for real M2.2 PyRosetta PPI execution.

The legacy workflow used the scheduler allocation as the upper bound for local
parallel work: PBS provided ``PBS_NP`` and the PyRosetta pipeline resolved that
into an effective worker count. This planner keeps the same contract for the
fresh workflow, while using subprocess chunks for stronger PyRosetta process
isolation and restartability.

This module never submits qsub and never imports PyRosetta.
"""

from __future__ import annotations

import csv
import json
import math
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from egfr_myo1d.core.logging_utils import append_job_status, append_phase_status
from egfr_myo1d.core.manifest import load_yaml_config, now_iso, write_json
from egfr_myo1d.core.run_context import RunContext, ensure_within
from egfr_myo1d.hpc.pbs import (
    DEFAULT_CONDA_ENV,
    DEFAULT_CONDA_SH,
    THREAD_ENV_KEYS,
    render_pbs_content,
)


DEFAULT_JOB_MANIFEST = "phase1_ppi/pyrosetta_adapter/pyrosetta_job_manifest.jsonl"
REAL_QSUB_ROOT = "phase1_ppi/pyrosetta_adapter/qsub"
REAL_CHUNK_MANIFEST = "phase1_ppi/pyrosetta_adapter/pyrosetta_real_chunk_manifest.csv"
REAL_PBS_MANIFEST = "phase1_ppi/pyrosetta_adapter/pyrosetta_real_pbs_manifest.csv"
REAL_PLAN_REPORT = "m2_2_pyrosetta_real_pbs_plan.md"
REAL_PLAN_JSON = "m2_2_pyrosetta_real_pbs_plan.json"
AUTO_MIN_MODELS_PER_CHUNK = 25

CHUNK_FIELDS = [
    "chunk_id",
    "run_id",
    "profile",
    "node",
    "pbs_job_name",
    "job_name",
    "state_id",
    "method",
    "seed",
    "models_per_seed",
    "model_start",
    "model_count",
    "model_end",
    "output_dir",
    "stdout_file",
    "stderr_file",
    "command",
    "chunk_status",
]

PBS_FIELDS = [
    "pbs_job_name",
    "run_id",
    "profile",
    "node",
    "ppn",
    "worker_count",
    "assigned_chunk_count",
    "assigned_model_count",
    "pbs_file",
    "pbs_stdout_file",
    "pbs_stderr_file",
    "queue",
    "walltime",
    "conda_env",
    "python_executable",
    "chunk_manifest",
    "qsub_command",
    "submit_by_default",
    "pbs_generation_status",
    "pbs_notes",
]


@dataclass
class PyRosettaRealPbsPlan:
    run_id: str
    profile: str
    status: str
    chunk_rows: list[dict[str, Any]] = field(default_factory=list)
    pbs_rows: list[dict[str, Any]] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    blockers: list[str] = field(default_factory=list)
    chunk_manifest_csv: Path | None = None
    pbs_manifest_csv: Path | None = None
    manifest_path: Path | None = None
    report_path: Path | None = None
    qsub_dir: Path | None = None
    counts: dict[str, int] = field(default_factory=dict)
    hpc: dict[str, Any] = field(default_factory=dict)


def _repo_path(ctx: RunContext, path: Path) -> str:
    return ctx.relative_to_repo(path)


def _resolve_repo_or_abs(ctx: RunContext, value: str | Path) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path.resolve()
    return (ctx.repo_root / path).resolve()


def _read_jsonl(path: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                rows.append(json.loads(line))
    return rows


def _write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str], ctx: RunContext) -> None:
    safe = ctx.require_within_run_dir(path)
    safe.parent.mkdir(parents=True, exist_ok=True)
    with safe.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _load_jobs(ctx: RunContext, job_manifest: Path | None) -> tuple[list[dict[str, Any]], Path, list[str]]:
    path = job_manifest or (ctx.run_dir / DEFAULT_JOB_MANIFEST)
    path = _resolve_repo_or_abs(ctx, path)
    ensure_within(path, ctx.run_dir)
    if not path.is_file():
        return [], path, ["missing_m2_2_job_manifest:{0}".format(path)]
    return _read_jsonl(path), path, []


def _hpc_values(
    ctx: RunContext,
    *,
    profile: str,
    nodes: list[str] | None = None,
    ppn: int | None = None,
    walltime: str | None = None,
    queue: str | None = None,
    conda_env: str | None = None,
    conda_sh: str | None = None,
    python_executable: str | None = None,
) -> dict[str, Any]:
    hpc = load_yaml_config(ctx.repo_root / "fresh" / "configs" / "hpc.yaml")
    configured_nodes = [str(item) for item in hpc.get("nodes", [])]
    selected_nodes = nodes or configured_nodes or ["node04", "node05", "node06"]
    selected_ppn = int(ppn or hpc.get("ppn", {}).get(profile, hpc.get("ppn", {}).get("production", 32)))
    selected_walltime = str(
        walltime
        or hpc.get("walltime", {}).get(profile, hpc.get("walltime", {}).get("production", "999:00:00"))
    )
    conda = hpc.get("conda", {})
    activate = str(conda.get("activate_command") or "")
    parsed_conda_sh = DEFAULT_CONDA_SH
    m_sh = re.search(r"source\s+(\S+)", activate)
    if m_sh:
        parsed_conda_sh = m_sh.group(1)
    python = hpc.get("python", {})
    return {
        "scheduler": hpc.get("scheduler", "PBS"),
        "queue": str(queue or hpc.get("queue", "workq")),
        "nodes": selected_nodes,
        "ppn": selected_ppn,
        "walltime": selected_walltime,
        "conda_env": str(conda_env or conda.get("env_name", DEFAULT_CONDA_ENV)),
        "conda_sh": str(conda_sh or parsed_conda_sh),
        "python_executable": str(python_executable or python.get("path_hint", "python")),
        "thread_limits": {
            key: int((hpc.get("thread_limits", {}) or {}).get(key, 1))
            for key in THREAD_ENV_KEYS
        },
    }


def _sanitize_id(text: str) -> str:
    keep = []
    for char in str(text):
        if char.isalnum() or char in "._-":
            keep.append(char)
        else:
            keep.append("_")
    return "".join(keep).strip("_") or "item"


def _build_chunk_rows(
    ctx: RunContext,
    *,
    jobs: list[dict[str, Any]],
    profile: str,
    nodes: list[str],
    models_per_chunk: int | None,
    target_workers_per_seed: int,
) -> list[dict[str, Any]]:
    chunk_rows: list[dict[str, Any]] = []
    chunk_index = 0
    if models_per_chunk is not None and models_per_chunk <= 0:
        raise ValueError("models_per_chunk must be positive")
    for job in sorted(jobs, key=lambda item: (str(item.get("state_id", "")), int(item.get("seed", 0)), str(item.get("job_name", "")))):
        models_per_seed = int(job.get("models_per_seed", 0))
        if models_per_seed <= 0:
            continue
        chunk_size = int(models_per_chunk) if models_per_chunk else min(
            models_per_seed,
            max(
                AUTO_MIN_MODELS_PER_CHUNK,
                math.ceil(models_per_seed / max(1, target_workers_per_seed)),
            ),
        )
        job_name = str(job["job_name"])
        for model_start in range(1, models_per_seed + 1, chunk_size):
            chunk_index += 1
            model_count = min(chunk_size, models_per_seed - model_start + 1)
            model_end = model_start + model_count - 1
            node = nodes[(chunk_index - 1) % len(nodes)]
            pbs_job_name = "m2_pyrosetta_{0}".format(_sanitize_id(node))
            chunk_id = "{0}_m{1:05d}_{2:05d}".format(
                _sanitize_id(job_name),
                model_start,
                model_end,
            )
            stdout_file = ctx.jobs_log_dir / "{0}.stdout".format(chunk_id)
            stderr_file = ctx.jobs_log_dir / "{0}.stderr".format(chunk_id)
            command = (
                "python -m egfr_myo1d.ppi.run_ppi_job "
                "--run-id {run_id} --job-name {job_name} --dry-run false --allow-real true "
                "--model-start {model_start} --max-models {model_count} --chunk-id {chunk_id}"
            ).format(
                run_id=ctx.run_id,
                job_name=job_name,
                model_start=model_start,
                model_count=model_count,
                chunk_id=chunk_id,
            )
            chunk_rows.append(
                {
                    "chunk_id": chunk_id,
                    "run_id": ctx.run_id,
                    "profile": profile,
                    "node": node,
                    "pbs_job_name": pbs_job_name,
                    "job_name": job_name,
                    "state_id": str(job.get("state_id", "")),
                    "method": str(job.get("method", "")),
                    "seed": str(job.get("seed", "")),
                    "models_per_seed": str(models_per_seed),
                    "model_start": str(model_start),
                    "model_count": str(model_count),
                    "model_end": str(model_end),
                    "output_dir": str(job.get("output_dir", "")),
                    "stdout_file": _repo_path(ctx, stdout_file),
                    "stderr_file": _repo_path(ctx, stderr_file),
                    "command": command,
                    "chunk_status": "PLANNED",
                }
            )
    return chunk_rows


def _shell_quote(text: str) -> str:
    return "'" + str(text).replace("'", "'\"'\"'") + "'"


def _pbs_body_for_rows(ctx: RunContext, rows: list[dict[str, Any]], ppn: int) -> list[str]:
    lines = [
        "# Real M2.2 PyRosetta execution. One PBS job uses PBS_NP as the local worker cap,",
        "# mirroring the legacy allocated-cpus/PBS_NP contract while isolating PyRosetta per chunk.",
        'MAX_PARALLEL="${PBS_NP:-%d}"' % ppn,
        "if [ \"$MAX_PARALLEL\" -gt %d ]; then MAX_PARALLEL=%d; fi" % (ppn, ppn),
        "if [ \"$MAX_PARALLEL\" -lt 1 ]; then MAX_PARALLEL=1; fi",
        "echo \"MAX_PARALLEL=$MAX_PARALLEL\"",
        "failures=0",
        "",
        "wait_for_slot() {",
        "  while [ \"$(jobs -rp | wc -l)\" -ge \"$MAX_PARALLEL\" ]; do",
        "    if ! wait -n; then failures=$((failures + 1)); fi",
        "  done",
        "}",
        "",
        "run_chunk() {",
        "  local job_name=\"$1\"",
        "  local chunk_id=\"$2\"",
        "  local model_start=\"$3\"",
        "  local model_count=\"$4\"",
        "  local stdout_file=\"$5\"",
        "  local stderr_file=\"$6\"",
        "  mkdir -p \"$(dirname \"$stdout_file\")\" \"$(dirname \"$stderr_file\")\"",
        "  python -m egfr_myo1d.ppi.run_ppi_job \\",
        "    --run-id %s \\" % _shell_quote(ctx.run_id),
        "    --job-name \"$job_name\" \\",
        "    --dry-run false \\",
        "    --allow-real true \\",
        "    --model-start \"$model_start\" \\",
        "    --max-models \"$model_count\" \\",
        "    --chunk-id \"$chunk_id\" \\",
        "    > \"$stdout_file\" 2> \"$stderr_file\"",
        "}",
        "",
    ]
    if not rows:
        lines.append("echo 'No assigned chunks for this compatibility PBS file.'")
        return lines
    for row in rows:
        lines.extend(
            [
                "wait_for_slot",
                "run_chunk {job_name} {chunk_id} {model_start} {model_count} {stdout_file} {stderr_file} &".format(
                    job_name=_shell_quote(row["job_name"]),
                    chunk_id=_shell_quote(row["chunk_id"]),
                    model_start=_shell_quote(row["model_start"]),
                    model_count=_shell_quote(row["model_count"]),
                    stdout_file=_shell_quote(str((ctx.repo_root / row["stdout_file"]).resolve())),
                    stderr_file=_shell_quote(str((ctx.repo_root / row["stderr_file"]).resolve())),
                ),
            ]
        )
    lines.extend(
        [
            "",
            "while [ \"$(jobs -rp | wc -l)\" -gt 0 ]; do",
            "  if ! wait -n; then failures=$((failures + 1)); fi",
            "done",
            "if [ \"$failures\" -ne 0 ]; then",
            "  echo \"FAILED_CHUNKS=$failures\" >&2",
            "  exit 1",
            "fi",
            "echo 'All assigned PyRosetta chunks completed.'",
        ]
    )
    return lines


def _write_summary(ctx: RunContext, plan: PyRosettaRealPbsPlan) -> Path:
    target = ctx.require_within_run_dir(ctx.reports_dir / REAL_PLAN_REPORT)
    lines = [
        "# M2.2 PyRosetta Real PBS Plan",
        "",
        "Status: {0}".format(plan.status),
        "Profile: {0}".format(plan.profile),
        "PBS files: {0}".format(len(plan.pbs_rows)),
        "Chunks: {0}".format(len(plan.chunk_rows)),
        "Models planned: {0}".format(plan.counts.get("planned_models", 0)),
        "",
        "This plan mirrors the legacy resource contract: PBS grants cores via `PBS_NP`, and the runner caps local concurrent work at that value.",
        "Fresh executes independent subprocess chunks instead of sharing PyRosetta objects across a Python multiprocessing pool.",
        "",
        "## Nodes",
    ]
    for row in plan.pbs_rows:
        lines.append(
            "- {node}: ppn={ppn}, chunks={chunks}, models={models}, file={pbs}".format(
                node=row["node"],
                ppn=row["ppn"],
                chunks=row["assigned_chunk_count"],
                models=row["assigned_model_count"],
                pbs=row["pbs_file"],
            )
        )
    if plan.warnings:
        lines.extend(["", "## Warnings"])
        for warning in plan.warnings:
            lines.append("- {0}".format(warning))
    if plan.blockers:
        lines.extend(["", "## Blockers"])
        for blocker in plan.blockers:
            lines.append("- {0}".format(blocker))
    lines.extend(
        [
            "",
            "## Submission",
            "",
            "No qsub command was invoked by this planner. Inspect the PBS files and submit manually on HPC.",
        ]
    )
    target.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return target


def plan_pyrosetta_real_pbs(
    ctx: RunContext,
    *,
    profile: str = "production",
    job_manifest: Path | None = None,
    nodes: list[str] | None = None,
    ppn: int | None = None,
    walltime: str | None = None,
    queue: str | None = None,
    conda_env: str | None = None,
    conda_sh: str | None = None,
    python_executable: str | None = None,
    models_per_chunk: int | None = None,
) -> PyRosettaRealPbsPlan:
    """Generate node-level PBS scripts for chunked real PyRosetta execution."""

    ctx.create_directories()
    jobs, job_manifest_path, blockers = _load_jobs(ctx, job_manifest)
    warnings: list[str] = []
    hpc = _hpc_values(
        ctx,
        profile=profile,
        nodes=nodes,
        ppn=ppn,
        walltime=walltime,
        queue=queue,
        conda_env=conda_env,
        conda_sh=conda_sh,
        python_executable=python_executable,
    )
    if not hpc["nodes"]:
        blockers.append("no_hpc_nodes_configured")
    if int(hpc["ppn"]) <= 0:
        blockers.append("ppn_must_be_positive")
    if models_per_chunk is not None and models_per_chunk <= 0:
        blockers.append("models_per_chunk_must_be_positive")

    real_jobs = [
        job for job in jobs if str(job.get("method", "")) == "pyrosetta_global_ppi"
    ]
    if jobs and not real_jobs:
        blockers.append("no_pyrosetta_global_ppi_jobs_in_manifest")
    if not jobs:
        blockers.append("zero_m2_2_jobs")

    chunk_rows: list[dict[str, Any]] = []
    if not blockers:
        chunk_rows = _build_chunk_rows(
            ctx,
            jobs=real_jobs,
            profile=profile,
            nodes=hpc["nodes"],
            models_per_chunk=models_per_chunk,
            target_workers_per_seed=int(hpc["ppn"]),
        )
        if not chunk_rows:
            blockers.append("zero_planned_chunks")

    qsub_dir = ctx.require_within_run_dir(ctx.run_dir / REAL_QSUB_ROOT)
    qsub_dir.mkdir(parents=True, exist_ok=True)
    chunk_manifest = ctx.require_within_run_dir(ctx.run_dir / REAL_CHUNK_MANIFEST)
    pbs_manifest = ctx.require_within_run_dir(ctx.run_dir / REAL_PBS_MANIFEST)

    pbs_rows: list[dict[str, Any]] = []
    grouped: dict[str, list[dict[str, Any]]] = {node: [] for node in hpc["nodes"]}
    for row in chunk_rows:
        grouped.setdefault(row["node"], []).append(row)

    for node in hpc["nodes"]:
        assigned = grouped.get(node, [])
        pbs_job_name = "m2_pyrosetta_{0}".format(_sanitize_id(node))
        pbs_file = qsub_dir / "{0}.pbs".format(pbs_job_name)
        stdout_file = ctx.jobs_log_dir / "{0}.stdout".format(pbs_job_name)
        stderr_file = ctx.jobs_log_dir / "{0}.stderr".format(pbs_job_name)
        assigned_model_count = sum(int(row["model_count"]) for row in assigned)
        body = _pbs_body_for_rows(ctx, assigned, int(hpc["ppn"]))
        content = render_pbs_content(
            job_name=pbs_job_name,
            node=node,
            ppn=int(hpc["ppn"]),
            walltime=str(hpc["walltime"]),
            queue=str(hpc["queue"]),
            repo_root=str(ctx.repo_root.resolve()),
            conda_sh=str(hpc["conda_sh"]),
            conda_env=str(hpc["conda_env"]),
            thread_limits=hpc["thread_limits"],
            stdout_path=str(stdout_file.resolve()),
            stderr_path=str(stderr_file.resolve()),
            run_id=ctx.run_id,
            command_lines=body,
        )
        with ctx.require_within_run_dir(pbs_file).open("w", encoding="utf-8", newline="\n") as handle:
            handle.write(content)
        pbs_rows.append(
            {
                "pbs_job_name": pbs_job_name,
                "run_id": ctx.run_id,
                "profile": profile,
                "node": node,
                "ppn": str(hpc["ppn"]),
                "worker_count": str(min(int(hpc["ppn"]), max(1, len(assigned)))),
                "assigned_chunk_count": str(len(assigned)),
                "assigned_model_count": str(assigned_model_count),
                "pbs_file": _repo_path(ctx, pbs_file),
                "pbs_stdout_file": _repo_path(ctx, stdout_file),
                "pbs_stderr_file": _repo_path(ctx, stderr_file),
                "queue": str(hpc["queue"]),
                "walltime": str(hpc["walltime"]),
                "conda_env": str(hpc["conda_env"]),
                "python_executable": str(hpc["python_executable"]),
                "chunk_manifest": _repo_path(ctx, chunk_manifest),
                "qsub_command": "qsub {0}".format(_repo_path(ctx, pbs_file)) if assigned else "",
                "submit_by_default": "false",
                "pbs_generation_status": "PASS" if assigned else "EMPTY_COMPATIBILITY",
                "pbs_notes": "submit manually; planner did not invoke qsub" if assigned else "no chunks assigned to this node",
            }
        )
        append_job_status(
            ctx,
            pbs_job_name,
            "PLANNED" if assigned else "EMPTY_COMPATIBILITY",
            node=node,
            ppn=int(hpc["ppn"]),
            stdout=_repo_path(ctx, stdout_file),
            stderr=_repo_path(ctx, stderr_file),
            details={
                "phase": "m2.2-pyrosetta-real-pbs-plan",
                "job_type": "pyrosetta_chunked_pbs",
                "assigned_chunk_count": len(assigned),
                "assigned_model_count": assigned_model_count,
                "pbs_file": _repo_path(ctx, pbs_file),
                "qsub_command": "qsub {0}".format(_repo_path(ctx, pbs_file)) if assigned else "",
                "submitted": False,
            },
        )

    _write_csv(chunk_manifest, chunk_rows, CHUNK_FIELDS, ctx)
    _write_csv(pbs_manifest, pbs_rows, PBS_FIELDS, ctx)

    planned_models = sum(int(row["model_count"]) for row in chunk_rows)
    expected_models = sum(int(job.get("models_per_seed", 0)) for job in real_jobs)
    if planned_models != expected_models:
        blockers.append("planned_model_count_mismatch:{0}!={1}".format(planned_models, expected_models))

    pbs_files_ok = all((ctx.repo_root / row["pbs_file"]).is_file() for row in pbs_rows)
    pbs_thread_ok = all(
        all("export {0}=1".format(key) in (ctx.repo_root / row["pbs_file"]).read_text(encoding="utf-8") for key in THREAD_ENV_KEYS)
        for row in pbs_rows
    )
    pbs_np_ok = all(
        "PBS_NP" in (ctx.repo_root / row["pbs_file"]).read_text(encoding="utf-8")
        for row in pbs_rows
    )
    pbs_ppn_ok = all(int(row["ppn"]) == int(hpc["ppn"]) for row in pbs_rows)
    for ok, msg in [
        (pbs_files_ok, "pbs_files_not_written"),
        (pbs_thread_ok, "pbs_thread_limits_missing"),
        (pbs_np_ok, "pbs_np_worker_cap_missing"),
        (pbs_ppn_ok, "pbs_ppn_mismatch"),
    ]:
        if not ok:
            blockers.append(msg)

    status = "FAIL" if blockers else ("PASS_WITH_WARNINGS" if warnings else "PASS")
    plan = PyRosettaRealPbsPlan(
        run_id=ctx.run_id,
        profile=profile,
        status=status,
        chunk_rows=chunk_rows,
        pbs_rows=pbs_rows,
        warnings=warnings,
        blockers=blockers,
        chunk_manifest_csv=chunk_manifest,
        pbs_manifest_csv=pbs_manifest,
        qsub_dir=qsub_dir,
        counts={
            "manifest_jobs": len(jobs),
            "pyrosetta_jobs": len(real_jobs),
            "planned_chunks": len(chunk_rows),
            "planned_models": planned_models,
            "pbs_files": len(pbs_rows),
            "models_per_chunk": int(models_per_chunk) if models_per_chunk else 0,
        },
        hpc={
            "nodes": hpc["nodes"],
            "ppn": int(hpc["ppn"]),
            "walltime": hpc["walltime"],
            "queue": hpc["queue"],
            "worker_contract": "PBS_NP capped at ppn; one subprocess per active chunk",
            "chunking": "automatic approximately one wave per seed when models_per_chunk=0",
            "legacy_reference": "run_production.py --allocated-cpus ${PBS_NP}; PipelineManager multiprocessing.Pool(n_cpus)",
        },
    )
    plan.report_path = _write_summary(ctx, plan)
    plan.manifest_path = ctx.require_within_run_dir(ctx.manifest_dir / REAL_PLAN_JSON)
    write_json(
        plan.manifest_path,
        {
            "phase": "m2.2-pyrosetta-real-pbs-plan",
            "run_id": ctx.run_id,
            "timestamp": now_iso(),
            "status": status,
            "profile": profile,
            "source_job_manifest": _repo_path(ctx, job_manifest_path),
            "chunk_manifest": _repo_path(ctx, chunk_manifest),
            "pbs_manifest": _repo_path(ctx, pbs_manifest),
            "qsub_dir": _repo_path(ctx, qsub_dir),
            "counts": plan.counts,
            "hpc": plan.hpc,
            "warnings": warnings,
            "blockers": blockers,
            "qsub_invoked": False,
            "pyrosetta_imported_by_planner": False,
            "docking_executed_by_planner": False,
        },
        ctx,
    )
    append_phase_status(
        ctx,
        "m2.2-pyrosetta-real-pbs-plan",
        status,
        "M2.2 real PyRosetta PBS plan generated without qsub submission",
        {
            "chunk_manifest": _repo_path(ctx, chunk_manifest),
            "pbs_manifest": _repo_path(ctx, pbs_manifest),
            "planned_chunks": len(chunk_rows),
            "planned_models": planned_models,
            "qsub_invoked": False,
        },
    )
    return plan
