"""M1 Phase 6: PBS job-file generation per handoff §10.

Generates concrete PBS files with absolute stdout/stderr paths, PYTHONPATH
exports, all five thread limits set to 1, and conda activation of the
``pyrosetta`` environment. Mode → ppn / walltime mapping is read from
``fresh/configs/hpc.yaml`` (no hardcoding).

This module never calls qsub. PBS submission is a user-side HPC step.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

from egfr_myo1d.core.logging_utils import append_job_status, append_phase_status
from egfr_myo1d.core.manifest import load_yaml_config
from egfr_myo1d.core.run_context import RunContext, RunContextError
from egfr_myo1d.io.hashing import describe_file


THREAD_ENV_KEYS = (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
)

DEFAULT_CONDA_SH = "/usr/local/anaconda/3/2023.09/etc/profile.d/conda.sh"
DEFAULT_CONDA_ENV = "pyrosetta"

WALLTIME_PATTERN = re.compile(r"^\d+:\d{2}:\d{2}$")


@dataclass
class PBSFile:
    job_name: str
    node: str
    ppn: int
    walltime: str
    queue: str
    mode: str
    output_path: Path
    content: str
    sha256: str | None = None
    run_id: str = ""
    profile: str = "codex_dev"
    warnings: list[str] = field(default_factory=list)
    status: str = "PASS"


# ---------------------------------------------------------------------------
# hpc.yaml loader
# ---------------------------------------------------------------------------


def load_hpc_gate(ctx: RunContext) -> dict[str, Any]:
    return load_yaml_config(ctx.repo_root / "fresh" / "configs" / "hpc.yaml")


def _parse_conda_paths(activate_command: str | None) -> tuple[str, str]:
    """Extract (conda_sh, env_name) from an activate_command string of the form
    ``"source <conda_sh> && conda activate <env>"``. Falls back to defaults
    when the string is absent or malformed."""
    if not activate_command:
        return DEFAULT_CONDA_SH, DEFAULT_CONDA_ENV
    sh = DEFAULT_CONDA_SH
    env = DEFAULT_CONDA_ENV
    m_sh = re.search(r"source\s+(\S+)", activate_command)
    if m_sh:
        sh = m_sh.group(1)
    m_env = re.search(r"conda\s+activate\s+(\S+)", activate_command)
    if m_env:
        env = m_env.group(1)
    return sh, env


# ---------------------------------------------------------------------------
# Pure rendering function (used by golden-file tests)
# ---------------------------------------------------------------------------


def render_pbs_content(
    job_name: str,
    node: str,
    ppn: int,
    walltime: str,
    queue: str,
    repo_root: str,
    conda_sh: str,
    conda_env: str,
    thread_limits: dict[str, int],
    stdout_path: str,
    stderr_path: str,
    run_id: str,
    command_lines: Iterable[str],
) -> str:
    thread_lines = "\n".join(
        "export {0}={1}".format(key, int(thread_limits.get(key, 1)))
        for key in THREAD_ENV_KEYS
    )
    body = "\n".join(line.rstrip() for line in command_lines)
    if not body:
        body = "# (no command lines provided for mode={0})".format("placeholder")
    return (
        "#!/bin/bash\n"
        "#PBS -q {queue}\n"
        "#PBS -N {job_name}\n"
        "#PBS -l nodes={node}:ppn={ppn}\n"
        "#PBS -l walltime={walltime}\n"
        "#PBS -V\n"
        "#PBS -o {stdout_path}\n"
        "#PBS -e {stderr_path}\n"
        "\n"
        "set -eo pipefail\n"
        "cd {repo_root}\n"
        "\n"
        "set +u\n"
        "source {conda_sh}\n"
        "conda activate {conda_env}\n"
        "set -u\n"
        "\n"
        "export LANG=C.UTF-8\n"
        "export LC_ALL=C.UTF-8\n"
        "export PYTHONIOENCODING=utf-8\n"
        "export REPO_ROOT=\"{repo_root}\"\n"
        "export PYTHONPATH=\"$REPO_ROOT/fresh/src:${{PYTHONPATH:-}}\"\n"
        "\n"
        "{thread_lines}\n"
        "\n"
        "{body}\n"
    ).format(
        queue=queue,
        job_name=job_name,
        node=node,
        ppn=ppn,
        walltime=walltime,
        stdout_path=stdout_path,
        stderr_path=stderr_path,
        repo_root=repo_root,
        conda_sh=conda_sh,
        conda_env=conda_env,
        thread_lines=thread_lines,
        body=body,
    )


# ---------------------------------------------------------------------------
# Mode-specific command lines
# ---------------------------------------------------------------------------


def derive_smoke_env_command_lines(run_id: str) -> list[str]:
    """Smoke A: env-only PBS body (handoff §10.2)."""
    return [
        "python -m egfr_myo1d.cli preflight --run-id {0} --mode smoke_env --profile hpc_strict".format(run_id),
        "python -m egfr_myo1d.cli cleanup --run-id {0} --mode test --dry-run false".format(run_id),
    ]


def derive_smoke_input_command_lines(run_id: str, input_root: str = "fresh/data/raw") -> list[str]:
    """Smoke B: input-prep PBS body (handoff §10.2 + Phase 8 prepare-inputs)."""
    return [
        "python -m egfr_myo1d.cli prepare-inputs --run-id {0} --mode smoke_input --profile hpc_strict --input-root {1}".format(
            run_id, input_root
        ),
        "python -m egfr_myo1d.cli cleanup --run-id {0} --mode test --dry-run false".format(run_id),
    ]


def derive_command_lines_for_mode(mode: str, run_id: str, input_root: str = "fresh/data/raw") -> list[str]:
    if mode == "smoke_env":
        return derive_smoke_env_command_lines(run_id)
    if mode == "smoke_input":
        return derive_smoke_input_command_lines(run_id, input_root)
    # mini / scaling / production: M1 does not prescribe specific commands;
    # the user fills these in for actual M2/M3 execution.
    return [
        "# placeholder for {0} mode; actual M2/M3 commands inserted by user".format(mode),
        "python -m egfr_myo1d.cli status --run-id {0}".format(run_id),
    ]


# ---------------------------------------------------------------------------
# generate_pbs entry point
# ---------------------------------------------------------------------------


def _validate_walltime(walltime: str) -> None:
    if not WALLTIME_PATTERN.match(walltime):
        raise ValueError("walltime must match HH:MM:SS pattern; got: {0}".format(walltime))


def generate_pbs(
    ctx: RunContext,
    job_name: str,
    mode: str,
    node: str | None = None,
    ppn: int | None = None,
    walltime: str | None = None,
    output_path: Path | None = None,
    profile: str = "codex_dev",
    input_root: str = "fresh/data/raw",
    command_lines: Iterable[str] | None = None,
) -> PBSFile:
    """Resolve defaults from hpc.yaml, render a concrete PBS file, and write
    it under ``ctx.run_dir`` (default ``runs/<run_id>/scripts/<job_name>.pbs``).

    Never calls qsub. Records a job_status entry with status=``GENERATED``.
    """
    if profile not in {"codex_dev", "hpc_strict"}:
        raise ValueError("profile must be codex_dev or hpc_strict; got: {0}".format(profile))
    if not job_name or not re.match(r"^[A-Za-z0-9_.\-]+$", job_name):
        raise ValueError("job_name must be alnum/underscore/dot/dash; got: {0}".format(job_name))

    config = load_hpc_gate(ctx)
    queue = config.get("queue", "workq")
    nodes = list(config.get("nodes", ["node04"]))
    ppn_map = config.get("ppn", {})
    walltime_map = config.get("walltime", {})
    thread_limits = dict(config.get("thread_limits", {}))
    for key in THREAD_ENV_KEYS:
        thread_limits.setdefault(key, 1)
    conda_sh, conda_env = _parse_conda_paths(
        config.get("conda", {}).get("activate_command")
    )
    if "env_name" in config.get("conda", {}):
        conda_env = config["conda"]["env_name"]

    # hpc.yaml lookup key: smoke_env/smoke_input both map to "smoke"
    lookup_mode = "smoke" if mode in {"smoke_env", "smoke_input"} else mode
    if lookup_mode not in ppn_map:
        raise ValueError(
            "mode '{0}' (lookup '{1}') has no ppn entry in hpc.yaml; known: {2}".format(
                mode, lookup_mode, sorted(ppn_map.keys())
            )
        )
    if lookup_mode not in walltime_map:
        raise ValueError(
            "mode '{0}' (lookup '{1}') has no walltime entry in hpc.yaml; known: {2}".format(
                mode, lookup_mode, sorted(walltime_map.keys())
            )
        )

    resolved_node = node if node else nodes[0]
    if resolved_node not in nodes:
        raise ValueError(
            "node '{0}' is not in hpc.yaml nodes={1}".format(resolved_node, nodes)
        )

    warnings: list[str] = []
    resolved_ppn = int(ppn) if ppn is not None else int(ppn_map[lookup_mode])
    if ppn is not None and int(ppn) != int(ppn_map[lookup_mode]):
        warnings.append(
            "ppn_override:requested={0} hpc_yaml_default={1}".format(ppn, ppn_map[lookup_mode])
        )

    resolved_walltime = walltime if walltime else str(walltime_map[lookup_mode])
    _validate_walltime(resolved_walltime)
    if walltime is not None and walltime != str(walltime_map[lookup_mode]):
        warnings.append(
            "walltime_override:requested={0} hpc_yaml_default={1}".format(
                walltime, walltime_map[lookup_mode]
            )
        )

    if output_path is None:
        target = ctx.run_dir / "scripts" / "{0}.pbs".format(job_name)
    else:
        candidate = Path(output_path)
        if not candidate.is_absolute():
            target = ctx.run_dir / candidate
        else:
            target = candidate

    safe_target = ctx.require_within_run_dir(target)
    safe_target.parent.mkdir(parents=True, exist_ok=True)

    repo_root = str(ctx.repo_root.resolve())
    stdout_path = "{0}/fresh/runs/{1}/logs/jobs/{2}.stdout".format(
        repo_root, ctx.run_id, job_name
    )
    stderr_path = stdout_path[:-7] + ".stderr"

    cmds = list(command_lines) if command_lines is not None else derive_command_lines_for_mode(
        mode, ctx.run_id, input_root=input_root
    )

    content = render_pbs_content(
        job_name=job_name,
        node=resolved_node,
        ppn=resolved_ppn,
        walltime=resolved_walltime,
        queue=queue,
        repo_root=repo_root,
        conda_sh=conda_sh,
        conda_env=conda_env,
        thread_limits=thread_limits,
        stdout_path=stdout_path,
        stderr_path=stderr_path,
        run_id=ctx.run_id,
        command_lines=cmds,
    )

    # Force Unix line endings; HPC is Linux and Path.write_text on Windows
    # would otherwise translate \n -> \r\n.
    with safe_target.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(content)
    descriptor = describe_file(safe_target)

    pbs = PBSFile(
        job_name=job_name,
        node=resolved_node,
        ppn=resolved_ppn,
        walltime=resolved_walltime,
        queue=queue,
        mode=mode,
        output_path=safe_target,
        content=content,
        sha256=descriptor.get("sha256"),
        run_id=ctx.run_id,
        profile=profile,
        warnings=warnings,
        status="WARN" if warnings else "PASS",
    )

    append_job_status(
        ctx,
        job_name=job_name,
        status="GENERATED",
        node=resolved_node,
        ppn=resolved_ppn,
        stdout=stdout_path,
        stderr=stderr_path,
        details={
            "mode": mode,
            "walltime": resolved_walltime,
            "queue": queue,
            "pbs_path": str(safe_target),
            "sha256": descriptor.get("sha256"),
            "profile": profile,
            "warnings": warnings,
        },
    )

    append_phase_status(
        ctx,
        phase="prepare-pbs",
        status=pbs.status,
        message="prepare-pbs job_name={0} mode={1} node={2} ppn={3}".format(
            job_name, mode, resolved_node, resolved_ppn
        ),
        details={
            "job_name": job_name,
            "mode": mode,
            "node": resolved_node,
            "ppn": resolved_ppn,
            "walltime": resolved_walltime,
            "queue": queue,
            "pbs_path": str(safe_target),
            "warnings": warnings,
            "profile": profile,
        },
    )

    return pbs
