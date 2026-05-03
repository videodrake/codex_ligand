"""HPC tool discovery and smoke checks for PPI-surface integrations.

The checks are intentionally non-installing. They record what is already
available on the active HPC environment before Milestone 3 relies on it.
"""

from __future__ import annotations

import importlib
import json
import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from egfr_myo1d.core.logging_utils import append_phase_status
from egfr_myo1d.core.manifest import load_yaml_config, now_iso, write_json
from egfr_myo1d.core.run_context import RunContext


DISCOVER_MODES = {"discover", "smoke"}
PROFILES = {"codex_dev", "hpc_strict"}

STATUS_OK = {"available", "installed_and_smoke_passed", "optional_disabled", "external_server_only"}


@dataclass
class ToolCheck:
    name: str
    required_level: str
    env: str
    status: str
    path: str | None = None
    smoke_test: str = "not_run"
    version: str | None = None
    details: dict[str, Any] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)
    blockers: list[str] = field(default_factory=list)


def _is_required(level: str) -> bool:
    return level in {"core", "core_addition"}


def _read_registry(ctx: RunContext, registry_path: Path | None = None) -> dict[str, Any]:
    path = registry_path or ctx.repo_root / "fresh" / "configs" / "tool_registry.yaml"
    return load_yaml_config(path).get("tools", {})


def _import_tool(name: str, module: str, required_level: str, profile: str) -> ToolCheck:
    try:
        imported = importlib.import_module(module)
    except Exception as exc:
        status = "not_installed"
        blockers = []
        warnings = [f"import_failed:{type(exc).__name__}:{exc}"]
        if _is_required(required_level) and profile == "hpc_strict":
            status = "installed_but_smoke_failed"
            blockers = warnings[:]
        return ToolCheck(
            name=name,
            required_level=required_level,
            env="",
            status=status,
            smoke_test="failed",
            details={"module": module, "error": f"{type(exc).__name__}: {exc}"},
            warnings=warnings,
            blockers=blockers,
        )
    return ToolCheck(
        name=name,
        required_level=required_level,
        env="",
        status="installed_and_smoke_passed",
        smoke_test="passed",
        version=getattr(imported, "__version__", None),
        details={"module": module},
    )


def _command_version_line(output: str) -> str | None:
    for line in output.splitlines():
        stripped = line.strip()
        if stripped:
            return stripped[:240]
    return None


def _run_command_smoke(command_path: str, args: list[str], timeout: int) -> tuple[str, str, int | None, str | None]:
    try:
        proc = subprocess.run(
            [command_path, *args],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            timeout=timeout,
            check=False,
        )
    except subprocess.TimeoutExpired as exc:
        output = (exc.stdout or "") + "\nTIMEOUT"
        return output, "timeout", None, None
    except OSError as exc:
        return f"{type(exc).__name__}: {exc}", "failed", None, None

    output = proc.stdout or ""
    # Several HPC tools print help/version and return non-zero. Treat any
    # non-empty output as a smoke pass for help-style probes.
    smoke = "passed" if output.strip() or proc.returncode == 0 else "failed"
    return output, smoke, proc.returncode, _command_version_line(output)


def _resolve_command(command: str) -> str | None:
    if "/" in command:
        path = Path(command)
        return str(path) if path.exists() else None
    return shutil.which(command)


def _command_tool(
    name: str,
    commands: list[str],
    required_level: str,
    profile: str,
    smoke_args: list[str],
    mode: str,
    timeout: int,
) -> ToolCheck:
    for command in commands:
        found = _resolve_command(command)
        if not found:
            continue
        if mode == "discover":
            return ToolCheck(
                name=name,
                required_level=required_level,
                env="",
                status="available",
                path=found,
                smoke_test="not_run",
                details={"command": command},
            )
        output, smoke, returncode, version = _run_command_smoke(found, smoke_args, timeout)
        status = "installed_and_smoke_passed" if smoke == "passed" else "installed_but_smoke_failed"
        warnings = [] if smoke == "passed" else [f"smoke_failed:returncode={returncode}"]
        blockers = warnings[:] if _is_required(required_level) and profile == "hpc_strict" else []
        return ToolCheck(
            name=name,
            required_level=required_level,
            env="",
            status=status,
            path=found,
            smoke_test=smoke,
            version=version,
            details={"command": command, "returncode": returncode, "smoke_args": smoke_args, "output_preview": output[:2000]},
            warnings=warnings,
            blockers=blockers,
        )

    status = "not_installed"
    message = f"not_found_on_path:{','.join(commands)}"
    warnings = [message]
    blockers: list[str] = []
    if _is_required(required_level) and profile == "hpc_strict":
        blockers = [message]
    return ToolCheck(
        name=name,
        required_level=required_level,
        env="",
        status=status,
        smoke_test="not_run",
        details={"commands": commands},
        warnings=warnings,
        blockers=blockers,
    )


def _check_external_path(path_value: str | None) -> dict[str, Any]:
    if not path_value:
        return {"path": None, "present": False}
    path = Path(path_value)
    return {"path": path_value, "present": path.exists(), "is_dir": path.is_dir()}


def _check_one_tool(name: str, spec: dict[str, Any], mode: str, profile: str, timeout: int) -> ToolCheck:
    required_level = str(spec.get("required_level", "optional"))
    external_path = _check_external_path(spec.get("external_path"))
    if spec.get("status_override"):
        result = ToolCheck(
            name=name,
            required_level=required_level,
            env=str(spec.get("env", "")),
            status=str(spec["status_override"]),
            smoke_test="not_run",
            details={"reason": "status_override", "external_path": external_path},
        )
        if "install_candidate" in spec:
            result.details["install_candidate"] = spec["install_candidate"]
        if "license_note" in spec:
            result.details["license_note"] = spec["license_note"]
        return result

    command_probes: list[str] = []
    if spec.get("binary"):
        command_probes.append(str(spec["binary"]))
    command_probes.extend(str(item) for item in spec.get("commands", []))

    if spec.get("python_import") and not command_probes:
        result = _import_tool(name, str(spec["python_import"]), required_level, profile)
    elif command_probes:
        result = _command_tool(
            name=name,
            commands=command_probes,
            required_level=required_level,
            profile=profile,
            smoke_args=[str(item) for item in spec.get("smoke_args", [])],
            mode=mode,
            timeout=timeout,
        )
    elif external_path["present"]:
        result = ToolCheck(
            name=name,
            required_level=required_level,
            env=str(spec.get("env", "")),
            status="available",
            path=str(external_path["path"]),
            smoke_test="not_run",
            details={"reason": "external_path_present"},
        )
    else:
        status = "optional_disabled" if not _is_required(required_level) else "not_installed"
        warning = "no import, command, binary, or present external_path configured"
        result = ToolCheck(
            name=name,
            required_level=required_level,
            env=str(spec.get("env", "")),
            status=status,
            smoke_test="not_run",
            details={"reason": warning},
            warnings=[warning],
            blockers=[warning] if _is_required(required_level) and profile == "hpc_strict" else [],
        )

    result.env = str(spec.get("env", ""))
    result.details["external_path"] = external_path
    if "install_candidate" in spec:
        result.details["install_candidate"] = spec["install_candidate"]
    if "license_note" in spec:
        result.details["license_note"] = spec["license_note"]
    return result


def _overall_status(checks: list[ToolCheck]) -> str:
    if any(item.blockers for item in checks):
        return "FAIL"
    if any(item.status not in STATUS_OK for item in checks):
        return "WARN"
    return "PASS"


def _tool_to_dict(tool: ToolCheck) -> dict[str, Any]:
    return {
        "required_level": tool.required_level,
        "env": tool.env,
        "status": tool.status,
        "path": tool.path,
        "smoke_test": tool.smoke_test,
        "version": tool.version,
        "details": tool.details,
        "warnings": tool.warnings,
        "blockers": tool.blockers,
    }


def _write_tool_log(ctx: RunContext, tool: ToolCheck) -> None:
    log_dir = ctx.logs_dir / "tools"
    ctx.require_within_run_dir(log_dir)
    log_dir.mkdir(parents=True, exist_ok=True)
    log_path = ctx.require_within_run_dir(log_dir / f"{tool.name}.smoke.log")
    log_path.write_text(json.dumps(_tool_to_dict(tool), indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _write_phase_log(ctx: RunContext, payload: dict[str, Any]) -> None:
    log_path = ctx.require_within_run_dir(ctx.logs_dir / "phase_tool_preflight.log")
    log_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _write_report(ctx: RunContext, payload: dict[str, Any]) -> Path:
    path = ctx.require_within_run_dir(ctx.reports_dir / "tool_installation_report.md")
    lines = [
        "# Tool Installation / Runtime Preflight Report",
        "",
        f"- run_id: `{payload['run_id']}`",
        f"- mode: `{payload['mode']}`",
        f"- profile: `{payload['profile']}`",
        f"- status: `{payload['status']}`",
        "",
        "| Tool | Required level | Env | Status | Smoke | Path/version |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for name, item in payload["tools"].items():
        path_version = item.get("path") or item.get("version") or ""
        lines.append(
            "| {0} | {1} | {2} | {3} | {4} | {5} |".format(
                name,
                item.get("required_level", ""),
                item.get("env", ""),
                item.get("status", ""),
                item.get("smoke_test", ""),
                str(path_version).replace("|", "/"),
            )
        )
    lines.extend(["", "## Blockers", ""])
    blockers = payload.get("blockers", [])
    if blockers:
        lines.extend(f"- {item}" for item in blockers)
    else:
        lines.append("- None")
    lines.extend(["", "## Notes", "", "- This command does not install packages or submit jobs."])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def run_tool_preflight(
    ctx: RunContext,
    mode: str = "discover",
    profile: str = "codex_dev",
    registry_path: Path | None = None,
    timeout: int = 15,
) -> dict[str, Any]:
    if mode not in DISCOVER_MODES:
        raise ValueError(f"unsupported tool preflight mode: {mode}")
    if profile not in PROFILES:
        raise ValueError(f"unsupported tool preflight profile: {profile}")

    registry = _read_registry(ctx, registry_path)
    tools = [
        _check_one_tool(name, dict(spec or {}), mode=mode, profile=profile, timeout=timeout)
        for name, spec in registry.items()
    ]
    for tool in tools:
        _write_tool_log(ctx, tool)

    blockers = [f"{tool.name}:{blocker}" for tool in tools for blocker in tool.blockers]
    payload = {
        "run_id": ctx.run_id,
        "created_at": now_iso(),
        "mode": mode,
        "profile": profile,
        "python_executable": sys.executable,
        "status": _overall_status(tools),
        "tools": {tool.name: _tool_to_dict(tool) for tool in tools},
        "counts": {
            "total": len(tools),
            "available_or_passed": sum(1 for tool in tools if tool.status in STATUS_OK),
            "not_installed": sum(1 for tool in tools if tool.status == "not_installed"),
            "smoke_failed": sum(1 for tool in tools if tool.status == "installed_but_smoke_failed"),
            "blockers": len(blockers),
        },
        "blockers": blockers,
    }
    write_json(ctx.manifest_dir / "tool_status.json", payload, ctx)
    _write_phase_log(ctx, payload)
    report_path = _write_report(ctx, payload)
    append_phase_status(
        ctx,
        phase="tool-preflight",
        status=payload["status"],
        message=f"tool preflight completed with {payload['status']}",
        details={
            "mode": mode,
            "profile": profile,
            "tool_count": len(tools),
            "blockers": blockers,
            "report": ctx.relative_to_repo(report_path),
        },
    )
    return payload
