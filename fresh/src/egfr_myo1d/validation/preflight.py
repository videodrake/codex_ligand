"""Environment and input preflight checks for fresh runs."""

from __future__ import annotations

import importlib
import shutil
import sys
from pathlib import Path
from typing import Any

from egfr_myo1d.core.logging_utils import append_phase_status
from egfr_myo1d.core.manifest import (
    initialize_manifests,
    now_iso,
    write_environment_report,
    write_input_manifest,
)
from egfr_myo1d.core.run_context import RunContext


PROFILES = {"codex_dev", "hpc_strict"}
MODES = {"smoke_env", "smoke_input"}

IMPORT_CHECKS = [
    {"name": "pyrosetta", "module": "pyrosetta", "required_kind": "heavy"},
    {"name": "rdkit", "module": "rdkit", "required_kind": "heavy"},
    {"name": "Bio", "module": "Bio", "required_kind": "heavy"},
    {"name": "numpy", "module": "numpy", "required_kind": "core"},
    {"name": "pandas", "module": "pandas", "required_kind": "core"},
]

COMMAND_CHECKS = [
    {"name": "vina", "commands": ["vina"], "required_kind": "heavy"},
    {"name": "fpocket", "commands": ["fpocket"], "required_kind": "heavy"},
    {"name": "obabel", "commands": ["obabel"], "required_kind": "heavy"},
    {"name": "qsub", "commands": ["qsub"], "required_kind": "heavy"},
    {"name": "qstat", "commands": ["qstat"], "required_kind": "heavy"},
    {"name": "qdel", "commands": ["qdel"], "required_kind": "heavy"},
    {"name": "P2Rank/prank", "commands": ["prank", "p2rank", "P2Rank"], "required_kind": "optional"},
]


def _status_for_missing(required_kind: str, profile: str) -> str:
    if required_kind == "optional":
        return "WARN"
    if required_kind == "core":
        return "FAIL"
    if profile == "hpc_strict":
        return "FAIL"
    return "WARN"


def check_import(name: str, module: str, required_kind: str, profile: str) -> dict[str, Any]:
    try:
        imported = importlib.import_module(module)
    except Exception as exc:
        return {
            "category": "python_import",
            "name": name,
            "status": _status_for_missing(required_kind, profile),
            "present": False,
            "details": {"module": module, "error": f"{type(exc).__name__}: {exc}"},
        }
    version = getattr(imported, "__version__", None)
    return {
        "category": "python_import",
        "name": name,
        "status": "PASS",
        "present": True,
        "details": {"module": module, "version": version},
    }


def check_command(name: str, commands: list[str], required_kind: str, profile: str) -> dict[str, Any]:
    for command in commands:
        found = shutil.which(command)
        if found:
            return {
                "category": "external_command",
                "name": name,
                "status": "PASS",
                "present": True,
                "details": {"command": command, "path": found},
            }
    return {
        "category": "external_command",
        "name": name,
        "status": _status_for_missing(required_kind, profile),
        "present": False,
        "details": {"commands": commands, "error": "not found on PATH"},
    }


def python_check() -> dict[str, Any]:
    major_minor = sys.version_info[:2]
    status = "PASS" if major_minor >= (3, 9) else "FAIL"
    return {
        "category": "python",
        "name": "python",
        "status": status,
        "present": True,
        "details": {
            "executable": sys.executable,
            "version": sys.version.replace("\n", " "),
        },
    }


def environment_checks(profile: str) -> list[dict[str, Any]]:
    checks = [python_check()]
    for item in IMPORT_CHECKS:
        checks.append(check_import(item["name"], item["module"], item["required_kind"], profile))
    for item in COMMAND_CHECKS:
        checks.append(check_command(item["name"], item["commands"], item["required_kind"], profile))
    return checks


def input_checks(ctx: RunContext, mode: str) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    input_manifest = write_input_manifest(ctx, mode)
    checks: list[dict[str, Any]] = []
    compound_stage_enabled = bool(input_manifest.get("compound_stage_enabled", False))

    for entry in input_manifest["inputs"]:
        input_id = entry["input_id"]
        present = bool(entry["present"])
        is_ligand = input_id.startswith("Cpd-")
        if present:
            status = "PASS"
        elif mode == "smoke_env":
            status = "WARN"
        elif mode == "smoke_input" and is_ligand and not compound_stage_enabled:
            status = "WARN"
        else:
            status = "FAIL"

        checks.append(
            {
                "category": "input_file",
                "name": input_id,
                "status": status,
                "present": present,
                "details": {
                    "expected_path": entry["expected_path"],
                    "required_for": entry["required_for"],
                    "sha256": entry["sha256"],
                    "size_bytes": entry["size_bytes"],
                },
            }
        )

    return checks, input_manifest


def summarize_status(checks: list[dict[str, Any]]) -> str:
    statuses = {check["status"] for check in checks}
    if "FAIL" in statuses:
        return "FAIL"
    if "WARN" in statuses:
        return "WARN"
    return "PASS"


def run_preflight(ctx: RunContext, mode: str, profile: str = "codex_dev") -> dict[str, Any]:
    if mode not in MODES:
        raise ValueError(f"unsupported preflight mode: {mode}")
    if profile not in PROFILES:
        raise ValueError(f"unsupported preflight profile: {profile}")

    initialize_manifests(ctx, mode)
    checks = environment_checks(profile)
    file_checks, input_manifest = input_checks(ctx, mode)
    checks.extend(file_checks)
    status = summarize_status(checks)
    counts = {
        "PASS": sum(1 for check in checks if check["status"] == "PASS"),
        "WARN": sum(1 for check in checks if check["status"] == "WARN"),
        "FAIL": sum(1 for check in checks if check["status"] == "FAIL"),
    }

    report = {
        "run_id": ctx.run_id,
        "created_at": now_iso(),
        "mode": mode,
        "profile": profile,
        "status": status,
        "python_executable": sys.executable,
        "python_version": sys.version.replace("\n", " "),
        "checks": checks,
        "counts": counts,
        "input_manifest": {
            "path": ctx.relative_to_repo(ctx.manifest_dir / "input_manifest.json"),
            "missing_required_inputs": input_manifest.get("missing_required_inputs", []),
        },
        "policy": {
            "codex_dev_heavy_missing_is_warn": profile == "codex_dev",
            "hpc_strict_heavy_missing_is_fail": profile == "hpc_strict",
            "p2rank_missing_is_warn": True,
            "docking_enabled": False,
        },
    }
    write_environment_report(ctx, report)
    append_phase_status(
        ctx,
        phase="preflight",
        status=status,
        message=f"preflight completed with {status}",
        details={"mode": mode, "profile": profile, "counts": counts},
    )
    return report
