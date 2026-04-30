"""Manifest creation for fresh workflow runs."""

from __future__ import annotations

import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from egfr_myo1d import WORKFLOW_STAGE
from egfr_myo1d.core.run_context import RunContext
from egfr_myo1d.io.hashing import describe_file


CONFIG_FILES = {
    "fresh_run": "fresh/configs/fresh_run.yaml",
    "hpc": "fresh/configs/hpc.yaml",
    "paths": "fresh/configs/paths.yaml",
    "gates": "fresh/configs/gates.yaml",
    "receptor_states": "fresh/configs/receptor_states.yaml",
    "pocket": "fresh/configs/pocket.yaml",
    "tool_registry": "fresh/configs/tool_registry.yaml",
    "tool_envs": "fresh/configs/tool_envs.yaml",
}


EXPECTED_INPUTS = [
    {
        "input_id": "EGFR_160-185",
        "expected_path": "fresh/data/raw/receptors/EGFR_160-185.pdb",
        "required_for": ["smoke_input", "production"],
        "notes": "primary membrane-validated receptor state",
    },
    {
        "input_id": "EGFR_170-200",
        "expected_path": "fresh/data/raw/receptors/EGFR_170-200.pdb",
        "required_for": ["smoke_input", "production"],
        "notes": "primary membrane-validated receptor state",
    },
    {
        "input_id": "3GT8_raw",
        "expected_path": "fresh/data/raw/receptors/3GT8_raw.pdb",
        "required_for": ["smoke_input", "production"],
        "notes": "crystallographic reference/control state",
    },
    {
        "input_id": "plus10_full_frame",
        "expected_path": "fresh/data/raw/receptors/plus10_full_frame.pdb",
        "required_for": ["smoke_input", "production"],
        "notes": "future membrane-frame source",
    },
    {
        "input_id": "AF-O94832-F1-model_v6",
        "expected_path": "fresh/data/raw/myo1d/AF-O94832-F1-model_v6.pdb",
        "required_for": ["smoke_input", "production"],
        "notes": "future MYO1D 955-1006 source model",
    },
    {
        "input_id": "Cpd-A",
        "expected_path": "fresh/data/raw/ligands/Cpd-A.sdf",
        "required_for": ["production"],
        "notes": "confidential ligand structure; public ID only",
    },
    {
        "input_id": "Cpd-B",
        "expected_path": "fresh/data/raw/ligands/Cpd-B.sdf",
        "required_for": ["production"],
        "notes": "confidential ligand structure; public ID only",
    },
    {
        "input_id": "Cpd-C",
        "expected_path": "fresh/data/raw/ligands/Cpd-C.sdf",
        "required_for": ["production"],
        "notes": "confidential ligand structure; public ID only",
    },
]


def now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def load_yaml_config(path: Path) -> dict[str, Any]:
    try:
        import yaml
    except ImportError as exc:
        raise RuntimeError(
            "PyYAML is required to load fresh workflow YAML configs. "
            "Install pyyaml in the active Python environment."
        ) from exc

    with Path(path).open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle)
    return data or {}


def load_configs(ctx: RunContext) -> dict[str, dict[str, Any]]:
    configs: dict[str, dict[str, Any]] = {}
    for name, rel_path in CONFIG_FILES.items():
        configs[name] = load_yaml_config(ctx.repo_root / rel_path)
    return configs


def config_file_descriptions(ctx: RunContext) -> dict[str, dict[str, Any]]:
    descriptions: dict[str, dict[str, Any]] = {}
    for name, rel_path in CONFIG_FILES.items():
        info = describe_file(ctx.repo_root / rel_path)
        descriptions[name] = {
            "path": rel_path,
            "sha256": info["sha256"],
            "present": info["present"],
            "size_bytes": info["size_bytes"],
        }
    return descriptions


def write_json(path: Path, data: dict[str, Any], ctx: RunContext) -> None:
    safe_path = ctx.require_within_run_dir(path)
    safe_path.parent.mkdir(parents=True, exist_ok=True)
    safe_path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def write_run_manifest(ctx: RunContext, mode: str) -> dict[str, Any]:
    configs = load_configs(ctx)
    fresh_run = configs.get("fresh_run", {})
    manifest = {
        "run_id": ctx.run_id,
        "created_at": now_iso(),
        "mode": mode,
        "workflow_stage": fresh_run.get("stage", "foundation"),
        "workflow_version": fresh_run.get("workflow_version", WORKFLOW_STAGE),
        "repo_root": str(ctx.repo_root),
        "fresh_root": str(ctx.fresh_root),
        "run_dir": str(ctx.run_dir),
        "python_executable": sys.executable,
        "python_version": sys.version.replace("\n", " "),
        "config_files": config_file_descriptions(ctx),
        "policy": {
            "old_workflow_modified": False,
            "outputs_restricted_to_run_dir": True,
            "docking_enabled": False,
            "compound_stage_enabled": bool(fresh_run.get("compound_stage1", False)),
        },
    }
    write_json(ctx.manifest_dir / "run_manifest.json", manifest, ctx)
    return manifest


def write_input_manifest(ctx: RunContext, mode: str) -> dict[str, Any]:
    configs = load_configs(ctx)
    fresh_run = configs.get("fresh_run", {})
    compound_stage_enabled = bool(fresh_run.get("compound_stage1", False))
    inputs = []
    missing_required_inputs = []

    for entry in EXPECTED_INPUTS:
        expected_path = entry["expected_path"]
        desc = describe_file(ctx.repo_root / expected_path)
        record = {
            "input_id": entry["input_id"],
            "expected_path": expected_path,
            "required_for": entry["required_for"],
            "present": desc["present"],
            "sha256": desc["sha256"],
            "size_bytes": desc["size_bytes"],
            "mtime": desc["mtime"],
            "notes": entry["notes"],
        }
        inputs.append(record)

        is_ligand = entry["input_id"].startswith("Cpd-")
        if mode == "smoke_input" and not desc["present"]:
            if not is_ligand or compound_stage_enabled:
                missing_required_inputs.append(entry["input_id"])

    manifest = {
        "run_id": ctx.run_id,
        "created_at": now_iso(),
        "mode": mode,
        "compound_stage_enabled": compound_stage_enabled,
        "inputs": inputs,
        "missing_required_inputs": missing_required_inputs,
    }
    write_json(ctx.manifest_dir / "input_manifest.json", manifest, ctx)
    return manifest


def _git_command(ctx: RunContext, args: list[str]) -> dict[str, Any]:
    try:
        proc = subprocess.run(
            ["git", *args],
            cwd=str(ctx.repo_root),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )
    except OSError as exc:
        return {"ok": False, "stdout": "", "stderr": str(exc), "returncode": None}
    return {
        "ok": proc.returncode == 0,
        "stdout": proc.stdout.strip(),
        "stderr": proc.stderr.strip(),
        "returncode": proc.returncode,
    }


def write_git_snapshot(ctx: RunContext) -> dict[str, Any]:
    commit = _git_command(ctx, ["rev-parse", "HEAD"])
    branch = _git_command(ctx, ["branch", "--show-current"])
    status = _git_command(ctx, ["status", "--short"])
    remote = _git_command(ctx, ["config", "--get", "remote.origin.url"])

    warnings = []
    for label, result in {
        "commit": commit,
        "branch": branch,
        "status": status,
        "remote": remote,
    }.items():
        if not result["ok"]:
            warnings.append({"check": label, "message": result["stderr"] or "git command failed"})

    snapshot = {
        "created_at": now_iso(),
        "commit": commit["stdout"] if commit["ok"] else None,
        "branch": branch["stdout"] if branch["ok"] else None,
        "dirty": bool(status["stdout"]) if status["ok"] else None,
        "status_short": status["stdout"] if status["ok"] else None,
        "remote_url": remote["stdout"] if remote["ok"] else None,
        "warnings": warnings,
    }
    write_json(ctx.manifest_dir / "git_snapshot.json", snapshot, ctx)
    return snapshot


def write_environment_report(ctx: RunContext, report: dict[str, Any]) -> dict[str, Any]:
    write_json(ctx.manifest_dir / "environment_report.json", report, ctx)
    return report


def initialize_manifests(ctx: RunContext, mode: str) -> dict[str, Any]:
    run_manifest = write_run_manifest(ctx, mode)
    input_manifest = write_input_manifest(ctx, mode)
    git_snapshot = write_git_snapshot(ctx)
    environment_report = write_environment_report(
        ctx,
        {
            "run_id": ctx.run_id,
            "created_at": now_iso(),
            "mode": mode,
            "profile": None,
            "status": "NOT_RUN",
            "checks": [],
            "message": "preflight has not been run yet",
        },
    )
    return {
        "run_manifest": run_manifest,
        "input_manifest": input_manifest,
        "git_snapshot": git_snapshot,
        "environment_report": environment_report,
    }
