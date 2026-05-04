import json
import os
import uuid
from pathlib import Path

import pytest

from egfr_myo1d.cli import main
from egfr_myo1d.core.logging_utils import initialize_logs
from egfr_myo1d.core.run_context import RunContext
from egfr_myo1d.tools import tool_preflight


REPO_ROOT = Path(__file__).resolve().parents[2]


def unique_run_id(prefix="tool_preflight"):
    return "{}_{}".format(prefix, uuid.uuid4().hex[:12])


def write_registry(path):
    path.write_text(
        "\n".join(
            [
                "tools:",
                "  numpy:",
                "    required_level: core",
                "    env: pyrosetta",
                "    python_import: numpy",
                "  vina:",
                "    required_level: core",
                "    env: pyrosetta",
                "    commands: [vina]",
                "    smoke_args: ['--version']",
                "  pyKVFinder:",
                "    required_level: optional",
                "    env: ppi_surface_or_pyrosetta",
                "    python_import: pyKVFinder",
                "  passer:",
                "    required_level: external_server_only",
                "    env: external_server",
                "    status_override: external_server_only",
                "",
            ]
        ),
        encoding="utf-8",
    )


def test_tool_preflight_writes_status_logs_and_report(tmp_path, monkeypatch):
    registry = tmp_path / "tool_registry.yaml"
    write_registry(registry)

    def fake_import(name):
        if name == "numpy":
            class Module:
                __version__ = "1.0"

            return Module()
        raise ImportError("mock missing")

    def fake_which(command):
        if command == "vina":
            return "/mock/bin/vina"
        return None

    def fake_run(args, stdout, stderr, text, timeout, check):
        class Proc:
            returncode = 0
            stdout = "AutoDock Vina mock\n"

        return Proc()

    monkeypatch.setattr(tool_preflight.importlib, "import_module", fake_import)
    monkeypatch.setattr(tool_preflight.shutil, "which", fake_which)
    monkeypatch.setattr(tool_preflight.subprocess, "run", fake_run)

    ctx = RunContext.create(unique_run_id(), "smoke_env", repo_root=REPO_ROOT)
    initialize_logs(ctx)
    report = tool_preflight.run_tool_preflight(
        ctx,
        mode="smoke",
        profile="codex_dev",
        registry_path=registry,
    )

    assert report["status"] == "WARN"
    assert report["tools"]["numpy"]["status"] == "installed_and_smoke_passed"
    assert report["tools"]["vina"]["smoke_test"] == "passed"
    assert report["tools"]["pyKVFinder"]["status"] == "not_installed"
    assert report["tools"]["passer"]["status"] == "external_server_only"
    assert (ctx.manifest_dir / "tool_status.json").is_file()
    assert (ctx.logs_dir / "phase_tool_preflight.log").is_file()
    assert (ctx.logs_dir / "tools" / "vina.smoke.log").is_file()
    assert (ctx.reports_dir / "tool_installation_report.md").is_file()


def test_tool_preflight_hpc_strict_blocks_missing_core_but_not_optional(tmp_path, monkeypatch):
    registry = tmp_path / "tool_registry.yaml"
    write_registry(registry)

    monkeypatch.setattr(tool_preflight.importlib, "import_module", lambda name: (_ for _ in ()).throw(ImportError("missing")))
    monkeypatch.setattr(tool_preflight.shutil, "which", lambda command: None)

    ctx = RunContext.create(unique_run_id(), "smoke_env", repo_root=REPO_ROOT)
    initialize_logs(ctx)
    report = tool_preflight.run_tool_preflight(
        ctx,
        mode="discover",
        profile="hpc_strict",
        registry_path=registry,
    )

    assert report["status"] == "FAIL"
    assert any(item.startswith("numpy:") for item in report["blockers"])
    assert not any(item.startswith("pyKVFinder:") for item in report["blockers"])
    assert report["tools"]["pyKVFinder"]["required_level"] == "optional"


def test_cli_tool_preflight_help_includes_modes(capsys):
    code = main(["tool-preflight", "--help"])
    out = capsys.readouterr().out

    assert code == 0
    assert "discover" in out
    assert "smoke" in out


def test_preflight_tools_script_wraps_tool_preflight_cli():
    script = REPO_ROOT / "fresh" / "scripts" / "preflight_tools.sh"

    assert script.is_file()
    content = script.read_text(encoding="utf-8")
    assert "egfr_myo1d.cli tool-preflight" in content
    assert "init-run" in content


def test_default_registry_surfaces_ppi_surface_handoff_tools(monkeypatch):
    monkeypatch.setattr(tool_preflight.importlib, "import_module", lambda name: (_ for _ in ()).throw(ImportError("missing")))
    monkeypatch.setattr(tool_preflight.shutil, "which", lambda command: None)

    ctx = RunContext.create(unique_run_id(), "smoke_env", repo_root=REPO_ROOT)
    initialize_logs(ctx)
    report = tool_preflight.run_tool_preflight(ctx, mode="discover", profile="codex_dev")

    expected = {
        "fpocket",
        "mdpocket",
        "vina",
        "obabel",
        "rdkit",
        "pyrosetta",
        "gromacs",
        "pyKVFinder",
        "mini_ftmap",
        "indeep",
        "pesto",
        "masif",
        "msmd",
        "pocketminer",
        "passer",
    }
    assert expected.issubset(set(report["tools"]))
    assert report["tools"]["indeep"]["required_level"] == "optional_ai"
    assert report["tools"]["mini_ftmap"]["required_level"] == "core_custom"


def test_external_path_only_tool_reports_available_when_path_exists(tmp_path, monkeypatch):
    external_tool = tmp_path / "InDeep"
    external_tool.mkdir()
    registry = tmp_path / "tool_registry.yaml"
    registry.write_text(
        "\n".join(
            [
                "tools:",
                "  indeep:",
                "    required_level: optional_ai",
                "    env: indeep",
                f"    external_path: {external_tool}",
                "    license_note: AGPLv3; do not vendor",
                "",
            ]
        ),
        encoding="utf-8",
    )

    monkeypatch.setattr(tool_preflight.importlib, "import_module", lambda name: (_ for _ in ()).throw(ImportError("missing")))
    monkeypatch.setattr(tool_preflight.shutil, "which", lambda command: None)

    ctx = RunContext.create(unique_run_id(), "smoke_env", repo_root=REPO_ROOT)
    initialize_logs(ctx)
    report = tool_preflight.run_tool_preflight(
        ctx,
        mode="discover",
        profile="codex_dev",
        registry_path=registry,
    )

    tool = report["tools"]["indeep"]
    assert tool["status"] == "available"
    assert tool["path"] == str(external_tool)
    assert tool["details"]["external_path"]["present"] is True


def test_binary_field_is_used_as_command_probe(tmp_path, monkeypatch):
    binary = tmp_path / ("fpocket.bat" if os.name == "nt" else "fpocket")
    if os.name == "nt":
        binary.write_text("@echo off\necho fpocket mock\n", encoding="utf-8")
    else:
        binary.write_text("#!/bin/sh\necho fpocket mock\n", encoding="utf-8")
    binary.chmod(0o755)
    registry = tmp_path / "tool_registry.yaml"
    registry.write_text(
        "\n".join(
            [
                "tools:",
                "  fpocket:",
                "    required_level: core",
                "    env: pyrosetta",
                f"    binary: {binary}",
                "    smoke_args: ['--help']",
                "",
            ]
        ),
        encoding="utf-8",
    )

    monkeypatch.setattr(tool_preflight.shutil, "which", lambda command: None)

    ctx = RunContext.create(unique_run_id(), "smoke_env", repo_root=REPO_ROOT)
    initialize_logs(ctx)
    report = tool_preflight.run_tool_preflight(
        ctx,
        mode="smoke",
        profile="hpc_strict",
        registry_path=registry,
    )

    tool = report["tools"]["fpocket"]
    assert report["status"] == "PASS"
    assert tool["status"] == "installed_and_smoke_passed"
    assert tool["path"] == str(binary)
