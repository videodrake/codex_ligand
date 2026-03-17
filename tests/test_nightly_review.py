import json
import subprocess
from pathlib import Path

from scripts.nightly_review import (
    build_unit_entries,
    collect_git_summary,
    create_review_bundle,
    list_repo_files,
    review_units,
)


def _write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def _make_minimal_repo(root: Path) -> None:
    _write(
        root / "config" / "example-project.yaml",
        "\n".join(
            [
                "project_name: demo_project",
                "output_root: ./output",
                "mode: vina",
                "max_workers: 4",
                "receptors:",
                "  - id: 3GT8_raw",
                "  - id: EGFR_160-185",
                "ligands:",
                "  - id: 173940",
            ]
        ),
    )
    _write(root / "main.py", "def main():\n    return 0\n")
    _write(root / "run_production.py", "def main():\n    return 0\n")
    _write(root / "run_full_test.py", "def main():\n    return 0\n")
    _write(root / "README.md", "# Demo\n")
    _write(root / "CLAUDE.md", "# Notes\n")
    _write(root / "pytest.ini", "[pytest]\n")
    _write(root / "requirements-test.txt", "pytest\npyyaml\n")
    _write(root / "output" / "README.md", "# Output\n")
    _write(root / "docs" / "runbook.md", "# Runbook\n")
    _write(root / "docs" / "manual_execution.md", "# Manual\n")
    _write(root / "docs" / "architecture.md", "# Architecture\n")
    _write(root / "docs" / "output_artifact_map.md", "# Output Map\n")
    _write(root / "docs" / "data_inventory.md", "# Data Inventory\n")
    _write(root / "docs" / "test_suite_triage.md", "# Tests\n")
    _write(root / "config" / "README.md", "# Config\n")
    _write(root / "egfr_pipeline" / "config.py", "VALUE = 1\n")
    _write(root / "egfr_pipeline" / "runtime.py", "VALUE = 1\n")
    _write(root / "egfr_pipeline" / "residue_utils.py", "VALUE = 1\n")
    _write(root / "egfr_pipeline" / "step_view.py", "VALUE = 1\n")
    _write(root / "egfr_pipeline" / "vina" / "vina_executor.py", "VALUE = 1\n")
    _write(root / "egfr_pipeline" / "vina" / "parse_poses.py", "VALUE = 1\n")
    _write(root / "egfr_pipeline" / "verdict.py", "VALUE = 1\n")
    _write(root / "egfr_pipeline" / "report.py", "VALUE = 1\n")
    _write(root / "egfr_pipeline" / "validate.py", "VALUE = 1\n")
    _write(root / "scripts" / "precheck_guard.py", "VALUE = 1\n")
    _write(root / "scripts" / "rebuild_step_views.py", "VALUE = 1\n")
    _write(root / "tests" / "test_smoke_cli.py", "def test_ok():\n    assert True\n")
    _write(root / "tests" / "test_run_production.py", "def test_ok():\n    assert True\n")
    _write(root / "tests" / "test_validation_smoke.py", "def test_ok():\n    assert True\n")


def test_create_review_bundle_writes_expected_files(tmp_path: Path):
    _make_minimal_repo(tmp_path)

    result = create_review_bundle(
        repo_root=tmp_path,
        output_root=tmp_path / "output" / "review_automation",
        config_path=tmp_path / "config" / "example-project.yaml",
        label="nightly-test",
    )

    assert result["bundle_dir"].name == "nightly-test"
    assert result["manifest_path"].exists()
    assert result["checklist_path"].exists()
    assert result["prompt_path"].exists()
    assert result["report_path"].exists()
    assert result["latest_path"].exists()

    manifest = json.loads(result["manifest_path"].read_text(encoding="utf-8"))
    assert manifest["config_summary"]["loaded"] is True
    assert manifest["config_summary"]["project_name"] == "demo_project"
    assert manifest["inventory_summary"]["review_unit_count"] == len(review_units())

    control_plane = next(
        item for item in manifest["review_units"] if item["unit_id"] == "control_plane"
    )
    assert control_plane["status"] == "present"
    assert "main.py" in control_plane["matched_files"]
    assert "tests/test_smoke_cli.py" in control_plane["related_tests"]


def test_collect_git_summary_detects_dirty_files(tmp_path: Path):
    _make_minimal_repo(tmp_path)
    subprocess.run(("git", "init"), cwd=tmp_path, check=True, capture_output=True)
    subprocess.run(
        ("git", "config", "user.email", "codex@example.com"),
        cwd=tmp_path,
        check=True,
        capture_output=True,
    )
    subprocess.run(
        ("git", "config", "user.name", "Codex"),
        cwd=tmp_path,
        check=True,
        capture_output=True,
    )
    subprocess.run(("git", "add", "."), cwd=tmp_path, check=True, capture_output=True)
    subprocess.run(
        ("git", "commit", "-m", "init"),
        cwd=tmp_path,
        check=True,
        capture_output=True,
    )

    (tmp_path / "main.py").write_text("def main():\n    return 1\n", encoding="utf-8")
    _write(tmp_path / "notes.txt", "new\n")

    summary = collect_git_summary(tmp_path)

    assert summary["available"] is True
    assert "main.py" in summary["changed_files"]
    assert "notes.txt" in summary["changed_files"]


def test_build_unit_entries_prioritizes_changed_critical_units(tmp_path: Path):
    _make_minimal_repo(tmp_path)
    all_files = list_repo_files(tmp_path)

    entries = build_unit_entries(
        repo_root=tmp_path,
        all_files=all_files,
        changed_files=["main.py"],
    )

    assert entries[0]["unit_id"] == "control_plane"
    assert "main.py" in entries[0]["changed_files"]
