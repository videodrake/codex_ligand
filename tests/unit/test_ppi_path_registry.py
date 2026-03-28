import pytest
from pathlib import Path

from egfr_pipeline.ppi.path_registry import (
    resolve_metadata_path,
    resolve_phase1_interface_report,
    resolve_ranking_path,
)

pytestmark = pytest.mark.unit


def test_resolve_ranking_path_prefers_canonical(tmp_path: Path) -> None:
    run_dir = tmp_path / "run"
    canonical = run_dir / "final_result" / "final_ranking.csv"
    legacy = run_dir / "final_ranking.csv"
    canonical.parent.mkdir(parents=True)
    canonical.write_text("rank\n", encoding="utf-8")
    legacy.write_text("rank\n", encoding="utf-8")

    result = resolve_ranking_path(run_dir)

    assert result.path == canonical
    assert result.source_type == "canonical"


def test_resolve_metadata_path_legacy_restored(tmp_path: Path) -> None:
    run_dir = tmp_path / "run"
    restored = run_dir / "restored" / "pyrosetta_run_metadata.json"
    restored.parent.mkdir(parents=True)
    restored.write_text("{}", encoding="utf-8")

    result = resolve_metadata_path(run_dir)

    assert result.path == restored
    assert result.source_type == "legacy"


def test_resolve_phase1_interface_report_legacy_flag(tmp_path: Path) -> None:
    legacy = tmp_path / "output" / "phase1_ppi" / "phase1_interface_report.md"
    legacy.parent.mkdir(parents=True)
    legacy.write_text("legacy", encoding="utf-8")

    blocked = resolve_phase1_interface_report(tmp_path, allow_legacy_paths=False)
    assert blocked.path is None
    assert blocked.source_type == "missing"

    allowed = resolve_phase1_interface_report(tmp_path, allow_legacy_paths=True)
    assert allowed.path == legacy
    assert allowed.source_type == "legacy"
