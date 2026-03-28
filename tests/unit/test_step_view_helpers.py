import pytest
from pathlib import Path

from egfr_pipeline.step_view import _ppi_target_ranking_path, _write_step_view_migration_notice

pytestmark = pytest.mark.unit


def test_ppi_target_ranking_path_returns_legacy_candidate_when_missing(tmp_path: Path) -> None:
    run_dir = tmp_path / "run"
    run_dir.mkdir(parents=True)

    assert _ppi_target_ranking_path(run_dir) == run_dir / "final_ranking.csv"


def test_write_step_view_migration_notice_is_idempotent(tmp_path: Path) -> None:
    legacy = tmp_path / "legacy_project"
    current = tmp_path / "output" / "workflow_a" / "step_views" / "project"
    legacy.mkdir(parents=True)
    current.mkdir(parents=True)

    _write_step_view_migration_notice(legacy_project_root=legacy, current_project_root=current)
    marker = legacy / "MOVED_TO_WORKFLOW_A_STEP_VIEWS.txt"
    first = marker.read_text(encoding="utf-8")

    _write_step_view_migration_notice(legacy_project_root=legacy, current_project_root=tmp_path / "ignored")
    second = marker.read_text(encoding="utf-8")

    assert first == second
    assert "workflow_a" in first
