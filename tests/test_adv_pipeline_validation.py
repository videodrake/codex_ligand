"""Tests for Workflow B (Advanced PPI-First Pipeline) validation guards."""
from __future__ import annotations

from pathlib import Path

import pytest

import run_production
from egfr_pipeline import paths as paths_mod


def _patch_output_root(tmp_path: Path, monkeypatch):
    """Patch paths module to use tmp_path as output root."""
    monkeypatch.setattr(paths_mod, "REPO_ROOT", tmp_path)
    monkeypatch.setattr(run_production, "REPO_ROOT", tmp_path)

    def _fake_output_root(config=None):
        return tmp_path / "output"

    monkeypatch.setattr(paths_mod, "output_root", _fake_output_root)


# ---------------------------------------------------------------------------
# _validate_adv_handoff tests
# ---------------------------------------------------------------------------

def test_validate_adv_handoff_phase2_raises_without_phase1_handoff(
    tmp_path: Path, monkeypatch,
) -> None:
    """Phase 2 requires phase1_downstream_patch_reference.csv from Phase 1."""
    _patch_output_root(tmp_path, monkeypatch)
    wb_p1 = tmp_path / "output" / "workflow_b" / "phase1_ppi_analysis"
    wb_p1.mkdir(parents=True)

    with pytest.raises(FileNotFoundError, match="Phase 1 핸드오프 없음"):
        run_production._validate_adv_handoff("adv-phase2")


def test_validate_adv_handoff_phase2_passes_with_handoff(
    tmp_path: Path, monkeypatch,
) -> None:
    """Phase 2 passes when handoff file exists."""
    _patch_output_root(tmp_path, monkeypatch)
    wb_p1 = tmp_path / "output" / "workflow_b" / "phase1_ppi_analysis"
    wb_p1.mkdir(parents=True)
    (wb_p1 / "phase1_downstream_patch_reference.csv").write_text(
        "residue_id,is_hotspot_any_state\nALA700,True\nALA701,True\n",
        encoding="utf-8",
    )

    # Should not raise
    run_production._validate_adv_handoff("adv-phase2")


def test_validate_adv_handoff_phase1_checks_ppi_completion(
    tmp_path: Path, monkeypatch,
) -> None:
    """Phase 1 (adv-phase1) fails when PPI docking is incomplete."""
    _patch_output_root(tmp_path, monkeypatch)
    monkeypatch.setattr(
        run_production, "check_phase2", lambda: ["3GT8_raw_seed0"],
    )

    with pytest.raises(FileNotFoundError, match="PPI 도킹 미완료"):
        run_production._validate_adv_handoff("adv-phase1")


def test_validate_adv_handoff_phase3_post_requires_round_log(
    tmp_path: Path, monkeypatch,
) -> None:
    """Phase 3 post mode requires a non-empty round log."""
    _patch_output_root(tmp_path, monkeypatch)
    wb_p3 = tmp_path / "output" / "workflow_b" / "phase3_focused_docking"
    wb_p3.mkdir(parents=True)

    with pytest.raises(FileNotFoundError, match="round log"):
        run_production._validate_adv_handoff("adv-phase3-post")


# ---------------------------------------------------------------------------
# Phase 3 cascade pre-checks
# ---------------------------------------------------------------------------

def test_phase3_execute_raises_without_setup(tmp_path: Path, monkeypatch) -> None:
    """Phase 3 cascade execute mode raises when job table is missing."""
    from egfr_pipeline.phase3 import rerun_cascade

    def _fake_output_root(config=None):
        return tmp_path

    monkeypatch.setattr(paths_mod, "output_root", _fake_output_root)

    wb_p2 = tmp_path / "workflow_b" / "phase2_pocket_analysis"
    wb_p2.mkdir(parents=True, exist_ok=True)
    (wb_p2 / "phase3_candidate_pocket_reference.csv").write_text(
        "col\nval\n", encoding="utf-8",
    )
    wb_p3 = tmp_path / "workflow_b" / "phase3_focused_docking"
    wb_p3.mkdir(parents=True, exist_ok=True)

    with pytest.raises(FileNotFoundError, match="job table"):
        rerun_cascade.run_phase3_cascade(mode="execute", round_id=0)


# ---------------------------------------------------------------------------
# Vina availability guard
# ---------------------------------------------------------------------------

def test_execute_round_raises_when_vina_unavailable(monkeypatch) -> None:
    """execute_round raises RuntimeError when Vina import fails."""
    from egfr_pipeline.phase3 import run_diverse_docking

    jobs = [{"seed_index": "0", "cpu_per_job": "4", "job_id": "test_001"}]

    import builtins
    original_import = builtins.__import__

    def mock_import(name, *args, **kwargs):
        if name == "egfr_pipeline.vina.vina_executor":
            raise ImportError("No vina")
        return original_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", mock_import)

    with pytest.raises(RuntimeError, match="Vina 모듈"):
        run_diverse_docking.execute_round(jobs, round_id=0)
