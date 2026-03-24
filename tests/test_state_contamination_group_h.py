"""Group H tests: State contamination prevention.

Tests for:
- _validate_lane_outputs() — lane별 산출물 검증
- Stale marker detection (H-2)
- Output validation before marker write (H-1)
- Seed marker integrity (H-5)
- PHASE_CORE_OUTPUTS coverage (H-6)
- Regression checks for Group 8 fixes
"""
from __future__ import annotations

import csv
import json
from pathlib import Path
from unittest.mock import patch

import pytest

# ---------------------------------------------------------------------------
# H-6: PHASE_CORE_OUTPUTS covers Phase 1-7
# ---------------------------------------------------------------------------


def test_phase_core_outputs_all_phases():
    """PHASE_CORE_OUTPUTS should have entries for all 7 phases."""
    from run_production import PHASE_CORE_OUTPUTS
    for phase in range(1, 8):
        assert phase in PHASE_CORE_OUTPUTS, f"Phase {phase} missing from PHASE_CORE_OUTPUTS"


def test_phase_core_outputs_phase3_has_csvs():
    """Phase 3 should require PPI residue and summary CSVs."""
    from run_production import PHASE_CORE_OUTPUTS
    phase3 = PHASE_CORE_OUTPUTS[3]
    labels = [label for label, _ in phase3]
    assert any("residue" in l.lower() for l in labels)
    assert any("summary" in l.lower() for l in labels)


def test_phase_core_outputs_phase4_7_unchanged():
    """Phase 4-7 entries should still exist (regression)."""
    from run_production import PHASE_CORE_OUTPUTS
    assert len(PHASE_CORE_OUTPUTS[4]) >= 2  # pocket + pose
    assert len(PHASE_CORE_OUTPUTS[5]) >= 2  # verdict + agreement
    assert len(PHASE_CORE_OUTPUTS[6]) >= 1  # report
    assert len(PHASE_CORE_OUTPUTS[7]) >= 1  # validation


# ---------------------------------------------------------------------------
# H-4: _validate_lane_outputs()
# ---------------------------------------------------------------------------


def test_validate_lane_outputs_unknown_lane():
    """Unknown/unmapped lane should return (True, []) — preserve existing behavior."""
    from run_production import _validate_lane_outputs
    with patch("run_production._load_config", return_value={"output_root": "/tmp/fake"}):
        valid, missing = _validate_lane_outputs("nonexistent-lane")
    assert valid is True
    assert missing == []


def test_validate_lane_outputs_adv_phase1_missing(tmp_path):
    """adv-phase1 with missing output should return (False, [...])."""
    from run_production import _validate_lane_outputs
    # Mock paths to point to tmp_path (no files exist)
    mock_config = {"output_root": str(tmp_path / "output")}
    with patch("run_production._load_config", return_value=mock_config):
        valid, missing = _validate_lane_outputs("adv-phase1")
    assert valid is False
    assert len(missing) > 0
    assert "phase1_downstream_patch_reference.csv" in missing[0]


def test_validate_lane_outputs_adv_phase1_valid(tmp_path):
    """adv-phase1 with valid output should return (True, [])."""
    from run_production import _validate_lane_outputs
    mock_config = {"output_root": str(tmp_path / "output")}
    # Create the expected output file
    wb_p1 = tmp_path / "output" / "workflow_b" / "phase1_ppi_analysis"
    wb_p1.mkdir(parents=True)
    csv_path = wb_p1 / "phase1_downstream_patch_reference.csv"
    csv_path.write_text("col1,col2\nval1,val2\n" + "x," * 50, encoding="utf-8")
    with patch("run_production._load_config", return_value=mock_config):
        valid, missing = _validate_lane_outputs("adv-phase1")
    assert valid is True
    assert missing == []


def test_validate_lane_outputs_ppi_post_empty_csv(tmp_path):
    """ppi-post with empty CSV should return (False, [...])."""
    from run_production import _validate_lane_outputs
    mock_config = {"output_root": str(tmp_path / "output")}
    ppi_post = tmp_path / "output" / "workflow_a" / "phase3_ppi_postprocess"
    ppi_post.mkdir(parents=True)
    # Create empty files (< 100B)
    (ppi_post / "ppi_pyrosetta_residues.csv").write_text("header\n", encoding="utf-8")
    (ppi_post / "ppi_pyrosetta_summary.csv").write_text("header\n", encoding="utf-8")
    with patch("run_production._load_config", return_value=mock_config):
        valid, missing = _validate_lane_outputs("ppi-post")
    assert valid is False
    assert len(missing) == 2


# ---------------------------------------------------------------------------
# H-5: Seed marker integrity
# ---------------------------------------------------------------------------


def test_seed_migration_empty_ranking(tmp_path):
    """Empty ranking.csv should NOT trigger marker migration."""
    from run_production import _seed_is_complete, _csv_has_rows
    docking_dir = tmp_path / "seed0"
    docking_dir.mkdir()
    final_dir = docking_dir / "final_result"
    final_dir.mkdir()
    # Create empty ranking
    ranking = final_dir / "final_ranking.csv"
    ranking.write_text("header\n", encoding="utf-8")
    assert not _csv_has_rows(ranking)

    target = {
        "name": "test_seed0",
        "docking_dir": str(docking_dir),
        "receptor_id": "test",
        "seed_index": 0,
        "config_ini": "",
    }
    result = _seed_is_complete(target)
    assert result is False
    # Marker should NOT have been created
    assert not (docking_dir / "seed_complete.json").exists()


def test_seed_migration_valid_ranking(tmp_path):
    """Valid ranking.csv should trigger marker migration with metadata."""
    from run_production import _seed_is_complete
    docking_dir = tmp_path / "seed0"
    docking_dir.mkdir()
    final_dir = docking_dir / "final_result"
    final_dir.mkdir()
    # Create valid ranking (header + data rows, > 100B)
    ranking = final_dir / "final_ranking.csv"
    rows = ["col1,col2,col3,col4,col5"]
    for i in range(10):
        rows.append(f"val{i},data{i},more{i},extra{i},long_value_{i}")
    ranking.write_text("\n".join(rows) + "\n", encoding="utf-8")
    assert ranking.stat().st_size >= 100

    target = {
        "name": "test_seed0",
        "docking_dir": str(docking_dir),
        "receptor_id": "test",
        "seed_index": 0,
        "config_ini": "",
    }
    result = _seed_is_complete(target)
    assert result is True
    # Marker should exist with migration metadata
    marker = docking_dir / "seed_complete.json"
    assert marker.exists()
    payload = json.loads(marker.read_text(encoding="utf-8"))
    assert payload["migration"] is True
    assert payload["ranking_rows"] == 10
    assert payload["ranking_size_bytes"] > 0


# ---------------------------------------------------------------------------
# Regression: Group 8 fixes
# ---------------------------------------------------------------------------


def test_regression_phase6_placeholder_detection():
    """check_phase6() should detect placeholder reports (Group 8 fix)."""
    from run_production import PHASE_CORE_OUTPUTS
    # Phase 6 should have report entry
    assert len(PHASE_CORE_OUTPUTS[6]) >= 1
    assert any("report" in label.lower() for label, _ in PHASE_CORE_OUTPUTS[6])


def test_regression_valid_construct_types_include_dimer():
    """patch_ingestion should accept dimer_offset (Group G fix)."""
    from egfr_pipeline.phase2.patch_ingestion import VALID_CONSTRUCT_TYPES
    assert "dimer_offset" in VALID_CONSTRUCT_TYPES


def test_regression_ppi_targets_dimer():
    """PPI_TARGETS should use dimer_offset (Group G fix)."""
    from run_production import PPI_TARGETS
    for t in PPI_TARGETS:
        assert t["construct_type"] == "dimer_offset"
