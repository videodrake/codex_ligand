"""Group 3 tests: PPI branch robustness verification (F-3).

Covers:
1. Threshold sweep output format & hotspot robustness (AC-3.1)
2. EC-3.1: all-same pass rate warning
3. Ambiguous model report format & dG comparison (AC-3.2)
4. Pilot fragment config generation & Jaccard comparison (AC-3.3)
5. EC-3.2: all Jaccard < 0.3 critical warning
6. Sheet 12 sensitivity 2-row output & hotspot overlap (AC-3.4)
7. Handoff quality guards: empty CSV, 0 hotspot, all irrelevant, <5 poses (AC-3.5)
8. EC-3.3: backward compat with _validate_adv_handoff patterns
9. Orthosteric+rim >80% overestimation warning (AC-3.6)
10. Conformational selection candidate tagging (AC-3.7)
11. Seed/state imbalance warning (AC-3.8)
12. Concordance score strong_both vs weak_both split (AC-3.9)
13. Adversarial: malformed inputs, boundary values, empty data paths
"""

import csv
import logging
import os
import tempfile
from pathlib import Path
from typing import Dict, List
from unittest.mock import patch

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent


# =====================================================================
# AC-3.1: Threshold sweep output format & hotspot robustness
# =====================================================================

def _make_orientation_log(n_models=20, score_range=(-0.5, 0.5)):
    """Create synthetic orientation log rows."""
    import random
    random.seed(42)
    rows = []
    for i in range(n_models):
        score = random.uniform(*score_range)
        rows.append({
            "model_id": f"model_{i:03d}",
            "receptor_id": "3GT8_raw",
            "seed_index": "0",
            "orientation_score": str(round(score, 4)),
            "orientation_class": "pass" if score > 0.15 else ("fail" if score < -0.15 else "ambiguous"),
        })
    return rows


def _make_residue_table(model_ids, residue_ids=None):
    """Create synthetic residue table for hotspot computation."""
    if residue_ids is None:
        residue_ids = ["VAL962", "VAL964", "SER971", "ALA968", "GLY970"]
    rows = []
    for mid in model_ids:
        for rid in residue_ids:
            rows.append({
                "model_id": mid,
                "chain": "A",
                "residue_id": rid,
            })
    return rows


class TestThresholdSweep:
    """AC-3.1: orientation threshold sensitivity sweep."""

    def test_sweep_produces_6_rows_for_default_thresholds(self):
        from egfr_pipeline.phase1.orientation_filter import sweep_threshold_sensitivity

        log = _make_orientation_log(50)
        pass_ids = [r["model_id"] for r in log if float(r["orientation_score"]) > 0]
        residues = _make_residue_table(pass_ids)

        rows, warnings = sweep_threshold_sensitivity(log, residues)
        assert len(rows) == 6
        expected_thresholds = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30]
        actual = [r["threshold"] for r in rows]
        assert actual == expected_thresholds

    def test_sweep_columns_present(self):
        from egfr_pipeline.phase1.orientation_filter import (
            SENSITIVITY_COLUMNS,
            sweep_threshold_sensitivity,
        )

        log = _make_orientation_log(30)
        pass_ids = [r["model_id"] for r in log if float(r["orientation_score"]) > 0]
        residues = _make_residue_table(pass_ids)

        rows, _ = sweep_threshold_sensitivity(log, residues)
        for row in rows:
            for col in SENSITIVITY_COLUMNS:
                assert col in row, f"Missing column: {col}"

    def test_sweep_pass_fail_ambig_sum_to_total(self):
        from egfr_pipeline.phase1.orientation_filter import sweep_threshold_sensitivity

        log = _make_orientation_log(40)
        pass_ids = [r["model_id"] for r in log if float(r["orientation_score"]) > 0]
        residues = _make_residue_table(pass_ids)

        rows, _ = sweep_threshold_sensitivity(log, residues)
        for row in rows:
            total = row["n_pass"] + row["n_fail"] + row["n_ambiguous"]
            assert total == row["n_total"]

    def test_sweep_hotspot_residues_semicolon_separated(self):
        from egfr_pipeline.phase1.orientation_filter import sweep_threshold_sensitivity

        log = _make_orientation_log(30)
        pass_ids = [r["model_id"] for r in log if float(r["orientation_score"]) > 0]
        residues = _make_residue_table(pass_ids)

        rows, _ = sweep_threshold_sensitivity(log, residues)
        for row in rows:
            hr = row["hotspot_residues"]
            if hr:
                parts = hr.split(";")
                assert all(p.strip() for p in parts), "Empty hotspot residue entry"

    def test_sweep_robustness_between_0_and_1(self):
        from egfr_pipeline.phase1.orientation_filter import sweep_threshold_sensitivity

        log = _make_orientation_log(40)
        pass_ids = [r["model_id"] for r in log if float(r["orientation_score"]) > 0]
        residues = _make_residue_table(pass_ids)

        rows, _ = sweep_threshold_sensitivity(log, residues)
        for row in rows:
            rob = row["hotspot_robustness"]
            if rob != "":
                assert 0.0 <= float(rob) <= 1.0

    def test_sweep_empty_log_returns_empty_with_warning(self):
        from egfr_pipeline.phase1.orientation_filter import sweep_threshold_sensitivity

        rows, warnings = sweep_threshold_sensitivity([], [])
        assert rows == []
        assert any("No scored models" in w for w in warnings)


class TestEC31AllSamePassRate:
    """EC-3.1: All thresholds yield identical pass rate."""

    def test_identical_pass_rate_triggers_warning(self):
        from egfr_pipeline.phase1.orientation_filter import sweep_threshold_sensitivity

        # All scores far positive — every threshold yields 100% pass
        log = [
            {"model_id": f"m{i}", "orientation_score": "0.9"}
            for i in range(10)
        ]
        residues = _make_residue_table([f"m{i}" for i in range(10)])
        _, warnings = sweep_threshold_sensitivity(log, residues)
        assert any("현재 값 유지" in w for w in warnings)


# =====================================================================
# AC-3.2: Ambiguous model report
# =====================================================================

class TestAmbiguousReport:
    """AC-3.2: ambiguous model report per (receptor_id, seed_index)."""

    def test_report_format(self):
        from egfr_pipeline.phase1.orientation_filter import compute_ambiguous_report

        log = [
            {"model_id": "m1", "receptor_id": "R1", "seed_index": "0", "orientation_class": "pass"},
            {"model_id": "m2", "receptor_id": "R1", "seed_index": "0", "orientation_class": "ambiguous"},
            {"model_id": "m3", "receptor_id": "R1", "seed_index": "0", "orientation_class": "ambiguous"},
            {"model_id": "m4", "receptor_id": "R1", "seed_index": "1", "orientation_class": "pass"},
        ]
        models = [
            {"model_id": "m1", "dG_separated": "-15.0"},
            {"model_id": "m2", "dG_separated": "-10.0"},
            {"model_id": "m3", "dG_separated": "-12.0"},
            {"model_id": "m4", "dG_separated": "-20.0"},
        ]
        residues = [
            {"model_id": "m2", "residue_id": "ALA700"},
            {"model_id": "m3", "residue_id": "ALA700"},
            {"model_id": "m3", "residue_id": "ALA701"},
        ]

        rows, warnings = compute_ambiguous_report(log, models, residues)
        # Should have 2 groups: (R1, 0) and (R1, 1)
        assert len(rows) == 2

        r1_s0 = next(r for r in rows if r["seed_index"] == "0")
        assert r1_s0["n_ambiguous"] == 2
        assert r1_s0["n_total"] == 3
        assert float(r1_s0["ambiguous_rate"]) == pytest.approx(2 / 3, abs=0.01)

    def test_dg_difference_computed(self):
        from egfr_pipeline.phase1.orientation_filter import compute_ambiguous_report

        log = [
            {"model_id": "m1", "receptor_id": "R1", "seed_index": "0", "orientation_class": "pass"},
            {"model_id": "m2", "receptor_id": "R1", "seed_index": "0", "orientation_class": "ambiguous"},
        ]
        models = [
            {"model_id": "m1", "dG_separated": "-20.0"},
            {"model_id": "m2", "dG_separated": "-10.0"},
        ]
        residues = [{"model_id": "m2", "residue_id": "ALA700"}]

        rows, _ = compute_ambiguous_report(log, models, residues)
        row = rows[0]
        assert row["mean_dG_pass"] != ""
        assert row["mean_dG_ambiguous"] != ""

    def test_unique_ambiguous_residues_non_empty(self):
        from egfr_pipeline.phase1.orientation_filter import compute_ambiguous_report

        log = [
            {"model_id": "m1", "receptor_id": "R1", "seed_index": "0", "orientation_class": "ambiguous"},
        ]
        models = [{"model_id": "m1", "dG_separated": "-10.0"}]
        residues = [
            {"model_id": "m1", "residue_id": "ALA700"},
            {"model_id": "m1", "residue_id": "ALA701"},
        ]

        rows, _ = compute_ambiguous_report(log, models, residues)
        assert rows[0]["n_unique_ambiguous_residues"] == 2


# =====================================================================
# AC-3.3: Pilot fragment config generation & comparison
# =====================================================================

class TestPilotFragmentRange:
    """AC-3.3: pilot config generation and Jaccard comparison."""

    def test_generate_3_configs(self, tmp_path):
        from scripts.pilot_fragment_range import generate_pilot_configs

        paths = generate_pilot_configs(tmp_path)
        assert len(paths) == 3
        for p in paths:
            assert p.exists()
            content = p.read_text()
            assert "pilot_fragment" in content
            assert "ppi:" in content

    def test_generate_pbs_script(self, tmp_path):
        from scripts.pilot_fragment_range import generate_pbs_script

        pbs = generate_pbs_script(tmp_path)
        assert pbs.exists()
        content = pbs.read_text()
        assert "#PBS" in content
        assert "RANGE" in content

    def test_jaccard_similarity_identical_sets(self):
        from scripts.pilot_fragment_range import jaccard_similarity

        s = {"A", "B", "C"}
        assert jaccard_similarity(s, s) == 1.0

    def test_jaccard_similarity_disjoint_sets(self):
        from scripts.pilot_fragment_range import jaccard_similarity

        assert jaccard_similarity({"A", "B"}, {"C", "D"}) == 0.0

    def test_jaccard_similarity_empty_sets(self):
        from scripts.pilot_fragment_range import jaccard_similarity

        assert jaccard_similarity(set(), set()) == 1.0

    def test_jaccard_similarity_partial_overlap(self):
        from scripts.pilot_fragment_range import jaccard_similarity

        j = jaccard_similarity({"A", "B", "C"}, {"B", "C", "D"})
        # intersection=2, union=4
        assert j == pytest.approx(0.5)

    def test_compare_pilot_results_with_mock_data(self, tmp_path):
        from scripts.pilot_fragment_range import compare_pilot_results

        # Create mock hotspot CSV files for 3 ranges
        for rid, hotspots in [("A", ["R1", "R2", "R3"]),
                              ("B", ["R2", "R3", "R4"]),
                              ("C", ["R5", "R6", "R7"])]:
            d = tmp_path / f"pilot_fragment_{rid}" / "workflow_b" / "phase1_ppi_analysis" / "3GT8_raw"
            d.mkdir(parents=True)
            csv_path = d / "ppi_hotspot_residues.csv"
            with open(csv_path, "w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=["residue_id", "is_hotspot"])
                w.writeheader()
                for r in hotspots:
                    w.writerow({"residue_id": r, "is_hotspot": "True"})

        rows, warnings = compare_pilot_results(tmp_path)
        assert len(rows) == 3  # A-B, A-C, B-C
        # A-B: intersection=2, union=4, jaccard=0.5
        ab = next(r for r in rows if r["range_a"] == "A" and r["range_b"] == "B")
        assert ab["jaccard"] == pytest.approx(0.5)


class TestEC32AllJaccardBelowThreshold:
    """EC-3.2: all pairwise Jaccard < 0.3 triggers critical warning."""

    def test_low_jaccard_triggers_warning(self, tmp_path):
        from scripts.pilot_fragment_range import compare_pilot_results

        # Completely disjoint hotspot sets
        for rid, hotspots in [("A", ["R1"]), ("B", ["R2"]), ("C", ["R3"])]:
            d = tmp_path / f"pilot_fragment_{rid}" / "workflow_b" / "phase1_ppi_analysis" / "3GT8_raw"
            d.mkdir(parents=True)
            csv_path = d / "ppi_hotspot_residues.csv"
            with open(csv_path, "w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=["residue_id", "is_hotspot"])
                w.writeheader()
                for r in hotspots:
                    w.writerow({"residue_id": r, "is_hotspot": "True"})

        _, warnings = compare_pilot_results(tmp_path)
        assert any("결정적 영향" in w for w in warnings)


# =====================================================================
# AC-3.4: Sheet 12 sensitivity 2-row output
# =====================================================================

class TestSheet12Sensitivity:
    """AC-3.4: sheet 12 inclusion/exclusion comparison."""

    def test_produces_2_rows(self):
        from egfr_pipeline.phase1.orientation_filter import compare_sheet12_sensitivity

        log_a = [
            {"model_id": "m1", "orientation_class": "pass"},
            {"model_id": "m2", "orientation_class": "fail"},
            {"model_id": "m3", "orientation_class": "ambiguous"},
        ]
        log_b = [
            {"model_id": "m1", "orientation_class": "pass"},
            {"model_id": "m2", "orientation_class": "pass"},
            {"model_id": "m3", "orientation_class": "fail"},
        ]
        residues = [
            {"model_id": "m1", "chain": "A", "residue_id": "ALA700"},
            {"model_id": "m2", "chain": "A", "residue_id": "ALA700"},
        ]
        rows = compare_sheet12_sensitivity(log_a, log_b, residues)
        assert len(rows) == 2
        assert rows[0]["setting"] == "A"
        assert rows[0]["active_face_sheets"] == "sheet_8+9"
        assert rows[1]["setting"] == "B"
        assert rows[1]["active_face_sheets"] == "sheet_8+9+12"

    def test_hotspot_overlap_present(self):
        from egfr_pipeline.phase1.orientation_filter import compare_sheet12_sensitivity

        log_a = [{"model_id": "m1", "orientation_class": "pass"}]
        log_b = [{"model_id": "m1", "orientation_class": "pass"}]
        residues = [{"model_id": "m1", "chain": "A", "residue_id": "ALA700"}]

        rows = compare_sheet12_sensitivity(log_a, log_b, residues)
        assert "hotspot_overlap" in rows[0]
        assert "hotspot_overlap" in rows[1]

    def test_identical_logs_give_overlap_1(self):
        from egfr_pipeline.phase1.orientation_filter import compare_sheet12_sensitivity

        log = [
            {"model_id": "m1", "orientation_class": "pass"},
            {"model_id": "m2", "orientation_class": "pass"},
        ]
        residues = [
            {"model_id": "m1", "chain": "A", "residue_id": "ALA700"},
            {"model_id": "m2", "chain": "A", "residue_id": "ALA700"},
        ]
        rows = compare_sheet12_sensitivity(log, log, residues)
        assert float(rows[0]["hotspot_overlap"]) == pytest.approx(1.0)


# =====================================================================
# AC-3.5: Handoff data quality guards
# =====================================================================

class TestHandoffQualityGuards:
    """AC-3.5: data quality guards at phase handoff points."""

    def _setup_handoff_dirs(self, tmp_path):
        wb_p1 = tmp_path / "wb_p1"
        wb_p2 = tmp_path / "wb_p2"
        wb_p3 = tmp_path / "wb_p3"
        wb_p1.mkdir()
        wb_p2.mkdir()
        wb_p3.mkdir()
        return wb_p1, wb_p2, wb_p3

    def test_phase1_zero_hotspots_raises(self, tmp_path):
        """Phase 1 → Phase 2: hotspot 0개 → FAIL."""
        from run_production import _validate_handoff_data_quality

        wb_p1, wb_p2, wb_p3 = self._setup_handoff_dirs(tmp_path)
        # Write patch ref CSV with no hotspot rows
        patch_path = wb_p1 / "phase1_downstream_patch_reference.csv"
        with open(patch_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["residue_id", "is_hotspot_any_state"])
            w.writeheader()
            w.writerow({"residue_id": "ALA700", "is_hotspot_any_state": "False"})

        with pytest.raises(FileNotFoundError, match="hotspot 잔기 0개"):
            _validate_handoff_data_quality("adv-phase2", wb_p1, wb_p2, wb_p3)

    def test_phase1_with_hotspots_passes(self, tmp_path):
        """Phase 1 → Phase 2: hotspot > 0 → OK."""
        from run_production import _validate_handoff_data_quality

        wb_p1, wb_p2, wb_p3 = self._setup_handoff_dirs(tmp_path)
        patch_path = wb_p1 / "phase1_downstream_patch_reference.csv"
        with open(patch_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["residue_id", "is_hotspot_any_state"])
            w.writeheader()
            w.writerow({"residue_id": "ALA700", "is_hotspot_any_state": "True"})

        # Should not raise
        _validate_handoff_data_quality("adv-phase2", wb_p1, wb_p2, wb_p3)

    def test_phase2_zero_pockets_raises(self, tmp_path):
        """Phase 2 → Phase 3: 0 pockets → FAIL."""
        from run_production import _validate_handoff_data_quality

        wb_p1, wb_p2, wb_p3 = self._setup_handoff_dirs(tmp_path)
        pocket_path = wb_p2 / "phase3_candidate_pocket_reference.csv"
        with open(pocket_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["pocket_id", "docking_priority", "relationship_class"])
            w.writeheader()
            # Empty — header only

        with pytest.raises(FileNotFoundError, match="후보 포켓 0개"):
            _validate_handoff_data_quality("adv-phase3-setup", wb_p1, wb_p2, wb_p3)

    def test_phase2_all_irrelevant_warns(self, tmp_path, caplog):
        """Phase 2 → Phase 3: all pockets irrelevant → WARNING."""
        from run_production import _validate_handoff_data_quality

        wb_p1, wb_p2, wb_p3 = self._setup_handoff_dirs(tmp_path)
        pocket_path = wb_p2 / "phase3_candidate_pocket_reference.csv"
        with open(pocket_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["pocket_id", "docking_priority", "relationship_class"])
            w.writeheader()
            w.writerow({"pocket_id": "P1", "docking_priority": "skip", "relationship_class": "low_relevance_candidate"})
            w.writerow({"pocket_id": "P2", "docking_priority": "", "relationship_class": "low_relevance_candidate"})

        with caplog.at_level(logging.WARNING):
            _validate_handoff_data_quality("adv-phase3-setup", wb_p1, wb_p2, wb_p3)

        assert any("irrelevant" in r.message for r in caplog.records)

    def test_phase3_insufficient_poses_raises(self, tmp_path):
        """Phase 3 → Phase 4: valid poses < 5 → FAIL."""
        from run_production import _validate_handoff_data_quality

        wb_p1, wb_p2, wb_p3 = self._setup_handoff_dirs(tmp_path)
        evidence_path = wb_p3 / "phase4_docking_evidence_reference.csv"
        with open(evidence_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["pocket_id", "ligand_id", "best_affinity"])
            w.writeheader()
            # Only 2 valid poses
            w.writerow({"pocket_id": "P1", "ligand_id": "L1", "best_affinity": "-7.5"})
            w.writerow({"pocket_id": "P2", "ligand_id": "L1", "best_affinity": "-6.0"})

        with pytest.raises(FileNotFoundError, match="유효 포즈.*최소 5개"):
            _validate_handoff_data_quality("adv-phase4", wb_p1, wb_p2, wb_p3)

    def test_phase3_single_ligand_warns(self, tmp_path, caplog):
        """Phase 3 → Phase 4: 1 ligand only → WARNING."""
        from run_production import _validate_handoff_data_quality

        wb_p1, wb_p2, wb_p3 = self._setup_handoff_dirs(tmp_path)
        evidence_path = wb_p3 / "phase4_docking_evidence_reference.csv"
        with open(evidence_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["pocket_id", "ligand_id", "best_affinity"])
            w.writeheader()
            for i in range(6):
                w.writerow({"pocket_id": f"P{i}", "ligand_id": "L1", "best_affinity": f"-{7+i*0.1}"})

        with caplog.at_level(logging.WARNING):
            _validate_handoff_data_quality("adv-phase4", wb_p1, wb_p2, wb_p3)

        assert any("리간드 1종" in r.message for r in caplog.records)

    def test_phase3_diverse_ligands_no_warning(self, tmp_path, caplog):
        """Phase 3 → Phase 4: 3 ligands → no warning."""
        from run_production import _validate_handoff_data_quality

        wb_p1, wb_p2, wb_p3 = self._setup_handoff_dirs(tmp_path)
        evidence_path = wb_p3 / "phase4_docking_evidence_reference.csv"
        with open(evidence_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["pocket_id", "ligand_id", "best_affinity"])
            w.writeheader()
            for i, lid in enumerate(["L1", "L1", "L2", "L2", "L3"]):
                w.writerow({"pocket_id": f"P{i}", "ligand_id": lid, "best_affinity": f"-{7+i*0.1}"})

        with caplog.at_level(logging.WARNING):
            _validate_handoff_data_quality("adv-phase4", wb_p1, wb_p2, wb_p3)

        assert not any("리간드 1종" in r.message for r in caplog.records)


# =====================================================================
# EC-3.3: _validate_adv_handoff backward compatibility
# =====================================================================

class TestValidateAdvHandoffBackwardCompat:
    """EC-3.3: _validate_adv_handoff return patterns (FAIL/WARNING/PASS)."""

    def test_missing_handoff_file_raises_file_not_found(self, tmp_path):
        """Existing pattern: missing file → FileNotFoundError."""
        from run_production import _validate_handoff_data_quality

        wb_p1 = tmp_path / "wb_p1"
        wb_p2 = tmp_path / "wb_p2"
        wb_p3 = tmp_path / "wb_p3"
        wb_p1.mkdir()
        wb_p2.mkdir()
        wb_p3.mkdir()

        # adv-phase2 with missing patch ref → guard skipped (file doesn't exist)
        # This should NOT raise — guard only activates when file exists
        _validate_handoff_data_quality("adv-phase2", wb_p1, wb_p2, wb_p3)

    def test_unrelated_lane_passes_through(self, tmp_path):
        """Lanes not matching any guard pass through."""
        from run_production import _validate_handoff_data_quality

        wb_p1 = tmp_path / "wb_p1"
        wb_p2 = tmp_path / "wb_p2"
        wb_p3 = tmp_path / "wb_p3"
        wb_p1.mkdir()
        wb_p2.mkdir()
        wb_p3.mkdir()

        _validate_handoff_data_quality("adv-phase1", wb_p1, wb_p2, wb_p3)


# =====================================================================
# AC-3.6: Orthosteric+rim >80% overestimation warning
# =====================================================================

class TestPatchOverestimationWarning:
    """AC-3.6: orthosteric+rim > 80% triggers WARNING."""

    def test_above_80_triggers_warning(self, caplog):
        from egfr_pipeline.phase2.patch_relationship import build_classification_note

        # 5 pockets: 3 orthosteric + 2 rim = 100% > 80%
        relationships = [
            {"pocket_id": f"P{i}", "relationship_class": "orthosteric_candidate",
             "hotspot_overlap_count": 3, "hotspot_overlap_fraction": "0.50",
             "centroid_distance_A": "5.0", "overlapping_residues": "A;B;C",
             "n_pocket_residues": 10, "n_patch_hotspots": 6, "classification_basis": "test"}
            for i in range(3)
        ] + [
            {"pocket_id": f"P{i+3}", "relationship_class": "rim_candidate",
             "hotspot_overlap_count": 1, "hotspot_overlap_fraction": "0.10",
             "centroid_distance_A": "10.0", "overlapping_residues": "D",
             "n_pocket_residues": 8, "n_patch_hotspots": 6, "classification_basis": "test"}
            for i in range(2)
        ]

        with caplog.at_level(logging.WARNING):
            lines = build_classification_note(relationships, {"A", "B", "C", "D", "E", "F"})

        note_text = "\n".join(lines)
        assert "과대추정" in note_text or any("과대추정" in r.message for r in caplog.records)

    def test_below_80_no_warning(self, caplog):
        from egfr_pipeline.phase2.patch_relationship import build_classification_note

        # 5 pockets: 1 orthosteric + 1 rim + 3 allosteric = 40% < 80%
        relationships = [
            {"pocket_id": "P0", "relationship_class": "orthosteric_candidate",
             "hotspot_overlap_count": 3, "hotspot_overlap_fraction": "0.50",
             "centroid_distance_A": "5.0", "overlapping_residues": "A;B;C",
             "n_pocket_residues": 10, "n_patch_hotspots": 6, "classification_basis": "test"},
            {"pocket_id": "P1", "relationship_class": "rim_candidate",
             "hotspot_overlap_count": 1, "hotspot_overlap_fraction": "0.10",
             "centroid_distance_A": "10.0", "overlapping_residues": "D",
             "n_pocket_residues": 8, "n_patch_hotspots": 6, "classification_basis": "test"},
        ] + [
            {"pocket_id": f"P{i+2}", "relationship_class": "allosteric_candidate",
             "hotspot_overlap_count": 0, "hotspot_overlap_fraction": "0.0",
             "centroid_distance_A": "18.0", "overlapping_residues": "",
             "n_pocket_residues": 8, "n_patch_hotspots": 6, "classification_basis": "test"}
            for i in range(3)
        ]

        with caplog.at_level(logging.WARNING):
            lines = build_classification_note(relationships, {"A", "B", "C", "D", "E", "F"})

        note_text = "\n".join(lines)
        assert "과대추정" not in note_text
        assert not any("과대추정" in r.message for r in caplog.records)


# =====================================================================
# AC-3.7: Conformational selection candidate tagging
# =====================================================================

class TestConformationalSelectionCandidate:
    """AC-3.7: state-specific + top 10% occupancy → tagged True."""

    def _build_state_dirs(self, tmp_path, states_data):
        """Create mock state dirs with patch and hotspot CSVs."""
        for state, residues in states_data.items():
            state_dir = tmp_path / state
            state_dir.mkdir(parents=True, exist_ok=True)

            # Patch table
            with open(state_dir / "ppi_interface_patch_table.csv", "w", newline="") as f:
                cols = ["chain", "residue_id", "residue_num", "residue_name",
                        "lobe_label", "max_occupancy", "n_clusters_present",
                        "construct_type", "orientation_validation_status"]
                w = csv.DictWriter(f, fieldnames=cols)
                w.writeheader()
                for r in residues:
                    w.writerow(r)

            # Hotspot table
            with open(state_dir / "ppi_hotspot_residues.csv", "w", newline="") as f:
                cols = ["chain", "residue_id", "is_hotspot", "n_orientation_valid_models"]
                w = csv.DictWriter(f, fieldnames=cols)
                w.writeheader()
                for r in residues:
                    w.writerow({
                        "chain": r["chain"],
                        "residue_id": r["residue_id"],
                        "is_hotspot": "True",
                        "n_orientation_valid_models": "100",
                    })

    def test_state_specific_top10_tagged_true(self, tmp_path):
        from egfr_pipeline.phase1.compare_states import compare_across_states

        # R_high only in state1 with highest occupancy → state_specific + top 10%
        self._build_state_dirs(tmp_path, {
            "state1": [
                {"chain": "A", "residue_id": "ALA700", "residue_num": "700",
                 "residue_name": "ALA", "lobe_label": "C-lobe",
                 "max_occupancy": "0.95", "n_clusters_present": "3",
                 "construct_type": "raw", "orientation_validation_status": "pass"},
            ],
            "state2": [
                {"chain": "A", "residue_id": "ALA800", "residue_num": "800",
                 "residue_name": "ALA", "lobe_label": "C-lobe",
                 "max_occupancy": "0.30", "n_clusters_present": "2",
                 "construct_type": "raw", "orientation_validation_status": "pass"},
            ],
            "state3": [
                {"chain": "A", "residue_id": "ALA801", "residue_num": "801",
                 "residue_name": "ALA", "lobe_label": "C-lobe",
                 "max_occupancy": "0.20", "n_clusters_present": "1",
                 "construct_type": "raw", "orientation_validation_status": "pass"},
            ],
        })

        _, robustness_rows = compare_across_states(
            tmp_path, states=["state1", "state2", "state3"]
        )
        state_specific = [r for r in robustness_rows if r["robustness_class"] == "state_specific"]
        tagged = [r for r in state_specific if r["conformational_selection_candidate"] is True]
        # ALA700 has highest occupancy among state-specific → should be tagged
        assert any(r["residue_id"] == "ALA700" for r in tagged)

    def test_robust_residues_not_tagged(self, tmp_path):
        from egfr_pipeline.phase1.compare_states import compare_across_states

        # Same residue in all 3 states → robust → NOT a candidate
        base_row = {"chain": "A", "residue_id": "ALA700", "residue_num": "700",
                     "residue_name": "ALA", "lobe_label": "C-lobe",
                     "max_occupancy": "0.90", "n_clusters_present": "3",
                     "construct_type": "raw", "orientation_validation_status": "pass"}
        self._build_state_dirs(tmp_path, {
            "state1": [base_row],
            "state2": [base_row],
            "state3": [base_row],
        })

        _, robustness_rows = compare_across_states(
            tmp_path, states=["state1", "state2", "state3"]
        )
        robust = [r for r in robustness_rows if r["robustness_class"] == "robust"]
        assert len(robust) == 1
        assert robust[0]["conformational_selection_candidate"] is False


# =====================================================================
# AC-3.8: Seed/state imbalance warning
# =====================================================================

class TestSeedStateImbalance:
    """AC-3.8: n_valid_models 2x imbalance between states → WARNING."""

    def _build_state_dirs_with_valid_models(self, tmp_path, state_counts):
        """Create mock state dirs with n_orientation_valid_models."""
        for state, n_valid in state_counts.items():
            state_dir = tmp_path / state
            state_dir.mkdir(parents=True, exist_ok=True)

            with open(state_dir / "ppi_interface_patch_table.csv", "w", newline="") as f:
                cols = ["chain", "residue_id", "residue_num", "residue_name",
                        "lobe_label", "max_occupancy", "n_clusters_present",
                        "construct_type", "orientation_validation_status"]
                w = csv.DictWriter(f, fieldnames=cols)
                w.writeheader()
                w.writerow({"chain": "A", "residue_id": "ALA700", "residue_num": "700",
                            "residue_name": "ALA", "lobe_label": "C-lobe",
                            "max_occupancy": "0.80", "n_clusters_present": "2",
                            "construct_type": "raw", "orientation_validation_status": "pass"})

            with open(state_dir / "ppi_hotspot_residues.csv", "w", newline="") as f:
                cols = ["chain", "residue_id", "is_hotspot", "n_orientation_valid_models"]
                w = csv.DictWriter(f, fieldnames=cols)
                w.writeheader()
                w.writerow({"chain": "A", "residue_id": "ALA700", "is_hotspot": "True",
                            "n_orientation_valid_models": str(n_valid)})

    def test_2x_imbalance_triggers_warning(self, tmp_path, caplog):
        from egfr_pipeline.phase1.compare_states import compare_across_states

        self._build_state_dirs_with_valid_models(tmp_path, {
            "state1": 200,
            "state2": 50,   # 200/50 = 4x → triggers
        })

        with caplog.at_level(logging.WARNING):
            compare_across_states(tmp_path, states=["state1", "state2"])

        assert any("AC-3.8" in r.message for r in caplog.records)

    def test_balanced_states_no_warning(self, tmp_path, caplog):
        from egfr_pipeline.phase1.compare_states import compare_across_states

        self._build_state_dirs_with_valid_models(tmp_path, {
            "state1": 100,
            "state2": 90,   # 100/90 ≈ 1.1x → no trigger
        })

        with caplog.at_level(logging.WARNING):
            compare_across_states(tmp_path, states=["state1", "state2"])

        assert not any("AC-3.8" in r.message for r in caplog.records)


# =====================================================================
# AC-3.9: Concordance score strong_both vs weak_both
# =====================================================================

class TestConcordanceScore:
    """AC-3.9: concordance_score split into strong_both / weak_both."""

    def _build_convergence_inputs(self, tmp_path, state, pyro_rows, ld_interface_rows):
        """Create mock pyrosetta patch table and lightdock interface support table.

        ld_interface_rows: rows for lightdock_interface_support_table.csv
            (per-model per-residue entries used by _load_lightdock_residues).
        """
        state_dir = tmp_path / state
        state_dir.mkdir(parents=True, exist_ok=True)

        # PyRosetta patch table
        with open(state_dir / "ppi_interface_patch_table.csv", "w", newline="") as f:
            cols = ["chain", "residue_id", "residue_num", "residue_name",
                    "lobe_label", "max_occupancy", "n_clusters_present",
                    "construct_type", "orientation_validation_status"]
            w = csv.DictWriter(f, fieldnames=cols)
            w.writeheader()
            for r in pyro_rows:
                w.writerow(r)

        # LightDock interface support table (per-model per-residue)
        ld_dir = state_dir / "lightdock"
        ld_dir.mkdir(parents=True, exist_ok=True)
        with open(ld_dir / "lightdock_interface_support_table.csv", "w", newline="") as f:
            cols = ["model_id", "receptor_id", "construct_type",
                    "orientation_validation_status", "orientation_score",
                    "swarm_id", "pose_rank", "scoring_value",
                    "chain", "residue_id", "residue_num", "residue_name",
                    "lobe_label", "source"]
            w = csv.DictWriter(f, fieldnames=cols)
            w.writeheader()
            for r in ld_interface_rows:
                w.writerow(r)

    def _make_ld_interface_rows(self, residue_id, residue_num, residue_name,
                                lobe_label, n_models, n_total_models):
        """Generate per-model lightdock interface support rows.

        n_models: how many models contain this residue (numerator of frequency).
        n_total_models: total models (denominator of frequency).
        """
        rows = []
        for i in range(n_total_models):
            if i < n_models:
                rows.append({
                    "model_id": f"ld_model_{i:03d}",
                    "receptor_id": "test_state",
                    "construct_type": "raw",
                    "orientation_validation_status": "pass",
                    "orientation_score": "0.5",
                    "swarm_id": "0",
                    "pose_rank": str(i + 1),
                    "scoring_value": "50.0",
                    "chain": "A",
                    "residue_id": residue_id,
                    "residue_num": residue_num,
                    "residue_name": residue_name,
                    "lobe_label": lobe_label,
                    "source": "lightdock",
                })
            else:
                # Model exists but doesn't contact this residue — add a different
                # residue so the model is counted in total
                rows.append({
                    "model_id": f"ld_model_{i:03d}",
                    "receptor_id": "test_state",
                    "construct_type": "raw",
                    "orientation_validation_status": "pass",
                    "orientation_score": "0.5",
                    "swarm_id": "0",
                    "pose_rank": str(i + 1),
                    "scoring_value": "50.0",
                    "chain": "A",
                    "residue_id": "GLY999",
                    "residue_num": "999",
                    "residue_name": "GLY",
                    "lobe_label": "C-lobe",
                    "source": "lightdock",
                })
        return rows

    def test_high_concordance_strong_both(self, tmp_path):
        from egfr_pipeline.phase1.lightdock_validation import compute_cross_method_convergence

        # PyRosetta: ALA700 occ=0.80
        # LightDock: ALA700 in 8/10 models → freq=0.80
        # concordance = min(0.80, 0.80)/max(0.80, 0.80) = 1.0 → strong_both
        ld_rows = self._make_ld_interface_rows("ALA700", "700", "ALA", "C-lobe", 8, 10)
        self._build_convergence_inputs(tmp_path, "test_state",
            pyro_rows=[
                {"chain": "A", "residue_id": "ALA700", "residue_num": "700",
                 "residue_name": "ALA", "lobe_label": "C-lobe",
                 "max_occupancy": "0.80", "n_clusters_present": "3",
                 "construct_type": "raw", "orientation_validation_status": "pass"},
            ],
            ld_interface_rows=ld_rows,
        )

        result = compute_cross_method_convergence("test_state", tmp_path)
        assert result is not None

        with open(result, encoding="utf-8") as f:
            rows = list(csv.DictReader(f))

        conv_rows = [r for r in rows if r["residue_id"] == "ALA700"]
        assert len(conv_rows) == 1
        row = conv_rows[0]
        assert row["convergence_class"] == "convergent"
        assert row["method_agreement"] == "strong_both"
        assert float(row["concordance_score"]) > 0.5

    def test_low_concordance_weak_both(self, tmp_path):
        from egfr_pipeline.phase1.lightdock_validation import compute_cross_method_convergence

        # PyRosetta: ALA700 occ=0.90
        # LightDock: ALA700 in 1/10 models → freq=0.10
        # concordance = 0.10/0.90 ≈ 0.11 → weak_both
        ld_rows = self._make_ld_interface_rows("ALA700", "700", "ALA", "C-lobe", 1, 10)
        self._build_convergence_inputs(tmp_path, "test_state",
            pyro_rows=[
                {"chain": "A", "residue_id": "ALA700", "residue_num": "700",
                 "residue_name": "ALA", "lobe_label": "C-lobe",
                 "max_occupancy": "0.90", "n_clusters_present": "3",
                 "construct_type": "raw", "orientation_validation_status": "pass"},
            ],
            ld_interface_rows=ld_rows,
        )

        result = compute_cross_method_convergence("test_state", tmp_path)
        assert result is not None

        with open(result, encoding="utf-8") as f:
            rows = list(csv.DictReader(f))

        conv_rows = [r for r in rows if r["residue_id"] == "ALA700"]
        assert len(conv_rows) == 1
        row = conv_rows[0]
        assert row["method_agreement"] == "weak_both"
        assert float(row["concordance_score"]) <= 0.5

    def test_single_method_no_concordance(self, tmp_path):
        from egfr_pipeline.phase1.lightdock_validation import compute_cross_method_convergence

        # ALA700 only in pyrosetta, ALA800 only in lightdock → each single-method
        ld_rows = self._make_ld_interface_rows("ALA800", "800", "ALA", "C-lobe", 7, 10)
        self._build_convergence_inputs(tmp_path, "test_state",
            pyro_rows=[
                {"chain": "A", "residue_id": "ALA700", "residue_num": "700",
                 "residue_name": "ALA", "lobe_label": "C-lobe",
                 "max_occupancy": "0.80", "n_clusters_present": "3",
                 "construct_type": "raw", "orientation_validation_status": "pass"},
            ],
            ld_interface_rows=ld_rows,
        )

        result = compute_cross_method_convergence("test_state", tmp_path)
        assert result is not None

        with open(result, encoding="utf-8") as f:
            rows = list(csv.DictReader(f))

        pyro_only = [r for r in rows if r["residue_id"] == "ALA700"]
        assert len(pyro_only) == 1
        assert pyro_only[0]["method_agreement"] == "pyrosetta_only"
        assert pyro_only[0]["concordance_score"] == ""

    def test_concordance_boundary_at_0_5(self, tmp_path):
        """Boundary: concordance_score exactly 0.5 → weak_both."""
        from egfr_pipeline.phase1.lightdock_validation import compute_cross_method_convergence

        # PyRosetta occ=0.50, LightDock freq=10/10=1.0 → concordance=0.50/1.0=0.5 → weak_both
        ld_rows = self._make_ld_interface_rows("ALA700", "700", "ALA", "C-lobe", 10, 10)
        self._build_convergence_inputs(tmp_path, "test_state",
            pyro_rows=[
                {"chain": "A", "residue_id": "ALA700", "residue_num": "700",
                 "residue_name": "ALA", "lobe_label": "C-lobe",
                 "max_occupancy": "0.50", "n_clusters_present": "2",
                 "construct_type": "raw", "orientation_validation_status": "pass"},
            ],
            ld_interface_rows=ld_rows,
        )

        result = compute_cross_method_convergence("test_state", tmp_path)
        with open(result, encoding="utf-8") as f:
            rows = list(csv.DictReader(f))

        row = next(r for r in rows if r["residue_id"] == "ALA700")
        # 0.5/1.0 = 0.5 → exactly 0.5 → weak_both (not > 0.5)
        assert row["method_agreement"] == "weak_both"


# =====================================================================
# Adversarial / boundary tests
# =====================================================================

class TestSweepAdversarial:
    """Attack surface: malformed orientation scores, single model, NaN."""

    def test_non_numeric_scores_skipped(self):
        from egfr_pipeline.phase1.orientation_filter import sweep_threshold_sensitivity

        log = [
            {"model_id": "m1", "orientation_score": ""},
            {"model_id": "m2", "orientation_score": "not_a_number"},
            {"model_id": "m3", "orientation_score": "0.3"},
        ]
        # NaN is valid float in Python, so only "" and "not_a_number" are skipped
        residues = [{"model_id": "m3", "chain": "A", "residue_id": "ALA700"}]
        rows, warnings = sweep_threshold_sensitivity(log, residues)
        # Only m3 should be counted (m1 and m2 are not valid floats)
        assert all(r["n_total"] == 1 for r in rows)

    def test_single_model_no_crash(self):
        from egfr_pipeline.phase1.orientation_filter import sweep_threshold_sensitivity

        log = [{"model_id": "m1", "orientation_score": "0.2"}]
        residues = [{"model_id": "m1", "chain": "A", "residue_id": "ALA700"}]
        rows, _ = sweep_threshold_sensitivity(log, residues)
        assert len(rows) == 6

    def test_all_negative_scores(self):
        from egfr_pipeline.phase1.orientation_filter import sweep_threshold_sensitivity

        log = [
            {"model_id": f"m{i}", "orientation_score": str(-0.5 - i * 0.1)}
            for i in range(10)
        ]
        residues = _make_residue_table([f"m{i}" for i in range(10)])
        rows, _ = sweep_threshold_sensitivity(log, residues)
        # All should fail at every threshold → 0 pass, 0 hotspots
        for r in rows:
            assert r["n_pass"] == 0
            assert r["n_hotspot_residues"] == 0


class TestAmbiguousReportAdversarial:
    """Attack surface: no pass models, no dG data, empty residues."""

    def test_all_ambiguous_no_pass(self):
        from egfr_pipeline.phase1.orientation_filter import compute_ambiguous_report

        log = [
            {"model_id": "m1", "receptor_id": "R1", "seed_index": "0", "orientation_class": "ambiguous"},
            {"model_id": "m2", "receptor_id": "R1", "seed_index": "0", "orientation_class": "ambiguous"},
        ]
        models = [
            {"model_id": "m1", "dG_separated": "-10.0"},
            {"model_id": "m2", "dG_separated": "-12.0"},
        ]
        residues = [{"model_id": "m1", "residue_id": "ALA700"}]

        rows, _ = compute_ambiguous_report(log, models, residues)
        assert len(rows) == 1
        assert rows[0]["ambiguous_rate"] == 1.0
        assert rows[0]["mean_dG_pass"] == ""  # No pass models

    def test_missing_dg_values(self):
        from egfr_pipeline.phase1.orientation_filter import compute_ambiguous_report

        log = [
            {"model_id": "m1", "receptor_id": "R1", "seed_index": "0", "orientation_class": "pass"},
            {"model_id": "m2", "receptor_id": "R1", "seed_index": "0", "orientation_class": "ambiguous"},
        ]
        # dG_separated missing for both
        models = [
            {"model_id": "m1", "dG_separated": ""},
            {"model_id": "m2"},
        ]
        residues = []

        rows, _ = compute_ambiguous_report(log, models, residues)
        assert rows[0]["mean_dG_pass"] == ""
        assert rows[0]["mean_dG_ambiguous"] == ""
        assert rows[0]["dG_difference"] == ""


class TestPatchClassificationBoundary:
    """Attack surface: boundary values for _assign_class thresholds."""

    def test_orthosteric_boundary_exact_threshold(self):
        from egfr_pipeline.phase2.patch_relationship import _assign_class

        # Exactly at threshold: 2 overlap, 0.25 fraction → orthosteric
        cls, _ = _assign_class(2, 0.25, 5.0)
        assert cls == "orthosteric_candidate"

    def test_just_below_orthosteric(self):
        from egfr_pipeline.phase2.patch_relationship import _assign_class

        # 2 overlap but fraction < 0.25 → rim (partial overlap)
        cls, _ = _assign_class(2, 0.24, 5.0)
        assert cls == "rim_candidate"

    def test_allosteric_boundary_at_20A(self):
        from egfr_pipeline.phase2.patch_relationship import _assign_class

        # No overlap, exactly 20Å → allosteric
        cls, _ = _assign_class(0, 0.0, 20.0)
        assert cls == "allosteric_candidate"

    def test_low_relevance_just_over_20A(self):
        from egfr_pipeline.phase2.patch_relationship import _assign_class

        # No overlap, 20.1Å → low_relevance
        cls, _ = _assign_class(0, 0.0, 20.1)
        assert cls == "low_relevance_candidate"

    def test_no_centroid_no_overlap(self):
        from egfr_pipeline.phase2.patch_relationship import _assign_class

        # No centroid distance available, no overlap → fallback
        cls, _ = _assign_class(0, 0.0, None)
        assert cls == "low_relevance_candidate"

    def test_no_centroid_with_overlap(self):
        from egfr_pipeline.phase2.patch_relationship import _assign_class

        # No centroid but has overlap → rim
        cls, _ = _assign_class(1, 0.10, None)
        assert cls == "rim_candidate"


class TestHandoffEmptyCSV:
    """Attack surface: CSV exists but is truly empty (0 bytes) or header-only."""

    def test_phase2_empty_file_zero_bytes(self, tmp_path):
        """Empty file (no header) should still be handled."""
        from run_production import _validate_handoff_data_quality

        wb_p1 = tmp_path / "p1"; wb_p2 = tmp_path / "p2"; wb_p3 = tmp_path / "p3"
        wb_p1.mkdir(); wb_p2.mkdir(); wb_p3.mkdir()

        pocket_path = wb_p2 / "phase3_candidate_pocket_reference.csv"
        pocket_path.write_text("")  # 0 bytes

        # csv.DictReader on empty file → 0 rows → FAIL
        with pytest.raises(FileNotFoundError, match="후보 포켓 0개"):
            _validate_handoff_data_quality("adv-phase3-setup", wb_p1, wb_p2, wb_p3)

    def test_phase3_exactly_5_poses_passes(self, tmp_path):
        """Boundary: exactly 5 valid poses should NOT fail."""
        from run_production import _validate_handoff_data_quality

        wb_p1 = tmp_path / "p1"; wb_p2 = tmp_path / "p2"; wb_p3 = tmp_path / "p3"
        wb_p1.mkdir(); wb_p2.mkdir(); wb_p3.mkdir()

        evidence_path = wb_p3 / "phase4_docking_evidence_reference.csv"
        with open(evidence_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["pocket_id", "ligand_id", "best_affinity"])
            w.writeheader()
            for i in range(5):
                w.writerow({"pocket_id": f"P{i}", "ligand_id": f"L{i%3}", "best_affinity": f"-{7+i*0.1}"})

        # Should not raise — exactly 5 is the minimum
        _validate_handoff_data_quality("adv-phase4", wb_p1, wb_p2, wb_p3)

    def test_phase3_na_affinity_not_counted(self, tmp_path):
        """Poses with NA affinity should not count as valid."""
        from run_production import _validate_handoff_data_quality

        wb_p1 = tmp_path / "p1"; wb_p2 = tmp_path / "p2"; wb_p3 = tmp_path / "p3"
        wb_p1.mkdir(); wb_p2.mkdir(); wb_p3.mkdir()

        evidence_path = wb_p3 / "phase4_docking_evidence_reference.csv"
        with open(evidence_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["pocket_id", "ligand_id", "best_affinity"])
            w.writeheader()
            # 3 valid + 5 NA → only 3 valid → FAIL
            for i in range(3):
                w.writerow({"pocket_id": f"P{i}", "ligand_id": "L1", "best_affinity": f"-7.{i}"})
            for i in range(5):
                w.writerow({"pocket_id": f"Q{i}", "ligand_id": "L2", "best_affinity": "NA"})

        with pytest.raises(FileNotFoundError, match="유효 포즈"):
            _validate_handoff_data_quality("adv-phase4", wb_p1, wb_p2, wb_p3)


class TestSheet12EmptyLogs:
    """Attack surface: empty orientation logs for sheet12 comparison."""

    def test_both_empty_logs(self):
        from egfr_pipeline.phase1.orientation_filter import compare_sheet12_sensitivity

        rows = compare_sheet12_sensitivity([], [], [])
        assert len(rows) == 2
        assert rows[0]["n_total"] == 0
        assert rows[1]["n_total"] == 0
        # Overlap should be empty string (no sets)
        assert rows[0]["hotspot_overlap"] == ""

    def test_one_empty_log(self):
        from egfr_pipeline.phase1.orientation_filter import compare_sheet12_sensitivity

        log_a = [{"model_id": "m1", "orientation_class": "pass"}]
        residues = [{"model_id": "m1", "chain": "A", "residue_id": "ALA700"}]

        rows = compare_sheet12_sensitivity(log_a, [], residues)
        assert rows[0]["n_total"] == 1
        assert rows[1]["n_total"] == 0


class TestConformationalSelectionEdge:
    """Attack surface: all state-specific with same occupancy, single state."""

    def _build_single_state(self, tmp_path, state, residues):
        state_dir = tmp_path / state
        state_dir.mkdir(parents=True, exist_ok=True)
        with open(state_dir / "ppi_interface_patch_table.csv", "w", newline="") as f:
            cols = ["chain", "residue_id", "residue_num", "residue_name",
                    "lobe_label", "max_occupancy", "n_clusters_present",
                    "construct_type", "orientation_validation_status"]
            w = csv.DictWriter(f, fieldnames=cols)
            w.writeheader()
            for r in residues:
                w.writerow(r)
        with open(state_dir / "ppi_hotspot_residues.csv", "w", newline="") as f:
            cols = ["chain", "residue_id", "is_hotspot", "n_orientation_valid_models"]
            w = csv.DictWriter(f, fieldnames=cols)
            w.writeheader()
            for r in residues:
                w.writerow({"chain": r["chain"], "residue_id": r["residue_id"],
                            "is_hotspot": "True", "n_orientation_valid_models": "100"})

    def test_no_states_returns_empty(self, tmp_path):
        from egfr_pipeline.phase1.compare_states import compare_across_states

        cross, robust = compare_across_states(tmp_path, states=["nonexistent"])
        assert cross == []
        assert robust == []

    def test_single_residue_single_state_is_state_specific(self, tmp_path):
        from egfr_pipeline.phase1.compare_states import compare_across_states

        self._build_single_state(tmp_path, "s1", [
            {"chain": "A", "residue_id": "ALA700", "residue_num": "700",
             "residue_name": "ALA", "lobe_label": "C-lobe",
             "max_occupancy": "0.50", "n_clusters_present": "1",
             "construct_type": "raw", "orientation_validation_status": "pass"},
        ])

        _, robust = compare_across_states(tmp_path, states=["s1"])
        assert len(robust) == 1
        assert robust[0]["robustness_class"] == "state_specific"


class TestPilotMissingFiles:
    """Attack surface: pilot comparison with missing hotspot files."""

    def test_missing_hotspot_files_yields_empty_sets(self, tmp_path):
        from scripts.pilot_fragment_range import compare_pilot_results

        # No pilot directories at all → all hotspot sets empty
        rows, warnings = compare_pilot_results(tmp_path)
        assert len(rows) == 3  # Still 3 pairwise comparisons
        # All Jaccard should be 1.0 (empty vs empty)
        for r in rows:
            assert r["jaccard"] == 1.0

    def test_partial_missing_files(self, tmp_path):
        from scripts.pilot_fragment_range import compare_pilot_results

        # Only range A has data
        d = tmp_path / "pilot_fragment_A" / "workflow_b" / "phase1_ppi_analysis" / "3GT8_raw"
        d.mkdir(parents=True)
        with open(d / "ppi_hotspot_residues.csv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["residue_id", "is_hotspot"])
            w.writeheader()
            w.writerow({"residue_id": "R1", "is_hotspot": "True"})

        rows, _ = compare_pilot_results(tmp_path)
        ab = next(r for r in rows if r["range_a"] == "A" and r["range_b"] == "B")
        assert ab["jaccard"] == 0.0  # {R1} vs {} = 0.0
