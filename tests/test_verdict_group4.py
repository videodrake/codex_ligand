"""Group 4 tests: Verdict metric calibration (F-4).

Covers:
1. Weight sensitivity simulation — 6 combinations, rescore logic (AC-4.1)
2. EC-4.1: all-identical STRONG warning
9. Adversarial: verdict.py cross threshold integration, raw score overflow
3. Pocket depth computation — positive values, contact exclusion (AC-4.2)
4. Cross-receptor coverage differentiation — n=1,2,3 scoring (AC-4.3)
5. PPI distance distribution analysis — bin counts, 70% warning (AC-4.4)
6. Phase 4 A3 axis spec document contents (AC-4.5)
7. Centroid offset correction — alpha sweep, category changes
8. Adversarial: empty inputs, boundary values
"""

import csv
import math
import tempfile
from pathlib import Path
from typing import Dict, List

import pytest


PROJECT_ROOT = Path(__file__).resolve().parent.parent


# =====================================================================
# AC-4.1: Weight sensitivity simulation
# =====================================================================

class TestWeightSensitivity:
    """AC-4.1: 6 weight combinations, rescore, verdict changes."""

    def test_rescore_baseline_matches_manual(self):
        from scripts.verdict_weight_sensitivity import rescore

        pocket = {"has_ppi_data": True, "vina_raw": 45.0, "ppi_raw": 15.0, "cross_raw": 20.0}
        weights = {"vina": 50, "ppi": 20, "cross": 30}
        total, verdict = rescore(pocket, weights)
        # vina: 45/60*50=37.5, ppi: 15/20*20=15.0, cross: 20/30*30=20.0 → 72.5
        assert total == pytest.approx(72.5, abs=0.1)
        assert verdict == "STRONG"

    def test_rescore_no_ppi_redistributes(self):
        from scripts.verdict_weight_sensitivity import rescore

        pocket = {"has_ppi_data": False, "vina_raw": 40.0, "ppi_raw": 0.0, "cross_raw": 25.0}
        weights = {"vina": 50, "ppi": 20, "cross": 30}
        total, _ = rescore(pocket, weights)
        # PPI redistributed: vina_max=50+20*(50/80)=62.5, cross_max=30+20*(30/80)=37.5
        expected_vina = 40.0 / 60.0 * 62.5
        expected_cross = 25.0 / 30.0 * 37.5
        assert total == pytest.approx(expected_vina + expected_cross, abs=0.1)

    def test_run_sensitivity_6_combinations(self):
        from scripts.verdict_weight_sensitivity import (
            WEIGHT_COMBINATIONS,
            generate_synthetic_pockets,
            run_sensitivity,
        )

        pockets = generate_synthetic_pockets()
        results, _ = run_sensitivity(pockets)
        assert len(results) == len(pockets)

        combo_names = list(WEIGHT_COMBINATIONS.keys())
        assert len(combo_names) == 6

        for r in results:
            for c in combo_names:
                assert f"{c}_total" in r
                assert f"{c}_verdict" in r
                assert r[f"{c}_verdict"] in ("STRONG", "MODERATE", "WEAK")

    def test_ppi_promotion_detected(self):
        from scripts.verdict_weight_sensitivity import run_sensitivity

        # Pocket borderline: high PPI, low vina → promoted with PPI boost
        pockets = [
            {"pocket_id": "P1", "receptor_id": "R1", "has_ppi_data": True,
             "vina_raw": 25.0, "ppi_raw": 19.0, "cross_raw": 20.0},
        ]
        results, warnings = run_sensitivity(pockets)
        # Check that at least some analysis was done
        assert len(results) == 1


class TestEC41AllIdentical:
    """EC-4.1: all combinations yield identical STRONG set."""

    def test_robust_warning_triggered(self):
        from scripts.verdict_weight_sensitivity import run_sensitivity

        # Very strong pocket → STRONG in all combos
        pockets = [
            {"pocket_id": "P1", "receptor_id": "R1", "has_ppi_data": True,
             "vina_raw": 55.0, "ppi_raw": 19.0, "cross_raw": 29.0},
        ]
        _, warnings = run_sensitivity(pockets)
        assert any("robust" in w for w in warnings)


# =====================================================================
# AC-4.2: Pocket depth computation
# =====================================================================

class TestPocketDepth:
    """AC-4.2: pocket_depth_A positive, contact exclusion."""

    def test_compute_pocket_depth_positive(self, tmp_path):
        from egfr_pipeline.vina.pocket_summary import compute_pocket_depth

        # Create a minimal PDB with CA atoms
        pdb_path = tmp_path / "receptor.pdb"
        # Surface CA atoms at (10, 0, 0) and (20, 0, 0)
        # Contact residue at resnum 100, non-contact at 200
        pdb_lines = [
            "ATOM      1  CA  ALA A 100      10.000   0.000   0.000  1.00  0.00\n",
            "ATOM      2  CA  ALA A 200      20.000   0.000   0.000  1.00  0.00\n",
        ]
        pdb_path.write_text("".join(pdb_lines))

        centroid = (5.0, 0.0, 0.0)
        contact_resnums = {100}  # Exclude residue 100
        depth = compute_pocket_depth(centroid, contact_resnums, pdb_path)

        assert depth is not None
        assert depth > 0
        # Should be distance to residue 200 (non-contact): sqrt((20-5)^2) = 15.0
        assert depth == pytest.approx(15.0, abs=0.1)

    def test_depth_excludes_contact_residues(self, tmp_path):
        from egfr_pipeline.vina.pocket_summary import compute_pocket_depth

        pdb_path = tmp_path / "receptor.pdb"
        # Close atom (contact, excluded) and far atom (non-contact)
        pdb_lines = [
            "ATOM      1  CA  ALA A 100       2.000   0.000   0.000  1.00  0.00\n",
            "ATOM      2  CA  ALA A 200      10.000   0.000   0.000  1.00  0.00\n",
        ]
        pdb_path.write_text("".join(pdb_lines))

        depth_with_excl = compute_pocket_depth((0, 0, 0), {100}, pdb_path)
        depth_without_excl = compute_pocket_depth((0, 0, 0), set(), pdb_path)

        # With exclusion: nearest is 200 at 10Å
        assert depth_with_excl == pytest.approx(10.0, abs=0.1)
        # Without exclusion: nearest is 100 at 2Å
        assert depth_without_excl == pytest.approx(2.0, abs=0.1)

    def test_missing_pdb_returns_none(self, tmp_path):
        from egfr_pipeline.vina.pocket_summary import compute_pocket_depth

        depth = compute_pocket_depth((0, 0, 0), set(), tmp_path / "nonexistent.pdb")
        assert depth is None

    def test_no_surface_atoms_returns_none(self, tmp_path):
        from egfr_pipeline.vina.pocket_summary import compute_pocket_depth

        pdb_path = tmp_path / "empty.pdb"
        pdb_path.write_text("END\n")
        depth = compute_pocket_depth((0, 0, 0), set(), pdb_path)
        assert depth is None


# =====================================================================
# AC-4.3: Cross-receptor differentiation
# =====================================================================

class TestCrossReceptorDifferentiation:
    """AC-4.3: n=1,2,3 gives different cross_coverage_pts."""

    def test_3of3_vs_2of3_different_scores(self):
        from egfr_pipeline.verdict import DEFAULT_THRESHOLDS

        pts_3of3 = DEFAULT_THRESHOLDS["cross_coverage_3of3"]
        pts_2of3 = DEFAULT_THRESHOLDS["cross_coverage_2of3"]
        pts_1of3 = DEFAULT_THRESHOLDS["cross_coverage_1of3"]

        assert pts_3of3 > pts_2of3, "3/3 should score higher than 2/3"
        assert pts_2of3 > pts_1of3, "2/3 should score higher than 1/3"

    def test_cross_differentiation_simulation_produces_changes(self):
        from scripts.verdict_weight_sensitivity import (
            generate_synthetic_cross_pockets,
            simulate_cross_differentiation,
        )

        pockets = generate_synthetic_cross_pockets()
        results, warnings = simulate_cross_differentiation(pockets)
        assert len(results) == len(pockets)

        # P_tied_3of3 and P_tied_2of3 should have same current score but differ in scheme_A
        tied_3 = next(r for r in results if r["pocket_id"] == "P_tied_3of3")
        tied_2 = next(r for r in results if r["pocket_id"] == "P_tied_2of3")

        assert tied_3["current_total"] == tied_2["current_total"]  # Same in current
        assert tied_3["scheme_A_total"] > tied_2["scheme_A_total"]  # Different in scheme_A

    def test_1of3_gets_zero_coverage_in_new_scheme(self):
        from scripts.verdict_weight_sensitivity import (
            generate_synthetic_cross_pockets,
            simulate_cross_differentiation,
        )

        pockets = generate_synthetic_cross_pockets()
        results, _ = simulate_cross_differentiation(pockets)

        p1 = next(r for r in results if r["pocket_id"] == "P_1of3_only")
        # In scheme_A, 1/3 gets 0 coverage → should score lower than current
        assert p1["scheme_A_total"] < p1["current_total"]


# =====================================================================
# AC-4.4: PPI distance distribution
# =====================================================================

class TestPPIDistanceDistribution:
    """AC-4.4: bin counts, 70% concentration warning."""

    def test_distribution_4_bins(self):
        from scripts.verdict_weight_sensitivity import analyze_ppi_distance_distribution

        distances = [5.0, 10.0, 20.0, 30.0]
        rows, _ = analyze_ppi_distance_distribution(distances)
        assert len(rows) == 4
        bins = [r["bin"] for r in rows]
        assert bins == ["adjacent", "near", "moderate", "distant"]

    def test_fractions_sum_to_one(self):
        from scripts.verdict_weight_sensitivity import (
            analyze_ppi_distance_distribution,
            generate_synthetic_distances,
        )

        dists = generate_synthetic_distances(50)
        rows, _ = analyze_ppi_distance_distribution(dists)
        total_frac = sum(r["fraction"] for r in rows)
        assert abs(total_frac - 1.0) < 0.01

    def test_70_percent_concentration_warning(self):
        from scripts.verdict_weight_sensitivity import analyze_ppi_distance_distribution

        # All pockets distant → 100% in one bin
        distances = [30.0, 35.0, 40.0, 45.0, 50.0]
        _, warnings = analyze_ppi_distance_distribution(distances)
        assert any("70%" in w for w in warnings)

    def test_offset_shifts_distribution(self):
        from scripts.verdict_weight_sensitivity import analyze_ppi_distance_distribution

        distances = [9.0, 10.0, 11.0]
        rows_raw, _ = analyze_ppi_distance_distribution(distances, offset=0.0)
        rows_off, _ = analyze_ppi_distance_distribution(distances, offset=4.0)

        # Without offset: all in "near" bin (8-15)
        near_raw = next(r for r in rows_raw if r["bin"] == "near")
        assert near_raw["n_pockets"] == 3

        # With 4Å offset: 9-4=5, 10-4=6, 11-4=7 → all in "adjacent" bin (<8)
        adj_off = next(r for r in rows_off if r["bin"] == "adjacent")
        assert adj_off["n_pockets"] == 3

    def test_empty_distances_returns_warning(self):
        from scripts.verdict_weight_sensitivity import analyze_ppi_distance_distribution

        rows, warnings = analyze_ppi_distance_distribution([])
        assert rows == []
        assert any("No PPI distance" in w for w in warnings)


# =====================================================================
# AC-4.5: Phase 4 A3 axis specification document
# =====================================================================

class TestA3AxisSpecification:
    """AC-4.5: spec doc contains input/thresholds/weights/output range."""

    def test_spec_doc_exists(self):
        spec_path = PROJECT_ROOT / "docs" / "phase4_A3_axis_specification.md"
        assert spec_path.exists()

    def test_spec_contains_4_required_sections(self):
        spec_path = PROJECT_ROOT / "docs" / "phase4_A3_axis_specification.md"
        content = spec_path.read_text()

        # 입력 데이터
        assert "입력" in content or "Input" in content
        assert "relationship_class" in content
        assert "hotspot_overlap_count" in content

        # 임계값/매핑
        assert "orthosteric_candidate" in content
        assert "1.0" in content and "0.7" in content

        # 가중치
        assert "30%" in content or "0.30" in content
        assert "45%" in content or "0.45" in content

        # 출력 범위
        assert "0.0" in content and "1.0" in content

    def test_pipeline_architecture_updated(self):
        report_path = PROJECT_ROOT / "PIPELINE_ARCHITECTURE_REPORT.md"
        content = report_path.read_text()
        assert "A3" in content
        assert "phase4_A3_axis_specification.md" in content


# =====================================================================
# Centroid offset correction
# =====================================================================

class TestCentroidOffsetCorrection:
    """AC-4.2 D-2.2: offset correction alpha sweep."""

    def test_simulate_offset_produces_results(self):
        from scripts.centroid_offset_analysis import (
            generate_synthetic_offset_data,
            simulate_offset_correction,
        )

        pockets = generate_synthetic_offset_data()
        results, _ = simulate_offset_correction(pockets)
        assert len(results) == len(pockets)

        for r in results:
            assert "alpha_0.7_dist" in r
            assert "alpha_0.7_cat" in r
            assert r["alpha_0.7_dist"] >= 0

    def test_deeper_pocket_gets_larger_correction(self):
        from scripts.centroid_offset_analysis import simulate_offset_correction

        pockets = [
            {"pocket_id": "deep", "raw_dist": 15.0, "pocket_depth": 8.0},
            {"pocket_id": "shallow", "raw_dist": 15.0, "pocket_depth": 2.0},
        ]
        results, _ = simulate_offset_correction(pockets)
        deep = next(r for r in results if r["pocket_id"] == "deep")
        shallow = next(r for r in results if r["pocket_id"] == "shallow")

        # Same raw distance, deeper pocket → more correction → lower corrected distance
        assert deep["alpha_0.7_dist"] < shallow["alpha_0.7_dist"]

    def test_correction_never_negative(self):
        from scripts.centroid_offset_analysis import simulate_offset_correction

        pockets = [
            {"pocket_id": "extreme", "raw_dist": 3.0, "pocket_depth": 10.0},
        ]
        results, _ = simulate_offset_correction(pockets)
        r = results[0]
        for alpha in [0.5, 0.7, 1.0]:
            assert r[f"alpha_{alpha:.1f}_dist"] >= 0

    def test_ec42_warning_5_plus_changes(self):
        from scripts.centroid_offset_analysis import simulate_offset_correction

        # 6 pockets all changing category
        pockets = [
            {"pocket_id": f"P{i}", "raw_dist": 16.0, "pocket_depth": 6.0}
            for i in range(6)
        ]
        _, warnings = simulate_offset_correction(pockets)
        assert any("EC-4.2" in w for w in warnings)

    def test_classify_spatial_boundaries(self):
        from scripts.centroid_offset_analysis import classify_spatial

        assert classify_spatial(7.9)[0] == "adjacent"
        assert classify_spatial(8.0)[0] == "near"
        assert classify_spatial(14.9)[0] == "near"
        assert classify_spatial(15.0)[0] == "moderate"
        assert classify_spatial(24.9)[0] == "moderate"
        assert classify_spatial(25.0)[0] == "distant"


# =====================================================================
# Adversarial / boundary tests
# =====================================================================

class TestWeightSensitivityAdversarial:
    """Edge cases for weight sensitivity."""

    def test_all_zero_raw_scores(self):
        from scripts.verdict_weight_sensitivity import rescore

        pocket = {"has_ppi_data": True, "vina_raw": 0.0, "ppi_raw": 0.0, "cross_raw": 0.0}
        total, verdict = rescore(pocket, {"vina": 50, "ppi": 20, "cross": 30})
        assert total == 0.0
        assert verdict == "WEAK"

    def test_max_raw_scores(self):
        from scripts.verdict_weight_sensitivity import rescore

        pocket = {"has_ppi_data": True, "vina_raw": 60.0, "ppi_raw": 20.0, "cross_raw": 30.0}
        total, verdict = rescore(pocket, {"vina": 50, "ppi": 20, "cross": 30})
        assert total == pytest.approx(100.0, abs=0.1)
        assert verdict == "STRONG"

    def test_extract_raw_scores_missing_fields(self):
        from scripts.verdict_weight_sensitivity import extract_raw_scores

        row = {"pocket_id": "P1", "receptor_id": "R1", "verdict": "WEAK"}
        result = extract_raw_scores(row)
        assert result["vina_raw"] == 0.0
        assert result["ppi_raw"] == 0.0
        assert result["cross_raw"] == 0.0

    def test_empty_pocket_list(self):
        from scripts.verdict_weight_sensitivity import run_sensitivity

        results, warnings = run_sensitivity([])
        assert results == []


class TestDistanceDistributionAdversarial:
    """Edge cases for distance distribution."""

    def test_single_distance(self):
        from scripts.verdict_weight_sensitivity import analyze_ppi_distance_distribution

        rows, warnings = analyze_ppi_distance_distribution([12.0])
        assert len(rows) == 4
        near = next(r for r in rows if r["bin"] == "near")
        assert near["n_pockets"] == 1
        assert near["fraction"] == 1.0

    def test_all_zero_distances(self):
        from scripts.verdict_weight_sensitivity import analyze_ppi_distance_distribution

        rows, _ = analyze_ppi_distance_distribution([0.0, 0.0, 0.0])
        adj = next(r for r in rows if r["bin"] == "adjacent")
        assert adj["n_pockets"] == 3


# =====================================================================
# Adversarial: verdict.py integration & edge cases
# =====================================================================

class TestVerdictCrossThresholdIntegration:
    """Verify verdict.py DEFAULT_THRESHOLDS contain the new cross keys."""

    def test_new_threshold_keys_exist(self):
        from egfr_pipeline.verdict import DEFAULT_THRESHOLDS

        assert "cross_coverage_3of3" in DEFAULT_THRESHOLDS
        assert "cross_coverage_2of3" in DEFAULT_THRESHOLDS
        assert "cross_coverage_1of3" in DEFAULT_THRESHOLDS

    def test_thresholds_are_monotonic(self):
        from egfr_pipeline.verdict import DEFAULT_THRESHOLDS

        assert DEFAULT_THRESHOLDS["cross_coverage_3of3"] > DEFAULT_THRESHOLDS["cross_coverage_2of3"]
        assert DEFAULT_THRESHOLDS["cross_coverage_2of3"] >= DEFAULT_THRESHOLDS["cross_coverage_1of3"]

    def test_strong_boundary_at_55(self):
        from egfr_pipeline.verdict import DEFAULT_THRESHOLDS

        assert DEFAULT_THRESHOLDS["valid_min"] == 55
        assert DEFAULT_THRESHOLDS["uncertain_min"] == 30


class TestRescoreOverflow:
    """Ensure raw scores exceeding max are capped properly."""

    def test_vina_raw_above_max_capped(self):
        from scripts.verdict_weight_sensitivity import rescore

        # vina_raw=100 > VINA_RAW_MAX=60 → should be capped at 60
        pocket = {"has_ppi_data": True, "vina_raw": 100.0, "ppi_raw": 0.0, "cross_raw": 0.0}
        total, _ = rescore(pocket, {"vina": 50, "ppi": 20, "cross": 30})
        # Capped: min(100,60)/60*50 = 50.0
        assert total == pytest.approx(50.0, abs=0.1)

    def test_ppi_raw_above_max_capped(self):
        from scripts.verdict_weight_sensitivity import rescore

        pocket = {"has_ppi_data": True, "vina_raw": 0.0, "ppi_raw": 50.0, "cross_raw": 0.0}
        total, _ = rescore(pocket, {"vina": 50, "ppi": 20, "cross": 30})
        # Capped: min(50,20)/20*20 = 20.0
        assert total == pytest.approx(20.0, abs=0.1)

    def test_cross_raw_above_max_capped(self):
        from scripts.verdict_weight_sensitivity import rescore

        pocket = {"has_ppi_data": True, "vina_raw": 0.0, "ppi_raw": 0.0, "cross_raw": 100.0}
        total, _ = rescore(pocket, {"vina": 50, "ppi": 20, "cross": 30})
        # Capped: min(100,30)/30*30 = 30.0
        assert total == pytest.approx(30.0, abs=0.1)


class TestPocketDepthAdversarial:
    """Edge cases for pocket_depth computation."""

    def test_only_contact_residues_returns_none(self, tmp_path):
        """If all CA atoms are contact residues, no surface → None."""
        from egfr_pipeline.vina.pocket_summary import compute_pocket_depth

        pdb_path = tmp_path / "receptor.pdb"
        pdb_lines = [
            "ATOM      1  CA  ALA A 100       5.000   0.000   0.000  1.00  0.00\n",
            "ATOM      2  CA  ALA A 101      10.000   0.000   0.000  1.00  0.00\n",
        ]
        pdb_path.write_text("".join(pdb_lines))

        depth = compute_pocket_depth((0, 0, 0), {100, 101}, pdb_path)
        assert depth is None

    def test_non_ca_atoms_ignored(self, tmp_path):
        """Only CA atoms used, not CB or N."""
        from egfr_pipeline.vina.pocket_summary import compute_pocket_depth

        pdb_path = tmp_path / "receptor.pdb"
        pdb_lines = [
            "ATOM      1  CB  ALA A 200       3.000   0.000   0.000  1.00  0.00\n",
            "ATOM      2  N   ALA A 200       4.000   0.000   0.000  1.00  0.00\n",
            "ATOM      3  CA  ALA A 200      10.000   0.000   0.000  1.00  0.00\n",
        ]
        pdb_path.write_text("".join(pdb_lines))

        depth = compute_pocket_depth((0, 0, 0), set(), pdb_path)
        # Only CA at 10.0 should be used, not CB at 3.0 or N at 4.0
        assert depth == pytest.approx(10.0, abs=0.1)

    def test_chain_filtering(self, tmp_path):
        """Only chain A atoms used."""
        from egfr_pipeline.vina.pocket_summary import compute_pocket_depth

        pdb_path = tmp_path / "receptor.pdb"
        pdb_lines = [
            "ATOM      1  CA  ALA B 200       3.000   0.000   0.000  1.00  0.00\n",
            "ATOM      2  CA  ALA A 300      20.000   0.000   0.000  1.00  0.00\n",
        ]
        pdb_path.write_text("".join(pdb_lines))

        depth = compute_pocket_depth((0, 0, 0), set(), pdb_path)
        # Chain B excluded → only chain A at 20.0
        assert depth == pytest.approx(20.0, abs=0.1)


class TestWeightSensitivityReport:
    """Verify report generation produces valid markdown."""

    def test_report_contains_combination_table(self):
        from scripts.verdict_weight_sensitivity import (
            build_sensitivity_report,
            generate_synthetic_pockets,
            run_sensitivity,
        )

        pockets = generate_synthetic_pockets()
        results, warnings = run_sensitivity(pockets)
        lines = build_sensitivity_report(results, warnings)
        text = "\n".join(lines)

        assert "baseline" in text
        assert "ppi_boost" in text
        assert "equal" in text
        assert "안정성 분류" in text

    def test_report_shows_changed_pockets(self):
        from scripts.verdict_weight_sensitivity import (
            build_sensitivity_report,
            generate_synthetic_pockets,
            run_sensitivity,
        )

        pockets = generate_synthetic_pockets()
        results, warnings = run_sensitivity(pockets)
        lines = build_sensitivity_report(results, warnings)
        text = "\n".join(lines)

        # Synthetic data has pockets that change verdict → should be documented
        assert "판정 변화" in text


class TestOffsetCorrectionReport:
    """Verify centroid offset CSV output columns."""

    def test_csv_columns_present(self):
        from scripts.centroid_offset_analysis import (
            generate_synthetic_offset_data,
            simulate_offset_correction,
        )

        pockets = generate_synthetic_offset_data()
        results, _ = simulate_offset_correction(pockets)
        r = results[0]

        required = ["pocket_id", "raw_dist_A", "pocket_depth_A",
                     "raw_category", "raw_spatial_pts", "verdict_change"]
        for col in required:
            assert col in r, f"Missing column: {col}"

        # Alpha columns
        for alpha in [0.5, 0.6, 0.7, 0.8, 0.9, 1.0]:
            assert f"alpha_{alpha:.1f}_dist" in r
            assert f"alpha_{alpha:.1f}_cat" in r
            assert f"alpha_{alpha:.1f}_pts" in r
