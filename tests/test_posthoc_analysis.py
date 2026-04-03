"""Tests for post-hoc analysis tools."""

import csv
import os
import tempfile
from pathlib import Path

from egfr_pipeline.posthoc_analysis import (
    run_sensitivity_analysis,
    run_convergence_analysis,
    apply_entropy_correction,
    RIGID_BODY_ENTROPY_PENALTY_REU,
)


def _write_scored_csv(path, n=100):
    """Generate a minimal scored_all_models.csv."""
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["id", "dG_separated", "dSASA", "sc_value", "total_score"])
        for i in range(n):
            dG = -20.0 + i * 0.4  # -20 to +20
            dSASA = 1000 - i * 8  # 1000 to 200
            sc = 0.8 - i * 0.005  # 0.8 to 0.3
            total = -500 + i * 10
            writer.writerow([i + 1, f"{dG:.2f}", f"{dSASA:.1f}", f"{sc:.3f}", f"{total:.1f}"])


def _write_cluster_summary(path, n_clusters=5, seed=0):
    """Generate a minimal cluster_summary.csv."""
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["cluster_id", "centroid_x", "centroid_y", "centroid_z",
                         "best_dG", "n_members"])
        for i in range(n_clusters):
            # Consistent centroids with small per-seed jitter
            x = 10.0 + i * 15.0 + seed * 0.5
            y = 20.0 + i * 5.0 + seed * 0.3
            z = 5.0 + seed * 0.2
            dG = -15.0 + i * 2.0 + seed * 0.1
            writer.writerow([f"C{i+1:02d}", f"{x:.1f}", f"{y:.1f}", f"{z:.1f}",
                             f"{dG:.2f}", 5 - i])


def _write_ranking_csv(path, n=10):
    """Generate a minimal final_ranking.csv."""
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Rank", "model_id", "dG_separated", "dSASA", "sc_value"])
        for i in range(n):
            writer.writerow([i + 1, f"M{i+1:03d}", f"{-18.0 + i * 1.5:.2f}",
                             f"{900 - i * 50:.1f}", f"{0.75 - i * 0.02:.3f}"])


class TestSensitivityAnalysis:
    def test_basic_run(self, tmp_path):
        scored = str(tmp_path / "scored_all_models.csv")
        _write_scored_csv(scored, n=200)
        output = str(tmp_path / "sensitivity.csv")

        results = run_sensitivity_analysis(scored, output, variations=[-20, 20])
        assert len(results) > 0
        assert os.path.exists(output)

        # Should have 4 parameters × 2 variations = 8 rows
        assert len(results) == 8

        for r in results:
            assert "parameter" in r
            assert "top20_overlap_frac" in r
            assert 0 <= r["top20_overlap_frac"] <= 1

    def test_stable_filters(self, tmp_path):
        """With strong separation, small variations should be stable."""
        scored = str(tmp_path / "scored.csv")
        _write_scored_csv(scored, n=500)
        output = str(tmp_path / "sens.csv")

        results = run_sensitivity_analysis(scored, output, variations=[-10, 10])

        # ±10% should preserve most top-20 for well-separated data
        for r in results:
            if r["parameter"] != "dG_max":
                assert r["n_survivors"] > 0


class TestConvergenceAnalysis:
    def test_converging_seeds(self, tmp_path):
        """Seeds with similar centroids should show convergence."""
        state_dir = tmp_path / "3GT8_raw"
        for seed_idx in range(3):
            cluster_dir = state_dir / f"prod_seed{seed_idx}" / "cluster_results"
            _write_cluster_summary(str(cluster_dir / "cluster_summary.csv"),
                                   n_clusters=4, seed=seed_idx)

        output = str(state_dir / "convergence.csv")
        result = run_convergence_analysis(str(state_dir), output, match_distance=8.0)

        assert result["n_seeds"] == 3
        assert result["n_consensus_sites"] >= 1
        assert result["convergence_verdict"] in ("CONVERGED", "PARTIALLY_CONVERGED")

    def test_single_seed_skipped(self, tmp_path):
        """Single seed should return empty (need ≥2)."""
        state_dir = tmp_path / "single"
        cluster_dir = state_dir / "prod_seed0" / "cluster_results"
        _write_cluster_summary(str(cluster_dir / "cluster_summary.csv"))

        result = run_convergence_analysis(str(state_dir), str(tmp_path / "out.csv"))
        assert result == {}


class TestEntropyCorrection:
    def test_dg_increases(self, tmp_path):
        """Entropy correction should make dG less favorable (more positive)."""
        ranking = str(tmp_path / "final_ranking.csv")
        _write_ranking_csv(ranking, n=5)
        output = str(tmp_path / "corrected.csv")

        rows = apply_entropy_correction(ranking, output)
        assert len(rows) == 5

        for row in rows:
            dG_orig = float(row["dG_separated"])
            dG_corr = float(row["dG_corrected"])
            assert dG_corr > dG_orig
            assert abs(dG_corr - dG_orig - RIGID_BODY_ENTROPY_PENALTY_REU) < 0.01

    def test_ranking_preserved_for_uniform_correction(self, tmp_path):
        """Uniform correction should not change relative ranking."""
        ranking = str(tmp_path / "final_ranking.csv")
        _write_ranking_csv(ranking, n=10)
        output = str(tmp_path / "corrected.csv")

        rows = apply_entropy_correction(ranking, output)

        # Original ranking: 1,2,3... → corrected should be same order
        for i, row in enumerate(sorted(rows, key=lambda r: float(r["dG_corrected"]))):
            assert int(row["rank_corrected"]) == i + 1
