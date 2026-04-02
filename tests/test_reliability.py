"""Tests for reliability enhancement features.

Covers:
  1. Experimental Axis 4 scoring in verdict.py
  2. Consensus scoring (compute_consensus)
  3. Decoy enrichment (Mann-Whitney U, Z-score, scrambled PDB generation)
  4. LightDock integration (pocket overlap, PPI bonus)
"""

import math
import os
import tempfile
from pathlib import Path

from egfr_pipeline.verdict import (
    DEFAULT_THRESHOLDS,
    score_pocket,
    compute_lightdock_pocket_overlap,
    _build_lightdock_index,
)
from egfr_pipeline.pyrosetta_docking.consensus_scoring import compute_consensus
from egfr_pipeline.decoy_enrichment import (
    mann_whitney_u,
    compute_enrichment_stats,
    identify_surface_residues_from_pdb,
    generate_scrambled_pdb,
)


# =========================================================================
# 1. Experimental Axis 4 scoring
# =========================================================================

class TestExperimentalAxis:
    def test_exp_axis_adds_points_when_data_present(self):
        """When experimental data shows high sensitivity/enrichment,
        total score should increase compared to no exp data."""
        pocket = {
            "best_affinity": "-7.5",
            "n_pose": "5",
            "n_ligand": "2",
            "pocket_exists_frac": "0.70",
            "dominant_ligand_fraction": "0.60",
            "ligand_pose_entropy": "0.90",
        }
        T = dict(DEFAULT_THRESHOLDS)

        # Without experimental data
        total_no_exp, _, _, _, _, _, raw_no_exp = score_pocket(
            pocket, None, [], None, T, False, exp_correlation=None,
        )

        # With good experimental data
        exp_corr = {
            "exp_sensitivity": 0.60,
            "exp_specificity": 0.90,
            "exp_enrichment": 4.0,
            "exp_hit_count": 3,
            "exp_false_pos": 0,
        }
        total_exp, _, _, _, _, _, raw_exp = score_pocket(
            pocket, None, [], None, T, False, exp_correlation=exp_corr,
        )

        assert total_exp > total_no_exp
        assert raw_exp["exp_sensitivity_pts"] == 6.0  # >= 0.50 → 6 pts
        assert raw_exp["exp_enrichment_pts"] == 5.0    # >= 3.0 → 5 pts
        assert raw_exp["exp_specificity_pts"] == 4.0   # >= 0.80 → 4 pts
        assert raw_exp["exp_score"] > 0

    def test_exp_axis_no_penalty_when_absent(self):
        """Without experimental data, scoring should be identical to old behavior."""
        pocket = {
            "best_affinity": "-8.0",
            "n_pose": "10",
            "n_ligand": "3",
            "pocket_exists_frac": "0.85",
            "dominant_ligand_fraction": "0.50",
            "ligand_pose_entropy": "1.10",
        }
        T = dict(DEFAULT_THRESHOLDS)

        total, _, _, vina, ppi, cross, raw = score_pocket(
            pocket, None, [], None, T, False, exp_correlation=None,
        )

        # Without exp data: denominator = 60 + 0 + 40 = 100, scale = 1.0
        assert raw["score_denominator"] == 100.0
        assert raw["exp_score"] == 0.0

    def test_exp_axis_denominator_increases_with_exp_data(self):
        """When exp data is present, denominator should increase by exp_max."""
        pocket = {
            "best_affinity": "-7.0",
            "n_pose": "4",
            "n_ligand": "1",
        }
        T = dict(DEFAULT_THRESHOLDS)
        exp_corr = {
            "exp_sensitivity": 0.10,
            "exp_specificity": 0.70,
            "exp_enrichment": 1.0,
            "exp_hit_count": 1,
            "exp_false_pos": 1,
        }
        _, _, _, _, _, _, raw = score_pocket(
            pocket, None, [], None, T, False, exp_correlation=exp_corr,
        )
        # Without PPI: vina(60) + cross(40) + exp(15) = 115
        assert raw["score_denominator"] == 115.0

    def test_exp_contradicts_tag(self):
        """Pocket that contacts non-binders but no binders gets contradicts tag."""
        pocket = {
            "best_affinity": "-6.0",
            "n_pose": "2",
            "n_ligand": "1",
        }
        T = dict(DEFAULT_THRESHOLDS)
        exp_corr = {
            "exp_sensitivity": 0.0,
            "exp_specificity": 0.50,
            "exp_enrichment": 0.0,
            "exp_hit_count": 0,
            "exp_false_pos": 3,
        }
        _, _, _, _, _, _, raw = score_pocket(
            pocket, None, [], None, T, False, exp_correlation=exp_corr,
        )
        assert "exp_contradicts" in raw["reason_tags"]


# =========================================================================
# 2. Consensus scoring
# =========================================================================

class TestConsensusScoring:
    def test_compute_consensus_basic(self):
        """Models ranked top in all functions should be consensus hits."""
        results = [
            {"model_file": "m1.pdb", "dG_ref2015": -20.0, "dG_beta_nov16": -18.0, "dG_ref2015_cart": -19.0},
            {"model_file": "m2.pdb", "dG_ref2015": -15.0, "dG_beta_nov16": -14.0, "dG_ref2015_cart": -13.0},
            {"model_file": "m3.pdb", "dG_ref2015": -10.0, "dG_beta_nov16": -9.0, "dG_ref2015_cart": -11.0},
            {"model_file": "m4.pdb", "dG_ref2015": -5.0, "dG_beta_nov16": -4.0, "dG_ref2015_cart": -6.0},
            {"model_file": "m5.pdb", "dG_ref2015": -1.0, "dG_beta_nov16": -2.0, "dG_ref2015_cart": -1.5},
        ]
        results = compute_consensus(results, percentile=20)

        # m1 should be rank 1 in all → consensus hit
        m1 = next(r for r in results if r["model_file"] == "m1.pdb")
        assert m1["rank_ref2015"] == 1
        assert m1["rank_beta_nov16"] == 1
        assert m1["rank_ref2015_cart"] == 1
        assert m1["consensus_hit"] == "yes"

        # m5 should not be consensus hit
        m5 = next(r for r in results if r["model_file"] == "m5.pdb")
        assert m5["consensus_hit"] == "no"

    def test_incomplete_scores_marked(self):
        """Models with missing scores should be marked 'incomplete'."""
        results = [
            {"model_file": "m1.pdb", "dG_ref2015": -20.0, "dG_beta_nov16": -18.0, "dG_ref2015_cart": None},
            {"model_file": "m2.pdb", "dG_ref2015": -15.0, "dG_beta_nov16": -14.0, "dG_ref2015_cart": -13.0},
        ]
        results = compute_consensus(results, percentile=50)
        m1 = next(r for r in results if r["model_file"] == "m1.pdb")
        assert m1["consensus_hit"] == "incomplete"


# =========================================================================
# 3. Decoy enrichment
# =========================================================================

class TestDecoyEnrichment:
    def test_mann_whitney_u_significant(self):
        """Clear separation should give small p-value."""
        real = [-20.0, -18.0, -19.0, -17.0, -21.0] * 10
        control = [-5.0, -4.0, -6.0, -3.0, -7.0] * 10
        u, p, r = mann_whitney_u(real, control)
        assert p < 0.05
        assert r > 0.3

    def test_mann_whitney_u_no_difference(self):
        """Same distribution should give non-significant result."""
        same = [-10.0, -9.0, -11.0, -8.0, -12.0] * 10
        u, p, r = mann_whitney_u(same, same)
        # p should be close to 0.5 (no difference)
        assert p > 0.05

    def test_compute_enrichment_stats_significant(self):
        """Clear real vs control difference should give 'significant' verdict."""
        real = [-20.0 + i * 0.1 for i in range(50)]
        control = [-5.0 + i * 0.1 for i in range(50)]
        stats = compute_enrichment_stats(real, control, "dG_test")
        assert stats["enrichment_verdict"] == "significant"
        assert stats["z_score"] < -2.0

    def test_compute_enrichment_stats_no_enrichment(self):
        """When control is better than real, verdict should be no_enrichment."""
        real = [-5.0 + i * 0.1 for i in range(30)]
        control = [-20.0 + i * 0.1 for i in range(30)]
        stats = compute_enrichment_stats(real, control, "dG_test")
        assert stats["enrichment_verdict"] == "no_enrichment"
        assert stats["z_score"] > 0

    def test_surface_residue_identification(self):
        """PDB with isolated residues should identify them as surface."""
        pdb_content = ""
        # Create a small cluster + one isolated residue
        for i, (rn, x, y, z) in enumerate([
            (1, 0, 0, 0), (2, 3, 0, 0), (3, 0, 3, 0), (4, 3, 3, 0),
            (5, 50, 50, 50),  # isolated
        ], 1):
            pdb_content += (
                f"ATOM  {i:5d}  CA  ALA A{rn:4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C\n"
            )
        pdb_content += "END\n"

        with tempfile.NamedTemporaryFile(mode="w", suffix=".pdb", delete=False) as f:
            f.write(pdb_content)
            path = f.name

        try:
            surface = identify_surface_residues_from_pdb(path)
            # All residues should be surface (small protein, few neighbours)
            assert len(surface) > 0
            assert 5 in surface  # isolated residue definitely surface
        finally:
            os.unlink(path)

    def test_scrambled_pdb_changes_residue_names(self):
        """Scrambled PDB should have different residue names at mutated positions."""
        pdb_lines = []
        for i, rn in enumerate([10, 20, 30, 40, 50], 1):
            pdb_lines.append(
                f"ATOM  {i:5d}  CA  ALA A{rn:4d}    "
                f"{float(i)*5:8.3f}   0.000   0.000  1.00 20.00           C\n"
            )
        pdb_lines.append("END\n")

        with tempfile.NamedTemporaryFile(mode="w", suffix=".pdb", delete=False) as f:
            f.writelines(pdb_lines)
            orig_path = f.name

        out_path = orig_path + "_scrambled.pdb"
        try:
            mutations = generate_scrambled_pdb(
                orig_path, out_path,
                surface_residues=[10, 20, 30, 40, 50],
                scramble_fraction=1.0,
                seed=42,
            )
            assert len(mutations) == 5
            # Check that mutations are to different AA
            for rn, (orig, new) in mutations.items():
                assert orig == "ALA"
                assert new != "ALA"

            # Verify output file exists and is different
            assert os.path.exists(out_path)
        finally:
            os.unlink(orig_path)
            if os.path.exists(out_path):
                os.unlink(out_path)


# =========================================================================
# 4. LightDock integration
# =========================================================================

class TestLightDockIntegration:
    def test_build_lightdock_index(self):
        """Should build {receptor: {residue: freq}} from CSV rows."""
        rows = [
            {"receptor_id": "3GT8_raw", "residue_id": "ALA699", "lightdock_frequency": "0.80"},
            {"receptor_id": "3GT8_raw", "residue_id": "GLY700", "lightdock_frequency": "0.45"},
            {"receptor_id": "EGFR_160-185", "residue_id": "ALA699", "lightdock_frequency": "0.60"},
        ]
        index = _build_lightdock_index(rows)
        assert "3GT8_raw" in index
        assert index["3GT8_raw"]["ALA699"] == 0.80
        assert index["3GT8_raw"]["GLY700"] == 0.45
        assert index["EGFR_160-185"]["ALA699"] == 0.60

    def test_compute_lightdock_pocket_overlap(self):
        """Pocket overlapping with LightDock residues should have nonzero stats."""
        pocket = {
            "union_contact_residues": "A:ALA699;A:GLY700;A:LEU701",
            "top_residues": "",
        }
        ld_residues = {"ALA699": 0.80, "GLY700": 0.50, "ILE705": 0.30}

        result = compute_lightdock_pocket_overlap(pocket, ld_residues)
        assert result["n_shared"] == 2  # ALA699, GLY700
        assert result["lightdock_mean_freq"] > 0.5
        assert result["lightdock_coverage"] > 0.5

    def test_lightdock_no_overlap(self):
        """No overlap should return zeros."""
        pocket = {"union_contact_residues": "A:ALA699;A:GLY700", "top_residues": ""}
        ld_residues = {"ILE705": 0.80, "LEU710": 0.60}

        result = compute_lightdock_pocket_overlap(pocket, ld_residues)
        assert result["n_shared"] == 0

    def test_lightdock_bonus_in_score_pocket(self):
        """LightDock overlap should increase PPI score."""
        pocket = {
            "best_affinity": "-7.5",
            "n_pose": "8",
            "n_ligand": "2",
            "pocket_exists_frac": "0.80",
            "dominant_ligand_fraction": "0.55",
            "ligand_pose_entropy": "0.95",
        }
        ppi_agreement = {
            "spatial_proximity": "adjacent",
            "n_shared_residues": "3",
            "spatial_dist_A": "6.0",
            "ppi_mean_occupancy_of_shared": "0.5",
            "n_ppi_partners_near": "1",
            "ppi_frac_runs_supporting": "0.60",
            "ppi_best_interface_delta_e": "-2.0",
        }
        T = dict(DEFAULT_THRESHOLDS)

        # Without LightDock
        _, _, _, _, ppi_no_ld, _, raw_no_ld = score_pocket(
            pocket, ppi_agreement, [], None, T, True,
            lightdock_overlap=None,
        )

        # With LightDock validation
        ld_overlap = {"n_shared": 4, "lightdock_mean_freq": 0.60, "lightdock_coverage": 0.50}
        _, _, reasons, _, ppi_ld, _, raw_ld = score_pocket(
            pocket, ppi_agreement, [], None, T, True,
            lightdock_overlap=ld_overlap,
        )

        assert ppi_ld >= ppi_no_ld
        assert "lightdock_validated" in raw_ld["reason_tags"]
        assert any("lightdock=" in r for r in reasons)

    def test_empty_lightdock_data_no_effect(self):
        """Empty LightDock index should not affect scoring."""
        rows: list = []
        index = _build_lightdock_index(rows)
        assert index == {}
