import math

import pytest

from egfr_pipeline.phase1.cluster_consensus import (
    compute_cluster_consensus,
    _compute_centroid,
    _compute_centroid_spread,
)

pytestmark = pytest.mark.smoke


def test_cluster_consensus_uses_orientation_pass_models_and_propagates_metadata():
    models = [
        {
            "model_id": "m_pass",
            "cluster_id": "cluster_01",
            "orientation_class": "pass",
            "construct_type": "full_kinase_domain",
            "dG_separated": "-12.5",
            "dSASA": "100.0",
            "sc_value": "0.7",
            "packstat": "0.5",
            "nres_int": "6",
        },
        {
            "model_id": "m_fail",
            "cluster_id": "cluster_01",
            "orientation_class": "fail",
            "construct_type": "full_kinase_domain",
            "dG_separated": "-1.0",
            "dSASA": "40.0",
            "sc_value": "0.2",
            "packstat": "0.1",
            "nres_int": "2",
        },
    ]
    residues = [
        {
            "model_id": "m_pass",
            "construct_type": "full_kinase_domain",
            "chain": "A",
            "residue_id": "LEU819",
            "residue_num": "819",
            "residue_name": "LEU",
            "lobe_label": "N-lobe",
            "delta_e_total": "-1.5",
        },
        {
            "model_id": "m_pass",
            "construct_type": "full_kinase_domain",
            "chain": "B",
            "residue_id": "VAL962",
            "residue_num": "962",
            "residue_name": "VAL",
            "lobe_label": "partner",
            "delta_e_total": "-0.5",
        },
        {
            "model_id": "m_fail",
            "construct_type": "full_kinase_domain",
            "chain": "A",
            "residue_id": "ASP855",
            "residue_num": "855",
            "residue_name": "ASP",
            "lobe_label": "C-lobe",
            "delta_e_total": "-0.2",
        },
    ]

    summaries, hotspots, patches = compute_cluster_consensus(models, residues, "3GT8_raw")

    assert len(summaries) == 1
    summary = summaries[0]
    assert summary["construct_type"] == "full_kinase_domain"
    assert summary["orientation_validation_status"] == "orientation_validated"
    assert summary["n_members_total"] == 2
    assert summary["n_members_orientation_valid"] == 1
    assert summary["representative_model"] == "m_pass"
    assert summary["receptor_hotspot_residues"] == "LEU819"

    receptor_hotspots = [row for row in hotspots if row["chain"] == "A"]
    assert len(receptor_hotspots) == 1
    assert receptor_hotspots[0]["residue_id"] == "LEU819"
    assert receptor_hotspots[0]["construct_type"] == "full_kinase_domain"
    assert receptor_hotspots[0]["orientation_validation_status"] == "orientation_validated"
    assert receptor_hotspots[0]["occupancy"] == 1.0

    patch_row = next(row for row in patches if row["residue_id"] == "LEU819")
    assert patch_row["construct_type"] == "full_kinase_domain"
    assert patch_row["orientation_validation_status"] == "orientation_validated"
    assert patch_row["global_model_count"] == 1
    assert patch_row["global_model_fraction"] == 1.0


def test_cluster_consensus_falls_back_to_all_models_without_orientation_columns():
    models = [
        {
            "model_id": "m1",
            "cluster_id": "cluster_01",
            "construct_type": "full_kinase_domain",
            "dG_separated": "-10.0",
        },
        {
            "model_id": "m2",
            "cluster_id": "cluster_01",
            "construct_type": "full_kinase_domain",
            "dG_separated": "-5.0",
        },
    ]
    residues = [
        {
            "model_id": "m1",
            "construct_type": "full_kinase_domain",
            "chain": "A",
            "residue_id": "LEU819",
            "residue_num": "819",
            "residue_name": "LEU",
            "lobe_label": "N-lobe",
            "delta_e_total": "-1.0",
        },
        {
            "model_id": "m2",
            "construct_type": "full_kinase_domain",
            "chain": "A",
            "residue_id": "LEU819",
            "residue_num": "819",
            "residue_name": "LEU",
            "lobe_label": "N-lobe",
            "delta_e_total": "-0.8",
        },
    ]

    summaries, hotspots, patches = compute_cluster_consensus(models, residues, "3GT8_raw")

    assert len(summaries) == 1
    summary = summaries[0]
    assert summary["orientation_validation_status"] == "not_available"
    assert summary["n_members_total"] == 2
    assert summary["n_members_orientation_valid"] == 2
    assert summary["representative_model"] == "m1"
    assert summary["receptor_hotspot_residues"] == "LEU819"

    hotspot = hotspots[0]
    assert hotspot["orientation_validation_status"] == "not_available"
    assert hotspot["n_orientation_valid_models"] == 2
    assert hotspot["occupancy"] == 1.0

    patch_row = patches[0]
    assert patch_row["orientation_validation_status"] == "not_available"
    assert patch_row["global_model_count"] == 2
    assert patch_row["global_model_fraction"] == 1.0


def test_cluster_consensus_falls_back_when_orientation_labels_are_empty():
    models = [
        {
            "model_id": "m1",
            "cluster_id": "cluster_01",
            "construct_type": "full_kinase_domain",
            "orientation_class": "",
            "dG_separated": "-10.0",
        },
        {
            "model_id": "m2",
            "cluster_id": "cluster_01",
            "construct_type": "full_kinase_domain",
            "orientation_class": "   ",
            "dG_separated": "-5.0",
        },
    ]
    residues = [
        {
            "model_id": "m1",
            "construct_type": "full_kinase_domain",
            "chain": "A",
            "residue_id": "LEU819",
            "residue_num": "819",
            "residue_name": "LEU",
            "lobe_label": "N-lobe",
            "delta_e_total": "-1.0",
        },
        {
            "model_id": "m2",
            "construct_type": "full_kinase_domain",
            "chain": "A",
            "residue_id": "ASP855",
            "residue_num": "855",
            "residue_name": "ASP",
            "lobe_label": "C-lobe",
            "delta_e_total": "-0.8",
        },
    ]

    summaries, hotspots, patches = compute_cluster_consensus(models, residues, "3GT8_raw")

    assert len(summaries) == 1
    summary = summaries[0]
    assert summary["orientation_validation_status"] == "not_available"
    assert summary["n_members_total"] == 2
    assert summary["n_members_orientation_valid"] == 2
    assert summary["representative_model"] == "m1"

    hotspot_ids = {row["residue_id"] for row in hotspots}
    assert hotspot_ids == {"LEU819", "ASP855"}

    for row in patches:
        assert row["orientation_validation_status"] == "not_available"
        assert row["global_model_fraction"] == 0.5


def test_centroid_computed_from_ca_coords():
    """When ca_coords are provided, cluster summary includes centroid_x/y/z."""
    models = [
        {
            "model_id": "m1",
            "cluster_id": "cluster_01",
            "construct_type": "full_kinase_domain",
            "orientation_class": "pass",
            "dG_separated": "-10.0",
        },
        {
            "model_id": "m2",
            "cluster_id": "cluster_01",
            "construct_type": "full_kinase_domain",
            "orientation_class": "pass",
            "dG_separated": "-8.0",
        },
    ]
    residues = [
        # m1 has residues 819 and 822
        {
            "model_id": "m1", "construct_type": "full_kinase_domain",
            "chain": "A", "residue_id": "LEU819", "residue_num": "819",
            "residue_name": "LEU", "lobe_label": "N-lobe", "delta_e_total": "-1.0",
        },
        {
            "model_id": "m1", "construct_type": "full_kinase_domain",
            "chain": "A", "residue_id": "ALA822", "residue_num": "822",
            "residue_name": "ALA", "lobe_label": "N-lobe", "delta_e_total": "-0.5",
        },
        # m2 has only residue 819
        {
            "model_id": "m2", "construct_type": "full_kinase_domain",
            "chain": "A", "residue_id": "LEU819", "residue_num": "819",
            "residue_name": "LEU", "lobe_label": "N-lobe", "delta_e_total": "-0.8",
        },
    ]
    # Fake CA coordinates
    ca_coords = {
        819: (10.0, 20.0, 30.0),
        822: (14.0, 24.0, 34.0),
    }

    summaries, hotspots, patches = compute_cluster_consensus(
        models, residues, "3GT8_raw", ca_coords=ca_coords
    )

    s = summaries[0]
    # Cluster centroid = mean of residues 819 and 822 (union across models)
    assert s["centroid_x"] == 12.0
    assert s["centroid_y"] == 22.0
    assert s["centroid_z"] == 32.0

    # Spread: m1 centroid = mean(819,822) = (12,22,32), m2 centroid = (10,20,30)
    # Cluster centroid = (12,22,32) (from union of 819+822)
    # m1 dist^2 = 0, m2 dist^2 = 4+4+4 = 12
    # RMSD = sqrt((0+12)/2) = sqrt(6)
    assert abs(s["centroid_spread_A"] - math.sqrt(6)) < 0.01

    # Patch table should have CA coordinates for receptor residues
    patch_819 = next(p for p in patches if p["residue_id"] == "LEU819")
    assert patch_819["ca_x"] == 10.0
    assert patch_819["ca_y"] == 20.0
    assert patch_819["ca_z"] == 30.0


def test_centroid_empty_without_ca_coords():
    """When ca_coords are not provided, centroid fields are empty strings."""
    models = [
        {
            "model_id": "m1", "cluster_id": "cluster_01",
            "construct_type": "full_kinase_domain",
            "dG_separated": "-10.0",
        },
    ]
    residues = [
        {
            "model_id": "m1", "construct_type": "full_kinase_domain",
            "chain": "A", "residue_id": "LEU819", "residue_num": "819",
            "residue_name": "LEU", "lobe_label": "N-lobe", "delta_e_total": "-1.0",
        },
    ]

    summaries, _, patches = compute_cluster_consensus(models, residues, "3GT8_raw")

    s = summaries[0]
    assert s["centroid_x"] == ""
    assert s["centroid_y"] == ""
    assert s["centroid_z"] == ""
    assert s["centroid_spread_A"] == ""

    patch_819 = next(p for p in patches if p["residue_id"] == "LEU819")
    assert patch_819["ca_x"] == ""
    assert patch_819["ca_y"] == ""
    assert patch_819["ca_z"] == ""


def test_centroid_spread_none_for_single_member():
    """centroid_spread_A should be empty for clusters with a single model."""
    models = [
        {
            "model_id": "m1", "cluster_id": "cluster_01",
            "construct_type": "full_kinase_domain",
            "orientation_class": "pass",
            "dG_separated": "-10.0",
        },
    ]
    residues = [
        {
            "model_id": "m1", "construct_type": "full_kinase_domain",
            "chain": "A", "residue_id": "LEU819", "residue_num": "819",
            "residue_name": "LEU", "lobe_label": "N-lobe", "delta_e_total": "-1.0",
        },
    ]
    ca_coords = {819: (10.0, 20.0, 30.0)}

    summaries, _, _ = compute_cluster_consensus(
        models, residues, "3GT8_raw", ca_coords=ca_coords
    )

    s = summaries[0]
    assert s["centroid_x"] == 10.0
    assert s["centroid_y"] == 20.0
    assert s["centroid_z"] == 30.0
    # Single model -> spread is undefined
    assert s["centroid_spread_A"] == ""


def test_compute_centroid_basic():
    """Unit test for _compute_centroid."""
    assert _compute_centroid([]) is None
    result = _compute_centroid([(1.0, 2.0, 3.0), (3.0, 4.0, 5.0)])
    assert result == (2.0, 3.0, 4.0)


def test_compute_centroid_spread_basic():
    """Unit test for _compute_centroid_spread."""
    assert _compute_centroid_spread([(1.0, 0.0, 0.0)], (0.0, 0.0, 0.0)) is None
    spread = _compute_centroid_spread(
        [(1.0, 0.0, 0.0), (-1.0, 0.0, 0.0)], (0.0, 0.0, 0.0)
    )
    assert abs(spread - 1.0) < 1e-10
