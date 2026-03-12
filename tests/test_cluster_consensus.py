from egfr_pipeline.phase1.cluster_consensus import compute_cluster_consensus


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
