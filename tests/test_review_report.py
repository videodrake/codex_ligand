from egfr_pipeline.phase1.review_report import _build_patch_reference, _build_report


def test_build_patch_reference_preserves_metadata_from_robustness_rows():
    robustness = [
        {
            "chain": "A",
            "residue_id": "ASP855",
            "residue_num": "855",
            "residue_name": "ASP",
            "lobe_label": "C-lobe",
            "robustness_class": "moderate",
            "n_states_present": "2",
            "states_present": "3GT8_raw;3GT8_cl38_48",
            "global_max_occupancy": "0.75",
            "is_hotspot_any_state": "True",
            "construct_type": "full_kinase_domain",
            "orientation_validation_status": "not_available",
        }
    ]

    rows = _build_patch_reference(robustness, state_summaries={})

    assert len(rows) == 1
    row = rows[0]
    assert row["construct_type"] == "full_kinase_domain"
    assert row["orientation_validation_status"] == "not_available"


def test_build_report_surfaces_mixed_orientation_status_in_handoff_text():
    robustness = [
        {
            "chain": "A",
            "residue_id": "ASP855",
            "residue_num": "855",
            "residue_name": "ASP",
            "lobe_label": "C-lobe",
            "robustness_class": "robust",
            "n_states_present": "3",
            "states_present": "3GT8_raw;3GT8_cl38_48;3GT8_cl85_100",
            "global_max_occupancy": "0.90",
            "is_hotspot_any_state": "True",
            "is_hotspot_all_present_states": "True",
            "construct_type": "full_kinase_domain",
            "orientation_validation_status": "orientation_validated",
        },
        {
            "chain": "A",
            "residue_id": "GLU866",
            "residue_num": "866",
            "residue_name": "GLU",
            "lobe_label": "C-lobe",
            "robustness_class": "moderate",
            "n_states_present": "2",
            "states_present": "3GT8_raw;3GT8_cl38_48",
            "global_max_occupancy": "0.70",
            "is_hotspot_any_state": "True",
            "is_hotspot_all_present_states": "False",
            "construct_type": "full_kinase_domain",
            "orientation_validation_status": "not_available",
        },
    ]

    report_lines = _build_report(
        state_summaries={},
        robustness=robustness,
        cross_state=[],
        receptor_meta=[],
        partner_meta=[],
    )
    report_text = "\n".join(report_lines)

    assert "`construct_type` preserved from upstream consensus: `full_kinase_domain`" in report_text
    assert "`orientation_validation_status` preserved from upstream consensus: `mixed_orientation_status`" in report_text
    assert "orientation status: mixed across contributing states" in report_text
    assert "LightDock remains secondary evidence only." in report_text
    assert "orientation_validation_status = not_available" in report_text
