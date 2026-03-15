import csv

from egfr_pipeline.phase1.compare_states import compare_across_states, generate_comparison_report


def _write_csv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def test_compare_states_propagates_metadata_from_patch_tables(tmp_path):
    output_dir = tmp_path / "phase1_ppi"
    raw_dir = output_dir / "3GT8_raw"
    cl38_dir = output_dir / "EGFR_160-185"

    _write_csv(
        raw_dir / "ppi_interface_patch_table.csv",
        [
            {
                "chain": "A",
                "residue_id": "LEU819",
                "residue_num": "819",
                "residue_name": "LEU",
                "lobe_label": "N-lobe",
                "n_clusters_present": "1",
                "max_occupancy": "0.80",
                "construct_type": "full_kinase_domain",
                "orientation_validation_status": "orientation_validated",
            }
        ],
    )
    _write_csv(
        raw_dir / "ppi_hotspot_residues.csv",
        [{"chain": "A", "residue_id": "LEU819", "is_hotspot": "True"}],
    )

    _write_csv(
        cl38_dir / "ppi_interface_patch_table.csv",
        [
            {
                "chain": "A",
                "residue_id": "LEU819",
                "residue_num": "819",
                "residue_name": "LEU",
                "lobe_label": "N-lobe",
                "n_clusters_present": "2",
                "max_occupancy": "0.60",
                "construct_type": "full_kinase_domain",
                "orientation_validation_status": "orientation_validated",
            }
        ],
    )
    _write_csv(
        cl38_dir / "ppi_hotspot_residues.csv",
        [{"chain": "A", "residue_id": "LEU819", "is_hotspot": "True"}],
    )

    cross_rows, robustness_rows = compare_across_states(
        output_dir, ["3GT8_raw", "EGFR_160-185"]
    )

    assert len(cross_rows) == 1
    assert len(robustness_rows) == 1
    row = robustness_rows[0]
    assert row["construct_type"] == "full_kinase_domain"
    assert row["orientation_validation_status"] == "orientation_validated"
    assert row["robustness_class"] == "moderate"
    assert row["states_present"] == "3GT8_raw;EGFR_160-185"


def test_compare_states_marks_mixed_orientation_status_in_report(tmp_path):
    output_dir = tmp_path / "phase1_ppi"
    raw_dir = output_dir / "3GT8_raw"
    cl38_dir = output_dir / "EGFR_160-185"

    for state_dir, status in (
        (raw_dir, "orientation_validated"),
        (cl38_dir, "not_available"),
    ):
        _write_csv(
            state_dir / "ppi_interface_patch_table.csv",
            [
                {
                    "chain": "A",
                    "residue_id": "ASP855",
                    "residue_num": "855",
                    "residue_name": "ASP",
                    "lobe_label": "C-lobe",
                    "n_clusters_present": "1",
                    "max_occupancy": "0.75",
                    "construct_type": "full_kinase_domain",
                    "orientation_validation_status": status,
                }
            ],
        )
        _write_csv(
            state_dir / "ppi_hotspot_residues.csv",
            [{"chain": "A", "residue_id": "ASP855", "is_hotspot": "True"}],
        )
        _write_csv(
            state_dir / "ppi_cluster_summary.csv",
            [
                {
                    "n_members_orientation_valid": "1",
                    "receptor_hotspot_residues": "ASP855",
                }
            ],
        )

    cross_rows, robustness_rows = compare_across_states(
        output_dir, ["3GT8_raw", "EGFR_160-185"]
    )
    report_path = generate_comparison_report(
        output_dir, cross_rows, robustness_rows, ["3GT8_raw", "EGFR_160-185"]
    )
    report_text = report_path.read_text(encoding="utf-8")

    assert robustness_rows[0]["orientation_validation_status"] == "mixed_orientation_status"
    assert "Construct type: full_kinase_domain" in report_text
    assert "mixed orientation validation status across contributing states" in report_text
    assert "Use the Phase 1 input validation outputs as the first authority on" in report_text
    assert "potentially ambiguous until resolved" in report_text
