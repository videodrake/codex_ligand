import csv

from egfr_pipeline.phase2.patch_ingestion import run_patch_ingestion


PATCH_REFERENCE_COLUMNS = [
    "residue_id",
    "residue_num",
    "residue_name",
    "chain",
    "lobe_label",
    "evidence_source",
    "construct_type",
    "orientation_validation_status",
    "robustness_class",
    "n_states_present",
    "states_present",
    "global_max_occupancy",
    "is_hotspot_any_state",
    "pyrosetta_evidence",
    "lightdock_evidence",
    "method_agreement",
    "confidence",
    "notes",
]


def test_blocked_md_gate_decision_fails_patch_ingestion(tmp_path):
    phase1_dir = tmp_path / "phase1_ppi"
    output_dir = tmp_path / "phase2_pockets"
    phase1_dir.mkdir()

    with open(
        phase1_dir / "phase1_downstream_patch_reference.csv",
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=PATCH_REFERENCE_COLUMNS)
        writer.writeheader()
        writer.writerow(
            {
                "residue_id": "LEU838",
                "residue_num": "838",
                "residue_name": "LEU",
                "chain": "A",
                "lobe_label": "N-lobe",
                "evidence_source": "primary",
                "construct_type": "full_kinase_domain",
                "orientation_validation_status": "orientation_validated",
                "robustness_class": "robust",
                "n_states_present": "3",
                "states_present": "3GT8_raw;3GT8_cl38_48;3GT8_cl85_100",
                "global_max_occupancy": "0.82",
                "is_hotspot_any_state": "True",
                "pyrosetta_evidence": "True",
                "lightdock_evidence": "True",
                "method_agreement": "both",
                "confidence": "high",
                "notes": "",
            }
        )

    (phase1_dir / "md_gate_decision.md").write_text(
        "# MD Gate Decision\n\nPhase 2 blocked due to unstable interface persistence across replicas.\n",
        encoding="utf-8",
    )
    (phase1_dir / "phase1_md_validation_report.md").write_text(
        "# Phase 1 MD Validation Report\n",
        encoding="utf-8",
    )
    (phase1_dir / "md_stability_metrics.csv").write_text(
        "complex_id,replica,status\npose_1,1,unstable\n",
        encoding="utf-8",
    )

    _, report_path = run_patch_ingestion(phase1_dir, output_dir)
    report_text = report_path.read_text(encoding="utf-8")

    assert "## MD Gate Traceability" in report_text
    assert "Decision status: `blocked`" in report_text
    assert "MD gate decision blocks Phase 2 progression" in report_text
