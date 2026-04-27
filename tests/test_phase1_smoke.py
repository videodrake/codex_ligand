import csv
import json

import pytest

from egfr_pipeline.phase1.cluster_consensus import process_state
from egfr_pipeline.phase1.compare_states import (
    CROSS_STATE_COLUMNS,
    ROBUSTNESS_COLUMNS,
    compare_across_states,
    generate_comparison_report,
)
from egfr_pipeline.phase1.lightdock_validation import compute_cross_method_convergence
from egfr_pipeline.phase1 import review_report as review_report_module

pytestmark = pytest.mark.smoke


def _write_csv(path, rows, fieldnames=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        fieldnames = list(rows[0].keys())
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def test_phase1_smoke_chain_runs_through_review_handoff(tmp_path, monkeypatch):
    output_dir = tmp_path / "phase1_ppi"

    for state_name, dg in (("3GT8_raw", "-12.5"), ("EGFR_160-185", "-11.0")):
        state_dir = output_dir / state_name
        _write_csv(
            state_dir / "pyrosetta_interface_models.csv",
            [
                {
                    "model_id": f"{state_name}_m1",
                    "cluster_id": "cluster_01",
                    "orientation_class": "pass",
                    "construct_type": "full_kinase_domain",
                    "dG_separated": dg,
                    "dSASA": "100.0",
                    "sc_value": "0.70",
                    "packstat": "0.50",
                    "nres_int": "6",
                }
            ],
        )
        _write_csv(
            state_dir / "pyrosetta_interface_residue_table.csv",
            [
                {
                    "model_id": f"{state_name}_m1",
                    "chain": "A",
                    "residue_id": "LEU819",
                    "residue_num": "819",
                    "residue_name": "LEU",
                    "lobe_label": "N-lobe",
                    "delta_e_total": "-1.5",
                    "construct_type": "full_kinase_domain",
                },
                {
                    "model_id": f"{state_name}_m1",
                    "chain": "B",
                    "residue_id": "VAL962",
                    "residue_num": "962",
                    "residue_name": "VAL",
                    "lobe_label": "partner",
                    "delta_e_total": "-0.4",
                    "construct_type": "full_kinase_domain",
                },
            ],
        )
        assert process_state(state_name, output_dir) is True

    lightdock_dir = output_dir / "3GT8_raw" / "lightdock"
    lightdock_dir.mkdir(parents=True, exist_ok=True)
    _write_csv(
        lightdock_dir / "lightdock_interface_support_table.csv",
        [
            {
                "model_id": "swarm1_rank1",
                "receptor_id": "3GT8_raw",
                "construct_type": "full_kinase_domain",
                "orientation_validation_status": "not_available",
                "swarm_id": "1",
                "pose_rank": "1",
                "scoring_value": "12.0",
                "chain": "A",
                "residue_id": "LEU819",
                "residue_num": "819",
                "residue_name": "LEU",
                "lobe_label": "N-lobe",
                "source": "lightdock",
            }
        ],
    )
    with open(lightdock_dir / "lightdock_run_metadata.json", "w", encoding="utf-8") as handle:
        json.dump({"construct_type": "full_kinase_domain"}, handle)

    convergence_path = compute_cross_method_convergence("3GT8_raw", output_dir)
    assert convergence_path is not None
    assert convergence_path.exists()

    states = ["3GT8_raw", "EGFR_160-185"]
    cross_rows, robustness_rows = compare_across_states(output_dir, states)
    assert cross_rows
    assert robustness_rows

    _write_csv(output_dir / "ppi_patch_cross_state_comparison.csv", cross_rows, CROSS_STATE_COLUMNS)
    _write_csv(output_dir / "ppi_patch_state_robustness.csv", robustness_rows, ROBUSTNESS_COLUMNS)
    comparison_report = generate_comparison_report(output_dir, cross_rows, robustness_rows, states)
    assert comparison_report.exists()

    input_phase1_dir = tmp_path / "input" / "PPI" / "phase1"
    _write_csv(
        input_phase1_dir / "receptor_metadata.csv",
        [
            {
                "state_name": "3GT8_raw",
                "residue_start": "699",
                "residue_end": "1007",
                "n_residues": "309",
                "source_pdb": "3gt8_raw.pdb",
            },
            {
                "state_name": "EGFR_160-185",
                "residue_start": "699",
                "residue_end": "1007",
                "n_residues": "309",
                "source_pdb": "3gt8_cl38_48.pdb",
            },
        ],
    )
    _write_csv(
        input_phase1_dir / "partner_metadata.csv",
        [
            {
                "partner_id": "extended_beta_meander_955_1006",
                "residue_start": "955",
                "residue_end": "1006",
                "n_residues": "52",
                "active_face_residues": "961;962;963;964;968;969;970;971;972",
            }
        ],
    )

    monkeypatch.setattr(review_report_module, "PROJECT_ROOT", tmp_path)
    report_path, patch_path = review_report_module.generate_review_report(output_dir)

    assert report_path.exists()
    assert patch_path.exists()

    with open(patch_path, encoding="utf-8") as handle:
        patch_rows = list(csv.DictReader(handle))

    assert len(patch_rows) == 1
    assert patch_rows[0]["residue_id"] == "LEU819"
    assert patch_rows[0]["construct_type"] == "full_kinase_domain"
    assert patch_rows[0]["orientation_validation_status"] == "orientation_validated"

    report_text = report_path.read_text(encoding="utf-8")
    assert "LightDock Secondary Validation" in report_text
    assert "`construct_type` preserved from upstream consensus: `full_kinase_domain`" in report_text
