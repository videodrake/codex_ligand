import json

from egfr_pipeline.ppi.afm_extract import extract_afm_batch


def test_extract_afm_batch_defaults_to_workflow_a_phase3_output(tmp_path):
    model_path = tmp_path / "afm_model.pdb"
    model_path.write_text(
        "\n".join(
            [
                "ATOM      1  CA  ALA A 699      0.000   0.000   0.000  1.00 20.00           C",
                "ATOM      2  CA  GLY B 961      0.000   0.000   6.000  1.00 20.00           C",
                "END",
                "",
            ]
        ),
        encoding="utf-8",
    )

    config_path = tmp_path / "config.json"
    config_path.write_text(
        json.dumps(
            {
                "output_root": str(tmp_path / "output"),
                "project_name": "egfr_myo1d_vina",
                "receptors": [{"id": "3GT8_raw"}],
                "ppi": {
                    "afm_models": {"3GT8_raw": str(model_path)},
                    "afm_settings": {
                        "receptor_chain": "A",
                        "partner_chain": "B",
                        "contact_cutoff": 8.0,
                    },
                },
            }
        ),
        encoding="utf-8",
    )

    residue_csv, summary_csv = extract_afm_batch(str(config_path))

    expected_dir = tmp_path / "output" / "workflow_a" / "phase3_ppi_postprocess"
    assert residue_csv == expected_dir / "ppi_afm_residues.csv"
    assert summary_csv == expected_dir / "ppi_afm_summary.csv"
    assert residue_csv.exists()
    assert summary_csv.exists()
