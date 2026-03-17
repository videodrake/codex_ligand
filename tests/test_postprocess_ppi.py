"""Tests for PPI postprocess config registration."""

from __future__ import annotations

import json
from pathlib import Path

from egfr_pipeline.config import load_config
from egfr_pipeline.ppi.postprocess_ppi import _update_config_ppi_dir


def test_update_config_ppi_dir_adds_construct_and_orientation_fields(tmp_path: Path) -> None:
    config_path = tmp_path / "config.json"
    config_path.write_text(
        json.dumps(
            {
                "project_name": "test_project",
                "ppi": {"pyrosetta_result_dirs": {}},
            }
        ),
        encoding="utf-8",
    )

    restored_dir = tmp_path / "restored"
    restored_dir.mkdir()

    _update_config_ppi_dir(
        str(config_path),
        receptor_id="3GT8_raw",
        restored_dir=restored_dir,
        partner_name="MYO1D_beta_meander",
        construct_type="full_kinase_domain",
        orientation_validation_status="orientation_validated",
    )

    config = load_config(str(config_path))
    entry = config["ppi"]["pyrosetta_result_dirs"]["3GT8_raw"][0]
    assert entry["path"] == str(restored_dir)
    assert entry["partner"] == "MYO1D_beta_meander"
    assert entry["construct_type"] == "full_kinase_domain"
    assert entry["orientation_validation_status"] == "orientation_validated"


def test_update_config_ppi_dir_upserts_existing_path_partner(tmp_path: Path) -> None:
    restored_dir = tmp_path / "restored"
    restored_dir.mkdir()
    config_path = tmp_path / "config.json"
    config_path.write_text(
        json.dumps(
            {
                "ppi": {
                    "pyrosetta_result_dirs": {
                        "3GT8_raw": [
                            {
                                "path": str(restored_dir),
                                "partner": "MYO1D_TH1",
                                "construct_type": "unknown_construct_type",
                                "orientation_validation_status": "not_available",
                            }
                        ]
                    }
                }
            }
        ),
        encoding="utf-8",
    )

    _update_config_ppi_dir(
        str(config_path),
        receptor_id="3GT8_raw",
        restored_dir=restored_dir,
        partner_name="MYO1D_TH1",
        construct_type="full_kinase_domain",
        orientation_validation_status="orientation_validated",
    )

    config = load_config(str(config_path))
    entries = config["ppi"]["pyrosetta_result_dirs"]["3GT8_raw"]
    assert len(entries) == 1
    assert entries[0]["construct_type"] == "full_kinase_domain"
    assert entries[0]["orientation_validation_status"] == "orientation_validated"
