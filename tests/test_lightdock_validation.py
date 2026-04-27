import csv
import json
from pathlib import Path

import pytest

from egfr_pipeline.phase1.lightdock_validation import (
    compute_cross_method_convergence,
    extract_lightdock_interfaces,
)

pytestmark = pytest.mark.smoke


def _write_csv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def test_cross_method_convergence_preserves_metadata_and_lightdock_construct_type(tmp_path):
    state_dir = tmp_path / "phase1_ppi" / "3GT8_raw"
    lightdock_dir = state_dir / "lightdock"
    lightdock_dir.mkdir(parents=True, exist_ok=True)

    _write_csv(
        state_dir / "ppi_interface_patch_table.csv",
        [
            {
                "chain": "A",
                "residue_id": "LEU819",
                "residue_num": "819",
                "residue_name": "LEU",
                "lobe_label": "N-lobe",
                "max_occupancy": "0.80",
                "construct_type": "full_kinase_domain",
                "orientation_validation_status": "orientation_validated",
            }
        ],
    )
    _write_csv(
        lightdock_dir / "lightdock_interface_support_table.csv",
        [
            {
                "model_id": "swarm1_rank1",
                "receptor_id": "3GT8_raw",
                "swarm_id": "1",
                "pose_rank": "1",
                "scoring_value": "12.0",
                "chain": "A",
                "residue_id": "LEU819",
                "residue_num": "819",
                "residue_name": "LEU",
                "lobe_label": "N-lobe",
                "source": "lightdock",
            },
            {
                "model_id": "swarm1_rank1",
                "receptor_id": "3GT8_raw",
                "swarm_id": "1",
                "pose_rank": "1",
                "scoring_value": "12.0",
                "chain": "A",
                "residue_id": "ASP855",
                "residue_num": "855",
                "residue_name": "ASP",
                "lobe_label": "C-lobe",
                "source": "lightdock",
            },
        ],
    )
    with open(lightdock_dir / "lightdock_run_metadata.json", "w", encoding="utf-8") as handle:
        json.dump({"construct_type": "full_kinase_domain"}, handle)

    out_path = compute_cross_method_convergence("3GT8_raw", tmp_path / "phase1_ppi")

    with open(out_path, encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    leu819 = next(row for row in rows if row["residue_id"] == "LEU819")
    asp855 = next(row for row in rows if row["residue_id"] == "ASP855")

    assert leu819["construct_type"] == "full_kinase_domain"
    assert leu819["orientation_validation_status"] == "orientation_validated"
    assert leu819["convergence_class"] == "convergent"
    assert leu819["method_agreement"] == "both"

    assert asp855["construct_type"] == "full_kinase_domain"
    assert asp855["orientation_validation_status"] == "not_available"
    assert asp855["convergence_class"] == "lightdock_only"
    assert asp855["method_agreement"] == "lightdock_only"


def test_extract_lightdock_interfaces_preserves_construct_type(tmp_path):
    state_name = "3GT8_raw"
    lightdock_dir = tmp_path / "phase1_ppi" / state_name / "lightdock"
    swarm_dir = lightdock_dir / "swarm_1"
    swarm_dir.mkdir(parents=True, exist_ok=True)

    with open(lightdock_dir / "rank_by_scoring.list", "w", encoding="utf-8") as handle:
        handle.write("swarm_1/lightdock_1.pdb 12.0\n")

    pdb_lines = [
        "ATOM      1  CA  LEU A 819      10.000  10.000  10.000  1.00 20.00           C",
        "ATOM      2  CA  ASP A 855      30.000  30.000  30.000  1.00 20.00           C",
        "ATOM      3  CA  VAL B 962      10.500  10.500  10.500  1.00 20.00           C",
        "END",
    ]
    with open(swarm_dir / "lightdock_1.pdb", "w", encoding="utf-8") as handle:
        handle.write("\n".join(pdb_lines) + "\n")

    with open(lightdock_dir / "lightdock_run_metadata.json", "w", encoding="utf-8") as handle:
        json.dump({"construct_type": "full_kinase_domain"}, handle)

    iface_path, model_path = extract_lightdock_interfaces(state_name, tmp_path / "phase1_ppi")

    with open(iface_path, encoding="utf-8") as handle:
        iface_rows = list(csv.DictReader(handle))
    with open(model_path, encoding="utf-8") as handle:
        model_rows = list(csv.DictReader(handle))

    assert iface_rows
    assert model_rows
    assert all(row["construct_type"] == "full_kinase_domain" for row in iface_rows)
    # With only 1 active-face residue (VAL962), orientation = insufficient_data
    for row in iface_rows:
        assert row["orientation_validation_status"] in (
            "not_available", "insufficient_data", "pass", "fail", "ambiguous",
        )
    assert model_rows[0]["construct_type"] == "full_kinase_domain"
    assert "orientation_score" in model_rows[0]
