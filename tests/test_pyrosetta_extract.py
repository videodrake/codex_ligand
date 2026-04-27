"""Tests for PyRosetta PPI residue extraction handoff."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from egfr_pipeline.ppi.pyrosetta_extract import (
    extract_pyrosetta_batch,
    extract_pyrosetta_interface_residues,
)

pytestmark = pytest.mark.smoke


def _write_csv(path: Path, header: list[str], rows: list[list[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)


def test_extract_pyrosetta_interface_residues_preserves_restored_receptor_chain() -> None:
    from tempfile import TemporaryDirectory

    with TemporaryDirectory() as tmp:
        result_dir = Path(tmp)
        _write_csv(
            result_dir / "final_result" / "final_ranking.csv",
            ["Binding_Residues_A", "Binding_Residues_B", "dG_separated", "dSASA", "File_CSV"],
            [
                [
                    "A:LEU819,B:ASP855",
                    "B:VAL962,B:ASN963",
                    "-12.5",
                    "820.0",
                    "Rank01_C01_M01_Energies.csv",
                ]
            ],
        )
        _write_csv(
            result_dir / "cluster_results" / "cluster_summary.csv",
            ["Binding_Residues_A"],
            [["A:LEU819,B:ASP855"]],
        )
        _write_csv(
            result_dir / "final_result" / "Rank01_C01_M01_InterfaceEnergies.csv",
            ["Residue_ID", "Residue_Name", "Chain", "DeltaE_total"],
            [
                ["819", "LEU", "A", "-2.1"],
                ["855", "ASP", "B", "-3.4"],
            ],
        )

        data = extract_pyrosetta_interface_residues(
            result_dir,
            receptor_id="3GT8_raw",
            partner_id="MYO1D_beta_meander",
            construct_type="full_kinase_domain",
        )

        rows = data["residue_rows"]
        assert len(rows) == 4
        residue_map = {row["residue_id"]: row for row in rows}

        assert "LEU819" in residue_map
        assert residue_map["LEU819"]["chain"] == "A"
        assert residue_map["LEU819"]["lobe_label"] == "N-lobe"

        assert "ASP855" in residue_map
        assert residue_map["ASP855"]["chain"] == "B"
        assert residue_map["ASP855"]["lobe_label"] == "C-lobe"
        assert residue_map["ASP855"]["partner_id"] == "MYO1D_beta_meander"
        assert residue_map["ASP855"]["construct_type"] == "full_kinase_domain"
        assert residue_map["ASP855"]["orientation_validation_status"] == "not_available"
        assert residue_map["ASP855"]["mean_interface_delta_e"] == -3.4
        assert residue_map["ASP855"]["best_interface_delta_e"] == -3.4

        # Partner-side residues (chain B) are now included for binding site analysis
        assert "VAL962" in residue_map
        assert residue_map["VAL962"]["lobe_label"] == "partner"
        assert data["summary"]["n_nlobe_interface_residues"] == 1
        assert data["summary"]["n_clobe_interface_residues"] == 1
        assert data["summary"]["orientation_validation_status"] == "not_available"


def test_extract_pyrosetta_batch_writes_extended_residue_schema(tmp_path: Path) -> None:
    result_dir = tmp_path / "ppi_result"
    _write_csv(
        result_dir / "final_result" / "final_ranking.csv",
        ["Binding_Residues_A", "Binding_Residues_B", "dG_separated", "dSASA"],
        [["A:LEU819,B:ASP855", "B:VAL962", "-10.0", "700.0"]],
    )
    _write_csv(
        result_dir / "cluster_results" / "cluster_summary.csv",
        ["Binding_Residues_A"],
        [["A:LEU819,B:ASP855"]],
    )

    config_path = tmp_path / "config.json"
    config_path.write_text(
        json.dumps(
            {
                "project_name": "test_project",
                "output_root": str(tmp_path / "out"),
                "receptors": [{"id": "3GT8_raw"}],
                "ppi": {
                    "pyrosetta_result_dirs": {
                        "3GT8_raw": [
                            {
                                "path": str(result_dir),
                                "partner": "MYO1D_beta_meander",
                                "construct_type": "full_kinase_domain",
                                "orientation_validation_status": "orientation_validated",
                            }
                        ]
                    }
                },
            }
        ),
        encoding="utf-8",
    )

    residue_csv, summary_csv = extract_pyrosetta_batch(str(config_path))

    with open(residue_csv, encoding="utf-8") as handle:
        residue_header = next(csv.reader(handle))
    with open(summary_csv, encoding="utf-8") as handle:
        summary_header = next(csv.reader(handle))

    assert "partner_id" in residue_header
    assert "chain" in residue_header
    assert "residue_name" in residue_header
    assert "lobe_label" in residue_header
    assert "construct_type" in residue_header
    assert "orientation_validation_status" in residue_header

    assert "partner_id" in summary_header
    assert "construct_type" in summary_header
    assert "orientation_validation_status" in summary_header
    assert "n_nlobe_interface_residues" in summary_header
    assert "n_clobe_interface_residues" in summary_header

    with open(residue_csv, encoding="utf-8") as handle:
        residue_rows = list(csv.DictReader(handle))
    with open(summary_csv, encoding="utf-8") as handle:
        summary_rows = list(csv.DictReader(handle))

    assert residue_rows[0]["orientation_validation_status"] == "orientation_validated"
    assert summary_rows[0]["orientation_validation_status"] == "orientation_validated"
    assert (residue_csv.parent / "ppi_pyrosetta_residue_long.csv").exists()
    assert (residue_csv.parent / "ppi_pyrosetta_model_table.csv").exists()


def test_extract_pyrosetta_batch_accepts_runtime_override_dirs(tmp_path: Path) -> None:
    result_dir = tmp_path / "ppi_result"
    _write_csv(
        result_dir / "final_result" / "final_ranking.csv",
        ["Binding_Residues_A", "Binding_Residues_B", "dG_separated", "dSASA"],
        [["A:LEU819,B:ASP855", "B:VAL962", "-10.0", "700.0"]],
    )
    _write_csv(
        result_dir / "cluster_results" / "cluster_summary.csv",
        ["Binding_Residues_A"],
        [["A:LEU819,B:ASP855"]],
    )

    config_path = tmp_path / "config.json"
    config_path.write_text(
        json.dumps(
            {
                "project_name": "test_project",
                "output_root": str(tmp_path / "out"),
                "receptors": [{"id": "3GT8_raw"}],
                "ppi": {},
            }
        ),
        encoding="utf-8",
    )

    residue_csv, summary_csv = extract_pyrosetta_batch(
        str(config_path),
        pyrosetta_result_dirs={
            "3GT8_raw": [
                {
                    "path": str(result_dir),
                    "partner": "MYO1D_beta_meander",
                    "construct_type": "full_kinase_domain",
                    "orientation_validation_status": "orientation_validated",
                }
            ]
        },
    )

    with open(residue_csv, encoding="utf-8") as handle:
        residue_rows = list(csv.DictReader(handle))
    with open(summary_csv, encoding="utf-8") as handle:
        summary_rows = list(csv.DictReader(handle))

    assert residue_rows[0]["partner_id"] == "MYO1D_beta_meander"
    assert residue_rows[0]["orientation_validation_status"] == "orientation_validated"
    assert summary_rows[0]["construct_type"] == "full_kinase_domain"


def test_extract_pyrosetta_interface_residues_resolves_relative_file_csv_energy_paths(
    tmp_path: Path,
) -> None:
    result_dir = tmp_path / "ppi_result"
    _write_csv(
        result_dir / "final_ranking.csv",
        [
            "Binding_Residues_A",
            "Binding_Residues_B",
            "dG_separated",
            "dSASA",
            "File_CSV",
        ],
        [
            [
                "A:THR1940,A:ILE1941",
                "B:LEU972",
                "-12.5",
                "820.0",
                "Rank01_C01_M01_Energies.csv",
            ]
        ],
    )
    _write_csv(
        result_dir / "final_result" / "Rank01_C01_M01_InterfaceEnergies.csv",
        ["Residue_ID", "Residue_Name", "Chain", "DeltaE_total"],
        [
            ["A1940", "THR", "A", "-1.1442"],
            ["A1941", "ILE", "A", "-2.3430"],
            ["B972", "LEU", "B", "-1.1336"],
        ],
    )

    data = extract_pyrosetta_interface_residues(
        result_dir,
        receptor_id="3GT8_raw",
        partner_id="MYO1D_beta_meander",
    )

    residue_map = {row["residue_id"]: row for row in data["residue_rows"]}
    assert residue_map["THR1940"]["mean_interface_delta_e"] == -1.144
    assert residue_map["THR1940"]["best_interface_delta_e"] == -1.144
    assert residue_map["ILE1941"]["mean_interface_delta_e"] == -2.343


def test_extract_pyrosetta_batch_aggregates_run_support_and_seed_provenance(tmp_path: Path) -> None:
    run0 = tmp_path / "phase1_ppi" / "3GT8_raw" / "prod_seed0"
    run1 = tmp_path / "phase1_ppi" / "3GT8_raw" / "prod_seed1"

    for run_dir, seed_index, binding in (
        (run0, 0, "A:LEU819,B:ASP855"),
        (run1, 1, "A:LEU819"),
    ):
        _write_csv(
            run_dir / "final_result" / "final_ranking.csv",
            ["Rank", "Parent", "File_PDB", "Binding_Residues_A", "Binding_Residues_B", "dG_separated", "dSASA"],
            [[
                "1",
                "C01_M01",
                "Rank01_C01_M01.pdb",
                binding,
                "B:VAL962",
                "-10.0",
                "700.0",
            ]],
        )
        _write_csv(
            run_dir / "cluster_results" / "cluster_summary.csv",
            ["Binding_Residues_A"],
            [[binding]],
        )
        (run_dir / "pyrosetta_run_metadata.json").write_text(
            json.dumps(
                {
                    "receptor_id": "3GT8_raw",
                    "partner_id": "MYO1D_beta_meander",
                    "construct_type": "full_kinase_domain",
                    "seed_index": seed_index,
                }
            ),
            encoding="utf-8",
        )

    config_path = tmp_path / "config.json"
    config_path.write_text(
        json.dumps(
            {
                "project_name": "test_project",
                "output_root": str(tmp_path / "out"),
                "receptors": [{"id": "3GT8_raw"}],
                "ppi": {},
            }
        ),
        encoding="utf-8",
    )

    residue_csv, summary_csv = extract_pyrosetta_batch(
        str(config_path),
        pyrosetta_result_dirs={
            "3GT8_raw": [
                {
                    "path": str(run0),
                    "partner": "MYO1D_beta_meander",
                    "construct_type": "full_kinase_domain",
                    "orientation_validation_status": "orientation_validated",
                },
                {
                    "path": str(run1),
                    "partner": "MYO1D_beta_meander",
                    "construct_type": "full_kinase_domain",
                    "orientation_validation_status": "orientation_validated",
                },
            ]
        },
    )

    residue_rows = list(csv.DictReader(residue_csv.open(encoding="utf-8")))
    summary_rows = list(csv.DictReader(summary_csv.open(encoding="utf-8")))
    long_rows = list(csv.DictReader((residue_csv.parent / "ppi_pyrosetta_residue_long.csv").open(encoding="utf-8")))
    model_rows = list(csv.DictReader((residue_csv.parent / "ppi_pyrosetta_model_table.csv").open(encoding="utf-8")))

    residue_map = {row["residue_id"]: row for row in residue_rows}
    assert residue_map["LEU819"]["n_runs_total"] == "2"
    assert residue_map["LEU819"]["n_runs_supporting"] == "2"
    assert residue_map["LEU819"]["frac_runs_supporting"] == "1.0"
    assert residue_map["LEU819"]["supporting_seed_indices"] == "0;1"
    assert residue_map["ASP855"]["n_runs_supporting"] == "1"
    assert residue_map["ASP855"]["frac_runs_supporting"] == "0.5"

    assert summary_rows[0]["n_runs_total"] == "2"
    assert summary_rows[0]["seed_indices"] == "0;1"
    assert summary_rows[0]["orientation_validation_status"] == "orientation_validated"

    assert len(long_rows) >= 3
    assert {row["seed_index"] for row in long_rows} == {"0", "1"}
    assert len(model_rows) == 2
