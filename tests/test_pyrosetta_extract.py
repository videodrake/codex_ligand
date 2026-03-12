"""Tests for PyRosetta PPI residue extraction handoff."""

from __future__ import annotations

import csv
import json
from pathlib import Path

from egfr_pipeline.ppi.pyrosetta_extract import (
    extract_pyrosetta_batch,
    extract_pyrosetta_interface_residues,
)


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
            ["Binding_Residues_A", "Binding_Residues_B", "dG_separated", "dSASA"],
            [
                [
                    "A:LEU819,B:ASP855",
                    "B:VAL962,B:ASN963",
                    "-12.5",
                    "820.0",
                ]
            ],
        )
        _write_csv(
            result_dir / "cluster_results" / "cluster_summary.csv",
            ["Binding_Residues_A"],
            [["A:LEU819,B:ASP855"]],
        )

        data = extract_pyrosetta_interface_residues(
            result_dir,
            receptor_id="3GT8_raw",
            partner_id="MYO1D_beta_meander",
            construct_type="full_kinase_domain",
        )

        rows = data["residue_rows"]
        assert len(rows) == 2
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

        # Partner-side residues from Binding_Residues_B must not be mixed in.
        assert "VAL962" not in residue_map
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
