"""Tests for project-config input preparation in Vina docking."""

import os
import sys
import types
from pathlib import Path

from egfr_pipeline.vina import vina_executor


def test_project_config_prepares_missing_pdbqt_inputs(tmp_path: Path, monkeypatch) -> None:
    receptors_dir = tmp_path / "input" / "receptors"
    ligands_dir = tmp_path / "input" / "ligands"
    receptors_dir.mkdir(parents=True)
    ligands_dir.mkdir(parents=True)

    receptor_pdb = receptors_dir / "3GT8_raw.pdb"
    ligand_sdf = ligands_dir / "173940_ligand.sdf"
    receptor_pdb.write_text("ATOM\n", encoding="utf-8")
    ligand_sdf.write_text("mock sdf\n", encoding="utf-8")

    config = {
        "receptors": [
            {
                "id": "3GT8_raw",
                "pdb": str(receptor_pdb),
                "pdbqt": str(receptors_dir / "3GT8_raw_receptor.pdbqt"),
            },
            {
                "id": "EGFR_160-185",
                "pdb": str(receptors_dir / "3GT8_raw.pdb"),
                "pdbqt": str(receptors_dir / "EGFR_160-185_receptor.pdbqt"),
            },
            {
                "id": "EGFR_170-200",
                "pdb": str(receptors_dir / "3GT8_raw.pdb"),
                "pdbqt": str(receptors_dir / "EGFR_170-200_receptor.pdbqt"),
            },
        ],
        "ligands": [
            {
                "id": "173940",
                "sdf": str(ligand_sdf),
                "pdbqt": str(ligands_dir / "173940_ligand.pdbqt"),
            }
        ],
    }

    prepared = {"receptors": [], "ligands": []}

    def fake_prepare_receptor(source: Path, output: Path) -> bool:
        prepared["receptors"].append((Path(source), Path(output)))
        output.write_text("RECEPTOR\n", encoding="utf-8")
        return True

    def fake_prepare_ligand(source: Path, output: Path) -> bool:
        prepared["ligands"].append((Path(source), Path(output)))
        output.write_text("LIGAND\n", encoding="utf-8")
        return True

    monkeypatch.setattr(vina_executor, "prepare_receptor", fake_prepare_receptor)
    monkeypatch.setattr(vina_executor, "prepare_ligand", fake_prepare_ligand)

    vina_executor.validate_project_config(config)
    vina_executor.ensure_project_config_inputs_ready(config)

    assert len(prepared["receptors"]) == 3
    assert len(prepared["ligands"]) == 1
    assert Path(config["receptors"][0]["pdbqt"]).exists()
    assert Path(config["ligands"][0]["pdbqt"]).exists()


def test_project_config_still_fails_without_any_source(tmp_path: Path) -> None:
    config = {
        "receptors": [
            {"id": "3GT8_raw", "pdbqt": str(tmp_path / "3GT8_raw_receptor.pdbqt")},
            {"id": "EGFR_160-185", "pdbqt": str(tmp_path / "EGFR_160-185_receptor.pdbqt")},
            {"id": "EGFR_170-200", "pdbqt": str(tmp_path / "EGFR_170-200_receptor.pdbqt")},
        ],
        "ligands": [
            {"id": "173940", "pdbqt": str(tmp_path / "173940_ligand.pdbqt")},
        ],
    }

    try:
        vina_executor.validate_project_config(config)
    except FileNotFoundError as exc:
        assert "no receptor source file is available" in str(exc)
    else:
        raise AssertionError("validate_project_config should fail when no source file exists")


def test_project_config_rebuilds_stale_pdbqt_inputs(tmp_path: Path, monkeypatch) -> None:
    receptors_dir = tmp_path / "input" / "receptors"
    ligands_dir = tmp_path / "input" / "ligands"
    receptors_dir.mkdir(parents=True)
    ligands_dir.mkdir(parents=True)

    receptor_pdb = receptors_dir / "3GT8_raw.pdb"
    receptor_pdb.write_text("ATOM\n", encoding="utf-8")
    receptor_pdbqt = receptors_dir / "3GT8_raw_receptor.pdbqt"
    receptor_pdbqt.write_text("OLD_RECEPTOR\n", encoding="utf-8")

    ligand_sdf = ligands_dir / "173940_ligand.sdf"
    ligand_sdf.write_text("old sdf\n", encoding="utf-8")
    ligand_pdbqt = ligands_dir / "173940_ligand.pdbqt"
    ligand_pdbqt.write_text("OLD_LIGAND\n", encoding="utf-8")

    os.utime(receptor_pdbqt, (100, 100))
    os.utime(ligand_pdbqt, (100, 100))
    os.utime(receptor_pdb, (200, 200))
    os.utime(ligand_sdf, (200, 200))

    config = {
        "receptors": [
            {
                "id": "3GT8_raw",
                "pdb": str(receptor_pdb),
                "pdbqt": str(receptor_pdbqt),
            },
        ],
        "ligands": [
            {
                "id": "173940",
                "sdf": str(ligand_sdf),
                "pdbqt": str(ligand_pdbqt),
            },
        ],
    }

    rebuilt = {"receptors": 0, "ligands": 0}

    def fake_prepare_receptor(source: Path, output: Path) -> bool:
        rebuilt["receptors"] += 1
        output.write_text("NEW_RECEPTOR\n", encoding="utf-8")
        return True

    def fake_prepare_ligand(source: Path, output: Path) -> bool:
        rebuilt["ligands"] += 1
        output.write_text("NEW_LIGAND\n", encoding="utf-8")
        return True

    monkeypatch.setattr(vina_executor, "prepare_receptor", fake_prepare_receptor)
    monkeypatch.setattr(vina_executor, "prepare_ligand", fake_prepare_ligand)

    vina_executor.ensure_project_config_inputs_ready(config)

    assert rebuilt == {"receptors": 1, "ligands": 1}
    assert receptor_pdbqt.read_text(encoding="utf-8") == "NEW_RECEPTOR\n"
    assert ligand_pdbqt.read_text(encoding="utf-8") == "NEW_LIGAND\n"


def test_derive_docking_seed_is_stable_and_job_specific() -> None:
    seed_a = vina_executor.derive_docking_seed(20260316, "3GT8_raw", "173940")
    seed_b = vina_executor.derive_docking_seed(20260316, "3GT8_raw", "173940")
    seed_c = vina_executor.derive_docking_seed(20260316, "3GT8_raw", "97806")

    assert seed_a == seed_b
    assert seed_a != seed_c
    assert seed_a is not None


def test_print_results_does_not_label_energy_terms_as_rmsd(caplog) -> None:
    import logging

    with caplog.at_level(logging.INFO, logger="pipeline.vina"):
        vina_executor.print_results(
            "173940",
            "blind",
            [[-8.5, -10.1, -0.7], [-8.2, -9.9, -0.4]],
            [0.0, 1.0, 2.0],
            [30.0, 30.0, 30.0],
            display_n=2,
        )

    out = caplog.text
    assert "RMSD_lb" not in out
    assert "Mode" in out
    assert "-8.50" in out
