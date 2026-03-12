"""Tests for project-config input preparation in Vina docking."""

from pathlib import Path

from egfr_pipeline.vina import dock


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
                "id": "3GT8_cl38_48",
                "pdb": str(receptors_dir / "3GT8_raw.pdb"),
                "pdbqt": str(receptors_dir / "3GT8_cl38_48_receptor.pdbqt"),
            },
            {
                "id": "3GT8_cl85_100",
                "pdb": str(receptors_dir / "3GT8_raw.pdb"),
                "pdbqt": str(receptors_dir / "3GT8_cl85_100_receptor.pdbqt"),
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

    monkeypatch.setattr(dock, "prepare_receptor", fake_prepare_receptor)
    monkeypatch.setattr(dock, "prepare_ligand", fake_prepare_ligand)

    dock.validate_project_config(config)
    dock.ensure_project_config_inputs_ready(config)

    assert len(prepared["receptors"]) == 3
    assert len(prepared["ligands"]) == 1
    assert Path(config["receptors"][0]["pdbqt"]).exists()
    assert Path(config["ligands"][0]["pdbqt"]).exists()


def test_project_config_still_fails_without_any_source(tmp_path: Path) -> None:
    config = {
        "receptors": [
            {"id": "3GT8_raw", "pdbqt": str(tmp_path / "3GT8_raw_receptor.pdbqt")},
            {"id": "3GT8_cl38_48", "pdbqt": str(tmp_path / "3GT8_cl38_48_receptor.pdbqt")},
            {"id": "3GT8_cl85_100", "pdbqt": str(tmp_path / "3GT8_cl85_100_receptor.pdbqt")},
        ],
        "ligands": [
            {"id": "173940", "pdbqt": str(tmp_path / "173940_ligand.pdbqt")},
        ],
    }

    try:
        dock.validate_project_config(config)
    except FileNotFoundError as exc:
        assert "no receptor source file is available" in str(exc)
    else:
        raise AssertionError("validate_project_config should fail when no source file exists")
