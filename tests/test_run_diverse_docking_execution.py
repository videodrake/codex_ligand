import sys
import types
from unittest.mock import patch

from egfr_pipeline.phase3.run_diverse_docking import execute_single_job


def test_execute_single_job_fails_fast_when_ligand_prep_returns_false(tmp_path):
    ligand_file = tmp_path / "ligand.sdf"
    receptor_file = tmp_path / "receptor.pdbqt"
    ligand_file.write_text("dummy ligand\n", encoding="utf-8")
    receptor_file.write_text("dummy receptor\n", encoding="utf-8")

    calls = {"run_docking": 0}
    fake_dock = types.ModuleType("egfr_pipeline.vina.dock")

    def fake_prepare_ligand(input_path, output_path):
        return False

    def fake_run_docking(*args, **kwargs):
        calls["run_docking"] += 1
        raise AssertionError("run_docking should not be called when ligand prep fails")

    fake_dock.prepare_ligand = fake_prepare_ligand
    fake_dock.run_docking = fake_run_docking

    job = {
        "job_id": "job-001",
        "receptor_pdb": str(receptor_file),
        "ligand_file": str(ligand_file),
        "output_dir": str(tmp_path / "out"),
        "output_prefix": "dock_result",
        "center_x": "0.0",
        "center_y": "0.0",
        "center_z": "0.0",
        "box_size_x": "20.0",
        "box_size_y": "20.0",
        "box_size_z": "20.0",
        "exhaustiveness": "8",
        "n_poses": "10",
    }

    with patch.dict(sys.modules, {"egfr_pipeline.vina.dock": fake_dock}):
        result = execute_single_job(job)

    assert result["status"] == "error"
    assert "no PDBQT produced" in result["error_message"]
    assert calls["run_docking"] == 0
