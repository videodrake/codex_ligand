import sys
import types
from unittest.mock import patch

from egfr_pipeline.phase3 import run_diverse_docking as diverse_module
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


def test_execute_single_job_passes_cpu_per_job_to_vina(tmp_path):
    ligand_file = tmp_path / "ligand.pdbqt"
    receptor_file = tmp_path / "receptor.pdbqt"
    ligand_file.write_text("dummy ligand\n", encoding="utf-8")
    receptor_file.write_text("dummy receptor\n", encoding="utf-8")

    calls = {}
    fake_dock = types.ModuleType("egfr_pipeline.vina.dock")

    def fake_run_docking(*args, **kwargs):
        calls["cpu"] = kwargs.get("cpu")
        return [[-7.2]]

    fake_dock.prepare_ligand = lambda *args, **kwargs: True
    fake_dock.run_docking = fake_run_docking

    job = {
        "job_id": "job-002",
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
        "cpu_per_job": "4",
    }

    with patch.dict(sys.modules, {"egfr_pipeline.vina.dock": fake_dock}):
        result = execute_single_job(job)

    assert result["status"] == "ok"
    assert calls["cpu"] == 4


def test_execute_round_caps_workers_to_allocated_cpu_budget(monkeypatch):
    jobs = [
        {"job_id": "j1", "seed_index": "1", "cpu_per_job": "4"},
        {"job_id": "j2", "seed_index": "1", "cpu_per_job": "4"},
        {"job_id": "j3", "seed_index": "1", "cpu_per_job": "4"},
    ]
    seen = {"workers": None}

    class _FakeFuture:
        def __init__(self, value):
            self._value = value

        def result(self):
            return self._value

    class _FakeExecutor:
        def __init__(self, max_workers):
            seen["workers"] = max_workers

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            return False

        def submit(self, fn, *args, **kwargs):
            return _FakeFuture(fn(*args, **kwargs))

    monkeypatch.setattr(
        diverse_module,
        "resolve_runtime_resources",
        lambda **kwargs: types.SimpleNamespace(effective_cpus=8),
    )
    monkeypatch.setattr(diverse_module, "execute_single_job", lambda job: {"job_id": job["job_id"], "status": "ok"})
    monkeypatch.setattr("concurrent.futures.ThreadPoolExecutor", _FakeExecutor)
    monkeypatch.setattr("concurrent.futures.as_completed", lambda futures: list(futures))

    results = diverse_module.execute_round(jobs, 1, max_workers=16, allocated_cpus=8)

    assert seen["workers"] == 2
    assert len(results) == 3
