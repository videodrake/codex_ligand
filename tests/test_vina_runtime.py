from __future__ import annotations

import sys
import types
from pathlib import Path

from egfr_pipeline.vina import vina_executor


def test_run_docking_applies_energy_range_to_returned_energies(
    tmp_path: Path,
    monkeypatch,
) -> None:
    calls = {}

    class FakeVina:
        def __init__(self, **kwargs):
            calls["constructor"] = kwargs

        def set_receptor(self, path):
            calls["receptor"] = path

        def set_ligand_from_file(self, path):
            calls["ligand"] = path

        def compute_vina_maps(self, center, box_size):
            calls["maps"] = {"center": center, "box_size": box_size}

        def dock(self, exhaustiveness, n_poses):
            calls["dock"] = {"exhaustiveness": exhaustiveness, "n_poses": n_poses}

        def write_poses(self, output_path, n_poses, energy_range, overwrite):
            calls["write_poses"] = {
                "output_path": output_path,
                "n_poses": n_poses,
                "energy_range": energy_range,
                "overwrite": overwrite,
            }

        def energies(self, n_poses, energy_range):
            calls["energies"] = {"n_poses": n_poses, "energy_range": energy_range}
            return [[-7.5, -0.1, -0.2]]

    monkeypatch.setitem(sys.modules, "vina", types.SimpleNamespace(Vina=FakeVina))

    receptor = tmp_path / "receptor.pdbqt"
    ligand = tmp_path / "ligand.pdbqt"
    output = tmp_path / "out.pdbqt"
    receptor.write_text("RECEPTOR\n", encoding="utf-8")
    ligand.write_text("LIGAND\n", encoding="utf-8")

    energies = vina_executor.run_docking(
        receptor,
        ligand,
        output,
        [1.0, 2.0, 3.0],
        [20.0, 20.0, 20.0],
        exhaustiveness=384,
        n_poses=100,
        cpu=8,
        seed=123,
        energy_range=5.0,
    )

    assert energies == [[-7.5, -0.1, -0.2]]
    assert calls["write_poses"]["energy_range"] == 5.0
    assert calls["energies"]["energy_range"] == 5.0
    assert calls["energies"]["n_poses"] == 100


def test_resolve_ligand_worker_count_caps_cpu_budget() -> None:
    assert vina_executor.resolve_ligand_worker_count(16, 3, cpu_per_job=8, available_cores=32) == 3
    assert vina_executor.resolve_ligand_worker_count(16, 3, cpu_per_job=8, available_cores=16) == 2
    assert vina_executor.resolve_ligand_worker_count(16, 3, cpu_per_job=0, available_cores=32) == 1
