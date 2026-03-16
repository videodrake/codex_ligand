from __future__ import annotations

import concurrent.futures
from pathlib import Path

from egfr_pipeline.phase1 import launch_docking as launch_module
from egfr_pipeline.phase1 import prepare_inputs as prepare_inputs_module
from egfr_pipeline.vina import dock as dock_module
import run_production


class _FakeFuture:
    def __init__(self, result):
        self._result = result

    def result(self):
        return self._result


class _FakeExecutor:
    def __init__(self, max_workers):
        self.max_workers = max_workers

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        return False

    def submit(self, fn, *args, **kwargs):
        return _FakeFuture(fn(*args, **kwargs))


def test_phase1_vina_prepares_project_inputs_before_dispatch(monkeypatch) -> None:
    config = {
        "receptors": [
            {
                "id": "3GT8_raw",
                "pdb": "input/receptors/3GT8_raw.pdb",
                "pdbqt": "input/receptors/3GT8_raw_receptor.pdbqt",
            },
            {
                "id": "EGFR_160-185",
                "pdb": "input/receptors/EGFR_160-185.pdb",
                "pdbqt": "input/receptors/EGFR_160-185_receptor.pdbqt",
            },
            {
                "id": "EGFR_170-200",
                "pdb": "input/receptors/EGFR_170-200.pdb",
                "pdbqt": "input/receptors/EGFR_170-200_receptor.pdbqt",
            },
        ],
        "ligands": [
            {
                "id": "173940",
                "pdbqt": "input/ligands/173940_ligand.pdbqt",
            },
        ],
        "vina": {
            "exhaustiveness": 384,
            "n_poses": 100,
        },
    }
    prepared = []
    dispatched = []

    monkeypatch.setattr(run_production, "_load_config", lambda: config)

    def fake_prepare_inputs(payload):
        prepared.append(payload)

    def fake_dock_one_receptor(receptor_entry, ligand_entries, payload):
        dispatched.append((receptor_entry["id"], len(ligand_entries), payload))
        return receptor_entry["id"], {}

    monkeypatch.setattr(
        dock_module,
        "ensure_project_config_inputs_ready",
        fake_prepare_inputs,
    )
    monkeypatch.setattr(dock_module, "dock_one_receptor", fake_dock_one_receptor)
    monkeypatch.setattr(concurrent.futures, "ProcessPoolExecutor", _FakeExecutor)
    monkeypatch.setattr(concurrent.futures, "as_completed", lambda futures: list(futures))

    run_production.phase1_vina()

    assert prepared == [config]
    assert [item[0] for item in dispatched] == [
        "3GT8_raw",
        "EGFR_160-185",
        "EGFR_170-200",
    ]
    assert all(item[1] == 1 for item in dispatched)
    assert all(item[2] is config for item in dispatched)


def test_ppi_docking_dir_prefers_explicit_directory(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr(run_production, "REPO_ROOT", tmp_path)

    target = {
        "name": "3GT8_raw_seed0",
        "docking_dir": "output/phase1_ppi/3GT8_raw/prod_seed0",
    }

    assert run_production._ppi_docking_dir(target) == (
        tmp_path / "output" / "phase1_ppi" / "3GT8_raw" / "prod_seed0"
    )


def test_phase2_ppi_prepares_runtime_inputs_and_dispatches_targets(
    tmp_path,
    monkeypatch,
) -> None:
    config_dir = tmp_path / "config" / "phase1"
    config_dir.mkdir(parents=True)

    targets = []
    for state_name in ["3GT8_raw", "EGFR_160-185", "EGFR_170-200"]:
        config_path = config_dir / f"phase1_prod_{state_name}_seed0.ini"
        config_path.write_text("[Path]\ninput_pdb_name = ignored\n", encoding="utf-8")
        targets.append(
            {
                "name": f"{state_name}_seed0",
                "config_ini": str(Path("config") / "phase1" / config_path.name),
                "docking_dir": str(Path("output") / "phase1_ppi" / state_name / "prod_seed0"),
                "mapping_csv": "",
                "receptor_id": state_name,
                "partner_name": "ext_beta_meander",
                "construct_type": "full_kinase_domain",
                "orientation_validation_status": "not_available",
                "seed_index": 0,
                "is_production": True,
            }
        )

    prepared = []
    dispatched = []

    def fake_prepare_phase1_inputs():
        prepared.append(True)
        return {}

    def fake_run_single(config_path, state_name, seed_index=0, is_production=False, dry_run=False):
        dispatched.append(
            {
                "config_path": config_path,
                "state_name": state_name,
                "seed_index": seed_index,
                "is_production": is_production,
                "dry_run": dry_run,
            }
        )
        return {"status": "completed"}

    monkeypatch.setattr(run_production, "REPO_ROOT", tmp_path)
    monkeypatch.setattr(run_production, "PPI_TARGETS", targets)
    monkeypatch.setattr(run_production, "_FORCE_MODE", True)
    monkeypatch.setattr(prepare_inputs_module, "prepare_phase1_inputs", fake_prepare_phase1_inputs)
    monkeypatch.setattr(launch_module, "run_single", fake_run_single)

    run_production.phase2_ppi()

    assert prepared == [True]
    assert [item["state_name"] for item in dispatched] == [
        "3GT8_raw",
        "EGFR_160-185",
        "EGFR_170-200",
    ]
    assert all(item["seed_index"] == 0 for item in dispatched)
    assert all(item["is_production"] is True for item in dispatched)
    assert all(item["dry_run"] is False for item in dispatched)
    assert all(str(item["config_path"]).startswith(str(config_dir)) for item in dispatched)
