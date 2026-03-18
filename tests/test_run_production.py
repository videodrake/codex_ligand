from __future__ import annotations

import concurrent.futures
import csv
import json
import sys
from pathlib import Path

from egfr_pipeline.phase1 import launch_docking as launch_module
from egfr_pipeline.phase1 import prepare_inputs as prepare_inputs_module
from egfr_pipeline.vina import vina_executor as vina_executor_module
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


def _write_csv(path: Path, fieldnames: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _write_json(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2), encoding="utf-8")


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
            "cpu": 8,
            "energy_range": 5.0,
            "exhaustiveness": 384,
            "n_poses": 100,
            "parallel_receptors": 2,
            "seed": 20260316,
        },
    }
    prepared = []
    dispatched = []
    executor_workers = []

    monkeypatch.setattr(run_production, "_load_config", lambda: config)
    monkeypatch.setattr(run_production.os, "cpu_count", lambda: 32)

    def fake_prepare_inputs(payload):
        prepared.append(payload)

    def fake_dock_one_receptor(receptor_entry, ligand_entries, payload):
        dispatched.append((receptor_entry["id"], len(ligand_entries), payload))
        return receptor_entry["id"], {}

    monkeypatch.setattr(
        vina_executor_module,
        "ensure_project_config_inputs_ready",
        fake_prepare_inputs,
    )
    monkeypatch.setattr(vina_executor_module, "dock_one_receptor", fake_dock_one_receptor)

    def fake_executor_factory(max_workers):
        executor_workers.append(max_workers)
        return _FakeExecutor(max_workers)

    monkeypatch.setattr(concurrent.futures, "ProcessPoolExecutor", fake_executor_factory)
    monkeypatch.setattr(concurrent.futures, "as_completed", lambda futures: list(futures))

    run_production.phase1_vina()

    assert prepared == [config]
    assert [item[0] for item in dispatched] == [
        "3GT8_raw",
        "EGFR_160-185",
        "EGFR_170-200",
    ]
    assert executor_workers == [2]
    assert all(item[1] == 1 for item in dispatched)
    assert all(item[2] is config for item in dispatched)


def test_phase1_vina_caps_parallelism_to_visible_cores(monkeypatch) -> None:
    config = {
        "receptors": [
            {"id": "3GT8_raw", "pdb": "a.pdb", "pdbqt": "a_receptor.pdbqt"},
            {"id": "EGFR_160-185", "pdb": "b.pdb", "pdbqt": "b_receptor.pdbqt"},
            {"id": "EGFR_170-200", "pdb": "c.pdb", "pdbqt": "c_receptor.pdbqt"},
        ],
        "ligands": [{"id": "173940", "pdbqt": "lig.pdbqt"}],
        "vina": {
            "cpu": 8,
            "energy_range": 5.0,
            "exhaustiveness": 384,
            "n_poses": 100,
            "parallel_receptors": 3,
            "seed": 20260316,
        },
    }
    executor_workers = []

    monkeypatch.setattr(run_production, "_load_config", lambda: config)
    monkeypatch.setattr(run_production.os, "cpu_count", lambda: 16)
    monkeypatch.setattr(vina_executor_module, "ensure_project_config_inputs_ready", lambda payload: None)
    monkeypatch.setattr(vina_executor_module, "dock_one_receptor", lambda receptor_entry, ligand_entries, payload: (receptor_entry["id"], {}))

    def fake_executor_factory(max_workers):
        executor_workers.append(max_workers)
        return _FakeExecutor(max_workers)

    monkeypatch.setattr(concurrent.futures, "ProcessPoolExecutor", fake_executor_factory)
    monkeypatch.setattr(concurrent.futures, "as_completed", lambda futures: list(futures))

    run_production.phase1_vina()

    assert executor_workers == [2]


def test_ppi_docking_dir_prefers_explicit_directory(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr(run_production, "REPO_ROOT", tmp_path)

    target = {
        "name": "3GT8_raw_seed0",
        "docking_dir": "output/workflow_a/phase2_ppi_docking/3GT8_raw/prod_seed0",
    }

    assert run_production._ppi_docking_dir(target) == (
        tmp_path / "output" / "workflow_a" / "phase2_ppi_docking" / "3GT8_raw" / "prod_seed0"
    )


def test_launch_docking_run_single_handles_relative_phase2_output_dir(
    tmp_path: Path,
    monkeypatch,
) -> None:
    project_root = tmp_path
    config_path = project_root / "config" / "phase1" / "phase1_prod_3GT8_raw_seed0.ini"
    config_path.parent.mkdir(parents=True)
    config_path.write_text(
        "\n".join(
            [
                "[Path]",
                "input_pdb_name = input/receptors/3GT8_raw.pdb",
                "[Docking]",
                "total_global_models = 20000",
                "[System]",
                "n_cpus = 16",
                "random_seed = 123",
            ]
        ),
        encoding="utf-8",
    )

    monkeypatch.setattr(launch_module, "PROJECT_ROOT", project_root)
    monkeypatch.setattr(
        launch_module,
        "PHASE1_OUTPUT_DIR",
        Path("output") / "workflow_a" / "phase2_ppi_docking",
    )

    metadata = launch_module.run_single(
        config_path,
        "3GT8_raw",
        seed_index=0,
        is_production=True,
        dry_run=True,
    )

    output_dir = (
        project_root
        / "output"
        / "workflow_a"
        / "phase2_ppi_docking"
        / "3GT8_raw"
        / "prod_seed0"
    )

    assert metadata["config_file"] == "config/phase1/phase1_prod_3GT8_raw_seed0.ini"
    assert metadata["output_dir"] == "output/workflow_a/phase2_ppi_docking/3GT8_raw/prod_seed0"
    assert (output_dir / "pyrosetta_run_metadata.json").exists()


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
                "docking_dir": str(
                    Path("output") / "workflow_a" / "phase2_ppi_docking" / state_name / "prod_seed0"
                ),
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


def test_check_phase4_uses_coverage_to_flag_partial_vina_postprocess(
    tmp_path: Path,
    monkeypatch,
) -> None:
    project_root = tmp_path / "output" / "phase4_project"
    config = {
        "project_name": "phase4_project",
        "output_root": str(tmp_path / "output"),
        "receptors": [{"id": "3GT8_raw"}],
        "ligands": [{"id": "173940"}, {"id": "97806"}],
    }

    for name in (
        "vina_pose_table.csv",
        "vina_pocket_table.csv",
        "vina_drug_pocket_map.csv",
        "vina_pocket_comparison.csv",
        "vina_pocket_bootstrap.csv",
    ):
        _write_csv(project_root / name, ["stub"], [{"stub": "ok"}])

    _write_csv(
        project_root / "vina_postprocess_coverage.csv",
        [
            "receptor_id",
            "ligand_id",
            "raw_pose_file",
            "status",
            "parsed_n_poses",
            "requested_n_poses",
            "coverage_fraction",
            "docking_mode",
            "exhaustiveness",
            "energy_range",
            "cpu_per_job",
            "base_seed",
            "vina_seed",
        ],
        [
            {
                "receptor_id": "3GT8_raw",
                "ligand_id": "173940",
                "raw_pose_file": "3GT8_raw/173940_blind.pdbqt",
                "status": "parsed",
                "parsed_n_poses": "100",
                "requested_n_poses": "100",
                "coverage_fraction": "1.0",
                "docking_mode": "blind",
                "exhaustiveness": "384",
                "energy_range": "5.0",
                "cpu_per_job": "8",
                "base_seed": "20260316",
                "vina_seed": "1",
            },
            {
                "receptor_id": "3GT8_raw",
                "ligand_id": "97806",
                "raw_pose_file": "3GT8_raw/97806_blind.pdbqt",
                "status": "missing",
                "parsed_n_poses": "0",
                "requested_n_poses": "100",
                "coverage_fraction": "0.0",
                "docking_mode": "blind",
                "exhaustiveness": "384",
                "energy_range": "5.0",
                "cpu_per_job": "8",
                "base_seed": "20260316",
                "vina_seed": "2",
            },
        ],
    )

    monkeypatch.setattr(run_production, "_load_config", lambda: config)
    monkeypatch.setattr(run_production, "_project_root", lambda: project_root)

    missing = run_production.check_phase4()
    assert "3GT8_raw/97806=missing" in missing


def test_main_writes_run_status_and_overviews_when_step_view_enabled(
    tmp_path: Path,
    monkeypatch,
) -> None:
    config = {
        "project_name": "telemetry_project",
        "output_root": "./output",
        "step_output_view": {"enabled": True},
        "receptors": [{"id": "3GT8_raw"}],
        "ligands": [{"id": "173940", "pdbqt": "input/ligands/173940_ligand.pdbqt"}],
        "vina": {"exhaustiveness": 384, "n_poses": 100},
    }
    config_path = tmp_path / "config.json"
    project_root = tmp_path / "output" / "telemetry_project"
    call_log: list[str] = []

    _write_json(config_path, config)

    def phase5() -> None:
        call_log.append("phase5")
        _write_csv(
            project_root / "valid_sites.csv",
            ["receptor_id", "pocket_id", "verdict", "confidence_score", "best_affinity"],
            [
                {
                    "receptor_id": "3GT8_raw",
                    "pocket_id": "3GT8_raw_PKT01",
                    "verdict": "STRONG",
                    "confidence_score": "88.0",
                    "best_affinity": "-7.2",
                }
            ],
        )
        _write_csv(
            project_root / "cross_method_agreement.csv",
            ["receptor_id", "pocket_id", "agreement_level"],
            [{"receptor_id": "3GT8_raw", "pocket_id": "3GT8_raw_PKT01", "agreement_level": "HIGH"}],
        )

    monkeypatch.setattr(run_production, "CONFIG_PATH", config_path)
    monkeypatch.setattr(run_production, "REPO_ROOT", tmp_path)
    monkeypatch.setattr(run_production, "PPI_TARGETS", [])
    monkeypatch.setattr(run_production, "_load_config", lambda: config)
    monkeypatch.setattr(run_production, "_project_root", lambda: project_root)
    monkeypatch.setattr(
        run_production,
        "PHASES",
        [(5, "Phase 5", phase5)],
    )
    monkeypatch.setattr(
        run_production,
        "PHASE_CHECKS",
        {
            5: (
                "step5",
                lambda: [] if (project_root / "valid_sites.csv").exists() else ["valid_sites.csv"],
            )
        },
    )
    monkeypatch.setattr(run_production, "banner", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(sys, "argv", ["run_production.py"])

    run_production.main()

    assert call_log == ["phase5"]
    run_status = json.loads((project_root / "run_status.json").read_text(encoding="utf-8"))
    assert run_status["overall_status"] == "completed"
    assert run_status["phase_states"][0]["status"] == "completed"
    assert run_status["summary"]["progress_percent"] == 100
    assert (project_root / "run_overview.md").exists()
    assert (project_root / "run_overview.html").exists()
    assert (project_root / "report_digest.md").exists()
    assert (project_root / "operational_recovery_playbook.md").exists()
    assert (project_root / "operational_recovery_plan.json").exists()
    assert (project_root / "step_index.md").exists()


def test_main_records_skip_reason_when_phase_is_already_complete(
    tmp_path: Path,
    monkeypatch,
) -> None:
    config = {
        "project_name": "skip_project",
        "output_root": "./output",
        "step_output_view": {"enabled": True},
        "receptors": [{"id": "3GT8_raw"}],
        "ligands": [{"id": "173940", "pdbqt": "input/ligands/173940_ligand.pdbqt"}],
        "vina": {"exhaustiveness": 384, "n_poses": 100},
    }
    config_path = tmp_path / "config.json"
    project_root = tmp_path / "output" / "skip_project"
    _write_json(config_path, config)
    _write_csv(
        project_root / "valid_sites.csv",
        ["receptor_id", "pocket_id", "verdict"],
        [{"receptor_id": "3GT8_raw", "pocket_id": "3GT8_raw_PKT01", "verdict": "STRONG"}],
    )
    _write_csv(
        project_root / "cross_method_agreement.csv",
        ["receptor_id", "pocket_id", "agreement_level"],
        [{"receptor_id": "3GT8_raw", "pocket_id": "3GT8_raw_PKT01", "agreement_level": "HIGH"}],
    )

    monkeypatch.setattr(run_production, "CONFIG_PATH", config_path)
    monkeypatch.setattr(run_production, "REPO_ROOT", tmp_path)
    monkeypatch.setattr(run_production, "PPI_TARGETS", [])
    monkeypatch.setattr(run_production, "_load_config", lambda: config)
    monkeypatch.setattr(run_production, "_project_root", lambda: project_root)
    monkeypatch.setattr(run_production, "PHASES", [(5, "Phase 5", lambda: None)])
    monkeypatch.setattr(
        run_production,
        "PHASE_CHECKS",
        {5: ("step5", lambda: [])},
    )
    monkeypatch.setattr(run_production, "banner", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(sys, "argv", ["run_production.py"])

    run_production.main()

    run_status = json.loads((project_root / "run_status.json").read_text(encoding="utf-8"))
    assert run_status["overall_status"] == "completed"
    assert run_status["phase_states"][0]["status"] == "skipped"
    assert run_status["phase_states"][0]["skip_reason"] == "Canonical outputs already exist"
