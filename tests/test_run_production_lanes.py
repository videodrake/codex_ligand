from __future__ import annotations

import json
import sys
from pathlib import Path

import run_production


def test_phase2_ppi_selector_runs_single_target(tmp_path: Path, monkeypatch) -> None:
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
                "receptor_id": state_name,
                "seed_index": 0,
                "is_production": True,
            }
        )

    dispatched = []
    monkeypatch.setattr(run_production, "REPO_ROOT", tmp_path)
    monkeypatch.setattr(run_production, "PPI_TARGETS", targets)
    monkeypatch.setattr(run_production, "_FORCE_MODE", True)
    monkeypatch.setattr(
        "egfr_pipeline.phase1.prepare_inputs.prepare_phase1_inputs",
        lambda *args, **kwargs: None,
    )

    def fake_run_single(config_path, state_name, **kwargs):
        dispatched.append((Path(config_path).name, state_name, kwargs))
        return {"status": "completed"}

    monkeypatch.setattr("egfr_pipeline.phase1.launch_docking.run_single", fake_run_single)

    run_production.phase2_ppi(state="EGFR_160-185", seed=0, allocated_cpus=16)

    assert dispatched == [
        (
            "phase1_prod_EGFR_160-185_seed0.ini",
            "EGFR_160-185",
                {
                    "seed_index": 0,
                    "is_production": True,
                    "dry_run": False,
                    "allocated_cpus": 16,
                },
            )
        ]


def test_main_lane_run_writes_marker_and_manifest(tmp_path: Path, monkeypatch) -> None:
    config = {
        "project_name": "lane_project",
        "output_root": str(tmp_path / "output"),
        "step_output_view": {"enabled": False},
        "resources": {
            "vina": {"cpu_per_job": 8, "parallel_receptors": 3},
            "ppi": {"cpu_per_job": 16},
            "lightdock": {"cpu_per_job": 16},
            "phase3": {"cpu_per_job": 4},
            "gpu": {"enabled": False, "engine": "gnina", "cpu_per_job": 4},
        },
    }
    config_path = tmp_path / "config.json"
    config_path.write_text(json.dumps(config), encoding="utf-8")
    project_root = tmp_path / "output" / "lane_project"
    calls = []

    monkeypatch.setattr(run_production, "CONFIG_PATH", config_path)
    monkeypatch.setattr(run_production, "REPO_ROOT", tmp_path)
    monkeypatch.setattr(run_production, "_load_config", lambda: config)
    monkeypatch.setattr(run_production, "_project_root", lambda: project_root)
    monkeypatch.setattr(
        run_production,
        "phase3_ppi_postprocess",
        lambda: calls.append("ppi-post"),
    )
    monkeypatch.setattr(
        run_production,
        "_refresh_lane_step_views",
        lambda lane, cfg: calls.append(f"refresh:{lane}"),
    )
    monkeypatch.setattr(sys, "argv", ["run_production.py", "--lane", "ppi-post", "--allocated-cpus", "4"])

    run_production.main()

    lane_dir = project_root / "logs" / "lanes" / "ppi-post__cpu-vina"
    marker = json.loads((lane_dir / "lane_complete.json").read_text(encoding="utf-8"))
    manifest = json.loads((lane_dir / "runtime_manifest.json").read_text(encoding="utf-8"))

    assert calls == ["ppi-post", "refresh:ppi-post"]
    assert marker["status"] == "completed"
    assert manifest["status"] == "completed"
    assert manifest["resources"]["effective_cpus"] == 4
