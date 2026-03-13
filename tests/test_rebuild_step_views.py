"""Tests for rebuilding derived step views from existing canonical outputs."""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

from egfr_pipeline.validate import ValidationResult


RECEPTOR_IDS = ["3GT8_raw", "3GT8_cl38_48", "3GT8_cl85_100"]
LIGAND_IDS = ["173940", "97806", "VAX-C12_0"]
SCRIPT_PATH = Path(__file__).resolve().parents[1] / "scripts" / "rebuild_step_views.py"


def _write_json(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2), encoding="utf-8")


def _write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def _read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _make_config(tmp_path: Path) -> tuple[Path, Path, Path, Path]:
    output_root = tmp_path / "output"
    th1_dir = tmp_path / "ppi_raw" / "TH1_run"
    beta_dir = tmp_path / "ppi_raw" / "beta_run"
    config = {
        "project_name": "rebuild_project",
        "output_root": str(output_root),
        "mode": "blind",
        "vina": {"mode": "blind", "exhaustiveness": 384, "n_poses": 100},
        "receptors": [{"id": receptor_id} for receptor_id in RECEPTOR_IDS],
        "ligands": [
            {"id": "173940", "pdbqt": "input/ligands/173940_ligand.pdbqt"},
            {"id": "97806", "pdbqt": "input/ligands/97806_ligand.pdbqt"},
            {"id": "VAX-C12_0", "pdbqt": "input/ligands/VAX-C12_0_ligand.pdbqt"},
        ],
        "ppi": {
            "pyrosetta_result_dirs": {
                "3GT8_raw": [
                    {"partner": "MYO1D_TH1", "path": str(th1_dir)},
                    {"partner": "MYO1D_beta", "path": str(beta_dir)},
                ]
            }
        },
    }
    config_path = tmp_path / "config.json"
    _write_json(config_path, config)
    return config_path, output_root / config["project_name"], th1_dir, beta_dir


def _seed_step1_pose_tree(project_root: Path) -> None:
    for receptor_id in RECEPTOR_IDS:
        for ligand_id in LIGAND_IDS:
            _write_text(
                project_root / receptor_id / f"{ligand_id}_blind.pdbqt",
                "MODEL 1\nENDMDL\nMODEL 2\nENDMDL\n",
            )


def _seed_step2_run_dir(run_dir: Path, *, target: str) -> None:
    _write_text(
        run_dir / "final_result" / "final_ranking.csv",
        "rank,score\n1,-20.5\n",
    )
    _write_json(
        run_dir / "pyrosetta_run_metadata.json",
        {"status": "complete", "n_models": 20000, "target": target},
    )


def _seed_project_root(project_root: Path) -> None:
    _write_text(
        project_root / "vina_pose_table.csv",
        "receptor_id,ligand_id\n3GT8_raw,173940\n",
    )
    _write_text(
        project_root / "vina_pocket_table.csv",
        "receptor_id,pocket_id\n3GT8_raw,3GT8_raw_PKT01\n",
    )
    _write_text(
        project_root / "vina_drug_pocket_map.csv",
        "receptor_id,ligand_id,dominant_pocket_id\n3GT8_raw,173940,3GT8_raw_PKT01\n",
    )
    _write_text(
        project_root / "vina_pocket_comparison.csv",
        "receptor_a,pocket_a,receptor_b,pocket_b\n3GT8_raw,3GT8_raw_PKT01,3GT8_cl38_48,3GT8_cl38_48_PKT01\n",
    )
    _write_text(
        project_root / "vina_pocket_bootstrap.csv",
        "receptor_id,pocket_id,pocket_exists_frac\n3GT8_raw,3GT8_raw_PKT01,0.9\n",
    )
    _write_text(
        project_root / "ppi_pyrosetta_residues.csv",
        "receptor_id,residue_id,occupancy\n3GT8_raw,A:ALA699,0.8\n",
    )
    _write_text(
        project_root / "ppi_pyrosetta_summary.csv",
        "receptor_id,n_final_models\n3GT8_raw,20\n",
    )
    _write_text(
        project_root / "valid_sites.csv",
        "receptor_id,pocket_id,verdict\n3GT8_raw,3GT8_raw_PKT01,STRONG\n",
    )
    _write_text(
        project_root / "cross_method_agreement.csv",
        "receptor_id,pocket_id,agreement_level\n3GT8_raw,3GT8_raw_PKT01,HIGH\n",
    )
    _write_text(
        project_root / "vina_consensus_sites.csv",
        "consensus_site_id,n_receptors\nCS01,3\n",
    )
    _write_text(project_root / "project_report.txt", "Step 6 narrative report.\n")
    _write_text(
        project_root / "combined_residue_evidence.csv",
        "receptor_id,residue_id,evidence\n3GT8_raw,A:ALA699,vina\n",
    )


def _load_script_module():
    module_name = "test_rebuild_step_views_module"
    sys.modules.pop(module_name, None)
    spec = importlib.util.spec_from_file_location(module_name, SCRIPT_PATH)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


def _ok_validation_result() -> ValidationResult:
    result = ValidationResult()
    result.ok("validation ok")
    return result


def test_rebuild_step_views_backfills_existing_outputs(monkeypatch, tmp_path: Path) -> None:
    config_path, project_root, th1_dir, beta_dir = _make_config(tmp_path)
    _seed_step1_pose_tree(project_root)
    _seed_step2_run_dir(th1_dir, target="TH1")
    _seed_step2_run_dir(beta_dir, target="beta-meander")
    _seed_project_root(project_root)
    _write_text(
        tmp_path / "output" / "phase1_ppi" / "phase1_interface_report.md",
        "# Phase 1 Interface Report\n",
    )

    module = _load_script_module()
    monkeypatch.setattr(module, "REPO_ROOT", tmp_path)
    monkeypatch.setattr(
        module,
        "PPI_TARGETS",
        [
            {"name": "TH1", "config_ini": "config/ppi_prod_TH1.ini", "docking_dir": th1_dir},
            {
                "name": "beta-meander",
                "config_ini": "config/ppi_prod_beta_meander.ini",
                "docking_dir": beta_dir,
            },
        ],
    )
    monkeypatch.setattr(module, "run_validation", lambda *_args, **_kwargs: _ok_validation_result())
    monkeypatch.setattr(
        sys,
        "argv",
        ["rebuild_step_views.py", "--config", str(config_path)],
    )

    exit_code = module.main()

    manifest = _read_json(project_root / "current_run_manifest.json")
    assert exit_code == 0
    assert manifest["execution_mode"] == "backfill:1,2,3,4,5,6,7"
    assert manifest["step_status"] == {
        "step1": "complete",
        "step2": "complete",
        "step3": "complete",
        "step4": "complete",
        "step5": "complete",
        "step6": "complete",
        "step7": "passed",
    }
    assert (project_root / "step1_vina_raw" / "raw_pose_index.csv").exists()
    assert (project_root / "step2_ppi_raw" / "raw_run_paths.tsv").exists()
    assert (project_root / "step6_report" / "project_report.txt").exists()
    assert (project_root / "step7_validate" / "validation_status.json").exists()
    assert (project_root / "step_index.md").exists()


def test_rebuild_step_views_can_limit_steps(monkeypatch, tmp_path: Path) -> None:
    config_path, project_root, th1_dir, beta_dir = _make_config(tmp_path)
    _seed_project_root(project_root)

    module = _load_script_module()
    monkeypatch.setattr(module, "REPO_ROOT", tmp_path)
    monkeypatch.setattr(
        module,
        "PPI_TARGETS",
        [
            {"name": "TH1", "config_ini": "config/ppi_prod_TH1.ini", "docking_dir": th1_dir},
            {
                "name": "beta-meander",
                "config_ini": "config/ppi_prod_beta_meander.ini",
                "docking_dir": beta_dir,
            },
        ],
    )
    monkeypatch.setattr(
        sys,
        "argv",
        ["rebuild_step_views.py", "--config", str(config_path), "--steps", "4,5,6"],
    )

    exit_code = module.main()

    manifest = _read_json(project_root / "current_run_manifest.json")
    assert exit_code == 0
    assert manifest["execution_mode"] == "backfill:4,5,6"
    assert manifest["step_status"]["step4"] == "complete"
    assert manifest["step_status"]["step5"] == "complete"
    assert manifest["step_status"]["step6"] == "complete"
    assert manifest["step_status"]["step1"] == "not_generated"
    assert not (project_root / "step1_vina_raw").exists()
    assert not (project_root / "step2_ppi_raw").exists()
    assert (project_root / "step4_vina_postprocess").exists()
    assert (project_root / "step5_verdict").exists()
    assert (project_root / "step6_report").exists()
