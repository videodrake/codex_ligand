"""Synthetic end-to-end tests for the additive step output layer."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import run_production
from egfr_pipeline.output_steps import STEP_SPECS
from egfr_pipeline.validate import ValidationResult


RECEPTOR_IDS = ["3GT8_raw", "EGFR_160-185", "EGFR_170-200"]
LIGAND_NAMES = ["173940", "97806", "VAX-C12_0"]


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
        "project_name": "step_e2e_project",
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


def _ppi_targets(th1_dir: Path, beta_dir: Path) -> list[dict]:
    return [
        {
            "name": "TH1",
            "config_ini": "config/ppi_prod_TH1.ini",
            "docking_dir": th1_dir,
        },
        {
            "name": "beta-meander",
            "config_ini": "config/ppi_prod_beta_meander.ini",
            "docking_dir": beta_dir,
        },
    ]


def _phase1_writer(project_root: Path, call_log: list[str]):
    def _phase1() -> None:
        call_log.append("phase1")
        for receptor_id in RECEPTOR_IDS:
            for ligand_name in LIGAND_NAMES:
                _write_text(
                    project_root / receptor_id / f"{ligand_name}_blind.pdbqt",
                    "MODEL 1\nENDMDL\nMODEL 2\nENDMDL\n",
                )

    return _phase1


def _phase2_writer(th1_dir: Path, beta_dir: Path, call_log: list[str]):
    def _phase2() -> None:
        call_log.append("phase2")
        _write_text(
            th1_dir / "final_result" / "final_ranking.csv",
            "rank,score\n1,-20.5\n",
        )
        _write_json(
            th1_dir / "pyrosetta_run_metadata.json",
            {"status": "complete", "n_models": 20000, "target": "TH1"},
        )
        _write_text(
            beta_dir / "final_result" / "final_ranking.csv",
            "rank,score\n1,-19.8\n",
        )
        _write_json(
            beta_dir / "pyrosetta_run_metadata.json",
            {"status": "complete", "n_models": 20000, "target": "beta-meander"},
        )

    return _phase2


def _phase3_writer(project_root: Path, repo_root: Path, call_log: list[str]):
    def _phase3() -> None:
        call_log.append("phase3")
        _write_text(
            project_root / "ppi_pyrosetta_residues.csv",
            "receptor_id,residue_id,occupancy\n3GT8_raw,A:ALA699,0.8\n",
        )
        _write_text(
            project_root / "ppi_pyrosetta_summary.csv",
            "receptor_id,n_final_models\n3GT8_raw,20\n",
        )
        _write_text(
            repo_root / "output" / "phase1_ppi" / "phase1_interface_report.md",
            "# Phase 1 Interface Report\n",
        )

    return _phase3


def _phase4_writer(project_root: Path, call_log: list[str], *, token: str):
    def _phase4() -> None:
        call_log.append("phase4")
        _write_text(
            project_root / "vina_pose_table.csv",
            f"receptor_id,ligand_id,token\n3GT8_raw,173940,{token}\n",
        )
        _write_text(
            project_root / "vina_pocket_table.csv",
            f"receptor_id,pocket_id,token\n3GT8_raw,3GT8_raw_PKT01_{token},{token}\n",
        )
        _write_text(
            project_root / "vina_drug_pocket_map.csv",
            f"receptor_id,ligand_id,dominant_pocket_id\n3GT8_raw,173940,3GT8_raw_PKT01_{token}\n",
        )
        _write_text(
            project_root / "vina_pocket_comparison.csv",
            f"receptor_a,pocket_a,receptor_b,pocket_b\n3GT8_raw,3GT8_raw_PKT01_{token},EGFR_160-185,EGFR_160-185_PKT01_{token}\n",
        )
        _write_text(
            project_root / "vina_pocket_bootstrap.csv",
            f"receptor_id,pocket_id,pocket_exists_frac\n3GT8_raw,3GT8_raw_PKT01_{token},0.9\n",
        )

    return _phase4


def _phase5_writer(project_root: Path, call_log: list[str], *, token: str):
    def _phase5() -> None:
        call_log.append("phase5")
        _write_text(
            project_root / "valid_sites.csv",
            f"receptor_id,pocket_id,verdict\n3GT8_raw,3GT8_raw_PKT01_{token},STRONG\n",
        )
        _write_text(
            project_root / "cross_method_agreement.csv",
            f"receptor_id,pocket_id,agreement_level\n3GT8_raw,3GT8_raw_PKT01_{token},HIGH\n",
        )
        _write_text(
            project_root / "vina_consensus_sites.csv",
            f"consensus_site_id,n_receptors\nCS_{token},3\n",
        )

    return _phase5


def _phase6_writer(project_root: Path, call_log: list[str], *, token: str):
    def _phase6() -> None:
        call_log.append("phase6")
        _write_text(
            project_root / "project_report.txt",
            (f"project report {token} " * 20).strip() + "\n",
        )
        _write_text(
            project_root / "combined_residue_evidence.csv",
            f"receptor_id,residue_id,evidence\n3GT8_raw,A:ALA699,{token}\n",
        )

    return _phase6


def _phase7_writer(call_log: list[str], *, token: str):
    def _phase7() -> ValidationResult:
        call_log.append("phase7")
        result = ValidationResult()
        result.ok(f"validation ok {token}")
        return result

    return _phase7


def _phase1_check(project_root: Path):
    def _check() -> list[str]:
        missing = []
        for receptor_id in RECEPTOR_IDS:
            for ligand_name in LIGAND_NAMES:
                pose_path = project_root / receptor_id / f"{ligand_name}_blind.pdbqt"
                if not pose_path.exists():
                    missing.append(f"{receptor_id}/{ligand_name}_blind.pdbqt")
        return missing

    return _check


def _paths_exist_check(*paths: Path):
    def _check() -> list[str]:
        return [path.name for path in paths if not path.exists()]

    return _check


def _report_check(path: Path):
    def _check() -> list[str]:
        if not path.exists() or path.stat().st_size < 100:
            return [path.name]
        return []

    return _check


def _install_run_context(
    monkeypatch,
    *,
    config_path: Path,
    repo_root: Path,
    targets: list[dict],
    phases: list[tuple[int, str, object]],
    checks: dict[int, tuple[str, object]],
) -> None:
    monkeypatch.setattr(run_production, "CONFIG_PATH", config_path)
    monkeypatch.setattr(run_production, "REPO_ROOT", repo_root)
    monkeypatch.setattr(run_production, "PPI_TARGETS", targets)
    monkeypatch.setattr(run_production, "PHASES", phases)
    monkeypatch.setattr(run_production, "PHASE_CHECKS", checks)
    monkeypatch.setattr(run_production, "banner", lambda *_args, **_kwargs: None)


def _run_main(monkeypatch, *args: str) -> None:
    monkeypatch.setattr(sys, "argv", ["run_production.py", *args])
    run_production.main()


def _build_phases(
    *,
    project_root: Path,
    repo_root: Path,
    th1_dir: Path,
    beta_dir: Path,
    call_log: list[str],
    phase4_token: str = "phase4",
    phase5_token: str = "phase5",
    phase6_token: str = "phase6",
    phase7_token: str = "phase7",
) -> list[tuple[int, str, object]]:
    return [
        (1, "Phase 1", _phase1_writer(project_root, call_log)),
        (2, "Phase 2", _phase2_writer(th1_dir, beta_dir, call_log)),
        (3, "Phase 3", _phase3_writer(project_root, repo_root, call_log)),
        (4, "Phase 4", _phase4_writer(project_root, call_log, token=phase4_token)),
        (5, "Phase 5", _phase5_writer(project_root, call_log, token=phase5_token)),
        (6, "Phase 6", _phase6_writer(project_root, call_log, token=phase6_token)),
        (7, "Phase 7", _phase7_writer(call_log, token=phase7_token)),
    ]


def _build_checks(project_root: Path, th1_dir: Path, beta_dir: Path) -> dict[int, tuple[str, object]]:
    return {
        1: ("step1", _phase1_check(project_root)),
        2: (
            "step2",
            _paths_exist_check(
                th1_dir / "final_result" / "final_ranking.csv",
                beta_dir / "final_result" / "final_ranking.csv",
            ),
        ),
        3: (
            "step3",
            _paths_exist_check(
                project_root / "ppi_pyrosetta_residues.csv",
                project_root / "ppi_pyrosetta_summary.csv",
            ),
        ),
        4: (
            "step4",
            _paths_exist_check(
                project_root / "vina_pose_table.csv",
                project_root / "vina_pocket_table.csv",
                project_root / "vina_drug_pocket_map.csv",
                project_root / "vina_pocket_comparison.csv",
                project_root / "vina_pocket_bootstrap.csv",
            ),
        ),
        5: (
            "step5",
            _paths_exist_check(
                project_root / "valid_sites.csv",
                project_root / "cross_method_agreement.csv",
            ),
        ),
        6: ("step6", _report_check(project_root / "project_report.txt")),
    }


def _run_full_pipeline(monkeypatch, tmp_path: Path) -> tuple[Path, Path, Path, Path, list[str]]:
    config_path, project_root, th1_dir, beta_dir = _make_config(tmp_path)
    call_log: list[str] = []
    _install_run_context(
        monkeypatch,
        config_path=config_path,
        repo_root=tmp_path,
        targets=_ppi_targets(th1_dir, beta_dir),
        phases=_build_phases(
            project_root=project_root,
            repo_root=tmp_path,
            th1_dir=th1_dir,
            beta_dir=beta_dir,
            call_log=call_log,
        ),
        checks=_build_checks(project_root, th1_dir, beta_dir),
    )
    _run_main(monkeypatch)
    return config_path, project_root, th1_dir, beta_dir, call_log


def test_full_fresh_run_e2e_creates_step1_to_step7_additively(
    monkeypatch,
    tmp_path: Path,
) -> None:
    config_path, project_root, th1_dir, beta_dir, call_log = _run_full_pipeline(monkeypatch, tmp_path)

    assert call_log == [
        "phase1",
        "phase2",
        "phase3",
        "phase4",
        "phase5",
        "phase6",
        "phase7",
    ]
    assert config_path.exists()
    assert project_root.exists()

    for step_num, spec in STEP_SPECS.items():
        step_dir = project_root / spec.folder_name
        assert step_dir.exists(), spec.folder_name
        assert (step_dir / "step_manifest.json").exists()
        assert (step_dir / "summary.md").exists()
        for artifact_name in spec.required_artifacts:
            assert (step_dir / artifact_name).exists(), f"Missing {artifact_name} in {spec.folder_name}"

    run_manifest = _read_json(project_root / "current_run_manifest.json")
    assert run_manifest["execution_mode"] == "full"
    assert run_manifest["fresh_run"] is True
    assert run_manifest["stale_steps"] == []
    assert run_manifest["step_status"] == {
        "step1": "complete",
        "step2": "complete",
        "step3": "complete",
        "step4": "complete",
        "step5": "complete",
        "step6": "complete",
        "step7": "passed",
    }

    index_text = (project_root / "step_index.md").read_text(encoding="utf-8")
    assert "## Step Overview Table" in index_text
    assert "## Where To Read First" in index_text
    assert "`step7_validate/validation_status.json`" in index_text

    assert (project_root / "step1_vina_raw" / "raw_pose_index.csv").read_text(encoding="utf-8").count("\n") == 10
    assert (project_root / "step2_ppi_raw" / "raw_run_paths.tsv").exists()
    assert not (project_root / "step2_ppi_raw" / "final_result").exists()
    assert (th1_dir / "final_result" / "final_ranking.csv").exists()
    assert (beta_dir / "final_result" / "final_ranking.csv").exists()

    for receptor_id in RECEPTOR_IDS:
        for ligand_name in LIGAND_NAMES:
            assert (project_root / receptor_id / f"{ligand_name}_blind.pdbqt").exists()

    for name in (
        "vina_pose_table.csv",
        "vina_pocket_table.csv",
        "vina_drug_pocket_map.csv",
        "ppi_pyrosetta_residues.csv",
        "ppi_pyrosetta_summary.csv",
        "valid_sites.csv",
        "cross_method_agreement.csv",
        "combined_residue_evidence.csv",
        "project_report.txt",
    ):
        assert (project_root / name).exists(), name


def test_partial_rerun_e2e_keeps_step1_to_step4_and_marks_unrebuilt_step6_stale(
    monkeypatch,
    tmp_path: Path,
) -> None:
    config_path, project_root, th1_dir, beta_dir, _ = _run_full_pipeline(monkeypatch, tmp_path)

    step14_before = {
        step_num: (project_root / STEP_SPECS[step_num].folder_name / "step_manifest.json").read_text(encoding="utf-8")
        for step_num in (1, 2, 3, 4)
    }
    step6_before = (project_root / "step6_report" / "step_manifest.json").read_text(encoding="utf-8")
    pocket_table_before = (project_root / "vina_pocket_table.csv").read_text(encoding="utf-8")
    valid_sites_before = (project_root / "valid_sites.csv").read_text(encoding="utf-8")

    rerun_log: list[str] = []
    _install_run_context(
        monkeypatch,
        config_path=config_path,
        repo_root=tmp_path,
        targets=_ppi_targets(th1_dir, beta_dir),
        phases=_build_phases(
            project_root=project_root,
            repo_root=tmp_path,
            th1_dir=th1_dir,
            beta_dir=beta_dir,
            call_log=rerun_log,
            phase4_token="rerun4",
            phase5_token="rerun5",
            phase6_token="rerun6",
            phase7_token="rerun7",
        ),
        checks=_build_checks(project_root, th1_dir, beta_dir),
    )
    _run_main(monkeypatch, "--from", "5")

    assert rerun_log == ["phase5", "phase7"]
    for step_num, original_text in step14_before.items():
        manifest_path = project_root / STEP_SPECS[step_num].folder_name / "step_manifest.json"
        assert manifest_path.read_text(encoding="utf-8") == original_text

    assert (project_root / "step6_report" / "step_manifest.json").read_text(encoding="utf-8") == step6_before
    assert (project_root / "vina_pocket_table.csv").read_text(encoding="utf-8") == pocket_table_before
    assert (project_root / "valid_sites.csv").read_text(encoding="utf-8") != valid_sites_before

    step5_manifest = _read_json(project_root / "step5_verdict" / "step_manifest.json")
    assert step5_manifest["status"] == "complete"
    assert "regenerated_at" in step5_manifest

    run_manifest = _read_json(project_root / "current_run_manifest.json")
    assert run_manifest["execution_mode"] == "from:5"
    assert run_manifest["step_status"]["step1"] == "complete"
    assert run_manifest["step_status"]["step2"] == "complete"
    assert run_manifest["step_status"]["step3"] == "complete"
    assert run_manifest["step_status"]["step4"] == "complete"
    assert run_manifest["step_status"]["step5"] == "complete"
    assert run_manifest["step_status"]["step6"] == "stale"
    assert run_manifest["step_status"]["step7"] == "passed"
    assert run_manifest["stale_steps"] == [6]


def test_status_only_e2e_does_not_mutate_existing_step_outputs(
    monkeypatch,
    tmp_path: Path,
) -> None:
    config_path, project_root, th1_dir, beta_dir, _ = _run_full_pipeline(monkeypatch, tmp_path)

    tracked_files = [
        project_root / "step1_vina_raw" / "step_manifest.json",
        project_root / "step5_verdict" / "step_manifest.json",
        project_root / "step7_validate" / "validation_status.json",
        project_root / "current_run_manifest.json",
        project_root / "step_index.md",
    ]
    snapshots = {
        path: (path.read_text(encoding="utf-8"), path.stat().st_mtime_ns)
        for path in tracked_files
    }

    _install_run_context(
        monkeypatch,
        config_path=config_path,
        repo_root=tmp_path,
        targets=_ppi_targets(th1_dir, beta_dir),
        phases=_build_phases(
            project_root=project_root,
            repo_root=tmp_path,
            th1_dir=th1_dir,
            beta_dir=beta_dir,
            call_log=[],
        ),
        checks=_build_checks(project_root, th1_dir, beta_dir),
    )
    _run_main(monkeypatch, "--status")

    for path, (text_before, mtime_before) in snapshots.items():
        assert path.read_text(encoding="utf-8") == text_before
        assert path.stat().st_mtime_ns == mtime_before


def test_root_output_compatibility_e2e_keeps_canonical_files_in_place(
    monkeypatch,
    tmp_path: Path,
) -> None:
    _, project_root, th1_dir, beta_dir, _ = _run_full_pipeline(monkeypatch, tmp_path)

    canonical_root_files = [
        project_root / "vina_pose_table.csv",
        project_root / "vina_pocket_table.csv",
        project_root / "vina_drug_pocket_map.csv",
        project_root / "vina_pocket_comparison.csv",
        project_root / "vina_pocket_bootstrap.csv",
        project_root / "ppi_pyrosetta_residues.csv",
        project_root / "ppi_pyrosetta_summary.csv",
        project_root / "valid_sites.csv",
        project_root / "cross_method_agreement.csv",
        project_root / "combined_residue_evidence.csv",
        project_root / "project_report.txt",
    ]
    for path in canonical_root_files:
        assert path.exists(), path.name

    assert not any(path.suffix == ".pdbqt" for path in (project_root / "step1_vina_raw").iterdir())
    assert not (project_root / "step2_ppi_raw" / th1_dir.name).exists()
    assert not (project_root / "step2_ppi_raw" / beta_dir.name).exists()
    assert (th1_dir / "final_result" / "final_ranking.csv").exists()
    assert (beta_dir / "final_result" / "final_ranking.csv").exists()

    raw_run_text = (project_root / "step2_ppi_raw" / "raw_run_paths.tsv").read_text(encoding="utf-8")
    assert "ppi_raw/TH1_run" in raw_run_text
    assert "ppi_raw/beta_run" in raw_run_text

    status = _read_json(project_root / "step7_validate" / "validation_status.json")
    assert status["project_root"] == "output/step_e2e_project"
