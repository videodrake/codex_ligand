"""Focused mode-compatibility tests for the additive step view layer."""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import run_production
from egfr_pipeline.step_view import (
    record_step4_outputs,
    record_step5_outputs,
    record_step6_outputs,
    record_step7_outputs,
)
from egfr_pipeline.validate import ValidationResult


def _write_json(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2), encoding="utf-8")


def _write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def _read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _make_config(tmp_path: Path) -> tuple[Path, Path]:
    output_root = tmp_path / "output"
    config = {
        "project_name": "step_modes_project",
        "output_root": str(output_root),
        "mode": "blind",
        "vina": {"mode": "blind", "exhaustiveness": 384, "n_poses": 100},
        "receptors": [{"id": "3GT8_raw"}],
        "ligands": [{"id": "173940", "pdbqt": "input/ligands/173940_ligand.pdbqt"}],
        "ppi": {"pyrosetta_result_dirs": {}},
    }
    config_path = tmp_path / "config.json"
    _write_json(config_path, config)
    return config_path, output_root / config["project_name"]


def _seed_phase456_outputs(project_root: Path, *, token: str = "seed") -> None:
    _write_text(
        project_root / "vina_pose_table.csv",
        f"receptor_id,ligand_id,token\n3GT8_raw,173940,{token}\n",
    )
    _write_text(
        project_root / "vina_pocket_table.csv",
        f"receptor_id,pocket_id,token\n3GT8_raw,3GT8_raw_PKT01,{token}\n",
    )
    _write_text(
        project_root / "vina_drug_pocket_map.csv",
        f"receptor_id,ligand_id,dominant_pocket_id\n3GT8_raw,173940,3GT8_raw_PKT01_{token}\n",
    )
    _write_text(
        project_root / "vina_pocket_comparison.csv",
        f"receptor_a,pocket_a,receptor_b,pocket_b\n3GT8_raw,3GT8_raw_PKT01_{token},3GT8_raw,3GT8_raw_PKT01_{token}\n",
    )
    _write_text(
        project_root / "vina_pocket_bootstrap.csv",
        f"receptor_id,pocket_id,pocket_exists_frac\n3GT8_raw,3GT8_raw_PKT01_{token},0.9\n",
    )
    _write_text(
        project_root / "vina_postprocess_coverage.csv",
        "receptor_id,ligand_id,status\n3GT8_raw,173940,parsed\n",
    )
    _write_text(
        project_root / "valid_sites.csv",
        f"receptor_id,pocket_id,verdict\n3GT8_raw,3GT8_raw_PKT01_{token},STRONG\n",
    )
    _write_text(
        project_root / "cross_method_agreement.csv",
        f"receptor_id,pocket_id,agreement_level\n3GT8_raw,3GT8_raw_PKT01_{token},HIGH\n",
    )
    _write_text(
        project_root / "combined_residue_evidence.csv",
        f"receptor_id,residue_id,evidence\n3GT8_raw,A:ALA699,{token}\n",
    )
    _write_text(
        project_root / "project_report.txt",
        (f"report {token} " * 20).strip() + "\n",
    )


def _seed_minimal_step(project_root: Path, folder_name: str, *, sentinel: str) -> None:
    step_dir = project_root / folder_name
    step_dir.mkdir(parents=True, exist_ok=True)
    _write_json(
        step_dir / "step_manifest.json",
        {
            "status": "complete",
            "generated_at": "2026-03-13T00:00:00Z",
        },
    )
    _write_text(step_dir / "keep.txt", sentinel)


def _seed_existing_step_views(config_path: Path, repo_root: Path) -> None:
    record_step4_outputs(config_path, repo_root=repo_root)
    record_step5_outputs(config_path, repo_root=repo_root)
    record_step6_outputs(config_path, repo_root=repo_root)

    validation_result = ValidationResult()
    validation_result.ok("preexisting validation output")
    record_step7_outputs(
        config_path,
        repo_root=repo_root,
        validation_result=validation_result,
    )


def _phase4_writer(project_root: Path, call_log: list[str], *, token: str):
    def _phase4() -> None:
        call_log.append("phase4")
        _seed_phase456_outputs(project_root, token=token)

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

    return _phase5


def _phase6_writer(project_root: Path, call_log: list[str], *, token: str):
    def _phase6() -> None:
        call_log.append("phase6")
        _write_text(
            project_root / "project_report.txt",
            (f"phase6 {token} " * 20).strip() + "\n",
        )
        _write_text(
            project_root / "combined_residue_evidence.csv",
            f"receptor_id,residue_id,evidence\n3GT8_raw,A:ALA699,{token}\n",
        )

    return _phase6


def _phase7_result(call_log: list[str], *, warning: str | None = None):
    def _phase7() -> ValidationResult:
        call_log.append("phase7")
        result = ValidationResult()
        result.ok("validation ok")
        if warning:
            result.warn(warning)
        return result

    return _phase7


def _phase_check(path: Path):
    def _check() -> list[str]:
        return [] if path.exists() else [path.name]

    return _check


def _install_run_context(
    monkeypatch,
    *,
    config_path: Path,
    repo_root: Path,
    phases: list[tuple[int, str, object]],
    checks: dict[int, tuple[str, object]],
) -> None:
    monkeypatch.setattr(run_production, "CONFIG_PATH", config_path)
    monkeypatch.setattr(run_production, "REPO_ROOT", repo_root)
    monkeypatch.setattr(run_production, "PPI_TARGETS", [])
    monkeypatch.setattr(run_production, "PHASES", phases)
    monkeypatch.setattr(run_production, "PHASE_CHECKS", checks)
    monkeypatch.setattr(run_production, "banner", lambda *_args, **_kwargs: None)


def _run_main(monkeypatch, *args: str) -> None:
    monkeypatch.setattr(sys, "argv", ["run_production.py", *args])
    run_production.main()


def test_status_mode_does_not_write_step_views(monkeypatch, tmp_path: Path) -> None:
    config_path, project_root = _make_config(tmp_path)
    _seed_phase456_outputs(project_root)
    record_step4_outputs(config_path, repo_root=tmp_path)

    step4_manifest = project_root / "step4_vina_postprocess" / "step_manifest.json"
    original_manifest_text = step4_manifest.read_text(encoding="utf-8")

    _install_run_context(
        monkeypatch,
        config_path=config_path,
        repo_root=tmp_path,
        phases=[],
        checks={
            4: ("step4", _phase_check(project_root / "vina_pocket_table.csv")),
            5: ("step5", _phase_check(project_root / "valid_sites.csv")),
            6: ("step6", _phase_check(project_root / "project_report.txt")),
        },
    )

    _run_main(monkeypatch, "--status")

    assert step4_manifest.read_text(encoding="utf-8") == original_manifest_text
    assert not (project_root / "current_run_manifest.json").exists()
    assert not (project_root / "step_index.md").exists()


def test_from_mode_marks_downstream_steps_stale_without_touching_earlier_steps(
    monkeypatch,
    tmp_path: Path,
) -> None:
    config_path, project_root = _make_config(tmp_path)
    _seed_phase456_outputs(project_root, token="before")
    _seed_minimal_step(project_root, "step1_vina_raw", sentinel="keep-step1")
    _seed_minimal_step(project_root, "step2_ppi_raw", sentinel="keep-step2")
    _seed_minimal_step(project_root, "step3_ppi_postprocess", sentinel="keep-step3")
    _seed_existing_step_views(config_path, tmp_path)

    call_log: list[str] = []
    _install_run_context(
        monkeypatch,
        config_path=config_path,
        repo_root=tmp_path,
        phases=[
            (4, "Phase 4", _phase4_writer(project_root, call_log, token="from4")),
            (5, "Phase 5", _phase5_writer(project_root, call_log, token="from5")),
            (6, "Phase 6", _phase6_writer(project_root, call_log, token="from6")),
            (7, "Phase 7", _phase7_result(call_log)),
        ],
        checks={
            4: ("step4", _phase_check(project_root / "vina_pocket_table.csv")),
            5: ("step5", _phase_check(project_root / "valid_sites.csv")),
            6: ("step6", _phase_check(project_root / "project_report.txt")),
        },
    )

    _run_main(monkeypatch, "--from", "4")

    assert call_log == ["phase4", "phase7"]
    assert (project_root / "step1_vina_raw" / "keep.txt").read_text(encoding="utf-8") == "keep-step1"

    step4_manifest = _read_json(project_root / "step4_vina_postprocess" / "step_manifest.json")
    assert step4_manifest["status"] == "complete"
    assert "regenerated_at" in step4_manifest

    run_manifest = _read_json(project_root / "current_run_manifest.json")
    assert run_manifest["execution_mode"] == "from:4"
    assert run_manifest["step_status"]["step1"] == "complete"
    assert run_manifest["step_status"]["step2"] == "complete"
    assert run_manifest["step_status"]["step3"] == "complete"
    assert run_manifest["step_status"]["step4"] == "complete"
    assert run_manifest["step_status"]["step5"] == "stale"
    assert run_manifest["step_status"]["step6"] == "stale"
    assert run_manifest["step_status"]["step7"] == "passed"
    assert run_manifest["stale_steps"] == [5, 6]

    index_text = (project_root / "step_index.md").read_text(encoding="utf-8")
    assert "Some derived steps are stale for the current rerun scope." in index_text


def test_only_mode_rebuilds_only_selected_step(monkeypatch, tmp_path: Path) -> None:
    config_path, project_root = _make_config(tmp_path)
    _seed_phase456_outputs(project_root, token="before")
    _seed_existing_step_views(config_path, tmp_path)

    step4_manifest_path = project_root / "step4_vina_postprocess" / "step_manifest.json"
    step5_manifest_path = project_root / "step5_verdict" / "step_manifest.json"
    step6_manifest_path = project_root / "step6_report" / "step_manifest.json"
    step4_before = step4_manifest_path.read_text(encoding="utf-8")
    step6_before = step6_manifest_path.read_text(encoding="utf-8")

    call_log: list[str] = []
    _install_run_context(
        monkeypatch,
        config_path=config_path,
        repo_root=tmp_path,
        phases=[
            (4, "Phase 4", _phase4_writer(project_root, call_log, token="only4")),
            (5, "Phase 5", _phase5_writer(project_root, call_log, token="only5")),
            (6, "Phase 6", _phase6_writer(project_root, call_log, token="only6")),
        ],
        checks={
            4: ("step4", _phase_check(project_root / "vina_pocket_table.csv")),
            5: ("step5", _phase_check(project_root / "valid_sites.csv")),
            6: ("step6", _phase_check(project_root / "project_report.txt")),
        },
    )

    _run_main(monkeypatch, "--only", "5")

    assert call_log == ["phase5"]
    assert step4_manifest_path.read_text(encoding="utf-8") == step4_before
    assert step6_manifest_path.read_text(encoding="utf-8") == step6_before

    step5_manifest = _read_json(step5_manifest_path)
    assert step5_manifest["status"] == "complete"
    assert "regenerated_at" in step5_manifest

    run_manifest = _read_json(project_root / "current_run_manifest.json")
    assert run_manifest["execution_mode"] == "only:5"
    assert run_manifest["step_status"]["step4"] == "complete"
    assert run_manifest["step_status"]["step5"] == "complete"
    assert run_manifest["step_status"]["step6"] == "complete"
    assert run_manifest["stale_steps"] == []


def test_force_mode_rebuilds_even_when_outputs_exist(monkeypatch, tmp_path: Path) -> None:
    config_path, project_root = _make_config(tmp_path)
    _seed_phase456_outputs(project_root, token="before")
    record_step5_outputs(config_path, repo_root=tmp_path)

    step5_manifest_path = project_root / "step5_verdict" / "step_manifest.json"
    call_log: list[str] = []
    phases = [(5, "Phase 5", _phase5_writer(project_root, call_log, token="forced"))]
    checks = {5: ("step5", _phase_check(project_root / "valid_sites.csv"))}

    _install_run_context(
        monkeypatch,
        config_path=config_path,
        repo_root=tmp_path,
        phases=phases,
        checks=checks,
    )
    _run_main(monkeypatch)

    step5_manifest = _read_json(step5_manifest_path)
    assert call_log == []
    assert "regenerated_at" not in step5_manifest

    _install_run_context(
        monkeypatch,
        config_path=config_path,
        repo_root=tmp_path,
        phases=phases,
        checks=checks,
    )
    _run_main(monkeypatch, "--force")

    step5_manifest = _read_json(step5_manifest_path)
    assert call_log == ["phase5"]
    assert step5_manifest["status"] == "complete"
    assert "regenerated_at" in step5_manifest


def test_fresh_run_creates_new_manifest_for_selected_step(monkeypatch, tmp_path: Path) -> None:
    config_path, project_root = _make_config(tmp_path)
    call_log: list[str] = []

    _install_run_context(
        monkeypatch,
        config_path=config_path,
        repo_root=tmp_path,
        phases=[(5, "Phase 5", _phase5_writer(project_root, call_log, token="fresh"))],
        checks={5: ("step5", _phase_check(project_root / "valid_sites.csv"))},
    )

    _run_main(monkeypatch, "--only", "5")

    assert call_log == ["phase5"]
    assert (project_root / "step5_verdict" / "step_manifest.json").exists()

    run_manifest = _read_json(project_root / "current_run_manifest.json")
    assert run_manifest["fresh_run"] is True
    assert run_manifest["execution_mode"] == "only:5"
    assert run_manifest["step_status"]["step5"] == "complete"


def test_reset_cleanup_targets_include_step_view_artifacts(monkeypatch, tmp_path: Path) -> None:
    script_path = Path(__file__).resolve().parents[1] / "scripts" / "reset_production_outputs.py"
    spec = importlib.util.spec_from_file_location("reset_production_outputs_test", script_path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)

    project_root = tmp_path / "output" / "step_modes_project"
    monkeypatch.setattr(module, "_project_output_dir", lambda: project_root)
    monkeypatch.setattr(module, "_ppi_output_dirs", lambda: set())
    monkeypatch.setattr(module, "PHASE_OUTPUT_DIRS", ())

    targets = module.collect_cleanup_targets()

    assert project_root / "step1_vina_raw" in targets
    assert project_root / "step7_validate" in targets
    assert project_root / "current_run_manifest.json" in targets
    assert project_root / "step_index.md" in targets
    assert project_root in targets


def test_step_view_docs_cover_reading_order_and_canonical_vs_derived_contract() -> None:
    repo_root = Path(__file__).resolve().parents[1]
    expectations = {
        "README.md": [
            "step_index.md",
            "canonical runtime outputs",
            "derived interpretation views",
        ],
        "output/README.md": [
            "step_index.md",
            "canonical",
            "derived interpretation views",
        ],
        "docs/data_flow_guide.md": [
            "step_index.md",
            "canonical runtime outputs",
            "derived interpretation view",
            "step1_vina_raw/",
            "step7_validate/",
        ],
        "docs/output_artifact_map.md": [
            "step_index.md",
            "canonical",
            "derived interpretation view",
            "step1_vina_raw/",
            "step7_validate/",
        ],
        "docs/current_pipeline_status.md": [
            "step_index.md",
            "canonical runtime outputs",
            "derived view",
        ],
        "docs/architecture.md": [
            "step_index.md",
            "canonical",
            "derived interpretation view",
            "step1_vina_raw/",
        ],
        "docs/runbook.md": [
            "step_index.md",
            "--status",
            "--from",
            "--only",
            "--force",
            "reset_production_outputs.py",
            "derived interpretation view",
        ],
    }

    for relative_path, snippets in expectations.items():
        text = (repo_root / relative_path).read_text(encoding="utf-8").lower()
        for snippet in snippets:
            assert snippet.lower() in text, f"Missing '{snippet}' in {relative_path}"
