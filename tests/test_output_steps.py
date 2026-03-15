"""Tests for the additive step output view helpers."""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import run_production
from egfr_pipeline.output_steps import (
    STEP_SPECS,
    build_current_run_manifest,
    copy_artifact_if_exists,
    ensure_step_dir,
    record_step1_outputs,
    record_step2_outputs,
    record_step3_outputs,
    record_step4_outputs,
    record_step5_outputs,
    record_step6_outputs,
    record_step7_outputs,
    refresh_root_step_views,
    step_output_view_enabled,
    write_artifact_index,
    write_step_manifest,
    write_step_summary,
)
from egfr_pipeline.validate import ValidationResult


def _write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def _write_json(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2), encoding="utf-8")


def _make_config(tmp_path: Path, *, step_view_enabled_flag: bool = True) -> Tuple[Path, Path]:
    config = {
        "project_name": "step_view_project",
        "output_root": "./output",
        "mode": "blind",
        "vina": {"mode": "blind", "exhaustiveness": 384, "n_poses": 100},
        "step_output_view": {"enabled": step_view_enabled_flag},
        "receptors": [
            {"id": "3GT8_raw"},
            {"id": "EGFR_160-185"},
            {"id": "EGFR_170-200"},
        ],
        "ligands": [
            {"id": "173940", "pdbqt": "input/ligands/173940_ligand.pdbqt"},
            {"id": "97806", "pdbqt": "input/ligands/97806_ligand.pdbqt"},
            {"id": "VAX-C12_0", "pdbqt": "input/ligands/VAX-C12_0_ligand.pdbqt"},
        ],
        "ppi": {
            "pyrosetta_result_dirs": {
                "3GT8_raw": [
                    {"partner": "MYO1D_TH1", "path": "/ppi/raw/TH1"},
                    {"partner": "MYO1D_beta", "path": "/ppi/raw/beta"},
                ]
            }
        },
    }
    config_path = tmp_path / "config.json"
    _write_json(config_path, config)
    return config_path, tmp_path / "output" / config["project_name"]


def _canonical_artifact_contents() -> Dict[str, str]:
    return {
        "vina_pose_table.csv": "receptor_id,ligand_id\n3GT8_raw,173940\n",
        "vina_pocket_table.csv": "receptor_id,pocket_id\n3GT8_raw,3GT8_raw_PKT01\n",
        "vina_drug_pocket_map.csv": "receptor_id,ligand_id,dominant_pocket_id\n3GT8_raw,173940,3GT8_raw_PKT01\n",
        "vina_pocket_comparison.csv": "receptor_a,pocket_a,receptor_b,pocket_b\n3GT8_raw,3GT8_raw_PKT01,EGFR_160-185,EGFR_160-185_PKT01\n",
        "vina_pocket_bootstrap.csv": "receptor_id,pocket_id,pocket_exists_frac\n3GT8_raw,3GT8_raw_PKT01,0.9\n",
        "valid_sites.csv": "receptor_id,pocket_id,verdict\n3GT8_raw,3GT8_raw_PKT01,STRONG\n",
        "cross_method_agreement.csv": "receptor_id,pocket_id,agreement_level\n3GT8_raw,3GT8_raw_PKT01,HIGH\n",
        "project_report.txt": "Step 6 narrative report.\n",
        "combined_residue_evidence.csv": "receptor_id,residue_id,evidence\n3GT8_raw,A:ALA699,vina\n",
    }


def _seed_project_root(project_root: Path, *, include_optional_step5: bool = False) -> Dict[str, str]:
    contents = _canonical_artifact_contents()
    for name, text in contents.items():
        _write_text(project_root / name, text)
    if include_optional_step5:
        _write_text(
            project_root / "vina_consensus_sites.csv",
            "consensus_site_id,n_receptors\nCS01,3\n",
        )
    (project_root / "3GT8_raw").mkdir(parents=True, exist_ok=True)
    _write_text(project_root / "3GT8_raw" / "173940_blind.pdbqt", "MODEL 1\nENDMDL\n")
    return contents


def _seed_step1_pose_tree(project_root: Path, *, missing_pairs: Tuple = ()) -> None:
    receptors = ["3GT8_raw", "EGFR_160-185", "EGFR_170-200"]
    ligands = ["173940", "97806", "VAX-C12_0"]
    for receptor_id in receptors:
        for ligand_id in ligands:
            if (receptor_id, ligand_id) in missing_pairs:
                continue
            pose_path = project_root / receptor_id / f"{ligand_id}_blind.pdbqt"
            _write_text(pose_path, "MODEL 1\nENDMDL\nMODEL 2\nENDMDL\n")


def _seed_step2_run_dir(
    run_dir: Path,
    *,
    include_ranking: bool = True,
    include_metadata: bool = False,
) -> None:
    if include_ranking:
        _write_text(
            run_dir / "final_result" / "final_ranking.csv",
            "rank,score\n1,-20.5\n",
        )
    _write_text(run_dir / "config_snapshot.ini", "[Metadata]\npartner=demo\n")
    if include_metadata:
        _write_json(
            run_dir / "pyrosetta_run_metadata.json",
            {"status": "complete", "n_models": 20000},
        )


def _seed_step3_outputs(
    project_root: Path,
    repo_root: Path,
    *,
    use_fallback_ppi_dir: bool = False,
    include_phase1_report: bool = True,
) -> None:
    base = project_root / "ppi" if use_fallback_ppi_dir else project_root
    _write_text(
        base / "ppi_pyrosetta_residues.csv",
        "receptor_id,residue_id,occupancy\n3GT8_raw,A:ALA699,0.8\n",
    )
    _write_text(
        base / "ppi_pyrosetta_summary.csv",
        "receptor_id,n_final_models\n3GT8_raw,20\n",
    )
    if include_phase1_report:
        _write_text(
            repo_root / "output" / "phase1_ppi" / "phase1_interface_report.md",
            "# Phase 1 Interface Report\n",
        )


def _read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _phase5_only_writer(project_root: Path, call_log: List[str]) -> object:
    def _phase5() -> None:
        call_log.append("phase5")
        _write_text(
            project_root / "valid_sites.csv",
            "receptor_id,pocket_id,verdict\n3GT8_raw,3GT8_raw_PKT01,STRONG\n",
        )
        _write_text(
            project_root / "cross_method_agreement.csv",
            "receptor_id,pocket_id,agreement_level\n3GT8_raw,3GT8_raw_PKT01,HIGH\n",
        )

    return _phase5


def _phase5_check(project_root: Path) -> object:
    def _check() -> List[str]:
        missing = []
        for name in ("valid_sites.csv", "cross_method_agreement.csv"):
            if not (project_root / name).exists():
                missing.append(name)
        return missing

    return _check


def test_helper_writers_cover_minimum_contracts(tmp_path: Path) -> None:
    project_root = tmp_path / "output" / "helper_project"
    step_dir = ensure_step_dir(project_root, 4)
    assert step_dir == project_root / "step4_vina_postprocess"
    assert step_dir.exists()

    source = tmp_path / "canonical.csv"
    _write_text(source, "a,b\n1,2\n")
    missing: List[str] = []
    copied = copy_artifact_if_exists(source, step_dir / "copied.csv", missing)
    assert copied is True
    assert missing == []
    assert (step_dir / "copied.csv").read_text(encoding="utf-8") == "a,b\n1,2\n"

    copied_missing = copy_artifact_if_exists(
        tmp_path / "absent.csv",
        step_dir / "absent.csv",
        missing,
        "absent.csv",
    )
    assert copied_missing is False
    assert missing == ["absent.csv"]

    index_path = write_artifact_index(
        step_dir,
        "raw_run_paths.tsv",
        rows=[{"target": "TH1", "path": "ppi/TH1"}],
        fieldnames=["target", "path"],
    )
    assert "\t" in index_path.read_text(encoding="utf-8")

    manifest_path = write_step_manifest(
        step_dir,
        {
            "step_number": 4,
            "step_name": "vina_postprocess",
            "phase_number": 4,
            "generated_at": "2026-03-13T00:00:00Z",
            "project_name": "helper_project",
            "project_root": "output/helper_project",
            "source_config": "config.json",
            "receptor_ids": ["3GT8_raw"],
            "ligand_ids": ["173940"],
            "upstream_steps": [1],
            "artifact_paths": ["copied.csv"],
            "notes": "",
        },
    )
    summary_path = write_step_summary(
        step_dir,
        step_num=4,
        description="Pocket-level interpretation view.",
        key_files=[{"path": "copied.csv", "description": "Inspect this first."}],
        next_step_reads=[{"path": "step5_verdict/valid_sites.csv", "description": "Next summary."}],
        warnings=["No warnings."],
    )

    manifest = _read_json(manifest_path)
    assert manifest["step_number"] == 4
    summary_text = summary_path.read_text(encoding="utf-8")
    assert "## Scientific Meaning" in summary_text
    assert "## Key Files" in summary_text
    assert "## Next Step Reads" in summary_text


def test_record_step1_outputs_indexes_raw_pose_paths_without_copying_pose_files(tmp_path: Path) -> None:
    config_path, project_root = _make_config(tmp_path)
    _seed_step1_pose_tree(
        project_root,
        missing_pairs=(("EGFR_170-200", "VAX-C12_0"),),
    )

    step1_dir = record_step1_outputs(config_path, repo_root=tmp_path)
    manifest_path, index_path = refresh_root_step_views(
        config_path,
        repo_root=tmp_path,
        execution_mode="from:1",
        ppi_config_paths=["config/ppi_prod_TH1.ini"],
    )

    assert step1_dir == project_root / "step1_vina_raw"
    assert (step1_dir / "raw_pose_index.csv").exists()
    assert not any(path.suffix == ".pdbqt" for path in step1_dir.iterdir())

    rows = (step1_dir / "raw_pose_index.csv").read_text(encoding="utf-8").strip().splitlines()
    assert rows[0] == "receptor_id,ligand_id,raw_pose_file,n_models,source,docking_mode,exhaustiveness,n_poses"
    assert len(rows) == 10
    assert "3GT8_raw,173940,3GT8_raw/173940_blind.pdbqt,2,canonical_output,blind,384,100" in rows
    assert "EGFR_170-200,VAX-C12_0,,0,missing,blind,384,100" in rows

    manifest = _read_json(step1_dir / "step_manifest.json")
    assert manifest["status"] == "partial"
    assert "EGFR_170-200/VAX-C12_0_blind.pdbqt" in manifest["missing_files"]
    assert manifest["artifact_paths"] == ["raw_pose_index.csv"]
    assert manifest["source_artifacts"][0]["canonical_path"] == "3GT8_raw/173940_blind.pdbqt"

    summary_text = (step1_dir / "summary.md").read_text(encoding="utf-8")
    assert "Indexed 8 of 9 expected docking pairs" in summary_text
    assert "step4_vina_postprocess/vina_pocket_table.csv" in summary_text

    run_manifest = _read_json(manifest_path)
    assert run_manifest["step_status"]["step1"] == "partial"
    assert run_manifest["step_status"]["step4"] == "not_generated"
    index_text = index_path.read_text(encoding="utf-8")
    assert "`step1_vina_raw/raw_pose_index.csv`" in index_text
    assert "Partial rerun detected" in index_text


def test_record_step2_outputs_copies_rankings_and_indexes_raw_run_paths(tmp_path: Path) -> None:
    config_path, project_root = _make_config(tmp_path)
    th1_dir = tmp_path / "ppi_raw" / "TH1_run"
    beta_dir = tmp_path / "ppi_raw" / "beta_run"
    _seed_step2_run_dir(th1_dir, include_ranking=True, include_metadata=True)
    _seed_step2_run_dir(beta_dir, include_ranking=False, include_metadata=False)

    step2_dir = record_step2_outputs(
        config_path,
        repo_root=tmp_path,
        ppi_targets=[
            {"name": "TH1", "docking_dir": th1_dir},
            {"name": "beta-meander", "docking_dir": beta_dir},
        ],
    )
    manifest_path, index_path = refresh_root_step_views(
        config_path,
        repo_root=tmp_path,
        execution_mode="from:2",
        ppi_config_paths=["config/ppi_prod_TH1.ini", "config/ppi_prod_beta_meander.ini"],
    )

    assert step2_dir == project_root / "step2_ppi_raw"
    assert (step2_dir / "TH1_final_ranking.csv").exists()
    assert not (step2_dir / "beta_meander_final_ranking.csv").exists()
    assert (step2_dir / "raw_run_paths.tsv").exists()
    assert (step2_dir / "pyrosetta_run_metadata.json").exists()
    assert not (step2_dir / "final_result").exists()

    raw_run_text = (step2_dir / "raw_run_paths.tsv").read_text(encoding="utf-8")
    assert "TH1\tppi_raw/TH1_run\tppi_raw/TH1_run/final_result/final_ranking.csv\tppi_raw/TH1_run/pyrosetta_run_metadata.json\tcanonical_output" in raw_run_text
    assert "beta-meander\tppi_raw/beta_run\t\t\tmissing" in raw_run_text

    step2_manifest = _read_json(step2_dir / "step_manifest.json")
    assert step2_manifest["status"] == "partial"
    assert step2_manifest["missing_files"] == ["beta_meander_final_ranking.csv"]
    assert any("Raw run directory intentionally referenced, not duplicated" in warning for warning in step2_manifest["warnings"])
    assert any("Optional metadata JSON missing for beta-meander." == warning for warning in step2_manifest["warnings"])

    metadata = _read_json(step2_dir / "pyrosetta_run_metadata.json")
    assert metadata["TH1"]["data"]["n_models"] == 20000

    run_manifest = _read_json(manifest_path)
    assert run_manifest["step_status"]["step2"] == "partial"
    index_text = index_path.read_text(encoding="utf-8")
    assert "`step2_ppi_raw`" in index_text


def test_record_step3_outputs_copies_root_level_ppi_postprocess_artifacts(tmp_path: Path) -> None:
    config_path, project_root = _make_config(tmp_path)
    _seed_step3_outputs(project_root, tmp_path, use_fallback_ppi_dir=False, include_phase1_report=True)

    step3_dir = record_step3_outputs(config_path, repo_root=tmp_path)
    manifest_path, index_path = refresh_root_step_views(
        config_path,
        repo_root=tmp_path,
        execution_mode="from:3",
        ppi_config_paths=["config/ppi_prod_TH1.ini"],
    )

    assert step3_dir == project_root / "step3_ppi_postprocess"
    assert (step3_dir / "ppi_pyrosetta_residues.csv").exists()
    assert (step3_dir / "ppi_pyrosetta_summary.csv").exists()
    assert (step3_dir / "phase1_interface_report.md").exists()

    manifest = _read_json(step3_dir / "step_manifest.json")
    assert manifest["status"] == "complete"
    assert manifest["missing_files"] == []
    assert "phase1_interface_report.md" in manifest["artifact_paths"]
    assert manifest["source_artifacts"][0]["canonical_path"] == "output/step_view_project/ppi_pyrosetta_residues.csv"

    run_manifest = _read_json(manifest_path)
    assert run_manifest["step_status"]["step3"] == "complete"
    index_text = index_path.read_text(encoding="utf-8")
    assert "`step3_ppi_postprocess`" in index_text


def test_record_step3_outputs_falls_back_to_ppi_subdir_and_marks_missing_optional_report(tmp_path: Path) -> None:
    config_path, project_root = _make_config(tmp_path)
    _seed_step3_outputs(project_root, tmp_path, use_fallback_ppi_dir=True, include_phase1_report=False)

    step3_dir = record_step3_outputs(config_path, repo_root=tmp_path)

    manifest = _read_json(step3_dir / "step_manifest.json")
    assert manifest["status"] == "complete"
    assert any(
        "Using fallback canonical source for ppi_pyrosetta_residues.csv" in warning
        for warning in manifest["warnings"]
    )
    assert "Optional Phase 1 interface report is not available." in manifest["warnings"]
    assert (step3_dir / "ppi_pyrosetta_residues.csv").exists()
    assert not (step3_dir / "phase1_interface_report.md").exists()


def test_step_collectors_refresh_root_views_without_changing_canonical_outputs(tmp_path: Path) -> None:
    config_path, project_root = _make_config(tmp_path)
    canonical_contents = _seed_project_root(project_root)
    original_root = {
        name: (project_root / name).read_text(encoding="utf-8")
        for name in canonical_contents
    }

    step4_dir = record_step4_outputs(config_path, repo_root=tmp_path)
    step5_dir = record_step5_outputs(config_path, repo_root=tmp_path)
    step6_dir = record_step6_outputs(config_path, repo_root=tmp_path)
    manifest_path, index_path = refresh_root_step_views(
        config_path,
        repo_root=tmp_path,
        execution_mode="full",
        ppi_config_paths=["config/ppi_prod_TH1.ini", "config/ppi_prod_beta_meander.ini"],
    )

    assert step4_dir == project_root / "step4_vina_postprocess"
    assert step5_dir == project_root / "step5_verdict"
    assert step6_dir == project_root / "step6_report"
    assert (step4_dir / "vina_pocket_table.csv").exists()
    assert (step5_dir / "valid_sites.csv").exists()
    assert (step6_dir / "project_report.txt").exists()

    step4_manifest = _read_json(step4_dir / "step_manifest.json")
    step5_manifest = _read_json(step5_dir / "step_manifest.json")
    step6_manifest = _read_json(step6_dir / "step_manifest.json")
    assert step4_manifest["status"] == "complete"
    assert step5_manifest["status"] == "complete"
    assert step6_manifest["status"] == "complete"
    assert "Optional artifact missing: vina_consensus_sites.csv" in step5_manifest["warnings"]

    run_manifest = _read_json(manifest_path)
    assert run_manifest["step_status"]["step4"] == "complete"
    assert run_manifest["step_status"]["step5"] == "complete"
    assert run_manifest["step_status"]["step6"] == "complete"
    assert run_manifest["ppi_config_paths"] == [
        "config/ppi_prod_TH1.ini",
        "config/ppi_prod_beta_meander.ini",
    ]

    index_text = index_path.read_text(encoding="utf-8")
    assert "## Run Summary" in index_text
    assert "## Step Overview Table" in index_text
    assert "## Where To Read First" in index_text
    assert "## Raw Debug Paths" in index_text
    assert "## Notes and Warnings" in index_text

    for name, original_text in original_root.items():
        assert (project_root / name).read_text(encoding="utf-8") == original_text


def test_collectors_record_missing_required_and_optional_artifacts(tmp_path: Path) -> None:
    config_path, project_root = _make_config(tmp_path)
    _seed_project_root(project_root, include_optional_step5=True)

    (project_root / "vina_pocket_bootstrap.csv").unlink()
    (project_root / "combined_residue_evidence.csv").unlink()

    record_step4_outputs(config_path, repo_root=tmp_path)
    record_step5_outputs(config_path, repo_root=tmp_path)
    record_step6_outputs(config_path, repo_root=tmp_path)

    step4_manifest = _read_json(project_root / "step4_vina_postprocess" / "step_manifest.json")
    step5_manifest = _read_json(project_root / "step5_verdict" / "step_manifest.json")
    step6_manifest = _read_json(project_root / "step6_report" / "step_manifest.json")

    assert step4_manifest["status"] == "partial"
    assert step4_manifest["missing_files"] == ["vina_pocket_bootstrap.csv"]
    assert (project_root / "step4_vina_postprocess" / "vina_pose_table.csv").exists()
    assert not (project_root / "step4_vina_postprocess" / "vina_pocket_bootstrap.csv").exists()

    assert step5_manifest["status"] == "complete"
    assert step5_manifest["warnings"] == []
    assert (project_root / "step5_verdict" / "vina_consensus_sites.csv").exists()

    assert step6_manifest["status"] == "partial"
    assert step6_manifest["missing_files"] == ["combined_residue_evidence.csv"]
    assert (project_root / "step6_report" / "project_report.txt").exists()
    assert not (project_root / "step6_report" / "combined_residue_evidence.csv").exists()


def test_build_current_run_manifest_reports_step_statuses(tmp_path: Path) -> None:
    config_path, project_root = _make_config(tmp_path)
    _seed_project_root(project_root)
    record_step4_outputs(config_path, repo_root=tmp_path)

    manifest = build_current_run_manifest(
        config_path,
        repo_root=tmp_path,
        execution_mode="from:4",
        ppi_config_paths=["config/ppi_prod_TH1.ini"],
        pyrosetta_raw_run_paths=["/ppi/raw/TH1"],
    )

    assert manifest["project_name"] == "step_view_project"
    assert manifest["execution_mode"] == "from:4"
    assert manifest["step_status"]["step4"] == "complete"
    assert manifest["step_status"]["step5"] == "not_generated"
    assert manifest["pyrosetta_raw_run_paths"] == ["/ppi/raw/TH1"]


def test_record_step7_outputs_persists_validation_result_files(tmp_path: Path) -> None:
    config_path, project_root = _make_config(tmp_path)
    _seed_project_root(project_root)

    validation_result = ValidationResult()
    validation_result.ok("vina_pose_table.csv exists (1 rows)")
    validation_result.warn("project_report.txt exists but is empty")
    validation_result.fail("Missing required artifact: vina_drug_pocket_map.csv")
    validation_result.fail("Schema mismatch: valid_sites.csv missing required column verdict")

    step7_dir = record_step7_outputs(
        config_path,
        repo_root=tmp_path,
        validation_result=validation_result,
    )
    manifest_path, index_path = refresh_root_step_views(
        config_path,
        repo_root=tmp_path,
        execution_mode="from:7",
    )

    assert step7_dir == project_root / "step7_validate"
    assert (step7_dir / "validation_status.json").exists()
    assert (step7_dir / "validation_summary.txt").exists()

    status = _read_json(step7_dir / "validation_status.json")
    assert status["status"] == "failed"
    assert status["project_root"] == "output/step_view_project"
    assert status["missing_files"] == ["Missing required artifact: vina_drug_pocket_map.csv"]
    assert status["schema_errors"] == [
        "Schema mismatch: valid_sites.csv missing required column verdict"
    ]
    assert status["warnings"] == ["project_report.txt exists but is empty"]
    assert status["validated_steps"] == [1, 2, 3, 4, 5, 6, 7]
    assert status["pass_count"] == 1
    assert status["warning_count"] == 1
    assert status["failure_count"] == 2

    summary_text = (step7_dir / "validation_summary.txt").read_text(encoding="utf-8")
    assert "Overall status: failed" in summary_text
    assert "Missing artifacts:" in summary_text
    assert "Schema mismatches:" in summary_text
    assert "Warning summary:" in summary_text
    assert "Validation failures:" in summary_text
    assert "Next action:" in summary_text

    manifest = _read_json(step7_dir / "step_manifest.json")
    assert manifest["status"] == "failed"
    assert manifest["artifact_paths"] == [
        "validation_status.json",
        "validation_summary.txt",
    ]

    run_manifest = _read_json(manifest_path)
    assert run_manifest["step_status"]["step7"] == "failed"
    index_text = index_path.read_text(encoding="utf-8")
    assert "`step7_validate/validation_status.json`" in index_text


def test_record_step7_outputs_persists_validation_error_state(tmp_path: Path) -> None:
    config_path, project_root = _make_config(tmp_path)
    _seed_project_root(project_root)

    step7_dir = record_step7_outputs(
        config_path,
        repo_root=tmp_path,
        error_message="validation crashed",
    )

    status = _read_json(step7_dir / "validation_status.json")
    assert status["status"] == "error"
    assert status["error_message"] == "validation crashed"
    assert status["failure_count"] == 0
    assert status["warning_count"] == 0

    manifest = _read_json(step7_dir / "step_manifest.json")
    assert manifest["status"] == "error"
    assert manifest["error_message"] == "validation crashed"

    summary_text = (step7_dir / "validation_summary.txt").read_text(encoding="utf-8")
    assert "Overall status: error" in summary_text
    assert "Fix the validation exception and rerun Phase 7." in summary_text


def test_step_output_view_enabled_honors_config_and_cli_override() -> None:
    assert step_output_view_enabled({}) is True
    assert step_output_view_enabled({"step_output_view": {"enabled": False}}) is False
    assert step_output_view_enabled({"step_output_view_enabled": False}) is False
    assert step_output_view_enabled({"step_output_view": {"enabled": True}}, cli_disabled=True) is False


def test_disabled_step_layer_preserves_canonical_outputs_without_creating_step_views(
    monkeypatch,
    tmp_path: Path,
) -> None:
    config_path, project_root = _make_config(tmp_path, step_view_enabled_flag=False)
    call_log: List[str] = []

    monkeypatch.setattr(run_production, "CONFIG_PATH", config_path)
    monkeypatch.setattr(run_production, "REPO_ROOT", tmp_path)
    monkeypatch.setattr(run_production, "PPI_TARGETS", [])
    monkeypatch.setattr(
        run_production,
        "PHASES",
        [(5, "Phase 5", _phase5_only_writer(project_root, call_log))],
    )
    monkeypatch.setattr(
        run_production,
        "PHASE_CHECKS",
        {5: ("step5", _phase5_check(project_root))},
    )
    monkeypatch.setattr(run_production, "banner", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(sys, "argv", ["run_production.py", "--only", "5"])

    run_production.main()

    assert call_log == ["phase5"]
    assert (project_root / "valid_sites.csv").exists()
    assert (project_root / "cross_method_agreement.csv").exists()
    assert not (project_root / "current_run_manifest.json").exists()
    assert not (project_root / "step_index.md").exists()
    for spec in STEP_SPECS.values():
        assert not (project_root / spec.folder_name).exists()
