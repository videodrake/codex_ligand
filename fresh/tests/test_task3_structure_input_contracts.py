import json
import os
import subprocess
import sys
import uuid
from pathlib import Path

import pytest

from egfr_myo1d.core.logging_utils import initialize_logs
from egfr_myo1d.core.run_context import RunContext, RunContextError
from egfr_myo1d.structure.contracts import validate_receptor_contract
from egfr_myo1d.structure.geometry import validate_membrane_frame
from egfr_myo1d.structure.myo1d_annotation import validate_myo1d_file
from egfr_myo1d.structure.pdb_parser import parse_pdb
from egfr_myo1d.validation.structure_inputs import validate_structure_inputs


REPO_ROOT = Path(__file__).resolve().parents[2]
TASK3_INPUTS = REPO_ROOT / "fresh" / "tests" / "fixtures" / "task3_inputs"


def unique_run_id(prefix="pytest_task3"):
    return "{}_{}".format(prefix, uuid.uuid4().hex[:12])


def receptor_contract(path, allow_single=True):
    return {
        "receptor_id": Path(path).stem,
        "path": path,
        "expected_assembly": "dimer",
        "expected_protomers": [
            {"protomer_id": "A", "chain_ids": ["A"], "required_ranges": [[699, 700]]},
            {"protomer_id": "B", "chain_ids": ["B"], "required_ranges": [[699, 700]]},
        ],
        "numbering_system": "project_uniprot_or_documented_mapping",
        "allow_single_chain_dimer_with_mapping": allow_single,
        "known_mutation_checks": [
            {"label": "3GT8_V924R", "residue_number": 924, "expected_wt": "VAL", "warn_if_resname": "ARG"}
        ],
    }


def test_pdb_parser_preserves_chain_ids_and_residue_numbers():
    structure = parse_pdb(TASK3_INPUTS / "egfr_valid_dimer_AB.pdb")

    assert structure.chain_ids == ("A", "B")
    assert structure.residue_numbers_for_chain("A") == {699, 700}
    assert structure.residue_numbers_for_chain("B") == {699, 700}


def test_pdb_parser_includes_hetatm_records():
    structure = parse_pdb(REPO_ROOT / "fresh" / "tests" / "fixtures" / "mini_atom_hetatm_cap.pdb")

    record_types = {atom.record_type for atom in structure.atoms}
    assert {"ATOM", "HETATM"} <= record_types
    assert any(residue.resname == "ACE" for residue in structure.residues)


def test_valid_egfr_dimer_fixture_returns_not_fail():
    report, _chain_rows, _residue_rows, _manifest = validate_receptor_contract(
        receptor_contract("egfr_valid_dimer_AB.pdb"), TASK3_INPUTS, "smoke_env", "codex_dev"
    )

    assert report["status"] in {"PASS", "WARN"}
    assert report["status"] != "FAIL"


def test_single_chain_dimer_ambiguity_returns_warn():
    report, _chain_rows, _residue_rows, _manifest = validate_receptor_contract(
        receptor_contract("egfr_single_chain_dimer_ambiguous_X.pdb"), TASK3_INPUTS, "smoke_env", "codex_dev"
    )

    assert report["status"] == "WARN"
    assert any(message["code"] == "single_chain_dimer_ambiguity" for message in report["messages"])


def test_residue_reset_fixture_fails_in_strict_validation():
    report, _chain_rows, _residue_rows, _manifest = validate_receptor_contract(
        receptor_contract("egfr_residue_reset_bad.pdb"), TASK3_INPUTS, "smoke_input", "hpc_strict"
    )

    assert report["status"] == "FAIL"
    assert any(message["code"] == "residue_number_reset_risk" for message in report["messages"])


def test_v924r_fixture_returns_warn():
    contract = receptor_contract("egfr_v924r_warn.pdb")
    contract["expected_protomers"] = [
        {"protomer_id": "A", "chain_ids": ["A"], "required_ranges": [[923, 925]]},
        {"protomer_id": "B", "chain_ids": ["B"], "required_ranges": [[923, 925]]},
    ]
    report, _chain_rows, _residue_rows, _manifest = validate_receptor_contract(
        contract, TASK3_INPUTS, "smoke_env", "codex_dev"
    )

    assert report["status"] == "WARN"
    assert any(message["code"] == "3GT8_V924R" for message in report["messages"])


def test_valid_myo1d_955_1001_fixture_returns_pass():
    report = validate_myo1d_file(TASK3_INPUTS / "myo1d_955_1001_valid.pdb", "smoke_env", "codex_dev")

    assert report["status"] == "PASS"
    assert report["construct_class"] == "955-1001"
    assert report["key_residue_bonus_weight"] == 0.0


def test_myo1d_962_start_fixture_reports_terminal_artifact():
    report = validate_myo1d_file(
        TASK3_INPUTS / "myo1d_962_1006_bad_terminal.pdb", "smoke_env", "codex_dev"
    )

    assert report["status"] in {"WARN", "FAIL"}
    assert any(message["code"] == "terminal_artifact_start_962" for message in report["messages"])


def test_myo1d_955_1006_fixture_warns_tail_noise_zone():
    report = validate_myo1d_file(
        TASK3_INPUTS / "myo1d_955_1006_tail_masked_warn.pdb", "smoke_env", "codex_dev"
    )

    assert report["status"] == "WARN"
    assert report["construct_class"] == "955-1006"
    assert any(message["code"] == "tail_noise_zone_present" for message in report["messages"])


def test_valid_membrane_frame_returns_pass_and_normalized_vectors():
    data = json.loads((TASK3_INPUTS / "membrane_frame_valid.json").read_text(encoding="utf-8"))
    report = validate_membrane_frame(data)

    assert report["status"] == "PASS"
    assert report["membrane_normal_unit"] == [0.0, 0.0, 1.0]
    assert report["dimer_axis_unit"] == [1.0, 0.0, 0.0]


def test_zero_membrane_normal_returns_fail():
    data = json.loads((TASK3_INPUTS / "membrane_frame_zero_normal_bad.json").read_text(encoding="utf-8"))
    report = validate_membrane_frame(data)

    assert report["status"] == "FAIL"
    assert any(message["code"] == "invalid_membrane_vector" for message in report["messages"])


def make_tmp_run_context(tmp_path):
    run_id = unique_run_id()
    run_dir = tmp_path / run_id
    return RunContext(
        repo_root=REPO_ROOT,
        fresh_root=REPO_ROOT / "fresh",
        run_id=run_id,
        run_dir=run_dir,
        manifest_dir=run_dir / "manifest",
        logs_dir=run_dir / "logs",
        jobs_log_dir=run_dir / "logs" / "jobs",
        errors_dir=run_dir / "logs" / "errors",
        qc_dir=run_dir / "qc",
        reports_dir=run_dir / "reports",
        scratch_dir=run_dir / "scratch",
        tmp_dir=run_dir / "tmp",
    )


def test_cli_writes_required_task3_outputs_inside_run_dir(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    report = validate_structure_inputs(ctx, "smoke_env", "codex_dev", TASK3_INPUTS)

    assert report["status"] == "WARN"
    required_paths = [
        ctx.manifest_dir / "structure_input_manifest.json",
        ctx.manifest_dir / "structure_qc_report.json",
        ctx.qc_dir / "chain_protomer_audit.csv",
        ctx.qc_dir / "residue_mapping_audit.csv",
        ctx.qc_dir / "membrane_frame_report.json",
        ctx.qc_dir / "myo1d_annotation_report.json",
        ctx.logs_dir / "master.log",
        ctx.logs_dir / "phase_status.jsonl",
        ctx.logs_dir / "job_status.jsonl",
    ]
    for path in required_paths:
        assert path.is_file()
        assert path.resolve().is_relative_to(ctx.run_dir.resolve())


def test_validate_structures_rejects_path_traversal_run_id():
    with pytest.raises(RunContextError):
        RunContext.create("../bad_run", "smoke_env", repo_root=REPO_ROOT)


def test_legacy_workflow_paths_not_modified_by_task3_diff_check():
    proc = subprocess.run(
        [
            "git",
            "diff",
            "--name-only",
            "--",
            "run_production.py",
            "main.py",
            "egfr_pipeline/",
            "config/",
            "docs/runbook.md",
            "output/",
            "results_export/",
        ],
        cwd=str(REPO_ROOT),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )

    assert proc.returncode == 0
    assert proc.stdout.strip() == ""


def test_validate_structures_cli_help_includes_command():
    env = os.environ.copy()
    env["PYTHONPATH"] = "{}{}{}".format(REPO_ROOT / "fresh" / "src", os.pathsep, env.get("PYTHONPATH", ""))
    proc = subprocess.run(
        [sys.executable, "-m", "egfr_myo1d.cli", "--help"],
        cwd=str(REPO_ROOT),
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )

    assert proc.returncode == 0
    assert "validate-structures" in proc.stdout
