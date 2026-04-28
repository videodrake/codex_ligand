import json
import os
import subprocess
import sys
import uuid
from pathlib import Path

import pytest

from egfr_myo1d.core.logging_utils import initialize_logs
from egfr_myo1d.core.run_context import RunContext, RunContextError
from egfr_myo1d.preparation.constructs import prepare_myo1d_construct
from egfr_myo1d.structure.pdb_parser import parse_pdb
from egfr_myo1d.validation.prepared_inputs import (
    audit_receptor,
    prepare_ppi_inputs,
    read_contract,
)


REPO_ROOT = Path(__file__).resolve().parents[2]
TASK4_INPUTS = REPO_ROOT / "fresh" / "tests" / "fixtures" / "task4_inputs"


def unique_run_id(prefix="pytest_task4"):
    return "{}_{}".format(prefix, uuid.uuid4().hex[:12])


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


def load_contract():
    return read_contract(TASK4_INPUTS, None)[0]


def run_prepare(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    report = prepare_ppi_inputs(ctx, "smoke_env", "codex_dev", TASK4_INPUTS)
    return ctx, report


def test_cli_help_lists_prepare_ppi_inputs():
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
    assert "prepare-ppi-inputs" in proc.stdout


def test_safe_run_id_passes_and_path_traversal_fails(tmp_path):
    ctx = RunContext.create(unique_run_id(), "smoke_env", repo_root=REPO_ROOT)
    assert ctx.run_dir.is_dir()
    with pytest.raises(RunContextError):
        RunContext.create("../bad_run", "smoke_env", repo_root=REPO_ROOT)


def test_prepared_egfr_preserves_ab_chains_residue_numbers_and_dimer_context(tmp_path):
    ctx, report = run_prepare(tmp_path)
    output = ctx.run_dir / "prepared" / "egfr" / "egfr_receptor_normalized.pdb"
    structure = parse_pdb(output)
    manifest = json.loads((ctx.manifest_dir / "prepared_input_manifest.json").read_text(encoding="utf-8"))

    assert report["status"] == "WARN"
    assert structure.chain_ids == ("A", "B")
    assert 924 in structure.residue_numbers_for_chain("A")
    assert manifest["egfr"]["dimer_context_preserved"] is True
    assert manifest["egfr"]["residue_numbering_preserved"] is True


def test_single_chain_ambiguous_receptor_warns_only_when_whitelisted():
    contract = load_contract()
    structure = parse_pdb(TASK4_INPUTS / "egfr_single_chain_X_ambiguous.pdb")
    status, _rows, messages = audit_receptor(
        structure,
        contract["receptor"],
        contract,
        "smoke_env",
        "codex_dev",
        fixture_id="egfr_single_chain_X_ambiguous",
    )

    assert status == "WARN"
    assert {message["classification"] for message in messages} == {"fixture_warning"}


def test_v924r_fixture_warns_and_is_not_mutated_when_whitelisted():
    contract = load_contract()
    structure = parse_pdb(TASK4_INPUTS / "egfr_v924r_warn.pdb")
    status, _rows, messages = audit_receptor(
        structure,
        contract["receptor"],
        contract,
        "smoke_env",
        "codex_dev",
        fixture_id="egfr_v924r_warn",
    )

    assert status == "WARN"
    assert any(message["code"] == "v924r_like_mutation" for message in messages)
    assert "ARG" in structure.residue_resnames(924)
    assert "VAL" not in structure.residue_resnames(924)


def test_non_whitelisted_warning_classes_fail_or_quarantine():
    contract = load_contract()
    contract["task3_closeout"]["fixture_warning_whitelist"] = []
    ambiguous = parse_pdb(TASK4_INPUTS / "egfr_single_chain_X_ambiguous.pdb")
    v924r = parse_pdb(TASK4_INPUTS / "egfr_v924r_warn.pdb")

    ambiguous_status, _rows, ambiguous_messages = audit_receptor(
        ambiguous, contract["receptor"], contract, "smoke_input", "hpc_strict", strict=True
    )
    v924r_status, _rows, v924r_messages = audit_receptor(
        v924r, contract["receptor"], contract, "smoke_input", "hpc_strict", strict=True
    )

    assert ambiguous_status == "FAIL"
    assert v924r_status == "FAIL"
    assert all(message["classification"] == "production_blocker" for message in ambiguous_messages + v924r_messages)


def test_task3_closeout_warning_whitelist_serialized(tmp_path):
    ctx, report = run_prepare(tmp_path)
    manifest = json.loads((ctx.manifest_dir / "prepared_input_manifest.json").read_text(encoding="utf-8"))
    qc = json.loads((ctx.manifest_dir / "preparation_qc_report.json").read_text(encoding="utf-8"))

    assert manifest["task3_closeout"]["warning_whitelist_applied"] is True
    assert "egfr_v924r_warn" in manifest["task3_closeout"]["fixture_warning_whitelist"]
    assert qc["task3_closeout"]["production_same_warning_policy"] == "fail_or_quarantine"


def test_myo1d_primary_construct_emits_955_1001_original_numbering(tmp_path):
    ctx, _report = run_prepare(tmp_path)
    output = ctx.run_dir / "prepared" / "myo1d" / "MYO1D_sheet8_9_12_core_955_1001.pdb"
    residues = set()
    for residue in parse_pdb(output).residues:
        if residue.resname not in {"ACE", "NME"}:
            residues.add(residue.residue_number)

    assert min(residues) == 955
    assert max(residues) == 1001
    assert 1002 not in residues


def test_myo1d_comparator_emits_955_1006_and_marks_998_plus_noise(tmp_path):
    ctx, _report = run_prepare(tmp_path)
    output = ctx.run_dir / "prepared" / "myo1d" / "MYO1D_ext_beta_meander_955_1006_tail_masked.pdb"
    residues = {residue.residue_number for residue in parse_pdb(output).residues}
    audit = (ctx.qc_dir / "myo1d_construct_audit.csv").read_text(encoding="utf-8")

    assert min(residues) == 955
    assert max(residues) == 1006
    assert "contact_monitoring_cap/noise_monitor" in audit


def test_bad_962_1006_terminal_artifact_detected_and_not_emitted_as_partner(tmp_path):
    ctx, _report = run_prepare(tmp_path)
    terminal_audit = (ctx.qc_dir / "terminal_artifact_audit.csv").read_text(encoding="utf-8")

    assert "myo1d_962_1006_terminal_artifact_bad" in terminal_audit
    assert "starts at 962" in terminal_audit
    assert not (ctx.run_dir / "prepared" / "myo1d" / "myo1d_962_1006_terminal_artifact_bad.pdb").exists()


def test_hetatm_cap_regression_keeps_ile1000_and_caps(tmp_path):
    contract = load_contract()
    source = TASK4_INPUTS / "myo1d_capped_hetatm_955_1001.pdb"
    output = tmp_path / "capped_prepared.pdb"
    record, audit_rows, terminal = prepare_myo1d_construct(
        source,
        output,
        "MYO1D_capped_regression",
        (955, 1001),
        contract["myo1d"],
        "regression_capped",
    )
    prepared = parse_pdb(output)

    assert any(atom.record_type == "HETATM" and atom.resname == "ILE" and atom.residue_number == 1000 for atom in prepared.atoms)
    assert any(atom.record_type == "HETATM" and atom.resname == "ACE" for atom in prepared.atoms)
    assert any(atom.record_type == "HETATM" and atom.resname == "NME" for atom in prepared.atoms)
    assert any(row["residue_name"] == "ILE" and row["record_type"] == "HETATM" for row in audit_rows)


def test_active_face_contract_contents_and_no_score_bonus(tmp_path):
    ctx, _report = run_prepare(tmp_path)
    data = json.loads((ctx.run_dir / "prepared" / "restraints" / "myo1d_active_face_contract.json").read_text(encoding="utf-8"))

    assert {961, 962, 963, 964, 968, 969, 970, 971, 972} <= set(data["primary_active_face_residues"])
    assert {993, 994, 995, 996, 997} <= set(data["sheet12_support_residues"])
    assert data["score_bonus_allowed"] is False


def test_membrane_exclusion_mask_written_and_references_frame(tmp_path):
    ctx, _report = run_prepare(tmp_path)
    data = json.loads((ctx.run_dir / "prepared" / "restraints" / "egfr_membrane_exclusion_mask.json").read_text(encoding="utf-8"))

    assert data["membrane_frame_source"] == "membrane_frame_valid.json"
    assert data["membrane_proximal_policy"] == "exclude_or_warn_for_future_ppi_sampling"


def test_manifest_contains_output_hashes_and_warning_classification(tmp_path):
    ctx, _report = run_prepare(tmp_path)
    manifest = json.loads((ctx.manifest_dir / "prepared_input_manifest.json").read_text(encoding="utf-8"))

    assert manifest["outputs"]["egfr_receptor"]["sha256"]
    assert manifest["outputs"]["myo1d_primary"]["sha256"]
    assert manifest["qc_summary"]["fixture_warning"] >= 1
    assert manifest["qc_summary"]["explicitly_whitelisted_comparator"] >= 1


def test_status_reports_prepare_phase(tmp_path, capsys):
    from egfr_myo1d.core.logging_utils import append_phase_status
    from egfr_myo1d.core.manifest import initialize_manifests
    from egfr_myo1d.cli import main

    run_id = unique_run_id()
    ctx = RunContext.create(run_id, "smoke_env", repo_root=REPO_ROOT)
    initialize_logs(ctx)
    initialize_manifests(ctx, "smoke_env")
    append_phase_status(ctx, "prepare-ppi-inputs", "WARN", "PPI input preparation completed with WARN", {})

    assert main(["status", "--run-id", run_id]) == 0
    output = capsys.readouterr().out
    assert "last_phase=prepare-ppi-inputs" in output


def test_old_workflow_diff_check_prints_nothing():
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
