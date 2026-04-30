import csv
import json
import os
import subprocess
import sys
import uuid
from pathlib import Path

import pytest

from egfr_myo1d.core.logging_utils import initialize_logs
from egfr_myo1d.core.run_context import RunContext, RunContextError
from egfr_myo1d.io.residue_mapping import read_residue_mapping
from egfr_myo1d.m2.ppi_inputs import generate_m2_1_ppi_inputs
from egfr_myo1d.orchestrator.prepare_inputs import run_prepare_inputs
from egfr_myo1d.ppi.collect_ppi_outputs import collect_ppi_outputs
from egfr_myo1d.ppi.consensus_patch import (
    CONSENSUS_FIELDS,
    build_m2_4_consensus_patch,
)
from egfr_myo1d.ppi.pyrosetta_adapter import generate_pyrosetta_harness
from egfr_myo1d.ppi.run_ppi_job import dry_run_job


REPO_ROOT = Path(__file__).resolve().parents[2]
INTEGRATION_FIXTURE = (
    REPO_ROOT / "fresh" / "tests" / "fixtures" / "m1_phase8_integration"
)
RAW_CONTACT_FIELDS = [
    "run_id",
    "state_id",
    "state_role",
    "method",
    "seed",
    "job_name",
    "pose_id",
    "score",
    "protomer_id",
    "egfr_runtime_residue",
    "egfr_resname",
    "myo1d_residue",
    "myo1d_resname",
    "contact_type",
    "min_distance_angstrom",
    "egfr_atom",
    "myo1d_atom",
    "egfr_contact_centroid_x",
    "egfr_contact_centroid_y",
    "egfr_contact_centroid_z",
    "pose_acceptance_class",
    "accepted_pose_flag",
    "membrane_side_class",
]


def unique_run_id(prefix="pytest_m2_phase4"):
    return "{0}_{1}".format(prefix, uuid.uuid4().hex[:12])


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


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def read_csv(path):
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_rows(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    return path


@pytest.fixture
def ctx_with_m2_2(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    aggregate = run_prepare_inputs(
        ctx,
        mode="smoke_input",
        profile="codex_dev",
        input_root=INTEGRATION_FIXTURE,
    )
    assert aggregate.status in {"PASS", "PASS_WITH_WARNINGS"}
    m2_1 = generate_m2_1_ppi_inputs(ctx)
    assert m2_1.status in {"PASS", "PASS_WITH_WARNINGS"}
    m2_2 = generate_pyrosetta_harness(ctx, mode="smoke_input")
    assert m2_2.status in {"PASS", "PASS_WITH_WARNINGS"}
    dry_run_job(ctx, m2_2.jobs[0].job_name)
    return ctx, m2_2


def _same_residue_a_b_mapping(ctx, state_id="EGFR_160-185"):
    rows = read_residue_mapping(ctx.qc_dir / "{0}_receptor_mapping.csv".format(state_id))
    a_rows = {row["source_resseq"]: row for row in rows if row["protomer_id"] == "A"}
    for row in rows:
        if row["protomer_id"] == "B" and row["source_resseq"] in a_rows:
            return a_rows[row["source_resseq"]], row
    raise AssertionError("fixture does not contain A/B paired mapping rows")


def _synthetic_raw_contact_rows(ctx, job_name):
    a_row, b_row = _same_residue_a_b_mapping(ctx)
    base = {
        "run_id": ctx.run_id,
        "state_id": "EGFR_160-185",
        "state_role": "primary_membrane_validated_state",
        "method": "pyrosetta_global_ppi",
        "seed": "0",
        "job_name": job_name,
        "score": "-1.5",
        "egfr_resname": a_row["source_resname"],
        "myo1d_resname": "GLY",
        "contact_type": "ca_contact",
        "min_distance_angstrom": "7.2",
        "egfr_atom": "CA",
        "myo1d_atom": "CA",
        "pose_acceptance_class": "accepted_active_8_9_supported_12",
        "accepted_pose_flag": "true",
        "membrane_side_class": "solvent_exposed",
    }
    a = dict(base)
    a.update(
        {
            "pose_id": "pose_a",
            "protomer_id": "A",
            "egfr_runtime_residue": a_row["runtime_resseq"],
            "myo1d_residue": "961",
            "egfr_contact_centroid_x": "1.0",
            "egfr_contact_centroid_y": "2.0",
            "egfr_contact_centroid_z": "3.0",
        }
    )
    b = dict(base)
    b.update(
        {
            "pose_id": "pose_b",
            "protomer_id": "B",
            "egfr_runtime_residue": b_row["runtime_resseq"],
            "myo1d_residue": "993",
            "egfr_contact_centroid_x": "4.0",
            "egfr_contact_centroid_y": "5.0",
            "egfr_contact_centroid_z": "6.0",
        }
    )
    return [a, b]


def test_m2_4_cli_help_includes_build_consensus_patch():
    proc = subprocess.run(
        [sys.executable, "-m", "egfr_myo1d.cli", "--help"],
        cwd=str(REPO_ROOT),
        env=_cli_env(),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    assert proc.returncode == 0
    out = proc.stdout.decode("utf-8", errors="replace") + proc.stderr.decode(
        "utf-8", errors="replace"
    )
    assert "build-m2-ppi-consensus-patch" in out


def test_m2_4_dry_run_no_contacts_writes_schema_valid_warning_output(ctx_with_m2_2):
    ctx, _m2_2 = ctx_with_m2_2
    m2_3 = collect_ppi_outputs(ctx)
    assert m2_3.status == "PASS_WITH_WARNINGS"

    report = build_m2_4_consensus_patch(ctx)

    assert report.status == "PASS_WITH_WARNINGS"
    assert report.patch_count == 0
    assert report.consensus_patch_csv.is_file()
    with report.consensus_patch_csv.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        assert header == CONSENSUS_FIELDS
        assert list(reader) == []
    payload = read_json(report.manifest_path)
    assert payload["docking_executed"] is False
    assert payload["pocket_discovery_executed"] is False


def test_m2_4_aggregates_restored_contacts_and_preserves_chain_protomer(ctx_with_m2_2):
    ctx, m2_2 = ctx_with_m2_2
    raw_table = ctx.run_dir / "phase1_ppi" / "tables" / "synthetic_m2_4_contacts.csv"
    write_rows(raw_table, _synthetic_raw_contact_rows(ctx, m2_2.jobs[0].job_name), RAW_CONTACT_FIELDS)
    m2_3 = collect_ppi_outputs(ctx, raw_contact_table=raw_table)
    assert m2_3.status == "PASS"

    report = build_m2_4_consensus_patch(ctx)

    assert report.status == "PASS"
    assert report.patch_count == 2
    rows = read_csv(report.consensus_patch_csv)
    assert [row["egfr_protomer_id"] for row in rows] == ["A", "B"]
    assert [row["egfr_chain_id"] for row in rows] == ["A", "B"]
    assert {row["evidence_status"] for row in rows} == {"accepted_contact_evidence"}
    assert {row["support_fraction"] for row in rows} == {"0.200"}
    assert rows[0]["myo1d_sheet_support"] == "sheet8"
    assert rows[1]["myo1d_sheet_support"] == "sheet12"
    assert (ctx.run_dir / "phase1_ppi" / "tables" / "ppi_consensus_patch.csv").is_file()
    state_support = read_csv(report.state_support_csv)
    assert state_support[0]["unique_pose_count"] == "5"
    assert state_support[0]["unique_seed_count"] == "1"

    symmetry = read_csv(report.symmetry_support_csv)
    assert any(
        row["protomer_a_supported"] == "True"
        and row["protomer_b_supported"] == "True"
        and row["symmetry_support_count"] == "1"
        for row in symmetry
    )


def test_m2_4_excludes_pending_pose_qc_contacts_from_consensus(ctx_with_m2_2):
    ctx, m2_2 = ctx_with_m2_2
    raw_rows = _synthetic_raw_contact_rows(ctx, m2_2.jobs[0].job_name)
    for row in raw_rows:
        row["pose_acceptance_class"] = ""
        row["accepted_pose_flag"] = ""
    raw_table = ctx.run_dir / "phase1_ppi" / "tables" / "synthetic_m2_4_pending_contacts.csv"
    write_rows(raw_table, raw_rows, RAW_CONTACT_FIELDS)
    m2_3 = collect_ppi_outputs(ctx, raw_contact_table=raw_table)
    assert m2_3.status == "PASS"

    report = build_m2_4_consensus_patch(ctx)

    assert report.status == "PASS_WITH_WARNINGS"
    assert report.patch_count == 0
    assert report.accepted_pose_count == 0
    assert report.pending_pose_qc_count == 2
    assert "contacts_lack_pose_qc_acceptance_class" in report.warnings
    assert "no_accepted_contact_evidence_after_pose_qc_filter" in report.warnings
    assert read_csv(report.consensus_patch_csv) == []


def test_m2_4_detects_residue_number_drift(ctx_with_m2_2):
    ctx, m2_2 = ctx_with_m2_2
    raw_table = ctx.run_dir / "phase1_ppi" / "tables" / "synthetic_m2_4_contacts.csv"
    write_rows(raw_table, _synthetic_raw_contact_rows(ctx, m2_2.jobs[0].job_name), RAW_CONTACT_FIELDS)
    m2_3 = collect_ppi_outputs(ctx, raw_contact_table=raw_table)
    assert m2_3.status == "PASS"
    contacts = read_csv(m2_3.pose_contacts_csv)
    contacts[0]["egfr_original_residue"] = "1"
    write_rows(m2_3.pose_contacts_csv, contacts, contacts[0].keys())

    report = build_m2_4_consensus_patch(ctx)

    assert report.status == "FAIL"
    assert any("residue_number_drift" in blocker for blocker in report.blockers)
    qc_rows = read_csv(report.qc_csv_path)
    assert any(row["check"] == "chain_protomer_identity_preserved" and row["status"] == "FAIL" for row in qc_rows)


def test_m2_4_detects_chain_and_protomer_identity_drift(ctx_with_m2_2):
    ctx, m2_2 = ctx_with_m2_2
    raw_table = ctx.run_dir / "phase1_ppi" / "tables" / "synthetic_m2_4_contacts.csv"
    write_rows(raw_table, _synthetic_raw_contact_rows(ctx, m2_2.jobs[0].job_name), RAW_CONTACT_FIELDS)
    m2_3 = collect_ppi_outputs(ctx, raw_contact_table=raw_table)
    assert m2_3.status == "PASS"
    base_contacts = read_csv(m2_3.pose_contacts_csv)

    chain_drift_contacts = [dict(row) for row in base_contacts]
    chain_drift_contacts[0]["egfr_chain_id"] = "Z"
    write_rows(m2_3.pose_contacts_csv, chain_drift_contacts, chain_drift_contacts[0].keys())
    chain_report = build_m2_4_consensus_patch(ctx)
    assert chain_report.status == "FAIL"
    assert any("chain_id_drift" in blocker for blocker in chain_report.blockers)

    protomer_drift_contacts = [dict(row) for row in base_contacts]
    protomer_drift_contacts[0]["protomer_id"] = "Z"
    write_rows(m2_3.pose_contacts_csv, protomer_drift_contacts, protomer_drift_contacts[0].keys())
    protomer_report = build_m2_4_consensus_patch(ctx)
    assert protomer_report.status == "FAIL"
    assert any("protomer_id_drift" in blocker for blocker in protomer_report.blockers)


def test_m2_4_rejects_contact_table_outside_run_dir(ctx_with_m2_2, tmp_path):
    ctx, _m2_2 = ctx_with_m2_2
    outside = tmp_path / "outside_contacts.csv"
    write_rows(outside, [], CONSENSUS_FIELDS)

    with pytest.raises(RunContextError):
        build_m2_4_consensus_patch(ctx, contact_table=outside)


def test_m2_4_outputs_stay_under_run_dir(ctx_with_m2_2):
    ctx, m2_2 = ctx_with_m2_2
    raw_table = ctx.run_dir / "phase1_ppi" / "tables" / "synthetic_m2_4_contacts.csv"
    write_rows(raw_table, _synthetic_raw_contact_rows(ctx, m2_2.jobs[0].job_name), RAW_CONTACT_FIELDS)
    collect_ppi_outputs(ctx, raw_contact_table=raw_table)

    report = build_m2_4_consensus_patch(ctx)

    for path in [
        report.consensus_patch_csv,
        report.m2_table_patch_csv,
        report.residue_csv,
        report.state_support_csv,
        report.symmetry_support_csv,
        report.qc_csv_path,
        report.summary_report_path,
        report.manifest_path,
    ]:
        path.resolve().relative_to(ctx.run_dir.resolve())


def _cli_env():
    env = os.environ.copy()
    env["PYTHONPATH"] = (
        str(REPO_ROOT / "fresh" / "src") + os.pathsep + env.get("PYTHONPATH", "")
    )
    return env
