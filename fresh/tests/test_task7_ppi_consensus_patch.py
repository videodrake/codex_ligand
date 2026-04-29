import csv
import json
import os
import shutil
import subprocess
import sys
import uuid
from pathlib import Path

from egfr_myo1d.core.logging_utils import initialize_logs
from egfr_myo1d.core.run_context import RunContext
from egfr_myo1d.validation.ppi_consensus import summarize_ppi_consensus


REPO_ROOT = Path(__file__).resolve().parents[2]
TASK7_INPUTS = REPO_ROOT / "fresh" / "tests" / "fixtures" / "task7_ppi_consensus"
HEADER = [
    "receptor_id",
    "receptor_state",
    "method",
    "seed",
    "replicate",
    "pose_id",
    "protomer_id",
    "egfr_contact_residues",
    "egfr_contact_centroid_x",
    "egfr_contact_centroid_y",
    "egfr_contact_centroid_z",
    "myo1d_contact_residues",
    "myo1d_active_face_contacts",
    "sheet12_support_contacts",
    "tail_noise_contacts",
    "membrane_side_class",
    "terminal_artifact_flag",
    "atp_overlap_flag",
    "membrane_proximal_flag",
    "pose_acceptance_class",
]


def unique_run_id(prefix="pytest_task7"):
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


def copy_task7_inputs(tmp_path):
    root = tmp_path / "task7_inputs"
    shutil.copytree(TASK7_INPUTS, root)
    return root


def write_contacts(root, rows, header=None):
    header = header or HEADER
    path = root / "accepted_ppi_contacts.csv"
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=header, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    return path


def base_row(pose_id, protomer_id="A", centroid=None):
    centroid = centroid or [1.0, 2.0, 3.0]
    return {
        "receptor_id": "EGFR_synthetic_inactive_dimer_AB",
        "receptor_state": "EGFR_fixture_state",
        "method": "pyrosetta_global_ppi",
        "seed": "0",
        "replicate": "1",
        "pose_id": pose_id,
        "protomer_id": protomer_id,
        "egfr_contact_residues": "{0}:742;{0}:746;{0}:923".format(protomer_id),
        "egfr_contact_centroid_x": str(centroid[0]),
        "egfr_contact_centroid_y": str(centroid[1]),
        "egfr_contact_centroid_z": str(centroid[2]),
        "myo1d_contact_residues": "961;962;968;993",
        "myo1d_active_face_contacts": "961;962;968",
        "sheet12_support_contacts": "993",
        "tail_noise_contacts": "",
        "membrane_side_class": "solvent_exposed",
        "terminal_artifact_flag": "false",
        "atp_overlap_flag": "false",
        "membrane_proximal_flag": "false",
        "pose_acceptance_class": "accepted_active_8_9_supported_12",
    }


def run_consensus(tmp_path, root):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    report = summarize_ppi_consensus(ctx, "smoke_env", "codex_dev", root)
    return ctx, report


def read_csv(path):
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def assert_under_run_dir(ctx, path):
    try:
        path.resolve().relative_to(ctx.run_dir.resolve())
    except ValueError:
        raise AssertionError("{0} escaped run dir".format(path))


def test_cli_help_includes_summarize_ppi_consensus():
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
    assert "summarize-ppi-consensus" in proc.stdout


def test_valid_synthetic_contact_table_writes_outputs(tmp_path):
    root = copy_task7_inputs(tmp_path)
    ctx, report = run_consensus(tmp_path, root)

    required = [
        ctx.manifest_dir / "ppi_consensus_input_manifest.json",
        ctx.manifest_dir / "ppi_consensus_qc_report.json",
        ctx.qc_dir / "ppi_consensus_schema_audit.csv",
        ctx.qc_dir / "ppi_residue_parsing_audit.csv",
        ctx.qc_dir / "ppi_tail_noise_audit.csv",
        ctx.qc_dir / "ppi_active_face_audit.csv",
        ctx.qc_dir / "ppi_membrane_atp_audit.csv",
        ctx.qc_dir / "ppi_convergence_audit.csv",
        ctx.run_dir / "output" / "ppi" / "ppi_consensus_patch.csv",
        ctx.reports_dir / "ppi_consensus_report.txt",
    ]
    assert report["status"] == "WARN"
    for path in required:
        assert path.is_file()
        assert_under_run_dir(ctx, path)


def test_missing_required_schema_fails_without_consensus_csv(tmp_path):
    root = tmp_path / "bad_schema"
    root.mkdir()
    write_contacts(root, [base_row("bad")], header=[column for column in HEADER if column != "pose_id"])
    ctx, report = run_consensus(tmp_path, root)

    assert report["status"] == "FAIL"
    assert not (ctx.run_dir / "output" / "ppi" / "ppi_consensus_patch.csv").exists()
    schema_rows = read_csv(ctx.qc_dir / "ppi_consensus_schema_audit.csv")
    assert any(row["column"] == "pose_id" and row["status"] == "FAIL" for row in schema_rows)


def test_tail_dominated_contacts_are_not_promoted(tmp_path):
    root = tmp_path / "tail"
    root.mkdir()
    rows = [base_row("tail_1"), base_row("tail_2")]
    for row in rows:
        row["myo1d_contact_residues"] = "998;1002;1003"
        row["myo1d_active_face_contacts"] = ""
        row["tail_noise_contacts"] = "1002;1003"
    write_contacts(root, rows)
    ctx, report = run_consensus(tmp_path, root)
    patches = read_csv(ctx.run_dir / "output" / "ppi" / "ppi_consensus_patch.csv")

    assert report["status"] == "WARN"
    assert patches[0]["evidence_class"] == "TAIL_DOMINATED"
    assert "tail_or_terminal_artifact_visible" in patches[0]["warnings"]


def test_terminal_artifact_records_are_quarantined_in_audit(tmp_path):
    root = tmp_path / "terminal"
    root.mkdir()
    rows = [base_row("terminal_1"), base_row("terminal_2")]
    rows[0]["terminal_artifact_flag"] = "true"
    rows[0]["myo1d_contact_residues"] = "962;1006"
    rows[0]["tail_noise_contacts"] = "1006"
    write_contacts(root, rows)
    ctx, _report = run_consensus(tmp_path, root)
    tail_rows = read_csv(ctx.qc_dir / "ppi_tail_noise_audit.csv")
    patches = read_csv(ctx.run_dir / "output" / "ppi" / "ppi_consensus_patch.csv")

    assert any(row["terminal_artifact_flag"] == "True" and row["quarantine_recommended"] == "True" for row in tail_rows)
    assert patches[0]["evidence_class"] == "TAIL_DOMINATED"


def test_active_face_sheet12_annotation_has_no_score_bonus(tmp_path):
    root = copy_task7_inputs(tmp_path)
    ctx, _report = run_consensus(tmp_path, root)
    active_rows = read_csv(ctx.qc_dir / "ppi_active_face_audit.csv")
    report = read_json(ctx.manifest_dir / "ppi_consensus_qc_report.json")

    assert active_rows
    assert all(row["score_bonus_allowed"] == "False" for row in active_rows)
    assert all(row["key_residue_bonus_weight"] == "0.0" for row in active_rows)
    assert report["active_face_and_sheet12_are_annotation_only"] is True


def test_atp_overlap_flags_propagate(tmp_path):
    root = tmp_path / "atp"
    root.mkdir()
    rows = [base_row("atp_1"), base_row("atp_2")]
    for row in rows:
        row["atp_overlap_flag"] = "true"
    write_contacts(root, rows)
    ctx, _report = run_consensus(tmp_path, root)
    membrane_rows = read_csv(ctx.qc_dir / "ppi_membrane_atp_audit.csv")
    patches = read_csv(ctx.run_dir / "output" / "ppi" / "ppi_consensus_patch.csv")

    assert any(row["atp_overlap_flag"] == "True" for row in membrane_rows)
    assert patches[0]["evidence_class"] == "ATP_OVERLAPPING"
    assert patches[0]["atp_overlap_fraction"] == "1.000"


def test_membrane_proximal_flags_propagate(tmp_path):
    root = tmp_path / "membrane"
    root.mkdir()
    rows = [base_row("mem_1"), base_row("mem_2")]
    for row in rows:
        row["membrane_proximal_flag"] = "true"
        row["membrane_side_class"] = "membrane_proximal"
    write_contacts(root, rows)
    ctx, _report = run_consensus(tmp_path, root)
    patches = read_csv(ctx.run_dir / "output" / "ppi" / "ppi_consensus_patch.csv")

    assert patches[0]["evidence_class"] == "MEMBRANE_PROXIMAL"
    assert patches[0]["membrane_proximal_fraction"] == "1.000"
    assert patches[0]["membrane_side_classes"] == "membrane_proximal"


def test_dispersed_centroids_are_classified(tmp_path):
    root = tmp_path / "dispersed"
    root.mkdir()
    rows = [
        base_row("disp_1", centroid=[0.0, 0.0, 0.0]),
        base_row("disp_2", centroid=[50.0, 0.0, 0.0]),
    ]
    write_contacts(root, rows)
    ctx, _report = run_consensus(tmp_path, root)
    patches = read_csv(ctx.run_dir / "output" / "ppi" / "ppi_consensus_patch.csv")

    assert patches[0]["evidence_class"] == "DISPERSED"
    assert float(patches[0]["centroid_spread_angstrom"]) > 20.0


def test_consensus_csv_is_deterministic(tmp_path):
    root = copy_task7_inputs(tmp_path)
    ctx_a, _report_a = run_consensus(tmp_path, root)
    ctx_b, _report_b = run_consensus(tmp_path, root)
    bytes_a = (ctx_a.run_dir / "output" / "ppi" / "ppi_consensus_patch.csv").read_bytes()
    bytes_b = (ctx_b.run_dir / "output" / "ppi" / "ppi_consensus_patch.csv").read_bytes()

    assert bytes_a == bytes_b


def test_no_output_escapes_run_dir(tmp_path):
    root = copy_task7_inputs(tmp_path)
    ctx, _report = run_consensus(tmp_path, root)
    for path in [
        ctx.manifest_dir / "ppi_consensus_input_manifest.json",
        ctx.manifest_dir / "ppi_consensus_qc_report.json",
        ctx.qc_dir / "ppi_consensus_schema_audit.csv",
        ctx.qc_dir / "ppi_residue_parsing_audit.csv",
        ctx.qc_dir / "ppi_tail_noise_audit.csv",
        ctx.qc_dir / "ppi_active_face_audit.csv",
        ctx.qc_dir / "ppi_membrane_atp_audit.csv",
        ctx.qc_dir / "ppi_convergence_audit.csv",
        ctx.run_dir / "output" / "ppi" / "ppi_consensus_patch.csv",
        ctx.reports_dir / "ppi_consensus_report.txt",
    ]:
        assert_under_run_dir(ctx, path)

