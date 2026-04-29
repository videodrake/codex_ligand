import builtins
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
from egfr_myo1d.validation.pocket_discovery import plan_pocket_discovery


REPO_ROOT = Path(__file__).resolve().parents[2]
TASK8_INPUTS = REPO_ROOT / "fresh" / "tests" / "fixtures" / "task8_pocket_planning"
HEADER = [
    "ppi_patch_id",
    "receptor_id",
    "receptor_state",
    "protomer_id",
    "support_methods",
    "support_seeds",
    "support_replicates",
    "n_accepted_poses",
    "egfr_consensus_residues",
    "egfr_contact_centroid_x",
    "egfr_contact_centroid_y",
    "egfr_contact_centroid_z",
    "centroid_spread_angstrom",
    "residue_occupancy_summary",
    "active_face_support_fraction",
    "sheet12_support_fraction",
    "tail_noise_fraction",
    "terminal_artifact_fraction",
    "atp_overlap_fraction",
    "membrane_proximal_fraction",
    "membrane_side_classes",
    "evidence_class",
    "warnings",
]


def unique_run_id(prefix="pytest_task8"):
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


def copy_task8_inputs(tmp_path):
    root = tmp_path / "task8_inputs"
    shutil.copytree(TASK8_INPUTS, root)
    return root


def base_patch(patch_id="ppi_patch_001", evidence_class="CONVERGENT_PATCH"):
    return {
        "ppi_patch_id": patch_id,
        "receptor_id": "EGFR_synthetic_inactive_dimer_AB",
        "receptor_state": "EGFR_fixture_state",
        "protomer_id": "A",
        "support_methods": "lightdock_gso_ppi;pyrosetta_global_ppi",
        "support_seeds": "0;1",
        "support_replicates": "1;2",
        "n_accepted_poses": "4",
        "egfr_consensus_residues": "A:742;A:746;A:923",
        "egfr_contact_centroid_x": "1.100",
        "egfr_contact_centroid_y": "2.050",
        "egfr_contact_centroid_z": "3.050",
        "centroid_spread_angstrom": "0.430",
        "residue_occupancy_summary": "A:742=4;A:746=4;A:923=3",
        "active_face_support_fraction": "1.000",
        "sheet12_support_fraction": "1.000",
        "tail_noise_fraction": "0.000",
        "terminal_artifact_fraction": "0.000",
        "atp_overlap_fraction": "0.000",
        "membrane_proximal_fraction": "0.000",
        "membrane_side_classes": "solvent_exposed",
        "evidence_class": evidence_class,
        "warnings": "",
    }


def write_consensus(root, rows, header=None):
    header = header or HEADER
    path = root / "ppi_consensus_patch.csv"
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=header, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    return path


def run_plan(tmp_path, root):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    report = plan_pocket_discovery(ctx, "smoke_env", "codex_dev", root, input_kind="synthetic_fixture")
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


def test_cli_help_includes_plan_pocket_discovery():
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
    assert "plan-pocket-discovery" in proc.stdout


def test_accepted_ppi_patch_produces_pocket_selection_plan(tmp_path):
    root = copy_task8_inputs(tmp_path)
    ctx, report = run_plan(tmp_path, root)
    pocket_csv = ctx.run_dir / "pockets" / "egfr_myo1d_ppi_adjacent_pockets.csv"
    alias_csv = ctx.run_dir / "pockets" / "egfr_myo1d_ppi_guided_pocket_plan_records.csv"
    plan_json = ctx.run_dir / "pockets" / "pocket_discovery_plan.json"
    rows = read_csv(pocket_csv)
    alias_rows = read_csv(alias_csv)
    plan = read_json(plan_json)

    assert report["status"] == "PASS"
    assert report["accepted_ppi_patch_count"] == 1
    assert pocket_csv.is_file()
    assert alias_csv.is_file()
    assert plan_json.is_file()
    assert rows == alias_rows
    assert plan["record_semantics"] == "ppi_guided_pocket_plan"
    assert plan["pocket_detector_runtime_executed"] is False
    assert plan["detected_pocket_record"] is False
    assert plan["compound_docking_or_scoring_executed"] is False
    assert any(row["ppi_patch_id"] == "ppi_patch_001" and row["future_discovery_allowed"] == "True" for row in rows)


def test_no_accepted_ppi_patch_produces_no_go(tmp_path):
    root = tmp_path / "no_go"
    root.mkdir()
    row = base_patch(evidence_class="TAIL_DOMINATED")
    row["tail_noise_fraction"] = "1.000"
    row["terminal_artifact_fraction"] = "0.500"
    write_consensus(root, [row])
    ctx, report = run_plan(tmp_path, root)
    rows = read_csv(ctx.run_dir / "pockets" / "egfr_myo1d_ppi_adjacent_pockets.csv")

    assert report["status"] == "NO_GO"
    assert report["planned_docking_job_count"] == 0
    assert rows == []


def test_schema_audit_written_for_pass_and_missing_schema(tmp_path):
    root = copy_task8_inputs(tmp_path)
    ctx, report = run_plan(tmp_path, root)
    schema_rows = read_csv(ctx.qc_dir / "task7_consensus_schema_audit.csv")

    assert report["status"] == "PASS"
    assert schema_rows
    assert all(row["status"] == "PASS" for row in schema_rows)

    bad_root = tmp_path / "bad_schema"
    bad_root.mkdir()
    row = base_patch()
    bad_header = [column for column in HEADER if column != "evidence_class"]
    write_consensus(bad_root, [row], header=bad_header)
    bad_ctx, bad_report = run_plan(tmp_path, bad_root)
    bad_schema_rows = read_csv(bad_ctx.qc_dir / "task7_consensus_schema_audit.csv")

    assert bad_report["status"] == "NO_GO"
    assert any(row["field"] == "evidence_class" and row["status"] == "FAIL" for row in bad_schema_rows)
    assert any(blocker["code"] == "invalid_task7_consensus_schema" for blocker in bad_report["blockers"])


def test_atp_overlap_pocket_is_flagged(tmp_path):
    root = tmp_path / "atp"
    root.mkdir()
    row = base_patch()
    row["atp_overlap_fraction"] = "1.000"
    write_consensus(root, [row])
    ctx, report = run_plan(tmp_path, root)
    atp_rows = read_csv(ctx.qc_dir / "atp_overlap_audit.csv")
    pocket_rows = read_csv(ctx.run_dir / "pockets" / "egfr_myo1d_ppi_adjacent_pockets.csv")

    assert report["status"] == "NO_GO"
    assert atp_rows[0]["atp_relation_class"] == "ATP-overlap"
    assert pocket_rows == []


def test_membrane_proximal_pocket_is_flagged(tmp_path):
    root = tmp_path / "membrane"
    root.mkdir()
    row = base_patch()
    row["membrane_proximal_fraction"] = "1.000"
    row["membrane_side_classes"] = "membrane_proximal"
    write_consensus(root, [row])
    ctx, report = run_plan(tmp_path, root)
    membrane_rows = read_csv(ctx.qc_dir / "membrane_accessibility_audit.csv")

    assert report["status"] == "NO_GO"
    assert membrane_rows[0]["membrane_accessibility_class"] == "membrane-proximal-risk"


def test_protomer_and_residue_numbering_are_preserved(tmp_path):
    root = copy_task8_inputs(tmp_path)
    ctx, _report = run_plan(tmp_path, root)
    rows = read_csv(ctx.run_dir / "pockets" / "egfr_myo1d_ppi_adjacent_pockets.csv")
    dimer_rows = read_csv(ctx.qc_dir / "dimer_accessibility_audit.csv")

    allowed = [row for row in rows if row["future_discovery_allowed"] == "True"][0]
    assert allowed["protomer_id"] == "A"
    assert allowed["egfr_residues"] == "A:742;A:746;A:923"
    assert any(row["protomer_id"] == "A" and row["dimer_accessibility_class"] == "dimer-accessible" for row in dimer_rows)


def test_all_outputs_remain_under_run_dir(tmp_path):
    root = copy_task8_inputs(tmp_path)
    ctx, _report = run_plan(tmp_path, root)
    required = [
        ctx.run_dir / "pockets" / "egfr_myo1d_ppi_adjacent_pockets.csv",
        ctx.run_dir / "pockets" / "egfr_myo1d_ppi_guided_pocket_plan_records.csv",
        ctx.run_dir / "pockets" / "pocket_discovery_plan.json",
        ctx.qc_dir / "task7_consensus_schema_audit.csv",
        ctx.qc_dir / "pocket_selection_audit.csv",
        ctx.qc_dir / "atp_overlap_audit.csv",
        ctx.qc_dir / "membrane_accessibility_audit.csv",
        ctx.qc_dir / "dimer_accessibility_audit.csv",
        ctx.manifest_dir / "pocket_discovery_manifest.json",
        ctx.reports_dir / "pocket_discovery_summary.md",
    ]
    for path in required:
        assert path.is_file()
        assert_under_run_dir(ctx, path)


def test_no_external_scientific_runtime_tool_is_executed(tmp_path, monkeypatch):
    root = copy_task8_inputs(tmp_path)
    blocked = {"pyrosetta", "lightdock", "vina", "fpocket", "p2rank", "prank"}
    original_import = builtins.__import__

    def guarded_import(name, globals=None, locals=None, fromlist=(), level=0):
        if name.split(".")[0].lower() in blocked:
            raise AssertionError("Task 8 must not import {}".format(name))
        return original_import(name, globals, locals, fromlist, level)

    monkeypatch.setattr(builtins, "__import__", guarded_import)
    ctx, report = run_plan(tmp_path, root)
    plan = read_json(ctx.run_dir / "pockets" / "pocket_discovery_plan.json")

    assert report["pocket_detector_runtime_executed"] is False
    assert plan["detected_pocket_record"] is False
    assert plan["record_semantics"] == "ppi_guided_pocket_plan"
    assert plan["compound_docking_or_scoring_executed"] is False
    assert plan["candidate_nomination_executed"] is False
