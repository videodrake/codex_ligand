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
from egfr_myo1d.pocket.m2_export import (
    ACCEPTED_M3_FIELDS,
    BOX_FIELDS,
    EVIDENCE_FIELDS,
    STATE_MAP_FIELDS,
    TRANSITION_GATE_FIELDS,
    export_m2_results,
    validate_box_row,
)
from egfr_myo1d.pocket.pocket_gate_apply import apply_pocket_gates

from test_m2_phase7_pocket_gates import make_complete_gate_fixture


REPO_ROOT = Path(__file__).resolve().parents[2]


def unique_run_id(prefix="pytest_m2_phase8"):
    return "{0}_{1}".format(prefix, uuid.uuid4().hex[:12])


def make_tmp_run_context(tmp_path):
    run_id = unique_run_id()
    run_dir = tmp_path / run_id
    ctx = RunContext(
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
    initialize_logs(ctx)
    return ctx


def read_csv(path):
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def read_header(path):
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        return next(csv.reader(handle))


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def write_csv_file(path, fieldnames, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def make_export_fixture(ctx):
    make_complete_gate_fixture(ctx)
    gate_report = apply_pocket_gates(ctx, synthetic_fixture=True)
    assert gate_report.status == "PASS"
    return export_m2_results(ctx, synthetic_fixture=True)


def test_m2_8_cli_smoke_runs_without_science_engines():
    run_id = unique_run_id("m2_8_cli_smoke")
    ctx = RunContext.create(run_id, "smoke_input", repo_root=REPO_ROOT)
    make_complete_gate_fixture(ctx)
    apply_pocket_gates(ctx, synthetic_fixture=True)

    proc = subprocess.run(
        [
            sys.executable,
            "-m",
            "egfr_myo1d.cli",
            "export-m2-results",
            "--run-id",
            run_id,
            "--synthetic-fixture",
            "true",
        ],
        cwd=str(REPO_ROOT),
        env=_cli_env(),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    assert proc.returncode == 0
    assert b"export-m2-results" in proc.stdout
    assert (ctx.run_dir / "phase2_pockets" / "export_for_m3" / "accepted_pockets_for_m3.csv").is_file()
    assert not (ctx.run_dir / "phase3_compounds").exists()


def test_m2_8_export_schema_regression(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    report = make_export_fixture(ctx)

    assert read_header(report.export_accepted_csv) == ACCEPTED_M3_FIELDS
    assert read_header(report.export_boxes_csv) == BOX_FIELDS
    assert read_header(report.export_state_map_csv) == STATE_MAP_FIELDS
    assert read_header(report.final_evidence_csv) == EVIDENCE_FIELDS
    assert read_header(report.transition_gate_csv) == TRANSITION_GATE_FIELDS


def test_m2_8_residue_set_json_preserves_state_chain_protomer(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    report = make_export_fixture(ctx)
    payload = read_json(report.export_residue_sets_json)
    pockets = {row["pocket_family_id"]: row for row in payload["pockets"]}

    primary = pockets["fam_primary"]
    assert primary["export_tier"] == "primary"
    state = primary["residue_sets_by_state"]["EGFR_160-185"]
    residue = state["protomer_residue_sets"]["A"][0]
    assert residue["chain_id"] == "A"
    assert residue["uniprot_residue"] == "900"
    assert residue["mapping_status"] == "PASS"


def test_m2_8_primary_secondary_and_reference_export_rules(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    report = make_export_fixture(ctx)
    rows = {row["pocket_family_id"]: row for row in read_csv(report.export_accepted_csv)}
    gate_rows = {row["pocket_family_id"]: row for row in read_csv(ctx.run_dir / "phase2_pockets" / "gated" / "pocket_gate_qc.csv")}

    assert rows["fam_primary"]["export_tier"] == "primary"
    assert rows["fam_primary"]["m3_primary_allowed"] == "True"
    assert rows["fam_secondary"]["export_tier"] == "secondary_review"
    assert rows["fam_secondary"]["m3_primary_allowed"] == "False"
    assert "fam_reference" not in rows
    assert gate_rows["fam_reference"]["gate_class"] == "reference_only"


def test_m2_8_rejects_are_not_rescued_by_soft_scores(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    report = make_export_fixture(ctx)
    exported_ids = {row["pocket_family_id"] for row in read_csv(report.export_accepted_csv)}

    assert "fam_atp_overlap" not in exported_ids
    assert "fam_atp_centroid" not in exported_ids
    assert "fam_low" not in exported_ids
    assert "fam_membrane" not in exported_ids
    assert "fam_dimer" not in exported_ids
    assert "fam_mapping" not in exported_ids
    evidence = {row["pocket_family_id"]: row for row in read_csv(report.final_evidence_csv)}
    assert evidence["fam_atp_overlap"]["export_tier"] == "not_exported"


def test_m2_8_zero_accepted_pockets_writes_empty_exports_and_blocks_m3(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    report = export_m2_results(ctx, synthetic_fixture=True)

    assert report.status == "PASS_WITH_WARNINGS"
    assert report.m3_docking_allowed is False
    assert read_csv(report.export_accepted_csv) == []
    text = report.final_summary_md.read_text(encoding="utf-8")
    assert "Milestone 3 docking should not proceed without manual review" in text


def test_m2_8_box_validation_rules():
    row = {
        "pocket_family_id": "fam",
        "state_id": "EGFR_160-185",
        "protomer_id": "A",
        "box_center_x": "1.0",
        "box_center_y": "2.0",
        "box_center_z": "3.0",
        "box_size_x": "12.0",
        "box_size_y": "12.0",
        "box_size_z": "12.0",
        "non_atp_pass": "True",
        "reference_only": "False",
        "m3_primary_allowed": "True",
    }
    assert validate_box_row(row)[0] == "PASS"
    bad_center = dict(row, box_center_x="")
    assert validate_box_row(bad_center)[0] == "FAIL"
    bad_size = dict(row, box_size_x="0")
    assert validate_box_row(bad_size)[0] == "FAIL"


def test_m2_8_transition_gate_fails_missing_inputs_in_production(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    report = export_m2_results(ctx, mode="production", synthetic_fixture=False)
    rows = {row["gate_id"]: row for row in read_csv(report.transition_gate_csv)}

    assert report.status == "FAIL"
    assert rows["G3"]["status"] == "FAIL"
    assert rows["G4"]["status"] == "FAIL"
    assert rows["G6"]["status"] == "FAIL"
    assert report.m3_docking_allowed is False


def test_m2_8_path_safety_rejects_inputs_outside_fresh_or_run_dir(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    outside = tmp_path / "outside_ppi.csv"
    outside.write_text("x\n", encoding="utf-8")

    with pytest.raises(RunContextError):
        export_m2_results(ctx, ppi_consensus=outside, synthetic_fixture=True)


def test_m2_8_non_goals_no_m3_or_engine_outputs(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    report = make_export_fixture(ctx)
    payload = read_json(report.final_status_json)

    assert payload["vina_executed"] is False
    assert payload["ligand_preparation_executed"] is False
    assert payload["receptor_pdbqt_preparation_executed"] is False
    assert payload["pyrosetta_docking_or_relaxation_executed"] is False
    assert payload["lightdock_executed"] is False
    assert payload["p2rank_executed"] is False
    assert payload["compound_scoring_executed"] is False
    assert payload["candidate_nomination_executed"] is False
    assert payload["cleanup_mode"] == "report_only"
    assert not (ctx.run_dir / "phase3_compounds").exists()
    assert not (ctx.run_dir / "phase4_scoring").exists()
    cleanup = read_json(report.cleanup_report_json)
    assert cleanup["deleted_files"] == []


def _cli_env():
    env = os.environ.copy()
    env["PYTHONPATH"] = str(REPO_ROOT / "fresh" / "src") + os.pathsep + env.get("PYTHONPATH", "")
    return env
