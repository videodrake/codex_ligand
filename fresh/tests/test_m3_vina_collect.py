import csv
import json
from pathlib import Path

from egfr_myo1d.compound.pdbqt_parse import parse_vina_log, parse_vina_pdbqt
from egfr_myo1d.compound.vina_collect import run_m3_vina_collect
from egfr_myo1d.compound.vina_jobs import JOB_FIELDS
from egfr_myo1d.core.logging_utils import initialize_logs
from egfr_myo1d.core.run_context import RunContext


POSE_MULTI = """MODEL 1
REMARK VINA RESULT: -7.500 0.000 0.000
ATOM      1  C   LIG A   1       0.000   0.000   0.000  0.00  0.00     0.000 C
ATOM      2  N   LIG A   1       2.000   0.000   0.000  0.00  0.00     0.000 N
ENDMDL
MODEL 2
REMARK VINA RESULT: -6.250 0.000 0.000
ATOM      1  C   LIG A   1       0.000   2.000   0.000  0.00  0.00     0.000 C
ATOM      2  H   LIG A   1       0.000   4.000   0.000  0.00  0.00     0.000 H
ENDMDL
"""

POSE_SINGLE = """REMARK VINA RESULT: -5.000 0.000 0.000
ATOM      1  C   LIG A   1       1.000   1.000   1.000  0.00  0.00     0.000 C
ATOM      2  O   LIG A   1       3.000   1.000   1.000  0.00  0.00     0.000 O
END
"""

LOG_TEXT = """-----+------------+----------+----------
mode | affinity | dist from best mode
     | (kcal/mol) | rmsd l.b.| rmsd u.b.
-----+------------+----------+----------
   1       -7.500      0.000      0.000
   2       -6.250      1.000      2.000
"""


def make_ctx(tmp_path, run_id="m3_vina_collect"):
    repo = tmp_path / "repo"
    fresh = repo / "fresh"
    run_dir = fresh / "runs" / run_id
    ctx = RunContext(
        repo_root=repo,
        fresh_root=fresh,
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
    (repo / ".gitignore").write_text("", encoding="utf-8")
    private = fresh / "data" / "private"
    private.mkdir(parents=True, exist_ok=True)
    (private / "compound_id_map.csv").write_text("public_id,internal_id,notes\nCpd-A,SECRET_A,test\n", encoding="utf-8")
    return ctx


def write_csv(path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def read_csv(path):
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def job_row(ctx, job_id="job_1", profile="mini", cid="Cpd-A", state="EGFR_160-185", box="box_1", repeat="1", allowed=True, collect_allowed=True, output=None, log=None, stdout=None, stderr=None):
    output = output or (ctx.run_dir / "phase3_compounds" / "docking_outputs" / "focused_pocket_first" / profile / job_id / "pose_out.pdbqt")
    log = log or output.with_name("vina.log")
    stdout = stdout or (ctx.jobs_log_dir / f"{job_id}.vina.stdout")
    stderr = stderr or (ctx.jobs_log_dir / f"{job_id}.vina.stderr")
    row = {field: "" for field in JOB_FIELDS}
    row.update(
        {
            "job_id": job_id,
            "run_id": ctx.run_id,
            "m2_run_id": "m2_run",
            "profile": profile,
            "mode": "generate",
            "job_scope": "focused_pocket_first",
            "chunk_id": "mini_EGFR_160-185_fam_1_chunk01",
            "node": "node04",
            "compound_public_id": cid,
            "state_id": state,
            "state_role": "reference" if state == "3GT8_raw" else "primary",
            "pocket_family_id": "fam_1",
            "box_id": box,
            "protomer_id": "A",
            "vina_repeat_id": repeat,
            "num_modes": "2",
            "output_pdbqt_file": ctx.relative_to_repo(output) if isinstance(output, Path) else output,
            "vina_log_file": ctx.relative_to_repo(log) if isinstance(log, Path) else log,
            "stdout_file": ctx.relative_to_repo(stdout) if isinstance(stdout, Path) else stdout,
            "stderr_file": ctx.relative_to_repo(stderr) if isinstance(stderr, Path) else stderr,
            "docking_status": "PLANNED",
            "allowed_for_execution": "true" if allowed else "false",
            "allowed_for_collection": "true" if collect_allowed else "false",
        }
    )
    return row, output, log, stdout, stderr


def write_manifest(ctx, rows, plan_status="PASS"):
    manifest = ctx.run_dir / "phase3_compounds" / "manifests" / "vina_job_manifest.csv"
    write_csv(manifest, JOB_FIELDS, rows)
    qc = ctx.run_dir / "phase3_compounds" / "qc"
    qc.mkdir(parents=True, exist_ok=True)
    (qc / "vina_job_plan_qc.json").write_text(json.dumps({"overall_status": plan_status, "counts": {"planned_vina_commands": len(rows), "planned_chunks": 1}}), encoding="utf-8")
    return manifest


def write_output(output, log, stdout, stderr, pose_text=POSE_MULTI, log_text=LOG_TEXT):
    output.parent.mkdir(parents=True, exist_ok=True)
    log.parent.mkdir(parents=True, exist_ok=True)
    stdout.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(pose_text, encoding="utf-8")
    log.write_text(log_text, encoding="utf-8")
    stdout.write_text("done\n", encoding="utf-8")
    stderr.write_text("", encoding="utf-8")


def test_pdbqt_parser_multimodel_single_model_and_log_affinity(tmp_path):
    multi = tmp_path / "multi.pdbqt"
    single = tmp_path / "single.pdbqt"
    log = tmp_path / "vina.log"
    multi.write_text(POSE_MULTI, encoding="utf-8")
    single.write_text(POSE_SINGLE, encoding="utf-8")
    log.write_text(LOG_TEXT, encoding="utf-8")

    poses = parse_vina_pdbqt(multi)
    one = parse_vina_pdbqt(single)
    affinities = parse_vina_log(log)

    assert len(poses) == 2
    assert poses[0].pose_rank == 1
    assert poses[0].vina_affinity_kcal_mol == -7.5
    assert poses[0].pose_center_x == 1.0
    assert poses[1].pose_heavy_atom_count == 1
    assert len(one) == 1
    assert affinities == {1: -7.5, 2: -6.25}


def test_dry_run_missing_manifest_warns_and_collect_missing_manifest_fails(tmp_path):
    ctx = make_ctx(tmp_path)

    dry = run_m3_vina_collect(ctx, m2_run_id="m2_run", mode="dry-run", profile="mini", write_empty_raw_table=True)
    fail = run_m3_vina_collect(ctx, m2_run_id="m2_run", mode="collect", profile="mini")

    assert dry.status == "WARN"
    assert dry.qc_json.is_file()
    assert dry.raw_pose_csv.is_file()
    assert fail.status == "FAIL"
    assert any("vina_job_manifest" in item for item in fail.blockers)


def test_collect_synthetic_manifest_writes_raw_and_completion_tables(tmp_path):
    ctx = make_ctx(tmp_path)
    row, output, log, stdout, stderr = job_row(ctx)
    write_output(output, log, stdout, stderr)
    manifest = write_manifest(ctx, [row])

    result = run_m3_vina_collect(ctx, m2_run_id="m2_run", job_manifest=manifest, profile="mini")
    raw = read_csv(result.raw_pose_csv)
    completion = read_csv(result.completion_csv)

    assert result.status == "PASS"
    assert result.m3_t7_attribution_ready is True
    assert len(raw) == 2
    assert len(completion) == 1
    assert completion[0]["completion_status"] == "COMPLETED"
    assert completion[0]["output_pdbqt_sha256"]
    assert completion[0]["vina_log_sha256"]
    assert raw[0]["compound_public_id"] == "Cpd-A"
    assert raw[0]["output_pdbqt_sha256"] == completion[0]["output_pdbqt_sha256"]
    assert raw[0]["vina_log_sha256"] == completion[0]["vina_log_sha256"]
    assert raw[0]["allowed_for_attribution"] == "true"


def test_schema_duplicate_nonpublic_and_path_validation_fail(tmp_path):
    ctx = make_ctx(tmp_path / "schema")
    row, output, log, stdout, stderr = job_row(ctx)
    write_output(output, log, stdout, stderr)
    bad_manifest = ctx.run_dir / "phase3_compounds" / "manifests" / "vina_job_manifest.csv"
    write_csv(bad_manifest, ["job_id"], [{"job_id": "x"}])
    assert run_m3_vina_collect(ctx, m2_run_id="m2_run", profile="mini").status == "FAIL"

    ctx2 = make_ctx(tmp_path / "dupe")
    row1, output1, log1, stdout1, stderr1 = job_row(ctx2, job_id="same")
    row2, output2, log2, stdout2, stderr2 = job_row(ctx2, job_id="same", box="box_2")
    write_output(output1, log1, stdout1, stderr1)
    write_output(output2, log2, stdout2, stderr2)
    write_manifest(ctx2, [row1, row2])
    assert run_m3_vina_collect(ctx2, m2_run_id="m2_run", profile="mini").status == "FAIL"

    ctx3 = make_ctx(tmp_path / "badid")
    row, output, log, stdout, stderr = job_row(ctx3, cid="SECRET")
    write_output(output, log, stdout, stderr)
    write_manifest(ctx3, [row])
    assert run_m3_vina_collect(ctx3, m2_run_id="m2_run", profile="mini").status == "FAIL"

    ctx4 = make_ctx(tmp_path / "path")
    row, output, log, stdout, stderr = job_row(ctx4, output=ctx4.repo_root / "outside.pdbqt")
    write_output(output, log, stdout, stderr)
    write_manifest(ctx4, [row])
    assert run_m3_vina_collect(ctx4, m2_run_id="m2_run", profile="mini").status == "FAIL"


def test_missing_empty_and_parse_failed_jobs_are_recorded(tmp_path):
    ctx = make_ctx(tmp_path)
    rows = []
    row1, output1, log1, stdout1, stderr1 = job_row(ctx, job_id="missing_output")
    log1.parent.mkdir(parents=True, exist_ok=True)
    log1.write_text(LOG_TEXT, encoding="utf-8")
    stdout1.parent.mkdir(parents=True, exist_ok=True)
    stdout1.write_text("", encoding="utf-8")
    stderr1.write_text("", encoding="utf-8")
    rows.append(row1)
    row2, output2, log2, stdout2, stderr2 = job_row(ctx, job_id="missing_log", box="box_2")
    output2.parent.mkdir(parents=True, exist_ok=True)
    output2.write_text(POSE_SINGLE, encoding="utf-8")
    rows.append(row2)
    row3, output3, log3, stdout3, stderr3 = job_row(ctx, job_id="empty_output", box="box_3")
    output3.parent.mkdir(parents=True, exist_ok=True)
    output3.write_text("", encoding="utf-8")
    log3.write_text(LOG_TEXT, encoding="utf-8")
    rows.append(row3)
    row4, output4, log4, stdout4, stderr4 = job_row(ctx, job_id="parse_fail", box="box_4")
    write_output(output4, log4, stdout4, stderr4, pose_text="REMARK VINA RESULT: -1.0\n")
    rows.append(row4)
    manifest = write_manifest(ctx, rows)

    result = run_m3_vina_collect(ctx, m2_run_id="m2_run", job_manifest=manifest, profile="mini", allow_partial=True)
    completion = {row["job_id"]: row["completion_status"] for row in read_csv(result.completion_csv)}

    assert result.status == "FAIL"
    assert completion["missing_output"] == "MISSING_OUTPUT"
    assert completion["missing_log"] == "MISSING_LOG"
    assert completion["empty_output"] == "OUTPUT_EMPTY"
    assert completion["parse_fail"] == "PARSE_FAIL"


def test_allow_partial_warns_when_some_poses_parse(tmp_path):
    ctx = make_ctx(tmp_path)
    good, output, log, stdout, stderr = job_row(ctx, job_id="good")
    write_output(output, log, stdout, stderr)
    missing, _, missing_log, _, _ = job_row(ctx, job_id="missing", box="box_2")
    missing_log.parent.mkdir(parents=True, exist_ok=True)
    missing_log.write_text(LOG_TEXT, encoding="utf-8")
    manifest = write_manifest(ctx, [good, missing])

    blocked = run_m3_vina_collect(ctx, m2_run_id="m2_run", job_manifest=manifest, profile="mini")
    partial = run_m3_vina_collect(ctx, m2_run_id="m2_run", job_manifest=manifest, profile="mini", allow_partial=True, force=True)

    assert blocked.status == "FAIL"
    assert partial.status == "WARN"
    assert partial.collection["parsed_pose_rows"] == 2


def test_log_error_classified_failed_by_log(tmp_path):
    ctx = make_ctx(tmp_path)
    row, output, log, stdout, stderr = job_row(ctx)
    log.parent.mkdir(parents=True, exist_ok=True)
    log.write_text("ERROR Cannot open receptor\n", encoding="utf-8")
    manifest = write_manifest(ctx, [row])

    result = run_m3_vina_collect(ctx, m2_run_id="m2_run", job_manifest=manifest, profile="mini", allow_partial=True)
    completion = read_csv(result.completion_csv)[0]

    assert completion["completion_status"] == "FAILED_BY_LOG"
    assert completion["failure_reason_code"] == "ERROR"


def test_affinity_disagreement_is_parse_warn(tmp_path):
    ctx = make_ctx(tmp_path)
    row, output, log, stdout, stderr = job_row(ctx)
    write_output(output, log, stdout, stderr, log_text=LOG_TEXT.replace("-7.500", "-8.000", 1))
    manifest = write_manifest(ctx, [row])

    result = run_m3_vina_collect(ctx, m2_run_id="m2_run", job_manifest=manifest, profile="mini")
    completion = read_csv(result.completion_csv)[0]
    raw = read_csv(result.raw_pose_csv)

    assert result.status == "WARN"
    assert completion["completion_status"] == "COMPLETED_PARSE_WARN"
    assert "disagreement" in raw[0]["parse_notes"]


def test_reference_only_collection_warns_and_can_be_allowed(tmp_path):
    ctx = make_ctx(tmp_path)
    row, output, log, stdout, stderr = job_row(ctx, state="3GT8_raw")
    write_output(output, log, stdout, stderr)
    manifest = write_manifest(ctx, [row])

    warn = run_m3_vina_collect(ctx, m2_run_id="m2_run", job_manifest=manifest, profile="mini")
    allowed = run_m3_vina_collect(ctx, m2_run_id="m2_run", job_manifest=manifest, profile="mini", allow_reference_only=True, force=True)

    assert warn.status == "WARN"
    assert allowed.status == "WARN"
    assert allowed.m3_t7_attribution_ready is False
    assert read_csv(allowed.raw_pose_csv)[0]["state_role"] == "reference"


def test_not_expected_rows_are_recorded_not_parsed(tmp_path):
    ctx = make_ctx(tmp_path)
    row, output, log, stdout, stderr = job_row(ctx, allowed=False)
    write_output(output, log, stdout, stderr)
    manifest = write_manifest(ctx, [row])

    result = run_m3_vina_collect(ctx, m2_run_id="m2_run", job_manifest=manifest, profile="mini", mode="dry-run")
    completion = read_csv(result.completion_csv)[0]

    assert completion["completion_status"] == "SKIPPED_NOT_ALLOWED"
    assert result.collection["not_expected_jobs"] == 1


def test_collection_requires_collection_allow_flag(tmp_path):
    ctx = make_ctx(tmp_path)
    row, output, log, stdout, stderr = job_row(ctx, collect_allowed=False)
    write_output(output, log, stdout, stderr)
    manifest = write_manifest(ctx, [row])

    result = run_m3_vina_collect(ctx, m2_run_id="m2_run", job_manifest=manifest, profile="mini")
    completion = read_csv(result.completion_csv)[0]

    assert result.status == "FAIL"
    assert completion["completion_status"] == "SKIPPED_COLLECTION_NOT_ALLOWED"
    assert result.collection["expected_jobs"] == 0
    assert result.collection["parsed_pose_rows"] == 0


def test_collect_ignores_header_only_downstream_and_prior_smiles_warning_text(tmp_path):
    ctx = make_ctx(tmp_path)
    row, output, log, stdout, stderr = job_row(ctx)
    write_output(output, log, stdout, stderr)
    manifest = write_manifest(ctx, [row])
    write_csv(ctx.run_dir / "phase3_compounds" / "tables" / "compound_pose_attribution.csv", ["job_id", "pose_rank"], [])
    (ctx.errors_dir / "error_summary.txt").write_text("previous blocker: SMILES printed in public outputs/logs\n", encoding="utf-8")

    result = run_m3_vina_collect(ctx, m2_run_id="m2_run", job_manifest=manifest, profile="mini")

    assert result.status == "PASS"
    assert result.collection["expected_jobs"] == 1
    assert result.collection["parsed_pose_rows"] == 2
    assert "forbidden downstream output created" not in result.blockers


def test_collect_requires_passed_m3_t5_plan(tmp_path):
    ctx = make_ctx(tmp_path)
    row, output, log, stdout, stderr = job_row(ctx)
    write_output(output, log, stdout, stderr)
    manifest = write_manifest(ctx, [row], plan_status="WARN")

    result = run_m3_vina_collect(ctx, m2_run_id="m2_run", job_manifest=manifest, profile="mini")

    assert result.status == "FAIL"
    assert any("M3-T5 job plan QC did not PASS" in item for item in result.blockers)
    assert result.m3_t7_attribution_ready is False


def test_no_coordinates_or_sensitive_tokens_in_public_outputs(tmp_path):
    ctx = make_ctx(tmp_path)
    row, output, log, stdout, stderr = job_row(ctx)
    write_output(output, log, stdout, stderr)
    manifest = write_manifest(ctx, [row])

    result = run_m3_vina_collect(ctx, m2_run_id="m2_run", job_manifest=manifest, profile="mini")
    public_text = "\n".join(
        path.read_text(encoding="utf-8", errors="ignore")
        for path in [result.raw_pose_csv, result.completion_csv, result.qc_csv, result.qc_json, result.report_md, result.phase3_log]
    )

    assert result.status == "PASS"
    assert "SECRET_A" not in public_text
    assert "SMILES" not in public_text
    assert "ATOM      1" not in public_text
    assert "final_m3_candidate_hypotheses.csv" not in public_text
    assert not (ctx.run_dir / "phase3_compounds" / "tables" / "compound_pose_attribution.csv").exists()
    assert not (ctx.run_dir / "phase3_compounds" / "qc" / "atp_migration_qc.csv").exists()


def test_deterministic_collection_and_status_records(tmp_path):
    ctx = make_ctx(tmp_path)
    row1, output1, log1, stdout1, stderr1 = job_row(ctx, job_id="job_a")
    row2, output2, log2, stdout2, stderr2 = job_row(ctx, job_id="job_b", box="box_2")
    write_output(output1, log1, stdout1, stderr1)
    write_output(output2, log2, stdout2, stderr2, pose_text=POSE_SINGLE)
    manifest = write_manifest(ctx, [row2, row1])

    first = run_m3_vina_collect(ctx, m2_run_id="m2_run", job_manifest=manifest, profile="mini")
    first_jobs = [row["job_id"] for row in read_csv(first.raw_pose_csv)]
    second = run_m3_vina_collect(ctx, m2_run_id="m2_run", job_manifest=manifest, profile="mini", force=True)
    second_jobs = [row["job_id"] for row in read_csv(second.raw_pose_csv)]
    summary = json.loads(second.qc_json.read_text(encoding="utf-8"))

    assert first_jobs == second_jobs
    assert summary["collection"]["parsed_pose_rows"] == len(read_csv(second.raw_pose_csv))
    assert "M3-T6" in (ctx.logs_dir / "phase_status.jsonl").read_text(encoding="utf-8")
    assert "vina_output_collection" in (ctx.logs_dir / "job_status.jsonl").read_text(encoding="utf-8")


def test_failed_jobs_csv_updated_for_missing_output(tmp_path):
    ctx = make_ctx(tmp_path)
    row, _, log, stdout, stderr = job_row(ctx)
    log.parent.mkdir(parents=True, exist_ok=True)
    log.write_text(LOG_TEXT, encoding="utf-8")
    stdout.parent.mkdir(parents=True, exist_ok=True)
    stdout.write_text("", encoding="utf-8")
    stderr.write_text("", encoding="utf-8")
    manifest = write_manifest(ctx, [row])

    result = run_m3_vina_collect(ctx, m2_run_id="m2_run", job_manifest=manifest, profile="mini", allow_partial=True)
    failed_jobs = (ctx.errors_dir / "failed_jobs.csv").read_text(encoding="utf-8")

    assert result.status == "FAIL"
    assert "MISSING_OUTPUT" in failed_jobs
