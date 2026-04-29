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
from egfr_myo1d.validation.pocket_candidate_prioritization import prioritize_pocket_candidates


REPO_ROOT = Path(__file__).resolve().parents[2]
TASK9_INPUTS = REPO_ROOT / "fresh" / "tests" / "fixtures" / "task9_pocket_candidates"


def unique_run_id(prefix="pytest_task9"):
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


def copy_task9_inputs(tmp_path):
    root = tmp_path / "task9_inputs"
    shutil.copytree(TASK9_INPUTS, root)
    return root


def run_prioritization(tmp_path, root):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    report = prioritize_pocket_candidates(
        ctx,
        "smoke_env",
        "codex_dev",
        root,
        root / "pocket_discovery_plan.json",
        root / "pocket_detector_candidates.csv",
        root / "ppi_consensus_patch.csv",
    )
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


def test_cli_help_includes_prioritize_pocket_candidates():
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
    assert "prioritize-pocket-candidates" in proc.stdout


def test_valid_ppi_rim_non_atp_candidate_carries_forward(tmp_path):
    root = copy_task9_inputs(tmp_path)
    ctx, report = run_prioritization(tmp_path, root)
    rows = read_csv(ctx.run_dir / "pockets" / "egfr_myo1d_prioritized_pocket_candidates.csv")
    strong = [row for row in rows if row["candidate_pocket_id"] == "pocket_strong_001"][0]

    assert report["status"] == "PASS_WITH_WARNINGS"
    assert strong["final_priority_class"] in ("CARRY_FORWARD_STRONG", "CARRY_FORWARD_WEAK")
    assert strong["ppi_relationship_class"] in ("orthosteric_ppi_overlap", "ppi_rim")
    assert strong["atp_relationship_class"] == "non_ATP"
    assert strong["membrane_compatibility_class"] == "membrane_compatible"
    assert strong["dimer_accessibility_class"] == "dimer_accessible"


def test_atp_overlap_candidate_is_blocked_even_with_high_detector_score(tmp_path):
    root = copy_task9_inputs(tmp_path)
    ctx, _report = run_prioritization(tmp_path, root)
    rows = read_csv(ctx.run_dir / "pockets" / "egfr_myo1d_prioritized_pocket_candidates.csv")
    atp = [row for row in rows if row["candidate_pocket_id"] == "pocket_atp_001"][0]
    audit = read_csv(ctx.qc_dir / "atp_exclusion_audit.csv")

    assert atp["final_priority_class"] == "BLOCKED_ATP"
    assert atp["blocker"] == "BLOCKED_ATP"
    assert any(row["candidate_pocket_id"] == "pocket_atp_001" and row["status"] == "FAIL" for row in audit)


def test_membrane_and_dimer_blockers_are_preserved(tmp_path):
    root = copy_task9_inputs(tmp_path)
    ctx, _report = run_prioritization(tmp_path, root)
    rows = read_csv(ctx.run_dir / "pockets" / "egfr_myo1d_prioritized_pocket_candidates.csv")
    membrane = [row for row in rows if row["candidate_pocket_id"] == "pocket_membrane_001"][0]
    dimer = [row for row in rows if row["candidate_pocket_id"] == "pocket_dimer_001"][0]

    assert membrane["final_priority_class"] == "BLOCKED_MEMBRANE"
    assert membrane["membrane_compatibility_class"] == "membrane_blocked"
    assert dimer["final_priority_class"] == "BLOCKED_DIMER"
    assert dimer["dimer_accessibility_class"] == "dimer_buried"


def test_missing_candidate_schema_produces_no_go_and_schema_fail(tmp_path):
    root = copy_task9_inputs(tmp_path)
    candidate_path = root / "pocket_detector_candidates.csv"
    rows = read_csv(candidate_path)
    bad_header = [key for key in rows[0].keys() if key != "nearest_ppi_residue_distance_angstrom"]
    with candidate_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=bad_header, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

    ctx, report = run_prioritization(tmp_path, root)
    schema_rows = read_csv(ctx.qc_dir / "pocket_candidate_schema_audit.csv")

    assert report["status"] == "NO_GO"
    assert report["candidate_pocket_count"] == 0
    assert any(row["field"] == "nearest_ppi_residue_distance_angstrom" and row["status"] == "FAIL" for row in schema_rows)


def test_centroid_only_distance_does_not_promote_to_carry_forward(tmp_path):
    root = copy_task9_inputs(tmp_path)
    ctx, _report = run_prioritization(tmp_path, root)
    rows = read_csv(ctx.run_dir / "pockets" / "egfr_myo1d_prioritized_pocket_candidates.csv")
    centroid_only = [row for row in rows if row["candidate_pocket_id"] == "pocket_centroid_only_001"][0]

    assert centroid_only["ppi_relationship_class"] == "allosteric_far"
    assert centroid_only["final_priority_class"] == "EXPLORATORY_ONLY"
    assert centroid_only["carry_forward"] == "False"


def test_protomer_id_and_egfr_residue_numbering_are_preserved(tmp_path):
    root = copy_task9_inputs(tmp_path)
    ctx, _report = run_prioritization(tmp_path, root)
    rows = read_csv(ctx.run_dir / "pockets" / "egfr_myo1d_prioritized_pocket_candidates.csv")
    strong = [row for row in rows if row["candidate_pocket_id"] == "pocket_strong_001"][0]

    assert strong["protomer_id"] == "A"
    assert strong["pocket_residues"] == "A:742;A:746;A:923"
    assert strong["pocket_surface_residues"] == "A:742;A:746"


def test_task9_outputs_remain_under_run_dir_and_manifest_has_no_runtime_execution(tmp_path):
    root = copy_task9_inputs(tmp_path)
    ctx, report = run_prioritization(tmp_path, root)
    required = [
        ctx.run_dir / "pockets" / "egfr_myo1d_prioritized_pocket_candidates.csv",
        ctx.run_dir / "pockets" / "egfr_myo1d_pocket_candidate_families.csv",
        ctx.run_dir / "pockets" / "pocket_candidate_prioritization.json",
        ctx.qc_dir / "pocket_candidate_schema_audit.csv",
        ctx.qc_dir / "pocket_candidate_input_audit.csv",
        ctx.qc_dir / "ppi_pocket_proximity_audit.csv",
        ctx.qc_dir / "atp_exclusion_audit.csv",
        ctx.qc_dir / "membrane_compatibility_audit.csv",
        ctx.qc_dir / "dimer_accessibility_audit.csv",
        ctx.qc_dir / "pocket_candidate_blockers.csv",
        ctx.manifest_dir / "pocket_candidate_prioritization_manifest.json",
        ctx.reports_dir / "task9_pocket_candidate_prioritization_summary.md",
    ]
    manifest = read_json(ctx.manifest_dir / "pocket_candidate_prioritization_manifest.json")

    for path in required:
        assert path.is_file()
        assert_under_run_dir(ctx, path)
    assert report["planned_compound_docking_jobs"] == 0
    assert manifest["pocket_detector_runtime_executed"] is False
    assert manifest["compound_docking_or_scoring_executed"] is False
    assert manifest["candidate_nomination_executed"] is False
    assert manifest["source_file_sha256_entries"]["candidate_pockets"]["sha256"]


def test_no_external_scientific_runtime_tool_is_imported_or_executed(tmp_path, monkeypatch):
    root = copy_task9_inputs(tmp_path)
    blocked = {"pyrosetta", "lightdock", "vina", "fpocket", "p2rank", "prank"}
    original_import = builtins.__import__

    def guarded_import(name, globals=None, locals=None, fromlist=(), level=0):
        if name.split(".")[0].lower() in blocked:
            raise AssertionError("Task 9 must not import {}".format(name))
        return original_import(name, globals, locals, fromlist, level)

    monkeypatch.setattr(builtins, "__import__", guarded_import)
    ctx, report = run_prioritization(tmp_path, root)
    manifest = read_json(ctx.manifest_dir / "pocket_candidate_prioritization_manifest.json")

    assert report["pocket_detector_runtime_executed"] is False
    assert manifest["external_detector_outputs_ingested"] is True


def test_old_workflow_protected_paths_have_no_diff():
    proc = subprocess.run(
        ["git", "diff", "--name-only", "--", "run_production.py", "main.py", "egfr_pipeline/", "config/", "docs/runbook.md", "output/", "results_export/"],
        cwd=str(REPO_ROOT),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )

    assert proc.returncode == 0
    assert proc.stdout.strip() == ""

