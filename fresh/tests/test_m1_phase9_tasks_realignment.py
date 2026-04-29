"""Tests for M1 Phase 9 — Tasks 4-9 schema realignment.

Phase 9 takes an additive approach: it does NOT modify Task 4-9 module logic.
Instead it adds an alignment helper (`validation/m1_alignment.py`) that
records the cross-reference between M1 canonical outputs and the Task 4-9
modules that would consume them in a future M2 actual-execution phase.

The realignment tests verify:
  1. M1 outputs and Task 4-9 outputs coexist in a single run directory
  2. The m1_alignment.json manifest captures the consumption map
  3. The 79 prior Task 4-9 tests are unmodified (verified by full suite)
"""

import json
import os
import shutil
import subprocess
import sys
import uuid
from pathlib import Path

import pytest

from egfr_myo1d.core.logging_utils import initialize_logs
from egfr_myo1d.core.run_context import RunContext
from egfr_myo1d.orchestrator.prepare_inputs import run_prepare_inputs
from egfr_myo1d.validation.m1_alignment import (
    M1AlignmentReport,
    M1_TO_TASK_CONSUMPTION,
    RECEPTOR_TASK_CONSUMERS,
    record_m1_alignment,
)


REPO_ROOT = Path(__file__).resolve().parents[2]
INTEGRATION_FIXTURE = (
    REPO_ROOT / "fresh" / "tests" / "fixtures" / "m1_phase8_integration"
)
TASK4_FIXTURE = REPO_ROOT / "fresh" / "tests" / "fixtures" / "task4_inputs"
TASK7_FIXTURE = REPO_ROOT / "fresh" / "tests" / "fixtures" / "task7_ppi_consensus"


def unique_run_id(prefix="pytest_m1_phase9"):
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


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def assert_under_run_dir(ctx, path):
    Path(path).resolve().relative_to(ctx.run_dir.resolve())


# ---------------------------------------------------------------------------
# m1_alignment helper unit tests
# ---------------------------------------------------------------------------


def test_record_m1_alignment_returns_report_after_prepare_inputs(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    run_prepare_inputs(
        ctx,
        mode="smoke_input",
        profile="codex_dev",
        input_root=INTEGRATION_FIXTURE,
    )
    report = record_m1_alignment(ctx)
    assert isinstance(report, M1AlignmentReport)
    assert report.aligned_count > 0
    assert report.status in {"PASS", "WARN"}


def test_record_m1_alignment_records_membrane_frame(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    run_prepare_inputs(
        ctx,
        mode="smoke_input",
        profile="codex_dev",
        input_root=INTEGRATION_FIXTURE,
    )
    report = record_m1_alignment(ctx)
    assert report.m1_outputs_present.get("membrane_frame_json") is True
    assert (
        "manifest/membrane_frame.json"
        in report.m1_canonical_paths.get("membrane_frame_json", "")
    )


def test_record_m1_alignment_records_myo1d_normalized_pdb(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    run_prepare_inputs(
        ctx,
        mode="smoke_input",
        profile="codex_dev",
        input_root=INTEGRATION_FIXTURE,
    )
    report = record_m1_alignment(ctx)
    assert report.m1_outputs_present.get("myo1d_normalized_pdb") is True


def test_record_m1_alignment_records_per_state_receptor_artifacts(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    run_prepare_inputs(
        ctx,
        mode="smoke_input",
        profile="codex_dev",
        input_root=INTEGRATION_FIXTURE,
    )
    report = record_m1_alignment(ctx)
    # Look for at least one state's normalized receptor PDB
    receptor_keys = [
        k for k in report.m1_outputs_present if "normalized/receptors" in k
    ]
    assert receptor_keys, "no normalized receptor PDBs recorded"


def test_record_m1_alignment_consumption_map_includes_task_ids(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    run_prepare_inputs(
        ctx,
        mode="smoke_input",
        profile="codex_dev",
        input_root=INTEGRATION_FIXTURE,
    )
    report = record_m1_alignment(ctx)
    # All listed receptor task consumers should appear at least once
    for task_id in RECEPTOR_TASK_CONSUMERS:
        assert task_id in report.task_consumption_map, (
            "missing task_id in consumption map: {0}".format(task_id)
        )


def test_record_m1_alignment_writes_manifest_under_run_dir(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    run_prepare_inputs(
        ctx,
        mode="smoke_input",
        profile="codex_dev",
        input_root=INTEGRATION_FIXTURE,
    )
    record_m1_alignment(ctx)
    manifest = ctx.manifest_dir / "m1_alignment.json"
    assert manifest.is_file()
    assert_under_run_dir(ctx, manifest)
    payload = read_json(manifest)
    assert payload["run_id"] == ctx.run_id
    assert payload["score_bonus_allowed"] is False
    assert payload["approach"] == "additive_phase9_no_task_module_logic_changed"


def test_record_m1_alignment_warns_when_run_dir_empty_of_m1_artifacts(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    # Do NOT run prepare_inputs; run dir has no M1 outputs
    report = record_m1_alignment(ctx)
    assert report.status == "WARN"
    assert report.aligned_count == 0
    assert any("no_m1_canonical_artifacts_found" in n for n in report.notes)


def test_record_m1_alignment_appends_phase_status(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    run_prepare_inputs(
        ctx,
        mode="smoke_input",
        profile="codex_dev",
        input_root=INTEGRATION_FIXTURE,
    )
    record_m1_alignment(ctx)
    phase_log = ctx.logs_dir / "phase_status.jsonl"
    lines = [
        json.loads(l)
        for l in phase_log.read_text(encoding="utf-8").splitlines()
        if l.strip()
    ]
    records = [r for r in lines if r.get("phase") == "m1-alignment"]
    assert len(records) == 1
    assert records[0]["details"]["approach"] == "additive"
    assert records[0]["details"]["score_bonus_allowed"] is False


# ---------------------------------------------------------------------------
# End-to-end M1 + Tasks 4 / 6 / 7 / 8 / 9 coexistence
# ---------------------------------------------------------------------------


def test_m1_outputs_and_task4_outputs_coexist_in_same_run_dir(tmp_path):
    """Run M1 prepare-inputs first, then Task 4 prepare-ppi-inputs, both
    inside the same RunContext. Verify both layers' outputs are present."""
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)

    # M1 layer
    agg = run_prepare_inputs(
        ctx,
        mode="smoke_input",
        profile="codex_dev",
        input_root=INTEGRATION_FIXTURE,
    )
    assert agg.status in {"PASS", "PASS_WITH_WARNINGS"}

    # Task 4 layer (using existing Task 4 fixture)
    from egfr_myo1d.validation.prepared_inputs import prepare_ppi_inputs

    prepare_ppi_inputs(
        ctx, "smoke_env", "codex_dev", TASK4_FIXTURE, None, strict=False
    )

    # Both layers' outputs coexist
    # M1
    assert (ctx.manifest_dir / "membrane_frame.json").is_file()
    assert (ctx.run_dir / "normalized" / "myo1d" / "MYO1D_955_1006.pdb").is_file()
    # Task 4 (its existing output paths under prepared/)
    prepared_dir = ctx.run_dir / "prepared"
    assert prepared_dir.is_dir()
    assert (ctx.manifest_dir / "preparation_qc_report.json").is_file()


def test_m1_alignment_records_task_consumption_after_full_pipeline(tmp_path):
    """After M1 + Task 4 + Task 6 in one run dir, the alignment helper records
    the consumption map covering both layers."""
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)

    # M1
    run_prepare_inputs(
        ctx,
        mode="smoke_input",
        profile="codex_dev",
        input_root=INTEGRATION_FIXTURE,
    )

    # Task 4
    from egfr_myo1d.validation.prepared_inputs import prepare_ppi_inputs

    prepare_ppi_inputs(
        ctx, "smoke_env", "codex_dev", TASK4_FIXTURE, None, strict=False
    )

    # Task 6
    from egfr_myo1d.validation.ppi_sampling_plan import plan_ppi_sampling

    plan_ppi_sampling(
        ctx,
        "smoke_input",
        "codex_dev",
        TASK4_FIXTURE,
        None,
        strict=False,
        command_line="pytest m1_phase9_realignment",
    )

    # Alignment record
    report = record_m1_alignment(ctx)
    assert "task4_prepare_ppi_inputs" in report.task_consumption_map
    assert "task6_plan_ppi_sampling" in report.task_consumption_map
    # M1 alignment manifest exists alongside Task 4/6 outputs
    assert (ctx.manifest_dir / "m1_alignment.json").is_file()
    assert (ctx.manifest_dir / "preparation_qc_report.json").is_file()
    assert (ctx.manifest_dir / "ppi_sampling_plan_report.json").is_file()


def test_legacy_task4_outputs_still_emitted_unchanged(tmp_path):
    """Phase 9 explicitly preserves Task 4's existing `prepared/` emissions.
    This test asserts the additive approach: Task 4's output paths are still
    created exactly as before."""
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    from egfr_myo1d.validation.prepared_inputs import prepare_ppi_inputs

    prepare_ppi_inputs(
        ctx, "smoke_env", "codex_dev", TASK4_FIXTURE, None, strict=False
    )
    # The Task 4 output paths exist
    assert (ctx.manifest_dir / "preparation_qc_report.json").is_file()
    assert (ctx.manifest_dir / "prepared_input_manifest.json").is_file()
    prepared_root = ctx.run_dir / "prepared"
    assert prepared_root.is_dir()
    # Phase 9 made NO changes here — unchanged behaviour preserved
