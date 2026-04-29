"""Tests for M1 Phase 3 — MYO1D construct slicing + QC (handoff §16).

Closes M1 §23 #13: MYO1D 955-1006 construct QC works on synthetic fixture.
"""

import csv
import json
import os
import subprocess
import sys
import uuid
from pathlib import Path

import pytest

from egfr_myo1d.core.logging_utils import initialize_logs
from egfr_myo1d.core.run_context import RunContext
from egfr_myo1d.myo1d.construct import (
    CAP_RESNAMES,
    emit_myo1d_construct_pdb,
    slice_myo1d_construct,
)
from egfr_myo1d.myo1d.qc import (
    MYO1D_CONSTRUCT_QC_COLUMNS,
    Myo1dQcReport,
    expand_residue_set,
    parse_construct_range,
    run_myo1d_qc,
)
from egfr_myo1d.structure.pdb_parser import parse_pdb


REPO_ROOT = Path(__file__).resolve().parents[2]
M1_PHASE3_FIXTURES = REPO_ROOT / "fresh" / "tests" / "fixtures" / "m1_phase3_myo1d"

VALID_955_1006 = M1_PHASE3_FIXTURES / "myo1d_955_1006_valid.pdb"
TERMINAL_BAD_962_1006 = M1_PHASE3_FIXTURES / "myo1d_962_1006_terminal_bad.pdb"
WITH_CAPS = M1_PHASE3_FIXTURES / "myo1d_with_ace_nme_caps.pdb"
SHORT_955_1001 = M1_PHASE3_FIXTURES / "myo1d_955_1001_short.pdb"


def unique_run_id(prefix="pytest_m1_phase3"):
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


def read_csv_rows(path):
    with Path(path).open("r", encoding="utf-8", newline="") as h:
        return list(csv.DictReader(h))


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def assert_under_run_dir(ctx, path):
    try:
        Path(path).resolve().relative_to(ctx.run_dir.resolve())
    except ValueError:
        raise AssertionError("{0} escaped run dir".format(path))


# ---------------------------------------------------------------------------
# Helper unit tests
# ---------------------------------------------------------------------------


def test_parse_construct_range_string():
    assert parse_construct_range("955-1006") == (955, 1006)


def test_parse_construct_range_tuple():
    assert parse_construct_range((955, 1001)) == (955, 1001)


def test_parse_construct_range_invalid_end_lt_start():
    with pytest.raises(ValueError):
        parse_construct_range("1006-955")


def test_expand_residue_set_single():
    assert expand_residue_set("961-964") == [961, 962, 963, 964]


def test_expand_residue_set_multi():
    assert expand_residue_set("961-964,968-972") == [
        961, 962, 963, 964, 968, 969, 970, 971, 972
    ]


# ---------------------------------------------------------------------------
# Slice / emit primitives
# ---------------------------------------------------------------------------


def test_slice_myo1d_construct_preserves_numbering():
    structure = parse_pdb(VALID_955_1006)
    sliced = slice_myo1d_construct(structure, 955, 1006, include_caps=True)
    residues = sorted({r.residue_number for r in sliced.residues if r.resname not in CAP_RESNAMES})
    assert residues[0] >= 955 and residues[-1] <= 1006


def test_slice_myo1d_construct_excludes_outside_range():
    structure = parse_pdb(VALID_955_1006)
    sliced = slice_myo1d_construct(structure, 968, 972, include_caps=False)
    residues = sorted({r.residue_number for r in sliced.residues})
    assert residues == [968, 969, 970, 971, 972]


def test_emit_myo1d_construct_pdb_writes_under_run_dir(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    structure = parse_pdb(VALID_955_1006)
    sliced = slice_myo1d_construct(structure, 955, 1006)
    out = ctx.run_dir / "normalized" / "myo1d" / "MYO1D_955_1006.pdb"
    emit_myo1d_construct_pdb(ctx, sliced, out)
    assert out.is_file()
    assert_under_run_dir(ctx, out)


# ---------------------------------------------------------------------------
# run_myo1d_qc — primary cases
# ---------------------------------------------------------------------------


def test_myo1d_955_1006_construct_emitted_with_original_numbering(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)

    report = run_myo1d_qc(ctx, VALID_955_1006, "955-1006", profile="codex_dev")

    assert isinstance(report, Myo1dQcReport)
    assert report.construct_id == "MYO1D_955_1006"
    assert report.start_residue == 955
    assert report.end_residue == 1006
    out_pdb = ctx.run_dir / "normalized" / "myo1d" / "MYO1D_955_1006.pdb"
    assert out_pdb.is_file()
    # Re-parse output and confirm residues span 955-1006
    parsed = parse_pdb(out_pdb)
    biological = sorted({r.residue_number for r in parsed.residues if r.resname not in CAP_RESNAMES})
    assert biological[0] == 955
    assert biological[-1] == 1006


def test_myo1d_key_residues_present_assertions(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    report = run_myo1d_qc(ctx, VALID_955_1006, "955-1006", profile="codex_dev")
    assert report.key_sheet8_present is True
    assert report.key_sheet9_present is True
    assert report.key_sheet12_present is True


def test_myo1d_n_watch_recorded(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    report = run_myo1d_qc(ctx, VALID_955_1006, "955-1006", profile="codex_dev")
    assert report.n_watch_present is True


def test_myo1d_qc_csv_columns_match_spec(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    run_myo1d_qc(ctx, VALID_955_1006, "955-1006", profile="codex_dev")

    qc_csv = ctx.qc_dir / "myo1d_construct_qc.csv"
    assert qc_csv.is_file()
    rows = read_csv_rows(qc_csv)
    assert len(rows) == 1
    row = rows[0]
    for col in MYO1D_CONSTRUCT_QC_COLUMNS:
        assert col in row, "missing column: {0}".format(col)


def test_myo1d_manifest_includes_sha256(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    run_myo1d_qc(ctx, VALID_955_1006, "955-1006", profile="codex_dev")

    manifest = read_json(ctx.manifest_dir / "myo1d_construct_manifest.json")
    assert manifest["construct_id"] == "MYO1D_955_1006"
    assert manifest["source_sha256"] is not None
    assert manifest["output_pdb_sha256"] is not None
    assert manifest["score_bonus_allowed"] is False
    assert manifest["key_residue_bonus_weight"] == 0.0


def test_myo1d_phase_status_appended(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    run_myo1d_qc(ctx, VALID_955_1006, "955-1006", profile="codex_dev")

    phase_log = ctx.logs_dir / "phase_status.jsonl"
    lines = [json.loads(line) for line in phase_log.read_text(encoding="utf-8").splitlines() if line.strip()]
    cleanup_records = [r for r in lines if r.get("phase") == "prepare-myo1d"]
    assert len(cleanup_records) == 1
    rec = cleanup_records[0]
    assert rec["status"] in ("PASS", "WARN")
    assert rec["details"]["score_bonus_allowed"] is False


# ---------------------------------------------------------------------------
# run_myo1d_qc — variants and negative regression
# ---------------------------------------------------------------------------


def test_myo1d_955_1001_comparator_pass(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    report = run_myo1d_qc(ctx, SHORT_955_1001, "955-1001", profile="codex_dev")
    assert report.status in ("PASS", "WARN")
    assert report.start_residue == 955
    assert report.end_residue == 1001
    # Sheet 8/9/12 should still be present
    assert report.key_sheet8_present is True
    assert report.key_sheet9_present is True
    assert report.key_sheet12_present is True


def test_myo1d_962_start_terminal_artifact_warned_in_codex_dev(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    report = run_myo1d_qc(ctx, TERMINAL_BAD_962_1006, "962-1006", profile="codex_dev")
    assert report.status == "WARN"
    assert any("962_start_terminal_artifact" in w for w in report.warnings)


def test_myo1d_962_start_terminal_artifact_fails_in_hpc_strict(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    report = run_myo1d_qc(ctx, TERMINAL_BAD_962_1006, "962-1006", profile="hpc_strict")
    assert report.status == "FAIL"


def test_myo1d_ace_nme_cap_preserved_in_output_pdb(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    report = run_myo1d_qc(ctx, WITH_CAPS, "955-1001", profile="codex_dev")
    assert report.ace_nme_caps_present is True

    parsed = parse_pdb(report.output_pdb)
    cap_resnames = {r.resname for r in parsed.residues if r.resname in CAP_RESNAMES}
    assert "ACE" in cap_resnames or "NME" in cap_resnames


def test_myo1d_missing_source_warns_in_codex_dev(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    bogus = tmp_path / "does_not_exist.pdb"
    report = run_myo1d_qc(ctx, bogus, "955-1006", profile="codex_dev")
    assert report.status == "WARN"
    assert any("source_missing" in w for w in report.warnings)


def test_myo1d_missing_source_fails_in_hpc_strict(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    bogus = tmp_path / "does_not_exist.pdb"
    report = run_myo1d_qc(ctx, bogus, "955-1006", profile="hpc_strict")
    assert report.status == "FAIL"


def test_myo1d_invalid_construct_range_fails(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    with pytest.raises(ValueError):
        run_myo1d_qc(ctx, VALID_955_1006, "1006-955", profile="codex_dev")


def test_myo1d_writes_under_run_dir_only(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    report = run_myo1d_qc(ctx, VALID_955_1006, "955-1006", profile="codex_dev")
    assert_under_run_dir(ctx, report.output_pdb)
    assert_under_run_dir(ctx, report.output_qc_csv)
    assert_under_run_dir(ctx, report.output_manifest)


def test_myo1d_score_bonus_zero_enforced_via_gates_yaml(tmp_path):
    """run_myo1d_qc must read key_residue_bonus_weight from gates.yaml and
    refuse to proceed if it is non-zero (no hardcoding in module)."""
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    # Default gates.yaml has 0.0; report.key_residue_bonus_weight should equal it.
    report = run_myo1d_qc(ctx, VALID_955_1006, "955-1006", profile="codex_dev")
    assert report.key_residue_bonus_weight == 0.0
    manifest = read_json(report.output_manifest)
    assert manifest["key_residue_bonus_weight"] == 0.0
    assert manifest["score_bonus_allowed"] is False


def test_myo1d_default_construct_from_gates_yaml(tmp_path):
    """Passing construct_range=None resolves to gates.yaml myo1d.construct (955-1006)."""
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    report = run_myo1d_qc(ctx, VALID_955_1006, construct_range=None, profile="codex_dev")
    assert report.start_residue == 955
    assert report.end_residue == 1006


# ---------------------------------------------------------------------------
# CLI tests
# ---------------------------------------------------------------------------


def test_cli_help_includes_prepare_myo1d():
    env = os.environ.copy()
    env["PYTHONPATH"] = (
        str(REPO_ROOT / "fresh" / "src") + os.pathsep + env.get("PYTHONPATH", "")
    )
    proc = subprocess.run(
        [sys.executable, "-m", "egfr_myo1d.cli", "--help"],
        cwd=str(REPO_ROOT),
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    assert proc.returncode == 0
    out = (proc.stdout or b"").decode("utf-8", errors="replace") + (
        proc.stderr or b""
    ).decode("utf-8", errors="replace")
    assert "prepare-myo1d" in out


def test_cli_prepare_myo1d_subcommand_help():
    env = os.environ.copy()
    env["PYTHONPATH"] = (
        str(REPO_ROOT / "fresh" / "src") + os.pathsep + env.get("PYTHONPATH", "")
    )
    proc = subprocess.run(
        [sys.executable, "-m", "egfr_myo1d.cli", "prepare-myo1d", "--help"],
        cwd=str(REPO_ROOT),
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    assert proc.returncode == 0
    out = (proc.stdout or b"").decode("utf-8", errors="replace")
    assert "--run-id" in out
    assert "--source" in out
    assert "--construct" in out
    assert "--profile" in out


def test_cli_prepare_myo1d_path_traversal_rejected():
    env = os.environ.copy()
    env["PYTHONPATH"] = (
        str(REPO_ROOT / "fresh" / "src") + os.pathsep + env.get("PYTHONPATH", "")
    )
    proc = subprocess.run(
        [
            sys.executable,
            "-m",
            "egfr_myo1d.cli",
            "prepare-myo1d",
            "--run-id",
            "../bad_run",
            "--source",
            str(VALID_955_1006),
            "--construct",
            "955-1006",
        ],
        cwd=str(REPO_ROOT),
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    assert proc.returncode != 0
