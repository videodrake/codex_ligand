"""Tests for M1 Phase 5 — membrane frame generation (handoff §15).

Closes M1 §23 #12: state-aware membrane_frame.json schema is implemented and
populated from coordinates (no hardcoded fallback vectors).
"""

import csv
import json
import math
import os
import re
import subprocess
import sys
import uuid
from pathlib import Path

import pytest

from egfr_myo1d.core.logging_utils import initialize_logs
from egfr_myo1d.core.run_context import RunContext
from egfr_myo1d.model.membrane_frame import (
    ALL_STATES,
    COORDINATE_CONVENTION,
    MEMBRANE_FRAME_QC_COLUMNS,
    PRIMARY_STATES,
    REFERENCE_CONTROL_STATES,
    StateMembraneFrame,
    compute_membrane_frame,
    run_membrane_frame_computation,
    write_membrane_frame_qc_csv,
    write_state_aware_membrane_frame_json,
)


REPO_ROOT = Path(__file__).resolve().parents[2]
PHASE5_FIXTURES = REPO_ROOT / "fresh" / "tests" / "fixtures" / "m1_phase5_membrane_frame"

DIMER_TM_JM = PHASE5_FIXTURES / "synthetic_dimer_with_TM_JM.pdb"
DIMER_KINASE_ONLY = PHASE5_FIXTURES / "synthetic_dimer_kinase_only.pdb"
DIMER_3GT8 = PHASE5_FIXTURES / "synthetic_3gt8_raw_kinase.pdb"


def unique_run_id(prefix="pytest_m1_phase5"):
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


def read_csv_rows(path):
    with Path(path).open("r", encoding="utf-8", newline="") as h:
        return list(csv.DictReader(h))


# ---------------------------------------------------------------------------
# Per-state computation tests
# ---------------------------------------------------------------------------


def test_membrane_normal_computed_from_TM_JM_coords_when_available():
    frame = compute_membrane_frame(DIMER_TM_JM, None, "EGFR_160-185", profile="codex_dev")
    assert frame.frame_source == "state_full_frame"
    assert frame.n_membrane is not None
    # Fixture is built so principal axis points along +z (positive Z component)
    assert frame.n_membrane[2] > 0
    assert frame.n_membrane_norm == pytest.approx(1.0, abs=1e-6)
    # Z component should dominate (|z| > |x| + |y|) in the synthetic fixture
    nx, ny, nz = frame.n_membrane
    assert abs(nz) > abs(nx) + abs(ny)


def test_x_dimer_axis_computed_from_chain_A_to_B_centroids():
    frame = compute_membrane_frame(DIMER_TM_JM, None, "EGFR_160-185", profile="codex_dev")
    assert frame.x_dimer_axis is not None
    # X component dominant (chain B at x=15..18, chain A at x=0..-3)
    xx, xy, xz = frame.x_dimer_axis
    assert xx > 0.9  # almost pure +x
    assert frame.x_dimer_axis_norm == pytest.approx(1.0, abs=1e-6)
    assert frame.protomer_a_centroid is not None
    assert frame.protomer_b_centroid is not None
    # A on left, B on right
    assert frame.protomer_a_centroid[0] < frame.protomer_b_centroid[0]


def test_protomer_centroids_recorded():
    frame = compute_membrane_frame(DIMER_TM_JM, None, "EGFR_160-185", profile="codex_dev")
    assert frame.protomer_a_centroid is not None and len(frame.protomer_a_centroid) == 3
    assert frame.protomer_b_centroid is not None and len(frame.protomer_b_centroid) == 3
    assert frame.centroid_distance is not None and frame.centroid_distance > 0


def test_plus10_fallback_used_when_state_lacks_TM_JM():
    """state_full_frame_pdb=None, plus10_full_frame_pdb=DIMER_TM_JM -> plus10_inherited."""
    frame = compute_membrane_frame(None, DIMER_TM_JM, "EGFR_160-185", profile="codex_dev")
    assert frame.frame_source == "plus10_inherited"
    # Even when the fallback succeeds, the result is WARN-tagged
    assert frame.status == "WARN"
    assert any("plus10_inherited_frame_used_as_fallback" in w for w in frame.warnings)


def test_plus10_blank_chain_duplicate_ca_split_in_hpc_strict(tmp_path):
    blank_chain = tmp_path / "blank_chain_plus10.pdb"
    lines = []
    for line in DIMER_TM_JM.read_text(encoding="utf-8").splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            line = line[:21] + " " + line[22:]
        lines.append(line)
    blank_chain.write_text("\n".join(lines) + "\n", encoding="utf-8")

    frame = compute_membrane_frame(
        None,
        blank_chain,
        "EGFR_160-185",
        profile="hpc_strict",
    )

    assert frame.frame_source == "plus10_inherited"
    assert frame.status == "WARN"
    assert frame.n_membrane is not None
    assert frame.x_dimer_axis is not None
    assert "duplicate_chain_ca_split_into_A_B_for_membrane_frame" in frame.notes


def test_3gt8_raw_marked_reference_control_not_primary_source():
    frame = compute_membrane_frame(DIMER_3GT8, None, "3GT8_raw", profile="codex_dev")
    assert frame.role == "crystallographic_reference_control"
    assert frame.status == "reference_control"
    # No vectors emitted for reference_control
    assert frame.n_membrane is None
    assert frame.x_dimer_axis is None
    assert "not_primary" in frame.frame_source


def test_missing_frame_source_reported_cleanly_no_invented_vectors_codex_dev():
    bogus = REPO_ROOT / "fresh" / "tests" / "fixtures" / "_does_not_exist.pdb"
    frame = compute_membrane_frame(bogus, None, "EGFR_160-185", profile="codex_dev")
    assert frame.status == "WARN"
    assert frame.n_membrane is None
    assert frame.x_dimer_axis is None
    assert frame.frame_source == "missing"
    assert any("missing_frame_source" in w for w in frame.warnings)


def test_missing_frame_source_fails_in_hpc_strict():
    bogus = REPO_ROOT / "fresh" / "tests" / "fixtures" / "_does_not_exist.pdb"
    frame = compute_membrane_frame(bogus, None, "EGFR_160-185", profile="hpc_strict")
    assert frame.status == "FAIL"
    assert frame.n_membrane is None


def test_kinase_only_input_lacks_TM_JM_warning():
    frame = compute_membrane_frame(DIMER_KINASE_ONLY, None, "EGFR_160-185", profile="codex_dev")
    # No 634-674 atoms; n_membrane cannot be computed
    assert frame.n_membrane is None
    assert any("missing_TM_JM_residues_634_674" in w for w in frame.warnings)
    assert frame.n_tm_jm_residues == 0


# ---------------------------------------------------------------------------
# State-aware JSON + CSV writers
# ---------------------------------------------------------------------------


def test_membrane_frame_json_state_aware_schema(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    frames = [
        compute_membrane_frame(DIMER_TM_JM, None, "EGFR_160-185", profile="codex_dev"),
        compute_membrane_frame(DIMER_TM_JM, None, "EGFR_170-200", profile="codex_dev"),
        compute_membrane_frame(DIMER_3GT8, None, "3GT8_raw", profile="codex_dev"),
    ]
    json_path = write_state_aware_membrane_frame_json(ctx, frames)
    payload = read_json(json_path)

    assert payload["coordinate_convention"] == COORDINATE_CONVENTION
    assert "frame_source_policy" in payload
    assert set(payload["states"].keys()) == set(ALL_STATES)
    s_primary = payload["states"]["EGFR_160-185"]
    assert s_primary["role"] == "primary_membrane_validated_state"
    assert s_primary["n_membrane"] is not None
    assert len(s_primary["n_membrane"]) == 3
    s_ref = payload["states"]["3GT8_raw"]
    assert s_ref["role"] == "crystallographic_reference_control"
    assert s_ref["n_membrane"] is None


def test_membrane_frame_qc_csv_columns_match_spec(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    frames = [
        compute_membrane_frame(DIMER_TM_JM, None, "EGFR_160-185", profile="codex_dev"),
        compute_membrane_frame(DIMER_3GT8, None, "3GT8_raw", profile="codex_dev"),
    ]
    csv_path = write_membrane_frame_qc_csv(ctx, frames)
    rows = read_csv_rows(csv_path)
    assert len(rows) == 2
    for col in MEMBRANE_FRAME_QC_COLUMNS:
        assert col in rows[0]


# ---------------------------------------------------------------------------
# Top-level orchestrator
# ---------------------------------------------------------------------------


def test_run_membrane_frame_computation_all_states(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    frames, overall = run_membrane_frame_computation(
        ctx, full_frame_source=DIMER_TM_JM, profile="codex_dev"
    )
    assert {f.state_id for f in frames} == set(ALL_STATES)
    # 3GT8_raw is reference_control; primary states use the explicit override.
    assert overall == "PASS"
    assert (ctx.manifest_dir / "membrane_frame.json").is_file()
    assert (ctx.qc_dir / "membrane_frame_qc.csv").is_file()
    # Phase status appended
    phase_log = ctx.logs_dir / "phase_status.jsonl"
    lines = [
        json.loads(l) for l in phase_log.read_text(encoding="utf-8").splitlines() if l.strip()
    ]
    assert any(r.get("phase") == "compute-membrane-frame" for r in lines)


def test_full_frame_source_override_takes_priority_over_raw_state_pdb(tmp_path):
    repo_root = tmp_path / "repo"
    raw_dir = repo_root / "fresh" / "data" / "raw" / "receptors"
    config_dir = repo_root / "fresh" / "configs"
    raw_dir.mkdir(parents=True)
    config_dir.mkdir(parents=True)
    (raw_dir / "EGFR_160-185.pdb").write_text("not a usable pdb\n", encoding="utf-8")
    (config_dir / "receptor_states.yaml").write_text(
        "\n".join(
            [
                "input_files:",
                "  EGFR_160-185: fresh/data/raw/receptors/EGFR_160-185.pdb",
                "  EGFR_170-200: fresh/data/raw/receptors/EGFR_170-200.pdb",
                "  3GT8_raw: fresh/data/raw/receptors/3GT8_raw.pdb",
                "  plus10_full_frame: fresh/data/raw/receptors/plus10_full_frame.pdb",
                "",
            ]
        ),
        encoding="utf-8",
    )
    ctx = RunContext(
        repo_root=repo_root,
        fresh_root=repo_root / "fresh",
        run_id=unique_run_id(),
        run_dir=tmp_path / "run",
        manifest_dir=tmp_path / "run" / "manifest",
        logs_dir=tmp_path / "run" / "logs",
        jobs_log_dir=tmp_path / "run" / "logs" / "jobs",
        errors_dir=tmp_path / "run" / "logs" / "errors",
        qc_dir=tmp_path / "run" / "qc",
        reports_dir=tmp_path / "run" / "reports",
        scratch_dir=tmp_path / "run" / "scratch",
        tmp_dir=tmp_path / "run" / "tmp",
    )
    initialize_logs(ctx)

    frames, overall = run_membrane_frame_computation(
        ctx,
        state_ids=["EGFR_160-185"],
        full_frame_source=DIMER_TM_JM,
        profile="hpc_strict",
    )

    assert overall == "PASS"
    assert frames[0].frame_source == "state_full_frame"
    assert frames[0].status == "PASS"
    assert not any("parse_error" in warning for warning in frames[0].warnings)


def test_run_membrane_frame_computation_writes_under_run_dir_only(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    frames, _ = run_membrane_frame_computation(
        ctx, full_frame_source=DIMER_TM_JM, profile="codex_dev"
    )
    json_path = ctx.manifest_dir / "membrane_frame.json"
    csv_path = ctx.qc_dir / "membrane_frame_qc.csv"
    json_path.resolve().relative_to(ctx.run_dir.resolve())
    csv_path.resolve().relative_to(ctx.run_dir.resolve())


def test_run_membrane_frame_computation_unknown_state_raises(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    with pytest.raises(ValueError):
        run_membrane_frame_computation(
            ctx, state_ids=["bogus_state"], profile="codex_dev"
        )


def test_run_membrane_frame_computation_invalid_profile(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    with pytest.raises(ValueError):
        run_membrane_frame_computation(ctx, profile="rogue")


# ---------------------------------------------------------------------------
# Anti-hardcoded-vector check (handoff §15.1)
# ---------------------------------------------------------------------------


def test_no_hardcoded_001_or_100_vectors_in_module_source():
    """The module must not contain literal axis vectors (e.g., [0, 0, 1]) as
    fallback computation results. References inside docstrings or comments are
    allowed and stripped out before the search."""
    src_path = REPO_ROOT / "fresh" / "src" / "egfr_myo1d" / "model" / "membrane_frame.py"
    text = src_path.read_text(encoding="utf-8")
    # Strip triple-quoted strings (docstrings) so docs may mention literals
    cleaned = re.sub(r'"""[\s\S]*?"""', "", text)
    # Strip line comments
    cleaned = re.sub(r"#.*", "", cleaned)
    # Search for [0,0,1] / [1,0,0] tolerant of whitespace
    bad_z = re.search(r"\[\s*0\s*,\s*0\s*,\s*1\s*\]", cleaned)
    bad_x = re.search(r"\[\s*1\s*,\s*0\s*,\s*0\s*\]", cleaned)
    bad_neg_z = re.search(r"\[\s*0\s*,\s*0\s*,\s*-\s*1\s*\]", cleaned)
    assert bad_z is None, "literal [0, 0, 1] found in module source outside docstrings"
    assert bad_x is None, "literal [1, 0, 0] found in module source outside docstrings"
    assert bad_neg_z is None, "literal [0, 0, -1] found in module source outside docstrings"


# ---------------------------------------------------------------------------
# CLI tests
# ---------------------------------------------------------------------------


def _cli_env():
    env = os.environ.copy()
    env["PYTHONPATH"] = (
        str(REPO_ROOT / "fresh" / "src") + os.pathsep + env.get("PYTHONPATH", "")
    )
    return env


def test_cli_help_includes_compute_membrane_frame():
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
    assert "compute-membrane-frame" in out


def test_cli_compute_membrane_frame_subcommand_help():
    proc = subprocess.run(
        [sys.executable, "-m", "egfr_myo1d.cli", "compute-membrane-frame", "--help"],
        cwd=str(REPO_ROOT),
        env=_cli_env(),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    assert proc.returncode == 0
    out = proc.stdout.decode("utf-8", errors="replace")
    for flag in ("--run-id", "--state", "--full-frame-source", "--profile"):
        assert flag in out


def test_cli_compute_membrane_frame_path_traversal_rejected():
    proc = subprocess.run(
        [
            sys.executable,
            "-m",
            "egfr_myo1d.cli",
            "compute-membrane-frame",
            "--run-id",
            "../bad_run",
            "--state",
            "all",
            "--full-frame-source",
            str(DIMER_TM_JM),
        ],
        cwd=str(REPO_ROOT),
        env=_cli_env(),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    assert proc.returncode != 0
