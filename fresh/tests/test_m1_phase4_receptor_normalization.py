"""Tests for M1 Phase 4 — receptor normalization (handoff §14).

Closes M1 §23 #10 (explicit A/B normalization) and #11 (+1000 runtime offset).
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
from egfr_myo1d.io.residue_mapping import (
    MAPPING_CSV_COLUMNS,
    read_residue_mapping,
)
from egfr_myo1d.model.receptor_normalize import (
    NormalizedReceptor,
    PRIMARY_STATES,
    REFERENCE_CONTROL_STATES,
    normalize_receptor,
    resolve_role,
)
from egfr_myo1d.model.receptor_qc import (
    NON_RECEPTOR_HETATM,
    detect_normalization_case,
    detect_warn_mutations,
    split_duplicate_chain,
)
from egfr_myo1d.structure.pdb_parser import parse_pdb


REPO_ROOT = Path(__file__).resolve().parents[2]
M1_PHASE4_FIXTURES = REPO_ROOT / "fresh" / "tests" / "fixtures" / "m1_phase4_receptor"

EXPLICIT_AB = M1_PHASE4_FIXTURES / "explicit_AB_dimer.pdb"
DUPLICATE_X = M1_PHASE4_FIXTURES / "duplicate_chain_X_dimer.pdb"
MONOMER = M1_PHASE4_FIXTURES / "single_chain_monomer.pdb"
V924R = M1_PHASE4_FIXTURES / "v924r_warn.pdb"
DIMER_WITH_TM = M1_PHASE4_FIXTURES / "dimer_with_TM_excluded_range.pdb"


def unique_run_id(prefix="pytest_m1_phase4"):
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
    if path is None:
        return
    try:
        Path(path).resolve().relative_to(ctx.run_dir.resolve())
    except ValueError:
        raise AssertionError("{0} escaped run dir".format(path))


# ---------------------------------------------------------------------------
# Helper unit tests
# ---------------------------------------------------------------------------


def test_resolve_role_primary():
    assert resolve_role("EGFR_160-185") == "primary_membrane_validated_state"
    assert resolve_role("EGFR_170-200") == "primary_membrane_validated_state"


def test_resolve_role_3gt8_raw_reference_control():
    assert resolve_role("3GT8_raw") == "crystallographic_reference_control_not_primary_membrane_state"


def test_detect_normalization_case_explicit_AB():
    structure = parse_pdb(EXPLICIT_AB)
    assert detect_normalization_case(structure) == "A_explicit_AB"


def test_detect_normalization_case_duplicate_X():
    structure = parse_pdb(DUPLICATE_X)
    assert detect_normalization_case(structure) == "B_duplicate_chain"


def test_detect_normalization_case_monomer():
    structure = parse_pdb(MONOMER)
    assert detect_normalization_case(structure) == "C_monomer"


def test_split_duplicate_chain_separates_two_protomers():
    structure = parse_pdb(DUPLICATE_X)
    a, b = split_duplicate_chain(list(structure.atoms))
    assert len(a) == 5  # 5 atoms in first half (per fixture)
    assert len(b) == 5  # 5 in second half
    # Each half should contain the same residues 669, 670
    a_residues = sorted({atom.residue_number for atom in a})
    b_residues = sorted({atom.residue_number for atom in b})
    assert a_residues == [669, 670]
    assert b_residues == [669, 670]


def test_detect_warn_mutations_v924r():
    structure = parse_pdb(V924R)
    hits = detect_warn_mutations(structure.atoms)
    assert any(h["residue_number"] == 924 and h["observed_resname"] == "ARG" for h in hits)


# ---------------------------------------------------------------------------
# Case A — explicit AB
# ---------------------------------------------------------------------------


def test_explicit_AB_passthrough_preserves_chains_and_numbers(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    report = normalize_receptor(ctx, EXPLICIT_AB, "EGFR_160-185", profile="codex_dev")
    assert isinstance(report, NormalizedReceptor)
    assert report.case == "A_explicit_AB"
    assert report.protomer_count == 2
    assert report.status in ("PASS", "WARN")

    # Full frame contains both chains
    parsed = parse_pdb(report.full_frame_pdb)
    assert "A" in parsed.chain_ids and "B" in parsed.chain_ids
    # Residue numbers preserved (669, 670 in both chains)
    a_residues = parsed.residue_numbers_for_chain("A")
    b_residues = parsed.residue_numbers_for_chain("B")
    assert 669 in a_residues and 670 in a_residues
    assert 669 in b_residues and 670 in b_residues


# ---------------------------------------------------------------------------
# Case B — duplicate chain X
# ---------------------------------------------------------------------------


def test_duplicate_chain_X_detected_and_split_into_A_B(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    report = normalize_receptor(ctx, DUPLICATE_X, "EGFR_170-200", profile="codex_dev")
    assert report.case == "B_duplicate_chain"
    assert report.protomer_count == 2
    assert report.status == "WARN"
    assert any("case_B_duplicate_chain_split_into_A_B" in w for w in report.warnings)

    parsed = parse_pdb(report.full_frame_pdb)
    chains = set(parsed.chain_ids)
    assert "A" in chains and "B" in chains
    # Original X chain should not appear in normalized output
    assert "X" not in chains


# ---------------------------------------------------------------------------
# Case C — single-chain monomer
# ---------------------------------------------------------------------------


def test_monomer_only_input_warns_in_codex_dev(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    report = normalize_receptor(ctx, MONOMER, "EGFR_160-185", profile="codex_dev")
    assert report.case == "C_monomer"
    assert report.status == "WARN"
    assert any("not_valid_for_dimer_only_main_analysis" in w for w in report.warnings)


def test_monomer_only_input_fails_in_hpc_strict(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    report = normalize_receptor(ctx, MONOMER, "EGFR_160-185", profile="hpc_strict")
    assert report.case == "C_monomer"
    assert report.status == "FAIL"


# ---------------------------------------------------------------------------
# Dockable crop + lipid/water/ion drop
# ---------------------------------------------------------------------------


def test_dockable_crop_excludes_634_668_keeps_669_1014(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    report = normalize_receptor(ctx, DIMER_WITH_TM, "EGFR_160-185", profile="codex_dev")
    parsed = parse_pdb(report.dockable_pdb)
    residues = {atom.residue_number for atom in parsed.atoms}
    # TM-core residues 640, 660 must NOT be present
    assert 640 not in residues and 660 not in residues
    # JM/KD residues 670-1014 must be present
    assert {670, 700, 800, 1000}.issubset(residues)
    # Outside-crop 1100 must NOT be present
    assert 1100 not in residues


def test_lipid_hetatm_dropped_from_dockable(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    report = normalize_receptor(ctx, DIMER_WITH_TM, "EGFR_160-185", profile="codex_dev")
    parsed = parse_pdb(report.dockable_pdb)
    resnames = {atom.resname for atom in parsed.atoms}
    assert "CHL" not in resnames
    assert "HOH" not in resnames


# ---------------------------------------------------------------------------
# +1000 runtime offset
# ---------------------------------------------------------------------------


def test_runtime_offset_protomer_B_plus_1000_only(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    report = normalize_receptor(ctx, EXPLICIT_AB, "EGFR_160-185", profile="codex_dev")
    parsed = parse_pdb(report.runtime_offset_pdb)
    a_residues = parsed.residue_numbers_for_chain("A")
    b_residues = parsed.residue_numbers_for_chain("B")
    # Protomer A unchanged
    assert 669 in a_residues and 670 in a_residues
    # Protomer B residue numbers + 1000
    assert 1669 in b_residues and 1670 in b_residues
    # Original B residue numbers absent from offset receptor
    assert 669 not in b_residues and 670 not in b_residues


def test_runtime_offset_applied_to_case_B_split(tmp_path):
    """Case B: after split, protomer B (originally chain X) gets +1000 offset."""
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    report = normalize_receptor(ctx, DUPLICATE_X, "EGFR_170-200", profile="codex_dev")
    parsed = parse_pdb(report.runtime_offset_pdb)
    a_residues = parsed.residue_numbers_for_chain("A")
    b_residues = parsed.residue_numbers_for_chain("B")
    assert 669 in a_residues and 670 in a_residues
    assert 1669 in b_residues and 1670 in b_residues


# ---------------------------------------------------------------------------
# Mapping CSV
# ---------------------------------------------------------------------------


def test_mapping_csv_columns_match_spec(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    report = normalize_receptor(ctx, EXPLICIT_AB, "EGFR_160-185", profile="codex_dev")
    rows = read_residue_mapping(report.mapping_csv)
    assert rows
    for col in MAPPING_CSV_COLUMNS:
        assert col in rows[0]


def test_mapping_csv_round_trip_residue_identity(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    report = normalize_receptor(ctx, EXPLICIT_AB, "EGFR_160-185", profile="codex_dev")
    rows = read_residue_mapping(report.mapping_csv)
    a_rows = [r for r in rows if r["protomer_id"] == "A"]
    b_rows = [r for r in rows if r["protomer_id"] == "B"]
    assert a_rows and b_rows
    # Protomer A: source_resseq == runtime_resseq
    for r in a_rows:
        assert int(r["source_resseq"]) == int(r["runtime_resseq"])
    # Protomer B: runtime_resseq == source_resseq + 1000
    for r in b_rows:
        assert int(r["runtime_resseq"]) == int(r["source_resseq"]) + 1000


# ---------------------------------------------------------------------------
# 3GT8_raw role
# ---------------------------------------------------------------------------


def test_3gt8_raw_marked_reference_control_not_primary(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    report = normalize_receptor(ctx, EXPLICIT_AB, "3GT8_raw", profile="codex_dev")
    assert report.role == "crystallographic_reference_control_not_primary_membrane_state"
    # 3GT8_raw should still emit explicit_AB and runtime_offset, but not the
    # dockable_669_1014 crop file (per handoff §14.2: kinase-only crystal)
    assert report.full_frame_pdb is not None
    assert report.runtime_offset_pdb is not None
    assert report.dockable_pdb is None


# ---------------------------------------------------------------------------
# V924R warning
# ---------------------------------------------------------------------------


def test_v924r_residue_warned_not_mutated(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    report = normalize_receptor(ctx, V924R, "EGFR_160-185", profile="codex_dev")
    assert report.v924r_warn is True
    assert report.v924r_handled == "warn_only_not_mutated"
    assert any(h["residue_number"] == 924 and h["observed_resname"] == "ARG" for h in report.warn_mutation_hits)
    # Confirm output PDB still contains ARG924 (not mutated)
    if report.full_frame_pdb:
        parsed = parse_pdb(report.full_frame_pdb)
        residue_924_resnames = {a.resname for a in parsed.atoms if a.residue_number == 924}
        assert "ARG" in residue_924_resnames
        assert "VAL" not in residue_924_resnames


# ---------------------------------------------------------------------------
# Output containment + manifest
# ---------------------------------------------------------------------------


def test_normalization_writes_under_run_dir_only(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    report = normalize_receptor(ctx, EXPLICIT_AB, "EGFR_160-185", profile="codex_dev")
    assert_under_run_dir(ctx, report.full_frame_pdb)
    assert_under_run_dir(ctx, report.dockable_pdb)
    assert_under_run_dir(ctx, report.runtime_offset_pdb)
    assert_under_run_dir(ctx, report.mapping_csv)
    assert_under_run_dir(ctx, report.audit_csv)
    assert_under_run_dir(ctx, report.manifest_json)


def test_normalization_audit_csv_records_records(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    report = normalize_receptor(ctx, EXPLICIT_AB, "EGFR_160-185", profile="codex_dev")
    rows = read_csv_rows(report.audit_csv)
    assert rows
    assert any(row["protomer_id"] == "A" for row in rows)
    assert any(row["protomer_id"] == "B" for row in rows)


def test_manifest_includes_all_outputs_and_sha256(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    report = normalize_receptor(ctx, EXPLICIT_AB, "EGFR_160-185", profile="codex_dev")
    manifest = read_json(report.manifest_json)
    assert manifest["state_id"] == "EGFR_160-185"
    assert manifest["role"] == "primary_membrane_validated_state"
    assert manifest["case"] == "A_explicit_AB"
    assert manifest["source_sha256"] is not None
    outputs = manifest["outputs"]
    assert outputs["full_frame_sha256"] is not None
    assert outputs["dockable_sha256"] is not None
    assert outputs["runtime_offset_sha256"] is not None
    assert outputs["mapping_csv"] is not None
    assert outputs["audit_csv"] is not None
    assert manifest["score_bonus_allowed"] is False


def test_phase_status_appended(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    normalize_receptor(ctx, EXPLICIT_AB, "EGFR_160-185", profile="codex_dev")
    phase_log = ctx.logs_dir / "phase_status.jsonl"
    lines = [
        json.loads(line)
        for line in phase_log.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    records = [r for r in lines if r.get("phase") == "prepare-receptor"]
    assert len(records) == 1
    assert records[0]["details"]["score_bonus_allowed"] is False


# ---------------------------------------------------------------------------
# CLI tests
# ---------------------------------------------------------------------------


def _cli_env():
    env = os.environ.copy()
    env["PYTHONPATH"] = (
        str(REPO_ROOT / "fresh" / "src") + os.pathsep + env.get("PYTHONPATH", "")
    )
    return env


def test_cli_help_includes_prepare_receptor():
    proc = subprocess.run(
        [sys.executable, "-m", "egfr_myo1d.cli", "--help"],
        cwd=str(REPO_ROOT),
        env=_cli_env(),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    assert proc.returncode == 0
    out = proc.stdout.decode("utf-8", errors="replace") + proc.stderr.decode("utf-8", errors="replace")
    assert "prepare-receptor" in out


def test_cli_prepare_receptor_subcommand_help():
    proc = subprocess.run(
        [sys.executable, "-m", "egfr_myo1d.cli", "prepare-receptor", "--help"],
        cwd=str(REPO_ROOT),
        env=_cli_env(),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    assert proc.returncode == 0
    out = proc.stdout.decode("utf-8", errors="replace")
    for flag in ("--run-id", "--state", "--source", "--profile", "--strict"):
        assert flag in out


def test_cli_prepare_receptor_path_traversal_rejected():
    proc = subprocess.run(
        [
            sys.executable,
            "-m",
            "egfr_myo1d.cli",
            "prepare-receptor",
            "--run-id",
            "../bad_run",
            "--state",
            "EGFR_160-185",
            "--source",
            str(EXPLICIT_AB),
        ],
        cwd=str(REPO_ROOT),
        env=_cli_env(),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    assert proc.returncode != 0
