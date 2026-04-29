"""Tests for M1 Phase 7 — ligand manifest shell (handoff §17).

Closes M1 §23 #14: ligand manifest shell exists and never leaks internal IDs.
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
from egfr_myo1d.ligand.manifest import (
    LIGAND_MANIFEST_QC_COLUMNS,
    LigandManifest,
    PRIVATE_MAPPING_REQUIRED_COLUMNS,
    SUPPORTED_FORMATS,
    build_ligand_manifest,
    load_default_paths,
    load_public_ids,
)


REPO_ROOT = Path(__file__).resolve().parents[2]
PHASE7_FIXTURES = REPO_ROOT / "fresh" / "tests" / "fixtures" / "m1_phase7_ligand"
LIGANDS_DIR = PHASE7_FIXTURES
PRIVATE_MAPPING = PHASE7_FIXTURES / "compound_id_map.csv"

INTERNAL_TEST_IDS = (
    "INTERNAL_TEST_PLACEHOLDER_A",
    "INTERNAL_TEST_PLACEHOLDER_B",
    "INTERNAL_TEST_PLACEHOLDER_C",
)


def unique_run_id(prefix="pytest_m1_phase7"):
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
    Path(path).resolve().relative_to(ctx.run_dir.resolve())


def assert_no_internal_id_leak(*paths):
    for path in paths:
        if path is None or not Path(path).is_file():
            continue
        text = Path(path).read_text(encoding="utf-8")
        for sentinel in INTERNAL_TEST_IDS:
            assert sentinel not in text, "internal ID '{0}' leaked into {1}".format(
                sentinel, path
            )


# ---------------------------------------------------------------------------
# Config loaders
# ---------------------------------------------------------------------------


def test_load_public_ids_from_fresh_run_yaml(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    ids = load_public_ids(ctx)
    assert "Cpd-A" in ids
    assert "Cpd-B" in ids
    assert "Cpd-C" in ids


def test_load_default_paths_resolves_to_fresh_data(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    raw_ligands, private_mapping = load_default_paths(ctx)
    assert raw_ligands.parts[-2:] == ("raw", "ligands")
    assert private_mapping.name == "compound_id_map.csv"


# ---------------------------------------------------------------------------
# Primary cases
# ---------------------------------------------------------------------------


def test_ligand_manifest_uses_only_public_ids_in_qc_csv(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    manifest = build_ligand_manifest(
        ctx,
        ligands_dir=LIGANDS_DIR,
        private_mapping_path=PRIVATE_MAPPING,
        profile="codex_dev",
    )
    rows = read_csv_rows(manifest.output_qc_csv)
    public_ids = {row["public_id"] for row in rows}
    assert public_ids == {"Cpd-A", "Cpd-B", "Cpd-C"}
    # Required column set matches spec
    for col in LIGAND_MANIFEST_QC_COLUMNS:
        assert col in rows[0]


def test_ligand_internal_ids_never_appear_in_public_outputs(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    manifest = build_ligand_manifest(
        ctx,
        ligands_dir=LIGANDS_DIR,
        private_mapping_path=PRIVATE_MAPPING,
        profile="codex_dev",
    )
    assert manifest.internal_ids_leaked_into_outputs is False
    assert manifest.leaked_internal_ids == []
    assert_no_internal_id_leak(manifest.output_qc_csv, manifest.output_manifest_json)


def test_ligand_files_present_records_sha256(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    manifest = build_ligand_manifest(
        ctx,
        ligands_dir=LIGANDS_DIR,
        private_mapping_path=PRIVATE_MAPPING,
        profile="codex_dev",
    )
    for record in manifest.ligand_files:
        assert record.exists is True
        assert record.format == "sdf"
        assert record.sha256 is not None
        assert record.size_bytes is not None and record.size_bytes > 0
    assert manifest.present_count == 3
    assert manifest.missing_count == 0


def test_ligand_manifest_status_pass_when_complete(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    manifest = build_ligand_manifest(
        ctx,
        ligands_dir=LIGANDS_DIR,
        private_mapping_path=PRIVATE_MAPPING,
        profile="codex_dev",
    )
    assert manifest.status == "PASS"


# ---------------------------------------------------------------------------
# Severity matrix (codex_dev/hpc_strict × stage true/false)
# ---------------------------------------------------------------------------


def test_ligand_files_missing_codex_dev_no_stage_warns(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    bogus_dir = tmp_path / "no_ligands"
    bogus_dir.mkdir()
    manifest = build_ligand_manifest(
        ctx,
        ligands_dir=bogus_dir,
        private_mapping_path=PRIVATE_MAPPING,
        profile="codex_dev",
        compound_stage_enabled=False,
    )
    assert manifest.status == "WARN"
    assert manifest.missing_count == 3


def test_ligand_files_missing_codex_dev_with_stage_warns(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    bogus_dir = tmp_path / "no_ligands"
    bogus_dir.mkdir()
    manifest = build_ligand_manifest(
        ctx,
        ligands_dir=bogus_dir,
        private_mapping_path=PRIVATE_MAPPING,
        profile="codex_dev",
        compound_stage_enabled=True,
    )
    assert manifest.status == "WARN"


def test_ligand_files_missing_hpc_strict_no_stage_warns(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    bogus_dir = tmp_path / "no_ligands"
    bogus_dir.mkdir()
    manifest = build_ligand_manifest(
        ctx,
        ligands_dir=bogus_dir,
        private_mapping_path=PRIVATE_MAPPING,
        profile="hpc_strict",
        compound_stage_enabled=False,
    )
    assert manifest.status == "WARN"


def test_ligand_files_missing_hpc_strict_with_stage_fails(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    bogus_dir = tmp_path / "no_ligands"
    bogus_dir.mkdir()
    manifest = build_ligand_manifest(
        ctx,
        ligands_dir=bogus_dir,
        private_mapping_path=PRIVATE_MAPPING,
        profile="hpc_strict",
        compound_stage_enabled=True,
    )
    assert manifest.status == "FAIL"


# ---------------------------------------------------------------------------
# Private mapping
# ---------------------------------------------------------------------------


def test_private_mapping_schema_validated_when_present(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    manifest = build_ligand_manifest(
        ctx,
        ligands_dir=LIGANDS_DIR,
        private_mapping_path=PRIVATE_MAPPING,
        profile="codex_dev",
    )
    assert manifest.private_mapping_present is True
    assert manifest.private_mapping_schema_valid is True


def test_private_mapping_absent_warns(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    manifest = build_ligand_manifest(
        ctx,
        ligands_dir=LIGANDS_DIR,
        private_mapping_path=tmp_path / "nonexistent_mapping.csv",
        profile="codex_dev",
    )
    assert manifest.private_mapping_present is False
    assert manifest.status == "WARN"
    assert any("private_mapping_not_present_or_unreadable" in w for w in manifest.warnings)


def test_private_mapping_invalid_schema_warns(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    bad_mapping = tmp_path / "bad_mapping.csv"
    bad_mapping.write_text("wrong_col,another\nfoo,bar\n", encoding="utf-8")
    manifest = build_ligand_manifest(
        ctx,
        ligands_dir=LIGANDS_DIR,
        private_mapping_path=bad_mapping,
        profile="codex_dev",
    )
    assert manifest.private_mapping_schema_valid is False
    assert any("private_mapping_schema_invalid" in w for w in manifest.warnings)


# ---------------------------------------------------------------------------
# Internal-ID leak detection (FAIL path)
# ---------------------------------------------------------------------------


def test_internal_id_leak_detection_fails_status(tmp_path):
    """If something puts an internal_id into the QC outputs, status FAILs.

    We force a leak by giving a ligands_dir that contains a file whose name
    matches a public_id but whose content contains the internal_id sentinel.
    Since the file's sha256 ends up in outputs (not its content), the leak
    needs another vector — we test by having the public_id match an internal
    sentinel directly (constructed test).
    """
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    # Make a tampered private mapping where public_id == sentinel
    leaked_dir = tmp_path / "leaked_ligands"
    leaked_dir.mkdir()
    (leaked_dir / "INTERNAL_TEST_LEAK.sdf").write_text(
        (PHASE7_FIXTURES / "Cpd-A.sdf").read_text(encoding="utf-8"),
        encoding="utf-8",
    )
    leak_mapping = tmp_path / "leak_mapping.csv"
    leak_mapping.write_text(
        "public_id,internal_id,notes\n"
        "INTERNAL_TEST_LEAK,INTERNAL_TEST_LEAK,placeholder\n",
        encoding="utf-8",
    )
    manifest = build_ligand_manifest(
        ctx,
        public_ids=["INTERNAL_TEST_LEAK"],
        ligands_dir=leaked_dir,
        private_mapping_path=leak_mapping,
        profile="codex_dev",
    )
    assert manifest.internal_ids_leaked_into_outputs is True
    assert manifest.status == "FAIL"
    assert manifest.leaked_internal_id_count == 1
    assert "INTERNAL_TEST_LEAK" not in "\n".join(manifest.leaked_internal_ids)
    assert manifest.leaked_internal_ids[0].startswith("<redacted_internal_id_sha256:")
    assert all("INTERNAL_TEST_LEAK" not in warning for warning in manifest.warnings)
    payload = read_json(manifest.output_manifest_json)
    assert payload["leaked_internal_id_count"] == 1
    assert payload["leaked_internal_ids"][0].startswith(
        "<redacted_internal_id_sha256:"
    )
    assert "INTERNAL_TEST_LEAK" not in "\n".join(payload["leaked_internal_ids"])
    assert all("INTERNAL_TEST_LEAK" not in warning for warning in payload["warnings"])


# ---------------------------------------------------------------------------
# Output containment + invariants
# ---------------------------------------------------------------------------


def test_ligand_manifest_writes_under_run_dir_only(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    manifest = build_ligand_manifest(
        ctx,
        ligands_dir=LIGANDS_DIR,
        private_mapping_path=PRIVATE_MAPPING,
        profile="codex_dev",
    )
    assert_under_run_dir(ctx, manifest.output_qc_csv)
    assert_under_run_dir(ctx, manifest.output_manifest_json)


def test_phase_status_appended(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    build_ligand_manifest(
        ctx,
        ligands_dir=LIGANDS_DIR,
        private_mapping_path=PRIVATE_MAPPING,
        profile="codex_dev",
    )
    phase_log = ctx.logs_dir / "phase_status.jsonl"
    lines = [
        json.loads(l)
        for l in phase_log.read_text(encoding="utf-8").splitlines()
        if l.strip()
    ]
    records = [r for r in lines if r.get("phase") == "manifest-ligands"]
    assert len(records) == 1
    assert records[0]["details"]["score_bonus_allowed"] is False


def test_invalid_profile_raises(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    with pytest.raises(ValueError):
        build_ligand_manifest(
            ctx,
            ligands_dir=LIGANDS_DIR,
            private_mapping_path=PRIVATE_MAPPING,
            profile="rogue",
        )


# ---------------------------------------------------------------------------
# .gitignore protection of private mapping
# ---------------------------------------------------------------------------


def test_gitignore_protects_private_mapping():
    gitignore = (REPO_ROOT / ".gitignore").read_text(encoding="utf-8")
    # Either the directory or the specific file pattern must be ignored
    assert any(
        pattern in gitignore
        for pattern in (
            "fresh/data/private/*",
            "fresh/data/private/compound_id_map.csv",
            "fresh/data/private",
        )
    ), ".gitignore must protect fresh/data/private (compound_id_map.csv)"


# ---------------------------------------------------------------------------
# CLI tests
# ---------------------------------------------------------------------------


def _cli_env():
    env = os.environ.copy()
    env["PYTHONPATH"] = (
        str(REPO_ROOT / "fresh" / "src") + os.pathsep + env.get("PYTHONPATH", "")
    )
    return env


def test_cli_help_includes_manifest_ligands():
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
    assert "manifest-ligands" in out


def test_cli_manifest_ligands_subcommand_help():
    proc = subprocess.run(
        [sys.executable, "-m", "egfr_myo1d.cli", "manifest-ligands", "--help"],
        cwd=str(REPO_ROOT),
        env=_cli_env(),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    assert proc.returncode == 0
    out = proc.stdout.decode("utf-8", errors="replace")
    for flag in ("--run-id", "--ligands-dir", "--private-mapping", "--profile"):
        assert flag in out


def test_cli_manifest_ligands_path_traversal_rejected():
    proc = subprocess.run(
        [
            sys.executable,
            "-m",
            "egfr_myo1d.cli",
            "manifest-ligands",
            "--run-id",
            "../bad_run",
            "--ligands-dir",
            str(LIGANDS_DIR),
        ],
        cwd=str(REPO_ROOT),
        env=_cli_env(),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    assert proc.returncode != 0
