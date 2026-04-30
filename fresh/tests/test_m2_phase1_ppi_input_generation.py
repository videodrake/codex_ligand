import json
import os
import subprocess
import sys
import uuid
from pathlib import Path

import pytest

from egfr_myo1d.core.logging_utils import initialize_logs
from egfr_myo1d.core.run_context import RunContext
from egfr_myo1d.io.residue_mapping import read_residue_mapping
from egfr_myo1d.m2.ppi_inputs import generate_m2_1_ppi_inputs
from egfr_myo1d.orchestrator.prepare_inputs import run_prepare_inputs
from egfr_myo1d.structure.pdb_parser import parse_pdb


REPO_ROOT = Path(__file__).resolve().parents[2]
INTEGRATION_FIXTURE = (
    REPO_ROOT / "fresh" / "tests" / "fixtures" / "m1_phase8_integration"
)


def unique_run_id(prefix="pytest_m2_phase1"):
    return "{0}_{1}".format(prefix, uuid.uuid4().hex[:12])


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


@pytest.fixture
def ctx_with_m1(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    aggregate = run_prepare_inputs(
        ctx,
        mode="smoke_input",
        profile="codex_dev",
        input_root=INTEGRATION_FIXTURE,
    )
    assert aggregate.status in {"PASS", "PASS_WITH_WARNINGS"}
    return ctx


def test_m2_1_generates_specs_from_m1_canonical_outputs(ctx_with_m1):
    report = generate_m2_1_ppi_inputs(ctx_with_m1)
    assert report.status in {"PASS", "PASS_WITH_WARNINGS"}
    assert report.manifest_path.is_file()
    assert report.qc_csv_path.is_file()
    assert report.summary_report_path.is_file()
    assert {pack.state_id for pack in report.packs} == {"EGFR_160-185", "EGFR_170-200"}

    payload = read_json(report.manifest_path)
    assert payload["execution_allowed"] is False
    assert payload["docking_executed"] is False
    assert payload["external_scientific_execution"] is False
    assert payload["ligand_or_compound_inputs_used"] is False
    assert payload["score_bonus_allowed"] is False

    first_pack = ctx_with_m1.run_dir / "prepared" / "m2_1_ppi_inputs" / "EGFR_160-185"
    pyrosetta_spec = read_json(
        first_pack / "specs" / "EGFR_160-185_pyrosetta-global-ppi_input_spec.json"
    )
    lightdock_spec = read_json(
        first_pack / "specs" / "EGFR_160-185_lightdock-gso-ppi_input_spec.json"
    )
    assert pyrosetta_spec["receptor"]["numbering_policy"] == "use_runtime_offset_receptor_and_mapping_csv"
    assert lightdock_spec["receptor"]["numbering_policy"] == "use_explicit_AB_dockable_receptor_and_mapping_csv"
    assert pyrosetta_spec["engine_command"] is None
    assert lightdock_spec["engine_command"] is None


def test_m2_1_preserves_chain_identity_and_runtime_mapping(ctx_with_m1):
    report = generate_m2_1_ppi_inputs(ctx_with_m1, states=["EGFR_160-185"])
    assert report.status in {"PASS", "PASS_WITH_WARNINGS"}
    pack_dir = ctx_with_m1.run_dir / "prepared" / "m2_1_ppi_inputs" / "EGFR_160-185"
    dockable_pack = pack_dir / "receptor" / "EGFR_160-185_dockable_669_1014_explicit_AB.pdb"
    assert parse_pdb(dockable_pack).chain_ids == ("A", "B")

    mapping = read_residue_mapping(ctx_with_m1.qc_dir / "EGFR_160-185_receptor_mapping.csv")
    b_rows = [row for row in mapping if row["protomer_id"] == "B"]
    assert b_rows
    assert {
        int(row["runtime_resseq"]) - int(row["source_resseq"])
        for row in b_rows
    } == {1000}


def test_m2_1_myo1d_primary_and_comparator_roles(ctx_with_m1):
    report = generate_m2_1_ppi_inputs(ctx_with_m1, states=["EGFR_160-185"])
    assert report.status in {"PASS", "PASS_WITH_WARNINGS"}
    pack_dir = ctx_with_m1.run_dir / "prepared" / "m2_1_ppi_inputs" / "EGFR_160-185"
    primary = parse_pdb(pack_dir / "myo1d" / "MYO1D_955_1001_primary.pdb")
    comparator = parse_pdb(pack_dir / "myo1d" / "MYO1D_955_1006_comparator.pdb")
    primary_residues = {residue.residue_number for residue in primary.residues}
    comparator_residues = {residue.residue_number for residue in comparator.residues}
    assert min(primary_residues) == 955
    assert max(primary_residues) <= 1001
    assert max(comparator_residues) == 1006

    spec = read_json(
        pack_dir / "specs" / "EGFR_160-185_pyrosetta-global-ppi_input_spec.json"
    )
    assert spec["partner"]["primary_role"] == "production_oriented_ppi_partner"
    assert spec["partner"]["comparator_role"] == "noise_monitor_only"
    assert spec["partner"]["terminal_artifact_962_start_allowed"] is False
    assert spec["partner"]["annotation"]["score_bonus_allowed"] is False


def test_m2_1_missing_m1_outputs_fails_without_external_execution(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    report = generate_m2_1_ppi_inputs(ctx, states=["EGFR_160-185"])
    assert report.status == "FAIL"
    assert report.manifest_path.is_file()
    payload = read_json(report.manifest_path)
    assert payload["docking_executed"] is False
    assert payload["external_scientific_execution"] is False
    assert payload["packs"][0]["blockers"]


def test_m2_1_outputs_stay_under_run_dir(ctx_with_m1):
    report = generate_m2_1_ppi_inputs(ctx_with_m1)
    paths = [report.manifest_path, report.qc_csv_path, report.summary_report_path]
    for pack in report.packs:
        for record in pack.output_files.values():
            raw = Path(record["path"])
            paths.append(raw if raw.is_absolute() else REPO_ROOT / raw)
        for spec_path in pack.method_specs.values():
            raw = Path(spec_path)
            paths.append(raw if raw.is_absolute() else REPO_ROOT / raw)

    for path in paths:
        path.resolve().relative_to(ctx_with_m1.run_dir.resolve())


def test_m2_1_module_does_not_import_or_execute_science_engines():
    source = (REPO_ROOT / "fresh" / "src" / "egfr_myo1d" / "m2" / "ppi_inputs.py").read_text(
        encoding="utf-8"
    )
    forbidden_snippets = [
        "import pyrosetta",
        "import lightdock",
        "import subprocess",
        "os.system",
        "Popen",
        ".run(",
        "qsub ",
    ]
    for snippet in forbidden_snippets:
        assert snippet not in source


def _cli_env():
    env = os.environ.copy()
    env["PYTHONPATH"] = (
        str(REPO_ROOT / "fresh" / "src") + os.pathsep + env.get("PYTHONPATH", "")
    )
    return env


def test_cli_help_includes_generate_m2_ppi_inputs():
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
    assert "generate-m2-ppi-inputs" in out

