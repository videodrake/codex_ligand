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
from egfr_myo1d.io.residue_mapping import read_residue_mapping
from egfr_myo1d.m2.ppi_inputs import generate_m2_1_ppi_inputs
from egfr_myo1d.orchestrator.prepare_inputs import run_prepare_inputs
from egfr_myo1d.ppi.collect_ppi_outputs import collect_ppi_outputs
from egfr_myo1d.ppi.pyrosetta_adapter import generate_pyrosetta_harness
from egfr_myo1d.ppi.run_ppi_job import dry_run_job


REPO_ROOT = Path(__file__).resolve().parents[2]
INTEGRATION_FIXTURE = (
    REPO_ROOT / "fresh" / "tests" / "fixtures" / "m1_phase8_integration"
)
RAW_CONTACT_FIELDS = [
    "run_id",
    "state_id",
    "state_role",
    "method",
    "seed",
    "job_name",
    "pose_id",
    "score",
    "protomer_id",
    "egfr_runtime_residue",
    "egfr_resname",
    "myo1d_residue",
    "myo1d_resname",
    "contact_type",
    "min_distance_angstrom",
    "egfr_atom",
    "myo1d_atom",
]


def unique_run_id(prefix="pytest_m2_phase3"):
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


def read_csv(path):
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_raw_contacts(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=RAW_CONTACT_FIELDS, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    return path


@pytest.fixture
def ctx_with_m2_2(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    aggregate = run_prepare_inputs(
        ctx,
        mode="smoke_input",
        profile="codex_dev",
        input_root=INTEGRATION_FIXTURE,
    )
    assert aggregate.status in {"PASS", "PASS_WITH_WARNINGS"}
    m2_1 = generate_m2_1_ppi_inputs(ctx)
    assert m2_1.status in {"PASS", "PASS_WITH_WARNINGS"}
    m2_2 = generate_pyrosetta_harness(ctx, mode="smoke_input")
    assert m2_2.status in {"PASS", "PASS_WITH_WARNINGS"}
    return ctx, m2_2


def _b_mapping_row(ctx, state_id="EGFR_160-185"):
    rows = read_residue_mapping(ctx.qc_dir / "{0}_receptor_mapping.csv".format(state_id))
    return next(row for row in rows if row["protomer_id"] == "B")


def test_m2_3_collects_dry_run_status_without_claiming_poses(ctx_with_m2_2):
    ctx, m2_2 = ctx_with_m2_2
    dry_run_job(ctx, m2_2.jobs[0].job_name)

    report = collect_ppi_outputs(ctx)

    assert report.status == "PASS_WITH_WARNINGS"
    assert report.raw_pose_table.is_file()
    assert report.pose_contacts_csv.is_file()
    raw_rows = read_csv(report.raw_pose_table)
    contact_rows = read_csv(report.pose_contacts_csv)
    assert raw_rows[0]["status"] == "DRY_RUN_ONLY"
    assert raw_rows[0]["pose_id"] == ""
    assert contact_rows == []

    payload = read_json(report.manifest_path)
    assert payload["docking_executed"] is False
    assert payload["pyrosetta_imported"] is False
    assert payload["contact_extraction_from_pose_pdb_executed"] is False
    assert payload["pose_acceptance_classification_executed"] is False


def test_m2_3_restores_runtime_b_chain_residue_to_original_numbering(ctx_with_m2_2):
    ctx, m2_2 = ctx_with_m2_2
    dry_run_job(ctx, m2_2.jobs[0].job_name)
    mapping = _b_mapping_row(ctx)
    raw_table = ctx.run_dir / "phase1_ppi" / "tables" / "synthetic_raw_contacts.csv"
    write_raw_contacts(
        raw_table,
        [
            {
                "run_id": ctx.run_id,
                "state_id": "EGFR_160-185",
                "state_role": "primary_membrane_validated_state",
                "method": "pyrosetta_global_ppi",
                "seed": "0",
                "job_name": m2_2.jobs[0].job_name,
                "pose_id": "synthetic_pose_001",
                "score": "-1.5",
                "protomer_id": "B",
                "egfr_runtime_residue": mapping["runtime_resseq"],
                "egfr_resname": mapping["source_resname"],
                "myo1d_residue": "961",
                "myo1d_resname": "GLY",
                "contact_type": "ca_contact",
                "min_distance_angstrom": "7.2",
                "egfr_atom": "CA",
                "myo1d_atom": "CA",
            }
        ],
    )

    report = collect_ppi_outputs(ctx, raw_contact_table=raw_table)

    assert report.status == "PASS"
    contacts = read_csv(report.pose_contacts_csv)
    assert len(contacts) == 1
    assert contacts[0]["protomer_id"] == "B"
    assert contacts[0]["egfr_runtime_residue"] == mapping["runtime_resseq"]
    assert contacts[0]["egfr_original_residue"] == mapping["source_resseq"]
    assert contacts[0]["egfr_uniprot_residue"] == mapping["source_resseq"]
    assert contacts[0]["state_role"] == "primary_membrane_validated_state"

    qc_rows = read_csv(report.mapping_qc_csv)
    assert qc_rows[0]["status"] == "PASS"
    for path in [
        report.raw_pose_table,
        report.pose_contacts_csv,
        report.mapping_qc_csv,
        report.unmapped_contacts_csv,
        report.summary_report_path,
        report.manifest_path,
    ]:
        path.resolve().relative_to(ctx.run_dir.resolve())


def test_m2_3_deduplicates_real_score_and_contact_pose_rows(ctx_with_m2_2):
    ctx, m2_2 = ctx_with_m2_2
    job = m2_2.jobs[0]
    output_dir = ctx.run_dir / "phase1_ppi" / "pyrosetta_adapter" / "outputs" / job.job_name
    score_csv = output_dir / "chunks" / "chunk_001" / "pose_scores.csv"
    score_csv.parent.mkdir(parents=True, exist_ok=True)
    pose_id = "{0}_model_00001".format(job.job_name)
    score_csv.write_text(
        "run_id,job_name,state_id,seed,model_index,pose_id,score,pose_path,status\n"
        "{0},{1},EGFR_160-185,0,1,{2},-10.0,path,PASS\n".format(
            ctx.run_id,
            job.job_name,
            pose_id,
        ),
        encoding="utf-8",
    )
    mapping = _b_mapping_row(ctx)
    raw_table = ctx.run_dir / "phase1_ppi" / "tables" / "synthetic_dedupe_contacts.csv"
    write_raw_contacts(
        raw_table,
        [
            {
                "run_id": ctx.run_id,
                "state_id": "EGFR_160-185",
                "state_role": "primary_membrane_validated_state",
                "method": "pyrosetta_global_ppi",
                "seed": "0",
                "job_name": job.job_name,
                "pose_id": pose_id,
                "score": "-10.0",
                "protomer_id": "B",
                "egfr_runtime_residue": mapping["runtime_resseq"],
                "egfr_resname": mapping["source_resname"],
                "myo1d_residue": "961",
                "myo1d_resname": "GLY",
                "contact_type": "ca_contact",
                "min_distance_angstrom": "7.2",
                "egfr_atom": "CA",
                "myo1d_atom": "CA",
            }
        ],
    )

    report = collect_ppi_outputs(ctx, raw_contact_table=raw_table)
    raw_rows = read_csv(report.raw_pose_table)

    assert report.status == "PASS"
    assert len([row for row in raw_rows if row["pose_id"] == pose_id]) == 1
    assert raw_rows[0]["notes"] == "real_run_pose_score"


def test_m2_3_unmapped_runtime_residue_fails_and_writes_unmapped_table(ctx_with_m2_2):
    ctx, m2_2 = ctx_with_m2_2
    raw_table = ctx.run_dir / "phase1_ppi" / "tables" / "synthetic_unmapped_contacts.csv"
    write_raw_contacts(
        raw_table,
        [
            {
                "run_id": ctx.run_id,
                "state_id": "EGFR_160-185",
                "state_role": "primary_membrane_validated_state",
                "method": "pyrosetta_global_ppi",
                "seed": "0",
                "job_name": m2_2.jobs[0].job_name,
                "pose_id": "synthetic_pose_unmapped",
                "score": "",
                "protomer_id": "B",
                "egfr_runtime_residue": "99999",
                "egfr_resname": "UNK",
                "myo1d_residue": "961",
                "myo1d_resname": "GLY",
                "contact_type": "ca_contact",
                "min_distance_angstrom": "8.0",
                "egfr_atom": "CA",
                "myo1d_atom": "CA",
            }
        ],
    )

    report = collect_ppi_outputs(ctx, raw_contact_table=raw_table)

    assert report.status == "FAIL"
    unmapped = read_csv(report.unmapped_contacts_csv)
    assert len(unmapped) == 1
    assert unmapped[0]["reason"] == "runtime_residue_not_in_m1_mapping"
    assert read_json(report.manifest_path)["docking_executed"] is False


def test_m2_3_sources_do_not_execute_science_engines():
    sources = [
        REPO_ROOT / "fresh/src/egfr_myo1d/ppi/collect_ppi_outputs.py",
        REPO_ROOT / "fresh/src/egfr_myo1d/ppi/restore_residue_mapping.py",
    ]
    forbidden_snippets = [
        "import pyrosetta",
        "import lightdock",
        "import fpocket",
        "import p2rank",
        "import subprocess",
        "os.system",
        "qsub ",
        "sbatch ",
        "vina ",
    ]
    for path in sources:
        source = path.read_text(encoding="utf-8")
        for snippet in forbidden_snippets:
            assert snippet not in source


def _cli_env():
    env = os.environ.copy()
    env["PYTHONPATH"] = (
        str(REPO_ROOT / "fresh" / "src") + os.pathsep + env.get("PYTHONPATH", "")
    )
    return env


def test_cli_help_includes_collect_m2_ppi_outputs():
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
    assert "collect-m2-ppi-outputs" in out
