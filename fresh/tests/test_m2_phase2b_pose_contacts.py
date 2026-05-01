import csv
import json
import os
import subprocess
import sys
import uuid
from pathlib import Path

from egfr_myo1d.core.logging_utils import initialize_logs
from egfr_myo1d.core.run_context import RunContext
from egfr_myo1d.io.residue_mapping import write_residue_mapping
from egfr_myo1d.ppi.collect_ppi_outputs import collect_ppi_outputs
from egfr_myo1d.ppi.pose_contacts import extract_m2_pose_contacts


REPO_ROOT = Path(__file__).resolve().parents[2]


def unique_run_id(prefix="pytest_m2_phase2b"):
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


def write_jsonl(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(row, sort_keys=True) + "\n")


def pdb_atom(serial, atom, resname, chain, resseq, x, y, z, element):
    return (
        "{record:<6}{serial:>5d} {atom:<4}{altloc:1}{resname:>3} {chain:1}"
        "{resseq:>4d}{icode:1}   {x:>8.3f}{y:>8.3f}{z:>8.3f}"
        "{occupancy:>6.2f}{bfactor:>6.2f}          {element:>2}\n"
    ).format(
        record="ATOM",
        serial=serial,
        atom=atom,
        altloc="",
        resname=resname,
        chain=chain,
        resseq=resseq,
        icode="",
        x=x,
        y=y,
        z=z,
        occupancy=1.0,
        bfactor=20.0,
        element=element,
    )


def make_pose_pdb(path):
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        pdb_atom(1, "CA", "GLY", "A", 669, 0.0, 0.0, 0.0, "C"),
        pdb_atom(2, "CB", "GLY", "A", 669, 0.0, 0.0, 1.0, "C"),
        pdb_atom(3, "CA", "LEU", "B", 1669, 10.0, 0.0, 0.0, "C"),
        pdb_atom(4, "CB", "LEU", "B", 1669, 10.0, 0.0, 1.0, "C"),
        pdb_atom(5, "CA", "SER", "C", 961, 0.0, 0.0, 4.0, "C"),
        pdb_atom(6, "CB", "SER", "C", 961, 0.0, 0.0, 3.5, "C"),
        pdb_atom(7, "CA", "THR", "C", 968, 10.0, 0.0, 4.0, "C"),
        pdb_atom(8, "CB", "THR", "C", 968, 10.0, 0.0, 3.5, "C"),
        "END\n",
    ]
    path.write_text("".join(lines), encoding="utf-8")
    return path


def setup_run_with_pose(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    ctx.create_directories()
    state_id = "EGFR_160-185"
    output_dir = ctx.run_dir / "phase1_ppi" / "pyrosetta_adapter" / "outputs" / "job_001"
    pose = make_pose_pdb(output_dir / "poses" / "pose_00001.pdb")
    mapping_csv = ctx.qc_dir / "{0}_receptor_mapping.csv".format(state_id)
    write_residue_mapping(
        mapping_csv,
        [
            {
                "state": state_id,
                "source_file": "synthetic",
                "protomer_id": "A",
                "source_chain": "A",
                "source_resseq": "669",
                "source_resname": "GLY",
                "runtime_chain": "A",
                "runtime_resseq": "669",
                "atom_count": "2",
                "role": "primary_membrane_validated_state",
            },
            {
                "state": state_id,
                "source_file": "synthetic",
                "protomer_id": "B",
                "source_chain": "B",
                "source_resseq": "669",
                "source_resname": "LEU",
                "runtime_chain": "B",
                "runtime_resseq": "1669",
                "atom_count": "2",
                "role": "primary_membrane_validated_state",
            },
        ],
        ctx,
    )
    job_manifest = ctx.run_dir / "phase1_ppi" / "pyrosetta_adapter" / "pyrosetta_job_manifest.jsonl"
    write_jsonl(
        job_manifest,
        [
            {
                "job_name": "job_001",
                "state_id": state_id,
                "state_role": "primary_membrane_validated_state",
                "method": "pyrosetta_global_ppi",
                "seed": 0,
                "models_per_seed": 1,
                "output_dir": str(output_dir),
                "mapping_csv": str(mapping_csv),
            }
        ],
    )
    return ctx, pose


def test_m2_2b_extracts_ab_to_c_contacts_for_mapping_restoration(tmp_path):
    ctx, pose = setup_run_with_pose(tmp_path)

    report = extract_m2_pose_contacts(ctx)

    assert report.status == "PASS"
    assert report.pose_count == 1
    assert report.raw_contact_count == 2
    rows = read_csv(report.raw_contact_table)
    assert {row["protomer_id"] for row in rows} == {"A", "B"}
    assert {row["egfr_chain_id"] for row in rows} == {"A", "B"}
    assert {row["myo1d_residue"] for row in rows} == {"961", "968"}
    assert all(row["pose_id"] == pose.stem for row in rows)
    assert all(row["accepted_pose_flag"] == "true" for row in rows)

    payload = read_json(report.manifest_path)
    assert payload["phase"] == "m2.2b-ppi-pose-contact-extraction"
    assert payload["raw_contact_count"] == 2
    assert "no_pyrosetta_import" in payload["non_goals_confirmed"]


def test_m2_2b_raw_contacts_feed_m2_3_mapping_restoration(tmp_path):
    ctx, _pose = setup_run_with_pose(tmp_path)
    extraction = extract_m2_pose_contacts(ctx)

    collection = collect_ppi_outputs(ctx, raw_contact_table=extraction.raw_contact_table)

    assert collection.status == "PASS"
    contacts = read_csv(collection.pose_contacts_csv)
    assert len(contacts) == 2
    by_protomer = {row["protomer_id"]: row for row in contacts}
    assert by_protomer["A"]["egfr_runtime_residue"] == "669"
    assert by_protomer["A"]["egfr_original_residue"] == "669"
    assert by_protomer["B"]["egfr_runtime_residue"] == "1669"
    assert by_protomer["B"]["egfr_original_residue"] == "669"
    assert by_protomer["B"]["egfr_chain_id"] == "B"
    assert by_protomer["B"]["accepted_pose_flag"] == "true"


def test_m2_2b_reads_scores_from_chunked_real_run_outputs(tmp_path):
    ctx, _pose = setup_run_with_pose(tmp_path)
    score_csv = (
        ctx.run_dir
        / "phase1_ppi"
        / "pyrosetta_adapter"
        / "outputs"
        / "job_001"
        / "chunks"
        / "chunk_001"
        / "pose_scores.csv"
    )
    score_csv.parent.mkdir(parents=True, exist_ok=True)
    score_csv.write_text(
        "run_id,job_name,state_id,seed,model_index,pose_id,score,pose_path,status\n"
        "{0},job_001,EGFR_160-185,0,1,pose_00001,-12.345,path,PASS\n".format(
            ctx.run_id
        ),
        encoding="utf-8",
    )

    report = extract_m2_pose_contacts(ctx)
    rows = read_csv(report.raw_contact_table)

    assert report.status == "PASS"
    assert rows
    assert {row["score"] for row in rows} == {"-12.345"}


def test_m2_2b_missing_pose_pdbs_warn_without_fabricating_contacts(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    ctx.create_directories()
    job_manifest = ctx.run_dir / "phase1_ppi" / "pyrosetta_adapter" / "pyrosetta_job_manifest.jsonl"
    write_jsonl(
        job_manifest,
        [
            {
                "job_name": "job_without_pose",
                "state_id": "EGFR_160-185",
                "method": "pyrosetta_global_ppi",
                "seed": 0,
                "output_dir": str(ctx.run_dir / "phase1_ppi" / "pyrosetta_adapter" / "outputs" / "job_without_pose"),
            }
        ],
    )

    report = extract_m2_pose_contacts(ctx)

    assert report.status == "PASS_WITH_WARNINGS"
    assert report.pose_count == 0
    assert report.raw_contact_count == 0
    assert "no_pose_pdbs_for_job:job_without_pose" in report.warnings
    assert read_csv(report.raw_contact_table) == []


def test_m2_2b_sources_do_not_execute_science_engines():
    source = (REPO_ROOT / "fresh/src/egfr_myo1d/ppi/pose_contacts.py").read_text(
        encoding="utf-8"
    )
    forbidden_snippets = [
        "import pyrosetta",
        "import lightdock",
        "import fpocket",
        "import p2rank",
        "subprocess",
        "os.system",
        "qsub ",
        "sbatch ",
        "vina ",
    ]
    for snippet in forbidden_snippets:
        assert snippet not in source


def test_cli_help_includes_extract_m2_ppi_contacts():
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
    out = proc.stdout.decode("utf-8", errors="replace") + proc.stderr.decode(
        "utf-8", errors="replace"
    )
    assert "extract-m2-ppi-contacts" in out
