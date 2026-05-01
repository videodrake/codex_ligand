import csv
import json
import subprocess
import sys
import uuid
from pathlib import Path

from egfr_myo1d.core.logging_utils import initialize_logs
from egfr_myo1d.core.run_context import RunContext
from egfr_myo1d.ppi.pyrosetta_real_jobs import plan_pyrosetta_real_pbs


REPO_ROOT = Path(__file__).resolve().parents[2]


def unique_run_id(prefix="pytest_m2_phase2c"):
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


def read_csv(path):
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_jsonl(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(row, sort_keys=True) + "\n")


def setup_m2_2_job_manifest(ctx):
    manifest = ctx.run_dir / "phase1_ppi" / "pyrosetta_adapter" / "pyrosetta_job_manifest.jsonl"
    rows = []
    for state in ["EGFR_160-185", "EGFR_170-200"]:
        for seed in [0]:
            job_name = "m2_2_pyrosetta_{0}_seed{1:03d}".format(state, seed)
            rows.append(
                {
                    "job_name": job_name,
                    "state_id": state,
                    "method": "pyrosetta_global_ppi",
                    "seed": seed,
                    "models_per_seed": 100,
                    "output_dir": str(
                        ctx.run_dir
                        / "phase1_ppi"
                        / "pyrosetta_adapter"
                        / "outputs"
                        / job_name
                    ),
                    "receptor_pdb": str(ctx.run_dir / "prepared" / state / "receptor.pdb"),
                    "partner_pdb": str(ctx.run_dir / "prepared" / state / "myo1d.pdb"),
                    "mapping_csv": str(ctx.qc_dir / "{0}_receptor_mapping.csv".format(state)),
                }
            )
    write_jsonl(manifest, rows)
    return manifest


def test_m2_2c_plans_three_node_ppn32_chunked_pyrosetta_pbs(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    ctx.create_directories()
    manifest = setup_m2_2_job_manifest(ctx)

    plan = plan_pyrosetta_real_pbs(
        ctx,
        profile="production",
        job_manifest=manifest,
        nodes=["node04", "node05", "node06"],
        ppn=32,
        models_per_chunk=10,
    )

    assert plan.status == "PASS"
    assert plan.counts["planned_models"] == 200
    assert plan.counts["planned_chunks"] == 20
    assert len(plan.pbs_rows) == 3
    assert {row["node"] for row in plan.pbs_rows} == {"node04", "node05", "node06"}
    assert {row["ppn"] for row in plan.pbs_rows} == {"32"}
    assert all(row["qsub_command"].startswith("qsub ") for row in plan.pbs_rows)
    assert all(row["qsub_command"].endswith(".pbs") for row in plan.pbs_rows)

    chunks = read_csv(plan.chunk_manifest_csv)
    assert len(chunks) == 20
    assert {row["model_start"] for row in chunks if row["seed"] == "0"} >= {"1", "11", "91"}
    assert all("--model-start" in row["command"] for row in chunks)
    assert all("--chunk-id" in row["command"] for row in chunks)

    pbs_text = "\n".join(
        (ctx.repo_root / row["pbs_file"]).read_text(encoding="utf-8")
        for row in plan.pbs_rows
    )
    assert "#PBS -l nodes=node04:ppn=32" in pbs_text
    assert "#PBS -l nodes=node05:ppn=32" in pbs_text
    assert "#PBS -l nodes=node06:ppn=32" in pbs_text
    assert 'MAX_PARALLEL="${PBS_NP:-32}"' in pbs_text
    assert "wait_for_slot" in pbs_text
    assert "--dry-run false" in pbs_text
    assert "--allow-real true" in pbs_text
    assert "qsub " not in pbs_text


def test_m2_2c_auto_chunk_size_uses_one_worker_wave_per_seed(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    ctx.create_directories()
    manifest = setup_m2_2_job_manifest(ctx)

    plan = plan_pyrosetta_real_pbs(
        ctx,
        profile="production",
        job_manifest=manifest,
        nodes=["node04", "node05", "node06"],
        ppn=32,
    )

    chunks = read_csv(plan.chunk_manifest_csv)
    assert plan.status == "PASS"
    assert plan.counts["models_per_chunk"] == 0
    assert len(chunks) == 8
    assert {row["model_count"] for row in chunks} <= {"25"}


def test_m2_2c_safe_submit_helper_is_explicit_opt_in():
    script = REPO_ROOT / "fresh" / "scripts" / "submit_m2_pyrosetta_real_jobs.py"
    text = script.read_text(encoding="utf-8")

    assert "--i-understand-this-submits-hpc-jobs" in text
    assert "default=False" in text

    proc = subprocess.run(
        [sys.executable, str(script), "--help"],
        cwd=str(REPO_ROOT),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    assert proc.returncode == 0
    out = proc.stdout.decode("utf-8", errors="replace")
    assert "--models-per-chunk" in out
