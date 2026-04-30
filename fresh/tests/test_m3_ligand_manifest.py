import csv
import json
import subprocess
from pathlib import Path

from egfr_myo1d.compound.confidentiality import PRIVATE_PATH_REDACTED
from egfr_myo1d.compound.ligand_manifest import run_m3_ligand_qc
from egfr_myo1d.core.logging_utils import initialize_logs
from egfr_myo1d.core.run_context import RunContext


def make_ctx(tmp_path, run_id="m3_ligand_qc"):
    repo = tmp_path / "repo"
    fresh = repo / "fresh"
    run_dir = fresh / "runs" / run_id
    ctx = RunContext(
        repo_root=repo,
        fresh_root=fresh,
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
    initialize_logs(ctx)
    return ctx


def write_smiles(ctx, public_id, text="C\n", suffix=".smi"):
    ligand_dir = ctx.fresh_root / "data" / "raw" / "ligands"
    ligand_dir.mkdir(parents=True, exist_ok=True)
    path = ligand_dir / f"{public_id}{suffix}"
    path.write_text(text, encoding="utf-8")
    return path


def write_all_public_ligands(ctx):
    for public_id in ["Cpd-A", "Cpd-B", "Cpd-C"]:
        write_smiles(ctx, public_id)


def write_private_map(ctx):
    private_dir = ctx.fresh_root / "data" / "private"
    private_dir.mkdir(parents=True, exist_ok=True)
    path = private_dir / "compound_id_map.csv"
    path.write_text(
        "public_id,internal_id,notes\nCpd-A,INTERNAL_ALPHA,synthetic\nCpd-B,INTERNAL_BETA,synthetic\nCpd-C,INTERNAL_GAMMA,synthetic\n",
        encoding="utf-8",
    )
    return path


def write_m3_t0_summary(ctx):
    qc_dir = ctx.run_dir / "phase3_compounds" / "qc"
    qc_dir.mkdir(parents=True, exist_ok=True)
    (qc_dir / "m3_readiness_summary.json").write_text(
        '{"overall_status":"PASS","m3_docking_allowed":true}\n',
        encoding="utf-8",
    )


def write_m3_t0_fail_summary(ctx):
    qc_dir = ctx.run_dir / "phase3_compounds" / "qc"
    qc_dir.mkdir(parents=True, exist_ok=True)
    (qc_dir / "m3_readiness_summary.json").write_text(
        '{"overall_status":"FAIL","m3_docking_allowed":false}\n',
        encoding="utf-8",
    )


def read_csv(path):
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def test_m3_ligand_qc_valid_public_files_pass_and_allow_t2(tmp_path):
    ctx = make_ctx(tmp_path)
    write_all_public_ligands(ctx)
    write_private_map(ctx)
    write_m3_t0_summary(ctx)

    result = run_m3_ligand_qc(ctx, mode="docking")
    summary = read_json(result.summary_json)
    rows = read_csv(result.manifest_csv)

    assert result.status == "PASS"
    assert summary["m3_t2_allowed"] is True
    assert {row["compound_public_id"] for row in rows} == {"Cpd-A", "Cpd-B", "Cpd-C"}
    assert all(row["file_sha256"] for row in rows)
    assert all(int(row["size_bytes"]) > 0 for row in rows)


def test_m3_ligand_qc_missing_files_warn_in_dry_run_and_fail_in_docking(tmp_path):
    dry_ctx = make_ctx(tmp_path / "dry")
    dry = run_m3_ligand_qc(dry_ctx, mode="dry-run")
    assert dry.status == "WARN"
    assert dry.counts["ligands_missing"] == 3

    docking_ctx = make_ctx(tmp_path / "docking")
    docking = run_m3_ligand_qc(docking_ctx, mode="docking")
    assert docking.status == "FAIL"
    assert docking.counts["ligands_missing"] == 3


def test_m3_ligand_qc_zero_byte_file_fails(tmp_path):
    ctx = make_ctx(tmp_path)
    write_all_public_ligands(ctx)
    (ctx.fresh_root / "data" / "raw" / "ligands" / "Cpd-B.smi").write_text("", encoding="utf-8")

    result = run_m3_ligand_qc(ctx, mode="docking")

    assert result.status == "FAIL"
    assert any("zero-byte" in blocker for blocker in result.blockers)


def test_m3_ligand_qc_unsupported_extension_fails(tmp_path):
    ctx = make_ctx(tmp_path)
    write_smiles(ctx, "Cpd-A")
    write_smiles(ctx, "Cpd-B")
    ligand_dir = ctx.fresh_root / "data" / "raw" / "ligands"
    (ligand_dir / "Cpd-C.txt").write_text("C\n", encoding="utf-8")

    result = run_m3_ligand_qc(ctx, mode="docking")

    assert result.status == "FAIL"
    assert any("unsupported ligand format" in blocker for blocker in result.blockers)


def test_m3_ligand_qc_represents_all_public_ids_when_missing(tmp_path):
    ctx = make_ctx(tmp_path)

    result = run_m3_ligand_qc(ctx, mode="dry-run")
    rows = read_csv(result.manifest_csv)

    assert [row["compound_public_id"] for row in rows] == ["Cpd-A", "Cpd-B", "Cpd-C"]


def test_m3_ligand_qc_multiple_files_warn_and_selects_preferred_then_lexicographic(tmp_path):
    ctx = make_ctx(tmp_path)
    write_all_public_ligands(ctx)
    write_m3_t0_summary(ctx)
    ligand_dir = ctx.fresh_root / "data" / "raw" / "ligands"
    (ligand_dir / "Cpd-A.sdf").write_text("mol\n$$$$\n", encoding="utf-8")
    (ligand_dir / "Cpd-A.2.sdf").write_text("mol\n$$$$\n", encoding="utf-8")

    result = run_m3_ligand_qc(ctx, mode="docking")
    rows = {row["compound_public_id"]: row for row in read_csv(result.manifest_csv)}

    assert result.status == "WARN"
    assert rows["Cpd-A"]["input_format"] == "sdf"
    assert rows["Cpd-A"]["alternative_file_count"] == "2"
    assert "multiple matching ligand files" in rows["Cpd-A"]["qc_notes"]


def test_m3_ligand_qc_m3_t0_failure_blocks_t2_even_with_valid_ligands(tmp_path):
    ctx = make_ctx(tmp_path)
    write_all_public_ligands(ctx)
    write_private_map(ctx)
    write_m3_t0_fail_summary(ctx)

    result = run_m3_ligand_qc(ctx, mode="docking")
    summary = read_json(result.summary_json)

    assert result.status == "FAIL"
    assert summary["m3_t2_allowed"] is False
    assert any("M3-T0 readiness must PASS" in blocker for blocker in result.blockers)


def test_m3_ligand_qc_tracked_alternative_candidate_fails(tmp_path):
    ctx = make_ctx(tmp_path)
    write_all_public_ligands(ctx)
    write_m3_t0_summary(ctx)
    ligand_dir = ctx.fresh_root / "data" / "raw" / "ligands"
    (ligand_dir / "Cpd-A.sdf").write_text("mol\n$$$$\n", encoding="utf-8")
    subprocess.run(["git", "init"], cwd=ctx.repo_root, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    subprocess.run(["git", "add", "-f", "fresh/data/raw/ligands/Cpd-A.smi"], cwd=ctx.repo_root, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    result = run_m3_ligand_qc(ctx, mode="docking")

    assert result.status == "FAIL"
    assert any("alternative or unsupported raw ligand candidate is tracked by git" in blocker for blocker in result.blockers)


def test_m3_ligand_qc_prebuilt_pdbqt_input_warns_until_m3_t2(tmp_path):
    ctx = make_ctx(tmp_path)
    write_smiles(ctx, "Cpd-A")
    write_smiles(ctx, "Cpd-B")
    write_m3_t0_summary(ctx)
    ligand_dir = ctx.fresh_root / "data" / "raw" / "ligands"
    (ligand_dir / "Cpd-C.pdbqt").write_text("REMARK prebuilt pdbqt\n", encoding="utf-8")

    result = run_m3_ligand_qc(ctx, mode="docking")

    assert result.status == "WARN"
    assert result.summary_json.read_text(encoding="utf-8").count("PDBQT supplied as input") >= 1


def test_m3_ligand_qc_private_mapping_redacts_internal_path_and_token(tmp_path):
    ctx = make_ctx(tmp_path)
    write_m3_t0_summary(ctx)
    private_dir = ctx.fresh_root / "data" / "private"
    private_ligands = private_dir / "ligands"
    private_ligands.mkdir(parents=True, exist_ok=True)
    private_map = private_dir / "compound_id_map.csv"
    private_map.write_text(
        "public_id,internal_id,notes\nCpd-A,INTERNAL_ALPHA,synthetic\nCpd-B,INTERNAL_BETA,synthetic\nCpd-C,INTERNAL_GAMMA,synthetic\n",
        encoding="utf-8",
    )
    for token in ["INTERNAL_ALPHA", "INTERNAL_BETA", "INTERNAL_GAMMA"]:
        (private_ligands / f"{token}.smi").write_text("C\n", encoding="utf-8")

    result = run_m3_ligand_qc(ctx, mode="docking", private_map=private_map)
    manifest_text = result.manifest_csv.read_text(encoding="utf-8")
    summary_text = result.summary_json.read_text(encoding="utf-8")

    assert result.status == "PASS"
    assert PRIVATE_PATH_REDACTED in manifest_text
    assert "INTERNAL_ALPHA" not in manifest_text
    assert "INTERNAL_ALPHA" not in summary_text
    rows = read_csv(result.manifest_csv)
    assert all(row["resolved_by_private_mapping"] == "true" for row in rows)


def test_m3_ligand_qc_detects_internal_id_leak_without_printing_token(tmp_path):
    ctx = make_ctx(tmp_path)
    write_all_public_ligands(ctx)
    private_dir = ctx.fresh_root / "data" / "private"
    private_dir.mkdir(parents=True, exist_ok=True)
    private_map = private_dir / "compound_id_map.csv"
    private_map.write_text("public_id,internal_id,notes\nCpd-A,SECRET123,synthetic\n", encoding="utf-8")
    report_dir = ctx.run_dir / "phase3_compounds" / "reports"
    report_dir.mkdir(parents=True, exist_ok=True)
    (report_dir / "leaky.md").write_text("SECRET123\n", encoding="utf-8")

    result = run_m3_ligand_qc(ctx, mode="docking", private_map=private_map)
    audit_text = result.audit_csv.read_text(encoding="utf-8")

    assert result.status == "FAIL"
    assert "INTERNAL_ID_FOR_Cpd-A" in audit_text
    assert "SECRET123" not in audit_text


def test_m3_ligand_qc_gitignore_patterns_are_added_when_missing(tmp_path):
    ctx = make_ctx(tmp_path)
    write_all_public_ligands(ctx)
    write_private_map(ctx)
    write_m3_t0_summary(ctx)
    (ctx.repo_root / ".gitignore").write_text("__pycache__/\n", encoding="utf-8")

    result = run_m3_ligand_qc(ctx, mode="docking", update_ignore=True)
    gitignore_text = (ctx.repo_root / ".gitignore").read_text(encoding="utf-8")

    assert result.status == "PASS"
    assert "fresh/data/raw/ligands/*" in gitignore_text
    assert "fresh/data/private/**" in gitignore_text


def test_m3_ligand_qc_tracked_raw_ligand_and_private_map_fail(tmp_path):
    ctx = make_ctx(tmp_path)
    write_all_public_ligands(ctx)
    private_dir = ctx.fresh_root / "data" / "private"
    private_dir.mkdir(parents=True, exist_ok=True)
    private_map = private_dir / "compound_id_map.csv"
    private_map.write_text("public_id,internal_id,notes\nCpd-A,INTERNAL_ALPHA,synthetic\n", encoding="utf-8")
    subprocess.run(["git", "init"], cwd=ctx.repo_root, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    subprocess.run(["git", "add", "-f", "fresh/data/raw/ligands/Cpd-A.smi", "fresh/data/private/compound_id_map.csv"], cwd=ctx.repo_root, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    result = run_m3_ligand_qc(ctx, mode="docking", private_map=private_map)

    assert result.status == "FAIL"
    assert any("raw ligand file is tracked by git" in blocker for blocker in result.blockers)
    assert any("private compound_id_map.csv is tracked by git" in blocker for blocker in result.blockers)


def test_m3_ligand_qc_non_goals_no_preparation_outputs(tmp_path):
    ctx = make_ctx(tmp_path)
    write_all_public_ligands(ctx)

    result = run_m3_ligand_qc(ctx, mode="docking")
    summary = read_json(result.summary_json)

    assert all(summary["non_goals_preserved"].values())
    assert not list((ctx.run_dir / "phase3_compounds" / "prepared_ligands").rglob("*.pdbqt"))
    assert not (ctx.run_dir / "phase3_compounds" / "receptor_pdbqt").exists()
    assert not (ctx.run_dir / "phase3_compounds" / "docking_inputs").exists()
    assert not (ctx.run_dir / "phase3_compounds" / "docking_outputs").exists()
