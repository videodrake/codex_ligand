import csv
import hashlib
import json
import subprocess
from pathlib import Path

from egfr_myo1d.compound.confidentiality import PRIVATE_PATH_REDACTED
from egfr_myo1d.compound.ligand_prepare import ToolStatus, run_m3_ligand_prepare
from egfr_myo1d.core.logging_utils import initialize_logs
from egfr_myo1d.core.run_context import RunContext


PDBQT_TEXT = "ATOM      1  C   LIG A   1       0.000   0.000   0.000  0.00  0.00     0.000 C\nEND\n"


def make_ctx(tmp_path, run_id="m3_ligand_prepare"):
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
    (ctx.repo_root / ".gitignore").write_text(
        "\n".join(
            [
                "fresh/data/raw/ligands/*",
                "fresh/data/private/*",
                "fresh/data/private/**",
                "fresh/runs/*/phase3_compounds/prepared_ligands/*",
                "fresh/runs/*/phase3_compounds/prepared_ligands/**",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    return ctx


def write_csv(path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def read_csv(path):
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def ligand_path(ctx, public_id, suffix, text):
    root = ctx.fresh_root / "data" / "raw" / "ligands"
    root.mkdir(parents=True, exist_ok=True)
    path = root / f"{public_id}{suffix}"
    path.write_text(text, encoding="utf-8")
    return path


def private_map(ctx):
    path = ctx.fresh_root / "data" / "private" / "compound_id_map.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "public_id,internal_id,notes\nCpd-A,INTERNAL_ALPHA,synthetic\nCpd-B,INTERNAL_BETA,synthetic\nCpd-C,INTERNAL_GAMMA,synthetic\n",
        encoding="utf-8",
    )
    return path


def write_m3_t1(ctx, rows, status="PASS", allowed=True, leaks=0):
    manifest = ctx.run_dir / "phase3_compounds" / "manifests" / "ligand_manifest_public.csv"
    fields = [
        "compound_public_id",
        "input_file",
        "input_format",
        "file_exists",
        "file_sha256",
        "size_bytes",
        "has_private_mapping",
        "resolved_by_private_mapping",
        "allowed_for_public_output",
        "selected_for_preparation",
        "alternative_file_count",
        "tracked_by_git",
        "gitignore_protected",
        "qc_status",
        "qc_notes",
        "supplied_pdbqt_manual_review",
    ]
    write_csv(manifest, fields, rows)
    qc = ctx.run_dir / "phase3_compounds" / "qc"
    qc.mkdir(parents=True, exist_ok=True)
    (qc / "m3_ligand_qc_summary.json").write_text(
        json.dumps(
            {
                "overall_status": status,
                "m3_t2_allowed": allowed,
                "counts": {"confidentiality_leaks": leaks},
            }
        ),
        encoding="utf-8",
    )


def row_for(ctx, public_id, path, private=False):
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    return {
        "compound_public_id": public_id,
        "input_file": PRIVATE_PATH_REDACTED if private else ctx.relative_to_repo(path),
        "input_format": path.suffix.lower().lstrip("."),
        "file_exists": "true",
        "file_sha256": digest,
        "size_bytes": str(path.stat().st_size),
        "has_private_mapping": "true" if private else "false",
        "resolved_by_private_mapping": "true" if private else "false",
        "allowed_for_public_output": "false" if private else "true",
        "selected_for_preparation": "true",
        "alternative_file_count": "0",
        "tracked_by_git": "false",
        "gitignore_protected": "true",
        "qc_status": "PASS",
        "qc_notes": "",
        "supplied_pdbqt_manual_review": "false",
    }


def fake_tools(openbabel=True, rdkit=False):
    return ToolStatus(
        rdkit_available=rdkit,
        rdkit_version="test-rdkit" if rdkit else None,
        openbabel_available=openbabel,
        openbabel_version="Open Babel test",
        meeko_available=False,
        meeko_version=None,
        obabel_binary="obabel",
    )


def patch_fake_conversion(monkeypatch):
    def fake_convert(source, output, public_id, ctx, tools, source_display):
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(PDBQT_TEXT, encoding="utf-8")
        return True, "fake Open Babel conversion wrote valid PDBQT"

    monkeypatch.setattr("egfr_myo1d.compound.ligand_prepare._convert_with_obabel", fake_convert)
    monkeypatch.setattr("egfr_myo1d.compound.ligand_prepare.detect_tools", lambda obabel_binary="obabel": fake_tools(openbabel=True))


def test_m3_ligand_prepare_supplied_pdbqt_requires_manual_review(tmp_path):
    ctx = make_ctx(tmp_path)
    private_map(ctx)
    paths = [ligand_path(ctx, public_id, ".pdbqt", PDBQT_TEXT) for public_id in ["Cpd-A", "Cpd-B", "Cpd-C"]]
    write_m3_t1(ctx, [row_for(ctx, public_id, path) for public_id, path in zip(["Cpd-A", "Cpd-B", "Cpd-C"], paths)])

    result = run_m3_ligand_prepare(ctx, mode="prepare")

    assert result.status == "FAIL"
    assert any("supplied PDBQT requires explicit manual-review approval" in blocker for blocker in result.blockers)


def test_m3_ligand_prepare_supplied_pdbqt_is_reused_after_manual_review(tmp_path):
    ctx = make_ctx(tmp_path)
    private_map(ctx)
    paths = [ligand_path(ctx, public_id, ".pdbqt", PDBQT_TEXT) for public_id in ["Cpd-A", "Cpd-B", "Cpd-C"]]
    rows = [row_for(ctx, public_id, path) for public_id, path in zip(["Cpd-A", "Cpd-B", "Cpd-C"], paths)]
    for row in rows:
        row["supplied_pdbqt_manual_review"] = "true"
    write_m3_t1(ctx, rows)

    result = run_m3_ligand_prepare(ctx, mode="prepare")
    output_rows = read_csv(result.preparation_manifest_csv)

    assert result.status == "WARN"
    assert result.counts["ligands_reused_supplied_pdbqt"] == 3
    assert all(row["preparation_method"] == "supplied_pdbqt_validated" for row in output_rows)
    for public_id in ["Cpd-A", "Cpd-B", "Cpd-C"]:
        assert (ctx.run_dir / "phase3_compounds" / "prepared_ligands" / public_id / f"{public_id}.pdbqt").is_file()


def test_m3_ligand_prepare_conversion_route_records_without_modifying_raw(tmp_path, monkeypatch):
    ctx = make_ctx(tmp_path)
    private_map(ctx)
    patch_fake_conversion(monkeypatch)
    paths = [ligand_path(ctx, public_id, ".sdf", "mol\n$$$$\n") for public_id in ["Cpd-A", "Cpd-B", "Cpd-C"]]
    before = {path: path.read_text(encoding="utf-8") for path in paths}
    write_m3_t1(ctx, [row_for(ctx, public_id, path) for public_id, path in zip(["Cpd-A", "Cpd-B", "Cpd-C"], paths)])

    result = run_m3_ligand_prepare(ctx, mode="prepare")
    rows = read_csv(result.preparation_manifest_csv)

    assert result.status == "PASS"
    assert all(row["preparation_method"] == "openbabel_conversion" for row in rows)
    assert all(path.read_text(encoding="utf-8") == before[path] for path in paths)


def test_m3_ligand_prepare_blocks_source_hash_drift_after_m3_t1(tmp_path, monkeypatch):
    ctx = make_ctx(tmp_path)
    private_map(ctx)
    patch_fake_conversion(monkeypatch)
    paths = [ligand_path(ctx, public_id, ".sdf", "mol\n$$$$\n") for public_id in ["Cpd-A", "Cpd-B", "Cpd-C"]]
    write_m3_t1(ctx, [row_for(ctx, public_id, path) for public_id, path in zip(["Cpd-A", "Cpd-B", "Cpd-C"], paths)])
    paths[0].write_text("changed\n$$$$\n", encoding="utf-8")

    result = run_m3_ligand_prepare(ctx, mode="prepare")

    assert result.status == "FAIL"
    assert any("source ligand SHA256 differs from M3-T1 manifest" in blocker for blocker in result.blockers)


def test_m3_ligand_prepare_smiles_input_never_appears_in_outputs(tmp_path, monkeypatch):
    ctx = make_ctx(tmp_path)
    private_map(ctx)
    patch_fake_conversion(monkeypatch)
    paths = [ligand_path(ctx, public_id, ".smi", "CCO\n") for public_id in ["Cpd-A", "Cpd-B", "Cpd-C"]]
    write_m3_t1(ctx, [row_for(ctx, public_id, path) for public_id, path in zip(["Cpd-A", "Cpd-B", "Cpd-C"], paths)])

    result = run_m3_ligand_prepare(ctx, mode="prepare")
    combined = "\n".join(
        path.read_text(encoding="utf-8", errors="ignore")
        for path in [
            result.preparation_manifest_csv,
            result.prepared_hashes_csv,
            result.prep_qc_csv,
            result.summary_json,
            result.report_md,
            result.phase3_log,
        ]
    )

    assert result.status == "WARN"
    assert "CCO" not in combined


def test_m3_ligand_prepare_missing_inputs_fail_prepare_and_warn_dry_run(tmp_path):
    ctx = make_ctx(tmp_path)
    private_map(ctx)
    write_m3_t1(ctx, [])

    prepare = run_m3_ligand_prepare(ctx, mode="prepare")
    dry_ctx = make_ctx(tmp_path / "dry")
    private_map(dry_ctx)
    write_m3_t1(dry_ctx, [])
    dry = run_m3_ligand_prepare(dry_ctx, mode="dry-run")

    assert prepare.status == "FAIL"
    assert dry.status == "WARN"
    assert dry.counts["ligands_prepared"] == 0
    assert any("ligand input file is missing" in warning for warning in dry.warnings)


def test_m3_ligand_prepare_unsupported_and_zero_byte_inputs_fail(tmp_path):
    ctx = make_ctx(tmp_path)
    private_map(ctx)
    a = ligand_path(ctx, "Cpd-A", ".xyz", "x\n")
    b = ligand_path(ctx, "Cpd-B", ".sdf", "")
    c = ligand_path(ctx, "Cpd-C", ".pdbqt", PDBQT_TEXT)
    write_m3_t1(ctx, [row_for(ctx, "Cpd-A", a), row_for(ctx, "Cpd-B", b), row_for(ctx, "Cpd-C", c)])

    result = run_m3_ligand_prepare(ctx, mode="prepare")

    assert result.status == "FAIL"
    assert any("unsupported input format" in blocker for blocker in result.blockers)
    assert any("zero-byte ligand input" in blocker for blocker in result.blockers)


def test_m3_ligand_prepare_openbabel_unavailable_blocks_conversion_not_pdbqt(tmp_path, monkeypatch):
    ctx = make_ctx(tmp_path)
    private_map(ctx)
    monkeypatch.setattr("egfr_myo1d.compound.ligand_prepare.detect_tools", lambda obabel_binary="obabel": fake_tools(openbabel=False))
    a = ligand_path(ctx, "Cpd-A", ".sdf", "mol\n$$$$\n")
    b = ligand_path(ctx, "Cpd-B", ".pdbqt", PDBQT_TEXT)
    c = ligand_path(ctx, "Cpd-C", ".pdbqt", PDBQT_TEXT)
    rows = [row_for(ctx, "Cpd-A", a), row_for(ctx, "Cpd-B", b), row_for(ctx, "Cpd-C", c)]
    for row in rows[1:]:
        row["supplied_pdbqt_manual_review"] = "true"
    write_m3_t1(ctx, rows)

    result = run_m3_ligand_prepare(ctx, mode="prepare")

    assert result.status == "FAIL"
    assert any("Open Babel unavailable" in blocker for blocker in result.blockers)
    assert (ctx.run_dir / "phase3_compounds" / "prepared_ligands" / "Cpd-B" / "Cpd-B.pdbqt").is_file()


def test_m3_ligand_prepare_multiple_sdf_records_warn(tmp_path, monkeypatch):
    ctx = make_ctx(tmp_path)
    private_map(ctx)
    patch_fake_conversion(monkeypatch)
    paths = [ligand_path(ctx, public_id, ".sdf", "mol1\n$$$$\nmol2\n$$$$\n") for public_id in ["Cpd-A", "Cpd-B", "Cpd-C"]]
    write_m3_t1(ctx, [row_for(ctx, public_id, path) for public_id, path in zip(["Cpd-A", "Cpd-B", "Cpd-C"], paths)])

    result = run_m3_ligand_prepare(ctx, mode="prepare")
    rows = read_csv(result.preparation_manifest_csv)

    assert result.status == "WARN"
    assert all(row["multi_molecule_input"] == "true" for row in rows)


def test_m3_ligand_prepare_private_path_redacted(tmp_path, monkeypatch):
    ctx = make_ctx(tmp_path)
    patch_fake_conversion(monkeypatch)
    pmap = private_map(ctx)
    private_ligands = ctx.fresh_root / "data" / "private" / "ligands"
    private_ligands.mkdir(parents=True, exist_ok=True)
    for token in ["INTERNAL_ALPHA", "INTERNAL_BETA", "INTERNAL_GAMMA"]:
        (private_ligands / f"{token}.sdf").write_text("mol\n$$$$\n", encoding="utf-8")
    rows = [
        {
            **row_for(ctx, public_id, private_ligands / f"{token}.sdf", private=True),
            "input_file": PRIVATE_PATH_REDACTED,
        }
        for public_id, token in zip(["Cpd-A", "Cpd-B", "Cpd-C"], ["INTERNAL_ALPHA", "INTERNAL_BETA", "INTERNAL_GAMMA"])
    ]
    write_m3_t1(ctx, rows)

    result = run_m3_ligand_prepare(ctx, mode="prepare", private_map=pmap)
    text = result.preparation_manifest_csv.read_text(encoding="utf-8") + result.summary_json.read_text(encoding="utf-8")

    assert result.status == "PASS"
    assert PRIVATE_PATH_REDACTED in text
    assert "INTERNAL_ALPHA" not in text


def test_m3_ligand_prepare_internal_id_leak_fails_without_printing_token(tmp_path, monkeypatch):
    ctx = make_ctx(tmp_path)
    private_map(ctx)
    patch_fake_conversion(monkeypatch)
    paths = [ligand_path(ctx, public_id, ".sdf", "mol\n$$$$\n") for public_id in ["Cpd-A", "Cpd-B", "Cpd-C"]]
    write_m3_t1(ctx, [row_for(ctx, public_id, path) for public_id, path in zip(["Cpd-A", "Cpd-B", "Cpd-C"], paths)])
    reports = ctx.run_dir / "phase3_compounds" / "reports"
    reports.mkdir(parents=True, exist_ok=True)
    (reports / "leaky.md").write_text("INTERNAL_ALPHA\n", encoding="utf-8")

    result = run_m3_ligand_prepare(ctx, mode="prepare")
    public_text = result.preparation_manifest_csv.read_text(encoding="utf-8") + result.summary_json.read_text(encoding="utf-8")

    assert result.status == "FAIL"
    assert "INTERNAL_ALPHA" not in public_text
    assert any("internal ID leakage" in blocker for blocker in result.blockers)


def test_m3_ligand_prepare_prepared_outputs_only_under_run_dir_and_no_non_goals(tmp_path, monkeypatch):
    ctx = make_ctx(tmp_path)
    private_map(ctx)
    patch_fake_conversion(monkeypatch)
    paths = [ligand_path(ctx, public_id, ".sdf", "mol\n$$$$\n") for public_id in ["Cpd-A", "Cpd-B", "Cpd-C"]]
    write_m3_t1(ctx, [row_for(ctx, public_id, path) for public_id, path in zip(["Cpd-A", "Cpd-B", "Cpd-C"], paths)])

    result = run_m3_ligand_prepare(ctx, mode="prepare")
    summary = read_json(result.summary_json)

    for public_id in ["Cpd-A", "Cpd-B", "Cpd-C"]:
        prepared = ctx.run_dir / "phase3_compounds" / "prepared_ligands" / public_id / f"{public_id}.pdbqt"
        prepared.resolve().relative_to(ctx.run_dir.resolve())
        assert prepared.is_file()
    assert all(summary["non_goals_preserved"].values())
    assert not (ctx.run_dir / "phase3_compounds" / "receptor_pdbqt").exists()
    assert not (ctx.run_dir / "phase3_compounds" / "docking_boxes").exists()
    assert not (ctx.run_dir / "phase3_compounds" / "vina_raw").exists()
    assert not (ctx.run_dir / "phase3_compounds" / "tables" / "final_m3_candidate_hypotheses.csv").exists()


def test_m3_ligand_prepare_tracked_sensitive_files_fail(tmp_path):
    ctx = make_ctx(tmp_path)
    pmap = private_map(ctx)
    a = ligand_path(ctx, "Cpd-A", ".pdbqt", PDBQT_TEXT)
    b = ligand_path(ctx, "Cpd-B", ".pdbqt", PDBQT_TEXT)
    c = ligand_path(ctx, "Cpd-C", ".pdbqt", PDBQT_TEXT)
    write_m3_t1(ctx, [row_for(ctx, "Cpd-A", a), row_for(ctx, "Cpd-B", b), row_for(ctx, "Cpd-C", c)])
    subprocess.run(["git", "init"], cwd=ctx.repo_root, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    subprocess.run(["git", "add", "-f", "fresh/data/raw/ligands/Cpd-A.pdbqt", "fresh/data/private/compound_id_map.csv"], cwd=ctx.repo_root, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)

    result = run_m3_ligand_prepare(ctx, mode="prepare", private_map=pmap)

    assert result.status == "FAIL"
    assert any("raw ligand file is tracked by git" in blocker for blocker in result.blockers)
    assert any("private compound_id_map.csv is tracked by git" in blocker for blocker in result.blockers)


def test_m3_ligand_prepare_m3_t1_blocks_prepare_mode(tmp_path):
    ctx = make_ctx(tmp_path)
    private_map(ctx)
    paths = [ligand_path(ctx, public_id, ".pdbqt", PDBQT_TEXT) for public_id in ["Cpd-A", "Cpd-B", "Cpd-C"]]
    write_m3_t1(ctx, [row_for(ctx, public_id, path) for public_id, path in zip(["Cpd-A", "Cpd-B", "Cpd-C"], paths)], status="FAIL", allowed=False)

    result = run_m3_ligand_prepare(ctx, mode="prepare")

    assert result.status == "FAIL"
    assert any("M3-T1 overall_status is FAIL" in blocker for blocker in result.blockers)
