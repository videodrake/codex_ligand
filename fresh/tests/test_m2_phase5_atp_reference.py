import csv
import json
import os
import subprocess
import sys
import uuid
from pathlib import Path

import pytest

from egfr_myo1d.core.logging_utils import initialize_logs
from egfr_myo1d.core.run_context import RunContext, RunContextError
from egfr_myo1d.io.residue_mapping import write_residue_mapping
from egfr_myo1d.pocket.atp_reference import (
    CENTROID_FIELDS,
    MAPPING_FIELDS,
    REFERENCE_FIELDS,
    build_atp_reference,
)


REPO_ROOT = Path(__file__).resolve().parents[2]


def unique_run_id(prefix="pytest_m2_phase5"):
    return "{0}_{1}".format(prefix, uuid.uuid4().hex[:12])


def make_tmp_run_context(tmp_path):
    run_id = unique_run_id()
    run_dir = tmp_path / run_id
    ctx = RunContext(
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
    initialize_logs(ctx)
    return ctx


def read_csv(path):
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def pdb_line(record, serial, atom_name, resname, chain, resseq, x, y, z, element="C"):
    return (
        "{record:<6}{serial:5d} {atom_name:<4} {resname:>3} {chain:1}"
        "{resseq:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 10.00          {element:>2}\n"
    ).format(
        record=record,
        serial=serial,
        atom_name=atom_name,
        resname=resname,
        chain=chain,
        resseq=resseq,
        x=x,
        y=y,
        z=z,
        element=element,
    )


def write_receptor_pdb(path, residues=(745, 855), resnames=None):
    resnames = resnames or {745: "LEU", 855: "VAL"}
    lines = []
    serial = 1
    for chain, y_offset in [("A", 0.0), ("B", 10.0)]:
        for residue in residues:
            resname = resnames.get(residue, "ALA")
            x = float(residue - 740)
            lines.append(pdb_line("ATOM", serial, "N", resname, chain, residue, x, y_offset, 0.0, "N"))
            serial += 1
            lines.append(pdb_line("ATOM", serial, "CA", resname, chain, residue, x + 1.0, y_offset, 0.0, "C"))
            serial += 1
            lines.append(pdb_line("ATOM", serial, "C", resname, chain, residue, x + 2.0, y_offset, 0.0, "C"))
            serial += 1
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(lines) + "END\n", encoding="utf-8")
    return path


def write_reference_ligand_pdb(path):
    lines = [
        pdb_line("ATOM", 1, "CA", "LEU", "A", 745, 1.0, 0.0, 0.0, "C"),
        pdb_line("ATOM", 2, "CB", "LEU", "A", 745, 2.0, 0.0, 0.0, "C"),
        pdb_line("ATOM", 3, "CA", "VAL", "A", 855, 40.0, 0.0, 0.0, "C"),
        pdb_line("HETATM", 4, "P", "ANP", "X", 1, 1.5, 0.0, 0.0, "P"),
        pdb_line("HETATM", 5, "O1", "ANP", "X", 1, 2.5, 0.0, 0.0, "O"),
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(lines) + "END\n", encoding="utf-8")
    return path


def write_runtime_numbered_reference_ligand_pdb(path):
    lines = [
        pdb_line("ATOM", 1, "CA", "LEU", "B", 1745, 1.0, 0.0, 0.0, "C"),
        pdb_line("ATOM", 2, "CB", "LEU", "B", 1745, 2.0, 0.0, 0.0, "C"),
        pdb_line("HETATM", 3, "P", "ANP", "X", 1, 1.5, 0.0, 0.0, "P"),
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(lines) + "END\n", encoding="utf-8")
    return path


def write_mapping(ctx, state_id="EGFR_160-185", residues=(745, 855), resnames=None):
    resnames = resnames or {745: "LEU", 855: "VAL"}
    rows = []
    for protomer, chain, offset in [("A", "A", 0), ("B", "B", 1000)]:
        for residue in residues:
            rows.append(
                {
                    "state": state_id,
                    "source_file": "{0}.pdb".format(state_id),
                    "protomer_id": protomer,
                    "source_chain": chain,
                    "source_resseq": residue,
                    "source_icode": "",
                    "source_resname": resnames.get(residue, "ALA"),
                    "runtime_chain": chain,
                    "runtime_resseq": residue + offset,
                    "atom_count": 3,
                    "role": "dockable",
                }
            )
    write_residue_mapping(ctx.qc_dir / "{0}_receptor_mapping.csv".format(state_id), rows, ctx)


def make_synthetic_state(ctx, state_id="EGFR_160-185", residues=(745, 855), resnames=None):
    ctx.create_directories()
    write_mapping(ctx, state_id, residues=residues, resnames=resnames)
    return write_receptor_pdb(
        ctx.run_dir
        / "normalized"
        / "receptors"
        / "{0}_dockable_669_1014_explicit_AB.pdb".format(state_id),
        residues=residues,
        resnames=resnames,
    )


def write_config(path, residues=None, allowed=True, minimum_residue_count=1):
    residues = residues if residues is not None else [
        {"uniprot_residue_number": 745, "residue_name": "LEU", "atp_site_region": "hinge"},
        {"uniprot_residue_number": 855, "residue_name": "VAL", "atp_site_region": "dfg_or_activation_loop"},
    ]
    ligand = "allowed_ligand_resnames:\n  - ANP\n" if allowed else "allowed_ligand_resnames: []\n"
    residue_text = "\n".join(
        [
            "    - uniprot_residue_number: {0}\n      residue_name: {1}\n      atp_site_region: {2}".format(
                item["uniprot_residue_number"],
                item.get("residue_name", ""),
                item.get("atp_site_region", "unknown_or_configured"),
            )
            for item in residues
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        (
            "schema_version: test_atp_reference\n"
            "{ligand}"
            "ligand_neighbor_distance_angstrom: 5.0\n"
            "reference_radius_angstrom: 10.0\n"
            "configured_residue_set:\n"
            "  source_id: synthetic_test_residue_set\n"
            "  evidence_source: fresh/tests/test_m2_phase5_atp_reference.py\n"
            "  minimum_residue_count: {minimum_residue_count}\n"
            "  residues:\n"
            "{residues}\n"
            "gate_policy:\n"
            "  atp_overlap_fraction_reject_threshold: 0.15\n"
            "  atp_centroid_distance_reject_threshold_angstrom: 10.0\n"
        ).format(
            ligand=ligand,
            minimum_residue_count=minimum_residue_count,
            residues=residue_text,
        ),
        encoding="utf-8",
    )
    return path


def test_m2_5_cli_smoke_runs_without_science_engines():
    run_id = unique_run_id("m2_5_cli_smoke")
    ctx = RunContext.create(run_id, "smoke_input", repo_root=REPO_ROOT)
    make_synthetic_state(ctx)
    config = write_config(ctx.run_dir / "synthetic_atp_reference.yaml")

    proc = subprocess.run(
        [
            sys.executable,
            "-m",
            "egfr_myo1d.cli",
            "build-m2-atp-reference",
            "--run-id",
            run_id,
            "--config",
            str(config),
            "--states",
            "EGFR_160-185",
            "--include-reference-states",
            "false",
            "--synthetic-fixture",
            "true",
        ],
        cwd=str(REPO_ROOT),
        env=_cli_env(),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    assert proc.returncode == 0
    assert (ctx.run_dir / "phase2_pockets" / "atp_reference" / "atp_site_reference.csv").is_file()
    assert b"build-m2-atp-reference" in proc.stdout


def test_m2_5_schema_regression_for_required_outputs(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    make_synthetic_state(ctx)
    config = write_config(ctx.run_dir / "atp_reference.yaml")
    report = build_atp_reference(
        ctx,
        states=["EGFR_160-185"],
        include_reference_states=False,
        config_path=config,
        synthetic_fixture=True,
    )

    assert report.status == "PASS_WITH_WARNINGS"
    for path, fields in [
        (report.reference_csv, REFERENCE_FIELDS),
        (report.residue_mapping_csv, MAPPING_FIELDS),
        (report.centroid_csv, CENTROID_FIELDS),
    ]:
        with path.open("r", encoding="utf-8", newline="") as handle:
            assert next(csv.reader(handle)) == fields


def test_m2_5_ligand_bearing_synthetic_reference_detects_anp_and_centroid(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    make_synthetic_state(ctx)
    config = write_config(ctx.run_dir / "ligand_only_atp_reference.yaml", residues=[])
    reference_pdb = write_reference_ligand_pdb(ctx.run_dir / "reference" / "anp_reference.pdb")

    report = build_atp_reference(
        ctx,
        states=["EGFR_160-185"],
        include_reference_states=False,
        config_path=config,
        reference_pdb=reference_pdb,
        synthetic_fixture=True,
    )

    assert report.status == "PASS_WITH_WARNINGS"
    refs = read_csv(report.reference_csv)
    assert len(refs) == 1
    assert refs[0]["ligand_resname"] == "ANP"
    assert refs[0]["reference_mode"] == "synthetic_fixture"
    centroids = read_csv(report.centroid_csv)
    assert any(row["centroid_method"] == "ligand_heavy_atom_centroid" for row in centroids)


def test_m2_5_ligand_reference_restores_runtime_residue_through_mapping(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    make_synthetic_state(ctx)
    config = write_config(ctx.run_dir / "ligand_only_atp_reference.yaml", residues=[])
    reference_pdb = write_runtime_numbered_reference_ligand_pdb(
        ctx.run_dir / "reference" / "runtime_numbered_anp_reference.pdb"
    )

    report = build_atp_reference(
        ctx,
        states=["EGFR_160-185"],
        include_reference_states=False,
        config_path=config,
        reference_pdb=reference_pdb,
        synthetic_fixture=True,
    )

    assert report.status == "PASS_WITH_WARNINGS"
    refs = read_csv(report.reference_csv)
    assert len(refs) == 1
    assert refs[0]["uniprot_residue_number"] == "745"
    assert refs[0]["residue_name"] == "LEU"
    assert "m1_mapping_restored_ligand_neighbor_identity" in refs[0]["warnings"]


def test_m2_5_configured_residue_set_maps_chain_a_and_b_separately(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    make_synthetic_state(ctx)
    config = write_config(ctx.run_dir / "configured_atp_reference.yaml")

    report = build_atp_reference(
        ctx,
        states=["EGFR_160-185"],
        include_reference_states=False,
        config_path=config,
        synthetic_fixture=True,
    )

    assert report.status == "PASS_WITH_WARNINGS"
    mappings = read_csv(report.residue_mapping_csv)
    assert {row["egfr_chain_id"] for row in mappings} == {"A", "B"}
    assert {row["egfr_protomer_id"] for row in mappings} == {"A", "B"}
    assert len(read_csv(report.centroid_csv)) == 2


def test_m2_5_missing_atp_evidence_fails_without_fabricated_rows(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    make_synthetic_state(ctx)
    config = write_config(ctx.run_dir / "empty_atp_reference.yaml", residues=[], allowed=False)

    report = build_atp_reference(
        ctx,
        states=["EGFR_160-185"],
        include_reference_states=False,
        config_path=config,
    )

    assert report.status == "FAIL"
    assert report.reference_count == 0
    assert read_csv(report.reference_csv) == []
    assert "no_atp_reference_evidence_ligand_or_configured_residue_set" in report.blockers


def test_m2_5_residue_drift_or_missing_mapping_is_reported(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    make_synthetic_state(ctx, residues=(746, 855))
    config = write_config(ctx.run_dir / "atp_reference.yaml")

    report = build_atp_reference(
        ctx,
        states=["EGFR_160-185"],
        include_reference_states=False,
        config_path=config,
    )

    assert report.status == "PASS_WITH_WARNINGS"
    mappings = read_csv(report.residue_mapping_csv)
    assert any(row["mapping_status"] == "MISSING" for row in mappings)
    assert any("state_missing_atp_residue_mappings" in warning for warning in report.warnings)


def test_m2_5_missing_primary_state_centroid_is_blocking(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    make_synthetic_state(ctx, residues=(746, 856))
    config = write_config(ctx.run_dir / "atp_reference.yaml")

    report = build_atp_reference(
        ctx,
        states=["EGFR_160-185"],
        include_reference_states=False,
        config_path=config,
    )

    assert report.status == "FAIL"
    assert "no_primary_state_atp_centroid_available" in report.blockers
    assert read_csv(report.centroid_csv) == []


def test_m2_5_sparse_configured_fallback_blocks_production_progression(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    make_synthetic_state(ctx)
    config = write_config(
        ctx.run_dir / "sparse_atp_reference.yaml",
        residues=[
            {"uniprot_residue_number": 745, "residue_name": "LEU", "atp_site_region": "hinge"},
            {"uniprot_residue_number": 855, "residue_name": "VAL", "atp_site_region": "dfg_or_activation_loop"},
        ],
        minimum_residue_count=3,
    )

    report = build_atp_reference(
        ctx,
        states=["EGFR_160-185"],
        include_reference_states=False,
        config_path=config,
    )

    assert report.status == "FAIL"
    assert any("configured_atp_residue_set_insufficient_for_hard_gate" in item for item in report.blockers)


def test_m2_5_residue_name_mismatch_is_reported(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    make_synthetic_state(ctx, resnames={745: "ALA", 855: "VAL"})
    config = write_config(ctx.run_dir / "atp_reference.yaml")

    report = build_atp_reference(
        ctx,
        states=["EGFR_160-185"],
        include_reference_states=False,
        config_path=config,
    )

    assert report.status == "PASS_WITH_WARNINGS"
    mappings = read_csv(report.residue_mapping_csv)
    assert any(row["warnings"] == "residue_name_mismatch" for row in mappings)
    assert any("state_residue_name_mismatches" in warning for warning in report.warnings)


def test_m2_5_rejects_reference_pdb_outside_fresh_or_run_dir(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    make_synthetic_state(ctx)
    config = write_config(ctx.run_dir / "atp_reference.yaml")
    outside = tmp_path / "outside_reference.pdb"
    write_reference_ligand_pdb(outside)

    with pytest.raises(RunContextError):
        build_atp_reference(
            ctx,
            states=["EGFR_160-185"],
            include_reference_states=False,
            config_path=config,
            reference_pdb=outside,
        )


def test_m2_5_non_goal_no_pocket_discovery_outputs_or_engine_calls(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    make_synthetic_state(ctx)
    config = write_config(ctx.run_dir / "atp_reference.yaml")
    report = build_atp_reference(
        ctx,
        states=["EGFR_160-185"],
        include_reference_states=False,
        config_path=config,
    )
    assert report.status == "PASS_WITH_WARNINGS"
    assert not (ctx.run_dir / "phase2_pockets" / "fpocket").exists()
    assert not (ctx.run_dir / "phase2_pockets" / "pocket_candidates").exists()

    source = (REPO_ROOT / "fresh/src/egfr_myo1d/pocket/atp_reference.py").read_text(
        encoding="utf-8"
    )
    for forbidden in [
        "import pyrosetta",
        "import subprocess",
        "os.system",
        "qsub ",
        "sbatch ",
        "vina ",
    ]:
        assert forbidden not in source
    payload = read_json(report.status_json)
    assert payload["fpocket_executed"] is False
    assert payload["pocket_discovery_executed"] is False


def _cli_env():
    env = os.environ.copy()
    env["PYTHONPATH"] = (
        str(REPO_ROOT / "fresh" / "src") + os.pathsep + env.get("PYTHONPATH", "")
    )
    return env
