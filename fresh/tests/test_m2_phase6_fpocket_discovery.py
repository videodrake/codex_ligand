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
from egfr_myo1d.pocket.family_merge import MERGED_FIELDS, merge_pocket_families
from egfr_myo1d.pocket.fpocket_adapter import (
    FORBIDDEN_M2_7_OUTPUTS,
    RAW_FIELDS,
    discover_and_normalize_pockets,
    parse_fpocket_output,
)
from egfr_myo1d.pocket.normalization import RAW_CANDIDATE_FIELDS


REPO_ROOT = Path(__file__).resolve().parents[2]


def unique_run_id(prefix="pytest_m2_phase6"):
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


def read_header(path):
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        return next(csv.reader(handle))


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


def write_receptor_pdb(path, residues=(745, 855, 900), resnames=None):
    resnames = resnames or {745: "LEU", 855: "VAL", 900: "ASP"}
    lines = []
    serial = 1
    for chain, y_offset in [("A", 0.0), ("B", 40.0)]:
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


def write_mapping(ctx, state_id="EGFR_160-185", residues=(745, 855, 900), resnames=None):
    resnames = resnames or {745: "LEU", 855: "VAL", 900: "ASP"}
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


def make_synthetic_state(ctx, state_id="EGFR_160-185", residues=(745, 855, 900), resnames=None):
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


def write_preconditions(ctx):
    consensus = ctx.run_dir / "output" / "ppi" / "ppi_consensus_patch.csv"
    consensus.parent.mkdir(parents=True, exist_ok=True)
    consensus.write_text("ppi_patch_id,state_id,egfr_chain_id\np1,EGFR_160-185,A\n", encoding="utf-8")
    atp = ctx.run_dir / "phase2_pockets" / "atp_reference" / "atp_site_reference.csv"
    atp.parent.mkdir(parents=True, exist_ok=True)
    atp.write_text("atp_site_id,uniprot_residue_number\natp1,745\n", encoding="utf-8")


def write_fpocket_output(root, state_id, pockets):
    out_dir = root / "{0}_out".format(state_id)
    pockets_dir = out_dir / "pockets"
    pockets_dir.mkdir(parents=True, exist_ok=True)
    info_lines = []
    for index, pocket in enumerate(pockets, start=1):
        info_lines.extend(
            [
                "Pocket {0} :\n".format(index),
                "Score : {0}\n".format(pocket.get("score", 10.0 + index)),
                "Druggability Score : {0}\n".format(pocket.get("druggability", 0.5)),
                "Pocket volume (Monte Carlo) : {0}\n".format(pocket.get("volume", 100.0 + index)),
                "Number of Alpha Spheres : {0}\n".format(pocket.get("alpha_spheres", 20 + index)),
            ]
        )
        lines = []
        serial = 1
        for atom in pocket["atoms"]:
            lines.append(
                pdb_line(
                    "ATOM",
                    serial,
                    atom.get("atom_name", "CA"),
                    atom.get("resname", "LEU"),
                    atom.get("chain", "A"),
                    atom.get("resseq", 745),
                    atom.get("x", 0.0),
                    atom.get("y", 0.0),
                    atom.get("z", 0.0),
                    atom.get("element", "C"),
                )
            )
            serial += 1
        (pockets_dir / "pocket{0}_atm.pdb".format(index)).write_text("".join(lines) + "END\n", encoding="utf-8")
    (out_dir / "pockets_info.txt").write_text("".join(info_lines), encoding="utf-8")
    return out_dir


def single_pocket(chain="A", residue=745, resname="LEU", x=0.0, y=0.0, z=0.0):
    return {
        "score": 42.0,
        "druggability": 0.77,
        "volume": 123.4,
        "alpha_spheres": 31,
        "atoms": [
            {"chain": chain, "resseq": residue, "resname": resname, "atom_name": "CA", "x": x, "y": y, "z": z},
            {"chain": chain, "resseq": residue, "resname": resname, "atom_name": "CB", "x": x + 2.0, "y": y, "z": z},
        ],
    }


def test_m2_6_cli_smoke_parser_only_runs_without_science_engines():
    run_id = unique_run_id("m2_6_cli_smoke")
    ctx = RunContext.create(run_id, "smoke_input", repo_root=REPO_ROOT)
    make_synthetic_state(ctx)
    write_preconditions(ctx)
    output_root = ctx.run_dir / "synthetic_fpocket"
    write_fpocket_output(output_root, "EGFR_160-185", [single_pocket()])

    proc = subprocess.run(
        [
            sys.executable,
            "-m",
            "egfr_myo1d.cli",
            "run-m2-fpocket-discovery",
            "--run-id",
            run_id,
            "--states",
            "EGFR_160-185",
            "--include-reference-states",
            "false",
            "--fpocket-output-root",
            str(output_root),
            "--synthetic-fixture",
            "true",
        ],
        cwd=str(REPO_ROOT),
        env=_cli_env(),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    assert proc.returncode == 0
    assert b"run-m2-fpocket-discovery" in proc.stdout
    assert (ctx.run_dir / "phase2_pockets" / "merged" / "pocket_candidates_raw.csv").is_file()


def test_m2_6_fpocket_parser_extracts_metrics_centroid_and_residues(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    receptor = make_synthetic_state(ctx)
    output_dir = write_fpocket_output(
        ctx.run_dir / "synthetic_fpocket",
        "EGFR_160-185",
        [
            single_pocket(chain="A", residue=745, x=0.0),
            single_pocket(chain="B", residue=855, resname="VAL", x=10.0, y=40.0),
        ],
    )

    rows, warnings = parse_fpocket_output(
        ctx,
        "EGFR_160-185",
        "primary",
        receptor,
        output_dir,
        ctx.qc_dir / "EGFR_160-185_receptor_mapping.csv",
    )

    assert warnings == []
    assert len(rows) == 2
    assert rows[0]["fpocket_score"] == "42.0"
    assert rows[0]["fpocket_druggability_score"] == "0.77"
    assert rows[0]["fpocket_volume"] == "123.4"
    assert rows[0]["num_alpha_spheres"] == "31"
    assert rows[0]["centroid_x"] == "1.000"
    assert rows[0]["raw_residue_ids"] == "A:745"
    assert rows[0]["mapped_uniprot_residues"] == "745"


def test_m2_6_schema_regression_for_required_outputs(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    make_synthetic_state(ctx)
    write_preconditions(ctx)
    output_root = ctx.run_dir / "synthetic_fpocket"
    write_fpocket_output(output_root, "EGFR_160-185", [single_pocket()])

    report = discover_and_normalize_pockets(
        ctx,
        states=["EGFR_160-185"],
        include_reference_states=False,
        fpocket_output_root=output_root,
        synthetic_fixture=True,
    )

    assert report.status == "PASS"
    assert read_header(report.raw_csv_paths[0]) == RAW_FIELDS
    assert read_header(report.raw_candidate_csv) == RAW_CANDIDATE_FIELDS
    assert read_header(report.merged_candidate_csv) == MERGED_FIELDS


def test_m2_6_dimer_receptor_is_staged_inside_run_dir_without_modifying_source(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    receptor = make_synthetic_state(ctx)
    before = receptor.read_text(encoding="utf-8")
    write_preconditions(ctx)
    output_root = ctx.run_dir / "synthetic_fpocket"
    write_fpocket_output(output_root, "EGFR_160-185", [single_pocket()])

    report = discover_and_normalize_pockets(
        ctx,
        states=["EGFR_160-185"],
        include_reference_states=False,
        fpocket_output_root=output_root,
        synthetic_fixture=True,
    )

    staged = ctx.run_dir / "phase2_pockets" / "raw" / "EGFR_160-185_fpocket_pockets" / "EGFR_160-185_fpocket_input.pdb"
    assert report.status == "PASS"
    assert staged.is_file()
    assert staged.resolve().is_relative_to(ctx.run_dir.resolve())
    assert receptor.read_text(encoding="utf-8") == before


def test_m2_6_chain_and_protomer_identity_are_preserved(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    make_synthetic_state(ctx)
    write_preconditions(ctx)
    output_root = ctx.run_dir / "synthetic_fpocket"
    write_fpocket_output(
        output_root,
        "EGFR_160-185",
        [
            single_pocket(chain="A", residue=745, x=0.0, y=0.0),
            single_pocket(chain="B", residue=855, resname="VAL", x=100.0, y=40.0),
        ],
    )

    report = discover_and_normalize_pockets(
        ctx,
        states=["EGFR_160-185"],
        include_reference_states=False,
        fpocket_output_root=output_root,
        synthetic_fixture=True,
    )

    candidates = read_csv(report.raw_candidate_csv)
    assert {row["egfr_chain_ids"] for row in candidates} == {"A", "B"}
    assert {row["egfr_protomer_ids"] for row in candidates} == {"A", "B"}
    families = read_csv(report.merged_candidate_csv)
    assert {row["member_protomer_ids"] for row in families} == {"A", "B"}


def test_m2_6_multi_protomer_pocket_is_marked(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    make_synthetic_state(ctx)
    write_preconditions(ctx)
    output_root = ctx.run_dir / "synthetic_fpocket"
    write_fpocket_output(
        output_root,
        "EGFR_160-185",
        [
            {
                "atoms": [
                    {"chain": "A", "resseq": 745, "resname": "LEU", "x": 0.0, "y": 0.0, "z": 0.0},
                    {"chain": "B", "resseq": 855, "resname": "VAL", "x": 1.0, "y": 1.0, "z": 0.0},
                ]
            }
        ],
    )

    report = discover_and_normalize_pockets(
        ctx,
        states=["EGFR_160-185"],
        include_reference_states=False,
        fpocket_output_root=output_root,
        synthetic_fixture=True,
    )

    candidates = read_csv(report.raw_candidate_csv)
    assert candidates[0]["multi_protomer"] == "True"
    assert candidates[0]["egfr_protomer_ids"] == "A;B"


def test_m2_6_mapping_restoration_reports_residue_name_mismatch(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    make_synthetic_state(ctx)
    write_preconditions(ctx)
    output_root = ctx.run_dir / "synthetic_fpocket"
    write_fpocket_output(output_root, "EGFR_160-185", [single_pocket(chain="A", residue=745, resname="ALA")])

    report = discover_and_normalize_pockets(
        ctx,
        states=["EGFR_160-185"],
        include_reference_states=False,
        fpocket_output_root=output_root,
        synthetic_fixture=True,
    )

    raw_rows = read_csv(report.raw_csv_paths[0])
    assert raw_rows[0]["mapping_status"] == "WARN"
    assert "residue_name_mismatch" in raw_rows[0]["warnings"]


def test_m2_6_pocket_merge_uses_centroid_or_residue_jaccard_deterministically():
    candidates = [
        _candidate("p1", "EGFR_160-185", "primary", "A", "0", "0", "0", "745"),
        _candidate("p2", "EGFR_170-200", "primary", "A", "3", "0", "0", "900"),
        _candidate("p3", "EGFR_160-185", "primary", "B", "50", "0", "0", "855"),
        _candidate("p4", "EGFR_170-200", "primary", "B", "80", "0", "0", "855"),
        _candidate("p5", "EGFR_170-200", "primary", "A", "200", "0", "0", "999"),
    ]

    families, _members = merge_pocket_families(candidates, 6.0, 0.30)
    memberships = [row["member_candidate_pocket_ids"] for row in families]

    assert memberships == ["p1;p2", "p3;p4", "p5"]
    assert families[0]["merge_threshold_centroid_angstrom"] == "6.0" or families[0]["merge_threshold_centroid_angstrom"] == 6.0
    assert families[1]["max_pairwise_residue_jaccard"] == "1.000"


def test_m2_6_3gt8_raw_is_reference_only_and_not_primary_support(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    make_synthetic_state(ctx, state_id="3GT8_raw")
    write_preconditions(ctx)
    output_root = ctx.run_dir / "synthetic_fpocket"
    write_fpocket_output(output_root, "3GT8_raw", [single_pocket()])

    report = discover_and_normalize_pockets(
        ctx,
        states=["3GT8_raw"],
        include_reference_states=False,
        fpocket_output_root=output_root,
        synthetic_fixture=True,
    )

    families = read_csv(report.merged_candidate_csv)
    assert families[0]["representative_state_role"] == "reference"
    assert families[0]["primary_state_count"] == "0"
    assert families[0]["reference_state_count"] == "1"
    assert "reference_only_family_not_promotable_in_m2_6" in families[0]["warnings"]


def test_m2_6_missing_fpocket_binary_fails_production_mode(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    make_synthetic_state(ctx)
    write_preconditions(ctx)

    report = discover_and_normalize_pockets(
        ctx,
        states=["EGFR_160-185"],
        include_reference_states=False,
        execution_mode="production",
        fpocket_binary="fpocket_definitely_missing_for_test",
    )

    assert report.status == "FAIL"
    assert any("fpocket_binary_unavailable" in item for item in report.blockers)
    assert read_json(report.raw_status_json)["fpocket_executed"] is False


def test_m2_6_missing_m2_4_and_m2_5_evidence_warns_in_parser_only(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    make_synthetic_state(ctx)
    output_root = ctx.run_dir / "synthetic_fpocket"
    write_fpocket_output(output_root, "EGFR_160-185", [single_pocket()])

    report = discover_and_normalize_pockets(
        ctx,
        states=["EGFR_160-185"],
        include_reference_states=False,
        fpocket_output_root=output_root,
        synthetic_fixture=True,
    )

    assert report.status == "PASS_WITH_WARNINGS"
    assert any("m2_4_consensus_patch_missing" in item for item in report.warnings)
    assert any("m2_5_atp_reference_missing" in item for item in report.warnings)
    assert not (ctx.run_dir / "output" / "ppi" / "ppi_consensus_patch.csv").exists()
    assert not (ctx.run_dir / "phase2_pockets" / "atp_reference" / "atp_site_reference.csv").exists()


def test_m2_6_rejects_parser_output_root_outside_fresh_or_run_dir(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    make_synthetic_state(ctx)
    write_preconditions(ctx)
    outside = tmp_path / "outside_fpocket"
    write_fpocket_output(outside, "EGFR_160-185", [single_pocket()])

    with pytest.raises(RunContextError):
        discover_and_normalize_pockets(
            ctx,
            states=["EGFR_160-185"],
            include_reference_states=False,
            fpocket_output_root=outside,
            synthetic_fixture=True,
        )


def test_m2_6_non_goals_no_gates_exports_or_other_engines(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    make_synthetic_state(ctx)
    write_preconditions(ctx)
    output_root = ctx.run_dir / "synthetic_fpocket"
    write_fpocket_output(output_root, "EGFR_160-185", [single_pocket()])

    report = discover_and_normalize_pockets(
        ctx,
        states=["EGFR_160-185"],
        include_reference_states=False,
        fpocket_output_root=output_root,
        synthetic_fixture=True,
    )
    payload = read_json(report.raw_status_json)

    assert report.status == "PASS"
    assert payload["fpocket_executed"] is False
    assert payload["p2rank_executed"] is False
    assert payload["vina_executed"] is False
    assert payload["pyrosetta_docking_or_relaxation_executed"] is False
    assert payload["m2_7_gates_executed"] is False
    assert payload["accepted_pocket_export_created"] is False
    for relative in FORBIDDEN_M2_7_OUTPUTS:
        assert not (ctx.run_dir / relative).exists()


def _candidate(candidate_id, state_id, state_role, protomer, x, y, z, residues):
    return {
        "candidate_pocket_id": candidate_id,
        "state_id": state_id,
        "state_role": state_role,
        "receptor_id": "EGFR",
        "receptor_path": "",
        "source_tool": "fpocket",
        "source_pocket_id": candidate_id,
        "egfr_chain_ids": protomer,
        "egfr_protomer_ids": protomer,
        "multi_protomer": False,
        "pocket_center_x": x,
        "pocket_center_y": y,
        "pocket_center_z": z,
        "pocket_extent_x": "1",
        "pocket_extent_y": "1",
        "pocket_extent_z": "1",
        "pocket_volume": "",
        "fpocket_score": "",
        "fpocket_druggability_score": "",
        "pocket_residue_count": 1,
        "receptor_residues": residues,
        "mapped_uniprot_residues": residues,
        "mapping_status": "PASS",
        "raw_candidate_status": "raw_ungated",
        "evidence_status": "RAW_UNGATED",
        "warnings": "",
    }


def _cli_env():
    env = os.environ.copy()
    env["PYTHONPATH"] = str(REPO_ROOT / "fresh" / "src") + os.pathsep + env.get("PYTHONPATH", "")
    return env
