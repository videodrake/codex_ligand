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
from egfr_myo1d.pocket.dimer_accessibility import DIMER_ACCESSIBILITY_FIELDS, evaluate_dimer_accessibility
from egfr_myo1d.pocket.membrane_gate import MEMBRANE_GEOMETRY_FIELDS, evaluate_membrane_geometry
from egfr_myo1d.pocket.pocket_gate_apply import (
    ACCEPTED_FIELDS,
    FORBIDDEN_EXPORT_OUTPUTS,
    GATE_QC_FIELDS,
    apply_pocket_gates,
)
from egfr_myo1d.pocket.ppi_relationship import (
    PPI_RELATIONSHIP_FIELDS,
    classify_ppi_relationship,
)


REPO_ROOT = Path(__file__).resolve().parents[2]


def unique_run_id(prefix="pytest_m2_phase7"):
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


RESIDUE_COORDS = {
    745: (40.0, 0.0, 10.0, "LYS"),
    900: (-20.0, 0.0, 10.0, "GLU"),
    901: (-14.0, 0.0, 10.0, "ASP"),
    902: (-4.0, 0.0, 10.0, "ASN"),
    930: (-1.0, 0.0, 10.0, "VAL"),
    950: (80.0, 0.0, 10.0, "PHE"),
}


def write_receptor_pdb(path):
    lines = []
    serial = 1
    for chain in ["A", "B"]:
        sign = -1.0 if chain == "A" else 1.0
        for residue, (x, y, z, resname) in RESIDUE_COORDS.items():
            x_value = x if chain == "A" else -x
            if residue == 930:
                x_value = sign * 1.0
            for atom_name, dx, element in [("N", -0.5, "N"), ("CA", 0.0, "C"), ("C", 0.5, "C")]:
                lines.append(pdb_line("ATOM", serial, atom_name, resname, chain, residue, x_value + dx, y, z, element))
                serial += 1
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(lines) + "END\n", encoding="utf-8")


def write_mapping(ctx, state_id):
    rows = []
    for protomer, chain, offset in [("A", "A", 0), ("B", "B", 1000)]:
        for residue, (_x, _y, _z, resname) in RESIDUE_COORDS.items():
            rows.append(
                {
                    "state": state_id,
                    "source_file": "{0}.pdb".format(state_id),
                    "protomer_id": protomer,
                    "source_chain": chain,
                    "source_resseq": residue,
                    "source_icode": "",
                    "source_resname": resname,
                    "runtime_chain": chain,
                    "runtime_resseq": residue + offset,
                    "atom_count": 3,
                    "role": "dockable",
                }
            )
    write_residue_mapping(ctx.qc_dir / "{0}_receptor_mapping.csv".format(state_id), rows, ctx)


def make_state(ctx, state_id):
    ctx.create_directories()
    write_mapping(ctx, state_id)
    write_receptor_pdb(
        ctx.run_dir
        / "normalized"
        / "receptors"
        / "{0}_dockable_669_1014_explicit_AB.pdb".format(state_id)
    )


def write_csv_file(path, fieldnames, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_membrane_frame(ctx):
    payload = {
        "coordinate_convention": "Z+ is membrane normal from extracellular to intracellular/cytosolic; X+ is chain A to chain B",
        "states": {},
    }
    for state in ["EGFR_160-185", "EGFR_170-200", "3GT8_raw"]:
        payload["states"][state] = {
            "role": "primary_membrane_validated_state" if state != "3GT8_raw" else "crystallographic_reference_control",
            "status": "PASS",
            "n_membrane": [0.0, 0.0, 1.0],
            "x_dimer_axis": [1.0, 0.0, 0.0],
            "protomer_a_centroid": [-10.0, 0.0, 0.0],
            "protomer_b_centroid": [10.0, 0.0, 0.0],
            "p_jm_anchor": [0.0, 0.0, 0.0],
            "lower_direction_sign": 1,
        }
    target = ctx.manifest_dir / "membrane_frame.json"
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return target


def ppi_row(state, protomer, residue, x):
    return {
        "ppi_patch_id": "ppi_{0}_{1}_{2}".format(state, protomer, residue),
        "receptor_id": "EGFR",
        "state_id": state,
        "state_role": "primary",
        "egfr_chain_id": protomer,
        "egfr_protomer_id": protomer,
        "egfr_residue_number": residue,
        "egfr_residue_name": RESIDUE_COORDS[residue][3],
        "egfr_contact_centroid_x": x,
        "egfr_contact_centroid_y": 0.0,
        "egfr_contact_centroid_z": 10.0,
    }


def family_row(family_id, members, states, roles, protomers, residues, center, primary_count, reference_count):
    return {
        "pocket_family_id": family_id,
        "representative_candidate_pocket_id": members.split(";")[0],
        "representative_state_id": states.split(";")[0],
        "representative_state_role": roles.split(";")[0],
        "representative_protomer_id": protomers.split(";")[0],
        "member_candidate_pocket_ids": members,
        "member_state_ids": states,
        "member_state_roles": roles,
        "member_protomer_ids": protomers,
        "primary_state_count": primary_count,
        "reference_state_count": reference_count,
        "raw_candidate_count": len(members.split(";")),
        "family_center_x": center[0],
        "family_center_y": center[1],
        "family_center_z": center[2],
        "centroid_spread_angstrom": "0.000",
        "residue_union_uniprot": residues,
        "residue_intersection_uniprot": residues,
        "max_pairwise_residue_jaccard": "1.000",
        "median_pairwise_residue_jaccard": "1.000",
        "merge_method": "synthetic_fixture",
        "merge_threshold_centroid_angstrom": "6.0",
        "merge_threshold_residue_jaccard": "0.30",
        "merge_status": "raw_family_ungated",
        "evidence_status": "RAW_UNGATED",
        "warnings": "",
    }


def raw_candidate(candidate_id, state, protomer, residues, mapping_status="PASS"):
    return {
        "candidate_pocket_id": candidate_id,
        "state_id": state,
        "state_role": "reference" if state == "3GT8_raw" else "primary",
        "receptor_id": "EGFR",
        "receptor_path": "",
        "source_tool": "fpocket",
        "source_pocket_id": candidate_id,
        "egfr_chain_ids": protomer,
        "egfr_protomer_ids": protomer,
        "multi_protomer": "False",
        "pocket_center_x": "",
        "pocket_center_y": "",
        "pocket_center_z": "",
        "pocket_extent_x": "",
        "pocket_extent_y": "",
        "pocket_extent_z": "",
        "pocket_volume": "",
        "fpocket_score": "",
        "fpocket_druggability_score": "",
        "pocket_residue_count": len(residues.split(";")) if residues else 0,
        "receptor_residues": residues,
        "mapped_uniprot_residues": residues,
        "mapping_status": mapping_status,
        "raw_candidate_status": "raw_ungated",
        "evidence_status": "RAW_UNGATED",
        "warnings": "",
    }


def make_complete_gate_fixture(ctx):
    for state in ["EGFR_160-185", "EGFR_170-200", "3GT8_raw"]:
        make_state(ctx, state)
    write_membrane_frame(ctx)
    ppi_fields = [
        "ppi_patch_id", "receptor_id", "state_id", "state_role", "egfr_chain_id",
        "egfr_protomer_id", "egfr_residue_number", "egfr_residue_name",
        "egfr_contact_centroid_x", "egfr_contact_centroid_y", "egfr_contact_centroid_z",
    ]
    ppi_rows = [
        ppi_row("EGFR_160-185", "A", 900, -20.0),
        ppi_row("EGFR_170-200", "A", 900, -20.0),
        ppi_row("EGFR_160-185", "A", 930, -1.0),
    ]
    write_csv_file(ctx.run_dir / "output" / "ppi" / "ppi_consensus_patch.csv", ppi_fields, ppi_rows)
    atp_fields = ["atp_site_id", "uniprot_residue_number", "residue_name"]
    write_csv_file(ctx.run_dir / "phase2_pockets" / "atp_reference" / "atp_site_reference.csv", atp_fields, [{"atp_site_id": "atp1", "uniprot_residue_number": "745", "residue_name": "LYS"}])
    centroid_fields = [
        "state_id", "receptor_id", "egfr_chain_id", "egfr_protomer_id", "centroid_method",
        "atp_centroid_x", "atp_centroid_y", "atp_centroid_z", "residue_count", "atom_count",
        "reference_radius_angstrom", "evidence_status", "warnings",
    ]
    centroid_rows = [
        {"state_id": "EGFR_160-185", "receptor_id": "EGFR", "egfr_chain_id": "A", "egfr_protomer_id": "A", "centroid_method": "synthetic", "atp_centroid_x": "40", "atp_centroid_y": "0", "atp_centroid_z": "10", "residue_count": "1", "atom_count": "3", "reference_radius_angstrom": "10", "evidence_status": "PASS", "warnings": ""},
        {"state_id": "EGFR_170-200", "receptor_id": "EGFR", "egfr_chain_id": "A", "egfr_protomer_id": "A", "centroid_method": "synthetic", "atp_centroid_x": "40", "atp_centroid_y": "0", "atp_centroid_z": "10", "residue_count": "1", "atom_count": "3", "reference_radius_angstrom": "10", "evidence_status": "PASS", "warnings": ""},
    ]
    write_csv_file(ctx.run_dir / "phase2_pockets" / "atp_reference" / "atp_site_centroid_by_state.csv", centroid_fields, centroid_rows)
    family_fields = [
        "pocket_family_id", "representative_candidate_pocket_id", "representative_state_id",
        "representative_state_role", "representative_protomer_id", "member_candidate_pocket_ids",
        "member_state_ids", "member_state_roles", "member_protomer_ids", "primary_state_count",
        "reference_state_count", "raw_candidate_count", "family_center_x", "family_center_y",
        "family_center_z", "centroid_spread_angstrom", "residue_union_uniprot",
        "residue_intersection_uniprot", "max_pairwise_residue_jaccard",
        "median_pairwise_residue_jaccard", "merge_method", "merge_threshold_centroid_angstrom",
        "merge_threshold_residue_jaccard", "merge_status", "evidence_status", "warnings",
    ]
    families = [
        family_row("fam_primary", "p1;p2", "EGFR_160-185;EGFR_170-200", "primary;primary", "A", "900", (-20, 0, 10), 2, 0),
        family_row("fam_secondary", "p3", "EGFR_160-185", "primary", "A", "901", (-14, 0, 10), 1, 0),
        family_row("fam_reference", "p4", "3GT8_raw", "reference", "A", "900", (-20, 0, 10), 0, 1),
        family_row("fam_atp_overlap", "p5", "EGFR_160-185", "primary", "A", "745", (-20, 0, 10), 1, 0),
        family_row("fam_atp_centroid", "p6", "EGFR_160-185", "primary", "A", "900", (40, 0, 10), 1, 0),
        family_row("fam_low", "p7", "EGFR_160-185", "primary", "A", "950", (-18, 0, 10), 1, 0),
        family_row("fam_membrane", "p8", "EGFR_160-185", "primary", "A", "900", (-5, 0, 10), 1, 0),
        family_row("fam_dimer", "p9", "EGFR_160-185", "primary", "A", "930", (-20, 0, 10), 1, 0),
        family_row("fam_mapping", "p10", "EGFR_160-185", "primary", "A", "999", (-20, 0, 10), 1, 0),
    ]
    write_csv_file(ctx.run_dir / "phase2_pockets" / "merged" / "pocket_candidates_merged.csv", family_fields, families)
    raw_fields = [
        "candidate_pocket_id", "state_id", "state_role", "receptor_id", "receptor_path",
        "source_tool", "source_pocket_id", "egfr_chain_ids", "egfr_protomer_ids",
        "multi_protomer", "pocket_center_x", "pocket_center_y", "pocket_center_z",
        "pocket_extent_x", "pocket_extent_y", "pocket_extent_z", "pocket_volume",
        "fpocket_score", "fpocket_druggability_score", "pocket_residue_count",
        "receptor_residues", "mapped_uniprot_residues", "mapping_status",
        "raw_candidate_status", "evidence_status", "warnings",
    ]
    raw_rows = [
        raw_candidate("p1", "EGFR_160-185", "A", "900"),
        raw_candidate("p2", "EGFR_170-200", "A", "900"),
        raw_candidate("p3", "EGFR_160-185", "A", "901"),
        raw_candidate("p4", "3GT8_raw", "A", "900"),
        raw_candidate("p5", "EGFR_160-185", "A", "745"),
        raw_candidate("p6", "EGFR_160-185", "A", "900"),
        raw_candidate("p7", "EGFR_160-185", "A", "950"),
        raw_candidate("p8", "EGFR_160-185", "A", "900"),
        raw_candidate("p9", "EGFR_160-185", "A", "930"),
        raw_candidate("p10", "EGFR_160-185", "A", "999", mapping_status="MISSING"),
    ]
    write_csv_file(ctx.run_dir / "phase2_pockets" / "merged" / "pocket_candidates_raw.csv", raw_fields, raw_rows)


def test_m2_7_cli_smoke_runs_without_science_engines():
    run_id = unique_run_id("m2_7_cli_smoke")
    ctx = RunContext.create(run_id, "smoke_input", repo_root=REPO_ROOT)
    make_complete_gate_fixture(ctx)

    proc = subprocess.run(
        [
            sys.executable,
            "-m",
            "egfr_myo1d.cli",
            "gate-m2-pockets",
            "--run-id",
            run_id,
            "--synthetic-fixture",
            "true",
        ],
        cwd=str(REPO_ROOT),
        env=_cli_env(),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    assert proc.returncode == 0
    assert b"gate-m2-pockets" in proc.stdout
    assert (ctx.run_dir / "phase2_pockets" / "gated" / "pocket_gate_qc.csv").is_file()


def test_m2_7_ppi_relationship_classifier_classes_and_not_centroid_only():
    ppi_rows = [{"state_id": "s", "egfr_protomer_id": "A", "egfr_residue_number": "900", "egfr_contact_centroid_x": "0", "egfr_contact_centroid_y": "0", "egfr_contact_centroid_z": "0"}]
    geometry = {
        ("s", "A", "900"): {"ca": [(0, 0, 0)], "heavy": [(0, 0, 0)]},
        ("s", "A", "901"): {"ca": [(6, 0, 0)], "heavy": [(6, 0, 0)]},
        ("s", "A", "902"): {"ca": [(14, 0, 0)], "heavy": [(14, 0, 0)]},
        ("s", "A", "950"): {"ca": [(90, 0, 0)], "heavy": [(90, 0, 0)]},
    }
    thresholds = {
        "ppi_rim_heavy_atom_cutoff_angstrom": 6.0,
        "ppi_rim_ca_cutoff_angstrom": 10.0,
        "ppi_allosteric_heavy_atom_cutoff_angstrom": 12.0,
        "ppi_allosteric_ca_cutoff_angstrom": 18.0,
        "ppi_allosteric_centroid_cutoff_angstrom": 18.0,
    }

    assert _classify("900", (50, 0, 0), ppi_rows, geometry, thresholds) == "orthosteric"
    assert _classify("901", (50, 0, 0), ppi_rows, geometry, thresholds) == "rim"
    assert _classify("902", (50, 0, 0), ppi_rows, geometry, thresholds, same_lower=True) == "allosteric_near"
    assert _classify("950", (1, 0, 0), ppi_rows, geometry, thresholds, same_lower=True) == "low_relevance"


def test_m2_7_ppi_overlap_requires_matching_state_and_protomer():
    thresholds = {
        "ppi_rim_heavy_atom_cutoff_angstrom": 6.0,
        "ppi_rim_ca_cutoff_angstrom": 10.0,
        "ppi_allosteric_heavy_atom_cutoff_angstrom": 12.0,
        "ppi_allosteric_ca_cutoff_angstrom": 18.0,
        "ppi_allosteric_centroid_cutoff_angstrom": 18.0,
    }
    ppi_rows = [
        {
            "state_id": "other_state",
            "egfr_protomer_id": "A",
            "egfr_residue_number": "900",
            "egfr_contact_centroid_x": "0",
            "egfr_contact_centroid_y": "0",
            "egfr_contact_centroid_z": "0",
        },
        {
            "state_id": "s",
            "egfr_protomer_id": "B",
            "egfr_residue_number": "900",
            "egfr_contact_centroid_x": "0",
            "egfr_contact_centroid_y": "0",
            "egfr_contact_centroid_z": "0",
        },
    ]
    family = {
        "pocket_family_id": "f",
        "representative_candidate_pocket_id": "p",
        "member_candidate_pocket_ids": "p",
        "member_state_ids": "s",
        "member_state_roles": "primary",
        "member_protomer_ids": "A",
        "residue_union_uniprot": "900",
        "family_center_x": 0,
        "family_center_y": 0,
        "family_center_z": 0,
    }

    row = classify_ppi_relationship(family, ppi_rows, {}, thresholds)

    assert row["ppi_relationship_class"] == "not_evaluable"
    assert row["ppi_relationship_pass"] is False
    assert row["residue_overlap_count"] == 0
    assert row["ppi_relationship_reason"] == "no_state_protomer_matched_ppi_patch_evidence"


def test_m2_7_gate_classes_cover_atp_ppi_membrane_dimer_and_mapping(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    make_complete_gate_fixture(ctx)
    report = apply_pocket_gates(ctx, synthetic_fixture=True)
    rows = {row["pocket_family_id"]: row for row in read_csv(report.gate_qc_csv)}

    assert rows["fam_primary"]["gate_class"] == "accepted_primary_pocket"
    assert rows["fam_secondary"]["gate_class"] == "accepted_secondary_cryptic"
    assert rows["fam_reference"]["gate_class"] == "reference_only"
    assert rows["fam_atp_overlap"]["gate_class"] == "atp_reject"
    assert rows["fam_atp_centroid"]["gate_class"] == "atp_reject"
    assert rows["fam_low"]["gate_class"] == "ppi_low_relevance_reject"
    assert rows["fam_membrane"]["gate_class"] == "membrane_geometry_reject"
    assert rows["fam_dimer"]["gate_class"] == "dimer_buried_reject"
    assert rows["fam_mapping"]["gate_class"] == "mapping_reject"
    assert rows["fam_atp_overlap"]["soft_score_status"] == "not_applied_hard_gate_failed"
    accepted_rows = read_csv(ctx.run_dir / "phase2_pockets" / "tables" / "accepted_pockets_for_m3.csv")
    accepted_ids = {row["pocket_family_id"] for row in accepted_rows}
    assert "fam_reference" not in accepted_ids
    assert accepted_ids == {"fam_primary", "fam_secondary"}


def test_m2_7_membrane_frame_and_lateral_score_behavior(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    make_complete_gate_fixture(ctx)
    report = apply_pocket_gates(ctx, synthetic_fixture=True)
    rows = {row["pocket_family_id"]: row for row in read_csv(report.membrane_geometry_csv)}

    assert rows["fam_primary"]["membrane_frame_status"] == "PASS"
    assert float(rows["fam_primary"]["lateral_score"]) >= 0.15
    assert rows["fam_primary"]["lower_lateral_pass"] == "True"
    assert float(rows["fam_membrane"]["lateral_score"]) < 0.15
    assert rows["fam_membrane"]["lower_lateral_pass"] == "False"
    source = (REPO_ROOT / "fresh" / "src" / "egfr_myo1d" / "pocket" / "membrane_gate.py").read_text(encoding="utf-8")
    assert "[0, 0, 1]" not in source
    assert "[0,0,1]" not in source


def test_m2_7_membrane_gate_skips_reference_frame_in_merged_family():
    membrane_frame = {
        "coordinate_convention": "Z+ is membrane normal from extracellular to intracellular/cytosolic; X+ is chain A to chain B",
        "states": {
            "3GT8_raw": {
                "role": "crystallographic_reference_control",
                "status": "reference_control",
                "n_membrane": [0.0, 0.0, 1.0],
                "p_jm_anchor": [0.0, 0.0, 0.0],
                "protomer_a_centroid": [-10.0, 0.0, 0.0],
                "protomer_b_centroid": [10.0, 0.0, 0.0],
            },
            "EGFR_160-185": {
                "role": "primary_membrane_validated_state",
                "status": "WARN",
                "n_membrane": [0.0, 0.0, 1.0],
                "p_jm_anchor": [0.0, 0.0, 0.0],
                "protomer_a_centroid": [-10.0, 0.0, 0.0],
                "protomer_b_centroid": [10.0, 0.0, 0.0],
            },
        },
    }
    family = {
        "pocket_family_id": "fam_reference_first",
        "representative_candidate_pocket_id": "p1",
        "member_state_ids": "3GT8_raw;EGFR_160-185",
        "member_state_roles": "primary;reference",
        "member_protomer_ids": "A;B",
        "family_center_x": -20.0,
        "family_center_y": 0.0,
        "family_center_z": 10.0,
    }

    row = evaluate_membrane_geometry(
        family,
        membrane_frame,
        {"EGFR_160-185": (0.0, 20.0)},
        {"lower_z_minimum": -5.0, "lateral_score_pass_threshold": 0.15},
    )

    assert row["membrane_frame_id"] == "EGFR_160-185"
    assert row["lower_lateral_pass"] is True
    assert row["membrane_geometry_reason"] == "lower_lateral_pass"
    assert "membrane_normal_missing" not in row["warnings"]
    assert "reference_state_skipped_for_membrane_frame:3GT8_raw" in row["warnings"]
    assert "multi_protomer_lateral_score_max_of_A_B" in row["warnings"]


def test_m2_7_membrane_gate_infers_lower_direction_from_receptor_centroids():
    membrane_frame = {
        "coordinate_convention": "Z+ is membrane normal from extracellular to intracellular/cytosolic; X+ is chain A to chain B",
        "states": {
            "EGFR_160-185": {
                "role": "primary_membrane_validated_state",
                "status": "PASS",
                "n_membrane": [0.0, 0.0, 1.0],
                "p_jm_anchor": [0.0, 0.0, 0.0],
                "protomer_a_centroid": [-10.0, 0.0, -50.0],
                "protomer_b_centroid": [10.0, 0.0, -50.0],
            },
        },
    }
    family = {
        "pocket_family_id": "fam_negative_z_lower",
        "representative_candidate_pocket_id": "p1",
        "member_state_ids": "EGFR_160-185",
        "member_state_roles": "primary",
        "member_protomer_ids": "A",
        "family_center_x": -20.0,
        "family_center_y": 0.0,
        "family_center_z": -45.0,
    }

    row = evaluate_membrane_geometry(
        family,
        membrane_frame,
        {"EGFR_160-185": (-80.0, 20.0)},
        {"lower_z_minimum": -5.0, "lateral_score_pass_threshold": 0.15},
    )

    assert row["pocket_membrane_z"] == "-45.000"
    assert row["lower_gate_pass"] is True
    assert row["lateral_gate_pass"] is True
    assert row["lower_lateral_pass"] is True
    assert "lower_direction_inferred_from_receptor_centroids" in row["warnings"]


def test_m2_7_dimer_accessibility_sasa_optional_warning(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    make_complete_gate_fixture(ctx)
    report = apply_pocket_gates(ctx, synthetic_fixture=True)
    rows = {row["pocket_family_id"]: row for row in read_csv(report.dimer_accessibility_csv)}

    assert rows["fam_dimer"]["dimer_accessibility_pass"] == "False"
    assert float(rows["fam_dimer"]["interface_overlap_fraction"]) >= 0.25
    assert rows["fam_primary"]["dimer_accessibility_pass"] == "True"
    assert rows["fam_primary"]["sasa_status"] == "not_computed_m2_7_optional"


def test_m2_7_dimer_accessibility_fails_when_interface_index_unavailable():
    family = {
        "pocket_family_id": "fam",
        "representative_candidate_pocket_id": "p",
        "member_state_ids": "EGFR_160-185",
        "member_state_roles": "primary",
        "member_protomer_ids": "A",
        "residue_union_uniprot": "900;901",
    }
    row = evaluate_dimer_accessibility(family, set(), threshold=0.25)

    assert row["dimer_accessibility_pass"] is False
    assert row["dimer_accessibility_reason"] == "interface_residue_index_unavailable"
    assert "interface_residue_index_unavailable" in row["warnings"]


def test_m2_7_schema_regression_for_gate_outputs(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    make_complete_gate_fixture(ctx)
    report = apply_pocket_gates(ctx, synthetic_fixture=True)

    assert read_header(report.gate_qc_csv) == GATE_QC_FIELDS
    assert read_header(report.accepted_csv) == ACCEPTED_FIELDS
    assert read_header(report.rejected_csv) == ACCEPTED_FIELDS
    assert read_header(report.ppi_relationship_csv) == PPI_RELATIONSHIP_FIELDS
    assert read_header(report.membrane_geometry_csv) == MEMBRANE_GEOMETRY_FIELDS
    assert read_header(report.dimer_accessibility_csv) == DIMER_ACCESSIBILITY_FIELDS
    assert (ctx.run_dir / "phase2_pockets" / "tables" / "accepted_pockets_for_m3.csv").is_file()


def test_m2_7_chain_protomer_and_reference_support_are_preserved(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    make_complete_gate_fixture(ctx)
    report = apply_pocket_gates(ctx, synthetic_fixture=True)
    rows = {row["pocket_family_id"]: row for row in read_csv(report.gate_qc_csv)}

    assert rows["fam_primary"]["member_protomer_ids"] == "A"
    assert rows["fam_primary"]["primary_state_count"] == "2"
    assert rows["fam_reference"]["reference_state_count"] == "1"
    assert rows["fam_reference"]["primary_state_count"] == "0"
    assert rows["fam_reference"]["gate_class"] == "reference_only"


def test_m2_7_missing_inputs_fail_production_but_warn_in_synthetic(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    report = apply_pocket_gates(ctx, mode="production", synthetic_fixture=False)
    assert report.status == "FAIL"
    assert any("missing_m2_4_ppi_consensus_patch" in item for item in report.blockers)
    assert any("missing_m2_5_atp_reference" in item for item in report.blockers)
    assert any("missing_m2_6_pocket_candidates_merged" in item for item in report.blockers)

    ctx2 = make_tmp_run_context(tmp_path)
    report2 = apply_pocket_gates(ctx2, mode="smoke_input", synthetic_fixture=True)
    assert report2.status == "PASS_WITH_WARNINGS"
    assert any("missing_m2_4_ppi_consensus_patch" in item for item in report2.warnings)


def test_m2_7_missing_atp_centroid_fails_production(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    make_complete_gate_fixture(ctx)
    centroid_path = ctx.run_dir / "phase2_pockets" / "atp_reference" / "atp_site_centroid_by_state.csv"
    centroid_path.unlink()

    report = apply_pocket_gates(ctx, mode="production", synthetic_fixture=False)

    assert report.status == "FAIL"
    assert any("missing_m2_5_atp_centroid_by_state" in item for item in report.blockers)


def test_m2_7_rejects_input_path_outside_fresh_or_run_dir(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    outside = tmp_path / "outside_ppi.csv"
    outside.write_text("x\n", encoding="utf-8")
    with pytest.raises(RunContextError):
        apply_pocket_gates(ctx, ppi_consensus=outside, synthetic_fixture=True)


def test_m2_7_non_goals_no_export_package_or_engine_calls(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    make_complete_gate_fixture(ctx)
    report = apply_pocket_gates(ctx, synthetic_fixture=True)
    payload = read_json(report.status_json)

    assert payload["vina_executed"] is False
    assert payload["ligand_preparation_executed"] is False
    assert payload["pyrosetta_docking_or_relaxation_executed"] is False
    assert payload["lightdock_executed"] is False
    assert payload["p2rank_executed"] is False
    assert payload["compound_scoring_executed"] is False
    assert payload["candidate_nomination_executed"] is False
    assert payload["m2_8_export_package_created"] is False
    for relative in FORBIDDEN_EXPORT_OUTPUTS:
        assert not (ctx.run_dir / relative).exists()


def _classify(residues, center, ppi_rows, geometry, thresholds, same_lower=False):
    family = {
        "pocket_family_id": "f",
        "representative_candidate_pocket_id": "p",
        "member_candidate_pocket_ids": "p",
        "member_state_ids": "s",
        "member_state_roles": "primary",
        "member_protomer_ids": "A",
        "residue_union_uniprot": residues,
        "family_center_x": center[0],
        "family_center_y": center[1],
        "family_center_z": center[2],
    }
    row = classify_ppi_relationship(
        family,
        ppi_rows,
        geometry,
        thresholds,
        same_lower_lateral_neighborhood_flag=same_lower,
    )
    return row["ppi_relationship_class"]


def _cli_env():
    env = os.environ.copy()
    env["PYTHONPATH"] = str(REPO_ROOT / "fresh" / "src") + os.pathsep + env.get("PYTHONPATH", "")
    return env
