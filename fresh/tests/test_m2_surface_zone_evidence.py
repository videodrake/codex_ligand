import csv
import json
import uuid
from pathlib import Path

from egfr_myo1d.cli import main
from egfr_myo1d.core.logging_utils import initialize_logs
from egfr_myo1d.core.manifest import load_configs
from egfr_myo1d.core.run_context import RunContext
from egfr_myo1d.pocket.zone_evidence import build_ppi_surface_zone_evidence


REPO_ROOT = Path(__file__).resolve().parents[2]


def make_ctx(tmp_path):
    run_id = "zone_evidence_{}".format(uuid.uuid4().hex[:12])
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


def write_csv(path, fieldnames, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def read_csv(path):
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def test_surface_zone_evidence_integrates_optional_tool_support_and_exports_m3_boxes(tmp_path):
    ctx = make_ctx(tmp_path)
    export = ctx.run_dir / "phase2_pockets" / "export_for_m3"
    accepted_fields = [
        "pocket_family_id",
        "export_tier",
        "acceptance_tier",
        "member_state_ids",
        "member_protomer_ids",
        "primary_state_count",
        "reference_state_count",
        "hard_gate_pass",
        "non_atp_pass",
        "ppi_relationship_pass",
        "lower_lateral_pass",
        "dimer_accessibility_pass",
        "primary_state_support_or_review_flag",
        "ppi_relationship_class",
        "interface_overlap_fraction",
        "family_center_x",
        "family_center_y",
        "family_center_z",
        "acceptance_reason",
        "warnings",
    ]
    write_csv(
        export / "accepted_pockets_for_m3.csv",
        accepted_fields,
        [
            {
                "pocket_family_id": "fam_primary",
                "export_tier": "primary",
                "acceptance_tier": "tier_1",
                "member_state_ids": "EGFR_160-185;EGFR_170-200",
                "member_protomer_ids": "A;B",
                "primary_state_count": "2",
                "reference_state_count": "0",
                "hard_gate_pass": "True",
                "non_atp_pass": "True",
                "ppi_relationship_pass": "True",
                "lower_lateral_pass": "True",
                "dimer_accessibility_pass": "True",
                "primary_state_support_or_review_flag": "True",
                "ppi_relationship_class": "rim_candidate",
                "interface_overlap_fraction": "0.35",
                "family_center_x": "1.0",
                "family_center_y": "2.0",
                "family_center_z": "3.0",
                "acceptance_reason": "M2 gates pass",
                "warnings": "",
            },
            {
                "pocket_family_id": "fam_atp",
                "export_tier": "primary",
                "acceptance_tier": "tier_1",
                "member_state_ids": "EGFR_160-185",
                "member_protomer_ids": "A",
                "primary_state_count": "1",
                "reference_state_count": "0",
                "hard_gate_pass": "False",
                "non_atp_pass": "False",
                "ppi_relationship_pass": "True",
                "lower_lateral_pass": "True",
                "dimer_accessibility_pass": "True",
                "primary_state_support_or_review_flag": "True",
                "ppi_relationship_class": "rim_candidate",
                "interface_overlap_fraction": "0.30",
                "family_center_x": "4.0",
                "family_center_y": "5.0",
                "family_center_z": "6.0",
                "acceptance_reason": "ATP rejected",
                "warnings": "ATP overlap",
            },
        ],
    )
    box_fields = [
        "pocket_family_id",
        "box_id",
        "state_id",
        "protomer_id",
        "box_center_x",
        "box_center_y",
        "box_center_z",
        "box_size_x",
        "box_size_y",
        "box_size_z",
        "box_source",
    ]
    write_csv(
        export / "accepted_pocket_boxes.csv",
        box_fields,
        [
            {
                "pocket_family_id": "fam_primary",
                "box_id": "box_primary_A",
                "state_id": "EGFR_160-185",
                "protomer_id": "A",
                "box_center_x": "1.0",
                "box_center_y": "2.0",
                "box_center_z": "3.0",
                "box_size_x": "16.0",
                "box_size_y": "16.0",
                "box_size_z": "16.0",
                "box_source": "accepted_pocket_boxes",
            },
            {
                "pocket_family_id": "fam_atp",
                "box_id": "box_atp",
                "state_id": "EGFR_160-185",
                "protomer_id": "A",
                "box_center_x": "4.0",
                "box_center_y": "5.0",
                "box_center_z": "6.0",
                "box_size_x": "16.0",
                "box_size_y": "16.0",
                "box_size_z": "16.0",
                "box_source": "accepted_pocket_boxes",
            },
        ],
    )
    evidence_dir = ctx.run_dir / "phase2_pockets" / "tool_evidence"
    write_csv(
        evidence_dir / "pykvfinder_cavity.csv",
        ["pocket_family_id", "volume_A3", "depth_A", "hydropathy_score"],
        [{"pocket_family_id": "fam_primary", "volume_A3": "420", "depth_A": "9.5", "hydropathy_score": "0.2"}],
    )
    write_csv(
        evidence_dir / "mini_ftmap_probe_hotspots.csv",
        ["pocket_family_id", "probe_cluster_count", "probe_hotspot_support"],
        [{"pocket_family_id": "fam_primary", "probe_cluster_count": "5", "probe_hotspot_support": "strong"}],
    )
    write_csv(
        evidence_dir / "indeep_ligandability.csv",
        ["pocket_family_id", "indeep_ligandability_score", "overlap_with_ppi_patch"],
        [{"pocket_family_id": "fam_primary", "indeep_ligandability_score": "0.82", "overlap_with_ppi_patch": "0.44"}],
    )

    result = build_ppi_surface_zone_evidence(ctx)

    assert result.status == "PASS"
    evidence_rows = {row["pocket_family_id"]: row for row in read_csv(result.zone_evidence_csv)}
    assert evidence_rows["fam_primary"]["zone_id"] == "zone_fam_primary"
    assert evidence_rows["fam_primary"]["cavity_geometry_support"] == "pyKVFinder:volume_A3=420;depth_A=9.5"
    assert evidence_rows["fam_primary"]["probe_hotspot_support"] == "mini_ftmap:strong;clusters=5"
    assert evidence_rows["fam_primary"]["surface_interface_support"] == "InDeep:ligandability=0.82;ppi_overlap=0.44"
    assert evidence_rows["fam_atp"]["final_zone_tier"] == "blocked"

    accepted_rows = read_csv(result.accepted_zones_csv)
    assert [row["pocket_family_id"] for row in accepted_rows] == ["fam_primary"]
    box_rows = read_csv(result.zone_boxes_csv)
    assert [row["pocket_family_id"] for row in box_rows] == ["fam_primary"]
    assert box_rows[0]["zone_id"] == "zone_fam_primary"
    manifest = json.loads(result.manifest_json.read_text(encoding="utf-8"))
    assert manifest["accepted_zone_count"] == 1


def test_cli_integrate_surface_zones_writes_public_zone_targets(capsys):
    run_id = "zone_cli_{}".format(uuid.uuid4().hex[:12])
    ctx = RunContext.create(run_id, "smoke_input", repo_root=REPO_ROOT)
    initialize_logs(ctx)
    export = ctx.run_dir / "phase2_pockets" / "export_for_m3"
    write_csv(
        export / "accepted_pockets_for_m3.csv",
        [
            "pocket_family_id",
            "export_tier",
            "hard_gate_pass",
            "non_atp_pass",
            "ppi_relationship_pass",
            "lower_lateral_pass",
            "dimer_accessibility_pass",
            "primary_state_support_or_review_flag",
        ],
        [
            {
                "pocket_family_id": "fam_cli",
                "export_tier": "primary",
                "hard_gate_pass": "True",
                "non_atp_pass": "True",
                "ppi_relationship_pass": "True",
                "lower_lateral_pass": "True",
                "dimer_accessibility_pass": "True",
                "primary_state_support_or_review_flag": "True",
            }
        ],
    )
    write_csv(
        export / "accepted_pocket_boxes.csv",
        ["pocket_family_id", "box_id", "state_id", "protomer_id", "box_center_x", "box_center_y", "box_center_z", "box_size_x", "box_size_y", "box_size_z"],
        [
            {
                "pocket_family_id": "fam_cli",
                "box_id": "box_cli",
                "state_id": "EGFR_160-185",
                "protomer_id": "A",
                "box_center_x": "1",
                "box_center_y": "2",
                "box_center_z": "3",
                "box_size_x": "14",
                "box_size_y": "14",
                "box_size_z": "14",
            }
        ],
    )

    code = main(["integrate-surface-zones", "--run-id", run_id])
    out = capsys.readouterr().out

    assert code == 0
    assert "integrate-surface-zones PASS" in out
    assert (ctx.run_dir / "phase2_pockets" / "tables" / "accepted_ppi_surface_zones_for_m3.csv").is_file()
    assert (ctx.run_dir / "phase2_pockets" / "export_for_m3" / "zone_composite_boxes.csv").is_file()


def test_pocket_zone_scoring_config_is_loaded_by_manifest(tmp_path):
    ctx = make_ctx(tmp_path)
    configs = load_configs(ctx)

    assert "pocket_zone_scoring" in configs
    assert configs["pocket_zone_scoring"]["schema_version"] == "pocket_zone_scoring_v1"
    assert "hard_gates" in configs["pocket_zone_scoring"]
