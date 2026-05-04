"""PPI-surface zone evidence integration for M2 to M3 handoff.

This module does not execute external tools. It consumes M2 accepted pocket
exports plus optional detector/support tables produced by fpocket-family tools,
pyKVFinder, mdpocket, mini-FTMap, InDeep, PeSTo, MaSIF, PocketMiner, or PASSer.
"""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from egfr_myo1d.core.logging_utils import append_phase_status
from egfr_myo1d.core.manifest import now_iso, write_json
from egfr_myo1d.core.run_context import RunContext


TABLES_ROOT = "phase2_pockets/tables"
EXPORT_ROOT = "phase2_pockets/export_for_m3"
TOOL_EVIDENCE_ROOT = "phase2_pockets/tool_evidence"

ZONE_EVIDENCE_FIELDS = [
    "zone_id",
    "pocket_family_id",
    "state_support",
    "symmetry_support",
    "ppi_consensus_support",
    "ppi_hotspot_support",
    "surface_interface_support",
    "cavity_geometry_support",
    "dynamic_persistence_support",
    "probe_hotspot_support",
    "cryptic_support",
    "allosteric_support",
    "non_atp_pass",
    "lower_lateral_pass",
    "dimer_accessibility_pass",
    "primary_state_support_pass",
    "final_zone_tier",
    "decision_reason",
]

ZONE_BOX_FIELDS = [
    "zone_id",
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
    "source_box_id",
]


@dataclass
class ZoneEvidenceResult:
    status: str
    blockers: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    zone_count: int = 0
    accepted_zone_count: int = 0
    zone_evidence_csv: Path | None = None
    accepted_zones_csv: Path | None = None
    zone_boxes_csv: Path | None = None
    manifest_json: Path | None = None
    report_md: Path | None = None


def _read_csv(path: Path | None) -> tuple[list[dict[str, str]], list[str]]:
    if path is None or not path.is_file():
        return [], []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        return list(reader), list(reader.fieldnames or [])


def _write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, Any]], ctx: RunContext) -> None:
    safe = ctx.require_within_run_dir(path)
    safe.parent.mkdir(parents=True, exist_ok=True)
    with safe.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _first_existing(paths: list[Path]) -> Path | None:
    for path in paths:
        if path.is_file():
            return path
    return None


def _truthy(value: str | None) -> bool:
    if value is None:
        return False
    return str(value).strip().lower() in {"1", "true", "yes", "y", "pass", "passed"}


def _truthy_or_missing(row: dict[str, str], key: str) -> bool:
    if key not in row or row.get(key, "") == "":
        return True
    return _truthy(row.get(key))


def _zone_id(pocket_family_id: str) -> str:
    safe = "".join(ch if ch.isalnum() or ch in {"_", "-"} else "_" for ch in pocket_family_id)
    return "zone_{0}".format(safe or "unknown")


def _index_optional(path: Path, fields: list[str]) -> dict[str, dict[str, str]]:
    rows, _ = _read_csv(path)
    indexed: dict[str, dict[str, str]] = {}
    for row in rows:
        key = row.get("pocket_family_id") or row.get("zone_id", "").removeprefix("zone_")
        if key:
            indexed[key] = {field: row.get(field, "") for field in fields if field in row}
    return indexed


def _support_indexes(evidence_dir: Path) -> dict[str, dict[str, dict[str, str]]]:
    return {
        "pykvfinder": _index_optional(evidence_dir / "pykvfinder_cavity.csv", ["volume_A3", "depth_A", "hydropathy_score"]),
        "mdpocket": _index_optional(evidence_dir / "mdpocket_persistence.csv", ["persistence_fraction", "frame_count", "mean_volume_A3"]),
        "mini_ftmap": _index_optional(evidence_dir / "mini_ftmap_probe_hotspots.csv", ["probe_cluster_count", "probe_hotspot_support"]),
        "indeep": _index_optional(evidence_dir / "indeep_ligandability.csv", ["indeep_ligandability_score", "indeep_epitope_score", "overlap_with_ppi_patch"]),
        "pesto": _index_optional(evidence_dir / "pesto_interface_support.csv", ["pesto_interface_score", "overlap_with_ppi_patch"]),
        "masif": _index_optional(evidence_dir / "masif_site_support.csv", ["masif_site_score", "surface_pattern_support"]),
        "pocketminer": _index_optional(evidence_dir / "pocketminer_cryptic_support.csv", ["cryptic_score", "cryptic_support"]),
        "passer": _index_optional(evidence_dir / "passer_allosteric_support.csv", ["allosteric_score", "allosteric_support"]),
    }


def _state_support(row: dict[str, str]) -> str:
    states = row.get("member_state_ids") or row.get("state_ids") or row.get("representative_state_id") or row.get("state_id") or ""
    primary = row.get("primary_state_count", "")
    reference = row.get("reference_state_count", "")
    pieces = []
    if states:
        pieces.append("states={0}".format(states))
    if primary:
        pieces.append("primary={0}".format(primary))
    if reference:
        pieces.append("reference={0}".format(reference))
    return ";".join(pieces) if pieces else "not_available"


def _symmetry_support(row: dict[str, str]) -> str:
    protomers = row.get("member_protomer_ids") or row.get("protomer_id") or ""
    unique = sorted({item.strip() for item in protomers.replace(",", ";").split(";") if item.strip()})
    if len(unique) >= 2:
        return "multi_protomer:{0}".format(";".join(unique))
    if unique:
        return "single_protomer:{0}".format(unique[0])
    return "not_available"


def _surface_support(pid: str, supports: dict[str, dict[str, dict[str, str]]]) -> str:
    indeep = supports["indeep"].get(pid)
    if indeep:
        return "InDeep:ligandability={0};ppi_overlap={1}".format(
            indeep.get("indeep_ligandability_score", ""),
            indeep.get("overlap_with_ppi_patch", ""),
        )
    pesto = supports["pesto"].get(pid)
    if pesto:
        return "PeSTo:interface={0};ppi_overlap={1}".format(
            pesto.get("pesto_interface_score", ""),
            pesto.get("overlap_with_ppi_patch", ""),
        )
    masif = supports["masif"].get(pid)
    if masif:
        return "MaSIF:site={0};support={1}".format(
            masif.get("masif_site_score", ""),
            masif.get("surface_pattern_support", ""),
        )
    return "not_available"


def _cavity_support(pid: str, supports: dict[str, dict[str, dict[str, str]]]) -> str:
    pykv = supports["pykvfinder"].get(pid)
    if pykv:
        return "pyKVFinder:volume_A3={0};depth_A={1}".format(pykv.get("volume_A3", ""), pykv.get("depth_A", ""))
    return "fallback_fpocket_or_m2_geometry"


def _dynamic_support(pid: str, supports: dict[str, dict[str, dict[str, str]]]) -> str:
    md = supports["mdpocket"].get(pid)
    if md:
        return "mdpocket:persistence={0};frames={1}".format(md.get("persistence_fraction", ""), md.get("frame_count", ""))
    return "not_available"


def _probe_support(pid: str, supports: dict[str, dict[str, dict[str, str]]]) -> str:
    probe = supports["mini_ftmap"].get(pid)
    if probe:
        return "mini_ftmap:{0};clusters={1}".format(
            probe.get("probe_hotspot_support", ""),
            probe.get("probe_cluster_count", ""),
        )
    return "not_available"


def _cryptic_support(pid: str, supports: dict[str, dict[str, dict[str, str]]]) -> str:
    cryptic = supports["pocketminer"].get(pid)
    if cryptic:
        return "PocketMiner:{0};score={1}".format(cryptic.get("cryptic_support", ""), cryptic.get("cryptic_score", ""))
    return "not_available"


def _allosteric_support(pid: str, row: dict[str, str], supports: dict[str, dict[str, dict[str, str]]]) -> str:
    passer = supports["passer"].get(pid)
    if passer:
        return "PASSer:{0};score={1}".format(passer.get("allosteric_support", ""), passer.get("allosteric_score", ""))
    relationship = row.get("ppi_relationship_class", "")
    if "allosteric" in relationship:
        return "m2_relationship:allosteric"
    return "not_available"


def _hard_gate_failures(row: dict[str, str]) -> list[str]:
    gates = {
        "non_atp_pass": "G1_non_atp",
        "ppi_relationship_pass": "G2_ppi_related",
        "lower_lateral_pass": "G3_lower_lateral",
        "dimer_accessibility_pass": "G4_dimer_accessible",
        "primary_state_support_or_review_flag": "G6_primary_state_or_review",
    }
    failures = []
    if "hard_gate_pass" in row and row.get("hard_gate_pass", "") and not _truthy(row.get("hard_gate_pass")):
        failures.append("hard_gate_pass")
    for key, label in gates.items():
        if not _truthy_or_missing(row, key):
            failures.append(label)
    return failures


def _evidence_row(row: dict[str, str], supports: dict[str, dict[str, dict[str, str]]]) -> dict[str, str]:
    pid = row.get("pocket_family_id", "")
    failures = _hard_gate_failures(row)
    final_tier = "blocked" if failures else (row.get("export_tier") or row.get("acceptance_tier") or "accepted")
    reasons = []
    if failures:
        reasons.append("blocked:{0}".format(",".join(failures)))
    else:
        reasons.append(row.get("acceptance_reason") or "hard gates pass")
    support_labels = [
        _surface_support(pid, supports),
        _cavity_support(pid, supports),
        _dynamic_support(pid, supports),
        _probe_support(pid, supports),
        _cryptic_support(pid, supports),
        _allosteric_support(pid, row, supports),
    ]
    present_support = [item for item in support_labels if item != "not_available"]
    if present_support:
        reasons.append("support={0}".format("|".join(present_support)))
    return {
        "zone_id": _zone_id(pid),
        "pocket_family_id": pid,
        "state_support": _state_support(row),
        "symmetry_support": _symmetry_support(row),
        "ppi_consensus_support": "pass" if _truthy_or_missing(row, "ppi_relationship_pass") else "fail",
        "ppi_hotspot_support": row.get("interface_overlap_fraction") or row.get("residue_overlap_count") or "not_available",
        "surface_interface_support": _surface_support(pid, supports),
        "cavity_geometry_support": _cavity_support(pid, supports),
        "dynamic_persistence_support": _dynamic_support(pid, supports),
        "probe_hotspot_support": _probe_support(pid, supports),
        "cryptic_support": _cryptic_support(pid, supports),
        "allosteric_support": _allosteric_support(pid, row, supports),
        "non_atp_pass": str(_truthy_or_missing(row, "non_atp_pass")),
        "lower_lateral_pass": str(_truthy_or_missing(row, "lower_lateral_pass")),
        "dimer_accessibility_pass": str(_truthy_or_missing(row, "dimer_accessibility_pass")),
        "primary_state_support_pass": str(_truthy_or_missing(row, "primary_state_support_or_review_flag")),
        "final_zone_tier": final_tier,
        "decision_reason": ";".join(reasons),
    }


def _box_rows(boxes: list[dict[str, str]], accepted_zone_ids: dict[str, str]) -> list[dict[str, str]]:
    rows = []
    for row in boxes:
        pid = row.get("pocket_family_id", "")
        if pid not in accepted_zone_ids:
            continue
        box_id = row.get("box_id") or "{0}_{1}_{2}".format(pid, row.get("state_id", "state"), row.get("protomer_id", "protomer"))
        rows.append(
            {
                "zone_id": accepted_zone_ids[pid],
                "pocket_family_id": pid,
                "box_id": box_id,
                "state_id": row.get("state_id", ""),
                "protomer_id": row.get("protomer_id", ""),
                "box_center_x": row.get("box_center_x", ""),
                "box_center_y": row.get("box_center_y", ""),
                "box_center_z": row.get("box_center_z", ""),
                "box_size_x": row.get("box_size_x", ""),
                "box_size_y": row.get("box_size_y", ""),
                "box_size_z": row.get("box_size_z", ""),
                "box_source": row.get("box_source") or "accepted_pocket_boxes",
                "source_box_id": row.get("box_id", ""),
            }
        )
    return rows


def _write_report(path: Path, ctx: RunContext, result: ZoneEvidenceResult, tool_files: dict[str, bool]) -> None:
    lines = [
        "# Milestone 2 Surface Zone Filtering Summary",
        "",
        f"- run_id: `{ctx.run_id}`",
        f"- status: `{result.status}`",
        f"- zones evaluated: `{result.zone_count}`",
        f"- accepted zones for M3: `{result.accepted_zone_count}`",
        "",
        "## Optional Tool Evidence",
        "",
    ]
    for name, present in sorted(tool_files.items()):
        lines.append(f"- {name}: {'present' if present else 'missing'}")
    lines.extend(["", "## Blockers", ""])
    lines.extend([f"- {item}" for item in result.blockers] or ["- None"])
    lines.extend(["", "## Warnings", ""])
    lines.extend([f"- {item}" for item in result.warnings] or ["- None"])
    safe = ctx.require_within_run_dir(path)
    safe.parent.mkdir(parents=True, exist_ok=True)
    safe.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_ppi_surface_zone_evidence(ctx: RunContext) -> ZoneEvidenceResult:
    ctx.create_directories()
    tables_dir = ctx.require_within_run_dir(ctx.run_dir / TABLES_ROOT)
    export_dir = ctx.require_within_run_dir(ctx.run_dir / EXPORT_ROOT)
    evidence_dir = ctx.require_within_run_dir(ctx.run_dir / TOOL_EVIDENCE_ROOT)
    tables_dir.mkdir(parents=True, exist_ok=True)
    export_dir.mkdir(parents=True, exist_ok=True)

    accepted_path = _first_existing(
        [
            export_dir / "accepted_pockets_for_m3.csv",
            ctx.run_dir / "phase2_pockets" / "final" / "accepted_pockets_for_m3.csv",
            tables_dir / "accepted_pockets_for_m3.csv",
        ]
    )
    boxes_path = _first_existing(
        [
            export_dir / "accepted_pocket_boxes.csv",
            ctx.run_dir / "phase2_pockets" / "final" / "accepted_pocket_boxes.csv",
            tables_dir / "accepted_pocket_boxes.csv",
        ]
    )

    accepted_rows, _ = _read_csv(accepted_path)
    box_input_rows, _ = _read_csv(boxes_path)
    blockers = []
    warnings = []
    if accepted_path is None:
        blockers.append("accepted_pockets_for_m3.csv not found")
    if boxes_path is None:
        warnings.append("accepted_pocket_boxes.csv not found; zone_composite_boxes.csv will be empty")

    supports = _support_indexes(evidence_dir)
    tool_files = {
        "pykvfinder_cavity.csv": (evidence_dir / "pykvfinder_cavity.csv").is_file(),
        "mdpocket_persistence.csv": (evidence_dir / "mdpocket_persistence.csv").is_file(),
        "mini_ftmap_probe_hotspots.csv": (evidence_dir / "mini_ftmap_probe_hotspots.csv").is_file(),
        "indeep_ligandability.csv": (evidence_dir / "indeep_ligandability.csv").is_file(),
        "pesto_interface_support.csv": (evidence_dir / "pesto_interface_support.csv").is_file(),
        "masif_site_support.csv": (evidence_dir / "masif_site_support.csv").is_file(),
        "pocketminer_cryptic_support.csv": (evidence_dir / "pocketminer_cryptic_support.csv").is_file(),
        "passer_allosteric_support.csv": (evidence_dir / "passer_allosteric_support.csv").is_file(),
    }

    zone_rows = [_evidence_row(row, supports) for row in accepted_rows]
    accepted_zone_rows = [row for row in zone_rows if row["final_zone_tier"] != "blocked"]
    accepted_zone_ids = {row["pocket_family_id"]: row["zone_id"] for row in accepted_zone_rows}
    box_rows = _box_rows(box_input_rows, accepted_zone_ids)

    zone_evidence_csv = tables_dir / "ppi_surface_zone_evidence.csv"
    accepted_zones_csv = tables_dir / "accepted_ppi_surface_zones_for_m3.csv"
    zone_boxes_csv = export_dir / "zone_composite_boxes.csv"
    _write_csv(zone_evidence_csv, ZONE_EVIDENCE_FIELDS, zone_rows, ctx)
    _write_csv(accepted_zones_csv, ZONE_EVIDENCE_FIELDS, accepted_zone_rows, ctx)
    _write_csv(zone_boxes_csv, ZONE_BOX_FIELDS, box_rows, ctx)

    status = "FAIL" if blockers else ("WARN" if warnings else "PASS")
    manifest_json = ctx.require_within_run_dir(ctx.run_dir / "manifest" / "ppi_surface_zone_evidence_manifest.json")
    report_md = ctx.require_within_run_dir(ctx.reports_dir / "milestone2_surface_zone_filtering_summary.md")
    result = ZoneEvidenceResult(
        status=status,
        blockers=blockers,
        warnings=warnings,
        zone_count=len(zone_rows),
        accepted_zone_count=len(accepted_zone_rows),
        zone_evidence_csv=zone_evidence_csv,
        accepted_zones_csv=accepted_zones_csv,
        zone_boxes_csv=zone_boxes_csv,
        manifest_json=manifest_json,
        report_md=report_md,
    )
    write_json(
        manifest_json,
        {
            "run_id": ctx.run_id,
            "created_at": now_iso(),
            "status": status,
            "zone_count": len(zone_rows),
            "accepted_zone_count": len(accepted_zone_rows),
            "inputs": {
                "accepted_pockets_for_m3": str(accepted_path) if accepted_path else None,
                "accepted_pocket_boxes": str(boxes_path) if boxes_path else None,
                "tool_evidence_dir": str(evidence_dir),
            },
            "optional_tool_files": tool_files,
            "outputs": {
                "ppi_surface_zone_evidence": str(zone_evidence_csv),
                "accepted_ppi_surface_zones_for_m3": str(accepted_zones_csv),
                "zone_composite_boxes": str(zone_boxes_csv),
            },
            "blockers": blockers,
            "warnings": warnings,
        },
        ctx,
    )
    _write_report(report_md, ctx, result, tool_files)
    append_phase_status(
        ctx,
        phase="m2-surface-zone-evidence",
        status=status,
        message="PPI-surface zone evidence integration completed with {0}".format(status),
        details={
            "zone_count": len(zone_rows),
            "accepted_zone_count": len(accepted_zone_rows),
            "zone_evidence": ctx.relative_to_repo(zone_evidence_csv),
            "zone_boxes": ctx.relative_to_repo(zone_boxes_csv),
            "blockers": blockers,
            "warnings": warnings,
        },
    )
    return result
