"""Task 9 pocket candidate prioritization writer."""

import csv
from pathlib import Path

from egfr_myo1d.analysis.pocket_candidate_prioritization import (
    CANDIDATE_REQUIRED_COLUMNS,
    FAMILY_FIELDS,
    PRIORITIZED_FIELDS,
    build_candidate_families,
    load_pocket_plan,
    ppi_patch_residue_map_from_plan,
    prioritize_candidate_rows,
    read_csv_rows,
    validate_candidate_schema,
)
from egfr_myo1d.core.logging_utils import append_job_status, append_phase_status
from egfr_myo1d.core.manifest import initialize_manifests, now_iso, write_json
from egfr_myo1d.core.run_context import RunContextError
from egfr_myo1d.io.hashing import describe_file


SCHEMA_AUDIT_FIELDS = ["check", "field", "status", "details"]
INPUT_AUDIT_FIELDS = ["candidate_pocket_id", "source_detector", "status", "notes"]
PPI_AUDIT_FIELDS = [
    "candidate_pocket_id",
    "nearest_ppi_patch_id",
    "ppi_relationship_class",
    "ppi_residue_overlap_count",
    "nearest_ppi_centroid_distance_angstrom",
    "nearest_ppi_residue_distance_angstrom",
    "status",
    "notes",
]
ATP_AUDIT_FIELDS = ["candidate_pocket_id", "atp_relationship_class", "atp_overlap_fraction", "status", "notes"]
MEMBRANE_AUDIT_FIELDS = ["candidate_pocket_id", "membrane_compatibility_class", "membrane_proximal_fraction", "status", "notes"]
DIMER_AUDIT_FIELDS = ["candidate_pocket_id", "protomer_id", "dimer_accessibility_class", "dimer_buried_fraction", "status", "notes"]
BLOCKER_FIELDS = ["candidate_pocket_id", "blocker", "details"]

NON_GOALS_CONFIRMED = [
    "no_vina_execution",
    "no_pyrosetta_execution",
    "no_lightdock_execution",
    "no_fpocket_or_p2rank_runtime_execution",
    "no_qsub_pbs_sbatch_generation",
    "no_ligand_preparation",
    "no_compound_docking",
    "no_compound_scoring",
    "no_candidate_nomination",
    "no_egfr_mutation_or_residue_renumbering",
]


def _safe_input_path(ctx, input_root, selected_path, default_name, label):
    root = input_root if input_root.is_absolute() else (ctx.repo_root / input_root)
    root = root.resolve()
    if selected_path is None:
        path = root / default_name
    else:
        candidate = Path(selected_path)
        if candidate.is_absolute():
            path = candidate
        else:
            repo_relative = (ctx.repo_root / candidate)
            path = repo_relative if repo_relative.exists() else root / candidate
    path = path.resolve()
    try:
        path.relative_to(root)
    except ValueError:
        raise RunContextError("{0} path escapes input root: {1}".format(label, path))
    return root, path


def _write_csv(path, rows, fieldnames, ctx):
    safe_path = ctx.require_within_run_dir(path)
    safe_path.parent.mkdir(parents=True, exist_ok=True)
    with safe_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _write_text(path, text, ctx):
    safe_path = ctx.require_within_run_dir(path)
    safe_path.parent.mkdir(parents=True, exist_ok=True)
    safe_path.write_text(text, encoding="utf-8")


def _input_record(ctx, path):
    record = {"path": ctx.relative_to_repo(path)}
    record.update(describe_file(path))
    return record


def _priority_counts(rows):
    counts = {}
    for row in rows:
        label = row.get("final_priority_class", "")
        counts[label] = counts.get(label, 0) + 1
    return counts


def _summary_text(report):
    lines = [
        "# Task 9 Pocket Candidate Prioritization Summary",
        "",
        "- status: {0}".format(report["status"]),
        "- candidate pockets ingested: {0}".format(report["candidate_pocket_count"]),
        "- carry-forward candidates: {0}".format(report["carry_forward_count"]),
        "- blockers: {0}".format(report["blocker_count"]),
        "- planned compound docking jobs: 0",
        "- pocket detector runtime executed: false",
        "",
        "Task 9 ingests provided detector-style pocket records and prioritizes them against Task 8 PPI-guided planning evidence.",
        "It does not run fpocket, P2Rank, Vina, PyRosetta, LightDock, compound docking, scoring, or candidate nomination.",
        "",
        "Priority counts:",
    ]
    for key in sorted(report.get("priority_counts", {}).keys()):
        lines.append("- {0}: {1}".format(key, report["priority_counts"][key]))
    if report.get("blockers"):
        lines.append("")
        lines.append("Blockers:")
        for item in report["blockers"]:
            lines.append("- {0}: {1}".format(item.get("candidate_pocket_id", "schema"), item.get("blocker", "")))
    return "\n".join(lines) + "\n"


def prioritize_pocket_candidates(ctx, mode, profile, input_root, pocket_plan=None, candidate_pockets=None, ppi_consensus=None, strict=False):
    initialize_manifests(ctx, mode)
    root, pocket_plan_path = _safe_input_path(ctx, input_root, pocket_plan, "pocket_discovery_plan.json", "pocket plan")
    _root2, candidate_path = _safe_input_path(ctx, input_root, candidate_pockets, "pocket_detector_candidates.csv", "candidate pockets")

    pockets_dir = ctx.require_within_run_dir(ctx.run_dir / "pockets")
    pockets_dir.mkdir(parents=True, exist_ok=True)
    prioritized_path = pockets_dir / "egfr_myo1d_prioritized_pocket_candidates.csv"
    family_path = pockets_dir / "egfr_myo1d_pocket_candidate_families.csv"
    json_path = pockets_dir / "pocket_candidate_prioritization.json"
    schema_audit_path = ctx.qc_dir / "pocket_candidate_schema_audit.csv"
    input_audit_path = ctx.qc_dir / "pocket_candidate_input_audit.csv"
    ppi_audit_path = ctx.qc_dir / "ppi_pocket_proximity_audit.csv"
    atp_audit_path = ctx.qc_dir / "atp_exclusion_audit.csv"
    membrane_audit_path = ctx.qc_dir / "membrane_compatibility_audit.csv"
    dimer_audit_path = ctx.qc_dir / "dimer_accessibility_audit.csv"
    blocker_path = ctx.qc_dir / "pocket_candidate_blockers.csv"
    manifest_path = ctx.manifest_dir / "pocket_candidate_prioritization_manifest.json"
    summary_path = ctx.reports_dir / "task9_pocket_candidate_prioritization_summary.md"

    blockers = []
    warnings = []
    schema_audit = []
    candidate_rows = []
    fieldnames = []
    plan = {}
    ppi_residue_map = {}
    status = "PASS"

    if not pocket_plan_path.exists():
        blockers.append({"candidate_pocket_id": "plan", "blocker": "missing_pocket_plan", "details": "Task 8 pocket_discovery_plan.json is missing"})
        status = "NO_GO"
    else:
        plan = load_pocket_plan(pocket_plan_path)
        ppi_residue_map = ppi_patch_residue_map_from_plan(plan, pocket_plan_path, ctx.repo_root)

    if not candidate_path.exists():
        blockers.append({"candidate_pocket_id": "schema", "blocker": "missing_candidate_pockets", "details": "pocket_detector_candidates.csv is missing"})
        schema_audit = validate_candidate_schema([])[1]
        status = "NO_GO"
    else:
        candidate_rows, fieldnames = read_csv_rows(candidate_path)
        missing, schema_audit = validate_candidate_schema(fieldnames)
        if missing:
            blockers.append({"candidate_pocket_id": "schema", "blocker": "BLOCKED_SCHEMA", "details": ";".join(missing)})
            status = "NO_GO"

    prioritized = []
    families = []
    input_audit = []
    ppi_audit = []
    atp_audit = []
    membrane_audit = []
    dimer_audit = []
    if status != "NO_GO":
        prioritized, input_audit, ppi_audit, atp_audit, membrane_audit, dimer_audit, row_blockers = prioritize_candidate_rows(candidate_rows, ppi_residue_map)
        blockers.extend(row_blockers)
        families = build_candidate_families(prioritized)
        if blockers:
            status = "PASS_WITH_WARNINGS"
        for row in prioritized:
            if row.get("warnings"):
                warnings.append({"candidate_pocket_id": row.get("candidate_pocket_id", ""), "message": row.get("warnings", "")})
    else:
        input_audit = [{
            "candidate_pocket_id": "",
            "source_detector": "",
            "status": "FAIL",
            "notes": ";".join([item.get("details", "") for item in blockers]),
        }]

    _write_csv(prioritized_path, prioritized, PRIORITIZED_FIELDS, ctx)
    _write_csv(family_path, families, FAMILY_FIELDS, ctx)
    _write_csv(schema_audit_path, schema_audit, SCHEMA_AUDIT_FIELDS, ctx)
    _write_csv(input_audit_path, input_audit, INPUT_AUDIT_FIELDS, ctx)
    _write_csv(ppi_audit_path, ppi_audit, PPI_AUDIT_FIELDS, ctx)
    _write_csv(atp_audit_path, atp_audit, ATP_AUDIT_FIELDS, ctx)
    _write_csv(membrane_audit_path, membrane_audit, MEMBRANE_AUDIT_FIELDS, ctx)
    _write_csv(dimer_audit_path, dimer_audit, DIMER_AUDIT_FIELDS, ctx)
    _write_csv(blocker_path, blockers, BLOCKER_FIELDS, ctx)

    outputs = {
        "prioritized_pocket_candidates": ctx.relative_to_repo(prioritized_path),
        "pocket_candidate_families": ctx.relative_to_repo(family_path),
        "pocket_candidate_prioritization": ctx.relative_to_repo(json_path),
        "pocket_candidate_schema_audit": ctx.relative_to_repo(schema_audit_path),
        "pocket_candidate_input_audit": ctx.relative_to_repo(input_audit_path),
        "ppi_pocket_proximity_audit": ctx.relative_to_repo(ppi_audit_path),
        "atp_exclusion_audit": ctx.relative_to_repo(atp_audit_path),
        "membrane_compatibility_audit": ctx.relative_to_repo(membrane_audit_path),
        "dimer_accessibility_audit": ctx.relative_to_repo(dimer_audit_path),
        "pocket_candidate_blockers": ctx.relative_to_repo(blocker_path),
        "manifest": ctx.relative_to_repo(manifest_path),
        "summary_report": ctx.relative_to_repo(summary_path),
    }
    priority_counts = _priority_counts(prioritized)
    carry_forward_count = sum(1 for row in prioritized if row.get("carry_forward") is True)
    report = {
        "run_id": ctx.run_id,
        "created_at": now_iso(),
        "phase": "prioritize-pocket-candidates",
        "status": status,
        "mode": mode,
        "profile": profile,
        "input_root": str(root),
        "pocket_plan": _input_record(ctx, pocket_plan_path),
        "candidate_pockets": _input_record(ctx, candidate_path),
        "ppi_consensus": ctx.relative_to_repo(ppi_consensus) if ppi_consensus is not None else None,
        "candidate_pocket_count": len(prioritized),
        "candidate_family_count": len(families),
        "carry_forward_count": carry_forward_count,
        "blocker_count": len(blockers),
        "priority_counts": priority_counts,
        "external_detector_outputs_ingested": bool(candidate_path.exists() and status != "NO_GO"),
        "pocket_detector_runtime_executed": False,
        "planned_compound_docking_jobs": 0,
        "compound_docking_or_scoring_executed": False,
        "candidate_nomination_executed": False,
        "outputs": outputs,
        "warnings": warnings,
        "blockers": blockers,
        "non_goals_confirmed": NON_GOALS_CONFIRMED,
    }
    manifest = dict(report)
    manifest["source_file_sha256_entries"] = {
        "pocket_plan": _input_record(ctx, pocket_plan_path),
        "candidate_pockets": _input_record(ctx, candidate_path),
    }
    write_json(json_path, report, ctx)
    write_json(manifest_path, manifest, ctx)
    _write_text(summary_path, _summary_text(report), ctx)

    phase_status = "WARN" if status != "PASS" or warnings else "PASS"
    append_phase_status(
        ctx,
        "prioritize-pocket-candidates",
        phase_status,
        "Pocket candidate prioritization completed with {0}".format(status),
        {"mode": mode, "profile": profile, "candidate_pocket_count": len(prioritized)},
    )
    append_job_status(
        ctx,
        "prioritize_pocket_candidates_local",
        phase_status,
        details={"message": "local pure-Python pocket candidate intake/prioritization only"},
    )
    return report
