"""M2.4 EGFR-side PPI consensus patch builder.

The builder consumes restored M2.3 contact evidence and creates deterministic,
chain-resolved EGFR-side consensus patch rows. It does not run docking, score
poses, discover pockets, or nominate compounds.
"""

from __future__ import annotations

import csv
import json
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from egfr_myo1d.core.logging_utils import append_job_status, append_phase_status
from egfr_myo1d.core.manifest import now_iso, write_json
from egfr_myo1d.core.run_context import RunContext, RunContextError, ensure_within
from egfr_myo1d.ppi.restore_residue_mapping import load_runtime_residue_map


M2_4_PHASE = "m2.4-symmetry-aware-ppi-consensus-patch"
DEFAULT_CONTACT_TABLE = "phase1_ppi/tables/ppi_pose_contacts.csv"
DEFAULT_JOB_MANIFEST = "phase1_ppi/pyrosetta_adapter/pyrosetta_job_manifest.jsonl"
OUTPUT_PATCH_TABLE = "output/ppi/ppi_consensus_patch.csv"
M2_TABLE_PATCH_TABLE = "phase1_ppi/tables/ppi_consensus_patch.csv"
CONSENSUS_FIELDS = [
    "ppi_patch_id",
    "receptor_id",
    "state_id",
    "state_role",
    "egfr_chain_id",
    "egfr_protomer_id",
    "egfr_residue_number",
    "egfr_residue_name",
    "egfr_contact_count",
    "seed_support_count",
    "pose_support_count",
    "support_fraction",
    "egfr_contact_centroid_x",
    "egfr_contact_centroid_y",
    "egfr_contact_centroid_z",
    "myo1d_contact_residues",
    "myo1d_sheet_support",
    "myo1d_active_face_score",
    "sheet8_support_score",
    "sheet9_support_score",
    "sheet12_support_score",
    "membrane_side_class",
    "evidence_status",
    "warnings",
]
RESIDUE_FIELDS = [
    "state_id",
    "state_role",
    "egfr_chain_id",
    "egfr_protomer_id",
    "egfr_residue_number",
    "egfr_residue_name",
    "contact_count",
    "pose_support_count",
    "seed_support_count",
    "myo1d_contact_residues",
    "evidence_status",
    "warnings",
]
STATE_SUPPORT_FIELDS = [
    "state_id",
    "state_role",
    "patch_count",
    "unique_pose_count",
    "unique_seed_count",
    "status",
    "notes",
]
SYMMETRY_SUPPORT_FIELDS = [
    "state_id",
    "egfr_residue_number",
    "protomer_a_supported",
    "protomer_b_supported",
    "symmetry_support_count",
    "notes",
]
QC_FIELDS = ["check", "status", "details"]
NON_GOALS_CONFIRMED = [
    "no_pyrosetta_import",
    "no_pyrosetta_docking_or_relaxation",
    "no_lightdock_execution",
    "no_vina_execution",
    "no_fpocket_or_p2rank_runtime",
    "no_compound_docking_or_scoring",
    "no_pocket_discovery",
    "no_candidate_nomination",
    "no_scheduler_submission",
    "no_cleanup_deletion",
]

SHEET8 = set(range(961, 965))
SHEET9 = set(range(968, 973))
SHEET12 = set(range(993, 998))
ACTIVE_FACE = SHEET8 | SHEET9
ACCEPTED_CLASSES = {
    "accepted_active_8_9_supported_12",
    "accepted_active_8_9_only",
    "key_8_9_12_engaged",
    "partial_key_engaged",
}
REJECT_PREFIXES = ("reject", "rejected")


@dataclass
class M2ConsensusPatchReport:
    run_id: str
    mode: str
    profile: str
    status: str
    contact_count: int = 0
    patch_count: int = 0
    accepted_pose_count: int = 0
    pending_pose_qc_count: int = 0
    warnings: list[str] = field(default_factory=list)
    blockers: list[str] = field(default_factory=list)
    manifest_path: Path | None = None
    consensus_patch_csv: Path | None = None
    m2_table_patch_csv: Path | None = None
    residue_csv: Path | None = None
    state_support_csv: Path | None = None
    symmetry_support_csv: Path | None = None
    qc_csv_path: Path | None = None
    summary_report_path: Path | None = None


@dataclass
class SupportDenominators:
    pose_by_state_protomer: dict[tuple[str, str], set[str]] = field(default_factory=dict)
    seed_by_state_protomer: dict[tuple[str, str], set[str]] = field(default_factory=dict)
    pose_by_state: dict[str, set[str]] = field(default_factory=dict)
    seed_by_state: dict[str, set[str]] = field(default_factory=dict)


def _run_relative_path(ctx: RunContext, text: str | Path | None, default_text: str) -> Path:
    raw = Path(text or default_text)
    if raw.is_absolute():
        path = raw.resolve()
    elif raw.parts and raw.parts[0] == "fresh":
        path = (ctx.repo_root / raw).resolve()
    else:
        path = (ctx.run_dir / raw).resolve()
    return ensure_within(path, ctx.run_dir)


def _repo_path(ctx: RunContext, path: Path) -> str:
    return ctx.relative_to_repo(path)


def _read_csv(path: Path) -> list[dict[str, Any]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _read_jsonl(path: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                rows.append(json.loads(line))
    return rows


def _write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str], ctx: RunContext) -> None:
    safe = ctx.require_within_run_dir(path)
    safe.parent.mkdir(parents=True, exist_ok=True)
    with safe.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _state_roles(ctx: RunContext) -> dict[str, str]:
    path = ctx.manifest_dir / "m2_1_ppi_input_manifest.json"
    if not path.is_file():
        return {}
    payload = json.loads(path.read_text(encoding="utf-8"))
    return {
        str(pack.get("state_id", "")): str(pack.get("role", ""))
        for pack in payload.get("packs", [])
    }


def _mapping_paths(ctx: RunContext, job_manifest_path: Path) -> tuple[dict[str, Path], list[dict[str, Any]], list[str]]:
    warnings: list[str] = []
    mapping_by_state: dict[str, Path] = {}
    if not job_manifest_path.is_file():
        return mapping_by_state, [], ["missing_m2_2_job_manifest:{0}".format(job_manifest_path)]
    job_rows = _read_jsonl(job_manifest_path)
    for row in job_rows:
        state_id = str(row.get("state_id", ""))
        mapping_text = str(row.get("mapping_csv", ""))
        if not state_id or not mapping_text:
            warnings.append("job_manifest_row_missing_state_or_mapping")
            continue
        path = Path(mapping_text)
        if path.is_absolute():
            mapping_path = path.resolve()
        else:
            mapping_path = (ctx.repo_root / path).resolve()
        try:
            ensure_within(mapping_path, ctx.run_dir)
        except RunContextError:
            warnings.append("mapping_path_escapes_run_dir:{0}".format(mapping_path))
            continue
        if not mapping_path.is_file():
            warnings.append("mapping_path_missing:{0}".format(mapping_path))
            continue
        mapping_by_state[state_id] = mapping_path
    return mapping_by_state, job_rows, warnings


def _validated_contacts(
    ctx: RunContext,
    contacts: list[dict[str, Any]],
    mapping_by_state: dict[str, Path],
    state_roles: dict[str, str],
) -> tuple[list[dict[str, Any]], list[str], list[str]]:
    warnings: list[str] = []
    blockers: list[str] = []
    validated: list[dict[str, Any]] = []
    mapping_cache: dict[str, Any] = {}
    for index, row in enumerate(contacts, start=1):
        state_id = str(row.get("state_id", ""))
        protomer_id = str(row.get("protomer_id", ""))
        runtime_residue = str(row.get("egfr_runtime_residue", ""))
        if not state_id or not protomer_id or not runtime_residue:
            blockers.append("contact_missing_state_protomer_or_runtime_residue:row{0}".format(index))
            continue
        mapping_path = mapping_by_state.get(state_id)
        if mapping_path is None:
            blockers.append("contact_missing_mapping_for_state:{0}:row{1}".format(state_id, index))
            continue
        cache_key = str(mapping_path)
        if cache_key not in mapping_cache:
            mapping_cache[cache_key] = load_runtime_residue_map(mapping_path)
        record, reason = mapping_cache[cache_key].restore(protomer_id, runtime_residue)
        if record is None:
            blockers.append("contact_unmapped:{0}:row{1}:{2}".format(state_id, index, reason))
            continue
        original = str(row.get("egfr_original_residue", ""))
        if original and original != record.original_residue:
            blockers.append(
                "residue_number_drift:{0}:{1}:{2}:row{3}".format(
                    state_id, original, record.original_residue, index
                )
            )
            continue
        expected_protomer = record.protomer_id
        if expected_protomer and protomer_id != expected_protomer:
            blockers.append(
                "protomer_id_drift:{0}:{1}:{2}:row{3}".format(
                    state_id, protomer_id, expected_protomer, index
                )
            )
            continue
        expected_chain = record.source_chain or record.runtime_chain
        observed_chain = str(row.get("egfr_chain_id", ""))
        if observed_chain and expected_chain and observed_chain != expected_chain:
            blockers.append(
                "chain_id_drift:{0}:{1}:{2}:row{3}".format(
                    state_id, observed_chain, expected_chain, index
                )
            )
            continue
        chain_id = observed_chain or expected_chain
        if not observed_chain:
            warnings.append("chain_id_restored_from_mapping:{0}:row{1}".format(state_id, index))
        enriched = dict(row)
        enriched["state_role"] = str(row.get("state_role") or state_roles.get(state_id, ""))
        enriched["protomer_id"] = expected_protomer or protomer_id
        enriched["egfr_chain_id"] = chain_id
        enriched["egfr_original_residue"] = record.original_residue
        enriched["egfr_uniprot_residue"] = str(row.get("egfr_uniprot_residue") or record.original_residue)
        enriched["egfr_resname"] = str(row.get("egfr_resname") or record.source_resname)
        validated.append(enriched)
    return validated, warnings, blockers


def _int_or_none(value: Any) -> int | None:
    try:
        return int(str(value).strip())
    except (TypeError, ValueError):
        return None


def _float_or_none(value: Any) -> float | None:
    try:
        text = str(value).strip()
        if not text:
            return None
        return float(text)
    except (TypeError, ValueError):
        return None


def _is_rejected(row: dict[str, Any]) -> bool:
    text = str(row.get("pose_acceptance_class") or row.get("evidence_status") or "").lower()
    accepted_flag = str(row.get("accepted_pose_flag", "")).strip().lower()
    if accepted_flag in {"0", "false", "no"}:
        return True
    return text.startswith(REJECT_PREFIXES)


def _is_accepted(row: dict[str, Any]) -> bool:
    text = str(row.get("pose_acceptance_class") or row.get("evidence_status") or "").lower()
    accepted_flag = str(row.get("accepted_pose_flag", "")).strip().lower()
    if accepted_flag in {"1", "true", "yes"}:
        return True
    if text in ACCEPTED_CLASSES:
        return True
    return text.startswith("accepted")


def _eligible_contacts(rows: list[dict[str, Any]]) -> tuple[list[dict[str, Any]], int, int]:
    eligible: list[dict[str, Any]] = []
    accepted_pose_keys: set[tuple[str, str, str]] = set()
    pending_pose_keys: set[tuple[str, str, str]] = set()
    for row in rows:
        if _is_rejected(row):
            continue
        key = (
            str(row.get("state_id", "")),
            str(row.get("seed", "")),
            str(row.get("pose_id", "")),
        )
        if _is_accepted(row):
            accepted_pose_keys.add(key)
            eligible.append(row)
        else:
            pending_pose_keys.add(key)
    return eligible, len(accepted_pose_keys), len(pending_pose_keys)


def _protomer_ids_by_state(mapping_by_state: dict[str, Path], warnings: list[str]) -> dict[str, list[str]]:
    protomers_by_state: dict[str, list[str]] = {}
    for state_id, mapping_path in mapping_by_state.items():
        residue_map = load_runtime_residue_map(mapping_path)
        protomers = sorted(
            {
                protomer_id
                for protomer_id, _runtime_residue in residue_map.by_protomer_runtime
                if protomer_id
            }
        )
        if protomers:
            protomers_by_state[state_id] = protomers
        else:
            warnings.append("mapping_has_no_protomer_ids:{0}".format(state_id))
    return protomers_by_state


def _support_denominators(
    job_rows: list[dict[str, Any]],
    mapping_by_state: dict[str, Path],
    accepted_rows: list[dict[str, Any]],
    warnings: list[str],
) -> SupportDenominators:
    denominators = SupportDenominators(
        pose_by_state_protomer=defaultdict(set),
        seed_by_state_protomer=defaultdict(set),
        pose_by_state=defaultdict(set),
        seed_by_state=defaultdict(set),
    )
    protomers_by_state = _protomer_ids_by_state(mapping_by_state, warnings)
    for job in job_rows:
        state_id = str(job.get("state_id", ""))
        seed = str(job.get("seed", ""))
        if not state_id or seed == "":
            warnings.append("job_manifest_row_missing_state_or_seed")
            continue
        protomers = protomers_by_state.get(state_id, [])
        if not protomers:
            protomers = sorted(
                {
                    str(row.get("protomer_id", ""))
                    for row in accepted_rows
                    if str(row.get("state_id", "")) == state_id and row.get("protomer_id")
                }
            )
        if not protomers:
            continue
        models_per_seed = _int_or_none(job.get("models_per_seed"))
        if models_per_seed is None or models_per_seed <= 0:
            warnings.append("job_manifest_models_per_seed_missing_or_invalid:{0}:{1}".format(state_id, seed))
            models_per_seed = 1
        for protomer_id in protomers:
            denominators.seed_by_state_protomer[(state_id, protomer_id)].add(seed)
            denominators.seed_by_state[state_id].add(seed)
            for model_index in range(1, models_per_seed + 1):
                pose_key = "{0}:planned_model_{1:05d}".format(seed, model_index)
                denominators.pose_by_state_protomer[(state_id, protomer_id)].add(pose_key)
                denominators.pose_by_state[state_id].add(pose_key)

    fallback_keys: set[tuple[str, str]] = set()
    for row in accepted_rows:
        state_id = str(row.get("state_id", ""))
        protomer_id = str(row.get("protomer_id", ""))
        seed = str(row.get("seed", ""))
        pose_key = "{0}:{1}".format(seed, row.get("pose_id", ""))
        key = (state_id, protomer_id)
        if not denominators.pose_by_state_protomer[key]:
            fallback_keys.add(key)
            denominators.pose_by_state_protomer[key].add(pose_key)
            denominators.pose_by_state[state_id].add(pose_key)
        if not denominators.seed_by_state_protomer[key]:
            denominators.seed_by_state_protomer[key].add(seed)
            denominators.seed_by_state[state_id].add(seed)
    for state_id, protomer_id in sorted(fallback_keys):
        warnings.append(
            "support_denominator_fell_back_to_contact_evidence:{0}:{1}".format(
                state_id, protomer_id
            )
        )
    return denominators


def _sheet_labels(myo_residues: set[int]) -> str:
    labels = []
    if myo_residues & SHEET8:
        labels.append("sheet8")
    if myo_residues & SHEET9:
        labels.append("sheet9")
    if myo_residues & SHEET12:
        labels.append("sheet12")
    return ";".join(labels)


def _mean(values: list[float]) -> str:
    if not values:
        return ""
    return "{0:.3f}".format(sum(values) / float(len(values)))


def _build_patch_rows(
    rows: list[dict[str, Any]],
    warnings: list[str],
    denominators: SupportDenominators,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    groups: dict[tuple[str, str, str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        state_id = str(row.get("state_id", ""))
        protomer_id = str(row.get("protomer_id", ""))
        original = str(row.get("egfr_original_residue", ""))
        chain_id = str(row.get("egfr_chain_id", ""))
        groups[(state_id, chain_id, protomer_id, original)].append(row)

    symmetry_index: dict[tuple[str, str], set[str]] = defaultdict(set)
    for state_id, _chain_id, protomer_id, original in groups:
        symmetry_index[(state_id, original)].add(protomer_id)

    patch_rows: list[dict[str, Any]] = []
    residue_rows: list[dict[str, Any]] = []
    for idx, key in enumerate(sorted(groups), start=1):
        state_id, chain_id, protomer_id, original = key
        group = groups[key]
        seeds = sorted({str(row.get("seed", "")) for row in group})
        poses = sorted({"{0}:{1}".format(row.get("seed", ""), row.get("pose_id", "")) for row in group})
        methods = sorted({str(row.get("method", "")) for row in group if row.get("method")})
        myo_residues = {
            value
            for value in (_int_or_none(row.get("myo1d_residue")) for row in group)
            if value is not None
        }
        centroids_x = [_float_or_none(row.get("egfr_contact_centroid_x")) for row in group]
        centroids_y = [_float_or_none(row.get("egfr_contact_centroid_y")) for row in group]
        centroids_z = [_float_or_none(row.get("egfr_contact_centroid_z")) for row in group]
        centroid_missing = not all(value is not None for value in centroids_x + centroids_y + centroids_z)
        row_warnings = []
        if centroid_missing:
            row_warnings.append("contact_centroid_not_available")
        if warnings:
            row_warnings.extend(sorted(set(warnings)))
        support_denom = max(
            1,
            len(denominators.pose_by_state_protomer.get((state_id, protomer_id), set())),
        )
        pose_support_count = len(poses)
        support_fraction = float(pose_support_count) / float(support_denom)
        symmetry_supported = symmetry_index[(state_id, original)]
        evidence_status = "accepted_contact_evidence"
        patch_row = {
            "ppi_patch_id": "m2_4_ppi_patch_{0:03d}".format(idx),
            "receptor_id": "EGFR_MYO1D_membrane_compatible_inactive_dimer",
            "state_id": state_id,
            "state_role": str(group[0].get("state_role", "")),
            "egfr_chain_id": chain_id,
            "egfr_protomer_id": protomer_id,
            "egfr_residue_number": original,
            "egfr_residue_name": str(group[0].get("egfr_resname", "")),
            "egfr_contact_count": len(group),
            "seed_support_count": len(seeds),
            "pose_support_count": pose_support_count,
            "support_fraction": "{0:.3f}".format(support_fraction),
            "egfr_contact_centroid_x": _mean([v for v in centroids_x if v is not None]),
            "egfr_contact_centroid_y": _mean([v for v in centroids_y if v is not None]),
            "egfr_contact_centroid_z": _mean([v for v in centroids_z if v is not None]),
            "myo1d_contact_residues": ";".join(str(value) for value in sorted(myo_residues)),
            "myo1d_sheet_support": _sheet_labels(myo_residues),
            "myo1d_active_face_score": len(myo_residues & ACTIVE_FACE),
            "sheet8_support_score": len(myo_residues & SHEET8),
            "sheet9_support_score": len(myo_residues & SHEET9),
            "sheet12_support_score": len(myo_residues & SHEET12),
            "membrane_side_class": ";".join(sorted({str(row.get("membrane_side_class", "")) for row in group if row.get("membrane_side_class")})) or "unknown",
            "evidence_status": evidence_status,
            "warnings": ";".join(sorted(set(row_warnings))),
        }
        patch_rows.append(patch_row)
        residue_rows.append(
            {
                "state_id": state_id,
                "state_role": patch_row["state_role"],
                "egfr_chain_id": chain_id,
                "egfr_protomer_id": protomer_id,
                "egfr_residue_number": original,
                "egfr_residue_name": patch_row["egfr_residue_name"],
                "contact_count": len(group),
                "pose_support_count": pose_support_count,
                "seed_support_count": len(seeds),
                "myo1d_contact_residues": patch_row["myo1d_contact_residues"],
                "evidence_status": evidence_status,
                "warnings": patch_row["warnings"],
            }
        )

    state_rows = []
    for state_id in sorted({row["state_id"] for row in patch_rows}):
        state_patches = [row for row in patch_rows if row["state_id"] == state_id]
        state_rows.append(
            {
                "state_id": state_id,
                "state_role": state_patches[0]["state_role"] if state_patches else "",
                "patch_count": len(state_patches),
                "unique_pose_count": len(denominators.pose_by_state.get(state_id, set())),
                "unique_seed_count": len(denominators.seed_by_state.get(state_id, set())),
                "status": "PASS",
                "notes": "state_support_separate_from_symmetry_support",
            }
        )

    symmetry_rows = []
    for (state_id, original), protomers in sorted(symmetry_index.items()):
        symmetry_rows.append(
            {
                "state_id": state_id,
                "egfr_residue_number": original,
                "protomer_a_supported": "A" in protomers,
                "protomer_b_supported": "B" in protomers,
                "symmetry_support_count": 1 if {"A", "B"}.issubset(protomers) else 0,
                "notes": "symmetry_not_counted_as_independent_state_support",
            }
        )
    return patch_rows, residue_rows, state_rows, symmetry_rows


def _write_summary(ctx: RunContext, report: M2ConsensusPatchReport) -> Path:
    target = ctx.require_within_run_dir(ctx.reports_dir / "m2_4_ppi_consensus_patch.md")
    lines = [
        "# M2.4 Symmetry-aware PPI Consensus Patch",
        "",
        "Status: {0}".format(report.status),
        "Input contacts: {0}".format(report.contact_count),
        "Consensus patch rows: {0}".format(report.patch_count),
        "Accepted pose count: {0}".format(report.accepted_pose_count),
        "Pending pose-QC pose count: {0}".format(report.pending_pose_qc_count),
        "",
        "No single best pose is selected. Chain/protomer identity remains explicit.",
        "",
        "## Non-goals Confirmed",
    ]
    for item in NON_GOALS_CONFIRMED:
        lines.append("- {0}".format(item))
    if report.warnings:
        lines.extend(["", "## Warnings"])
        for warning in report.warnings:
            lines.append("- {0}".format(warning))
    if report.blockers:
        lines.extend(["", "## Blockers"])
        for blocker in report.blockers:
            lines.append("- {0}".format(blocker))
    target.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return target


def _status(blockers: list[str], warnings: list[str], patch_count: int) -> str:
    if blockers:
        return "FAIL"
    if warnings or patch_count == 0:
        return "PASS_WITH_WARNINGS"
    return "PASS"


def build_m2_4_consensus_patch(
    ctx: RunContext,
    mode: str = "smoke_input",
    profile: str = "codex_dev",
    contact_table: Path | None = None,
    job_manifest: Path | None = None,
) -> M2ConsensusPatchReport:
    """Build deterministic M2.4 consensus patch tables from M2.3 contacts."""

    ctx.create_directories()
    contact_path = _run_relative_path(ctx, contact_table, DEFAULT_CONTACT_TABLE)
    job_manifest_path = _run_relative_path(ctx, job_manifest, DEFAULT_JOB_MANIFEST)
    warnings: list[str] = []
    blockers: list[str] = []
    contacts: list[dict[str, Any]] = []

    if not contact_path.is_file():
        blockers.append("missing_m2_3_pose_contacts:{0}".format(contact_path))
    else:
        contacts = _read_csv(contact_path)
        if not contacts:
            warnings.append("no_restored_contact_evidence_available")

    mapping_by_state, job_rows, mapping_warnings = _mapping_paths(ctx, job_manifest_path)
    warnings.extend(mapping_warnings)
    state_roles = _state_roles(ctx)
    validated: list[dict[str, Any]] = []
    if contacts:
        validated, validation_warnings, validation_blockers = _validated_contacts(
            ctx, contacts, mapping_by_state, state_roles
        )
        warnings.extend(validation_warnings)
        blockers.extend(validation_blockers)

    eligible, accepted_pose_count, pending_pose_qc_count = _eligible_contacts(validated)
    if contacts and not eligible and not blockers:
        warnings.append("no_accepted_contact_evidence_after_pose_qc_filter")
    if pending_pose_qc_count:
        warnings.append("contacts_lack_pose_qc_acceptance_class")
    patch_rows: list[dict[str, Any]] = []
    residue_rows: list[dict[str, Any]] = []
    state_rows: list[dict[str, Any]] = []
    symmetry_rows: list[dict[str, Any]] = []
    if not blockers:
        denominators = _support_denominators(job_rows, mapping_by_state, eligible, warnings)
        patch_rows, residue_rows, state_rows, symmetry_rows = _build_patch_rows(
            eligible, warnings, denominators
        )

    output_patch = ctx.require_within_run_dir(ctx.run_dir / OUTPUT_PATCH_TABLE)
    m2_table_patch = ctx.require_within_run_dir(ctx.run_dir / M2_TABLE_PATCH_TABLE)
    residue_csv = ctx.require_within_run_dir(ctx.run_dir / "phase1_ppi/tables/ppi_consensus_residues.csv")
    state_support_csv = ctx.require_within_run_dir(ctx.run_dir / "phase1_ppi/tables/ppi_patch_state_support.csv")
    symmetry_support_csv = ctx.require_within_run_dir(ctx.run_dir / "phase1_ppi/tables/ppi_patch_symmetry_support.csv")
    qc_csv = ctx.require_within_run_dir(ctx.qc_dir / "m2_4_ppi_consensus_patch_qc.csv")

    _write_csv(output_patch, patch_rows, CONSENSUS_FIELDS, ctx)
    _write_csv(m2_table_patch, patch_rows, CONSENSUS_FIELDS, ctx)
    _write_csv(residue_csv, residue_rows, RESIDUE_FIELDS, ctx)
    _write_csv(state_support_csv, state_rows, STATE_SUPPORT_FIELDS, ctx)
    _write_csv(symmetry_support_csv, symmetry_rows, SYMMETRY_SUPPORT_FIELDS, ctx)
    qc_rows = [
        {"check": "contact_table_exists", "status": "PASS" if contact_path.is_file() else "FAIL", "details": _repo_path(ctx, contact_path)},
        {"check": "mapping_manifest_exists", "status": "PASS" if job_manifest_path.is_file() else "WARN", "details": _repo_path(ctx, job_manifest_path)},
        {"check": "chain_protomer_identity_preserved", "status": "FAIL" if blockers else "PASS", "details": ";".join(blockers)},
        {"check": "real_contact_evidence_present", "status": "PASS" if patch_rows else "WARN", "details": "patch_count={0}".format(len(patch_rows))},
    ]
    _write_csv(qc_csv, qc_rows, QC_FIELDS, ctx)

    status = _status(blockers, warnings, len(patch_rows))
    report = M2ConsensusPatchReport(
        run_id=ctx.run_id,
        mode=mode,
        profile=profile,
        status=status,
        contact_count=len(contacts),
        patch_count=len(patch_rows),
        accepted_pose_count=accepted_pose_count,
        pending_pose_qc_count=pending_pose_qc_count,
        warnings=sorted(set(warnings)),
        blockers=blockers,
        consensus_patch_csv=output_patch,
        m2_table_patch_csv=m2_table_patch,
        residue_csv=residue_csv,
        state_support_csv=state_support_csv,
        symmetry_support_csv=symmetry_support_csv,
        qc_csv_path=qc_csv,
    )
    summary = _write_summary(ctx, report)
    report.summary_report_path = summary
    manifest = ctx.require_within_run_dir(ctx.manifest_dir / "m2_4_ppi_consensus_patch_manifest.json")
    payload = {
        "run_id": ctx.run_id,
        "created_at": now_iso(),
        "phase": M2_4_PHASE,
        "mode": mode,
        "profile": profile,
        "status": status,
        "source_contact_table": _repo_path(ctx, contact_path),
        "source_job_manifest": _repo_path(ctx, job_manifest_path),
        "contact_count": len(contacts),
        "patch_count": len(patch_rows),
        "accepted_pose_count": accepted_pose_count,
        "pending_pose_qc_count": pending_pose_qc_count,
        "execution_allowed": False,
        "pyrosetta_imported": False,
        "docking_executed": False,
        "relaxation_executed": False,
        "pocket_discovery_executed": False,
        "compound_scoring_executed": False,
        "candidate_nomination_executed": False,
        "score_bonus_allowed": False,
        "outputs": {
            "consensus_patch_csv": _repo_path(ctx, output_patch),
            "phase1_ppi_consensus_patch_csv": _repo_path(ctx, m2_table_patch),
            "consensus_residues_csv": _repo_path(ctx, residue_csv),
            "state_support_csv": _repo_path(ctx, state_support_csv),
            "symmetry_support_csv": _repo_path(ctx, symmetry_support_csv),
            "qc_csv": _repo_path(ctx, qc_csv),
            "summary_report": _repo_path(ctx, summary),
        },
        "warnings": sorted(set(warnings)),
        "blockers": blockers,
        "non_goals_confirmed": NON_GOALS_CONFIRMED,
    }
    write_json(manifest, payload, ctx)
    report.manifest_path = manifest

    phase_status = "FAIL" if status == "FAIL" else (
        "WARN" if status == "PASS_WITH_WARNINGS" else "PASS"
    )
    append_phase_status(
        ctx,
        M2_4_PHASE,
        phase_status,
        "M2.4 consensus patch build {0} with {1} patch row(s)".format(status, len(patch_rows)),
        {"patch_count": len(patch_rows), "docking_executed": False},
    )
    append_job_status(
        ctx,
        "m2_4_ppi_consensus_patch_local",
        phase_status,
        details={"message": "local consensus patch post-processing only"},
    )
    return report
