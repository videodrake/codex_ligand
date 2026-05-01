"""M2 PPI pose contact extraction.

This module turns run-local posed EGFR-MYO1D complex PDB files into the raw
contact table consumed by M2.3 mapping restoration. It does not run docking or
score poses. The expected scientific convention is EGFR dimer chains A/B and
MYO1D partner chain C.
"""

from __future__ import annotations

import csv
import json
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

from egfr_myo1d.core.logging_utils import append_job_status, append_phase_status
from egfr_myo1d.core.manifest import load_yaml_config, now_iso, write_json
from egfr_myo1d.core.run_context import RunContext, RunContextError, ensure_within
from egfr_myo1d.ppi.collect_ppi_outputs import RAW_CONTACT_REQUIRED_FIELDS
from egfr_myo1d.ppi.pyrosetta_adapter import load_job_by_name
from egfr_myo1d.structure.pdb_parser import AtomRecord, PDBParseError, parse_pdb


M2_2_CONTACT_PHASE = "m2.2b-ppi-pose-contact-extraction"
DEFAULT_JOB_MANIFEST = "phase1_ppi/pyrosetta_adapter/pyrosetta_job_manifest.jsonl"
DEFAULT_RAW_CONTACT_TABLE = "phase1_ppi/tables/pyrosetta_ab_c_raw_contacts.csv"
RAW_CONTACT_FIELDS = [
    "run_id",
    "state_id",
    "state_role",
    "method",
    "seed",
    "job_name",
    "pose_id",
    "score",
    "pose_path",
    "protomer_id",
    "egfr_chain_id",
    "egfr_runtime_residue",
    "egfr_resname",
    "myo1d_residue",
    "myo1d_resname",
    "contact_type",
    "min_distance_angstrom",
    "egfr_atom",
    "myo1d_atom",
    "egfr_contact_centroid_x",
    "egfr_contact_centroid_y",
    "egfr_contact_centroid_z",
    "pose_acceptance_class",
    "accepted_pose_flag",
    "membrane_side_class",
]
QC_FIELDS = ["check", "status", "details"]
NON_GOALS_CONFIRMED = [
    "no_pyrosetta_import",
    "no_pyrosetta_docking_or_relaxation",
    "no_lightdock_execution",
    "no_vina_execution",
    "no_fpocket_or_p2rank_runtime",
    "no_compound_docking_or_scoring",
    "no_candidate_nomination",
    "no_scheduler_submission",
    "no_cleanup_deletion",
]
CHAIN_B_RUNTIME_OFFSET_THRESHOLD = 1500
CHAIN_B_RUNTIME_OFFSET = 1000


@dataclass(frozen=True)
class EgfrContactFilter:
    dockable_range: tuple[int, int] | None = None
    excluded_ranges: tuple[tuple[int, int], ...] = ()
    warnings: tuple[str, ...] = ()


@dataclass
class PoseContactExtractionReport:
    run_id: str
    mode: str
    profile: str
    status: str
    pose_count: int = 0
    raw_contact_count: int = 0
    accepted_pose_count: int = 0
    warnings: list[str] = field(default_factory=list)
    blockers: list[str] = field(default_factory=list)
    raw_contact_table: Path | None = None
    qc_csv_path: Path | None = None
    manifest_path: Path | None = None
    summary_report_path: Path | None = None


def _repo_path(ctx: RunContext, path: Path) -> str:
    return ctx.relative_to_repo(path)


def _resolve_run_path(ctx: RunContext, value: str | Path | None, default_text: str) -> Path:
    raw = Path(value or default_text)
    if raw.is_absolute():
        path = raw.resolve()
    else:
        repo_candidate = (ctx.repo_root / raw).resolve()
        try:
            path = ensure_within(repo_candidate, ctx.run_dir)
        except RunContextError:
            path = (ctx.run_dir / raw).resolve()
    return ensure_within(path, ctx.run_dir)


def _resolve_repo_or_abs(ctx: RunContext, text: str | Path) -> Path:
    path = Path(text)
    if path.is_absolute():
        return path.resolve()
    return (ctx.repo_root / path).resolve()


def _read_jsonl(path: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                rows.append(json.loads(line))
    return rows


def _read_csv(path: Path) -> list[dict[str, Any]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str], ctx: RunContext) -> None:
    safe = ctx.require_within_run_dir(path)
    safe.parent.mkdir(parents=True, exist_ok=True)
    with safe.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _is_heavy(atom: AtomRecord) -> bool:
    element = (atom.element or atom.atom_name[:1]).strip().upper()
    return element != "H"


def _distance(a: AtomRecord, b: AtomRecord) -> float:
    return math.sqrt((a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z - b.z) ** 2)


def _grid_key(atom: AtomRecord, cell: float) -> tuple[int, int, int]:
    return (
        math.floor(atom.x / cell),
        math.floor(atom.y / cell),
        math.floor(atom.z / cell),
    )


def _nearby_keys(key: tuple[int, int, int]) -> Iterable[tuple[int, int, int]]:
    x0, y0, z0 = key
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for dz in (-1, 0, 1):
                yield (x0 + dx, y0 + dy, z0 + dz)


def _ca_centroids(atoms: Iterable[AtomRecord]) -> dict[tuple[str, int], tuple[float, float, float]]:
    centroids: dict[tuple[str, int], tuple[float, float, float]] = {}
    for atom in atoms:
        if atom.atom_name == "CA":
            centroids[(atom.chain_id, atom.residue_number)] = (atom.x, atom.y, atom.z)
    return centroids


def _format_float(value: float) -> str:
    return "{0:.3f}".format(value)


def _pose_accepted(contact_count: int, min_contacts_for_acceptance: int) -> bool:
    return contact_count >= min_contacts_for_acceptance


def _parse_range(value: Any, label: str) -> tuple[tuple[int, int] | None, str | None]:
    if value in (None, ""):
        return None, None
    if isinstance(value, str):
        text = value.strip()
        if "-" in text:
            left, right = text.split("-", 1)
            try:
                start = int(left.strip())
                end = int(right.strip())
            except ValueError:
                return None, "invalid_residue_range:{0}:{1}".format(label, value)
        else:
            try:
                start = end = int(text)
            except ValueError:
                return None, "invalid_residue_range:{0}:{1}".format(label, value)
    elif isinstance(value, (list, tuple)) and len(value) == 2:
        try:
            start = int(value[0])
            end = int(value[1])
        except (TypeError, ValueError):
            return None, "invalid_residue_range:{0}:{1}".format(label, value)
    else:
        return None, "invalid_residue_range:{0}:{1}".format(label, value)
    if start > end:
        start, end = end, start
    return (start, end), None


def _load_egfr_contact_filter(ctx: RunContext) -> EgfrContactFilter:
    warnings: list[str] = []
    config_path = ctx.repo_root / "fresh" / "configs" / "gates.yaml"
    try:
        gates = load_yaml_config(config_path)
    except FileNotFoundError:
        return EgfrContactFilter(warnings=("missing_gates_yaml_for_contact_filter",))
    receptor = gates.get("receptor", {}) if isinstance(gates, dict) else {}
    dockable_range, warning = _parse_range(
        receptor.get("dockable_crop_default"), "dockable_crop_default"
    )
    if warning:
        warnings.append(warning)
    excluded_range, warning = _parse_range(
        receptor.get("excluded_tm_core_default"), "excluded_tm_core_default"
    )
    if warning:
        warnings.append(warning)
    excluded_ranges = (excluded_range,) if excluded_range else ()
    return EgfrContactFilter(
        dockable_range=dockable_range,
        excluded_ranges=excluded_ranges,
        warnings=tuple(warnings),
    )


def _egfr_original_like_residue(atom: AtomRecord) -> int:
    if atom.chain_id == "B" and atom.residue_number >= CHAIN_B_RUNTIME_OFFSET_THRESHOLD:
        return atom.residue_number - CHAIN_B_RUNTIME_OFFSET
    return atom.residue_number


def _egfr_contact_allowed(atom: AtomRecord, contact_filter: EgfrContactFilter) -> bool:
    residue = _egfr_original_like_residue(atom)
    if contact_filter.dockable_range is not None:
        start, end = contact_filter.dockable_range
        if residue < start or residue > end:
            return False
    for start, end in contact_filter.excluded_ranges:
        if start <= residue <= end:
            return False
    return True


def extract_pose_contacts(
    *,
    ctx: RunContext,
    pose_pdb: Path,
    job: dict[str, Any],
    contact_cutoff_angstrom: float = 6.0,
    partner_chain: str = "C",
    min_contacts_for_acceptance: int = 1,
    score: str = "",
    egfr_contact_filter: EgfrContactFilter | None = None,
) -> tuple[list[dict[str, Any]], list[str]]:
    """Extract EGFR A/B to MYO1D partner-chain contacts from one pose PDB."""

    warnings: list[str] = []
    try:
        structure = parse_pdb(pose_pdb)
    except PDBParseError as exc:
        return [], ["pose_parse_failed:{0}:{1}".format(pose_pdb, exc)]

    receptor_chains = {"A", "B"}
    atoms = [atom for atom in structure.atoms if atom.record_type in {"ATOM", "HETATM"}]
    contact_filter = egfr_contact_filter or EgfrContactFilter()
    receptor_atoms = [
        atom
        for atom in atoms
        if atom.chain_id in receptor_chains
        and _is_heavy(atom)
        and _egfr_contact_allowed(atom, contact_filter)
    ]
    partner_atoms = [
        atom for atom in atoms if atom.chain_id == partner_chain and _is_heavy(atom)
    ]
    if not receptor_atoms:
        warnings.append("pose_has_no_egfr_receptor_atoms:{0}".format(pose_pdb))
    if not partner_atoms:
        warnings.append("pose_has_no_partner_chain_{0}:{1}".format(partner_chain, pose_pdb))
    if not receptor_atoms or not partner_atoms:
        return [], warnings

    cell = max(contact_cutoff_angstrom, 0.1)
    grid: dict[tuple[int, int, int], list[AtomRecord]] = {}
    for atom in receptor_atoms:
        grid.setdefault(_grid_key(atom, cell), []).append(atom)

    best_by_pair: dict[tuple[str, int, int], tuple[float, AtomRecord, AtomRecord]] = {}
    cutoff_sq = contact_cutoff_angstrom * contact_cutoff_angstrom
    for partner_atom in partner_atoms:
        for key in _nearby_keys(_grid_key(partner_atom, cell)):
            for egfr_atom in grid.get(key, []):
                dx = egfr_atom.x - partner_atom.x
                dy = egfr_atom.y - partner_atom.y
                dz = egfr_atom.z - partner_atom.z
                dist_sq = dx * dx + dy * dy + dz * dz
                if dist_sq > cutoff_sq:
                    continue
                pair_key = (
                    egfr_atom.chain_id,
                    egfr_atom.residue_number,
                    partner_atom.residue_number,
                )
                dist = math.sqrt(dist_sq)
                previous = best_by_pair.get(pair_key)
                if previous is None or dist < previous[0]:
                    best_by_pair[pair_key] = (dist, egfr_atom, partner_atom)

    accepted = _pose_accepted(len(best_by_pair), min_contacts_for_acceptance)
    acceptance_class = (
        "accepted_contact_geometry_unscored"
        if accepted
        else "pending_pose_qc_geometry_insufficient_contacts"
    )
    ca_centroids = _ca_centroids(receptor_atoms)
    rows: list[dict[str, Any]] = []
    state_id = str(job.get("state_id", ""))
    pose_id = pose_pdb.stem
    for _pair_key, (dist, egfr_atom, partner_atom) in sorted(best_by_pair.items()):
        protomer_id = egfr_atom.chain_id
        centroid = ca_centroids.get(
            (egfr_atom.chain_id, egfr_atom.residue_number),
            (egfr_atom.x, egfr_atom.y, egfr_atom.z),
        )
        rows.append(
            {
                "run_id": ctx.run_id,
                "state_id": state_id,
                "state_role": str(job.get("state_role", "primary_membrane_validated_state")),
                "method": str(job.get("method", "pyrosetta_global_ppi")),
                "seed": str(job.get("seed", "")),
                "job_name": str(job.get("job_name", "")),
                "pose_id": pose_id,
                "score": str(score),
                "pose_path": _repo_path(ctx, pose_pdb),
                "protomer_id": protomer_id,
                "egfr_chain_id": egfr_atom.chain_id,
                "egfr_runtime_residue": str(egfr_atom.residue_number),
                "egfr_resname": egfr_atom.resname,
                "myo1d_residue": str(partner_atom.residue_number),
                "myo1d_resname": partner_atom.resname,
                "contact_type": "heavy_atom_contact",
                "min_distance_angstrom": _format_float(dist),
                "egfr_atom": egfr_atom.atom_name,
                "myo1d_atom": partner_atom.atom_name,
                "egfr_contact_centroid_x": _format_float(centroid[0]),
                "egfr_contact_centroid_y": _format_float(centroid[1]),
                "egfr_contact_centroid_z": _format_float(centroid[2]),
                "pose_acceptance_class": acceptance_class,
                "accepted_pose_flag": "true" if accepted else "false",
                "membrane_side_class": "not_classified",
            }
        )
    return rows, warnings


def _load_jobs(ctx: RunContext, job_manifest: Path | None) -> tuple[list[dict[str, Any]], Path, list[str]]:
    path = _resolve_run_path(ctx, job_manifest, DEFAULT_JOB_MANIFEST)
    if not path.is_file():
        return [], path, ["missing_m2_2_job_manifest:{0}".format(path)]
    return _read_jsonl(path), path, []


def _pose_files_for_job(ctx: RunContext, job: dict[str, Any], pose_root: Path | None) -> list[Path]:
    candidates: list[Path] = []
    if pose_root is not None:
        root = ensure_within(pose_root.resolve(), ctx.run_dir)
        candidates.extend(
            [
                root / str(job.get("job_name", "")) / "poses",
                root / str(job.get("job_name", "")),
            ]
        )
    output_text = str(job.get("output_dir", ""))
    if output_text:
        output_dir = ensure_within(_resolve_repo_or_abs(ctx, output_text), ctx.run_dir)
        candidates.extend([output_dir / "poses", output_dir])

    seen: set[Path] = set()
    files: list[Path] = []
    for candidate in candidates:
        if not candidate.is_dir():
            continue
        for pdb in sorted(candidate.glob("*.pdb")):
            resolved = ensure_within(pdb.resolve(), ctx.run_dir)
            if resolved not in seen:
                seen.add(resolved)
                files.append(resolved)
    return files


def _score_by_pose_id_for_job(ctx: RunContext, job: dict[str, Any]) -> dict[str, str]:
    output_text = str(job.get("output_dir", ""))
    if not output_text:
        return {}
    output_dir = ensure_within(_resolve_repo_or_abs(ctx, output_text), ctx.run_dir)
    candidates = [output_dir / "pose_scores.csv"]
    chunks_dir = output_dir / "chunks"
    if chunks_dir.is_dir():
        candidates.extend(sorted(chunks_dir.glob("*/pose_scores.csv")))

    scores: dict[str, str] = {}
    for candidate in candidates:
        if not candidate.is_file():
            continue
        for row in _read_csv(candidate):
            pose_id = str(row.get("pose_id", ""))
            if not pose_id:
                continue
            scores[pose_id] = str(row.get("score", ""))
    return scores


def _write_summary(ctx: RunContext, report: PoseContactExtractionReport) -> Path:
    target = ctx.require_within_run_dir(ctx.reports_dir / "m2_2b_ppi_pose_contact_extraction.md")
    lines = [
        "# M2.2b PPI Pose Contact Extraction",
        "",
        "Status: {0}".format(report.status),
        "Pose PDBs parsed: {0}".format(report.pose_count),
        "Raw contact rows: {0}".format(report.raw_contact_count),
        "Accepted pose count: {0}".format(report.accepted_pose_count),
        "",
        "This phase parses existing posed complex PDB files only. It does not run docking, relaxation, or scoring.",
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


def extract_m2_pose_contacts(
    ctx: RunContext,
    *,
    mode: str = "production",
    profile: str = "hpc_strict",
    job_manifest: Path | None = None,
    pose_root: Path | None = None,
    output_csv: Path | None = None,
    contact_cutoff_angstrom: float = 6.0,
    partner_chain: str = "C",
    min_contacts_for_acceptance: int = 1,
) -> PoseContactExtractionReport:
    """Build an M2.3-compatible raw contact table from run-local pose PDBs."""

    ctx.create_directories()
    jobs, job_manifest_path, blockers = _load_jobs(ctx, job_manifest)
    resolved_pose_root = (
        ensure_within(_resolve_repo_or_abs(ctx, pose_root), ctx.run_dir)
        if pose_root is not None
        else None
    )
    output_path = _resolve_run_path(ctx, output_csv, DEFAULT_RAW_CONTACT_TABLE)
    egfr_contact_filter = _load_egfr_contact_filter(ctx)

    rows: list[dict[str, Any]] = []
    warnings: list[str] = list(egfr_contact_filter.warnings)
    pose_count = 0
    accepted_pose_ids: set[tuple[str, str]] = set()
    for job in jobs:
        pose_files = _pose_files_for_job(ctx, job, resolved_pose_root)
        if not pose_files:
            warnings.append("no_pose_pdbs_for_job:{0}".format(job.get("job_name", "")))
            continue
        score_by_pose_id = _score_by_pose_id_for_job(ctx, job)
        for pose_file in pose_files:
            pose_count += 1
            pose_rows, pose_warnings = extract_pose_contacts(
                ctx=ctx,
                pose_pdb=pose_file,
                job=job,
                contact_cutoff_angstrom=contact_cutoff_angstrom,
                partner_chain=partner_chain,
                min_contacts_for_acceptance=min_contacts_for_acceptance,
                score=score_by_pose_id.get(pose_file.stem, ""),
                egfr_contact_filter=egfr_contact_filter,
            )
            rows.extend(pose_rows)
            warnings.extend(pose_warnings)
            if any(row.get("accepted_pose_flag") == "true" for row in pose_rows):
                accepted_pose_ids.add((str(job.get("job_name", "")), pose_file.stem))

    for required in RAW_CONTACT_REQUIRED_FIELDS:
        if required not in RAW_CONTACT_FIELDS:
            blockers.append("extractor_missing_required_raw_contact_field:{0}".format(required))

    _write_csv(output_path, rows, RAW_CONTACT_FIELDS, ctx)
    qc_path = ctx.qc_dir / "m2_2b_ppi_pose_contact_extraction_qc.csv"
    qc_rows = [
        {
            "check": "job_manifest_exists",
            "status": "PASS" if job_manifest_path.is_file() else "FAIL",
            "details": _repo_path(ctx, job_manifest_path),
        },
        {
            "check": "pose_pdbs_present",
            "status": "PASS" if pose_count else "WARN",
            "details": "pose_count={0}".format(pose_count),
        },
        {
            "check": "raw_contacts_present",
            "status": "PASS" if rows else "WARN",
            "details": "raw_contact_count={0}".format(len(rows)),
        },
        {
            "check": "partner_chain",
            "status": "PASS",
            "details": partner_chain,
        },
        {
            "check": "egfr_contact_filter",
            "status": "PASS" if not egfr_contact_filter.warnings else "WARN",
            "details": "dockable={0}; excluded={1}".format(
                egfr_contact_filter.dockable_range,
                list(egfr_contact_filter.excluded_ranges),
            ),
        },
    ]
    _write_csv(qc_path, qc_rows, QC_FIELDS, ctx)

    status = "FAIL" if blockers else ("PASS_WITH_WARNINGS" if warnings or not rows else "PASS")
    report = PoseContactExtractionReport(
        run_id=ctx.run_id,
        mode=mode,
        profile=profile,
        status=status,
        pose_count=pose_count,
        raw_contact_count=len(rows),
        accepted_pose_count=len(accepted_pose_ids),
        warnings=sorted(set(warnings)),
        blockers=blockers,
        raw_contact_table=output_path,
        qc_csv_path=qc_path,
    )
    summary = _write_summary(ctx, report)
    report.summary_report_path = summary
    manifest = ctx.manifest_dir / "m2_2b_ppi_pose_contact_extraction_manifest.json"
    payload = {
        "run_id": ctx.run_id,
        "created_at": now_iso(),
        "phase": M2_2_CONTACT_PHASE,
        "mode": mode,
        "profile": profile,
        "status": status,
        "source_job_manifest": _repo_path(ctx, job_manifest_path),
        "pose_root": _repo_path(ctx, resolved_pose_root) if resolved_pose_root else None,
        "contact_cutoff_angstrom": contact_cutoff_angstrom,
        "partner_chain": partner_chain,
        "min_contacts_for_acceptance": min_contacts_for_acceptance,
        "egfr_contact_filter": {
            "dockable_range": list(egfr_contact_filter.dockable_range)
            if egfr_contact_filter.dockable_range
            else None,
            "excluded_ranges": [list(item) for item in egfr_contact_filter.excluded_ranges],
            "chain_b_runtime_offset_threshold": CHAIN_B_RUNTIME_OFFSET_THRESHOLD,
            "chain_b_runtime_offset": CHAIN_B_RUNTIME_OFFSET,
        },
        "pose_count": pose_count,
        "raw_contact_count": len(rows),
        "accepted_pose_count": len(accepted_pose_ids),
        "outputs": {
            "raw_contact_table": _repo_path(ctx, output_path),
            "qc_csv": _repo_path(ctx, qc_path),
            "summary_report": _repo_path(ctx, summary),
        },
        "warnings": sorted(set(warnings)),
        "blockers": blockers,
        "non_goals_confirmed": NON_GOALS_CONFIRMED,
    }
    write_json(manifest, payload, ctx)
    report.manifest_path = manifest

    phase_status = "FAIL" if status == "FAIL" else ("WARN" if status == "PASS_WITH_WARNINGS" else "PASS")
    append_phase_status(
        ctx,
        M2_2_CONTACT_PHASE,
        phase_status,
        "M2.2b PPI pose contact extraction {0} with {1} raw contact row(s)".format(
            status, len(rows)
        ),
        {"raw_contact_table": _repo_path(ctx, output_path)},
    )
    append_job_status(
        ctx,
        "m2_2b_ppi_pose_contact_extraction_local",
        phase_status,
        details={"raw_contact_count": len(rows), "docking_executed": False},
    )
    return report


def load_job_for_contact_runner(ctx: RunContext, job_name: str) -> dict[str, Any]:
    """Wrapper for the real runner so imports stay local to the PPI package."""
    return load_job_by_name(ctx, job_name)
