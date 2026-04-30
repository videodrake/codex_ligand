"""M2.6 fpocket adapter and raw pocket parser.

The adapter can run fpocket in production mode, but the default test/smoke path
is parser-only. No M2.7 gate, compound docking, or candidate nomination is
performed here.
"""

from __future__ import annotations

import csv
import json
import math
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

from egfr_myo1d.core.logging_utils import append_job_status, append_phase_status
from egfr_myo1d.core.manifest import load_yaml_config, now_iso, write_json
from egfr_myo1d.core.run_context import RunContext, RunContextError, ensure_within
from egfr_myo1d.io.residue_mapping import read_residue_mapping
from egfr_myo1d.pocket.family_merge import MEMBER_FIELDS, MERGED_FIELDS, merge_pocket_families
from egfr_myo1d.pocket.normalization import RAW_CANDIDATE_FIELDS, normalize_raw_pockets
from egfr_myo1d.structure.pdb_parser import AtomRecord, PDBParseError, parse_pdb


M2_6_PHASE = "m2.6-fpocket-adapter-and-pocket-normalization"
RAW_ROOT = "phase2_pockets/raw"
MERGED_ROOT = "phase2_pockets/merged"
TABLE_ROOT = "phase2_pockets/tables"
RAW_FIELDS = [
    "state_id",
    "state_role",
    "receptor_id",
    "receptor_path",
    "source_tool",
    "source_tool_version",
    "source_pocket_id",
    "source_pocket_dir",
    "fpocket_score",
    "fpocket_druggability_score",
    "fpocket_volume",
    "num_alpha_spheres",
    "num_pocket_atoms",
    "centroid_x",
    "centroid_y",
    "centroid_z",
    "raw_chain_ids",
    "raw_residue_ids",
    "raw_residue_names",
    "egfr_chain_ids",
    "egfr_protomer_ids",
    "mapped_uniprot_residues",
    "mapping_status",
    "parser_status",
    "evidence_status",
    "warnings",
]
QC_FIELDS = ["check", "status", "details"]
NON_GOALS_CONFIRMED = [
    "no_pyrosetta_import",
    "no_pyrosetta_docking_or_relaxation",
    "no_lightdock_execution",
    "no_p2rank_execution",
    "no_vina_execution",
    "no_compound_docking_or_scoring",
    "no_candidate_nomination",
    "no_m2_7_gate_application",
    "no_accepted_pocket_export",
    "no_scheduler_submission",
    "no_cleanup_deletion",
]
FORBIDDEN_M2_7_OUTPUTS = [
    "phase2_pockets/tables/pocket_gate_qc.csv",
    "phase2_pockets/tables/pocket_rejects.csv",
    "phase2_pockets/tables/accepted_pockets_for_m3.csv",
    "phase2_pockets/export_for_m3/accepted_pocket_boxes.csv",
]


@dataclass
class FpocketDiscoveryReport:
    run_id: str
    mode: str
    profile: str
    execution_mode: str
    status: str
    raw_pocket_count: int = 0
    raw_candidate_count: int = 0
    merged_family_count: int = 0
    warnings: list[str] = field(default_factory=list)
    blockers: list[str] = field(default_factory=list)
    raw_status_json: Path | None = None
    normalization_status_json: Path | None = None
    manifest_path: Path | None = None
    summary_report_path: Path | None = None
    raw_csv_paths: list[Path] = field(default_factory=list)
    raw_candidate_csv: Path | None = None
    merged_candidate_csv: Path | None = None
    member_csv: Path | None = None


def _repo_path(ctx: RunContext, path: Path) -> str:
    return ctx.relative_to_repo(path)


def _is_within(child: Path, parent: Path) -> bool:
    try:
        Path(child).resolve().relative_to(Path(parent).resolve())
        return True
    except ValueError:
        return False


def _resolve_fresh_or_run_path(ctx: RunContext, value: str | Path) -> Path:
    raw = Path(value)
    path = raw.resolve() if raw.is_absolute() else (ctx.repo_root / raw).resolve()
    if not (_is_within(path, ctx.fresh_root) or _is_within(path, ctx.run_dir)):
        raise RunContextError("input path escapes fresh/ and run dir: {0}".format(path))
    return path


def _run_relative_path(ctx: RunContext, value: str | Path) -> Path:
    raw = Path(value)
    path = raw.resolve() if raw.is_absolute() else (ctx.run_dir / raw).resolve()
    return ensure_within(path, ctx.run_dir)


def _write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str], ctx: RunContext) -> None:
    safe = ctx.require_within_run_dir(path)
    safe.parent.mkdir(parents=True, exist_ok=True)
    with safe.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _load_config(ctx: RunContext, config_path: Path | None) -> dict[str, Any]:
    default = ctx.repo_root / "fresh" / "configs" / "pocket.yaml"
    path = _resolve_fresh_or_run_path(ctx, config_path or default)
    if path.is_file():
        return load_yaml_config(path)
    return {
        "fpocket": {"binary": "fpocket"},
        "merge": {
            "centroid_merge_threshold_angstrom": 6.0,
            "residue_jaccard_merge_threshold": 0.30,
        },
    }


def _load_receptor_states(ctx: RunContext) -> dict[str, Any]:
    return load_yaml_config(ctx.repo_root / "fresh" / "configs" / "receptor_states.yaml")


def _default_states(ctx: RunContext, include_reference_states: bool) -> list[str]:
    data = _load_receptor_states(ctx)
    states = list(data.get("primary_membrane_validated_states", []))
    if include_reference_states:
        states.extend(data.get("reference_control_states", []))
    return list(dict.fromkeys(states))


def _state_role_map(ctx: RunContext) -> dict[str, str]:
    data = _load_receptor_states(ctx)
    roles: dict[str, str] = {}
    for state_id, role in (data.get("state_policy", {}) or {}).items():
        role_text = str(role)
        if "reference" in role_text or "control" in role_text:
            roles[str(state_id)] = "reference"
        elif "primary" in role_text:
            roles[str(state_id)] = "primary"
        else:
            roles[str(state_id)] = role_text or "unknown"
    return roles


def _state_receptor_path(ctx: RunContext, state_id: str) -> Path | None:
    candidates = [
        ctx.run_dir / "normalized" / "receptors" / "{0}_dockable_669_1014_explicit_AB.pdb".format(state_id),
        ctx.run_dir / "normalized" / "receptors" / "{0}_dockable_reference_explicit_AB.pdb".format(state_id),
        ctx.run_dir / "normalized" / "receptors" / "{0}_full_frame_explicit_AB.pdb".format(state_id),
        ctx.run_dir / "normalized" / "receptors" / "{0}_runtime_offset_receptor_only.pdb".format(state_id),
    ]
    for candidate in candidates:
        if candidate.is_file():
            return ctx.require_within_run_dir(candidate)
    return None


def _mapping_path(ctx: RunContext, state_id: str) -> Path | None:
    path = ctx.qc_dir / "{0}_receptor_mapping.csv".format(state_id)
    if path.is_file():
        return ctx.require_within_run_dir(path)
    return None


def _is_heavy(atom: AtomRecord) -> bool:
    element = (atom.element or atom.atom_name[:1]).upper()
    return element != "H"


def _format_float(value: float | None) -> str:
    return "" if value is None else "{0:.3f}".format(value)


def _safe_float(value: Any) -> float | None:
    try:
        text = str(value).strip()
        if not text:
            return None
        return float(text)
    except (TypeError, ValueError):
        return None


def _mean(values: Iterable[float]) -> float | None:
    items = list(values)
    if not items:
        return None
    return sum(items) / float(len(items))


def _stage_receptor(ctx: RunContext, state_id: str, target_dir: Path) -> tuple[Path | None, list[str]]:
    warnings: list[str] = []
    receptor = _state_receptor_path(ctx, state_id)
    if receptor is None:
        return None, ["missing_normalized_receptor_for_state:{0}".format(state_id)]
    safe_target_dir = ctx.require_within_run_dir(target_dir)
    safe_target_dir.mkdir(parents=True, exist_ok=True)
    staged = safe_target_dir / "{0}_fpocket_input.pdb".format(state_id)
    shutil.copy2(receptor, staged)
    return staged, warnings


def _find_fpocket_output_dir(root: Path, state_id: str) -> Path | None:
    candidates = [
        root,
        root / "{0}_out".format(state_id),
        root / "{0}_fpocket_out".format(state_id),
        root / "{0}_fpocket_input_out".format(state_id),
        root / state_id,
        root / state_id / "{0}_out".format(state_id),
        root / state_id / "{0}_fpocket_out".format(state_id),
        root / state_id / "{0}_fpocket_input_out".format(state_id),
    ]
    for candidate in candidates:
        if candidate.is_dir() and ((candidate / "pockets_info.txt").is_file() or (candidate / "pockets").is_dir()):
            return candidate.resolve()
    return None


def _copy_parser_output(ctx: RunContext, source_dir: Path, state_output_dir: Path) -> Path:
    source = _resolve_fresh_or_run_path(ctx, source_dir)
    target_root = ctx.require_within_run_dir(state_output_dir)
    target_root.mkdir(parents=True, exist_ok=True)
    if _is_within(source, target_root):
        return source
    target = target_root / source.name
    if target.exists():
        shutil.rmtree(target)
    shutil.copytree(source, target)
    return target


def _fpocket_version(binary: str) -> str:
    resolved = shutil.which(binary)
    if not resolved:
        return "unavailable"
    try:
        proc = subprocess.run(
            [resolved, "-h"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            timeout=10,
            check=False,
        )
    except (OSError, subprocess.TimeoutExpired):
        return "unknown"
    text = (proc.stdout or proc.stderr or "").strip().splitlines()
    return text[0][:120] if text else "unknown"


def _run_fpocket(
    ctx: RunContext,
    state_id: str,
    state_dir: Path,
    staged_receptor: Path,
    binary: str,
) -> tuple[Path | None, dict[str, Any], list[str], list[str]]:
    warnings: list[str] = []
    blockers: list[str] = []
    resolved = shutil.which(binary)
    metadata = {
        "state_id": state_id,
        "fpocket_binary": binary,
        "fpocket_binary_path": resolved,
        "fpocket_version": _fpocket_version(binary),
        "command": [binary, "-f", staged_receptor.name],
        "return_code": None,
        "stdout_log": None,
        "stderr_log": None,
    }
    if not resolved:
        blockers.append("fpocket_binary_unavailable:{0}".format(binary))
        return None, metadata, warnings, blockers
    stdout_path = ctx.require_within_run_dir(ctx.jobs_log_dir / "m2_6_fpocket_{0}_stdout.log".format(state_id))
    stderr_path = ctx.require_within_run_dir(ctx.jobs_log_dir / "m2_6_fpocket_{0}_stderr.log".format(state_id))
    proc = subprocess.run(
        [resolved, "-f", staged_receptor.name],
        cwd=str(ctx.require_within_run_dir(state_dir)),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )
    stdout_path.write_text(proc.stdout or "", encoding="utf-8")
    stderr_path.write_text(proc.stderr or "", encoding="utf-8")
    metadata.update(
        {
            "return_code": proc.returncode,
            "stdout_log": _repo_path(ctx, stdout_path),
            "stderr_log": _repo_path(ctx, stderr_path),
        }
    )
    if proc.returncode != 0:
        blockers.append("fpocket_failed:{0}:return_code_{1}".format(state_id, proc.returncode))
        return None, metadata, warnings, blockers
    output_dir = _find_fpocket_output_dir(state_dir, state_id)
    if output_dir is None:
        blockers.append("fpocket_output_missing_after_execution:{0}".format(state_id))
    return output_dir, metadata, warnings, blockers


def _parse_pockets_info(path: Path) -> dict[str, dict[str, str]]:
    if not path.is_file():
        return {}
    info: dict[str, dict[str, str]] = {}
    current: str | None = None
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        lower = stripped.lower()
        if lower.startswith("pocket"):
            parts = stripped.replace(":", " ").split()
            if len(parts) >= 2 and parts[1].isdigit():
                current = "pocket{0}".format(parts[1])
                info.setdefault(current, {})
                continue
        if current is None or ":" not in stripped:
            continue
        key, value = [item.strip() for item in stripped.split(":", 1)]
        key_lower = key.lower()
        if "druggability" in key_lower and "score" in key_lower:
            info[current]["fpocket_druggability_score"] = value
        elif "score" in key_lower:
            info[current]["fpocket_score"] = value
        elif "volume" in key_lower:
            info[current]["fpocket_volume"] = value
        elif "alpha" in key_lower and ("sphere" in key_lower or "spheres" in key_lower):
            info[current]["num_alpha_spheres"] = value
    return info


def _pocket_number_from_file(path: Path) -> str:
    name = path.stem
    digits = ""
    for char in name:
        if char.isdigit():
            digits += char
        elif digits:
            break
    return "pocket{0}".format(digits) if digits else name


def _residue_key(atom: AtomRecord) -> tuple[str, int, str]:
    return (atom.chain_id, atom.residue_number, atom.resname)


def _mapping_indexes(mapping_path: Path | None) -> tuple[dict[tuple[str, int], dict[str, Any]], dict[tuple[str, int], dict[str, Any]]]:
    by_source: dict[tuple[str, int], dict[str, Any]] = {}
    by_runtime: dict[tuple[str, int], dict[str, Any]] = {}
    if mapping_path is None:
        return by_source, by_runtime
    for row in read_residue_mapping(mapping_path):
        source_resseq = _safe_int(row.get("source_resseq"))
        runtime_resseq = _safe_int(row.get("runtime_resseq"))
        source_chain = str(row.get("source_chain", ""))
        runtime_chain = str(row.get("runtime_chain", ""))
        if source_chain and source_resseq is not None:
            by_source[(source_chain, source_resseq)] = row
        if runtime_chain and runtime_resseq is not None:
            by_runtime[(runtime_chain, runtime_resseq)] = row
    return by_source, by_runtime


def _safe_int(value: Any) -> int | None:
    try:
        return int(str(value).strip())
    except (TypeError, ValueError):
        return None


def _map_residues(
    residues: set[tuple[str, int, str]],
    mapping_path: Path | None,
) -> tuple[dict[str, Any], list[str]]:
    warnings: list[str] = []
    if mapping_path is None:
        return {
            "egfr_chain_ids": "",
            "egfr_protomer_ids": "",
            "mapped_uniprot_residues": "",
            "mapping_status": "MISSING",
        }, ["missing_receptor_mapping"]
    by_source, by_runtime = _mapping_indexes(mapping_path)
    chains: set[str] = set()
    protomers: set[str] = set()
    uniprot_residues: set[str] = set()
    mapped = 0
    mismatches = 0
    drift = 0
    for chain, residue_number, resname in sorted(residues):
        row = by_source.get((chain, residue_number))
        used_runtime = False
        if row is None:
            row = by_runtime.get((chain, residue_number))
            used_runtime = row is not None
        if row is None:
            warnings.append("unmapped_pocket_residue:{0}:{1}".format(chain, residue_number))
            continue
        mapped += 1
        source_chain = str(row.get("source_chain", "")) or chain
        protomer_id = str(row.get("protomer_id", ""))
        source_resseq = str(row.get("source_resseq", ""))
        chains.add(source_chain)
        if protomer_id:
            protomers.add(protomer_id)
        if source_resseq:
            uniprot_residues.add(source_resseq)
        source_resname = str(row.get("source_resname", ""))
        if source_resname and resname and source_resname.upper() != resname.upper():
            mismatches += 1
            warnings.append(
                "residue_name_mismatch:{0}:{1}:{2}:{3}".format(
                    chain, residue_number, resname, source_resname
                )
            )
        if used_runtime:
            drift += 1
            warnings.append(
                "residue_number_drift_runtime_to_source:{0}:{1}->{2}".format(
                    chain, residue_number, source_resseq
                )
            )
    if not residues:
        status = "NO_RESIDUES"
    elif mapped == 0:
        status = "MISSING"
    elif mapped < len(residues):
        status = "PARTIAL"
    elif mismatches or drift:
        status = "WARN"
    else:
        status = "PASS"
    return {
        "egfr_chain_ids": ";".join(sorted(chains)),
        "egfr_protomer_ids": ";".join(sorted(protomers)),
        "mapped_uniprot_residues": ";".join(sorted(uniprot_residues, key=lambda value: (not value.isdigit(), value))),
        "mapping_status": status,
    }, warnings


def parse_fpocket_output(
    ctx: RunContext,
    state_id: str,
    state_role: str,
    receptor_path: Path | None,
    output_dir: Path,
    mapping_path: Path | None,
    source_tool_version: str = "not_run_parser_only",
) -> tuple[list[dict[str, Any]], list[str]]:
    """Parse a single fpocket output directory into raw M2.6 rows."""

    safe_output_dir = _resolve_fresh_or_run_path(ctx, output_dir)
    info = _parse_pockets_info(safe_output_dir / "pockets_info.txt")
    pockets_dir = safe_output_dir / "pockets"
    warnings: list[str] = []
    if not pockets_dir.is_dir():
        return [], ["missing_fpocket_pockets_dir:{0}".format(safe_output_dir)]
    pocket_files = sorted(pockets_dir.glob("pocket*_atm.pdb"))
    if not pocket_files:
        return [], ["no_fpocket_pocket_atom_files:{0}".format(pockets_dir)]
    rows: list[dict[str, Any]] = []
    for pocket_file in pocket_files:
        source_pocket_id = _pocket_number_from_file(pocket_file)
        row_warnings: list[str] = []
        try:
            structure = parse_pdb(pocket_file)
            atoms = [atom for atom in structure.atoms if _is_heavy(atom)]
        except PDBParseError as exc:
            row_warnings.append("pocket_atom_parse_failed:{0}".format(exc))
            atoms = []
        residues = {_residue_key(atom) for atom in atoms if atom.record_type in {"ATOM", "HETATM"}}
        mapping, mapping_warnings = _map_residues(residues, mapping_path)
        row_warnings.extend(mapping_warnings)
        xs = [atom.x for atom in atoms]
        ys = [atom.y for atom in atoms]
        zs = [atom.z for atom in atoms]
        centroid_x = _mean(xs)
        centroid_y = _mean(ys)
        centroid_z = _mean(zs)
        if not atoms:
            row_warnings.append("pocket_atoms_unavailable_for_centroid")
        raw_residue_ids = ["{0}:{1}".format(chain, residue_number) for chain, residue_number, _ in sorted(residues)]
        raw_residue_names = ["{0}:{1}:{2}".format(chain, residue_number, resname) for chain, residue_number, resname in sorted(residues)]
        metrics = info.get(source_pocket_id, {})
        parser_status = "PASS" if not row_warnings else "WARN"
        if mapping.get("mapping_status") in {"MISSING", "NO_RESIDUES"}:
            parser_status = "WARN"
        evidence_status = "REFERENCE_ONLY_RAW_UNGATED" if state_role == "reference" else "RAW_UNGATED"
        rows.append(
            {
                "state_id": state_id,
                "state_role": state_role,
                "receptor_id": "EGFR_MYO1D_membrane_compatible_inactive_dimer",
                "receptor_path": _repo_path(ctx, receptor_path) if receptor_path else "",
                "source_tool": "fpocket",
                "source_tool_version": source_tool_version,
                "source_pocket_id": source_pocket_id,
                "source_pocket_dir": _repo_path(ctx, safe_output_dir),
                "fpocket_score": metrics.get("fpocket_score", ""),
                "fpocket_druggability_score": metrics.get("fpocket_druggability_score", ""),
                "fpocket_volume": metrics.get("fpocket_volume", ""),
                "num_alpha_spheres": metrics.get("num_alpha_spheres", ""),
                "num_pocket_atoms": len(atoms),
                "centroid_x": _format_float(centroid_x),
                "centroid_y": _format_float(centroid_y),
                "centroid_z": _format_float(centroid_z),
                "raw_chain_ids": ";".join(sorted({chain for chain, _, _ in residues})),
                "raw_residue_ids": ";".join(raw_residue_ids),
                "raw_residue_names": ";".join(raw_residue_names),
                "egfr_chain_ids": mapping.get("egfr_chain_ids", ""),
                "egfr_protomer_ids": mapping.get("egfr_protomer_ids", ""),
                "mapped_uniprot_residues": mapping.get("mapped_uniprot_residues", ""),
                "mapping_status": mapping.get("mapping_status", ""),
                "parser_status": parser_status,
                "evidence_status": evidence_status,
                "warnings": ";".join(sorted(set(row_warnings))),
                "_extent_x": _format_float(max(xs) - min(xs) if xs else None),
                "_extent_y": _format_float(max(ys) - min(ys) if ys else None),
                "_extent_z": _format_float(max(zs) - min(zs) if zs else None),
            }
        )
    return rows, warnings


def _precondition_warnings_and_blockers(ctx: RunContext, execution_mode: str) -> tuple[list[str], list[str], dict[str, bool]]:
    warnings: list[str] = []
    blockers: list[str] = []
    consensus = ctx.run_dir / "output" / "ppi" / "ppi_consensus_patch.csv"
    atp = ctx.run_dir / "phase2_pockets" / "atp_reference" / "atp_site_reference.csv"
    evidence = {
        "m2_4_consensus_patch_exists": consensus.is_file(),
        "m2_5_atp_reference_exists": atp.is_file(),
    }
    if not consensus.is_file():
        message = "m2_4_consensus_patch_missing"
        if execution_mode == "production":
            blockers.append(message + "_production_blocker")
        else:
            warnings.append(message + "_parser_only_allowed")
    if not atp.is_file():
        warnings.append("m2_5_atp_reference_missing_m2_7_atp_gate_blocked")
    return warnings, blockers, evidence


def _status(blockers: list[str], warnings: list[str], raw_candidate_count: int) -> str:
    if blockers:
        return "FAIL"
    if warnings or raw_candidate_count == 0:
        return "PASS_WITH_WARNINGS"
    return "PASS"


def _write_summary(ctx: RunContext, report: FpocketDiscoveryReport) -> Path:
    target = ctx.require_within_run_dir(ctx.reports_dir / "m2_6_fpocket_adapter_and_pocket_normalization.md")
    lines = [
        "# M2.6 fpocket Adapter and Pocket Normalization",
        "",
        "Status: {0}".format(report.status),
        "Execution mode: {0}".format(report.execution_mode),
        "Raw fpocket pockets: {0}".format(report.raw_pocket_count),
        "Raw candidate rows: {0}".format(report.raw_candidate_count),
        "Merged raw families: {0}".format(report.merged_family_count),
        "",
        "All candidates remain raw_ungated. M2.7 gates and accepted pocket exports were not run.",
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


def discover_and_normalize_pockets(
    ctx: RunContext,
    mode: str = "smoke_input",
    profile: str = "codex_dev",
    states: Iterable[str] | None = None,
    include_reference_states: bool = True,
    execution_mode: str = "parser_only",
    fpocket_output_root: Path | None = None,
    config_path: Path | None = None,
    fpocket_binary: str | None = None,
    synthetic_fixture: bool = False,
) -> FpocketDiscoveryReport:
    """Run or parse fpocket outputs and create M2.6 raw/merged tables."""

    ctx.create_directories()
    if execution_mode not in {"parser_only", "production"}:
        raise ValueError("execution_mode must be parser_only or production")
    config = _load_config(ctx, config_path)
    merge_config = config.get("merge", {}) or {}
    fpocket_config = config.get("fpocket", {}) or {}
    binary = fpocket_binary or str(fpocket_config.get("binary", "fpocket"))
    centroid_threshold = float(merge_config.get("centroid_merge_threshold_angstrom", 6.0))
    jaccard_threshold = float(merge_config.get("residue_jaccard_merge_threshold", 0.30))

    state_list = list(states) if states else _default_states(ctx, include_reference_states)
    roles = _state_role_map(ctx)
    raw_root = ctx.require_within_run_dir(ctx.run_dir / RAW_ROOT)
    merged_root = ctx.require_within_run_dir(ctx.run_dir / MERGED_ROOT)
    raw_root.mkdir(parents=True, exist_ok=True)
    merged_root.mkdir(parents=True, exist_ok=True)

    warnings, blockers, evidence = _precondition_warnings_and_blockers(ctx, execution_mode)
    fpocket_runs: list[dict[str, Any]] = []
    all_raw_rows: list[dict[str, Any]] = []
    raw_csv_paths: list[Path] = []
    fpocket_executed = False
    resolved_output_root = (
        _resolve_fresh_or_run_path(ctx, fpocket_output_root)
        if fpocket_output_root is not None
        else None
    )

    for state_id in state_list:
        state_role = roles.get(state_id, "unknown")
        state_dir = ctx.require_within_run_dir(raw_root / "{0}_fpocket_pockets".format(state_id))
        state_dir.mkdir(parents=True, exist_ok=True)
        staged_receptor, stage_warnings = _stage_receptor(ctx, state_id, state_dir)
        if stage_warnings:
            warnings.extend(stage_warnings)
            if execution_mode == "production":
                blockers.extend("{0}_production_blocker".format(item) for item in stage_warnings)
        output_dir: Path | None = None
        source_tool_version = "synthetic_fixture" if synthetic_fixture else "not_run_parser_only"
        if execution_mode == "production" and not blockers:
            output_dir, run_metadata, run_warnings, run_blockers = _run_fpocket(
                ctx, state_id, state_dir, staged_receptor, binary
            )
            fpocket_executed = fpocket_executed or run_metadata.get("return_code") is not None
            source_tool_version = str(run_metadata.get("fpocket_version", "unknown"))
            fpocket_runs.append(run_metadata)
            warnings.extend(run_warnings)
            blockers.extend(run_blockers)
        elif resolved_output_root is not None:
            found = _find_fpocket_output_dir(resolved_output_root, state_id)
            if found is None:
                warnings.append("parser_only_fpocket_output_missing:{0}".format(state_id))
            else:
                output_dir = _copy_parser_output(ctx, found, state_dir)
        elif execution_mode == "parser_only":
            warnings.append("parser_only_no_fpocket_output_root:{0}".format(state_id))

        state_rows: list[dict[str, Any]] = []
        if output_dir is not None:
            mapping = _mapping_path(ctx, state_id)
            parsed_rows, parse_warnings = parse_fpocket_output(
                ctx,
                state_id,
                state_role,
                staged_receptor,
                output_dir,
                mapping,
                source_tool_version=source_tool_version,
            )
            state_rows.extend(parsed_rows)
            warnings.extend(parse_warnings)
        raw_csv = raw_root / "{0}_fpocket_raw.csv".format(state_id)
        _write_csv(raw_csv, state_rows, RAW_FIELDS, ctx)
        raw_csv_paths.append(raw_csv)
        all_raw_rows.extend(state_rows)

    raw_candidates = normalize_raw_pockets(all_raw_rows)
    merged_families, member_rows = merge_pocket_families(
        raw_candidates,
        centroid_threshold_angstrom=centroid_threshold,
        residue_jaccard_threshold=jaccard_threshold,
    )
    raw_candidates_csv = merged_root / "pocket_candidates_raw.csv"
    merged_csv = merged_root / "pocket_candidates_merged.csv"
    member_csv = merged_root / "pocket_candidate_members.csv"
    tables_dir = ctx.require_within_run_dir(ctx.run_dir / TABLE_ROOT)
    tables_dir.mkdir(parents=True, exist_ok=True)
    tables_merged_csv = tables_dir / "pocket_candidates_merged.csv"
    _write_csv(raw_candidates_csv, raw_candidates, RAW_CANDIDATE_FIELDS, ctx)
    _write_csv(merged_csv, merged_families, MERGED_FIELDS, ctx)
    _write_csv(member_csv, member_rows, MEMBER_FIELDS, ctx)
    _write_csv(tables_merged_csv, merged_families, MERGED_FIELDS, ctx)

    for forbidden in FORBIDDEN_M2_7_OUTPUTS:
        if (ctx.run_dir / forbidden).exists():
            warnings.append("preexisting_m2_7_output_not_modified:{0}".format(forbidden))

    status = _status(blockers, warnings, len(raw_candidates))
    report = FpocketDiscoveryReport(
        run_id=ctx.run_id,
        mode=mode,
        profile=profile,
        execution_mode=execution_mode,
        status=status,
        raw_pocket_count=len(all_raw_rows),
        raw_candidate_count=len(raw_candidates),
        merged_family_count=len(merged_families),
        warnings=sorted(set(warnings)),
        blockers=blockers,
        raw_csv_paths=raw_csv_paths,
        raw_candidate_csv=raw_candidates_csv,
        merged_candidate_csv=merged_csv,
        member_csv=member_csv,
    )
    summary = _write_summary(ctx, report)
    report.summary_report_path = summary
    raw_status = raw_root / "fpocket_discovery_status.json"
    normalization_status = merged_root / "pocket_normalization_status.json"
    manifest = ctx.manifest_dir / "m2_6_fpocket_adapter_manifest.json"
    payload = {
        "run_id": ctx.run_id,
        "created_at": now_iso(),
        "phase": M2_6_PHASE,
        "mode": mode,
        "profile": profile,
        "execution_mode": execution_mode,
        "status": status,
        "states": state_list,
        "include_reference_states": include_reference_states,
        "synthetic_fixture": synthetic_fixture,
        "preconditions": evidence,
        "fpocket_binary": binary,
        "fpocket_executed": fpocket_executed,
        "fpocket_runs": fpocket_runs,
        "p2rank_executed": False,
        "vina_executed": False,
        "pyrosetta_docking_or_relaxation_executed": False,
        "lightdock_executed": False,
        "m2_7_gates_executed": False,
        "accepted_pocket_export_created": False,
        "compound_scoring_executed": False,
        "candidate_nomination_executed": False,
        "raw_pocket_count": len(all_raw_rows),
        "raw_candidate_count": len(raw_candidates),
        "merged_family_count": len(merged_families),
        "merge_policy": {
            "centroid_merge_threshold_angstrom": centroid_threshold,
            "residue_jaccard_merge_threshold": jaccard_threshold,
            "merge_rule": "centroid_distance <= threshold OR residue_jaccard >= threshold",
        },
        "outputs": {
            "raw_csvs": [_repo_path(ctx, path) for path in raw_csv_paths],
            "raw_candidate_csv": _repo_path(ctx, raw_candidates_csv),
            "merged_candidate_csv": _repo_path(ctx, merged_csv),
            "tables_merged_candidate_csv": _repo_path(ctx, tables_merged_csv),
            "member_csv": _repo_path(ctx, member_csv),
            "summary_report": _repo_path(ctx, summary),
            "raw_status_json": _repo_path(ctx, raw_status),
            "normalization_status_json": _repo_path(ctx, normalization_status),
        },
        "warnings": sorted(set(warnings)),
        "blockers": blockers,
        "non_goals_confirmed": NON_GOALS_CONFIRMED,
    }
    write_json(raw_status, payload, ctx)
    write_json(normalization_status, payload, ctx)
    write_json(manifest, payload, ctx)
    report.raw_status_json = raw_status
    report.normalization_status_json = normalization_status
    report.manifest_path = manifest

    qc_rows = [
        {"check": "m2_4_consensus_patch_present", "status": "PASS" if evidence["m2_4_consensus_patch_exists"] else ("FAIL" if execution_mode == "production" else "WARN"), "details": str(evidence["m2_4_consensus_patch_exists"])},
        {"check": "m2_5_atp_reference_present", "status": "PASS" if evidence["m2_5_atp_reference_exists"] else "WARN", "details": str(evidence["m2_5_atp_reference_exists"])},
        {"check": "raw_candidate_status_ungated", "status": "PASS", "details": "M2.6 does not apply M2.7 gates"},
        {"check": "accepted_pocket_export_not_created", "status": "PASS", "details": "accepted_pockets_for_m3.csv not written by M2.6"},
    ]
    _write_csv(ctx.qc_dir / "m2_6_fpocket_adapter_qc.csv", qc_rows, QC_FIELDS, ctx)

    phase_status = "FAIL" if status == "FAIL" else ("WARN" if status == "PASS_WITH_WARNINGS" else "PASS")
    append_phase_status(
        ctx,
        M2_6_PHASE,
        phase_status,
        "M2.6 pocket discovery/normalization {0} with {1} raw candidate(s)".format(status, len(raw_candidates)),
        {"raw_candidate_count": len(raw_candidates), "m2_7_gates_executed": False},
    )
    append_job_status(
        ctx,
        "m2_6_fpocket_adapter_local",
        phase_status,
        details={"message": "fpocket adapter/parser and raw normalization only"},
    )
    return report
