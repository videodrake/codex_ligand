"""M2.1 PPI input generation from M1 canonical outputs.

This module prepares auditable EGFR-MYO1D input packs and method specification
files for future PyRosetta/LightDock work. It deliberately does not import,
submit, or execute any scientific runtime.
"""

from __future__ import annotations

import csv
import json
import shutil
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

from egfr_myo1d.core.logging_utils import append_job_status, append_phase_status
from egfr_myo1d.core.manifest import load_yaml_config, now_iso, write_json
from egfr_myo1d.core.run_context import RunContext
from egfr_myo1d.io.hashing import describe_file
from egfr_myo1d.io.residue_mapping import read_residue_mapping
from egfr_myo1d.myo1d.construct import emit_myo1d_construct_pdb, slice_myo1d_construct
from egfr_myo1d.structure.myo1d_annotation import expand_range
from egfr_myo1d.structure.pdb_parser import PDBParseError, PDBStructure, parse_pdb


M2_1_PHASE = "m2.1-ppi-input-generation"
PACK_ROOT = "prepared/m2_1_ppi_inputs"
METHODS = ("pyrosetta_global_ppi", "lightdock_gso_ppi")
NON_GOALS_CONFIRMED = [
    "no_pyrosetta_docking_or_relaxation",
    "no_lightdock_execution",
    "no_vina_execution",
    "no_fpocket_or_p2rank_runtime",
    "no_compound_docking",
    "no_compound_scoring",
    "no_candidate_nomination",
    "no_scheduler_submission",
    "no_cleanup_deletion",
    "no_mutation_repair",
    "no_receptor_renumbering",
]


@dataclass
class M2PpiPack:
    state_id: str
    role: str
    status: str
    input_files: dict[str, dict[str, Any]] = field(default_factory=dict)
    output_files: dict[str, dict[str, Any]] = field(default_factory=dict)
    method_specs: dict[str, str] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)
    blockers: list[str] = field(default_factory=list)


@dataclass
class M2PpiInputReport:
    run_id: str
    mode: str
    profile: str
    states_requested: list[str]
    status: str
    packs: list[M2PpiPack] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    blockers: list[str] = field(default_factory=list)
    manifest_path: Path | None = None
    qc_csv_path: Path | None = None
    summary_report_path: Path | None = None


def _load_gates(ctx: RunContext) -> dict[str, Any]:
    return load_yaml_config(ctx.repo_root / "fresh" / "configs" / "gates.yaml")


def _load_receptor_states(ctx: RunContext) -> dict[str, Any]:
    return load_yaml_config(
        ctx.repo_root / "fresh" / "configs" / "receptor_states.yaml"
    )


def _primary_states(ctx: RunContext) -> list[str]:
    data = _load_receptor_states(ctx)
    return list(data.get("primary_membrane_validated_states", []))


def _state_role(ctx: RunContext, state_id: str) -> str:
    data = _load_receptor_states(ctx)
    return data.get("state_policy", {}).get(state_id, "unknown")


def _resolve_states(ctx: RunContext, states: Iterable[str] | None) -> list[str]:
    if states:
        requested = [state.strip() for state in states if state and state.strip()]
        if requested:
            return requested
    return _primary_states(ctx)


def _file_entry(ctx: RunContext, path: Path) -> dict[str, Any]:
    entry = {"path": ctx.relative_to_repo(path)}
    entry.update(describe_file(path))
    return entry


def _copy_into_run(ctx: RunContext, source: Path, target: Path) -> Path:
    safe_target = ctx.require_within_run_dir(target)
    safe_target.parent.mkdir(parents=True, exist_ok=True)
    if not source.is_file():
        raise FileNotFoundError(str(source))
    shutil.copy2(str(source), str(safe_target))
    return safe_target


def _write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str], ctx: RunContext) -> None:
    safe_path = ctx.require_within_run_dir(path)
    safe_path.parent.mkdir(parents=True, exist_ok=True)
    with safe_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            encoded = {}
            for key, value in row.items():
                if isinstance(value, (list, dict)):
                    encoded[key] = json.dumps(value, sort_keys=True)
                else:
                    encoded[key] = value
            writer.writerow(encoded)


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def _m1_paths(ctx: RunContext, state_id: str) -> dict[str, Path]:
    return {
        "dockable_receptor": ctx.run_dir
        / "normalized"
        / "receptors"
        / "{0}_dockable_669_1014_explicit_AB.pdb".format(state_id),
        "runtime_offset_receptor": ctx.run_dir
        / "normalized"
        / "receptors"
        / "{0}_runtime_offset_receptor_only.pdb".format(state_id),
        "full_frame_receptor": ctx.run_dir
        / "normalized"
        / "receptors"
        / "{0}_full_frame_explicit_AB.pdb".format(state_id),
        "receptor_mapping": ctx.qc_dir / "{0}_receptor_mapping.csv".format(state_id),
        "receptor_audit": ctx.qc_dir
        / "{0}_receptor_normalization_audit.csv".format(state_id),
        "receptor_manifest": ctx.manifest_dir
        / "{0}_receptor_manifest.json".format(state_id),
        "membrane_frame": ctx.manifest_dir / "membrane_frame.json",
        "myo1d_normalized": ctx.run_dir / "normalized" / "myo1d" / "MYO1D_955_1006.pdb",
        "myo1d_qc": ctx.qc_dir / "myo1d_construct_qc.csv",
        "myo1d_manifest": ctx.manifest_dir / "myo1d_construct_manifest.json",
    }


def _gate_annotation(ctx: RunContext) -> dict[str, Any]:
    gates = _load_gates(ctx)
    myo = gates.get("myo1d", {})
    key = myo.get("key_residues", {})
    watch = myo.get("watch_residues", {})
    sheet8 = sorted(expand_range(str(key.get("sheet8", "961-964"))))
    sheet9 = sorted(expand_range(str(key.get("sheet9", "968-972"))))
    sheet12 = sorted(expand_range(str(key.get("sheet12", "993-997"))))
    n_watch = sorted(expand_range(str(watch.get("n_terminal", "955-957"))))
    c_watch = sorted(expand_range(str(watch.get("c_terminal", "1001-1006"))))
    return {
        "primary_construct_id": "MYO1D_955_1001_primary",
        "primary_range": [955, 1001],
        "comparator_construct_id": "MYO1D_955_1006_comparator",
        "comparator_range": [955, 1006],
        "active_face_residues": sorted(set(sheet8 + sheet9)),
        "sheet8_residues": sheet8,
        "sheet9_residues": sheet9,
        "sheet12_support_residues": sheet12,
        "n_terminal_watch_residues": n_watch,
        "tail_noise_monitor_residues": c_watch,
        "score_bonus_allowed": False,
        "key_residue_bonus_weight": float(myo.get("key_residue_bonus_weight", 0.0)),
        "annotation_use": "qc_evidence_only_no_score_bonus",
    }


def _parse_chain_ids(path: Path) -> tuple[str, ...]:
    try:
        return parse_pdb(path).chain_ids
    except PDBParseError:
        return tuple()


def _mapping_summary(path: Path) -> dict[str, Any]:
    rows = read_residue_mapping(path)
    protomers = sorted({row.get("protomer_id", "") for row in rows if row.get("protomer_id")})
    b_offsets = []
    for row in rows:
        if row.get("protomer_id") != "B":
            continue
        try:
            source = int(row.get("source_resseq", ""))
            runtime = int(row.get("runtime_resseq", ""))
        except ValueError:
            continue
        b_offsets.append(runtime - source)
    return {
        "row_count": len(rows),
        "protomers": protomers,
        "protomer_b_offsets": sorted(set(b_offsets)),
        "runtime_offset_second_protomer_preserved": bool(b_offsets)
        and sorted(set(b_offsets)) == [1000],
    }


def _write_myo1d_pack(
    ctx: RunContext,
    source_myo1d: Path,
    pack_dir: Path,
) -> dict[str, Path]:
    myo_dir = ctx.require_within_run_dir(pack_dir / "myo1d")
    myo_dir.mkdir(parents=True, exist_ok=True)
    source_structure = parse_pdb(source_myo1d)
    primary_structure = slice_myo1d_construct(source_structure, 955, 1001)
    primary_path = myo_dir / "MYO1D_955_1001_primary.pdb"
    comparator_path = myo_dir / "MYO1D_955_1006_comparator.pdb"
    emit_myo1d_construct_pdb(ctx, primary_structure, primary_path)
    _copy_into_run(ctx, source_myo1d, comparator_path)
    return {"primary": primary_path, "comparator": comparator_path}


def _method_spec(
    ctx: RunContext,
    state_id: str,
    method: str,
    pack_dir: Path,
    packed: dict[str, Path],
    inputs: dict[str, Path],
    annotation: dict[str, Any],
    mapping_summary: dict[str, Any],
) -> dict[str, Any]:
    if method == "pyrosetta_global_ppi":
        receptor_key = "runtime_offset_receptor_pack"
        numbering_policy = "use_runtime_offset_receptor_and_mapping_csv"
    else:
        receptor_key = "dockable_receptor_pack"
        numbering_policy = "use_explicit_AB_dockable_receptor_and_mapping_csv"

    return {
        "schema_version": "m2_1_ppi_input_spec_v0_1",
        "state_id": state_id,
        "method": method,
        "execution_mode": "input_spec_only",
        "execution_allowed": False,
        "docking_executed": False,
        "relaxation_executed": False,
        "engine_command": None,
        "receptor": {
            "assembly": "explicit_AB_dimer",
            "pack_path": ctx.relative_to_repo(packed[receptor_key]),
            "m1_source_path": ctx.relative_to_repo(
                inputs[
                    "runtime_offset_receptor"
                    if method == "pyrosetta_global_ppi"
                    else "dockable_receptor"
                ]
            ),
            "chain_ids": ["A", "B"],
            "target_protomers": ["A", "B"],
            "numbering_policy": numbering_policy,
            "mapping_csv": ctx.relative_to_repo(inputs["receptor_mapping"]),
            "mapping_summary": mapping_summary,
            "insertion_codes_preserved": True,
            "mutation_repair_allowed": False,
            "v924r_policy": "warn_only_not_mutated_if_present_in_m1_manifest",
        },
        "partner": {
            "primary_pack_path": ctx.relative_to_repo(packed["myo1d_primary_pack"]),
            "comparator_pack_path": ctx.relative_to_repo(packed["myo1d_comparator_pack"]),
            "primary_construct_id": annotation["primary_construct_id"],
            "primary_role": "production_oriented_ppi_partner",
            "primary_range": annotation["primary_range"],
            "comparator_construct_id": annotation["comparator_construct_id"],
            "comparator_role": "noise_monitor_only",
            "comparator_range": annotation["comparator_range"],
            "terminal_artifact_962_start_allowed": False,
            "annotation": annotation,
        },
        "membrane_context": {
            "membrane_frame_json": ctx.relative_to_repo(inputs["membrane_frame"]),
            "membrane_proximal_policy": "flag_or_block_in_later_pose_qc",
        },
        "non_target_context": {
            "atp_site_overlap_policy": "blocking_non_target_for_ppi_disruptive_objective",
            "compound_or_ligand_inputs_used": False,
        },
        "non_goals_confirmed": NON_GOALS_CONFIRMED,
    }


def _write_report_markdown(ctx: RunContext, report: M2PpiInputReport) -> Path:
    target = ctx.require_within_run_dir(ctx.reports_dir / "m2_1_ppi_input_generation.md")
    lines = [
        "# M2.1 PPI Input Generation",
        "",
        "Status: {0}".format(report.status),
        "",
        "This report records input packs and method specs only. No docking, relaxation, pocket discovery, scheduler submission, compound scoring, or candidate nomination was run.",
        "",
        "## Packs",
    ]
    for pack in report.packs:
        lines.append(
            "- {0}: {1} (warnings={2}, blockers={3})".format(
                pack.state_id, pack.status, len(pack.warnings), len(pack.blockers)
            )
        )
    lines.extend(["", "## Non-goals Confirmed"])
    for item in NON_GOALS_CONFIRMED:
        lines.append("- {0}".format(item))
    target.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return target


def _missing_blockers(paths: dict[str, Path], required_keys: Iterable[str]) -> list[str]:
    blockers = []
    for key in required_keys:
        if not paths[key].is_file():
            blockers.append("missing_m1_artifact:{0}:{1}".format(key, paths[key]))
    return blockers


def _generate_pack_for_state(
    ctx: RunContext,
    state_id: str,
    mode: str,
    profile: str,
    strict: bool,
) -> M2PpiPack:
    role = _state_role(ctx, state_id)
    pack = M2PpiPack(state_id=state_id, role=role, status="PASS")
    if role != "primary_membrane_validated_state":
        message = "state_not_primary_membrane_validated:{0}:{1}".format(state_id, role)
        if strict or profile == "hpc_strict":
            pack.blockers.append(message)
            pack.status = "FAIL"
        else:
            pack.warnings.append(message)
            pack.status = "WARN"
        return pack

    paths = _m1_paths(ctx, state_id)
    required = [
        "dockable_receptor",
        "runtime_offset_receptor",
        "receptor_mapping",
        "membrane_frame",
        "myo1d_normalized",
        "myo1d_qc",
    ]
    pack.blockers.extend(_missing_blockers(paths, required))
    if pack.blockers:
        pack.status = "FAIL"
        return pack

    pack.input_files = {key: _file_entry(ctx, value) for key, value in paths.items()}
    pack_dir = ctx.require_within_run_dir(ctx.run_dir / PACK_ROOT / state_id)
    receptor_dir = ctx.require_within_run_dir(pack_dir / "receptor")
    specs_dir = ctx.require_within_run_dir(pack_dir / "specs")
    receptor_dir.mkdir(parents=True, exist_ok=True)
    specs_dir.mkdir(parents=True, exist_ok=True)

    dockable_pack = _copy_into_run(
        ctx,
        paths["dockable_receptor"],
        receptor_dir / paths["dockable_receptor"].name,
    )
    runtime_pack = _copy_into_run(
        ctx,
        paths["runtime_offset_receptor"],
        receptor_dir / paths["runtime_offset_receptor"].name,
    )
    mapping_pack = _copy_into_run(
        ctx,
        paths["receptor_mapping"],
        receptor_dir / paths["receptor_mapping"].name,
    )
    myo_paths = _write_myo1d_pack(ctx, paths["myo1d_normalized"], pack_dir)

    chain_ids = _parse_chain_ids(dockable_pack)
    if tuple(chain_ids) != ("A", "B"):
        pack.blockers.append("dockable_receptor_not_explicit_AB:{0}".format(chain_ids))

    mapping_summary = _mapping_summary(paths["receptor_mapping"])
    if not mapping_summary["runtime_offset_second_protomer_preserved"]:
        pack.blockers.append("protomer_B_runtime_offset_not_preserved")

    membrane = _read_json(paths["membrane_frame"])
    state_frame = membrane.get("states", {}).get(state_id, {})
    if state_frame.get("status") not in ("PASS", "WARN"):
        pack.warnings.append(
            "membrane_frame_status:{0}".format(state_frame.get("status", "missing"))
        )

    annotation = _gate_annotation(ctx)
    if annotation["key_residue_bonus_weight"] != 0.0:
        pack.blockers.append("key_residue_bonus_weight_not_zero")

    packed = {
        "dockable_receptor_pack": dockable_pack,
        "runtime_offset_receptor_pack": runtime_pack,
        "receptor_mapping_pack": mapping_pack,
        "myo1d_primary_pack": myo_paths["primary"],
        "myo1d_comparator_pack": myo_paths["comparator"],
    }
    inputs = {
        "dockable_receptor": paths["dockable_receptor"],
        "runtime_offset_receptor": paths["runtime_offset_receptor"],
        "receptor_mapping": paths["receptor_mapping"],
        "membrane_frame": paths["membrane_frame"],
    }

    for method in METHODS:
        spec = _method_spec(
            ctx,
            state_id,
            method,
            pack_dir,
            packed,
            inputs,
            annotation,
            mapping_summary,
        )
        spec_path = specs_dir / "{0}_{1}_input_spec.json".format(
            state_id, method.replace("_", "-")
        )
        write_json(spec_path, spec, ctx)
        pack.method_specs[method] = ctx.relative_to_repo(spec_path)

    pack.output_files = {key: _file_entry(ctx, value) for key, value in packed.items()}
    pack.output_files.update(
        {
            "pyrosetta_spec": _file_entry(
                ctx, specs_dir / "{0}_pyrosetta-global-ppi_input_spec.json".format(state_id)
            ),
            "lightdock_spec": _file_entry(
                ctx, specs_dir / "{0}_lightdock-gso-ppi_input_spec.json".format(state_id)
            ),
        }
    )
    pack.status = "FAIL" if pack.blockers else ("WARN" if pack.warnings else "PASS")
    return pack


def _qc_rows(report: M2PpiInputReport) -> list[dict[str, Any]]:
    rows = []
    for pack in report.packs:
        rows.append(
            {
                "state_id": pack.state_id,
                "role": pack.role,
                "status": pack.status,
                "dockable_receptor_pack": pack.output_files.get(
                    "dockable_receptor_pack", {}
                ).get("path", ""),
                "runtime_offset_receptor_pack": pack.output_files.get(
                    "runtime_offset_receptor_pack", {}
                ).get("path", ""),
                "myo1d_primary_pack": pack.output_files.get("myo1d_primary_pack", {}).get(
                    "path", ""
                ),
                "myo1d_comparator_pack": pack.output_files.get(
                    "myo1d_comparator_pack", {}
                ).get("path", ""),
                "score_bonus_allowed": False,
                "external_execution_allowed": False,
                "warnings": pack.warnings,
                "blockers": pack.blockers,
            }
        )
    return rows


def generate_m2_1_ppi_inputs(
    ctx: RunContext,
    mode: str = "smoke_input",
    profile: str = "codex_dev",
    states: Iterable[str] | None = None,
    strict: bool = False,
) -> M2PpiInputReport:
    """Generate M2.1 PPI input packs/specs from M1 canonical outputs.

    The function only copies M1 artifacts, derives MYO1D 955-1001 from the M1
    955-1006 construct, and writes JSON/CSV/Markdown reports under
    ``ctx.run_dir``.
    """
    ctx.create_directories()
    state_list = _resolve_states(ctx, states)
    packs = [
        _generate_pack_for_state(ctx, state, mode, profile, strict)
        for state in state_list
    ]

    blockers = [
        "{0}:{1}".format(pack.state_id, blocker)
        for pack in packs
        for blocker in pack.blockers
    ]
    warnings = [
        "{0}:{1}".format(pack.state_id, warning)
        for pack in packs
        for warning in pack.warnings
    ]
    if blockers:
        status = "FAIL"
    elif warnings or any(pack.status == "WARN" for pack in packs):
        status = "PASS_WITH_WARNINGS"
    else:
        status = "PASS"

    report = M2PpiInputReport(
        run_id=ctx.run_id,
        mode=mode,
        profile=profile,
        states_requested=state_list,
        status=status,
        packs=packs,
        warnings=warnings,
        blockers=blockers,
    )

    qc_path = ctx.qc_dir / "m2_1_ppi_input_qc.csv"
    _write_csv(
        qc_path,
        _qc_rows(report),
        [
            "state_id",
            "role",
            "status",
            "dockable_receptor_pack",
            "runtime_offset_receptor_pack",
            "myo1d_primary_pack",
            "myo1d_comparator_pack",
            "score_bonus_allowed",
            "external_execution_allowed",
            "warnings",
            "blockers",
        ],
        ctx,
    )
    report.qc_csv_path = qc_path

    summary_path = _write_report_markdown(ctx, report)
    report.summary_report_path = summary_path

    payload = {
        "run_id": ctx.run_id,
        "created_at": now_iso(),
        "phase": M2_1_PHASE,
        "mode": mode,
        "profile": profile,
        "states_requested": state_list,
        "status": status,
        "packs": [
            {
                "state_id": pack.state_id,
                "role": pack.role,
                "status": pack.status,
                "inputs": pack.input_files,
                "outputs": pack.output_files,
                "method_specs": pack.method_specs,
                "warnings": pack.warnings,
                "blockers": pack.blockers,
            }
            for pack in packs
        ],
        "execution_allowed": False,
        "docking_executed": False,
        "external_scientific_execution": False,
        "score_bonus_allowed": False,
        "ligand_or_compound_inputs_used": False,
        "non_goals_confirmed": NON_GOALS_CONFIRMED,
        "qc_csv": ctx.relative_to_repo(qc_path),
        "summary_report": ctx.relative_to_repo(summary_path),
    }
    manifest_path = ctx.manifest_dir / "m2_1_ppi_input_manifest.json"
    write_json(manifest_path, payload, ctx)
    report.manifest_path = manifest_path

    phase_status = "FAIL" if status == "FAIL" else (
        "WARN" if status == "PASS_WITH_WARNINGS" else "PASS"
    )
    append_phase_status(
        ctx,
        M2_1_PHASE,
        phase_status,
        "M2.1 PPI input generation {0} for {1} state(s)".format(
            status, len(packs)
        ),
        {
            "states": state_list,
            "manifest": ctx.relative_to_repo(manifest_path),
            "execution_allowed": False,
            "score_bonus_allowed": False,
        },
    )
    append_job_status(
        ctx,
        "m2_1_ppi_input_generation_local",
        phase_status,
        details={
            "message": "local input pack/spec generation only",
            "external_execution_allowed": False,
        },
    )
    return report
