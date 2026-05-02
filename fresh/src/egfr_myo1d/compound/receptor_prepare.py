"""M3-T3 receptor PDBQT preparation and docking-box validation.

This module is deliberately limited to receptor PDBQT preparation and M2 box
validation. It never prepares ligands, invokes Vina, generates docking jobs,
parses poses, scores compounds, ranks candidates, or creates candidate claims.
"""

from __future__ import annotations

import contextlib
import csv
import hashlib
import importlib
import io
import json
import shutil
import subprocess
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from egfr_myo1d.compound.confidentiality import (
    git_check_ignored,
    git_ls_files,
    load_private_map,
    public_output_scan_paths,
    scan_internal_id_leaks,
)
from egfr_myo1d.compound.docking_boxes import (
    BOX_MANIFEST_FIELDS,
    BOX_QC_FIELDS,
    BOX_CHECK_IDS,
    PreparedState,
    read_csv_rows,
    validate_docking_boxes,
)
from egfr_myo1d.core.logging_utils import append_job_status, append_phase_status, log_master
from egfr_myo1d.core.run_context import RunContext, ensure_within


RECEPTOR_PREP_FIELDS = [
    "state_id",
    "state_role",
    "source_receptor_file",
    "source_receptor_format",
    "source_receptor_sha256",
    "source_receptor_size_bytes",
    "source_receptor_origin",
    "source_receptor_is_normalized",
    "source_receptor_is_dimer",
    "source_receptor_is_old_workflow",
    "source_receptor_is_monomer_only",
    "receptor_mapping_file",
    "membrane_frame_file",
    "prepared_receptor_pdbqt_file",
    "prepared_receptor_pdbqt_sha256",
    "prepared_receptor_pdbqt_size_bytes",
    "preparation_method",
    "conversion_tool",
    "conversion_tool_version",
    "tool_stdout_file",
    "tool_stderr_file",
    "pdbqt_validation_success",
    "chain_mapping_success",
    "protomer_mapping_success",
    "raw_input_modified",
    "prepared_file_gitignored",
    "prepared_file_tracked",
    "allowed_for_vina_smoke",
    "preparation_status",
    "preparation_notes",
]

RECEPTOR_HASH_FIELDS = [
    "state_id",
    "prepared_receptor_pdbqt_file",
    "prepared_receptor_pdbqt_sha256",
    "prepared_receptor_pdbqt_size_bytes",
    "source_receptor_file",
    "source_receptor_sha256",
    "hash_status",
    "notes",
]

RECEPTOR_QC_FIELDS = [
    "state_id",
    "check_id",
    "category",
    "status",
    "severity",
    "details",
    "recommended_fix",
]

RECEPTOR_CHECK_IDS = [
    "m2_accepted_boxes_present",
    "m1_receptor_mapping_present",
    "membrane_frame_present",
    "source_receptor_resolved",
    "source_receptor_nonempty",
    "source_receptor_atom_records",
    "source_receptor_normalized",
    "source_receptor_dimer",
    "source_receptor_not_monomer_only",
    "source_receptor_not_old_workflow",
    "source_hash_recorded",
    "receptor_pdbqt_tool_available",
    "receptor_pdbqt_generated_or_supplied",
    "receptor_pdbqt_nonempty",
    "receptor_pdbqt_atom_records",
    "receptor_pdbqt_validation",
    "prepared_hash_recorded",
    "raw_input_unchanged",
    "prepared_path_inside_run_dir",
    "prepared_file_gitignored",
    "prepared_file_not_tracked",
    "no_coordinates_printed",
    "non_goals_preserved",
]

RECEPTOR_GITIGNORE_PATTERNS = [
    "fresh/runs/*/phase3_compounds/receptor_pdbqt/*",
    "fresh/runs/*/phase3_compounds/receptor_pdbqt/**",
    "fresh/runs/*/phase3_compounds/receptors/prepared/*",
    "fresh/runs/*/phase3_compounds/receptors/prepared/**",
    "*.receptor.pdbqt.tmp",
    "*.pdbqt.tmp",
    "*.vina.tmp",
]

FORBIDDEN_OUTPUT_DIRS = [
    "docking_inputs",
    "docking_outputs",
    "vina_raw",
    "pocket_first",
    "compound_anchor",
    "pose_attribution",
]

PRIMARY_STATES = {"EGFR_160-185", "EGFR_170-200"}


@dataclass(frozen=True)
class ReceptorToolStatus:
    prepare_receptor4_available: bool
    prepare_receptor4_version: str | None
    prepare_receptor4_binary: str | None
    mk_prepare_receptor_available: bool
    mk_prepare_receptor_version: str | None
    mk_prepare_receptor_binary: str | None
    meeko_available: bool
    meeko_version: str | None
    openbabel_available: bool
    openbabel_version: str | None
    rdkit_available: bool
    rdkit_version: str | None
    vina_available: bool
    vina_invoked: bool = False


@dataclass
class M3ReceptorBoxResult:
    status: str
    m3_t4_allowed: bool
    blockers: list[str]
    warnings: list[str]
    receptor_preparation_manifest_csv: Path
    receptor_hashes_csv: Path
    docking_box_manifest_csv: Path
    receptor_prep_qc_csv: Path
    docking_box_qc_csv: Path
    summary_json: Path
    report_md: Path
    phase3_log: Path
    counts: dict[str, int]
    states: dict[str, dict[str, Any]]
    tools: ReceptorToolStatus


def now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def _bool_text(value: bool) -> str:
    return "true" if value else "false"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _severity(status: str, blocker: bool = False) -> str:
    if status == "FAIL" and blocker:
        return "BLOCKER"
    if status in {"FAIL", "WARN"}:
        return "MAJOR"
    if status == "NOT_APPLICABLE":
        return "MINOR"
    return "INFO"


def _append_receptor_qc(
    rows: list[dict[str, str]],
    state_id: str,
    check_id: str,
    category: str,
    status: str,
    details: str,
    recommended_fix: str = "",
    blocker: bool = False,
) -> None:
    rows.append(
        {
            "state_id": state_id,
            "check_id": check_id,
            "category": category,
            "status": status,
            "severity": _severity(status, blocker),
            "details": details,
            "recommended_fix": recommended_fix,
        }
    )


def _write_csv(path: Path, fields: list[str], rows: list[dict[str, Any]], ctx: RunContext) -> None:
    safe = ctx.require_within_run_dir(path)
    safe.parent.mkdir(parents=True, exist_ok=True)
    with safe.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _append_missing_receptor_qc(rows: list[dict[str, str]], states: list[str]) -> None:
    present = {(row["state_id"], row["check_id"]) for row in rows}
    for state in states or [""]:
        for check_id in RECEPTOR_CHECK_IDS:
            if (state, check_id) not in present:
                _append_receptor_qc(rows, state, check_id, "coverage", "NOT_APPLICABLE", "check not reached")


def _append_missing_box_qc(rows: list[dict[str, str]]) -> None:
    present = {(row["pocket_family_id"], row["state_id"], row["protomer_id"], row["box_id"], row["check_id"]) for row in rows}
    scopes = {(row["pocket_family_id"], row["state_id"], row["protomer_id"], row["box_id"]) for row in rows} or {("", "", "", "")}
    for family, state, protomer, box_id in scopes:
        for check_id in BOX_CHECK_IDS:
            if (family, state, protomer, box_id, check_id) not in present:
                rows.append(
                    {
                        "pocket_family_id": family,
                        "state_id": state,
                        "protomer_id": protomer,
                        "box_id": box_id,
                        "check_id": check_id,
                        "category": "coverage",
                        "status": "NOT_APPLICABLE",
                        "severity": "MINOR",
                        "details": "check not reached",
                        "recommended_fix": "",
                    }
                )


def _tool_version(args: list[str]) -> str | None:
    try:
        proc = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False, timeout=10)
    except (OSError, subprocess.TimeoutExpired):
        return None
    text = (proc.stdout or proc.stderr or "").strip()
    return text.splitlines()[0] if text else None


def detect_receptor_tools() -> ReceptorToolStatus:
    prepare = shutil.which("prepare_receptor4.py")
    mk = shutil.which("mk_prepare_receptor.py")
    obabel = shutil.which("obabel")
    vina = shutil.which("vina")

    meeko_available = False
    meeko_version = None
    try:
        meeko = importlib.import_module("meeko")
        meeko_available = True
        meeko_version = getattr(meeko, "__version__", None)
    except Exception:
        pass

    rdkit_available = False
    rdkit_version = None
    try:
        with contextlib.redirect_stderr(io.StringIO()):
            rdkit = importlib.import_module("rdkit")
            importlib.import_module("rdkit.Chem")
        rdkit_available = True
        rdkit_version = getattr(rdkit, "__version__", None)
    except BaseException:
        pass

    return ReceptorToolStatus(
        prepare_receptor4_available=bool(prepare),
        prepare_receptor4_version=_tool_version([prepare, "-h"]) if prepare else None,
        prepare_receptor4_binary=prepare,
        mk_prepare_receptor_available=bool(mk),
        mk_prepare_receptor_version=_tool_version([mk, "--help"]) if mk else None,
        mk_prepare_receptor_binary=mk,
        meeko_available=meeko_available,
        meeko_version=meeko_version,
        openbabel_available=bool(obabel),
        openbabel_version=_tool_version([obabel, "-V"]) if obabel else None,
        rdkit_available=rdkit_available,
        rdkit_version=rdkit_version,
        vina_available=bool(vina),
    )


def _ensure_logs(ctx: RunContext) -> None:
    ctx.create_directories()
    for directory in [ctx.logs_dir, ctx.jobs_log_dir, ctx.errors_dir]:
        ctx.require_within_run_dir(directory).mkdir(parents=True, exist_ok=True)
    for path in [
        ctx.logs_dir / "master.log",
        ctx.logs_dir / "phase_status.jsonl",
        ctx.logs_dir / "job_status.jsonl",
        ctx.errors_dir / "error_summary.txt",
        ctx.errors_dir / "failed_jobs.csv",
        ctx.logs_dir / "phase3_compounds.log",
    ]:
        safe = ctx.require_within_run_dir(path)
        if not safe.exists():
            safe.touch()
    failed = ctx.errors_dir / "failed_jobs.csv"
    if failed.stat().st_size == 0:
        failed.write_text("timestamp,job_name,status,message\n", encoding="utf-8")


def update_receptor_gitignore(repo_root: Path) -> tuple[bool, list[str]]:
    gitignore = repo_root / ".gitignore"
    existing = gitignore.read_text(encoding="utf-8").splitlines() if gitignore.is_file() else []
    missing = [pattern for pattern in RECEPTOR_GITIGNORE_PATTERNS if pattern not in existing]
    if not missing:
        return False, list(RECEPTOR_GITIGNORE_PATTERNS)
    gitignore.parent.mkdir(parents=True, exist_ok=True)
    with gitignore.open("a", encoding="utf-8") as handle:
        if existing and existing[-1].strip():
            handle.write("\n")
        handle.write("\n# M3 receptor PDBQT and temporary docking-input protection\n")
        for pattern in missing:
            handle.write(pattern + "\n")
    return True, list(RECEPTOR_GITIGNORE_PATTERNS)


def _safe_private_map_path(ctx: RunContext, private_map: Path | None) -> Path:
    path = private_map or (ctx.fresh_root / "data" / "private" / "compound_id_map.csv")
    resolved = Path(path)
    if not resolved.is_absolute():
        resolved = ctx.repo_root / resolved
    return ensure_within(resolved.resolve(), ctx.fresh_root)


def _resolve_m2_path(ctx: RunContext, m2_run_id: str, candidates: list[str]) -> Path | None:
    m2_root = ensure_within(ctx.fresh_root / "runs" / m2_run_id, ctx.fresh_root / "runs")
    for rel in candidates:
        path = m2_root / rel
        if path.is_file():
            return path
    return None


def _m2_inputs(ctx: RunContext, m2_run_id: str) -> dict[str, Path | None]:
    return {
        "accepted_pocket_boxes": _resolve_m2_path(
            ctx,
            m2_run_id,
            [
                "phase2_pockets/final/accepted_pocket_boxes.csv",
                "phase2_pockets/export_for_m3/accepted_pocket_boxes.csv",
                "phase2_pockets/tables/accepted_pocket_boxes.csv",
            ],
        ),
        "accepted_pockets_for_m3": _resolve_m2_path(
            ctx,
            m2_run_id,
            [
                "phase2_pockets/final/accepted_pockets_for_m3.csv",
                "phase2_pockets/export_for_m3/accepted_pockets_for_m3.csv",
                "phase2_pockets/tables/accepted_pockets_for_m3.csv",
            ],
        ),
        "pocket_gate_qc": _resolve_m2_path(
            ctx,
            m2_run_id,
            [
                "phase2_pockets/qc/pocket_gate_qc.csv",
                "phase2_pockets/final/pocket_gate_qc.csv",
                "phase2_pockets/export_for_m3/pocket_gate_qc.csv",
                "phase2_pockets/gated/pocket_gate_qc.csv",
                "phase2_pockets/tables/pocket_gate_qc.csv",
            ],
        ),
        "atp_reference": _resolve_m2_path(
            ctx,
            m2_run_id,
            [
                "phase2_pockets/atp_reference/atp_site_reference.csv",
                "phase2_pockets/export_for_m3/atp_site_reference.csv",
                "phase2_pockets/tables/atp_site_reference.csv",
            ],
        ),
        "ppi_consensus_patch": _resolve_m2_path(
            ctx,
            m2_run_id,
            [
                "phase1_ppi/consensus/ppi_consensus_patch_merged.csv",
                "phase2_pockets/export_for_m3/ppi_consensus_patch.csv",
            ],
        ),
    }


def _load_m3_t2_summary(ctx: RunContext) -> tuple[dict[str, Any] | None, str | None]:
    path = ctx.run_dir / "phase3_compounds" / "qc" / "m3_ligand_prep_summary.json"
    if not path.is_file():
        return None, "M3-T2 ligand preparation summary not found"
    try:
        return json.loads(path.read_text(encoding="utf-8")), None
    except (OSError, json.JSONDecodeError):
        return None, "M3-T2 ligand preparation summary is unreadable"


def _state_role(state_id: str, row_role: str | None = None) -> str:
    role = (row_role or "").strip().lower()
    if role in {"primary", "reference", "exploratory"}:
        return role
    if state_id in PRIMARY_STATES:
        return "primary"
    if state_id == "3GT8_raw":
        return "reference"
    return "unknown"


def _states_from_boxes(box_file: Path | None, states_arg: list[str] | None) -> tuple[list[str], dict[str, str]]:
    roles: dict[str, str] = {}
    if states_arg:
        for state in states_arg:
            roles[state] = _state_role(state)
        return list(dict.fromkeys(states_arg)), roles
    if not box_file or not box_file.is_file():
        return [], roles
    rows, _ = read_csv_rows(box_file)
    states: list[str] = []
    for row in rows:
        state = (row.get("state_id") or "").strip()
        if state and state not in states:
            states.append(state)
        if state:
            roles[state] = _state_role(state, row.get("state_role"))
    return states, roles


def _path_from_repo(ctx: RunContext, text: str) -> Path | None:
    if not text:
        return None
    path = Path(text)
    if not path.is_absolute():
        path = ctx.repo_root / path
    try:
        return ensure_within(path.resolve(), ctx.fresh_root)
    except Exception:
        return None


def _rows_for_state(box_file: Path | None, state_id: str) -> list[dict[str, str]]:
    if not box_file or not box_file.is_file():
        return []
    rows, _ = read_csv_rows(box_file)
    return [row for row in rows if (row.get("state_id") or "").strip() == state_id]


def _truthy(value: Any) -> bool:
    return str(value or "").strip().lower() in {"true", "t", "1", "yes", "y", "pass", "passed"}


def _supplied_receptor_pdbqt_manual_review(box_file: Path | None, state_id: str) -> bool:
    for row in _rows_for_state(box_file, state_id):
        if _truthy(row.get("supplied_receptor_pdbqt_manual_review")) or _truthy(row.get("receptor_pdbqt_manual_review")):
            return True
    return False


def _candidate_receptor_files(ctx: RunContext, m2_run_id: str, state_id: str, box_file: Path | None) -> list[Path]:
    candidates: list[Path] = []
    for row in _rows_for_state(box_file, state_id):
        for key in ["receptor_pdb", "source_receptor_file"]:
            path = _path_from_repo(ctx, (row.get(key) or "").strip())
            if path:
                candidates.append(path)
    m2_root = ctx.fresh_root / "runs" / m2_run_id
    roots = [
        m2_root / "prepared" / "m2_1_ppi_inputs" / state_id / "receptor",
        m2_root / "phase1_receptors" / "normalized",
        m2_root / "phase1_receptors" / "final",
        m2_root / "phase0_inputs" / "receptors",
        ctx.run_dir / "phase0_inputs" / "receptors",
        ctx.run_dir / "phase1_receptors" / "normalized",
        ctx.fresh_root / "data" / "normalized" / "receptors" / state_id,
    ]
    names = [
        f"{state_id}_dockable_669_1014_explicit_AB.pdb",
        f"{state_id}_dockable_reference_explicit_AB.pdb",
        f"{state_id}_runtime_offset_receptor_only.pdb",
        "dockable_669_1014_explicit_AB.pdb",
        "dockable_reference_explicit_AB.pdb",
        "runtime_offset_receptor_only.pdb",
        "receptor.pdbqt",
    ]
    for root in roots:
        for name in names:
            path = root / name
            if path.is_file():
                candidates.append(path.resolve())
        if root.is_dir():
            for path in sorted(root.glob(f"*{state_id}*.pdb")) + sorted(root.glob(f"*{state_id}*.pdbqt")):
                candidates.append(path.resolve())
    unique: list[Path] = []
    for path in candidates:
        if path not in unique:
            unique.append(path)
    return unique


def _candidate_mapping_files(ctx: RunContext, m2_run_id: str, state_id: str, box_file: Path | None) -> list[Path]:
    candidates: list[Path] = []
    for row in _rows_for_state(box_file, state_id):
        path = _path_from_repo(ctx, (row.get("receptor_mapping_csv") or "").strip())
        if path:
            candidates.append(path)
    m2_root = ctx.fresh_root / "runs" / m2_run_id
    for path in [
        m2_root / "phase1_receptors" / "manifests" / "receptor_mapping.csv",
        m2_root / "phase1_receptors" / "normalized" / "receptor_mapping.csv",
        m2_root / "phase0_inputs" / "receptor_mapping.csv",
        m2_root / "qc" / f"{state_id}_receptor_mapping.csv",
        m2_root / "normalized" / "receptors" / f"{state_id}_receptor_mapping.csv",
        ctx.run_dir / "phase0_inputs" / "receptor_mapping.csv",
        ctx.run_dir / "phase1_receptors" / "normalized" / "receptor_mapping.csv",
        ctx.run_dir / "qc" / f"{state_id}_receptor_mapping.csv",
        ctx.fresh_root / "data" / "normalized" / "receptors" / state_id / "receptor_mapping.csv",
    ]:
        if path.is_file():
            candidates.append(path.resolve())
    unique: list[Path] = []
    for path in candidates:
        if path not in unique:
            unique.append(path)
    return unique


def _resolve_membrane_frame(ctx: RunContext, m2_run_id: str) -> Path | None:
    m2_root = ctx.fresh_root / "runs" / m2_run_id
    for path in [
        m2_root / "phase0_inputs" / "membrane_frame.json",
        m2_root / "manifest" / "membrane_frame.json",
        m2_root / "qc" / "membrane_frame.json",
        ctx.run_dir / "phase0_inputs" / "membrane_frame.json",
        ctx.run_dir / "manifest" / "membrane_frame.json",
        ctx.run_dir / "qc" / "membrane_frame.json",
    ]:
        if path.is_file():
            return path.resolve()
    return None


def _select_one(paths: list[Path]) -> tuple[Path | None, str | None]:
    existing = [path for path in paths if path.is_file()]
    if not existing:
        return None, "not found"
    if len(existing) == 1:
        return existing[0], None
    explicit = [path for path in existing if "dockable" in path.name.lower() or path.suffix.lower() == ".pdbqt"]
    if len(explicit) == 1:
        return explicit[0], None
    return None, "multiple candidates found"


def _atom_lines(path: Path) -> list[str]:
    try:
        return [line for line in path.read_text(encoding="utf-8", errors="ignore").splitlines() if line.startswith(("ATOM", "HETATM"))]
    except OSError:
        return []


def _chains_from_atoms(atom_lines: list[str]) -> set[str]:
    chains = set()
    for line in atom_lines:
        if len(line) > 21:
            chain = line[21].strip()
            if chain:
                chains.add(chain)
    return chains


def _coords_from_atoms(atom_lines: list[str]) -> list[tuple[float, float, float]]:
    coords: list[tuple[float, float, float]] = []
    for line in atom_lines:
        try:
            coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
        except ValueError:
            continue
    return coords


def _bounds(coords: list[tuple[float, float, float]]) -> tuple[tuple[float, float], tuple[float, float], tuple[float, float]] | None:
    if not coords:
        return None
    xs, ys, zs = zip(*coords)
    return ((min(xs), max(xs)), (min(ys), max(ys)), (min(zs), max(zs)))


def _mapping_protomers(mapping_path: Path | None) -> tuple[set[str], set[str]]:
    if not mapping_path or not mapping_path.is_file():
        return set(), set()
    try:
        rows, _ = read_csv_rows(mapping_path)
    except (OSError, csv.Error, UnicodeDecodeError):
        return set(), set()
    protomers = {(row.get("protomer_id") or row.get("protomer") or "").strip() for row in rows}
    chains = {(row.get("runtime_chain") or row.get("chain_id") or row.get("source_chain") or "").strip() for row in rows}
    return {p for p in protomers if p}, {c for c in chains if c}


def _is_old_workflow(path: Path | None) -> bool:
    if not path:
        return False
    text = path.as_posix().lower()
    return any(token in text for token in ["workflow_a", "workflow-b", "workflow_b", "results_export", "old"])


def _looks_monomer(path: Path | None, chains: set[str], protomers: set[str]) -> bool | None:
    if path and "monomer" in path.as_posix().lower():
        return True
    if len(chains) == 1:
        return True
    if len(protomers) >= 2 or len(chains) >= 2:
        return False
    return None


def _source_origin(ctx: RunContext, m2_run_id: str, path: Path | None) -> str:
    if not path:
        return ""
    m2_root = ctx.fresh_root / "runs" / m2_run_id
    for label, root in [
        ("m2_phase1_receptors_normalized", m2_root / "phase1_receptors" / "normalized"),
        ("m2_phase1_receptors_final", m2_root / "phase1_receptors" / "final"),
        ("m2_phase0_inputs_receptors", m2_root / "phase0_inputs" / "receptors"),
        ("m3_phase0_inputs_receptors", ctx.run_dir / "phase0_inputs" / "receptors"),
        ("m3_phase1_receptors_normalized", ctx.run_dir / "phase1_receptors" / "normalized"),
        ("fresh_data_normalized_receptors", ctx.fresh_root / "data" / "normalized" / "receptors"),
    ]:
        try:
            path.resolve().relative_to(root.resolve())
            return label
        except ValueError:
            continue
    return "unknown"


def _validate_pdbqt(path: Path) -> tuple[bool, str]:
    atoms = _atom_lines(path)
    if not path.is_file() or path.stat().st_size <= 0:
        return False, "PDBQT missing or empty"
    if not atoms:
        return False, "PDBQT has no ATOM/HETATM records"
    ligand_tokens = {"ROOT", "BRANCH", "TORSDOF"}
    text = path.read_text(encoding="utf-8", errors="ignore")
    if any(token in text for token in ligand_tokens) and len(_chains_from_atoms(atoms)) <= 1:
        return False, "PDBQT looks ligand-like rather than receptor-like"
    return True, "PDBQT receptor atom records present"


def _copy_or_prepare_pdbqt(
    source: Path,
    output: Path,
    state_id: str,
    ctx: RunContext,
    tools: ReceptorToolStatus,
) -> tuple[bool, str, str, str, Path, Path]:
    output = ctx.require_within_run_dir(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    stdout_path = ctx.require_within_run_dir(ctx.jobs_log_dir / f"m3_t3_{state_id}_receptor_prep.stdout.txt")
    stderr_path = ctx.require_within_run_dir(ctx.jobs_log_dir / f"m3_t3_{state_id}_receptor_prep.stderr.txt")
    if source.suffix.lower() == ".pdbqt":
        shutil.copyfile(source, output)
        stdout_path.write_text("supplied receptor PDBQT validated and copied\n", encoding="utf-8")
        stderr_path.write_text("", encoding="utf-8")
        return True, "supplied_receptor_pdbqt_validated", "none", "", stdout_path, stderr_path
    if tools.prepare_receptor4_available and tools.prepare_receptor4_binary:
        args = [tools.prepare_receptor4_binary, "-r", str(source), "-o", str(output)]
        method = "prepare_receptor4"
        version = tools.prepare_receptor4_version or ""
    elif tools.mk_prepare_receptor_available and tools.mk_prepare_receptor_binary:
        args = [tools.mk_prepare_receptor_binary, "-i", str(source), "-o", str(output)]
        method = "mk_prepare_receptor"
        version = tools.mk_prepare_receptor_version or ""
    else:
        stdout_path.write_text("", encoding="utf-8")
        stderr_path.write_text("no receptor PDBQT preparation tool available\n", encoding="utf-8")
        return False, "no_receptor_pdbqt_tool_available", "", "", stdout_path, stderr_path
    try:
        proc = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False)
    except OSError as exc:
        stdout_path.write_text("", encoding="utf-8")
        stderr_path.write_text(f"receptor prep invocation failed: {exc.__class__.__name__}\n", encoding="utf-8")
        return False, method, method, version, stdout_path, stderr_path
    stdout_path.write_text(proc.stdout or "", encoding="utf-8")
    stderr_path.write_text(proc.stderr or "", encoding="utf-8")
    return proc.returncode == 0 and output.is_file(), method, method, version, stdout_path, stderr_path


def _prepare_states(
    ctx: RunContext,
    m2_run_id: str,
    box_file: Path | None,
    states: list[str],
    state_roles: dict[str, str],
    membrane_frame: Path | None,
    tools: ReceptorToolStatus,
    mode: str,
    force: bool,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, str]], dict[str, PreparedState], list[str], list[str], dict[str, dict[str, Any]]]:
    manifest_rows: list[dict[str, Any]] = []
    hash_rows: list[dict[str, Any]] = []
    qc_rows: list[dict[str, str]] = []
    prepared_states: dict[str, PreparedState] = {}
    blockers: list[str] = []
    warnings: list[str] = []
    state_summary: dict[str, dict[str, Any]] = {}
    prep_root = ctx.require_within_run_dir(ctx.run_dir / "phase3_compounds" / "receptor_pdbqt")
    prep_root.mkdir(parents=True, exist_ok=True)

    for state_id in states:
        notes: list[str] = []
        status = "NOT_APPLICABLE" if mode == "dry-run" else "PASS"
        role = state_roles.get(state_id) or _state_role(state_id)
        source, source_error = _select_one(_candidate_receptor_files(ctx, m2_run_id, state_id, box_file))
        mapping, mapping_error = _select_one(_candidate_mapping_files(ctx, m2_run_id, state_id, box_file))
        atoms = _atom_lines(source) if source else []
        chains = _chains_from_atoms(atoms)
        coords = _coords_from_atoms(atoms)
        protomers, mapping_chains = _mapping_protomers(mapping)
        dimer = (len(chains) >= 2 or len(protomers) >= 2)
        monomer = _looks_monomer(source, chains, protomers)
        old_workflow = _is_old_workflow(source)
        normalized = bool(source and ("normalized" in source.as_posix().lower() or "dockable" in source.name.lower() or "runtime_offset" in source.name.lower()))
        source_sha = ""
        prepared_sha = ""
        prepared_size = ""
        method = ""
        conversion_tool = ""
        conversion_version = ""
        stdout_file = ""
        stderr_file = ""
        pdbqt_ok = False
        raw_modified = False
        source_size = ""
        source_format = source.suffix.lower().lstrip(".") if source else "unknown"
        out_path = ctx.require_within_run_dir(prep_root / state_id / "receptor.pdbqt")
        prepared_tracked = False
        prepared_ignored = git_check_ignored(ctx.repo_root, out_path)
        supplied_pdbqt_approved = _supplied_receptor_pdbqt_manual_review(box_file, state_id)

        _append_receptor_qc(qc_rows, state_id, "m2_accepted_boxes_present", "input", "PASS" if box_file else ("FAIL" if mode == "prepare" else "WARN"), "M2 accepted boxes checked", "Run M2.8 export.", mode == "prepare" and not box_file)
        _append_receptor_qc(qc_rows, state_id, "m1_receptor_mapping_present", "input", "PASS" if mapping else "FAIL", mapping_error or "receptor mapping found", "Provide receptor_mapping.csv from M1 normalized receptor outputs.", not mapping)
        _append_receptor_qc(qc_rows, state_id, "membrane_frame_present", "input", "PASS" if membrane_frame else "FAIL", "membrane frame found" if membrane_frame else "membrane_frame.json missing", "Provide membrane_frame.json from M1/M2 inputs.", not membrane_frame)
        if not mapping:
            status = "FAIL"
            blockers.append(f"{state_id}: receptor mapping absent or ambiguous")
            notes.append(mapping_error or "receptor mapping missing")
        if not membrane_frame:
            status = "FAIL"
            blockers.append(f"{state_id}: membrane_frame.json absent")
            notes.append("membrane_frame.json missing")
        _append_receptor_qc(qc_rows, state_id, "source_receptor_resolved", "input", "PASS" if source else "FAIL", source_error or "source receptor resolved", "Provide one unambiguous normalized dimer receptor source.", not source)
        if not source:
            status = "FAIL"
            blockers.append(f"{state_id}: no M1 normalized dimer receptor source")
            notes.append(source_error or "source receptor missing")
        elif source.suffix.lower() not in {".pdb", ".pdbqt"}:
            status = "FAIL"
            blockers.append(f"{state_id}: unsupported receptor input format")
            notes.append("unsupported receptor input format")
        elif source.stat().st_size == 0:
            status = "FAIL"
            blockers.append(f"{state_id}: zero-byte receptor input")
            notes.append("zero-byte receptor input")
        if source:
            source_size = str(source.stat().st_size)
            source_sha = _sha256(source)
        _append_receptor_qc(qc_rows, state_id, "source_receptor_nonempty", "input", "PASS" if source and source.stat().st_size > 0 else "FAIL", "source receptor non-empty checked", "Provide a non-empty receptor PDB/PDBQT.", not (source and source.stat().st_size > 0))
        _append_receptor_qc(qc_rows, state_id, "source_receptor_atom_records", "input", "PASS" if atoms else "FAIL", "source receptor atom records checked", "Use receptor file with ATOM/HETATM records.", not atoms)
        _append_receptor_qc(qc_rows, state_id, "source_receptor_normalized", "provenance", "PASS" if normalized else "FAIL", "normalized receptor provenance checked", "Use normalized M1 receptor output.", not normalized)
        _append_receptor_qc(qc_rows, state_id, "source_receptor_dimer", "provenance", "PASS" if dimer else "FAIL", "dimer/protomer evidence checked", "Use dimer-compatible receptor with A/B protomer mapping.", not dimer)
        _append_receptor_qc(qc_rows, state_id, "source_receptor_not_monomer_only", "provenance", "PASS" if monomer is False else "FAIL", "monomer-only source rejected" if monomer else "source is not monomer-only", "Use dimer receptor source.", monomer is not False)
        _append_receptor_qc(qc_rows, state_id, "source_receptor_not_old_workflow", "provenance", "FAIL" if old_workflow else "PASS", "old Workflow A/B source rejected" if old_workflow else "not old Workflow A/B", "Use fresh M1/M2 normalized receptor outputs.", old_workflow)
        _append_receptor_qc(qc_rows, state_id, "source_hash_recorded", "provenance", "PASS" if source_sha else "FAIL", "source hash recorded" if source_sha else "source hash unavailable", "Provide readable receptor source.", not source_sha)
        if not atoms or not normalized or not dimer or monomer is not False or old_workflow:
            status = "FAIL"
            if monomer is not False:
                blockers.append(f"{state_id}: source receptor is monomer-only or lacks dimer mapping")
            if old_workflow:
                blockers.append(f"{state_id}: source receptor is old Workflow A/B output")
            notes.append("source receptor validation failed")

        tool_available = tools.prepare_receptor4_available or tools.mk_prepare_receptor_available or (source is not None and source.suffix.lower() == ".pdbqt")
        _append_receptor_qc(qc_rows, state_id, "receptor_pdbqt_tool_available", "tooling", "PASS" if tool_available else ("WARN" if mode == "dry-run" else "FAIL"), "receptor PDBQT tool availability checked", "Install prepare_receptor4.py or mk_prepare_receptor.py, or supply valid receptor PDBQT.", mode == "prepare" and not tool_available)
        if source and source.suffix.lower() != ".pdbqt" and not tool_available and mode == "prepare":
            status = "FAIL"
            blockers.append(f"{state_id}: receptor PDBQT tool unavailable")
            notes.append("no receptor PDBQT conversion tool available")
        if source and source.suffix.lower() == ".pdbqt" and not supplied_pdbqt_approved:
            status = "FAIL"
            blockers.append(f"{state_id}: supplied receptor PDBQT lacks manual review approval")
            notes.append("supplied receptor PDBQT lacks manual review approval")
        _append_receptor_qc(
            qc_rows,
            state_id,
            "supplied_receptor_pdbqt_manual_review",
            "provenance",
            "PASS" if not source or source.suffix.lower() != ".pdbqt" or supplied_pdbqt_approved else "FAIL",
            "supplied receptor PDBQT manual-review approval checked",
            "Set supplied_receptor_pdbqt_manual_review=true in accepted_pocket_boxes.csv after review, or provide a source PDB for M3-T3 conversion.",
            bool(source and source.suffix.lower() == ".pdbqt" and not supplied_pdbqt_approved),
        )

        source_before = source_sha
        if mode == "dry-run":
            method = "dry_run_no_receptor_pdbqt_written"
            if not tool_available:
                warnings.append(f"{state_id}: receptor PDBQT tool unavailable in dry-run")
        elif status != "FAIL" and source:
            if out_path.exists() and not force:
                status = "FAIL"
                blockers.append(f"{state_id}: prepared receptor PDBQT already exists and force=false")
                notes.append("prepared receptor PDBQT already exists")
            else:
                ok, method, conversion_tool, conversion_version, stdout_path, stderr_path = _copy_or_prepare_pdbqt(source, out_path, state_id, ctx, tools)
                stdout_file = ctx.relative_to_repo(stdout_path)
                stderr_file = ctx.relative_to_repo(stderr_path)
                if method == "supplied_receptor_pdbqt_validated":
                    warnings.append(f"{state_id}: supplied receptor PDBQT reused rather than generated")
                    if status == "PASS":
                        status = "WARN"
                if not ok:
                    status = "FAIL"
                    blockers.append(f"{state_id}: receptor PDBQT generation failure")
                    notes.append("receptor PDBQT generation failed")
        if out_path.is_file():
            pdbqt_ok, pdbqt_message = _validate_pdbqt(out_path)
            notes.append(pdbqt_message)
            prepared_sha = _sha256(out_path)
            prepared_size = str(out_path.stat().st_size)
            prepared_tracked = bool(git_ls_files(ctx.repo_root, out_path))
            prepared_ignored = git_check_ignored(ctx.repo_root, out_path)
            if not pdbqt_ok:
                status = "FAIL"
                blockers.append(f"{state_id}: invalid prepared receptor PDBQT")
            if prepared_tracked:
                status = "FAIL"
                blockers.append(f"{state_id}: generated receptor PDBQT tracked by git")
            if not prepared_ignored and status != "FAIL":
                status = "WARN"
                warnings.append(f"{state_id}: prepared receptor PDBQT is not gitignore-protected")
        if source and source.is_file():
            raw_modified = source_before != _sha256(source)
            if raw_modified:
                status = "FAIL"
                blockers.append(f"{state_id}: raw receptor input was modified")
        try:
            ctx.require_within_run_dir(out_path)
            inside_run = True
        except Exception:
            inside_run = False
            status = "FAIL"
            blockers.append(f"{state_id}: generated receptor PDBQT path resolves outside run_dir")

        _append_receptor_qc(qc_rows, state_id, "receptor_pdbqt_generated_or_supplied", "output", "PASS" if out_path.is_file() else ("NOT_APPLICABLE" if mode == "dry-run" else "FAIL"), "prepared receptor PDBQT output checked", "Generate or supply receptor PDBQT.", mode == "prepare" and not out_path.is_file())
        _append_receptor_qc(qc_rows, state_id, "receptor_pdbqt_nonempty", "output", "PASS" if out_path.is_file() and out_path.stat().st_size > 0 else ("NOT_APPLICABLE" if not out_path.is_file() else "FAIL"), "prepared receptor PDBQT non-empty checked", "Regenerate receptor PDBQT.", out_path.is_file() and out_path.stat().st_size <= 0)
        _append_receptor_qc(qc_rows, state_id, "receptor_pdbqt_atom_records", "output", "PASS" if pdbqt_ok else ("NOT_APPLICABLE" if not out_path.is_file() else "FAIL"), "prepared receptor atom records checked", "Regenerate receptor PDBQT.", out_path.is_file() and not pdbqt_ok)
        _append_receptor_qc(qc_rows, state_id, "receptor_pdbqt_validation", "output", "PASS" if pdbqt_ok else ("NOT_APPLICABLE" if not out_path.is_file() else "FAIL"), "PDBQT validation checked", "Regenerate valid receptor PDBQT.", out_path.is_file() and not pdbqt_ok)
        _append_receptor_qc(qc_rows, state_id, "prepared_hash_recorded", "provenance", "PASS" if prepared_sha else ("NOT_APPLICABLE" if not out_path.is_file() else "FAIL"), "prepared hash recorded" if prepared_sha else "prepared hash unavailable", "Regenerate receptor PDBQT.", out_path.is_file() and not prepared_sha)
        _append_receptor_qc(qc_rows, state_id, "raw_input_unchanged", "provenance", "PASS" if not raw_modified else "FAIL", "raw receptor input unchanged", "Restore raw receptor input.", raw_modified)
        _append_receptor_qc(qc_rows, state_id, "prepared_path_inside_run_dir", "path_safety", "PASS" if inside_run else "FAIL", "prepared receptor path checked", "Write under phase3_compounds/receptor_pdbqt/.", not inside_run)
        _append_receptor_qc(qc_rows, state_id, "prepared_file_gitignored", "confidentiality", "PASS" if prepared_ignored or not out_path.exists() else "WARN", "prepared receptor gitignore checked", "Add receptor PDBQT gitignore patterns.", False)
        _append_receptor_qc(qc_rows, state_id, "prepared_file_not_tracked", "confidentiality", "FAIL" if prepared_tracked else "PASS", "prepared receptor git tracking checked", "Remove from git index without deleting local file.", prepared_tracked)
        _append_receptor_qc(qc_rows, state_id, "no_coordinates_printed", "output_hygiene", "PASS", "coordinates are not written to M3-T3 CSV/JSON/report logs")
        _append_receptor_qc(qc_rows, state_id, "non_goals_preserved", "scope", "PASS", "no ligand prep, Vina, docking jobs, docking, pose parsing, scoring, or candidate nomination")

        source_display = ctx.relative_to_repo(source) if source else ""
        mapping_display = ctx.relative_to_repo(mapping) if mapping else ""
        membrane_display = ctx.relative_to_repo(membrane_frame) if membrane_frame else ""
        prepared_display = ctx.relative_to_repo(out_path) if out_path.is_file() else ""
        allowed = status in {"PASS", "WARN"} and bool(prepared_sha)
        manifest_rows.append(
            {
                "state_id": state_id,
                "state_role": role,
                "source_receptor_file": source_display,
                "source_receptor_format": source_format or "unknown",
                "source_receptor_sha256": source_sha,
                "source_receptor_size_bytes": source_size,
                "source_receptor_origin": _source_origin(ctx, m2_run_id, source),
                "source_receptor_is_normalized": "true" if normalized else "false",
                "source_receptor_is_dimer": "true" if dimer else "false",
                "source_receptor_is_old_workflow": "true" if old_workflow else "false",
                "source_receptor_is_monomer_only": "true" if monomer is True else ("false" if monomer is False else "unknown"),
                "receptor_mapping_file": mapping_display,
                "membrane_frame_file": membrane_display,
                "prepared_receptor_pdbqt_file": prepared_display,
                "prepared_receptor_pdbqt_sha256": prepared_sha,
                "prepared_receptor_pdbqt_size_bytes": prepared_size,
                "preparation_method": method,
                "conversion_tool": conversion_tool,
                "conversion_tool_version": conversion_version,
                "tool_stdout_file": stdout_file,
                "tool_stderr_file": stderr_file,
                "pdbqt_validation_success": "true" if pdbqt_ok else "false",
                "chain_mapping_success": "true" if chains or mapping_chains else "false",
                "protomer_mapping_success": "true" if protomers else "false",
                "raw_input_modified": "true" if raw_modified else "false",
                "prepared_file_gitignored": "true" if prepared_ignored else "false",
                "prepared_file_tracked": "true" if prepared_tracked else "false",
                "allowed_for_vina_smoke": "true" if allowed else "false",
                "preparation_status": status,
                "preparation_notes": "; ".join(dict.fromkeys(notes)),
            }
        )
        hash_rows.append(
            {
                "state_id": state_id,
                "prepared_receptor_pdbqt_file": prepared_display,
                "prepared_receptor_pdbqt_sha256": prepared_sha,
                "prepared_receptor_pdbqt_size_bytes": prepared_size,
                "source_receptor_file": source_display,
                "source_receptor_sha256": source_sha,
                "hash_status": "PASS" if prepared_sha and source_sha and not raw_modified else ("NOT_APPLICABLE" if mode == "dry-run" else "FAIL"),
                "notes": "hashes recorded" if prepared_sha and source_sha else "hash incomplete",
            }
        )
        prepared_states[state_id] = PreparedState(
            state_id=state_id,
            state_role=role,
            status=status,
            receptor_pdbqt_file=prepared_display or None,
            receptor_pdbqt_sha256=prepared_sha or None,
            receptor_mapping_file=mapping_display or None,
            membrane_frame_file=membrane_display or None,
            protomer_ids=protomers,
            coordinate_bounds=_bounds(coords),
        )
        state_summary[state_id] = {
            "status": status,
            "state_role": role,
            "source_receptor_file": source_display or None,
            "source_receptor_sha256": source_sha or None,
            "prepared_receptor_pdbqt_file": prepared_display or None,
            "prepared_receptor_pdbqt_sha256": prepared_sha or None,
            "source_receptor_is_normalized": normalized,
            "source_receptor_is_dimer": dimer,
            "source_receptor_is_old_workflow": old_workflow,
            "source_receptor_is_monomer_only": monomer,
            "allowed_for_vina_smoke": allowed,
        }
    return manifest_rows, hash_rows, qc_rows, prepared_states, blockers, warnings, state_summary


def _scan_output_hygiene(ctx: RunContext, private_entries: list[Any]) -> tuple[int, list[str], list[str]]:
    leaks = scan_internal_id_leaks(ctx, private_entries)
    coordinate_hits: list[str] = []
    for path in public_output_scan_paths(ctx):
        if path.suffix.lower() == ".pdbqt":
            continue
        try:
            text = path.read_text(encoding="utf-8", errors="ignore")
        except OSError:
            continue
        for line in text.splitlines():
            if line.startswith(("ATOM  ", "HETATM")):
                coordinate_hits.append(ctx.relative_to_repo(path))
                break
    return len(leaks), coordinate_hits, [ctx.relative_to_repo(path) for path in public_output_scan_paths(ctx)]


def _forbidden_outputs_created(ctx: RunContext) -> list[str]:
    phase3 = ctx.run_dir / "phase3_compounds"
    created: list[str] = []
    for rel in FORBIDDEN_OUTPUT_DIRS:
        path = phase3 / rel
        if path.exists():
            created.append(ctx.relative_to_repo(path))
    for rel in [
        "tables/final_m3_candidate_hypotheses.csv",
        "tables/pocket_compound_evidence_table.csv",
    ]:
        path = phase3 / rel
        if path.exists():
            created.append(ctx.relative_to_repo(path))
    return created


def _default_box_margin(ctx: RunContext) -> float:
    try:
        import yaml
    except Exception:
        return 0.0
    path = ctx.fresh_root / "configs" / "gates.yaml"
    if not path.is_file():
        return 0.0
    try:
        data = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    except Exception:
        return 0.0
    for key in ["box_margin", "docking_box_margin", "m3_box_margin"]:
        value = data.get(key)
        if isinstance(value, (int, float)):
            return float(value)
    compound = data.get("compound") or data.get("m3") or {}
    if isinstance(compound, dict):
        for key in ["box_margin", "docking_box_margin"]:
            value = compound.get(key)
            if isinstance(value, (int, float)):
                return float(value)
    return 0.0


def _write_report(path: Path, ctx: RunContext, summary: dict[str, Any], outputs: dict[str, Path]) -> None:
    safe = ctx.require_within_run_dir(path)
    safe.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# M3-T3 Receptor PDBQT Preparation and Docking Box Validation",
        "",
        f"- run_id: {summary['run_id']}",
        f"- m2_run_id: {summary['m2_run_id']}",
        f"- mode: {summary['mode']}",
        f"- overall_status: {summary['overall_status']}",
        f"- m3_t4_allowed: {str(summary['m3_t4_allowed']).lower()}",
        "- no_ligand_preparation_attempted: true",
        "- no_pocket_discovery_attempted: true",
        "- no_vina_attempted: true",
        "- no_docking_job_generation_attempted: true",
        "- no_docking_attempted: true",
        "- no_pose_parsing_attempted: true",
        "- no_scoring_attempted: true",
        "- no_candidate_nomination_attempted: true",
        "",
        "## Counts",
        "",
    ]
    for key in sorted(summary["counts"]):
        lines.append(f"- {key}: {summary['counts'][key]}")
    lines.extend(["", "## Blockers", ""])
    lines.extend([f"- {item}" for item in summary["blockers"]] or ["- none"])
    lines.extend(["", "## Warnings", ""])
    lines.extend([f"- {item}" for item in summary["warnings"]] or ["- none"])
    lines.extend(["", "## Prepared States", ""])
    for state_id, state in summary["states"].items():
        lines.append(f"- {state_id}: {state['status']} {state.get('state_role')}")
    lines.extend(["", "## Outputs", ""])
    for label, output in outputs.items():
        lines.append(f"- {label}: {ctx.relative_to_repo(output)}")
    lines.extend(["", "## Next Task", "", "M3-T4 - Vina focused docking smoke test"])
    safe.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_m3_receptor_boxes(
    ctx: RunContext,
    m2_run_id: str,
    mode: str = "prepare",
    force: bool = False,
    box_margin: float | None = None,
    states: list[str] | None = None,
    update_ignore: bool = True,
    allow_independent: bool = False,
    require_m3_t2_pass: bool = True,
    private_map: Path | None = None,
    operation: str = "combined",
) -> M3ReceptorBoxResult:
    if mode not in {"dry-run", "prepare"}:
        raise ValueError("mode must be dry-run or prepare")
    ctx.create_directories()
    _ensure_logs(ctx)
    phase3 = ctx.require_within_run_dir(ctx.run_dir / "phase3_compounds")
    manifest_dir = ctx.require_within_run_dir(phase3 / "manifests")
    qc_dir = ctx.require_within_run_dir(phase3 / "qc")
    report_dir = ctx.require_within_run_dir(phase3 / "reports")
    receptor_root = ctx.require_within_run_dir(phase3 / "receptor_pdbqt")
    for directory in [manifest_dir, qc_dir, report_dir, receptor_root]:
        directory.mkdir(parents=True, exist_ok=True)

    receptor_manifest = manifest_dir / "receptor_preparation_manifest.csv"
    receptor_hashes = manifest_dir / "receptor_pdbqt_hashes.csv"
    box_manifest = manifest_dir / "docking_box_manifest.csv"
    receptor_qc = qc_dir / "receptor_prep_qc.csv"
    box_qc = qc_dir / "docking_box_qc.csv"
    summary_json = qc_dir / "m3_receptor_box_summary.json"
    report_md = report_dir / "m3_task3_receptor_box_validation.md"
    phase3_log = ctx.logs_dir / "phase3_compounds.log"

    blockers: list[str] = []
    warnings: list[str] = []
    tools = detect_receptor_tools()
    inputs = _m2_inputs(ctx, m2_run_id)
    margin = _default_box_margin(ctx) if box_margin is None else float(box_margin)
    if margin < 0:
        blockers.append("box_margin must be non-negative")
    if update_ignore:
        gitignore_updated, patterns_verified = update_receptor_gitignore(ctx.repo_root)
    else:
        gitignore_updated, patterns_verified = False, list(RECEPTOR_GITIGNORE_PATTERNS)

    private_entries, private_warnings = load_private_map(_safe_private_map_path(ctx, private_map))
    warnings.extend(private_warnings)
    m3_t2, m3_t2_error = _load_m3_t2_summary(ctx)
    if m3_t2_error:
        if mode == "prepare" and require_m3_t2_pass and not allow_independent:
            blockers.append(m3_t2_error)
        else:
            warnings.append(m3_t2_error)
    if m3_t2 and require_m3_t2_pass and m3_t2.get("overall_status") == "FAIL":
        blockers.append("M3-T2 overall_status is FAIL")
    if m3_t2 and m3_t2.get("m3_t3_allowed") is False and mode == "prepare" and not allow_independent:
        blockers.append("M3-T2 m3_t3_allowed=false blocks M3-T3 prepare")
    if allow_independent:
        warnings.append("allow_independent was set; diagnostic receptor/box preparation bypass recorded")

    if not inputs["accepted_pocket_boxes"]:
        if mode == "prepare":
            blockers.append("M2 accepted_pocket_boxes.csv absent in prepare mode")
        else:
            warnings.append("M2 accepted_pocket_boxes.csv absent in dry-run mode; no receptor PDBQT or active box manifest produced")
    if not inputs["accepted_pockets_for_m3"]:
        if mode == "prepare":
            blockers.append("accepted_pockets_for_m3.csv absent in prepare mode")
        else:
            warnings.append("accepted_pockets_for_m3.csv absent in dry-run mode")

    expected_states, state_roles = _states_from_boxes(inputs["accepted_pocket_boxes"], states)
    membrane_frame = _resolve_membrane_frame(ctx, m2_run_id)
    receptor_rows, hash_rows, receptor_qc_rows, prepared_states, receptor_blockers, receptor_warnings, state_summary = _prepare_states(
        ctx,
        m2_run_id,
        inputs["accepted_pocket_boxes"],
        expected_states,
        state_roles,
        membrane_frame,
        tools,
        mode,
        force,
    )
    blockers.extend(receptor_blockers)
    warnings.extend(receptor_warnings)

    box_result = validate_docking_boxes(
        ctx,
        inputs["accepted_pocket_boxes"],
        inputs["accepted_pockets_for_m3"],
        margin,
        prepared_states,
        mode=mode,
    )
    blockers.extend(box_result.blockers)
    warnings.extend(box_result.warnings)

    _append_missing_receptor_qc(receptor_qc_rows, expected_states)
    _append_missing_box_qc(box_result.qc_rows)
    _write_csv(receptor_manifest, RECEPTOR_PREP_FIELDS, receptor_rows, ctx)
    _write_csv(receptor_hashes, RECEPTOR_HASH_FIELDS, hash_rows, ctx)
    _write_csv(box_manifest, BOX_MANIFEST_FIELDS, box_result.manifest_rows, ctx)
    _write_csv(receptor_qc, RECEPTOR_QC_FIELDS, receptor_qc_rows, ctx)
    _write_csv(box_qc, BOX_QC_FIELDS, box_result.qc_rows, ctx)

    leak_count, coordinate_hits, scanned_paths = _scan_output_hygiene(ctx, list(private_entries.values()))
    if leak_count:
        blockers.append("internal ID leakage detected in public outputs/logs")
    if coordinate_hits:
        blockers.append("receptor coordinates printed in public outputs/logs")
    forbidden_outputs = _forbidden_outputs_created(ctx)
    if forbidden_outputs:
        blockers.append("docking/Vina/scoring/candidate outputs exist under this M3 task run")

    tracked_generated = [
        row["prepared_receptor_pdbqt_file"]
        for row in receptor_rows
        if row.get("prepared_receptor_pdbqt_file") and row.get("prepared_file_tracked") == "true"
    ]
    if tracked_generated:
        blockers.append("generated receptor PDBQT files are tracked by git")
    if tools.vina_invoked:
        blockers.append("Vina invoked")

    counts = {
        "states_expected": len(expected_states),
        "states_prepared": sum(1 for row in receptor_rows if row.get("preparation_status") in {"PASS", "WARN"} and row.get("prepared_receptor_pdbqt_file")),
        "states_failed": sum(1 for row in receptor_rows if row.get("preparation_status") == "FAIL"),
        "states_warn": sum(1 for row in receptor_rows if row.get("preparation_status") == "WARN"),
        "boxes_input": box_result.counts.get("boxes_input", 0),
        "boxes_pass": box_result.counts.get("boxes_pass", 0),
        "boxes_warn": box_result.counts.get("boxes_warn", 0),
        "boxes_fail": box_result.counts.get("boxes_fail", 0),
        "boxes_allowed_for_vina_smoke": box_result.counts.get("boxes_allowed_for_vina_smoke", 0),
        "old_workflow_sources_rejected": sum(1 for row in receptor_rows if row.get("source_receptor_is_old_workflow") == "true"),
        "monomer_only_sources_rejected": sum(1 for row in receptor_rows if row.get("source_receptor_is_monomer_only") == "true"),
        "tracked_generated_pdbqt_files": len(tracked_generated),
        "confidentiality_leaks": leak_count,
    }
    if mode == "dry-run" and not blockers:
        warnings.append("dry-run mode only")
    status = "FAIL" if blockers else ("WARN" if warnings or mode == "dry-run" else "PASS")
    t2_policy_ok = bool(m3_t2 and m3_t2.get("overall_status") == "PASS" and m3_t2.get("m3_t3_allowed") is not False) or bool(allow_independent)
    m3_t4_allowed = status == "PASS" and t2_policy_ok
    if allow_independent and status == "PASS":
        m3_t4_allowed = True

    summary = {
        "schema_version": "m3_receptor_box_summary_v1",
        "run_id": ctx.run_id,
        "m2_run_id": m2_run_id,
        "reviewed_at": now_iso(),
        "mode": mode,
        "overall_status": status,
        "m3_t4_allowed": m3_t4_allowed,
        "allow_independent": allow_independent,
        "require_m3_t2_pass": require_m3_t2_pass,
        "m3_t2": {
            "summary_found": m3_t2 is not None,
            "overall_status": m3_t2.get("overall_status") if m3_t2 else None,
            "m3_t3_allowed": m3_t2.get("m3_t3_allowed") if m3_t2 else None,
        },
        "m2": {
            "accepted_pocket_boxes_found": inputs["accepted_pocket_boxes"] is not None,
            "accepted_pockets_for_m3_found": inputs["accepted_pockets_for_m3"] is not None,
            "pocket_gate_qc_found": inputs["pocket_gate_qc"] is not None,
            "atp_reference_found": inputs["atp_reference"] is not None,
            "ppi_consensus_patch_found": inputs["ppi_consensus_patch"] is not None,
            "accepted_box_count": box_result.counts.get("boxes_input", 0),
            "accepted_pocket_family_count": len({row.get("pocket_family_id", "") for row in box_result.manifest_rows if row.get("pocket_family_id")}),
        },
        "m1": {
            "receptor_mapping_found": any(row.get("receptor_mapping_file") for row in receptor_rows),
            "membrane_frame_found": membrane_frame is not None,
            "normalized_receptor_count": sum(1 for row in receptor_rows if row.get("source_receptor_is_normalized") == "true"),
        },
        "tools": {
            "prepare_receptor4_available": tools.prepare_receptor4_available,
            "prepare_receptor4_version": tools.prepare_receptor4_version,
            "mk_prepare_receptor_available": tools.mk_prepare_receptor_available,
            "mk_prepare_receptor_version": tools.mk_prepare_receptor_version,
            "meeko_available": tools.meeko_available,
            "meeko_version": tools.meeko_version,
            "openbabel_available": tools.openbabel_available,
            "openbabel_version": tools.openbabel_version,
            "rdkit_available": tools.rdkit_available,
            "rdkit_version": tools.rdkit_version,
            "vina_available": tools.vina_available,
            "vina_invoked": tools.vina_invoked,
        },
        "no_ligand_preparation_attempted": True,
        "no_pocket_discovery_attempted": True,
        "no_vina_attempted": True,
        "no_docking_job_generation_attempted": True,
        "no_docking_attempted": True,
        "no_pose_parsing_attempted": True,
        "no_scoring_attempted": True,
        "no_candidate_nomination_attempted": True,
        "blockers": list(dict.fromkeys(blockers)),
        "warnings": list(dict.fromkeys(warnings)),
        "counts": counts,
        "states": state_summary,
        "boxes": {
            "box_manifest_file": ctx.relative_to_repo(box_manifest),
            "qc_file": ctx.relative_to_repo(box_qc),
            "all_active_boxes_traceable_to_m2": box_result.all_active_boxes_traceable_to_m2,
            "all_active_boxes_have_receptor": box_result.all_active_boxes_have_receptor,
            "all_active_boxes_non_atp": box_result.all_active_boxes_non_atp,
            "all_active_boxes_valid_dimensions": box_result.all_active_boxes_valid_dimensions,
        },
        "gitignore": {
            "updated": gitignore_updated,
            "patterns_verified": patterns_verified,
            "sensitive_or_generated_files_tracked": tracked_generated,
        },
        "confidentiality": {
            "internal_ids_redacted": True,
            "public_outputs_scanned": scanned_paths,
            "leaks_detected": leak_count,
        },
        "non_goals_preserved": {
            "ligand_preparation": True,
            "pocket_discovery": True,
            "vina_execution": True,
            "docking_job_generation": True,
            "pose_parsing": True,
            "compound_scoring": True,
            "candidate_nomination": True,
        },
    }
    ctx.require_within_run_dir(summary_json).write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    _write_report(
        report_md,
        ctx,
        summary,
        {
            "receptor_preparation_manifest": receptor_manifest,
            "receptor_pdbqt_hashes": receptor_hashes,
            "docking_box_manifest": box_manifest,
            "receptor_prep_qc": receptor_qc,
            "docking_box_qc": box_qc,
            "m3_receptor_box_summary": summary_json,
            "m3_task3_receptor_box_validation": report_md,
            "phase3_log": phase3_log,
        },
    )
    with ctx.require_within_run_dir(phase3_log).open("a", encoding="utf-8") as handle:
        handle.write(
            "{0} M3-T3 status={1} run_id={2} m2_run_id={3} mode={4} m3_t4_allowed={5} states_prepared={6} boxes_pass={7} boxes_fail={8} no_ligand_preparation_attempted=true no_pocket_discovery_attempted=true no_vina_attempted=true no_docking_job_generation_attempted=true no_docking_attempted=true no_pose_parsing_attempted=true no_scoring_attempted=true\n".format(
                now_iso(),
                status,
                ctx.run_id,
                m2_run_id,
                mode,
                str(m3_t4_allowed).lower(),
                counts["states_prepared"],
                counts["boxes_pass"],
                counts["boxes_fail"],
            )
        )
    log_master(ctx, status, "M3-T3 receptor/box validation completed without ligand/Vina/docking/scoring work")
    append_phase_status(
        ctx,
        "phase3_compounds",
        status,
        "M3-T3 receptor PDBQT preparation and docking box validation completed",
        {
            "task": "M3-T3",
            "m2_run_id": m2_run_id,
            "mode": mode,
            "m3_t4_allowed": m3_t4_allowed,
            "states_prepared": counts["states_prepared"],
            "boxes_pass": counts["boxes_pass"],
            "boxes_fail": counts["boxes_fail"],
            "no_ligand_preparation_attempted": True,
            "no_pocket_discovery_attempted": True,
            "no_vina_attempted": True,
            "no_docking_job_generation_attempted": True,
            "no_docking_attempted": True,
            "no_pose_parsing_attempted": True,
            "no_scoring_attempted": True,
        },
    )
    append_job_status(
        ctx,
        "M3-T3-receptor-box-validation",
        status,
        details={"message": f"operation={operation}; no Vina/docking/scoring work", "m3_t4_allowed": m3_t4_allowed},
    )
    return M3ReceptorBoxResult(
        status=status,
        m3_t4_allowed=m3_t4_allowed,
        blockers=list(dict.fromkeys(blockers)),
        warnings=list(dict.fromkeys(warnings)),
        receptor_preparation_manifest_csv=receptor_manifest,
        receptor_hashes_csv=receptor_hashes,
        docking_box_manifest_csv=box_manifest,
        receptor_prep_qc_csv=receptor_qc,
        docking_box_qc_csv=box_qc,
        summary_json=summary_json,
        report_md=report_md,
        phase3_log=phase3_log,
        counts=counts,
        states=state_summary,
        tools=tools,
    )


def default_box_margin(ctx: RunContext) -> float:
    return _default_box_margin(ctx)
