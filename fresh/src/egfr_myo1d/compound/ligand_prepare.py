"""M3-T2 ligand preparation and PDBQT generation.

This module prepares Cpd-A/B/C ligand PDBQT files only. It does not prepare
receptors, generate docking boxes, run Vina, create docking jobs, parse poses,
score compounds, rank candidates, or create candidate claims.
"""

from __future__ import annotations

import csv
import contextlib
import io
import importlib
import json
import shutil
import subprocess
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from egfr_myo1d.compound.confidentiality import (
    PRIVATE_PATH_REDACTED,
    PUBLIC_COMPOUND_IDS,
    RECOMMENDED_GITIGNORE_PATTERNS,
    PrivateMapEntry,
    assert_fresh_input_path,
    git_check_ignored,
    git_ls_files,
    load_private_map,
    public_output_scan_paths,
    scan_internal_id_leaks,
    update_gitignore,
)
from egfr_myo1d.compound.ligand_manifest import (
    ALLOWED_EXTENSIONS,
    basic_parse_check,
    resolve_ligand,
)
from egfr_myo1d.core.logging_utils import append_job_status, append_phase_status, log_master
from egfr_myo1d.core.run_context import RunContext


PREP_MANIFEST_FIELDS = [
    "compound_public_id",
    "source_input_file",
    "source_input_format",
    "source_input_sha256",
    "source_input_size_bytes",
    "has_private_mapping",
    "resolved_by_private_mapping",
    "prepared_pdbqt_file",
    "prepared_pdbqt_sha256",
    "prepared_pdbqt_size_bytes",
    "preparation_method",
    "conversion_tool",
    "conversion_tool_version",
    "rdkit_available",
    "openbabel_available",
    "rdkit_parse_success",
    "obabel_conversion_success",
    "pdbqt_validation_success",
    "multi_molecule_input",
    "stereochemistry_defined",
    "formal_charge",
    "mol_weight",
    "num_heavy_atoms",
    "num_rotatable_bonds",
    "num_hbd",
    "num_hba",
    "protonation_status_known",
    "tautomer_enumeration_performed",
    "protonation_enumeration_performed",
    "raw_input_modified",
    "allowed_for_public_output",
    "preparation_status",
    "preparation_notes",
]
HASH_FIELDS = [
    "compound_public_id",
    "prepared_pdbqt_file",
    "prepared_pdbqt_sha256",
    "prepared_pdbqt_size_bytes",
    "source_input_file",
    "source_input_sha256",
    "hash_status",
    "notes",
]
QC_FIELDS = [
    "compound_public_id",
    "check_id",
    "category",
    "status",
    "severity",
    "details",
    "recommended_fix",
]


@dataclass(frozen=True)
class ToolStatus:
    rdkit_available: bool
    rdkit_version: str | None
    openbabel_available: bool
    openbabel_version: str | None
    meeko_available: bool
    meeko_version: str | None
    obabel_binary: str = "obabel"


@dataclass
class M3LigandPrepResult:
    status: str
    m3_t3_allowed: bool
    blockers: list[str]
    warnings: list[str]
    preparation_manifest_csv: Path
    prepared_hashes_csv: Path
    prep_qc_csv: Path
    summary_json: Path
    report_md: Path
    phase3_log: Path
    counts: dict[str, int]


def now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def _bool_text(value: bool) -> str:
    return "true" if value else "false"


def _boolish(value: Any) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "y", "pass"}


def _sha256(path: Path) -> str:
    import hashlib

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _safe_write_csv(path: Path, fields: list[str], rows: list[dict[str, Any]], ctx: RunContext) -> None:
    safe = ctx.require_within_run_dir(path)
    safe.parent.mkdir(parents=True, exist_ok=True)
    with safe.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _read_csv_rows(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        return list(reader), list(reader.fieldnames or [])


def _severity(status: str, blocker: bool = False) -> str:
    if status == "FAIL" and blocker:
        return "BLOCKER"
    if status == "FAIL":
        return "MAJOR"
    if status == "WARN":
        return "MAJOR"
    if status == "NOT_APPLICABLE":
        return "MINOR"
    return "INFO"


def _append_qc(
    rows: list[dict[str, str]],
    public_id: str,
    check_id: str,
    category: str,
    status: str,
    details: str,
    recommended_fix: str = "",
    blocker: bool = False,
) -> None:
    rows.append(
        {
            "compound_public_id": public_id,
            "check_id": check_id,
            "category": category,
            "status": status,
            "severity": _severity(status, blocker),
            "details": details,
            "recommended_fix": recommended_fix,
        }
    )


def detect_tools(obabel_binary: str = "obabel") -> ToolStatus:
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

    meeko_available = False
    meeko_version = None
    try:
        meeko = importlib.import_module("meeko")
        meeko_available = True
        meeko_version = getattr(meeko, "__version__", None)
    except Exception:
        pass

    openbabel_available = False
    openbabel_version = None
    try:
        proc = subprocess.run(
            [obabel_binary, "-V"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )
        openbabel_available = proc.returncode == 0
        text = (proc.stdout or proc.stderr or "").strip()
        openbabel_version = text.splitlines()[0] if text else None
    except OSError:
        pass

    return ToolStatus(
        rdkit_available=rdkit_available,
        rdkit_version=rdkit_version,
        openbabel_available=openbabel_available,
        openbabel_version=openbabel_version,
        meeko_available=meeko_available,
        meeko_version=meeko_version,
        obabel_binary=obabel_binary,
    )


def _safe_private_map_path(ctx: RunContext, private_map: Path | None) -> Path:
    path = private_map or (ctx.fresh_root / "data" / "private" / "compound_id_map.csv")
    resolved = Path(path)
    if not resolved.is_absolute():
        resolved = ctx.repo_root / resolved
    return assert_fresh_input_path(ctx, resolved.resolve())


def _load_m3_t1_summary(ctx: RunContext) -> tuple[dict[str, Any] | None, str | None]:
    path = ctx.run_dir / "phase3_compounds" / "qc" / "m3_ligand_qc_summary.json"
    if not path.is_file():
        return None, "M3-T1 ligand QC summary not found"
    try:
        return json.loads(path.read_text(encoding="utf-8")), None
    except (OSError, json.JSONDecodeError):
        return None, "M3-T1 ligand QC summary is unreadable"


def _load_manifest(ctx: RunContext) -> tuple[list[dict[str, str]], str | None]:
    path = ctx.run_dir / "phase3_compounds" / "manifests" / "ligand_manifest_public.csv"
    if not path.is_file():
        return [], "M3-T1 ligand manifest not found"
    try:
        rows, _ = _read_csv_rows(path)
        return rows, None
    except (OSError, csv.Error, UnicodeDecodeError):
        return [], "M3-T1 ligand manifest is unreadable"


def _validate_pdbqt(path: Path) -> tuple[bool, str]:
    try:
        text = path.read_text(encoding="utf-8", errors="ignore")
    except OSError:
        return False, "prepared PDBQT could not be read"
    if not text.strip():
        return False, "prepared PDBQT is empty"
    atom_lines = [line for line in text.splitlines() if line.startswith("ATOM") or line.startswith("HETATM")]
    if not atom_lines:
        return False, "prepared PDBQT has no ATOM/HETATM records"
    return True, "prepared PDBQT has atom records"


def _obabel_input_flag(suffix: str) -> str:
    return {
        ".sdf": "-isdf",
        ".mol2": "-imol2",
        ".pdb": "-ipdb",
        ".smi": "-ismi",
        ".smiles": "-ismi",
    }.get(suffix, "-i" + suffix.lstrip("."))


def _convert_with_obabel(
    source: Path,
    output: Path,
    public_id: str,
    ctx: RunContext,
    tools: ToolStatus,
    source_display: str,
) -> tuple[bool, str]:
    ctx.require_within_run_dir(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    stdout_path = ctx.require_within_run_dir(ctx.jobs_log_dir / f"m3_t2_{public_id}_obabel.stdout.txt")
    stderr_path = ctx.require_within_run_dir(ctx.jobs_log_dir / f"m3_t2_{public_id}_obabel.stderr.txt")
    args = [
        tools.obabel_binary,
        _obabel_input_flag(source.suffix.lower()),
        str(source),
        "-opdbqt",
        "-O",
        str(output),
    ]
    if source.suffix.lower() in {".smi", ".smiles"}:
        args.insert(-2, "--gen3d")
    try:
        proc = subprocess.run(
            args,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )
    except OSError as exc:
        stdout_path.write_text("", encoding="utf-8")
        stderr_path.write_text("Open Babel invocation failed: {0}\n".format(exc.__class__.__name__), encoding="utf-8")
        return False, "Open Babel invocation failed"
    stdout_path.write_text(proc.stdout or "", encoding="utf-8")
    stderr_path.write_text(proc.stderr or "", encoding="utf-8")
    if proc.returncode != 0:
        return False, "Open Babel conversion failed; stdout/stderr captured in logs/jobs"
    if not output.is_file():
        return False, "Open Babel did not create a PDBQT output"
    ok, message = _validate_pdbqt(output)
    return ok, message


def _rdkit_metrics(path: Path, tools: ToolStatus) -> dict[str, Any]:
    metrics: dict[str, Any] = {
        "rdkit_parse_success": "NOT_APPLICABLE",
        "stereochemistry_defined": "",
        "formal_charge": "",
        "mol_weight": "",
        "num_heavy_atoms": "",
        "num_rotatable_bonds": "",
        "num_hbd": "",
        "num_hba": "",
        "protonation_status_known": "false",
    }
    if not tools.rdkit_available:
        return metrics
    try:
        with contextlib.redirect_stderr(io.StringIO()):
            from rdkit import Chem
            from rdkit.Chem import Descriptors, Lipinski
    except BaseException:
        return metrics
    suffix = path.suffix.lower()
    mol = None
    try:
        if suffix == ".sdf":
            supplier = Chem.SDMolSupplier(str(path), removeHs=False)
            mol = next((m for m in supplier if m is not None), None)
        elif suffix in {".smi", ".smiles"}:
            first = None
            for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
                if line.strip() and not line.lstrip().startswith("#"):
                    first = line.split()[0]
                    break
            mol = Chem.MolFromSmiles(first) if first else None
        elif suffix == ".mol2":
            mol = Chem.MolFromMol2File(str(path), sanitize=True, removeHs=False)
        elif suffix == ".pdb":
            mol = Chem.MolFromPDBFile(str(path), sanitize=True, removeHs=False)
    except BaseException:
        mol = None
    if mol is None:
        metrics["rdkit_parse_success"] = "false"
        return metrics
    metrics["rdkit_parse_success"] = "true"
    metrics["formal_charge"] = str(sum(atom.GetFormalCharge() for atom in mol.GetAtoms()))
    metrics["mol_weight"] = "{0:.3f}".format(Descriptors.MolWt(mol))
    metrics["num_heavy_atoms"] = str(mol.GetNumHeavyAtoms())
    metrics["num_rotatable_bonds"] = str(Lipinski.NumRotatableBonds(mol))
    metrics["num_hbd"] = str(Lipinski.NumHDonors(mol))
    metrics["num_hba"] = str(Lipinski.NumHAcceptors(mol))
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    metrics["stereochemistry_defined"] = "false" if any(center[1] == "?" for center in chiral_centers) else "true"
    return metrics


def _multi_molecule_flag(path: Path) -> tuple[bool, str | None]:
    suffix = path.suffix.lower()
    try:
        text = path.read_text(encoding="utf-8", errors="ignore")
    except OSError:
        return False, None
    if suffix == ".sdf" and text.count("$$$$") > 1:
        return True, "multiple SDF records present; first molecule selected by converter/tool default"
    if suffix in {".smi", ".smiles"}:
        lines = [line for line in text.splitlines() if line.strip() and not line.lstrip().startswith("#")]
        if len(lines) > 1:
            return True, "multiple SMILES records present; first record policy applies"
    return False, None


def _ensure_logs(ctx: RunContext) -> None:
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


def _write_report(
    path: Path,
    ctx: RunContext,
    status: str,
    m3_t3_allowed: bool,
    blockers: list[str],
    warnings: list[str],
    counts: dict[str, int],
    outputs: dict[str, Path],
) -> None:
    safe = ctx.require_within_run_dir(path)
    safe.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# M3-T2 Ligand Preparation and PDBQT Generation",
        "",
        f"- run_id: {ctx.run_id}",
        f"- overall_status: {status}",
        f"- m3_t3_allowed: {str(m3_t3_allowed).lower()}",
        "- compounds: confidential MYO1D activity-associated compounds",
        "- no_receptor_pdbqt_generation_attempted: true",
        "- no_docking_box_generation_attempted: true",
        "- no_vina_attempted: true",
        "- no_docking_attempted: true",
        "",
        "## Counts",
        "",
    ]
    for key in sorted(counts):
        lines.append(f"- {key}: {counts[key]}")
    lines.extend(["", "## Blockers", ""])
    lines.extend([f"- {item}" for item in blockers] or ["- none"])
    lines.extend(["", "## Warnings", ""])
    lines.extend([f"- {item}" for item in warnings] or ["- none"])
    lines.extend(["", "## Outputs", ""])
    for label, output in outputs.items():
        lines.append(f"- {label}: {ctx.relative_to_repo(output)}")
    lines.extend(["", "## Next Task", "", "M3-T3 - receptor PDBQT preparation and docking box validation"])
    safe.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _manifest_public_ids(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    return {row.get("compound_public_id", ""): row for row in rows}


def _source_from_manifest(
    ctx: RunContext,
    public_id: str,
    manifest_row: dict[str, str] | None,
    private_entries: dict[str, PrivateMapEntry],
) -> tuple[Path | None, str, str, bool, bool, str | None]:
    if manifest_row is None:
        return None, "", "", False, False, "ligand manifest row missing"
    input_file = (manifest_row.get("input_file") or "").strip()
    resolved_by_private = input_file == PRIVATE_PATH_REDACTED or manifest_row.get("resolved_by_private_mapping") == "true"
    has_private = manifest_row.get("has_private_mapping") == "true" or public_id in private_entries
    if resolved_by_private:
        resolution = resolve_ligand(ctx, public_id, private_entries.get(public_id))
        if resolution.selected_path is None:
            return None, PRIVATE_PATH_REDACTED, manifest_row.get("input_format", ""), has_private, True, "private-mapped ligand input could not be resolved without exposing private ID"
        return resolution.selected_path, PRIVATE_PATH_REDACTED, resolution.input_format, has_private, True, None
    if not input_file:
        return None, "", manifest_row.get("input_format", ""), has_private, False, "ligand input file is missing"
    path = Path(input_file)
    if not path.is_absolute():
        path = ctx.repo_root / path
    path = assert_fresh_input_path(ctx, path.resolve())
    return path, ctx.relative_to_repo(path), path.suffix.lower().lstrip("."), has_private, False, None


def run_m3_ligand_prepare(
    ctx: RunContext,
    mode: str = "prepare",
    private_map: Path | None = None,
    force: bool = False,
    keep_intermediates: bool = False,
    update_ignore: bool = True,
    allow_independent: bool = False,
    require_m3_t1_pass: bool = True,
    obabel_binary: str = "obabel",
) -> M3LigandPrepResult:
    if mode not in {"dry-run", "prepare"}:
        raise ValueError("mode must be dry-run or prepare")

    ctx.create_directories()
    _ensure_logs(ctx)
    phase3 = ctx.run_dir / "phase3_compounds"
    manifest_dir = ctx.require_within_run_dir(phase3 / "manifests")
    qc_dir = ctx.require_within_run_dir(phase3 / "qc")
    reports_dir = ctx.require_within_run_dir(phase3 / "reports")
    prep_root = ctx.require_within_run_dir(phase3 / "prepared_ligands")
    for directory in [manifest_dir, qc_dir, reports_dir, prep_root]:
        directory.mkdir(parents=True, exist_ok=True)

    prep_manifest = manifest_dir / "ligand_preparation_manifest.csv"
    prepared_hashes = manifest_dir / "prepared_ligand_hashes.csv"
    prep_qc = qc_dir / "ligand_prep_qc.csv"
    summary_json = qc_dir / "m3_ligand_prep_summary.json"
    report_md = reports_dir / "m3_task2_ligand_preparation.md"
    phase3_log = ctx.logs_dir / "phase3_compounds.log"

    blockers: list[str] = []
    warnings: list[str] = []
    manifest_rows: list[dict[str, Any]] = []
    hash_rows: list[dict[str, Any]] = []
    qc_rows: list[dict[str, str]] = []
    ligand_summary: dict[str, dict[str, Any]] = {}

    tools = detect_tools(obabel_binary=obabel_binary)
    if update_ignore:
        gitignore_updated, patterns_verified = update_gitignore(ctx.repo_root)
    else:
        gitignore_updated, patterns_verified = False, list(RECOMMENDED_GITIGNORE_PATTERNS)

    private_map_path = _safe_private_map_path(ctx, private_map)
    private_entries, private_warnings = load_private_map(private_map_path)
    if private_warnings:
        warnings.extend(private_warnings)
    private_map_tracked = bool(git_ls_files(ctx.repo_root, private_map_path)) if private_map_path.is_file() else False
    tracked_sensitive: list[str] = []
    if private_map_tracked:
        blockers.append("private compound_id_map.csv is tracked by git")
        tracked_sensitive.append("PRIVATE_MAPPING")

    m3_t1_summary, m3_t1_error = _load_m3_t1_summary(ctx)
    if m3_t1_error:
        if mode == "prepare":
            blockers.append(m3_t1_error)
        else:
            warnings.append(m3_t1_error)
    if m3_t1_summary:
        if require_m3_t1_pass and m3_t1_summary.get("overall_status") == "FAIL":
            blockers.append("M3-T1 overall_status is FAIL")
        if m3_t1_summary.get("m3_t2_allowed") is False and mode == "prepare" and not allow_independent:
            blockers.append("M3-T1 m3_t2_allowed=false blocks ligand preparation")
        if (m3_t1_summary.get("counts") or {}).get("confidentiality_leaks", 0):
            blockers.append("M3-T1 confidentiality leaks block ligand preparation")
    if allow_independent:
        warnings.append("allow_independent was set; M3-T2 is running outside normal M3-T1 pass gate")

    t1_manifest_rows, manifest_error = _load_manifest(ctx)
    if manifest_error:
        if mode == "prepare":
            blockers.append(manifest_error)
        else:
            warnings.append(manifest_error)
    manifest_by_id = _manifest_public_ids(t1_manifest_rows)
    _append_qc(qc_rows, "", "m3_t1_manifest_present", "precondition", "PASS" if not manifest_error else ("FAIL" if mode == "prepare" else "WARN"), manifest_error or "M3-T1 ligand manifest present", "Run M3-T1 ligand QC first.", mode == "prepare" and bool(manifest_error))
    _append_qc(qc_rows, "", "m3_t1_confidentiality_pass", "precondition", "PASS" if not (m3_t1_summary and (m3_t1_summary.get("counts") or {}).get("confidentiality_leaks", 0)) else "FAIL", "M3-T1 confidentiality gate checked", "Fix M3-T1 confidentiality audit.", bool(m3_t1_summary and (m3_t1_summary.get("counts") or {}).get("confidentiality_leaks", 0)))

    prepared_count = 0
    reused_count = 0
    failed_count = 0
    warn_count = 0
    private_used_count = 0

    for public_id in PUBLIC_COMPOUND_IDS:
        notes: list[str] = []
        status = "NOT_APPLICABLE" if mode == "dry-run" else "PASS"
        method = ""
        conversion_tool = ""
        conversion_version = ""
        manifest_row = manifest_by_id.get(public_id, {})
        source_path, source_display, source_format, has_private, resolved_private, source_error = _source_from_manifest(
            ctx,
            public_id,
            manifest_row,
            private_entries,
        )
        if resolved_private:
            private_used_count += 1
        prepared_path = ctx.require_within_run_dir(prep_root / public_id / f"{public_id}.pdbqt")
        prepared_display = ctx.relative_to_repo(prepared_path)
        source_sha = ""
        source_size = ""
        prepared_sha = ""
        prepared_size = ""
        raw_before = raw_after = None
        metrics = {
            "rdkit_parse_success": "NOT_APPLICABLE",
            "stereochemistry_defined": "",
            "formal_charge": "",
            "mol_weight": "",
            "num_heavy_atoms": "",
            "num_rotatable_bonds": "",
            "num_hbd": "",
            "num_hba": "",
            "protonation_status_known": "false",
        }
        multi_molecule = False
        pdbqt_valid = False
        obabel_success = "NOT_APPLICABLE"
        rdkit_parse_success = "NOT_APPLICABLE"
        raw_modified = False
        source_tracked = False
        prepared_tracked = False
        prepared_ignored = git_check_ignored(ctx.repo_root, prepared_path)

        if source_error:
            status = "WARN" if mode == "dry-run" else "FAIL"
            notes.append(source_error)
            if mode == "prepare":
                blockers.append(f"{public_id}: {source_error}")
            else:
                warnings.append(f"{public_id}: {source_error}")
        elif source_path is not None:
            source_format = source_path.suffix.lower().lstrip(".")
            source_sha = _sha256(source_path)
            source_size = str(source_path.stat().st_size)
            manifest_sha = str(manifest_row.get("file_sha256", "")).strip()
            manifest_size = str(manifest_row.get("size_bytes", "")).strip()
            if manifest_sha and manifest_sha != source_sha:
                status = "FAIL"
                blockers.append(f"{public_id}: source ligand SHA256 differs from M3-T1 manifest")
                notes.append("source hash mismatch with M3-T1 manifest")
            if manifest_size and manifest_size != source_size:
                status = "FAIL"
                blockers.append(f"{public_id}: source ligand size differs from M3-T1 manifest")
                notes.append("source size mismatch with M3-T1 manifest")
            raw_before = source_sha
            source_tracked = bool(git_ls_files(ctx.repo_root, source_path))
            if source_tracked:
                status = "FAIL"
                blockers.append(f"{public_id}: raw ligand file is tracked by git")
                tracked_sensitive.append(f"RAW_LIGAND_{public_id}")
            if source_path.suffix.lower() not in ALLOWED_EXTENSIONS:
                status = "FAIL"
                blockers.append(f"{public_id}: unsupported input format")
                notes.append("unsupported input format")
            if source_path.stat().st_size == 0:
                status = "FAIL"
                blockers.append(f"{public_id}: zero-byte ligand input")
                notes.append("zero-byte ligand input")
            parse_status, parse_notes = basic_parse_check(source_path)
            if parse_status == "FAIL":
                status = "FAIL"
                blockers.append(f"{public_id}: input parse failed")
            elif parse_status == "WARN" and status != "FAIL":
                status = "WARN" if mode == "prepare" else "NOT_APPLICABLE"
                warnings.append(f"{public_id}: input parse warning")
            notes.append(parse_notes)
            multi_molecule, multi_note = _multi_molecule_flag(source_path)
            if multi_note:
                if status == "PASS":
                    status = "WARN"
                warnings.append(f"{public_id}: multiple molecule records detected")
                notes.append(multi_note)
            if source_path.suffix.lower() in {".smi", ".smiles"} and not tools.rdkit_available and tools.openbabel_available:
                if status == "PASS":
                    status = "WARN"
                warnings.append(f"{public_id}: Open Babel-only SMILES 3D generation used")
                notes.append("Open Babel-only SMILES 3D generation used")
            metrics = _rdkit_metrics(source_path, tools)
            rdkit_parse_success = metrics["rdkit_parse_success"]

            if mode == "dry-run":
                method = "dry_run_no_pdbqt_written"
            elif status != "FAIL":
                prepared_path.parent.mkdir(parents=True, exist_ok=True)
                if prepared_path.exists() and not force:
                    status = "FAIL"
                    blockers.append(f"{public_id}: prepared PDBQT already exists and force=false")
                    notes.append("prepared output already exists")
                elif source_path.suffix.lower() == ".pdbqt":
                    if not _boolish(manifest_row.get("supplied_pdbqt_manual_review")):
                        status = "FAIL"
                        blockers.append(f"{public_id}: supplied PDBQT requires explicit manual-review approval")
                        notes.append("supplied PDBQT reuse not approved")
                    else:
                        shutil.copyfile(source_path, prepared_path)
                        method = "supplied_pdbqt_validated"
                        conversion_tool = "none"
                        conversion_version = ""
                        warnings.append(f"{public_id}: supplied PDBQT reused rather than generated")
                        if status == "PASS":
                            status = "WARN"
                        reused_count += 1
                else:
                    if not tools.openbabel_available:
                        status = "FAIL"
                        blockers.append(f"{public_id}: Open Babel unavailable for required conversion route")
                        notes.append("Open Babel unavailable for conversion")
                    else:
                        ok, message = _convert_with_obabel(
                            source_path,
                            prepared_path,
                            public_id,
                            ctx,
                            tools,
                            source_display,
                        )
                        obabel_success = "true" if ok else "false"
                        method = "openbabel_conversion"
                        conversion_tool = "obabel"
                        conversion_version = tools.openbabel_version or ""
                        notes.append(message)
                        if not ok:
                            status = "FAIL"
                            blockers.append(f"{public_id}: conversion failed")
                if prepared_path.is_file():
                    pdbqt_valid, pdbqt_message = _validate_pdbqt(prepared_path)
                    notes.append(pdbqt_message)
                    if not pdbqt_valid:
                        status = "FAIL"
                        blockers.append(f"{public_id}: invalid prepared PDBQT")
                    prepared_sha = _sha256(prepared_path)
                    prepared_size = str(prepared_path.stat().st_size)
                    prepared_tracked = bool(git_ls_files(ctx.repo_root, prepared_path))
                    if prepared_tracked:
                        status = "FAIL"
                        blockers.append(f"{public_id}: prepared PDBQT is tracked by git")
                        tracked_sensitive.append(f"PREPARED_PDBQT_{public_id}")
                    prepared_ignored = git_check_ignored(ctx.repo_root, prepared_path)
                    if not prepared_ignored and status != "FAIL":
                        status = "WARN"
                        warnings.append(f"{public_id}: prepared PDBQT is not gitignore-protected")
                if source_path.is_file():
                    raw_after = _sha256(source_path)
                    raw_modified = raw_before != raw_after
                    if raw_modified:
                        status = "FAIL"
                        blockers.append(f"{public_id}: raw ligand input was modified")

        if status == "FAIL":
            failed_count += 1
        elif status == "WARN":
            warn_count += 1
            if prepared_path.is_file():
                prepared_count += 1
        elif status == "PASS":
            prepared_count += 1
        elif status == "NOT_APPLICABLE" and mode == "dry-run":
            warnings.append(f"{public_id}: dry-run mode did not prepare PDBQT")

        _append_qc(qc_rows, public_id, "input_file_resolved", "input", "PASS" if source_path else ("FAIL" if mode == "prepare" else "WARN"), "input path resolved" if source_path else "input path missing", "Run M3-T1 with valid ligand inputs.", mode == "prepare" and not source_path)
        _append_qc(qc_rows, public_id, "input_format_supported", "input", "PASS" if source_path and source_path.suffix.lower() in ALLOWED_EXTENSIONS else ("NOT_APPLICABLE" if not source_path else "FAIL"), "input format checked", "Use .sdf/.mol2/.pdb/.pdbqt/.smi/.smiles.", bool(source_path and source_path.suffix.lower() not in ALLOWED_EXTENSIONS))
        _append_qc(qc_rows, public_id, "input_file_nonempty", "input", "PASS" if source_path and source_path.is_file() and source_path.stat().st_size > 0 else ("NOT_APPLICABLE" if not source_path else "FAIL"), "input file non-empty" if source_path else "no input file", "Provide non-empty ligand input.", bool(source_path and source_path.stat().st_size <= 0))
        _append_qc(qc_rows, public_id, "source_hash_recorded", "provenance", "PASS" if source_sha else ("NOT_APPLICABLE" if mode == "dry-run" and not source_path else "FAIL"), "source hash recorded" if source_sha else "source hash unavailable", "Provide a readable ligand input.", mode == "prepare" and not source_sha)
        _append_qc(qc_rows, public_id, "rdkit_parse", "chemistry_qc", "PASS" if rdkit_parse_success == "true" else ("FAIL" if rdkit_parse_success == "false" else "NOT_APPLICABLE"), "RDKit parse status recorded", "Install RDKit for chemistry descriptors if required.", False)
        _append_qc(qc_rows, public_id, "openbabel_conversion", "conversion", "PASS" if obabel_success == "true" else ("FAIL" if obabel_success == "false" else "NOT_APPLICABLE"), "Open Babel conversion status recorded", "Install Open Babel or provide valid PDBQT input.", obabel_success == "false")
        _append_qc(qc_rows, public_id, "pdbqt_generated_or_supplied", "output", "PASS" if prepared_path.is_file() else ("NOT_APPLICABLE" if mode == "dry-run" else "FAIL"), "prepared PDBQT output checked", "Prepare or supply PDBQT.", mode == "prepare" and not prepared_path.is_file())
        _append_qc(qc_rows, public_id, "pdbqt_nonempty", "output", "PASS" if prepared_path.is_file() and prepared_path.stat().st_size > 0 else ("NOT_APPLICABLE" if not prepared_path.is_file() else "FAIL"), "prepared PDBQT non-empty check", "Regenerate PDBQT.", prepared_path.is_file() and prepared_path.stat().st_size <= 0)
        _append_qc(qc_rows, public_id, "pdbqt_atom_records", "output", "PASS" if pdbqt_valid else ("NOT_APPLICABLE" if not prepared_path.is_file() else "FAIL"), "prepared PDBQT atom-record validation", "Regenerate valid PDBQT.", prepared_path.is_file() and not pdbqt_valid)
        _append_qc(qc_rows, public_id, "prepared_hash_recorded", "provenance", "PASS" if prepared_sha else ("NOT_APPLICABLE" if not prepared_path.is_file() else "FAIL"), "prepared hash recorded" if prepared_sha else "prepared hash unavailable", "Regenerate PDBQT.", prepared_path.is_file() and not prepared_sha)
        _append_qc(qc_rows, public_id, "raw_input_unchanged", "provenance", "PASS" if not raw_modified else "FAIL", "raw ligand input unchanged", "Restore raw input and rerun.", raw_modified)
        _append_qc(qc_rows, public_id, "prepared_file_gitignored", "confidentiality", "PASS" if prepared_ignored or not prepared_path.exists() else "WARN", "prepared file gitignore protection checked", "Add prepared ligand ignore pattern.", False)
        _append_qc(qc_rows, public_id, "prepared_file_not_tracked", "confidentiality", "FAIL" if prepared_tracked else "PASS", "prepared PDBQT git tracking checked", "Remove from git index without deleting local file.", prepared_tracked)

        allowed_public = not resolved_private
        manifest_rows.append(
            {
                "compound_public_id": public_id,
                "source_input_file": source_display,
                "source_input_format": source_format,
                "source_input_sha256": source_sha,
                "source_input_size_bytes": source_size,
                "has_private_mapping": _bool_text(has_private),
                "resolved_by_private_mapping": _bool_text(resolved_private),
                "prepared_pdbqt_file": prepared_display if prepared_path.is_file() else "",
                "prepared_pdbqt_sha256": prepared_sha,
                "prepared_pdbqt_size_bytes": prepared_size,
                "preparation_method": method,
                "conversion_tool": conversion_tool,
                "conversion_tool_version": conversion_version,
                "rdkit_available": _bool_text(tools.rdkit_available),
                "openbabel_available": _bool_text(tools.openbabel_available),
                "rdkit_parse_success": rdkit_parse_success,
                "obabel_conversion_success": obabel_success,
                "pdbqt_validation_success": _bool_text(pdbqt_valid),
                "multi_molecule_input": _bool_text(multi_molecule),
                "stereochemistry_defined": metrics["stereochemistry_defined"],
                "formal_charge": metrics["formal_charge"],
                "mol_weight": metrics["mol_weight"],
                "num_heavy_atoms": metrics["num_heavy_atoms"],
                "num_rotatable_bonds": metrics["num_rotatable_bonds"],
                "num_hbd": metrics["num_hbd"],
                "num_hba": metrics["num_hba"],
                "protonation_status_known": metrics["protonation_status_known"],
                "tautomer_enumeration_performed": "false",
                "protonation_enumeration_performed": "false",
                "raw_input_modified": _bool_text(raw_modified),
                "allowed_for_public_output": _bool_text(allowed_public),
                "preparation_status": status,
                "preparation_notes": "; ".join(notes),
            }
        )
        hash_rows.append(
            {
                "compound_public_id": public_id,
                "prepared_pdbqt_file": prepared_display if prepared_path.is_file() else "",
                "prepared_pdbqt_sha256": prepared_sha,
                "prepared_pdbqt_size_bytes": prepared_size,
                "source_input_file": source_display,
                "source_input_sha256": source_sha,
                "hash_status": "PASS" if prepared_sha and source_sha else ("NOT_APPLICABLE" if mode == "dry-run" else "FAIL"),
                "notes": "hashes recorded" if prepared_sha and source_sha else "hash incomplete",
            }
        )
        ligand_summary[public_id] = {
            "status": status,
            "source_input_file": source_display or None,
            "source_input_format": source_format or None,
            "prepared_pdbqt_file": prepared_display if prepared_path.is_file() else None,
            "prepared_pdbqt_sha256": prepared_sha or None,
            "has_private_mapping": has_private,
            "resolved_by_private_mapping": resolved_private,
            "preparation_method": method or None,
            "formal_charge": metrics["formal_charge"] or None,
            "mol_weight": metrics["mol_weight"] or None,
            "num_heavy_atoms": metrics["num_heavy_atoms"] or None,
            "num_rotatable_bonds": metrics["num_rotatable_bonds"] or None,
            "stereochemistry_defined": metrics["stereochemistry_defined"] or None,
            "protonation_status_known": False,
        }

    _append_qc(qc_rows, "", "internal_id_not_leaked", "confidentiality", "PASS", "confidentiality scan pending final write", "Regenerate outputs after removing private tokens.", False)
    _append_qc(qc_rows, "", "non_goals_preserved", "scope", "PASS", "no receptor PDBQT, docking boxes, Vina, docking jobs, pose parsing, scoring, or candidate nomination", "", False)
    _safe_write_csv(prep_manifest, PREP_MANIFEST_FIELDS, manifest_rows, ctx)
    _safe_write_csv(prepared_hashes, HASH_FIELDS, hash_rows, ctx)
    _safe_write_csv(prep_qc, QC_FIELDS, qc_rows, ctx)

    leak_rows = scan_internal_id_leaks(ctx, private_entries.values())
    leak_count = len(leak_rows)
    if leak_count:
        blockers.append("internal ID leakage detected in public outputs/logs")

    counts = {
        "ligands_expected": len(PUBLIC_COMPOUND_IDS),
        "ligands_prepared": prepared_count,
        "ligands_reused_supplied_pdbqt": reused_count,
        "ligands_failed": failed_count,
        "ligands_warn": warn_count,
        "private_mappings_used": private_used_count,
        "internal_id_leaks": leak_count,
        "tracked_sensitive_files": len(tracked_sensitive),
    }
    if mode == "dry-run" and not blockers:
        warnings.append("dry-run mode only; no PDBQT files were created")
    status = "FAIL" if blockers else ("WARN" if warnings or mode == "dry-run" else "PASS")
    m3_t3_allowed = status == "PASS"

    summary = {
        "schema_version": "m3_ligand_prep_summary_v1",
        "run_id": ctx.run_id,
        "reviewed_at": now_iso(),
        "mode": mode,
        "overall_status": status,
        "m3_t3_allowed": m3_t3_allowed,
        "m3_t1": {
            "summary_found": m3_t1_summary is not None,
            "overall_status": m3_t1_summary.get("overall_status") if m3_t1_summary else None,
            "m3_t2_allowed": m3_t1_summary.get("m3_t2_allowed") if m3_t1_summary else None,
            "confidentiality_leaks": (m3_t1_summary.get("counts") or {}).get("confidentiality_leaks") if m3_t1_summary else None,
        },
        "tools": {
            "rdkit_available": tools.rdkit_available,
            "rdkit_version": tools.rdkit_version,
            "openbabel_available": tools.openbabel_available,
            "openbabel_version": tools.openbabel_version,
            "meeko_available": tools.meeko_available,
            "meeko_version": tools.meeko_version,
            "vina_invoked": False,
        },
        "no_receptor_pdbqt_generation_attempted": True,
        "no_docking_box_generation_attempted": True,
        "no_vina_attempted": True,
        "no_docking_attempted": True,
        "no_scoring_attempted": True,
        "no_candidate_nomination_attempted": True,
        "blockers": blockers,
        "warnings": warnings,
        "counts": counts,
        "ligands": ligand_summary,
        "confidentiality": {
            "internal_ids_redacted": True,
            "raw_structures_not_printed": True,
            "public_outputs_scanned": [ctx.relative_to_repo(path) for path in public_output_scan_paths(ctx)],
            "leaks_detected": leak_count,
        },
        "gitignore": {
            "updated": gitignore_updated,
            "patterns_verified": patterns_verified,
            "sensitive_files_tracked": tracked_sensitive,
        },
        "non_goals_preserved": {
            "receptor_pdbqt_generation": True,
            "docking_box_generation": True,
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
        status,
        m3_t3_allowed,
        blockers,
        warnings,
        counts,
        {
            "ligand_preparation_manifest": prep_manifest,
            "prepared_ligand_hashes": prepared_hashes,
            "ligand_prep_qc": prep_qc,
            "m3_ligand_prep_summary": summary_json,
            "m3_task2_ligand_preparation": report_md,
            "phase3_log": phase3_log,
        },
    )
    with ctx.require_within_run_dir(phase3_log).open("a", encoding="utf-8") as handle:
        handle.write(
            "{0} M3-T2 status={1} run_id={2} mode={3} m3_t3_allowed={4} ligands_prepared={5} no_receptor_pdbqt_generation_attempted=true no_docking_box_generation_attempted=true no_vina_attempted=true no_docking_attempted=true no_scoring_attempted=true\n".format(
                now_iso(),
                status,
                ctx.run_id,
                mode,
                str(m3_t3_allowed).lower(),
                prepared_count,
            )
        )
    log_master(ctx, status, "M3-T2 ligand preparation completed without receptor/Vina/docking work")
    append_phase_status(
        ctx,
        "phase3_compounds",
        status if status != "WARN" else "WARN",
        "M3-T2 ligand preparation and PDBQT generation completed",
        {
            "task": "M3-T2",
            "mode": mode,
            "m3_t3_allowed": m3_t3_allowed,
            "ligands_prepared": prepared_count,
            "no_receptor_pdbqt_generation_attempted": True,
            "no_docking_box_generation_attempted": True,
            "no_vina_attempted": True,
            "no_docking_attempted": True,
            "no_scoring_attempted": True,
        },
    )
    append_job_status(
        ctx,
        "M3-T2-ligand-preparation",
        status if status != "WARN" else "WARN",
        details={
            "message": "ligand preparation only; no receptor/Vina/docking/scoring work",
            "m3_t3_allowed": m3_t3_allowed,
        },
    )

    return M3LigandPrepResult(
        status=status,
        m3_t3_allowed=m3_t3_allowed,
        blockers=blockers,
        warnings=warnings,
        preparation_manifest_csv=prep_manifest,
        prepared_hashes_csv=prepared_hashes,
        prep_qc_csv=prep_qc,
        summary_json=summary_json,
        report_md=report_md,
        phase3_log=phase3_log,
        counts=counts,
    )
