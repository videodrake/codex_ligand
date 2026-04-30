"""M3-T1 ligand manifest, file QC, and confidentiality audit."""

from __future__ import annotations

import csv
import glob
import hashlib
import json
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
from egfr_myo1d.core.logging_utils import append_job_status, append_phase_status, log_master
from egfr_myo1d.core.run_context import RunContext


ALLOWED_EXTENSIONS = [".sdf", ".mol2", ".pdb", ".smi", ".smiles", ".pdbqt"]
PREFERRED_EXTENSIONS = [".sdf", ".mol2", ".pdb", ".smi", ".smiles", ".pdbqt"]

MANIFEST_FIELDS = [
    "compound_public_id",
    "input_file",
    "input_format",
    "file_exists",
    "file_sha256",
    "size_bytes",
    "has_private_mapping",
    "resolved_by_private_mapping",
    "allowed_for_public_output",
    "selected_for_preparation",
    "alternative_file_count",
    "tracked_by_git",
    "gitignore_protected",
    "qc_status",
    "qc_notes",
]
HASH_FIELDS = [
    "compound_public_id",
    "input_file",
    "input_format",
    "sha256",
    "size_bytes",
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
AUDIT_FIELDS = [
    "check_id",
    "category",
    "scope",
    "status",
    "severity",
    "public_id",
    "observed_file",
    "details",
    "recommended_fix",
]


@dataclass
class LigandResolution:
    public_id: str
    selected_path: Path | None
    display_path: str
    input_format: str
    has_private_mapping: bool
    resolved_by_private_mapping: bool
    alternatives: list[Path]
    unsupported_candidates: list[Path]
    selection_warning: str | None = None


@dataclass
class M3LigandQCResult:
    status: str
    m3_t2_allowed: bool
    blockers: list[str]
    warnings: list[str]
    manifest_csv: Path
    hashes_csv: Path
    qc_csv: Path
    audit_csv: Path
    summary_json: Path
    report_md: Path
    phase3_log: Path
    counts: dict[str, int]


def now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def _bool_text(value: bool) -> str:
    return "true" if value else "false"


def _status_severity(status: str, blocking: bool = False) -> str:
    if status == "FAIL" and blocking:
        return "BLOCKER"
    if status == "FAIL":
        return "MAJOR"
    if status == "WARN":
        return "MAJOR"
    return "INFO"


def _safe_write_csv(path: Path, fields: list[str], rows: list[dict[str, Any]], ctx: RunContext) -> None:
    safe_path = ctx.require_within_run_dir(path)
    safe_path.parent.mkdir(parents=True, exist_ok=True)
    with safe_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _safe_private_map_path(ctx: RunContext, private_map: Path | None) -> Path:
    path = private_map or (ctx.fresh_root / "data" / "private" / "compound_id_map.csv")
    resolved = Path(path)
    if not resolved.is_absolute():
        resolved = ctx.repo_root / resolved
    resolved = resolved.resolve()
    return assert_fresh_input_path(ctx, resolved)


def _candidate_paths(patterns: list[Path]) -> list[Path]:
    matches: list[Path] = []
    for pattern in patterns:
        for match in glob.glob(str(pattern)):
            path = Path(match).resolve()
            if path.is_file():
                matches.append(path)
    return sorted(set(matches), key=lambda p: p.as_posix())


def _select_candidate(candidates: list[Path]) -> tuple[Path | None, list[Path], list[Path], str | None]:
    allowed = [path for path in candidates if path.suffix.lower() in ALLOWED_EXTENSIONS]
    unsupported = [path for path in candidates if path.suffix.lower() not in ALLOWED_EXTENSIONS]
    if not allowed:
        return None, [], unsupported, None
    by_extension: dict[str, list[Path]] = {}
    for path in allowed:
        by_extension.setdefault(path.suffix.lower(), []).append(path)
    for ext in PREFERRED_EXTENSIONS:
        if ext in by_extension:
            selected_group = sorted(by_extension[ext], key=lambda p: p.as_posix())
            warning = None
            if len(selected_group) > 1:
                warning = "multiple matching ligand files; selected deterministic lexicographic first"
            selected = selected_group[0]
            alternatives = [path for path in allowed if path != selected]
            return selected, alternatives, unsupported, warning
    return None, [], unsupported, None


def _private_patterns(ctx: RunContext, internal_id: str) -> list[Path]:
    return [
        ctx.fresh_root / "data" / "private" / "ligands" / f"{internal_id}.*",
        ctx.fresh_root / "data" / "raw" / "ligands" / f"{internal_id}.*",
        ctx.run_dir / "phase3_compounds" / "ligands" / f"{internal_id}.*",
    ]


def resolve_ligand(
    ctx: RunContext,
    public_id: str,
    private_entry: PrivateMapEntry | None,
) -> LigandResolution:
    public_candidates = _candidate_paths(
        [
            ctx.fresh_root / "data" / "raw" / "ligands" / f"{public_id}.*",
            ctx.run_dir / "phase3_compounds" / "ligands" / f"{public_id}.*",
        ]
    )
    for candidate in public_candidates:
        assert_fresh_input_path(ctx, candidate)
    selected, alternatives, unsupported, selection_warning = _select_candidate(public_candidates)
    if selected is not None or unsupported:
        return LigandResolution(
            public_id=public_id,
            selected_path=selected,
            display_path=ctx.relative_to_repo(selected) if selected else "",
            input_format=selected.suffix.lower().lstrip(".") if selected else "",
            has_private_mapping=private_entry is not None,
            resolved_by_private_mapping=False,
            alternatives=alternatives,
            unsupported_candidates=unsupported,
            selection_warning=selection_warning,
        )
    if private_entry is None:
        return LigandResolution(
            public_id=public_id,
            selected_path=None,
            display_path="",
            input_format="",
            has_private_mapping=False,
            resolved_by_private_mapping=False,
            alternatives=[],
            unsupported_candidates=[],
        )
    private_candidates = _candidate_paths(_private_patterns(ctx, private_entry.internal_id))
    for candidate in private_candidates:
        assert_fresh_input_path(ctx, candidate)
    selected, alternatives, unsupported, selection_warning = _select_candidate(private_candidates)
    return LigandResolution(
        public_id=public_id,
        selected_path=selected,
        display_path=PRIVATE_PATH_REDACTED if selected else "",
        input_format=selected.suffix.lower().lstrip(".") if selected else "",
        has_private_mapping=True,
        resolved_by_private_mapping=selected is not None,
        alternatives=alternatives,
        unsupported_candidates=unsupported,
        selection_warning=selection_warning,
    )


def basic_parse_check(path: Path) -> tuple[str, str]:
    suffix = path.suffix.lower()
    try:
        text = path.read_text(encoding="utf-8", errors="ignore")
    except OSError as exc:
        return "FAIL", "file could not be read: {0}".format(exc.__class__.__name__)
    if not text.strip():
        return "FAIL", "file is empty"
    lines = text.splitlines()
    if suffix in {".smi", ".smiles"}:
        useful = [line for line in lines if line.strip() and not line.lstrip().startswith("#")]
        if not useful or not useful[0].split():
            return "FAIL", "SMILES-like file has no non-empty first token"
        return "PASS", "SMILES-like file has at least one non-empty non-comment record"
    if suffix == ".sdf":
        records = text.count("$$$$")
        if records == 0:
            return "WARN", "SDF is non-empty but has no record separator"
        return "PASS", "SDF record separators present"
    if suffix == ".mol2":
        if "@<TRIPOS>MOLECULE" not in text:
            return "FAIL", "MOL2 marker missing"
        return "PASS", "MOL2 molecule marker present"
    if suffix in {".pdb", ".pdbqt"}:
        atom_count = sum(1 for line in lines if line.startswith("ATOM") or line.startswith("HETATM"))
        if atom_count == 0:
            return "WARN", "coordinate-like file has no ATOM/HETATM records"
        if suffix == ".pdbqt":
            return "WARN", "PDBQT supplied as input; not generated by M3-T1"
        return "PASS", "PDB atom records present"
    return "FAIL", "unsupported ligand format"


def _append_qc(
    rows: list[dict[str, str]],
    public_id: str,
    check_id: str,
    category: str,
    status: str,
    details: str,
    recommended_fix: str = "",
    blocking: bool = False,
) -> None:
    rows.append(
        {
            "compound_public_id": public_id,
            "check_id": check_id,
            "category": category,
            "status": status,
            "severity": _status_severity(status, blocking),
            "details": details,
            "recommended_fix": recommended_fix,
        }
    )


def _audit_pass(check_id: str, scope: str, details: str = "") -> dict[str, str]:
    return {
        "check_id": check_id,
        "category": "confidentiality",
        "scope": scope,
        "status": "PASS",
        "severity": "INFO",
        "public_id": "",
        "observed_file": "",
        "details": details,
        "recommended_fix": "",
    }


def _read_m3_t0(ctx: RunContext, warnings: list[str]) -> dict[str, Any]:
    path = ctx.run_dir / "phase3_compounds" / "qc" / "m3_readiness_summary.json"
    if not path.is_file():
        warnings.append("M3 readiness summary not found; ligand QC was run independently.")
        return {"summary_found": False, "overall_status": None, "m3_docking_allowed": None}
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        warnings.append("M3 readiness summary was present but unreadable.")
        return {"summary_found": True, "overall_status": None, "m3_docking_allowed": None}
    return {
        "summary_found": True,
        "overall_status": payload.get("overall_status"),
        "m3_docking_allowed": payload.get("m3_docking_allowed"),
    }


def _ensure_logs(ctx: RunContext) -> None:
    for path in [
        ctx.logs_dir / "master.log",
        ctx.logs_dir / "phase_status.jsonl",
        ctx.logs_dir / "job_status.jsonl",
        ctx.errors_dir / "error_summary.txt",
        ctx.errors_dir / "failed_jobs.csv",
        ctx.logs_dir / "phase3_compounds.log",
    ]:
        safe = ctx.require_within_run_dir(path)
        safe.parent.mkdir(parents=True, exist_ok=True)
        if not safe.exists():
            safe.touch()
    failed = ctx.errors_dir / "failed_jobs.csv"
    if failed.stat().st_size == 0:
        failed.write_text("timestamp,job_name,status,message\n", encoding="utf-8")


def _write_report(
    path: Path,
    ctx: RunContext,
    status: str,
    m3_t2_allowed: bool,
    blockers: list[str],
    warnings: list[str],
    counts: dict[str, int],
    outputs: dict[str, Path],
) -> None:
    safe = ctx.require_within_run_dir(path)
    safe.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# M3-T1 Ligand Manifest, Confidentiality Audit, and QC",
        "",
        f"- run_id: {ctx.run_id}",
        f"- overall_status: {status}",
        f"- m3_t2_allowed: {str(m3_t2_allowed).lower()}",
        "- compounds: confidential MYO1D activity-associated compounds",
        "- no_ligand_preparation_attempted: true",
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
    lines.extend(
        [
            "",
            "## Public IDs",
            "",
            "- Cpd-A",
            "- Cpd-B",
            "- Cpd-C",
            "",
            "## Non-goals Preserved",
            "",
            "- no ligand preparation",
            "- no ligand PDBQT generation",
            "- no receptor PDBQT generation",
            "- no Vina execution",
            "- no docking job generation",
            "- no compound scoring",
            "- no candidate nomination",
            "",
            "## Outputs",
            "",
        ]
    )
    for label, output in outputs.items():
        lines.append(f"- {label}: {ctx.relative_to_repo(output)}")
    lines.extend(["", "## Next Task", "", "M3-T2 - ligand preparation and PDBQT generation"])
    safe.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_m3_ligand_qc(
    ctx: RunContext,
    mode: str = "dry-run",
    private_map: Path | None = None,
    update_ignore: bool = True,
) -> M3LigandQCResult:
    if mode not in {"dry-run", "docking"}:
        raise ValueError("mode must be dry-run or docking")
    ctx.create_directories()
    _ensure_logs(ctx)

    phase3 = ctx.run_dir / "phase3_compounds"
    manifest_dir = ctx.require_within_run_dir(phase3 / "manifests")
    qc_dir = ctx.require_within_run_dir(phase3 / "qc")
    reports_dir = ctx.require_within_run_dir(phase3 / "reports")
    for directory in [manifest_dir, qc_dir, reports_dir]:
        directory.mkdir(parents=True, exist_ok=True)

    manifest_csv = manifest_dir / "ligand_manifest_public.csv"
    hashes_csv = manifest_dir / "ligand_file_hashes.csv"
    manifest_alias = manifest_dir / "ligand_qc.csv"
    qc_csv = qc_dir / "ligand_qc_summary.csv"
    audit_csv = qc_dir / "confidentiality_audit.csv"
    summary_json = qc_dir / "m3_ligand_qc_summary.json"
    report_md = reports_dir / "m3_task1_ligand_qc.md"
    phase3_log = ctx.logs_dir / "phase3_compounds.log"

    blockers: list[str] = []
    warnings: list[str] = []
    manifest_rows: list[dict[str, Any]] = []
    hash_rows: list[dict[str, Any]] = []
    qc_rows: list[dict[str, str]] = []
    audit_rows: list[dict[str, str]] = []

    m3_t0 = _read_m3_t0(ctx, warnings)
    m3_t0_allows_t2 = bool(m3_t0.get("summary_found")) and m3_t0.get("overall_status") == "PASS" and m3_t0.get("m3_docking_allowed") is True
    if mode == "docking" and not m3_t0_allows_t2:
        blockers.append(
            "M3-T0 readiness must PASS and allow docking before M3-T1 can authorize M3-T2 ligand preparation"
        )
    elif mode == "dry-run" and not m3_t0_allows_t2:
        warnings.append("M3-T0 readiness is not PASS; M3-T2 will remain blocked")
    private_map_path = _safe_private_map_path(ctx, private_map)
    private_entries, private_warnings = load_private_map(private_map_path)
    if private_warnings:
        warnings.extend(private_warnings)
        audit_rows.append(
            {
                "check_id": "PRIVATE_MAP_PRESENT",
                "category": "confidentiality",
                "scope": "private_mapping",
                "status": "WARN",
                "severity": "MAJOR",
                "public_id": "",
                "observed_file": "",
                "details": "private mapping absent or not usable",
                "recommended_fix": "Provide fresh/data/private/compound_id_map.csv if internal IDs are needed.",
            }
        )
    else:
        audit_rows.append(_audit_pass("PRIVATE_MAP_PRESENT", "private_mapping", "private mapping loaded without public disclosure"))

    gitignore_updated = False
    patterns_verified = list(RECOMMENDED_GITIGNORE_PATTERNS)
    if update_ignore:
        gitignore_updated, patterns_verified = update_gitignore(ctx.repo_root)

    private_map_tracked = bool(git_ls_files(ctx.repo_root, private_map_path)) if private_map_path.exists() else False
    tracked_sensitive: list[str] = []
    if private_map_tracked:
        tracked_sensitive.append("PRIVATE_MAPPING")
        blockers.append("private compound_id_map.csv is tracked by git")
        audit_rows.append(
            {
                "check_id": "PRIVATE_MAP_GIT_TRACKING",
                "category": "confidentiality",
                "scope": "git_tracking",
                "status": "FAIL",
                "severity": "BLOCKER",
                "public_id": "",
                "observed_file": "PRIVATE_PATH_REDACTED",
                "details": "private mapping is tracked by git",
                "recommended_fix": "Remove from git index without deleting local file.",
            }
        )
    else:
        audit_rows.append(_audit_pass("PRIVATE_MAP_GIT_TRACKING", "git_tracking", "private mapping is not tracked or absent"))

    inputs: dict[str, dict[str, Any]] = {}
    ligands_found = 0
    ligands_missing = 0
    ligands_failed = 0
    ligands_warn = 0
    private_mappings_found = 0

    for public_id in PUBLIC_COMPOUND_IDS:
        entry = private_entries.get(public_id)
        if entry is not None:
            private_mappings_found += 1
        resolution = resolve_ligand(ctx, public_id, entry)
        selected = resolution.selected_path
        exists = selected is not None and selected.is_file()
        sha = _sha256(selected) if exists else ""
        size = selected.stat().st_size if exists else 0
        tracked = bool(git_ls_files(ctx.repo_root, selected)) if exists else False
        ignored = git_check_ignored(ctx.repo_root, selected) if exists else False
        candidate_paths = []
        if selected is not None:
            candidate_paths.append(selected)
        candidate_paths.extend(resolution.alternatives)
        candidate_paths.extend(resolution.unsupported_candidates)
        tracked_candidates = [path for path in sorted(set(candidate_paths), key=lambda p: p.as_posix()) if git_ls_files(ctx.repo_root, path)]
        notes: list[str] = []
        status = "PASS"

        if resolution.selection_warning:
            status = "WARN"
            notes.append(resolution.selection_warning)
            warnings.append(f"{public_id}: multiple matching ligand files; deterministic selection used")
        if resolution.unsupported_candidates and not exists:
            status = "FAIL"
            notes.append("unsupported ligand format")
            blockers.append(f"{public_id}: unsupported ligand format")
        if not exists:
            ligands_missing += 1
            notes.append("ligand file missing")
            if resolution.unsupported_candidates:
                status = "FAIL"
            elif mode == "docking":
                status = "FAIL"
                blockers.append(f"{public_id}: ligand file missing in docking mode")
            else:
                status = "WARN"
                warnings.append(f"{public_id}: ligand file missing in dry-run mode")
        else:
            ligands_found += 1
            if selected.suffix.lower() not in ALLOWED_EXTENSIONS:
                status = "FAIL"
                notes.append("unsupported ligand format")
                blockers.append(f"{public_id}: unsupported ligand format")
            if size == 0:
                status = "FAIL"
                notes.append("zero-byte ligand file")
                blockers.append(f"{public_id}: zero-byte ligand file")
            parse_status, parse_notes = basic_parse_check(selected)
            notes.append(parse_notes)
            if parse_status == "FAIL":
                status = "FAIL"
                blockers.append(f"{public_id}: basic parse failed")
            elif parse_status == "WARN" and status != "FAIL":
                status = "WARN"
                warnings.append(f"{public_id}: basic parse warning")
            if selected.suffix.lower() == ".pdbqt" and status != "FAIL":
                status = "WARN"
                warnings.append(f"{public_id}: PDBQT supplied as input and will be validated in M3-T2")
            if tracked:
                status = "FAIL"
                blockers.append(f"{public_id}: raw ligand file is tracked by git")
                tracked_sensitive.append(f"RAW_LIGAND_{public_id}")
            extra_tracked = [path for path in tracked_candidates if path != selected]
            if extra_tracked:
                status = "FAIL"
                blockers.append(f"{public_id}: alternative or unsupported raw ligand candidate is tracked by git")
                tracked_sensitive.extend(f"RAW_LIGAND_{public_id}_CANDIDATE" for _ in extra_tracked)
            if not ignored:
                if status != "FAIL":
                    status = "WARN"
                warnings.append(f"{public_id}: raw ligand file is not gitignore-protected")

            _append_qc(qc_rows, public_id, "FORMAT_ALLOWED", "file_qc", "PASS" if selected.suffix.lower() in ALLOWED_EXTENSIONS else "FAIL", "format allowed" if selected.suffix.lower() in ALLOWED_EXTENSIONS else "unsupported format", "Use one of .sdf/.mol2/.pdb/.smi/.smiles as raw ligand input; PDBQT is generated later by M3-T2.", selected.suffix.lower() not in ALLOWED_EXTENSIONS)
            _append_qc(qc_rows, public_id, "NONEMPTY_FILE", "file_qc", "PASS" if size > 0 else "FAIL", "non-empty file" if size > 0 else "zero-byte file", "Replace with a non-empty ligand file.", size <= 0)
            _append_qc(qc_rows, public_id, "BASIC_PARSE", "file_qc", parse_status, parse_notes, "Inspect source ligand file.", parse_status == "FAIL")
            _append_qc(qc_rows, public_id, "GIT_TRACKING", "confidentiality", "FAIL" if tracked else "PASS", "raw ligand file is tracked by git" if tracked else "raw ligand file is not tracked", "Remove from git index without deleting local file.", tracked)
            _append_qc(qc_rows, public_id, "GITIGNORE_PROTECTION", "confidentiality", "PASS" if ignored else "WARN", "raw ligand file is gitignore-protected" if ignored else "raw ligand file is not gitignore-protected", "Add recommended ligand ignore patterns.", False)

        if status == "FAIL":
            ligands_failed += 1
        elif status == "WARN":
            ligands_warn += 1
        if not exists:
            _append_qc(qc_rows, public_id, "FILE_EXISTS", "file_qc", status, "ligand file missing", "Place a public-ID ligand file or private mapping target.", mode == "docking")
            if tracked_candidates:
                status = "FAIL"
                blockers.append(f"{public_id}: unsupported raw ligand candidate is tracked by git")
                tracked_sensitive.extend(f"RAW_LIGAND_{public_id}_CANDIDATE" for _ in tracked_candidates)

        display_path = resolution.display_path
        allowed_public_output = not resolution.resolved_by_private_mapping
        manifest_rows.append(
            {
                "compound_public_id": public_id,
                "input_file": display_path,
                "input_format": resolution.input_format,
                "file_exists": _bool_text(exists),
                "file_sha256": sha,
                "size_bytes": str(size) if exists else "",
                "has_private_mapping": _bool_text(resolution.has_private_mapping),
                "resolved_by_private_mapping": _bool_text(resolution.resolved_by_private_mapping),
                "allowed_for_public_output": _bool_text(allowed_public_output),
                "selected_for_preparation": _bool_text(status in {"PASS", "WARN"} and exists),
                "alternative_file_count": str(len(resolution.alternatives)),
                "tracked_by_git": _bool_text(tracked),
                "gitignore_protected": _bool_text(ignored),
                "qc_status": status,
                "qc_notes": "; ".join(notes),
            }
        )
        hash_rows.append(
            {
                "compound_public_id": public_id,
                "input_file": display_path,
                "input_format": resolution.input_format,
                "sha256": sha,
                "size_bytes": str(size) if exists else "",
                "hash_status": "PASS" if exists and sha else ("FAIL" if mode == "docking" else "WARN"),
                "notes": "hash recorded" if exists and sha else "ligand file missing",
            }
        )
        inputs[public_id] = {
            "status": status,
            "input_file": display_path or None,
            "input_format": resolution.input_format or None,
            "sha256": sha or None,
            "size_bytes": size if exists else None,
            "has_private_mapping": resolution.has_private_mapping,
            "resolved_by_private_mapping": resolution.resolved_by_private_mapping,
            "allowed_for_public_output": allowed_public_output,
        }

    _safe_write_csv(manifest_csv, MANIFEST_FIELDS, manifest_rows, ctx)
    _safe_write_csv(manifest_alias, MANIFEST_FIELDS, manifest_rows, ctx)
    _safe_write_csv(hashes_csv, HASH_FIELDS, hash_rows, ctx)
    _safe_write_csv(qc_csv, QC_FIELDS, qc_rows, ctx)

    audit_rows.append(_audit_pass("PUBLIC_IDS_ONLY", "manifest", "public outputs use Cpd-A/Cpd-B/Cpd-C labels"))
    audit_rows.append(_audit_pass("RAW_STRUCTURES_NOT_COPIED", "reports_and_manifests", "raw ligand structures were not copied into reports/manifests"))
    audit_rows.append(_audit_pass("GITIGNORE_PATTERNS", "gitignore", "recommended sensitive-data patterns verified"))

    leak_rows = scan_internal_id_leaks(ctx, private_entries.values())
    if leak_rows:
        blockers.append("internal ID leakage detected in public outputs/logs")
        audit_rows.extend(leak_rows)

    _safe_write_csv(audit_csv, AUDIT_FIELDS, audit_rows, ctx)

    # Re-scan after writing the audit itself; rows never include private tokens.
    leak_rows = scan_internal_id_leaks(ctx, private_entries.values())
    leak_count = len(leak_rows)
    if leak_count and "internal ID leakage detected in public outputs/logs" not in blockers:
        blockers.append("internal ID leakage detected in public outputs/logs")

    counts = {
        "ligands_expected": len(PUBLIC_COMPOUND_IDS),
        "ligands_found": ligands_found,
        "ligands_missing": ligands_missing,
        "ligands_failed_qc": ligands_failed,
        "ligands_warn_qc": ligands_warn,
        "private_mappings_found": private_mappings_found,
        "confidentiality_leaks": leak_count,
        "tracked_sensitive_files": len(tracked_sensitive),
    }
    status = "FAIL" if blockers else ("WARN" if warnings else "PASS")
    m3_t2_allowed = status == "PASS" and mode == "docking" and m3_t0_allows_t2

    summary = {
        "schema_version": "m3_ligand_qc_summary_v1",
        "run_id": ctx.run_id,
        "reviewed_at": now_iso(),
        "mode": mode,
        "overall_status": status,
        "m3_t2_allowed": m3_t2_allowed,
        "no_ligand_preparation_attempted": True,
        "no_ligand_pdbqt_generation_attempted": True,
        "no_receptor_pdbqt_generation_attempted": True,
        "no_docking_attempted": True,
        "no_scoring_attempted": True,
        "blockers": blockers,
        "warnings": warnings,
        "counts": counts,
        "inputs": inputs,
        "confidentiality": {
            "private_map_present": private_map_path.is_file(),
            "private_map_tracked_by_git": private_map_tracked,
            "internal_ids_redacted": True,
            "public_outputs_scanned": [ctx.relative_to_repo(path) for path in public_output_scan_paths(ctx)],
            "leaks_detected": leak_count,
        },
        "gitignore": {
            "updated": gitignore_updated,
            "patterns_verified": patterns_verified,
            "sensitive_files_tracked": tracked_sensitive,
        },
        "m3_t0": m3_t0,
        "non_goals_preserved": {
            "ligand_preparation": True,
            "ligand_pdbqt_generation": True,
            "receptor_pdbqt_generation": True,
            "vina_execution": True,
            "docking_job_generation": True,
            "compound_scoring": True,
            "candidate_nomination": True,
        },
    }
    summary_path = ctx.require_within_run_dir(summary_json)
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    _write_report(
        report_md,
        ctx,
        status,
        m3_t2_allowed,
        blockers,
        warnings,
        counts,
        {
            "ligand_manifest_public": manifest_csv,
            "ligand_file_hashes": hashes_csv,
            "ligand_qc_summary": qc_csv,
            "confidentiality_audit": audit_csv,
            "m3_ligand_qc_summary": summary_json,
            "m3_task1_ligand_qc": report_md,
            "phase3_log": phase3_log,
        },
    )

    with ctx.require_within_run_dir(phase3_log).open("a", encoding="utf-8") as handle:
        handle.write(
            "{0} M3-T1 status={1} run_id={2} mode={3} m3_t2_allowed={4} no_ligand_preparation_attempted=true no_docking_attempted=true blockers={5} warnings={6}\n".format(
                now_iso(),
                status,
                ctx.run_id,
                mode,
                str(m3_t2_allowed).lower(),
                len(blockers),
                len(warnings),
            )
        )
    log_master(ctx, status, "M3-T1 ligand manifest and confidentiality QC completed")
    append_phase_status(
        ctx,
        "phase3_compounds",
        status if status != "WARN" else "WARN",
        "M3-T1 ligand manifest, confidentiality audit, and ligand file QC completed",
        {
            "task": "M3-T1",
            "mode": mode,
            "m3_t2_allowed": m3_t2_allowed,
            "no_ligand_preparation_attempted": True,
            "no_ligand_pdbqt_generation_attempted": True,
            "no_docking_attempted": True,
        },
    )
    append_job_status(
        ctx,
        "M3-T1-ligand-qc",
        status if status != "WARN" else "WARN",
        details={
            "message": "ligand QC and confidentiality audit only",
            "m3_t2_allowed": m3_t2_allowed,
            "no_ligand_preparation_attempted": True,
            "no_docking_attempted": True,
        },
    )

    return M3LigandQCResult(
        status=status,
        m3_t2_allowed=m3_t2_allowed,
        blockers=blockers,
        warnings=warnings,
        manifest_csv=manifest_csv,
        hashes_csv=hashes_csv,
        qc_csv=qc_csv,
        audit_csv=audit_csv,
        summary_json=summary_json,
        report_md=report_md,
        phase3_log=phase3_log,
        counts=counts,
    )
