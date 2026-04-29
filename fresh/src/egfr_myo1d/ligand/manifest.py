"""M1 Phase 7: Ligand manifest shell per handoff §17.

For each public ligand ID (Cpd-A/B/C from fresh_run.yaml):
  - Look for <ligands_dir>/<public_id>.{sdf,mol,mol2,pdb}
  - Record sha256 + format if present, or missing-status if absent
  - Write qc/ligand_manifest_qc.csv (public IDs only)
  - Read private_mapping (fresh/data/private/compound_id_map.csv) if present;
    validate the 3-column schema (public_id, internal_id, notes)
  - Verify NO internal IDs from the private mapping leak into any run output
    file. A leak immediately escalates status to FAIL.

M1 does NOT run ligand docking. This module is shell-only.
"""

from __future__ import annotations

import csv
import hashlib
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

from egfr_myo1d.core.logging_utils import append_phase_status
from egfr_myo1d.core.manifest import load_yaml_config, now_iso, write_json
from egfr_myo1d.core.run_context import RunContext, RunContextError
from egfr_myo1d.io.hashing import describe_file


SUPPORTED_FORMATS = ("sdf", "mol", "mol2", "pdb")
PRIVATE_MAPPING_REQUIRED_COLUMNS = ("public_id", "internal_id", "notes")

LIGAND_MANIFEST_QC_COLUMNS = [
    "public_id",
    "expected_path",
    "exists",
    "format",
    "sha256_or_missing",
    "size_bytes_or_missing",
    "status",
    "notes",
]


@dataclass
class LigandFileRecord:
    public_id: str
    expected_path: Path
    exists: bool = False
    format: str = ""
    sha256: str | None = None
    size_bytes: int | None = None
    alternative_formats: list[str] = field(default_factory=list)
    status: str = "PASS"
    notes: str = ""

    def to_csv_row(self) -> dict[str, Any]:
        return {
            "public_id": self.public_id,
            "expected_path": str(self.expected_path),
            "exists": str(self.exists).lower(),
            "format": self.format,
            "sha256_or_missing": self.sha256 if self.sha256 else "MISSING",
            "size_bytes_or_missing": (
                self.size_bytes if self.size_bytes is not None else "MISSING"
            ),
            "status": self.status,
            "notes": self.notes,
        }


@dataclass
class LigandManifest:
    public_ids: list[str]
    ligand_files: list[LigandFileRecord] = field(default_factory=list)
    present_count: int = 0
    missing_count: int = 0
    private_mapping_path: Path | None = None
    private_mapping_present: bool = False
    private_mapping_schema_valid: bool | None = None
    internal_ids_leaked_into_outputs: bool = False
    leaked_internal_id_count: int = 0
    leaked_internal_id_hashes: list[str] = field(default_factory=list)
    leaked_internal_ids: list[str] = field(default_factory=list)
    compound_stage_enabled: bool = False
    profile: str = "codex_dev"
    status: str = "PASS"
    warnings: list[str] = field(default_factory=list)
    output_qc_csv: Path | None = None
    output_manifest_json: Path | None = None


# ---------------------------------------------------------------------------
# Config + private-mapping helpers
# ---------------------------------------------------------------------------


def load_public_ids(ctx: RunContext) -> list[str]:
    config = load_yaml_config(ctx.repo_root / "fresh" / "configs" / "fresh_run.yaml")
    ligands = config.get("ligands", {})
    ids = ligands.get("public_ids") or config.get("compound_public_ids") or []
    return [str(x) for x in ids]


def load_default_paths(ctx: RunContext) -> tuple[Path, Path]:
    config = load_yaml_config(ctx.repo_root / "fresh" / "configs" / "paths.yaml")
    raw_ligands = ctx.repo_root / config.get("raw_ligands", "fresh/data/raw/ligands")
    private_dir = ctx.repo_root / config.get("private_data", "fresh/data/private")
    private_mapping = private_dir / "compound_id_map.csv"
    return raw_ligands, private_mapping


def _read_private_mapping(path: Path) -> tuple[list[dict[str, str]], bool, list[str]]:
    """Read and validate the private mapping CSV. Returns (rows, schema_valid,
    warnings). The mapping content (internal_ids) is used only for leak
    detection — it is never copied into any run output."""
    warnings: list[str] = []
    if not path.is_file():
        return [], False, ["private_mapping_not_present_or_unreadable"]
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            cols = reader.fieldnames or []
            schema_valid = all(c in cols for c in PRIVATE_MAPPING_REQUIRED_COLUMNS)
            rows = list(reader)
    except (OSError, csv.Error) as exc:
        return [], False, ["private_mapping_read_error: {0}".format(exc)]
    if not schema_valid:
        warnings.append(
            "private_mapping_schema_invalid: required={0} got={1}".format(
                list(PRIVATE_MAPPING_REQUIRED_COLUMNS), cols
            )
        )
    return rows, schema_valid, warnings


def _detect_format(public_id: str, ligands_dir: Path) -> tuple[Path | None, str, list[str]]:
    """Return (path, format, alternative_formats). ``path`` is None if no
    file matches any supported format. Format priority: sdf > mol > mol2 > pdb."""
    found: list[tuple[str, Path]] = []
    for fmt in SUPPORTED_FORMATS:
        candidate = ligands_dir / "{0}.{1}".format(public_id, fmt)
        if candidate.is_file():
            found.append((fmt, candidate))
    if not found:
        return None, "", []
    primary_fmt, primary_path = found[0]
    alts = [f for f, _ in found[1:]]
    return primary_path, primary_fmt, alts


# ---------------------------------------------------------------------------
# Leak detection
# ---------------------------------------------------------------------------


def _scan_for_internal_id_leaks(
    output_files: Iterable[Path], internal_ids: Iterable[str]
) -> list[str]:
    """Read each output file and report any internal_id substring matches.

    The mapping is read once and pulled out of the comparison if empty. We
    treat any substring match as a leak (conservative — better to false-flag
    than miss)."""
    sentinel_ids = [s for s in internal_ids if s]
    if not sentinel_ids:
        return []
    leaks: set[str] = set()
    for path in output_files:
        if not Path(path).is_file():
            continue
        try:
            text = Path(path).read_text(encoding="utf-8")
        except (OSError, UnicodeDecodeError):
            continue
        for sentinel in sentinel_ids:
            if sentinel in text:
                leaks.add(sentinel)
    return sorted(leaks)


def _redact_internal_id_leaks(leaks: Iterable[str]) -> tuple[list[str], list[str]]:
    """Return redacted placeholders and stable hashes for leaked IDs.

    Raw internal IDs must not be copied back into run outputs, including the
    failure report that records a leak. Hashes make repeated leaks auditable
    without exposing the confidential identifier itself.
    """
    hashes = []
    redacted = []
    for value in sorted(set(str(x) for x in leaks if x)):
        digest = hashlib.sha256(value.encode("utf-8")).hexdigest()
        hashes.append(digest)
        redacted.append("<redacted_internal_id_sha256:{0}>".format(digest[:12]))
    return redacted, hashes


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------


def _write_qc_csv(ctx: RunContext, records: list[LigandFileRecord]) -> Path:
    path = ctx.qc_dir / "ligand_manifest_qc.csv"
    safe = ctx.require_within_run_dir(path)
    safe.parent.mkdir(parents=True, exist_ok=True)
    with safe.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=LIGAND_MANIFEST_QC_COLUMNS, extrasaction="ignore"
        )
        writer.writeheader()
        for record in records:
            writer.writerow(record.to_csv_row())
    return safe


def _write_manifest_json(ctx: RunContext, manifest: LigandManifest) -> Path:
    payload: dict[str, Any] = {
        "public_ids": manifest.public_ids,
        "ligands_dir": str(manifest.ligand_files[0].expected_path.parent)
        if manifest.ligand_files
        else None,
        "private_mapping_path": str(manifest.private_mapping_path)
        if manifest.private_mapping_path
        else None,
        "private_mapping_present": manifest.private_mapping_present,
        "private_mapping_schema_valid": manifest.private_mapping_schema_valid,
        "internal_ids_leaked_into_outputs": manifest.internal_ids_leaked_into_outputs,
        "leaked_internal_id_count": manifest.leaked_internal_id_count,
        "leaked_internal_id_hashes": manifest.leaked_internal_id_hashes,
        "leaked_internal_ids": manifest.leaked_internal_ids,
        "compound_stage_enabled": manifest.compound_stage_enabled,
        "ligands": [
            {
                "public_id": r.public_id,
                "expected_path": str(r.expected_path),
                "exists": r.exists,
                "format": r.format,
                "sha256": r.sha256,
                "size_bytes": r.size_bytes,
                "alternative_formats": r.alternative_formats,
                "status": r.status,
                "notes": r.notes,
            }
            for r in manifest.ligand_files
        ],
        "present_count": manifest.present_count,
        "missing_count": manifest.missing_count,
        "status": manifest.status,
        "warnings": manifest.warnings,
        "profile": manifest.profile,
        "score_bonus_allowed": False,
        "timestamp": now_iso(),
    }
    target = ctx.manifest_dir / "ligand_manifest_report.json"
    write_json(target, payload, ctx)
    return target


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------


def build_ligand_manifest(
    ctx: RunContext,
    public_ids: Iterable[str] | None = None,
    ligands_dir: Path | None = None,
    private_mapping_path: Path | None = None,
    profile: str = "codex_dev",
    compound_stage_enabled: bool = False,
) -> LigandManifest:
    """Build the ligand manifest shell. M1 does NOT run docking; this only
    inventories what's expected and detects internal-ID leaks."""
    if profile not in {"codex_dev", "hpc_strict"}:
        raise ValueError("profile must be codex_dev or hpc_strict; got: {0}".format(profile))

    default_ligands_dir, default_private_mapping = load_default_paths(ctx)
    if ligands_dir is None:
        ligands_dir = default_ligands_dir
    else:
        ligands_dir = Path(ligands_dir)
    if private_mapping_path is None:
        private_mapping_path = default_private_mapping
    else:
        private_mapping_path = Path(private_mapping_path)

    public_ids_resolved = list(public_ids) if public_ids is not None else load_public_ids(ctx)
    if not public_ids_resolved:
        raise ValueError("no public_ids resolved from arguments or fresh_run.yaml")

    manifest = LigandManifest(
        public_ids=public_ids_resolved,
        private_mapping_path=private_mapping_path,
        compound_stage_enabled=compound_stage_enabled,
        profile=profile,
    )

    # Per-public-id file record
    present_count = 0
    missing_count = 0
    for public_id in public_ids_resolved:
        path, fmt, alts = _detect_format(public_id, ligands_dir)
        record = LigandFileRecord(
            public_id=public_id,
            expected_path=ligands_dir / "{0}.sdf".format(public_id),
            alternative_formats=alts,
        )
        if path is not None:
            descriptor = describe_file(path)
            record.exists = True
            record.format = fmt
            record.sha256 = descriptor.get("sha256")
            record.size_bytes = descriptor.get("size_bytes")
            record.expected_path = path
            record.status = "PASS"
            present_count += 1
            if alts:
                record.notes = "alternative_formats={0}".format(",".join(alts))
        else:
            record.exists = False
            record.format = ""
            record.notes = "missing"
            missing_count += 1
            # Severity per matrix:
            #   codex_dev / stage=false → WARN
            #   codex_dev / stage=true  → WARN
            #   hpc_strict / stage=false → WARN
            #   hpc_strict / stage=true  → FAIL
            if profile == "hpc_strict" and compound_stage_enabled:
                record.status = "FAIL"
            else:
                record.status = "WARN"
        manifest.ligand_files.append(record)

    manifest.present_count = present_count
    manifest.missing_count = missing_count

    # Aggregate ligand-file status to manifest status
    if any(r.status == "FAIL" for r in manifest.ligand_files):
        manifest.status = "FAIL"
    elif any(r.status == "WARN" for r in manifest.ligand_files):
        manifest.status = "WARN"

    # Read private mapping (only for leak detection; never copied to outputs)
    rows, schema_valid, mapping_warnings = _read_private_mapping(private_mapping_path)
    manifest.private_mapping_present = bool(rows) or private_mapping_path.is_file()
    manifest.private_mapping_schema_valid = (
        schema_valid if private_mapping_path.is_file() else None
    )
    manifest.warnings.extend(mapping_warnings)
    if mapping_warnings and manifest.status == "PASS":
        manifest.status = "WARN"

    internal_ids = [
        str(r.get("internal_id", "")).strip() for r in rows if r.get("internal_id")
    ]

    # Write outputs first so leak detection can scan them
    qc_path = _write_qc_csv(ctx, manifest.ligand_files)
    json_path = _write_manifest_json(ctx, manifest)
    manifest.output_qc_csv = qc_path
    manifest.output_manifest_json = json_path

    # Leak detection scans the outputs JUST written
    leaks = _scan_for_internal_id_leaks([qc_path, json_path], internal_ids)
    if leaks:
        redacted_leaks, leak_hashes = _redact_internal_id_leaks(leaks)
        manifest.internal_ids_leaked_into_outputs = True
        manifest.leaked_internal_id_count = len(leaks)
        manifest.leaked_internal_id_hashes = leak_hashes
        manifest.leaked_internal_ids = redacted_leaks
        manifest.status = "FAIL"
        manifest.warnings.append(
            "internal_id_leak_detected: count={0} hashes={1}".format(
                len(leaks), ",".join(leak_hashes)
            )
        )
        # Re-write JSON to record the leak in the saved report
        json_path = _write_manifest_json(ctx, manifest)
        manifest.output_manifest_json = json_path

    append_phase_status(
        ctx,
        phase="manifest-ligands",
        status=manifest.status,
        message=(
            "manifest-ligands status={0} present={1} missing={2} leak={3}".format(
                manifest.status,
                manifest.present_count,
                manifest.missing_count,
                manifest.internal_ids_leaked_into_outputs,
            )
        ),
        details={
            "public_ids": manifest.public_ids,
            "private_mapping_present": manifest.private_mapping_present,
            "private_mapping_schema_valid": manifest.private_mapping_schema_valid,
            "internal_ids_leaked_into_outputs": manifest.internal_ids_leaked_into_outputs,
            "compound_stage_enabled": manifest.compound_stage_enabled,
            "profile": manifest.profile,
            "score_bonus_allowed": False,
            "ligand_manifest_qc_csv": str(qc_path),
            "ligand_manifest_report_json": str(manifest.output_manifest_json),
        },
    )

    return manifest
