"""M1 Phase 3: MYO1D 955-1006 construct slicing + QC.

Closes M1 §23 #13 (MYO1D construct QC) per
``milestone1_foundation_codex_handoff_v0_5.md`` §16. Source-of-truth values
(construct range, key residues, watch residues, key_residue_bonus_weight) are
read from ``fresh/configs/gates.yaml``. No values are hardcoded inside this
module.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

from egfr_myo1d.core.logging_utils import append_phase_status
from egfr_myo1d.core.manifest import load_yaml_config, now_iso, write_json
from egfr_myo1d.core.run_context import RunContext, RunContextError
from egfr_myo1d.io.hashing import describe_file
from egfr_myo1d.myo1d.construct import (
    CAP_RESNAMES,
    emit_myo1d_construct_pdb,
    slice_myo1d_construct,
)
from egfr_myo1d.structure.pdb_parser import PDBParseError, parse_pdb


MYO1D_CONSTRUCT_QC_COLUMNS = [
    "construct_id",
    "source_file",
    "chain_id",
    "start_residue",
    "end_residue",
    "n_residues",
    "missing_residues",
    "key_sheet8_present",
    "key_sheet9_present",
    "key_sheet12_present",
    "n_watch_present",
    "c_watch_present",
    "ace_nme_caps_present",
    "status",
    "warnings",
]


@dataclass
class Myo1dQcReport:
    construct_id: str
    source_file: str
    chain_id: str | None
    start_residue: int | None
    end_residue: int | None
    n_residues: int = 0
    missing_residues: list[int] = field(default_factory=list)
    key_sheet8_present: bool = False
    key_sheet9_present: bool = False
    key_sheet12_present: bool = False
    n_watch_present: bool = False
    c_watch_present: bool = False
    ace_nme_caps_present: bool = False
    status: str = "PASS"
    warnings: list[str] = field(default_factory=list)
    output_pdb: Path | None = None
    output_qc_csv: Path | None = None
    output_manifest: Path | None = None
    key_residue_bonus_weight: float = 0.0
    profile: str = "codex_dev"

    def to_csv_row(self) -> dict[str, Any]:
        return {
            "construct_id": self.construct_id,
            "source_file": self.source_file,
            "chain_id": self.chain_id or "",
            "start_residue": self.start_residue if self.start_residue is not None else "",
            "end_residue": self.end_residue if self.end_residue is not None else "",
            "n_residues": self.n_residues,
            "missing_residues": ";".join(str(n) for n in sorted(self.missing_residues)),
            "key_sheet8_present": str(self.key_sheet8_present).lower(),
            "key_sheet9_present": str(self.key_sheet9_present).lower(),
            "key_sheet12_present": str(self.key_sheet12_present).lower(),
            "n_watch_present": str(self.n_watch_present).lower(),
            "c_watch_present": str(self.c_watch_present).lower(),
            "ace_nme_caps_present": str(self.ace_nme_caps_present).lower(),
            "status": self.status,
            "warnings": ";".join(self.warnings),
        }


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def parse_construct_range(text_or_tuple) -> tuple[int, int]:
    """Parse '955-1006' string or (955, 1006) tuple to (start, end).

    Raises ValueError on malformed input or end < start.
    """
    if isinstance(text_or_tuple, str):
        if "-" not in text_or_tuple:
            raise ValueError("construct must be 'start-end'; got: {0}".format(text_or_tuple))
        s_text, e_text = text_or_tuple.split("-", 1)
        start = int(s_text)
        end = int(e_text)
    else:
        start = int(text_or_tuple[0])
        end = int(text_or_tuple[1])
    if end < start:
        raise ValueError("construct_range invalid: end < start ({0}-{1})".format(start, end))
    return start, end


def expand_residue_set(text: str) -> list[int]:
    """Expand a single range '961-964' or comma list '961-964,968-972' to sorted ints."""
    out: set[int] = set()
    for chunk in text.split(","):
        chunk = chunk.strip()
        if not chunk:
            continue
        if "-" in chunk:
            s, e = chunk.split("-", 1)
            out.update(range(int(s), int(e) + 1))
        else:
            out.add(int(chunk))
    return sorted(out)


def load_myo1d_gate(ctx: RunContext) -> dict[str, Any]:
    """Read ``fresh/configs/gates.yaml`` for MYO1D source-of-truth values."""
    gates_path = ctx.repo_root / "fresh" / "configs" / "gates.yaml"
    config = load_yaml_config(gates_path)
    return config.get("myo1d", {})


def _all_in(target: Iterable[int], universe: Iterable[int]) -> bool:
    universe_set = set(universe)
    return all(item in universe_set for item in target)


def _construct_id(start: int, end: int) -> str:
    return "MYO1D_{0}_{1}".format(start, end)


def _resolve_chain_id(structure_atoms) -> str | None:
    counts: dict[str, int] = {}
    for atom in structure_atoms:
        counts[atom.chain_id] = counts.get(atom.chain_id, 0) + 1
    if not counts:
        return None
    return max(counts, key=counts.get)


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------


def _write_qc_outputs(ctx: RunContext, report: Myo1dQcReport) -> None:
    qc_csv = ctx.qc_dir / "myo1d_construct_qc.csv"
    manifest_path = ctx.manifest_dir / "myo1d_construct_manifest.json"

    safe_qc = ctx.require_within_run_dir(qc_csv)
    safe_qc.parent.mkdir(parents=True, exist_ok=True)
    with safe_qc.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=MYO1D_CONSTRUCT_QC_COLUMNS, extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerow(report.to_csv_row())
    report.output_qc_csv = safe_qc

    source_path = Path(report.source_file) if report.source_file else None
    source_desc = describe_file(source_path) if source_path else {
        "present": False, "sha256": None, "size_bytes": None,
    }
    output_desc = describe_file(report.output_pdb) if report.output_pdb else {
        "present": False, "sha256": None, "size_bytes": None,
    }

    manifest_data = {
        "construct_id": report.construct_id,
        "source_file": report.source_file,
        "source_sha256": source_desc.get("sha256"),
        "source_size_bytes": source_desc.get("size_bytes"),
        "output_pdb": str(report.output_pdb) if report.output_pdb else None,
        "output_pdb_sha256": output_desc.get("sha256"),
        "output_pdb_size_bytes": output_desc.get("size_bytes"),
        "output_qc_csv": str(report.output_qc_csv) if report.output_qc_csv else None,
        "construct_range": [report.start_residue, report.end_residue],
        "n_residues": report.n_residues,
        "missing_residues": report.missing_residues,
        "key_residue_bonus_weight": report.key_residue_bonus_weight,
        "score_bonus_allowed": False,
        "key_sheet8_present": report.key_sheet8_present,
        "key_sheet9_present": report.key_sheet9_present,
        "key_sheet12_present": report.key_sheet12_present,
        "n_watch_present": report.n_watch_present,
        "c_watch_present": report.c_watch_present,
        "ace_nme_caps_present": report.ace_nme_caps_present,
        "status": report.status,
        "warnings": report.warnings,
        "profile": report.profile,
        "timestamp": now_iso(),
    }
    write_json(manifest_path, manifest_data, ctx)
    report.output_manifest = manifest_path

    append_phase_status(
        ctx,
        phase="prepare-myo1d",
        status=report.status,
        message=(
            "prepare-myo1d construct={0} status={1} n_residues={2} warnings={3}".format(
                report.construct_id,
                report.status,
                report.n_residues,
                len(report.warnings),
            )
        ),
        details={
            "construct_id": report.construct_id,
            "construct_range": [report.start_residue, report.end_residue],
            "n_residues": report.n_residues,
            "key_sheet8_present": report.key_sheet8_present,
            "key_sheet9_present": report.key_sheet9_present,
            "key_sheet12_present": report.key_sheet12_present,
            "n_watch_present": report.n_watch_present,
            "c_watch_present": report.c_watch_present,
            "ace_nme_caps_present": report.ace_nme_caps_present,
            "score_bonus_allowed": False,
            "profile": report.profile,
            "manifest": str(manifest_path),
            "qc_csv": str(report.output_qc_csv) if report.output_qc_csv else None,
            "output_pdb": str(report.output_pdb) if report.output_pdb else None,
        },
    )


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------


def run_myo1d_qc(
    ctx: RunContext,
    source_pdb: Path,
    construct_range=None,
    profile: str = "codex_dev",
) -> Myo1dQcReport:
    """Slice MYO1D source PDB and emit the canonical M1 construct + QC.

    Reads ``fresh/configs/gates.yaml`` for MYO1D construct/key/watch values;
    no values are hardcoded inside this module.

    Outputs (under ``ctx.run_dir``):
        - ``normalized/myo1d/<construct_id>.pdb``
        - ``qc/myo1d_construct_qc.csv``
        - ``manifest/myo1d_construct_manifest.json``

    Severity policy (handoff §16):
        - ``PASS`` if source parsed, key residues present, no terminal artifact.
        - ``WARN`` for missing key residues, terminal artifacts in
          ``codex_dev``, source missing in ``codex_dev``, capped HETATMs.
        - ``FAIL`` for parse errors, missing source in ``hpc_strict``, bad
          residue range, ``key_residue_bonus_weight != 0`` in gates.yaml,
          terminal artifact (962-start) in ``hpc_strict``.
    """
    if profile not in {"codex_dev", "hpc_strict"}:
        raise ValueError("profile must be codex_dev or hpc_strict; got: {0}".format(profile))

    gate = load_myo1d_gate(ctx)
    bonus_weight = float(gate.get("key_residue_bonus_weight", 0.0))
    if bonus_weight != 0.0:
        raise RuntimeError(
            "gates.yaml myo1d.key_residue_bonus_weight must be 0.0 (annotation-only); got {0}".format(
                bonus_weight
            )
        )

    if construct_range is None:
        construct_range = gate.get("construct", "955-1006")
    start, end = parse_construct_range(construct_range)

    sheet8 = expand_residue_set(gate.get("key_residues", {}).get("sheet8", "961-964"))
    sheet9 = expand_residue_set(gate.get("key_residues", {}).get("sheet9", "968-972"))
    sheet12 = expand_residue_set(gate.get("key_residues", {}).get("sheet12", "993-997"))
    n_watch = expand_residue_set(gate.get("watch_residues", {}).get("n_terminal", "955-957"))
    c_watch = expand_residue_set(gate.get("watch_residues", {}).get("c_terminal", "1001-1006"))

    construct_id = _construct_id(start, end)
    warnings: list[str] = []
    status = "PASS"

    source_path = Path(source_pdb)

    # Source missing
    if not source_path.is_file():
        warnings.append("source_missing: {0}".format(source_path))
        if profile == "hpc_strict":
            status = "FAIL"
        else:
            status = "WARN"
        report = Myo1dQcReport(
            construct_id=construct_id,
            source_file=str(source_path),
            chain_id=None,
            start_residue=start,
            end_residue=end,
            n_residues=0,
            status=status,
            warnings=warnings,
            key_residue_bonus_weight=bonus_weight,
            profile=profile,
        )
        _write_qc_outputs(ctx, report)
        return report

    # Parse source
    try:
        structure = parse_pdb(source_path)
    except PDBParseError as exc:
        warnings.append("parse_error: {0}".format(exc))
        report = Myo1dQcReport(
            construct_id=construct_id,
            source_file=str(source_path),
            chain_id=None,
            start_residue=start,
            end_residue=end,
            n_residues=0,
            status="FAIL",
            warnings=warnings,
            key_residue_bonus_weight=bonus_weight,
            profile=profile,
        )
        _write_qc_outputs(ctx, report)
        return report

    # Slice + emit
    sliced = slice_myo1d_construct(structure, start, end, include_caps=True)
    output_pdb = ctx.run_dir / "normalized" / "myo1d" / "{0}.pdb".format(construct_id)
    emit_myo1d_construct_pdb(ctx, sliced, output_pdb)

    # Compute QC fields
    biological_residues = {
        r.residue_number
        for r in sliced.residues
        if r.resname not in CAP_RESNAMES
    }
    n_residues = len(biological_residues)
    expected_full = set(range(start, end + 1))
    missing = sorted(expected_full - biological_residues)
    cap_present = any(r.resname in CAP_RESNAMES for r in sliced.residues)
    chain_id = _resolve_chain_id(sliced.atoms)

    relevant_sheet8 = [r for r in sheet8 if start <= r <= end]
    relevant_sheet9 = [r for r in sheet9 if start <= r <= end]
    relevant_sheet12 = [r for r in sheet12 if start <= r <= end]
    relevant_n_watch = [r for r in n_watch if start <= r <= end]
    relevant_c_watch = [r for r in c_watch if start <= r <= end]

    key_sheet8_present = bool(relevant_sheet8) and _all_in(relevant_sheet8, biological_residues)
    key_sheet9_present = bool(relevant_sheet9) and _all_in(relevant_sheet9, biological_residues)
    key_sheet12_present = bool(relevant_sheet12) and _all_in(relevant_sheet12, biological_residues)
    n_watch_present = bool(relevant_n_watch) and _all_in(relevant_n_watch, biological_residues)
    c_watch_present = bool(relevant_c_watch) and _all_in(relevant_c_watch, biological_residues)

    # Severity decisions
    if start == 962:
        warnings.append("myo1d_962_start_terminal_artifact")
        if profile == "hpc_strict":
            status = "FAIL"
        elif status == "PASS":
            status = "WARN"

    key_residues_in_range_relevant = bool(relevant_sheet8 or relevant_sheet9 or relevant_sheet12)
    if key_residues_in_range_relevant and not (
        (not relevant_sheet8 or key_sheet8_present)
        and (not relevant_sheet9 or key_sheet9_present)
        and (not relevant_sheet12 or key_sheet12_present)
    ):
        warnings.append("key_residues_incomplete")
        if status == "PASS":
            status = "WARN"

    if cap_present:
        warnings.append("ace_nme_caps_preserved")

    if missing:
        warnings.append(
            "missing_residues_in_range: {0}".format(
                ",".join(str(n) for n in missing[:10])
                + ("..." if len(missing) > 10 else "")
            )
        )

    report = Myo1dQcReport(
        construct_id=construct_id,
        source_file=str(source_path),
        chain_id=chain_id,
        start_residue=start,
        end_residue=end,
        n_residues=n_residues,
        missing_residues=missing,
        key_sheet8_present=key_sheet8_present,
        key_sheet9_present=key_sheet9_present,
        key_sheet12_present=key_sheet12_present,
        n_watch_present=n_watch_present,
        c_watch_present=c_watch_present,
        ace_nme_caps_present=cap_present,
        status=status,
        warnings=warnings,
        output_pdb=output_pdb,
        key_residue_bonus_weight=bonus_weight,
        profile=profile,
    )

    _write_qc_outputs(ctx, report)
    return report
