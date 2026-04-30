"""Residue mapping restoration for M2.3 PPI contact tables.

PyRosetta pose numbering is not a source of truth for EGFR identity. This
module restores EGFR runtime residue numbers through the M1 receptor mapping
CSV before downstream pose QC or consensus work can interpret contacts.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

from egfr_myo1d.io.residue_mapping import read_residue_mapping


POSE_CONTACT_FIELDS = [
    "run_id",
    "state_id",
    "state_role",
    "method",
    "job_name",
    "seed",
    "pose_id",
    "score",
    "egfr_chain_id",
    "protomer_id",
    "egfr_runtime_residue",
    "egfr_original_residue",
    "egfr_uniprot_residue",
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

UNMAPPED_CONTACT_FIELDS = [
    "run_id",
    "state_id",
    "seed",
    "pose_id",
    "protomer_id",
    "egfr_runtime_residue",
    "reason",
    "mapping_csv",
]

MAPPING_QC_FIELDS = [
    "state_id",
    "mapping_csv",
    "input_contact_count",
    "restored_contact_count",
    "unmapped_contact_count",
    "unmapped_fraction",
    "status",
    "notes",
]


@dataclass(frozen=True)
class RuntimeMappingRecord:
    state_id: str
    protomer_id: str
    source_chain: str
    source_resseq: str
    source_icode: str
    source_resname: str
    runtime_chain: str
    runtime_resseq: int
    mapping_csv: str

    @property
    def original_residue(self) -> str:
        if self.source_icode:
            return "{0}{1}".format(self.source_resseq, self.source_icode)
        return self.source_resseq


@dataclass
class RuntimeResidueMap:
    mapping_csv: Path
    by_runtime: dict[int, list[RuntimeMappingRecord]] = field(default_factory=dict)
    by_protomer_runtime: dict[tuple[str, int], RuntimeMappingRecord] = field(default_factory=dict)

    def restore(self, protomer_id: str, runtime_residue: str) -> tuple[RuntimeMappingRecord | None, str]:
        try:
            runtime = int(str(runtime_residue).strip())
        except ValueError:
            return None, "invalid_runtime_residue"

        if protomer_id:
            by_pair = self.by_protomer_runtime.get((protomer_id, runtime))
            if by_pair is not None:
                return by_pair, ""

        candidates = self.by_runtime.get(runtime, [])
        if len(candidates) == 1:
            return candidates[0], ""
        if len(candidates) > 1:
            return None, "ambiguous_runtime_residue_requires_protomer"
        return None, "runtime_residue_not_in_m1_mapping"


@dataclass
class MappingRestorationResult:
    status: str
    restored_rows: list[dict[str, Any]]
    unmapped_rows: list[dict[str, Any]]
    qc_rows: list[dict[str, Any]]
    warnings: list[str] = field(default_factory=list)
    blockers: list[str] = field(default_factory=list)


def load_runtime_residue_map(mapping_csv: Path) -> RuntimeResidueMap:
    residue_map = RuntimeResidueMap(mapping_csv=Path(mapping_csv).resolve())
    for row in read_residue_mapping(mapping_csv):
        runtime_text = str(row.get("runtime_resseq", "")).strip()
        try:
            runtime = int(runtime_text)
        except ValueError:
            continue
        record = RuntimeMappingRecord(
            state_id=str(row.get("state", "")),
            protomer_id=str(row.get("protomer_id", "")),
            source_chain=str(row.get("source_chain", "")),
            source_resseq=str(row.get("source_resseq", "")),
            source_icode=str(row.get("source_icode", "")),
            source_resname=str(row.get("source_resname", "")),
            runtime_chain=str(row.get("runtime_chain", "")),
            runtime_resseq=runtime,
            mapping_csv=str(Path(mapping_csv).resolve()),
        )
        residue_map.by_runtime.setdefault(runtime, []).append(record)
        if record.protomer_id:
            residue_map.by_protomer_runtime[(record.protomer_id, runtime)] = record
    return residue_map


def _group_by_state(rows: Iterable[dict[str, Any]]) -> dict[str, list[dict[str, Any]]]:
    grouped: dict[str, list[dict[str, Any]]] = {}
    for row in rows:
        grouped.setdefault(str(row.get("state_id", "")), []).append(row)
    return grouped


def _raw_contact_to_output(
    run_id: str,
    row: dict[str, Any],
    mapping_record: RuntimeMappingRecord,
) -> dict[str, Any]:
    return {
        "run_id": str(row.get("run_id") or run_id),
        "state_id": str(row.get("state_id", "")),
        "state_role": str(row.get("state_role") or row.get("role") or ""),
        "method": str(row.get("method", "")),
        "job_name": str(row.get("job_name", "")),
        "seed": str(row.get("seed", "")),
        "pose_id": str(row.get("pose_id", "")),
        "score": str(row.get("score", "")),
        "egfr_chain_id": mapping_record.source_chain or mapping_record.runtime_chain,
        "protomer_id": mapping_record.protomer_id,
        "egfr_runtime_residue": str(mapping_record.runtime_resseq),
        "egfr_original_residue": mapping_record.original_residue,
        "egfr_uniprot_residue": mapping_record.original_residue,
        "egfr_resname": str(row.get("egfr_resname") or mapping_record.source_resname),
        "myo1d_residue": str(row.get("myo1d_residue", "")),
        "myo1d_resname": str(row.get("myo1d_resname", "")),
        "contact_type": str(row.get("contact_type", "")),
        "min_distance_angstrom": str(row.get("min_distance_angstrom", "")),
        "egfr_atom": str(row.get("egfr_atom", "")),
        "myo1d_atom": str(row.get("myo1d_atom", "")),
        "egfr_contact_centroid_x": str(row.get("egfr_contact_centroid_x", "")),
        "egfr_contact_centroid_y": str(row.get("egfr_contact_centroid_y", "")),
        "egfr_contact_centroid_z": str(row.get("egfr_contact_centroid_z", "")),
        "pose_acceptance_class": str(row.get("pose_acceptance_class", "")),
        "accepted_pose_flag": str(row.get("accepted_pose_flag", "")),
        "membrane_side_class": str(row.get("membrane_side_class", "")),
    }


def restore_contact_table_mapping(
    *,
    run_id: str,
    raw_contact_rows: list[dict[str, Any]],
    mapping_csv_by_state: dict[str, Path],
    unmapped_fraction_fail_threshold: float = 0.0,
) -> MappingRestorationResult:
    """Restore raw EGFR runtime contacts through M1 receptor mappings.

    ``raw_contact_rows`` are assumed to be already collected under the current
    run directory. The function does not parse pose PDB files, score contacts,
    or classify pose acceptance.
    """

    restored_rows: list[dict[str, Any]] = []
    unmapped_rows: list[dict[str, Any]] = []
    qc_rows: list[dict[str, Any]] = []
    warnings: list[str] = []
    blockers: list[str] = []
    mapping_cache: dict[str, RuntimeResidueMap] = {}

    grouped = _group_by_state(raw_contact_rows)
    if not raw_contact_rows:
        return MappingRestorationResult(
            status="WARN",
            restored_rows=[],
            unmapped_rows=[],
            qc_rows=[
                {
                    "state_id": "",
                    "mapping_csv": "",
                    "input_contact_count": 0,
                    "restored_contact_count": 0,
                    "unmapped_contact_count": 0,
                    "unmapped_fraction": "0.000000",
                    "status": "WARN",
                    "notes": "no_raw_contact_rows_to_restore",
                }
            ],
            warnings=["no_raw_contact_rows_to_restore"],
        )

    for state_id, rows in sorted(grouped.items()):
        mapping_csv = mapping_csv_by_state.get(state_id)
        state_restored = 0
        state_unmapped = 0
        if mapping_csv is None:
            for row in rows:
                unmapped_rows.append(
                    {
                        "run_id": str(row.get("run_id") or run_id),
                        "state_id": state_id,
                        "seed": str(row.get("seed", "")),
                        "pose_id": str(row.get("pose_id", "")),
                        "protomer_id": str(row.get("protomer_id", "")),
                        "egfr_runtime_residue": str(row.get("egfr_runtime_residue", "")),
                        "reason": "missing_mapping_csv_for_state",
                        "mapping_csv": "",
                    }
                )
            state_unmapped = len(rows)
        else:
            cache_key = str(Path(mapping_csv).resolve())
            if cache_key not in mapping_cache:
                mapping_cache[cache_key] = load_runtime_residue_map(mapping_csv)
            residue_map = mapping_cache[cache_key]
            for row in rows:
                record, reason = residue_map.restore(
                    str(row.get("protomer_id", "")),
                    str(row.get("egfr_runtime_residue", "")),
                )
                if record is None:
                    state_unmapped += 1
                    unmapped_rows.append(
                        {
                            "run_id": str(row.get("run_id") or run_id),
                            "state_id": state_id,
                            "seed": str(row.get("seed", "")),
                            "pose_id": str(row.get("pose_id", "")),
                            "protomer_id": str(row.get("protomer_id", "")),
                            "egfr_runtime_residue": str(row.get("egfr_runtime_residue", "")),
                            "reason": reason,
                            "mapping_csv": str(Path(mapping_csv).resolve()),
                        }
                    )
                else:
                    state_restored += 1
                    restored_rows.append(_raw_contact_to_output(run_id, row, record))

        total = len(rows)
        fraction = float(state_unmapped) / float(total) if total else 0.0
        if fraction > unmapped_fraction_fail_threshold:
            status = "FAIL"
            blockers.append(
                "unmapped_fraction_exceeds_threshold:{0}:{1:.6f}".format(
                    state_id, fraction
                )
            )
        elif state_unmapped:
            status = "WARN"
            warnings.append("unmapped_contacts_present:{0}:{1}".format(state_id, state_unmapped))
        else:
            status = "PASS"
        qc_rows.append(
            {
                "state_id": state_id,
                "mapping_csv": str(Path(mapping_csv).resolve()) if mapping_csv else "",
                "input_contact_count": total,
                "restored_contact_count": state_restored,
                "unmapped_contact_count": state_unmapped,
                "unmapped_fraction": "{0:.6f}".format(fraction),
                "status": status,
                "notes": "egfr_uniprot_residue_uses_m1_source_resseq_until_uniprot_column_exists",
            }
        )

    if blockers:
        status = "FAIL"
    elif warnings:
        status = "WARN"
    else:
        status = "PASS"
    return MappingRestorationResult(
        status=status,
        restored_rows=restored_rows,
        unmapped_rows=unmapped_rows,
        qc_rows=qc_rows,
        warnings=warnings,
        blockers=blockers,
    )
