"""MYO1D beta-meander annotation validation for Task 3."""

from __future__ import annotations

from pathlib import Path

from egfr_myo1d.structure.pdb_parser import PDBStructure, parse_pdb


DEFAULT_ANNOTATION = {
    "primary_construct": "MYO1D_sheet8_9_12_core_955_1001",
    "primary_range": "955-1001",
    "comparator_construct": "MYO1D_ext_beta_meander_955_1006_tail_masked",
    "comparator_range": "955-1006",
    "structural_buffer": "955-960",
    "primary_active_face": "961-964,968-972",
    "support_region": "993-997",
    "short_c_terminal_cap": "998-1001",
    "extended_tail_noise_zone": "998-1006",
    "key_residue_bonus_weight": 0.0,
    "key_residue_contact_is_annotation_only": True,
}


def expand_range(text: str) -> set[int]:
    start_text, end_text = text.split("-", 1)
    start = int(start_text)
    end = int(end_text)
    return set(range(start, end + 1))


def expand_multi_range(text: str) -> set[int]:
    residues: set[int] = set()
    for chunk in text.split(","):
        residues.update(expand_range(chunk.strip()))
    return residues


def classify_construct(residue_numbers: set[int]) -> str:
    if not residue_numbers:
        return "unknown"
    min_residue = min(residue_numbers)
    max_residue = max(residue_numbers)
    if min_residue == 955 and max_residue == 1001:
        return "955-1001"
    if min_residue == 955 and max_residue == 1006:
        return "955-1006"
    return "unknown"


def _severity_for_missing(mode: str, profile: str) -> str:
    if mode == "smoke_env" and profile == "codex_dev":
        return "WARN"
    return "FAIL"


def validate_myo1d_structure(
    structure: PDBStructure,
    source_path: Path,
    mode: str,
    profile: str,
) -> dict[str, object]:
    residue_numbers = {residue.residue_number for residue in structure.residues}
    messages: list[dict[str, object]] = []
    active_face = expand_multi_range(DEFAULT_ANNOTATION["primary_active_face"])
    support_region = expand_range(DEFAULT_ANNOTATION["support_region"])
    tail_noise_zone = expand_range("1002-1006")

    missing_active = sorted(active_face - residue_numbers)
    missing_support = sorted(support_region - residue_numbers)
    if missing_active:
        messages.append(
            {
                "severity": _severity_for_missing(mode, profile),
                "code": "missing_active_face_residues",
                "message": f"missing MYO1D active-face residues: {missing_active}",
            }
        )
    if missing_support:
        messages.append(
            {
                "severity": _severity_for_missing(mode, profile),
                "code": "missing_sheet12_support_residues",
                "message": f"missing MYO1D sheet-12 support residues: {missing_support}",
            }
        )

    if residue_numbers:
        first_residue = min(residue_numbers)
        if first_residue == 962:
            messages.append(
                {
                    "severity": _severity_for_missing(mode, profile),
                    "code": "terminal_artifact_start_962",
                    "message": "main MYO1D partner starts at 962; expected 955 buffer is absent",
                }
            )
        elif first_residue > 955:
            messages.append(
                {
                    "severity": "WARN",
                    "code": "missing_n_terminal_buffer",
                    "message": f"MYO1D first residue is {first_residue}; expected 955 buffer",
                }
            )

    tail_present = sorted(tail_noise_zone & residue_numbers)
    if tail_present:
        messages.append(
            {
                "severity": "WARN",
                "code": "tail_noise_zone_present",
                "message": f"residues {tail_present} are tail/noise-zone residues for downstream accounting",
            }
        )

    if residue_numbers:
        terminal_residues = {number for number in residue_numbers if number >= 998}
        if len(terminal_residues) >= max(3, len(residue_numbers) // 3):
            messages.append(
                {
                    "severity": "WARN",
                    "code": "terminal_records_dominate",
                    "message": "terminal 998+ residues are prominent; inspect for terminal artifact effects",
                }
            )

    status = "PASS"
    if any(message["severity"] == "FAIL" for message in messages):
        status = "FAIL"
    elif messages:
        status = "WARN"

    return {
        "input_path": str(source_path),
        "status": status,
        "construct_class": classify_construct(residue_numbers),
        "residue_count": len(residue_numbers),
        "min_residue": min(residue_numbers) if residue_numbers else None,
        "max_residue": max(residue_numbers) if residue_numbers else None,
        "active_face_present": sorted(active_face & residue_numbers),
        "support_region_present": sorted(support_region & residue_numbers),
        "tail_noise_zone_present": tail_present,
        "annotation": DEFAULT_ANNOTATION,
        "key_residue_bonus_weight": 0.0,
        "key_residue_contact_is_annotation_only": True,
        "messages": messages,
    }


def validate_myo1d_file(path: Path, mode: str, profile: str) -> dict[str, object]:
    return validate_myo1d_structure(parse_pdb(path), path, mode, profile)

