"""Shared residue normalization utilities."""
import re
from typing import Optional, Set

_RESNAME_NORMALIZE = {
    "HSD": "HIS", "HSE": "HIS", "HSP": "HIS",
    "CYX": "CYS",
    "HIE": "HIS", "HID": "HIS", "HIP": "HIS",
}


def normalize_residue_id(residue_id: str) -> str:
    """Normalize 'A:MET971' -> 'MET971' (strip chain, normalize resname)."""
    if ":" in residue_id:
        residue_id = residue_id.split(":", 1)[1]
    i = 0
    while i < len(residue_id) and not residue_id[i].isdigit() and residue_id[i] != '-':
        i += 1
    resname = residue_id[:i]
    resnum = residue_id[i:]
    resname = _RESNAME_NORMALIZE.get(resname, resname)
    return f"{resname}{resnum}"


def extract_resnum(normalized_id: str) -> Optional[int]:
    """Extract residue number from normalized id like 'MET971' -> 971."""
    m = re.search(r'(-?\d+)$', normalized_id)
    return int(m.group(1)) if m else None


def parse_residue_set(raw: str) -> Set[str]:
    """Parse semicolon-separated residue string into a normalized set."""
    if not raw:
        return set()
    return {normalize_residue_id(r) for r in raw.split(";") if r.strip()}
