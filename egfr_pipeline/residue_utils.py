"""Shared residue normalization utilities."""
import re
from typing import Optional, Set

_RESNAME_NORMALIZE = {
    "HSD": "HIS", "HSE": "HIS", "HSP": "HIS",
    "CYX": "CYS",
    "HIE": "HIS", "HID": "HIS", "HIP": "HIS",
}


def normalize_residue_id(residue_id: str, keep_chain: bool = False) -> str:
    """Normalize residue identifier with optional chain preservation.

    Args:
        residue_id: Raw residue id, e.g. ``'A:MET971'`` or ``'MET971'``.
        keep_chain: When *True*, preserve the chain prefix (``'A:MET971'``
            with normalized resname).  When *False* (default), strip the
            chain to produce ``'MET971'``.
    """
    chain_prefix = ""
    if ":" in residue_id:
        chain, residue_id = residue_id.split(":", 1)
        if keep_chain:
            chain_prefix = f"{chain}:"
    i = 0
    while i < len(residue_id) and not residue_id[i].isdigit() and residue_id[i] != '-':
        i += 1
    resname = residue_id[:i]
    resnum = residue_id[i:]
    resname = _RESNAME_NORMALIZE.get(resname, resname)
    return f"{chain_prefix}{resname}{resnum}"


def extract_resnum(normalized_id: str) -> Optional[int]:
    """Extract residue number from normalized id like 'MET971' -> 971."""
    m = re.search(r'(-?\d+)$', normalized_id)
    return int(m.group(1)) if m else None


def parse_residue_set(raw: str, keep_chain: bool = False) -> Set[str]:
    """Parse semicolon-separated residue string into a normalized set.

    Args:
        raw: Semicolon-separated residue identifiers.
        keep_chain: Passed through to :func:`normalize_residue_id`.
    """
    if not raw:
        return set()
    return {normalize_residue_id(r, keep_chain=keep_chain)
            for r in raw.split(";") if r.strip()}
